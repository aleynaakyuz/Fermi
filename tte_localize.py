"""
Localize bursts from TTE instead of trigdat, and see whether it helps.

The hypothesis being tested: trigdat's 64 ms grid over-constrains the best-S/N
window search for short bursts. Their median window is 4 bins wide, so the
search picks among roughly 8 candidate intervals. TTE can be binned as finely
as wanted, giving the same search hundreds of candidates.

Everything else is held fixed -- same background windows, same 50-300 keV
band, same legacy_DoL with scattering on -- so any difference is attributable
to the time binning alone.

Two details the TTE path has to get right:

  energy    legacy_DoL's response tables are built for the 8 CTIME channel
            boundaries. TTE has 128 channels, so it is rebinned onto those
            edges before the fit. They line up to better than a keV at the
            low end and a few keV at 300.

  deadtime  TTE carries per-event deadtime, which gdt exposes through
            get_exposure(). Rates are counts/exposure, matching what the fit
            expects, rather than counts/bin-width.

Usage:
    python tte_localize.py --dt 0.008
    python tte_localize.py --triggers bn180402406 --dt 0.004
"""
import argparse
import csv
import glob
import json
import os
import warnings

warnings.filterwarnings("ignore")

import numpy as np
from astropy.coordinates import SkyCoord

import analyze_grbs as A
from gdt.missions.fermi.gbm.trigdat import Trigdat
from gdt.missions.fermi.gbm.tte import GbmTte
from gdt.missions.fermi.time import Time
from gdt.core.binning.unbinned import bin_by_time
from gdt.core.background.fitter import BackgroundFitter
from gdt.core.background.binned import Polynomial
from gdt.missions.fermi.gbm.localization.dol.legacy_dol import legacy_DoL

TTEDIR = "data/tte"
DATADIR = "data/gbm"
MANIFEST = "catalog/fermigbrst.csv"
OUT = "site/tte_results.json"

# The 8 CTIME/trigdat channel edges the DoL response tables assume.
CTIME_EDGES = np.array([3.4, 10.0, 22.0, 44.0, 95.0, 300.0, 500.0, 800.0, 2000.0])
NAI = [f"n{d}" for d in "0123456789ab"]


def load_events(trigger, det, grid, chan_edges):
    """Raw events for one detector, with each event's CTIME channel index.

    Counting events directly is exact. Binning first and then summing whole
    bins over a window is not: any bin overlapping the window edge is taken
    whole, which on a 0.064 s window at 8 ms resolution adds up to two bins
    and 32% spurious counts.

    `chan_edges` are THAT detector's own trigdat channel boundaries. Detector
    gains differ, so a shared set of edges biases each detector differently,
    and the fit reads precisely that relative pattern across detectors.
    """
    hits = glob.glob(os.path.join(TTEDIR, trigger, f"glg_tte_{det}_*.fit"))
    if not hits:
        return None
    tte = GbmTte.open(sorted(hits)[-1])
    tte = tte.slice_time((grid[0], grid[1]))
    data = tte.data
    # data.times is already relative to the trigger -- gdt subtracts trigtime
    # on read, so subtracting it again sends every event to -5e8 s.
    times = np.asarray(data.times, dtype=float)
    emin = np.asarray(tte.ebounds.low_edges(), dtype=float)
    emax = np.asarray(tte.ebounds.high_edges(), dtype=float)
    centres = 0.5 * (emin + emax)
    chan_of_pha = np.clip(np.digitize(centres, chan_edges) - 1, 0, 7)
    pha = np.asarray(data.channels, dtype=int)
    chan = chan_of_pha[np.clip(pha, 0, len(chan_of_pha) - 1)]
    deadtime = float(tte.event_deadtime), float(tte.overflow_deadtime)
    tte.close()
    order = np.argsort(times)
    return times[order], chan[order], deadtime


def count_exact(times, chan, window, deadtime):
    """Counts per CTIME channel in [t0, t1), and the live exposure."""
    t0, t1 = window
    i0, i1 = np.searchsorted(times, [t0, t1])
    sel = chan[i0:i1]
    counts = np.bincount(sel, minlength=8)[:8].astype(float)
    # Deadtime: every event costs event_deadtime, overflow events cost more.
    n_over = int((sel == 7).sum())
    dead = (len(sel) - n_over) * deadtime[0] + n_over * deadtime[1]
    return counts, max((t1 - t0) - dead, 1e-6)


def load_tte_phaii(trigger, det, dt, grid):
    """TTE for one detector, binned at dt and rebinned onto the CTIME edges.

    `grid` is a shared (tstart, tstop) so every detector lands on identical
    bin edges -- detectors cover slightly different time ranges, and binning
    each on its own span gives arrays that differ by a bin and cannot be
    summed. Restricting to the span actually used also avoids binning the full
    610 s at millisecond resolution when only ~160 s is needed.
    """
    hits = glob.glob(os.path.join(TTEDIR, trigger, f"glg_tte_{det}_*.fit"))
    if not hits:
        return None
    tte = GbmTte.open(sorted(hits)[-1])
    tte = tte.slice_time((grid[0], grid[1]))
    phaii = tte.to_phaii(bin_by_time, dt, tstart=grid[0], tstop=grid[1])

    # Map the 128 TTE channels onto the 8 the response tables use.
    emin = np.asarray(phaii.data.emin, dtype=float)
    emax = np.asarray(phaii.data.emax, dtype=float)
    centres = 0.5 * (emin + emax)
    idx = np.clip(np.digitize(centres, CTIME_EDGES) - 1, 0, 7)

    counts = np.asarray(phaii.data.counts, dtype=float)
    rebinned = np.zeros((counts.shape[0], 8))
    for ch in range(8):
        sel = idx == ch
        if sel.any():
            rebinned[:, ch] = counts[:, sel].sum(axis=1)

    tstart = np.asarray(phaii.data.tstart, dtype=float)
    tstop = np.asarray(phaii.data.tstop, dtype=float)
    exposure = np.asarray(phaii.data.exposure, dtype=float)
    tte.close()
    return tstart, tstop, exposure, rebinned


def integrate(tstart, tstop, exposure, counts, window):
    sel = (tstop > window[0]) & (tstart < window[1])
    if not sel.any():
        return np.zeros(8), 0.0
    return counts[sel].sum(axis=0), float(exposure[sel].sum())


def fit_background(tstart, tstop, counts, exposure, bg_times, order=1):
    """Per-channel linear background, fitted on the background windows only."""
    mid = 0.5 * (tstart + tstop)
    inbg = np.zeros(len(mid), dtype=bool)
    for lo, hi in bg_times:
        inbg |= (mid >= lo) & (mid <= hi)
    if inbg.sum() < 4:
        raise ValueError("too few bins in the background windows")
    exp = np.clip(exposure, 1e-9, None)
    model = np.zeros_like(counts)
    for ch in range(counts.shape[1]):
        rate = counts[inbg, ch] / exp[inbg]
        coef = np.polyfit(mid[inbg], rate, order)
        model[:, ch] = np.clip(np.polyval(coef, mid), 0.0, None) * exp
    return model


def best_sn_window(tstart, tstop, counts, bkg, search_span, min_dur):
    """Same criterion as the trigdat path, on a much finer grid.

    Candidate intervals are capped so the scan stays O(n * k) rather than
    O(n^2): at 8 ms over a 2 s search span that is already 250 bins.
    """
    mid = 0.5 * (tstart + tstop)
    inside = np.where((mid >= search_span[0]) & (mid <= search_span[1]))[0]
    if inside.size < 3:
        return None
    lo, hi = int(inside[0]), int(inside[-1]) + 1
    # Channels 3 and 4 are the 50-300 keV band the fit uses.
    net = (counts[:, 3] + counts[:, 4]) - (bkg[:, 3] + bkg[:, 4])
    bg = np.clip(bkg[:, 3] + bkg[:, 4], 0.0, None)
    cn = np.concatenate([[0.0], np.cumsum(net)])
    cb = np.concatenate([[0.0], np.cumsum(bg)])
    dt = float(np.median(tstop - tstart))
    min_bins = max(1, int(round(min_dur / dt)))

    best, best_snr = None, -np.inf
    for i in range(lo, hi):
        for j in range(i + min_bins, hi + 1):
            s = cn[j] - cn[i]
            n = cb[j] - cb[i]
            if s <= 0 or n <= 0:
                continue
            snr = s / np.sqrt(n)
            if snr > best_snr:
                best_snr, best = snr, (i, j)
    if best is None:
        return None
    i, j = best
    return float(tstart[i]), float(tstop[j - 1]), float(best_snr)


def run(trigger, row, dt):
    """Localize one burst from TTE. Returns a comparison record."""
    # Geometry, window definitions and the per-detector channel edges all come
    # from trigdat, which also supplies the position history TTE lacks.
    td = Trigdat.open(glob.glob(f"{DATADIR}/{trigger}/glg_trigdat_*")[0])
    t90 = float(row["t90"])
    t90_start = float(row["t90_start"])
    probe = td.to_phaii(td._detectors[0], timescale=64)
    span = (float(probe.data.tstart[0]), float(probe.data.tstop[-1]))
    _, t90_win, bg_times = A.choose_windows(t90, t90_start, span)

    lo = min(w[0] for w in bg_times) - 1.0
    hi = max(w[1] for w in bg_times) + 1.0
    grid = (lo, hi)

    # Each detector gets ITS OWN trigdat channel edges, because gains differ.
    edges, events = {}, {}
    for det in NAI:
        ph = td.to_phaii(det, timescale=1024)
        e = np.concatenate([np.asarray(ph.data.emin, dtype=float),
                            [float(ph.data.emax[-1])]])
        edges[det] = e
        got = load_events(trigger, det, grid, e)
        if got is None:
            raise FileNotFoundError(f"missing TTE for {det}")
        events[det] = got

    # --- window search, on the summed triggered detectors -------------------
    trig_dets = [d for d in td.triggered_detectors if d in events] or NAI[:3]
    bin_edges = np.arange(grid[0], grid[1] + dt, dt)
    centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    summed = np.zeros((len(centres), 8))
    for det in trig_dets:
        times, chan, _ = events[det]
        for ch in (3, 4):
            sel = times[chan == ch]
            summed[:, ch] += np.histogram(sel, bins=bin_edges)[0]

    inbg = np.zeros(len(centres), dtype=bool)
    for a, b in bg_times:
        inbg |= (centres >= a) & (centres <= b)
    bkg_sum = np.zeros_like(summed)
    for ch in (3, 4):
        coef = np.polyfit(centres[inbg], summed[inbg, ch] / dt, 1)
        bkg_sum[:, ch] = np.clip(np.polyval(coef, centres), 0.0, None) * dt

    margin = 0.15 * max(t90_win[1] - t90_win[0], 1.0)
    found = best_sn_window(bin_edges[:-1], bin_edges[1:], summed, bkg_sum,
                           (t90_win[0] - margin, t90_win[1] + margin),
                           min_dur=dt)
    sn_win = (found[0], found[1]) if found else t90_win

    # --- counts for the fit: exact over the window, per detector ------------
    mrates = np.zeros((14, 8))
    brates = np.zeros((14, 8))
    exposures = []
    for k, det in enumerate(NAI):
        times, chan, dead = events[det]
        src, exp = count_exact(times, chan, sn_win, dead)
        mrates[k] = src
        exposures.append(exp)

        # Background rate per channel from the background windows, then scaled
        # to the source window -- same construction as the trigdat path.
        bb = np.zeros(8)
        for ch in range(8):
            hist = np.histogram(times[chan == ch], bins=bin_edges)[0]
            if inbg.sum() < 4:
                raise ValueError("too few bins in the background windows")
            coef = np.polyfit(centres[inbg], hist[inbg] / dt, 1)
            mid = 0.5 * (sn_win[0] + sn_win[1])
            bb[ch] = max(float(np.polyval(coef, mid)), 0.0) * (sn_win[1] - sn_win[0])
        brates[k] = bb
    sdur = float(np.mean(exposures))

    crange = [int(np.digitize(e, edges[NAI[0]], right=True) - 1) for e in (50.0, 300.0)]
    tc = td.trigtime + 0.5 * sum(sn_win)
    frame = td.poshist.at(Time(tc, format="fermi"))
    loc = legacy_DoL().eval(
        crange,
        np.rint(mrates / sdur).astype(np.int32),
        np.rint(brates / sdur).astype(np.int32),
        sdur, sdur,
        frame.obsgeoloc.xyz.to_value("km"), frame.quaternion.scalar_last,
        edges[NAI[0]], td.headers["PRIMARY"]["RA_OBJ"],
        td.headers["PRIMARY"]["DEC_OBJ"], int(tc), scat_opt=1)
    ra, dec = np.degrees([loc["best"]["ra"], loc["best"]["dec"]])

    prev = json.load(open(f"site/all/results/{trigger}.json"))
    here = SkyCoord(ra, dec, unit="deg")
    official = SkyCoord(prev["official_ra"], prev["official_dec"], unit="deg")
    rec = {
        "trigger": trigger, "dt": dt, "t90": t90,
        "tte_ra": float(ra), "tte_dec": float(dec),
        "tte_window": [round(sn_win[0], 4), round(sn_win[1], 4)],
        "tte_snr": round(found[2], 2) if found else None,
        "trigdat_window": [round(prev["peak_mode"]["src_duration"], 3)],
        "sep_tte_official": float(here.separation(official).deg),
        "sep_trigdat_official": prev["peak_mode"]["sep"],
    }
    if "external" in prev:
        truth = SkyCoord(prev["external"]["ra"], prev["external"]["dec"], unit="deg")
        rec["sep_tte_truth"] = float(here.separation(truth).deg)
        rec["sep_trigdat_truth"] = prev["external"]["sep_peak"]
        rec["sep_official_truth"] = prev["external"]["sep_official"]
    return rec


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dt", type=float, default=0.008,
                        help="TTE bin width in seconds")
    parser.add_argument("--triggers", type=str, default="")
    args = parser.parse_args()

    rows = {r["trigger_name"].strip(): r
            for r in csv.DictReader(open(MANIFEST, newline=""))}
    if args.triggers:
        triggers = [t.strip() for t in args.triggers.split(",")]
    else:
        triggers = sorted(d for d in os.listdir(TTEDIR)
                          if os.path.isdir(os.path.join(TTEDIR, d)))

    print(f"{len(triggers)} bursts with TTE, binning at {args.dt*1000:.0f} ms\n")
    out, fails = [], []
    for trigger in triggers:
        try:
            rec = run(trigger, rows[trigger], args.dt)
            out.append(rec)
            line = (f"  {trigger}  win {rec['tte_window'][1]-rec['tte_window'][0]:6.3f}s "
                    f"(trigdat {rec['trigdat_window'][0]:5.3f}s)  "
                    f"vs official: TTE {rec['sep_tte_official']:6.2f} "
                    f"trigdat {rec['sep_trigdat_official']:6.2f}")
            if "sep_tte_truth" in rec:
                line += (f"  | vs truth: TTE {rec['sep_tte_truth']:6.2f} "
                         f"trigdat {rec['sep_trigdat_truth']:6.2f} "
                         f"official {rec['sep_official_truth']:6.2f}")
            print(line)
        except Exception as err:
            fails.append((trigger, f"{type(err).__name__}: {err}"))
            print(f"  {trigger}  FAILED {type(err).__name__}: {err}")

    with open(OUT, "w") as handle:
        json.dump({"dt": args.dt, "results": out, "failures": fails}, handle, indent=2)
    print(f"\n{len(out)} ok, {len(fails)} failed -> {OUT}")

    if out:
        tt = np.array([r["sep_tte_official"] for r in out])
        tg = np.array([r["sep_trigdat_official"] for r in out])
        print(f"\nvs official map: TTE median {np.median(tt):.2f} deg, "
              f"trigdat {np.median(tg):.2f} deg, TTE closer on {(tt<tg).sum()}/{len(tt)}")
        wt = [r for r in out if "sep_tte_truth" in r]
        if wt:
            a = np.array([r["sep_tte_truth"] for r in wt])
            b = np.array([r["sep_trigdat_truth"] for r in wt])
            c = np.array([r["sep_official_truth"] for r in wt])
            print(f"vs truth ({len(wt)} bursts): TTE {np.median(a):.2f}, "
                  f"trigdat {np.median(b):.2f}, official {np.median(c):.2f} deg")
            print(f"   TTE closer to truth than trigdat on {(a<b).sum()}/{len(a)}, "
                  f"than official on {(a<c).sum()}/{len(a)}")


if __name__ == "__main__":
    main()
