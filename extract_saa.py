"""
Fold three weeks of continuous GBM data on SAA exit, for plotting.

Reduces 5 GB of daily CTIME to a few hundred kilobytes of curves:

  fold    per detector, median count rate against minutes since Fermi last
          left the South Atlantic Anomaly, 0 to 90 minutes in half-minute
          bins, stacked over every exit in the window
  chan    the same fold split by CTIME's eight energy channels for one NaI,
          which shows whether the activation is harder or softer than the
          background it sits on
  fit     A exp(-t/tau) + B fitted per detector over the first hour
  day     one 24-hour lightcurve at 10 s resolution with the SAA passages
          marked, so the sawtooth is visible directly rather than only in
          the stacked form

Bins with no exposure are dropped: the detectors are switched off inside the
SAA, so those are real gaps rather than zeros.

Usage:
    python extract_saa.py --workers 7
"""
import argparse
import glob
import json
import os
import warnings

warnings.filterwarnings("ignore")

import numpy as np

DAILY = "data/daily"
OUT = "site/saa_profiles.json"

DETECTORS = [f"n{d}" for d in "0123456789ab"] + ["b0", "b1"]
CHAN_DET = "n0"
EXAMPLE_DAY = "220827"

FOLD_MAX = 90.0          # minutes after an exit
FOLD_STEP = 0.5
BASE_LO, BASE_HI = 60.0, 90.0     # quiet reference window


def saa_exits(day):
    """Exit times, and the (entry, exit) spans, for one day."""
    from astropy.io import fits
    paths = glob.glob(os.path.join(DAILY, day, "glg_poshist_all_*.fit"))
    if not paths:
        return None, None
    with fits.open(paths[0]) as handle:
        data = handle[1].data
        met = np.asarray(data["SCLK_UTC"], dtype=float)
        saa = ((np.asarray(data["FLAGS"]).astype(np.uint16) >> 1) & 1).astype(bool)
    edge = np.diff(saa.astype(np.int8))
    exits = met[1:][edge == -1]
    entries = met[1:][edge == 1]
    spans = []
    for e in entries:
        after = exits[exits > e]
        if after.size:
            spans.append((float(e), float(after[0])))
    return exits, spans


def read_ctime(day, det):
    from astropy.io import fits
    paths = glob.glob(os.path.join(DAILY, day, f"glg_ctime_{det}_*.pha"))
    if not paths:
        return None
    with fits.open(paths[0]) as handle:
        spec = handle["SPECTRUM"].data
        t0 = np.asarray(spec["TIME"], dtype=float)
        t1 = np.asarray(spec["ENDTIME"], dtype=float)
        counts = np.asarray(spec["COUNTS"], dtype=float)
        expo = np.asarray(spec["EXPOSURE"], dtype=float)
        eb = handle["EBOUNDS"].data
        edges = np.append(np.asarray(eb["E_MIN"], float),
                          float(eb["E_MAX"][-1]))
    good = expo > 0.1
    return (0.5 * (t0 + t1))[good], counts[good], expo[good], edges


#   Per fold bin, keep at most this many samples before taking the median.
#   About 22,000 land in each bin over 21 days, so this barely subsamples.
MAX_SAMPLES = 20000
MAX_SAMPLES_CHAN = 5000


def fold_one(args):
    """Folded MEDIAN rate for one detector across every day.

    Medians, not an exposure-weighted mean, and the difference is not
    cosmetic. The mean is dominated by outliers -- GRBs, solar flares,
    particle events -- and the 60-90 minute reference window has by far the
    most exposure, so it catches the most of them. Taking means put that
    'quiet' baseline at 1405 counts/s on n0 while the profile had already
    decayed to 1200 by six minutes, i.e. the reference was higher than the
    thing it was supposed to be a floor for, which dragged the fitted decay
    constant from 5.6 minutes down to 0.6. The median is unmoved by a handful
    of bright minutes in three weeks.
    """
    det, days = args
    nb = int(round(FOLD_MAX / FOLD_STEP))
    buckets = [[] for _ in range(nb)]
    chan_buckets = None
    ebounds = None
    want_chan = det == CHAN_DET

    for day in days:
        exits, _ = saa_exits(day)
        if exits is None or not exits.size:
            continue
        got = read_ctime(day, det)
        if got is None:
            continue
        mid, counts, expo, eb = got
        ebounds = eb if ebounds is None else ebounds
        if want_chan and chan_buckets is None:
            chan_buckets = [[[] for _ in range(counts.shape[1])]
                            for _ in range(nb)]
        idx = np.searchsorted(exits, mid) - 1
        ok = idx >= 0
        since = np.full(mid.shape, np.nan)
        since[ok] = (mid[ok] - exits[idx[ok]]) / 60.0
        sel = ok & (since >= 0) & (since < FOLD_MAX)
        if not sel.any():
            continue
        b = np.clip((since[sel] / FOLD_STEP).astype(int), 0, nb - 1)
        rate = counts[sel].sum(axis=1) / np.maximum(expo[sel], 1e-9)
        order = np.argsort(b, kind="stable")
        b_s, rate_s = b[order], rate[order]
        starts = np.searchsorted(b_s, np.arange(nb), side="left")
        stops = np.searchsorted(b_s, np.arange(nb), side="right")
        for i in range(nb):
            if stops[i] > starts[i]:
                buckets[i].append(rate_s[starts[i]:stops[i]])
        if want_chan:
            crate = counts[sel] / np.maximum(expo[sel], 1e-9)[:, None]
            crate_s = crate[order]
            for i in range(nb):
                if stops[i] > starts[i]:
                    chunk = crate_s[starts[i]:stops[i]]
                    for c in range(chunk.shape[1]):
                        chan_buckets[i][c].append(chunk[:, c])

    def collapse(parts, cap):
        if not parts:
            return float("nan")
        v = np.concatenate(parts)
        if v.size > cap:
            v = v[:: max(1, v.size // cap)]
        return float(np.median(v))

    rate = np.array([collapse(buckets[i], MAX_SAMPLES) for i in range(nb)])
    out = {"rate": rate.tolist(),
           "ebounds": (ebounds.tolist() if ebounds is not None else None)}
    if want_chan and chan_buckets is not None:
        nchan = len(chan_buckets[0])
        out["chan"] = [[collapse(chan_buckets[i][c], MAX_SAMPLES_CHAN)
                        for i in range(nb)] for c in range(nchan)]
    return det, out


def fit_decay(t, rate):
    """A exp(-t/tau) + B over the first hour."""
    from scipy.optimize import curve_fit
    model = lambda x, A, tau, B: A * np.exp(-x / tau) + B
    k = np.isfinite(rate) & (t < 60.0)
    if k.sum() < 10:
        return None
    try:
        p, _ = curve_fit(model, t[k], rate[k],
                         p0=[rate[k][0] - rate[k][-1], 5.0, rate[k][-1]],
                         maxfev=40000)
    except Exception:
        return None
    if not np.all(np.isfinite(p)) or p[1] <= 0:
        return None
    return {"A": round(float(p[0]), 1), "tau": round(float(p[1]), 2),
            "B": round(float(p[2]), 1)}


def example_day(day, det, step=10.0):
    """A 24 h lightcurve at coarse resolution, plus the SAA spans."""
    got = read_ctime(day, det)
    exits, spans = saa_exits(day)
    if got is None or exits is None:
        return None
    mid, counts, expo, _ = got
    t0 = mid.min()
    rel = mid - t0
    nb = int(np.ceil(rel.max() / step))
    b = np.clip((rel / step).astype(int), 0, nb - 1)
    c = np.bincount(b, weights=counts.sum(axis=1), minlength=nb)
    e = np.bincount(b, weights=expo, minlength=nb)
    rate = np.where(e > 1.0, c / np.maximum(e, 1e-9), np.nan)
    hours = (np.arange(nb) * step + step / 2) / 3600.0
    keep = np.isfinite(rate)
    return {
        "det": det, "day": day, "step": step,
        "t": [round(float(v), 4) for v in hours[keep]],
        "r": [int(round(float(v))) for v in rate[keep]],
        "saa": [[round((a - t0) / 3600.0, 4), round((b2 - t0) / 3600.0, 4)]
                for a, b2 in spans],
        "exits": [round(float(x - t0) / 3600.0, 4) for x in exits],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=7)
    args = parser.parse_args()

    days = sorted(d for d in os.listdir(DAILY)
                  if os.path.isdir(os.path.join(DAILY, d)))
    if not days:
        raise SystemExit("no daily data; run download_daily.py first")
    n_exits = 0
    for day in days:
        exits, _ = saa_exits(day)
        n_exits += 0 if exits is None else int(exits.size)
    print(f"{len(days)} days, {n_exits} SAA exits, {args.workers} workers",
          flush=True)

    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(args.workers) as pool:
        folded = dict(pool.map(fold_one, [(d, days) for d in DETECTORS]))

    centres = np.arange(0.0, FOLD_MAX, FOLD_STEP) + FOLD_STEP / 2
    out = {
        "days": len(days), "exits": n_exits,
        "first": days[0], "last": days[-1],
        "t": [round(float(v), 2) for v in centres],
        "fold": {}, "fit": {}, "chan": None, "day": None,
    }
    quiet = (centres >= BASE_LO) & (centres < BASE_HI)
    for det, got in folded.items():
        if got is None:
            continue
        rate = np.asarray(got["rate"], dtype=float)
        out["fold"][det] = [None if not np.isfinite(v) else round(float(v), 1)
                            for v in rate]
        fit = fit_decay(centres, rate)
        if fit:
            base = np.nanmedian(rate[quiet])
            fit["peak_excess"] = round(float(100 * (rate[0] / base - 1)), 1)
            fit["base_quiet"] = round(float(base), 1)
            out["fit"][det] = fit

        if det == CHAN_DET and got.get("chan"):
            eb = got["ebounds"]
            rows = []
            for c, series in enumerate(got["chan"]):
                r = np.asarray(series, dtype=float)
                base = np.nanmedian(r[quiet])
                rows.append({
                    "lo": round(float(eb[c]), 1), "hi": round(float(eb[c + 1]), 1),
                    "base": round(float(base), 2),
                    # Relative to each channel's own quiet level, because the
                    # bottom channel outcounts the top one by orders of
                    # magnitude and would otherwise be the only one visible.
                    "rel": [None if not np.isfinite(v) or base <= 0
                            else round(float(v / base), 4) for v in r],
                })
            out["chan"] = {"det": det, "rows": rows}

    day = EXAMPLE_DAY if EXAMPLE_DAY in days else days[len(days) // 2]
    out["day"] = example_day(day, CHAN_DET)

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as handle:
        json.dump(out, handle, separators=(",", ":"))
    print(f"wrote {OUT} ({os.path.getsize(OUT)/1e3:.0f} KB)")
    print(f"  folds for {len(out['fold'])} detectors, fits for {len(out['fit'])}")
    for det in ("n0", "nb", "b0"):
        if det in out["fit"]:
            f = out["fit"][det]
            print(f"    {det}: tau {f['tau']:5.2f} min   peak excess "
                  f"{f['peak_excess']:+.1f}%")


if __name__ == "__main__":
    main()
