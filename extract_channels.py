"""
Energy-resolved views of each burst, from TTE: a spectrogram and a set of
band lightcurves.

The detector panel answers "where was it". This answers "what was it doing
spectrally". GRBs are not static in energy -- most decay hard to soft, pulses
peak earlier at high energy than low, and some have a hard precursor that the
summed lightcurve hides entirely. None of that is visible in a single curve.

TTE carries 128 energy channels per detector, which is far more than can be
drawn as separate panels and more than the detectors can actually resolve: a
NaI has roughly 10% FWHM energy resolution, so 128 pseudo-logarithmic
channels over 4.5 keV to 2 MeV heavily oversample it. They are rebinned here
onto 24 logarithmic channels, which loses nothing physical and costs a
sixteenth of the page budget.

Two products per burst, both summed over the detectors that triggered:

  spectrogram   24 log channels x 48 time bins, in units of significance,
                (counts - background) / sqrt(background). Significance rather
                than raw rate because the background spans more than two
                orders of magnitude from 8 keV to 900 keV; a raw heatmap is
                simply a picture of the background, with the burst invisible.

  bands         8 energy bands as separate lightcurves, each scaled to its
                OWN range. This is the opposite of the choice made for the
                detector panel, and deliberately so: there the amplitude
                comparison between detectors is the localization signal, so a
                shared scale is essential. Here the interesting thing is how
                the pulse SHAPE changes with energy -- lags, asymmetry,
                hard-to-soft evolution -- and a shared scale would flatten
                the high bands into invisibility, since the 8-15 keV band
                outcounts 600-900 keV by a factor of a hundred.

Stored per burst:
    lo, hi   display window, seconds from trigger (same as the detector panel)
    nt, nc   time bins and energy channels in the spectrogram
    sg       base64 uint8, nc*nt bytes, channel-major
    sr       [lo, hi] significance mapped onto 0..255
    bd       base64 uint8, 8*nt bytes, one band per row
    bs       per-band [min, max] counts/s
    dets     detectors summed
    tot      net counts per band over the T90 window

Resumable: a burst whose JSON exists is skipped.

Usage:
    python extract_channels.py --workers 8
    python extract_channels.py --limit 10
"""
import argparse
import base64
import csv
import glob
import json
import os
import re
import warnings

warnings.filterwarnings("ignore")

import numpy as np

TTEDIR = "data/tte"
MANIFEST = "catalog/fermigbrst.csv"
RESULTS = "site/all/results"
OUTDIR = "site/vec/chan"

NAI = [f"n{d}" for d in "0123456789ab"]

NTIME = 48
NCHAN = 24
EMIN, EMAX = 8.0, 900.0

# Round-numbered rather than exactly logarithmic, so the labels read cleanly.
BANDS = [8, 15, 25, 50, 100, 200, 350, 600, 900]

# Bottom of the colour scale. Zero, not a negative number: anchoring at
# -5 sigma put zero significance in the middle of the ramp, so background
# noise rendered as mid-grey and the burst had only the top half of the range
# to stand out in. Clipping negative fluctuations away costs nothing for a
# display and puts the quiet sky at the dark end where it belongs.
SIG_FLOOR = 0.0

# Floor on the top of the scale too, so a burst with no real signal does not
# get its noise stretched across the full ramp and made to look structured.
SIG_CEIL_MIN = 5.0


def fnum(value, default=None):
    try:
        out = float(value)
        return out if np.isfinite(out) else default
    except (TypeError, ValueError):
        return default


def window_for(t90, t90_start):
    """Identical to the detector panel's window, so the two line up."""
    margin = max(2.0 * t90, 1.0)
    return t90_start - margin, t90_start + t90 + margin


def chosen_detectors(trigger, available):
    """Sum the detectors that triggered; fall back to all NaI."""
    path = f"{RESULTS}/{trigger}.json"
    if os.path.exists(path):
        try:
            with open(path) as handle:
                trg = json.load(handle).get("triggered_detectors") or []
        except Exception:
            trg = []
        keep = [d for d in trg if d in available and d.startswith("n")]
        if keep:
            return keep
    return [d for d in NAI if d in available]


def pack(arr, lo, hi):
    """uint8 over [lo, hi]."""
    span = max(hi - lo, 1e-9)
    q = np.clip(np.rint((arr - lo) / span * 255.0), 0, 255).astype(np.uint8)
    return base64.b64encode(q.tobytes()).decode("ascii")


def process(trigger, row):
    from gdt.missions.fermi.gbm.tte import GbmTte

    t90 = fnum(row["t90"])
    t90_start = fnum(row["t90_start"])
    if t90 is None or t90_start is None or t90 <= 0:
        raise ValueError("no usable t90")
    lo, hi = window_for(t90, t90_start)
    tedges = np.linspace(lo, hi, NTIME + 1)
    width = (hi - lo) / NTIME

    bdir = os.path.join(TTEDIR, trigger)
    files = {}
    for name in os.listdir(bdir):
        m = re.match(r"glg_tte_(n[0-9ab])_", name)
        if m:
            files[m.group(1)] = os.path.join(bdir, name)
    dets = chosen_detectors(trigger, files)
    if not dets:
        raise ValueError("no NaI detectors available")

    eedges = np.logspace(np.log10(EMIN), np.log10(EMAX), NCHAN + 1)
    spec = np.zeros((NCHAN, NTIME))
    band = np.zeros((len(BANDS) - 1, NTIME))

    for det in dets:
        tte = GbmTte.open(files[det])
        data = tte.data
        eb = data.ebounds
        low = np.asarray(eb.low_edges(), dtype=float)
        high = np.asarray(eb.high_edges(), dtype=float)
        # Geometric centre, since the channels are quasi-logarithmic.
        centre = np.sqrt(np.maximum(low, 1e-3) * np.maximum(high, 1e-3))
        # data.times is ALREADY trigger-relative; subtracting trigtime again
        # empties the window.
        times = np.asarray(data.times, dtype=float)
        chan = np.asarray(data.channels)
        energy = centre[np.clip(chan, 0, centre.size - 1)]
        keep = (energy >= EMIN) & (energy < EMAX) & (times >= lo) & (times < hi)
        if not keep.any():
            continue
        t_sel, e_sel = times[keep], energy[keep]
        spec += np.histogram2d(e_sel, t_sel, bins=[eedges, tedges])[0]
        band += np.histogram2d(e_sel, t_sel, bins=[np.array(BANDS, float),
                                                   tedges])[0]

    if spec.sum() == 0:
        raise ValueError("no counts in the window")

    # Background from the outer thirds, which the margin keeps clear of the
    # burst. Per channel, because the spectrum of the background is steep.
    centres = 0.5 * (tedges[:-1] + tedges[1:])
    outside = (centres < lo + (hi - lo) / 3.0) | (centres > hi - (hi - lo) / 3.0)
    if not outside.any():
        outside = np.ones(NTIME, dtype=bool)

    bg = np.median(spec[:, outside], axis=1)[:, None]
    sig = (spec - bg) / np.sqrt(np.maximum(bg, 1.0))
    shi = float(max(np.percentile(sig, 99.9), SIG_CEIL_MIN))
    sg = pack(sig, SIG_FLOOR, shi)

    brate = band / width
    bmin = brate.min(axis=1)
    bmax = brate.max(axis=1)
    rows = []
    for i in range(brate.shape[0]):
        span = max(float(bmax[i] - bmin[i]), 1e-9)
        rows.append(np.clip(np.rint((brate[i] - bmin[i]) / span * 255.0),
                            0, 255).astype(np.uint8))
    bd = base64.b64encode(np.concatenate(rows).tobytes()).decode("ascii")

    bbg = np.median(brate[:, outside], axis=1)
    inside = (centres >= t90_start) & (centres <= t90_start + t90)
    tot = [int(round(float(((brate[i] - bbg[i])[inside]).sum() * width)))
           for i in range(brate.shape[0])]

    return {
        "t": trigger,
        "lo": round(lo, 4), "hi": round(hi, 4),
        "nt": NTIME, "nc": NCHAN,
        "sg": sg, "sr": [SIG_FLOOR, round(shi, 2)],
        "bd": bd,
        "bs": [[round(float(bmin[i]), 1), round(float(bmax[i]), 1)]
               for i in range(brate.shape[0])],
        "dets": dets,
        "tot": tot,
    }


def worker(args):
    trigger, row = args
    out = os.path.join(OUTDIR, f"{trigger}.json")
    if os.path.exists(out):
        return trigger, "cached"
    try:
        rec = process(trigger, row)
        tmp = out + ".part"
        with open(tmp, "w") as handle:
            json.dump(rec, handle, separators=(",", ":"))
        os.replace(tmp, out)
        return trigger, "ok"
    except Exception as err:
        return trigger, f"FAIL {type(err).__name__}: {err}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()

    os.makedirs(OUTDIR, exist_ok=True)
    rows = {r["trigger_name"].strip(): r
            for r in csv.DictReader(open(MANIFEST, newline=""))}
    have = {os.path.splitext(f)[0]
            for f in os.listdir("site/vec/bursts") if f.endswith(".json")}
    eligible = sorted(t for t in have
                      if t in rows and os.path.isdir(os.path.join(TTEDIR, t)))
    if args.limit:
        eligible = eligible[: args.limit]

    print(f"{len(eligible)} bursts, {args.workers} workers", flush=True)
    tasks = [(t, rows[t]) for t in eligible]

    from concurrent.futures import ProcessPoolExecutor, as_completed
    import time
    done = ok = cached = failed = 0
    failures = []
    started = time.time()
    with ProcessPoolExecutor(args.workers) as pool:
        futures = {pool.submit(worker, task): task[0] for task in tasks}
        for future in as_completed(futures):
            trigger, status = future.result()
            done += 1
            if status == "ok":
                ok += 1
            elif status == "cached":
                cached += 1
            else:
                failed += 1
                failures.append((trigger, status))
            if done % 100 == 0 or done == len(tasks):
                rate = done / max(time.time() - started, 1e-9)
                print(f"  {done}/{len(tasks)}  ok={ok} cached={cached} "
                      f"failed={failed}  ETA {(len(tasks)-done)/rate/60:.0f} min",
                      flush=True)

    print(f"\n{ok} new, {cached} cached, {failed} failed")
    files = [f for f in os.listdir(OUTDIR) if f.endswith(".json")]
    total = sum(os.path.getsize(os.path.join(OUTDIR, f)) for f in files)
    print(f"channel data: {total/1e6:.1f} MB across {len(files)} files "
          f"({total/max(len(files),1)/1024:.1f} KB each)")
    for trigger, status in failures[:15]:
        print(f"  {trigger}: {status}")
    if failures:
        with open("site/vec/chan_failures.json", "w") as handle:
            json.dump(failures, handle, indent=2)


if __name__ == "__main__":
    main()
