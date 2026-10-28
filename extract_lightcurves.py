"""
Extract lightcurves for every eligible burst, as compact arrays the browser
can draw.

Deliberately does NOT run the DoL fit. A lightcurve needs only the summed
triggered detectors, the background polynomial and the two candidate windows,
so this pass is far cheaper than the localization run -- roughly a second per
burst rather than eleven.

Stored per burst, all rounded hard because the data is drawn into a panel a
few hundred pixels wide:
    e   bin edges, seconds from trigger
    r   count rate per bin, summed triggered detectors, 50-300 keV
    b   fitted background rate per bin
    w90 the full T90 window
    wsn the best signal-to-noise window
    wbg the background windows
    ts  binning in ms

Resumable: a burst whose JSON exists is skipped.

Usage:
    python extract_lightcurves.py --workers 9
"""
import argparse
import csv
import json
import os
import warnings

warnings.filterwarnings("ignore")

import numpy as np

DATADIR = "data/gbm"
MANIFEST = "catalog/fermigbrst.csv"
OUTDIR = "site/vec/lc"
SHORT_LONG_CUT = 2.0


def fnum(value, default=None):
    try:
        out = float(value)
        return out if np.isfinite(out) else default
    except (TypeError, ValueError):
        return default


def process(trigger, row):
    import analyze_grbs as A
    from gdt.missions.fermi.gbm.trigdat import Trigdat
    from gdt.core.background.fitter import BackgroundFitter
    from gdt.core.background.binned import Polynomial

    bdir = os.path.join(DATADIR, trigger)
    trigdat_file = next(f for f in sorted(os.listdir(bdir))
                        if f.startswith("glg_trigdat_"))
    trigdat = Trigdat.open(os.path.join(bdir, trigdat_file))

    t90 = float(row["t90"])
    t90_start = float(row["t90_start"])
    probe = trigdat.to_phaii(trigdat._detectors[0], timescale=64)
    span = (float(probe.data.tstart[0]), float(probe.data.tstop[-1]))
    timescale, t90_window, bg_times = A.choose_windows(t90, t90_start, span)

    found = A.peak_window(trigdat, timescale, t90_window, bg_times)
    sn_window = [found[0], found[1]] if found else list(t90_window)

    summed = trigdat.sum_detectors(trigdat.triggered_detectors,
                                   timescale=timescale)
    summed = summed.slice_energy(*(A.LOC_ERANGE,))
    fitter = BackgroundFitter.from_phaii(summed, Polynomial,
                                         time_ranges=bg_times)
    fitter.fit(order=1)
    background = fitter.interpolate_bins(summed.data.tstart, summed.data.tstop)

    tstart = np.asarray(summed.data.tstart, dtype=float)
    tstop = np.asarray(summed.data.tstop, dtype=float)

    # Divide by exposure, not bin width. trigdat's wide pre/post-trigger bins
    # carry real deadtime, and in a high-particle environment the livetime can
    # fall to half: bn250223692 runs at 45-100% livetime across its span.
    # Using width there understates the rate by up to 2x and, worse, plots the
    # data in different units from the fitted background, which is already
    # deadtime-corrected -- so the model looks wrong against a flat-looking
    # curve that is itself an artefact.
    exposure = np.clip(np.asarray(summed.data.exposure, dtype=float), 1e-6, None)
    rate = np.asarray(summed.data.counts).sum(axis=1) / exposure
    bkg = np.asarray(background.rates).sum(axis=1)

    # Livetime stability. A polynomial background assumes the detectors are
    # counting under roughly constant conditions. In a high-particle
    # environment that fails: bn250223692 swings between 45% and 100% livetime
    # as the spacecraft leaves the region, so the deadtime-corrected background
    # is genuinely non-stationary and no polynomial order recovers it. This is
    # the cause, measured directly, rather than a downstream symptom such as
    # the model extrapolating negative far from any window that is used.
    width = np.clip(tstop - tstart, 1e-9, None)
    live = exposure / width
    full = trigdat.sum_detectors(trigdat.triggered_detectors, timescale=timescale)
    allrate = float(np.mean(np.asarray(full.data.counts).sum(axis=1) / width))

    # Edges rather than centres, so the step plot is unambiguous where trigdat
    # switches resolution and bin widths change.
    edges = np.concatenate([tstart, [tstop[-1]]])

    return {
        "t": trigger,
        "ts": timescale,
        "e": [round(float(v), 3) for v in edges],
        "r": [int(round(float(v))) for v in rate],
        "b": [int(round(float(v))) for v in bkg],
        "w90": [round(float(t90_window[0]), 3), round(float(t90_window[1]), 3)],
        "wsn": [round(float(sn_window[0]), 3), round(float(sn_window[1]), 3)],
        "wbg": [[round(float(a), 1), round(float(b), 1)] for a, b in bg_times],
        "lv": [round(float(live.min()), 3), round(float(np.ptp(live)), 3),
               round(allrate)],
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
    parser.add_argument("--workers", type=int, default=9)
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()

    os.makedirs(OUTDIR, exist_ok=True)
    rows = {r["trigger_name"].strip(): r
            for r in csv.DictReader(open(MANIFEST, newline=""))}

    # Only bursts that already have contour data: the viewer pairs the two, so
    # a lightcurve without a sky plot would have nothing to attach to.
    have = {os.path.splitext(f)[0]
            for f in os.listdir("site/vec/bursts") if f.endswith(".json")}
    eligible = sorted(t for t in have if t in rows)
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
            if done % 200 == 0 or done == len(tasks):
                rate = done / max(time.time() - started, 1e-9)
                print(f"  {done}/{len(tasks)}  ok={ok} failed={failed}  "
                      f"{rate:.1f}/s  ETA {(len(tasks)-done)/rate/60:.1f} min",
                      flush=True)

    print(f"\n{ok} new, {cached} cached, {failed} failed")
    total = sum(os.path.getsize(os.path.join(OUTDIR, f))
                for f in os.listdir(OUTDIR) if f.endswith(".json"))
    print(f"lightcurve data: {total/1e6:.1f} MB across {len(os.listdir(OUTDIR))} files")
    for trigger, status in failures[:15]:
        print(f"  {trigger}: {status}")
    if failures:
        with open("site/vec/lc_failures.json", "w") as handle:
            json.dump(failures, handle, indent=2)


if __name__ == "__main__":
    main()
