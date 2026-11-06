"""
Per-detector lightcurves from TTE, as compact arrays the browser can draw.

GBM has no optics. It localizes by comparing how much each of its 14 detectors
saw, so the per-detector view is not a decoration -- it is the actual input to
the fit. A burst in front of n6 and n7 and nothing in n0-n3 is a burst on that
side of the spacecraft, and the panel shows that directly.

Resolution is the reason this reads TTE rather than trigdat. trigdat's bins
are a fixed template keyed to the trigger -- 64 ms only from -0.3 to +0.5 s,
then 256 ms, then 1.024 s, then 8.192 s -- so a 0.1 s short burst lands in one
or two bins and a long burst's tail is sampled eight seconds at a time. TTE
timestamps events to ~4 us, so the window can be chosen to suit the burst.

Sizing. The viewer embeds everything in one page against a 16 MB cap, and the
sky contours and summed lightcurve already account for most of it. So the 14
curves become one base64 uint8 blob rather than JSON numbers: 14 x 80 bins is
1120 bytes raw, about 1.5 KB encoded, against roughly 16 KB as text. 80 bins
is not a compromise at the size these are drawn -- each panel renders about
180 px wide, so the screen is the limit, not the data. It still beats trigdat
comfortably where it matters: a 2 s window round a short burst comes out at
25 ms bins, against trigdat's 64 ms floor which in any case only covers
-0.3 to +0.5 s.

Scaling is per GROUP, not per detector. All 12 NaIs share one scale, because
comparing them is the entire point and normalizing each to its own maximum
would make every detector look equally bright. The two BGOs get their own,
because their band is twenty times wider and their background correspondingly
higher -- on a shared scale they peg the top and squash all twelve NaI curves
into the bottom third. Each group's byte range also spans [min, max] rather
than [0, max]: a burst that adds 100 counts/s on a 1000 counts/s pedestal
would otherwise get about 25 of the 256 available levels.

Stored per burst:
    lo, hi  display window, seconds from trigger
    nb      bins per detector
    sn      [base, span] counts/s for the NaI group: rate = base + byte/255*span
    sb      [base, span] the same for the two BGOs
    dets    detector names, in the order they appear in the blob
    d       one base64 uint8 blob, len(dets) * nb bytes, detectors concatenated
    net     {detector: net counts in the T90 window, background subtracted}
    dead    detectors the run found dead
    trg     detectors that triggered
    eb      {nai: [lo, hi], bgo: [lo, hi]} display energy bands, keV

Resumable: a burst whose JSON exists is skipped.

Usage:
    python extract_detectors.py --workers 6
    python extract_detectors.py --limit 20
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
OUTDIR = "site/vec/det"

NAI = [f"n{d}" for d in "0123456789ab"]
BGO = ["b0", "b1"]
DETS = NAI + BGO

# Display bands. The NaIs are useful from the 8 keV electronics threshold to
# about 900 keV, above which their efficiency has collapsed; the BGOs cover
# the decade above that. These are for looking at the burst, and are wider
# than the 50-300 keV band the localization fit uses.
NAI_BAND = (8.0, 900.0)
BGO_BAND = (200.0, 40000.0)

NBINS = 80

# Set from --allow-partial; a list so complete_bursts() can see it.
ALLOW_PARTIAL = [False]


def fnum(value, default=None):
    try:
        out = float(value)
        return out if np.isfinite(out) else default
    except (TypeError, ValueError):
        return default


def window_for(t90, t90_start):
    """Display span: the T90 window with a margin that scales with it.

    A fixed margin cannot serve both ends of the catalog. Ten seconds either
    side of a 0.05 s short burst is 400 bins of pure background; the same ten
    seconds around a 200 s burst shows almost no background at all.
    """
    margin = max(2.0 * t90, 1.0)
    return t90_start - margin, t90_start + t90 + margin


def counts_in(times, lo, hi, nbins):
    """Exact event counts per bin.

    np.histogram with explicit edges, rather than any rebinning helper: a
    helper that lands a partial bin at the edge inflates the last bin, which
    on an earlier pass here showed up as a spurious +32% in the final bin.
    """
    edges = np.linspace(lo, hi, nbins + 1)
    counts, _ = np.histogram(times, bins=edges)
    return counts, edges


def process(trigger, row):
    from gdt.missions.fermi.gbm.tte import GbmTte

    bdir = os.path.join(TTEDIR, trigger)
    t90 = fnum(row["t90"])
    t90_start = fnum(row["t90_start"])
    if t90 is None or t90_start is None or t90 <= 0:
        raise ValueError("no usable t90")
    lo, hi = window_for(t90, t90_start)
    width = (hi - lo) / NBINS

    files = {}
    for name in os.listdir(bdir):
        m = re.match(r"glg_tte_([nb][0-9ab])_", name)
        if m:
            files[m.group(1)] = os.path.join(bdir, name)
    if not files:
        raise ValueError("no TTE files")

    rates, nets = {}, {}
    for det, path in sorted(files.items()):
        tte = GbmTte.open(path)
        data = tte.data
        band = BGO_BAND if det.startswith("b") else NAI_BAND
        eb = data.ebounds
        low = np.asarray(eb.low_edges(), dtype=float)
        high = np.asarray(eb.high_edges(), dtype=float)
        chans = np.flatnonzero((low >= band[0]) & (high <= band[1]))
        if chans.size == 0:
            continue
        # data.times is ALREADY relative to the trigger. Subtracting trigtime
        # again shifts every event by ~5e8 s and empties the window; that cost
        # a debugging session once, hence the note.
        sel = np.isin(data.channels, chans)
        times = np.asarray(data.times, dtype=float)[sel]
        counts, edges = counts_in(times, lo, hi, NBINS)
        rate = counts / width

        # Net counts in the T90 window, over a background taken from the
        # outer thirds of the display span, which the margin above guarantees
        # are outside the burst.
        centres = 0.5 * (edges[:-1] + edges[1:])
        inside = (centres >= t90_start) & (centres <= t90_start + t90)
        outside = (centres < lo + (hi - lo) / 3.0) | (centres > hi - (hi - lo) / 3.0)
        base = float(np.median(rate[outside])) if outside.any() else 0.0
        nets[det] = int(round(float((rate[inside] - base).sum() * width)))
        rates[det] = rate

    if not rates:
        raise ValueError("no detector produced a curve")

    def group_range(dets):
        members = [rates[d] for d in dets if d in rates]
        if not members:
            return 0.0, 1.0
        lo_r = float(min(r.min() for r in members))
        hi_r = float(max(r.max() for r in members))
        return lo_r, max(hi_r - lo_r, 1e-6)

    nai_base, nai_span = group_range(NAI)
    bgo_base, bgo_span = group_range(BGO)

    order = [d for d in DETS if d in rates]
    blob = []
    for det in order:
        base, span = (bgo_base, bgo_span) if det.startswith("b") else (nai_base, nai_span)
        blob.append(np.clip(np.rint((rates[det] - base) / span * 255.0),
                            0, 255).astype(np.uint8))

    return {
        "t": trigger,
        "lo": round(lo, 4), "hi": round(hi, 4), "nb": NBINS,
        "sn": [round(nai_base, 1), round(nai_span, 1)],
        "sb": [round(bgo_base, 1), round(bgo_span, 1)],
        "dets": order,
        "d": base64.b64encode(np.concatenate(blob).tobytes()).decode("ascii"),
        "net": [nets.get(d, 0) for d in order],
        "eb": {"nai": list(NAI_BAND), "bgo": list(BGO_BAND)},
    }


def add_fit_detectors(rec, trigger):
    """Which detectors the localization actually used, from the saved run."""
    path = f"site/all/results/{trigger}.json"
    if not os.path.exists(path):
        return
    try:
        with open(path) as handle:
            res = json.load(handle)
    except Exception:
        return
    trg = res.get("triggered_detectors")
    if trg:
        rec["trg"] = list(trg)
    dead = res.get("dead_detectors")
    if dead:
        rec["dead"] = list(dead)


def worker(args):
    trigger, row = args
    out = os.path.join(OUTDIR, f"{trigger}.json")
    if os.path.exists(out):
        return trigger, "cached"
    try:
        rec = process(trigger, row)
        add_fit_detectors(rec, trigger)
        tmp = out + ".part"
        with open(tmp, "w") as handle:
            json.dump(rec, handle, separators=(",", ":"))
        os.replace(tmp, out)
        return trigger, "ok"
    except Exception as err:
        return trigger, f"FAIL {type(err).__name__}: {err}"


def complete_bursts():
    """Bursts whose TTE download has finished.

    All 14 detectors are required, not merely most of them. This runs while
    the download is still going, and a burst caught mid-flight would otherwise
    be written out missing whichever files had not landed yet -- and then
    skipped as cached on the next pass, leaving a permanently incomplete
    panel. Bursts genuinely missing a detector in the archive are rare, and
    are picked up by the final --allow-partial pass.
    """
    out = []
    for d in sorted(glob.glob(f"{TTEDIR}/*")):
        if not os.path.isdir(d):
            continue
        dets = set()
        for name in os.listdir(d):
            m = re.match(r"glg_tte_([nb][0-9ab])_", name)
            if m and not name.endswith(".part"):
                dets.add(m.group(1))
        # Eight NaIs is the --allow-partial floor. HEASARC really does omit
        # TTE for some detectors on some triggers -- 17 bursts here, the worst
        # being bn181221565 with n0, n1 and n2 absent -- and a panel showing
        # the nine that exist is better than no panel at all.
        if dets >= set(DETS) or (ALLOW_PARTIAL[0] and len(dets & set(NAI)) >= 8):
            out.append(os.path.basename(d))
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--allow-partial", action="store_true",
                        help="accept bursts the archive is missing a detector "
                             "for; run only once the download has finished")
    args = parser.parse_args()

    ALLOW_PARTIAL[0] = args.allow_partial
    os.makedirs(OUTDIR, exist_ok=True)
    rows = {r["trigger_name"].strip(): r
            for r in csv.DictReader(open(MANIFEST, newline=""))}
    have = {os.path.splitext(f)[0]
            for f in os.listdir("site/vec/bursts") if f.endswith(".json")}

    eligible = [t for t in complete_bursts() if t in rows and t in have]
    if args.limit:
        eligible = eligible[: args.limit]
    print(f"{len(eligible)} bursts with complete TTE, {args.workers} workers",
          flush=True)
    if not eligible:
        return

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
            if done % 50 == 0 or done == len(tasks):
                rate = done / max(time.time() - started, 1e-9)
                print(f"  {done}/{len(tasks)}  ok={ok} cached={cached} "
                      f"failed={failed}  ETA {(len(tasks)-done)/rate/60:.0f} min",
                      flush=True)

    print(f"\n{ok} new, {cached} cached, {failed} failed")
    files = [f for f in os.listdir(OUTDIR) if f.endswith(".json")]
    total = sum(os.path.getsize(os.path.join(OUTDIR, f)) for f in files)
    print(f"detector data: {total/1e6:.1f} MB across {len(files)} files "
          f"({total/max(len(files),1)/1024:.1f} KB each)")
    for trigger, status in failures[:15]:
        print(f"  {trigger}: {status}")
    if failures:
        with open("site/vec/det_failures.json", "w") as handle:
            json.dump(failures, handle, indent=2)


if __name__ == "__main__":
    main()
