"""
Extract 90% confidence contours for every eligible burst, as compact vector
data the browser can draw.

The PNG gallery cannot back a dropdown over 2,049 bursts: at ~40 KB each the
images total 81 MB against a 16 MB page budget, and artifacts allow only 64
supporting files. Contour polygons are far smaller -- about 2 KB per burst once
subsampled and rounded -- so the whole catalog fits in one page and every
burst gets a real plot rather than only the few hundred that fit as images.

Per burst, stores the 90% region for:
    official   the published Fermi GBM HEALPix map
    t90        this reproduction, full T90 source window
    peak       this reproduction, best signal-to-noise sub-interval
plus centroids, the external mission position where one exists, and the
numbers behind the comparison.

Resumable: a burst whose JSON exists is skipped.

Usage:
    python extract_contours.py --workers 9
    python extract_contours.py --limit 10
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
EXTERNAL = "catalog/external_localizations.csv"
OUTDIR = "site/vec/bursts"
SHORT_LONG_CUT = 2.0

# A contour drawn into a ~460 px panel needs nothing like its full resolution.
MAX_POINTS = 56
COORD_DP = 2


def fnum(value, default=None):
    try:
        out = float(value)
        return out if np.isfinite(out) else default
    except (TypeError, ValueError):
        return default


def pack(paths):
    """Subsample and round contour paths into flat [x,y,x,y,...] arrays."""
    out = []
    for path in paths:
        if len(path) < 3:
            continue
        step = max(1, len(path) // MAX_POINTS)
        pts = path[::step]
        flat = []
        for x, y in pts:
            flat.append(round(float(x), COORD_DP))
            flat.append(round(float(y), COORD_DP))
        out.append(flat)
    return out


def contours(hpx):
    try:
        return pack(hpx.confidence_region_path(0.9))
    except Exception:
        return []


def process(trigger, row, ext):
    import analyze_grbs as A
    A.MAKE_PLOTS = False

    art_t90, art_peak = {}, {}
    res_t90 = A.localize(trigger, row, window_mode="t90", artifacts=art_t90)
    res_peak = A.localize(trigger, row, window_mode="peak", artifacts=art_peak)

    official = art_peak["official"]
    rec = {
        "t": trigger,
        "n": res_peak["name"],
        "k": res_peak["grb_type"],
        "t90": round(res_peak["t90"], 2),
        "det": res_peak["triggered_detectors"],
        "oc": [round(float(official.centroid[0]), 3),
               round(float(official.centroid[1]), 3)],
        "oa": round(res_peak["official_area90"]),
        "ac": [round(res_t90["my_ra"], 3), round(res_t90["my_dec"], 3)],
        "aa": round(res_t90["my_area90"]),
        "as": round(res_t90["sep_mine_official"], 2),
        "aw": round(res_t90["src_duration"], 2),
        "pc": [round(res_peak["my_ra"], 3), round(res_peak["my_dec"], 3)],
        "pa": round(res_peak["my_area90"]),
        "ps": round(res_peak["sep_mine_official"], 2),
        "pw": round(res_peak["src_duration"], 2),
        "cat": [round(res_peak["catalog_ra"], 3), round(res_peak["catalog_dec"], 3)],
        "O": contours(official),
        "A": contours(art_t90["sys_map"]),
        "P": contours(art_peak["sys_map"]),
    }
    if ext is not None:
        from astropy.coordinates import SkyCoord
        target = SkyCoord(ext["ra"], ext["dec"], unit="deg")
        rec["x"] = {
            "s": ext["source"], "n": ext["name"],
            "c": [round(ext["ra"], 4), round(ext["dec"], 4)],
            "do": round(float(SkyCoord(*rec["oc"], unit="deg").separation(target).deg), 2),
            "da": round(float(SkyCoord(*rec["ac"], unit="deg").separation(target).deg), 2),
            "dp": round(float(SkyCoord(*rec["pc"], unit="deg").separation(target).deg), 2),
        }
    return rec


def worker(args):
    trigger, row, ext = args
    out = os.path.join(OUTDIR, f"{trigger}.json")
    if os.path.exists(out):
        return trigger, "cached"
    try:
        rec = process(trigger, row, ext)
        tmp = out + ".part"
        with open(tmp, "w") as handle:
            json.dump(rec, handle, separators=(",", ":"))
        os.replace(tmp, out)
        return trigger, "ok"
    except Exception as err:
        return trigger, f"FAIL {type(err).__name__}: {err}"


def load_external():
    if not os.path.exists(EXTERNAL):
        return {}
    out = {}
    for row in csv.DictReader(open(EXTERNAL, newline="")):
        if row.get("match_quality") == "suspect":
            continue
        ra, dec = fnum(row["swift_ra"]), fnum(row["swift_dec"])
        source, name = "Swift", row["swift_name"]
        if ra is None or dec is None:
            ra, dec = fnum(row["lat_ra"]), fnum(row["lat_dec"])
            source, name = "Fermi-LAT", row["lat_name"]
        if ra is None or dec is None:
            continue
        out[row["trigger_name"]] = {"ra": ra, "dec": dec,
                                    "source": source, "name": name}
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=9)
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()

    os.makedirs(OUTDIR, exist_ok=True)
    rows = {r["trigger_name"].strip(): r
            for r in csv.DictReader(open(MANIFEST, newline=""))}
    ext_all = load_external()

    eligible = []
    for trigger, row in rows.items():
        bdir = os.path.join(DATADIR, trigger)
        if not os.path.isdir(bdir):
            continue
        names = os.listdir(bdir)
        if not any(n.startswith("glg_healpix_") for n in names):
            continue
        if not any(n.startswith("glg_trigdat_") for n in names):
            continue
        t90, t90_start = fnum(row["t90"]), fnum(row["t90_start"])
        if t90 is None or t90_start is None or t90 <= 0:
            continue
        eligible.append(trigger)
    eligible.sort()
    if args.limit:
        eligible = eligible[: args.limit]

    print(f"{len(eligible)} bursts, {args.workers} workers", flush=True)
    tasks = [(t, rows[t], ext_all.get(t)) for t in eligible]

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
    total = sum(os.path.getsize(os.path.join(OUTDIR, f))
                for f in os.listdir(OUTDIR) if f.endswith(".json"))
    print(f"contour data: {total/1e6:.1f} MB across "
          f"{len(os.listdir(OUTDIR))} files")
    if failures:
        with open("site/vec/failures.json", "w") as handle:
            json.dump(failures, handle, indent=2)


if __name__ == "__main__":
    main()
