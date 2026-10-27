"""
Localize every eligible GBM burst both ways and plot the result against the
official map and, where one exists, an external mission's position.

For each burst:
    - DoL with the full T90 source window
    - DoL with the best signal-to-noise sub-interval
    - one figure carrying both, the official Fermi HEALPix region, and the
      external position if another instrument pinned the burst down

Runs one process per burst across cores, and is resumable: a burst whose
result JSON already exists is skipped, so an interrupted run picks up where it
stopped.

Usage:
    python run_all.py                      # everything eligible
    python run_all.py --limit 20           # first 20, for a smoke test
    python run_all.py --workers 10
    python run_all.py --class short        # only short bursts
"""
import argparse
import csv
import json
import os
import traceback
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATADIR = "data/gbm"
MANIFEST = "catalog/fermigbrst.csv"
EXTERNAL = "catalog/external_localizations.csv"
RESULTDIR = "site/all/results"
FIGDIR_ALL = "site/all/figures"

OFFICIAL_C = "#2f6fd0"      # official Fermi map
T90_C = "#b8690a"           # my localization, full T90 window
PEAK_C = "#d9534f"          # my localization, best-S/N window
EXT_C = "#1f9e6e"           # another mission's position
SHORT_LONG_CUT = 2.0


def fnum(value, default=None):
    try:
        out = float(value)
        return out if np.isfinite(out) else default
    except (TypeError, ValueError):
        return default


def wrap_near(values, center):
    return center + (np.asarray(values, dtype=float) - center + 180.0) % 360.0 - 180.0


def four_way_figure(trigger, art_t90, art_peak, ext, res_t90, res_peak, path):
    """Both window modes, the official region, and any external position."""
    fig, ax = plt.subplots(figsize=(5.0, 4.4), dpi=96)
    official = art_peak["official"]
    center_ra = official.centroid[0]
    center_dec = official.centroid[1]

    # T90 is dashed and drawn wider: the two window modes often agree closely,
    # and a solid line underneath would simply vanish under the peak contour.
    layers = [
        (official, OFFICIAL_C, "official Fermi GBM", 1.9, "-"),
        (art_t90["sys_map"], T90_C, "mine, full T90", 1.9, "--"),
        (art_peak["sys_map"], PEAK_C, "mine, best S/N", 1.4, "-"),
    ]
    spans = []
    for hpx, color, label, lw, ls in layers:
        try:
            paths = hpx.confidence_region_path(0.9)
        except Exception:
            continue
        for k, path_xy in enumerate(paths):
            ra = wrap_near(path_xy[:, 0], center_ra)
            dec = path_xy[:, 1]
            spans.append((ra, dec))
            ax.plot(ra, dec, color=color, lw=lw, ls=ls,
                    label=label if k == 0 else None)

    ax.plot(wrap_near([official.centroid[0]], center_ra), [official.centroid[1]],
            marker="P", color=OFFICIAL_C, ms=8, mec="white", mew=.8, ls="none")
    ax.plot(wrap_near([res_peak["my_ra"]], center_ra), [res_peak["my_dec"]],
            marker="X", color=PEAK_C, ms=8, mec="white", mew=.8, ls="none")

    if ext is not None:
        era = wrap_near([ext["ra"]], center_ra)
        ax.plot(era, [ext["dec"]], marker="*", color=EXT_C, ms=15,
                mec="white", mew=.9, ls="none",
                label=f"{ext['source']} position")
        spans.append((np.array(era, dtype=float), np.array([ext["dec"]], dtype=float)))

    if spans:
        all_ra = np.concatenate([s[0] for s in spans])
        all_dec = np.concatenate([s[1] for s in spans])
        pad = 0.09 * max(np.ptp(all_ra), np.ptp(all_dec), 2.0)
        ax.set_xlim(all_ra.min() - pad, all_ra.max() + pad)
        ax.set_ylim(all_dec.min() - pad, all_dec.max() + pad)

    ax.invert_xaxis()
    ax.set_aspect(1.0 / max(np.cos(np.radians(center_dec)), 0.05))
    ax.set_xlabel("RA (deg)", fontsize=8)
    ax.set_ylabel("Dec (deg)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.grid(alpha=.22, lw=.5)
    ax.set_title(
        f"{trigger}   T90 {res_peak['t90']:.2f}s ({res_peak['grb_type']})\n"
        f"T90 win {res_t90['sep_mine_official']:.2f}deg   "
        f"best-S/N {res_peak['sep_mine_official']:.2f}deg",
        fontsize=8)
    ax.legend(fontsize=6.3, loc="best", framealpha=.9)
    plt.tight_layout()
    plt.savefig(path, bbox_inches="tight")
    plt.close("all")


def process(trigger, row, ext):
    """Both window modes plus the figure for one burst."""
    import analyze_grbs as A
    A.MAKE_PLOTS = False

    art_t90, art_peak = {}, {}
    res_t90 = A.localize(trigger, row, window_mode="t90", artifacts=art_t90)
    res_peak = A.localize(trigger, row, window_mode="peak", artifacts=art_peak)

    os.makedirs(FIGDIR_ALL, exist_ok=True)
    figpath = os.path.join(FIGDIR_ALL, f"cmp_{trigger}.png")
    four_way_figure(trigger, art_t90, art_peak, ext, res_t90, res_peak, figpath)

    combined = {
        "trigger": trigger,
        "name": res_peak["name"],
        "t90": res_peak["t90"],
        "grb_type": res_peak["grb_type"],
        "fluence": res_peak.get("fluence"),
        "triggered_detectors": res_peak["triggered_detectors"],
        "dead_detectors": res_peak["dead_detectors"],
        "official_ra": res_peak["official_ra"],
        "official_dec": res_peak["official_dec"],
        "official_area90": res_peak["official_area90"],
        "catalog_ra": res_peak["catalog_ra"],
        "catalog_dec": res_peak["catalog_dec"],
        "catalog_err": res_peak["catalog_err"],
        "figure": os.path.basename(figpath),
        "t90_mode": {
            "sep": res_t90["sep_mine_official"],
            "ra": res_t90["my_ra"], "dec": res_t90["my_dec"],
            "area90": res_t90["my_area90"],
            "src_duration": res_t90["src_duration"],
            "official_in_mine": res_t90["official_in_my_map"],
        },
        "peak_mode": {
            "sep": res_peak["sep_mine_official"],
            "ra": res_peak["my_ra"], "dec": res_peak["my_dec"],
            "area90": res_peak["my_area90"],
            "src_duration": res_peak["src_duration"],
            "snr": res_peak["peak_snr"],
            "net_fraction": res_peak["peak_net_fraction"],
            "official_in_mine": res_peak["official_in_my_map"],
        },
    }
    if ext is not None:
        from astropy.coordinates import SkyCoord
        target = SkyCoord(ext["ra"], ext["dec"], unit="deg")
        combined["external"] = {
            "source": ext["source"], "name": ext["name"],
            "ra": ext["ra"], "dec": ext["dec"],
            "sep_official": float(
                SkyCoord(res_peak["official_ra"], res_peak["official_dec"],
                         unit="deg").separation(target).deg),
            "sep_t90": float(
                SkyCoord(res_t90["my_ra"], res_t90["my_dec"],
                         unit="deg").separation(target).deg),
            "sep_peak": float(
                SkyCoord(res_peak["my_ra"], res_peak["my_dec"],
                         unit="deg").separation(target).deg),
        }
    return combined


def worker(args):
    trigger, row, ext = args
    out = os.path.join(RESULTDIR, f"{trigger}.json")
    if os.path.exists(out):
        return trigger, "cached", None
    try:
        combined = process(trigger, row, ext)
        tmp = out + ".part"
        with open(tmp, "w") as handle:
            json.dump(combined, handle)
        os.replace(tmp, out)
        return trigger, "ok", combined["peak_mode"]["sep"]
    except Exception as err:
        return trigger, f"FAIL {type(err).__name__}: {err}", None


def load_external():
    """Best external position per burst, preferring Swift over LAT precision."""
    if not os.path.exists(EXTERNAL):
        return {}
    out = {}
    for row in csv.DictReader(open(EXTERNAL, newline="")):
        # Skip matches already identified as the wrong burst.
        if row.get("match_quality") == "suspect":
            continue
        ra = fnum(row["swift_ra"])
        dec = fnum(row["swift_dec"])
        source, name = "Swift", row["swift_name"]
        if ra is None or dec is None:
            ra, dec = fnum(row["lat_ra"]), fnum(row["lat_dec"])
            source, name = "Fermi-LAT", row["lat_name"]
        if ra is None or dec is None:
            continue
        out[row["trigger_name"]] = {
            "ra": ra, "dec": dec, "source": source, "name": name,
        }
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--class", dest="cls", choices=("all", "short", "long"),
                        default="all")
    args = parser.parse_args()

    os.makedirs(RESULTDIR, exist_ok=True)
    os.makedirs(FIGDIR_ALL, exist_ok=True)

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
        kind = "short" if t90 < SHORT_LONG_CUT else "long"
        if args.cls != "all" and kind != args.cls:
            continue
        eligible.append(trigger)
    eligible.sort()
    if args.limit:
        eligible = eligible[: args.limit]

    n_short = sum(1 for t in eligible if fnum(rows[t]["t90"]) < SHORT_LONG_CUT)
    print(f"{len(eligible)} bursts to process "
          f"({n_short} short, {len(eligible) - n_short} long), "
          f"{len(ext_all)} external positions available, "
          f"{args.workers} workers", flush=True)

    tasks = [(t, rows[t], ext_all.get(t)) for t in eligible]

    from concurrent.futures import ProcessPoolExecutor, as_completed
    done = ok = cached = failed = 0
    failures = []
    import time
    started = time.time()

    with ProcessPoolExecutor(args.workers) as pool:
        futures = {pool.submit(worker, task): task[0] for task in tasks}
        for future in as_completed(futures):
            trigger, status, _ = future.result()
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
                eta = (len(tasks) - done) / rate / 60 if rate else 0
                print(f"  {done}/{len(tasks)}  ok={ok} cached={cached} "
                      f"failed={failed}  {rate:.2f}/s  ETA {eta:.0f} min",
                      flush=True)

    print(f"\nfinished: {ok} new, {cached} cached, {failed} failed")
    for trigger, status in failures[:25]:
        print(f"  {trigger}: {status}")
    if failures:
        with open("site/all/failures.json", "w") as handle:
            json.dump(failures, handle, indent=2)


if __name__ == "__main__":
    main()
