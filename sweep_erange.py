"""
Sweep the localization energy band over a fixed burst sample.

Trigdat carries only 8 coarse channels, so the band is not continuously
tunable: each keV pair collapses onto a channel-index pair, and only a handful
of distinct bands exist. The pairs below were chosen to land on distinct
cranges (see the printed mapping).

Selecting a band by minimising separation on the same bursts used to measure it
is fitting to the answer, so treat the winner here as a candidate only and
confirm it on an independent sample with --seed.

Usage:
    python sweep_erange.py --seed 20261022 --n 10
    python sweep_erange.py --seed 7 --n 12 --tag holdout
"""
import argparse
import json
import os
import statistics
import subprocess
import sys

PY = sys.executable
SITEDIR = "site"

# (lo_keV, hi_keV, label) chosen so each maps to a distinct channel pair.
BANDS = [
    (50, 300, "44-300 keV  (standard, ch3-4)"),
    (50, 95, "44-95 keV   (soft only, ch3)"),
    (100, 300, "95-300 keV  (hard only, ch4)"),
    (25, 300, "22-300 keV  (ch2-4)"),
    (50, 400, "44-500 keV  (ch3-5)"),
    (25, 95, "22-95 keV   (ch2-3)"),
    (100, 400, "95-500 keV  (ch4-5)"),
    (25, 400, "22-500 keV  (ch2-5)"),
]


def run(lo, hi, seed, n, tag):
    out = f"sweep_{tag}_{lo}_{hi}.json"
    path = os.path.join(SITEDIR, out)
    if not os.path.exists(path):
        cmd = [PY, "analyze_grbs.py", "--n", str(n), "--seed", str(seed),
               "--window", "peak", "--erange", f"{lo},{hi}",
               "--noplots", "--out", out]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            print(f"  band {lo}-{hi} FAILED\n{proc.stdout[-1500:]}\n{proc.stderr[-1500:]}")
            return None
    with open(path) as handle:
        return json.load(handle)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=20261022)
    parser.add_argument("--n", type=int, default=10)
    parser.add_argument("--tag", type=str, default="main")
    args = parser.parse_args()

    print(f"sweep over {len(BANDS)} bands, seed={args.seed}, n={args.n}\n")
    header = (f"{'band':30s}{'crange':>9s}{'median':>9s}{'mean':>8s}"
              f"{'max':>8s}{'<2deg':>7s}{'area':>8s}")
    print(header)
    print("-" * len(header))

    table = []
    for lo, hi, label in BANDS:
        data = run(lo, hi, args.seed, args.n, args.tag)
        if not data or not data["results"]:
            continue
        res = data["results"]
        seps = [r["sep_mine_official"] for r in res]
        ratios = [r["my_area90"] / r["official_area90"] for r in res]
        crange = res[0]["crange"]
        row = {
            "lo": lo, "hi": hi, "label": label, "crange": crange,
            "n": len(res),
            "median": statistics.median(seps),
            "mean": statistics.mean(seps),
            "max": max(seps),
            "under2": sum(1 for s in seps if s < 2.0),
            "area": statistics.median(ratios),
        }
        table.append(row)
        print(f"{label:30s}{str(crange):>9s}{row['median']:9.2f}{row['mean']:8.2f}"
              f"{row['max']:8.2f}{row['under2']:5d}/{len(res):<2d}{row['area']:8.2f}")

    if table:
        best = min(table, key=lambda r: r["median"])
        base = next((r for r in table if (r["lo"], r["hi"]) == (50, 300)), None)
        print()
        print(f"best median: {best['label'].strip()}  ->  {best['median']:.2f} deg")
        if base:
            print(f"standard band:                 ->  {base['median']:.2f} deg")
            print(f"difference: {best['median'] - base['median']:+.2f} deg")

    with open(os.path.join(SITEDIR, f"sweep_{args.tag}_summary.json"), "w") as handle:
        json.dump({"seed": args.seed, "n": args.n, "bands": table}, handle, indent=2)
    print(f"\nwrote {SITEDIR}/sweep_{args.tag}_summary.json")


if __name__ == "__main__":
    main()
