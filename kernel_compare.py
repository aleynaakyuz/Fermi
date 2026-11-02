"""
Which systematic kernel does the official GBM map look like it used?

Every official HEALPix file carries LMETHOD = 'Interactive', written
unconditionally, so the header cannot say whether a human or RoboBA produced
a given map. The two kernels differ enough in width that the published region
AREA may carry the answer:

    hitl_model     4.17 deg (91.8%) + 15.3 deg tail,
                   or 2.31 deg (88.4%) + 13.2 deg tail
    robo_ba_model  long  1.86 deg (57.9%) + 4.14 deg
                   short 2.55 deg (39.0%) + 4.43 deg

RoboBA's kernel has no 13-15 deg tail at all, so a RoboBA-produced map should
be markedly tighter.

Method. The SAME statistical map is convolved both ways, so the ratio between
my two areas is purely the kernel. For each burst, score each kernel by how
close it brings my area to the official one, in log space because areas span
orders of magnitude:

    score = |log(my_area / official_area)|

The confound, stated plainly. My statistical map is not GBM's -- different
source intervals give a different width before any kernel is applied. Over the
catalog my hitl areas run a median 1.30x the official ones, so some of any
apparent preference is that offset rather than the kernel. Two guards:

  - the per-burst winner is also reported after dividing out the median
    offset, so a uniform bias cannot by itself decide the outcome;
  - the result is broken down by year. RoboBA came into service around
    2017-2018. If the preference is real it should shift with epoch; if it is
    an artefact of my intervals it should not.

Usage:
    python kernel_compare.py
"""
import glob
import json
import collections

import numpy as np

SRC = "site/vec/bursts"


def main():
    rows = []
    for path in sorted(glob.glob(f"{SRC}/*.json")):
        with open(path) as handle:
            d = json.load(handle)
        oa, pa, ra = d.get("oa"), d.get("pa"), d.get("ra_area")
        if not all(isinstance(v, (int, float)) and v > 0 for v in (oa, pa, ra)):
            continue
        rows.append({
            "t": d["t"], "year": "20" + d["t"][2:4], "k": d["k"],
            "oa": float(oa), "hitl": float(pa), "robo": float(ra),
            "ps": d.get("ps"),
        })
    if not rows:
        raise SystemExit("no contour data with both kernel areas yet")

    hitl = np.array([r["hitl"] for r in rows])
    robo = np.array([r["robo"] for r in rows])
    off = np.array([r["oa"] for r in rows])
    print(f"{len(rows)} bursts with both kernels and an official area\n")

    print("Area relative to the official map:")
    for lab, arr in (("hitl", hitl), ("RoboBA", robo)):
        rt = arr / off
        print(f"  {lab:7s} median ratio {np.median(rt):5.2f}   "
              f"within 2x of official: {100*((rt > 0.5) & (rt < 2)).mean():4.0f}%")
    print(f"  RoboBA area / hitl area: median {np.median(robo/hitl):.2f} "
          f"(pure kernel effect, same statistical map)")
    print()

    def tally(scale_h=1.0, scale_r=1.0, label=""):
        sh = np.abs(np.log(hitl / scale_h / off))
        sr = np.abs(np.log(robo / scale_r / off))
        win_h = int((sh < sr).sum())
        print(f"  {label:34s} hitl {win_h:5d}  RoboBA {len(rows)-win_h:5d}  "
              f"({100*(len(rows)-win_h)/len(rows):.0f}% RoboBA)")
        return sh, sr

    print("Which kernel brings my area closer to the official one?")
    tally(1.0, 1.0, "raw")
    # Divide out each kernel's own median offset, so a uniform bias in my
    # statistical map cannot decide the winner on its own.
    tally(np.median(hitl / off), np.median(robo / off), "after removing median offset")
    print()

    print("By year (RoboBA entered service around 2017-2018):")
    print(f"  {'year':6s}{'n':>6s}{'% favouring RoboBA':>21s}{'median off area':>18s}")
    by = collections.defaultdict(list)
    for i, r in enumerate(rows):
        by[r["year"]].append(i)
    mh, mr = np.median(hitl / off), np.median(robo / off)
    for year in sorted(by):
        idx = np.array(by[year])
        if len(idx) < 10:
            continue
        sh = np.abs(np.log(hitl[idx] / mh / off[idx]))
        sr = np.abs(np.log(robo[idx] / mr / off[idx]))
        print(f"  {year:6s}{len(idx):6d}{100*(sr < sh).mean():20.0f}%"
              f"{np.median(off[idx]):18.0f}")
    print()

    print("By burst class:")
    for cls in ("long", "short"):
        idx = np.array([i for i, r in enumerate(rows) if r["k"] == cls])
        if not len(idx):
            continue
        sh = np.abs(np.log(hitl[idx] / mh / off[idx]))
        sr = np.abs(np.log(robo[idx] / mr / off[idx]))
        print(f"  {cls:6s} n={len(idx):5d}  favouring RoboBA {100*(sr < sh).mean():4.0f}%"
              f"   median official area {np.median(off[idx]):6.0f} sq deg")


if __name__ == "__main__":
    main()
