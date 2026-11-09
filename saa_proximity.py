"""
Is the livetime instability flagged on seven bursts caused by the SAA?

Seven of 2,049 bursts carry a livetime-swing flag: their deadtime-corrected
exposure varies by more than 20% across the trigdat span, which breaks the
polynomial background the localization depends on. Biltzinger et al. (2020)
model exactly this kind of background as activation decaying after each SAA
exit, so the obvious question is whether these seven sit soon after one.

There is a confound that has to be handled first. A very bright burst causes
its own deadtime, and six of the seven are above the 98th percentile in
count rate -- so "bright" and "flagged" are nearly the same set, and any
SAA signal would be invisible against it. The test therefore uses two
control groups:

    flagged       the seven
    rate_matched  equally bright bursts whose livetime is steady. If the SAA
                  is responsible, these should sit FURTHER from an SAA exit.
    random        ordinary bursts, for the baseline distribution

Fermi is inside the SAA about 13% of the time and the orbit is ~96 minutes,
so time-since-exit is bounded by roughly an hour and a half. A null result
looks like all three groups drawing from the same spread of that interval.

Needs the daily poshist files, which carry a once-a-second SAA flag;
trigdat has no SAA information at all.

Usage:
    python saa_proximity.py --sample sample.json
"""
import argparse
import glob
import json
import os
import warnings

warnings.filterwarnings("ignore")

import numpy as np

POSDIR = "data/poshist"
DATADIR = "data/gbm"


def trigtime_of(trigger):
    from gdt.missions.fermi.gbm.trigdat import Trigdat
    paths = sorted(glob.glob(os.path.join(DATADIR, trigger, "glg_trigdat_*.fit")))
    if not paths:
        return None
    return float(Trigdat.open(paths[0]).trigtime)


def saa_intervals(yymmdd, cache={}):
    """(entry, exit) mission-elapsed times for every SAA passage that day.

    Read straight from the FITS columns rather than through
    GbmPosHist.get_spacecraft_states(). The states table is fine, but its
    `time` column is 86,520 astropy Time objects per day and pulling a MET
    out of each one costs more than the rest of this script put together.
    SCLK_UTC is already MET, and bit 1 of FLAGS is the SAA flag -- checked
    against gdt's own decoding, 86,520 of 86,520 rows agree.
    """
    if yymmdd in cache:
        return cache[yymmdd]
    from astropy.io import fits
    paths = sorted(glob.glob(os.path.join(POSDIR, f"glg_poshist_all_{yymmdd}_v*.fit")))
    if not paths:
        cache[yymmdd] = None
        return None
    with fits.open(paths[-1]) as handle:
        data = handle[1].data
        met = np.asarray(data["SCLK_UTC"], dtype=float)
        saa = ((np.asarray(data["FLAGS"]).astype(np.uint16) >> 1) & 1).astype(bool)
    # Rising and falling edges of the flag.
    edge = np.diff(saa.astype(np.int8))
    entries = met[1:][edge == 1]
    exits = met[1:][edge == -1]
    cache[yymmdd] = (met, saa, entries, exits)
    return cache[yymmdd]


def measure(trigger):
    t0 = trigtime_of(trigger)
    if t0 is None:
        return None
    day = trigger[2:8]
    got = saa_intervals(day)
    if got is None:
        return None
    met, saa, entries, exits = got

    # The burst may fall near midnight; the day file still brackets it, but
    # guard against a trigger outside the file's span.
    if not (met[0] <= t0 <= met[-1]):
        return None

    in_saa = bool(saa[np.argmin(np.abs(met - t0))])
    prev_exits = exits[exits <= t0]
    next_entry = entries[entries > t0]

    # A burst early in the UTC day can have its last SAA exit on the previous
    # one. Dropping those silently would have been the wrong kind of missing
    # data: they are precisely the cases with the LONGEST time since an exit,
    # so losing them drags the control medians down and flatters whatever the
    # flagged group does. Reach back a day instead.
    spans_midnight = False
    if not prev_exits.size:
        import datetime
        prev_day = (datetime.datetime.strptime(day, "%y%m%d")
                    - datetime.timedelta(days=1)).strftime("%y%m%d")
        earlier = saa_intervals(prev_day)
        if earlier is not None:
            _, _, _, prev_exits_y = earlier
            if prev_exits_y.size:
                prev_exits = prev_exits_y
                spans_midnight = True

    return {
        "trigger": trigger,
        "in_saa": in_saa,
        "since_exit": float(t0 - prev_exits[-1]) if prev_exits.size else None,
        "until_entry": float(next_entry[0] - t0) if next_entry.size else None,
        "n_passages": int(exits.size),
        "spans_midnight": spans_midnight,
    }


def summarise(name, rows):
    vals = np.array([r["since_exit"] / 60.0 for r in rows
                     if r.get("since_exit") is not None])
    if not vals.size:
        print(f"  {name:14s} no measurements")
        return vals
    print(f"  {name:14s} n={len(vals):3d}   median {np.median(vals):6.1f} min   "
          f"min {vals.min():5.1f}   max {vals.max():6.1f}   "
          f"within 30 min: {100*(vals < 30).mean():3.0f}%")
    return vals


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    args = parser.parse_args()

    with open(args.sample) as handle:
        sample = json.load(handle)

    results, missing = {}, []
    for group, rows in sample.items():
        out = []
        for row in rows:
            got = measure(row["t"])
            if got is None:
                missing.append(row["t"])
                continue
            got.update(swing=row["swing"], rate=row["rate"], group=group)
            out.append(got)
        results[group] = out

    print("Time since Fermi last left the SAA, at the trigger instant\n")
    groups = {}
    for group in ("flagged", "rate_matched", "random"):
        if group in results:
            groups[group] = summarise(group, results[group])

    print("\nThe seven flagged bursts individually:")
    print(f"  {'trigger':14s}{'swing':>7}{'rate c/s':>11}{'since exit':>12}{'in SAA':>8}")
    for r in sorted(results.get("flagged", []), key=lambda x: -x["swing"]):
        se = r["since_exit"]
        print(f"  {r['trigger']:14s}{r['swing']:7.2f}{r['rate']:11,.0f}"
              f"{(f'{se/60:.1f} min' if se is not None else '--'):>12}"
              f"{str(r['in_saa']):>8}")

    if "flagged" in groups and "rate_matched" in groups:
        a, b = groups["flagged"], groups["rate_matched"]
        if a.size and b.size:
            from scipy.stats import mannwhitneyu
            p = mannwhitneyu(a, b, alternative="less").pvalue
            print(f"\n  flagged vs rate-matched, one-sided Mann-Whitney")
            print(f"  (is the flagged group CLOSER to an SAA exit?)  p = {p:.3f}")
    if missing:
        print(f"\n  {len(missing)} bursts had no usable poshist: "
              f"{', '.join(missing[:8])}")

    out = "site/saa_proximity.json"
    with open(out, "w") as handle:
        json.dump(results, handle, indent=1)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
