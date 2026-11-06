"""
Where each detector was pointing: at the burst, at the Earth, or at open sky.

This is what makes the detector panel readable rather than merely true. You
can already see that n7 is the bright one; this says n7 is bright *because it
was the detector pointing at the burst*, 26 degrees off, while n2 at 157
degrees saw nothing.

Nothing here comes from TTE -- TTE carries no pointing information at all.
It comes from trigdat's position history and quaternion, the same source the
localization already uses for the spacecraft attitude, so no new download.

Geometry, all at the trigger instant:

  nadir       straight down, toward the Earth's centre. From Fermi's 536 km
              the Earth is a disc of angular radius ~67 deg centred there,
              covering about 30% of the sky. A direction less than that from
              nadir is looking at the planet.
  Sun         orbital night when the Sun itself is inside the Earth disc.
  burst       the published Fermi GBM position, so the angles describe the
              real sky rather than this reproduction's answer.

TWO TRAPS, both of which produced confidently wrong numbers before being
caught:

  1. `frame.geocenter` is a GCRS position. Calling `.separation()` on it
     against an ICRS coordinate makes astropy transform it mid-comparison and
     silently returns the wrong angle -- it reported 39% of real triggered
     bursts as being behind the Earth, which is impossible.
  2. Round-tripping it through `.icrs` is also wrong: that re-refers a
     position with finite distance to the solar system barycentre and flips
     the direction.

     Taking the raw .ra/.dec and treating them as ICRS agrees with gdt's own
     `location_visible()` 54 times out of 54. The `.icrs` version agrees 33.

Detector pointings are the opposite case: `detector.skycoord(frame)` really
is in the spacecraft frame and really does need `.icrs`. That path is
confirmed by the physics -- angle to burst against net counts gives a
within-burst Spearman of -0.68, negative in 98% of bursts.

Stored per burst:
    er    Earth angular radius, degrees
    bn    burst angle from nadir (so bn - er is its clearance past the limb)
    sn    Sun angle from nadir (less than er means orbital night)
    a     base64 uint8, 3 x 14 bytes: angle to burst, to nadir, to Sun,
          for each detector in the standard order, rounded to whole degrees

Resumable: a burst whose JSON exists is skipped.

Usage:
    python extract_geometry.py --workers 8
"""
import argparse
import base64
import glob
import json
import os
import warnings

warnings.filterwarnings("ignore")

import numpy as np

DATADIR = "data/gbm"
RESULTS = "site/all/results"
DETSRC = "site/vec/det"
OUTDIR = "site/vec/geom"

ALL_DETS = [f"n{d}" for d in "0123456789ab"] + ["b0", "b1"]

# Beyond the limb by this much is called open sky rather than limb-grazing.
LIMB_PAD = 20.0


def raw(coord):
    """Take a coordinate's angles at face value, as ICRS.

    Deliberately not `.icrs`: see the module docstring. This is the form
    validated against gdt's own location_visible().
    """
    from astropy.coordinates import SkyCoord
    return SkyCoord(float(np.ravel(coord.ra.deg)[0]),
                    float(np.ravel(coord.dec.deg)[0]), unit="deg")


def deg(q):
    import astropy.units as u
    return float(np.ravel(u.Quantity(q).to_value(u.deg))[0])


def process(trigger):
    from gdt.missions.fermi.gbm.trigdat import Trigdat
    from gdt.missions.fermi.gbm.detectors import GbmDetectors
    from astropy.coordinates import SkyCoord, get_sun
    from astropy.time import Time

    path = sorted(glob.glob(os.path.join(DATADIR, trigger, "glg_trigdat_*.fit")))
    if not path:
        raise ValueError("no trigdat")
    trigdat = Trigdat.open(path[0])
    when = Time(trigdat.trigtime, format="fermi")
    frame = trigdat.poshist.at(when)

    nadir = raw(frame.geocenter)
    sun = raw(get_sun(when))
    radius = deg(frame.earth_angular_radius)

    res_path = f"{RESULTS}/{trigger}.json"
    with open(res_path) as handle:
        res = json.load(handle)
    burst = SkyCoord(res["official_ra"], res["official_dec"], unit="deg")

    to_burst, to_nadir, to_sun = [], [], []
    for det in ALL_DETS:
        # Detector pointings DO need .icrs; the geocenter does not.
        point = raw(GbmDetectors.from_str(det).skycoord(frame).icrs)
        to_burst.append(deg(point.separation(burst)))
        to_nadir.append(deg(point.separation(nadir)))
        to_sun.append(deg(point.separation(sun)))

    angles = np.clip(np.rint(np.concatenate([to_burst, to_nadir, to_sun])),
                     0, 255).astype(np.uint8)
    return {
        "t": trigger,
        "er": round(radius, 1),
        "bn": round(deg(burst.separation(nadir)), 1),
        "sn": round(deg(sun.separation(nadir)), 1),
        "a": base64.b64encode(angles.tobytes()).decode("ascii"),
    }


def worker(trigger):
    out = os.path.join(OUTDIR, f"{trigger}.json")
    if os.path.exists(out):
        return trigger, "cached"
    try:
        rec = process(trigger)
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
    eligible = sorted(os.path.splitext(f)[0] for f in os.listdir(DETSRC)
                      if f.endswith(".json"))
    eligible = [t for t in eligible if os.path.exists(f"{RESULTS}/{t}.json")]
    if args.limit:
        eligible = eligible[: args.limit]
    print(f"{len(eligible)} bursts, {args.workers} workers", flush=True)

    from concurrent.futures import ProcessPoolExecutor, as_completed
    import time
    done = ok = cached = failed = 0
    failures = []
    started = time.time()
    with ProcessPoolExecutor(args.workers) as pool:
        futures = {pool.submit(worker, t): t for t in eligible}
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
            if done % 200 == 0 or done == len(eligible):
                rate = done / max(time.time() - started, 1e-9)
                print(f"  {done}/{len(eligible)}  ok={ok} failed={failed}  "
                      f"ETA {(len(eligible)-done)/rate/60:.0f} min", flush=True)

    print(f"\n{ok} new, {cached} cached, {failed} failed")
    files = [f for f in os.listdir(OUTDIR) if f.endswith(".json")]
    total = sum(os.path.getsize(os.path.join(OUTDIR, f)) for f in files)
    print(f"geometry: {total/1e6:.2f} MB across {len(files)} files "
          f"({total/max(len(files),1):.0f} bytes each)")
    for trigger, status in failures[:10]:
        print(f"  {trigger}: {status}")


if __name__ == "__main__":
    main()
