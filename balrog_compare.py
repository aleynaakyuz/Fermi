"""
Score BALROG against the same independent positions used for everything else.

BALROG is a genuinely different method, not another error model: it fits sky
position and spectrum together, specifically to remove the systematics DoL
leaves behind. So for the bursts with an arcsecond position from Swift or the
LAT, four answers can be put side by side against the same truth:

    official   the published Fermi GBM map (DoL + a human's intervals)
    T90        this reproduction, full T90 window
    best S/N   this reproduction, brightest sub-interval
    BALROG     MPE's independent Bayesian fit

Two format differences from the Fermi maps are handled here: BALROG publishes
RING ordering at NSIDE 32 (~1.8 deg pixels) where Fermi uses NESTED at
NSIDE 128 (~0.46 deg). Centroids are therefore quantized about four times more
coarsely, which sets a floor of roughly a degree on any BALROG separation --
worth remembering before reading small differences as real.

Usage:
    python balrog_compare.py
"""
import glob
import json
import os

import numpy as np
import healpy as hp
from astropy.io import fits
from astropy.coordinates import SkyCoord

BALDIR = "data/balrog"
RESULTS = "site/all/results"
OUT = "site/balrog_comparison.json"


def read_balrog(path):
    """Centroid and 90% area from a BALROG HEALPix map."""
    with fits.open(path) as hdul:
        hdu = hdul[1]
        prob = np.asarray(hdu.data[hdu.columns[0].name], dtype=float).ravel()
        head = hdu.header
        nside = int(head.get("NSIDE", hp.npix2nside(prob.size)))
        nest = str(head.get("ORDERING", "RING")).strip().upper() == "NESTED"
    total = prob.sum()
    if not np.isfinite(total) or total <= 0:
        return None
    prob = prob / total
    theta, phi = hp.pix2ang(nside, int(np.argmax(prob)), nest=nest)
    lon = float(np.degrees(phi))
    lat = float(90.0 - np.degrees(theta))

    # The catalogue mixes frames: the SysErr maps are Galactic at NSIDE 512,
    # the plain ones celestial at NSIDE 32 or 256. Reading a Galactic map as
    # RA/Dec yields a position tens of degrees wrong and nothing complains,
    # so the frame is honoured explicitly.
    frame = str(head.get("COORDSYS", "C")).strip().upper()
    if frame.startswith("G"):
        c = SkyCoord(lon, lat, unit="deg", frame="galactic").icrs
        ra, dec = float(c.ra.deg), float(c.dec.deg)
    else:
        ra, dec = lon, lat
    order = np.argsort(prob)[::-1]
    inside = int(np.searchsorted(np.cumsum(prob[order]), 0.9) + 1)
    area = inside * hp.nside2pixarea(nside, degrees=True)
    return {"ra": ra, "dec": dec, "area90": float(area),
            "nside": nside, "nest": nest, "frame": frame}


def main():
    rows = []
    for d in sorted(glob.glob(f"{BALDIR}/*")):
        trigger = os.path.basename(d)
        maps = glob.glob(os.path.join(d, "*_healpix_*.fits"))
        res_path = f"{RESULTS}/{trigger}.json"
        if not maps or not os.path.exists(res_path):
            continue
        with open(res_path) as handle:
            prev = json.load(handle)
        if "external" not in prev:
            continue
        bal = read_balrog(sorted(maps)[-1])
        if bal is None:
            continue
        meta = {}
        mp = os.path.join(d, "balrog.json")
        if os.path.exists(mp):
            with open(mp) as handle:
                meta = json.load(handle)

        truth = SkyCoord(prev["external"]["ra"], prev["external"]["dec"], unit="deg")
        bpos = SkyCoord(bal["ra"], bal["dec"], unit="deg")
        row = {
            "trigger": trigger, "k": prev["grb_type"], "t90": prev["t90"],
            "source": prev["external"]["source"],
            "balrog_ra": bal["ra"], "balrog_dec": bal["dec"],
            "balrog_area": bal["area90"], "nside": bal["nside"],
            "frame": bal["frame"],
            "d_balrog": float(bpos.separation(truth).deg),
            "d_official": prev["external"]["sep_official"],
            "d_t90": prev["external"]["sep_t90"],
            "d_peak": prev["external"]["sep_peak"],
            "official_area": prev["official_area90"],
            "peak_area": prev["peak_mode"]["area90"],
        }
        if "ra" in meta and "dec" in meta:
            q = SkyCoord(meta["ra"], meta["dec"], unit="deg")
            row["quoted_ra"] = meta["ra"]
            row["quoted_dec"] = meta["dec"]
            row["d_quoted"] = float(q.separation(truth).deg)
            row["map_vs_quoted"] = float(q.separation(bpos).deg)
            row["sys_err"] = meta.get("sys_err")
        rows.append(row)

    if not rows:
        raise SystemExit("no BALROG maps paired with a truth position yet")

    with open(OUT, "w") as handle:
        json.dump(rows, handle, indent=2)

    n = len(rows)
    print(f"{n} bursts with BALROG and an independent position\n")

    cols = [("official", "d_official"), ("T90 window", "d_t90"),
            ("best S/N", "d_peak"), ("BALROG", "d_balrog")]
    print("Distance from the independent position:")
    print(f"  {'method':14s}{'median':>9s}{'mean':>8s}{'<5 deg':>9s}{'>10 deg':>9s}")
    arrs = {}
    for lab, key in cols:
        a = np.array([r[key] for r in rows], dtype=float)
        arrs[lab] = a
        print(f"  {lab:14s}{np.median(a):9.2f}{a.mean():8.2f}"
              f"{100*(a<5).mean():8.0f}%{100*(a>10).mean():8.0f}%")
    print()

    b = arrs["BALROG"]
    for lab in ("official", "best S/N"):
        a = arrs[lab]
        print(f"  BALROG closer than {lab:10s}: {(b<a).sum():3d}/{n} "
              f"({100*(b<a).mean():3.0f}%)   median paired change {np.median(b-a):+.2f} deg")
    try:
        from scipy.stats import wilcoxon
        print(f"  Wilcoxon BALROG vs official : p = {wilcoxon(b, arrs['official']).pvalue:.4f}")
        print(f"  Wilcoxon BALROG vs best S/N : p = {wilcoxon(b, arrs['best S/N']).pvalue:.4f}")
    except Exception:
        pass
    print()

    if any("d_quoted" in r for r in rows):
        q = np.array([r["d_quoted"] for r in rows if "d_quoted" in r])
        m = np.array([r["map_vs_quoted"] for r in rows if "map_vs_quoted" in r])
        print(f"BALROG's own quoted position vs its map's peak pixel: "
              f"median {np.median(m):.2f} deg apart")
        print(f"  quoted position vs truth: median {np.median(q):.2f} deg "
              f"({len(q)} bursts)")
        print(f"  (map centroid is quantized at NSIDE {rows[0]['nside']}, "
              f"{hp.nside2resol(rows[0]['nside'], arcmin=True)/60:.2f} deg)")
    print()

    print("By burst class:")
    for cls in ("long", "short"):
        idx = [i for i, r in enumerate(rows) if r["k"] == cls]
        if not idx:
            continue
        print(f"  {cls:6s} n={len(idx):4d}   " + "   ".join(
            f"{lab} {np.median(arrs[lab][idx]):5.2f}" for lab, _ in cols))
    print()
    print("Region size (median 90% area, sq deg):")
    print(f"  official {np.median([r['official_area'] for r in rows]):6.0f}   "
          f"best S/N {np.median([r['peak_area'] for r in rows]):6.0f}   "
          f"BALROG {np.median([r['balrog_area'] for r in rows]):6.0f}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
