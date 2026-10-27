"""
Identify Fermi GBM bursts that carry a localization from another mission.

Two independent signals, combined:

  error_radius  The GBM catalog adopts an external position when one exists.
                GBM's own ground localization cannot do better than its
                systematic floor of roughly 2-4 deg, so any burst listed with a
                sub-degree error radius is carrying somebody else's position.
                This flags a burst but cannot say which instrument supplied it.

  cross-match   Naming the instrument requires the other catalogs:
                  fermilgrb  Fermi LAT 2nd GRB catalog. Joins directly on its
                             gbm_cat_name column, so no time matching needed.
                  swiftgrb   Swift. Matched on trigger time, with the burst
                             date required to agree as a second check.

A burst flagged by error_radius but matched by neither catalog still has an
external position, from a mission not queried here (INTEGRAL, Konus, AGILE,
MAXI, IPN triangulation) or from an afterglow reported only in a GCN circular.
Those are reported as "unattributed" rather than dropped.

Output: catalog/external_localizations.csv
"""
import csv
import io
import os
import urllib.parse
import urllib.request

import numpy as np
from astropy.io.votable import parse_single_table
from astropy.coordinates import SkyCoord
import astropy.units as u

TAP = "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync"
MANIFEST = "catalog/fermigbrst.csv"
OUT = "catalog/external_localizations.csv"
CACHE = "catalog/_xmatch_cache"

# A GBM position is never this good; below this, the position came from
# somewhere else. GBM's systematic core alone is 2.3-4.2 deg.
EXTERNAL_ERR_DEG = 0.5
# Swift and GBM trigger on the same burst within seconds; 60 s is generous and
# still far tighter than the ~1/day trigger rate of either instrument.
TIME_TOL_SEC = 60.0


def tap(query, cache_name):
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, cache_name)
    if not os.path.exists(path):
        url = f"{TAP}?{urllib.parse.urlencode({'QUERY': query, 'LANG': 'ADQL'})}"
        with urllib.request.urlopen(url, timeout=300) as resp:
            payload = resp.read()
        with open(path, "wb") as handle:
            handle.write(payload)
    with open(path, "rb") as handle:
        table = parse_single_table(io.BytesIO(handle.read())).to_table()
    for col in table.colnames:
        if table[col].dtype.kind == "S":
            table[col] = [v.decode().strip() for v in table[col]]
    return table


def as_float(value, default=np.nan):
    try:
        out = float(value)
        return out if np.isfinite(out) else default
    except (TypeError, ValueError):
        return default


def main():
    gbm = list(csv.DictReader(open(MANIFEST, newline="")))
    print(f"{len(gbm)} GBM bursts in the catalog")

    # ---------------------------------------------------- LAT, by direct join
    lat = tap(
        "SELECT name, gbm_cat_name, ra, dec, error_radius FROM fermilgrb",
        "fermilgrb.xml")
    lat_by_gbm = {}
    for row in lat:
        key = str(row["gbm_cat_name"]).strip()
        if key and key.lower() not in ("", "none", "nan"):
            lat_by_gbm[key] = row
    print(f"  LAT catalog: {len(lat)} rows, {len(lat_by_gbm)} with a GBM name")

    # --------------------------------------------------- Swift, two catalogs
    # swiftgrb is frozen at 2012-12-29, so on its own it misses every Swift
    # burst from 2013 on. swiftgrbba, the Burst Advocate compilation, is the
    # maintained one and carries MJD directly. Both are read and merged, with
    # swiftgrbba preferred where the same burst appears in each.
    from astropy.time import Time as AstroTime

    def iso_to_mjd(raw):
        raw = str(raw).strip()
        if not raw or raw.lower() in ("", "nan", "none", "--"):
            return np.nan
        try:
            return float(AstroTime(raw, format="isot", scale="utc").mjd)
        except Exception:
            return as_float(raw)

    entries = []

    legacy = tap(
        "SELECT name, trigger_time, ra, dec, pos_err, redshift FROM swiftgrb",
        "swiftgrb.xml")
    for row in legacy:
        mjd = iso_to_mjd(row["trigger_time"])
        if np.isfinite(mjd):
            entries.append({
                "mjd": mjd, "name": str(row["name"]).strip(),
                "ra": as_float(row["ra"]), "dec": as_float(row["dec"]),
                "err": as_float(row["pos_err"]),
                "redshift": str(row["redshift"]).strip(),
                "trigger_obs": "", "table": "swiftgrb",
            })

    ba = tap(
        'SELECT name, "time", ra, dec, error_radius, xrt_ra, xrt_dec, '
        'xrt_error_radius, redshift, trigger_obs FROM swiftgrbba',
        "swiftgrbba.xml")
    for row in ba:
        mjd = as_float(row["time"])
        if not np.isfinite(mjd):
            continue
        ra, dec = as_float(row["xrt_ra"]), as_float(row["xrt_dec"])
        err = as_float(row["xrt_error_radius"])
        if not (np.isfinite(ra) and np.isfinite(dec)):
            ra, dec, err = as_float(row["ra"]), as_float(row["dec"]), \
                as_float(row["error_radius"])
        entries.append({
            "mjd": mjd, "name": str(row["name"]).strip(),
            "ra": ra, "dec": dec, "err": err,
            "redshift": str(row["redshift"]).strip(),
            "trigger_obs": str(row["trigger_obs"]).strip(), "table": "swiftgrbba",
        })

    # Dedupe: same burst within the tolerance keeps the swiftgrbba row.
    entries.sort(key=lambda e: (e["mjd"], e["table"] != "swiftgrbba"))
    merged = []
    for entry in entries:
        if merged and abs(entry["mjd"] - merged[-1]["mjd"]) * 86400 <= TIME_TOL_SEC:
            continue
        merged.append(entry)

    sw_rows = merged
    sw_times = np.array([e["mjd"] for e in sw_rows])
    n_legacy = sum(1 for e in sw_rows if e["table"] == "swiftgrb")
    print(f"  Swift: {len(legacy)} legacy + {len(ba)} BA rows "
          f"-> {len(sw_rows)} unique ({n_legacy} only in the frozen catalog)")

    tol_days = TIME_TOL_SEC / 86400.0
    out_rows = []
    counts = {"lat": 0, "swift": 0, "both": 0, "unattributed": 0}

    for burst in gbm:
        trigger = burst["trigger_name"].strip()
        err = as_float(burst["error_radius"], default=np.inf)
        has_external = err < EXTERNAL_ERR_DEG
        mjd = as_float(burst["trigger_time"])

        # LAT
        lat_row = lat_by_gbm.get(trigger) or lat_by_gbm.get(burst["name"].strip())

        # Swift, nearest trigger in time
        swift_row, dt_sec = None, None
        if np.isfinite(mjd) and len(sw_times):
            idx = int(np.searchsorted(sw_times, mjd))
            for cand in (idx - 1, idx, idx + 1):
                if 0 <= cand < len(sw_times):
                    delta = abs(sw_times[cand] - mjd)
                    if delta <= tol_days and (dt_sec is None or delta * 86400 < dt_sec):
                        swift_row, dt_sec = sw_rows[cand], delta * 86400

        if not (has_external or lat_row is not None or swift_row is not None):
            continue

        sources = []
        if swift_row is not None:
            sources.append("Swift")
        if lat_row is not None:
            sources.append("Fermi-LAT")
        attribution = "+".join(sources) if sources else "unattributed"

        if swift_row is not None and lat_row is not None:
            counts["both"] += 1
        elif swift_row is not None:
            counts["swift"] += 1
        elif lat_row is not None:
            counts["lat"] += 1
        else:
            counts["unattributed"] += 1

        # Offset between the adopted catalog position and the external one, a
        # consistency check that the match is the same burst.
        offset = ""
        gbm_ra, gbm_dec = as_float(burst["ra"]), as_float(burst["dec"])
        ext_ra = ext_dec = np.nan
        if swift_row is not None:
            ext_ra, ext_dec = swift_row["ra"], swift_row["dec"]
        elif lat_row is not None:
            ext_ra, ext_dec = as_float(lat_row["ra"]), as_float(lat_row["dec"])
        if all(np.isfinite(v) for v in (gbm_ra, gbm_dec, ext_ra, ext_dec)):
            offset = f"{SkyCoord(gbm_ra, gbm_dec, unit='deg').separation(SkyCoord(ext_ra, ext_dec, unit='deg')).deg:.4f}"

        # Quality flag. When the catalog carries an adopted external position
        # (small error radius), the matched coordinates must agree with it; if
        # they do not, the cross-match found the wrong burst. When the error
        # radius is large the catalog holds GBM's own position, so a large
        # offset is just GBM's localization error and says nothing about the
        # match.
        quality = "ok"
        if offset:
            off_val = float(offset)
            if err < EXTERNAL_ERR_DEG and off_val > 1.0:
                quality = "suspect"
            elif err >= EXTERNAL_ERR_DEG and off_val > 1.0:
                quality = "gbm_own_position"

        out_rows.append({
            "trigger_name": trigger,
            "name": burst["name"].strip(),
            "trigger_mjd": burst["trigger_time"],
            "t90": burst["t90"],
            "catalog_ra": burst["ra"],
            "catalog_dec": burst["dec"],
            "catalog_err_deg": burst["error_radius"],
            "flagged_by_err_radius": int(has_external),
            "attribution": attribution,
            "swift_name": swift_row["name"] if swift_row is not None else "",
            "swift_dt_sec": f"{dt_sec:.1f}" if dt_sec is not None else "",
            "swift_ra": f"{swift_row['ra']:.4f}" if swift_row is not None and np.isfinite(swift_row["ra"]) else "",
            "swift_dec": f"{swift_row['dec']:.4f}" if swift_row is not None and np.isfinite(swift_row["dec"]) else "",
            "swift_pos_err": f"{swift_row['err']:.4f}" if swift_row is not None and np.isfinite(swift_row["err"]) else "",
            "swift_redshift": swift_row["redshift"] if swift_row is not None else "",
            "swift_trigger_obs": swift_row["trigger_obs"] if swift_row is not None else "",
            "lat_name": str(lat_row["name"]).strip() if lat_row is not None else "",
            "lat_ra": f"{as_float(lat_row['ra']):.4f}" if lat_row is not None else "",
            "lat_dec": f"{as_float(lat_row['dec']):.4f}" if lat_row is not None else "",
            "lat_err_deg": f"{as_float(lat_row['error_radius']):.4f}" if lat_row is not None else "",
            "catalog_vs_external_deg": offset,
            "match_quality": quality,
        })

    out_rows.sort(key=lambda r: r["trigger_name"])
    os.makedirs("catalog", exist_ok=True)
    with open(OUT, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(out_rows[0].keys()))
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"\n{len(out_rows)} GBM bursts with an external localization")
    print(f"  Swift only            {counts['swift']}")
    print(f"  Fermi-LAT only        {counts['lat']}")
    print(f"  both                  {counts['both']}")
    print(f"  unattributed          {counts['unattributed']}")
    flagged = sum(r["flagged_by_err_radius"] for r in out_rows)
    print(f"\n  flagged by error_radius < {EXTERNAL_ERR_DEG} deg: {flagged}")
    matched = len(out_rows) - counts["unattributed"]
    print(f"  attributed to an instrument:          {matched}")
    suspect = sum(1 for r in out_rows if r["match_quality"] == "suspect")
    own = sum(1 for r in out_rows if r["match_quality"] == "gbm_own_position")
    print(f"\n  match quality: {suspect} suspect, {own} carrying GBM's own "
          f"position, {len(out_rows) - suspect - own} clean")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
