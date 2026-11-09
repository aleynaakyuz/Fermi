"""
Download continuous daily GBM data: CTIME, CSPEC and position history.

Everything else in this project is trigger-centred -- trigdat and TTE both
cover a few hundred seconds around a burst. That is the wrong shape for
background work. SAA activation decays over hours, so inside a trigger window
it is an almost constant offset, absorbed whole by the background
polynomial's constant term. Measuring it needs continuous days.

This is also the data gbmbkgpy (Biltzinger et al. 2020) expects: it ingests
daily CTIME or CSPEC plus poshist, never trigdat.

    https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/YYYY/MM/DD/current/
        glg_ctime_<det>_YYMMDD_vNN.pha    8 channels, 0.256 s
        glg_cspec_<det>_YYMMDD_vNN.pha  128 channels, 4.096 s
        glg_poshist_all_YYMMDD_vNN.fit   position, attitude, SAA flag, 1 s

Roughly 286 MB per day for all 14 detectors in both binnings plus poshist,
so a three-week pilot is about 6 GB and a year would be 105 GB.

Resumable: a file already on disk and intact is skipped.

Usage:
    python download_daily.py --start 2022-08-17 --days 21
    python download_daily.py --start 2022-08-17 --days 21 --types ctime,poshist
"""
import argparse
import datetime
import glob
import os
import re
import sys
import threading
import time
import urllib.error
import urllib.request

BASE = "https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily"
OUTDIR = "data/daily"

DETECTORS = [f"n{d}" for d in "0123456789ab"] + ["b0", "b1"]
USER_AGENT = "fermi-grb-compare/1.0 (academic use)"
RETRIES = 4
BACKOFF = 2.0
TIMEOUT = 600

_lock = threading.Lock()


def log(msg):
    with _lock:
        print(msg, flush=True)


def fetch(url):
    last = None
    for attempt in range(RETRIES):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
                return resp.read()
        except urllib.error.HTTPError as err:
            if err.code in (404, 410):
                return None
            last = err
        except Exception as err:
            last = err
        if attempt < RETRIES - 1:
            time.sleep(BACKOFF * (2 ** attempt))
    raise RuntimeError(f"{url}: {last}")


def is_intact(path):
    try:
        if os.path.getsize(path) < 2880:
            return False
        with open(path, "rb") as handle:
            return handle.read(6) == b"SIMPLE"
    except OSError:
        return False


def wanted_names(kind, yymmdd):
    """Filename patterns for a day, before version resolution."""
    if kind == "poshist":
        return [(f"glg_poshist_all_{yymmdd}", "fit")]
    return [(f"glg_{kind}_{det}_{yymmdd}", "pha") for det in DETECTORS]


def do_day(args):
    date, kinds = args
    yyyy, mm, dd = date.strftime("%Y"), date.strftime("%m"), date.strftime("%d")
    yymmdd = date.strftime("%y%m%d")
    dest = os.path.join(OUTDIR, yymmdd)
    rec = {"day": yymmdd, "ok": 0, "cached": 0, "missing": 0, "bytes": 0}

    stems = []
    for kind in kinds:
        stems += wanted_names(kind, yymmdd)

    # One directory listing serves every file for the day.
    need = []
    for stem, ext in stems:
        have = glob.glob(os.path.join(dest, f"{stem}_v*.{ext}"))
        if have and is_intact(have[0]):
            rec["cached"] += 1
        else:
            need.append((stem, ext))
    if not need:
        rec["status"] = "cached"
        return rec

    listing = fetch(f"{BASE}/{yyyy}/{mm}/{dd}/current/")
    if listing is None:
        rec["status"] = "no_directory"
        return rec
    text = listing.decode("utf-8", "replace")

    os.makedirs(dest, exist_ok=True)
    for stem, ext in need:
        names = sorted(set(re.findall(rf"{stem}_v(\d+)\.{ext}", text)))
        if not names:
            rec["missing"] += 1
            continue
        name = f"{stem}_v{names[-1]}.{ext}"
        payload = fetch(f"{BASE}/{yyyy}/{mm}/{dd}/current/{name}")
        if payload is None:
            rec["missing"] += 1
            continue
        tmp = os.path.join(dest, name + ".part")
        with open(tmp, "wb") as handle:
            handle.write(payload)
        os.replace(tmp, os.path.join(dest, name))
        rec["ok"] += 1
        rec["bytes"] += len(payload)

    rec["status"] = "ok" if rec["ok"] else "incomplete"
    return rec


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", required=True, help="YYYY-MM-DD")
    parser.add_argument("--days", type=int, default=21)
    parser.add_argument("--types", default="ctime,cspec,poshist")
    parser.add_argument("--workers", type=int, default=5)
    args = parser.parse_args()

    kinds = [k.strip() for k in args.types.split(",") if k.strip()]
    for kind in kinds:
        if kind not in ("ctime", "cspec", "poshist"):
            raise SystemExit(f"unknown type {kind!r}")

    start = datetime.datetime.strptime(args.start, "%Y-%m-%d").date()
    days = [start + datetime.timedelta(days=i) for i in range(args.days)]

    per_day = (190 if "ctime" in kinds else 0) + \
              (87 if "cspec" in kinds else 0) + (9 if "poshist" in kinds else 0)
    free = os.statvfs(".").f_bavail * os.statvfs(".").f_frsize / 1e9
    est = per_day * len(days) / 1000
    log(f"{len(days)} days from {days[0]} to {days[-1]}, types {','.join(kinds)}")
    log(f"  ~{est:.1f} GB estimated, {free:.0f} GB free, {args.workers} workers")
    if est > free * 0.8:
        log("  refusing: that would use more than 80% of free disk")
        return 1

    os.makedirs(OUTDIR, exist_ok=True)
    from concurrent.futures import ThreadPoolExecutor, as_completed
    done = total = ok = cached = missing = 0
    started = time.time()
    with ThreadPoolExecutor(args.workers) as pool:
        futures = {pool.submit(do_day, (d, kinds)): d for d in days}
        for future in as_completed(futures):
            try:
                rec = future.result()
            except Exception as err:
                log(f"  {futures[future]}: {err}")
                done += 1
                continue
            total += rec["bytes"]; ok += rec["ok"]
            cached += rec["cached"]; missing += rec["missing"]
            done += 1
            rate = done / max(time.time() - started, 1e-9)
            log(f"  {done}/{len(days)}  {total/1e9:.2f} GB  "
                f"ETA {(len(days)-done)/rate/60:.1f} min")

    log(f"\n{ok} files downloaded, {cached} already present, {missing} missing")
    log(f"  {total/1e9:.2f} GB this run")
    return 0


if __name__ == "__main__":
    sys.exit(main())
