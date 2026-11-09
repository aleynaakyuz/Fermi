"""
Download GBM daily position-history files for selected bursts.

trigdat carries a single attitude sample at the trigger and no SAA flag at
all, so it cannot say how long before the burst Fermi last left the South
Atlantic Anomaly. The daily poshist product can: it samples position,
attitude and an in-SAA flag once a second for a whole day.

One file covers a day, so bursts sharing a date share a download. About
12 MB each.

    https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/YYYY/MM/DD/current/
        glg_poshist_all_YYMMDD_vNN.fit

Resumable: a file already on disk and intact is skipped.

Usage:
    python download_poshist.py --triggers bn220827896,bn221009553
    python download_poshist.py --from-json sample.json --workers 4
"""
import argparse
import glob
import json
import os
import re
import sys
import threading
import time
import urllib.error
import urllib.request

BASE = "https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily"
OUTDIR = "data/poshist"

USER_AGENT = "fermi-grb-compare/1.0 (academic use)"
RETRIES = 4
BACKOFF = 2.0
TIMEOUT = 300

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


def day_of(trigger):
    """bnYYMMDDfff -> (YYYY, MM, DD, YYMMDD)."""
    yy, mm, dd = trigger[2:4], trigger[4:6], trigger[6:8]
    return f"20{yy}", mm, dd, f"{yy}{mm}{dd}"


def do_day(key):
    yyyy, mm, dd, yymmdd = key
    existing = glob.glob(os.path.join(OUTDIR, f"glg_poshist_all_{yymmdd}_v*.fit"))
    if existing and is_intact(existing[0]):
        return {"day": yymmdd, "status": "cached", "bytes": 0}

    listing = fetch(f"{BASE}/{yyyy}/{mm}/{dd}/current/")
    if listing is None:
        return {"day": yymmdd, "status": "no_directory", "bytes": 0}
    text = listing.decode("utf-8", "replace")
    names = sorted(set(re.findall(
        rf"glg_poshist_all_{yymmdd}_v\d+\.fit", text)))
    if not names:
        return {"day": yymmdd, "status": "no_file", "bytes": 0}

    name = names[-1]            # highest version
    payload = fetch(f"{BASE}/{yyyy}/{mm}/{dd}/current/{name}")
    if payload is None:
        return {"day": yymmdd, "status": "fetch_404", "bytes": 0}

    os.makedirs(OUTDIR, exist_ok=True)
    tmp = os.path.join(OUTDIR, name + ".part")
    with open(tmp, "wb") as handle:
        handle.write(payload)
    os.replace(tmp, os.path.join(OUTDIR, name))
    return {"day": yymmdd, "status": "ok", "bytes": len(payload)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--triggers", type=str, default="")
    parser.add_argument("--from-json", type=str, default="",
                        help="JSON of {group: [{t: trigger}, ...]}")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    triggers = []
    if args.triggers:
        triggers += [t.strip() for t in args.triggers.split(",") if t.strip()]
    if args.from_json:
        with open(args.from_json) as handle:
            for group in json.load(handle).values():
                triggers += [r["t"] for r in group]
    if not triggers:
        raise SystemExit("nothing to do: pass --triggers or --from-json")

    days = sorted({day_of(t) for t in triggers})
    os.makedirs(OUTDIR, exist_ok=True)
    log(f"{len(triggers)} bursts span {len(days)} days, {args.workers} workers")

    from concurrent.futures import ThreadPoolExecutor, as_completed
    counts, total, done = {}, 0, 0
    started = time.time()
    with ThreadPoolExecutor(args.workers) as pool:
        futures = {pool.submit(do_day, d): d for d in days}
        for future in as_completed(futures):
            try:
                rec = future.result()
            except Exception as err:
                rec = {"day": futures[future][3], "status": "error",
                       "bytes": 0, "error": str(err)}
            counts[rec["status"]] = counts.get(rec["status"], 0) + 1
            total += rec.get("bytes", 0)
            done += 1
            if done % 5 == 0 or done == len(days):
                rate = done / max(time.time() - started, 1e-9)
                log(f"  {done}/{len(days)}  {total/1e6:.0f} MB  "
                    f"ETA {(len(days)-done)/rate/60:.1f} min")

    log("\n" + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    log(f"  {total/1e6:.0f} MB downloaded")
    return 0


if __name__ == "__main__":
    sys.exit(main())
