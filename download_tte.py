"""
Download GBM time-tagged event (TTE) data for selected bursts.

Only the 12 NaI detectors are fetched. The two BGOs are deliberately skipped:
legacy_DoL's detector mask excludes them, and zeroing their counts was checked
against four bursts and changed the fitted position by 0.0000 deg in every
case. That saves 22 MB of 88 MB per burst.

TTE is ~66 MB per burst against 127 KB for trigdat, so this is selective by
design -- pass a trigger list or a burst class, never the whole catalog.

Resumable: a detector file already on disk and intact is skipped.

Usage:
    python download_tte.py --class short --limit 8
    python download_tte.py --triggers bn200826187,bn181120265
    python download_tte.py --class short --workers 6
"""
import argparse
import csv
import json
import os
import re
import sys
import threading
import time
import urllib.error
import urllib.request

BASE = "https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers"
MANIFEST = "catalog/fermigbrst.csv"
DATADIR = "data/tte"
LOGFILE = "catalog/tte_download_log.jsonl"

# The 12 NaI detectors. BGOs (b0, b1) are not fetched; see the module docstring.
NAI = [f"n{d}" for d in "0123456789ab"]

USER_AGENT = "fermi-grb-tte/1.0 (academic use)"
RETRIES = 4
BACKOFF = 2.0
TIMEOUT = 600
SHORT_LONG_CUT = 2.0

_print_lock = threading.Lock()


def log(msg):
    with _print_lock:
        print(msg, flush=True)


def fetch(url):
    last = None
    for attempt in range(RETRIES):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
                return resp.read()
        except urllib.error.HTTPError as err:
            if err.code == 404:
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


def save(path, payload):
    tmp = path + ".part"
    with open(tmp, "wb") as handle:
        handle.write(payload)
    os.replace(tmp, path)


def do_burst(trigger):
    year = f"20{trigger[2:4]}"
    dest = os.path.join(DATADIR, trigger)
    record = {"trigger": trigger, "got": [], "cached": [], "missing": [], "bytes": 0}

    have = {}
    if os.path.isdir(dest):
        for name in os.listdir(dest):
            m = re.match(r"glg_tte_(n[0-9ab])_", name)
            if m and is_intact(os.path.join(dest, name)):
                have[m.group(1)] = name

    need = [d for d in NAI if d not in have]
    record["cached"] = sorted(have)
    if not need:
        record["status"] = "cached"
        return record

    listing = fetch(f"{BASE}/{year}/{trigger}/current/")
    if listing is None:
        record["status"] = "no_directory"
        return record
    text = listing.decode("utf-8", "replace")

    os.makedirs(dest, exist_ok=True)
    for det in need:
        # Take the highest version present, as for trigdat.
        names = sorted(set(re.findall(
            rf"glg_tte_{det}_{trigger}_v(\d+)\.fit", text)))
        if not names:
            record["missing"].append(det)
            continue
        fname = f"glg_tte_{det}_{trigger}_v{names[-1]}.fit"
        payload = fetch(f"{BASE}/{year}/{trigger}/current/{fname}")
        if payload is None:
            record["missing"].append(det)
            continue
        save(os.path.join(dest, fname), payload)
        record["got"].append(det)
        record["bytes"] += len(payload)

    record["status"] = "ok" if record["got"] else (
        "incomplete" if record["missing"] else "cached")
    return record


def select(args):
    rows = {r["trigger_name"].strip(): r
            for r in csv.DictReader(open(MANIFEST, newline=""))}
    if args.triggers:
        return [t.strip() for t in args.triggers.split(",") if t.strip()]

    # Default pool: bursts that already have a localization result, so the TTE
    # run has something to be compared against.
    import glob
    done = {os.path.splitext(os.path.basename(f))[0]
            for f in glob.glob("site/all/results/*.json")}
    out = []
    for trigger in sorted(done):
        row = rows.get(trigger)
        if row is None:
            continue
        try:
            t90 = float(row["t90"])
        except (TypeError, ValueError):
            continue
        kind = "short" if t90 < SHORT_LONG_CUT else "long"
        if args.cls != "all" and kind != args.cls:
            continue
        if args.with_truth:
            res = json.load(open(f"site/all/results/{trigger}.json"))
            if "external" not in res:
                continue
        out.append(trigger)
    if args.limit:
        out = out[: args.limit]
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--triggers", type=str, default="")
    parser.add_argument("--class", dest="cls",
                        choices=("all", "short", "long"), default="short")
    parser.add_argument("--with-truth", action="store_true",
                        help="only bursts that have an external position")
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--workers", type=int, default=5)
    args = parser.parse_args()

    triggers = select(args)
    os.makedirs(DATADIR, exist_ok=True)
    os.makedirs(os.path.dirname(LOGFILE), exist_ok=True)

    est = len(triggers) * 66 / 1024
    free = os.statvfs(".").f_bavail * os.statvfs(".").f_frsize / 1e9
    log(f"{len(triggers)} bursts x 12 NaI detectors, ~{est:.1f} GB estimated, "
        f"{free:.0f} GB free, {args.workers} workers")
    if est > free * 0.8:
        log("  refusing: that would use more than 80% of free disk")
        return 1

    from concurrent.futures import ThreadPoolExecutor, as_completed
    counts = {}
    total = done = 0
    started = time.time()
    with open(LOGFILE, "a") as logfh, ThreadPoolExecutor(args.workers) as pool:
        futures = {pool.submit(do_burst, t): t for t in triggers}
        for future in as_completed(futures):
            trigger = futures[future]
            try:
                rec = future.result()
            except Exception as err:
                rec = {"trigger": trigger, "status": "error", "error": str(err)}
            counts[rec["status"]] = counts.get(rec["status"], 0) + 1
            total += rec.get("bytes", 0)
            done += 1
            logfh.write(json.dumps(rec) + "\n")
            logfh.flush()
            if done % 5 == 0 or done == len(triggers):
                rate = done / max(time.time() - started, 1e-9)
                log(f"  {done}/{len(triggers)}  {total/1e9:.2f} GB  "
                    f"{rate*60:.1f}/min  ETA {(len(triggers)-done)/rate/60:.0f} min")

    log("\nfinished: " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    log(f"  {total/1e9:.2f} GB downloaded this run")
    return 0


if __name__ == "__main__":
    sys.exit(main())
