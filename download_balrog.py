"""
Download BALROG localizations from the MPE catalogue.

BALROG (Burgess et al. 2018) is an independent Bayesian localization for GBM
that fits sky position and spectrum simultaneously, explicitly to remove the
systematics DoL leaves behind. MPE runs it automatically on GBM triggers and
publishes the result per burst:

    page   https://grb.mpe.mpg.de/grb/GRB<id>/
    map    /media/data/GRB<id>/data_files/GRB<id>_healpix_v<NN>.fits

The maps are HEALPix NSIDE 32 (~1.8 deg pixels, against NSIDE 128 for the
official Fermi maps), about 57 KB each. Not every GBM trigger has one --
roughly four in five do.

The page is fetched first and its healpix link read off it, rather than
guessing the version suffix.

Resumable: a burst whose map is already on disk and intact is skipped.

Usage:
    python download_balrog.py --with-truth
    python download_balrog.py --workers 4
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

BASE = "https://grb.mpe.mpg.de"
RESULTS = "site/all/results"
OUTDIR = "data/balrog"
LOGFILE = "catalog/balrog_download_log.jsonl"

USER_AGENT = "fermi-grb-compare/1.0 (academic use)"
RETRIES = 3
BACKOFF = 2.0
TIMEOUT = 120

_lock = threading.Lock()


def log(msg):
    with _lock:
        print(msg, flush=True)


def fetch(url, binary=True):
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


def quoted_position(text):
    """Pull BALROG's own quoted RA/Dec and errors off the burst page."""
    flat = re.sub(r"<[^>]+>", " ", text)
    flat = re.sub(r"\s+", " ", flat)
    out = {}
    m = re.search(r"Balrog Ra \[deg\]\s*([\d.+-]+)\s*([\d.]+)", flat)
    if m:
        out["ra"] = float(m.group(1))
        out["ra_err"] = float(m.group(2))
    m = re.search(r"Balrog Dec \[deg\]\s*([-\d.]+)\s*([\d.]+)", flat)
    if m:
        out["dec"] = float(m.group(1))
        out["dec_err"] = float(m.group(2))
    m = re.search(r"Estimated systematic error \[deg\]\s*([\d.]+)", flat)
    if m:
        out["sys_err"] = float(m.group(1))
    m = re.search(r"Fitted Model\s*(\w+)", flat)
    if m:
        out["model"] = m.group(1)
    return out


def do_burst(trigger):
    grb = "GRB" + trigger[2:]
    dest = os.path.join(OUTDIR, trigger)
    rec = {"trigger": trigger, "grb": grb}

    existing = glob.glob(os.path.join(dest, "*_healpix_*.fits"))
    meta_path = os.path.join(dest, "balrog.json")
    if existing and is_intact(existing[0]) and os.path.exists(meta_path):
        rec["status"] = "cached"
        return rec

    page = fetch(f"{BASE}/grb/{grb}/")
    if page is None:
        rec["status"] = "no_balrog"
        return rec
    text = page.decode("utf-8", "replace")

    hrefs = re.findall(r'href="(/media/data/[^"]*healpix[^"]*\.fits)"', text)
    if not hrefs:
        rec["status"] = "no_map"
        return rec
    href = sorted(set(hrefs))[-1]

    payload = fetch(BASE + href)
    if payload is None:
        rec["status"] = "map_404"
        return rec

    os.makedirs(dest, exist_ok=True)
    name = os.path.basename(href)
    tmp = os.path.join(dest, name + ".part")
    with open(tmp, "wb") as handle:
        handle.write(payload)
    os.replace(tmp, os.path.join(dest, name))

    meta = quoted_position(text)
    meta["map"] = name
    with open(meta_path, "w") as handle:
        json.dump(meta, handle)

    rec["status"] = "ok"
    rec["bytes"] = len(payload)
    rec["quoted"] = meta
    return rec


def select(args):
    triggers = sorted(os.path.splitext(os.path.basename(f))[0]
                      for f in glob.glob(f"{RESULTS}/*.json"))
    if args.triggers:
        want = {t.strip() for t in args.triggers.split(",")}
        return [t for t in triggers if t in want]
    if args.with_truth:
        keep = []
        for t in triggers:
            with open(f"{RESULTS}/{t}.json") as handle:
                if "external" in json.load(handle):
                    keep.append(t)
        triggers = keep
    if args.limit:
        triggers = triggers[: args.limit]
    return triggers


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--triggers", type=str, default="")
    parser.add_argument("--with-truth", action="store_true",
                        help="only bursts that have an independent position")
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    triggers = select(args)
    os.makedirs(OUTDIR, exist_ok=True)
    os.makedirs(os.path.dirname(LOGFILE), exist_ok=True)
    log(f"{len(triggers)} bursts, ~{len(triggers)*57/1024:.0f} MB if all present, "
        f"{args.workers} workers")

    from concurrent.futures import ThreadPoolExecutor, as_completed
    counts, total, done = {}, 0, 0
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
            if done % 25 == 0 or done == len(triggers):
                rate = done / max(time.time() - started, 1e-9)
                log(f"  {done}/{len(triggers)}  {total/1e6:.0f} MB  "
                    f"ETA {(len(triggers)-done)/rate/60:.0f} min")

    log("\n" + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    log(f"  {total/1e6:.0f} MB downloaded")
    return 0


if __name__ == "__main__":
    sys.exit(main())
