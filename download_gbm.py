"""
Bulk-download Fermi GBM trigdat and localization maps for every GRB in the
burst catalog.

For each burst the script reads the HEASARC `current/` directory listing once,
then picks the highest-version file of each wanted product. Reading the listing
rather than guessing version suffixes handles the naming shifts across the
mission:

    trigdat       glg_trigdat_all_bnYYMMDDfff_vNN.fit   all years
    healpix       glg_healpix_all_bnYYMMDDfff_vNN.fit   later bursts
    locprob       glg_locprob_all_bnYYMMDDfff_vNN.fit   earlier bursts (legacy
                                                        name for the same kind
                                                        of probability map)

Both map products are requested; a burst that has neither is recorded as such
rather than treated as an error.

Layout:
    data/gbm/bnYYMMDDfff/glg_trigdat_all_bnYYMMDDfff_v01.fit
    data/gbm/bnYYMMDDfff/glg_healpix_all_bnYYMMDDfff_v00.fit

The run is resumable: a burst whose files are already on disk and intact is
skipped without touching the network.

Usage:
    python download_gbm.py                     # everything in the manifest
    python download_gbm.py --limit 6           # first 6 bursts
    python download_gbm.py --triggers bn2008261871,bn170817529
    python download_gbm.py --workers 8
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
DATADIR = "data/gbm"
LOGFILE = "catalog/download_log.jsonl"

# Products to fetch, in priority order within each group. The map products are
# alternatives: healpix supersedes locprob, but older bursts only have locprob.
TRIGDAT = "trigdat"
MAPS = ("healpix", "locprob")

USER_AGENT = "fermi-grb-bulk/1.0 (academic use; astropy/urllib)"
RETRIES = 4
BACKOFF = 2.0
TIMEOUT = 180

_print_lock = threading.Lock()
_log_lock = threading.Lock()


def log(msg):
    with _print_lock:
        print(msg, flush=True)


def fetch(url, binary=True):
    """GET a URL with retries. Returns bytes, or None on a hard 404."""
    last = None
    for attempt in range(RETRIES):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
                return resp.read()
        except urllib.error.HTTPError as err:
            if err.code == 404:
                return None          # genuinely absent, no point retrying
            last = err
        except Exception as err:     # timeouts, connection resets, DNS blips
            last = err
        if attempt < RETRIES - 1:
            time.sleep(BACKOFF * (2 ** attempt))
    raise RuntimeError(f"{url}: {last}")


def pick_versions(listing, trigger):
    """Map product name -> highest-version filename found in a dir listing."""
    found = {}
    pattern = re.compile(
        rf"glg_(\w+)_all_{trigger}(?:_v(\d+))?\.fit\b"
    )
    for match in pattern.finditer(listing):
        product, version = match.group(1), match.group(2)
        if version is None:
            continue             # unversioned quicklook alias; prefer versioned
        version = int(version)
        if product not in found or version > found[product][0]:
            found[product] = (version, match.group(0))
    return {k: v[1] for k, v in found.items()}


def is_intact(path):
    """Cheap validity check: non-empty and starts with a FITS SIMPLE card."""
    try:
        if os.path.getsize(path) < 2880:
            return False
        with open(path, "rb") as handle:
            return handle.read(6) == b"SIMPLE"
    except OSError:
        return False


def save(path, payload):
    """Write atomically so an interrupted run never leaves a half file."""
    tmp = path + ".part"
    with open(tmp, "wb") as handle:
        handle.write(payload)
    os.replace(tmp, path)


def existing(dest):
    """Products already present and intact in a burst directory."""
    have = {}
    if not os.path.isdir(dest):
        return have
    for filename in os.listdir(dest):
        match = re.match(r"glg_(\w+)_all_\w+_v\d+\.fit$", filename)
        if match and is_intact(os.path.join(dest, filename)):
            have[match.group(1)] = filename
    return have


def nomap_marker(dest):
    """Path of the flag noting that the archive holds no map for this burst.

    Without it a re-run would re-request the directory listing for every
    pre-healpix burst, since the map would always look still-missing.
    """
    return os.path.join(dest, ".no_map")


def do_burst(trigger):
    """Fetch trigdat + a localization map for one burst. Returns a log record."""
    year = f"20{trigger[2:4]}"
    dest = os.path.join(DATADIR, trigger)
    record = {"trigger": trigger, "files": {}, "skipped": [], "missing": []}

    have = existing(dest)
    need_trigdat = TRIGDAT not in have
    need_map = not any(m in have for m in MAPS)
    if need_map and os.path.exists(nomap_marker(dest)):
        need_map = False
        record["missing"].append("map")

    if not need_trigdat and not need_map:
        record["status"] = "cached"
        record["files"] = have
        return record

    listing = fetch(f"{BASE}/{year}/{trigger}/current/")
    if listing is None:
        record["status"] = "no_directory"
        return record

    available = pick_versions(listing.decode("utf-8", "replace"), trigger)

    wanted = []
    if need_trigdat:
        if TRIGDAT in available:
            wanted.append(available[TRIGDAT])
        else:
            record["missing"].append(TRIGDAT)
    else:
        record["skipped"].append(TRIGDAT)

    os.makedirs(dest, exist_ok=True)

    if need_map:
        chosen = next((available[m] for m in MAPS if m in available), None)
        if chosen:
            wanted.append(chosen)
        else:
            record["missing"].append("map")
            open(nomap_marker(dest), "w").close()
    else:
        record["skipped"].append("map")

    for filename in wanted:
        url = f"{BASE}/{year}/{trigger}/current/{filename}"
        payload = fetch(url)
        if payload is None:
            record["missing"].append(filename)
            continue
        path = os.path.join(dest, filename)
        save(path, payload)
        product = re.match(r"glg_(\w+)_all_", filename).group(1)
        record["files"][product] = filename
        record.setdefault("bytes", 0)
        record["bytes"] += len(payload)

    record["files"].update(have)
    record["status"] = "ok" if record["files"] else "empty"
    return record


def load_triggers(args):
    with open(MANIFEST, newline="") as handle:
        rows = list(csv.DictReader(handle))
    triggers = [r["trigger_name"].strip() for r in rows]
    triggers = [t for t in triggers if re.fullmatch(r"bn\d{9}", t)]

    if args.triggers:
        wanted = {t.strip() for t in args.triggers.split(",")}
        triggers = [t for t in triggers if t in wanted]
    if args.limit:
        triggers = triggers[: args.limit]
    return triggers


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--triggers", type=str, default="")
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()

    triggers = load_triggers(args)
    os.makedirs(DATADIR, exist_ok=True)
    os.makedirs(os.path.dirname(LOGFILE), exist_ok=True)

    total = len(triggers)
    log(f"{total} bursts to process, {args.workers} workers")

    counts = {"ok": 0, "cached": 0, "no_directory": 0, "empty": 0, "error": 0}
    downloaded = 0
    started = time.time()
    done = 0

    from concurrent.futures import ThreadPoolExecutor, as_completed

    with open(LOGFILE, "a") as logfh, ThreadPoolExecutor(args.workers) as pool:
        futures = {pool.submit(do_burst, t): t for t in triggers}
        for future in as_completed(futures):
            trigger = futures[future]
            try:
                record = future.result()
            except Exception as err:
                record = {"trigger": trigger, "status": "error", "error": str(err)}
            counts[record["status"]] = counts.get(record["status"], 0) + 1
            downloaded += record.get("bytes", 0)
            done += 1

            with _log_lock:
                logfh.write(json.dumps(record) + "\n")
                logfh.flush()

            if done % 50 == 0 or done == total:
                rate = done / max(time.time() - started, 1e-9)
                eta = (total - done) / rate / 60 if rate else 0
                log(
                    f"  {done}/{total}  ok={counts['ok']} cached={counts['cached']} "
                    f"nodir={counts['no_directory']} err={counts['error']}  "
                    f"{downloaded / 1e9:.2f} GB  {rate:.1f}/s  ETA {eta:.0f} min"
                )

    log("\nFinished:")
    for key, value in sorted(counts.items()):
        if value:
            log(f"  {key:14s} {value}")
    log(f"  downloaded     {downloaded / 1e9:.2f} GB this run")
    return 0 if counts["error"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
