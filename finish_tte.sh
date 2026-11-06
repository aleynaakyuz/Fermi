#!/usr/bin/env bash
# Wait for the TTE download to finish, then extract every detector panel.
#
# Two passes. The first takes only bursts with all 14 detectors present. The
# second adds --allow-partial, which accepts a burst the archive is genuinely
# missing a detector for -- safe only once nothing is still downloading.
set -u
PY=/home/aakyuz/miniconda3/envs/gdt/bin/python
cd /home/aakyuz/Fermi

while pgrep -f "download_tte.py --class all" > /dev/null; do sleep 30; done
echo "download finished $(date -Is)"

while pgrep -f "extract_detectors.py" > /dev/null; do sleep 15; done

$PY extract_detectors.py --workers 6
$PY extract_detectors.py --workers 6 --allow-partial
echo "extraction finished $(date -Is)"
