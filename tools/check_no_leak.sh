#!/usr/bin/env bash
# Pre-push guard. Fails loudly if anything that must stay private has crept in.
set -u
fail=0
say(){ printf '%s\n' "$*"; }

# 1) No local absolute paths
if grep -rIl '/Volumes/' . --exclude-dir=.git --exclude=check_no_leak.sh 2>/dev/null | grep -q .; then
  say "FAIL: local absolute path (/Volumes/...) found in:"
  grep -rIl '/Volumes/' . --exclude-dir=.git --exclude=check_no_leak.sh; fail=1
fi

# 2) No reference to the private data root or to datasets 001/002
if grep -rIl -e '05data' -e '001_apple' -e '002_apple' . --exclude-dir=.git --exclude=check_no_leak.sh --exclude=.gitignore 2>/dev/null | grep -q .; then
  say "FAIL: reference to the private data root or to datasets 001/002:"
  grep -rIl -e '05data' -e '001_apple' -e '002_apple' . --exclude-dir=.git --exclude=check_no_leak.sh --exclude=.gitignore; fail=1
fi

# 3) No bulk/raw data blobs
if find . -path ./.git -prune -o \( -name '*.npz' -o -name '*.npy' -o -name '*.xlsx' -o -name '*.mat' -o -name '*.parquet' \) -print 2>/dev/null | grep -q .; then
  say "FAIL: raw/bulk data file present:"
  find . -path ./.git -prune -o \( -name '*.npz' -o -name '*.npy' -o -name '*.xlsx' -o -name '*.mat' -o -name '*.parquet' \) -print; fail=1
fi

[ "$fail" -eq 0 ] && say "OK: no private paths, no 001/002 references, no raw data blobs."
exit "$fail"
