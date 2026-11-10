#!/usr/bin/env bash
# stage_fcs_gx_db.sh DATABASE DIR [DIR ...]
# Copies the FCS-GX database into the first DIR with room for it (its size
# plus 20 GB) and prints the copy's path on stdout. FCS-GX is much faster
# reading its database from memory (/dev/shm) or a local disk than from a
# network file system. If no DIR has room, prints DATABASE itself (read in
# place) and says so on stderr. The caller deletes the copy when finished.
set -euo pipefail
db=${1:?usage: stage_fcs_gx_db.sh DATABASE DIR [DIR ...]}; shift
need=$(( $(du -sb "$db" | cut -f1) + 20000000000 ))
for d in "$@"; do
  [[ -d "$d" && -w "$d" ]] || continue
  if [[ "$(df -B1 --output=avail "$d" | tail -n 1)" -gt "$need" ]]; then
    dest=$d/boardercollie_gxdb_$$
    cp -r "$db" "$dest.partial"
    mv "$dest.partial" "$dest"
    echo "$dest"
    exit 0
  fi
done
echo "no staging directory had room for $db; reading it in place" >&2
echo "$db"
