#!/usr/bin/env bash
set -euo pipefail

log_file="/pscratch/sd/a/agett/e3sm_scratch/e3smv3_amip_test2/run/atm.log.49472970.260228-033054"

if [[ ! -f "$log_file" ]]; then
  echo "ERROR: File not found: $log_file" >&2
  exit 1
fi

# Keep only content after the last occurrence of "MASTER FIELD LIST".
section=$(awk '
  /MASTER FIELD LIST/ {buf=""; seen=1; next}
  seen {buf = buf $0 "\n"}
  END {
    if (!seen) exit 2
    printf "%s", buf
  }
' "$log_file") || {
  rc=$?
  if [[ $rc -eq 2 ]]; then
    echo "ERROR: 'MASTER FIELD LIST' not found in file." >&2
  else
    echo "ERROR: Failed while extracting MASTER FIELD LIST section." >&2
  fi
  exit 1
}

fields=(
  PHIS LANDFRAC T Q U V TS ICEFRAC OCNFRAC SST TREFHT
  FSUTOA FLUT FSNS FSDS FLDS FLNS PRECT LHFLX SHFLX LWC SWCF LWCF
  TH7001000 OMEGA700 OMEGA500 U10 TGCLDLWP TGCLDIWP CLOUD CDNUMC
  AEROD_v AODVIS CDNUMC CLDTOT EMIS TOT_CLD_VISTAU BURDENSO4 FSNTOA
  FLUTC FLNTC FLNSC FSNTC FSNTOAC FSDSC TAUX TAUY CLDLOW Z500
)

missing=0
for f in "${fields[@]}"; do
  if ! grep -Fqw -- "$f" <<< "$section"; then
    echo "$f"
    missing=1
  fi
done

if [[ $missing -eq 0 ]]; then
  echo "All requested strings were found in the MASTER FIELD LIST section."
fi
