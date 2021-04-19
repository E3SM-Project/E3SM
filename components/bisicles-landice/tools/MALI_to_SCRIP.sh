#!/bin/bash

# ================
# Usage Subroutine
# ================

usage() {
  echo 'MALI_to_SCRIP.sh'
  echo 'Creates a grid file in the SCRIP format from an MALI grid file.'
  echo ''
  echo 'usage:'
  echo './`basename $0` -mpas_file INPUT_FILENAME -scrip_file OUTPUT_FILENAME'
  echo '  -mpas_file   Input MALI grid file'
  echo '  -scrip_file  Output SCRIP grid file'
  echo '  [--help|-h]  Show usage directions'
}

# ========================
# Error / Exit Subroutines
# ========================

err_exit() {
  echo "ERROR: $1"
  exit 1
}

use_err_exit() {
  echo "ERROR: $1"
  echo "***********************"
  usage
  exit 1
}

# ============
# Main Program
# ============

MPAS_SPECIFIED=0
SCRIP_SPECIFIED=0

while [ $# -gt 0 ]; do
  case $1 in
    -h|--help)
      usage
      exit 0
    ;;
    -mpas_file)
      MPAS_SPECIFIED=1
      MPAS_FILE=$2
      if [ ! -e $MPAS_FILE ]; then
        err_exit "can not find MPAS grid file $MPAS_FILE"
      fi
      shift
    ;;
    -scrip_file)
      SCRIP_SPECIFIED=1
      SCRIP_FILE=$2
      if [ -e $SCRIP_FILE ]; then
        err_exit "SCRIP grid file $SCRIP_FILE already exists!"
      fi
      touch $SCRIP_FILE || err_exit "You do not have permission to create SCRIP file $SCRIP_FILE"
      rm -f $SCRIP_FILE
      shift
    ;;
    * )
      use_err_exit "ERROR: $1 is not a valid argument"
    ;;
  esac
  shift
done

if [ $MPAS_SPECIFIED -eq 0 ]; then
  use_err_exit "You must specify an input file!"
fi

if [ $SCRIP_SPECIFIED -eq 0 ]; then
  use_err_exit "You must specify an output file!"
fi

ncl MPAS_to_SCRIP.ncl in_file=\"$MPAS_FILE\" out_file=\"$SCRIP_FILE\"

