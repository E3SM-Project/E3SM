#!/bin/bash
#

################################################################################

main() {
  parse_args $@
  check_inputs
  if [ ! -z $VERBOSE ]; then
    print_summary
    echo "Running ncl script..."
    echo ""
  fi

  ncl ncl/merge_mapping_files.ncl
}

################################################################################

usage() {
  cat << EOF
$ ./run_merge_mapping_files --map_in_oo MAP_IN_OO_FNAME     \\
                            --map_in_ms MAP_IN_MS_FNAME     \\
                            --region_mask REGION_MASK_FNAME \\
                            --map_out MAP_OUT_FNAME         \\
                            [-h -v]

  -h,--help     show this message
  -v,--verbose  echo file names back to screen before running NCL script
  --map_in_oo MAP_IN_OO_FNAME
                mapping file containing map to open ocean points
  --map_in_ms MAP_IN_MS_FNAME
                mapping file containing map to marginal sea points
  --region_mask REGION_MASK_FNAME
                POP region mask file (to specify open ocean vs marginal sea)
  --map_out MAP_OUT_FNAME
                output file
EOF
}

################################################################################

parse_args() {

  # Clear pre-existing environment variables that will be set
  export MAP_IN_OO_FNAME=""
  export MAP_IN_MS_FNAME=""
  export REGION_MASK_FNAME=""
  export MAP_OUT_FNAME=""

  while [ $# -gt 0 ]; do
    case $1 in
      --map_in_oo )
        export MAP_IN_OO_FNAME=$2
        shift
      ;;
      --map_in_ms )
        export MAP_IN_MS_FNAME=$2
        shift
      ;;
      --region_mask )
        export REGION_MASK_FNAME=$2
        shift
      ;;
      --map_out )
        export MAP_OUT_FNAME=$2
        shift
      ;;
      -v | --verbose )
        VERBOSE=TRUE
      ;;
      -h | --help )
        usage
        exit 0
      ;;
      * )
        echo "ERROR: $1 is not a valid argument"
        echo ""
        echo "Run $0 -h for help"
        exit 1
    esac
    shift
  done
}

################################################################################

check_inputs() {

  ERROR=""

  if [ -z ${MAP_IN_OO_FNAME} ]; then
    echo "ERROR: Missing --map_in_oo option"
    ERROR="TRUE"
  fi

  if [ -z ${MAP_IN_MS_FNAME} ]; then
    echo "ERROR: Missing --map_in_ms option"
    ERROR="TRUE"
  fi

  if [ -z ${REGION_MASK_FNAME} ]; then
    echo "ERROR: Missing --region_mask option"
    ERROR="TRUE"
  fi

  if [ -z ${MAP_OUT_FNAME} ]; then
    echo "ERROR: Missing --map_out option"
    ERROR="TRUE"
  fi

  if [ ! -z ${ERROR} ]; then
    echo ""
    echo "Run $0 -h for help"
    exit 1
  fi
}

################################################################################

print_summary() {
  echo "Open ocean map: $MAP_IN_OO_FNAME"
  echo "Marginal seas map: $MAP_IN_MS_FNAME"
  echo "Region mask: $REGION_MASK_FNAME"
  echo "Output map: $MAP_OUT_FNAME"
}

################################################################################

main $@
