#!/bin/bash
#
# $ ./run_merge_mapping_files --map-in-oo-fname file1_name   \
#                             --map-in-ms-fname file2_name   \
#                             --region-mask-fname file3_name \
#                             --map-out-fname file4_name

usage() {
  echo "$ ./run_merge_mapping_files --map-in-oo-fname file1_name   \\"
  echo "                            --map-in-ms-fname file2_name   \\"
  echo "                            --region-mask-fname file3_name \\"
  echo "                            --map-out-fname file4_name     \\"
  echo "                            [-h -v]"
  echo ""
  echo "  -h,--help     show this message"
  echo "  -v,--verbose  show this message"
  echo "  file1_name    entries into the open ocean"
  echo "  file2_name    entries into marginal seas"
  echo "  file3_name    marginal seas designation"
  echo "  file4_name    output file"
}

parse_args() {

  export NCL_N_ARGS=4
  NCL_ARG_1=""
  NCL_ARG_2=""
  NCL_ARG_3=""
  NCL_ARG_4=""
  VERBOSE=""

  while [ $# -gt 0 ]; do
    case $1 in
      --map-in-oo-fname )
        export NCL_ARG_1=$2
        shift
      ;;
      --map-in-ms-fname )
        export NCL_ARG_2=$2
        shift
      ;;
      --region-mask-fname )
        export NCL_ARG_3=$2
        shift
      ;;
      --map-out-fname )
        export NCL_ARG_4=$2
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

  if [ -z ${NCL_ARG_1} ] || \
     [ -z ${NCL_ARG_2} ] || \
     [ -z ${NCL_ARG_3} ] || \
     [ -z ${NCL_ARG_4} ]; then
    echo "ERROR: Missing argument(s)!"
    echo ""
    echo "Run $0 -h for help"
    exit 1
  fi
}

print_summary() {
  echo "Open ocean map: $NCL_ARG_1"
  echo "Marginal seas map: $NCL_ARG_2"
  echo "Region mask: $NCL_ARG_3"
  echo "Output map: $NCL_ARG_4"
}

parse_args $@
if [ ! -z $VERBOSE ]; then
  print_summary
  echo "Running ncl script..."
  echo ""
fi

ncl ncl/merge_mapping_files.ncl
