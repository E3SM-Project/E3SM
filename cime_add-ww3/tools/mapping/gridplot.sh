#!/bin/bash

usage() {
  echo 'USAGE: gridplot.sh -gridfile SCRIP_GRID_FILE -gridname GRID_SHORT_NAME [-plottype FILE_FORMAT_OF_PLOT] [-orthographic] [-stereographic] [-center_lat LATITUDE] [-center_lon LONGITUDE] [-out_dir OUTPUT_DIRECTORY]'
  echo ''
  echo 'Required flags:'
  echo '-gridfile                 Full path of SCRIP grid file to plot'
  echo '-gridname                 Short name of grid (used for title in plot)'
  echo ''
  echo 'Optional flags:'
  echo '-plottype                 Type of plot to make (pdf, png, ps, x11); default is pdf'
  echo '-ortho|-orthographic      Only use orthographic projection'
  echo '-stereo|-stereographic    Only use stereographic projection'
  echo '-center_lon               Latitude at center of orthographic projection; default is 20 N'
  echo '-center_lat               Longitude at center of orthographic projection; default is 10 W'
  echo '-out_dir                  Directory where plot[s] will be saved'
}

SPECIFIED_PLOT=0
while [ $# -gt 0 ]; do
  case $1 in
    -gridfile )
      GRIDFILE=$2
      shift
      if [ ! -e $GRIDFILE ]; then
        echo "ERROR: can not find grid file $GRIDFILE"
        exit 1
      fi
    ;;
    -gridname )
      GRIDNAME=$2
      shift
    ;;
    -plottype )
      PLOTTYPE=$2
      shift
    ;;
    -ortho|-orthographic )
      NCL_OPTS="${NCL_OPTS} plot_ortho=True"
      SPECIFIED_PLOT=1
    ;;
    -stereo|-stereographic )
      NCL_OPTS="${NCL_OPTS} plot_stereo=True"
      SPECIFIED_PLOT=1
    ;;
    -center_lat )
      NCL_OPTS="${NCL_OPTS} center_lat=$2"
      shift
    ;;
    -center_lon )
      NCL_OPTS="${NCL_OPTS} center_lon=$2"
      shift
    ;;
    -out_dir|-plot_dir )
      if [ ! -d $2 ]; then
        echo "ERROR: $2 is not a valid directory"
        exit 5
      fi
      NCL_OPTS="${NCL_OPTS} out_dir=\"$2\""
      shift
    ;;
    -h|--help )
      usage
      exit 0
    ;;
    * )
      echo "ERROR: $1 is not a valid argument"
      echo ''
      usage
      exit 2
    ;;
  esac
  shift
done

if [ -z "$GRIDFILE" ]; then
  echo "ERROR: You must specify a SCRIP file to plot."
  echo ''
  usage
  exit 3
fi
if [ -z "$GRIDNAME" ]; then
  echo "ERROR: You must specify a short name for the grid file."
  echo ''
  usage
  exit 4
fi

# If not using -stereo or -ortho, generate both plots!
if [ ${SPECIFIED_PLOT} -eq 0 ]; then
  NCL_OPTS="${NCL_OPTS} plot_ortho=True plot_stereo=True"
fi
NCL_OPTS="${NCL_OPTS} gridfile=\"${GRIDFILE}\" gridname=\"${GRIDNAME}\""
if [ ! -z "$PLOTTYPE" ]; then
  NCL_OPTS="$NCL_OPTS out_type=\"${PLOTTYPE}\""
fi
ncl SCRIP2plot.ncl ${NCL_OPTS}
