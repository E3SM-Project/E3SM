#!/usr/bin/env bash
#===============================================================================
# SVN $Id: create_ESMF_map.sh 65202 2014-11-06 21:07:45Z mlevy@ucar.edu $
# SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/tools/mapping/trunk_tags/mapping_141106/gen_mapping_files/gen_ESMF_mapping_file/create_ESMF_map.sh $
#
# Create needed mapping files for gen_domain and coupler mapping
# Currently supported on cheyenne, geyser, caldera, and pronghorn
#
#===============================================================================
echo $0
date
SDIR=`dirname $0`

#===============================================================================
# Usage subroutine
#===============================================================================
usage() {
  echo ''
  echo '**********************************************************'
  echo 'usage:'
  echo './create_ESMF_map.sh'
  echo ' A wrapper for the ESMF mapping tool that creates a mapping file'
  echo ' from the source grid to the destination grid. Specify what type'
  echo ' of mapping to use with the -maptype flag (aave, blin, bilin,patc,'
  echo ' nearestdtos, or neareststod)'
  echo ''
  echo 'create_ESMF_map.sh'
  echo '  --filesrc|-fsrc  input source grid_filename (required)'
  echo '  --filedst|-fdst  input destination grid_filename (required)'
  echo '  --namesrc|-nsrc  output source name in mapping file (required)'
  echo '  --namedst|-ndst  output destination name in mapping file (required)'
  echo '  --maptype|-map   type of mapping [aave|blin|bilin|patc|nearestdtos|neareststod] (required)'
  echo '  [ --typesrc|tsrc ] [regional|global]'
  echo '  [ --typedst|tdst ] [regional|global]'
  echo '  [ --pass2esmf ]    ["options"]'
  echo '  [ --batch|-b ]'
  echo '  [ --clm_name ]'
  echo '  [ --serial ]'
  echo '  [ -mach|--machine ] [machine_name]'
  echo '  [ --large_file|-big ]'
  echo '  [ --help|-h ]'
  echo '  [ -v ]'
  echo ''
  echo 'where'
  echo ' --filesrc (or -fsrc)'
  echo '   SCRIP grid format source filename (full pathname)'
  echo ' --filedst (or -fdst)'
  echo '   SCRIP grid format destination filename (full pathname)'
  echo ' --namesrc (or -nsrc) and --namesrc (or -nsrc) will result in the'
  echo '   following mapping files'
  echo '     namesrc_TO_namedst_maptype.cdate.nc'
  echo ''
  echo ' --typesrc (or -tsrc)'
  echo '   source grid type,  valid values are regional or global'
  echo '   default is global'
  echo ' --typedst (or -tdst)'
  echo '   destination grid type, valid values are regional or global'
  echo '   default is global'
  echo ' --pass2esmf'
  echo '   pass options directly to the ESMF tool.'
  echo ' --clm_name'
  echo '   Use the CLM naming convention'
  echo ' --serial'
  echo '   Load the serial ESMF tools rather than the parallel tools'
  echo '   (necessary for mapping grids with a single point).'
  echo ' --machine (or -mach)'
  echo '   Name of the machine you are running on. Currently supports cheyenne,'
  echo '   geyser, caldera, and pronghorn. Note that this script will'
  echo '   determine the machine name automatically from the hostfile command.'
  echo ' -d'
  echo '   toggle debug-only'
  echo ' --help or -h'
  echo '   displays this help message'
  echo ''
  echo 'You can also set the following env variables:'
  echo '  ESMFBIN_PATH - Path to ESMF binaries'
  echo '                 (Leave unset on cheyenne/caldera/geyser/pronghorn'
  echo '                 and the tool will be loaded from modules)'
  echo '  MPIEXEC ------ Name of mpirun executable'
  echo '                 (ignored if --serial, which is default on cheyenne'
  echo '                 login nodes'
  echo '**********************************************************'
}

#===============================================================================
# runcmd subroutine
#===============================================================================
runcmd() {
  cmd=$@
  if [ -z "$cmd" ]; then
    echo "No command given to the runcmd function"
    exit 3
  fi

  if [ "$verbose" = "YES" ]; then
    echo "$cmd"
  fi

  if [ "$debug" != "YES" ]; then
    ${cmd}
    rc=$?
  else
    rc=0
  fi

  if [ $rc != 0 ]; then
    echo "Error status returned from create_ESMF_map script"
    exit 4
    undo
   fi

   return 0
}

#===============================================================================
# Main program
#===============================================================================

#-------------------------------------------------------------------------------
# Process input arguments
#-------------------------------------------------------------------------------

debug="no"
verbose="no"
type_src="global"
type_dst="global"
pass_thru=""
CLMNAME="FALSE"
serial="FALSE"
MACH="UNSET"
use_large="false"
CheckMapsFlag=""
use_rtm=0

while [ $# -gt 0 ]; do
  case $1 in
    -v )
      verbose="YES"
    ;;
    --clm_name )
      CLMNAME="TRUE"
    ;;
    --serial )
      serial="TRUE"
    ;;
    -mach|--machine )
      MACH=$2
      shift
    ;;
    -big|--large_file )
      use_large="true"
    ;;
    -fsrc|--filesrc )
      fsrc=$2
      shift
    ;;
    -fdst|--filedst )
      fdst=$2
      shift
    ;;
    -nsrc|--namesrc )
      nsrc=$2
      shift
    ;;
    -ndst|--namedst )
      ndst=$2
      shift
    ;;
    -map|--maptype )
      map_type=$2
      echo "map_type is $map_type"
      shift
    ;;
    -tsrc|--typesrc )
      type_src=$2
      echo "type_src is $type_src"
      shift
    ;;
      -tdst|--typedst )
      type_dst=$2
      echo "type_dst is $type_dst"
      shift
    ;;
    --pass2esmf )
      pass_thru=$2
      echo "Sending ESMF_RegridWeightGen $pass_thru"
      shift
    ;;
    -h|--help )
      usage
      exit 0
    ;;
    * )
      echo "****************************"
      echo "ERROR:: invalid argument $1"
      echo "****************************"
      usage
      exit 1
    ;;
  esac

  shift
done

#-------------------------------------------------------------------------------
# Determine machine to run on
#-------------------------------------------------------------------------------

shopt -s extglob
if [ $MACH == "UNSET" ]; then
  hostname=`hostname`
  case $hostname in
    cheyenne* )
      MACH="cheyenne"
      cheyenne_login=TRUE
    ;;
    r+([0-9])i+([0-9])n+([0-9]) )
      MACH="cheyenne"
    ;;
    geyser* )
      MACH="dav"
    ;;
    caldera* )
      MACH="dav"
    ;;
    pronghorn* )
      MACH="dav"
    ;;
    *)
      echo "Can not determine machine name from hostname '$hostname'"
      exit 1
    ;;
  esac
fi

# machine-specific restrictions
if [ "$cheyenne_login" == "TRUE" ] && [ "$serial" != "TRUE" ]; then
  echo "ERROR: You are trying to use parallal ESMF tools on the cheyenne login node."
  echo "       Either run with '--serial' or move to a compute node."
  exit 1
fi

# check for required arguments
echo "fsrc is $fsrc"
echo "fdst is $fdst"
if [ -z "$fsrc" ]; then
  echo "Must specfiy -fsrc or --filesrc argument "
  echo "Invoke create_ESMF_map.sh -h for usage"
  exit 1
fi

if [ -z "$fdst" ]; then
  echo "Must specfiy -fdst or --filedst argument "
  echo "Invoke create_ESMF_map.sh -h for usage"
  exit 2
fi

if [ -z "$nsrc" ]; then
  echo "Must specfiy -nsrc or --namesrc argument "
  echo "Invoke create_ESMF_map.sh -h for usage"
  exit 3
fi

if [ -z "$ndst" ]; then
  echo "Must specfiy -ndst or --namedst argument "
  echo "Invoke create_ESMF_map.sh -h for usage"
  exit 4
fi

if [ -z "$map_type" ]; then
  echo "Must specfiy -map or --maptype argument "
  echo "Invoke create_ESMF_map.sh -h for usage"
  exit 5
fi

# check for existence of files
if [ ! -f "${fsrc}" ]; then
  echo "Source grid file does NOT exist: $fsrc}"
  exit 6
fi

if [ ! -f "${fdst}" ]; then
  echo "Destination grid file does NOT exist: $fdst"
  exit 7
fi

# check for type of map
if [ $map_type != "aave" ] && [ $map_type != "blin" ] && [ $map_type != "bilin" ] && [ $map_type != "nearestdtos" ] && [ $map_type != "neareststod" ] && [ $map_type != "patc" ]; then
  echo "ERROR: $map_type is not a valid type of mapping."
  echo "(must be aave, blin, bilin, patc, nearestdtos, or neareststod)"
  exit 8
fi

#-------------------------------------------------------------------------------
# Machine specific env stuff
#-------------------------------------------------------------------------------

case $MACH in
  ## cheyenne
  "cheyenne" )
    module purge
    module load intel/17.0.1 esmf_libs/7.0.0
    if [ "$serial" == "TRUE" ]; then
      # No MPIEXEC
      if [ -z "$MPIEXEC" ]; then
        MPIEXEC=""
      fi
      module load esmf-7.1.0r-ncdfio-uni-O
    else
      # MPIEXEC should be mpirun -np
      if [ -z "$MPIEXEC" ]; then
        if [ -z "$NCPUS" ]; then
          NCPUS=1
        fi
        MPIEXEC="mpirun -np $NCPUS"
      fi
      module load esmf-7.1.0r-ncdfio-mpi-O
      module load mpt/2.15f
    fi
    # need to load module to access ncatted
    module load nco
  ;;
## geyser, caldera, or pronghorn
  "dav" )
    module purge
    module load intel/17.0.1 esmflibs/7.1.0r
    if [ "$serial" == "TRUE" ]; then
      # No MPIEXEC
      if [ -z "$MPIEXEC" ]; then
        MPIEXEC=""
      fi
      module load esmf-7.1.0r-ncdfio-uni-O
    else
      echo "ERROR: Parallel ESMF tools are not available on $MACH, use --serial"
      exit 1
    fi
    # need to load module to access ncatted
    module load nco
  ;;
  *)
    echo "Machine $MACH NOT recognized"
    exit
  ;;
esac

#-------------------------------------------------------------------------------
# run ESMF_RegridWeightGen
#-------------------------------------------------------------------------------

if [ ! -z $ESMFBIN_PATH ]; then
  ESMF_REGRID="$ESMFBIN_PATH/ESMF_RegridWeightGen"
else
  ESMF_REGRID="ESMF_RegridWeightGen"
fi

# Make sure $ESMF_REGRID is a valid command
command -v $ESMF_REGRID >/dev/null 2>&1 || { echo "Can not find ESMF_RegridWeightGen, make sure it is in your \$PATH or that you specify \$ESMFBIN_PATH." && exit 1; }

# Remove previous log files
rm -f PET*.Log

# Set output map name and create it
cdate=`date +%y%m%d`
if [ $CLMNAME == "TRUE" ]; then
  mapname=${nsrc}_to_${ndst}_${map_type}_da_c
else
  mapname=${nsrc}_TO_${ndst}_${map_type}.
fi
fmap=map_${mapname}${cdate}.nc

echo ""
echo "Creating $fmap..."

mapping="NULL"
case $map_type in
  "aave")
    mapping="conserve"
  ;;
  "blin" | "bilin")
    mapping="bilinear"
  ;;
  "patc")
    mapping="patch -p all"
  ;;
  "nearestdtos")
    mapping="nearestdtos"
  ;;
  "neareststod")
    mapping="neareststod"
  ;;
esac

if [ "$mapping" == "NULL" ]; then
  echo "ERROR: $map_type is not a valid option for --maptype"
  exit 9
fi
cmd="$MPIEXEC $ESMF_REGRID --ignore_unmapped -m $mapping -w $fmap -s $fsrc -d $fdst $pass_thru"

if [ $use_large == "true" ]; then
  cmd="$cmd --64bit_offset"
fi

if [ $type_src == "regional" ]; then
    cmd="$cmd --src_regional"
    echo "regional source"
fi
if [ $type_dst == "regional" ]; then
    cmd="$cmd --dst_regional"
    echo "regional destination"
fi
echo "cmd is $cmd"
if [ $serial == "TRUE" ]; then
  ESMF_RegridWeightGen --version | grep ESMF_VERSION_STRING
fi
echo ""
runcmd $cmd
if [ "$debug" != "YES" ] && [ ! -f "$fmap" ]; then
   echo "Output mapping file $fmap was NOT created: $fmap"
   exit 4
fi
HIST="$cmd"
HOST=`hostname`
ESMF_VER=`$MPIEXEC $ESMF_REGRID --version | grep VERSION_STRING | tr -d ' ' | sed s/:/": "/g`
if [ -z "$ESMF_VER" ]; then
  ESMF_VER="Version number unavailable, if CVS_revision att exists it contains version"
fi
ncatted -h -a creation_date,global,a,c,"`date`" -a created_by,global,a,c,"$LOGNAME" -a run_command,global,a,c,"$HIST" -a ESMF_ver,global,a,c,"$ESMF_VER" -a hostname,global,a,c,"$HOST" $fmap > /dev/null 2>&1 || { echo "WARNING: after creating $fmap, an error was encountered in ncatted. The map was generated correctly but its global attributes were not updated."; exit 0; }

echo "Successfully created mapping file $fmap "
