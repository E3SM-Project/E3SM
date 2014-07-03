#!/bin/bash
#===============================================================================
# SVN $Id: gen_cesm_maps.sh 46158 2013-04-19 18:41:34Z mlevy@ucar.edu $
# SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/tools/mapping/trunk_tags/mapping_130509/gen_mapping_files/gen_cesm_maps.sh $
#
# Create needed mapping files for gen_domain and coupler mapping
# Currently supported on yellowstone, geyser, caldera, and jaguarpf
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
  echo './gen_cesm_maps.sh  '
  echo ' Create a suite of mapping files to use in CESM. Depending on the'
  echo ' grid files specified, the result will be a subset of the following'
  echo ' mapping files: '
  echo ' * atm -> ocn: conservative, bilinear, patch'
  echo ' * ocn -> atm: conservative, bilinear'
  echo ' * atm -> lnd: conservative, bilinear'
  echo ' * lnd -> atm: conservative'
  echo ' * ocn -> lnd: conservative'
  echo ' * lnd -> rtm: conservative'
  echo ' * rtm -> lnd: conservative'
  echo ''
  echo 'gen_cesm_maps.sh '
  echo '  --fileatm|-fatm  input atm_grid_filename'
  echo '  --fileocn|-focn  input ocn_grid_filename'
  echo '  --filelnd|-flnd  input lnd_grid_filename'
  echo '  --filertm|-frtm  input rtm_grid_filename'
  echo '  --nameocn|-nocn  output ocn_name in mapping file' 
  echo '  --nameatm|-natm  output atm_name in mapping file'
  echo '  --namelnd|-nlnd  output lnd_name in mapping file'
  echo '  --namertm|-nrtm  output rtm_name in mapping file'
  echo '  [ --typeocn|tocn ] [regional|global]'
  echo '  [ --typeatm|tatm ] [regional|global]'
  echo '  [ --nogridcheck ]'
  echo '  [ --batch|-b ]'
  echo '  [ --help|-h ]'
  echo '  [ -v ]'
  echo ' '
  echo 'where '
  echo ' --fileatm (or -fatm) '
  echo '   SCRIP grid format atmosphere filename (full pathname)'
  echo ' --fileocn (or -focn) '
  echo '   SCRIP grid format ocean filename (full pathname)'
  echo ' --filelnd (or -flnd) '
  echo '   SCRIP grid format land filename (full pathname), must be global'
  echo ' --filertm (or -frtm) '
  echo '   SCRIP grid format runoff filename (full pathname)'
  echo ' --nameatm (or -natm) '
  echo '   Shortname to use for atm in mapping filename'
  echo ' --nameocn (or -nocn) '
  echo '   Shortname to use for ocn in mapping filename'
  echo ' --namelnd (or -nlnd) '
  echo '   Shortname to use for lnd in mapping filename'
  echo ' --namertm (or -nrtm) '
  echo '   Shortname to use for rtm in mapping filename'
  echo ' --typeocn (or -tocn) '
  echo '   ocean grid type,  valid values are regional or global'
  echo '   default is global'
  echo ' --typeatm (or -tatm) '
  echo '   atm grid type, valid values are regional or global'
  echo '   default is global'
  echo '   value must be global if -frtm and -nrtm are specified'
  echo ' --nogridcheck '
  echo '   By default, script will run consistency check on new'
  echo '   maps; this flag disables these checks'
  echo ' --batch (or -b) '
  echo '   Toggles batch mode usage. If you want to run in batch mode'
  echo '   you need to have a separate batch script for a supported machine'
  echo '   that calls this script interactively - you cannot submit this'
  echo '   script directly to the batch system'
  echo ' -rc '
  echo '   Pass the "--recompile" flag to the ESMF tool'
  echo '   (Only necessary if nothing has been built in ../check_maps/)'
  echo ' -d '
  echo '   toggle debug-only '
  echo ' --help or -h  '
  echo '   displays this help message'
  echo ''
  echo 'Note: if rtm is specified and lnd is not, then this tool will'
  echo '      assume lnd and atm are on the same grid.'
  echo ''
  echo 'You can also set the following env variables:'
  echo '  ESMFBIN_PATH - Path to ESMF binaries '
  echo '                 (Leave unset on yellowstone and caldera and the tool'
  echo '                 will be loaded from modules)'
  echo '  MPIEXEC ------ Name of mpirun executable'
  echo '                 (default is mpirun.lsf on yellowstone and caldera; if'
  echo '                 you run interactively on yellowstone, mpi is not used)'
  echo '  REGRID_PROC -- Number of MPI processors to use'
  echo '                 (default is 8)'
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
       echo "Error status returned from gen_cesm_maps script"
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

interactive="YES"
debug="no"
verbose="no"
type_atm="global"
type_ocn="global"
CheckMapsFlag=""
atm_ocn=0
atm_lnd=0
lnd_rtm=0
ocn_lnd=0


while [ $# -gt 0 ]; do
   case $1 in
       -v)
	   verbose="YES"
	   ;;
       -b|--batch)
	   interactive="NO"
	   ;;
       -focn|--fileocn )
	   focn=$2
	   shift
	   ;;
       -fatm|--fileatm )
	   fatm=$2
	   shift
	   ;;
       -flnd|--filelnd )
	   flnd=$2
	   shift
	   ;;
       -frtm|--filertm )
	   frtm=$2
	   shift
	   ;;
       -nocn|--nameocn )
	   nocn=$2
	   shift
	   ;;
       -natm|--nameatm )
	   natm=$2
	   shift
	   ;;
       -nlnd|--namelnd )
	   nlnd=$2
	   shift
	   ;;
       -nrtm|--namertm )
	   nrtm=$2
	   shift
	   ;;
       -tocn|--typeocn )
	   type_ocn=$2
	   shift
	   ;;
       --recompile|-rc )
	   CheckMapsFlag=-rc
	   echo "Will recompile ESMF gridcheck tool"
	   ;;
       --nogridcheck )
	   SkipGridCheck=TRUE
	   echo "Will not check quality of maps!"
	   ;;
       -tatm|--typeatm )
	   type_atm=$2
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

# Make file and name are specified for all desired grids
# (And that the files exist)
if [ ! -z "$fatm" ]; then
  if [ -z "$natm" ]; then
    echo "If you specify atm grid (-fatm), you must"\
         "also provide atm short name (-natm)."
    echo "Invoke gen_cesm_maps.sh -h for usage"
    exit 1
  fi
  if [ ! -f "${fatm}" ]; then
    echo "Atmosphere grid file does NOT exist: $fatm"
    exit 2
  fi
fi
if [ ! -z "$focn" ]; then
  if [ -z "$nocn" ]; then
    echo "If you specify ocn grid (-focn), you must"\
         "also provide ocn short name (-nocn)."
    echo "Invoke gen_cesm_maps.sh -h for usage"
    exit 3
  fi
  if [ ! -f "${focn}" ]; then
    echo "Ocean grid file does NOT exist: $focn"
    exit 4
  fi
fi
if [ ! -z "$flnd" ]; then
  if [ -z "$nlnd" ]; then
    echo "If you specify lnd grid (-flnd), you must"\
         "also provide lnd short name (-nlnd)."
    echo "Invoke gen_cesm_maps.sh -h for usage"
    exit 5
  fi
  if [ ! -f "${flnd}" ]; then
    echo "Land grid file does NOT exist: $flnd"
    exit 6
  fi
fi
if [ ! -z "$frtm" ]; then
  if [ -z "$nrtm" ]; then
    echo "If you specify rtm grid (-frtm), you must"\
         "also provide rtm short name (-nrtm)."
    echo "Invoke gen_cesm_maps.sh -h for usage"
    exit 7
  fi
  if [ ! -f "${frtm}" ]; then
    echo "Runoff grid file does NOT exist: $frtm"
    exit 8
  fi
fi

# Determine what maps to make
if [ ! -z "$fatm" ] && [ ! -z "$focn" ]; then
  atm_ocn=1
fi

if [ ! -z "$fatm" ] && [ ! -z "$flnd" ]; then
  atm_lnd=1
fi

if [ ! -z "$frtm" ]; then
  if [ ! -z "$flnd" ] || [ ! -z "$fatm" ]; then
    lnd_rtm=1
  fi
  if [ -z "$flnd" ]; then
    if [ "${type_atm}" == "regional" ]; then
      echo "WARNING: Can not use regional atmosphere to create atm2rtm map!"
      lnd_rtm=0
    else
      echo "Assuming atmosphere and land grid are identical, generating"\
           "atm2rof and rof2atm maps"
      nlnd=$natm
      flnd=$fatm
    fi
  fi
fi

if [ ! -z "$focn" ] && [ ! -z "$flnd" ]; then
  ocn_lnd=1
fi

# See if any maps are being made
if [ $((atm_ocn+atm_lnd+lnd_rtm+ocn_lnd)) == 0 ]; then
  echo "ERROR: can not generate any maps based on given input!"
  echo "Invoke gen_cesm_maps.sh -h for usage"
  exit 9
fi

# set some defaults
if [ -z "$REGRID_PROC" ]; then
   REGRID_PROC=8
fi

#-------------------------------------------------------------------------------
# run ESMF_RegridWeightGen
#-------------------------------------------------------------------------------

# Resolve interactive or batch mode command
# NOTE - if you want to run in batch mode - you need to have a separate
# batch file that calls this script interactively - you cannot submit
# this script to the batch system

if [ "$interactive" = "YES" ]; then
    batchrun=""
else
    batchrun="--batch"
fi

cdate=`date +%y%m%d`
file_list=""
make_map=$SDIR/gen_ESMF_mapping_file/create_ESMF_map.sh

#-------------------------------------------------------------------------------
# Make Maps
#-------------------------------------------------------------------------------

if [ ${atm_ocn} == 1 ]; then
  #--- ocn to atm conservative (area avg?) -------------------------------------
  $make_map -fsrc $focn -nsrc $nocn -fdst $fatm -ndst $natm \
            -tsrc $type_ocn -tdst $type_atm -map aave $batchrun
  file_list="$file_list map_${nocn}_TO_${natm}_aave.$cdate.nc"

  #--- ocn to atm bilinear (non-conservative) ----------------------------------
  echo ""
  echo "----------------------------------------------------------"
  echo ""
  $make_map -fsrc $focn -nsrc $nocn -fdst $fatm -ndst $natm \
            -tsrc $type_ocn -tdst $type_atm -map blin $batchrun
  file_list="$file_list map_${nocn}_TO_${natm}_blin.$cdate.nc"

  #--- atm to ocn conservative (area avg?) -------------------------------------
  echo ""
  echo "----------------------------------------------------------"
  echo ""
  $make_map -fsrc $fatm -nsrc $natm -fdst $focn -ndst $nocn \
            -tsrc $type_atm -tdst $type_ocn -map aave $batchrun
  file_list="$file_list map_${natm}_TO_${nocn}_aave.$cdate.nc"

  #--- atm to ocn bilinear (non-conservative) ----------------------------------
  echo ""
  echo "----------------------------------------------------------"
  echo ""
  $make_map -fsrc $fatm -nsrc $natm -fdst $focn -ndst $nocn \
            -tsrc $type_atm -tdst $type_ocn -map blin $batchrun
  file_list="$file_list map_${natm}_TO_${nocn}_blin.$cdate.nc"

  #--- atm to ocn patch mapping (non-conservative) -------------------------------
  echo ""
  echo "----------------------------------------------------------"
  echo ""
  $make_map -fsrc $fatm -nsrc $natm -fdst $focn -ndst $nocn \
            -tsrc $type_atm -tdst $type_ocn -map patc $batchrun
  file_list="$file_list map_${natm}_TO_${nocn}_patc.$cdate.nc"
fi

if [ $atm_lnd == 1 ]; then
  #--- atm to lnd conservative (area avg?) -------------------------------------
  if [ $atm_ocn == 1 ]; then
    echo ""
    echo "----------------------------------------------------------"
    echo ""
  fi
  $make_map -fsrc $fatm -nsrc $natm -fdst $flnd -ndst $nlnd \
            -tsrc $type_atm -tdst global -map aave $batchrun
  file_list="$file_list map_${natm}_TO_${nlnd}_aave.$cdate.nc"

  #--- atm to lnd bilinear (non-conservative) ----------------------------------
  echo ""
  echo "----------------------------------------------------------"
  echo ""
  $make_map -fsrc $fatm -nsrc $natm -fdst $flnd -ndst $nlnd \
            -tsrc $type_atm -tdst global -map blin $batchrun
  file_list="$file_list map_${natm}_TO_${nlnd}_blin.$cdate.nc"

  #--- lnd to atm conservative (area avg?) -------------------------------------
  echo ""
  echo "----------------------------------------------------------"
  echo ""
  $make_map -fsrc $flnd -nsrc $nlnd -fdst $fatm -ndst $natm \
            -tsrc global -tdst $type_atm -map aave $batchrun
  file_list="$file_list map_${nlnd}_TO_${natm}_aave.$cdate.nc"

fi

if [ $ocn_lnd == 1 ]; then
  #--- ocn to lnd conservative (area avg) -------------------------------------
  if [ $atm_ocn == 1 ] || [ $atm_lnd == 1 ]; then
    echo ""
    echo "----------------------------------------------------------"
    echo ""
  fi
  $make_map -fsrc $focn -nsrc $nocn -fdst $flnd -ndst $nlnd \
            -tsrc $type_ocn -tdst global -map aave $batchrun
  file_list="$file_list map_${nocn}_TO_${nlnd}_aave.$cdate.nc"
fi

if [ $lnd_rtm == 1 ]; then
  #--- lnd to rtm conservative (area avg) -------------------------------------
  if [ $atm_ocn == 1 ] || [ $atm_lnd == 1 ] || [$ocn_lnd == 1]; then
    echo ""
    echo "----------------------------------------------------------"
    echo ""
  fi
  $make_map -fsrc $flnd -nsrc $nlnd -fdst $frtm -ndst $nrtm \
            -tsrc global -tdst global -map aave $batchrun
  file_list="$file_list map_${nlnd}_TO_${nrtm}_aave.$cdate.nc"

  #--- rtm to lnd conservative (area avg) -------------------------------------
  echo ""
  echo "----------------------------------------------------------"
  echo ""
  $make_map -fsrc $frtm -nsrc $nrtm -fdst $flnd -ndst $nlnd \
            -tsrc global -tdst global -map aave $batchrun
  file_list="$file_list map_${nrtm}_TO_${nlnd}_aave.$cdate.nc"

fi

# Run ESMF Regrid Weight Check tool
if [ ! -z "$SkipGridCheck" ]; then
	echo "Skipping the consistency check"
	rm -f hostfile
	exit 0
fi
echo ""
echo "----------------------------------------------------------"
echo ""

CHECK_MAP="$SDIR/../check_maps/check_map.sh ${CheckMapsFlag}"

echo "Running ESMF Regrid Weight Check tool for each generated map"
echo CHECK_MAP = ${CHECK_MAP}
echo file_list = ${file_list}
echo ""
runcmd "${CHECK_MAP} ${file_list}"

rm -f hostfile

if [ ! -z "$MP_EUIDEVICE_tmp" ]; then
	# Re-enable MP_EUIDEVICE (if previously defined)
	export MP_EUIDEVICE=$MP_EUIDEVICE_tmp
fi

if [ ! -z "$MP_INSTANCES_tmp" ]; then
	# Re-enable MP_EUIDEVICE (if previously defined)
	export MP_INSTANCES=$MP_INSTANCES_tmp
fi

exit 0

#===============================================================================
