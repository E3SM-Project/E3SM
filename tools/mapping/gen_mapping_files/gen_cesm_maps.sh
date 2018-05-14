#!/bin/bash
#===============================================================================
# Create needed mapping files for gen_domain and coupler mapping
# Currently supported on cheyenne, geyser, caldera, and prognhorn
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
  echo './gen_cesm_maps.sh'
  echo ' Create a suite of mapping files to use in CESM. Depending on the'
  echo ' grid files specified, the result will be a subset of the following'
  echo ' mapping files:'
  echo ' * atm -> ocn: conservative, bilinear, patch'
  echo ' * ocn -> atm: conservative, bilinear'
  echo ' * atm -> lnd: conservative, bilinear'
  echo ' * lnd -> atm: conservative'
  echo ' * ocn -> lnd: conservative'
  echo ' * lnd -> rtm: conservative'
  echo ' * rtm -> lnd: conservative'
  echo ' * lnd -> glc: conservative, bilinear'
  echo ' * glc -> lnd: conservative, bilinear'
  echo ''
  echo 'gen_cesm_maps.sh'
  echo '  --fileatm|-fatm  input atm_grid_filename'
  echo '  --fileocn|-focn  input ocn_grid_filename'
  echo '  --filelnd|-flnd  input lnd_grid_filename'
  echo '  --filertm|-frtm  input rtm_grid_filename'
  echo '  --fileglc|-fglc  input glc_grid_filename'
  echo '  --nameocn|-nocn  output ocn_name in mapping file'
  echo '  --nameatm|-natm  output atm_name in mapping file'
  echo '  --namelnd|-nlnd  output lnd_name in mapping file'
  echo '  --namertm|-nrtm  output rtm_name in mapping file'
  echo '  --nameglc|-nglc  output glc_name in mapping file'
  echo '  [ --typeocn|tocn ] [regional|global]'
  echo '  [ --typeatm|tatm ] [regional|global]'
  echo '  [ --nogridcheck ]'
  echo '  [ --batch|-b ]'
  echo '  [ --help|-h ]'
  echo '  [ -v ]'
  echo ''
  echo 'where'
  echo ' --fileatm (or -fatm)'
  echo '   SCRIP grid format atmosphere filename (full pathname)'
  echo ' --fileocn (or -focn)'
  echo '   SCRIP grid format ocean filename (full pathname)'
  echo ' --filelnd (or -flnd)'
  echo '   SCRIP grid format land filename (full pathname), must be global'
  echo ' --filertm (or -frtm)'
  echo '   SCRIP grid format runoff filename (full pathname)'
  echo ' --fileglc (or -fglc)'
  echo '   SCRIP grid format glc filename (full pathname), assumed to be regional'
  echo ' --nameatm (or -natm)'
  echo '   Shortname to use for atm in mapping filename'
  echo ' --nameocn (or -nocn)'
  echo '   Shortname to use for ocn in mapping filename'
  echo ' --namelnd (or -nlnd)'
  echo '   Shortname to use for lnd in mapping filename'
  echo ' --namertm (or -nrtm)'
  echo '   Shortname to use for rtm in mapping filename'
  echo ' --nameglc (or -nglc)'
  echo '   Shortname to use for glc in mapping filename'
  echo ' --typeocn (or -tocn)'
  echo '   ocean grid type,  valid values are regional or global'
  echo '   default is global'
  echo ' --typeatm (or -tatm)'
  echo '   atm grid type, valid values are regional or global'
  echo '   default is global'
  echo '   value must be global if -frtm and -nrtm are specified'
  echo ' --nogridcheck'
  echo '   By default, script will run consistency check on new'
  echo '   maps; this flag disables these checks'
  echo ' --serial'
  echo '   Run the ESMF tools in serial rather than parallel'
  echo ' -rc'
  echo '   Pass the "--recompile" flag to the ESMF tool'
  echo '   (Only necessary if nothing has been built in ../check_maps/)'
  echo ' --help or -h'
  echo '   displays this help message'
  echo ''
  echo 'Note: if rtm or glc are specified and lnd is not, then this tool will'
  echo '      assume lnd and atm are on the same grid.'
  echo ''
  echo 'You can also set the following env variables:'
  echo '  ESMFBIN_PATH - Path to ESMF binaries'
  echo '                 (Known machines will load tools from modules)'
  echo '  MPIEXEC ------ Name of mpirun executable'
  echo '                 (currently tools only run in serial due to module issues)'
  echo '**********************************************************'
}

#===============================================================================
# make_map subroutine
#===============================================================================
make_map() {
  make_map_exe=$SDIR/gen_ESMF_mapping_file/create_ESMF_map.sh
  fsrc=$1
  nsrc=$2
  fdst=$3
  ndst=$4
  tsrc=$5
  tdst=$6
  map=$7
  if [ "$serial" == "TRUE" ]; then
    run_serial="--serial"
  else
    run_serial=""
  fi
  $make_map_exe -fsrc $fsrc -nsrc $nsrc -fdst $fdst -ndst $ndst  \
                -tsrc $tsrc -tdst $tdst -map $map $run_serial   || exit $?
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
   ${cmd}
   rc=$?
   if [ $rc != 0 ]; then
       echo "Error status $rc returned from gen_cesm_maps script"
       exit $rc
   fi
   return 0
}

#===============================================================================
# Main program
#===============================================================================

#-------------------------------------------------------------------------------
# Process input arguments
#-------------------------------------------------------------------------------

verbose="no"
type_atm="global"
type_ocn="global"
CheckMapsFlag=""
atm_ocn=0
atm_lnd=0
lnd_rtm=0
ocn_lnd=0
lnd_glc=0
serial="FALSE"


while [ $# -gt 0 ]; do
   case $1 in
     -v)
	     verbose="YES"
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
     -fglc|--fileglc )
  	   fglc=$2
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
     -nglc|--nameglc )
  	   nglc=$2
  	   shift
	   ;;
     -tocn|--typeocn )
  	   type_ocn=$2
  	   shift
	   ;;
     --serial )
       serial=TRUE
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
if [ ! -z "$fglc" ]; then
  if [ -z "$nglc" ]; then
    echo "If you specify glc grid (-fglc), you must"\
         "also provide glc short name (-nglc)."
    echo "Invoke gen_cesm_maps.sh -h for usage"
    exit 5
  fi
  if [ ! -f "${fglc}" ]; then
    echo "GLC grid file does NOT exist: $fglc"
    exit 6
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

if [ ! -z "$fglc" ]; then
  if [ ! -z "$flnd" ] || [ ! -z "$fatm" ]; then
    lnd_glc=1
  fi
  if [ -z "$flnd" ]; then
    if [ "${type_atm}" == "regional" ]; then
      echo "WARNING: Can not use regional atmosphere to create atm2glc map!"
      lnd_glc=0
    else
      echo "Assuming atmosphere and land grid are identical, generating"\
           "atm2glc and glc2atm maps"
      nlnd=$natm
      flnd=$fatm
    fi
  fi
fi

if [ ! -z "$focn" ] && [ ! -z "$flnd" ]; then
  ocn_lnd=1
fi

# See if any maps are being made
if [ $((atm_ocn+atm_lnd+lnd_rtm+lnd_glc+ocn_lnd)) == 0 ]; then
  echo "ERROR: can not generate any maps based on given input!"
  echo "Invoke gen_cesm_maps.sh -h for usage"
  exit 9
fi

#-------------------------------------------------------------------------------
# run ESMF_RegridWeightGen
#-------------------------------------------------------------------------------

cdate=`date +%y%m%d`
file_list=""

#-------------------------------------------------------------------------------
# Make Maps
#-------------------------------------------------------------------------------

if [ ${atm_ocn} == 1 ]; then
  #--- ocn to atm conservative (area avg?) -------------------------------------
  make_map "$focn" "$nocn" "$fatm" "$natm" "$type_ocn" "$type_atm" "aave"
  file_list="$file_list map_${nocn}_TO_${natm}_aave.$cdate.nc"

  #--- ocn to atm bilinear (non-conservative) ----------------------------------
  echo ""
  echo "----------------------------------------------------------"
  echo ""
  make_map "$focn" "$nocn" "$fatm" "$natm" "$type_ocn" "$type_atm" "blin"
  file_list="$file_list map_${nocn}_TO_${natm}_blin.$cdate.nc"

  #--- atm to ocn conservative (area avg?) -------------------------------------
  echo ""
  echo "----------------------------------------------------------"
  echo ""
  make_map "$fatm" "$natm" "$focn" "$nocn" "$type_atm" "$type_ocn" "aave"
  file_list="$file_list map_${natm}_TO_${nocn}_aave.$cdate.nc"

  #--- atm to ocn bilinear (non-conservative) ----------------------------------
  echo ""
  echo "----------------------------------------------------------"
  echo ""
  make_map "$fatm" "$natm" "$focn" "$nocn" "$type_atm" "$type_ocn" "blin"
  file_list="$file_list map_${natm}_TO_${nocn}_blin.$cdate.nc"

  #--- atm to ocn patch mapping (non-conservative) -------------------------------
  echo ""
  echo "----------------------------------------------------------"
  echo ""
  make_map "$fatm" "$natm" "$focn" "$nocn" "$type_atm" "$type_ocn" "patc"
  file_list="$file_list map_${natm}_TO_${nocn}_patc.$cdate.nc"
fi

if [ $atm_lnd == 1 ]; then
  #--- atm to lnd conservative (area avg?) -------------------------------------
  if [ ! -z "$file_list" ]; then
    echo ""
    echo "----------------------------------------------------------"
    echo ""
  fi
  make_map "$fatm" "$natm" "$flnd" "$nlnd" "$type_atm" "global" "aave"
  file_list="$file_list map_${natm}_TO_${nlnd}_aave.$cdate.nc"

  #--- atm to lnd bilinear (non-conservative) ----------------------------------
  echo ""
  echo "----------------------------------------------------------"
  echo ""
  make_map "$fatm" "$natm" "$flnd" "$nlnd" "$type_atm" "global" "blin"
  file_list="$file_list map_${natm}_TO_${nlnd}_blin.$cdate.nc"

  #--- lnd to atm conservative (area avg?) -------------------------------------
  echo ""
  echo "----------------------------------------------------------"
  echo ""
  make_map "$flnd" "$nlnd" "$fatm" "$natm" "global" "$type_atm" "aave"
  file_list="$file_list map_${nlnd}_TO_${natm}_aave.$cdate.nc"

fi

if [ $ocn_lnd == 1 ]; then
  #--- ocn to lnd conservative (area avg) -------------------------------------
  if [ ! -z "$file_list" ]; then
    echo ""
    echo "----------------------------------------------------------"
    echo ""
  fi
  make_map "$focn" "$nocn" "$flnd" "$nlnd" "$type_ocn" "global" "aave"
  file_list="$file_list map_${nocn}_TO_${nlnd}_aave.$cdate.nc"
fi

if [ $lnd_rtm == 1 ]; then
  #--- lnd to rtm conservative (area avg) -------------------------------------
  if [ ! -z "$file_list" ]; then
    echo ""
    echo "----------------------------------------------------------"
    echo ""
  fi
  make_map "$flnd" "$nlnd" "$frtm" "$nrtm" "global" "global" "aave"
  file_list="$file_list map_${nlnd}_TO_${nrtm}_aave.$cdate.nc"

  #--- rtm to lnd conservative (area avg) -------------------------------------
  echo ""
  echo "----------------------------------------------------------"
  echo ""
  make_map "$frtm" "$nrtm" "$flnd" "$nlnd" "global" "global" "aave"
  file_list="$file_list map_${nrtm}_TO_${nlnd}_aave.$cdate.nc"

fi

if [ $lnd_glc == 1 ]; then
  #--- glc to lnd conservative (area avg?) -------------------------------------
  if [ ! -z "$file_list" ]; then
    echo ""
    echo "----------------------------------------------------------"
    echo ""
  fi
  make_map "$fglc" "$nglc" "$flnd" "$nlnd" "regional" "global" "aave"
  file_list="$file_list map_${nglc}_TO_${nlnd}_aave.$cdate.nc"

  #--- lnd to glc conservative (area avg?) -------------------------------------
  echo ""
  echo "----------------------------------------------------------"
  echo ""
  make_map "$flnd" "$nlnd" "$fglc" "$nglc" "global" "regional" "aave"
  file_list="$file_list map_${nlnd}_TO_${nglc}_aave.$cdate.nc"

  #--- lnd to glc bilinear (non-conservative) ----------------------------------
  echo ""
  echo "----------------------------------------------------------"
  echo ""
  make_map "$flnd" "$nlnd" "$fglc" "$nglc" "global" "regional" "blin"
  file_list="$file_list map_${nlnd}_TO_${nglc}_blin.$cdate.nc"

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
