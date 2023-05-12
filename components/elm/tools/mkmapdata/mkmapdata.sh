#!/bin/bash
#----------------------------------------------------------------------
#
# mkmapdata.sh
#
# Create needed mapping files for mksurfdata_map and CLM.
# 
# Example to run for an output resolution of 4x5
#
# mkmapdata.sh -r 4x5
#
# valid arguments: 
# -f <scripfilename> Input grid filename 
# -t <type>          Output type, supported values are [regional, global]
# -r <res>           Output resolution
# -b                 use batch mode (not default)
# -l                 list mapping files required (so can use check_input_data to get them)
# -d                 debug usage -- display mkmapdata that will be run but don't execute them
# -v                 verbose usage -- log more information on what is happening
# -h                 displays this help message
#
# See usage (below) for additional arguments.
#
# You can also set the following env variables:
#
# ESMFBIN_PATH   -- Path to ESMF binaries
# INPUTDATA_PATH -- Path to root of input data directory
# MPIEXEC        -- Name of mpirun executable
# REGRID_PROC    -- Number of MPI processors to use; Note that this is only used
#                   in setting MPIEXEC for the supported machines. This should
#                   probably go away in the future, because the machines are
#                   likely out of date, and it is just easier to set this via
#                   setting MPIEXEC directly.
#
#----------------------------------------------------------------------
echo $0
dir=${0%/*}
if [ "$dir" = "$0" ];then
  dir="."
fi
outfilelist="clm.input_data_list"
default_res="10x15"

#----------------------------------------------------------------------
# SET SOME DEFAULTS -- if not set via env variables outside

if [ -z "$INPUTDATA_PATH" ]; then
   INPUTDATA_PATH=/glade/p/cesm/cseg/inputdata
fi
if [ -z "$REGRID_PROC" ]; then
   REGRID_PROC=8
fi
#----------------------------------------------------------------------
# Usage subroutine
usage() {
  echo ""
  echo "**********************"
  echo "usage: ./mkmapdata.sh <arguments>"
  echo ""
  echo "valid arguments: "
  echo "[-f|--gridfile <gridname>] "
  echo "    Full pathname of model SCRIP grid file to use "
  echo "    This variable should be set if this is not a supported grid" 
  echo "    This variable will override the automatic generation of the"
  echo "    filename generated from the -res argument "
  echo "    the filename is generated ASSUMING that this is a supported "
  echo "    grid that has entries in the file namelist_defaults_clm.xml"
  echo "    the -r|--res argument MUST be specied if this argument is specified" 
  echo "[-r|--res <res>]"
  echo "    Model output resolution (default is $default_res)"
  echo "[-t|--gridtype <type>]"
  echo "    Model output grid type"
  echo "    supported values are [regional,global], (default is global)"
  echo "[-i|--inputdata-path <inputdata_path>]"
  echo "    Full path to root of inputdata directory"
  echo "[-n|--ntasks <ntasks>]"
  echo "    Number of MPI tasks to use in regrid command. Note that this is only"
  echo "    used when using a recognized machine and MPIEXEC is not set. "
  echo "    Otherwise, setting MPIEXEC (either via environment variable or by "
  echo "    passing --mpiexec <command>) completely specifies the command and "
  echo "    number of tasks. I.e., mkmapdata.sh --mpiexec \"srun --ntasks=10\" "
  echo "    implies that ESMF_RegridWeightGen will be called with srun and 10 "
  echo "    MPI tasks."
  echo "[-m|--mpiexec <command>]"
  echo "    Command used to run mpi jobs (mpirun, aprun, etc)"
  echo "[-e|--esmf-path <path to esmf binaries>]"
  echo "    Path to ESMF executables to use. If not present, will default to"
  echo "    any executables in PATH environment variable"
  echo "[-o|--output-filetype <filetype>]"
  echo "    Specific type for output file [netcdf4, 64bit_offset]"
  echo "[-b|--batch]"
  echo "    Toggles batch mode usage (whether or not to use MPI). If this flag"
  echo "    is not set, then MPIEXEC is ignored and ESMF will be called directly."
  echo "[-l|--list]"
  echo "    List mapping files required (use check_input_data to get them)"
  echo "    also writes data to $outfilelist"
  echo "[-d|--debug]"
  echo "    Toggles debug-only (don't actually run mkmapdata just echo what would happen)"
  echo "[-h|--help]  "
  echo "    Displays this help message"
  echo "[-v|--verbose]"
  echo "    Toggle verbose usage -- log more information on what is happening "
  echo ""
  echo "You can also set the following env variables:"
  echo "    ESMFBIN_PATH   -- Path to ESMF binaries "
  echo "                      (default assume ESMF_RegridWeightGen is in PATH)"
  echo "    INPUTDATA_PATH -- Path to root of input data directory"
  echo "                      (default is /glade/p/cesm/cseg/inputdata)"
  echo "    MPIEXEC        -- Name of mpirun executable"
  echo "                      (default is mpirun.lsf)"
  echo "    REGRID_PROC    -- Number of MPI processors to use"
  echo "                      (default is 8)"
  echo ""
  echo "NOTE: environment variables can be passed by preceding command with"
  echo "      'env var1=setting var2=setting '"
  echo "**********************"
}
#----------------------------------------------------------------------
# runcmd subroutine
#----------------------------------------------------------------------

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
        echo "Error status returned from mkmapdata script"
        exit 4
    fi
    return 0
}

#----------------------------------------------------------------------
# Process input arguments
#----------------------------------------------------------------------

interactive="YES"
debug="no"
res="default"
gridtype="global"
verbose="no"
list="no"
outgrid=""
gridfile="default"
output_filetype=""

while [ $# -gt 0 ]; do
    case $1 in
        -v|--verbose)
            verbose="YES"
            ;;
        -b|--batch) 
            interactive="NO"
            ;;
        -d|--debug)
            debug="YES"
            ;;
        -l|--list)
            debug="YES"
            list="YES"
            ;;
        -r|--res)
            res=$2
            shift
            ;;
        -f|--gridfile)
            gridfile=$2
            shift
            ;;
        -t|--gridtype)
            gridtype=$2
            shift
            ;;
        -i|--inputdata-path)
            INPUTDATA_PATH=$2
            shift
            ;;
        -n|--ntasks)
            REGRID_PROC=$2
            shift
            ;;
        -m|--mpiexec)
            MPIEXEC=$2
            shift
            ;;
        -e|--esmf-path)
            ESMFBIN_PATH=$2
            shift
            ;;
        -o|--output-filetype)
            output_filetype=$2
            shift
            ;;
        -h|--help )
            usage
            exit 0
            ;;
        * )
            echo "ERROR:: invalid argument sent in: $2"
            usage
            exit 1
            ;;
    esac
    shift
done

echo "Script to create mapping files required by mksurfdata_map"

#----------------------------------------------------------------------
# Determine output scrip grid file
#----------------------------------------------------------------------

# Set general query command used below
QUERY="$dir/../../bld/queryDefaultNamelist.pl -silent -namelist elmexp "
QUERY="$QUERY -justvalue -options sim_year=2000 -csmdata $INPUTDATA_PATH"
echo "query command is $QUERY"

echo ""
DST_EXTRA_ARGS=""
if [ "$gridfile" != "default" ]; then
    GRIDFILE=$gridfile
    echo "Using user specified scrip grid file: $GRIDFILE" 
    if [ "$res" = "default" ]; then
        echo "When user specified grid file is given you MUST set the resolution (as the name of your grid)\n";
        exit 1
    fi
    
    # For now, make some assumptions about user-specified grids --
    # that they are SCRIP format, and small enough to not require
    # large file support for the output mapping file. In the future,
    # we may want to provide command-line options to allow the user to
    # override these defaults.
    DST_LRGFIL="none"
    DST_TYPE="SCRIP"
else
    echo "Error: no grid file specified"
    exit 1
    if [ "$res" = "default" ]; then
        res=$default_res
    fi

    QUERYARGS="-res $res -options lmask=nomask"

    # Find the output grid file for this resolution using the XML database
    QUERYFIL="$QUERY -var scripgriddata $QUERYARGS -onlyfiles"
    if [ "$verbose" = "YES" ]; then
        echo $QUERYFIL
    fi
    GRIDFILE=`$QUERYFIL`
    echo "Using default scrip grid file: $GRIDFILE" 

    # Determine extra information about the destination grid file
    DST_LRGFIL=`$QUERY -var scripgriddata_lrgfile_needed $QUERYARGS`
    DST_TYPE=`$QUERY -var scripgriddata_type $QUERYARGS`
    if [ "$DST_TYPE" = "UGRID" ]; then
        # For UGRID, we need extra information: the meshname variable
        dst_meshname=`$QUERY -var scripgriddata_meshname $QUERYARGS`
        DST_EXTRA_ARGS="$DST_EXTRA_ARGS --dst_meshname $dst_meshname"
    fi
fi

# Make sure grid type is consistent
if [ "$gridtype" = "global" ] && [ `echo "$res" | grep -c "1x1_"` = 1 ]; then
    echo "This is a regional resolution and yet it is being run as global, set type with '-t' option\n";
    exit 1
fi
echo "Output grid resolution is $res"

# Make sure we found a gridfile for given resolution
if [ -z "$GRIDFILE" ]; then
    echo "Output grid file was NOT found for this resolution: $res\n";
    exit 1
fi

# Make sure gridfile exists for given resolution
if [ "$list" = "YES" ]; then
    echo "outgrid = $GRIDFILE"
    echo "outgrid = $GRIDFILE" > $outfilelist
elif [ ! -f "$GRIDFILE" ]; then
    echo "Output grid file does NOT exist: $GRIDFILE\n";
    echo "Make sure INPUTDATA_PATH environment variable is set correctly"
    exit 1
fi

#----------------------------------------------------------------------
# Determine all input grid files and output file names 
#----------------------------------------------------------------------

  grids=(                                     \
      "0.5x0.5_AVHRR"                         \
      "0.5x0.5_MODIS"                         \
      "3x3min_LandScan2004"                   \
      "3x3min_MODIS"                          \
      "3x3min_USGS"                           \
      "5x5min_nomask"                         \
      "5x5min_IGBP-GSDP"                      \
      "5x5min_ISRIC-WISE"                     \
      "10x10min_nomask"                       \
      "10x10min_IGBPmergeICESatGIS"           \
      "3x3min_GLOBE-Gardner"                  \
      "3x3min_GLOBE-Gardner-mergeGIS"         \
      "0.9x1.25_GRDC"                         \
      "360x720cru_cruncep"                    \
      "1km-merge-10min_HYDRO1K-merge-nomask"  \
      "0.5x0.5_GSDTG2000"                     \
    )

# Set timestamp for names below 
CDATE="c"`date +%y%m%d`

# Set name of each output mapping file
# First determine the name of the input scrip grid file  
# for each of the above grids
declare -i nfile=1
for gridmask in ${grids[*]}
do
    grid=${gridmask%_*}
    lmask=${gridmask#*_}
 
    QUERYARGS="-res $grid -options lmask=$lmask,glc_nec=10 "
 
    QUERYFIL="$QUERY -var scripgriddata $QUERYARGS -onlyfiles"
    if [ "$verbose" = "YES" ]; then
        echo $QUERYFIL
    fi
    INGRID[nfile]=`$QUERYFIL`
    if [ "$list" = "YES" ]; then
        echo "ingrid = ${INGRID[nfile]}"
        echo "ingrid = ${INGRID[nfile]}" >> $outfilelist
    fi
 
    OUTFILE[nfile]=map_${grid}_${lmask}_to_${res}_nomask_aave_da_$CDATE.nc
 
    # Determine extra information about the source grid file
    SRC_EXTRA_ARGS[nfile]=""
    SRC_LRGFIL[nfile]=`$QUERY -var scripgriddata_lrgfile_needed $QUERYARGS`
    SRC_TYPE[nfile]=`$QUERY -var scripgriddata_type $QUERYARGS`
    if [ "${SRC_TYPE[nfile]}" = "UGRID" ]; then
        # For UGRID, we need extra information: the meshname variable
        src_meshname=`$QUERY -var scripgriddata_meshname $QUERYARGS`
        SRC_EXTRA_ARGS[nfile]="${SRC_EXTRA_ARGS[nfile]} --src_meshname $src_meshname"
    fi
 
    nfile=nfile+1
done

#----------------------------------------------------------------------
# Determine supported machine specific stuff
# TODO: machine-specific stuff should be handled by CIME now, so this
# should all be ripped out.
#----------------------------------------------------------------------

hostname=`hostname`
if [ -n "$NERSC_HOST" ]; then
    hostname=$NERSC_HOST
fi
case $hostname in
    ##yellowstone
    ys* | caldera* | geyser* )
    . /glade/apps/opt/lmod/lmod/init/bash
    module load esmf
    module load ncl
    module load nco
  
    if [ -z "$ESMFBIN_PATH" ]; then
        if [ "$gridtype" = "global" ]; then
            mpi=mpi
            mpitype="mpich2"
        else
            mpi=uni
            mpitype="mpiuni"
        fi
        ESMFBIN_PATH=/glade/apps/opt/esmf/6.3.0-ncdfio/intel/12.1.5/bin/binO/Linux.intel.64.${mpitype}.default
    fi
    if [ -z "$MPIEXEC" ]; then
        MPIEXEC="mpirun.lsf"
    fi
    ;;
  
    ## hopper
    hopper* )
    .  /opt/modules/default/init/bash
    module load ncl/6.1.2
    module load nco
    if [ -z "$ESMFBIN_PATH" ]; then
        module use -a /project/projectdirs/ccsm1/modulefiles/hopper
        if [ "$gridtype" = "global" ]; then
            mpi=mpi
            mpitype="mpi"
        else
            mpi=uni
            mpitype="mpiuni"
        fi
        module load esmf/6.3.0r-ncdfio-${mpitype}-O
        ESMFBIN_PATH=$ESMF_LIBDIR/../bin
    fi
    if [ -z "$MPIEXEC" ]; then
        MPIEXEC="aprun -n $REGRID_PROC"
    fi
  
    ;;
  
    ## edison
    edison* )
    .  /opt/modules/default/init/bash
    module load ncl/6.1.1
    module load nco
    if [ -z "$ESMFBIN_PATH" ]; then
        module use -a /project/projectdirs/ccsm1/modulefiles/edison
        if [ "$gridtype" = "global" ]; then
            mpi=mpi
            mpitype="mpi"
        else
            mpi=uni
            mpitype="mpiuni"
        fi
        module load esmf/6.3.0r-ncdfio-${mpitype}-O
        ESMFBIN_PATH=$ESMF_LIBDIR/../bin
    fi
    if [ -z "$MPIEXEC" ]; then
        MPIEXEC="aprun -n $REGRID_PROC"
    fi
  
    ;;
    ##no other machine currently supported    
    *)
    echo "Machine $hostname NOT recognized"
    ;;
 
esac

#----------------------------------------------------------------------
# Generate the mapping files needed for surface dataset generation
#----------------------------------------------------------------------
 
# Resolve interactive or batch mode command
# NOTE - we only need this here because the machine-specific settings above set
# MPIEXEC depending on if we are on a certain machine, so we need to override
# them if we want interactive mode. This bit of logic should go away in the
# future if we drop machine-specific support from this script.
if [ "$interactive" == "NO" ]; then
    if [ -z "$MPIEXEC" ]; then
        echo "Name of MPI exec to use was NOT set"
        echo "Set the environment variable: MPIEXEC"
        exit 1
    elif ! command -v $MPIEXEC >& /dev/null ; then
        echo "The MPIEXEC pathname given is NOT an executable: $MPIEXEC"
        echo "Set the environment variable: MPIEXEC or run in interactive mode without MPI"
        exit 1
    fi
    echo "Running with MPI via $MPIEXEC"
    mpirun=$MPIEXEC
else
    echo "Running without MPI"
    mpirun=""
fi

# Look for ESMF_RegridWeightGen. If ESMFBIN_PATH is set, then look for binary
# there. Otherwise, assume it exists in PATH and use which to find command.
if [ -z "${ESMFBIN_PATH}" ]; then
    ESMF_REGRID=`which ESMF_RegridWeightGen`
else
    ESMF_REGRID="${ESMFBIN_PATH}/ESMF_RegridWeightGen"
fi
if [ ! -x "$ESMF_REGRID" ]; then
    echo "ESMF_RegridWeightGen does NOT exist in ESMF binary directory: $ESMFBIN_PATH\n"
    echo "Upgrade to a newer version of ESMF with this utility included"
    echo "Set the environment variable: ESMFBIN_PATH"
    exit 1
fi

# Remove previous log files
rm -f PET*.Log

# Now run the mapping for each file, checking that input files exist
# and then afterwards that the output mapping file exists
declare -i nfile=1
until ((nfile>${#INGRID[*]})); do
    echo "Creating mapping file: ${OUTFILE[nfile]}"
    echo "From input grid: ${INGRID[nfile]}"
    echo "For output grid: $GRIDFILE"
    echo " "
 
    # Check that input and output grid files exist
    if [ -z "${INGRID[nfile]}" ] || [ -z "$GRIDFILE" ] || [ -z "${OUTFILE[nfile]}" ]; then
        echo "Either input or output grid or output mapping file is NOT set"
        exit 3
    fi
    if [ ! -f "${INGRID[nfile]}" ]; then
        echo "Input grid file does NOT exist: ${INGRID[nfile]}"
        if [ "$list" != "YES" ]; then
            echo "Set --list argument to output list of grid files to $outfilelist"
            echo "and then download separately using check_input_data (in CIME)."
            exit 2
        fi
    fi
    if [ ! -f "$GRIDFILE" ]; then
        echo "Output grid file does NOT exist: $GRIDFILE"
        exit 3
    fi
 
    # Determine what (if any) large file support is needed. Use the
    # most extreme large file support needed by either the source file
    # or the destination file.
    if [ "$DST_LRGFIL" = "netcdf4" ] || [ "${SRC_LRGFIL[nfile]}" = "netcdf4" ]; then
        lrgfil="--netcdf4"
    elif [ "$DST_LRGFIL" = "64bit_offset" ] || [ "${SRC_LRGFIL[nfile]}" = "64bit_offset" ]; then
        lrgfil="--64bit_offset"
    elif [ "$DST_LRGFIL" = "none" ] && [ "${SRC_LRGFIL[nfile]}" = "none" ]; then
        lrgfil=""
    else
        echo "Unknown LRGFIL type:"
        echo "DST_LRGFIL = $DST_LRGFIL"
        echo "SRC_LRGFIL = ${SRC_LRGFIL[nfile]}"
        exit 4
    fi
 
    # Override file type
    if [ "${output_filetype}" != "" ]; then
        lrgfil="--${output_filetype}"
    fi
 
    # Skip if file already exists
    if [ -f "${OUTFILE[nfile]}" ]; then
        echo "Skipping creation of ${OUTFILE[nfile]} as already exists"
    else
 
        # Build regrid command
        cmd="$mpirun $ESMF_REGRID --ignore_unmapped -s ${INGRID[nfile]} "
        cmd="$cmd -d $GRIDFILE -m conserve -w ${OUTFILE[nfile]}"
        if [ $gridtype = "regional" ]; then
            cmd="$cmd --dst_regional"
        fi
        cmd="$cmd --src_type ${SRC_TYPE[nfile]} ${SRC_EXTRA_ARGS[nfile]} --dst_type $DST_TYPE $DST_EXTRA_ARGS"
        cmd="$cmd $lrgfil"
  
        # Run regrid command
        runcmd $cmd
        if [ "$debug" != "YES" ] && [ ! -f "${OUTFILE[nfile]}" ]; then
            echo "Output mapping file was NOT created: ${OUTFILE[nfile]}"
            exit 6
        fi
  
        # Add some metadata to the file
        HOST=`hostname`
        history="$ESMF_REGRID"
        runcmd "ncatted -a history,global,a,c,"$history"  ${OUTFILE[nfile]}"
        runcmd "ncatted -a hostname,global,a,c,$HOST   -h ${OUTFILE[nfile]}"
        runcmd "ncatted -a logname,global,a,c,$LOGNAME -h ${OUTFILE[nfile]}"
  
        # Check for duplicate mapping weights
        newfile="rmdups_${OUTFILE[nfile]}"
        runcmd "rm -f $newfile"
        runcmd "env MAPFILE=${OUTFILE[nfile]} NEWMAPFILE=$newfile ncl $dir/rmdups.ncl"
        if [ -f "$newfile" ]; then
            runcmd "mv $newfile ${OUTFILE[nfile]}"
        fi
    fi
 
    nfile=nfile+1
done

if [ "${debug}" != "YES" ]; then
    echo ""
    echo "Successfully created needed mapping files for $res"
    echo ""
else
    echo ""
    echo "Script ran to completion; disable debug mode to generate mapping files."
    echo ""
fi

exit 0
