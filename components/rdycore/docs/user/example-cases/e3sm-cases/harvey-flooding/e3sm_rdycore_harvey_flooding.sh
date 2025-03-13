#!/bin/sh

res=MOS_USRDAT
compset=RMOSGPCC
mach=
rainfall_dataset=
project_id=
N=1
e3sm_dir=../../../../../../../../

# will be set autoamtically
compiler=
rdycore_dir=

display_help() {
  echo "Usage: $0 " >&2
  echo
  echo "   -h, --help                        Display this message"
  echo "   --e3sm-dir                        Path to E3SM-RDycore directory"
  echo "   --mach <pm-cpu|pm-gpu|frontier>   Supported machine name"
  echo "   --frontier-node-type <cpu|gpu>    To run on Frontier CPUs or GPUs"
  echo "   -N, --node  <N>                   Number of nodes (default = 1)"
  echo "   --project <project-id>            Project ID that will charged for the job"
  echo "   --rainfall_dataset <name>  Supported dataset name (i.e. daymet|imerg|mrms|mswep|nldas)"
  return 0
}  

while [ $# -gt 0 ]
do
  case "$1" in
    --mach ) mach="$2"; shift ;;
    --frontier-node-type) frontier_node_type="$2"; shift ;;
    --project-id) project_id="$2"; shift ;;
    --e3sm-dir) e3sm_dir="$2"; shift ;;
    -N | --node) N="$2"; shift ;;
    --rainfall_dataset) rainfall_dataset="$2"; shift ;;
    -h | --help)
      display_help
      exit 0
      ;;
    -*)
      echo "Unsupported argument: $1"
      display_help
      exit 0
      ;;
    *)  break;;    # terminate while loop
  esac
  shift
done

# Check if '--rainfall_dataset <name>' is supported
supported_dataset=0
if [ "$rainfall_dataset" = "daymet" ]; then
    supported_dataset=1
elif [ "$rainfall_dataset" = "imerg" ]; then
    supported_dataset=1
elif [ "$rainfall_dataset" = "mrms" ]; then
    supported_dataset=1
elif [ "$rainfall_dataset" = "mswep" ]; then
    rainfall_dataset=1
elif [ "$rainfall_dataset" = "nldas" ]; then
    supported_dataset=1
fi

if [ "$supported_dataset" -eq 0 ]; then
    echo "The following rainfall dataset is not supported: " $rainfall_dataset
    display_help
    exit 0
fi

#
# Determine the value of:
#  device
#  ntasks
#
device=""
if [ "$mach" == "pm-cpu" ]; then
  data_dir=/global/cfs/projectdirs/m4267/shared/data/harvey
  device="cpu"
  ntasks=$((N*128))
  macros_file_in=${PWD}/gnu_pm-cpu.cmake.pm-cpu-opt-32bit-gcc-11-2-0-fc2888174f5
  macros_file_out=gnu_pm-cpu.cmake
  compiler=gnu
elif [ "$mach" == "pm-gpu" ]; then
  data_dir=/global/cfs/projectdirs/m4267/shared/data/harvey
  device="gpu"
  ntasks=$((N*4))
  macros_file_in=${PWD}/gnugpu_pm-gpu.cmake.pm-gpu-opt-32bit-gcc-11-2-0-fc2888174f5
  macros_file_out=gnugpu_pm-gpu.cmake
  compiler=gnugpu
elif [ "$mach" == "frontier" ]; then

  data_dir=/lustre/orion/cli192/proj-shared/data/harvey
  macros_file_in=${PWD}/gnugpu_frontier.cmake.frontier-gpu-opt-32bit-gcc-11-2-0-fc288817
  macros_file_out=gnugpu_frontier.cmake
  compiler=gnugpu

  # Make sure both CPU and GPU options were not specified
  if [ "$frontier_node_type" == "cpu" ]
  then
    device="cpu"
    ntasks=$((N*56))
  elif [ "$frontier_node_type" == "gpu" ]
  then
    device="gpu"
    mach_prefix=".gpu"
    ntasks=$((N*8))
  else
    echo "Unknown -frontier-node-type $frontier_node_type"
    display_help
    exit 0
  fi

elif [ "$mach" == "" ]; then
  echo "No machine specified via --mach"
  display_help
  exit 0
else
  echo "Unsupported machine specified via --mach $mach"
  display_help
  exit 0
fi

domain_file=domain_Dlnd_2926x1_c240507.nc
domain_path=${data_dir}/e3sm/Turning_30m
RDYCORE_YAML_FILE=${PWD}/Turning_30m.critical_flow_bc.yaml
RDYCORE_IC_FILE=${data_dir}/Turning_30m/solution_219.int32.dat
RDYCORE_MESH_FILE=${data_dir}/Turning_30m/Turning_30m_with_z.updated.with_sidesets.exo
RDYCORE_BIN_MAP=${data_dir}/e3sm/Turning_30m/map_MOSART_to_RDycore_Turning_30_2926532x1.bin
frivinp_rtm=${data_dir}/e3sm/Turning_30m/MOSART_Dlnd_2926x1_c240507.nc
map_file=${data_dir}/e3sm/Turning_30m/map_${rainfall_dataset}_to_Dlnd_2926x1.nc

src_dir=${e3sm_dir}
case_dir=${src_dir}/cime/scripts

# Compile RDycore

rdycore_dir=$e3sm_dir/externals/rdycore/

cd $rdycore_dir

source config/set_petsc_settings.sh --mach $mach --config 3


if [ ! -d "$rdycore_dir/build-$PETSC_ARCH" ]
then
  echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
  echo "The following expected RDycore build directory not found:$rdycore_dir/build-$PETSC_ARCH "
  echo "So, attempting to build RDycore."
  echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

  cmake -S . -B build-$PETSC_ARCH -DCMAKE_INSTALL_PREFIX=$PWD/build-$PETSC_ARCH
fi

# Build 
cd build-$PETSC_ARCH
make -j4 install

# Determine the case name
cd $src_dir
git_hash=`git log -n 1 --format=%h`
case_name=Dlnd.${mach}.${device}.Turning_30m_${rainfall_dataset}_${git_hash}.NTASKS_${ntasks}.${compiler}.`date "+%Y-%m-%d"`

# Create the case
cd ${src_dir}/cime/scripts
case_dir=${src_dir}/cime/scripts
./create_newcase -case ${case_dir}/${case_name} \
--res ${res} --mach ${mach} --compiler ${compiler} --compset ${compset} --project ${project_id}

# Change to case dir
cd ${case_dir}/${case_name}

./xmlchange CALENDAR=NO_LEAP

./xmlchange LND_DOMAIN_FILE=$domain_file
./xmlchange ATM_DOMAIN_FILE=$domain_file
./xmlchange LND_DOMAIN_PATH=$domain_path
./xmlchange ATM_DOMAIN_PATH=$domain_path

./xmlchange DATM_CLMNCEP_YR_END=1979
./xmlchange DATM_CLMNCEP_YR_START=1979
./xmlchange DATM_CLMNCEP_YR_ALIGN=1979
./xmlchange DLND_CPLHIST_YR_START=1979
./xmlchange DLND_CPLHIST_YR_END=1979
./xmlchange DLND_CPLHIST_YR_ALIGN=1979
./xmlchange RUN_STARTDATE=1979-01-01

./xmlchange NTASKS=$ntasks
./xmlchange STOP_N=2,STOP_OPTION=nhours
./xmlchange PIO_NETCDF_FORMAT=classic
./xmlchange ROF_NCPL=24

cat >> user_nl_mosart << EOF
frivinp_rtm    = '$frivinp_rtm'
routingmethod  = 1
frivinp_mesh   = ''
rtmhist_nhtfrq = 1
rtmhist_mfilt  = 240
decomp_option  = 'rdycore'
EOF

cat >> user_nl_dlnd << EOF
dtlimit=2.0e0
mapalgo = "nn"
mapread = "$map_file"
EOF

cp ${data_dir}/e3sm/dlnd//user_dlnd.streams.txt.lnd.${rainfall_dataset} user_dlnd.streams.txt.lnd.gpcc

./case.setup

cp ${macros_file_in} cmake_macros/${macros_file_out}

#petsc_libs=`pkg-config --libs --static $PETSC_DIR/$PETSC_ARCH/lib/pkgconfig/petsc.pc`
#sed -i "s/PLACEHOLDER_PETSC_LIBS/${petsc_libs//\//\\/}/g" cmake_macros/${macros_file_out}
sed -i "s/PLACEHOLDER_E3SM_DIR/${e3sm_dir//\//\\/}/g" cmake_macros/${macros_file_out}
sed -i "s/PLACEHOLDER_PETSC_DIR/${PETSC_DIR//\//\\/}/g" cmake_macros/${macros_file_out}
sed -i "s/PLACEHOLDER_PETSC_ARCH/${PETSC_ARCH}/g" cmake_macros/${macros_file_out}

if [ "$mach" == "pm-cpu" ]; then
  ./xmlchange run_exe="\${EXEROOT}/e3sm.exe -ceed /cpu/self -log_view"
elif [ "$mach" == "pm-gpu" ]; then
  ./xmlchange run_exe="-G4 \${EXEROOT}/e3sm.exe -ceed /gpu/cuda -dm_vec_type cuda -use_gpu_aware_mpi 1 -log_view -log_view_gpu_time"
elif [ "$mach" == "frontier" ]; then
  # Make sure both CPU and GPU options were not specified
  if [ "$frontier_node_type" == "cpu" ]
  then
    ./xmlchange run_exe="\${EXEROOT}/e3sm.exe -ceed /cpu/self -log_view"
  elif [ "$frontier_node_type" == "gpu" ]
  then
    ./xmlchange run_exe="\${EXEROOT}/e3sm.exe -ceed /gpu/hip -dm_vec_type hip -use_gpu_aware_mpi 0 -log_view -log_view_gpu_time"
  fi
else
  echo "Unsupported machine specified via --mach $mach"
  display_help
  exit 0
fi

rundir=`./xmlquery RUNDIR --value`

cd $rundir

cp $RDYCORE_YAML_FILE    rdycore.yaml
ln -s $RDYCORE_IC_FILE   solution_219.dat
ln -s $RDYCORE_MESH_FILE .
ln -s $RDYCORE_BIN_MAP   map_MOSART_to_RDycore.bin

cd ${case_dir}/${case_name}

./case.build

