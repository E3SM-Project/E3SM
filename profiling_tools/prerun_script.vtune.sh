#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ $1 = "CPPFLAGS" ]; then
  echo "-g -I/global/project/projectdirs/acme/mjac/api_itt/ -DPROFILE_VTUNE_NEW"
  exit
elif [ $1 = "LDFLAGS" ]; then
  echo "-dynamic -L/global/project/projectdirs/acme/mjac/api_itt/ -litt -L\${VTUNE_AMPLIFIER_XE_2018_DIR}/lib64/ -littnotify"
  exit
elif [ $1 = "README" ]; then
  echo "This is a one step profiling tool."
  echo "Once data have been collected by launching a job using ./case.submit,"
  echo "the post_process action has to be called"
  echo "before the results can be visualized in amplxe-gui."
  exit
elif [ $1 = "POSTPROCESS" ]; then
  if [ "$#" -ne 2 ]; then
      echo "Illegal number of parameters"
      echo "usage: $0 POSTPROCESS /path/to/vtune_results_dir"
      exit
  elif [ ! -d "$2" ]; then
      echo "$2 is not a valid directory."
      exit
  fi
 
 
  module load vtune
  
  result_dir=$2
  case_dir=$(dirname $(dirname $result_dir))
  exe=$case_dir/bld/e3sm.exe
  cp $exe $exe\ \(deleted\)
  build_dir=$(dirname $exe)
  
  E3SM_SOURCE=$(pushd $case_dir 1>/dev/null; ./xmlquery --value SRCROOT; popd 1>/dev/null)
  
  amplxe-cl -finalize -r $result_dir \
  -search-dir=/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64 \
  -search-dir=/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64_lin \
  -search-dir=/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64_lin_mic \
  -search-dir=/opt/intel/compilers_and_libraries/linux/compiler/lib/intel64 \
  -search-dir=/opt/intel/compilers_and_libraries/linux/compiler/lib/intel64_lin \
  -search-dir=/opt/intel/compilers_and_libraries/linux/compiler/lib/intel64_lin_mic \
  -search-dir=/opt/intel/vtune_amplifier_xe/lib64/runtime \
  -search-dir=/opt/intel/vtune_amplifier_xe/bin64  \
  -search-dir=/opt/intel/vtune_amplifier_xe/lib64 \
  -search-dir=/usr/lib64 \
  -search-dir=/lib64 \
  -search-dir=/opt/cray/pe/hdf5-parallel/default/INTEL/15.0/lib \
  -search-dir=/opt/cray/pe/mpt/default/gni/mpich-intel/16.0/lib \
  -search-dir=/opt/cray/pe/pmi/default/lib64 \
  -search-dir=/opt/cray/xpmem/default/lib64 \
  -search-dir=/opt/cray/alps/default/lib64 \
  -search-dir=/opt/cray/ugni/default/lib64 \
  -search-dir=/opt/cray/udreg/default/lib64 \
  -search-dir=$build_dir \
  -search-dir=$build_dir/lib \
  -search-dir=$build_dir \
  -source-search-dir=$build_dir \
  -source-search-dir=$E3SM_SOURCE
  exit
elif [ $1 = "PROFILE" ]; then
  module load vtune
  analysis=advanced_hotspots
  tool=amplxe-cl
  resdir=vtune.${analysis}.${LID}
  echo "0 $tool -collect ${analysis} -finalization-mode=none -data-limit=0 -start-paused -r ${resdir} -- `./xmlquery --value EXEROOT`/e3sm.exe" > run/mpmd.${LID}.conf;
  echo "1-$( echo `./xmlquery --value TOTALPES`-1 | bc) `./xmlquery --value EXEROOT`/e3sm.exe" >> run/mpmd.${LID}.conf  
fi

