#!/bin/bash


############################################################################
## MAKE SURE WE HAVE WHAT WE NEED
############################################################################
function usage {
  printf "Usage: ./cmakescript.sh 2dfile.nc 3dfile.nc\n\n"
  printf "You must specify NCHOME and NFHOME environment variables to specify\n"
  printf "where the NetCDF libraries are located\n\n"
  printf "NetCDF binaries must include ncdump, nf-config, and nc-config\n\n"
  printf "You can also define FFLAGS to control optimizations and NCRMS \n"
  printf "to reduce the number of CRM samples and the runtime of the tests.\n\n"
  printf "./cmakescript.sh [-h|--help] for this message\n\n"
}
if [[ "$1" == "" || "$2" == "" ]]; then
  printf "Error: missing 2d and / or 3d NetCDF File parameters\n"
  usage
  exit -1
fi
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  usage
  exit 0
fi
if [[ "$NCHOME" == "" ]]; then
  printf "Error: NCHOME environment variable not set\n"
  printf "set NCHOME with the path to the NetCDF C installation\n\n"
  usage
  exit -1
fi
if [[ "$NFHOME" == "" ]]; then
  printf "Error: NFHOME environment variable not set\n"
  printf "set NFHOME with the path to the NetCDF Fortran installation\n\n"
  usage
  exit -1
fi


############################################################################
## GRAB 2D DATA FROM THE NETCDF FILE
############################################################################
NX=`$NCHOME/bin/ncdump -h $1  | grep "crm_nx =" | awk '{print $3}'`
NY=`$NCHOME/bin/ncdump -h $1  | grep "crm_ny =" | awk '{print $3}'`
NZ=`$NCHOME/bin/ncdump -h $1  | grep "crm_nz =" | awk '{print $3}'`
NX_RAD=`$NCHOME/bin/ncdump -h $1  | grep "crm_nx_rad =" | awk '{print $3}'`
NY_RAD=`$NCHOME/bin/ncdump -h $1  | grep "crm_ny_rad =" | awk '{print $3}'`
DX=1000
DT=5
NCRMS_FILE=`ncdump -h $1 | grep UNLIMITED | awk '{print $6}' | cut -d '(' -f 2`
if [[ $NY -eq 1 ]]; then
  YES3D=0
else
  echo "Error: 3D file specified as the 2D file\n\n"
  usage
fi
PLEV=`$NCHOME/bin/ncdump -h $1  | grep "nlev =" | awk '{print $3}'`
if [[ "$NCRMS" != "" ]]; then
  NCRMS2D=$NCRMS
  if [[ $NCRMS -gt $NCRMS_FILE ]]; then
    printf "WARNING: NCRMS environment variable is larger than the available samples in the 2D input NetCDF file\n\n"
    NCRMS2D=$NCRMS_FILE
  fi
else
  NCRMS2D=$NCRMS_FILE
fi

DEFS2D=" -DNCRMS=$NCRMS2D -DCRM -DCRM_NX=$NX -DCRM_NY=$NY -DCRM_NZ=$NZ -DCRM_NX_RAD=$NX_RAD -DCRM_NY_RAD=$NY_RAD -DCRM_DT=$DT -DCRM_DX=$DX -DYES3DVAL=$YES3D -DPLEV=$PLEV -Dsam1mom -DMMF_STANDALONE -DHAVE_MPI"
printf "2D Defs: $DEFS2D\n\n"


############################################################################
## GRAB 3D DATA FROM THE NETCDF FILE
############################################################################
NX=`$NCHOME/bin/ncdump -h $2  | grep "crm_nx =" | awk '{print $3}'`
NY=`$NCHOME/bin/ncdump -h $2  | grep "crm_ny =" | awk '{print $3}'`
NZ=`$NCHOME/bin/ncdump -h $2  | grep "crm_nz =" | awk '{print $3}'`
NX_RAD=`$NCHOME/bin/ncdump -h $2  | grep "crm_nx_rad =" | awk '{print $3}'`
NY_RAD=`$NCHOME/bin/ncdump -h $2  | grep "crm_ny_rad =" | awk '{print $3}'`
DX=1000
DT=5
NCRMS_FILE=`ncdump -h $2 | grep UNLIMITED | awk '{print $6}' | cut -d '(' -f 2`
if [[ $NY -eq 1 ]]; then
  echo "Error: 2D file specified as the 3D file\n\n"
  usage
else
  YES3D=1
fi
PLEV=`$NCHOME/bin/ncdump -h $2  | grep "nlev =" | awk '{print $3}'`
if [[ "$NCRMS" != "" ]]; then
  if [[ $NCRMS -gt $NCRMS_FILE ]]; then
    printf "ERROR: NCRMS environment variable is larger than the available samples in the 3D input NetCDF file\n\n"
    exit -1
  fi
  NCRMS3D=$NCRMS
else
  NCRMS3D=$NCRMS_FILE
fi

DEFS3D=" -DNCRMS=$NCRMS3D -DCRM -DCRM_NX=$NX -DCRM_NY=$NY -DCRM_NZ=$NZ -DCRM_NX_RAD=$NX_RAD -DCRM_NY_RAD=$NY_RAD -DCRM_DT=$DT -DCRM_DX=$DX -DYES3DVAL=$YES3D -DPLEV=$PLEV -Dsam1mom -DMMF_STANDALONE -DHAVE_MPI"
printf "3D Defs: $DEFS3D\n\n"


############################################################################
## CLEAN UP THE PREVIOUS BUILD
############################################################################
rm -rf CMakeCache.txt CMakeFiles cmake_install.cmake CTestTestfile.cmake Makefile fortran.exe cpp.exe cpp2d cpp3d fortran2d fortran3d


############################################################################
## SYMLINK INPUT FILES INTO EXECUTABLE DIRECTORIES
############################################################################
mkdir fortran2d
mkdir fortran3d
mkdir cpp2d    
mkdir cpp3d    
cd fortran2d   ; ln -s ../$1 ./input.nc
cd ../fortran3d; ln -s ../$2 ./input.nc
cd ../cpp2d    ; ln -s ../$1 ./input.nc
cd ../cpp3d    ; ln -s ../$2 ./input.nc
cd ..

### link non-standard data file
# rm *3d/input.nc
# ln -s /gpfs/alpine/cli115/proj-shared/hannah6/crm_standalone_data/crmdata_3d_bug_combined.nc /ccs/home/hannah6/E3SM/E3SM_SRC1/components/eam/src/physics/crm/samxx/test/build/fortran3d/input.nc
# ln -s /gpfs/alpine/cli115/proj-shared/hannah6/crm_standalone_data/crmdata_3d_bug_combined.nc /ccs/home/hannah6/E3SM/E3SM_SRC1/components/eam/src/physics/crm/samxx/test/build/cpp3d/input.nc


############################################################################
## GET THE NETCDF LINKING FLAGS
############################################################################
NCFLAGS="`$NFHOME/bin/nf-config --flibs` `$NCHOME/bin/nc-config --libs`"
printf "NetCDF Flags: $NCFLAGS\n\n"


############################################################################
## RUN THE CONFIGURE
############################################################################
FFLAGS="$FFLAGS -I$NCHOME/include -I$NFHOME/include"
CXXFLAGS="$CXXFLAGS -I$NCHOME/include -I$NFHOME/include"
CUDAFLAGS="$CUDAFLAGS ${CUDA_ARCH}"

printf "FFLAGS: $FFLAGS\n\n"
printf "CXXFLAGS: $CXXFLAGS\n\n"
printf "CUDAFLAGS: $CUDAFLAGS\n\n"

echo cmake                          \
  -DCMAKE_Fortran_FLAGS="$FFLAGS"   \
  -DCMAKE_CXX_FLAGS="$CXXFLAGS"     \
  -DNCFLAGS="$NCFLAGS"              \
  -DDEFS2D="$DEFS2D"                \
  -DDEFS3D="$DEFS3D"                \
  -DCUDA_FLAGS="$CUDAFLAGS"         \
  -DYAKL_HOME=${YAKL_HOME}          \
  -DYAKL_CUB_HOME=${YAKL_CUB_HOME}  \
  -DARCH="${ARCH}"                  \
  ..

cmake                               \
  -DCMAKE_Fortran_FLAGS="$FFLAGS"   \
  -DCMAKE_CXX_FLAGS="$CXXFLAGS"     \
  -DNCFLAGS="$NCFLAGS"              \
  -DDEFS2D="$DEFS2D"                \
  -DDEFS3D="$DEFS3D"                \
  -DCUDA_FLAGS="$CUDAFLAGS"         \
  -DYAKL_HOME=${YAKL_HOME}          \
  -DYAKL_CUB_HOME=${YAKL_CUB_HOME}  \
  -DARCH="${ARCH}"                  \
  ..


