#!/bin/bash

# get arguments
# Need --e3sm_root=

# Also need tempest in PATH for now...


# test mkdurfdat.pl to generate land surface data
# See step 7 in
# https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/872579110/Running+E3SM+on+New+Grids

e3sm_root="default"
test_root="default"
inputdata_root="default"
reference_files="default"

for arg in "$@"
do
case $arg in
    -e=*|--e3sm_root=*)
	e3sm_root="${arg#*=}"
	shift
	;;

    -i=*|--inputdata_root=*)
	inputdata_root="${arg#*=}"
	shift
	;;

    -r=*|--reference_files=*)
	reference_files="${arg#*=}"
	shift
	;;

esac
done

if [[ ${e3sm_root} == "default" ]]; then
    echo "Error: e3sm_root not set" >&2
    exit 1;
fi

if [[ ${inputdata_root} == "default" ]]; then
    echo "Error: inputdata_root not set" >&2
    exit 1;
fi

if [[ ${reference_files} == "default" ]]; then
    echo "Error: reference_files not set" >&2
    exit 1;
fi

output_root=$PWD
cime_root=${e3sm_root}/cime

# Add testing bin to path
PATH=${test_root}/bin:${PATH}

# We will redirect verbose test log output to a file; remove any existing
# versions of this file first
test_log=${PWD}/test.out
rm -f ${test_log}


target_grid=${reference_files}/ne30pg4_scrip.nc
input_topo=${inputdata_root}/atm/cam/topo/USGS-topo-cube3000.nc
output_topo=${PWD}/output.nc

echo "build cube_to_data in ${PWD}/builds ..." >> ${test_log}
mkdir -p builds
cd builds
${cime_root}/CIME/scripts/configure --mpilib mpich --macros-format Makefile >> ${test_log} 2>&1

if [ ! -f .env_mach_specific.sh ]; then
    if [ ! -f .env_mach_specific.sh ]; then
        echo "ERROR running ${cime_root}/CIME/scripts/configure" >&2
        echo "cat ${test_log} for more info" >&2
        exit 1
    fi
fi
cp ${e3sm_root}/components/eam/tools/topo_tool/cube_to_target/* .
sed  "s:^FFLAGS:#FFLAGS:g" Makefile | sed "s:^LDFLAGS:#LDFLAGS:g" > Makefile.tmp
echo "include Macros.make" > Makefile
echo 'FC=${MPIFC}' >> Makefile
echo 'LDFLAGS=${SLIBS}' >> Makefile
echo 'FFLAGS+=-I${NETCDF_PATH}/include' >> Makefile
cat  Makefile.tmp >> Makefile



#

(. .env_mach_specific.sh && export FC=gfortran && export LIB_NETCDF=${NETCDF_FORTRAN_PATH}/lib && export INC_NETCDF=${NETCDF_FORTRAN_PATH/include} && make) >> ${test_log} 2>&1

if [ ! -f cube_to_target ]; then
    echo "ERROR building cube_to_target" >&2
    echo "cat ${test_log} for more info" >&2
    exit 1
fi

echo "Running cube_to_target" >> ${test_log} 2>&1
(. .env_mach_specific.sh && ./cube_to_target --target-grid ${target_grid} --input-topography ${input_topo} --output-topography ${output_topo} )  >> ${test_log}

exit 0
