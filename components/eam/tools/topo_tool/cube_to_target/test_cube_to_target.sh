#!/bin/bash

display_help() {
    echo "Usage: $0 " >&2
    echo
    echo "  -e, --e3sm_root <e3sm_root_directory>   Specify location of E3SM"
    echo "  -h, --help                              Display this message"
    echo "  -i, --inputdata_root <data_directory>   Specify location of climate inputdata"
    echo
}    
 
    
# get arguments
# Need --e3sm_root=
#      --reference_files=
#      --inputdata_root=

e3sm_root="default"
test_root="default"
inputdata_root="default"

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

    -*)
	display_help
	exit 1;
	;;
    
    -h|--help)
	display_help
	exit 0;
	;;

esac
done

if [[ ${e3sm_root} == "default" ]]; then
    echo "Error: e3sm_root not set" >&2
    display_help
    exit 1;
fi

if [[ ${inputdata_root} == "default" ]]; then
    echo "Error: inputdata_root not set" >&2
    display_help
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


generatecsmesh=`which GenerateCSMesh`
generatevolumetricmesh=`which GenerateVolumetricMesh`
convertmeshtoscrip=`which ConvertMeshToSCRIP`

if [ "${generatecsmesh}x" == "x" ]; then
    echo "ERROR: tempestremap tool GenerateCSMesh not found in PATH" >&2
    echo "cat ${test_log} for more info" >&2
    exit 1
fi

if [ "${generatevolumetricmesh}x" == "x" ]; then
    echo "ERROR: tempestremap tool GenerateVolumetricMesh not found in PATH" >&2
    echo "cat ${test_log} for more info" >&2
    exit 1
fi

if [ "${convertmeshtoscrip}x" == "x" ]; then
    echo "ERROR: tempestremap tool ConvertMeshToScrip not found in PATH" >&2
    echo "cat ${test_log} for more info" >&2
    exit 1
fi


meshfile=ne30.g
gridfile=ne30pg4.g
scripfile=ne30pg4_scrip.nc
target_grid=${reference_files}/ne30pg4_scrip.nc
input_topo=${inputdata_root}/atm/cam/topo/USGS-topo-cube3000.nc
output_topo=${PWD}/output.nc

echo "Running ${generatecsmesh}" >> ${test_log} 2>&1
(${generatecsmesh} --alt --res 30 --file ${meshfile}) >> ${test_log} 2>&1
if [ ! -f ${meshfile} ]; then
    echo "ERROR: GenerateCSMesh: no ${meshfile} file created" >&2
    echo "cat ${test_log} for more info" >&2
    exit 1
fi

echo "Running ${generatevolumetricmesh}" >> ${test_log} 2>&1
(${generatevolumetricmesh} --in ${meshfile} --out ${gridfile}) >> ${test_log} 2>&1
if [ ! -f ${gridfile} ]; then
    echo "ERROR: GenerateVolumetricMesh: no ${gridfile} file created" >&2
    echo "cat ${test_log} for more info" >&2
    exit 1
fi

echo "Running ${convertmeshtoscrip}" >> ${test_log} 2>&1
(${convertmeshtoscrip} --in ${meshfile} --out ${scripfile}) >> ${test_log} 2>&1
if [ ! -f ${scripfile} ]; then
    echo "ERROR: ConvertMeshToSCRIP: no ${scripfile} file created" >&2
    echo "cat ${test_log} for more info" >&2
    exit 1
fi


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
#Edit Makefile to use macros defined by configure in Macros.make
sed  "s:^FFLAGS:#FFLAGS:g" Makefile | sed "s:^LDFLAGS:#LDFLAGS:g" > Makefile.tmp
echo "include Macros.make" > Makefile
echo 'FC=${MPIFC}' >> Makefile
echo 'LDFLAGS=${SLIBS}' >> Makefile
echo 'FFLAGS+=-I${NETCDF_PATH}/include' >> Makefile
cat  Makefile.tmp >> Makefile

# Compile
(. .env_mach_specific.sh && export FC=gfortran && export LIB_NETCDF=${NETCDF_FORTRAN_PATH}/lib && export INC_NETCDF=${NETCDF_FORTRAN_PATH/include} && make) >> ${test_log} 2>&1

if [ ! -f cube_to_target ]; then
    echo "ERROR building cube_to_target" >&2
    echo "cat ${test_log} for more info" >&2
    exit 1
fi

echo "Running cube_to_target" >> ${test_log} 2>&1
(. .env_mach_specific.sh && ./build/cube_to_target --target-grid ${target_grid} --input-topography ${input_topo} --output-topography ${output_topo} )  >> ${test_log}

exit 0
