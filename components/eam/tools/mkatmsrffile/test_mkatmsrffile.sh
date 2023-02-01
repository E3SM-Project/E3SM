#!/bin/bash

display_help() {
    echo "Usage: $0 " >&2
    echo
    echo "  -e, --e3sm_root <e3sm_root_directory>   Specify location of E3SM"
    echo "  -h, --help                              Display this message"
    echo "  -i, --inputdata_root <data_directory>   Specify location of climate inputdata"
    echo "  -r, --reference_files <ref_directory>   Specify location where files"
    
    echo "                                          1x1d.nc, ne30np4_pentagons.091226.nc,"
    echo "                                          map_1x1_to_ne30np4_aave.nc"
    echo "                                          are located" 
    echo
    echo "NOTE: requires tempestremap and ESMF tools to be in PATH environment variable"
}    

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

if [[ ${reference_files} == "default" ]]; then
    echo "Error: reference_files not set" >&2
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


srf_file=${reference_files}/1x1d.nc
atm_file=${reference_files}/ne30np4_pentagons.091226.nc
land_file=${inputdata_root}/atm/cam/chem/trop_mozart/dvel/regrid_vegetation.nc
soilw_file=${inputdata_root}/atm/cam/chem/trop_mozart/dvel/clim_soilw.nc
srf2atm_file=${reference_files}/map_1x1_to_ne30np4_aave.nc

for i in ${srf_file} ${atm_file} ${land_file} ${soilw_file} ${srf2atm_file}
do
    if [ ! -f $i ]; then
	echo "Error: file ${i} not found" >&2
	exit 1
    fi
done

output_file=${PWD}/atmsrf_ne30np4.nc

echo "build mkatmsrrfile in ${PWD}/builds ..." >> ${test_log}
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
cp ${e3sm_root}/components/eam/tools/mkatmsrffile/* .

# Edit Makefile to use variables created by configure
sed  "s:^FFLAGS:#FFLAGS:g" Makefile | sed "s:^INC:#INC^:g" | sed "s:^LIB:#LIB:g"  > Makefile.tmp
echo "include Macros.make" > Makefile
echo 'FC=${MPIFC}' >> Makefile
echo 'LIB=${SLIBS}' >> Makefile
echo 'FFLAGS+=-I${NETCDF_PATH}/include' >> Makefile
cat  Makefile.tmp >> Makefile

cat <<EOF > nml_atmsrf
&input
srfFileName =  '${srf_file}'
atmFileName = '${atm_file}'
landFileName = '${land_file}' 
soilwFileName = '${soilw_file}'
srf2atmFmapname = '${srf2atm_file}'
outputFileName = '${output_file}'
/

EOF



#

(. .env_mach_specific.sh && export FC=gfortran && export LIB_NETCDF=${NETCDF_FORTRAN_PATH}/lib && export INC_NETCDF=${NETCDF_FORTRAN_PATH/include} && make) >> ${test_log} 2>&1
if [ ! -f mkatmsrffile ]; then
    echo "ERROR building mkatmsrffile" >&2
    echo "cat ${test_log} for more info" >&2
    exit 1
fi

rm -f ${output_file}
echo "Running mkatmsrffile" >> ${test_log} 2>&1
(. .env_mach_specific.sh && ./mkatmsrffile )  >> ${test_log}

if [ ! -f ${output_file} ]; then
    echo "Error: file ${i} not found" >&2
    exit 1
else
    echo "output file ${i} created" >> ${test_log} 2>&1
fi


exit 0
