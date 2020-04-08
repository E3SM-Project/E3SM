#!/bin/bash

# get arguments
# Need --cime_root=
#      --inputdata_root=

cime_root="default"
inputdata_root="default"

for arg in "$@"
do
case $arg in
    -c=*|--cime_root=*)
    cime_root="${arg#*=}"
    shift
    ;;

    -i=*|--inputdata_root=*)
    inputdata_root="${arg#*=}"
    shift
    ;;
esac
done
if [[ ${inputdata_root} == "default" ]]; then
    echo "Error: inputdata_root not set" >&2
    exit 1;
fi
if [[ ${cime_root} == "default" ]]; then
    echo "Error: cime_root not set" >&2
    exit 1;
fi

# Setup the test case. gen_domain takes as inputs a mapping file and the names
# of the ocean and land grids. Note that we could probably get these by parsing
# the specified baseline filenames or global netcdf attributes, so that we
# could maybe simplify this to specifying just ocn_baseline and lnd_baseline,
# and then figuring out which mapping file and grid names we need.
ocn_name=oQU240
lnd_name=ne4np4
mapping_file=${inputdata_root}/cpl/gridmaps/oQU240/map_oQU240_to_ne4np4_aave.160614.nc
ocn_baseline=${inputdata_root}/share/domains/domain.ocn.ne4np4_oQU240.160614.nc
lnd_baseline=${inputdata_root}/share/domains/domain.lnd.ne4np4_oQU240.160614.nc

# We will redirect verbose test log output to a file; remove any existing
# versions of this file first
test_root=${PWD}
test_log=${test_root}/test.out
rm -f ${test_log}

# Build the gen_domain executable
echo "" >> ${test_log}
echo "Building gen_domain in ${PWD}/builds ..." >> ${test_log}
mkdir -p builds
cd builds
${cime_root}/tools/configure --macros-format Makefile --mpilib mpi-serial >> ${test_log} 2>&1
if [ ! -f .env_mach_specific.sh ]; then
    # try without mpi-serial flag
    echo "ERROR running ${cime_root}/tools/configure" >&2
    echo "It's possible mpi-serial doesn't work on this machine. Trying again with default" >&2
    ${cime_root}/tools/configure --clean --macros-format Makefile >> ${test_log} 2>&1
    if [ ! -f .env_mach_specific.sh ]; then
        echo "ERROR running ${cime_root}/tools/configure" >&2
        echo "cat ${test_log} for more info" >&2
        exit 1
    fi
fi

cp ${cime_root}/tools/mapping/gen_domain_files/src/* .
(. .env_mach_specific.sh && make) >> ${test_log} 2>&1
if [ ! -f ../gen_domain ]; then
    echo "ERROR building gen_domain" >&2
    echo "cat ${test_log} for more info" >&2
    exit 1
fi

# Build the cprnc executable (for comparison of netcdf files)
echo "" >> ${test_log}
echo "Building cprnc in ${PWD}/builds ..." >> ${test_log}
cp ${cime_root}/tools/cprnc/*.F90 .
cp ${cime_root}/tools/cprnc/Makefile .
cp ${cime_root}/tools/cprnc/Depends .
cp ${cime_root}/tools/cprnc/*.in .
(. .env_mach_specific.sh && make GENF90=${cime_root}/src/externals/genf90/genf90.pl) >> ${test_log} 2>&1
if [ ! -f cprnc ]; then
    echo "ERROR building cprnc" >&2
    echo "cat ${test_log} for more info" >&2
    exit 1
fi
cd ${test_root}

# Run example gen_domain code on test case
(. builds/.env_mach_specific.sh && ./gen_domain -m ${mapping_file} -o ${ocn_name} -l ${lnd_name}) >> ${test_log} 2>&1

# Compare outputs from test case against baselines
datestring=`date +'%y%m%d'`
for baseline in ${ocn_baseline} ${lnd_baseline}; do
    # Find file that matches prefix of specified baseline; do this by stripping
    # off last two tokens from baseline filename (file extension and datestring)
    # and adding in datestring for current day and .nc file extension.
    testfile=`basename ${baseline} | rev | cut -d. -f3- | rev`.${datestring}.nc
    if [ ! -f ${testfile} ]; then
	echo "ERROR: ${testfile} not generated" >&2
	echo "cat ${test_log} for more info" >&2
	exit 1
    fi
    # Compare against baseline and print report from cprnc comparison
    echo "Comparing $testfile against ${baseline}..."
    (. builds/.env_mach_specific.sh && ./builds/cprnc -m ${testfile} ${baseline}) >> ${test_log} 2>&1

    # Check results
    last=`tail -n3 ${test_log}`
    if [[ ${last} =~ "STOP" ]]; then
	echo ${last} >&2
	echo "Error running cprnc" >&2
	echo "cat ${test_log} for more info" >&2
	exit 1
    fi
    if [[ ${last} =~ "DIFFERENT" ]]; then
	echo ${last} >&2
	echo ${baseline} DIFFERENT FROM ${testfile} >&2
	echo "cat ${test_log} for more info" >&2
	exit 1
    fi
    if ! [[ ${last} =~ "IDENTICAL" ]]; then
	echo ${last} >&2
	echo "undetermined output from cprnc" >&2
	echo "cat ${test_log} for more info" >&2
	exit 1
    fi
done

# Exit gracefully. Alternatively, we could return a non-zero exit status if any
# of the cprnc comparisons failed.
exit 0
