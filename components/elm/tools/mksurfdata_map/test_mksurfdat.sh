#!/bin/bash

# get arguments
# Need --e3sm_root=
#      --ESMFBIN=
#      --inputdata_root=

# Also need tempest in PATH for now...


# test mkdurfdat.pl to generate land surface data
# See step 7 in
# https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/872579110/Running+E3SM+on+New+Grids

e3sm_root="default"
test_root="default"
esmfbin="default"
inputdata_root="default"

for arg in "$@"
do
case $arg in
    -e=*|--e3sm_root=*)
	e3sm_root="${arg#*=}"
	shift
	;;

    -E=*|--ESMFBIN=*)
	esmfbin="${arg#*=}"
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
if [[ ${e3sm_root} == "default" ]]; then
    echo "Error: e3sm_root not set" >&2
    exit 1;
fi
if [[ ${esmfbin} == "default" ]]; then
    echo "Error: ESMFBIN not set" >&2
    exit 1;
fi

export ESMFBIN_PATH=${esmfbin}

output_root=$PWD
cime_root=${e3sm_root}/cime

# Add testing bin to path
PATH=${test_root}/bin:${PATH}

# We will redirect verbose test log output to a file; remove any existing
# versions of this file first
test_log=${PWD}/test.out
rm -f ${test_log}

###########
###  1  ### Create mapping files for each land surface type if needed
###########

# a)
# Obtain or generate a target grid file in SCRIP format. For these example, we will use a ne1024pg2 grid file,
# which we will need to create (note that most np4 grid files can be found within the inputdata repository, for
# example, the ne1024np4 grid file is at
#   https://web.lcrc.anl.gov/public/e3sm/mapping/grids/ne1024np4_scrip_c20191023.nc)

generatecsmesh=`which GenerateCSMesh`
generatevolumetricmesh=`which GenerateVolumetricMesh`
convertexodustoscrip=`which ConvertMeshToSCRIP`
mksurfdata_map=`which mksurfdata_map` # TODO compile this

#ne1024g=${cime_root}/tools/mapping/tests/reference_files/ne1024.g
# These files will be created
meshfile=ne4.g
gridfile=ne4pg2.g
scripfile=ne4pg2_scrip.nc

echo "Running ${generatecsmesh}" >> ${test_log} 2>&1
(${generatecsmesh} --alt --res 30 --file ${meshfile}) >> ${test_log} 2>&1
if [ ! -f ${meshfile} ]; then
    echo "ERROR: GenerateVolumetricMesh: no ${meshfile} file created" >&2
    echo "cat ${test_log} for more info" >&2
    exit 1
fi

echo "Running ${generatevolumetricmesh}" >> ${test_log} 2>&1
(${generatevolumetricmesh} --in ${meshfile} --out ${gridfile} --np 2 --uniform) >> ${test_log} 2>&1
if [ ! -f ${gridfile} ]; then
    echo "ERROR: GenerateVolumetricMesh: no ${gridfile} file created" >&2
    echo "cat ${test_log} for more info" >&2
    exit 1
fi

echo "Running ${convertexodustoscrip}" >> ${test_log} 2>&1
(${convertexodustoscrip} --in ${meshfile} --out ${scripfile}) >> ${test_log} 2>&1
if [ ! -f ${scripfile} ]; then
    echo "ERROR: ConvertExodusToSCRIP: no ${scripfile} file created" >&2
    echo "cat ${test_log} for more info" >&2
    exit 1
fi



# b)
# Get list of input grid files for each land surface input data file. This is done by running the
# components/clm/tools/shared/mkmapdata/mkmapdata.sh script in debug mode to output a list of needed
# files (along with the commands that will be used to generate each map file; also make sure GRIDFILE
# is set to the SCRIP file from the above step): 
mkmapdata=${e3sm_root}/components/elm/tools/mkmapdata/mkmapdata.sh
if [ ! -f ${makemapdata} ]; then
   echo "ERROR: mkmapdata.sh not found"
   exit 1
fi


# Gen env_mach_specific
${cime_root}/CIME/scripts/configure --mpilib mpich --macros-format Makefile >> ${test_log} 2>&1

echo "Running ${mkmapdata}" >> ${test_log} 2>&1
(. .env_mach_specific.sh && set -x && ${mkmapdata} --gridfile ${scripfile} --inputdata-path ${inputdata_root} --res ne4pg2 --gridtype global --output-filetype 64bit_offset) >> ${test_log} 2>&1


set +x
#awk '/ingrid/ {print $3}' clm.input_data_list  > files.txt
missingfile="no"
#for filename in `cat files.txt`; do
#    if [ ! -f ${filename} ]; then
#	echo "${filename} not found"
#	missingfile="yes"
#    fi
#done

if [ ${missingfile} == "yes" ]; then
    echo "Error: missing files. Try downloading from"
    echo "https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata"
    echo "or"
    echo "https://web.lcrc.anl.gov/public/e3sm/inputdata"
    exit 1
fi

# d) Create mapping file

##############
####  3  #####
##############

echo "build mksurfdata_map" >> ${test_log}
if [ ! -f .env_mach_specific.sh ]; then
    # try without mpi-serial flag
    echo "ERROR running ${cime_root}/CIME/scripts/configure" >&2
    echo "It's possible mpi-serial doesn't work on this machine. Trying again with default" >&2
    ${cime_root}/CIME/scripts/configure --clean >> ${test_log} 2>&1
    (. .env_mach_specific.sh && ${cime_root}/CIME/scripts/configure --macros-format Makefile) >> ${test_log} 2>&1
    if [ ! -f .env_mach_specific.sh ]; then
        echo "ERROR running ${cime_root}/CIME/scripts/configure" >&2
        echo "cat ${test_log} for more info" >&2
        exit 1
    fi
else
    (. .env_mach_specific.sh && ${cime_root}/CIME/scripts/configure --macros-format Makefile --mpilib mpi-serial) >> ${test_log} 2>&1
fi

cp ${e3sm_root}/components/elm/tools/mksurfdata_map/src/* .
echo 'export USER_FC=$FC' >> .env_mach_specific.sh
echo 'export USER_CC=$CC' >> .env_mach_specific.sh
echo 'export LIB_NETCDF=$NETCDF_PATH/lib' >> .env_mach_specific.sh
echo 'export INC_NETCDF=$NETCDF_PATH/include' >> .env_mach_specific.sh
(. .env_mach_specific.sh && make) >> ${test_log} 2>&1
if [ ! -f ${mksurfdat_map} ]; then
    echo "ERROR finding/building mksurfdata_map" >&2
    echo "cat ${test_log} for more info" >&2
    exit 1
fi

echo "Running mksurfdata.pl" >> ${test_log} 2>&1
mapdir=${inputdat_root}/lnd/clm2/mappingdata/maps/ne4np4
(. .env_mach_specific.sh && ${e3sm_root}/components/elm/tools/mksurfdata_map/mksurfdata.pl -res ne4np4 -y 2010 -d -dinlc ${inputdata_root} -usr_mapdir ${mapdir}) >> ${test_log}




exit 0
