#!/bin/bash

display_help() {
    echo "Usage: $0 " >&2
    echo
    echo "  -e, --e3sm_root <e3sm_root_directory>   Specify location of E3SM"
    echo "  -h, --help                              Display this message"
    echo "  -i, --inputdata_root <data_directory>   Specify location of climate inputdata"
    echo
    echo "NOTE: requires tempestremap and ESMF tools to be in PATH environment variable"
}    

# get arguments
# Need --e3sm_root=
#      --inputdata_root=

# test mkdurfdat.pl to generate land surface data
# See step 7 in
# https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/872579110/Running+E3SM+on+New+Grids

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

if [[ ${inputdata_root} == "default" ]]; then
    echo "Error: inputdata_root not set" >&2
    display_help
    exit 1;
fi
if [[ ${e3sm_root} == "default" ]]; then
    echo "Error: e3sm_root not set" >&2
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
mkdir -p src
cd src

###########
###  1  ### Create mapping files for each land surface type if needed
###########

# a)
# Obtain or generate a target grid file in SCRIP format. For these example, we will use a ne1024pg2 grid file,
# which we will need to create (note that most np4 grid files can be found within the inputdata repository, for
# example, the ne1024np4 grid file is at
#   https://web.lcrc.anl.gov/public/e3sm/mapping/grids/ne1024np4_scrip_c20191023.nc)

generatecsmesh=$(which GenerateCSMesh)
generatevolumetricmesh=$(which GenerateVolumetricMesh)
convertmeshtoscrip=$(which ConvertMeshToSCRIP)
esmfregridweightgen=$(which ESMF_RegridWeightGen)

#test for tempestremap and ESMF tools
if [ "${esmfregridweightgen}x" == "x" ]; then
    echo "ERROR: ESMF tool ESMF_RegridWeightGen not found in PATH" >&2
    echo "cat ${test_log} for more info" >&2
    exit 1
fi

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



# These files will be created
meshfile=ne4.g
gridfile=ne4pg4.g
scripfile=ne4pg4_scrip.nc

echo "Running ${generatecsmesh}" >> ${test_log} 2>&1
(${generatecsmesh} --alt --res 4 --file ${meshfile}) >> ${test_log} 2>&1
if [ ! -f ${meshfile} ]; then
    echo "ERROR: GenerateCSMesh: no ${meshfile} file created" >&2
    echo "cat ${test_log} for more info" >&2
    exit 1
fi

echo "Running ${generatevolumetricmesh}" >> ${test_log} 2>&1
(${generatevolumetricmesh} --in ${meshfile} --out ${gridfile} --np 4 --uniform) >> ${test_log} 2>&1
if [ ! -f ${gridfile} ]; then
    echo "ERROR: GenerateVolumetricMesh: no ${gridfile} file created" >&2
    echo "cat ${test_log} for more info" >&2
    exit 1
fi

echo "Running ${convertmeshtoscrip}" >> ${test_log} 2>&1
(${convertmeshtoscrip} --in ${meshfile} --out ${scripfile}) >> ${test_log} 2>&1
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
(. .env_mach_specific.sh && set -x && ${mkmapdata} --gridfile ${scripfile} --inputdata-path ${inputdata_root} --res ne4pg4 --gridtype global --output-filetype 64bit_offset) >> ${test_log} 2>&1


# d) Create mapping file

##############
####  3  #####
##############
echo "build mksurfdata_map" >> ${test_log} 2>&1
today=$(date +%y%m%d)
cp ${e3sm_root}/components/elm/tools/mksurfdata_map/src/* .
cat <<EOF >> .env_mach_specific.sh
export USER_FC="$(awk '/MPIFC :=/ {$1=$2=""; print $0}' Macros.make)"
export USER_CPPDEFS="$(awk '/CPPDEFS :=/ {$1=$2=$3=""; print $0}' Macros.make)"
export USER_FFLAGS="$(awk '/FFLAGS :=/ {$1=$2=""; print $0}' Macros.make)"
export USER_LDFLAGS="$(awk '/SLIBS :=/ {$1=$2=$3=""; print $0}' Macros.make)"
export LIB_NETCDF=$NETCDF_PATH/lib
export INC_NETCDF=$NETCDF_PATH/include
EOF
sed -i  's|\.\./\.\./\.\.|..|' Makefile.common
(. .env_mach_specific.sh && make) >> ${test_log} 2>&1
if [ ! -f ${mksurfdat_map} ]; then
    echo "ERROR finding/building mksurfdata_map" >&2
    echo "cat ${test_log} for more info" >&2
    exit 1
fi

echo "Running mksurfdata.pl" >> ${test_log} 2>&1
echo "${e3sm_root}/components/elm/tools/mksurfdata_map/mksurfdata.pl -res usrspec -usr_gname ne4pg4 -usr_gdate ${today} -y 2010 -d -dinlc ${inputdata_root} -usr_mapdir ${PWD}"  >> ${test_log} 2>&1
(. .env_mach_specific.sh && ${e3sm_root}/components/elm/tools/mksurfdata_map/mksurfdata.pl -res usrspec -usr_gname ne4pg4 -usr_gdate ${today} -y 2010 -d -dinlc ${inputdata_root} -usr_mapdir ${PWD}) >> ${test_log} 2>&1


exit 0
