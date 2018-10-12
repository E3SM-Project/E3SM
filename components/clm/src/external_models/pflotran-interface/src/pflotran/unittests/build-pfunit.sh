#!/usr/bin/env bash
#
# script to build pfunit
#

PFUNIT_DIR=
COMPILER=
VENDOR=
DISTCLEAN=0
BUILD_STATUS=0
################################################################################
function determine-vendor() {
    _version=`${COMPILER} --version`
    if [[ ${_version} =~ "ifort" ]]; then
        VENDOR=Intel
    elif [[ ${_version} =~ "GNU" ]]; then
        VENDOR=GNU
    else
        echo "Warning: Could not determine the compiler vendor. Assuming GNU"
    fi
}

function build-pfunit() {

    echo "----------------------------------------------------------------------"
    echo "Building pFUnit :"
    make F90=${COMPILER} F90_VENDOR=${VENDOR} DEBUG=YES
    BUILD_STATUS=$?
}

function test-pfunit() {
    echo "----------------------------------------------------------------------"
    echo "Testing pFUnit :"
    ./tests/tests.x
    BUILD_STATUS=$?
}

function clean-pfunit() {
    echo "----------------------------------------------------------------------"
    echo "Cleaning pFUnit :"
    make distclean
    BUILD_STATUS=$?
}


################################################################################
#
# main program
#
################################################################################
function usage() {
     echo "
Usage: $0 [options]
    -c FORTRAN_COMPILER   path to fortran compiler
    -d PFUNIT_DIR         path to pfunit root directory
    -h                    print this help message
    -n                    run make distclean on pfunit source

Notes:

  * eventually add flags for mpi

"
}

# setup based on commandline args
while getopts "c:d:hn" FLAG
do
  case ${FLAG} in
    c) COMPILER=${OPTARG};;
    d) PFUNIT_DIR=${OPTARG};;
    h) usage; exit 0 ;;
    n) DISTCLEAN=1;;
  esac
done

# verify all required info is set
if [ -z "${PFUNIT_DIR}" ]; then
    echo "ERROR: The pFUnit root directory must be provided on the command line."
    exit 1
fi

if [ -z "${COMPILER}" ]; then
    echo "ERROR: The compiler must be provided on the command line."
    exit 1
fi

determine-vendor

echo "PFUNIT_DIR: ${PFUNIT_DIR}"
echo "FC: ${COMPILER}"
echo "VENDOR: ${VENDOR}"

pushd ${PFUNIT_DIR}
if [ "${DISTCLEAN}" -eq "1" ]; then
    clean-pfunit
else
    build-pfunit
    if [ "${BUILD_STATUS}" -eq "0" ]; then
	test-pfunit
    fi
fi
popd

exit ${BUILD_STATUS}
