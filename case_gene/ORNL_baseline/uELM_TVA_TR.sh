#!/bin/bash

set -e

# Create a test case uELM_TES_TVA_I1850uELMCNPRDCTCBC

#E3SM_DIN="/gpfs/wolf2/cades/cli185/proj-shared/pt-e3sm-inputdata"
E3SM_DIN="//gpfs/wolf2/cades/cli185/world-shared/e3sm/inputdata"
DATA_ROOT="/gpfs/wolf2/cades/cli185/proj-shared/wangd/kiloCraft/TES_cases_data/Daymet_ERA5_TESSFA2/"
E3SM_SRCROOT=$(git rev-parse --show-toplevel)
echo "E3SM_SRCROOT: $E3SM_SRCROOT"
echo "E3SM_DIN: $E3SM_DIN"

EXPID="TVA"
CASEDIR="$E3SM_SRCROOT/e3sm_cases/uELM_${EXPID}_I20TRuELMCNPRDCTCBC"
CASE_DATA="${DATA_ROOT}/${EXPID}"
DOMAIN_FILE="${EXPID}_domain.lnd.TES_SE.4km.1d.c241218.nc"
SURFDATA_FILE="${EXPID}_surfdata.TES_SE.4km.1d.NLCD.c241219.nc"

DOMAIN_SOFT_LINK="${CASE_DATA}/atm_forcing.datm7.km.1d/domain.lnd.Daymet.km.1d.nc"

# add a soft link to domain (lnd/atm)
# Check if the soft link exists
if [ ! -L "$DOMAIN_SOFT_LINK" ]; then
    # The soft link does not exist, create it
    ln -s "${CASE_DATA}/domain_surfdata/${DOMAIN_FILE}" "$SOFT_LINK"
    echo "Soft link created: $SOFT_LINK -> ${CASE_DATA}/domain_surfdata/${DOMAIN_FILE}"
fi

\rm -rf "${CASEDIR}"

#${E3SM_SRCROOT}/cime/scripts/create_newcase --case "${CASEDIR}" --mach summitPlus --compiler pgi --mpilib spectrum-mpi --compset I1850uELMCNPRDCTCBC --res ELM_USRDAT --pecount "${PECOUNT}" --handle-preexisting-dirs r --srcroot "${E3SM_SRCROOT}"

${E3SM_SRCROOT}/cime/scripts/create_newcase --case "${CASEDIR}" --mach cades-baseline --compiler gnu --mpilib openmpi --compset I20TRuELMTESCNPRDCTCBC --res ELM_USRDAT  --handle-preexisting-dirs r --srcroot "${E3SM_SRCROOT}"

cd "${CASEDIR}"

./xmlchange PIO_TYPENAME="pnetcdf"

./xmlchange PIO_NETCDF_FORMAT="64bit_data"

./xmlchange DIN_LOC_ROOT="${E3SM_DIN}"

./xmlchange DIN_LOC_ROOT_CLMFORC="${CASE_DATA}"

./xmlchange CIME_OUTPUT_ROOT="${E3SM_SRCROOT}/e3sm_runs/"

./xmlchange DATM_MODE="uELM_TES"

./xmlchange NTASKS="1"
./xmlchange NTASKS_LND="1280"
./xmlchange NTASKS_ATM="100"
./xmlchange NTASKS_CPL="100"

./xmlchange NTASKS_PER_INST="1"

./xmlchange MAX_MPITASKS_PER_NODE="128"

./xmlchange ATM_DOMAIN_PATH="${CASE_DATA}/domain_surfdata/"

./xmlchange ATM_DOMAIN_FILE="${DOMAIN_FILE}"

./xmlchange LND_DOMAIN_PATH="${CASE_DATA}/domain_surfdata/"

./xmlchange LND_DOMAIN_FILE="${DOMAIN_FILE}"

./xmlchange JOB_WALLCLOCK_TIME="24:00:00"

./xmlchange USER_REQUESTED_WALLTIME="24:00:00"

./xmlchange ATM_NCPL="24"

./xmlchange DATM_CLMNCEP_YR_START="1980"
./xmlchange DATM_CLMNCEP_YR_END="1999"
./xmlchange DATM_CLMNCEP_YR_ALIGN="1990"

./xmlchange STOP_N="150"
./xmlchange REST_N="50"
./xmlchange STOP_OPTION="nyears"

./xmlchange CONTINUE_RUN="FALSE"
./xmlchange RUN_TYPE="startup"
./xmlchange CCSM_BGC="CO2A"

echo "
finidat = '/gpfs/wolf2/cades/cli185/proj-shared/wangd/kmELM/e3sm_runs/uELM_TVA_I1850uELMCNPRDCTCBC/run_spinup/uELM_TVA_I1850uELMCNPRDCTCBC.elm.r.0801-01-01-00000.nc'
fsurdat = '${CASE_DATA}/domain_surfdata/${SURFDATA_FILE}'
      hist_mfilt=1
     flanduse_timeseries = ''
      " >> user_nl_elm

./case.setup --reset

./case.setup

./case.build --clean-all

./case.build

#./case.submit

