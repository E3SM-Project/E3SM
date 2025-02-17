#!/bin/bash

set -e

# Create a case of uELM_TES_TVA_I20TRuELMCNPRDCTCBC from 2000-2024 with ERA5 forcing

./xmlchange CONTINUE_RUN="TRUE"
./xmlchange STOP_N="24"
./xmlchange REST_N="24"

./xmlchange DATM_CLMNCEP_YR_ALIGN="2000"
./xmlchange DATM_CLMNCEP_YR_START="2000"
./xmlchange DATM_CLMNCEP_YR_END="2023"

: << 'FUTURE_PRJECTION'
./xmlchange STOP_N="76"
./xmlchange REST_N="38"

FUTURE_PRJECTION
