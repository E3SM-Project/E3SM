# kmELM Reference Simulation README

This document provides instructions for configuring and running the kmELM reference simulation on the several machines (Frontier, Perlmutter, and Baseline). This README outlines the steps necessary to set up a high-resolution reference case using the TESSFA branch of the kmELM code.

## Prerequisites

- Access to the CCSI/ORNL Linux environment.
- SSH client to log into the machine.
- Git installed to clone the repository.

## Setup Instructions

### 1. Log into your machine
Use SSH to log into the Baseline system:
```bash
ssh wangd@baseline.ccs.ornl.gov
```

### 2. Create a new directory
Create a directory for the kmELM source code:
```bash
mkdir -p /gpfs/wolf2/cades/cli185/proj-shared/wangd/kmELM
cd /gpfs/wolf2/cades/cli185/proj-shared/wangd/kmELM
```

### 3. Clone the kmELM repository
Checkout the kmELM code from GitHub:
```bash
git clone git@github.com:daliwang/kmELM.git
cd kmELM
git submodule update --init --recursive
```
> **Note:** Type `yes` when prompted to initialize submodules.

### 4. Checkout the required branch
Switch to the `TESSFA` branch:
```bash
git fetch origin
git checkout -b my-dev origin/TESSFA
```

### 5. Configure for your HPC system (Optional)
If you are running kmELM on a different HPC system, you will need to modify the configuration files in the `kmELM/cime_config/machines` directory. You may need to edit:
- `config_machines.xml`
- `config_batch.xml`
- `config_pio`
- `config_workflow`
- `cmake_macros`

Refer to the sections for CCSI-Baseline, SummitPlus, or Frontier for guidance on how to adjust these settings for your environment.

### 6. Navigate to case generation directory
Go to the appropriate case generation directory for your machine (e.g., ORNL_baseline or Frontier):
```bash
cd kmELM/case_gene/<machines>
```

### 7. Review case creation script
Examine the case generation script, `uELM_caseGEN_TVA.sh`, to understand its functionality.

### 8. Modify case generation script
You will need to modify the `uELM_caseGEN_TVA.sh` script with the appropriate paths:
```bash
E3SM_SRCROOT="/gpfs/wolf2/cades/cli185/proj-shared/wangd/kmELM"
CASEDIR="/gpfs/wolf2/cades/cli185/proj-shared/wangd/e3sm_cases/uELM_${EXPID}_I1850uELMCNPRDCTCBC"
```

### 9. Create the uELM case
Run the case generation script:
```bash
sh uELM_caseGEN_TVA.sh
```
This script will create the necessary case directories.

### 10. Build and run the case
Once the case has been created, you can build and run it:
```bash
./case.build
./case.submit
```

Following these steps will enable you to set up and run the kmELM reference simulation on the selected HPC systems. Adjust paths and configurations as needed based on your specific environment and requirements. 

For any further assistance, refer to the E3SM documentation or contact the kmELM development team.
