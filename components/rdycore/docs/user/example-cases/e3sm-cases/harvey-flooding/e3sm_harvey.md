# Overview

The E3SM-RDycore model has been tested on Perlmutter and Frontier for the RDycore's 5-day Hurricane Harvey benchmark. The E3SM-RDycore simulation uses a data-land configuration with an active river model. In an E3SM-RDycore run, RDycore can run on CPUs and GPUs. The overall workflow for an E3SM-RDycore run is as follows:

- Get the code:
    - Clone the E3SM fork from [https://github.com/rdycore/e3sm](https://github.com/rdycore/e3sm).
    - Switch to the E3SM-RDycore development branch and initialize submodules of E3SM **and** RDycore.
```bash
git clone git@github.com:rdycore/e3sm
cd e3sm
git checkout rdycore/mosart-rdycore/d4ca7d0606-2024-10-02
git submodule update --init
cd externals/rdycore
git submodule update --init
```

- Create, build, and run a case
    1. Compile RDycore to generate libraries (i.e. `librdycore.a`, `librdycore_f90.a`, `libcyaml.a`, `libyaml.a`, and `libcmocka.a`)
    2. Create an E3SM case. Currently, the coupled model has been tested for a case with `--comspet RMOSGPCC --res MOS_USRDAT`.
    3. Before building the case, make the following modifications:
        - Modify the Macros file to add settings for PETSc and RDycore
        - Update the DLND streamfile (i.e `user_dlnd.streams.txt.lnd.gpcc`)
        - In `user_nl_mosart`, specify a few settings for MOSART including providing a placeholder MOSART file via `frivinp_rtm`
        - In `user_nl_dlnd`, specify a few settings for DLND including a map from the DLND mesh to the placeholder MOSART mesh
    4. Build the case
    5. Before submitting the case, do the following
        - In the rundir (`./xmlquerry RUNDIR`), copy or add symbolic links to a RDycore input YAML (as `rdycore.yaml`),
 any files specified in the RDycore's YAML file (e.g. mesh, initial condition), and map file to exchange data
 from the placeholder MOSART mesh to RDycore mesh.
        - Change the value of `run_exe` in the `env_mach_specific.xml` to include commandline options for PETSc and libCEED.
    6. Submit the case

The steps a-e have been automated via the shell script via [`e3sm_rdycore_harvey_flooding.sh`](e3sm_rdycore_harvey_flooding.sh).

```bash
cd <e3sm-rdycore/externals/rdycore/docs/user/example-cases/e3sm-cases/harvey-flooding>

./e3sm_rdycore_harvey_flooding.sh -h
Usage: ./e3sm_rdycore_harvey_flooding.sh

   -h, --help                        Display this message
   --e3sm-dir                        Path to E3SM-RDycore directory
   --mach <pm-cpu|pm-gpu|frontier>   Supported machine name
   --frontier-node-type <cpu|gpu>    To run on Frontier CPUs or GPUs
   -N, --node  <N>                   Number of nodes (default = 1)
   --project <project-id>            Project ID that will charged for the job
   --rainfall_dataset <name>  Supported dataset name (i.e. daymet|imerg|mrms|mswep|nldas)
```

## Example for Perlmutter CPU nodes

```bash
./e3sm_rdycore_harvey_flooding.sh \
--e3sm-dir /global/cfs/projectdirs/m4267/gbisht/e3sm/ \
--mach pm-cpu  \
-N 1 \
--project-id m4267 \
--rainfall_dataset mrms
```

## Example for Perlmutter GPU nodes
```bash
./e3sm_rdycore_harvey_flooding.sh \
--e3sm-dir /global/cfs/projectdirs/m4267/gbisht/e3sm/ \
--mach pm-gpu  \
-N 1 \
--project-id m4267 \
--rainfall_dataset mrms
```

## Example for Frontier using CPUs

```bash
./e3sm_rdycore_harvey_flooding.sh \
--e3sm-dir /lustre/orion/cli192/proj-shared/gb9/e3sm \
--mach frontier --frontier-node-type cpu \
-N 1 \
--project-id cli192 \
--rainfall_dataset mrms
```

## Example for Frontier using GPUs

```bash
./e3sm_rdycore_harvey_flooding.sh \
--e3sm-dir /lustre/orion/cli192/proj-shared/gb9/e3sm \
--mach frontier --frontier-node-type gpu \
-N 1 \
--project-id cli192 \
--rainfall_dataset mrms
```
