# Hurricane Harvey Flooding Simulation using Spatially-homogeneous Rainfall and Time-varying Ocean BC dataset

Both datasets, the spatially-homogeneous rainfall (`share/conditions/Houston1km.rain.<int32|int64>.bin`) and
time-varying ocean boundary condition (`share/conditions/Houston1km.bc.<int32|int64>.bin`) are PETSc `Vec`
in binary format that contain the following information:

```text
// time_1 value_1
// time_2 value_2
// time_3 value_3
```

- `time_*` has the unit of `seconds` and should start from `0.0`.
- The unit of `value_*` is
    - `m/s` for the rainfall dataset, and
    - `m` for ocean water height boundary condition dataset.

The rainfall and boundary condition dataset can be specified to the RDycore driver through
the following command line options:

1. `-rain <binary-rainfall-dataset>`
2. `-homogeneous_bc_file <binary-bc-dataset>` and

## Script

`setup_harvey_flooding_ocean_bc.sh` Will create symbolic link to the mesh
and initial condition files locally, compile RDycore (if needed), and
create a batch script for DOE supercomputers that can be submitted via `sbatch`.

```bash
 ./setup_harvey_flooding_ocean_bc.sh -h
Usage: ./setup_harvey_flooding_ocean_bc.sh

   -h, --help                        Display this message
   --rdycore-dir                     Path to RDycore directory
   --mach <pm-cpu|pm-gpu|frontier>   Supported machine name
   --frontier-node-type <cpu|gpu>    To run on Frontier CPUs or GPUs
   -N --node  <N>                    Number of nodes (default = 1)
   --project <project-id>            Project ID that will charged for the job
```

- For Perlmutter:
  - `--mach pm-cpu` Run RDycore on CPU nodes
  - `--mach pm-gpu` Run RDycore GPU nodes using CUDA
- For Frontier (`--mach frontier`):
  - `--frontier-node-type cpu`: Run RDycore on CPUs
  - `--frontier-node-type gpu`: Run RDycore on GPUs using HIP

## Example for Perlmutter CPU nodes

```bash
./setup_harvey_flooding_ocean_bc.sh \
 --mach pm-cpu -N 1 --project m4267 \
--rdycore-dir /global/cfs/projectdirs/m4267/gbisht/rdycore
```

## Example for Perlmutter GPU nodes

```bash
./setup_harvey_flooding_ocean_bc.sh \
--mach pm-gpu -N 1 --project m4267_g \
--rdycore-dir /global/cfs/projectdirs/m4267/gbisht/rdycore 
```

## Example for Frontier using CPUs

```bash
./setup_harvey_flooding_ocean_bc.sh \
--mach frontier --frontier-node-type cpu -N 2 \
--project cli192 \
--rdycore-dir /lustre/orion/cli192/proj-shared/gb9/rdycore/rdycore 
```

## Example for Frontier using GPUs

```bash
./setup_harvey_flooding_ocean_bc.sh \
--mach frontier --frontier-node-type gpu -N 1 \
--project cli192 \
--rdycore-dir /lustre/orion/cli192/proj-shared/gb9/rdycore/rdycore 
```
