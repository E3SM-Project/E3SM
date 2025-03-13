# Hurricane Harvey Flooding using MRMS Dataset and Critical Outflow BC

The 1 km Multi-Radar Multi-Sensor System (MRMS) rainfall dataset ([Zhang et al., 2016](https://journals.ametsoc.org/view/journals/bams/97/4/bams-d-14-00174.1.xml)), is used
in this simulation. Each hourly MRMS dataset file is a PETSc Vec saved in
the binary format, named as `YYYY-MM-DD:HH-SS.<int32|int64>.bin`. These
binary files contains the following information:

- `ncols`    : number of columns in the rainfall dataset
- `nrows`    : number of rowns in the rainfall dataset
- `xlc`      : x coordinate of the lower left corner [m]
- `ylc`      : y coordinate of the lower left corner [m]
- `cellsize` : size of grid cells in the rainfall dataset [m]
- `data`     : rainfall rate for ncols * nrows cells [mm/hr]

The start date and the location of MRMS dataset is specified to the RDycore
driver through the following two command line options:

1. `-raster_rain_start_date YYYY,MM,DD,HH,SS`, and
2. `-raster_rain_dir <path/to/the/mrms/dataset>`

The MRMS dataset and mesh are available on Perlmutter and Frontier under
the RDycore's project directoy. Apart from the MRMS dataset, following the four
spatially-distributed rainfall datasets are also available:

1. Daymet,
2. North American Land Data Assimilation System (NLDAS),
3. Integrated multi-satellite retrievals for global precipitation measurement (IMERG), and
4. Multi-Source Weighted-Ensemble Precipitation (MSWEP)

A critical outflow boundary condition is applied to the 13 edges within the
`Turning_30m_with_z.updated.with_sidesets.exo` mesh file that are identified
with the labels `elem_ss1` and `side_ss1`.

## Script

`setup_harvey_flooding_critical_outflow_bc.sh` Will create symbolic link to the mesh
and initial condition files locally, compile RDycore (if needed), and
create a batch script for DOE supercomputers that can be submitted via `sbatch`.


```bash
 ./setup_harvey_flooding_critical_outflow_bc.sh -h
Usage: ./setup_harvey_flooding_critical_outflow_bc.sh

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
./setup_harvey_flooding_critical_outflow_bc.sh \
--mach pm-cpu -N 1 --project m4267 \
--rdycore-dir /global/cfs/projectdirs/m4267/gbisht/rdycore
```

## Example for Perlmutter GPU nodes

```bash
./setup_harvey_flooding_critical_outflow_bc.sh \
--mach pm-gpu -N 1 --project m4267_g \
--rdycore-dir /global/cfs/projectdirs/m4267/gbisht/rdycore
```

## Example for Frontier using CPUs

```bash
./setup_harvey_flooding_critical_outflow_bc.sh \
--mach frontier --frontier-node-type cpu -N 2 \
--project cli192 \
--rdycore-dir /lustre/orion/cli192/proj-shared/gb9/rdycore/rdycore 
```

## Example for Frontier using GPUs

```bash
./setup_harvey_flooding_critical_outflow_bc.sh \
--mach frontier --frontier-node-type gpu -N 1 \
--project cli192 \
--rdycore-dir /lustre/orion/cli192/proj-shared/gb9/rdycore/rdycore 
```