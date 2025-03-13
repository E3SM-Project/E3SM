# Running Examples on Supported DOE Machines

Here we provide a few examples cases for running RDycore on DOE's supercomputers.
Each case directory contains the following:

* `index.md` that describes the case
* RDycore input YAML file for the case
* Placeholder batch scripts for the DOE supercomputers on which the case has been previously run
* A bash script that:
    * Compiles RDycore, if needed, and
    * Generates a batch script that must submitted via `sbatch`

The files for the meshes, boundary conditions, and source-sink terms are not included
in the repository and are instead available in RDycore shared project directory on the
supported DOE's supercomputers.

The following cases are supported:

1. [Idealized dam break problem](dam-break/index.md)
2. [Flooding of Houston during Hurricane Harvey](harvey-flooding/harvey-flooding.md)
