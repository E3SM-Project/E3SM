The RDycore team has preinstalled PETSc on few of the DOE machinces
and this directory contains information about those PETSc installations.
Below is the information about files contained in this directory.

- `modules.<machine>.<compiler>`: Modules that should be loaded for a given machine and 
  compiler combination prior to installing PETSc or using preinstalled PETSc.
- `install.<machine>.<compiler>.<optimiztion-or-debug>.<32bit-or-64bit>.sh`: Script to install
  PETSc for a given machine, compiler, optimization level, and 32bit or 64bit integer support.
- `set_petsc_settings.sh`: A script that loads modules and sets PETSc environment variables
  (i.e., `PETSC_DIR` and `PETSC_ARCH`) for the preinstalled PETSc on the machine.

NOTE: The Perlmutter supercomputer has two types of nodes: (1) CPU-only, and (2) CPU-GPU.