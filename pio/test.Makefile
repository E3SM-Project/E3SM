MPICC=mpicc
MPIF90=mpif90

pio_spmdctest: pio_spmd.c
	$(MPICC) -std=c99 -g -DTESTSWAPM $< -o $@

pio_spmdftest: pio_spmd_utils.F90
	$(MPIF90) -g -DTESTSWAPM -ffree-line-length-none $< -o $@
