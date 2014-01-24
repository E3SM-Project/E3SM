MPICC=mpicc
MPIF90=mpif90

transposetest:  transposetest.c pio_spmd.o
	$(MPICC)  -std=c99 -g  $< -o $@ pio_spmd.o


pio_spmd.o: pio_spmd.c
	$(MPICC) -c -std=c99 -g -DTESTSWAPM $< -o $@

pio_spmdftest: pio_spmd_utils.F90
	$(MPIF90) -g -DTESTSWAPM -ffree-line-length-none $< -o $@
