############################################
# Makefile for HOMME
#############################################
F90 = mpif90
CC  = mpicc

FFLAGS = -w -fallow-argument-mismatch -ffree-line-length-none -O3 
CPPDEFS = -D_PRIM -DQSIZE_D=10 -DNP=4 -DPLEV=26  \
-DHOMME_WITHOUT_PIOLIBRARY -DPIO_INTERP -D_PRESTART \
-D_MPI -DTIMING -DHAVE_MPI \
-DSPMD  -DCPRGNU -DDARWIN -DFORTRANUNDERSCORE
 



LDFLAGS = -lblas -llapack

VPATH = ../src ../src/share ../src/test_src ../src/preqx ../src/preqx/share \
../../../cime/src/share/timing ../utils/csm_share


SRC = prim_main.F90 \
prim_driver_mod.F90 prim_driver_base.F90 \
element_ops.F90 element_state.F90 element_mod.F90 \
derivative_mod_base.F90 edge_mod_base.F90 bndry_mod_base.F90 viscosity_base.F90 \
viscosity_preqx_base.F90 prim_advection_base.F90 \
kinds.F90 parallel_mod.F90 dimensions_mod.F90 derivative_mod.F90  quadrature_mod.F90 \
physical_constants.F90 coordinate_systems_mod.F90 edgetype_mod.F90 gridgraph_mod.F90 \
control_mod.F90 params_mod.F90 hybrid_mod.F90 hybvcoord_mod.F90 physics_mod.F90 prim_si_mod.F90 \
time_mod.F90 global_norms_mod.F90 mesh_mod.F90 cube_mod.F90 spacecurve_mod.F90 repro_sum_mod.F90 \
zoltan_mod.F90 reduction_mod.F90 thread_mod.F90 edge_mod.F90 schedtype_mod.F90 \
metagraph_mod.F90 ll_mod.F90 bndry_mod.F90 prim_restart_mod.F90 restart_io_mod.F90 \
test_mod.F90 asp_tests.F90 baroclinic_inst_mod.F90 viscosity_mod.F90  \
dcmip12_wrapper.F90 dcmip2012_test1_2_3.F90 dcmip2012_test1_conv.F90 dcmip2012_test4.F90 \
mtests.F90 dcmip16_wrapper.F90 dcmip2016-baroclinic.F90 dcmip2016-supercell.F90 \
dcmip2016-tropical-cyclone.F90 dcmip2016-terminator.F90 held_suarez_mod.F90 \
dcmip2016-kessler.F90 dcmip2016-physics-z.F90 \
gllfvremap_mod.F90 gllfvremap_util_mod.F90 domain_mod.F90 namelist_mod.F90 cg_mod.F90  \
common_io_mod.F90 common_movie_mod.F90 prim_state_mod.F90 mass_matrix_mod.F90 schedule_mod.F90 \
metis_mod.F90 scalable_grid_init_mod.F90 dof_mod.F90 \
prim_advance_mod.F90 prim_advection_mod.F90 sl_advection.F90 compose_mod.F90 vertremap_base.F90 \
interpolate_mod.F90 model_init_mod.F90 vertremap_mod.F90 compose_test_mod.F90 \
interpolate_driver_mod.F90 interp_movie_mod.F90 \
shr_kind_mod.F90 shr_sys_mod.F90 shr_mpi_mod.F90 \
perf_mod.F90 perf_utils.F90

CC_SRC = gptl_papi.c gptl.c f_wrappers.c GPTLutil.c GPTLprint_memusage.c \
  GPTLget_memusage.c 



OBJ = $(SRC:.F90=.o)
CC_OBJ = $(CC_SRC:.c=.o)

preqx: $(OBJ) $(CC_OBJ)
	$(F90) $(FFLAGS) $(OBJ) $(CC_OBJ) $(LDFLAGS) -o preqx

clean:
	rm *.o *.mod

%.o : %.F90
	$(F90) $(FFLAGS) $(CPPDEFS) -c $< -o $@

%.o : %.c
	$(CC) $(CFLAGS) $(CPPDEFS) -I ../../../cime/src/share/timing  -c $< -o $@


SrcFiles:
	echo "$(SRC)" > SrcFiles

Depends: SrcFiles
	./mkDepends FilePath SrcFiles > Depends

dcmip16_wrapper.o: dcmip2016-baroclinic.o dcmip2016-supercell.o dcmip2016-tropical-cyclone.o \
dcmip2016-terminator.o

perf_mod.o: perf_utils.o

shr_sys_mod.o: shr_sys_mod.F90 shr_mpi_mod.o


include ./Depends

