.SUFFIXES: .F .o

OBJS = mpas_vector_operations.o \
       mpas_matrix_operations.o \
       mpas_tensor_operations.o \
       mpas_rbf_interpolation.o \
       mpas_vector_reconstruction.o \
       mpas_spline_interpolation.o \
       mpas_tracer_advection_helpers.o \
       mpas_tracer_advection_mono.o \
       mpas_tracer_advection_std.o \
	   mpas_geometry_utils.o

DEPS := $(shell find ../core_$(CORE)/ -type f -name "*.xml" ! -name "*processed.xml")

all: operators

operators: $(OBJS) $(DEPS)
	ar -ru libops.a $(OBJS)

mpas_vector_operations.o: $(DEPS)
mpas_matrix_operations.o: $(DEPS)
mpas_tensor_operations.o: mpas_vector_operations.o mpas_matrix_operations.o $(DEPS)
mpas_rbf_interpolation.o: mpas_vector_operations.o
mpas_vector_reconstruction.o: mpas_rbf_interpolation.o
mpas_spline_interpolation:
mpas_tracer_advection_helpers.o: mpas_geometry_utils.o $(DEPS)
mpas_tracer_advection_mono.o: mpas_tracer_advection_helpers.o
mpas_tracer_advection_std.o: mpas_tracer_advection_helpers.o
mpas_geometry_utils.o:

clean:
	$(RM) *.o *.mod *.f90 libops.a
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../framework -I../external/esmf_time_f90
endif
