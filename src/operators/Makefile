.SUFFIXES: .F .o

OBJS = mpas_rbf_interpolation.o mpas_vector_reconstruction.o mpas_spline_interpolation.o

all: operators

operators: $(OBJS)
	ar -ru libops.a $(OBJS)

mpas_vector_reconstruction.o: mpas_rbf_interpolation.o
mpas_rbf_interpolation.o:
mpas_spline_interpolation:

clean:
	$(RM) *.o *.mod *.f90 libops.a

.F.o:
	$(RM) $@ $*.mod
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework
