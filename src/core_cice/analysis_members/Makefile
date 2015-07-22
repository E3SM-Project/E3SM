.SUFFIXES: .F .c .o

OBJS = mpas_cice_analysis_driver.o

MEMBERS = mpas_cice_high_frequency_output.o

all: $(OBJS)

mpas_cice_analysis_driver.o: $(MEMBERS)

clean:
	$(RM) *.o *.i *.mod *.f90

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 -I../../external/esmf_time_f90 -I../../framework -I../../operators -I../forward_model -I../column $(FCINCLUDES) 
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F -I../../external/esmf_time_f90 -I../../framework -I../../operators -I../forward_model -I../column $(CPPINCLUDES) $(FCINCLUDES) 
endif

.c.o:
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CINCLUDES) -I../../external/esmf_time_f90 -I../../framework -I../../operators -I../forward_model -I../column -c $<
