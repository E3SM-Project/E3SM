.SUFFIXES: .F .c .o

OBJS = mpas_cice_analysis_driver.o

MEMBERS = 

all: $(OBJS)

mpas_cice_analysis_driver.o: $(MEMBERS)

#mpas_cice_TEMPLATE.o:

clean:
	$(RM) *.o *.i *.mod *.f90

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES)  -I../../framework -I../../operators -I../../external/esmf_time_f90 -I../column
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES)  -I../../framework -I../../operators -I../../external/esmf_time_f90 -I../column
endif

.c.o:
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CINCLUDES) -c $<
