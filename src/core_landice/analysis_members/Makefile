.SUFFIXES: .F .c .o

OBJS = mpas_li_analysis_driver.o

MEMBERS = mpas_li_global_stats.o \
          mpas_li_regional_stats.o 

all: $(OBJS)

mpas_li_analysis_driver.o: $(MEMBERS)

clean:
	$(RM) *.o *.i *.mod *.f90

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) 
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) 
endif

.c.o:
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CINCLUDES) -c $<
