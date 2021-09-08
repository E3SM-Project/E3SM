
.SUFFIXES: .F .o .cpp

OBJS =  mpas_li_constants.o \
        mpas_li_mask.o \
        mpas_li_setup.o

all: $(OBJS)

mpas_li_constants.o:

mpas_li_setup.o:

mpas_li_mask.o: mpas_li_setup.o



clean:
	$(RM) *.o *.mod *.f90
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) 
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) 
endif

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $*.cpp $(CXINCLUDES) $(CPPINCLUDES) -lmpi_cxx -lstdc++ $(CPPFLAGS)
