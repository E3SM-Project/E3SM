
.SUFFIXES: .F .o .cpp

OBJS =  shr_kind_mod.o \
        pcd_const.o \
        shr_const_mod.o \
        mpas_li_constants.o \
        mpas_li_mask.o \
        mpas_li_mesh.o \
        mpas_li_config.o \
        mpas_li_setup.o \
        mpas_li_time_average_coupled.o

all: $(OBJS)

mpas_li_constants.o: shr_const_mod.o

# Standalone builds compile the E3SM shared constants here; E3SM builds get them
# from the share library instead.
SHARE_UTIL = ../../../../share/util

shr_kind_mod.o: $(SHARE_UTIL)/shr_kind_mod.F90
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $< $(CPPINCLUDES) $(FCINCLUDES)

pcd_const.o: $(SHARE_UTIL)/pcd_const.F90
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $< $(CPPINCLUDES) $(FCINCLUDES)

shr_const_mod.o: $(SHARE_UTIL)/shr_const_mod.F90 shr_kind_mod.o pcd_const.o
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $< $(CPPINCLUDES) $(FCINCLUDES)

mpas_li_setup.o:

mpas_li_mask.o: mpas_li_setup.o

mpas_li_mesh.o:

mpas_li_config.o:

mpas_li_time_average_coupled.o: mpas_li_setup.o

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
