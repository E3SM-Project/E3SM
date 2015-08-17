
.SUFFIXES: .F .o .cpp

OBJS =  mpas_li_core.o \
    mpas_li_core_interface.o \
    mpas_li_time_integration.o \
    mpas_li_time_integration_fe.o \
    mpas_li_diagnostic_vars.o \
    mpas_li_tendency.o \
    mpas_li_setup.o \
    mpas_li_statistics.o \
    mpas_li_velocity.o \
    mpas_li_sia.o \
    mpas_li_mask.o \
    mpas_li_velocity_external.o

ifeq "$(BUILD_INTERFACE)" "true"
    OBJS += Interface_velocity_solver.o
endif

all: $(OBJS)

mpas_li_core_interface.o: mpas_li_core.o

mpas_li_core.o: mpas_li_time_integration.o \
                     mpas_li_setup.o \
                     mpas_li_velocity.o \
                     mpas_li_diagnostic_vars.o \
                     mpas_li_statistics.o \
                     mpas_li_mask.o

mpas_li_setup.o:

mpas_li_time_integration.o: mpas_li_time_integration_fe.o

mpas_li_time_integration_fe.o: mpas_li_velocity.o \
                               mpas_li_tendency.o \
                               mpas_li_diagnostic_vars.o \
                               mpas_li_setup.o

mpas_li_tendency.o: mpas_li_setup.o

mpas_li_diagnostic_vars.o: mpas_li_mask.o \
                           mpas_li_velocity.o \
                           mpas_li_constants.o

mpas_li_velocity.o: mpas_li_sia.o \
                    mpas_li_setup.o \
                    mpas_li_velocity_external.o

mpas_li_sia.o: mpas_li_mask.o \
               mpas_li_setup.o

mpas_li_statistics.o: mpas_li_mask.o \
                      mpas_li_setup.o \
                      mpas_li_constants.o

mpas_li_mask.o: mpas_li_setup.o

mpas_li_constants.o:

mpas_li_velocity_external.o:

Interface_velocity_solver.o:

clean:
	$(RM) *.o *.mod *.f90
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F.o:
	$(RM) $@ $*.mod
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../../framework -I../../operators -I../../external/esmf_time_f90 -I../shared

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $*.cpp $(CXINCLUDES) $(CPPINCLUDES) -lmpi_cxx -lstdc++ $(CPPFLAGS)
