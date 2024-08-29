
.SUFFIXES: .F .o .cpp

OBJS =  mpas_li_core.o \
    mpas_li_core_interface.o \
    mpas_li_time_integration.o \
    mpas_li_time_integration_fe_rk.o \
    mpas_li_diagnostic_vars.o \
    mpas_li_advection.o \
    mpas_li_advection_fct_shared.o \
    mpas_li_advection_fct.o \
    mpas_li_calving.o \
    mpas_li_statistics.o \
    mpas_li_velocity.o \
    mpas_li_thermal.o \
    mpas_li_iceshelf_melt.o \
    mpas_li_sia.o \
    mpas_li_velocity_simple.o \
    mpas_li_velocity_external.o \
    mpas_li_subglacial_hydro.o \
    mpas_li_bedtopo.o \
    mpas_li_ocean_extrap.o

ifneq (, $(findstring MPAS_LI_BUILD_INTERFACE,$(CPPFLAGS)))
    OBJS += Interface_velocity_solver.o
endif

all: $(OBJS)

mpas_li_core_interface.o: mpas_li_core.o

mpas_li_core.o: mpas_li_time_integration.o \
                     mpas_li_velocity.o \
                     mpas_li_thermal.o \
                     mpas_li_diagnostic_vars.o \
                     mpas_li_statistics.o \
                     mpas_li_calving.o

mpas_li_time_integration.o: mpas_li_time_integration_fe_rk.o

mpas_li_time_integration_fe_rk.o: mpas_li_advection.o \
	                       mpas_li_calving.o \
                               mpas_li_thermal.o \
			       mpas_li_iceshelf_melt.o \
                               mpas_li_diagnostic_vars.o \
			       mpas_li_velocity.o \
			       mpas_li_subglacial_hydro.o \
			       mpas_li_bedtopo.o \
                               mpas_li_ocean_extrap.o

mpas_li_advection.o: mpas_li_thermal.o \
                     mpas_li_advection_fct.o \
                     mpas_li_advection_fct_shared.o \
	             mpas_li_diagnostic_vars.o

mpas_li_advection_fct.o:  mpas_li_advection_fct_shared.o

mpas_li_calving.o: mpas_li_thermal.o \
		   mpas_li_advection.o

mpas_li_thermal.o:

mpas_li_iceshelf_melt.o: mpas_li_calving.o

mpas_li_diagnostic_vars.o: mpas_li_thermal.o

mpas_li_velocity.o: mpas_li_sia.o \
                    mpas_li_velocity_simple.o \
                    mpas_li_velocity_external.o \
		    mpas_li_advection.o \
		    mpas_li_thermal.o

mpas_li_sia.o: mpas_li_diagnostic_vars.o

mpas_li_subglacial_hydro.o: mpas_li_diagnostic_vars.o

mpas_li_velocity_simple.o:

mpas_li_statistics.o: mpas_li_diagnostic_vars.o

mpas_li_velocity_external.o: Interface_velocity_solver.o

mpas_li_bedtopo.o: mpas_li_advection.o

mpas_li_ocean_extrap.o: mpas_li_calving.o

Interface_velocity_solver.o:

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
