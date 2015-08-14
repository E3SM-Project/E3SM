
# ===================================
# Check if building with LifeV, Albany, and/or PHG external libraries

BUILD_INTERFACE=false  # This will become true if any of the external libraries are being used.

# LifeV can solve L1L2 or FO
ifeq "$(LIFEV)" "true"
    EXTERNAL_DYCORE_FLAG += -DUSE_EXTERNAL_L1L2
    EXTERNAL_DYCORE_FLAG += -DUSE_EXTERNAL_FIRSTORDER
    BUILD_INTERFACE = true
endif # LIFEV IF

# Albany can only solve FO at present
ifeq "$(ALBANY)" "true"
    EXTERNAL_DYCORE_FLAG += -DUSE_EXTERNAL_FIRSTORDER
    BUILD_INTERFACE = true
endif # ALBANY IF

# Currently LifeV AND Albany is not allowed
ifeq "$(LIFEV)" "true"
ifeq "$(ALBANY)" "true"
    $(error Compiling with both LifeV and Albany is not allowed at this time.)
endif
endif

# PHG currently requires LifeV
ifeq "$(PHG)" "true"
ifneq "$(LIFEV)" "true"
    $(error Compiling with PHG requires LifeV at this time.)
endif
endif
# PHG can only Stokes at present
ifeq "$(PHG)" "true"
    EXTERNAL_DYCORE_FLAG += -DUSE_EXTERNAL_STOKES
    BUILD_INTERFACE = true
endif # PHG IF

override CPPFLAGS += $(EXTERNAL_DYCORE_FLAG)
# ===================================

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
