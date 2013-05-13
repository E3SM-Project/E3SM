.SUFFIXES: .F .o

CVMIX_REPO_ADDRESS=http://cvmix.googlecode.com/svn/trunk/src/shared

OBJS = mpas_ocn_mpas_core.o \
       mpas_ocn_advection.o \
       mpas_ocn_thick_hadv.o \
       mpas_ocn_thick_vadv.o \
       mpas_ocn_gm.o \
       mpas_ocn_vel_coriolis.o \
       mpas_ocn_vel_vadv.o \
       mpas_ocn_vel_hmix.o \
       mpas_ocn_vel_hmix_del2.o \
       mpas_ocn_vel_hmix_leith.o \
       mpas_ocn_vel_hmix_del4.o \
       mpas_ocn_vel_forcing.o \
       mpas_ocn_vel_forcing_windstress.o \
       mpas_ocn_vel_forcing_rayleigh.o \
       mpas_ocn_vel_pressure_grad.o \
       mpas_ocn_vmix.o \
       mpas_ocn_vmix_coefs_const.o \
       mpas_ocn_vmix_coefs_rich.o \
       mpas_ocn_vmix_coefs_tanh.o \
       mpas_ocn_restoring.o \
       mpas_ocn_tendency.o \
       mpas_ocn_diagnostics.o \
       mpas_ocn_tracer_hmix.o \
       mpas_ocn_tracer_hmix_del2.o \
       mpas_ocn_tracer_hmix_del4.o \
       mpas_ocn_tracer_advection.o \
       mpas_ocn_tracer_advection_std.o \
       mpas_ocn_tracer_advection_std_hadv.o \
       mpas_ocn_tracer_advection_std_vadv.o \
       mpas_ocn_tracer_advection_std_vadv2.o \
       mpas_ocn_tracer_advection_std_vadv3.o \
       mpas_ocn_tracer_advection_std_vadv4.o \
       mpas_ocn_tracer_advection_mono.o \
       mpas_ocn_tracer_advection_helpers.o \
       mpas_ocn_time_integration.o \
       mpas_ocn_time_integration_rk4.o \
       mpas_ocn_time_integration_split.o \
       mpas_ocn_equation_of_state.o \
       mpas_ocn_equation_of_state_jm.o \
       mpas_ocn_equation_of_state_linear.o \
       mpas_ocn_diagnostics.o \
       mpas_ocn_global_diagnostics.o \
       mpas_ocn_time_average.o \
       mpas_ocn_monthly_forcing.o

all: libcvmix core_hyd

libcvmix:
	if [ ! -d cvmix ]; then \
	(svn checkout $(CVMIX_REPO_ADDRESS) cvmix) \
	fi
	if [ -d cvmix ]; then \
		(cd cvmix; svn update; make all FC="$(FC)" FFLAGS="$(FFLAGS)" FINCLUDES="$(FINCLUDES)") \
	fi

core_hyd: $(OBJS)
	ar -ru libdycore.a $(OBJS) cvmix/*.o

mpas_ocn_advection.o:

mpas_ocn_time_integration.o: mpas_ocn_time_integration_rk4.o mpas_ocn_time_integration_split.o

mpas_ocn_time_integration_rk4.o: mpas_ocn_tendency.o mpas_ocn_diagnostics.o

mpas_ocn_time_integration_split.o: mpas_ocn_tendency.o mpas_ocn_diagnostics.o

mpas_ocn_tendency.o: mpas_ocn_time_average.o

mpas_ocn_diagnostics.o: mpas_ocn_time_average.o

mpas_ocn_global_diagnostics.o: 

mpas_ocn_time_average.o:

mpas_ocn_thick_hadv.o:

mpas_ocn_thick_vadv.o:

mpas_ocn_gm.o: 

mpas_ocn_vel_pressure_grad.o:

mpas_ocn_vel_vadv.o:

mpas_ocn_vel_hmix.o: mpas_ocn_vel_hmix_del2.o mpas_ocn_vel_hmix_leith.o mpas_ocn_vel_hmix_del4.o

mpas_ocn_vel_hmix_del2.o:

mpas_ocn_vel_hmix_leith.o:

mpas_ocn_vel_hmix_del4.o:

mpas_ocn_vel_forcing.o: mpas_ocn_vel_forcing_windstress.o mpas_ocn_vel_forcing_rayleigh.o

mpas_ocn_vel_forcing_windstress.o:

mpas_ocn_vel_forcing_rayleigh.o:

mpas_ocn_vel_coriolis.o:

mpas_ocn_tracer_hmix.o: mpas_ocn_tracer_hmix_del2.o mpas_ocn_tracer_hmix_del4.o

mpas_ocn_tracer_hmix_del2.o:

mpas_ocn_tracer_hmix_del4.o:

mpas_ocn_tracer_advection.o: mpas_ocn_tracer_advection_std.o mpas_ocn_tracer_advection_mono.o

mpas_ocn_tracer_advection_std.o: mpas_ocn_tracer_advection_std_hadv.o mpas_ocn_tracer_advection_std_vadv.o

mpas_ocn_tracer_advection_std_hadv.o: mpas_ocn_tracer_advection_helpers.o

mpas_ocn_tracer_advection_std_vadv.o: mpas_ocn_tracer_advection_std_vadv2.o mpas_ocn_tracer_advection_std_vadv3.o mpas_ocn_tracer_advection_std_vadv4.o

mpas_ocn_tracer_advection_std_vadv2.o: mpas_ocn_tracer_advection_helpers.o

mpas_ocn_tracer_advection_std_vadv3.o: mpas_ocn_tracer_advection_helpers.o

mpas_ocn_tracer_advection_std_vadv4.o: mpas_ocn_tracer_advection_helpers.o

mpas_ocn_tracer_advection_mono.o: mpas_ocn_tracer_advection_helpers.o

mpas_ocn_tracer_advection_helpers.o:

mpas_ocn_restoring.o:

mpas_ocn_vmix.o: mpas_ocn_vmix_coefs_const.o mpas_ocn_vmix_coefs_rich.o mpas_ocn_vmix_coefs_tanh.o

mpas_ocn_vmix_coefs_const.o:

mpas_ocn_vmix_coefs_rich.o: mpas_ocn_equation_of_state.o

mpas_ocn_vmix_coefs_tanh.o:

mpas_ocn_equation_of_state.o: mpas_ocn_equation_of_state_jm.o mpas_ocn_equation_of_state_linear.o

mpas_ocn_equation_of_state_jm.o:

mpas_ocn_equation_of_state_linear.o:

mpas_ocn_monthly_forcing.o:

mpas_ocn_mpas_core.o: mpas_ocn_advection.o \
                      mpas_ocn_thick_hadv.o \
                      mpas_ocn_gm.o \
                      mpas_ocn_thick_vadv.o \
                      mpas_ocn_vel_coriolis.o \
                      mpas_ocn_vel_vadv.o \
                      mpas_ocn_vel_hmix.o \
                      mpas_ocn_vel_hmix_del2.o \
                      mpas_ocn_vel_hmix_leith.o \
                      mpas_ocn_vel_hmix_del4.o \
                      mpas_ocn_vel_forcing.o \
                      mpas_ocn_vel_forcing_windstress.o \
                      mpas_ocn_vel_pressure_grad.o \
                      mpas_ocn_tracer_hmix.o \
                      mpas_ocn_tracer_hmix_del2.o \
                      mpas_ocn_tracer_hmix_del4.o \
                      mpas_ocn_vmix.o \
                      mpas_ocn_vmix_coefs_const.o \
                      mpas_ocn_vmix_coefs_rich.o \
                      mpas_ocn_vmix_coefs_tanh.o \
                      mpas_ocn_restoring.o \
                      mpas_ocn_tracer_advection.o \
                      mpas_ocn_tracer_advection_std.o \
                      mpas_ocn_tracer_advection_std_hadv.o \
                      mpas_ocn_tracer_advection_std_vadv.o \
                      mpas_ocn_tracer_advection_std_vadv2.o \
                      mpas_ocn_tracer_advection_std_vadv3.o \
                      mpas_ocn_tracer_advection_std_vadv4.o \
                      mpas_ocn_tracer_advection_mono.o \
                      mpas_ocn_tracer_advection_helpers.o \
                      mpas_ocn_tendency.o \
                      mpas_ocn_diagnostics.o \
                      mpas_ocn_time_integration.o \
                      mpas_ocn_time_integration_rk4.o \
                      mpas_ocn_time_integration_split.o \
                      mpas_ocn_equation_of_state.o \
                      mpas_ocn_equation_of_state_jm.o \
                      mpas_ocn_equation_of_state_linear.o \
                      mpas_ocn_global_diagnostics.o \
                      mpas_ocn_time_average.o \
                      mpas_ocn_monthly_forcing.o

clean:
	if [ -d cvmix ]; then \
		(cd cvmix; make clean) \
	fi
	$(RM) *.o *.mod *.f90 libdycore.a

.F.o:
	$(RM) $@ $*.mod
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../operators -I../external/esmf_time_f90
