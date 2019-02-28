$(aPo)/scrip_kindsmod.o: SCRIP/scrip_kindsmod.f90
	@$(aPb)/ad3 scrip_kindsmod

$(aPo)/scrip_constants.o: SCRIP/scrip_constants.f \
	$(aPo)/scrip_kindsmod.o
	@$(aPb)/ad3 scrip_constants

$(aPo)/scrip_iounitsmod.o: SCRIP/scrip_iounitsmod.f90 \
	$(aPo)/scrip_kindsmod.o
	@$(aPb)/ad3 scrip_iounitsmod

$(aPo)/scrip_errormod.o: SCRIP/scrip_errormod.f90 \
	$(aPo)/scrip_kindsmod.o \
	$(aPo)/scrip_iounitsmod.o
	@$(aPb)/ad3 scrip_errormod

$(aPo)/scrip_netcdfmod.o: SCRIP/scrip_netcdfmod.f90 \
	$(aPo)/scrip_kindsmod.o \
	$(aPo)/scrip_errormod.o
	@$(aPb)/ad3 scrip_netcdfmod

$(aPo)/scrip_grids.o: SCRIP/scrip_grids.f \
	$(aPo)/scrip_kindsmod.o \
	$(aPo)/scrip_constants.o \
	$(aPo)/scrip_iounitsmod.o \
	$(aPo)/scrip_netcdfmod.o \
   $(aPo)/scrip_errormod.o
	@$(aPb)/ad3 scrip_grids

$(aPo)/scrip_remap_vars.o: SCRIP/scrip_remap_vars.f \
	$(aPo)/scrip_kindsmod.o \
	$(aPo)/scrip_constants.o \
	$(aPo)/scrip_grids.o \
	$(aPo)/scrip_errormod.o \
	$(aPo)/scrip_netcdfmod.o \
	$(aPo)/scrip_iounitsmod.o 
	@$(aPb)/ad3 scrip_remap_vars

$(aPo)/scrip_remap_conservative.o: SCRIP/scrip_remap_conservative.f \
	$(aPo)/scrip_kindsmod.o \
	$(aPo)/scrip_constants.o \
	$(aPo)/scrip_timers.o \
	$(aPo)/scrip_remap_vars.o \
	$(aPo)/scrip_grids.o \
	$(aPo)/scrip_errormod.o \
	$(aPo)/scrip_netcdfmod.o \
	$(aPo)/scrip_iounitsmod.o 
	@$(aPb)/ad3 scrip_remap_conservative

$(aPo)/scrip_timers.o: SCRIP/scrip_timers.f \
	$(aPo)/scrip_kindsmod.o
	@$(aPb)/ad3 scrip_timers

$(aPo)/scrip_remap_write.o: SCRIP/scrip_remap_write.f \
	$(aPo)/scrip_kindsmod.o \
	$(aPo)/scrip_constants.o \
	$(aPo)/scrip_netcdfmod.o \
	$(aPo)/scrip_remap_vars.o \
	$(aPo)/scrip_grids.o \
	$(aPo)/scrip_errormod.o \
	$(aPo)/scrip_iounitsmod.o
	@$(aPb)/ad3 scrip_remap_write

$(aPo)/scrip_remap_read.o: SCRIP/scrip_remap_read.f \
	$(aPo)/scrip_kindsmod.o \
	$(aPo)/scrip_constants.o \
	$(aPo)/scrip_netcdfmod.o \
	$(aPo)/scrip_remap_vars.o \
	$(aPo)/scrip_grids.o \
	$(aPo)/scrip_errormod.o \
	$(aPo)/scrip_iounitsmod.o
	@$(aPb)/ad3 scrip_remap_read

$(aPo)/scrip_interface.o: SCRIP/scrip_interface.ftn \
	$(aPo)/scrip_kindsmod.o \
	$(aPo)/scrip_constants.o \
	$(aPo)/scrip_timers.o \
	$(aPo)/scrip_remap_vars.o \
	$(aPo)/scrip_grids.o \
	$(aPo)/scrip_remap_conservative.o \
	$(aPo)/scrip_netcdfmod.o \
	$(aPo)/scrip_remap_write.o \
	$(aPo)/scrip_remap_read.o \
	$(aPo)/scrip_iounitsmod.o \
	$(aPo)/scrip_errormod.o 
	@$(aPb)/ad3 scrip_interface

