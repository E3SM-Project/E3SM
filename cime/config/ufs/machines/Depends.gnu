geopk.o:geopk.F90
	$(FC) -c $(INCLDIR) $(INCS) $(FFLAGS) $(FREEFLAGS) -fcray-pointer $<
