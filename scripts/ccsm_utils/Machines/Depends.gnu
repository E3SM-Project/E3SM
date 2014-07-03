shr_infnan_mod.o: shr_infnan_mod.F90
	$(FC) -c $(INCLDIR) $(INCS) $(FFLAGS) $(FREEFLAGS) -fno-range-check $<
shr_test_infnan_mod.o: shr_test_infnan_mod.F90
	$(FC) -c $(INCLDIR) $(INCS) $(FFLAGS) $(FREEFLAGS) -fno-range-check $<
seq_drydep_mod.o: seq_drydep_mod.F90
	$(FC) -c $(INCLDIR) $(INCS) $(FFLAGS) $(FREEFLAGS) -fno-range-check $<
infnan.o: infnan.F90
	$(FC) -c $(INCLDIR) $(INCS) $(FFLAGS) $(FREEFLAGS) -fno-range-check $<
nanMod.o: nanMod.F90
	$(FC) -c $(INCLDIR) $(INCS) $(FFLAGS) $(FREEFLAGS) -fno-range-check $<
geopk.o:geopk.F90
	$(FC) -c $(INCLDIR) $(INCS) $(FFLAGS) $(FREEFLAGS) -fcray-pointer $<
