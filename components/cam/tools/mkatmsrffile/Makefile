include $(CASEROOT)/Macros.bluefire

LIBDIR=/ptmp/jedwards/SMS_D.f09_g16.B1850CN.bluefire.C.083815/lib
LIBS = -L$(LIBDIR) -lcsm_share -lpio -lmct -lmpeu -L$(LIB_NETCDF) -lnetcdf -L$(LIB_PNETCDF) -lpnetcdf
MODDIR=$(LIBDIR)/include


mkatmsrffile: mkatmsrffile.F90 
	$(FC) $(FFLAGS) -I$(MODDIR) $< -o $@  $(LIBS) $(LDFLAGS) 


