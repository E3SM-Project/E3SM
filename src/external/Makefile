.SUFFIXES: .F .c .o

all: esmf_time ezxml

esmf_time:
	( cd esmf_time_f90; $(MAKE) FC="$(FC) $(FFLAGS)" CPP="$(CPP)" CPPFLAGS="$(CPPFLAGS) -DHIDE_MPI" )

ezxml:
	( cd ezxml; $(MAKE) )

clean:
	( cd esmf_time_f90; $(MAKE) clean )
	( cd ezxml; $(MAKE) clean )
