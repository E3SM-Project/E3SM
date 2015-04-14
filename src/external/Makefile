.SUFFIXES: .F .c .o

all: esmf_time ezxml-lib

esmf_time:
	( cd esmf_time_f90; $(MAKE) FC="$(FC) $(FFLAGS)" CPP="$(CPP)" CPPFLAGS="$(CPPFLAGS) -DHIDE_MPI" )

ezxml-lib:
	( cd ezxml; $(MAKE) )

clean:
	( cd esmf_time_f90; $(MAKE) clean )
	( cd ezxml; $(MAKE) clean )
