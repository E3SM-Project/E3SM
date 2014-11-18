.SUFFIXES: .F .c .o

all: esmf_time

esmf_time:
	( cd esmf_time_f90; $(MAKE) FC="$(FC) $(FFLAGS)" CPP="$(CPP)" CPPFLAGS="$(CPPFLAGS) -DHIDE_MPI" )

clean:
	( cd esmf_time_f90; $(MAKE) clean )
