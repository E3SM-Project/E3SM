# Makefile for pegasus.2.6.1
#
#----------------------------------------------------------------------
#
# includes
#
include ./make.in
include ./cpp.in

all: exe
	touch .success

#
# Suffix rules...
#
.SUFFIXES : .F .f .F90 .f90 .o

.F.o:
	$(CPP) $(CPPFLAGS) $*.F $*.f ;  $(FC) -c $(FCFLAGS) $*.f

.F.f:
	$(CPP) $(CPPFLAGS) $*.F $*.f

.f.o:
	$(CPP) $(CPPFLAGS) $*.F $*.f ;  $(FC) -c $(FCFLAGS) $*.f

.F90.o:
	$(FC) $(CPPFLAGS) -c $(FCFLAGS) $*.F90

.F90.f90:
	$(RM) $*.mod ; $(CPP) $(CPPFLAGS) $*.F90 

.f90.o:
	$(RM) $*.mod ; $(CPP) $(CPPFLAGS) $*.F90 ;  $(FC) -c $(FCFLAGS) $*.f90

#
# Object lists...
#
# OBJ1, OBJ2, and OBJ4 files 
#    These files are special versions of the equivalent CAM5 modules.
#    They provide the functionality needed by the cambox offline code
#    that is used for development and testing of the modal aerosol module (MAM),
#    but (in most cases) not all the functionality of the equivalent CAM5 module.
#    Also, they may have been taken from a version of CAM5 that was older
#    than ACME-V0 (i.e., pre 2014).
# OBJ3 and OBJ5 files
#    These files should be identical to the equivalent CAM5 module (or almost so).
# OBJ9 files 
#    These are specific to the offline test code and have no CAM equivalents.

CONFIGINFO = \
	./make.in  \
	./cpp.in

OBJ1 = \
	shr_kind_mod.o \
	shr_const_mod.o	\
        cam_logfile.o  \
        cam_abortutils.o  \
        spmd_utils.o  \

OBJ2 = \
        perf_mod.o \
	ppgrid.o  \
	physconst.o  \
	constituents.o  \
	physics_types.o  \
        glb_verif_smry.o \
	ascii_io.o \
	physpkg.o  \
        test_core.o 

OBJ3 = \
	driver.o  
#-----------------------------------------------------------------------------
# WARNING: Don't touch anything below this line unless you exactly know it !!!
#-----------------------------------------------------------------------------


#
# Dependancies...
#
exe : $(OBJ1) $(OBJ2) $(OBJ3) $(CONFIGINFO)
	$(FC) -o $@ $(FCFLAGS) $(LDFLAGS) $(OBJ1) $(OBJ2) $(OBJ3) 

$(OBJ1):  $(CONFIGINFO)

$(OBJ2):  $(OBJ1)

$(OBJ3):  $(OBJ2)


clean:
	/bin/rm  exe  *.o  *.mod unitTest_1
#/bin/rm  exe  *.o  *.f  *.f90  *.mod


cleanmod:
	/bin/rm  *.x  *.mod  module*.o  module*.f90  module*.f

