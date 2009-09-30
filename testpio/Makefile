#
# Defined externally (Makefile.conf and Makefile.arch)
#
#   INCLUDES
#   LIBS
#   MPICC
#   MPIF90
#   COPTS
#   FOPTS
#   CFLAGS
#   FFLAGS
#   AWK
#   AR
#

RM=/bin/rm -f
FDEPENDS=fdepends.awk
DEPSUF = .d
MODSUF = .mod
CPPSUF = .f90


ifeq (,$(PIOARCH))
  PIOARCH=conf
endif
include ../pio/Makefile.$(PIOARCH)

SRCS_C=
SRCS_F90= kinds_mod.F90 \
          check_mod.F90 \
          utils_mod.F90 \
          gdecomp_mod.F90 \
	  namelist_mod.F90 \
          testpio.F90

SRCS2_C=
SRCS2_F90= kinds_mod.F90 \
           gdecomp_mod.F90 \
	   testdecomp.F90 

OBJS=  $(SRCS_C:.c=.o) \
       $(SRCS_F90:.F90=.o)

OBJS2=  $(SRCS2_C:.c=.o) \
       $(SRCS2_F90:.F90=.o)

# File foo.DEPSUF will contain make dependency for foo.F90
DEPENDS := $(SRCS_F90:.F90=$(DEPSUF))


#
#
#

EXE=testpio
EXE2=testdecomp
LEVEL=0

all: depends
	$(MAKE) -f Makefile LEVEL=1 $(EXE)
	$(MAKE) -f Makefile LEVEL=1 $(EXE2)

depends: $(DEPENDS)
	@echo "Done updating dependencies"

PIODIR = ../pio
PIOINC = -I$(PIODIR)
ifeq ($(TIMING),yes)
  PIOLIB = -L$(PIODIR)/../timing -ltiming -L$(PIODIR) -lpio
  INCLUDES += -I$(PIODIR)/timing
else
  PIOLIB = -L$(PIODIR) -lpio
endif
PIOLIBA = $(PIODIR)/libpio.a

$(EXE): $(OBJS) $(PIOLIBA)
	$(MPIF90) $(FFLAGS) -o $@ $(OBJS) $(PIOLIB) $(LDLIBS)

$(EXE2): $(OBJS2) $(PIOLIBA)
	$(MPIF90) $(FFLAGS) -o $@ $(OBJS2) $(PIOLIB) $(LDLIBS)

#
# Suffix rules
#

.SUFFIXES:
.SUFFIXES: .o .c .F90 $(DEPSUF)


# For each file foo.F90, make the depency file foo.d

%$(DEPSUF): %.F90
	@echo 'Making dependencies for' $< '-->' $@
	@$(AWK) -f $(FDEPENDS) -v NAME=$(basename $<) -v SUF=$(suffix $<) $< > $@



ifeq ($(EXPLICIT_CPP),yes)
SRCS_CPP= $(SRCS_F90:.F90=$(CPPSUF))
.F90.o:
	@if [ -w $*.f90 ] ; then echo "ERROR: file $*.f90 is writable - the .f90 suffix is reserved for temporary cpp output" ; exit 1; fi
	$(RM) $*.f90
	$(CPP) $(CPPFLAGS) $(CFLAGS) $(COPTS) $(INCLUDES) -o $*.f90 $*.F90 
	chmod a-w $*.f90
	$(MPIF90) -c $(FFLAGS) $(FOPTS) $(INCLUDES) $(PIOINC) $*.f90
else
.F90.o:
	$(MPIF90) -c $(FFLAGS) $(FOPTS) $(INCLUDES) $(PIOINC) $*.F90
endif


.c.o:
	$(MPICC) -c $(CFLAGS) $(COPTS) $(INCLUDES) $(PIOINC) $*.c


DATA= foo.binr foo.i4.nc foo.r4.binr foo.r8.bin  foo.r8.nc \
      foo.i4.binr  foo.nc  foo.r4.nc foo.r8.binr
     
clean:
	$(RM) $(EXE) $(EXE2) $(OBJS) $(OBJS2) *.mod $(DEPENDS) $(DATA) $(SRCS_CPP)


#
# Automatically generated module dependencies
#

ifneq (0,$(LEVEL))
  include $(DEPENDS)
endif

