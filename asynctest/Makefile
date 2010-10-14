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
FDEPENDS=../pio/fdepends.awk
DEPSUF = .d
MODSUF = .mod
CPPSUF = .f90


ifeq (,$(PIOARCH))
  PIOARCH=conf
endif
include ../pio/Makefile.$(PIOARCH)

SRCS_C=
SRCS_F90= pioasync.F90


OBJS=  $(SRCS_C:.c=.o) \
       $(SRCS_F90:.F90=.o)

# File foo.DEPSUF will contain make dependency for foo.F90
DEPENDS := $(SRCS_F90:.F90=$(DEPSUF))


#
#
#

EXE=pioasync

LEVEL=0

default: depends
	$(MAKE) -f Makefile LEVEL=1 $(EXE)

all: depends
	$(MAKE) -f Makefile LEVEL=1 $(EXE)

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


DATA=
     
clean:
	$(RM) $(EXE)  $(OBJS)  *.mod $(DEPENDS) $(DATA) $(SRCS_CPP)


#
# Automatically generated module dependencies
#

ifneq (0,$(LEVEL))
  include $(DEPENDS)
endif

