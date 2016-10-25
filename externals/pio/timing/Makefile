#
# Defined externally in Makefile.arch for each platform
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

RM=/bin/rm -f
MODSUF = .mod
CPPSUF = .f90


ifeq (,$(PIOARCH))
  PIOARCH=conf
endif
include ../pio/Makefile.$(PIOARCH)
export PIOARCH

SRCS_C = GPTLget_memusage.c \
         GPTLprint_memusage.c \
         GPTLutil.c \
         f_wrappers.c \
         gptl.c \
         gptl_papi.c \
         threadutil.c \

SRCS_F90 =  perf_mod.F90 \
            perf_utils.F90  

OBJS=  $(SRCS_C:.c=.o) \
       $(SRCS_F90:.F90=.o)


MODFILES := $(SRCS_F90:.F90=$(MODSUF))

PERL = /usr/bin/perl

LIB= libtiming.a

all: $(LIB)


#
# Suffix rules
#

.SUFFIXES:
.SUFFIXES: .o .c .F90 $(DEPSUF)



ifeq ($(EXPLICIT_CPP),yes)
SRCS_CPP= $(SRCS_F90:.F90=$(CPPSUF))
.F90.o:
	@if [ -w $*.f90 ] ; then echo "ERROR: file $*.f90 is writable - the .f90 suffix is reserved for temporary cpp output" ; exit 1; fi
	$(RM) $*.f90
	$(CPP) $(CPPDEFS) $(CFLAGS) $(COPTS) $(INCLUDES) -o $*.f90 $*.F90 
	chmod a-w $*.f90
	$(MPIFC) -c $(FFLAGS) $(FOPTS) $(INCLUDES) $*.f90
else
SRCS_CPP=
.F90.o:
	$(MPIFC) -c $(FFLAGS) $(FOPTS) $(INCLUDES) $*.F90
endif

.c.o:
	$(MPICC) -c $(CFLAGS) $(CPPDEFS) $(COPTS) $(INCLUDES) $*.c


$(LIB): $(OBJS)
	$(RM) $@
	$(AR) $(ARFLAGS) $@ $(OBJS)

     
clean:
	$(RM) $(LIB) $(OBJS) $(MODFILES) $(SRCS_CPP)


perf_mod.o: perf_utils.o
perf_mod.o: gptl.inc
perf_utils.o: gptl.inc
