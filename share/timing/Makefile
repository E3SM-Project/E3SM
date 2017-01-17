.SUFFIXES:
.SUFFIXES: .F90 .o .c .f90
# name of macros file - but default this is generic

VPATH := $(GPTL_DIR)

# Determine whether to compile threaded or not
ifeq ($(strip $(BUILD_THREADED)),TRUE)
   compile_threaded = true
endif
ifeq ($(strip $(SMP)),TRUE)
   compile_threaded = true
endif


ifndef MOD_SUFFIX
   MOD_SUFFIX := mod
endif


ifeq ($(strip $(MACFILE)),)
   MACFILE  := Macros
endif

# Machine specific macros file 
# This must be included before any settings are overwritten 
# But must be AFTER any definitions it uses are defined.
# So be careful if moving this either earlier or later in the makefile!!!
include $(MACFILE)

INCLDIR += -I$(GPTL_DIR)

ifeq ($(strip $(MPILIB)), mpi-serial)
  CC      := $(SCC)
  FC      := $(SFC)
  MPIFC   := $(SFC)
  MPICC   := $(SCC)
  INCLDIR += -I$(GPTL_LIBDIR)/../mct/mpi-serial
else
  CC := $(MPICC)
  FC := $(MPIFC)
  CPPDEFS += -DHAVE_MPI
endif
ifdef CPRE
  FPPDEFS := $(patsubst -D%,$(CPRE)%,$(CPPDEFS)) 
else
  FPPDEFS := $(CPPDEFS)
endif



OBJS = gptl.o GPTLutil.o GPTLget_memusage.o GPTLprint_memusage.o \
       gptl_papi.o f_wrappers.o perf_mod.o perf_utils.o


libgptl.a: $(OBJS)
	$(AR) ruv $@ $(OBJS)



.c.o:
	$(CC) -c $(INCLDIR) $(INCS) $(CFLAGS) $(CPPDEFS) $<
.F.o:
	$(FC) -c $(INCLDIR) $(INCS) $(FFLAGS) $(FPPDEFS) $(FIXEDFLAGS) $<
.f90.o:
	$(FC) -c $(INCLDIR) $(INCS) $(FFLAGS) $(FREEFLAGS)  $<
.F90.o:
	$(FC) -c $(INCLDIR) $(INCS) $(FFLAGS) $(FPPDEFS) $(FREEFLAGS)  $<

mostlyclean:
	$(RM) -f *.f *.f90 

clean:
	$(RM) -f *.f *.f90 *.d *.$(MOD_SUFFIX) $(OBJS)


install: libgptl.a
	cp -p $(GPTL_DIR)/gptl.h $(SHAREDPATH)/include
	cp -p *.$(MOD_SUFFIX) $(SHAREDPATH)/include
	cp -p libgptl.a $(SHAREDPATH)/lib


perf_mod.o: perf_utils.o
f_wrappers.o: gptl.h private.h
f_wrappers_pmpi.o: gptl.h private.h
gptl.o: gptl.h private.h
util.o: gptl.h private.h
gptl_papi.o: gptl.h private.h
pmpi.o: gptl.h private.h
