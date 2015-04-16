include ../macros.make
# *xlf* (IBM) compilers use a different mechanism for #defines than other compilers
ifeq ($(findstring xlf, $(FC)),xlf)
  DEFINE = -WF,-D
else
  DEFINE = -D
endif

# May be reset below
THREADED = no

# Ensure that FFLAGS are set for low optimization
FFLAGS += -g -O1

# TESTS defines which of EXES have a built-in testing procedure.
TESTS = testbasics nlreader outoforder

ifeq ($(HAVE_PAPI),yes)
  TESTS += testinit testpapi
endif
EXES = $(TESTS) errtest utrtest overhead

LDFLAGS = -g -L.. -l$(LIBNAME) $(ABIFLAGS)

#For gptlf.mod
FFLAGS += -I..

ifeq ($(HAVE_MPI),yes)
  EXES += summary pmpi
  FFLAGS  += $(DEFINE)HAVE_MPI $(MPI_INCFLAGS)
  LDFLAGS += $(MPI_LIBFLAGS)
endif

ifeq ($(ENABLE_PMPI),yes)
  FFLAGS += $(DEFINE)ENABLE_PMPI
  LIBNAME = gptl_pmpi
  ifeq ($(HAVE_IARGCGETARG),yes)
    FFLAGS += $(DEFINE)HAVE_IARGCGETARG
  endif
else
  LIBNAME = gptl
endif

ifeq ($(OPENMP),yes)
  THREADED = yes
  EXES    += toomanythreads
  FFLAGS  += $(DEFINE)THREADED_OMP $(FOMPFLAG)
  LDFLAGS += $(FOMPFLAG)
else
  ifeq ($(PTHREADS),yes)
# Threaded tests use OpenMP
    THREADED = yes
    EXES    += toomanythreads
    FFLAGS  += $(DEFINE)THREADED_OMP $(FOMPFLAG)
    LDFLAGS += $(FOMPFLAG) -lpthread
  endif
endif

ifeq ($(HAVE_PAPI),yes)
  FFLAGS  += $(DEFINE)HAVE_PAPI
  LDFLAGS += $(PAPI_LIBFLAGS)
endif

OBJS = $(patsubst %,%.o,$(EXES))

%.o: %.F90
	$(FC) $(FFLAGS) -c $<

all: $(EXES)

clean:
	$(RM) $(OBJS) *.mod $(EXES) timing.*

errtest: errtest.o
	$(FC) $(FFLAGS) -o $@ $< $(LDFLAGS)

utrtest: utrtest.o
	$(FC) $(FFLAGS) -o $@ $< $(LDFLAGS)

summary: summary.o
	$(FC) $(FFLAGS) -o $@ $< $(LDFLAGS)

testinit: testinit.o
	$(FC) $(FFLAGS) -o $@ $< $(LDFLAGS)

testpapi: testpapi.o
	$(FC) $(FFLAGS) -o $@ $< $(LDFLAGS)

testbasics: testbasics.o
	$(FC) $(FFLAGS) -o $@ $< $(LDFLAGS)

nlreader: nlreader.o
	$(FC) $(FFLAGS) -o $@ $< $(LDFLAGS)

outoforder: outoforder.o
	$(FC) $(FFLAGS) -o $@ $< $(LDFLAGS)

pmpi: pmpi.o
	$(FC) $(FFLAGS) -o $@ $< $(LDFLAGS)
pmpi.o: pmpi.F90
	$(FC) -c $(FFLAGS) $<

overhead: overhead.F90 ../lib$(LIBNAME).a
	$(FC) $(FFLAGS) -o $@ $< $(LDFLAGS)

toomanythreads: toomanythreads.F90
	$(FC) $(FFLAGS) -o $@ $< $(LDFLAGS)

# Hidden executable testbacktrace requires some strange flags so don't build it unless asked
testbacktrace: testbacktrace.F90
	$(FC) $(FFLAGS) $(INSTRFLAG) -rdynamic -o $@ $< $(LDFLAGS)

# Built-in tests to ensure that GPTL is behaving as expected
# Invoke as "make test" from one dir above

test: $(TESTS)
	@$(RM) timing.0
	@echo Running ./testbasics...
	./testbasics
	@echo Success
	@echo
	@echo Checking contents of timing.0 for a specific region name...
	@grep -q testbasics timing.0 || (echo "Failed" && exit 1)
	@echo Success
	@echo
	@echo Running ./nlreader...
	./nlreader
	@echo Success
	@echo Running ./outoforder...
	./outoforder
	@echo Success
ifeq ($(THREADED),yes)
	@echo Running ./toomanythreads...
	./toomanythreads
endif
ifeq ($(HAVE_PAPI),yes)
	@echo
	@echo Running ./testinit...
	./testinit
	@echo
	@echo Running ./testpapi...
	./testpapi
endif
ifeq ($(ENABLE_PMPI),yes)
	@echo
	@echo Running $(MPICMD) ./pmpi
	env $(MPICMD) ./pmpi
	@echo Checking contents of timing.0 for region name MPI_Send...
	@grep -q MPI_Send timing.0 || (echo "Failed" && exit 1)
	@echo Success
endif
	@echo "All Fortran tests passed"
