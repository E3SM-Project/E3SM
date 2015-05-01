# TESTS defines which EXES to run when "make test" is done.
# Always run "summary". Add others depending on settings.

include ../macros.make

# Ensure that CFLAGS are set for no optimization
CFLAGS += -g -O0

TESTS = global
ifeq ($(HAVE_MPI),yes)
  TESTS += summary pmpi
endif

ifeq ($(HAVE_PAPI),yes)
  TESTS += testpapi
endif

ifeq ($(ENABLE_PMPI),yes)
  CFLAGS += -DENABLE_PMPI
  LIBNAME = gptl_pmpi
  ifeq ($(HAVE_IARGCGETARG),yes)
    CFLAGS += -DHAVE_IARGCGETARG
  endif
else
  LIBNAME = gptl
endif

LDFLAGS = $(CFLAGS) -L.. -l$(LIBNAME) $(ABIFLAGS) -lm

ifeq ($(OPENMP),yes)
  CFLAGS  += -DTHREADED_OMP $(COMPFLAG)
  LDFLAGS += $(COMPFLAG)
else
  ifeq ($(PTHREADS),yes)
# Threaded tests use OpenMP
    CFLAGS  += -DTHREADED_OMP $(COMPFLAG)
    LDFLAGS += -lpthread $(COMPFLAG)
  endif
endif

ifeq ($(HAVE_PAPI),yes)
  CFLAGS  += -DHAVE_PAPI $(PAPI_INCFLAGS)
  LDFLAGS += $(PAPI_LIBFLAGS)
endif

ifeq ($(HAVE_MPI),yes)
  CFLAGS  += -DHAVE_MPI $(MPI_INCFLAGS)
  LDFLAGS += $(MPI_LIBFLAGS)
endif

# Auto-profiling 
ifeq ($(TEST_AUTOPROFILE),yes)
  TESTS += cygprofile profcxx
# Initialize OBJS: cygprofilesubs.o is the one .o file needed that patsubst can't
# set automatically
  OBJS = cygprofilesubs.o
endif

# EXES is TESTS plus some other sample codes.
# Items in EXES but not TESTS are built by "make all" but
# not run by "make test"

EXES = $(TESTS) papiomptest gptl_gran_overhead printwhileon

UTILEXES = gptl_gran_overhead
ifeq ($(HAVE_PAPI),yes)
  UTILEXES += gptl_avail gptl_knownflopcount
endif
EXES += $(UTILEXES)

OBJS += $(patsubst %,%.o,$(EXES))

all: $(EXES)

clean:
	$(RM) *.o $(EXES) timing.*

gptl_avail: avail.o
	$(CC) -o $@ $< $(LDFLAGS)

gptl_knownflopcount: knownflopcount.o
	$(CC) -o $@ $< $(LDFLAGS)

gptl_gran_overhead: gran_overhead.o
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

papiomptest: papiomptest.o
	$(CC) $(THREADFLAGS) -o $@ $< $(LDFLAGS)

cygprofile: cygprofile.o cygprofilesubs.o
	$(CC) -o $@ cygprofile.o cygprofilesubs.o $(LDFLAGS) $(INSTRFLAG)
cygprofilesubs.o: cygprofilesubs.c
	$(CC) -c $(CFLAGS) $(INSTRFLAG) -o $@ $<

profcxx: profcxx.o
	$(CXX) -o $@ $< $(LDFLAGS) $(INSTRFLAG)
profcxx.o: profcxx.C myclasses.h
	$(CXX) -c $(CFLAGS) $(INSTRFLAG) -o $@ $<

summary: summary.o
	$(CC) -o $@ $< $(LDFLAGS)

testpapi: testpapi.o
	$(CC) -o $@ $< $(LDFLAGS)

pmpi: pmpi.o
	$(CC) -o $@ $< $(LDFLAGS)

global: global.o
	$(CC) -o $@ $< $(LDFLAGS)

printwhileon: printwhileon.o
	$(CC) -o $@ $< $(LDFLAGS)

# Built-in tests to ensure that GPTL is behaving as expected
# Invoke as "make test" from one dir above

test: $(TESTS)
	@echo

ifeq ($(HAVE_PAPI),yes)
	@./testpapi
	@echo
endif

ifeq ($(HAVE_MPI),yes)
	@echo Running env OMP_NUM_THREADS=2 $(MPICMD) ./summary...
	@env OMP_NUM_THREADS=2 $(MPICMD) ./summary || (echo "Failure to run" && exit 1)
	@echo Success
	@echo "Shell is checking contents of timing.summary produced by ./summary for a specific region name..."
	@grep -q sleep timing.summary || (echo "Failure" && exit 1)
	@echo Success
ifeq ($(ENABLE_PMPI),yes)
	@echo
	@echo Running $(MPICMD) ./pmpi...
	@env $(MPICMD) ./pmpi || (echo "Failure" && exit 1)
	@echo Success
endif
endif

ifeq ($(TEST_AUTOPROFILE),yes)
	@echo
	@echo "Testing auto-profiling..."
	./cygprofile || (echo "Failure to run" && exit 1)
	@echo Success at running
	@echo "Running hex2name.pl | grep callutil100times..."
	(../hex2name.pl ./cygprofile timing.0 | grep -q callutil100times) && echo "Success" || (echo "Failure" && exit 1)
endif
	@echo
	@echo "All C and C++ tests passed"
	@echo

install: $(UTILEXES)
ifneq ($(UTILEXES),)
	install -m 0755 $(UTILEXES) $(INSTALLDIR)/bin && exit 0
else
	@exit 0
endif
