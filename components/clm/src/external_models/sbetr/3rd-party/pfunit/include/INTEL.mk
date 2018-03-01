F90 ?=ifort

I=-I
M=-I
L=-L

ifneq ($(UNAME),Windows)
# Non-Windows (Linux) command line options for the intel compiler
version = $(shell $(F90) --version | grep -E '\(IFORT\) 13')

FFLAGS += -assume realloc_lhs
FFLAGS += -g -O0 -traceback -check uninit -check bounds -check stack -check uninit

ifeq ($(USEOPENMP),YES)
FFLAGS += -openmp
endif


else
# Windows command line options for the intel compiler
version = $(shell $(F90)  2>&1 | head -1 | grep 'Version 13')

# Suppress version information with each compile.
FFLAGS += /nologo
FFLAGS += /assume:realloc_lhs
FFLAGS += /Z7 /Od /traceback /check:uninit /check:bounds /check:stack /check:uninit
# Enable the Fortran preprocessor
FFLAGS += /fpp

# Remove the DEBUG_FLAGS -g option.
DEBUG_FLAGS = /Z7
endif


# Common command line options.

F90_PP_ONLY = -E
F90_PP_OUTPUT = >

CPPFLAGS +=-DIntel
FPPFLAGS +=-DIntel

# Check if the version of the compiler is 13
ifneq ($(version),)
  CPPFLAGS+=-DINTEL_13
  FPPFLAGS+=-DINTEL_13
endif

ifeq ($(USEOPENMP),YES)
FFLAGS += -openmp
LIBS += -openmp
endif
