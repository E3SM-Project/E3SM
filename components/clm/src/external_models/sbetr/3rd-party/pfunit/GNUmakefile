.PHONY: tests all install documentation config

TOP_DIR ?=$(shell pwd)

DOC_DIR  = $(TOP_DIR)/documentation
EXAMPLES_DIR  = $(TOP_DIR)/Examples
SOURCE_DIR  = $(TOP_DIR)/source
TESTS_DIR   = $(TOP_DIR)/tests
INCLUDE_DIR = $(TOP_DIR)/include
VPATH      += $(SOURCE_DIR) $(INCLUDE_DIR)

# Set DOXYGEN to the desired executable.
# NOTE: Doxygen Version 1.8.5 does not respect CamelCase names from
# Fortran source code by currently converting all to lowercase.  It
# does this to get HTML links correct for references in the source
# code that also do not respect the CamelCase convention.  The Fortran
# standard specifies case insensitivity.
#
# DOXYGEN = /opt/local/share/doxygen/doxygen-1.7.6/bin/doxygen
# DOXYGEN = /opt/local/share/doxygen/doxygen-1.7.5.1/bin/doxygen
DOXYGEN ?= doxygen

#
# Set PYTHON executable if needed.  Python scripts all use
# /usr/bin/env to call python. Explicit calls to python in the
# makefiles now use $(PYTHON).  If python defaults to version 3 on a
# system one can set PYTHON to python2 here. Or one can create a
# symbolic link to python2 from python in the PATH. For example:
#
# ln -s /usr/bin/python2 ~/bin/python
# ln -s /usr/bin/python2-config ~/bin/python-config
# export PATH=~/bin:${PATH}
#
# Then build pFUnit.
#
PYTHON = python

# Determine operating system, architecture and compiler
# automatically if possible

UNAME ?=$(shell uname)
ifeq ($(UNAME),)
  UNAME =UNKNOWN
else
# Check for Windows/CYGWIN compilation.
ifneq (,$(findstring CYGWIN,$(UNAME)))
UNAME =Windows
endif
endif

# ARCH  ?=$(shell arch)
ARCH  ?=$(shell uname -m)
ifeq ($(ARCH),)
  ARCH =UNKNOWN
endif

ifneq ($(UNAME),Windows)
# Also set the archiver and RANLIB options.
NULL :=
SPACE := ${NULL} ${NULL}
AR = ar -r$(SPACE)
RANLIB ?= ranlib
OUTPUT_FLAG ?= -o$(SPACE)
else
# Also set the archiver and RANLIB options.
AR = lib /out:
RANLIB ?= echo
OUTPUT_FLAG ?= /nologo /Fe
endif

# Set the relevant file extensions
include $(INCLUDE_DIR)/extensions.mk


# 32/64 ABI - almost all architectures are now 64 bit
ifeq ($(ARCH),i386)
  PFUNIT_ABI ?=32
else
  PFUNIT_ABI ?=64
endif

# Other defaults
MPI ?=NO # do not include MPI capabilities
OPENMP ?=NO # do not include OpenMP threading
# MAX_RANK ?=$(PFUNIT_MAX_RANK) # Maximum rank of arrays in assertions
MAX_RANK ?=5
# Remove the following at version 4.
ifneq ($(origin PFUNIT_MAX_RANK), undefined)
   MAX_RANK=$(PFUNIT_MAX_RANK)
endif
# The following supersedes the previous.
ifneq ($(origin PFUNIT_MAX_ARRAY_RANK), undefined)
   MAX_RANK=$(PFUNIT_MAX_ARRAY_RANK)
endif

ifneq ($(UNAME),Windows)
ROBUST ?=YES # for now include RobustRunner by default
else
ROBUST = NO
endif

# F90 Vendor common elements (override below)
FFLAGS ?=
D=-D
I=-I
MOD=-I
DEBUG_FLAGS =-g

# Default compiler by architecture - always gfortran for now:
F90 ?=gfortran
F90_VENDOR ?=GNU

# F90 Vendor specifics
# Possibly F90 defined - makes things simple:

ifneq (,$(findstring $(F90), ifort gfortran nag nagfor pgfortran xlf))
  ifeq ($(F90),ifort)
     COMPILER=Intel
  else ifeq ($(F90),gfortran)
     COMPILER=GNU
  else ifeq ($(F90),nagfor)
     COMPILER=NAG
  else ifeq ($(F90),pgfortran)
     COMPILER=PGI
  else ifneq (,$(findstring $(F90),xlf))
     COMPILER=IBM
	else
		COMPILER=UNKNOWN
  endif
# Override F90_VENDOR with COMPILER
	F90_VENDOR=$(COMPILER)
else # use F90_VENDOR to specify
  ifneq (,$(findstring $(F90_VENDOR),INTEL Intel intel ifort))
    COMPILER=Intel
  else ifneq (,$(findstring $(F90_VENDOR),GNU gnu gfortran GFortran GFORTRAN))
    COMPILER=GNU
  else ifneq (,$(findstring $(F90_VENDOR),nag NAG nagfor))
    COMPILER=NAG
  else ifneq (,$(findstring $(F90_VENDOR),pgi PGI pgfortran))
    COMPILER=PGI
  else ifneq (,$(findstring $(F90_VENDOR),ibm IBM xlf XLF))
    COMPILER=IBM
  endif
endif

# F90_VENDOR is no longer needed after this point.  We keep it around
# until we can verify that it's not needed in subdirectories or for
# recursive calls. TODO:  Check F90_VENDOR usage.

ifneq ($(findstring $(MPI),yes YES Yes),)
  USEMPI=YES
  MPIF90 ?= mpif90
  MPIRUN ?= mpirun
  FPPFLAGS += $DUSE_MPI
  CPPFLAGS += -DUSE_MPI
  ifeq ($(MPICH),YES)
     LIBMPI ?=-lmpich
  else
# The following may be redundant and better handled via an MPI's linking script. 2013-1104 MLR
#     LIBMPI ?=-lmpi
     LIBMPI ?=
  endif
  LDFLAGS += $(LIBMPI)
endif

ifneq ($(findstring $(OPENMP),yes YES Yes),)
  USEOPENMP=YES
endif

ifneq ($(findstring $(ROBUST),yes YES Yes),)
  BUILDROBUST=YES
  FPPFLAGS += $DBUILD_ROBUST
  CPPFLAGS += -DBUILD_ROBUST
endif

FPPFLAGS += $D$(COMPILER) $D$(UNAME)
CPPFLAGS += -D$(COMPILER) -D$(UNAME) -I$(INCLUDE_DIR)

# FPPFLAGS += $D$(F90_VENDOR) $D$(UNAME)
# CPPFLAGS += -D$(F90_VENDOR) -D$(UNAME) -I$(INCLUDE_DIR)

ifeq ($(PFUNIT_ABI),64)
  FPPFLAGS += $DLONG_PTR
  CPPFLAGS += -DLONG_PTR
endif

FFLAGS +=$I$(INCLUDE_DIR) $(MOD)$(SOURCE_DIR)
CFLAGS +=-I$(INCLUDE_DIR)

ifeq ($(DEBUG),YES)
        FFLAGS += $(DEBUG_FLAGS)
endif

all: include/configuration.mk
	$(MAKE) -C $(SOURCE_DIR) all
	$(MAKE) -C $(TESTS_DIR) all

documentation:
	$(DOXYGEN) documentation/doxygen.conf

documentation/pFUnit3-ReferenceManual.pdf: documentation
	$(MAKE) -C documentation/latex all
	mv -f documentation/latex/refman.pdf documentation/pFUnit3-ReferenceManual.pdf


clean: local-top1-clean local-top1-cleanExamples

local-top1-clean: local-top1-cleanExamples
	$(MAKE) -C $(SOURCE_DIR) clean
	$(MAKE) -C $(TESTS_DIR) clean
	\rm -f include/configuration.mk test.xml

local-top1-cleanExamples:
	\rm -f Examples/Simple/tests.xml
	tools/clean Examples

distclean: local-top1-distclean

local-top1-distclean: local-top1-cleanExamples
	$(MAKE) -C $(SOURCE_DIR) distclean
	$(MAKE) -C $(TESTS_DIR) distclean
	$(MAKE) -C $(DOC_DIR) distclean
	\rm -f include/configuration.mk test.xml

tests: all
ifeq ($(USEMPI),YES)
	$(MPIRUN) -np 4 ./tests/tests$(EXE_EXT)
else
	./tests/tests$(EXE_EXT)
endif

develop:
	cp -f $(TOP_DIR)/include/base-develop.mk $(TOP_DIR)/include/base.mk

install: libpfunit$(LIB_EXT)
ifndef INSTALL_DIR
	$(error Must specify INSTALL_DIR. Example: make install INSTALL_DIR=SOME_PATH, \
	where SOME_PATH is different than $(TOP_DIR))
else
ifeq ($(INSTALL_DIR),$(TOP_DIR))
        $(error INSTALL_DIR cannot be the same as TOP_DIR)
endif
ifeq ($(INSTALL_DIR),.)
        $(error INSTALL_DIR cannot be the same as TOP_DIR)
endif
	@echo Installing pFUnit in $(INSTALL_DIR)
	tools/install $(INSTALL_DIR)/lib source/lib*
	tools/install $(INSTALL_DIR)/mod source/*.mod
	tools/install $(INSTALL_DIR) include
	mv -f $(INSTALL_DIR)/include/base-install.mk $(INSTALL_DIR)/include/base.mk
	tools/install $(INSTALL_DIR) bin
	@echo +++
	@echo PFUNIT has been installed in $(INSTALL_DIR).
	@echo For normal usage please ensure PFUNIT is set to $(INSTALL_DIR).
	@echo For example, in bash:  export PFUNIT=$(INSTALL_DIR)
endif

include/configuration.mk:
	@echo "# include/configuration.mk generated automatically during build" \
		> include/configuration.mk
	@echo COMPILER  ?= $(COMPILER) >> include/configuration.mk
	@echo USEOPENMP ?= $(USEOPENMP) >> include/configuration.mk
	@echo USEMPI    ?= $(USEMPI) >> include/configuration.mk
	@echo MAX_RANK 	?= $(MAX_RANK) >> include/configuration.mk
	@echo BUILDROBUST ?= $(BUILDROBUST) >> include/configuration.mk
	@echo VERSION = \"`cat VERSION`\" >> include/configuration.mk

export UNAME
export OBJ_EXT
export EXE_EXT
export LIB_EXT
export AR
export RANLIB
export OUTPUT_FLAG
export F90
export F90_VENDOR
export FFLAGS
export FPPFLAGS
export CPP
export CFLAGS
export CPPFLAGS
export LDFLAGS
export SOURCE_DIR
export INCLUDE_DIR
export VPATH
export MPI
export USEMPI
export USEOPENMP
export MAX_RANK
export BUILDROBUST
export MPIF90
export LIBMPI
export COMPILER
export PYTHON

ifeq ($(DEBUG),YES)
  $(warning Compilation configuration is as follows:)
  $(warning     UNAME:  $(UNAME))
  $(warning     ARCH:   $(ARCH))
  $(warning     F90 vendor:     $(COMPILER))
  $(warning     F90 command:    $(F90))
  $(warning     USE MPI:        $(MPI))
  $(warning     ABI:            $(PFUNIT_ABI))
  $(warning File extensions:)
  $(warning     OBJ_EXT         $(OBJ_EXT))
  $(warning     LIB_EXT         $(LIB_EXT))
  $(warning     EXE_EXT         $(EXE_EXT))
endif
