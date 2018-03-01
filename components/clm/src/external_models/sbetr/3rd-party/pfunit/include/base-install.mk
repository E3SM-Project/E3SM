
BASEMK_INCLUDED=YES

INCLUDE_DIR =$(PFUNIT)/include
LIB_DIR     =$(PFUNIT)/lib
MOD_DIR     =$(PFUNIT)/mod

# These two are not part of an installation.
# SRC_DIR     =$(PFUNIT)/source
# TESTS_DIR   =$(PFUNIT)/tests

# Read in compile configuration to help set flags like -gomp for GNU.
include $(INCLUDE_DIR)/configuration.mk

# Set the required file extensions.
include $(INCLUDE_DIR)/extensions.mk

# F90 Vendor common elements (override below)
# FFLAGS ?=
D=-D
I=-I
MOD=-I
DEBUG_FLAGS =-g

# Include the compiler-specific options.
COMPILER ?= COMPILER_NOT_SET
COMPILER_ = $(shell echo $(COMPILER) | tr a-z A-Z )
include $(INCLUDE_DIR)/$(COMPILER_).mk

FFLAGS += $I$(INCLUDE_DIR)

ifeq ($(BUILDROBUST),YES)
  FPPFLAGS += $DBUILD_ROBUST
  CPPFLAGS += -DBUILD_ROBUST
endif

# include/driver.F90 needs both BUILD_ROBUST
ifneq ($(USEMPI),YES)
  FC=$(F90)
else
  FC=$(MPIF90)
endif

%$(OBJ_EXT): %.F90
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $@ $<

.PHONY: clean distclean echo

clean: local-base0-clean

local-base0-clean:
	$(RM) *$(OBJ_EXT) *.mod *.i90 *~ *.tmp *.s *.dbg
	$(RM) -r *.dSYM

distclean: local-base0-distclean

local-base0-distclean: clean
	$(RM) *$(LIB_EXT) *$(EXE_EXT)

echo:
	@echo COMPILER: $(COMPILER)
	@echo FC:	$(FC)
	@echo USEMPI:   $(USEMPI)
	@echo FFLAGS:   $(FFLAGS)
	@echo FPPFLAGS: $(FPPFLAGS)
	@echo CPPFLAGS: $(CPPFLAGS)

export FC
export BASEMK_INCLUDED
