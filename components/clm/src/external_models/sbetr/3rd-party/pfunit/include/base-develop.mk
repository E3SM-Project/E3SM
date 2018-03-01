BASEMK_INCLUDED=YES

SRC_DIR     =$(TOP)/source
TESTS_DIR   =$(TOP)/tests
INCLUDE_DIR =$(TOP)/include
LIB_DIR     =$(TOP)/source
MOD_DIR     =$(TOP)/source

# Set the required file extensions.
include $(INCLUDE_DIR)/extensions.mk

COMPILER_ = $(shell echo $(COMPILER) | tr a-z A-Z )

# Include the compiler-specific options.
include $(INCLUDE_DIR)/$(COMPILER_).mk

FFLAGS += $I$(INCLUDE_DIR)

ifeq ($(USEMPI),)
  FC=$(F90)
else
  override FC=$(MPIF90)
endif

%$(OBJ_EXT): %.F90
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $@ $<

.PHONY: clean distclean

clean:
	-$(RM) *$(OBJ_EXT) *.mod *.i90 *~ *.tmp *.dbg
	-$(RM) -r *.dSYM

distclean: clean
	-$(RM) *$(LIB_EXT) *$(EXE_EXT) dependencies.inc

export FC
export BASEMK_INCLUDED
