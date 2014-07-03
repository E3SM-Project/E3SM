#
# File:  preprocess.mk
#
#----------------------------------------------------------------------------
#
#  This makefile is called from the POP driver makefile and performs
#  only the preprocessing step.
#
#----------------------------------------------------------------------------

SHELL    = /bin/sh

#----------------------------------------------------------------------------
#
#  Set valid suffixes.
#
#----------------------------------------------------------------------------

#  First clean out current list of suffixes, then define them
.SUFFIXES: 
.SUFFIXES: .c .f .f90 .F .F90 .C

#----------------------------------------------------------------------------
#
#  Include architecture-specific flags and options. 
#
#----------------------------------------------------------------------------

ifneq (,$(POPARCH))
  include $(POPDIR)/build/$(POPARCH).gnu
  export POPARCH
else
  bogus:
	@echo "  Please set POPARCH environment variable"
endif

#----------------------------------------------------------------------------
#
#  Define paths to sources in variable SRCDIRS.
#
#----------------------------------------------------------------------------

SRCDIRS = $(POPEXEDIR)/
SRCDIRS := $(SRCDIRS) $(POPDIR)/source/
SRCDIRS := $(SRCDIRS) $(POPDIR)/$(COMMDIR)/
SRCDIRS := $(SRCDIRS) $(POPDIR)/drivers/cpl_none/

#----------------------------------------------------------------------------
#
#  VPATH is the built-in symbol whose value is the path that gmake will 
#  search for dependencies.
#
#----------------------------------------------------------------------------

VPATH = $(SRCDIRS)

#----------------------------------------------------------------------------
#
# Define .F sources that must be preprocessed into the build directory as .f
# and add .f version to list of target source files.
#
#----------------------------------------------------------------------------

SOURCES = 
FSRCS   = $(strip $(foreach dir,$(SRCDIRS),$(wildcard $(dir)*.F)))
ifneq (,$(FSRCS))
  SOURCES := $(addprefix $(POPEXEDIR)/compile/, $(notdir $(FSRCS:.F=.f))) \
             $(SOURCES)
endif

#----------------------------------------------------------------------------
#
# Define .F90 sources to be preprocessed into the build directory as .f90
#
#----------------------------------------------------------------------------

F90SRCS   = $(strip $(foreach dir,$(SRCDIRS),$(wildcard $(dir)*.F90)))
ifneq (,$(F90SRCS))
  SOURCES := $(addprefix $(POPEXEDIR)/compile/, $(notdir $(F90SRCS:.F90=.f90))) \
             $(SOURCES)
endif

#----------------------------------------------------------------------------
#
# Define .C sources that must be preprocessed into the build directory as .c
#
#----------------------------------------------------------------------------

CSRCS   = $(strip $(foreach dir,$(SRCDIRS),$(wildcard $(dir)*.C)))
ifneq (,$(CSRCS))
  SOURCES := $(addprefix $(POPEXEDIR)/compile/, $(notdir $(CSRCS:.C=.c))) \
             $(SOURCES)
endif

#----------------------------------------------------------------------------
#
# Define any .f sources that need to be copied into the build directory
#
#----------------------------------------------------------------------------

LFSRCS  = $(strip $(foreach dir,$(SRCDIRS),$(wildcard $(dir)*.f)))
ifneq (,$(LFSRCS))
  ifneq (,$(FSRCS))
    LFSRCS    := $(filter-out $(FSRCS:.F=.f),$(LFSRCS))
  endif
  ifneq (,$(LFSRCS))
    SOURCES := $(addprefix $(POPEXEDIR)/compile/, $(notdir $(LFSRCS))) \
               $(SOURCES)
  endif
endif

#----------------------------------------------------------------------------
#
# Define .f90 sources that need to be copied into the build directory
#
#----------------------------------------------------------------------------

LF90SRCS  = $(strip $(foreach dir,$(SRCDIRS),$(wildcard $(dir)*.f90)))
ifneq (,$(LF90SRCS))
  ifneq (,$(F90SRCS))
    LF90SRCS    := $(filter-out $(F90SRCS:.F90=.f90),$(LF90SRCS))
  endif
  ifneq (,$(LF90SRCS))
    SOURCES := $(addprefix $(POPEXEDIR)/compile/, $(notdir $(LF90SRCS))) \
               $(SOURCES)
  endif
endif

#----------------------------------------------------------------------------
#
# Define .c sources that need to be copied into the build directory
#
#----------------------------------------------------------------------------

LCSRCS   = $(strip $(foreach dir,$(SRCDIRS),$(wildcard $(dir)*.c)))
ifneq (,$(LCSRCS))
  ifneq (,$(CSRCS))
    LCSRCS    := $(filter-out $(CSRCS:.C=.c),$(LCSRCS))
  endif
  ifneq (,$(LCSRCS))
    SOURCES := $(addprefix $(POPEXEDIR)/compile/, $(notdir $(LCSRCS))) \
               $(SOURCES)
  endif
endif

#----------------------------------------------------------------------------
#
# Preprocess all source files.  Implicit rules should take care of all cases.
#
#----------------------------------------------------------------------------

.PHONY: preprocess

preprocess: $(SOURCES)

#----------------------------------------------------------------------------
#
# Implicit rules for preprocessing.
#
#----------------------------------------------------------------------------
 
# Cancel the implicit gmake rules for preprocessing

%.c : %.C
%.f90 : %.F90
%.f : %.F

# Preprocessing rules for Fortran (.F, F90) and C files

$(POPEXEDIR)/compile/%.f: %.F
	@echo '$(POPARCH) preprocessing ' $<
	@$(Cpp) $(Cpp_opts) $< > $(POPEXEDIR)/compile/$*.f

$(POPEXEDIR)/compile/%.f90: %.F90
	@echo '$(POPARCH) preprocessing ' $<
	@$(Cpp) $(Cpp_opts) $< > $(POPEXEDIR)/compile/$*.f90

#  For some reason, our current Cpp options are incorrect for C files
#  so let the C compiler take care of ifdefs and just copy.

$(POPEXEDIR)/compile/%.c: %.C
	@echo '$(POPARCH) preprocessing ' $<
	@$(Cp) $< $(POPEXEDIR)/compile/$*.c

# Preprocessing rules for Fortran f, f90 and c files
# Should only copy these files into the compile directory.

$(POPEXEDIR)/compile/%.f: %.f
	@echo '$(POPARCH) preprocessing ' $<
	@$(Cp) $< $(POPEXEDIR)/compile/$*.f

$(POPEXEDIR)/compile/%.f90: %.f90
	@echo '$(POPARCH) preprocessing ' $<
	@$(Cp) $< $(POPEXEDIR)/compile/$*.f90

$(POPEXEDIR)/compile/%.c: %.c
	@echo '$(POPARCH) preprocessing ' $<
	@$(Cp) $< $(POPEXEDIR)/compile/$*.c

#----------------------------------------------------------------------------
