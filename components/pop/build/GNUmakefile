#
# File:  GNUmakefile
#
#----------------------------------------------------------------------------
#
#  This is the main driver makefile for compiling POP.  It relies on
#  sub-makefiles to perform preprocessing, dependency analysis and compiling.
#
#  Several variables must be set either as environment variables or on the 
#  command line (e.g. gmake POPDIR=/your/path).  These variables are:
#
#  POPEXEDIR = the directory where you want the executable to end up.  A
#              common value might be . (the current working directory)
#  POPDIR    = the path to the POP distribution from which you want
#              to grab source files, etc.  This should be the path to the 
#              top level pop directory, not the individual source,mpi 
#              directories.
#  POPARCH   = the base name of a file in the build directory where 
#              architecture and location-specific information is defined
#
#  The optimize variable has a default value of yes and can not be set
#  as an environment variable.  It can only be changed to no via the
#  command line.
#  OPTIMIZE  = yes/no  whether you want full optimization or no optimization
#
#----------------------------------------------------------------------------
#
#  The bourne shell is safest for gmake
SHELL = /bin/sh
#
#
#  By default, you build optimized code.  To override, add "OPTIMIZE=no" 
#  to the gmake command line.
OPTIMIZE   = yes
#
#----------------------------------------------------------------------------
#
#  Check that path variables POPEXEDIR, POPDIR have been set.
#
#----------------------------------------------------------------------------

VARSDEFINED =

ifneq (,$(POPEXEDIR))
  VARSDEFINED := $(VARSDEFINED)POPEXEDIR
  export POPEXEDIR
else
  bogusexedir:
	@echo "  Please setenv POPEXEDIR"
endif

ifneq (,$(POPDIR))
  VARSDEFINED := $(VARSDEFINED)POPDIR
  export POPDIR
else
  boguspopdir:
	@echo "  Please setenv POPDIR"
endif

#  Set the directory where dependency info will reside.
DepDir = $(POPEXEDIR)/compile/Depends

#----------------------------------------------------------------------------
#
#  Include architecture-specific flags and options. 
#
#----------------------------------------------------------------------------

ifneq (,$(POPARCH))
  export POPARCH
  VARSDEFINED := $(VARSDEFINED)POPARCH
  ifneq (,$(POPDIR))
     include $(POPDIR)/build/$(POPARCH).gnu
  endif
else
  boguspoparch:
	@echo "  Please setenv POPARCH"
endif

#----------------------------------------------------------------------------
#
#  Define a symbol (TARGETX) for the executable name (pop)
#
#----------------------------------------------------------------------------

TARGETX = pop
ifeq ($(OPTIMIZE),no)
  #  If building a debug version, append "_db" to the executable name
  TARGETX := $(TARGETX)_db
endif

#----------------------------------------------------------------------------
#
#  Make the pop executable
#
#----------------------------------------------------------------------------

.PHONY: depends preprocess compile dirs

$(TARGETX): depends
	@echo "  GNUmakefile is making target '$(TARGETX)'"
	@gmake -f $(POPDIR)/build/compile.mk POPEXEDIR=$(POPEXEDIR) POPDIR=$(POPDIR) OPTIMIZE=$(OPTIMIZE) 

#----------------------------------------------------------------------------
#
#  Create the file dependencies for compiling.
#
#----------------------------------------------------------------------------

depends: preprocess
	@echo "$(POPEXEDIR) $(POPDIR)"
	@echo "$(POPARCH) automatically generating dependencies for compiling"
	@gmake -f $(POPDIR)/build/depends.mk POPEXEDIR=$(POPEXEDIR) OPTIMIZE=$(OPTIMIZE)

#----------------------------------------------------------------------------
#
#  Perform the preprocessing
#
#----------------------------------------------------------------------------

preprocess: dirs
	@echo "$(POPARCH) preprocess POP source"
	@gmake -f $(POPDIR)/build/preprocess.mk

#----------------------------------------------------------------------------
#
#  Create the directory structure if it does not exist.
#
#----------------------------------------------------------------------------

dirs:
	@echo "$(POPARCH) is creating the required directory structure"
	@cd $(POPEXEDIR) &&  if [ ! -d compile ]; then mkdir compile; fi
	@cd $(POPEXEDIR)/compile &&  \
            if [ ! -d Depends ]; then mkdir Depends; fi 
	@cd $(POPEXEDIR)

#----------------------------------------------------------------------------
#
#  Include utilities makefile with rules for clean, clobber, etc.
#
#----------------------------------------------------------------------------

clean:
ifeq ($(VARSDEFINED),POPEXEDIRPOPDIRPOPARCH)
	@cd $(POPEXEDIR) && $(RM) core pop pop_db 
	@cd $(POPEXEDIR) && $(RM) compile/Depends/*.* 
	@cd $(POPEXEDIR) && $(RM) compile/*.f compile/*.f90 compile/*.c
	@cd $(POPEXEDIR) && $(RM) compile/*.o compile/*.$(MODSUF) 
else
	@echo "  Please setenv POPEXEDIR POPDIR and POPARCH"
	@echo "  Only $(VARSDEFINED) have been defined."
endif

#----------------------------------------------------------------------------
