#-------------------------------------------------------------------------------
# SVN $Id: Makefile 49095 2013-07-16 21:51:37Z mlevy@ucar.edu $
# SVN $Source$
#-------------------------------------------------------------------------------
# Common Makefile: a framework for building all CCSM components and more
#
# Command-line variables
#   MACFILE=<file> ~ the macros definition file to use/include
#   EXEC=<name>    ~ name given to executable, default is a.out
#   VPATH=<vpath>  ~ VPATH            , default is . (cwd only)
#   SRCS=<files>   ~ list of src files, default is all .c .F .F90 files in VPATH
#   VPFILE=<file>  ~ file with list of dirs, used to create VPATH
#   SRCFILE=<file> ~ file with list of src files, used to create SRCS
#   DEPGEN=<exec>  ~ dependency generator utility, default is makdep
#
#   <macro defns>  ~ any macro definitions found in this file or the included
#                    MACFILE will be over-riden by cmd-line macro definitions
#   MODEL=<model>  ~ a standard macro definition, often found in the included
#                    MACFILE, used to trigger special compilation flags
#
# Usage examples:
#   % gmake MACFILE=Macros.AIX VPFILE=Filepath MODEL=ccm3 EXEC=atm
#   % gmake MACFILE=Macros.AIX VPFILE=Filepath SRCFILE=Srclist EXEC=pop
#   % gmake MACFILE=Macros.C90 VPATH="dir1 dir2" SRCS="file1.c file2.F90"
#   % gmake MACFILE=Macros.SUN SRCS="test.F"
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# parse cmd-line and establish values for EXEC, VPATH, SRCS, OBJS, etc
#-------------------------------------------------------------------------------

EXEC     := ../runoff_map
MACFILE  := Macros.make
MODEL    := NONE
VPFILE   := NONE
VPATH    := $(shell pwd)
SRCFILE  := NONE
SRCS     := NONE
DEPGEN   := ./makdep  # an externally provided dependency generator
#compile_threaded := TRUE
compile_threaded := FALSE

ifneq ($(VPATH),.)
  # this variable was specified on cmd line or in an env var
else
  ifneq ($(VPFILE),NONE)
    # explicit list of VPATH dirs is provided
    VPATH := $(wildcard . $(shell cat $(VPFILE) ) )
  endif
endif

ifneq ($(SRCS),NONE)
  # this variable was specified on cmd line or in an env var
else
  ifneq ($(SRCFILE),NONE)
    # explicit list of src files is provided
    SRCS := $(shell cat $(SRCFILE) )
  else
    # list of src files is all .F90 .F .c files in VPATH
    SRCS := $(wildcard $(addsuffix /*.F90 , $(VPATH)) \
		       $(addsuffix /*.[cF], $(VPATH)) )
  endif
endif

OBJS  := $(addsuffix .o, $(sort $(basename $(notdir $(SRCS)))))
DEPS  := $(addsuffix .d, $(sort $(basename $(notdir $(SRCS)))))
INCS  := $(patsubst %,-I%, $(VPATH) )
RM    := rm

.SUFFIXES:
.SUFFIXES: .F90 .F .c .o

all: $(EXEC)

#-------------------------------------------------------------------------------
# include the file that provides macro definitions required by build rules
# note: the MACFILE may not be needed for certain goals
#-------------------------------------------------------------------------------

ifneq ($(MAKECMDGOALS), db_files)
  -include $(MACFILE)
endif
FFLAGS   += -c
CFLAGS   += -c
INCLDIR  += -I$(NETCDF_PATH)/include

#-------------------------------------------------------------------------------
# echo file names, paths, compile flags, etc. used during build
#-------------------------------------------------------------------------------

.PHONY: db_files
db_files:
	@echo " "
	@echo "* EXEC    := $(EXEC)"
	@echo "* MACFILE := $(MACFILE)"
	@echo "* VPFILE  := $(VPFILE)"
	@echo "* VPATH   := $(VPATH)"
	@echo "* SRCFILE := $(SRCFILE)"
	@echo "* INCS    := $(INCS)"
	@echo "* SRCS    := $(SRCS)"
	@echo "* OBJS    := $(OBJS)"
	@echo "* DEPS    := $(DEPS)"

.PHONY: db_flags
db_flags:
	@echo " "
	@echo "* cc      := $(SCC)  $(CFLAGS) $(INCS) $(INCLDIR)"
	@echo "* .F.o    := $(SFC)  $(FFLAGS) $(FIXEDFLAGS) $(CPPDEFS) $(INCS) $(INCLDIR)"
	@echo "* .F90.o  := $(SFC)  $(FFLAGS) $(FREEFLAGS) $(CPPDEFS) $(INCS) $(INCLDIR)"

#-------------------------------------------------------------------------------
# build rules: MACFILE, cmd-line, or env vars must provide the needed macros
#-------------------------------------------------------------------------------

$(EXEC): $(OBJS)
	$(SFC) -o $(EXEC) $(OBJS) $(SLIBS) $(ULIBS) $(LDFLAGS)

.c.o:
	$(SCC) $(INCLDIR) $(INCS) $(CPPDEFS) $(CFLAGS)  $<

ifeq ($(CPP),NONE)
.F.o:
	$(SFC) $(FFLAGS) $(FIXEDFLAGS) $(CPPFLAGS) $(CPPDEFS) $(INCLDIR) $(INCS)  $<

.F90.o:
	$(SFC) $(FFLAGS) $(FREEFLAGS)  $(CPPFLAGS) $(CPPDEFS) $(INCLDIR) $(INCS)  $<
else
.F.o:
	$(SFC)  $(FFLAGS) $(FIXEDFLAGS) $(CPPDEFS) $(INCS) $(INCLDIR) $*.f

.F90.o:
	$(SFC)  $(FFLAGS) $(FREEFLAGS) $(CPPDEFS) $(INCS) $(INCLDIR) $*.F90
endif

.PHONY: clean
clean:
	$(RM) -f *.f *.f90 *.d *.mod *.optrpt $(OBJS) $(EXEC) $(DEPGEN)

.PHONY: distclean
distclean: clean
	../../../../configure --clean

#-------------------------------------------------------------------------------
# Build & include dependency files
#-------------------------------------------------------------------------------
# ASSUMPTIONS:
# o an externally provided dependency generator, $(DEPGEN), is available,
#   its cmd line syntax is compatible with the build rules below.   Eg, for
#   each .o file, there is a corresponding .d (dependency) file, and both
#   will be dependent on the same src file, eg.    foo.o foo.d : foo.F90
#   Also, the dependancy genorator's capabilities, limitations, and assumptions
#   are understood & accepted.
#-------------------------------------------------------------------------------

%.d : %.c $(DEPGEN)
	@ echo "Building dependency for $@"
	@ $(DEPGEN) -f $(INCS) $< | head -3  > $@
%.d : %.F $(DEPGEN)
	@ echo "Building dependency for $@"
	@ $(DEPGEN) -f $(INCS) $<  > $@
%.d : %.F90 $(DEPGEN)
	@ echo "Building dependency for $@"
	@ $(DEPGEN) -f $(INCS) $<  > $@
%.d : %.H $(DEPGEN)
	@ echo "Building dependency for $@"
	@ $(DEPGEN) -f $(INCS) $<  > $@

$(DEPGEN) :
	$(SCC) -o $(DEPGEN) ../tools/makdep.c

# the if-tests prevent DEPS files from being created when they're not needed
ifneq ($(MAKECMDGOALS), db_files)
ifneq ($(MAKECMDGOALS), db_flags)
ifneq ($(MAKECMDGOALS), clean)
ifneq ($(MAKECMDGOALS), distclean)
    -include $(DEPS)
endif
endif
endif
endif
