#-----------------------------------------------------------------------
# This Makefile is for building clm tools on AIX, Linux (with pgf90 or 
# lf95 compiler), Darwin or IRIX platforms.
#
# $Id: Makefile 35671 2012-03-22 21:23:14Z mlevy@ucar.edu $
# $URL$
#
# These macros can be changed by setting environment variables:
#
# NETCDF_PATH --- Diretory location of netcdf
# LIB_NETCDF --- Library directory location of netcdf. (defaults to /usr/local/lib)
# INC_NETCDF --- Include directory location of netcdf. (defaults to /usr/local/include)
# MOD_NETCDF --- Module directory location of netcdf.  (defaults to $LIB_NETCDF)
# USER_FC ------ Allow user to override the default Fortran compiler specified in Makefile.
# USER_FCTYP --- Allow user to override the default type of Fortran compiler (linux and USER_FC=ftn only).
# USER_CC ------ Allow user to override the default C compiler specified in Makefile (linux only).
# USER_LINKER -- Allow user to override the default linker specified in Makefile.
# USER_CPPDEFS - Additional CPP defines.
# USER_CFLAGS -- Additional C compiler flags that the user wishes to set.
# USER_FFLAGS -- Additional Fortran compiler flags that the user wishes to set.
# USER_LDLAGS -- Additional load flags that the user wishes to set.
# SMP ---------- Shared memory Multi-processing (TRUE or FALSE) [default is FALSE]
# OPT ---------- Use optimized options.
#
#------------------------------------------------------------------------

# Definitions
null :=
CURDIR = .
EXENAME = ../gen_domain
RM = rm

default: $(EXENAME)

Macros:
	@echo "use the configure script located in the Machines directory to create the Makefile Macros file"

-include Macros
# Check for the netcdf library and include directories 
ifdef NETCDF_PATH
  LIB_NETCDF:=$(NETCDF_PATH)/lib
  INC_NETCDF:=$(NETCDF_PATH)/include
  MOD_NETCDF:=$(NETCDF_PATH)/include
else
  ifeq ($(LIB_NETCDF),$(null))
    LIB_NETCDF := /usr/local/lib
  endif
  ifeq ($(INC_NETCDF),$(null))
    INC_NETCDF := /usr/local/include
  endif
 
  ifeq ($(MOD_NETCDF),$(null))
    MOD_NETCDF := $(LIB_NETCDF)
  endif
endif

# Set Fortran compiler
ifneq ($(SFC),$(null))
  FC := $(SFC)
endif

# Set user specified C compiler
ifneq ($(SCC),$(null))
  CC := $(SCC)
endif

# Set if Shared memory multi-processing will be used
ifeq ($(SMP),$(null))
  SMP := FALSE
endif

CPPDEF += $(USER_CPPDEFS)

# Set optimization on by default
ifeq ($(OPT),$(null))
  OPT := TRUE
endif

ifeq ($(OPT),TRUE)
  CPPDEF := -DOPT
endif

# Load dependency search path.
dirs := . $(shell cat Filepath)

# Set cpp search path, include netcdf
cpp_dirs := $(dirs) $(INC_NETCDF) $(MOD_NETCDF)
cpp_path := $(foreach dir,$(cpp_dirs),-I$(dir)) # format for command line

# Expand any tildes in directory names. Change spaces to colons.
VPATH := $(foreach dir,$(cpp_dirs),$(wildcard $(dir))) 
VPATH := $(subst $(space),:,$(VPATH))               
ifdef CPRE
  FPPDEFS := $(patsubst -D%,$(CPRE)%,$(CPPDEFS)) 
else
  FPPDEFS := $(CPPDEFS)
endif
#Primary Target: build the tool
all: $(EXENAME) Macros

OBJS := gen_domain.o 

#------------------------------------------------------------------------
# Default rules and macros
#------------------------------------------------------------------------

.SUFFIXES:
.SUFFIXES: .F90 .c .o

# Append user defined compiler and load flags to Makefile defaults
CFLAGS   += $(USER_CFLAGS) -I$(INC_NETCDF) 
FFLAGS   += $(USER_FFLAGS) -I$(MOD_NETCDF) -I$(INC_NETCDF) 
LDFLAGS  += $(USER_LDFLAGS)

# Set user specified linker
ifneq ($(USER_LINKER),$(null))
  LINKER := $(USER_LINKER)
else
  LINKER := $(FC)
endif

.F90.o:
	$(FC) -c $(FPPDEFS) $(FFLAGS) $<

.c.o:
	$(CC) -c $(CFLAGS) $<


$(EXENAME): $(OBJS)
	$(LINKER) -o $@ $(OBJS) $(SLIBS) $(LDFLAGS)

clean:
	$(RM) -f $(OBJS) $(EXENAME) *.mod 

veryclean:
	$(RM) -f $(OBJS) $(EXENAME) *.mod Macros env_mach* build*

gen_domain.o : gen_domain.F90
