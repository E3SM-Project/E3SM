#===============================================================================
#
# Common Makefile: a framework for building all CIME components and more
#
#===============================================================================

# Set up special characters
null  :=
comma := ,

# Load dependency search path.
dirs := .
dirs += $(shell cat Filepath)

cpp_dirs := $(dirs)
# Add INCROOT to path for Depends and Include
MINCROOT :=
ifdef INCROOT
  cpp_dirs += $(INCROOT)
  MINCROOT := $(INCROOT)
endif

# Expand any tildes in directory names. Change spaces to colons.
VPATH := $(foreach dir,$(cpp_dirs),$(wildcard $(dir)))
VPATH := $(subst $(space),:,$(VPATH))

RM    := rm
CP    := cp

exec_se: $(EXEC_SE)  $(CURDIR)/Depends
complib: $(COMPLIB)  $(CURDIR)/Depends

# Determine whether to compile threaded or not
# Set the THREADDIR for the shared build
# based on the threaded build status
compile_threaded = false
ifeq ($(strip $(SMP)),TRUE)
   compile_threaded = true
   THREADDIR = threads
else
   ifeq ($(strip $(BUILD_THREADED)),TRUE)
      compile_threaded = true
      THREADDIR = threads
   else
      THREADDIR = nothreads
   endif
endif

# set the debug directory based on the debug status
ifeq ($(strip $(DEBUG)),TRUE)
   DEBUGDIR = debug
else
   DEBUGDIR = nodebug
endif

ifeq ($(strip $(USE_ESMF_LIB)), TRUE)
   ESMFDIR = esmf
else
   ESMFDIR = noesmf
endif

# Determine whether any C++ code will be included in the build;
# currently, C++ code is included if and only if we're linking to the
# trilinos library or the Albany library.
USE_CXX = false
ifeq ($(strip $(USE_TRILINOS)), TRUE)
   USE_CXX = true
endif
ifeq ($(strip $(USE_ALBANY)), TRUE)
   USE_CXX = true
endif

ifndef MOD_SUFFIX
   MOD_SUFFIX := mod
endif

#===============================================================================
# set CPP options (must use this before any flags or cflags settings)
#===============================================================================

CPPDEFS := $(USER_CPPDEFS)

# Unless DEBUG mode is enabled, use NDEBUG to turn off assert statements.
ifneq ($(strip $(DEBUG)),TRUE)
   CPPDEFS += -DNDEBUG
endif

# USE_ESMF_LIB is currently only defined in env_build.xml
ifeq ($(USE_ESMF_LIB), TRUE)
   CPPDEFS += -DUSE_ESMF_LIB
endif

# ESMF_INTERFACE is currently only defined in env_build.xml
ifeq ($(COMP_INTERFACE), ESMF)
   CPPDEFS += -DESMF_INTERFACE
else
   CPPDEFS += -DMCT_INTERFACE
endif

ifeq ($(strip $(MPILIB)),mpi-serial)
  CPPDEFS += -DNO_MPI2
else
  CPPDEFS += -DHAVE_MPI
endif
ifeq ($(compile_threaded), true)
  CPPDEFS += -DTHREADED_OMP
endif

ifeq (,$(EXEROOT))
  EXEROOT = $(shell ./xmlquery EXEROOT -value)
endif
ifeq (,$(BUILD_THREADED))
  BUILD_THREADED = $(shell ./xmlquery BUILD_THREADED -value)
endif

ifeq (,$(LIBROOT))
  LIBROOT = $(shell ./xmlquery LIBROOT -value)
endif
ifeq (,$(SHAREDLIBROOT))
  SHAREDLIBROOT = $(shell ./xmlquery SHAREDLIBROOT -value)
endif
ifeq (,$(COMPILER))
  COMPILER = $(shell ./xmlquery COMPILER -value)
endif
ifeq (,$(NINST_VALUE))
  NINST_VALUE = $(shell ./xmlquery NINST_VALUE -value)
endif
ifeq (,$(MPILIB))
  MPILIB = $(shell ./xmlquery MPILIB -value)
endif
ifeq ($(strip $(PIO_VERSION)),1)
  CPPDEFS += -DPIO1
endif

ifeq (,$(SHAREDPATH))
  SHAREDPATH = $(COMPILER)/$(MPILIB)/$(DEBUGDIR)/$(THREADDIR)
  INSTALL_SHAREDPATH = $(EXEROOT)/$(SHAREDPATH)
endif

include $(CASEROOT)/Macros.make
# Decide whether to use a C++ or Fortran linker, based on whether we
# are using any C++ code and the compiler-dependent CXX_LINKER variable
ifeq ($(USE_CXX), true)
  # The following is essentially an "if... elseif... else", but gmake
  # 3.80 and earlier doesn't support elseif
  ifeq ($(CXX_LINKER), CXX)
    LD := $(MPICXX)
  endif
  ifeq ($(CXX_LINKER), FORTRAN)
    LD := $(MPIFC)
  endif
endif
# Use this if LD has not already been defined.
ifeq ($(origin LD), default)
  ifeq ($(strip $(MPILIB)),mpi-serial)
    LD := $(SFC)
  else
    LD := $(MPIFC)
  endif
endif



ifeq ($(USE_CXX), true)
  ifeq ($(SUPPORTS_CXX), FALSE)
    $(error Fatal attempt to include C++ code on a compiler/machine combo that has not been set up to support C++)
  endif
endif

# Not clear how to escape commas for libraries with their own configure
# script, and they don't need this defined anyway, so leave this out of
# FPPDEFS.
ifeq ($(HAS_F2008_CONTIGUOUS),TRUE)
  CONTIGUOUS_FLAG := -DUSE_CONTIGUOUS=contiguous,
else
  CONTIGUOUS_FLAG := -DUSE_CONTIGUOUS=
endif

ifdef CPRE
  CONTIGUOUS_FLAG := $(subst $(comma),\\$(comma),$(CONTIGUOUS_FLAG))
  CONTIGUOUS_FLAG := $(patsubst -D%,$(CPRE)%,$(CONTIGUOUS_FLAG))
endif

ifndef AR
   AR := ar
endif

ifdef NETCDF_C_PATH
  ifndef NETCDF_FORTRAN_PATH
    $(error "NETCDF_C_PATH specified without NETCDF_FORTRAN_PATH")
  endif
  NETCDF_SEPARATE:=true
  ifndef INC_NETCDF_C
    INC_NETCDF_C:=$(NETCDF_C_PATH)/include
  endif
  ifndef INC_NETCDF_FORTRAN
    INC_NETCDF_FORTRAN:=$(NETCDF_FORTRAN_PATH)/include
  endif
  ifndef LIB_NETCDF_C
    LIB_NETCDF_C:=$(NETCDF_C_PATH)/lib
  endif
  ifndef LIB_NETCDF_FORTRAN
    LIB_NETCDF_FORTRAN:=$(NETCDF_C_PATH)/lib
  endif
 else ifdef NETCDF_FORTRAN_PATH
  $(error "NETCDF_FORTRAN_PATH specified without NETCDF_C_PATH")
 else ifdef NETCDF_PATH
  NETCDF_SEPARATE:=false
  ifndef INC_NETCDF
    INC_NETCDF:=$(NETCDF_PATH)/include
  endif
  ifndef LIB_NETCDF
    LIB_NETCDF:=$(NETCDF_PATH)/lib
  endif
else
  # No Netcdf is an error unless target is clean or DEP
  ifneq ($(MAKECMDGOALS), db_files)
  ifneq ($(MAKECMDGOALS), db_flags)
  ifeq (,$(findstring clean,$(MAKECMDGOALS)))
      $(error NETCDF not found: Define NETCDF_PATH or NETCDF_C_PATH and NETCDF_FORTRAN_PATH in config_machines.xml or config_compilers.xml)
  endif
  endif
  endif
endif


ifeq ($(MPILIB),mpi-serial)
  undefine PNETCDF_PATH
else
  ifdef PNETCDF_PATH
    ifndef $(INC_PNETCDF)
      INC_PNETCDF:=$(PNETCDF_PATH)/include
    endif
    ifndef LIB_PNETCDF
      LIB_PNETCDF:=$(PNETCDF_PATH)/lib
    endif
  endif
endif
# Set PETSc info if it is being used
ifeq ($(strip $(USE_PETSC)), TRUE)
  ifdef PETSC_PATH
    ifndef INC_PETSC
      INC_PETSC:=$(PETSC_PATH)/include
    endif
    ifndef LIB_PETSC
      LIB_PETSC:=$(PETSC_PATH)/lib
    endif
  else
    $(error PETSC_PATH must be defined when USE_PETSC is TRUE)
  endif

  # Get the "PETSC_LIB" list an env var
  include $(PETSC_PATH)/conf/variables
endif

# Set Trilinos info if it is being used
ifeq ($(strip $(USE_TRILINOS)), TRUE)
  ifdef TRILINOS_PATH
    ifndef INC_TRILINOS
      INC_TRILINOS:=$(TRILINOS_PATH)/include
    endif
    ifndef LIB_TRILINOS
      LIB_TRILINOS:=$(TRILINOS_PATH)/lib
    endif
  else
    $(error TRILINOS_PATH must be defined when USE_TRILINOS is TRUE)
  endif

  # get a bunch of variables related to this trilinos installation;
  # these variables begin with "Trilinos_"
  include $(INC_TRILINOS)/Makefile.export.Trilinos
endif

# Set Albany info if it is being used
ifeq ($(strip $(USE_ALBANY)), TRUE)
  ifdef ALBANY_PATH
    ifndef INC_ALBANY
      INC_ALBANY:=$(ALBANY_PATH)/include
    endif
    ifndef LIB_ALBANY
      LIB_ALBANY:=$(ALBANY_PATH)/lib
    endif
  else
    $(error ALBANY_PATH must be defined when USE_ALBANY is TRUE)
  endif

  # get the "ALBANY_LINK_LIBS" list as an env var
  include $(ALBANY_PATH)/export_albany.in
endif

# Set HAVE_SLASHPROC on LINUX systems which are not bluegene or Darwin (OSx)

ifeq ($(findstring -DLINUX,$(CPPDEFS)),-DLINUX)
  ifneq ($(findstring DBG,$(CPPDEFS)),DBG)
    ifneq ($(findstring Darwin,$(CPPDEFS)),Darwin)
      CPPDEFS += -DHAVE_SLASHPROC
    endif
  endif
endif

ifdef CPRE
  FPPDEFS := $(subst $(comma),\\$(comma),$(CPPDEFS))
  FPPDEFS := $(patsubst -D%,$(CPRE)%,$(FPPDEFS))
else
  FPPDEFS := $(CPPDEFS)
endif

#===============================================================================
# Set config args for pio and mct to blank and then enable serial
#===============================================================================
ifndef CONFIG_ARGS
  CONFIG_ARGS :=
endif
ifeq  ($(findstring pio,$(MODEL)),pio)
  CONFIG_ARGS+= --enable-timing
  ifeq ($DEBUG,TRUE)
     CONFIG_ARGS+= --enable-debug
  endif
endif

#===============================================================================
# User-specified INCLDIR
#===============================================================================

INCLDIR := -I.
ifdef USER_INCLDIR
  INCLDIR += $(USER_INCLDIR)
endif

#===============================================================================
# MPI-serial library (part of MCT)
#===============================================================================

ifeq ($(strip $(MPILIB)), mpi-serial)
  CC      := $(SCC)
  FC      := $(SFC)
  CXX     := $(SCXX)
  MPIFC   := $(SFC)
  MPICC   := $(SCC)
  MPICXX  := $(SCXX)
  CONFIG_ARGS += MCT_PATH=$(SHAREDLIBROOT)/$(SHAREDPATH)/mct/mpi-serial
else
  CC  := $(MPICC)
  FC  := $(MPIFC)
  CXX := $(MPICXX)
  ifdef MPI_PATH
    INC_MPI := $(MPI_PATH)/include
    LIB_MPI := $(MPI_PATH)/lib
  endif
endif

#===============================================================================
# Set include paths (needed after override for any model specific builds below)
#===============================================================================
INCLDIR += -I$(INSTALL_SHAREDPATH)/include -I$(INSTALL_SHAREDPATH)/$(COMP_INTERFACE)/$(ESMFDIR)/$(NINST_VALUE)/include

ifeq ($(NETCDF_SEPARATE), false)
  INCLDIR += -I$(INC_NETCDF)
else ifeq ($(NETCDF_SEPARATE), true)
  INCLDIR += -I$(INC_NETCDF_C) -I$(INC_NETCDF_FORTRAN)
endif
ifdef MOD_NETCDF
  INCLDIR += -I$(MOD_NETCDF)
endif
ifdef INC_MPI
  INCLDIR += -I$(INC_MPI)
endif
ifdef INC_PNETCDF
  INCLDIR += -I$(INC_PNETCDF)
endif
ifdef INC_PETSC
  INCLDIR += -I$(INC_PETSC)
endif
ifdef INC_TRILINOS
  INCLDIR += -I$(INC_TRILINOS)
endif
ifdef INC_ALBANY
  INCLDIR += -I$(INC_ALBANY)
endif

ifeq ($(MODEL),driver)
  INCLDIR += -I$(EXEROOT)/atm/obj -I$(EXEROOT)/ice/obj -I$(EXEROOT)/ocn/obj -I$(EXEROOT)/glc/obj -I$(EXEROOT)/rof/obj -I$(EXEROOT)/wav/obj -I$(EXEROOT)/esp/obj
# nagfor and gcc have incompatible LDFLAGS.
# nagfor requires the weird "-Wl,-Wl,," syntax.
# If done in config_compilers.xml, we break MCT.
  ifeq ($(strip $(COMPILER)),nag)
     ifeq ($(NETCDF_SEPARATE), false)
       SLIBS += -Wl,-Wl,,-rpath=$(LIB_NETCDF)
     else ifeq ($(NETCDF_SEPARATE), true)
       SLIBS += -Wl,-Wl,,-rpath=$(LIB_NETCDF_C)
       SLIBS += -Wl,-Wl,,-rpath=$(LIB_NETCDF_FORTRAN)
     endif
  endif
else
  ifeq ($(strip $(COMPILER)),nag)
    ifeq ($(DEBUG), TRUE)
      ifeq ($(strip $(MACH)),hobart)
       # GCC needs to be able to link to
       # nagfor runtime to get autoconf
       # tests to work.
         CFLAGS += -Wl,--as-needed,--allow-shlib-undefined
         SLIBS += -L$(COMPILER_PATH)/lib/NAG_Fortran -lf60rts
       endif
    endif
  endif
endif


ifndef MCT_LIBDIR
  MCT_LIBDIR=$(INSTALL_SHAREDPATH)/lib
endif

ifndef PIO_LIBDIR
  PIO_LIBDIR=$(INSTALL_SHAREDPATH)/lib
endif

ifndef GPTL_LIBDIR
  GPTL_LIBDIR=$(INSTALL_SHAREDPATH)/lib
endif

ifndef GLC_DIR
  GLC_DIR=$(EXEROOT)/glc
endif
ifndef CISM_LIBDIR
  CISM_LIBDIR=$(GLC_DIR)/lib
endif

INCLDIR +=	-I$(INSTALL_SHAREDPATH)/include -I$(CIMEROOT)/share/csm_share/shr \
		-I$(CIMEROOT)/share/csm_share/include -I$(CIMEROOT)/share/shr_RandNum/include
#
# Use the MCT dir for the cache for all configure calls because it is the first one
#
CFLAGS+=$(CPPDEFS)
CXXFLAGS := $(CFLAGS)

CONFIG_ARGS +=  CC="$(CC)" FC="$(FC)" MPICC="$(MPICC)" \
                MPIFC="$(MPIFC)" FCFLAGS="$(FFLAGS) $(FREEFLAGS) $(INCLDIR)" \
                CPPDEFS="$(CPPDEFS)" CFLAGS="$(CFLAGS) -I.. $(INCLDIR)" \
                LDFLAGS="$(LDFLAGS)"

ifeq ($(NETCDF_SEPARATE), false)
  CONFIG_ARGS += NETCDF_PATH=$(NETCDF_PATH)
else ifeq ($(NETCDF_SEPARATE), true)
  # The mct library needs the NetCDF_C library
  CONFIG_ARGS += NETCDF_PATH=$(NETCDF_C_PATH)
endif

ifeq ($(COMPILER),nag)
  CONFIG_ARGS += LIBS="$(SLIBS)"
endif

FFLAGS += $(FPPDEFS)
FFLAGS_NOOPT += $(FPPDEFS)


ifeq ($(findstring -cosp,$(CAM_CONFIG_OPTS)),-cosp)
# The following is for the COSP simulator code:
COSP_LIBDIR:=$(EXEROOT)/atm/obj/cosp
endif

ifeq ($(MODEL),cam)
  # These RRTMG files take an extraordinarily long time to compile with optimization.
  # Until mods are made to read the data from files, just remove optimization from
  # their compilation.
rrtmg_lw_k_g.o: rrtmg_lw_k_g.f90
	$(FC) -c $(FPPFLAGS) $(INCLDIR) $(INCS) $(FREEFLAGS) $(FFLAGS_NOOPT) $<

rrtmg_sw_k_g.o: rrtmg_sw_k_g.f90
	$(FC) -c $(FPPFLAGS) $(INCLDIR) $(INCS) $(FREEFLAGS) $(FFLAGS_NOOPT) $<


  ifdef COSP_LIBDIR
    INCLDIR+=-I$(COSP_LIBDIR) -I$(COSP_LIBDIR)/../
    $(COSP_LIBDIR)/libcosp.a: cam_abortutils.o
    $(MAKE) -C $(COSP_LIBDIR) F90='$(FC)' F90FLAGS='$(INCLDIR) $(INCS) $(FREEFLAGS) $(FFLAGS) $(FC_AUTO_R8)' \
    F90FLAGS_noauto='$(INCLDIR) $(INCS) $(FREEFLAGS) $(FFLAGS)' \
    F90FLAGS_fixed='$(INCLDIR) $(INCS) $(FIXEDFLAGS) $(FFLAGS) $(FC_AUTO_R8)'

    cospsimulator_intr.o: $(COSP_LIBDIR)/libcosp.a
  endif

endif

# Set esmf.mk location with ESMF_LIBDIR having precedent over ESMFMKFILE
CCSM_ESMFMKFILE := undefined_CCSM_ESMFMKFILE
ifdef ESMFMKFILE
   CCSM_ESMFMKFILE := $(ESMFMKFILE)
endif
ifdef ESMF_LIBDIR
   CCSM_ESMFMKFILE := $(ESMF_LIBDIR)/esmf.mk
endif


# System libraries (netcdf, mpi, pnetcdf, esmf, trilinos, etc.)
ifndef SLIBS
  ifeq ($(NETCDF_SEPARATE), false)
    SLIBS := -L$(LIB_NETCDF) -lnetcdff -lnetcdf
  else ifeq ($(NETCDF_SEPARATE), true)
    SLIBS := -L$(LIB_NETCDF_FORTRAN) -L$(LIB_NETCDF_C) -lnetcdff -lnetcdf
  endif
endif

ifdef LIB_PNETCDF
   SLIBS += -L$(LIB_PNETCDF) -lpnetcdf
endif
ifdef LAPACK_LIBDIR
   SLIBS += -L$(LAPACK_LIBDIR) -llapack -lblas
endif
ifdef LIB_MPI
   ifndef MPI_LIB_NAME
      SLIBS += -L$(LIB_MPI) -lmpi
   else
      SLIBS += -L$(LIB_MPI) -l$(MPI_LIB_NAME)
   endif
endif

# For compiling and linking with external ESMF.
# If linking to external ESMF library then include esmf.mk
# ESMF_F90COMPILEPATHS
# ESMF_F90LINKPATHS
# ESMF_F90LINKRPATHS
# ESMF_F90ESMFLINKLIBS
ifeq ($(USE_ESMF_LIB), TRUE)
  -include $(CCSM_ESMFMKFILE)
  FFLAGS += $(ESMF_F90COMPILEPATHS)
  SLIBS  += $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) $(ESMF_F90ESMFLINKLIBS)
endif

# Add PETSc libraries
ifeq ($(strip $(USE_PETSC)), TRUE)
  SLIBS += -L${LIB_PETSC} ${PETSC_LIB}
endif

# Add trilinos libraries; too be safe, we include all libraries included in the trilinos build,
# as well as all necessary third-party libraries
ifeq ($(strip $(USE_TRILINOS)), TRUE)
  SLIBS += -L$(LIB_TRILINOS) $(Trilinos_LIBRARIES) $(Trilinos_TPL_LIBRARY_DIRS) $(Trilinos_TPL_LIBRARIES)
endif

# Add Albany libraries.  These are defined in the ALBANY_LINK_LIBS env var that was included above
ifeq ($(strip $(USE_ALBANY)), TRUE)
  SLIBS += $(ALBANY_LINK_LIBS)
endif

# Add libraries and flags that we need on the link line when C++ code is included
# We need to do these additions after CONFIG_ARGS is set, because they can sometimes break configure for mct, etc.,
# if they are added to LDFLAGS in CONFIG_ARGS.
ifeq ($(USE_CXX), true)
  ifdef CXX_LIBS
    SLIBS += $(CXX_LIBS)
  endif

  ifdef CXX_LDFLAGS
    LDFLAGS += $(CXX_LDFLAGS)
  endif
endif

# Machine stuff to appear last on the link step
ifndef MLIBS
     MLIBS  :=
endif

#------------------------------------------------------------------------------
# Drive configure scripts for support libraries (mct)
#------------------------------------------------------------------------------

$(SHAREDLIBROOT)/$(SHAREDPATH)/mct/Makefile.conf:
	@echo "SHAREDLIBROOT |$(SHAREDLIBROOT)| SHAREDPATH |$(SHAREDPATH)|"; \
	$(CONFIG_SHELL) $(CIMEROOT)/externals/mct/configure $(CONFIG_ARGS) --srcdir $(CIMEROOT)/externals/mct

$(SHAREDLIBROOT)/$(SHAREDPATH)/mct/mpi-serial/Makefile.conf:
	@echo "SHAREDLIBROOT |$(SHAREDLIBROOT)| SHAREDPATH |$(SHAREDPATH)|"; \
	$(CONFIG_SHELL) $(CIMEROOT)/externals/mct/mpi-serial/configure $(CONFIG_ARGS) --srcdir $(CIMEROOT)/externals/mct/mpi-serial

ifeq ($(PIO_VERSION),2)
# This is a pio2 library
  PIOLIB = $(PIO_LIBDIR)/libpiof.a $(PIO_LIBDIR)/libpioc.a
  PIOLIBNAME = -lpiof -lpioc
  PIO_SRC_DIR = $(CIMEROOT)/externals/pio2
else
# This is a pio1 library
  PIOLIB = $(PIO_LIBDIR)/libpio.a
  PIOLIBNAME = -lpio
  ifneq ("$(wildcard $(CIMEROOT)/externals/pio1/pio)", "")
    PIO_SRC_DIR = $(CIMEROOT)/externals/pio1
  else
    PIO_SRC_DIR = $(CIMEROOT)/externals/pio1/pio
  endif
endif
#endif

MCTLIBS = $(MCT_LIBDIR)/libmct.a $(MCT_LIBDIR)/libmpeu.a

GPTLLIB = $(GPTL_LIBDIR)/libgptl.a

ULIBS += -L$(INSTALL_SHAREDPATH)/$(COMP_INTERFACE)/$(ESMFDIR)/$(NINST_VALUE)/lib -lcsm_share -L$(INSTALL_SHAREDPATH)/lib $(PIOLIBNAME) -lgptl -lmct -lmpeu

#------------------------------------------------------------------------------
# Drive cmake script for cism and pio
#------------------------------------------------------------------------------

ifndef CMAKE_OPTS
  CMAKE_OPTS :=
endif
# note that the fortran flags include neither the FREEFLAGS nor the
# FIXEDFLAGS, so that both free & fixed code can be built (cmake
# doesn't seem to be able to differentiate between free & fixed
# fortran flags)
CMAKE_OPTS += -D CMAKE_Fortran_FLAGS:STRING="$(FFLAGS) $(INCLDIR)" \
              -D CMAKE_C_FLAGS:STRING="$(CFLAGS) $(INCLDIR)" \
              -D CMAKE_CXX_FLAGS:STRING="$(CXXFLAGS) $(INCLDIR)" \
              -D CMAKE_VERBOSE_MAKEFILE:BOOL=ON \
              -D GPTL_PATH:STRING=$(INSTALL_SHAREDPATH) \
              -D PIO_ENABLE_TESTS:BOOL=OFF \
              -D USER_CMAKE_MODULE_PATH:STRING=$(CIMEROOT)/externals/CMake

# Allow for separate installations of the NetCDF C and Fortran libraries
ifeq ($(NETCDF_SEPARATE), false)
  CMAKE_OPTS += -D NETCDF_DIR:PATH=$(NETCDF_PATH)
else ifeq ($(NETCDF_SEPARATE), true)
  # NETCDF_Fortran_DIR points to the separate
  # installation of Fortran NetCDF for PIO
  CMAKE_OPTS += -D NETCDF_DIR:PATH=$(NETCDF_C_PATH) \
                -D NETCDF_Fortran_DIR:PATH=$(NETCDF_FORTRAN_PATH)
endif

ifdef PNETCDF_PATH
	CMAKE_OPTS += -D PNETCDF_DIR:STRING="$(PNETCDF_PATH)"
else
        CMAKE_OPTS += -D WITH_PNETCDF:LOGICAL=FALSE -D PIO_USE_MPIIO:LOGICAL=FALSE
endif
ifdef PIO_FILESYSTEM_HINTS
	CMAKE_OPTS += -D PIO_FILESYSTEM_HINTS:STRING="$(PIO_FILESYSTEM_HINTS)"
endif

# This captures the many cism-specific options to cmake
CMAKE_OPTS += $(USER_CMAKE_OPTS)

# CMake doesn't seem to like it when you define compilers via -D
# CMAKE_C_COMPILER, etc., when you rerun cmake with an existing
# cache. So doing this via environment variables instead.
ifndef CMAKE_ENV_VARS
  CMAKE_ENV_VARS :=
endif
CMAKE_ENV_VARS += CC=$(CC) \
                  CXX=$(CXX) \
                  FC=$(FC) \
                  LDFLAGS="$(LDFLAGS)"


# We declare $(GLC_DIR)/Makefile to be a phony target so that cmake is
# always rerun whenever invoking 'make $(GLC_DIR)/Makefile'; this is
# desirable to pick up any new source files that may have been added
.PHONY: $(GLC_DIR)/Makefile
$(GLC_DIR)/Makefile:
	cd $(GLC_DIR); \
	$(CMAKE_ENV_VARS) cmake $(CMAKE_OPTS) $(CIMEROOT)/../components/cism/glimmer-cism

$(PIO_LIBDIR)/Makefile:
	cd $(PIO_LIBDIR); \
	$(CMAKE_ENV_VARS) cmake $(CMAKE_OPTS) $(PIO_SRC_DIR)

#-------------------------------------------------------------------------------
# Build & include dependency files
#-------------------------------------------------------------------------------

touch_filepath:
	touch $(CURDIR)/Filepath

# Get list of files and build dependency file for all .o files
#   using perl scripts mkSrcfiles and mkDepends
# if a source is of form .F90.in strip the .in before creating the list of objects
SOURCES := $(shell cat Srcfiles)
BASENAMES := $(basename $(basename $(SOURCES)))
OBJS    := $(addsuffix .o, $(BASENAMES))
INCS    := $(foreach dir,$(cpp_dirs),-I$(dir))

CURDIR := $(shell pwd)

$(CURDIR)/Depends: $(CURDIR)/Srcfiles $(CURDIR)/Deppath
	$(CASETOOLS)/mkDepends $(USER_MKDEPENDS_OPTS) Deppath Srcfiles > $@

$(CURDIR)/Deppath: $(CURDIR)/Filepath
	$(CP) -f $(CURDIR)/Filepath $@
	@echo "$(MINCROOT)" >> $@

$(CURDIR)/Srcfiles: $(CURDIR)/Filepath
	$(CASETOOLS)/mkSrcfiles

$(CURDIR)/Filepath:
	@echo "$(VPATH)" > $@


#-------------------------------------------------------------------------------
# echo file names, paths, compile flags, etc. used during build
#-------------------------------------------------------------------------------

db_files:
	@echo " "
	@echo "* MACFILE := $(MACFILE)"
	@echo "* VPATH   := $(VPATH)"
	@echo "* INCS    := $(INCS)"
	@echo "* OBJS    := $(OBJS)"
db_flags:
	@echo " "
	@echo "* cc      := $(CC)  $(CFLAGS) $(INCS) $(INCLDIR)"
	@echo "* .F.o    := $(FC)  $(FFLAGS) $(FIXEDFLAGS) $(INCS) $(INCLDIR)"
	@echo "* .F90.o  := $(FC)  $(FFLAGS) $(FREEFLAGS) $(INCS) $(INCLDIR)"
	ifeq ($(USE_CXX), true)
	  @echo "* .cpp.o  := $(CXX) $(CXXFLAGS) $(INCS) $(INCLDIR)"
	endif

#-------------------------------------------------------------------------------
# Rules used for the tests run by "configure -test"
#-------------------------------------------------------------------------------

test_fc: test_fc.o
	$(LD) -o $@ test_fc.o $(LDFLAGS)
ifeq ($(NETCDF_SEPARATE), false)
test_nc: test_nc.o
	$(LD) -o $@ test_nc.o -L$(LIB_NETCDF) -lnetcdff -lnetcdf $(LDFLAGS)
else ifeq ($(NETCDF_SEPARATE), true)
test_nc: test_nc.o
	$(LD) -o $@ test_nc.o -L$(LIB_NETCDF_FORTRAN) -L$(LIB_NETCDF_C) -lnetcdff -lnetcdf $(LDFLAGS)
endif
test_mpi: test_mpi.o
	$(LD) -o $@ test_mpi.o $(LDFLAGS)
test_esmf: test_esmf.o
	$(LD) -o $@ test_esmf.o $(LDFLAGS)

#-------------------------------------------------------------------------------
# create list of component libraries - hard-wired for current ccsm components
#-------------------------------------------------------------------------------

ifeq ($(ULIBDEP),$(null))
   ifneq ($(LIBROOT),$(null))
     ULIBDEP += $(LIBROOT)/libatm.a
     ULIBDEP += $(LIBROOT)/libice.a
     ifeq ($(findstring clm5_0,$(CLM_CONFIG_OPTS)),clm5_0)
          LNDLIB := libclm.a
     else
          ifeq ($(findstring clm4_5,$(CLM_CONFIG_OPTS)),clm4_5)
               LNDLIB := libclm.a
          else
               LNDLIB := liblnd.a
          endif
     endif
     ifeq ($(findstring libclm.a,$(LNDLIB)),libclm.a)
        ULIBDEP += $(INSTALL_SHAREDPATH)/$(COMP_INTERFACE)/$(ESMFDIR)/lib/$(LNDLIB)
     else
        ULIBDEP += $(LIBROOT)/$(LNDLIB)
     endif
     ifeq ($(MODEL),driver)
         ifeq ($(findstring libclm.a,$(LNDLIB)),libclm.a)
             INCLDIR += -I$(SHAREDLIBROOT)/$(SHAREDPATH)/$(COMP_INTERFACE)/$(ESMFDIR)/clm/obj
         else
             INCLDIR += -I$(EXEROOT)/lnd/obj
         endif
     endif


     ULIBDEP += $(LIBROOT)/libocn.a
     ULIBDEP += $(LIBROOT)/librof.a
     ULIBDEP += $(LIBROOT)/libglc.a
     ULIBDEP += $(LIBROOT)/libwav.a
     ULIBDEP += $(LIBROOT)/libesp.a
     ifeq ($(COMP_GLC), cism)
       ULIBDEP += $(CISM_LIBDIR)/libglimmercismfortran.a
       ifeq ($(CISM_USE_TRILINOS), TRUE)
         ULIBDEP += $(CISM_LIBDIR)/libglimmercismcpp.a
       endif
     endif
     ifeq ($(OCN_SUBMODEL),moby)
       ULIBDEP += $(LIBROOT)/libmoby.a
     endif
   endif
endif

ifdef COSP_LIBDIR
  ULIBDEP += $(COSP_LIBDIR)/libcosp.a
endif


ifndef CLIBS
   ifdef ULIBDEP
     # For each occurrence of something like /path/to/foo/libbar.a in ULIBDEP,
     # CLIBS will contain -L/path/to/foo -lbar
     CLIBS := $(foreach LIBDEP,$(strip $(ULIBDEP)), -L$(dir $(LIBDEP)) $(patsubst lib%.a,-l%,$(notdir $(LIBDEP))))
   endif
endif

# libcsm_share.a is in ULIBDEP, but -lcsm_share is in ULIBS rather than CLIBS,
# so this needs to be added after creating CLIBS above
CSMSHARELIB = $(INSTALL_SHAREDPATH)/$(COMP_INTERFACE)/$(ESMFDIR)/$(NINST_VALUE)/lib/libcsm_share.a
ULIBDEP += $(CSMSHARELIB)

#-------------------------------------------------------------------------------
# build rules:
#-------------------------------------------------------------------------------

.SUFFIXES:
.SUFFIXES: .F90 .F .f90 .f .c .cpp .o .in

ifeq ($(MPILIB),mpi-serial)
  MPISERIAL = $(INSTALL_SHAREDPATH)/lib/libmpi-serial.a
  MLIBS += $(MPISERIAL)
  CMAKE_OPTS += -DMPI_C_INCLUDE_PATH=$(INSTALL_SHAREDPATH)/include \
      -DMPI_Fortran_INCLUDE_PATH=$(INSTALL_SHAREDPATH)/include \
      -DMPI_C_LIBRARIES=$(INSTALL_SHAREDPATH)/lib/libmpi-serial.a \
      -DMPI_Fortran_LIBRARIES=$(INSTALL_SHAREDPATH)/lib/libmpi-serial.a
endif

$(MCTLIBS)  : $(MPISERIAL)

$(PIOLIB) : $(MPISERIAL) $(GPTLLIB)

$(CSMSHARELIB):  $(MCTLIBS) $(PIOLIB) $(GPTLLIB)
ifneq ($(MODEL),csm_share)
  $(OBJS):  $(CSMSHARELIB)
endif

$(EXEC_SE): $(OBJS) $(ULIBDEP) $(CSMSHARELIB) $(MCTLIBS) $(PIOLIB) $(GPTLLIB)
	$(LD) -o $(EXEC_SE) $(OBJS) $(CLIBS) $(ULIBS) $(SLIBS) $(MLIBS) $(LDFLAGS)

$(COMPLIB): $(OBJS)
	$(AR) -r $(COMPLIB) $(OBJS)


.c.o:
	$(CC) -c $(INCLDIR) $(INCS) $(CFLAGS)  $<
.F.o:
	$(FC) -c $(INCLDIR) $(INCS) $(FFLAGS) $(FIXEDFLAGS) $<
.f.o:
	$(FC) -c $(INCLDIR) $(INCS) $(FFLAGS) $(FIXEDFLAGS) $<
.f90.o:
	$(FC) -c $(INCLDIR) $(INCS) $(FFLAGS) $(FREEFLAGS)  $<
.F90.o:
	$(FC) -c $(INCLDIR) $(INCS) $(FFLAGS) $(FREEFLAGS) $(CONTIGUOUS_FLAG) $<
.cpp.o:
	$(CXX) -c $(INCLDIR) $(INCS) $(CXXFLAGS)  $<

%.F90: %.F90.in
	$(CIMEROOT)/externals/genf90/genf90.pl $< > $@

cleanatm:
	$(RM) -f $(LIBROOT)/libatm.a
	cd $(EXEROOT)/atm/obj;  $(RM) -f *.o *.mod

cleancpl:
	cd $(EXEROOT)/cpl/obj;  $(RM) -f *.o *.mod

cleanocn:
	$(RM) -f $(LIBROOT)/libocn.a
	cd $(EXEROOT)/ocn/obj ; $(RM) -f *.o *.mod

cleanwav:
	$(RM) -f $(LIBROOT)/libwav.a
	cd $(EXEROOT)/wav/obj ; $(RM) -f *.o *.mod

cleanesp:
	$(RM) -f $(LIBROOT)/libesp.a
	cd $(EXEROOT)/esp/obj ; $(RM) -f *.o *.mod

cleanglc:
	$(RM) -f $(LIBROOT)/libglc.a
	$(RM) -fr $(EXEROOT)/glc

cleanice:
	$(RM) -f $(LIBROOT)/libice.a
	cd $(EXEROOT)/ice/obj ; $(RM) -f *.o *.mod

cleanrof:
	$(RM) -f $(LIBROOT)/librof.a
	cd $(EXEROOT)/rof/obj ; $(RM) -f *.o *.mod

cleanlnd:
	$(RM) -f $(LIBROOT)/liblnd.a
	cd $(EXEROOT)/lnd/obj ; $(RM) -f *.o *.mod

cleancsmshare:
	$(RM) -f $(CSMSHARELIB)
	$(RM) -fr $(SHAREDLIBROOT)/$(SHAREDPATH)/$(COMP_INTERFACE)/$(ESMFDIR)/$(NINST_VALUE)/csm_share

cleanpio:
	$(RM) -f $(PIO_LIBDIR)/libpio*
	$(RM) -fr $(SHAREDLIBROOT)/$(SHAREDPATH)/pio

cleanmct:
	$(RM) -f $(MCTLIBS)
	$(RM) -fr $(SHAREDLIBROOT)/$(SHAREDPATH)/mct

cleangptl:
	$(RM) -f $(GPTLLIB)
	$(RM) -fr $(SHAREDLIBROOT)/$(SHAREDPATH)/gptl

clean: cleanatm cleanocn cleanwav cleanglc cleanice cleanrof cleanlnd cleanesp

realclean: clean cleancsmshare cleanpio cleanmct cleangptl

# the if-tests prevent DEPS files from being created when they're not needed
ifneq ($(MAKECMDGOALS), db_files)
ifneq ($(MAKECMDGOALS), db_flags)
ifeq (,$(findstring clean,$(MAKECMDGOALS)))
    -include $(CURDIR)/Depends $(CASEROOT)/Depends.$(COMPILER) $(CASEROOT)/Depends.$(MACH) $(CASEROOT)/Depends.$(MACH).$(COMPILER)
endif
endif
endif
