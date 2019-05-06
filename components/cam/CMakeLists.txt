#===============================================================================
#
# Common CMakeLists.txt: a framework for building all CIME components and more
#
# This is a port of cime/scripts/Tools/Makefile. As more components are ported to
# CMake, the directory level of this file will rise to the top-level directory.
#
# We will prefer space-separated strings over lists
#
#===============================================================================

cmake_minimum_required(VERSION 3.3)
cmake_policy(SET CMP0057 NEW)

project(E3SM C CXX Fortran)

# Load dependency search path.
file(STRINGS Filepath FILEPATH_DIRS)
string(REPLACE ";" " " DIRS "${FILEPATH_DIRS}")
set(DIRS ". ${DIRS}")

set(CPP_DIRS ${DIRS})
# Add INCROOT to path for Depends and Include
set(MINCROOT "")
if (INCROOT)
  set(CPP_DIRS "${CPP_DIRS} ${INCROOT}")
  set(MINCROOT ${INCROOT})
endif ()

# # Expand any tildes in directory names. Change spaces to colons. JGF
# VPATH := $(foreach dir,$(cpp_dirs),$(wildcard $(dir)))
# VPATH := $(subst $(space),:,$(VPATH))

# RM    := rm
# CP    := cp

# exec_se: $(EXEC_SE)  Depends
# complib: $(COMPLIB)  Depends

# Determine whether to compile threaded or not
# Set the THREADDIR for the shared build
# based on the threaded build status
if (SMP OR SMP_PRESENT)
  set(THREADDIR "threads")
  set(compile_threaded TRUE)
else()
  set(THREADDIR "nothreads")
  set(compile_threaded FALSE)
endif()

# set the debug directory based on the debug status
if (DEBUG)
  set(DEBUGDIR "debug")
else()
  set(DEBUGDIR "nodebug")
endif()

if (USE_ESMF_LIB)
  set(ESMFDIR "esmf")
else()
  set(ESMFDIR "noesmf")
endif()

# Determine whether any C++ code will be included in the build;
# currently, C++ code is included if and only if we're linking to the
# trilinos library or the Albany library.
set(USE_CXX FALSE)
if (USE_TRILINOS OR USE_ALBANY OR USE_KOKKOS)
  set(USE_CXX TRUE)
endif()

if (NOT MOD_SUFFIX)
  set(MOD_SUFFIX "mod")
endif()

#===============================================================================
# set CPP options (must use this before any flags or cflags settings)
#===============================================================================

set(CPPDEFS ${USER_CPPDEFS} "-D${OS}")

include(${CASEROOT}/Macros)

# SLIBS comes from Macros, so this append must come after Macros are included
if (USE_FMS)
  set(SLIBS "${SLIBS} -lfms")
endif()

if (DEBUG)
  if (CIME_MODEL STREQUAL "cesm")
    set(CPPDEFS "${CPPDEFS} -DDEBUG")
  endif()
else()
  set(CPPDEFS "${CPPDEFS} -DNDEBUG")
endif()

if (USE_ESMF_LIB)
  set(CPPDEFS "${CPPDEFS} -DUSE_ESMF_LIB")
endif()

if (COMP_INTERFACE STREQUAL "nuopc")
  set(CPPDEFS "${CPPDEFS} -DNUOPC_INTERFACE")
else()
  set(CPPDEFS "${CPPDEFS} -DMCT_INTERFACE")
endif()

if (COMPARE_TO_NUOPC)
  set(CPPDEFS "${CPPDEFS} -DCOMPARE_TO_NUOPC")
endif()

if (MPILIB STREQUAL mpi-serial)
  set(CPPDEFS "${CPPDEFS} -DNO_MPI2")
else()
  set(CPPDEFS "${CPPDEFS} -DHAVE_MPI")
endif()

if (PIO_VERSION STREQUAL "1")
  set(CPPDEFS "${CPPDEFS} -DPIO1")
else()
  set(USE_CXX TRUE)
endif()

if (NOT SHAREDPATH)
  set(SHAREDPATH ${COMPILER}/${MPILIB}/${DEBUGDIR}/${THREADDIR}/${COMP_INTERFACE})
  set(INSTALL_SHAREDPATH ${EXEROOT}/${SHAREDPATH})
endif()

# Decide whether to use a C++ or Fortran linker, based on whether we
# are using any C++ code and the compiler-dependent CXX_LINKER variable
if (USE_CXX AND CXX_LINKER STREQUAL "CXX")
  if (MPILIB STREQUAL mpi-serial)
    set(LD ${SCXX})
  else()
    set(LD ${MPICXX})
  endif()
else()
  if (MPILIB STREQUAL mpi-serial)
    set(LD ${SFC})
  else()
    set(LD ${MPIFC})
  endif()
endif()

if (USE_CXX AND NOT SUPPORTS_CXX)
  message(FATAL_ERROR "Fatal attempt to include C++ code on a compiler/machine combo that has not been set up to support C++")
endif()

# Not clear how to escape commas for libraries with their own configure
# script, and they don't need this defined anyway, so leave this out of
# FPPDEFS. JGF
if (HAS_F2008_CONTIGUOUS)
  set(CONTIGUOUS_FLAG "-DUSE_CONTIGUOUS=contiguous,")
else()
  set(CONTIGUOUS_FLAG "-DUSE_CONTIGUOUS")
endif()

if (CPRE)
  set(CONTIGUOUS_FLAG ${CPRE}) # JGF
endif()

if (NOT AR)
  set(AR "ar")
endif()
if (NOT ARFLAGS)
  set(ARFLAGS "-r")
endif()

if (NETCDF_C_PATH)
  if (NOT NETCDF_FORTRAN_PATH)
    message(FATAL_ERROR "NETCDF_C_PATH specified without NETCDF_FORTRAN_PATH")
  endif()
  set(NETCDF_SEPARATE TRUE)
  if (NOT INC_NETCDF_C)
    set(INC_NETCDF_C ${NETCDF_C_PATH}/include)
  endif()
  if (NOT INC_NETCDF_FORTRAN)
    set(INC_NETCDF_FORTRAN ${NETCDF_FORTRAN_PATH}/include)
  endif()
  if (NOT LIB_NETCDF_C)
    set(LIB_NETCDF_C ${NETCDF_C_PATH}/lib)
  endif()
  if (NOT LIB_NETCDF_FORTRAN)
    set(LIB_NETCDF_FORTRAN ${NETCDF_FORTRAN_PATH}/lib)
  endif()
elseif (NETCDF_FORTRAN_PATH)
  message(FATAL_ERROR "NETCDF_FORTRAN_PATH specified without NETCDF_C_PATH")
elseif (NETCDF_PATH)
  set(NETCDF_SEPARATE FALSE)
  if (NOT INC_NETCDF)
    set(INC_NETCDF ${NETCDF_PATH}/include)
  endif()
  if (NOT LIB_NETCDF)
    set(LIB_NETCDF ${NETCDF_PATH}/lib)
  endif()
else()
  message(FATAL_ERROR "NETCDF not found: Define NETCDF_PATH or NETCDF_C_PATH and NETCDF_FORTRAN_PATH in config_machines.xml or config_compilers.xml")
endif()

if (MPILIB STREQUAL mpi-serial)
  if (PNETCDF_PATH)
    unset(PNETCDF_PATH)
  endif()
else()
  if (PNETCDF_PATH)
    if (NOT INC_PNETCDF)
      set(INC_PNETCDF ${PNETCDF_PATH}/include)
    endif()
    if (NOT LIB_PNETCDF)
      set(LIB_PNETCDF ${PNETCDF_PATH}/lib)
    endif()
  endif()
endif()

# Set PETSc info if it is being used
if (USE_PETSC)
  if (PETSC_PATH)
    if (NOT INC_PETSC)
      set(INC_PETSC ${PETSC_PATH}/include)
    endif()
    if (NOT LIB_PETSC)
      set(LIB_PETSC ${PETSC_PATH}/lib)
    endif()
  else()
    message(FATAL_ERROR "PETSC_PATH must be defined when USE_PETSC is TRUE")
  endif()

  # # Get the "PETSC_LIB" list an env var  JGF
  #include $(PETSC_PATH)/lib/petsc/conf/variables
endif()

if (USE_TRILINOS)
  if (TRILINOS_PATH)
    if (NOT INC_TRILINOS)
      set(INC_TRILINOS ${TRILINOS_PATH}/include)
    endif()
    if (NOT LIB_TRILINOS)
      set(LIB_TRILINOS ${TRILINOS_PATH}/lib)
    endif()
  else()
    message(FATAL_ERROR "TRILINOS_PATH must be defined when USE_TRILINOS is TRUE")
  endif()

  # # get a bunch of variables related to this trilinos installation;
  # these variables begin with "Trilinos_" JGF
  #include $(INC_TRILINOS)/Makefile.export.Trilinos
endif()

if (USE_ALBANY)
  if (ALBANY_PATH)
    if (NOT INC_ALBANY)
      set(INC_ALBANY ${ALBANY_PATH}/include)
    endif()
    if (NOT LIB_ALBANY)
      set(LIB_ALBANY ${ALBANY_PATH}/lib)
    endif()
  else()
    message(FATAL_ERROR "ALBANY_PATH must be defined when USE_ALBANY is TRUE")
  endif()

  # # get the "ALBANY_LINK_LIBS" list as an env var JGF
  #include $(ALBANY_PATH)/export_albany.in
endif()

if (USE_KOKKOS)
  include(${INSTALL_SHAREDPATH}/Kokkos)
  # JGF - convert kokkos lists to strings
  set(SLIBS "${SLIBS} ${KOKKOS_LIBS}")
  set(CXXFLAGS "${CXXFLAGS} ${KOKKOS_CXXFLAGS}")
  set(CXX_LDFLAGS "${CXX_LDFLAGS} ${KOKKOS_LDFLAGS}")
endif()

if (USE_MOAB)
  if (MOAB_PATH)
    set(CPPDEFS "${CPPDEFS} -DHAVE_MOAB")
    if (NOT INC_MOAB)
      set(INC_MOAB ${MOAB_PATH}/include)
    endif()
    if (NOT LIB_MOAB)
      set(LIB_MOAB ${MOAB_PATH}/lib)
    endif()
  else()
    message(FATAL_ERROR "MOAB_PATH must be defined when USE_MOAB is TRUE")
  endif()

  # # get the "IMESH_LIBS" list as an env var JGF
  #include $(LIB_MOAB)/iMesh-Defs.inc
endif()

# Set HAVE_SLASHPROC on LINUX systems which are not bluegene or Darwin (OSx)
string(FIND CPPDEFS "-DLINUX" HAS_DLINUX)
string(FIND CPPDEFS "DBG" HAS_DBG)
string(FIND CPPDEFS "Darwin" HAS_DARWIN)
if (NOT HAS_DLINUX EQUAL -1 AND HAS_DBG EQUAL -1 AND HASDARWIN EQUAL -1)
  set(CPPDEFS "${CPPDEFS} -DHAVE_SLASHPROC")
endif()

# Atleast on Titan+cray mpi, MPI_Irsends() are buggy, causing hangs during I/O
# Force PIO to use MPI_Isends instead of the default, MPI_Irsends
if (PIO_VERSION STREQUAL 2)
  set(EXTRA_PIO_CPPDEFS "-DUSE_MPI_ISEND_FOR_FC")
else()
  set(EXTRA_PIO_CPPDEFS "-D_NO_MPI_RSEND")
endif

if (LIB_PNETCDF)
  set(CPPDEFS "${CPPDEFS} -D_PNETCDF")
  set(SLIBS "${SLIBS} -L${LIB_PNETCDF} -lpnetcdf")
endif()

# Set esmf.mk location with ESMF_LIBDIR having precedent over ESMFMKFILE
set(CIME_ESMFMKFILE "undefined_ESMFMKFILE")
if (ESMFMKFILE)
  set(CIME_ESMFMKFILE ${ESMFMKFILE})
endif()
if (ESMF_LIBDIR)
  set(CIME_ESMFMKFILE ${ESMF_LIBDIR}/esmf.mk)
endif()

# For compiling and linking with external ESMF.
# If linking to external ESMF library then include esmf.mk
# ESMF_F90COMPILEPATHS
# ESMF_F90LINKPATHS
# ESMF_F90LINKRPATHS
# ESMF_F90ESMFLINKLIBS
if (USE_ESMF_LIB)
  # include(${CIME_ESMFMKFILE}) # JGF
  set(CPPDEFS "${CPPDEFS} -DESMF_VERSION_MAJOR=${ESMF_VERSION_MAJOR} -DESMF_VERSION_MINOR=${ESMF_VERSION_MINOR}")
  set(FFLAGS "${FFLAGS} ${ESMF_F90COMPILEPATHS}")
  set(SLIBS "${SLIBS} ${ESMF_F90LINKPATHS} ${ESMF_F90LINKRPATHS} ${ESMF_F90ESMFLINKLIBS}")
endif()

# JGF
# ifdef CPRE
#   FPPDEFS := $(subst $(comma),\\$(comma),$(CPPDEFS))
#   FPPDEFS := $(patsubst -D%,$(CPRE)%,$(FPPDEFS))
#   EXTRA_PIO_FPPDEFS := $(subst $(comma),\\$(comma),$(EXTRA_PIO_CPPDEFS))
#   EXTRA_PIO_FPPDEFS := $(patsubst -D%,$(CPRE)%,$(EXTRA_PIO_FPPDEFS))
# else
#   FPPDEFS := $(CPPDEFS)
#   EXTRA_PIO_FPPDEFS := $(EXTRA_PIO_CPPDEFS)
# endif

#===============================================================================
# Set config args for pio and mct to blank and then enable serial
#===============================================================================
if (NOT CONFIG_ARGS)
  set(CONFIG_ARGS "")
endif()
string(FIND MODEL "pio" HAS_PIO)
if (NOT HAS_PIO EQUAL -1)
  set(CONFIG_ARGS "${CONFIG_ARGS} --enable-timing")
  if (DEBUG)
    set(CONFIG_ARGS "${CONFIG_ARGS} --enable-debug")
  endif()
endif()

#===============================================================================
# User-specified INCLDIR
#===============================================================================

set(INCLDIR "-I.")
if (USER_INCLDIR)
  set(INCLDIR "${INCLDIR} ${USER_INCLDIR}")
endif()

#===============================================================================
# MPI-serial library (part of MCT)
#===============================================================================

if (MPILIB STREQUAL "mpi-serial")
  set(CC ${SCC})
  set(FC ${SFC})
  set(CXX ${SCXX})
  set(MPIFC ${SFC})
  set(MPICC ${SCC})
  set(MPICXX ${SCXX})
  set(CONFIG_ARGS "${CONFIG_ARGS} MCT_PATH=${SHAREDLIBROOT}/${SHAREDPATH}/mct/mpi-serial")
else()
  set(CC ${MPICC})
  set(FC ${MPIFC})
  set(CXX ${MPICXX})
  if (MPI_PATH)
    set(INC_MPI ${MPI_PATH}/include)
    set(LIB_MPI ${MPI_PATH}/lib)
  endif()
endif()
set(CSM_SHR_INCLUDE ${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/${NINST_VALUE}/include)
# VPATH+=$(CSM_SHR_INCLUDE) # JGF

#===============================================================================
# Set include paths (needed after override for any model specific builds below)
#===============================================================================
set(INCLDIR "${INCLDIR} -I${INSTALL_SHAREDPATH}/include -I${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/${NINST_VALUE}/include")

if (NOT NETCDF_SEPARATE)
  set(INCLDIR "${INCLDIR} -I${INC_NETCDF}")
else()
  set(INCLDIR "${INCLDIR} -I${INC_NETCDF_C} -I${INC_NETCDF_FORTRAN}")
endif()

foreach(ITEM MOD_NETCDF INC_MPI INC_PNETCDF INC_PETSC INC_TRILINOS INC_ALBANY INC_MOAB)
  if (${ITEM})
    set(INCLDIR "${INCLDIR} -I${${ITEM}}")
  endif()
endforeach()

if (MODEL STREQUAL "driver")
  set(INCLDIR "${INCLDIR} -I${EXEROOT}/atm/obj -I${EXEROOT}/ice/obj -I${EXEROOT}/ocn/obj -I${EXEROOT}/glc/obj -I${EXEROOT}/rof/obj -I${EXEROOT}/wav/obj -I${EXEROOT}/esp/obj -I${EXEROOT}/iac/obj")
  # nagfor and gcc have incompatible LDFLAGS.
  # nagfor requires the weird "-Wl,-Wl,," syntax.
  # If done in config_compilers.xml, we break MCT.
  if (COMPILER STREQUAL "nag")
    if (NOT NETCDF_SEPARATE)
      set(SLIBS "${SLIBS} -Wl,-Wl,,-rpath=${LIB_NETCDF}")
    else()
      set(SLIBS "${SLIBS} -Wl,-Wl,,-rpath=${LIB_NETCDF_C}")
      set(SLIBS "${SLIBS} -Wl,-Wl,,-rpath=${LIB_NETCDF_FORTRAN}")
    endif()
  endif()
else()
  if (COMPILER STREQUAL "nag")
    if (DEBUG)
      if (MACH STREQUAL "hobart")
       # GCC needs to be able to link to
       # nagfor runtime to get autoconf
       # tests to work.
       set(CFLAGS "${CFLAGS} -Wl,--as-needed,--allow-shlib-undefined")
       set(SLIBS "${SLIBS} -L${COMPILER_PATH}/lib/NAG_Fortran -lf62rts")
     endif()
   endif()
 endif()
endif()

if (NOT MCT_LIBDIR)
  set(MCT_LIBDIR "${INSTALL_SHAREDPATH}/lib")
endif()

if (PIO_LIBDIR)
  if (PIO_VERSION STREQUAL ${PIO_VERSION_MAJOR})
    set(INCLDIR "${INCLDIR} -I${PIO_INCDIR}")
    set(SLIBS "${SLIBS} -L${PIO_LIBDIR}")
  else()
    # If PIO_VERSION_MAJOR doesnt match, build from source
    unset(PIO_LIBDIR)
  endif()
endif()
if (NOT PIO_LIBDIR)
  set(PIO_LIBDIR "${INSTALL_SHAREDPATH}/lib")
endif()

if (NOT GPTL_LIBDIR)
  set(GPTL_LIBDIR "${INSTALL_SHAREDPATH}/lib")
endif()

if (NOT GLC_DIR)
  set(GLC_DIR "${EXEROOT}/glc")
endif()

if (NOT CISM_LIBDIR)
  set(CISM_LIBDIR "${GLC_DIR}/lib")
endif()

if (NOT GLCROOT)
  # Backwards compatibility
  set(GLCROOT "${CIMEROOT}/../components/cism")
endif()

set(INCLDIR "-I${INSTALL_SHAREDPATH}/include")

#
# Use the MCT dir for the cache for all configure calls because it is the first one
#
set(CFLAGS "${CPPDEFS}")
set(CXXFLAGS "${CPPDEFS}")
set(CONFIG_ARGS "CC=\"${CC}\" FC=\"${FC}\" MPICC=\"${MPICC}\" MPIFC=\"${MPIFC}\" FCFLAGS=\"${FFLAGS} ${FREEFLAGS} ${INCLDIR}\" CPPDEFS=\"${CPPDEFS}\" CFLAGS=\"${CFLAGS} -I.. ${INCLDIR}\" LDFLAGS=\"${LDFLAGS}\"")

if (NOT NETCDF_SEPARATE)
  set(CONFIG_ARGS "${CONFIG_ARGS} NETCDF_PATH=${NETCDF_PATH}")
else()
  set(CONFIG_ARGS "${CONFIG_ARGS} NETCDF_PATH=${NETCDF_C_PATH}")
endif()

if (COMPILER STREQUAL "nag")
  set(CONFIG_ARGS "LIBS=${SLIBS}")
endif()

set(FFLAGS "${FPPDEFS}")
set(FFLAGS_NOOPT "${FPPDEFS}")

string(FIND CAM_CONFIG_OPTS "-cosp" HAS_COSP)
if (NOT HAS_COSP EQUAL -1)
  # The following is for the COSP simulator code:
  get_filename_component(COSP_LIBDIR ${EXEROOT}/atm/obj/cosp ABSOLUTE)
  if (MODEL STREQUAL "driver")
    set(INCLDIR "${INCLDIR} -I${COSP_LIBDIR}")
  endif()
endif()


if (MODEL STREQUAL "cam")
  # These RRTMG files take an extraordinarily long time to compile with optimization.
  # Until mods are made to read the data from files, just remove optimization from
  # their compilation.

  # JGF
# rrtmg_lw_k_g.o: rrtmg_lw_k_g.f90
# 	$(FC) -c $(FPPFLAGS) $(INCLDIR) $(INCS) $(FREEFLAGS) $(FFLAGS_NOOPT) $<

# rrtmg_sw_k_g.o: rrtmg_sw_k_g.f90
# 	$(FC) -c $(FPPFLAGS) $(INCLDIR) $(INCS) $(FREEFLAGS) $(FFLAGS_NOOPT) $<

  if (COSP_LIBDIR)
    set(INCLDIR "-I${COSP_LIBDIR} -I${COSP_LIBDIR}/../ -I../${INSTALL_SHAREDPATH}/include -I../${CSM_SHR_INCLUDE}")
    # JGF
    add_library(libcosp cam_abortutils.f90)
  endif()

endif()

# System libraries (netcdf, mpi, pnetcdf, esmf, trilinos, etc.)
if (NOT SLIBS)
  if (NOT NETCDF_SEPARATE)
    set(SLIBS "${SLIBS} -L${LIB_NETCDF} -lnetcdff -lnetcdf")
  else()
    set(SLIBS "${SLIBS} -L${LIB_NETCDF_FORTRAN} -L${LIB_NETCDF_C} -lnetcdff -lnetcdf")
  endif()
endif()

if (LAPACK_LIBDIR)
  set(SLIBS "${SLIBS} -L${LAPACK_LIBDIR} -llapack -lblas")
endif()

if (LIB_MPI)
  if (NOT MPI_LIB_NAME)
    set(SLIBS "${SLIBS} -L${LIB_MPI} -lmpi")
  else()
    set(SLIBS "${SLIBS} -L${LIB_MPI} -l${MPI_LIB_NAME}")
  endif()
endif()

# Add PETSc libraries
if (USE_PETSC)
  set(SLIBS "${SLIBS} ${PETSC_LIB}")
endif()

# Add trilinos libraries; too be safe, we include all libraries included in the trilinos build,
# as well as all necessary third-party libraries
if (USE_TRILINOS)
  set(SLIBS "${SLIBS} -L${LIB_TRILINOS} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARIES}")
endif()

# Add Albany libraries.  These are defined in the ALBANY_LINK_LIBS env var that was included above
if (USE_ALBANY)
  set(SLIBS "${SLIBS} ${ALBANY_LINK_LIBS}")
endif()

# Add MOAB libraries.  These are defined in the MOAB_LINK_LIBS env var that was included above
if (USE_MOAB)
  set(SLIBS "${SLIBS} ${IMESH_LIBS}")
endif()

# Add libraries and flags that we need on the link line when C++ code is included
# We need to do these additions after CONFIG_ARGS is set, because they can sometimes break configure for mct, etc.,
# if they are added to LDFLAGS in CONFIG_ARGS.
if (USE_CXX)
  if (CXX_LIBS)
    set(SLIBS "${SLIBS} ${CXX_LIBS}")
  endif()

  if (CXX_LDFLAGS)
    set(LDFLAGS "${SLIBS} ${CXX_LDFLAGS}")
  endif
endif

# Remove arch flag if it exists
string(REGEX REPLACE "-arch[^ ]+" "" F90_LDFLAGS "${LDFLAGS}")

#------------------------------------------------------------------------------
# Drive configure scripts for support libraries (mct)
#------------------------------------------------------------------------------

if (NOT EXISTS ${SHAREDLIBROOT}/${SHAREDPATH}/mct/Makefile.conf)
  execute_process(COMMAND ${CIMEROOT}/src/externals/mct/configure ${CONFIG_ARGS} --srcdir ${CIMEROOT}/src/externals/mct)
endif()

if (NOT EXISTS ${SHAREDLIBROOT}/${SHAREDPATH}/mct/mpi-serial/Makefile.conf)
  execute_process(COMMAND ${CIMEROOT}/src/externals/mct/mpi-serial/configure ${CONFIG_ARGS} --srcdir ${CIMEROOT}/src/externals/mct/mpi-serial)
endif()

if (PIO_VERSION STREQUAL 2)
  # This is a pio2 library
  set(PIOLIB "${PIO_LIBDIR}/libpiof.a ${PIO_LIBDIR}/libpioc.a")
  set(PIOLIBNAME "-lpiof -lpioc")
  set(PIO_SRC_DIR "${CIMEROOT}/src/externals/pio2")
else()
  # This is a pio1 library
  set(PIOLIB "${PIO_LIBDIR}/libpio.a")
  set(PIOLIBNAME "-lpio")
  if (NOT EXISTS "${CIMEROOT}/src/externals/pio1/pio")
    set(PIO_SRC_DIR "${CIMEROOT}/src/externals/pio1")
  else()
    set(PIO_SRC_DIR "${CIMEROOT}/src/externals/pio1/pio")
  endif()
endif()

set(MCTLIBS "${MCT_LIBDIR}/libmct.a ${MCT_LIBDIR}/libmpeu.a")

set(GPTLLIB "${GPTL_LIBDIR}/libgptl.a")

set(ULIBS "-L${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/${NINST_VALUE}/lib -lcsm_share -L${INSTALL_SHAREDPATH}/lib ${PIOLIBNAME} -lgptl -lmct -lmpeu")

#------------------------------------------------------------------------------
# Drive cmake script for cism and pio
#------------------------------------------------------------------------------

set(CMAKE_Fortran_FLAGS "${FFLAGS} ${EXTRA_PIO_FPPDEFS} ${INCLDIR}")
set(CMAKE_C_FLAGS "${CFLAGS} ${EXTRA_PIO_CPPDEFS} ${INCLDIR}")
set(CMAKE_CXX_FLAGS "${CXXFLAGS} ${EXTRA_PIO_CPPDEFS} ${INCLDIR}")
set(CMAKE_VERBOSE_MAKEFILE TRUE)
set(GPTL_PATH "${INSTALL_SHAREDPATH}")
set(PIO_ENABLE_TESTS FALSE)
set(PIO_USE_MALLOC TRUE)
set(USER_CMAKE_MODULE_PATH "${CIMEROOT}/src/CMake;${CIMEROOT}/src/externals/pio2/cmake")

if (NOT PNETCDF_PATH)
  set(WITH_PNETCDF FALSE)
endif()

if (MPILIB STREQUAL "mpi-serial")
  set(PIO_USE_MPISERIAL TRUE)
  set(MPISERIAL_PATH ${INSTALL_SHAREDPATH})
endif()

set(CMAKE_C_COMPILER ${CC})
set(CMAKE_CXX_COMPILER ${CXX})
set(CMAKE_FC_COMPILER ${FC})
set(CMAKE_EXE_LINKER_FLAGS ${LDFLAGS})

# JGF: We only want to do CAM for now
# add_directory(${GLCROOT}/source_cism)
# add_directory(${PIO_SRC_DIR})

#-------------------------------------------------------------------------------
# Build & include dependency files
#-------------------------------------------------------------------------------

# JGF, I think we can rely on CMake's built-in dependency analysus

# Get src files
execute_process(COMMAND ${CASETOOLS}/mkSrcfiles)
file(STRINGS Srcfiles SOURCES)

foreach(ITEM ${CPP_DIRS})
  set(INCS "${INCS} -I{}")
endforeach()

#-------------------------------------------------------------------------------
# create list of component libraries - hard-wired for current ccsm components
#-------------------------------------------------------------------------------

if (CIME_MODEL STREQUAL "cesm")
  if (COMP_LND STREQUAL "clm")
    set(USE_SHARED_CLM TRUE)
  else()
    set(USE_SHARED_CLM FALSE)
  endif()
else()
  set(USE_SHARED_CLM FALSE)
endif()

if (NOT USE_SHARED_CLM)
  set(LNDOBJDIR "${EXEROOT}/lnd/obj")
  set(LNDLIBDIR "${LIBROOT}")
  if (COMP_LND STREQUAL "clm")
    set(LNDLIB "libclm.a")
  else()
    set(LNDLIB "liblnd.a")
  endif()
  set(INCLDIR "${INCLDIR} -I${LNDOBJDIR}")
else()
  set(LNDLIB "libclm.a")
  set(LNDOBJDIR "${SHAREDLIBROOT}/${SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/clm/obj")
  set(LNDLIBDIR "${EXEROOT}/${SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/lib")
  set(INCLDIR "${INCLDIR} -I${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/include")
  if (MODEL STREQUAL "clm")
    set(INCLUDE_DIR "${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/include")
  endif()
endif()

if (NOT ULIBDEP)
  if (LIBROOT)
    set(ULIBDEP "${ULIBDEP} ${LIBROOT}/libatm.a")
    set(ULIBDEP "${ULIBDEP} ${LIBROOT}/libice.a")
    set(ULIBDEP "${ULIBDEP} ${LNDLIBDIR}/${LNDLIB}")
    set(INCLDIR "${INCLDIR} -I${LNDOBJDIR}")
    set(ULIBDEP "${ULIBDEP} ${LIBROOT}/libocn.a")
    set(ULIBDEP "${ULIBDEP} ${LIBROOT}/librof.a")
    set(ULIBDEP "${ULIBDEP} ${LIBROOT}/libglc.a")
    set(ULIBDEP "${ULIBDEP} ${LIBROOT}/libwav.a")
    set(ULIBDEP "${ULIBDEP} ${LIBROOT}/libiac.a")
    set(ULIBDEP "${ULIBDEP} ${LIBROOT}/libesp.a")
    if (COMP_GLC STREQUAL "cism")
      set(ULIBDEP "${ULIBDEP} ${CISM_LIBDIR}/libglimmercismfortran.a")
      if (CISM_USE_TRILINOS)
        set(ULIBDEP "${ULIBDEP} ${CISM_LIBDIR}/libglimmercismcpp.a")
      endif()
    endif()
    if (OCN_SUBMODEL STRQUAL "moby")
      set(ULIBDEP "${ULIBDEP} ${LIBROOT}/libmoby.a")
    endif()
  endif()
endif()

if (COSP_LIBDIR)
  set(ULIBDEP "${ULIBDEP} ${COSP_LIBDIR}/libcosp.a")
endif()

# JGF, this kind of thing should be done via target_link_library in cmake
if (NOT CLIBS)
  if (ULIBDEP)
    # For each occurrence of something like /path/to/foo/libbar.a in ULIBDEP,
    # CLIBS will contain -L/path/to/foo -lbar
    foreach (ITEM ${ULIBDEP})
      get_filename_component(ULIB_DIR ${ITEM} DIRECTORY)
      get_filename_component(ULIB_LIB ${ITEM} NAME)
      set(CLIBS "${CLIBS} -L${ULIB_DIR} -l${ULIB_LIB}")
    endforeach()
  endif()
endif()

# libcsm_share.a is in ULIBDEP, but -lcsm_share is in ULIBS rather than CLIBS,
# so this needs to be added after creating CLIBS above
set(CSMSHARELIB "${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/${NINST_VALUE}/lib/libcsm_share.a")
set(ULIBDEP "${ULIBDEP} ${CSMSHARELIB}")

#-------------------------------------------------------------------------------
# build rules:
#-------------------------------------------------------------------------------

if (MPILIB STREQUAL "mpi-serial")
  set(MPISERIAL "${INSTALL_SHAREDPATH}/lib/lib-mpi-serial.a")
  set(MLIBS "${MLIBS} ${MPISERIAL}")
endif()

# JGF: how to do genf90?