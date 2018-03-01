# This macro identifies compilers and third-party library needs 
# for particular hosts.
macro(set_up_platform)

  # Set library suffix based on whether we're building shared/static.
  # FIXME: We have to hack this together here, since CMAKE_SHARED_LIBRARY_SUFFIX
  # FIXME: isn't available before project() is called, which is when set_up_platform()
  # FIXME: is invoked. Gross.
  if (BUILD_SHARED_LIBS)
    if (APPLE)
      set(LIB_SUFFIX .dylib)
    elseif (WIN32)
      set(LIB_SUFFIX .dll)
    else()
      set(LIB_SUFFIX .so)
    endif()
  else()
    set(LIB_SUFFIX .a)
  endif()

  # Set defaults for the various third-party libraries. These defaults
  # are hardwired because the project can't have been defined before 
  # this macro is executed, and so PROJECT_BINARY_DIR is unavailable.
  set(Z_LIBRARY "${CMAKE_CURRENT_BINARY_DIR}/lib/libz.a")
  set(Z_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/include")
  get_filename_component(Z_LIBRARY_DIR ${Z_LIBRARY} DIRECTORY)
  set(HDF5_LIB_NAME hdf5)
  set(HDF5_HL_LIB_NAME hdf5_hl)
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(HDF5_LIB_NAME ${HDF5_LIB_NAME}_debug)
    set(HDF5_HL_LIB_NAME ${HDF5_HL_LIB_NAME}_debug)
  endif()
  set(HDF5_LIBRARY "${CMAKE_CURRENT_BINARY_DIR}/lib/lib${HDF5_LIB_NAME}${LIB_SUFFIX}")
  set(HDF5_HL_LIBRARY "${CMAKE_CURRENT_BINARY_DIR}/lib/lib${HDF5_HL_LIB_NAME}${LIB_SUFFIX}")
  set(HDF5_LIBRARIES ${HDF5_HL_LIB_NAME};${HDF5_LIB_NAME})
  set(HDF5_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/include")
  get_filename_component(HDF5_LIBRARY_DIR ${Z_LIBRARY} DIRECTORY)

  if (APPLE)
    set(NEED_LAPACK FALSE)
  else()
    set(NEED_LAPACK TRUE)
  endif()

  # Certain tools (e.g. patch) require TMPDIR to be defined. If it is not, 
  # we do so here.
  set(TMPDIR_VAR $ENV{TMPDIR})
  if (NOT TMPDIR_VAR)
    # FIXME: Does this exist everywhere?
    set(ENV{TMPDIR} "/tmp")
  endif()

  # Get the hostname for this machine. 
  site_name(HOSTNAME)

  if (HOSTNAME MATCHES "cori") # NERSC Cori phase1
    #ndk  make config debug=1 mpi=1 prefix=$SCRATCH/polymec
    # (Intel's compilers don't do C11.).
    set(CMAKE_C_COMPILER $ENV{CC})
    set(CMAKE_CXX_COMPILER $ENV{CXX})
    set(CMAKE_Fortran_COMPILER $ENV{FC})

    # We are cared for mathematically.
    set(NEED_LAPACK FALSE)

  elseif (HOSTNAME MATCHES "edison") # NERSC Edison
    # Edison likes Intel's compilers
    # (but Intel's compilers don't do C11.).
    set(CMAKE_C_COMPILER $ENV{CC})
    set(CMAKE_CXX_COMPILER $ENV{CXX})
    set(CMAKE_Fortran_COMPILER $ENV{FC})

    # We are cared for mathematically.
    set(NEED_LAPACK FALSE)

  elseif(HOSTNAME MATCHES "yslogin") # NCAR yellowstone
    message("-- Running on yellowstone.")
    set(NUM_BUILD_THREADS "4")
    if (FALSE) # gnu
      set(BLAS_INCLUDE_DIRS "/glade/apps/opt/lib")
      set(BLAS_LIBRARIES "${BLAS_INCLUDE_DIRS}/libblas.a")
      set(LAPACK_INCLUDE_DIRS "/glade/apps/opt/lib")
      set(LAPACK_LIBRARIES "${LAPACK_INCLUDE_DIRS}/liblapack.a")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${BLAS_LIBRARIES}")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LIBRARIES}")
      set(NEED_LAPACK FALSE)
    endif()

    if (TRUE) # intel
      set(NEED_LAPACK FALSE)
    endif()
    
  endif()

endmacro()
