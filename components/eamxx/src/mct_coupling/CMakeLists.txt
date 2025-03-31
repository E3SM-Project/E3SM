set (dynLibName)
if ("${SCREAM_DYNAMICS_DYCORE}" STREQUAL "HOMME")
  # Need these two to be linked, due to homme deps

  # Recall that src/dynamics/homme/CMakeLists.txt does not build a dyn lib.
  # It only defines a macro to build it, so that different tests can build
  # a dyn lib if necessary.
  # Here, we ask to create a dyn lib depending on the config options that
  # were set during e3sm configuration

  # This is a list of cmake vars whose values are then used when calling
  # CreateDynamicsLib, to correctly build the dynamics library within Homme.
  # We set them to a default, but each compset should set its values anyways.
  set (SCREAM_DYN_TARGET none CACHE STRING "The name of the desired Homme target.")
  set (SCREAM_NP 4 CACHE STRING "The number of Gauss points per element.")
  set (SCREAM_NUM_TRACERS 4 CACHE STRING "The max number of tracers.")
  set (SCREAM_USE_PIO FALSE CACHE STRING "Whether Homme can use PIO.")

  CreateDynamicsLib (${SCREAM_DYN_TARGET} ${SCREAM_NP} ${SCREAM_NUM_VERTICAL_LEV} ${SCREAM_NUM_TRACERS})
endif()

set(ATM_SRC
  eamxx_cxx_f90_interface.cpp
  atm_comp_mct.F90
  ekat_string_utils.F90
  eamxx_f2c_mod.F90
  eamxx_cpl_indices.F90
  ${SCREAM_BASE_DIR}/../eam/src/physics/cam/physics_utils.F90
)

set (SCREAM_LIBS
     scream_share
     scream_control
     ${dynLibName}
     p3
     shoc
     scream_rrtmgp
     eamxx_cosp
     cld_fraction
     spa
     iop_forcing
     nudging
     diagnostics
     tms
)
if (SCREAM_ENABLE_ML_CORRECTION)
  list (APPEND SCREAM_LIBS ml_correction)
endif ()
if (SCREAM_ENABLE_MAM)
  list (APPEND SCREAM_LIBS mam)
endif()

# Create atm lib
add_library(atm ${ATM_SRC})

# Link libraries
# NOTE: when ACE is enable, we can safely discard most of the model
if (EAMXX_ENABLE_ACE)
  # Set the default PyTorch paths if not already defined
  if(NOT DEFINED TORCH_CMAKE_PATH)
    if (Kokkos_ENABLE_CUDA)
      set(TORCH_CMAKE_PATH "/global/cfs/cdirs/e3sm/software/libtorch/libtorch-cxx11-abi-shared-with-deps-2.6.0+cu124" 
          CACHE PATH "Path to PyTorch CMake files")
      message(STATUS "MCT Using GPU version of PyTorch at ${TORCH_CMAKE_PATH}")
    else()
      set(TORCH_CMAKE_PATH "/global/cfs/cdirs/e3sm/software/libtorch/libtorch-cxx11-abi-shared-with-deps-2.6.0+cpu" 
          CACHE PATH "Path to PyTorch CMake files")
      message(STATUS "MCT Using CPU version of PyTorch at ${TORCH_CMAKE_PATH}")
    endif()
  endif()

  # Set CMake prefix path for PyTorch
  set(CMAKE_PREFIX_PATH ${TORCH_CMAKE_PATH})

  # Find PyTorch and require it
  find_package(Torch REQUIRED)

  if(Torch_FOUND)
    message(STATUS "MCT Found PyTorch: Building ACE component")
    set(SCREAM_FOUND_PYTORCH TRUE CACHE BOOL "Whether PyTorch was found" FORCE)

    # is this part needed with the torch includes?
    set_target_properties(atm PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/cmake/atm ${TORCH_INCLUDE_DIRS})

    # Explicitly add PyTorch include directories
    target_include_directories(atm PUBLIC ${TORCH_INCLUDE_DIRS})

    # Explicitly link PyTorch libraries
    # TODO: not sure why this must be PRIVATE?
    # TODO: Even though I don't include shoc, p3, etc., they still get built? Ugh!
    target_link_libraries(atm PRIVATE scream_share scream_control ${dynLibName} ace ${TORCH_LIBRARIES} diagnostics csm_share)
  else()
    message(FATAL_ERROR "PyTorch not found. Ensure TORCH_CMAKE_PATH is set correctly.")
  endif()
else()
  # Business as usual
  # Set f90 compiled modules dir
  set_target_properties(atm PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/cmake/atm)
  target_include_directories(atm PUBLIC ${CMAKE_BINARY_DIR}/cmake/atm)
  target_link_libraries(atm PRIVATE ${SCREAM_LIBS} csm_share)
endif()
