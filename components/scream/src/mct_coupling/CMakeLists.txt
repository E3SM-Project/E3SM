set (dynLibName)
if ("${SCREAM_DYNAMICS_DYCORE}" STREQUAL "HOMME")
  # Need these two to be linked, due to homme deps
  # Note: homme finds these packages, but it was configured
  #       outside this folder (and not in a parent folder),
  #       so the results of find_package have been lost
  #       as soon as we got out of that scope
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)

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
  set (SCREAM_USE_ENERGY FALSE CACHE STRING "Whether Homme has extra energy checks.")

  add_definitions(-DHOMME_WITHOUT_PIOLIBRARY)
  CreateDynamicsLib (${SCREAM_DYN_TARGET} ${SCREAM_NP} ${SCREAM_NUM_VERTICAL_LEV} ${SCREAM_NUM_TRACERS} ${SCREAM_USE_PIO} ${SCREAM_USE_ENERGY})

  target_include_directories (
    ${dynLibName}
    PUBLIC
    ${INSTALL_SHAREDPATH}/pio/pio${PIO_VERSION}/pio
    ${INSTALL_SHAREDPATH}/mct/noesmf/${NINST_VALUE}/csm_share
    ${NetCDF_Fortran_PATH}/include
  )

  target_link_libraries(${dynLibName} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
  # target_link_directories(${dynLibName} PUBLIC ${NetCDF_Fortran_PATH} ${INSTALL_SHAREDPATH}/pio/pio${PIO_VERSION}/pio)
  # target_link_libraries(${dynLibName} pio netcdff)
endif()

set(ATM_SRC
  scream_cxx_f90_interface.cpp
  surface_coupling.cpp
  atm_comp_mct.F90
  dead_mct_mod.F90
  dead_data_mod.F90
  dead_mod.F90
  scream_scorpio_interface.F90
  scream_scorpio_interface.cpp
  scream_scorpio_interface_iso_c2f.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/cam/physics_utils.F90
)

if (USE_ESMF_LIB)
  set(ESMFDIR "esmf")
else()
  set(ESMFDIR "noesmf")
endif()
set(CSM_SHARE "${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/${NINST_VALUE}/csm_share")

set (SCREAM_LIBS
     scream_share
     scream_control
     ${dynLibName}
     p3
     shoc
)

# In CIME builds, create imported targets for scorpio and csm_share
include (Scorpio)
CreateScorpioTarget(TRUE) # Note: the input is whether we want FLIB too. We want/need both, so look for both
include (CsmShare)
CreateCsmShareTarget()

set (SCREAM_CIME_LIBS
     scream_pioc
     scream_piof
     scream_csm_share
)

# Create atm lib
add_library(atm ${ATM_SRC})

# Set includes
target_include_directories(atm PUBLIC ${SCREAM_INCLUDE_DIRS})

# Set f90 compiled modules dir
set_target_properties(atm PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/cmake/atm)

# Link libraries
target_link_libraries(atm PRIVATE ${SCREAM_LIBS} ${SCREAM_CIME_LIBS})

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
