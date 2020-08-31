set(IO_SRC
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

# Create io lib
add_library(io ${IO_SRC})

# Set includes
target_include_directories(io PUBLIC ${SCREAM_INCLUDE_DIRS})

# Set f90 compiled modules dir
set_target_properties(io PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/cmake/atm)

# Link libraries
target_link_libraries(io PRIVATE ${SCREAM_LIBS} ${SCREAM_CIME_LIBS})

add_subdirectory(tests)
