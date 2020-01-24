set(ATM_SRC
  scream_cxx_f90_interface.cpp
  atm_comp_mct.F90
  dead_mct_mod.F90
  dead_data_mod.F90
  dead_mod.F90
)

if (USE_ESMF_LIB)
  set(ESMFDIR "esmf")
else()
  set(ESMFDIR "noesmf")
endif()

set(CSM_SHARE "${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/${NINST_VALUE}/csm_share")
link_directories(${SCREAM_TPL_LIBRARY_DIRS})

add_library(atm ${ATM_SRC})
target_include_directories(atm PUBLIC ${SCREAM_INCLUDE_DIRS} ${SCREAM_TPL_INCLUDE_DIRS} ${CSM_SHARE} ${INSTALL_SHAREDPATH}/include)
set_target_properties(atm PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/cmake/atm)
target_link_libraries(atm scream_share scream_control p3 shoc ${SCREAM_TPL_LIBRARIES})
