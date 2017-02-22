set(DRV_ROOT "${CIME_ROOT}/src/driver_cpl")

add_definitions(
  -DNUM_COMP_INST_ATM=1
  -DNUM_COMP_INST_LND=1
  -DNUM_COMP_INST_OCN=1
  -DNUM_COMP_INST_ICE=1
  -DNUM_COMP_INST_GLC=1
  -DNUM_COMP_INST_WAV=1
  -DNUM_COMP_INST_ROF=1
  -DNUM_COMP_INST_ESP=1
  )

# Add source directories from stubs. This should be done first, so that in the
# case of name collisions, the drv versions take precedence (when there are two
# files with the same name, the one added later wins).
add_subdirectory(${CIME_ROOT}/share/csm_share/unit_test_stubs/pio pio)

# Add drv source directories
add_subdirectory(${DRV_ROOT}/shr drv_shr)
add_subdirectory(${DRV_ROOT}/driver drv_driver)

# Add general unit test directories (stubbed out files, etc.)
add_subdirectory(utils drv_unit_test_utils)
add_subdirectory(stubs drv_unit_test_stubs)

# Build libraries containing stuff needed for the unit tests. 
# Eventually, these add_library calls should probably be distributed into the
# correct location, rather than being in this top-level CMakeLists.txt file.
# Note that we are including the stub pio in the csm_share library for simplicity.
add_library(csm_share ${share_sources} ${share_mct_sources} ${share_pio_sources}
  ${pio_sources})
declare_generated_dependencies(csm_share "${share_genf90_sources}")
declare_generated_dependencies(csm_share "${pio_genf90_sources}")
add_dependencies(csm_share mct_project)

add_library(esmf_wrf_timemgr ${esmf_wrf_timemgr_sources})
add_dependencies(esmf_wrf_timemgr csm_share)

add_library(drv ${drv_sources})
add_dependencies(drv csm_share esmf_wrf_timemgr)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Set the list of libraries needed for these unit tests. Note that not all unit
# tests need all of these libraries, but it's easiest just to set the same list
# for everyone.
set(DRV_UNIT_TEST_LIBS drv;csm_share;esmf_wrf_timemgr;mct;mpeu)
list(APPEND DRV_UNIT_TEST_LIBS ${NETCDF_LIBRARIES})

# Add the test directories
add_subdirectory(avect_wrapper_test)
add_subdirectory(seq_map_test)
add_subdirectory(glc_elevclass_test)
add_subdirectory(vertical_gradient_calculator_test)
add_subdirectory(map_lnd2glc_test)
add_subdirectory(map_glc2lnd_test)
add_subdirectory(map_lnd2rof_irrig_test)
