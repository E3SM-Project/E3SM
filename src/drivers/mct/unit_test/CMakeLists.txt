set(DRV_ROOT "${CIME_ROOT}/src/drivers/mct")

add_definitions(
  -DNUM_COMP_INST_ATM=1
  -DNUM_COMP_INST_LND=1
  -DNUM_COMP_INST_OCN=1
  -DNUM_COMP_INST_ICE=1
  -DNUM_COMP_INST_GLC=1
  -DNUM_COMP_INST_WAV=1
  -DNUM_COMP_INST_ROF=1
  -DNUM_COMP_INST_ESP=1
  -DNUM_COMP_INST_IAC=1
  )

# The following definitions are needed when building with the mpi-serial library
if (USE_MPI_SERIAL)
  add_definitions(-DNO_MPI2 -DNO_MPIMOD)
endif()

# Add source directories from stubs. This should be done first, so that in the
# case of name collisions, the drv versions take precedence (when there are two
# files with the same name, the one added later wins).
add_subdirectory(${CIME_ROOT}/src/share/unit_test_stubs/pio pio)

# Add drv source directories
add_subdirectory(${DRV_ROOT}/shr drv_shr)
add_subdirectory(${DRV_ROOT}/main drv_main)

# Add general unit test directories (stubbed out files, etc.)
add_subdirectory(utils drv_unit_test_utils)
add_subdirectory(stubs drv_unit_test_stubs)

# Build libraries containing stuff needed for the unit tests.
#
# Eventually, these add_library calls should probably be distributed into the
# correct location, rather than being in this top-level CMakeLists.txt file.
#
# Note that we are including the stub pio in the csm_share library for simplicity.
#
# We include the esmf_wrf_timemgr_sources in the csm_share library
# rather than building it into its own library because of circular
# dependencies: shr_cal_mod (in share_sources) depends on esmf, but esmf
# depends on shr_sys_mod (also in share_sources).
add_library(csm_share ${share_sources} ${share_mct_sources} ${share_pio_sources}
  ${pio_sources} ${esmf_wrf_timemgr_sources})
declare_generated_dependencies(csm_share "${share_genf90_sources}")
declare_generated_dependencies(csm_share "${pio_genf90_sources}")
add_dependencies(csm_share mct_project)

add_library(drv ${drv_sources})
add_dependencies(drv csm_share)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Set the list of libraries needed for these unit tests. Note that not all unit
# tests need all of these libraries, but it's easiest just to set the same list
# for everyone.
set(DRV_UNIT_TEST_LIBS drv;csm_share;mct;mpeu)
if (USE_MPI_SERIAL)
  list(APPEND DRV_UNIT_TEST_LIBS mpi-serial)
endif()
list(APPEND DRV_UNIT_TEST_LIBS ${NETCDF_LIBRARIES})

# Add the test directories
add_subdirectory(avect_wrapper_test)
add_subdirectory(seq_map_test)
add_subdirectory(glc_elevclass_test)
add_subdirectory(map_glc2lnd_test)
add_subdirectory(map_lnd2rof_irrig_test)
add_subdirectory(check_fields_test)
