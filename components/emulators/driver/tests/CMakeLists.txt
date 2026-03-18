cmake_minimum_required(VERSION 3.16)

enable_testing()
find_package(MPI REQUIRED COMPONENTS Fortran)

set(DRIVER_TEST_SOURCES
  test_emulator_factory.F90
)

add_executable(test_emulator_factory ${DRIVER_TEST_SOURCES})

# This test needs:
# - Fortran modules from emulator_common
# - the factory implementation from emulator_driver
# - MPI Fortran (if you use MPI in the test)
target_link_libraries(test_emulator_factory PRIVATE
  emulator_common
  emulator_driver
  MPI::MPI_Fortran
)

# Make sure it can find the Fortran modules from emulator_common
# (assuming emulator_common set Fortran_MODULE_DIRECTORY and added it to its include dirs)
# Linking to emulator_common is usually enough, but you can be explicit:
get_target_property(EMU_MOD_DIR emulator_common Fortran_MODULE_DIRECTORY)
if (EMU_MOD_DIR)
  target_include_directories(test_emulator_factory PRIVATE ${EMU_MOD_DIR})
endif()

add_test(NAME emulator_factory_basic
         COMMAND test_emulator_factory)

# On Cray systems GPU-aware MPI is enabled by default but requires GTL to
# be linked. Disable it for unit tests which don't need GPU MPI transfers.
set_tests_properties(emulator_factory_basic PROPERTIES
  ENVIRONMENT "MPICH_GPU_SUPPORT_ENABLED=0")
