# All the executables used by the Homme tests (and by the cases in test/cases).
#
# Each executable is built in <build-dir>/test/execs/<exec-name>. Since multiple tests can use the
# same executable, the executables are created here, rather than where the tests are.
# Use homme_add_model_exec (see cmake/HommeMacros.cmake) to add a new one.

# For quick manual runs (development/debugging only): configure a namelist (and copy the vcoord files
# it needs) in the folder of an executable. The KEY=VALUE pairs set the variables used in the namelist.
function (homme_add_exec_namelist exec namelist)
  foreach (kv IN LISTS ARGN)
    string (REGEX REPLACE "=.*" "" key "${kv}")
    string (REGEX REPLACE "^[^=]*=" "" val "${kv}")
    set (${key} "${val}")
  endforeach()
  set (dir ${CMAKE_CURRENT_BINARY_DIR}/${exec})
  configure_file (${namelist} ${dir}/namelist.nl)
  file (MAKE_DIRECTORY ${dir}/movies ${dir}/vcoord)
  file (COPY ${HOMME_TEST_DATA_DIR}/vcoord/${HOMME_TEST_VCOORD_INT_FILE}
             ${HOMME_TEST_DATA_DIR}/vcoord/${HOMME_TEST_VCOORD_MID_FILE}
        DESTINATION ${dir}/vcoord)
endfunction()

# Some executables that no test uses (e.g., with many levels). They can be useful for manual runs.
option (HOMME_BUILD_EXTRA_TEST_EXECS "Build additional executables that are not used by any test" OFF)

if (BUILD_HOMME_SWEQX)
  homme_add_model_exec (swtcA MODEL sweqx NP 8 NC 4 PLEV 4 QSIZE_D 4 OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/swtcA)
  homme_add_model_exec (swtcB MODEL sweqx PLEV 1 QSIZE_D 4 OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/swtcB)
  homme_add_model_exec (swtcC MODEL sweqx PLEV 1 QSIZE_D 1 PIO OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/swtcC)
endif()

if (BUILD_HOMME_PREQX)
  homme_add_model_exec (baroC MODEL preqx PLEV 26 QSIZE_D 25 ENERGY OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/baroC)
  homme_add_model_exec (baroCam MODEL preqx PLEV 26 QSIZE_D 25 PIO ENERGY OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/baroCam)
  homme_add_model_exec (preqx-nlev10-native MODEL preqx PLEV 10 QSIZE_D 1 PIO OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/preqx-nlev10-native)
  homme_add_model_exec (preqx-nlev20-interp MODEL preqx PLEV 20 QSIZE_D 1 OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/preqx-nlev20-interp)
  homme_add_model_exec (preqx-nlev20-native MODEL preqx PLEV 20 QSIZE_D 1 PIO OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/preqx-nlev20-native)
  homme_add_model_exec (preqx-nlev30-interp MODEL preqx PLEV 30 QSIZE_D 6 ENERGY OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/preqx-nlev30-interp)
  homme_add_model_exec (preqx-nlev30-native MODEL preqx PLEV 30 QSIZE_D 6 PIO OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/preqx-nlev30-native)
  homme_add_model_exec (preqx-nlev40-interp MODEL preqx PLEV 40 QSIZE_D 4 OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/preqx-nlev40-interp)
  homme_add_model_exec (preqx-nlev60-interp MODEL preqx PLEV 60 QSIZE_D 4 OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/preqx-nlev60-interp)
  homme_add_model_exec (preqx-nlev60-native MODEL preqx PLEV 60 QSIZE_D 4 PIO OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/preqx-nlev60-native)
  if (HOMME_BUILD_EXTRA_TEST_EXECS)
    homme_add_model_exec (preqx-nlev100-native MODEL preqx PLEV 100 QSIZE_D 1 PIO OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/preqx-nlev100-native)
    homme_add_model_exec (preqx-nlev140-native MODEL preqx PLEV 140 QSIZE_D 1 PIO OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/preqx-nlev140-native)
    homme_add_model_exec (preqx-nlev200-native MODEL preqx PLEV 200 QSIZE_D 1 PIO OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/preqx-nlev200-native)
  endif()
endif()

if (BUILD_HOMME_TOOL)
  homme_add_model_exec (tool-nlev26 MODEL tool PLEV 26 QSIZE_D 4 PIO OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/tool-nlev26)
endif()

if (PREQX_USE_PIO)
  set (theta_nlev128_pio PIO)
endif()

if (BUILD_HOMME_THETA)
  homme_add_model_exec (theta-l-nlev20 MODEL theta-l PLEV 20 QSIZE_D 1 ENERGY OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/theta-l-nlev20)
  homme_add_model_exec (theta-l-nlev30 MODEL theta-l PLEV 30 QSIZE_D 6 ENERGY OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/theta-l-nlev30)
  homme_add_model_exec (theta-l-nlev40 MODEL theta-l PLEV 40 QSIZE_D 6 ENERGY OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/theta-l-nlev40)
  homme_add_model_exec (theta-l-nlev60 MODEL theta-l PLEV 60 QSIZE_D 4 ENERGY OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/theta-l-nlev60)
  homme_add_model_exec (theta-l-nlev10-native MODEL theta-l PLEV 10 QSIZE_D 3 PIO ENERGY OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/theta-l-nlev10-native)
  homme_add_model_exec (theta-l-nlev20-native MODEL theta-l PLEV 20 QSIZE_D 3 PIO ENERGY OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/theta-l-nlev20-native)
  homme_add_model_exec (theta-l-nlev128-native MODEL theta-l PLEV 128 QSIZE_D 3 PIO ENERGY OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/theta-l-nlev128-native)
  # Special test for -DHOMME_WITHOUT_PIOLIBRARY (build only)
  homme_add_model_exec (stt MODEL theta-l PLEV 3 QSIZE_D 5 ENERGY COMPILE_DEFINITIONS HOMME_WITHOUT_PIOLIBRARY OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/stt)
  if (HOMME_BUILD_EXTRA_TEST_EXECS)
    homme_add_model_exec (theta-l-nlev30-native MODEL theta-l PLEV 30 QSIZE_D 3 PIO ENERGY OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/theta-l-nlev30-native)
    homme_add_model_exec (theta-l-nlev60-native MODEL theta-l PLEV 60 QSIZE_D 3 PIO ENERGY OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/theta-l-nlev60-native)
    homme_add_model_exec (theta-l-nlev64-native MODEL theta-l PLEV 64 QSIZE_D 3 PIO ENERGY OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/theta-l-nlev64-native)
    homme_add_model_exec (theta-l-nlev100-native MODEL theta-l PLEV 100 QSIZE_D 3 PIO ENERGY OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/theta-l-nlev100-native)
    homme_add_model_exec (theta-l-nlev140-native MODEL theta-l PLEV 140 QSIZE_D 3 PIO ENERGY OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/theta-l-nlev140-native)
    homme_add_model_exec (theta-l-nlev150-native MODEL theta-l PLEV 150 QSIZE_D 3 PIO ENERGY OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/theta-l-nlev150-native)
    homme_add_model_exec (theta-l-nlev200-native MODEL theta-l PLEV 200 QSIZE_D 3 PIO ENERGY OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/theta-l-nlev200-native)
    homme_add_model_exec (theta-l-nlev256-native MODEL theta-l PLEV 256 QSIZE_D 3 PIO ENERGY OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/theta-l-nlev256-native)
    homme_add_model_exec (theta-l-nlev300-native MODEL theta-l PLEV 300 QSIZE_D 3 PIO ENERGY OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/theta-l-nlev300-native)
    homme_add_model_exec (theta-nlev72 MODEL theta-l PLEV 72 QSIZE_D 10 OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/theta-nlev72)
    homme_add_exec_namelist (theta-nlev72 ${HOMME_SYSTEM_TESTS_SRC_DIR}/theta-l/baroclinic/theta.nl
      HOMME_TEST_VCOORD_INT_FILE=acme-72i.ascii
      HOMME_TEST_VCOORD_MID_FILE=acme-72m.ascii
      HOMME_TEST_NE=8
      HOMME_TEST_TIME_STEP=450
      HOMME_TEST_NDAYS=1
      HOMME_TEST_QSIZE=4
      HOMME_TEST_RSPLIT=2
      HOMME_TEST_NU=7e15
      HOMME_TEST_NUDIV=1e15
      HOMME_TEST_NUTOP=2.5e5
      HOMME_TEST_LIM=9
      HOMME_TEST_HVSCALING=2
      HOMME_TEST_HVS_TOM=0
      HOMME_THETA_HY_MODE=true
      HOMME_THETA_FORM=1
      HOMME_TTYPE=5)
  endif()
endif()

if (BUILD_HOMME_THETA_KOKKOS)
  homme_add_model_exec (theta-l-nlev30-kokkos MODEL theta-l_kokkos PLEV 30 QSIZE_D 6 ENERGY OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/theta-l-nlev30-kokkos)
  homme_add_exec_namelist (theta-l-nlev30-kokkos ${CMAKE_CURRENT_SOURCE_DIR}/thetanh.nl
    HOMME_TEST_VCOORD_INT_FILE=cami-30.ascii
    HOMME_TEST_VCOORD_MID_FILE=camm-30.ascii
    HOMME_TEST_NE=8
    HOMME_TEST_TIME_STEP=450
    HOMME_TEST_NDAYS=1
    HOMME_TEST_QSIZE=4
    HOMME_TEST_RSPLIT=2
    HOMME_TEST_QSPLIT=2
    HOMME_TEST_MOISTURE=dry
    HOMME_TEST_HYDROSTATIC_MODE=true)
  homme_add_model_exec (theta-l-nlev72-kokkos MODEL theta-l_kokkos PLEV 72 QSIZE_D 10 ENERGY OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/theta-l-nlev72-kokkos)
  homme_add_exec_namelist (theta-l-nlev72-kokkos ${CMAKE_CURRENT_SOURCE_DIR}/thetanh.nl
    HOMME_TEST_VCOORD_INT_FILE=acme-72i.ascii
    HOMME_TEST_VCOORD_MID_FILE=acme-72m.ascii
    HOMME_TEST_NE=8
    HOMME_TEST_TIME_STEP=450
    HOMME_TEST_NDAYS=1
    HOMME_TEST_QSIZE=10
    HOMME_TEST_RSPLIT=2
    HOMME_TEST_QSPLIT=2
    HOMME_TEST_MOISTURE=dry
    HOMME_TEST_HYDROSTATIC_MODE=true)
  homme_add_model_exec (theta-nlev128-kokkos MODEL theta-l_kokkos PLEV 128 QSIZE_D 10 ${theta_nlev128_pio} OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/theta-nlev128-kokkos)
  homme_add_exec_namelist (theta-nlev128-kokkos ${HOMME_SYSTEM_TESTS_SRC_DIR}/theta-l/baroclinic/theta.nl
    HOMME_TEST_VCOORD_INT_FILE=sabi-128.ascii
    HOMME_TEST_VCOORD_MID_FILE=sabm-128.ascii
    HOMME_TEST_NE=8
    HOMME_TEST_TIME_STEP=450
    HOMME_TEST_NDAYS=1
    HOMME_TEST_HVSCALING=0
    HOMME_TEST_QSIZE=4
    HOMME_TEST_RSPLIT=2
    HOMME_TEST_NU=7e15
    HOMME_TEST_NUDIV=1e15
    HOMME_TEST_NUTOP=2.5e5
    HOMME_TEST_LIM=9
    HOMME_TEST_HVS=1
    HOMME_TEST_HVS_TOM=0
    HOMME_THETA_HY_MODE=true
    HOMME_THETA_FORM=1
    HOMME_TTYPE=5)
  # For doubly periodic tests
  homme_add_model_exec (theta-l-nlev20-native-kokkos MODEL theta-l_kokkos PLEV 20 QSIZE_D 3 PIO OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/theta-l-nlev20-native-kokkos)
  homme_add_model_exec (theta-l-nlev10-native-kokkos MODEL theta-l_kokkos PLEV 10 QSIZE_D 3 PIO OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/theta-l-nlev10-native-kokkos)
  homme_add_model_exec (theta-l-nlev128-native-kokkos MODEL theta-l_kokkos PLEV 128 QSIZE_D 3 PIO OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/theta-l-nlev128-native-kokkos)
  # Only needed for perf/correctness comparisons against the C++ version
  if (BUILD_HOMME_THETA)
    homme_add_model_exec (theta-nlev128 MODEL theta-l PLEV 128 QSIZE_D 10 ${theta_nlev128_pio} OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/theta-nlev128)
    homme_add_exec_namelist (theta-nlev128 ${HOMME_SYSTEM_TESTS_SRC_DIR}/theta-l/baroclinic/theta.nl
      HOMME_TEST_VCOORD_INT_FILE=sabi-128.ascii
      HOMME_TEST_VCOORD_MID_FILE=sabm-128.ascii
      HOMME_TEST_NE=8
      HOMME_TEST_TIME_STEP=450
      HOMME_TEST_NDAYS=1
      HOMME_TEST_HVSCALING=0
      HOMME_TEST_QSIZE=4
      HOMME_TEST_RSPLIT=2
      HOMME_TEST_NU=7e15
      HOMME_TEST_NUDIV=1e15
      HOMME_TEST_NUTOP=2.5e5
      HOMME_TEST_LIM=9
      HOMME_TEST_HVS=1
      HOMME_TEST_HVS_TOM=0
      HOMME_THETA_HY_MODE=true
      HOMME_THETA_FORM=1
      HOMME_TTYPE=5)
  endif()
  if (HOMME_BUILD_EXTRA_TEST_EXECS AND BUILD_HOMME_THETA)
    homme_add_model_exec (theta-nlev72-kokkos MODEL theta-l_kokkos PLEV 72 QSIZE_D 10 OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/theta-nlev72-kokkos)
    homme_add_exec_namelist (theta-nlev72-kokkos ${HOMME_SYSTEM_TESTS_SRC_DIR}/theta-l/baroclinic/theta.nl
      HOMME_TEST_VCOORD_INT_FILE=acme-72i.ascii
      HOMME_TEST_VCOORD_MID_FILE=acme-72m.ascii
      HOMME_TEST_NE=8
      HOMME_TEST_TIME_STEP=450
      HOMME_TEST_NDAYS=1
      HOMME_TEST_QSIZE=4
      HOMME_TEST_RSPLIT=2
      HOMME_TEST_NU=7e15
      HOMME_TEST_NUDIV=1e15
      HOMME_TEST_NUTOP=2.5e5
      HOMME_TEST_LIM=9
      HOMME_TEST_HVSCALING=2
      HOMME_TEST_HVS_TOM=0
      HOMME_THETA_HY_MODE=true
      HOMME_THETA_FORM=1
      HOMME_TTYPE=5)
  endif()
endif()

if (BUILD_HOMME_PREQX_KOKKOS)
  homme_add_model_exec (preqx-nlev26-kokkos MODEL preqx_kokkos PLEV 26 QSIZE_D 4 OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/preqx-nlev26-kokkos)
  homme_add_exec_namelist (preqx-nlev26-kokkos ${HOMME_SYSTEM_TESTS_SRC_DIR}/preqx/preqx.nl
    HOMME_TEST_VCOORD_INT_FILE=cami-26.ascii
    HOMME_TEST_VCOORD_MID_FILE=camm-26.ascii
    HOMME_TEST_NE=2
    HOMME_TEST_TIME_STEP=300
    HOMME_TEST_NDAYS=1
    HOMME_TEST_QSIZE=4
    HOMME_TEST_RSPLIT=3
    HOMME_TEST_MOISTURE=notdry
    HOMME_TEST_TIME_STEP=600)
  homme_add_model_exec (preqx-nlev72-kokkos MODEL preqx_kokkos PLEV 72 QSIZE_D 40 OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/preqx-nlev72-kokkos)
  homme_add_exec_namelist (preqx-nlev72-kokkos ${HOMME_SYSTEM_TESTS_SRC_DIR}/preqx/preqx.nl
    HOMME_TEST_VCOORD_INT_FILE=acme-72i.ascii
    HOMME_TEST_VCOORD_MID_FILE=acme-72m.ascii
    HOMME_TEST_NE=2
    HOMME_TEST_TIME_STEP=300
    HOMME_TEST_NDAYS=1
    HOMME_TEST_QSIZE=4
    HOMME_TEST_RSPLIT=3
    HOMME_TEST_MOISTURE=notdry
    HOMME_TEST_TIME_STEP=600
    HOMME_TEST_HVSCALING=0
    HOMME_TEST_NU=7e15
    HOMME_TEST_NUDIV=1e15
    HOMME_TEST_NUTOP=0
    HOMME_TEST_LIM=8)
  if (HOMMEXX_BFB_TESTING)
    # F90 versions, to compare against the C++ ones
    homme_add_model_exec (preqx-nlev26 MODEL preqx PLEV 26 QSIZE_D 4 OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/preqx-nlev26)
    homme_add_exec_namelist (preqx-nlev26 ${HOMME_SYSTEM_TESTS_SRC_DIR}/preqx/preqx.nl
      HOMME_TEST_VCOORD_INT_FILE=cami-26.ascii
      HOMME_TEST_VCOORD_MID_FILE=camm-26.ascii
      HOMME_TEST_NE=2
      HOMME_TEST_TIME_STEP=300
      HOMME_TEST_NDAYS=1
      HOMME_TEST_QSIZE=4
      HOMME_TEST_RSPLIT=3
      HOMME_TEST_MOISTURE=notdry)
    homme_add_model_exec (preqx-nlev72 MODEL preqx PLEV 72 QSIZE_D 40 OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/preqx-nlev72)
    homme_add_exec_namelist (preqx-nlev72 ${HOMME_SYSTEM_TESTS_SRC_DIR}/preqx/preqx.nl
      HOMME_TEST_VCOORD_INT_FILE=acme-72i.ascii
      HOMME_TEST_VCOORD_MID_FILE=acme-72m.ascii
      HOMME_TEST_NE=2
      HOMME_TEST_TIME_STEP=300
      HOMME_TEST_NDAYS=1
      HOMME_TEST_QSIZE=4
      HOMME_TEST_RSPLIT=3
      HOMME_TEST_MOISTURE=notdry)
  endif()
  # For doubly periodic tests
  homme_add_model_exec (preqx-nlev10-native-kokkos MODEL preqx_kokkos PLEV 10 QSIZE_D 3 PIO OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/preqx-nlev10-native-kokkos)
endif()

# Executables used by the cases in test/cases (not by the ctest tests)
add_custom_target (cases-execs)
foreach (exec IN ITEMS preqx-nlev20-interp preqx-nlev30-interp preqx-nlev60-interp
                       theta-l-nlev20 theta-l-nlev30 theta-l-nlev40 theta-l-nlev60)
  if (TARGET ${exec})
    add_dependencies (cases-execs ${exec})
  endif()
endforeach()
