include (ScreamUtils)

# Standalone AD runs of the water_tracers process on a mesh-free point grid.
# All inputs are constant ICs (no input files). The ocean fluxes are
# diagnostic only, so tracers and qv stay at their initial values.

set (NUM_STEPS 2)
set (ATM_TIME_STEP 1800)
set (RUN_T0 2021-10-12-45000)

CreateADUnitTestExec(water_tracers_standalone LIBS water_tracers)

# Generic and property tracers only.
set (POSTFIX generic)
set (TRACERS_FILE tracers_generic.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input_${POSTFIX}.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output_${POSTFIX}.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tracers_generic.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/tracers_generic.yaml COPYONLY)
CreateUnitTestFromExec(water_tracers_standalone water_tracers_standalone
  EXE_ARGS "--args -ifile=input_${POSTFIX}.yaml"
  LABELS water_tracers physics)

if (EAMXX_ENABLE_WATER_ISOTOPES)
  # Same tracers, plus isotope entries in the same file/process.
  set (POSTFIX isotopes)
  set (TRACERS_FILE tracers_with_isotopes.yaml)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_${POSTFIX}.yaml)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/output_${POSTFIX}.yaml)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tracers_with_isotopes.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/tracers_with_isotopes.yaml COPYONLY)
  CreateUnitTestFromExec(water_tracers_standalone_isotopes water_tracers_standalone
    EXE_ARGS "--args -ifile=input_${POSTFIX}.yaml"
    LABELS water_tracers water_isotopes physics)
else()
  # Negative test: a build without isotope support must reject isotope
  # entries with a clear build-time-aware error, not silently ignore them
  # or require a different process name.
  set (POSTFIX isotopes_rejected)
  set (TRACERS_FILE tracers_with_isotopes.yaml)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_${POSTFIX}.yaml)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/output_${POSTFIX}.yaml)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tracers_with_isotopes.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/tracers_with_isotopes.yaml COPYONLY)
  CreateUnitTestFromExec(water_tracers_standalone_isotopes_rejected water_tracers_standalone
    EXE_ARGS "--args -ifile=input_${POSTFIX}.yaml"
    LABELS water_tracers physics)
  set_tests_properties(water_tracers_standalone_isotopes_rejected PROPERTIES
    WILL_FAIL TRUE)
endif()
