include (ScreamUtils)

# Standalone AD runs of the water tracer processes on a mesh-free point grid.
# All inputs are constant ICs (no input files). The ocean fluxes are
# diagnostic only, so tracers and qv stay at their initial values.

set (NUM_STEPS 2)
set (ATM_TIME_STEP 1800)
set (RUN_T0 2021-10-12-45000)

set (WTRC_LIBS water_tracers)
if (EAMXX_ENABLE_WATER_ISOTOPES)
  list (APPEND WTRC_LIBS water_isotopes)
endif()
CreateADUnitTestExec(water_tracers_standalone LIBS ${WTRC_LIBS})

# Generic and property tracers (water_tracers process)
set (POSTFIX generic)
set (WTRC_PROC water_tracers)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input_${POSTFIX}.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output_${POSTFIX}.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/water_tracers.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/water_tracers_${POSTFIX}.yaml COPYONLY)
CreateUnitTestFromExec(water_tracers_standalone water_tracers_standalone
  EXE_ARGS "--args -ifile=input_${POSTFIX}.yaml"
  LABELS water_tracers physics)

# Same tracers plus isotopes (water_isotopes process)
if (EAMXX_ENABLE_WATER_ISOTOPES)
  set (POSTFIX isotopes)
  set (WTRC_PROC water_isotopes)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/input_${POSTFIX}.yaml)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/output_${POSTFIX}.yaml)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/water_isotopes.yaml
                 ${CMAKE_CURRENT_BINARY_DIR}/water_tracers_${POSTFIX}.yaml COPYONLY)
  CreateUnitTestFromExec(water_tracers_standalone_isotopes water_tracers_standalone
    EXE_ARGS "--args -ifile=input_${POSTFIX}.yaml"
    LABELS water_tracers water_isotopes physics)
endif()
