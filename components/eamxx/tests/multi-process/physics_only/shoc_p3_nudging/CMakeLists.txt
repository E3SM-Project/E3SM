INCLUDE (ScreamUtils)

# Create the exec
CreateADUnitTestExec (shoc_p3_nudging
  LIBS shoc p3 nudging diagnostics)

# Ensure test input files are present in the data dir
GetInputFile(scream/init/${EAMxx_tests_IC_FILE_72lev})
GetInputFile(scream/init/${EAMxx_tests_IC_FILE_128lev})
GetInputFile(cam/topo/${EAMxx_tests_TOPO_FILE})

set (RUN_T0 2021-10-12-45000)

# Run a quick setup function to create vertical remap file:
CreateUnitTest(create_vert_remap_and_weights "create_vert_remap_and_weights.cpp"
  LIBS scream_share
  FIXTURES_SETUP shoc_p3_create_vertical_remap_and_weights_file)

# Run a test to setup nudging source data:
set (NUM_STEPS  5)
set (POSTFIX source_data)
set (ATM_TIME_STEP 300)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input_source_data.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input_source_data.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output_source_data.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output_remapped.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output_source_data_remapped.yaml)
CreateUnitTestFromExec (shoc_p3_source shoc_p3_nudging
      EXE_ARGS "--args -ifile=input_source_data.yaml"
      FIXTURES_SETUP shoc_p3_source_data
      FIXTURES_REQUIRED shoc_p3_create_vertical_remap_and_weights_file)

# Run a test with nudging turned on using raw source data for nudging:
set (NUM_STEPS  5)
set (ATM_TIME_STEP 300)
set (POSTFIX nudged)
set (VERT_TYPE TIME_DEPENDENT_3D_PROFILE)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input_nudging.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input_nudging.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output_nudged.yaml)
CreateUnitTestFromExec (shoc_p3_nudged shoc_p3_nudging
      EXE_ARGS "--args -ifile=input_nudging.yaml"
      FIXTURES_REQUIRED shoc_p3_source_data)

# Run a test with nudging turned on using remapped source data for nudging:
set (NUM_STEPS  5)
set (ATM_TIME_STEP 300)
set (POSTFIX nudged_remapped)
set (VERT_TYPE STATIC_1D_VERTICAL_PROFILE)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input_nudging.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input_nudging_remapped.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output_nudged_remapped.yaml)
CreateUnitTestFromExec (shoc_p3_nudged_remapped shoc_p3_nudging
      EXE_ARGS "--args -ifile=input_nudging_remapped.yaml"
      FIXTURES_REQUIRED shoc_p3_source_data)

# Run a test with nudging using data read in glob pattern and skip vertical interpolation:
set (NUM_STEPS  5)
set (ATM_TIME_STEP 300)
set (POSTFIX nudged_glob_novert)
set (VERT_TYPE TIME_DEPENDENT_3D_PROFILE)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input_nudging_glob_novert.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input_nudging_glob_novert.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output_nudged_glob_novert.yaml)
CreateUnitTestFromExec (shoc_p3_nudging_glob_novert shoc_p3_nudging
      EXE_ARGS "--args -ifile=input_nudging_glob_novert.yaml"
      FIXTURES_REQUIRED shoc_p3_source_data)

