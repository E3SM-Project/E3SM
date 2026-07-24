include(ScreamUtils)

if (NOT SCREAM_ONLY_GENERATE_BASELINES)
  CreateUnitTest(mam_photo_table_test
    SOURCES mam_photo_table_test.cpp
    LIBS mam eamxx_scorpio_interface eamxx_io yaml-cpp
    LABELS "mam;physics"
    MPI_RANKS 1
    THREADS 1
  )
  target_compile_definitions(mam_photo_table_test PRIVATE
    SCREAM_DATA_DIR="${SCREAM_DATA_DIR}"
  )
  # Ensure test input files are present in the data dir
  GetInputFile(scream/mam4xx/photolysis/table_photo_input_ts_355.yaml)
  GetInputFile(scream/mam4xx/photolysis/jlong_input_ts_355.yaml)
endif()
