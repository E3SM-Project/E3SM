set (NC_DUMP @NetCDF_C_PATHS@/bin/ncdump)
foreach (MPI_RANK RANGE 1 @SCREAM_TEST_MAX_RANKS@)
  set (SRC_FILE io_output_restart_check_np${MPI_RANK}.AVERAGE.Steps_x10.0000-01-01.000020.nc)
  set (TGT_FILE io_output_restart_np${MPI_RANK}.AVERAGE.Steps_x10.0000-01-01.000020.nc)

  # Dump vars that are expected to be equal and those expected to be different
  execute_process (
    COMMAND ${NC_DUMP} -v field_1,field_2,field_4 ${SRC_FILE}
    OUTPUT_VARIABLE ncdump_out_src
    ERROR_VARIABLE  ncdump_out_src)

  execute_process (
    COMMAND ${NC_DUMP} -v field_1,field_2,field_4 ${TGT_FILE}
    OUTPUT_VARIABLE ncdump_out_tgt
    ERROR_VARIABLE  ncdump_out_tgt)

  if (NOT "${ncdump_out1}" STREQUAL "${ncdump_out2}")
    message ("Restart test failed with ${MPI_RANK} mpi ranks.\n")
    message ("ncdump of original output:\n${ncdump_out_tgt}")
    message ("ncdump of restarted output:\n${ncdump_out_src}")
    message (FATAL_ERROR "Restart test FAILED")
  endif ()
endforeach()

message ("Restart test PASSED.")
