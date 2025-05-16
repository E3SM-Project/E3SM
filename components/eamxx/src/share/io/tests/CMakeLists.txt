#include(ScreamUtils)
include(ScreamUtils)

# The restart test requires cprnc
include (BuildCprnc)
BuildCprnc()

## Test scorpio interfaces
CreateUnitTest(scorpio_interface_tests "scorpio_interface_tests.cpp"
  LIBS eamxx_scorpio_interface LABELS "io"
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)

## Test io utils
CreateUnitTest(io_utils "io_utils.cpp"
  LIBS scream_io LABELS io
  PROPERTIES RESOURCE_LOCK rpointer_file
)

# Test creation of diagnostic from diag_field_name
CreateUnitTest(create_diag "create_diag.cpp"
  LIBS diagnostics scream_io
  LABELS io diagnostics
)

## Test basic output (no packs, no diags, all avg types, all freq units)
CreateUnitTest(io_basic "io_basic.cpp"
  LIBS scream_io LABELS io
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)

## Test output where we write one file per month
CreateUnitTest(io_monthly "io_monthly.cpp"
  LIBS scream_io LABELS io
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)

## Test basic output (no packs, no diags, all avg types, all freq units)
CreateUnitTest(io_filled "io_filled.cpp"
  LIBS scream_io LABELS io
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)

## Test packed I/O
CreateUnitTest(io_packed "io_packed.cpp"
  LIBS scream_io LABELS io
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)

## Test diagnostic output
CreateUnitTest(io_diags "io_diags.cpp"
  LIBS scream_io LABELS io
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)

# Test output on SE grid
CreateUnitTest(io_se_grid "io_se_grid.cpp"
  LIBS scream_io LABELS io
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)

## Test output restart
# NOTE: These tests cannot run in parallel due to contention of the rpointer file
CreateUnitTest(output_restart "output_restart.cpp"
  LIBS scream_io LABELS io
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
  FIXTURES_SETUP_INDIVIDUAL restart_check_setup
  PROPERTIES RESOURCE_LOCK rpointer_file
)

# For each avg_type and rank combination, compare the monolithic and restared run
include (CompareNCFiles)
foreach (AVG_TYPE IN ITEMS INSTANT AVERAGE)
  foreach (MPI_RANKS RANGE 1 ${SCREAM_TEST_MAX_RANKS})
    CompareNCFiles (
      TEST_NAME output_restart_check_${AVG_TYPE}_np${MPI_RANKS}
      SRC_FILE monolithic.${AVG_TYPE}.nsteps_x10.np${MPI_RANKS}.2000-01-01-00000.nc
      TGT_FILE restarted.${AVG_TYPE}.nsteps_x10.np${MPI_RANKS}.2000-01-01-00000.nc
      LABELS io
      FIXTURES_REQUIRED restart_check_setup_np${MPI_RANKS}_omp1
    )
  endforeach()
endforeach()

## Test remap output
CreateUnitTest(io_remap_test "io_remap_test.cpp"
  LIBS scream_io diagnostics LABELS io remap
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)

## Test remap output when map file is sub-sampling (ARM-style)
CreateUnitTest(io_horiz_sampling "io_horiz_sampling.cpp"
  LIBS scream_io LABELS io remap
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)

## Test single-column reader
CreateUnitTest(io_scm_reader "io_scm_reader.cpp"
  LIBS scream_io LABELS io
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)
