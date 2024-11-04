# - FindMPI4PY
# Find mpi4py includes
# This module defines:
# MPI4PY_INCLUDE_DIR, where to find mpi4py.h, etc.
# MPI4PY_FOUND

function (SetMpi4pyIncludeDir)
endfunction()

if (NOT TARGET mpi4py)
  # If user provided an include dir, we will use that, otherwise we'll ask python to find it
  if (NOT MPI4PY_INCLUDE_DIR)
    execute_process(COMMAND
      "${PYTHON_EXECUTABLE}" "-c" "import mpi4py; print (mpi4py.get_include())"
      OUTPUT_VARIABLE OUTPUT
      RESULT_VARIABLE RESULT
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    if (RESULT)
      set(MPI4PY_FOUND FALSE)
    else ()
      set (MPI4PY_INCLUDE_DIR ${OUTPUT} CACHE PATH "Path to mpi4py include directory" FORCE)
    endif()
  endif()

  # If we still don't have an include dir, it means we have no mpi4py installed
  if (NOT MPI4PY_INCLUDE_DIR)
    set(MPI4PY_FOUND FALSE)
  else ()
    add_library(mpi4py INTERFACE)
    target_include_directories(mpi4py INTERFACE SYSTEM ${MPI4PY_INCLUDE_DIR})
    set(MPI4PY_FOUND TRUE)
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(mpi4py DEFAULT_MSG MPI4PY_INCLUDE_DIR)
