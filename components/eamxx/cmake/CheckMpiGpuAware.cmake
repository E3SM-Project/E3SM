# Best-effort check that the MPI installation in use actually supports
# operating on GPU (device) memory directly, a.k.a. "GPU-aware" (or
# "CUDA-aware"/"ROCm-aware") MPI. This is what SCREAM_MPI_ON_DEVICE=ON
# requires: with that option on, EAMxx passes device pointers straight to
# MPI calls, so an MPI that cannot handle them will typically segfault, hang,
# or silently corrupt data at run time, often far from the actual root cause.
#
# There is no single portable API to query this at CMake-configure time, so
# this function relies on implementation/environment-specific signals:
#  - Open MPI: 'ompi_info' reports whether it was built with CUDA support.
#  - (Cray) MPICH: the MPICH_GPU_SUPPORT_ENABLED env var, set by Cray's
#    module system when a GPU-transport module is loaded (e.g.
#    craype-accel-nvidia80 or craype-accel-amd-gfx90a).
#  - MVAPICH2-GDR: the MV2_USE_CUDA env var.
#
# If none of these signals can be evaluated, GPU-awareness cannot be
# determined, and we only warn (rather than error), since erroring out on an
# unrecognized-but-possibly-fine MPI would be worse than a false negative.

function (eamxx_check_mpi_gpu_aware)
  set (found_signal FALSE)
  set (is_gpu_aware FALSE)
  set (signal_source "")

  # --- Open MPI: ask ompi_info whether CUDA support was built in --- #
  find_program (OMPI_INFO_EXE ompi_info)
  if (OMPI_INFO_EXE)
    execute_process (
      COMMAND ${OMPI_INFO_EXE} --parsable --all
      OUTPUT_VARIABLE OMPI_INFO_OUTPUT
      ERROR_QUIET
      RESULT_VARIABLE OMPI_INFO_RESULT
    )
    if (OMPI_INFO_RESULT EQUAL 0)
      if (Kokkos_ENABLE_CUDA)
        set (found_signal TRUE)
        set (signal_source "'${OMPI_INFO_EXE} --parsable --all' (mpi_built_with_cuda_support)")
        if (OMPI_INFO_OUTPUT MATCHES "mpi_built_with_cuda_support:value:true")
          set (is_gpu_aware TRUE)
        endif ()
      endif ()
      # Open MPI does not expose an analogous, well-documented query for
      # ROCm/HIP support, so we cannot confirm/deny HIP-awareness this way.
    endif ()
  endif ()

  # --- (Cray) MPICH: MPICH_GPU_SUPPORT_ENABLED env var --- #
  if (NOT found_signal AND DEFINED ENV{MPICH_GPU_SUPPORT_ENABLED})
    set (found_signal TRUE)
    set (signal_source "MPICH_GPU_SUPPORT_ENABLED=$ENV{MPICH_GPU_SUPPORT_ENABLED} env var")
    if ("$ENV{MPICH_GPU_SUPPORT_ENABLED}" STREQUAL "1")
      set (is_gpu_aware TRUE)
    endif ()
  endif ()

  # --- MVAPICH2-GDR: MV2_USE_CUDA env var --- #
  if (NOT found_signal AND DEFINED ENV{MV2_USE_CUDA})
    set (found_signal TRUE)
    set (signal_source "MV2_USE_CUDA=$ENV{MV2_USE_CUDA} env var")
    if ("$ENV{MV2_USE_CUDA}" STREQUAL "1")
      set (is_gpu_aware TRUE)
    endif ()
  endif ()

  if (found_signal AND NOT is_gpu_aware)
    message (FATAL_ERROR
      "SCREAM_MPI_ON_DEVICE is ON, but the MPI installation being used does not appear to "
      "support GPU-aware communication (i.e., operating directly on device pointers).\n"
      "  Detected via: ${signal_source}.\n"
      "  Fix options:\n"
      "   - Load/use an MPI (and, on Cray systems, the matching craype-accel-* module) "
      "that provides GPU-aware MPI, or\n"
      "   - Re-configure with -DSCREAM_MPI_ON_DEVICE=OFF, so EAMxx stages MPI data "
      "through host buffers instead."
    )
  elseif (NOT found_signal)
    message (WARNING
      "SCREAM_MPI_ON_DEVICE is ON, but EAMxx could not determine whether the underlying "
      "MPI installation supports GPU-aware communication (no known signal, such as "
      "ompi_info's CUDA support flag or the MPICH_GPU_SUPPORT_ENABLED env var, was found).\n"
      "  If this MPI does NOT support device pointers, expect crashes or silently wrong "
      "results at run time. If unsure, re-configure with -DSCREAM_MPI_ON_DEVICE=OFF."
    )
  else ()
    message (STATUS "SCREAM_MPI_ON_DEVICE=ON: MPI installation appears to be GPU-aware.")
  endif ()
endfunction ()
