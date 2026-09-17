# Check whether the MPI installation in use actually supports operating on
# GPU (device) memory directly, a.k.a. "GPU-aware" (or "CUDA-aware"/
# "ROCm-aware") MPI. This is what SCREAM_MPI_ON_DEVICE=ON requires: with that
# option on, EAMxx passes device pointers straight to MPI calls, so an MPI
# that cannot handle them will typically segfault, hang, or silently corrupt
# data at run time, often far from the actual root cause.
#
# There is no portable compile-time API to query this (and, as it turns out,
# even implementation-specific hints like the presence of MPIX_GPU_SUPPORT_*
# macros in mpi.h are NOT reliable: those macros can be present in the header
# while the MPI library itself was not actually built with device support).
# So this module's primary check is FUNCTIONAL: it builds and runs a tiny
# program that allocates a device buffer, ships it through MPI via a
# self-message (MPI_Sendrecv_replace, rank 0 to rank 0), and checks that the
# right value comes back. A non-GPU-aware MPI will typically crash (it will
# try to dereference the device pointer from the host) or, less commonly,
# return corrupted data; either way the diagnostic below detects it.
#
# If that functional test cannot be built/run in this environment (e.g. we
# cannot find the CUDA/HIP runtime to link against, or we are cross
# compiling, or the process cannot even reach MPI_Init here -- e.g. because
# launching an MPI program directly, without a batch scheduler, doesn't work
# on this machine), we fall back to a couple of well-known, implementation-
# specific signals (Open MPI's ompi_info, Cray MPICH's
# MPICH_GPU_SUPPORT_ENABLED). If nothing is conclusive, we only warn -- never
# silently claim support, but also never hard-fail on our own diagnostic's
# limitations, since that would be worse than an inconclusive check.

function (eamxx_check_mpi_gpu_aware)
  if (NOT Kokkos_ENABLE_CUDA AND NOT Kokkos_ENABLE_HIP AND NOT Kokkos_ENABLE_SYCL)
    return ()
  endif ()

  eamxx_run_mpi_gpu_aware_probe (probe_status backend_name)

  if (probe_status STREQUAL "AWARE")
    message (STATUS "SCREAM_MPI_ON_DEVICE=ON: diagnostic confirms MPI can operate on ${backend_name} device pointers.")
    return ()
  elseif (probe_status STREQUAL "NOT_AWARE")
    eamxx_report_mpi_not_gpu_aware ("a diagnostic program that ships a ${backend_name} device pointer through MPI (self MPI_Sendrecv_replace) failed")
    return ()
  endif ()

  # probe_status is "INCONCLUSIVE": the functional probe could not build or
  # run in this environment. Fall back to known implementation-specific
  # signals before giving up and just warning.
  eamxx_mpi_gpu_aware_signal_fallback (fallback_status fallback_source)
  if (fallback_status STREQUAL "AWARE")
    message (STATUS "SCREAM_MPI_ON_DEVICE=ON: ${fallback_source} indicates MPI is GPU-aware.")
  elseif (fallback_status STREQUAL "NOT_AWARE")
    eamxx_report_mpi_not_gpu_aware ("${fallback_source}")
  else ()
    message (WARNING
      "SCREAM_MPI_ON_DEVICE is ON, but EAMxx could not determine whether the underlying MPI "
      "installation supports GPU-aware communication: the functional diagnostic could not be "
      "built/run in this environment, and no known fallback signal (ompi_info's CUDA support "
      "flag, MPICH_GPU_SUPPORT_ENABLED, MV2_USE_CUDA) was found either.\n"
      "  If this MPI does NOT support device pointers, expect crashes or silently wrong "
      "results at run time. If unsure, re-configure with -DSCREAM_MPI_ON_DEVICE=OFF."
    )
  endif ()
endfunction ()

function (eamxx_report_mpi_not_gpu_aware reason)
  message (FATAL_ERROR
    "SCREAM_MPI_ON_DEVICE is ON, but the MPI installation being used does not appear to "
    "support GPU-aware communication (i.e., operating directly on device pointers): "
    "${reason}.\n"
    "  Fix options:\n"
    "   - Load/use an MPI (and, on Cray systems, the matching craype-accel-* module) "
    "that provides GPU-aware MPI, or\n"
    "   - Re-configure with -DSCREAM_MPI_ON_DEVICE=OFF, so EAMxx stages MPI data "
    "through host buffers instead."
  )
endfunction ()

# Builds and runs the functional GPU-aware-MPI probe. Sets <status_var> in
# the parent scope to one of AWARE / NOT_AWARE / INCONCLUSIVE, and
# <backend_var> to the GPU backend that was probed.
function (eamxx_run_mpi_gpu_aware_probe status_var backend_var)
  set (${status_var} "INCONCLUSIVE" PARENT_SCOPE)
  set (${backend_var} "" PARENT_SCOPE)

  if (CMAKE_CROSSCOMPILING)
    message (STATUS "EAMxx: cross-compiling, so the MPI GPU-awareness diagnostic cannot be run at configure time.")
    return ()
  endif ()

  set (probe_dir ${CMAKE_BINARY_DIR}/CMakeFiles/CheckMpiGpuAware)
  file (MAKE_DIRECTORY ${probe_dir})
  set (src_file ${probe_dir}/check_mpi_gpu_aware.cpp)
  set (exe_file ${probe_dir}/check_mpi_gpu_aware${CMAKE_EXECUTABLE_SUFFIX})

  set (link_libs "")
  set (include_dirs "")

  if (Kokkos_ENABLE_CUDA)
    set (backend "CUDA")
    set (runtime_header "cuda_runtime.h")
    set (malloc_fn "cudaMalloc")
    set (memcpy_fn "cudaMemcpy")
    set (free_fn "cudaFree")
    set (h2d "cudaMemcpyHostToDevice")
    set (d2h "cudaMemcpyDeviceToHost")

    if (CMAKE_CUDA_COMPILER)
      get_filename_component (cuda_bin_dir ${CMAKE_CUDA_COMPILER} DIRECTORY)
      get_filename_component (cuda_root_hint ${cuda_bin_dir} DIRECTORY)
    endif ()
    find_library (EAMXX_CUDART_LIBRARY NAMES cudart
      HINTS ${cuda_root_hint} ENV CUDA_HOME ENV CUDA_PATH
      PATH_SUFFIXES lib64 lib lib/x64)
    find_path (EAMXX_CUDA_INCLUDE_DIR NAMES cuda_runtime.h
      HINTS ${cuda_root_hint} ENV CUDA_HOME ENV CUDA_PATH
      PATH_SUFFIXES include)

    if (NOT EAMXX_CUDART_LIBRARY OR NOT EAMXX_CUDA_INCLUDE_DIR)
      message (STATUS "EAMxx: could not locate the CUDA runtime library/headers to build the MPI GPU-awareness diagnostic; skipping the functional check.")
      return ()
    endif ()
    set (link_libs ${EAMXX_CUDART_LIBRARY})
    set (include_dirs ${EAMXX_CUDA_INCLUDE_DIR})

  elseif (Kokkos_ENABLE_HIP)
    set (backend "HIP")
    set (runtime_header "hip/hip_runtime.h")
    set (malloc_fn "hipMalloc")
    set (memcpy_fn "hipMemcpy")
    set (free_fn "hipFree")
    set (h2d "hipMemcpyHostToDevice")
    set (d2h "hipMemcpyDeviceToHost")

    if (CMAKE_HIP_COMPILER)
      get_filename_component (hip_bin_dir ${CMAKE_HIP_COMPILER} DIRECTORY)
      get_filename_component (hip_root_hint ${hip_bin_dir} DIRECTORY)
    endif ()
    find_library (EAMXX_HIP_RUNTIME_LIBRARY NAMES amdhip64
      HINTS ${hip_root_hint} ENV ROCM_PATH
      PATH_SUFFIXES lib lib64)
    find_path (EAMXX_HIP_INCLUDE_DIR NAMES hip/hip_runtime.h
      HINTS ${hip_root_hint} ENV ROCM_PATH
      PATH_SUFFIXES include)

    if (NOT EAMXX_HIP_RUNTIME_LIBRARY OR NOT EAMXX_HIP_INCLUDE_DIR)
      message (STATUS "EAMxx: could not locate the HIP runtime library/headers to build the MPI GPU-awareness diagnostic; skipping the functional check.")
      return ()
    endif ()
    set (link_libs ${EAMXX_HIP_RUNTIME_LIBRARY})
    set (include_dirs ${EAMXX_HIP_INCLUDE_DIR})

  elseif (Kokkos_ENABLE_SYCL)
    set (backend "SYCL")
    # SYCL device allocation/queues don't map onto the malloc/memcpy/free
    # pattern used below for CUDA/HIP, so it gets its own source template.
  endif ()

  set (${backend_var} ${backend} PARENT_SCOPE)

  if (backend STREQUAL "SYCL")
    file (WRITE ${src_file} "
#include <sycl/sycl.hpp>
#include <mpi.h>
#include <cstdio>

int main (int argc, char** argv) {
  MPI_Init(&argc, &argv);
  std::fprintf(stdout, \"MPI_INIT_OK\\n\"); std::fflush(stdout);

  sycl::queue q;
  const int value = 42;
  int* dev_buf = sycl::malloc_device<int>(1, q);
  q.memcpy(dev_buf, &value, sizeof(int)).wait();

  MPI_Status status;
  MPI_Sendrecv_replace(dev_buf, 1, MPI_INT, 0, 0, 0, 0, MPI_COMM_WORLD, &status);

  int result = 0;
  q.memcpy(&result, dev_buf, sizeof(int)).wait();
  sycl::free(dev_buf, q);

  MPI_Finalize();
  if (result == value) {
    std::fprintf(stdout, \"GPU_AWARE_OK\\n\");
    return 0;
  }
  return 1;
}
")
    # Best-effort: forward whatever SYCL flags EAMxx's machine file uses, if any.
    string (APPEND CMAKE_CXX_FLAGS " ${SYCL_COMPILE_FLAGS} ${SYCL_FLAGS}")
  else ()
    file (WRITE ${src_file} "
#include <${runtime_header}>
#include <mpi.h>
#include <cstdio>

int main (int argc, char** argv) {
  MPI_Init(&argc, &argv);
  std::fprintf(stdout, \"MPI_INIT_OK\\n\"); std::fflush(stdout);

  const int value = 42;
  int host_buf = value;
  int* dev_buf = nullptr;
  ${malloc_fn}((void**)&dev_buf, sizeof(int));
  ${memcpy_fn}(dev_buf, &host_buf, sizeof(int), ${h2d});

  MPI_Status status;
  MPI_Sendrecv_replace(dev_buf, 1, MPI_INT, 0, 0, 0, 0, MPI_COMM_WORLD, &status);

  int result = 0;
  ${memcpy_fn}(&result, dev_buf, sizeof(int), ${d2h});
  ${free_fn}(dev_buf);

  MPI_Finalize();
  if (result == value) {
    std::fprintf(stdout, \"GPU_AWARE_OK\\n\");
    return 0;
  }
  return 1;
}
")
  endif ()

  set (try_compile_extra_args "")
  if (link_libs)
    list (APPEND try_compile_extra_args LINK_LIBRARIES ${link_libs})
  endif ()

  try_compile (probe_compiled
    ${probe_dir}/build
    ${src_file}
    ${try_compile_extra_args}
    CMAKE_FLAGS
      "-DINCLUDE_DIRECTORIES=${include_dirs}"
      "-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
    COPY_FILE ${exe_file}
    OUTPUT_VARIABLE probe_compile_output
  )

  if (NOT probe_compiled)
    message (STATUS
      "EAMxx: could not compile the MPI GPU-awareness diagnostic for ${backend} (this says "
      "nothing about whether MPI itself is GPU-aware); skipping the functional check.\n"
      "${probe_compile_output}")
    return ()
  endif ()

  execute_process (
    COMMAND ${exe_file}
    RESULT_VARIABLE probe_run_result
    OUTPUT_VARIABLE probe_run_output
    ERROR_VARIABLE probe_run_error
    TIMEOUT 60
  )

  if (NOT probe_run_output MATCHES "MPI_INIT_OK")
    message (STATUS
      "EAMxx: the MPI GPU-awareness diagnostic did not reach MPI_Init in this environment "
      "(e.g. running an MPI program directly may not work on a batch-scheduled login node); "
      "skipping the functional check.")
    return ()
  endif ()

  if (probe_run_result EQUAL 0 AND probe_run_output MATCHES "GPU_AWARE_OK")
    set (${status_var} "AWARE" PARENT_SCOPE)
  else ()
    set (${status_var} "NOT_AWARE" PARENT_SCOPE)
    message (STATUS "MPI GPU-awareness diagnostic (${backend}) failed. exit code: ${probe_run_result}")
    if (probe_run_output)
      message (STATUS "  stdout: ${probe_run_output}")
    endif ()
    if (probe_run_error)
      message (STATUS "  stderr: ${probe_run_error}")
    endif ()
  endif ()
endfunction ()

# Implementation-specific fallback signals, used only when the functional
# probe above was inconclusive. Sets <status_var> to AWARE / NOT_AWARE /
# INCONCLUSIVE, and <source_var> to a human-readable description of the
# signal that was used.
function (eamxx_mpi_gpu_aware_signal_fallback status_var source_var)
  set (${status_var} "INCONCLUSIVE" PARENT_SCOPE)
  set (${source_var} "" PARENT_SCOPE)

  # --- Open MPI: ask ompi_info whether CUDA support was built in --- #
  if (Kokkos_ENABLE_CUDA)
    find_program (OMPI_INFO_EXE ompi_info)
    if (OMPI_INFO_EXE)
      execute_process (
        COMMAND ${OMPI_INFO_EXE} --parsable --all
        OUTPUT_VARIABLE ompi_info_output
        ERROR_QUIET
        RESULT_VARIABLE ompi_info_result
      )
      if (ompi_info_result EQUAL 0)
        set (${status_var} "NOT_AWARE" PARENT_SCOPE)
        set (${source_var} "'${OMPI_INFO_EXE} --parsable --all' reports mpi_built_with_cuda_support:value:false" PARENT_SCOPE)
        if (ompi_info_output MATCHES "mpi_built_with_cuda_support:value:true")
          set (${status_var} "AWARE" PARENT_SCOPE)
          set (${source_var} "'${OMPI_INFO_EXE} --parsable --all' reports mpi_built_with_cuda_support:value:true" PARENT_SCOPE)
        endif ()
        return ()
      endif ()
    endif ()
  endif ()

  # --- (Cray) MPICH: MPICH_GPU_SUPPORT_ENABLED env var --- #
  if (DEFINED ENV{MPICH_GPU_SUPPORT_ENABLED})
    set (${source_var} "the MPICH_GPU_SUPPORT_ENABLED=$ENV{MPICH_GPU_SUPPORT_ENABLED} environment variable" PARENT_SCOPE)
    if ("$ENV{MPICH_GPU_SUPPORT_ENABLED}" STREQUAL "1")
      set (${status_var} "AWARE" PARENT_SCOPE)
    else ()
      set (${status_var} "NOT_AWARE" PARENT_SCOPE)
    endif ()
    return ()
  endif ()

  # --- MVAPICH2-GDR: MV2_USE_CUDA env var --- #
  if (DEFINED ENV{MV2_USE_CUDA})
    set (${source_var} "the MV2_USE_CUDA=$ENV{MV2_USE_CUDA} environment variable" PARENT_SCOPE)
    if ("$ENV{MV2_USE_CUDA}" STREQUAL "1")
      set (${status_var} "AWARE" PARENT_SCOPE)
    else ()
      set (${status_var} "NOT_AWARE" PARENT_SCOPE)
    endif ()
    return ()
  endif ()
endfunction ()
