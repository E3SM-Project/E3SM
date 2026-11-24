# LibTorch discovery for the emulator package.
#
# emulator_find_libtorch() defines the interface target emulator::torch. 
#  The search proceeds as follows:
#   1. Torch_DIR / CMAKE_PREFIX_PATH
#   2. the torch package of EMULATOR_TORCH_PYTHON, so the C++ and embedded
#      Python backends share one LibTorch (two in one process fail)
#   3. a prebuilt archive from download.pytorch.org, cached in FETCHCONTENT_BASE_DIR.


include(FetchContent)

# Sets <out_var> to the CMake prefix of the torch that EMULATOR_TORCH_PYTHON imports, or to "" if it can't import torch.
function(emulator_torch_python_prefix out_var)
  execute_process(COMMAND ${EMULATOR_TORCH_PYTHON} -c
      "import torch.utils; print(torch.utils.cmake_prefix_path)"
    RESULT_VARIABLE failed
    OUTPUT_VARIABLE prefix OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  if(failed)
    set(prefix "")
  endif()
  set(${out_var} "${prefix}" PARENT_SCOPE)
endfunction()

# Sets <out_var> to TRUE if an NVIDIA driver (i.e. a usable GPU) is present.
# A CUDA compiler alone doesn't count: nvcc from NVHPC on a CPU-only machine
# would select a CUDA LibTorch that can't configure or run!
function(emulator_torch_have_nvidia_gpu out_var)
  set(found FALSE)
  if(EXISTS "/proc/driver/nvidia/version")
    set(found TRUE)
  else()
    find_program(EMULATOR_NVIDIA_SMI nvidia-smi)
    # mark as advanced to not pollute cmake options
    mark_as_advanced(EMULATOR_NVIDIA_SMI)
    if(EMULATOR_NVIDIA_SMI)
      execute_process(COMMAND ${EMULATOR_NVIDIA_SMI} -L
        RESULT_VARIABLE failed OUTPUT_QUIET ERROR_QUIET)
      if(NOT failed)
        set(found TRUE)
      endif()
    endif()
  endif()
  set(${out_var} ${found} PARENT_SCOPE)
endfunction()

# Resolves EMULATOR_LIBTORCH_VARIANT ("auto" -> cpu or cu124) into <out_var>.
function(emulator_torch_resolve_variant out_var)
  set(variant "${EMULATOR_LIBTORCH_VARIANT}")
  if(variant STREQUAL "auto")
    emulator_torch_have_nvidia_gpu(have_gpu)
    if(have_gpu)
      set(variant cu124)
    else()
      set(variant cpu)
    endif()
    message(STATUS "LibTorch variant: ${variant} (auto-detected)")
  endif()
  set(${out_var} "${variant}" PARENT_SCOPE)
endfunction()

# A CUDA LibTorch runs find_package(CUDAToolkit REQUIRED) from its config, so
# it needs a toolkit with the same major version at configure time.
# this is the cause of the Caffe2 errors -- fail during configure instead!
function(emulator_torch_check_cuda variant)
  if(NOT variant MATCHES "^cu([0-9]+)[0-9]$")
    return()
  endif()
  set(major "${CMAKE_MATCH_1}")
  find_package(CUDAToolkit QUIET)
  if(CUDAToolkit_FOUND AND CUDAToolkit_VERSION_MAJOR VERSION_EQUAL major)
    return()
  endif()
  if(CUDAToolkit_FOUND)
    set(found "found CUDA ${CUDAToolkit_VERSION} in ${CUDAToolkit_BIN_DIR}")
  else()
    set(found "no CUDA toolkit found")
  endif()
  message(FATAL_ERROR
    "LibTorch variant '${variant}' needs a CUDA ${major}.x toolkit at "
    "configure time, but ${found}.\n"
    "Either set -DCUDAToolkit_ROOT=<matching toolkit> (on NVHPC also add "
    "its math_libs/<ver>/lib64 to CMAKE_LIBRARY_PATH), or use "
    "-DEMULATOR_LIBTORCH_VARIANT=cpu.")
endfunction()

# Sets <out_var> to the download URL: EMULATOR_LIBTORCH_URL if given, 
# else built from the version and resolved variant.
function(emulator_torch_url out_var)
  if(EMULATOR_LIBTORCH_URL)
    set(${out_var} "${EMULATOR_LIBTORCH_URL}" PARENT_SCOPE)
    return()
  endif()
  emulator_torch_resolve_variant(variant)
  emulator_torch_check_cuda(${variant})
  set(${out_var}
    "https://download.pytorch.org/libtorch/${variant}/libtorch-cxx11-abi-shared-with-deps-${EMULATOR_LIBTORCH_VERSION}%2B${variant}.zip"
    PARENT_SCOPE)
endfunction()

# Downloads <url> into the FetchContent cache and sets <out_var> to the
# extracted LibTorch root. FetchContent reuses an existing download even if
# the URL changed (e.g. cu124 -> cpu), so the URL is stamped next to it and a
# mismatch wipes the old copy.
function(emulator_torch_fetch url out_var)
  if(FETCHCONTENT_BASE_DIR)
    set(base "${FETCHCONTENT_BASE_DIR}")
  else()
    set(base "${CMAKE_BINARY_DIR}/_deps")
  endif()
  set(stamp "${base}/libtorch-url.txt")
  set(old_url "")
  if(EXISTS "${stamp}")
    file(READ "${stamp}" old_url)
  endif()
  if(NOT old_url STREQUAL url)
    if(EXISTS "${base}/libtorch-src")
      message(STATUS "LibTorch URL changed; removing previous download")
    endif()
    file(REMOVE_RECURSE "${base}/libtorch-src" "${base}/libtorch-build"
                        "${base}/libtorch-subbuild")
    file(WRITE "${stamp}" "${url}")
  endif()

  message(STATUS "LibTorch not found; fetching ${url}")
  FetchContent_Declare(libtorch URL ${url} DOWNLOAD_EXTRACT_TIMESTAMP OFF)
  FetchContent_GetProperties(libtorch)
  if(NOT libtorch_POPULATED)
    FetchContent_Populate(libtorch)
  endif()
  set(${out_var} "${libtorch_SOURCE_DIR}" PARENT_SCOPE)
endfunction()

# Ensures compatible version of LibTorch is installed.
# Defines an emulator::torch target since that persists across the whole build
# and we don't have to worry about cmake's function scope and Torch_LIBRARIES will
# be visible to the rest of the project
function(emulator_find_libtorch)
  if(TARGET emulator::torch)
    return()
  endif()

  if(NOT Torch_DIR AND NOT CMAKE_PREFIX_PATH)
    emulator_torch_python_prefix(py_prefix)
    if(py_prefix)
      message(STATUS "Using LibTorch from ${EMULATOR_TORCH_PYTHON}: ${py_prefix}")
      list(APPEND CMAKE_PREFIX_PATH "${py_prefix}")
    endif()
  endif()

  find_package(Torch QUIET)
  if(NOT Torch_FOUND)
    emulator_torch_url(url)
    emulator_torch_fetch("${url}" torch_root)
    list(APPEND CMAKE_PREFIX_PATH "${torch_root}")
    find_package(Torch REQUIRED)
  endif()

  # The torch target carries the C++ ABI flag LibTorch was built with,
  # so linking PUBLIC propagates it.
  add_library(emulator_torch INTERFACE)
  target_link_libraries(emulator_torch INTERFACE ${TORCH_LIBRARIES})
  add_library(emulator::torch ALIAS emulator_torch)
endfunction()
