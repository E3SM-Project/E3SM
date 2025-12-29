#!/bin/bash
#==============================================================================
# Build Functions for Emulator Test Scripts
#==============================================================================

do_clean() {
    if [[ -d "${BUILD_DIR}" ]]; then
        log_step "CLEAN" "Removing build directory..."
        rm -rf "${BUILD_DIR}"
        log_success "Build directory removed"
    else
        log_info "Build directory doesn't exist, nothing to clean"
    fi
}

setup_build_dir() {
    mkdir -p "${BUILD_DIR}"
    mkdir -p "${LOG_DIR}"
}

run_cmake() {
    log_step "1/3" "CMake configuration..."
    
    local cmake_log="${LOG_DIR}/cmake_configure.log"
    
    cd "${BUILD_DIR}"
    
    local skip_mpi_flag="OFF"
    if [[ "${SKIP_MPI_TESTS}" == "true" ]]; then
        skip_mpi_flag="ON"
    fi
    
    local cmake_cmd="cmake ${EMULATOR_COMPS_DIR} \
        -DEMULATOR_STANDALONE_BUILD=ON \
        -DEMULATOR_BUILD_TESTS=${BUILD_TESTS} \
        -DEMULATOR_TEST_LEVEL=${TEST_LEVEL} \
        -DEMULATOR_SKIP_MPI_TESTS=${skip_mpi_flag} \
        -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
        -DCMAKE_C_COMPILER=${CC} \
        -DCMAKE_CXX_COMPILER=${CXX} \
        -DCMAKE_Fortran_COMPILER=${FC} \
        -DEATM_INFERENCE_BACKEND=${BACKEND}"
    
    # Add LibTorch if backend is libtorch
    if [[ "${BACKEND}" == "libtorch" ]]; then
        if [[ -z "${Torch_DIR}" ]]; then
            log_error "Torch_DIR not set. LibTorch setup failed."
            exit 1
        fi
        cmake_cmd="${cmake_cmd} -DTorch_DIR=${Torch_DIR}"
        log_info "LibTorch: ${Torch_DIR}"
    fi
    
    if [[ "${VERBOSE}" == "true" ]]; then
        eval "${cmake_cmd}" 2>&1 | tee "${cmake_log}"
    elif [[ "${QUIET}" == "true" ]]; then
        eval "${cmake_cmd}" > "${cmake_log}" 2>&1
    else
        eval "${cmake_cmd}" 2>&1 | tee "${cmake_log}"
    fi
    
    local result=$?
    if [[ ${result} -eq 0 ]]; then
        log_success "CMake configuration complete"
    else
        log_error "CMake configuration failed (see ${cmake_log})"
        tail -20 "${cmake_log}"
        exit 1
    fi
}

run_build() {
    log_step "2/3" "Building (${JOBS} jobs)..."
    
    local build_log="${LOG_DIR}/build.log"
    
    cd "${BUILD_DIR}"
    
    if [[ "${VERBOSE}" == "true" ]]; then
        make -j${JOBS} 2>&1 | tee "${build_log}"
    elif [[ "${QUIET}" == "true" ]]; then
        make -j${JOBS} > "${build_log}" 2>&1
    else
        make -j${JOBS} 2>&1 | tee "${build_log}"
    fi
    
    local result=$?
    if [[ ${result} -eq 0 ]]; then
        log_success "Build complete"
    else
        log_error "Build failed (see ${build_log})"
        tail -30 "${build_log}"
        exit 1
    fi
}

run_tests() {
    if [[ "${RUN_TESTS}" != "true" ]]; then
        log_info "Skipping tests"
        return 0
    fi
    
    log_step "3/3" "Running tests..."
    
    local test_log="${LOG_DIR}/test.log"
    
    cd "${BUILD_DIR}"
    
    if [[ "${VERBOSE}" == "true" ]]; then
        ctest --output-on-failure 2>&1 | tee "${test_log}"
    elif [[ "${QUIET}" == "true" ]]; then
        ctest --output-on-failure > "${test_log}" 2>&1
    else
        ctest --output-on-failure 2>&1 | tee "${test_log}"
    fi
    local test_result=$?
    
    if [[ ${test_result} -eq 0 ]]; then
        log_success "All tests passed"
    else
        log_warning "Some tests failed (see ${test_log})"
        if [[ "${QUIET}" == "true" ]]; then
            grep -E "^\s*[0-9]+.*Failed" "${test_log}" 2>/dev/null || true
        fi
    fi
    
    return ${test_result}
}

generate_summary() {
    local summary_file="${LOG_DIR}/build_summary.txt"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    
    cat > "${summary_file}" << EOF
================================================================================
Emulator Components Build Summary
================================================================================
Timestamp:      ${timestamp}
Machine:        ${MACHINE}
Backend:        ${BACKEND}
Build Type:     ${BUILD_TYPE}
Build Dir:      ${BUILD_DIR}
Test Level:     ${TEST_LEVEL}
Skip MPI Tests: ${SKIP_MPI_TESTS}
Jobs:           ${JOBS}

Environment:
  CC:           ${CC}
  CXX:          ${CXX}
  FC:           ${FC}
  Torch_DIR:    ${Torch_DIR:-not set}

Logs:
  Configure:    ${LOG_DIR}/cmake_configure.log
  Build:        ${LOG_DIR}/build.log
  Test:         ${LOG_DIR}/test.log
================================================================================
EOF
}
