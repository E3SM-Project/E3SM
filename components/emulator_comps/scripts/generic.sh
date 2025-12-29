#!/bin/bash
#==============================================================================
# Generic Linux Environment Setup
#==============================================================================

load_modules_generic() {
    log_info "Using generic Linux configuration..."
    
    # Assume MPI compilers are in PATH
    export CC=${CC:-mpicc}
    export CXX=${CXX:-mpicxx}
    export FC=${FC:-mpif90}
    
    # Check for NetCDF
    if [[ -z "${NETCDF_DIR}" ]] && [[ -z "${NetCDF_C_PATH}" ]]; then
        log_warning "NETCDF_DIR not set. You may need to set it manually"
    fi
    
    log_success "Generic configuration applied"
}

#------------------------------------------------------------------------------
# LibTorch Detection (Generic)
#------------------------------------------------------------------------------
setup_libtorch_generic() {
    log_info "Checking LibTorch configuration..."
    
    # Check if Torch_DIR is already set
    if [[ -n "${Torch_DIR}" ]]; then
        if [[ -d "${Torch_DIR}" ]]; then
            log_success "LibTorch found: ${Torch_DIR}"
            return 0
        else
            log_error "Torch_DIR is set but directory does not exist: ${Torch_DIR}"
            exit 1
        fi
    fi
    
    # Try to detect from Python (if available)
    if command -v python &>/dev/null; then
        local torch_path
        torch_path=$(python -c "import torch; print(torch.__path__[0])" 2>/dev/null)
        
        if [[ -n "${torch_path}" ]]; then
            TORCH_DIR="${torch_path}/share/cmake/Torch"
            if [[ -d "${TORCH_DIR}" ]]; then
                export Torch_DIR="${TORCH_DIR}"
                log_success "LibTorch detected from Python: ${TORCH_DIR}"
                return 0
            fi
        fi
    fi
    
    # Not found
    log_error "LibTorch not found!"
    log_error ""
    log_error "To use the libtorch backend, either:"
    log_error "  1. Set Torch_DIR to your LibTorch cmake directory:"
    log_error "     export Torch_DIR=/path/to/libtorch/share/cmake/Torch"
    log_error ""
    log_error "  2. Install PyTorch in your Python environment:"
    log_error "     pip install torch"
    log_error ""
    log_error "  3. Use the stub backend instead:"
    log_error "     ./test rebuild --backend stub"
    exit 1
}
