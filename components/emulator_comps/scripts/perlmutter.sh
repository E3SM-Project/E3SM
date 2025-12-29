#!/bin/bash
#==============================================================================
# Perlmutter (NERSC) Environment Setup
#==============================================================================

load_modules_perlmutter() {
    log_info "Loading Perlmutter modules..."
    
    # Reset modules
    module purge 2>/dev/null || true
    
    # Load required modules (ORDER MATTERS!)
    module load PrgEnv-gnu
    module load craype-x86-milan  # CPU target for Perlmutter
    module load cmake
    module load cray-hdf5
    module load cray-netcdf
    module load cray-parallel-netcdf
    
    # Set compilers
    export CC=cc
    export CXX=CC
    export FC=ftn
    
    # Verify environment
    if [[ -z "${NETCDF_DIR}" ]]; then
        log_error "NETCDF_DIR not set after loading modules"
        exit 1
    fi
    
    log_success "Modules loaded"
}

#------------------------------------------------------------------------------
# LibTorch Detection via PyTorch Module
#------------------------------------------------------------------------------
setup_libtorch_perlmutter() {
    log_info "Detecting LibTorch from pytorch module..."
    
    # Load pytorch module (includes LibTorch)
    # Note: Use pytorch/2.6.0 for better NCCL compatibility on Perlmutter
    if ! module load pytorch/2.6.0 2>/dev/null; then
        log_error "Failed to load pytorch module"
        log_error "LibTorch backend requires pytorch module on Perlmutter"
        exit 1
    fi
    
    # Get Torch cmake directory from Python
    local torch_path
    torch_path=$(python -c "import torch; print(torch.__path__[0])" 2>/dev/null)
    
    if [[ -z "${torch_path}" ]]; then
        log_error "Could not detect torch installation path"
        exit 1
    fi
    
    TORCH_DIR="${torch_path}/share/cmake/Torch"
    
    if [[ ! -d "${TORCH_DIR}" ]]; then
        log_error "LibTorch cmake directory not found: ${TORCH_DIR}"
        exit 1
    fi
    
    export Torch_DIR="${TORCH_DIR}"
    log_success "LibTorch detected: ${TORCH_DIR}"
}
