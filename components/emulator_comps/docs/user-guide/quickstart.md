# Quick Start Guide

This guide covers setting up and running E3SM with the Atmosphere Emulator
(EATM) component.

## Prerequisites

- Working E3SM build environment
- For LibTorch backend: LibTorch C++ library (from PyTorch installation)

## Creating a Case with EATM

For initial testing, use the all-stub compset where only EATM is active:

```console
cd ${E3SM_ROOT}/cime/scripts
./create_newcase --case my_eatm_test \
    --compset 2000_EATM_SLND_SICE_SOCN_SROF_SGLC_SWAV \
    --res gauss180x360_IcoswISC30E3r5 \
    --mach pm-cpu
```

!!! tip
    The all-stub compset (`SLND`, `SICE`, etc.) isolates EATM behavior without
    coupling interactions, making it easier to debug emulator issues.

## Build Configuration

EATM build options are controlled via the `EATM_CMAKE_OPTIONS` XML variable.

### Selecting Inference Backend

By default, EATM uses the **STUB** backend (no-op, for testing). To use LibTorch:

```console
cd my_eatm_test

# Set inference backend to LibTorch
./xmlchange EATM_CMAKE_OPTIONS="EATM_INFERENCE_BACKEND libtorch"

# Add LibTorch path (append to existing options)
./xmlchange --append EATM_CMAKE_OPTIONS=" Torch_DIR /path/to/libtorch/share/cmake/Torch"
```

### Available Build Options

| Option | Values | Default | Description |
| ------ | ------ | ------- | ----------- |
| `EATM_INFERENCE_BACKEND` | `stub`, `libtorch` | `stub` | ML backend |
| `Torch_DIR` | path | — | Path to LibTorch CMake config |

!!! note
    Build options are set in `env_build.xml` and take effect at `./case.build` time.
    After changing `EATM_CMAKE_OPTIONS`, you must rebuild:
    `./case.build --clean-all && ./case.build`

## Runtime Configuration

EATM uses YAML configuration files for portability—parsable by both Python and C++.

### Configuration Files

| File | Location | Purpose |
| ---- | -------- | ------- |
| `defaults_yaml_eatm` | `eatm/cime_config/` | Shipped defaults |
| `user_yaml_eatm` | Case directory | User overrides |
| `atm_in` | Run directory | Merged runtime config (YAML format) |

### Configuring EATM

Create or modify `user_yaml_eatm` in your case directory:

```yaml
# User overrides for EATM configuration
eatm:
  # Grid file (SCRIP format)
  grid_file: "/path/to/grid_file.nc"

  # Build-time settings
  build:
    inference_backend: libtorch

  # Runtime settings
  runtime:
    model_path: "/path/to/ace_model.pt"
    ic_file: "/path/to/initial_conditions.nc"
    enabled: true

  # Model I/O configuration
  model_io:
    # Set to true for CNN models (ACE2), false for pointwise MLPs
    spatial_mode: true
    input_variables:
      - PRESsfc
      - TMP2m
      - SPFH2m
      - DSWRFtoa
      - global_mean_co2
    output_variables:
      - PRATEsfc
      - DLWRFsfc
      - USWRFsfc

  # Coupling options
  coupling:
    debug: false
```

### Inference Backends

| Backend | Use Case | Requirements |
| ------- | -------- | ------------ |
| `stub` | Testing, development | None (no-op) |
| `libtorch` | Production, performance | LibTorch C++ library |

!!! tip
    Use the `stub` backend for initial testing to verify coupling without
    requiring a trained AI model.

## Build and Run

```console
cd my_eatm_test
./case.setup
./case.build
./case.submit
```

## Spatial Mode

The `spatial_mode` configuration controls how data is formatted for inference:

| Mode | Use Case | Data Format |
| ---- | -------- | ----------- |
| `true` | CNN models (ACE2) | `[1, C, H, W]` - spatial grid |
| `false` | Pointwise MLPs | `[batch, C]` - per-point |

For CNN models like ACE2, the emulator:

1. Packs input fields into `[C, H, W]` tensor (channel-first)
2. Runs inference with `batch_size=1`
3. Unpacks output tensor back to field storage

## Troubleshooting

| Issue | Solution |
| ----- | -------- |
| Grid file not found | Verify `grid_file` path in `user_yaml_eatm` |
| YAML parse errors | Check YAML syntax (proper indentation, no tabs) |
| LibTorch linking errors | Ensure `LD_LIBRARY_PATH` includes LibTorch |
| Torch not found | Set `Torch_DIR` to LibTorch cmake directory |
| Grid mismatch | Ensure IC file matches grid dimensions |
| Missing input variables | Verify IC file contains required fields |

## Next Steps

- [Architecture Overview](../tech-guide/architecture.md)
- [Inference Backends](../tech-guide/inference-backends.md)
- [Standalone Development](../developer-guide/standalone-build.md)
