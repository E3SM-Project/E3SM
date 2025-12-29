# Inference Backend Examples

This directory contains examples and utilities for integrating AI models with the E3SM emulator inference backends.

## ACE2-ERA5 Integration

These scripts demonstrate integration with the [ai2cm/ace](https://github.com/ai2cm/ace) ACE2-ERA5 climate emulator.

### Quick Start

1. **Install dependencies**:

   ```bash
   pip install fme torch onnx onnxruntime
   ```

2. **Export ACE2-ERA5 model** (choose one):

   **For LibTorch backend**:

   ```bash
   python export_ace2_to_torchscript.py \
     --checkpoint /pscratch/sd/m/mahf708/ACE2-ERA5/checkpoint \
     --output ace2_era5_scripted.pt
   ```

   **For ONNX Runtime backend**:

   ```bash
   python export_ace2_to_onnx.py \
     --checkpoint /pscratch/sd/m/mahf708/ACE2-ERA5/checkpoint \
     --output ace2_era5.onnx
   ```

3. **Test the wrapper** (for PyTorch backend):

   ```bash
   python ace2_era5_wrapper.py \
     --checkpoint /pscratch/sd/m/mahf708/ACE2-ERA5/checkpoint \
     --device-id 0
   ```

## Files

- **`ace2_era5_wrapper.py`**: Python module for PyTorch backend
  - Loads ACE2-ERA5 checkpoint using fme library
  - Supports distributed inference with torchrun
  - Compatible with E3SM emulator PyTorch backend

- **`export_ace2_to_torchscript.py`**: Export script for LibTorch backend
  - Converts ACE2-ERA5 to TorchScript format
  - Validates exported model
  - No Python runtime needed for inference

- **`export_ace2_to_onnx.py`**: Export script for ONNX Runtime backend
  - Converts ACE2-ERA5 to ONNX format
  - Cross-framework compatibility
  - Optimized CPU inference

## Usage with E3SM Emulator

### PyTorch Backend

Configure in `atm_in`:

```yaml
eatm:
  inference_backend: "pytorch"
  model:
    path: "/pscratch/sd/m/mahf708/ACE2-ERA5/checkpoint"
  python_module: "ace2_era5_wrapper"
  python_executable: "/path/to/conda/envs/ace/bin/python"
```

### LibTorch Backend

Configure in `atm_in`:

```yaml
eatm:
  inference_backend: "libtorch"
  model:
    path: "/path/to/ace2_era5_scripted.pt"
```

### ONNX Runtime Backend

Configure in `atm_in`:

```yaml
eatm:
  inference_backend: "onnx"
  model:
    path: "/path/to/ace2_era5.onnx"
```

## Distributed Inference

For multi-GPU inference with PyTorch backend:

```bash
torchrun --nproc_per_node=4 ace2_era5_wrapper.py \
  --checkpoint /pscratch/sd/m/mahf708/ACE2-ERA5/checkpoint \
  --distributed
```

## NERSC Perlmutter

On NERSC Perlmutter:

```bash
# Load PyTorch module
module load pytorch/2.0.1

# For LibTorch: use system LibTorch
cd components/emulator_comps
./standalone_build.sh --configure \
  -DEMULATOR_HAS_LIBTORCH=ON \
  -DNERSC=ON

# For ONNX Runtime: install in conda environment
conda install -c conda-forge onnxruntime
```

## See Also

- [Inference Backends Documentation](../docs/tech-guide/inference-backends.md)
- [ACE2-ERA5 Hugging Face](https://huggingface.co/allenai/ACE2-ERA5)
- [ai2cm/ace GitHub](https://github.com/ai2cm/ace)
