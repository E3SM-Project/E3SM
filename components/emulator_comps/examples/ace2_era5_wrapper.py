"""
ACE2-ERA5 Inference Wrapper for E3SM Emulator

Integrates ai2cm/ace ACE2-ERA5 climate emulator with the PyTorch backend.
Handles checkpoint loading, inference, and distributed execution.

Usage with E3SM emulator PyTorch backend:
    - Set `python_module: ace2_era5_wrapper` in atm_in
    - Set `model_path` to ACE2-ERA5 checkpoint directory or tar file
    - For distributed inference, use torchrun or set distributed=True

Example atm_in configuration:
    inference_backend: "pytorch"
    model:
      path: "/path/to/ACE2-ERA5/ace2_era5_ckpt.tar"
    python_module: "ace2_era5_wrapper"
"""
import os
import sys
import torch
import numpy as np
import tarfile
from pathlib import Path

# Global model state
_model = None
_device = None
_config = None
_initialized = False


def extract_checkpoint(tar_path: str, extract_dir: str = None) -> Path:
    """
    Extract ACE2 checkpoint from tar file.
    
    Args:
        tar_path: Path to ace2_era5_ckpt.tar
        extract_dir: Directory to extract to (default: tar_path/../checkpoint)
    
    Returns:
        Path to extracted checkpoint directory
    """
    tar_path = Path(tar_path)
    if extract_dir is None:
        extract_dir = tar_path.parent / "checkpoint"
    else:
        extract_dir = Path(extract_dir)
    
    if extract_dir.exists():
        print(f"[ACE2] Checkpoint already extracted at {extract_dir}")
        return extract_dir
    
    print(f"[ACE2] Extracting checkpoint from {tar_path}")
    extract_dir.mkdir(parents=True, exist_ok=True)
    
    with tarfile.open(tar_path, "r") as tar:
        tar.extractall(extract_dir)
    
    print(f"[ACE2] Extraction complete: {extract_dir}")
    return extract_dir


def initialize(checkpoint_path: str, device_id: int = -1, distributed: bool = False):
    """
    Initialize ACE2-ERA5 model.
    
    Args:
        checkpoint_path: Path to ace2_era5_ckpt.tar or extracted directory
        device_id: GPU device ID (-1 for CPU)
        distributed: Enable distributed inference
    """
    global _model, _device, _config, _initialized
    
    if _initialized:
        print("[ACE2] Already initialized, skipping")
        return
    
    print(f"[ACE2] Initializing ACE2-ERA5...")
    print(f"[ACE2]   checkpoint={checkpoint_path}")
    print(f"[ACE2]   device_id={device_id}")
    print(f"[ACE2]   distributed={distributed}")
    
    # Set up device
    if distributed and torch.cuda.is_available():
        # Use local rank from environment (set by torchrun or MPI)
        local_rank = int(os.environ.get("LOCAL_RANK", device_id))
        _device = torch.device(f"cuda:{local_rank}")
        torch.cuda.set_device(_device)
        
        # Initialize distributed process group
        import torch.distributed as dist
        if not dist.is_initialized():
            # Try to initialize from torchrun environment
            if "RANK" in os.environ and "WORLD_SIZE" in os.environ:
                print(f"[ACE2] Initializing distributed (rank={os.environ['RANK']}, "
                      f"world_size={os.environ['WORLD_SIZE']})")
                dist.init_process_group(backend="nccl")
            else:
                print("[ACE2] Warning: distributed=True but no torchrun environment detected")
                distributed = False
    else:
        _device = torch.device(f"cuda:{device_id}" if device_id >= 0 and torch.cuda.is_available() else "cpu")
    
    print(f"[ACE2] Using device: {_device}")
    
    # Extract checkpoint if needed
    checkpoint_path = Path(checkpoint_path)
    if str(checkpoint_path).endswith(".tar"):
        checkpoint_path = extract_checkpoint(checkpoint_path)
    
    # Load ACE2 model
    try:
        # Try using fme library first
        from fme.ace.inference.inference import InferenceWrapper
        from fme.core.device import get_device
        
        print("[ACE2] Loading model using fme library...")
        wrapper = InferenceWrapper.from_checkpoint(str(checkpoint_path))
        _model = wrapper.model.to(_device)
        _model.eval()
        _config = wrapper.config
        
        print(f"[ACE2] Model loaded successfully from checkpoint")
        
    except ImportError as e:
        print(f"[ACE2] Warning: fme library not available ({e})")
        print("[ACE2] Attempting fallback: loading as TorchScript...")
        
        # Fallback: load as TorchScript if fme not available
        model_path = checkpoint_path / "model.pt"
        if not model_path.exists():
            raise FileNotFoundError(
                f"Neither fme library nor TorchScript model found. "
                f"Install fme (pip install fme) or export model to TorchScript."
            )
        
        _model = torch.jit.load(str(model_path), map_location=_device)
        _model.eval()
        print(f"[ACE2] TorchScript model loaded from {model_path}")
    
    _initialized = True
    print("[ACE2] Initialization complete")


def inference(inputs: np.ndarray) -> np.ndarray:
    """
    Run ACE2-ERA5 inference.
    
    Args:
        inputs: NumPy array [batch_size, input_channels]
                Input channels should match ACE2-ERA5 requirements
    
    Returns:
        NumPy array [batch_size, output_channels]
        Output channels match ACE2-ERA5 model outputs
    """
    global _model, _device, _initialized
    
    if not _initialized or _model is None:
        raise RuntimeError("Model not initialized. Call initialize() first.")
    
    with torch.no_grad():
        # Convert to tensor
        x = torch.from_numpy(inputs).to(_device, dtype=torch.float32)
        
        # Run inference
        y = _model(x)
        
        # Convert back to NumPy
        return y.cpu().numpy().astype(np.float64)


def finalize():
    """
    Cleanup distributed resources and model.
    """
    global _model, _device, _initialized
    
    print("[ACE2] Finalizing...")
    
    # Cleanup distributed
    try:
        import torch.distributed as dist
        if dist.is_initialized():
            print("[ACE2] Destroying process group")
            dist.destroy_process_group()
    except ImportError:
        pass
    
    _model = None
    _device = None
    _initialized = False
    
    print("[ACE2] Finalized")


if __name__ == "__main__":
    # Simple test
    import argparse
    
    parser = argparse.ArgumentParser(description="Test ACE2-ERA5 wrapper")
    parser.add_argument("--checkpoint", required=True, help="Path to ACE2-ERA5 checkpoint")
    parser.add_argument("--device-id", type=int, default=-1, help="GPU device ID")
    parser.add_argument("--batch-size", type=int, default=4, help="Test batch size")
    parser.add_argument("--distributed", action="store_true", help="Test distributed mode")
    args = parser.parse_args()
    
    # Initialize
    initialize(args.checkpoint, args.device_id, args.distributed)
    
    # Test inference (assuming 39 input channels for ACE2-ERA5)
    test_inputs = np.random.randn(args.batch_size, 39).astype(np.float64)
    print(f"[Test] Running inference with shape {test_inputs.shape}")
    
    outputs = inference(test_inputs)
    print(f"[Test] Output shape: {outputs.shape}")
    print(f"[Test] Output range: [{outputs.min():.4f}, {outputs.max():.4f}]")
    
    # Finalize
    finalize()
    print("[Test] Success!")
