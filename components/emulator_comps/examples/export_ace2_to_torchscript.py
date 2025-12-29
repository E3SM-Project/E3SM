"""
Export ACE2-ERA5 Checkpoint to TorchScript

Converts ACE2-ERA5 model to TorchScript format for use with LibTorch backend.

Usage:
    python export_ace2_to_torchscript.py \\
        --checkpoint /path/to/ACE2-ERA5/checkpoint \\
        --output ace2_era5_scripted.pt

Requirements:
    - torch
    - fme (pip install fme)
"""
import argparse
import torch
from pathlib import Path


def export_to_torchscript(checkpoint_path: str, output_path: str, method: str = "trace"):
    """
    Export ACE2 model to TorchScript format.
    
    Args:
        checkpoint_path: Path to ACE2-ERA5 checkpoint directory
        output_path: Path to save TorchScript model (.pt)
        method: 'trace' or 'script' (trace is usually more reliable)
    """
    print(f"Loading ACE2-ERA5 model from {checkpoint_path}")
    
    try:
        from fme.ace.inference.inference import InferenceWrapper
        wrapper = InferenceWrapper.from_checkpoint(checkpoint_path)
        model = wrapper.model
        model.eval()
        
        print(f"Model loaded successfully")
        print(f"  Model class: {type(model).__name__}")
        
    except ImportError:
        print("ERROR: fme library not found. Install with: pip install fme")
        return False
    except Exception as e:
        print(f"ERROR loading checkpoint: {e}")
        return False
    
    # Get input shape from model
    # ACE2-ERA5 typically expects [batch, channels] input
    # We'll create a dummy batch size of 1
    try:
        # Try to infer from model config
        if hasattr(wrapper, 'config') and hasattr(wrapper.config, 'in_channels'):
            in_channels = wrapper.config.in_channels
        else:
            # Default for ACE2-ERA5
            in_channels = 39
        
        print(f"Using input shape: [1, {in_channels}]")
        dummy_input = torch.randn(1, in_channels)
        
    except Exception as e:
        print(f"ERROR creating dummy input: {e}")
        return False
    
    # Export to TorchScript
    print(f"Exporting to TorchScript using method: {method}")
    
    try:
        if method == "trace":
            # Tracing: records operations on actual input
            scripted_model = torch.jit.trace(model, dummy_input)
        elif method == "script":
            # Scripting: analyzes code directly
            scripted_model = torch.jit.script(model)
        else:
            print(f"ERROR: Unknown method '{method}'. Use 'trace' or 'script'")
            return False
        
        # Save
        print(f"Saving TorchScript model to {output_path}")
        scripted_model.save(output_path)
        
        print("Export successful!")
        
        # Validate
        print("Validating exported model...")
        loaded_model = torch.jit.load(output_path)
        test_output = loaded_model(dummy_input)
        print(f"  Test output shape: {test_output.shape}")
        print("Validation passed!")
        
        return True
        
    except Exception as e:
        print(f"ERROR during export: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    parser = argparse.ArgumentParser(description="Export ACE2-ERA5 to TorchScript")
    parser.add_argument("--checkpoint", required=True,
                        help="Path to ACE2-ERA5 checkpoint directory or tar file")
    parser.add_argument("--output", required=True,
                        help="Output path for TorchScript model (.pt)")
    parser.add_argument("--method", choices=["trace", "script"], default="trace",
                        help="Export method: trace (default) or script")
    parser.add_argument("--extract", action="store_true",
                        help="Extract tar file if needed")
    args = parser.parse_args()
    
    checkpoint_path = Path(args.checkpoint)
    
    # Handle tar extraction if needed
    if str(checkpoint_path).endswith(".tar") and args.extract:
        print(f"Extracting checkpoint from tar file...")
        import tarfile
        extract_dir = checkpoint_path.parent / "checkpoint"
        with tarfile.open(checkpoint_path, "r") as tar:
            tar.extractall(extract_dir)
        checkpoint_path = extract_dir
        print(f"Extracted to {checkpoint_path}")
    
    # Export
    success = export_to_torchscript(str(checkpoint_path), args.output, args.method)
    
    if success:
        print("\n✓ Export complete!")
        print(f"  TorchScript model: {args.output}")
        print(f"\nTo use with E3SM emulator LibTorch backend:")
        print(f"  inference_backend: libtorch")
        print(f"  model:")
        print(f"    path: {args.output}")
    else:
        print("\n✗ Export failed")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())
