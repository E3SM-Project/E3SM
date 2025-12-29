"""
Export ACE2-ERA5 Checkpoint to ONNX

Converts ACE2-ERA5 model to ONNX format for use with ONNX Runtime backend.

Usage:
    python export_ace2_to_onnx.py \\
        --checkpoint /path/to/ACE2-ERA5/checkpoint \\
        --output ace2_era5.onnx

Requirements:
    - torch
    - fme (pip install fme)
    - onnx (pip install onnx)
"""
import argparse
import torch
from pathlib import Path


def export_to_onnx(checkpoint_path: str, output_path: str, opset_version: int = 14):
    """
    Export ACE2 model to ONNX format.
    
    Args:
        checkpoint_path: Path to ACE2-ERA5 checkpoint directory
        output_path: Path to save ONNX model (.onnx)
        opset_version: ONNX opset version (14 recommended)
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
    
    # Get input shape
    try:
        if hasattr(wrapper, 'config') and hasattr(wrapper.config, 'in_channels'):
            in_channels = wrapper.config.in_channels
        else:
            in_channels = 39  # Default for ACE2-ERA5
        
        print(f"Using input shape: [batch, {in_channels}]")
        dummy_input = torch.randn(1, in_channels)
        
    except Exception as e:
        print(f"ERROR creating dummy input: {e}")
        return False
    
    # Export to ONNX
    print(f"Exporting to ONNX (opset {opset_version})...")
    
    try:
        torch.onnx.export(
            model,
            dummy_input,
            output_path,
            export_params=True,
            opset_version=opset_version,
            do_constant_folding=True,
            input_names=['input'],
            output_names=['output'],
            dynamic_axes={
                'input': {0: 'batch_size'},
                'output': {0: 'batch_size'}
            }
        )
        
        print(f"ONNX model saved to {output_path}")
        
        # Validate with ONNX
        try:
            import onnx
            onnx_model = onnx.load(output_path)
            onnx.checker.check_model(onnx_model)
            print("✓ ONNX model validation passed")
            
            # Print model info
            print("\nModel Information:")
            print(f"  Inputs: {[i.name for i in onnx_model.graph.input]}")
            print(f"  Outputs: {[o.name for o in onnx_model.graph.output]}")
            
        except ImportError:
            print("Warning: onnx library not found, skipping validation")
            print("  Install with: pip install onnx")
        
        # Test with ONNX Runtime if available
        try:
            import onnxruntime as ort
            print("\nTesting with ONNX Runtime...")
            
            session = ort.InferenceSession(output_path)
            test_input = dummy_input.numpy()
            outputs = session.run(None, {'input': test_input})
            
            print(f"✓ ONNX Runtime inference successful")
            print(f"  Output shape: {outputs[0].shape}")
            
        except ImportError:
            print("Note: onnxruntime not found, skipping runtime test")
            print("  Install with: pip install onnxruntime")
        
        return True
        
    except Exception as e:
        print(f"ERROR during ONNX export: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    parser = argparse.ArgumentParser(description="Export ACE2-ERA5 to ONNX")
    parser.add_argument("--checkpoint", required=True,
                        help="Path to ACE2-ERA5 checkpoint directory or tar file")
    parser.add_argument("--output", required=True,
                        help="Output path for ONNX model (.onnx)")
    parser.add_argument("--opset", type=int, default=14,
                        help="ONNX opset version (default: 14)")
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
    success = export_to_onnx(str(checkpoint_path), args.output, args.opset)
    
    if success:
        print("\n✓ Export complete!")
        print(f"  ONNX model: {args.output}")
        print(f"\nTo use with E3SM emulator ONNX Runtime backend:")
        print(f"  inference_backend: onnx")
        print(f"  model:")
        print(f"    path: {args.output}")
    else:
        print("\n✗ Export failed")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())
