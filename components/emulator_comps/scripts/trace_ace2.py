#!/usr/bin/env python3
"""
Load an ACE model checkpoint and trace it for interoperability with LibTorch.

KNOWN ISSUE: The ACE2 SHT (Spherical Harmonic Transform) layer stores precomputed
weights as class attributes rather than registered buffers. These don't move when
calling model.to(device), causing device mismatch errors when tracing on GPU.

WORKAROUND: Trace on CPU (--device_type cpu), then use torch.jit.load(path, map_location)
in C++ to load the model onto CUDA. The LibTorch backend handles this automatically.

Alternatively, fix the SHT layer in fme/sht_fix.py to register its weights as buffers.
"""

import os
from pathlib import Path
from typing import Optional

import torch

from fme.ace.stepper import load_stepper


def trace_to_torchscript(
    model: torch.nn.Module,
    dummy_input: torch.Tensor,
    filename: Optional[str] = "traced_model.pt",
) -> None:
    """
    Save PyTorch model to TorchScript using tracing.

    Parameters
    ----------
    model : torch.NN.Module
        a PyTorch model
    dummy_input : torch.Tensor
        appropriate size Tensor to act as input to model
    filename : str
        name of file to save to
    """
    # FIXME: torch.jit.optimize_for_inference() when PyTorch issue #81085 is resolved
    traced_model = torch.jit.trace(model, dummy_input)
    frozen_model = torch.jit.freeze(traced_model)
    frozen_model.save(filename)


def load_torchscript(filename: Optional[str] = "saved_model.pt") -> torch.nn.Module:
    """
    Load a TorchScript from file.

    Parameters
    ----------
    filename : str
        name of file containing TorchScript model
    """
    model = torch.jit.load(filename)

    return model


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--device_type",
        help="Device type to run the inference on",
        type=str,
        choices=["cpu", "cuda", "hip", "xpu", "mps"],
        default="cpu",
    )
    parser.add_argument(
        "--checkpoint",
        help="Path to the file containing the PyTorch model",
        type=Path,
    )
    parsed_args = parser.parse_args()
    device_type = parsed_args.device_type
    checkpoint = parsed_args.checkpoint
    # add device signature to name as well
    traced_torch_script = Path(str(checkpoint).replace(".tar", f"_traced_{device_type}.pt"))

    # validate arguments
    if not checkpoint.exists():
        msg = f"Error: checkpoint file '{checkpoint}' does not exist."
        raise FileNotFoundError(msg)
    if not checkpoint.is_file():
        msg = f"Error: Path '{checkpoint}' is not a checkpoint file."
        raise IsADirectoryError(msg)

    # Load checkpoint and prepare for tracing
    stepper = load_stepper(checkpoint)

    # Prepare dummy input
    img_shape = stepper.get_state()['dataset_info']['img_shape']
    input_names = stepper.config.input_names
    output_names = stepper.config.output_names
    print(f"DEBUG: Model inputs ({len(input_names)}): {input_names}")
    print(f"DEBUG: Model outputs ({len(output_names)}): {output_names}")

    example_input = torch.rand(1, len(input_names), *img_shape)

    # Transfer the model and inputs to requested device
    if device_type == "hip":
        device = torch.device("cuda")  # NOTE: HIP is treated as CUDA in FTorch
    else:
        device = torch.device(device_type)
    
    # Warn about known ACE2 SHT layer issues on CUDA
    if device.type == "cuda":
        print("WARNING: Tracing on CUDA may fail due to ACE2 SHT layer device issues.")
        print("         If you encounter device mismatch errors, trace on CPU instead:")
        print("         python trace_ace2.py --device_type cpu --checkpoint ...")
        print("")
        # set_device requires an explicit index
        cuda_index = device.index if device.index is not None else 0
        torch.cuda.set_device(cuda_index)
    
    # Set default device so tensors created during forward pass are on the correct device
    torch.set_default_device(device)

    if len(stepper.modules) != 1:
        msg = (f"Expected a single module in the stepper, but found {len(stepper.modules)}.")
        raise RuntimeError(msg)

    trained_model = stepper.modules[0].to(device)
    trained_model.eval()
    
    # CRITICAL: Some model buffers (e.g., SHT weights) may not move with .to()
    # Recursively ensure ALL parameters and buffers are on the target device
    def move_all_tensors_to_device(module, device):
        for name, buf in module._buffers.items():
            if buf is not None and buf.device != device:
                module._buffers[name] = buf.to(device)
        for name, param in module._parameters.items():
            if param is not None and param.device != device:
                module._parameters[name] = param.to(device)
        for child in module.children():
            move_all_tensors_to_device(child, device)
    
    move_all_tensors_to_device(trained_model, device)
    
    example_input = example_input.to(device)

    # Run model for dummy inputs, should catch any high-level issues
    example_outputs = trained_model(example_input)

    # Save model
    trace_to_torchscript(
        trained_model, example_input, filename=traced_torch_script
    )

    # =====================================================
    # Check model saved OK
    # =====================================================

    # Load torchscript and run model as a test
    example_input = 2.0 * example_input
    trained_model_testing_output = trained_model(example_input)

    ts_model = load_torchscript(filename=traced_torch_script)
    ts_model_output = ts_model(example_input)

    if torch.all(ts_model_output.eq(trained_model_testing_output)):
        print("Saved TorchScript model working as expected in a basic test.")
        print("Users should perform further validation as appropriate.")
        print(f"Output saved to: {traced_torch_script}")
    else:
        model_error = (
            "Saved Torchscript model is not performing as expected.\n"
            "Consider using scripting if you used tracing, or investigate further."
        )
        raise RuntimeError(model_error)
