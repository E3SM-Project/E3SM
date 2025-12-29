#!/usr/bin/env python3
"""
Generate dummy TorchScript model for LibTorch backend testing.
Creates a simple linear model: output = W * input + b

Usage:
    python create_dummy_model.py [output_path]
    
If output_path is not provided, saves to same directory as script.
"""
import torch
import torch.nn as nn
import sys
import os

class DummyModel(nn.Module):
    """Simple linear model for testing."""
    def __init__(self, input_size=10, output_size=5):
        super().__init__()
        self.linear = nn.Linear(input_size, output_size)
        # Initialize with known weights for reproducibility
        # output = 0.1 * sum(inputs) + 0 = 1.0 for all-ones input
        nn.init.constant_(self.linear.weight, 0.1)
        nn.init.constant_(self.linear.bias, 0.0)
    
    def forward(self, x):
        return self.linear(x)

def main():
    input_size = 10
    output_size = 5
    
    # Get output path from args or use default
    if len(sys.argv) > 1:
        output_path = sys.argv[1]
    else:
        output_path = os.path.join(os.path.dirname(__file__), "dummy_model.pt")
    
    # Create model
    model = DummyModel(input_size, output_size)
    model.eval()
    
    # Convert to float32 precision (even though E3SM uses float64 by default)
    model = model.float()
    
    # Trace the model with float32 precision input
    dummy_input = torch.randn(1, input_size, dtype=torch.float32)
    traced_model = torch.jit.trace(model, dummy_input)
    
    # Save
    traced_model.save(output_path)
    
    print(f"Created: {output_path}")
    print(f"  Input:  (batch, {input_size})")
    print(f"  Output: (batch, {output_size})")
    
    # Verify with test input
    test_input = torch.ones(1, input_size, dtype=torch.float32)
    output = traced_model(test_input)
    expected = 1.0  # 0.1 * 10 = 1.0
    actual = output[0, 0].item()
    if abs(actual - expected) < 0.01:
        print(f"  Verified: output[0] = {actual:.3f} (expected {expected})")
    else:
        print(f"  WARNING: output[0] = {actual:.3f}, expected {expected}")
        sys.exit(1)

if __name__ == "__main__":
    main()
