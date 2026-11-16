#!/usr/bin/env python3
"""Write the TorchScript modules the LibTorch backend test loads.

Usage: make_libtorch_fixture.py OUTDIR

The arithmetic is exact in float32, so the test can compare with ==.
"""

import os
import sys
from typing import List

import torch


class Affine(torch.nn.Module):
    def forward(self, x):
        return 2.0 * x + 1.0


class TwoHeads(torch.nn.Module):
    """Returns a tuple."""

    def forward(self, x):
        return 2.0 * x + 1.0, x - 1.0


class TwoHeadsList(torch.nn.Module):
    """Returns a list."""

    def forward(self, x) -> List[torch.Tensor]:
        return [2.0 * x + 1.0, x - 1.0]


class ChannelWeights(torch.nn.Module):
    """Fixed to 3 input channels, as a real network is; 1 output channel."""

    def __init__(self):
        super().__init__()
        self.register_buffer("w", torch.tensor([1.0, 10.0, 100.0]))

    def forward(self, x):
        return (x * self.w.view(1, -1, 1, 1)).sum(dim=1, keepdim=True)


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: make_libtorch_fixture.py OUTDIR")
    outdir = sys.argv[1]
    os.makedirs(outdir, exist_ok=True)

    # [batch, channels, ny, nx]
    example = torch.zeros(1, 3, 2, 4)

    modules = {
        "affine.pt": torch.jit.trace(Affine().eval(), example),
        "tuple.pt": torch.jit.trace(TwoHeads().eval(), example),
        "list.pt": torch.jit.script(TwoHeadsList().eval()),
        "channels.pt": torch.jit.trace(ChannelWeights().eval(), example),
    }
    for name, module in modules.items():
        torch.jit.save(module, os.path.join(outdir, name))


if __name__ == "__main__":
    main()
