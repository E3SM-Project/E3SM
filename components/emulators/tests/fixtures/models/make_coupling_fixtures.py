import os
import sys
from pathlib import Path

import torch
from torch import Tensor, nn


class ForcingModel(nn.Module):
    def __init__(self, nlevels: int):
        super().__init__()
        self.nlevels = nlevels

    def forward(self, surface_flux: Tensor) -> tuple[Tensor, Tensor]:
        # Expected shapes: [ncols] and [ncols, nlevels].
        ncols = surface_flux.size(0)
        dtype = surface_flux.dtype
        device = surface_flux.device

        columns = torch.arange(ncols, dtype=dtype, device=device)
        levels = torch.arange(self.nlevels, dtype=dtype, device=device)

        temperature_out = 280.0 + columns
        moisture_out = 0.10 + 0.01 * columns.unsqueeze(1) + 0.001 * levels.unsqueeze(0)
        return temperature_out, moisture_out


class LandModel(nn.Module):
    def __init__(self, mask_val: float):
        super().__init__()
        self.mask_val = mask_val

    def forward(
        self,
        temperature: Tensor,
        soil_moisture: Tensor,
        filter_column: Tensor,
    ) -> Tensor:

        flux: Tensor = torch.full_like(temperature, self.mask_val)

        columns = filter_column[filter_column >= 0]

        flux[columns] = 2.0 * temperature[columns] + 100.0 * soil_moisture[
            columns
        ].mean(dim=1)

        return flux


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: make_libtorch_fixture.py OUTDIR")
    destination = Path(sys.argv[1])
    os.makedirs(destination, exist_ok=True)

    for filename, model in (
        ("forcing_emulator.pt", ForcingModel(nlevels=3)),
        ("land_emulator.pt", LandModel(mask_val=-9.99)),
    ):
        scripted = torch.jit.script(model.eval())
        scripted.save(str(destination / filename))


if __name__ == "__main__":
    main()
