# make_coupling_fixtures.py
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
    def forward(self, temperature: Tensor, soil_moisture: Tensor) -> Tensor:
        return 2.0 * temperature + 100.0 * soil_moisture.mean(dim=1)


def main():
    destination = Path(__file__).resolve().parent

    for filename, model in (
        ("forcing_emulator.pt", ForcingModel(nlevels=3)),
        ("land_emulator.pt", LandModel()),
    ):
        scripted = torch.jit.script(model.eval())
        scripted.save(str(destination / filename))


if __name__ == "__main__":
    main()
