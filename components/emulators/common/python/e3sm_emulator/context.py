"""What the coupler told us: which ranks we have, and which columns we own.

The C++ backend builds this dict from ``InferenceContext``, which in turn is
built from the *component* MPI communicator that MCT handed the emulator.
Everything downstream is derived from here and from nothing else.

That is the point.  A distributed model that discovers its rank from
``SLURM_PROCID`` / ``SLURM_NTASKS`` sees the *whole job* in a coupled run, and
a process group built from those blocks forever waiting for ocean and land
ranks that will never call in.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class Context:
    """Ranks and this rank's share of the grid."""

    rank: int = 0
    world_size: int = 1
    node_name: str = ""

    nx: int = 0
    ny: int = 0
    num_global_cols: int = 0
    col_gids: np.ndarray = field(default_factory=lambda: np.empty(0, np.int64))
    lat: np.ndarray = field(default_factory=lambda: np.empty(0, np.float64))
    lon: np.ndarray = field(default_factory=lambda: np.empty(0, np.float64))

    @classmethod
    def from_dict(cls, data: dict) -> "Context":
        """Build from the dict the C++ backend passes to the factory."""
        return cls(
            rank=int(data.get("rank", 0)),
            world_size=int(data.get("world_size", 1)),
            node_name=str(data.get("node_name", "")),
            nx=int(data.get("nx", 0)),
            ny=int(data.get("ny", 0)),
            num_global_cols=int(data.get("num_global_cols", 0)),
            col_gids=np.asarray(data.get("col_gids", []), dtype=np.int64),
            lat=np.asarray(data.get("lat", []), dtype=np.float64),
            lon=np.asarray(data.get("lon", []), dtype=np.float64),
        )

    @property
    def num_local_cols(self) -> int:
        return int(self.col_gids.size)

    @property
    def is_root(self) -> bool:
        return self.rank == 0

    def describe(self) -> str:
        return (
            f"rank {self.rank}/{self.world_size} on {self.node_name or '?'}, "
            f"{self.num_local_cols} of {self.num_global_cols} columns "
            f"on a {self.nx}x{self.ny} grid"
        )
