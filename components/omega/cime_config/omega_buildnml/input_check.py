"""
Check that input files set in ``user_nl_omega`` match the case's mesh.

A user may replace the ``InitialState`` and ``InitialVertCoord`` input files
(see ``USER_INPUT_STREAMS`` in ``validate.py``), but the mesh always comes
from ``input_files.yaml``. Each replacement file is compared against the
default mesh file for ``OCN_GRID``:

* errors: missing file, mesh dimensions that differ, variables the stream
  reads that are absent, and cell locations that differ at a sample of cells
* warnings: a different number of vertical layers, which may be intended

Mesh files carry no name or identifier, so two different meshes with the same
dimensions are told apart by the cell location spot check.
"""

import logging
import math
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Optional

from .validate import USER_INPUT_STREAMS, _raise

logger = logging.getLogger(__name__)

#: Variables Omega reads from each stream a user may replace
REQUIRED_VARIABLES = {
    "InitialState": (
        "NormalVelocity",
        "PseudoThickness",
        "Temperature",
        "Salinity",
    ),
    "InitialVertCoord": (
        "BottomGeomDepth",
        "RefPseudoThickness",
    ),
}

#: Dimensions that must match the mesh
MESH_DIMENSIONS = ("NCells", "NEdges", "NVertices", "MaxEdges")

#: Vertical dimension, which may differ from the default with a warning
LAYER_DIMENSION = "NVertLayers"

#: Cell location variables compared at a sample of cells
SPOT_CHECK_VARIABLES = ("LonCell", "LatCell")

#: Number of cells sampled by the spot check
SPOT_CHECK_POINTS = 16

#: Largest allowed difference in cell location, in radians (about 6 m)
SPOT_CHECK_TOLERANCE = 1.0e-6

SOURCE = "user_nl_omega"


def get_user_input_files(user_overrides: Mapping[str, Any]) -> dict[str, str]:
    """
    Find the input files a user has set in ``user_nl_omega``.

    Parameters:
    -----------
    user_overrides : dict[str, Any]
        Validated content of the case's ``user_nl_omega``.

    Returns:
    --------
    dict[str, str]
        Mapping of IOStream name to the user's input file.
    """
    io_streams = user_overrides.get("IOStreams") or {}
    user_files = {}

    for stream in USER_INPUT_STREAMS:
        config = io_streams.get(stream)
        if isinstance(config, Mapping) and "Filename" in config:
            user_files[stream] = config["Filename"]

    return user_files


def check_user_input_files(
    user_files: Mapping[str, str],
    default_files: Mapping[str, str],
    mesh_name: str,
) -> list[str]:
    """
    Check that user input files are consistent with the case's mesh.

    Missing files are always an error. The remaining checks need the
    ``netCDF4`` module and the default mesh file; if either is unavailable
    they are skipped with a warning. buildnml runs again at the start of
    ``case.run``, after CIME has downloaded the default files, so a check
    skipped at submission is normally done then.

    Parameters:
    -----------
    user_files : dict[str, str]
        Mapping of IOStream name to the user's input file.
    default_files : dict[str, str]
        Mapping of IOStream name to the default input file for the mesh, as
        returned by ``resolve_streams_files``.
    mesh_name : str
        CIME ocean grid name, used in messages.

    Returns:
    --------
    list[str]
        Warnings, which have also been logged.

    Raises:
    -------
    ValueError
        If any user input file is inconsistent with the mesh.
    """
    if not user_files:
        return []

    errors = [
        f"IOStreams.{stream}.Filename '{path}' does not exist or is not a "
        f"file."
        for stream, path in user_files.items()
        if not Path(path).is_file()
    ]
    _raise(errors, SOURCE)

    warnings: list[str] = []

    try:
        import netCDF4
    except ImportError:
        warnings.append(
            "The netCDF4 Python module is not available, so input files set "
            "in user_nl_omega were not checked against the mesh. Load the "
            "E3SM unified or CIME conda environment to enable the check."
        )
        return _log(warnings)

    mesh_file = default_files["HorzMeshIn"]
    if not Path(mesh_file).is_file():
        warnings.append(
            f"The default mesh file for {mesh_name}, '{mesh_file}', has not "
            f"been downloaded yet, so input files set in user_nl_omega were "
            f"not checked against it. They will be checked when buildnml "
            f"runs again at the start of case.run."
        )
        return _log(warnings)

    with netCDF4.Dataset(mesh_file) as mesh:
        n_cells = mesh.dimensions.get("NCells")
        samples = _sample_indices(len(n_cells) if n_cells is not None else 0)
        mesh_info = _read_info(mesh, sample=samples)

    effective_layers: dict[str, Optional[int]] = {}

    for stream in USER_INPUT_STREAMS:
        path = user_files.get(stream, default_files.get(stream))
        if path is None or not Path(path).is_file():
            continue
        with netCDF4.Dataset(path) as dataset:
            info = _read_info(dataset, sample=samples)
        effective_layers[stream] = info["dims"].get(LAYER_DIMENSION)

        if stream not in user_files:
            continue

        stream_errors, stream_warnings = _compare(
            stream, path, info, mesh_info, mesh_name
        )
        errors.extend(stream_errors)
        warnings.extend(stream_warnings)

        default_path = default_files.get(stream)
        if default_path and Path(default_path).is_file():
            with netCDF4.Dataset(default_path) as dataset:
                default_info = _read_info(dataset, sample=None)
            default_layers = default_info["dims"].get(LAYER_DIMENSION)
            layers = info["dims"].get(LAYER_DIMENSION)
            if (
                layers is not None
                and default_layers is not None
                and layers != default_layers
            ):
                warnings.append(
                    f"IOStreams.{stream}.Filename '{path}' has {layers} "
                    f"vertical layers, but the default {stream} file for "
                    f"{mesh_name} has {default_layers}. Make sure this is "
                    f"intended."
                )

    state_layers = effective_layers.get("InitialState")
    vert_layers = effective_layers.get("InitialVertCoord")
    if (
        state_layers is not None
        and vert_layers is not None
        and state_layers != vert_layers
    ):
        warnings.append(
            f"The InitialState file has {state_layers} vertical layers but "
            f"the InitialVertCoord file has {vert_layers}. Omega is unlikely "
            f"to run unless both are set to files with the same layers."
        )

    _raise(errors, SOURCE)

    return _log(warnings)


def _compare(
    stream: str,
    path: str,
    info: Mapping[str, Any],
    mesh_info: Mapping[str, Any],
    mesh_name: str,
) -> tuple[list[str], list[str]]:
    """
    Compare one user input file against the mesh.

    Returns:
    --------
    tuple[list[str], list[str]]
        Errors and warnings.
    """
    option = f"IOStreams.{stream}.Filename '{path}'"
    errors: list[str] = []
    warnings: list[str] = []

    if "NCells" not in info["dims"]:
        errors.append(f"{option} has no NCells dimension.")

    for dim in MESH_DIMENSIONS:
        size = info["dims"].get(dim)
        mesh_size = mesh_info["dims"].get(dim)
        if size is not None and mesh_size is not None and size != mesh_size:
            errors.append(
                f"{option} has {dim} = {size}, but mesh {mesh_name} has "
                f"{dim} = {mesh_size}."
            )

    missing = [
        name
        for name in REQUIRED_VARIABLES[stream]
        if name not in info["variables"]
    ]
    if missing:
        errors.append(
            f"{option} is missing variable(s) {', '.join(missing)} read by "
            f"the {stream} stream."
        )

    if errors:
        # the spot check is meaningless if the dimensions differ
        return errors, warnings

    for name in SPOT_CHECK_VARIABLES:
        values = info["samples"].get(name)
        mesh_values = mesh_info["samples"].get(name)
        if values is None or mesh_values is None:
            warnings.append(
                f"{option} or the mesh file has no {name}, so cell "
                f"locations were not compared."
            )
            continue
        largest = max(
            (
                _angle_difference(a, b)
                for a, b in zip(values, mesh_values, strict=True)
            ),
            default=0.0,
        )
        if largest > SPOT_CHECK_TOLERANCE:
            errors.append(
                f"{option} has {name} values that differ from mesh "
                f"{mesh_name} by up to {largest:.3g} radians, so it is on a "
                f"different mesh."
            )

    return errors, warnings


def _read_info(dataset: Any, sample: Optional[list[int]]) -> dict[str, Any]:
    """
    Read the dimensions, variable names and sampled cell locations of a file.
    """
    info: dict[str, Any] = {
        "dims": {name: len(dim) for name, dim in dataset.dimensions.items()},
        "variables": set(dataset.variables),
        "samples": {},
    }

    # skip sampling a file too small for the sample; _compare reports why
    if sample and info["dims"].get("NCells", 0) > max(sample):
        for name in SPOT_CHECK_VARIABLES:
            variable = dataset.variables.get(name)
            if variable is not None and variable.dimensions == ("NCells",):
                info["samples"][name] = [float(x) for x in variable[sample]]

    return info


def _sample_indices(n_cells: int) -> list[int]:
    """
    Evenly spaced cell indices, including the first and last cell.
    """
    if n_cells <= 0:
        return []
    if n_cells <= SPOT_CHECK_POINTS:
        return list(range(n_cells))
    step = (n_cells - 1) / (SPOT_CHECK_POINTS - 1)
    return sorted({round(i * step) for i in range(SPOT_CHECK_POINTS)})


def _angle_difference(a: float, b: float) -> float:
    """
    Absolute difference between two angles in radians, allowing for wrapping.
    """
    return abs((a - b + math.pi) % (2.0 * math.pi) - math.pi)


def _log(warnings: list[str]) -> list[str]:
    for warning in warnings:
        logger.warning("Omega: %s", warning)
    return warnings
