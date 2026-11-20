import math

import pytest

from omega_buildnml.input_check import (
    REQUIRED_VARIABLES,
    check_user_input_files,
    get_user_input_files,
)

netCDF4 = pytest.importorskip("netCDF4")

N_CELLS = 40
N_EDGES = 120
N_VERTICES = 80
N_LAYERS = 10


def _write_mesh_file(
    path,
    n_cells=N_CELLS,
    n_edges=N_EDGES,
    n_layers=N_LAYERS,
    lon_shift=0.0,
    variables=None,
    lon_lat=True,
):
    """Write a small stand in for an Omega mesh/initial condition file."""
    if variables is None:
        variables = (
            REQUIRED_VARIABLES["InitialState"]
            + REQUIRED_VARIABLES["InitialVertCoord"]
        )

    with netCDF4.Dataset(path, "w", format="NETCDF3_64BIT_DATA") as ds:
        ds.createDimension("NCells", n_cells)
        ds.createDimension("NEdges", n_edges)
        ds.createDimension("NVertices", N_VERTICES)
        ds.createDimension("MaxEdges", 7)
        ds.createDimension("NVertLayers", n_layers)
        ds.createDimension("time", 1)

        if lon_lat:
            lon = ds.createVariable("LonCell", "f8", ("NCells",))
            lat = ds.createVariable("LatCell", "f8", ("NCells",))
            lon[:] = [
                (2.0 * math.pi * i / n_cells + lon_shift) % (2.0 * math.pi)
                for i in range(n_cells)
            ]
            lat[:] = [
                -0.5 * math.pi + math.pi * i / n_cells for i in range(n_cells)
            ]

        for name in variables:
            ds.createVariable(name, "f8", ("NCells",))

    return str(path)


@pytest.fixture
def default_files(tmp_path):
    mesh = _write_mesh_file(tmp_path / "mesh.nc")
    return {
        "HorzMeshIn": mesh,
        "InitialVertCoord": mesh,
        "InitialState": mesh,
    }


def test_get_user_input_files():
    user_overrides = {
        "IOStreams": {
            "InitialState": {"Filename": "/a.nc"},
            "History": {"Filename": "hist"},
        }
    }

    assert get_user_input_files(user_overrides) == {"InitialState": "/a.nc"}
    assert get_user_input_files({}) == {}


def test_no_user_files_is_not_checked(default_files):
    assert check_user_input_files({}, default_files, "mesh") == []


def test_a_matching_file_passes(tmp_path, default_files):
    ic = _write_mesh_file(tmp_path / "ic.nc")

    warnings = check_user_input_files(
        {"InitialState": ic}, default_files, "mesh"
    )

    assert warnings == []


def test_a_missing_file_is_rejected(tmp_path, default_files):
    with pytest.raises(ValueError, match="does not exist"):
        check_user_input_files(
            {"InitialState": str(tmp_path / "nope.nc")}, default_files, "mesh"
        )


@pytest.mark.parametrize(
    ("kwargs", "dim"),
    [
        ({"n_cells": N_CELLS + 1}, "NCells"),
        ({"n_edges": N_EDGES - 1}, "NEdges"),
    ],
)
def test_a_file_on_a_different_mesh_is_rejected(
    tmp_path, default_files, kwargs, dim
):
    ic = _write_mesh_file(tmp_path / "ic.nc", **kwargs)

    with pytest.raises(ValueError, match=rf"has {dim} = "):
        check_user_input_files({"InitialState": ic}, default_files, "mesh")


def test_a_smaller_mesh_is_rejected_not_crashed(tmp_path, default_files):
    ic = _write_mesh_file(tmp_path / "ic.nc", n_cells=5)

    with pytest.raises(ValueError, match="has NCells = 5"):
        check_user_input_files({"InitialState": ic}, default_files, "mesh")


def test_different_cell_locations_are_rejected(tmp_path, default_files):
    ic = _write_mesh_file(tmp_path / "ic.nc", lon_shift=0.01)

    with pytest.raises(ValueError, match="LonCell values that differ"):
        check_user_input_files({"InitialState": ic}, default_files, "mesh")


def test_longitude_wrapping_is_not_a_difference(tmp_path, default_files):
    ic = _write_mesh_file(tmp_path / "ic.nc", lon_shift=2.0 * math.pi)

    assert (
        check_user_input_files({"InitialState": ic}, default_files, "mesh")
        == []
    )


def test_missing_variables_are_rejected(tmp_path, default_files):
    ic = _write_mesh_file(
        tmp_path / "ic.nc", variables=("Temperature", "Salinity")
    )

    with pytest.raises(ValueError, match="NormalVelocity, PseudoThickness"):
        check_user_input_files({"InitialState": ic}, default_files, "mesh")


def test_vertical_coordinate_optional_variables_are_not_required(
    tmp_path, default_files
):
    ic = _write_mesh_file(
        tmp_path / "ic.nc", variables=REQUIRED_VARIABLES["InitialVertCoord"]
    )

    assert (
        check_user_input_files({"InitialVertCoord": ic}, default_files, "mesh")
        == []
    )


def test_missing_required_vertical_coordinate_variable_is_rejected(
    tmp_path, default_files
):
    ic = _write_mesh_file(tmp_path / "ic.nc", variables=("BottomGeomDepth",))

    with pytest.raises(ValueError, match="RefPseudoThickness"):
        check_user_input_files({"InitialVertCoord": ic}, default_files, "mesh")


def test_missing_cell_locations_warn(tmp_path, default_files):
    ic = _write_mesh_file(tmp_path / "ic.nc", lon_lat=False)

    warnings = check_user_input_files(
        {"InitialState": ic}, default_files, "mesh"
    )

    assert any("not compared" in warning for warning in warnings)


def test_a_different_layer_count_warns(tmp_path, default_files):
    ic = _write_mesh_file(tmp_path / "ic.nc", n_layers=N_LAYERS + 5)

    warnings = check_user_input_files(
        {"InitialState": ic}, default_files, "mesh"
    )

    assert any("vertical layers" in warning for warning in warnings)
    assert any("InitialVertCoord file has" in warning for warning in warnings)


def test_consistent_new_layers_only_warn_once_per_stream(
    tmp_path, default_files
):
    ic = _write_mesh_file(tmp_path / "ic.nc", n_layers=N_LAYERS + 5)

    warnings = check_user_input_files(
        {"InitialState": ic, "InitialVertCoord": ic}, default_files, "mesh"
    )

    assert len(warnings) == 2
    assert not any("InitialVertCoord file has" in w for w in warnings)


def test_an_undownloaded_mesh_skips_the_check(tmp_path):
    ic = _write_mesh_file(tmp_path / "ic.nc")
    missing = str(tmp_path / "not_yet_downloaded.nc")
    default_files = {
        "HorzMeshIn": missing,
        "InitialVertCoord": missing,
        "InitialState": missing,
    }

    warnings = check_user_input_files(
        {"InitialState": ic}, default_files, "mesh"
    )

    assert any("has not been downloaded" in warning for warning in warnings)
