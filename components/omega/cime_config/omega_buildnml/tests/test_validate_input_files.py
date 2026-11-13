from copy import deepcopy

import pytest

from omega_buildnml.validate import validate_input_files_config


@pytest.fixture
def input_files():
    """A minimal, valid ``input_files.yaml`` configuration."""
    return {
        "meshes": {
            "Icos10": {
                "inputs": [
                    {
                        "file": "ocean.Icos10.nc",
                        "streams": [
                            "HorzMeshIn",
                            "InitialVertCoord",
                            "InitialState",
                        ],
                    }
                ]
            }
        }
    }


@pytest.fixture
def defaults():
    """A minimal stand-in for the parsed ``Default.yml`` configuration."""
    return {
        "IOStreams": {
            "HorzMeshIn": {},
            "InitialVertCoord": {},
            "InitialState": {},
            "Forcing": {},
        }
    }


@pytest.fixture
def mesh(input_files):
    """The single mesh entry of the valid configuration."""
    return input_files["meshes"]["Icos10"]


@pytest.fixture
def input_group(mesh):
    """The single input group of the valid configuration."""
    return mesh["inputs"][0]


def test_a_valid_configuration_is_returned(input_files, defaults):
    assert validate_input_files_config(input_files, defaults) == input_files


def test_a_valid_mesh_is_returned(input_files, defaults):
    validated = validate_input_files_config(
        input_files, defaults, mesh_name="Icos10"
    )

    assert validated == input_files


@pytest.mark.parametrize("config", [{}, [], "meshes", None])
def test_configurations_that_are_not_mappings_are_rejected(config, defaults):
    with pytest.raises(ValueError, match="empty or is not a mapping"):
        validate_input_files_config(config, defaults)


def test_unknown_top_level_keys_are_rejected(input_files, defaults):
    input_files["mesh"] = {}

    with pytest.raises(ValueError, match=r"Unknown top-level key\(s\): mesh"):
        validate_input_files_config(input_files, defaults)


@pytest.mark.parametrize("meshes", [{}, [], "Icos10", None])
def test_meshes_that_are_not_mappings_are_rejected(
    input_files, meshes, defaults
):
    input_files["meshes"] = meshes

    with pytest.raises(ValueError, match="`meshes` is missing, empty"):
        validate_input_files_config(input_files, defaults)


def test_an_unsupported_mesh_is_rejected(input_files, defaults):
    with pytest.raises(ValueError, match="Unsupported OCN_GRID"):
        validate_input_files_config(input_files, defaults, mesh_name="Icos30")


@pytest.mark.parametrize("entry", [{}, [], "inputs", None])
def test_mesh_entries_that_are_not_mappings_are_rejected(
    input_files, entry, defaults
):
    input_files["meshes"]["Icos10"] = entry

    with pytest.raises(ValueError, match="is empty or is not a mapping"):
        validate_input_files_config(input_files, defaults)


def test_unknown_keys_in_a_mesh_entry_are_rejected(
    input_files, mesh, defaults
):
    mesh["input"] = []

    with pytest.raises(ValueError, match=r"Unknown key\(s\) input"):
        validate_input_files_config(input_files, defaults)


@pytest.mark.parametrize("inputs", [{}, [], "file", None])
def test_inputs_that_are_not_lists_are_rejected(
    input_files, mesh, inputs, defaults
):
    mesh["inputs"] = inputs

    with pytest.raises(ValueError, match="Missing inputs"):
        validate_input_files_config(input_files, defaults)


def test_input_groups_that_are_not_mappings_are_rejected(
    input_files, mesh, defaults
):
    mesh["inputs"] = ["ocean.Icos10.nc"]

    with pytest.raises(ValueError, match="Input group 0 is not a mapping"):
        validate_input_files_config(input_files, defaults)


def test_unknown_keys_in_an_input_group_are_rejected(
    input_files, input_group, defaults
):
    input_group["stream"] = []

    with pytest.raises(ValueError, match=r"Unknown key\(s\) stream"):
        validate_input_files_config(input_files, defaults)


@pytest.mark.parametrize("file_name", ["", [], None])
def test_input_groups_without_a_file_are_rejected(
    input_files, input_group, file_name, defaults
):
    input_group["file"] = file_name

    with pytest.raises(ValueError, match="Missing file in input group 0"):
        validate_input_files_config(input_files, defaults)


@pytest.mark.parametrize("streams", [[], {}, "InitialState", None])
def test_input_groups_without_streams_are_rejected(
    input_files, input_group, streams, defaults
):
    input_group["streams"] = streams

    with pytest.raises(ValueError, match="Missing streams in input group 0"):
        validate_input_files_config(input_files, defaults)


@pytest.mark.parametrize("stream", ["", [], None])
def test_stream_names_that_are_not_strings_are_rejected(
    input_files, input_group, stream, defaults
):
    input_group["streams"].append(stream)

    with pytest.raises(ValueError, match="must be non-empty strings"):
        validate_input_files_config(input_files, defaults)


def test_unknown_stream_names_are_rejected(input_files, input_group, defaults):
    input_group["streams"].append("MyCustomStream")

    with pytest.raises(ValueError, match="Unknown IOStream 'MyCustomStream'"):
        validate_input_files_config(input_files, defaults)


def test_streams_assigned_more_than_once_are_rejected(
    input_files, mesh, defaults
):
    mesh["inputs"].append({"file": "state.nc", "streams": ["InitialState"]})

    with pytest.raises(ValueError, match="assigned more than once"):
        validate_input_files_config(input_files, defaults)


def test_streams_may_be_split_across_input_groups(input_files, mesh, defaults):
    mesh["inputs"] = [
        {"file": "mesh.nc", "streams": ["HorzMeshIn", "InitialVertCoord"]},
        {"file": "state.nc", "streams": ["InitialState"]},
    ]

    assert validate_input_files_config(input_files, defaults) == input_files


def test_missing_required_streams_are_rejected(
    input_files, input_group, defaults
):
    input_group["streams"].remove("InitialVertCoord")

    with pytest.raises(
        ValueError, match=r"Missing required IOStream\(s\) InitialVertCoord"
    ):
        validate_input_files_config(input_files, defaults)


def test_optional_streams_may_be_assigned(input_files, mesh, defaults):
    mesh["inputs"].append({"file": "forcing.nc", "streams": ["Forcing"]})

    assert validate_input_files_config(input_files, defaults) == input_files


def test_every_mesh_is_validated_when_no_mesh_is_given(input_files, defaults):
    broken = deepcopy(input_files["meshes"]["Icos10"])
    broken["inputs"][0]["streams"].remove("HorzMeshIn")
    input_files["meshes"]["Icos30"] = broken

    with pytest.raises(ValueError, match="mesh: Icos30"):
        validate_input_files_config(input_files, defaults)


def test_only_the_given_mesh_is_validated(input_files, defaults):
    input_files["meshes"]["Icos30"] = {"inputs": []}

    validated = validate_input_files_config(
        input_files, defaults, mesh_name="Icos10"
    )

    assert validated == input_files


def test_problems_are_reported_together(
    input_files, mesh, input_group, defaults
):
    mesh["input"] = []
    input_group["streams"].remove("InitialVertCoord")
    input_group["streams"].append(1)

    with pytest.raises(ValueError) as error:
        validate_input_files_config(input_files, defaults)

    reported = str(error.value)

    assert "Unknown key(s) input" in reported
    assert "must be non-empty strings" in reported
    assert "Missing required IOStream(s) InitialVertCoord" in reported
