import pytest

from omega_buildnml.validate import (
    validate_blocked_options,
    validate_config_overrides,
    validate_overrides,
)


@pytest.fixture
def defaults():
    """A minimal stand in for Omega's default configuration."""
    return {
        "TimeIntegration": {
            "TimeStepper": "Forward-Backward",
            "TimeStep": "0000_00:30:00",
        },
        "Tendencies": {"SurfaceTracerRestoringEnable": False},
        "IOStreams": {
            "InitialState": {"Filename": "ocean.nc"},
            "History": {"Freq": 1, "FreqUnits": "months"},
        },
    }


@pytest.fixture
def input_files():
    """An ``input_files.yaml`` defining the meshes Omega supports."""
    return {"meshes": {"Icos10": {}, "Icos30": {}}}


@pytest.fixture
def config_overrides():
    """A minimal, valid ``config_overrides.yaml`` configuration."""
    return {
        "coupled": {
            "TimeIntegration": {"TimeStepper": "RungeKutta4"},
        },
        "meshes": {
            "Icos10": {"TimeIntegration": {"TimeStep": "0000_00:05:00"}},
        },
    }


def test_a_valid_configuration_is_returned(
    config_overrides, input_files, defaults
):
    validated = validate_config_overrides(
        config_overrides, input_files, defaults
    )

    assert validated == config_overrides


@pytest.mark.parametrize("config", [{}, [], "coupled", None])
def test_configurations_that_are_not_mappings_are_rejected(
    config, input_files, defaults
):
    with pytest.raises(ValueError, match="empty or is not a mapping"):
        validate_config_overrides(config, input_files, defaults)


def test_unknown_top_level_keys_are_rejected(
    config_overrides, input_files, defaults
):
    config_overrides["mesh"] = {}

    with pytest.raises(ValueError, match=r"Unknown top-level key\(s\): mesh"):
        validate_config_overrides(config_overrides, input_files, defaults)


@pytest.mark.parametrize("coupled", [{}, [], "TimeIntegration", None])
def test_coupled_overrides_are_required(
    config_overrides, input_files, defaults, coupled
):
    config_overrides["coupled"] = coupled

    with pytest.raises(ValueError, match="`coupled` is missing, empty"):
        validate_config_overrides(config_overrides, input_files, defaults)


def test_mesh_overrides_are_optional(config_overrides, input_files, defaults):
    del config_overrides["meshes"]

    validated = validate_config_overrides(
        config_overrides, input_files, defaults
    )

    assert validated == config_overrides


@pytest.mark.parametrize("meshes", [[], "Icos10"])
def test_meshes_that_are_not_mappings_are_rejected(
    config_overrides, input_files, defaults, meshes
):
    config_overrides["meshes"] = meshes

    with pytest.raises(ValueError, match="`meshes` is not a mapping"):
        validate_config_overrides(config_overrides, input_files, defaults)


@pytest.mark.parametrize("overrides", [{}, [], "TimeIntegration", None])
def test_mesh_entries_that_are_not_mappings_are_rejected(
    config_overrides, input_files, defaults, overrides
):
    config_overrides["meshes"]["Icos10"] = overrides

    with pytest.raises(ValueError, match="is empty or is not a mapping"):
        validate_config_overrides(config_overrides, input_files, defaults)


def test_overrides_for_unsupported_meshes_are_rejected(
    config_overrides, input_files, defaults
):
    config_overrides["meshes"]["Icos120"] = {
        "TimeIntegration": {"TimeStep": "0000_01:00:00"}
    }

    with pytest.raises(ValueError, match="Unsupported mesh: Icos120"):
        validate_config_overrides(config_overrides, input_files, defaults)


def test_unknown_coupled_overrides_are_rejected(
    config_overrides, input_files, defaults
):
    config_overrides["coupled"]["TimeIntegration"]["TimeSteper"] = "RK4"

    with pytest.raises(ValueError, match="in coupled overrides"):
        validate_config_overrides(config_overrides, input_files, defaults)


def test_unknown_mesh_overrides_are_rejected(
    config_overrides, input_files, defaults
):
    config_overrides["meshes"]["Icos10"]["TimeIntegration"]["Step"] = "0"

    with pytest.raises(ValueError, match="in overrides for mesh: Icos10"):
        validate_config_overrides(config_overrides, input_files, defaults)


def test_iostreams_under_a_mesh_entry_are_rejected(
    config_overrides, input_files, defaults
):
    config_overrides["meshes"]["Icos10"]["IOStreams"] = {
        "MyStream": {"Freq": 1}
    }

    with pytest.raises(
        ValueError, match="`IOStreams` is not permitted under mesh: Icos10"
    ):
        validate_config_overrides(config_overrides, input_files, defaults)


def test_iostreams_are_allowed_under_coupled(
    config_overrides, input_files, defaults
):
    config_overrides["coupled"]["IOStreams"] = {"MyStream": {"Freq": 1}}

    validated = validate_config_overrides(
        config_overrides, input_files, defaults
    )

    assert validated == config_overrides


def test_every_mesh_is_validated_when_no_mesh_is_given(
    config_overrides, input_files, defaults
):
    config_overrides["meshes"]["Icos30"] = {"TimeIntegration": {"Step": "0"}}

    with pytest.raises(ValueError, match="in overrides for mesh: Icos30"):
        validate_config_overrides(config_overrides, input_files, defaults)


def test_only_the_given_mesh_is_validated(
    config_overrides, input_files, defaults
):
    config_overrides["meshes"]["Icos30"] = {"TimeIntegration": {"Step": "0"}}

    validated = validate_config_overrides(
        config_overrides, input_files, defaults, mesh_name="Icos10"
    )

    assert validated == config_overrides


def test_problems_are_reported_together(
    config_overrides, input_files, defaults
):
    config_overrides["mesh"] = {}
    config_overrides["coupled"]["Tendencies"] = {"Restoring": True}
    config_overrides["meshes"]["Icos120"] = {"TimeIntegration": {}}

    with pytest.raises(ValueError) as error:
        validate_config_overrides(config_overrides, input_files, defaults)

    reported = str(error.value)

    assert "Unknown top-level key(s): mesh" in reported
    assert "Tendencies.Restoring" in reported
    assert "Unsupported mesh: Icos120" in reported


def test_overrides_of_known_options_are_valid(defaults):
    overrides = {"TimeIntegration": {"TimeStep": "0000_00:05:00"}}

    assert validate_overrides(overrides, defaults, "test") == []


def test_overrides_of_unknown_options_are_reported(defaults):
    overrides = {"TimeIntegration": {"TimeSteps": "0000_00:05:00"}}

    assert validate_overrides(overrides, defaults, "test") == [
        "Unknown option(s) TimeIntegration.TimeSteps in test."
    ]


def test_unknown_sections_are_reported(defaults):
    overrides = {"Tendancies": {"SurfaceTracerRestoringEnable": True}}

    assert validate_overrides(overrides, defaults, "test") == [
        "Unknown option(s) Tendancies in test."
    ]


def test_options_that_are_not_sections_are_reported(defaults):
    overrides = {"TimeIntegration": "0000_00:05:00"}

    assert validate_overrides(overrides, defaults, "test") == [
        "Unknown option(s) TimeIntegration in test."
    ]


def test_sections_that_are_not_options_are_reported(defaults):
    overrides = {"TimeIntegration": {"TimeStep": {"Value": "0000_00:05:00"}}}

    assert validate_overrides(overrides, defaults, "test") == [
        "Unknown option(s) TimeIntegration.TimeStep in test."
    ]


def test_new_streams_are_allowed(defaults):
    overrides = {"IOStreams": {"MyStream": {"Freq": 1}}}

    assert validate_overrides(overrides, defaults, "test") == []


def test_unknown_options_are_reported_together(defaults):
    overrides = {
        "TimeIntegration": {"TimeSteps": "0000_00:05:00"},
        "Tendancies": {"SurfaceTracerRestoringEnable": True},
    }

    assert validate_overrides(overrides, defaults, "test") == [
        "Unknown option(s) TimeIntegration.TimeSteps, Tendancies in test."
    ]


def test_options_that_are_not_blocked_are_allowed():
    overrides = {"TimeIntegration": {"TimeStep": "0000_00:05:00"}}

    assert validate_blocked_options(overrides, "test") == []


def test_blocked_options_are_reported():
    overrides = {"TimeIntegration": {"StartTime": "0001-01-01_00:00:00"}}

    assert validate_blocked_options(overrides, "test") == [
        "Option(s) TimeIntegration.StartTime in test are set by CIME and "
        "cannot be overridden."
    ]


def test_blocked_sections_are_reported_rather_than_their_options():
    overrides = {"IOStreams": {"RestartWrite": {"Precision": "single"}}}

    assert validate_blocked_options(overrides, "test") == [
        "Option(s) IOStreams.RestartWrite in test are set by CIME and "
        "cannot be overridden."
    ]
