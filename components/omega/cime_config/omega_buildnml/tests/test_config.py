import pytest

from omega_buildnml.config import (
    build_omega_config,
    build_runtime_overrides,
)
from omega_buildnml.read_write import DEFAULT_CONFIG_PATH, read_default_config
from omega_buildnml.validate import BLOCKED_OPTIONS, validate_user_overrides


@pytest.fixture
def defaults():
    """A minimal stand in for Omega's default configuration."""
    return {
        "Omega": {
            "TimeIntegration": {
                "TimeStepper": "default",
                "TimeStep": "0000_00:30:00",
            },
            "IOStreams": {
                "History": {"Freq": 1, "FreqUnits": "months"},
            },
        }
    }


@pytest.fixture
def omega_defaults():
    """Omega's actual default configuration."""
    return read_default_config(DEFAULT_CONFIG_PATH)


def _time_stepper(source):
    """Build an override setting ``TimeStepper`` to ``source``."""
    return {"TimeIntegration": {"TimeStepper": source}}


def _nested_option(option, value):
    """Build a nested mapping from a dotted config path."""
    override = value

    for key in reversed(option.split(".")):
        override = {key: override}

    return override


def _build(defaults, **layers):
    """Build a config, defaulting any unspecified override layer to empty."""
    return build_omega_config(
        defaults=defaults,
        coupled_overrides=layers.get("coupled_overrides", {}),
        mesh_overrides=layers.get("mesh_overrides", {}),
        runtime_overrides=layers.get("runtime_overrides", {}),
        user_overrides=layers.get("user_overrides", {}),
    )


def test_defaults_are_used_when_nothing_is_overridden(defaults):
    config = _build(defaults)

    assert config["Omega"]["TimeIntegration"]["TimeStepper"] == "default"


def test_coupled_overrides_the_defaults(defaults):
    config = _build(defaults, coupled_overrides=_time_stepper("coupled"))

    assert config["Omega"]["TimeIntegration"]["TimeStepper"] == "coupled"


def test_mesh_overrides_the_coupled_overrides(defaults):
    config = _build(
        defaults,
        coupled_overrides=_time_stepper("coupled"),
        mesh_overrides=_time_stepper("mesh"),
    )

    assert config["Omega"]["TimeIntegration"]["TimeStepper"] == "mesh"


def test_runtime_overrides_the_mesh_overrides(defaults):
    config = _build(
        defaults,
        mesh_overrides=_time_stepper("mesh"),
        runtime_overrides=_time_stepper("runtime"),
    )

    assert config["Omega"]["TimeIntegration"]["TimeStepper"] == "runtime"


def test_user_overrides_the_runtime_overrides(defaults):
    config = _build(
        defaults,
        runtime_overrides=_time_stepper("runtime"),
        user_overrides=_time_stepper("user"),
    )

    assert config["Omega"]["TimeIntegration"]["TimeStepper"] == "user"


def test_user_overrides_are_applied_last(defaults):
    config = _build(
        defaults,
        coupled_overrides=_time_stepper("coupled"),
        mesh_overrides=_time_stepper("mesh"),
        runtime_overrides=_time_stepper("runtime"),
        user_overrides=_time_stepper("user"),
    )

    assert config["Omega"]["TimeIntegration"]["TimeStepper"] == "user"


def test_options_that_are_not_overridden_are_preserved(defaults):
    config = _build(defaults, user_overrides=_time_stepper("user"))

    time_integration = config["Omega"]["TimeIntegration"]

    assert time_integration["TimeStep"] == "0000_00:30:00"


def test_overrides_do_not_modify_their_inputs(defaults):
    user_overrides = _time_stepper("user")

    _build(defaults, user_overrides=user_overrides)

    assert defaults["Omega"]["TimeIntegration"]["TimeStepper"] == "default"
    assert user_overrides == _time_stepper("user")


def test_defaults_are_accepted_without_an_omega_section(defaults):
    config = _build(defaults["Omega"], user_overrides=_time_stepper("user"))

    assert config["Omega"]["TimeIntegration"]["TimeStepper"] == "user"


@pytest.mark.parametrize("blocked_option", sorted(BLOCKED_OPTIONS))
def test_blocked_options_are_rejected_before_they_are_merged(
    omega_defaults, blocked_option
):
    """
    Blocked options are enforced by validation, not by the merge order.

    User overrides are merged last, so a blocked option would win if it ever
    reached ``build_omega_config``.
    """
    user_overrides = _nested_option(blocked_option, "user")

    with pytest.raises(ValueError, match="cannot be overridden"):
        validate_user_overrides(user_overrides, omega_defaults)


@pytest.mark.parametrize("blocked_option", sorted(BLOCKED_OPTIONS))
def test_blocked_options_are_defined_in_the_defaults(
    omega_defaults, blocked_option
):
    """``BLOCKED_OPTIONS`` can drift from the options Omega defines."""
    option = omega_defaults

    for key in blocked_option.split("."):
        assert key in option, f"{blocked_option} is not in the defaults"
        option = option[key]


@pytest.fixture
def runtime_args():
    return {
        "calendar": "NO_LEAP",
        "continue_run": False,
        "case_name": "test_case",
        "streams_files": {"InitialState": "init.nc"},
    }


def test_runtime_overrides_prefixes_analysis_filenames(runtime_args):
    defaults = {
        "Analysis": {"GlobalStats": {"Filename": "global.stats.$Y.nc"}}
    }
    overrides = build_runtime_overrides(**runtime_args, defaults=defaults)
    assert (
        overrides["Analysis"]["GlobalStats"]["Filename"]
        == "test_case.omega.global.stats.$Y.nc"
    )


def test_runtime_overrides_omits_analysis_when_missing(runtime_args):
    overrides = build_runtime_overrides(**runtime_args, defaults={})
    assert "Analysis" not in overrides


def test_runtime_overrides_raises_missing_filename(runtime_args):
    defaults = {"Analysis": {"GlobalStats": {"Enable": True}}}
    with pytest.raises(ValueError, match="Missing 'Filename'"):
        build_runtime_overrides(**runtime_args, defaults=defaults)


def test_runtime_overrides_raises_non_mapping_member(runtime_args):
    defaults = {"Analysis": {"GlobalStats": "not_a_mapping"}}
    with pytest.raises(ValueError, match="expected a mapping"):
        build_runtime_overrides(**runtime_args, defaults=defaults)


def test_runtime_overrides_prefixes_write_iostreams_filenames(runtime_args):
    defaults = {
        "IOStreams": {
            "History": {"Mode": "write", "Filename": "history.$Y-$M"},
            "HorzMeshIn": {"Mode": "read", "Filename": "mesh.nc"},
        }
    }
    overrides = build_runtime_overrides(**runtime_args, defaults=defaults)
    assert (
        overrides["IOStreams"]["History"]["Filename"]
        == "test_case.omega.history.$Y-$M"
    )
    assert "HorzMeshIn" not in overrides["IOStreams"]


def test_runtime_overrides_raises_missing_filename_for_write_stream(
    runtime_args,
):
    defaults = {"IOStreams": {"History": {"Mode": "write"}}}
    with pytest.raises(ValueError, match="Missing 'Filename'"):
        build_runtime_overrides(**runtime_args, defaults=defaults)
