import pytest

from omega_buildnml.validate import validate_user_overrides


@pytest.fixture
def defaults():
    """A minimal stand in for Omega's default configuration."""
    return {
        "TimeIntegration": {
            "TimeStepper": "Forward-Backward",
            "TimeStep": "0000_00:30:00",
            "StartTime": "0001-01-01_00:00:00",
        },
        "Tendencies": {"SurfaceTracerRestoringEnable": False},
        "IO": {
            "IOTasks": 1,
            "IOStride": 1,
            "IOBaseTask": 0,
            "IORearranger": "box",
            "IODefaultFormat": "pnetcdf",
        },
        "IOStreams": {
            "InitialState": {"Filename": "ocean.nc"},
            "History": {"Freq": 1, "FreqUnits": "months"},
        },
    }


def test_an_empty_configuration_is_returned(defaults):
    assert validate_user_overrides({}, defaults) == {}


def test_a_valid_configuration_is_returned(defaults):
    user_overrides = {"Tendencies": {"SurfaceTracerRestoringEnable": True}}

    validated = validate_user_overrides(user_overrides, defaults)

    assert validated == user_overrides


@pytest.mark.parametrize("config", [[], "TimeIntegration", None])
def test_configurations_that_are_not_mappings_are_rejected(config, defaults):
    with pytest.raises(ValueError, match="is not a mapping"):
        validate_user_overrides(config, defaults)


def test_unknown_options_are_rejected(defaults):
    user_overrides = {"TimeIntegration": {"TimeSteps": "0000_00:05:00"}}

    with pytest.raises(ValueError, match="Unknown option"):
        validate_user_overrides(user_overrides, defaults)


def test_blocked_options_are_rejected(defaults):
    user_overrides = {"TimeIntegration": {"StartTime": "0002-01-01_00:00:00"}}

    with pytest.raises(ValueError, match="cannot be overridden"):
        validate_user_overrides(user_overrides, defaults)


def test_custom_iostreams_are_allowed(defaults):
    """
    A user can add a wholly new IOStream in ``user_nl_omega``, so long as
    the stream itself is a valid configuration (e.g. for debugging).

    Regression test for a bug where the now-removed ``KNOWN_STREAMS`` list
    was mistakenly thought to gate this; ``IOStreams`` is an
    ``OPEN_SECTIONS`` entry in :func:`validate_overrides`, so custom stream
    names were never actually blocked here.
    """
    user_overrides = {
        "IOStreams": {
            "MyCustomHiFreq": {
                "Filename": "ocn.hifreq.$Y-$M",
                "Mode": "write",
                "Freq": 1,
                "FreqUnits": "second",
                "Contents": ["State"],
            }
        }
    }

    validated = validate_user_overrides(user_overrides, defaults)

    assert validated == user_overrides


def test_overriding_an_existing_stream_is_allowed(defaults):
    user_overrides = {"IOStreams": {"History": {"Freq": 5}}}

    validated = validate_user_overrides(user_overrides, defaults)

    assert validated == user_overrides


_DRIVER_OWNED_IO_OPTIONS = {
    "IOTasks": 8,
    "IOStride": 2,
    "IOBaseTask": 4,
    "IORearranger": "subset",
    "IODefaultFormat": "netcdf3",
}


@pytest.mark.parametrize("option", list(_DRIVER_OWNED_IO_OPTIONS))
def test_driver_owned_io_options_are_rejected(option, defaults):
    """
    The parallel IO layout is owned by the driver (CIME/shr_pio), so a user
    may not override it in ``user_nl_omega``.
    """
    user_overrides = {"IO": {option: _DRIVER_OWNED_IO_OPTIONS[option]}}

    with pytest.raises(ValueError, match="cannot be overridden"):
        validate_user_overrides(user_overrides, defaults)
