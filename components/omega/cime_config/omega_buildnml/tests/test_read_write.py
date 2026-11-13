import pytest
import yaml

from omega_buildnml.read_write import (
    _read_yaml_file,
    _unwrap_omega_section,
    read_user_overrides,
)


@pytest.fixture
def user_nl(tmp_path):
    """Write a mapping to a user_nl_omega file, returning the path."""

    def _write(overrides):
        path = tmp_path / "user_nl_omega"

        with path.open("w", encoding="utf-8") as f:
            yaml.safe_dump(overrides, f)

        return path

    return _write


@pytest.fixture
def malformed_yaml(tmp_path):
    """
    Write a YAML file that a mapping cannot express.

    Duplicate keys, comments, and indentation mistakes are all lost when a
    mapping is dumped, so they have to be written out as text.
    """

    def _write(contents):
        path = tmp_path / "sample.yaml"
        path.write_text(contents, encoding="utf-8")

        return path

    return _write


def _read_yaml(path):
    """Read a YAML file the way the package reads its config files."""
    with path.open("r", encoding="utf-8") as f:
        return _read_yaml_file(f)


def test_mappings_without_duplicates_are_read(user_nl):
    overrides = {
        "TimeIntegration": {
            "TimeStep": "0000_00:30:00",
            "TimeStepper": "Forward-Backward",
        }
    }

    assert _read_yaml(user_nl(overrides)) == overrides


def test_duplicate_top_level_keys_are_rejected(malformed_yaml):
    path = malformed_yaml(
        "TimeIntegration:\n  TimeStep: a\nTimeIntegration:\n  TimeStep: b\n"
    )

    with pytest.raises(ValueError, match="Duplicate key 'TimeIntegration'"):
        _read_yaml(path)


def test_duplicate_nested_keys_are_rejected(malformed_yaml):
    path = malformed_yaml("TimeIntegration:\n  TimeStep: a\n  TimeStep: b\n")

    with pytest.raises(ValueError, match="Duplicate key 'TimeStep'"):
        _read_yaml(path)


def test_duplicate_keys_within_a_list_are_rejected(malformed_yaml):
    path = malformed_yaml("inputs:\n  - file: a.nc\n    file: b.nc\n")

    with pytest.raises(ValueError, match="Duplicate key 'file'"):
        _read_yaml(path)


def test_duplicate_key_errors_report_where_they_were_found(malformed_yaml):
    path = malformed_yaml("Tracers:\n  Base: a\n  Base: b\n")

    with pytest.raises(ValueError, match="sample.yaml"):
        _read_yaml(path)


def test_overrides_are_returned_when_there_is_no_omega_section():
    overrides = {"TimeIntegration": {"TimeStep": "0000_00:30:00"}}

    assert _unwrap_omega_section(overrides) == overrides


def test_the_omega_section_is_unwrapped():
    overrides = {"TimeIntegration": {"TimeStep": "0000_00:30:00"}}

    assert _unwrap_omega_section({"Omega": overrides}) == overrides


def test_an_empty_omega_section_is_not_an_error():
    assert _unwrap_omega_section({"Omega": None}) == {}


def test_keys_alongside_the_omega_section_are_rejected():
    overrides = {
        "Omega": {"IOStreams": {}},
        "TimeIntegration": {"TimeStep": "0000_00:30:00"},
    }

    with pytest.raises(ValueError, match="only top-level key"):
        _unwrap_omega_section(overrides)


def test_an_empty_user_nl_omega_is_not_an_error(malformed_yaml):
    assert read_user_overrides(malformed_yaml("")) == {}


def test_a_user_nl_omega_with_only_comments_is_not_an_error(malformed_yaml):
    assert read_user_overrides(malformed_yaml("# no overrides here\n")) == {}


def test_user_overrides_are_read_without_an_omega_section(user_nl):
    overrides = {"TimeIntegration": {"TimeStep": "0000_00:10:00"}}

    assert read_user_overrides(user_nl(overrides)) == overrides


def test_user_overrides_are_read_with_an_omega_section(user_nl):
    overrides = {"TimeIntegration": {"TimeStep": "0000_00:10:00"}}

    assert read_user_overrides(user_nl({"Omega": overrides})) == overrides


def test_unknown_user_overrides_are_rejected(user_nl):
    overrides = {"TimeIntegration": {"TimeStepp": "0000_00:10:00"}}

    with pytest.raises(ValueError, match="Unknown option"):
        read_user_overrides(user_nl(overrides))


def test_blocked_user_overrides_are_rejected(user_nl):
    overrides = {"IOStreams": {"RestartWrite": {"Precision": "single"}}}

    with pytest.raises(ValueError, match="cannot be overridden"):
        read_user_overrides(user_nl(overrides))


def test_user_overrides_that_are_not_a_mapping_are_rejected(user_nl):
    with pytest.raises(ValueError, match="is not a mapping"):
        read_user_overrides(user_nl(["TimeIntegration", "IOStreams"]))


def test_a_missing_user_nl_omega_is_reported(tmp_path):
    with pytest.raises(FileNotFoundError, match="does not exist"):
        read_user_overrides(tmp_path / "user_nl_omega")
