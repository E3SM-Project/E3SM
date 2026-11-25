import pytest

from omega_buildnml.validate import (
    check_user_input_files,
    get_user_input_files,
)


def test_get_user_input_files():
    user_overrides = {
        "IOStreams": {
            "InitialState": {"Filename": "/a.nc"},
            "InitialVertCoord": {"Filename": "/b.nc"},
            "History": {"Filename": "hist"},
        }
    }

    assert get_user_input_files(user_overrides) == {
        "InitialState": "/a.nc",
        "InitialVertCoord": "/b.nc",
    }
    assert get_user_input_files({}) == {}


def test_existing_user_files_pass(tmp_path):
    ic = tmp_path / "ic.nc"
    ic.touch()

    check_user_input_files({})
    check_user_input_files(
        {"InitialState": str(ic), "InitialVertCoord": str(ic)}
    )


def test_missing_files_and_directories_are_rejected(tmp_path):
    with pytest.raises(ValueError) as exc:
        check_user_input_files(
            {
                "InitialState": str(tmp_path / "missing.nc"),
                "InitialVertCoord": str(tmp_path),
            }
        )

    assert "IOStreams.InitialState.Filename" in str(exc.value)
    assert "IOStreams.InitialVertCoord.Filename" in str(exc.value)
    assert "does not exist or is not a file" in str(exc.value)
