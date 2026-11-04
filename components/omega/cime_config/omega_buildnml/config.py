from collections.abc import Mapping
from copy import deepcopy
from pathlib import Path
from typing import Any

from ._types import PathLike, YamlMapping


def build_omega_config(
    defaults: YamlMapping,
    coupled_overrides: YamlMapping,
    mesh_overrides: YamlMapping,
    runtime_overrides: YamlMapping,
    user_overrides: YamlMapping,
) -> YamlMapping:
    """
    Build the Omega configuration dictionary.

    Parameters:
    -----------
    defaults : dict[str, Any]
        Default configuration values, loaded from configs/Default.yml.
    coupled_overrides : dict[str, Any]
        Coupled model overrides, loaded from cime_config/config_overrides.yaml
    mesh_overrides: dict[str, Any]
        Mesh specific overrides, loaded from cime_config/mesh_overrides.yaml
    runtime_overrides : dict[str, Any]
        Runtime specific overrides, based on CIME case configuration.
    user_overrides : dict[str, Any]
        User specified overrides, loaded from user_nl_omega.

    Returns:
    --------
    dict[str, Any]:
        Final Omega configuration dictionary.
    """
    if "Omega" in defaults:
        defaults = defaults["Omega"]

    config = deepcopy(defaults)
    config = _deep_merge(config, coupled_overrides)
    config = _deep_merge(config, mesh_overrides)
    config = _deep_merge(config, runtime_overrides)
    config = _deep_merge(config, user_overrides)

    return {"Omega": config}


def resolve_streams_files(
    input_files: YamlMapping, mesh_name: str, din_loc_root: PathLike
) -> dict[str, str]:
    """
    Resolve input filenames for all Omega IOStreams

    Uses mesh name and input_files.yaml to determine the correct input files
    for each IOStream.

    Parameters:
    -----------
    input_files : dict[str, str]
        Parsed content of cime_config/input_files.yaml.
    mesh_name : str
        CIME ocean grid name
    din_loc_root : Path
        Path to root of E3SM input data directory.

    Returns:
    --------
    dict[str, str]
        Mapping of Omega IOStream names to resolved input filenames.
    """
    meshes: YamlMapping = input_files["meshes"]

    mesh_dir = Path(din_loc_root) / "ocn" / "omega" / mesh_name
    streams_files = {}

    for input_group in meshes[mesh_name]["inputs"]:
        resolved_file_path = mesh_dir / input_group["file"]

        for stream in input_group["streams"]:
            streams_files[stream] = str(resolved_file_path)

    return streams_files


def build_runtime_overrides(
    calendar: str,
    continue_run: bool,
    case_name: str,
    streams_files: dict[str, str],
    defaults: YamlMapping,
) -> YamlMapping:
    """
    Build overrides from CIME case config vars only available after case.setup

    Parameters:
    -----------
    calendar : str
        CIME calendar type (i.e., "NO_LEAP", "GREGORIAN").
    continue_run : bool
        Whether to continue a previous run.
    case_name : str
        CIME case name.
    streams_files : dict[str, str]
        Resolved input filenames keyed by Omega IOStream name.
    defaults : dict[str, Any]
        Default configuration values, loaded from configs/Default.yml.
        Needed to automatically prepend case_name to analysis files

    Returns:
    --------
    dict[str, Any]
        Runtime overrides dictionary.
    """
    io_overrides = {
        stream_name: {"Filename": filename}
        for stream_name, filename in streams_files.items()
    }

    if continue_run:
        io_overrides["InitialState"]["FreqUnits"] = "Never"
        io_overrides["RestartRead"] = {"FreqUnits": "OnStartup"}
    else:
        io_overrides["InitialState"]["FreqUnits"] = "OnStartup"
        io_overrides["RestartRead"] = {"FreqUnits": "Never"}

    analysis_overrides = _prefixed_output_streams(
        defaults.get("Analysis", {}), case_name, "Analysis"
    )
    io_overrides.update(
        _prefixed_output_streams(
            defaults.get("IOStreams", {}),
            case_name,
            "IOStream",
            write_only=True,
        )
    )

    # restart filename format is dictated by CIME's case_st_archive script
    io_overrides["RestartWrite"] = {
        "Filename": f"{case_name}.omega.r.$Y-$M-$D-$S.nc"
    }

    overrides = {
        "TimeIntegration": {"CalendarType": _to_omega_calendar(calendar)},
        "IOStreams": io_overrides,
    }

    if analysis_overrides:
        overrides["Analysis"] = analysis_overrides

    return overrides


def _deep_merge(
    base: Mapping[str, Any],
    override: Mapping[str, Any],
) -> YamlMapping:
    """
    Recursively merge two mappings without modifying either input.

    Nested mappings are merged. All other values, including lists, are replaced
    by the override value.

    Parameters:
    -----------
    base (dict[str, Any])
        The base mapping
    override (dict[str, Any]):
        The mapping with the override values

    Returns:
        dict[str, Any]: Merged mapping
    """
    merged: YamlMapping = deepcopy(dict(base))

    for key, override_value in override.items():
        base_value = merged.get(key)

        if isinstance(base_value, Mapping) and isinstance(
            override_value, Mapping
        ):
            merged[key] = _deep_merge(base_value, override_value)
        else:
            merged[key] = deepcopy(override_value)

    return merged


def _to_omega_calendar(calendar: str) -> str:
    """
    Convert a CIME calendar string to corresponding Omega calendar string

    Args:
        calendar (str): CIME calendar string.

    Returns:
        str: Omega calendar string.
    """
    if calendar == "NO_LEAP":
        return "No Leap"
    elif calendar == "GREGORIAN":
        return "Gregorian"
    else:
        msg = f"Unsupported calendar type: {calendar}"
        raise ValueError(msg)


def _prefixed_output_streams(
    entries: Mapping[str, Any],
    case_name: str,
    item_type: str,
    write_only: bool = False,
) -> dict[str, dict[str, str]]:
    """
    Build filename overrides prefixed with ${case_name}.omega.<filename>

    Parameters:
    -----------
    entries : Mapping[str, Any]
        Dictionary of IOStream or Analysis configs.
    case_name : str
        CIME case name.
    item_type : str
        Name of section (e.g., 'Analysis', 'IOStream') for error messages.
    write_only : bool
        If True, only entries with Mode == 'write' are processed.

    Returns:
    --------
    dict[str, dict[str, str]]
        Mapping of entry names to filename override dictionaries.
    """
    overrides = {}
    for name, config in entries.items():
        if not isinstance(config, Mapping):
            err_msg = (
                f"Invalid entry for {item_type} '{name}': "
                f"expected a mapping, got {type(config).__name__}"
            )
            raise ValueError(err_msg)

        if write_only and config.get("Mode", "").lower() != "write":
            continue

        if "Filename" not in config:
            err_msg = f"Missing 'Filename' for {item_type} '{name}'"
            raise ValueError(err_msg)

        overrides[name] = {
            "Filename": f"{case_name}.omega.{config['Filename']}"
        }

    return overrides
