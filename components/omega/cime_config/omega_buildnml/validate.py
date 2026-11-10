from typing import Optional

from ._types import YamlMapping

DATA_PATH = "components/omega/cime_config/omega_buildnml/data"

INPUT_FILES_PATH = f"{DATA_PATH}/input_files.yaml"

OVERRIDES_PATH = f"{DATA_PATH}/config_overrides.yaml"

DEFAULTS_PATH = "components/omega/configs/Default.yml"

#: IOStreams that every mesh must provide an input file for
REQUIRED_STREAMS = frozenset(
    {"HorzMeshIn", "InitialVertCoord", "InitialState"}
)

#: Keys allowed in a single mesh entry
MESH_KEYS = frozenset({"inputs"})

#: Keys allowed in a single entry of a mesh's ``inputs`` list
INPUT_GROUP_KEYS = frozenset({"file", "streams"})

#: Keys allowed at the top level of ``config_overrides.yaml``
OVERRIDES_KEYS = frozenset({"coupled", "meshes"})

#: Sections not strictly validated against the defaults. IOStreams entries may
#: define new streams, and have requirements that depend on the values passed,
#: so they need their own validation.
OPEN_SECTIONS = frozenset({"IOStreams"})

#: IOStreams wholly controlled by CIME and the coupler. The required streams
#: are staged through ``omega.input_data_list``, the restart streams are named
#: and scheduled by the coupler, and forcing is provided by the coupler.
BLOCKED_STREAMS = frozenset(
    REQUIRED_STREAMS | {"Forcing", "RestartRead", "RestartWrite"}
)

#: Config options set by CIME, which a user is not permitted to override.
#: Matched as prefixes, so naming a section blocks everything below it.
BLOCKED_OPTIONS = frozenset(
    {f"IOStreams.{stream}" for stream in BLOCKED_STREAMS}
    | {
        # start and stop types must be coupled
        "TimeIntegration.StartType",
        "TimeIntegration.StopType",
        # start time is provided by the coupler at runtime
        "TimeIntegration.StartTime",
        # stop criterion is ignored for coupled sims
        "TimeIntegration.StopCriterion",
        # calendar must agree with the CIME ``CALENDAR`` setting
        "TimeIntegration.CalendarType",
        # The parallel IO layout is owned by the driver (CIME/shr_pio) so it
        # stays consistent with the rest of the case. Set it with the
        # PIO_NUMTASKS_OCN, PIO_STRIDE_OCN, PIO_ROOT_OCN, PIO_REARRANGER_OCN
        # and PIO_TYPENAME_OCN case variables.
        "IO.IOTasks",
        "IO.IOStride",
        "IO.IOBaseTask",
        "IO.IORearranger",
        "IO.IODefaultFormat",
    }
)


def validate_input_files_config(
    input_files: YamlMapping,
    defaults: YamlMapping,
    mesh_name: Optional[str] = None,
) -> YamlMapping:
    """
    Validate the contents of the ``input_files.yaml`` configuration.

    All problems found are collected and reported together, rather than
    raising on the first one encountered.

    Parameters:
    -----------
    input_files : dict[str, Any]
        Parsed content of ``cime_config/omega_buildnml/data/input_files.yaml``
    defaults : dict[str, Any]
        Default configuration values, loaded from configs/Default.yml. Used
        to check that streams assigned an input file are streams Omega
        actually defines, since ``input_files.yaml`` only supplies a
        ``Filename`` for an existing ``IOStreams`` entry.
    mesh_name : str, optional
        The name of the mesh to validate. If not provided, all mesh entries
        will be validated.

    Returns:
    --------
    dict[str, Any]
        The validated configuration.

    Raises:
    -------
    ValueError
        If any required keys are missing or if any values are invalid.
    """
    if not isinstance(input_files, dict) or not input_files:
        err_msg = f"`{INPUT_FILES_PATH}` is empty or is not a mapping."
        raise ValueError(err_msg)

    unknown_keys = sorted(set(input_files) - {"meshes"})
    if unknown_keys:
        _raise(
            [f"Unknown top-level key(s): {', '.join(unknown_keys)}."],
            INPUT_FILES_PATH,
        )

    meshes: YamlMapping = input_files.get("meshes", {})
    if not isinstance(meshes, dict) or not meshes:
        _raise(
            ["`meshes` is missing, empty, or is not a mapping."],
            INPUT_FILES_PATH,
        )

    known_streams = frozenset(defaults.get("IOStreams", {}))

    if mesh_name is None:
        errors = []
        for name in meshes:
            errors.extend(
                _validate_input_files_entry(input_files, known_streams, name)
            )
        _raise(errors, INPUT_FILES_PATH)
        return input_files

    if mesh_name not in meshes:
        err_msg = (
            f"Unsupported OCN_GRID for Omega: {mesh_name}. \n"
            f"Could not find entry in `{INPUT_FILES_PATH}`"
        )
        raise ValueError(err_msg)

    _raise(
        _validate_input_files_entry(input_files, known_streams, mesh_name),
        INPUT_FILES_PATH,
    )

    return input_files


def validate_config_overrides(
    config_overrides: YamlMapping,
    input_files: YamlMapping,
    defaults: YamlMapping,
    mesh_name: Optional[str] = None,
) -> YamlMapping:
    """
    Validate the contents of the ``config_overrides.yaml`` configuration.

    All problems found are collected and reported together, rather than
    raising on the first one encountered.

    Mesh specific overrides are optional, so a mesh without an entry is not
    an error.

    Parameters:
    -----------
    config_overrides : dict[str, Any]
        Parsed content of
        ``cime_config/omega_buildnml/data/config_overrides.yaml``
    input_files : dict[str, Any]
        Parsed content of ``cime_config/omega_buildnml/data/input_files.yaml``,
        which defines the meshes Omega supports.
    defaults : dict[str, Any]
        Default configuration values, loaded from configs/Default.yml.
    mesh_name : str, optional
        The name of the mesh to validate. If not provided, all mesh entries
        will be validated.

    Returns:
    --------
    dict[str, Any]
        The validated configuration.

    Raises:
    -------
    ValueError
        If any required keys are missing or if any values are invalid.
    """
    if not isinstance(config_overrides, dict) or not config_overrides:
        err_msg = f"`{OVERRIDES_PATH}` is empty or is not a mapping."
        raise ValueError(err_msg)

    errors: list[str] = []

    unknown_keys = sorted(set(config_overrides) - OVERRIDES_KEYS)
    if unknown_keys:
        errors.append(f"Unknown top-level key(s): {', '.join(unknown_keys)}.")

    coupled = config_overrides.get("coupled")
    if not isinstance(coupled, dict) or not coupled:
        errors.append("`coupled` is missing, empty, or is not a mapping.")
    else:
        errors.extend(
            validate_overrides(coupled, defaults, "coupled overrides")
        )

    meshes: YamlMapping = config_overrides.get("meshes", {})
    if not isinstance(meshes, dict):
        errors.append("`meshes` is not a mapping.")
        _raise(errors, OVERRIDES_PATH)

    if mesh_name is None:
        for name in meshes:
            errors.extend(
                _validate_config_overrides_entry(
                    config_overrides, input_files, defaults, name
                )
            )
    elif mesh_name in meshes:
        errors.extend(
            _validate_config_overrides_entry(
                config_overrides, input_files, defaults, mesh_name
            )
        )

    _raise(errors, OVERRIDES_PATH)

    return config_overrides


def validate_user_overrides(
    user_overrides: YamlMapping, defaults: YamlMapping
) -> YamlMapping:
    """
    Validate the contents of a case's ``user_nl_omega`` file.

    All problems found are collected and reported together, rather than
    raising on the first one encountered.

    An empty ``user_nl_omega`` is not an error, it just means the user has
    not overridden anything.

    Parameters:
    -----------
    user_overrides : dict[str, Any]
        Parsed content of the case's ``user_nl_omega``.
    defaults : dict[str, Any]
        Default configuration values, loaded from configs/Default.yml.

    Returns:
    --------
    dict[str, Any]
        The validated user overrides.

    Raises:
    -------
    ValueError
        If any unknown or blocked options are set.
    """
    if not isinstance(user_overrides, dict):
        err_msg = "`user_nl_omega` is not a mapping."
        raise ValueError(err_msg)

    errors = validate_overrides(user_overrides, defaults, "user overrides")
    errors.extend(validate_blocked_options(user_overrides, "user overrides"))

    _raise(errors, "user_nl_omega")

    return user_overrides


def _raise(errors: list[str], config_path: str) -> None:
    """
    Raise a single ``ValueError`` describing all accumulated errors.

    Does nothing when ``errors`` is empty.

    Parameters:
    -----------
    errors : list[str]
        Error messages collected during validation.
    config_path : str
        Path of the configuration file the errors were found in.

    Raises:
    -------
    ValueError
        If ``errors`` is non-empty.
    """
    if not errors:
        return

    details = "\n".join(f"  - {error}" for error in errors)
    err_msg = (
        f"Invalid Omega configuration:\n{details}\n"
        f"Please check your setting in `{config_path}`"
    )
    raise ValueError(err_msg)


def _validate_input_files_entry(
    input_files: YamlMapping, known_streams: frozenset, mesh_name: str
) -> list[str]:
    """
    Validate that the specified mesh has a valid configuration in input_files.

    Parameters:
    -----------
    input_files : dict[str, Any]
        Parsed content of ``cime_config/omega_buildnml/data/input_files.yaml``
    known_streams : frozenset[str]
        Names of the ``IOStreams`` Omega defines in ``configs/Default.yml``.
    mesh_name : str
        The name of the mesh to validate.

    Returns:
    --------
    list[str]
        Error messages describing any problems found. Empty when the entry is
        valid.
    """
    meshes: YamlMapping = input_files["meshes"]
    mesh: YamlMapping = meshes[mesh_name]

    if not isinstance(mesh, dict) or not mesh:
        return [f"Mesh: {mesh_name} is empty or is not a mapping."]

    errors: list[str] = []
    streams_files = {}

    unknown_keys = sorted(set(mesh) - MESH_KEYS)
    if unknown_keys:
        errors.append(
            f"Unknown key(s) {', '.join(unknown_keys)} for mesh: {mesh_name}."
        )

    inputs = mesh.get("inputs")
    if not isinstance(inputs, list) or not inputs:
        errors.append(
            f"Missing inputs, or inputs is not a non-empty list, for "
            f"mesh: {mesh_name}."
        )
        return errors

    missing_msg = "Missing {key} in input group {index} for mesh: {mesh_name}."

    for index, input_group in enumerate(inputs):
        if not isinstance(input_group, dict):
            errors.append(
                f"Input group {index} is not a mapping for mesh: {mesh_name}."
            )
            continue

        unknown_keys = sorted(set(input_group) - INPUT_GROUP_KEYS)
        if unknown_keys:
            errors.append(
                f"Unknown key(s) {', '.join(unknown_keys)} in input group "
                f"{index} for mesh: {mesh_name}."
            )

        file_name = input_group.get("file")
        streams = input_group.get("streams")

        if not isinstance(file_name, str) or not file_name:
            errors.append(
                missing_msg.format(
                    key="file", index=index, mesh_name=mesh_name
                )
            )

        if not isinstance(streams, list) or not streams:
            errors.append(
                missing_msg.format(
                    key="streams", index=index, mesh_name=mesh_name
                )
            )
            continue

        for stream in streams:
            if not isinstance(stream, str) or not stream:
                errors.append(
                    f"Stream names must be non-empty strings, got "
                    f"'{stream}' in input group {index} for "
                    f"mesh: {mesh_name}."
                )
                continue

            if stream in streams_files:
                errors.append(
                    f"Stream '{stream}' is assigned more than once for "
                    f"mesh: {mesh_name}."
                )
                continue

            if stream not in known_streams:
                errors.append(
                    f"Unknown IOStream '{stream}' in input group {index} "
                    f"for mesh: {mesh_name}. Streams referenced in "
                    f"input_files.yaml must already be defined in "
                    f"Default.yml."
                )
                continue

            # just store filename; the full path will be resolved later
            # must store something to test for duplicates
            streams_files[stream] = str(file_name)

    missing_streams = REQUIRED_STREAMS - set(streams_files)
    if missing_streams:
        errors.append(
            f"Missing required IOStream(s) "
            f"{', '.join(sorted(missing_streams))} for mesh: {mesh_name}."
        )

    return errors


def _validate_config_overrides_entry(
    config_overrides: YamlMapping,
    input_files: YamlMapping,
    defaults: YamlMapping,
    mesh_name: str,
) -> list[str]:
    """
    Validate the overrides of a single mesh in config_overrides.

    Parameters:
    -----------
    config_overrides : dict[str, Any]
        Parsed content of
        ``cime_config/omega_buildnml/data/config_overrides.yaml``
    input_files : dict[str, Any]
        Parsed content of ``cime_config/omega_buildnml/data/input_files.yaml``,
        which defines the meshes Omega supports.
    defaults : dict[str, Any]
        Default configuration values, loaded from configs/Default.yml.
    mesh_name : str
        The name of the mesh to validate.

    Returns:
    --------
    list[str]
        Error messages describing any problems found. Empty when the entry is
        valid.
    """
    meshes: YamlMapping = config_overrides["meshes"]
    supported_meshes: YamlMapping = input_files.get("meshes", {})
    overrides: YamlMapping = meshes[mesh_name]

    errors: list[str] = []

    if not isinstance(overrides, dict) or not overrides:
        errors.append(f"Mesh: {mesh_name} is empty or is not a mapping.")
    else:
        errors.extend(
            validate_overrides(
                overrides, defaults, f"overrides for mesh: {mesh_name}"
            )
        )
        if "IOStreams" in overrides:
            errors.append(
                f"`IOStreams` is not permitted under mesh: {mesh_name}. "
                f"IOStreams shared by every mesh belong under `coupled`, "
                f"and case-specific IOStreams belong in `user_nl_omega`."
            )

    if mesh_name not in supported_meshes:
        errors.append(
            f"Unsupported mesh: {mesh_name}. Could not find entry in "
            f"`{INPUT_FILES_PATH}`."
        )

    return errors


def validate_overrides(
    overrides: YamlMapping, defaults: YamlMapping, source: str
) -> list[str]:
    """
    Validate that overrides only set options defined in Omega's defaults.

    Overrides that do not appear in the defaults would silently add new
    options, rather than overriding an existing one. Sections listed in
    ``OPEN_SECTIONS`` are not checked.

    Parameters:
    -----------
    overrides : dict[str, Any]
        Override options to validate.
    defaults : dict[str, Any]
        Default configuration values, loaded from configs/Default.yml.
    source : str
        Description of where the overrides came from, used in error messages.

    Returns:
    --------
    list[str]
        Error messages describing any problems found. Empty when the overrides
        are valid.
    """
    unknown_options = _unknown_override_options(overrides, defaults)

    if not unknown_options:
        return []

    return [f"Unknown option(s) {', '.join(unknown_options)} in {source}."]


def validate_blocked_options(overrides: YamlMapping, source: str) -> list[str]:
    """
    Validate that overrides do not set options controlled by CIME.

    Options listed in ``BLOCKED_OPTIONS`` are set from the case configuration,
    or by the coupler at runtime, so overriding them would either be silently
    discarded or leave the run inconsistent with the rest of the case.

    Parameters:
    -----------
    overrides : dict[str, Any]
        Override options to validate.
    source : str
        Description of where the overrides came from, used in error messages.

    Returns:
    --------
    list[str]
        Error messages describing any problems found. Empty when the overrides
        are valid.
    """
    blocked_options = _blocked_override_options(overrides)

    if not blocked_options:
        return []

    return [
        f"Option(s) {', '.join(blocked_options)} in {source} are set by CIME "
        f"and cannot be overridden."
    ]


def _unknown_override_options(
    overrides: YamlMapping, defaults: YamlMapping, prefix: str = ""
) -> list[str]:
    """
    Find override options that are not defined in Omega's defaults.

    Parameters:
    -----------
    overrides : dict[str, Any]
        Override options to validate.
    defaults : dict[str, Any]
        Default configuration values, loaded from configs/Default.yml.
    prefix : str, optional
        Dotted path of the parent section, used when recursing.

    Returns:
    --------
    list[str]
        Dotted paths of any options not defined in the defaults.
    """
    unknown_options: list[str] = []

    for key, value in overrides.items():
        option = f"{prefix}{key}"

        if option in OPEN_SECTIONS:
            continue

        if key not in defaults:
            unknown_options.append(option)
            continue

        default_value = defaults[key]

        if isinstance(value, dict) != isinstance(default_value, dict):
            unknown_options.append(option)
        elif isinstance(value, dict):
            unknown_options.extend(
                _unknown_override_options(value, default_value, f"{option}.")
            )

    return unknown_options


def _blocked_override_options(
    overrides: YamlMapping, prefix: str = ""
) -> list[str]:
    """
    Find override options that are controlled by CIME.

    Every level of the override tree is checked on the way down, so listing a
    section in ``BLOCKED_OPTIONS`` blocks all the options below it.

    Parameters:
    -----------
    overrides : dict[str, Any]
        Override options to validate.
    prefix : str, optional
        Dotted path of the parent section, used when recursing.

    Returns:
    --------
    list[str]
        Dotted paths of any options that are controlled by CIME.
    """
    blocked_options: list[str] = []

    for key, value in overrides.items():
        option = f"{prefix}{key}"

        if option in BLOCKED_OPTIONS:
            blocked_options.append(option)
        elif isinstance(value, dict):
            blocked_options.extend(
                _blocked_override_options(value, f"{option}.")
            )

    return blocked_options
