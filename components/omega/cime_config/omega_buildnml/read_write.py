from importlib import resources
from pathlib import Path
from typing import IO, Any, Hashable, Optional

import yaml

from ._types import PathLike, YamlMapping
from .validate import (
    validate_config_overrides,
    validate_input_files_config,
    validate_user_overrides,
)

#: Path to Omega's default configuration file, relative to this package
DEFAULT_CONFIG_PATH = Path(__file__).parents[2] / "configs" / "Default.yml"


def read_default_config(path: PathLike) -> YamlMapping:
    """
    Read Omega's default configuration from a YAML file on disk.

    Parameters:
    -----------
    path : PathLike
        Path to default config file (i.e. components/omega/config/Defaults.yml)

    Returns:
    --------
    dict[str, Any]
        The read configuration as a dictionary.
    """
    path = Path(path)

    if not path.is_file():
        err_msg = f"{path} does not exist or is not a file"
        raise FileNotFoundError(err_msg)

    with path.open("r", encoding="utf-8") as f:
        config = _read_yaml_file(f)

    if "Omega" not in config:
        err_msg = f"{path} does not contain a top-level 'Omega' section"
        raise ValueError(err_msg)

    # TODO: check that the config is valid (e.g., required keys are present)
    return config["Omega"]


def read_input_files_config(mesh_name: Optional[str] = None) -> YamlMapping:
    """
    Read the input_files.yaml configuration file from the package resources.

    Parameters:
    -----------
    mesh_name : str, optional
        The name of the mesh to validate. If not provided, all mesh entries
        will be validated.

    Returns:
    --------
    dict[str, Any]
        The read configuration as a dictionary.
    """
    input_files = _read_packaged_yaml_file("input_files.yaml")
    defaults = read_default_config(DEFAULT_CONFIG_PATH)

    return validate_input_files_config(
        input_files, defaults, mesh_name=mesh_name
    )


def read_config_overrides(mesh_name: Optional[str] = None) -> YamlMapping:
    """
    Read config_overrides.yaml configuration file from the package resources.

    Parameters:
    -----------
    mesh_name : str, optional
        The name of the mesh to validate. If not provided, all mesh entries
        will be validated.

    Returns:
    --------
    dict[str, Any]
        The read configuration as a dictionary.
    """
    config_overrides = _read_packaged_yaml_file("config_overrides.yaml")
    input_files = read_input_files_config(mesh_name=mesh_name)
    defaults = read_default_config(DEFAULT_CONFIG_PATH)

    return validate_config_overrides(
        config_overrides, input_files, defaults, mesh_name=mesh_name
    )


def read_user_overrides(path: PathLike) -> YamlMapping:
    """
    Read a case's user_nl_omega file.

    The overrides are expected to be a YAML snippet mirroring the structure
    of Omega's default configuration. An empty file is not an error, it just
    means the user has not overridden anything.

    Parameters:
    -----------
    path : PathLike
        Path to the case's ``user_nl_omega`` file.

    Returns:
    --------
    dict[str, Any]
        The read overrides as a dictionary.
    """
    path = Path(path)

    if not path.is_file():
        err_msg = f"{path} does not exist or is not a file"
        raise FileNotFoundError(err_msg)

    with path.open("r", encoding="utf-8") as f:
        user_overrides = _read_yaml_file(f)

    # a user_nl_omega with no overrides in it parses as ``None``
    if user_overrides is None:
        return {}

    if not isinstance(user_overrides, dict):
        err_msg = f"{path} is not a mapping."
        raise ValueError(err_msg)

    user_overrides = _unwrap_omega_section(user_overrides)
    defaults = read_default_config(DEFAULT_CONFIG_PATH)

    return validate_user_overrides(user_overrides, defaults)


def write_yaml_mapping(mapping: YamlMapping, file_path: PathLike) -> None:
    """
    Write a mapping to a YAML file.

    Parameters:
    -----------
    mapping : dict[str, Any]
        The mapping to write to the YAML file.
    file_path : PathLike
        The path to the output YAML file.

    Returns:
    --------
    None
    """
    with Path(file_path).open("w", encoding="utf-8") as f:
        yaml.safe_dump(mapping, f, sort_keys=False, default_flow_style=False)


def write_input_data_list(
    streams_files: dict[str, str], casebuild: PathLike
) -> None:
    """
    Build omega.input_data_list

    Enables automatic retrieval of missing input files

    Parameters
    ----------
    streams_files : dict[str, str]
        Dict mapping stream names to their corresponding input file paths.
    casebuild : Path
        Path to the case build directory (i.e. CASEBUILD)

    Returns
    -------
    None
    """

    unique_files = list(dict.fromkeys(streams_files.values()))

    input_data_list = [
        f"omega_input_{index} = {filename}"
        for index, filename in enumerate(unique_files, start=1)
    ]

    path = Path(casebuild) / "omega.input_data_list"
    with path.open("w", encoding="utf-8") as f:
        for input_file in input_data_list:
            f.write(f"{input_file}\n")


class _UniqueKeyLoader(yaml.SafeLoader):
    """
    YAML loader that rejects mappings containing duplicate keys.

    ``yaml.SafeLoader`` silently keeps the last value when a key is repeated,
    so a copy/pasted or misspelled entry would quietly override an earlier
    one instead of being reported.
    """

    def construct_mapping(
        self, node: yaml.MappingNode, deep: bool = False
    ) -> dict[Hashable, Any]:
        keys = set()

        for key_node, _ in node.value:
            key = self.construct_object(key_node, deep=deep)

            if key in keys:
                err_msg = (
                    f"Duplicate key {key!r} found:\n{key_node.start_mark}"
                )
                raise ValueError(err_msg)

            keys.add(key)

        return super().construct_mapping(node, deep=deep)


def _unwrap_omega_section(config: YamlMapping) -> YamlMapping:
    """
    Return the ``Omega`` section of a config, if one is present.

    Users are expected to wrap their overrides in a top-level ``Omega`` key,
    matching Omega's default configuration, but the key is optional. When it
    is present it must be the only top-level key, otherwise under indented
    sections would be silently dropped.

    Parameters
    ----------
    config : dict[str, Any]
        Parsed content of a YAML file.

    Returns:
    --------
    dict[str, Any]
        The ``Omega`` section, or the config itself when there isn't one.
    """
    if "Omega" not in config:
        return config

    siblings = sorted(set(config) - {"Omega"})
    if siblings:
        err_msg = (
            f"`Omega` must be the only top-level key. Found "
            f"{', '.join(siblings)} alongside it, please check the "
            f"indentation of your overrides."
        )
        raise ValueError(err_msg)

    omega_section: YamlMapping = config["Omega"]

    # an ``Omega`` key with nothing under it parses as ``None``
    if omega_section is None:
        return {}

    return omega_section


def _read_yaml_file(f: IO[str]) -> YamlMapping:
    """
    Read a YAML mapping from a file.

    Parameters
    ----------
    f : TextIO
        File-like object to read the YAML mapping from.
    Returns:
    -------
    dict[str, Any]
        Read YAML mapping.
    """

    return yaml.load(f, Loader=_UniqueKeyLoader)


def _read_packaged_yaml_file(file_name: str) -> YamlMapping:
    """
    Read a YAML mapping packaged within omega_buildnml/data

    Parameters
    ----------
    file_path : str
        Name of the packaged YAML file to read.

    Returns:
    --------
    dict[str, Any]
        Read YAML mapping.
    """
    resource = resources.files("omega_buildnml.data").joinpath(file_name)

    if not resource.is_file():
        err_msg = f"Packaged configuration file '{file_name}' not found."
        raise FileNotFoundError(err_msg)

    with resource.open("r", encoding="utf-8") as f:
        return _read_yaml_file(f)
