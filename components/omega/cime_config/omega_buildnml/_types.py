from os import PathLike as OsPathLike
from typing import Any, Union

YamlMapping = dict[str, Any]
PathLike = Union[str, OsPathLike[str]]
