# `omega_buildnml`

Python package that builds Omega's runtime configuration (`omega.yml`) during
`case.setup`/`case.submit`.

Special care is taken to validate the input configuration files against
`configs/Default.yml`, so mistakes (typos, unknown options, options CIME
controls exclusively) are caught early rather than surfacing as a confusing
runtime error. `buildnml` only validates the mesh requested by the case
being built; CI (`validate_config.py`) validates every mesh and override
entry in the packaged configuration files. See render documentation for
more details.

> [!CAUTION]
> **Required dependencies are limited to the standard library and PyYAML**
>
> `netCDF4` is optional and enables user input file checks against the mesh.
> If it is unavailable, the checks are skipped with a warning.
> CIME does not provide its own environment, it just assumes the case's python
> environment it is being called from already has PyYAML available. In practice
> this is often handled by loading a python module on an HPC system that
> includes `PyYAML`.

## Modules

```
omega_buildnml
├── __init__.py   # public API is defined by contents of __all__
├── _types.py     # shared type aliases
├── config.py     # merges configuration layers and resolves mesh input files
├── input_check.py # checks user input files against the case mesh
├── read_write.py # reads/writes YAML config files, including packaged data
├── validate.py   # validates configuration layers
├── data/         # packaged YAML data (see data/README.md)
└── tests/        # unit tests
```

## Development environment

For now, we use the same `conda` [environment](../../dev-conda.txt) used for
linting Omega source code. The development environment primarily provides
`pytest` needed for running the unit tests. The actual code defined in this
package requires the python standard library and PyYAML. Tests of input file
checks additionally require netCDF4.

## User input files

Users may set `IOStreams.InitialState.Filename` and
`IOStreams.InitialVertCoord.Filename` in `user_nl_omega` to absolute file
paths. Other settings in these streams remain controlled by CIME. On startup,
the files are checked against the default case mesh when it and netCDF4 are
available; missing prerequisites produce a warning. Different vertical layer
counts also produce warnings. Ensure the initial state and vertical coordinate
files use compatible layers.

```bash
conda create -n omega_dev --file ../../dev-conda.txt
conda activate omega_dev
pre-commit install
```

## Testing

To validate the configuration files and run unit tests, do:
```bash
cd components/omega/cime_config
./validate_config.py
pytest omega_buildnml/tests
```

## TODO

- [ ] Validate `IOStreams` conditionally on what's set within a stream
  - For example: if `UseStartEnd` is true, both `StartTime` and `EndTime`
    must be provided.
- [ ] Use `difflib` to suggest close matches for unknown/misspelled options.
- [ ] Check that mesh input files exist in the remote input data database
  - Gate failures on whether database is reachable at all.
    (Prevents false positives when LCRC is down for monthly maintenance)
