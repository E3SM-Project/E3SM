(omega-dev-buildnml)=

# CIME `buildnml` and Configuration Validation

`components/omega/cime_config/buildnml` generates a case's `omega.yml` by
layering Omega's defaults, coupled/mesh overrides, runtime overrides derived
from the case, and `user_nl_omega` (see {ref}`omega-user-config`). This
page covers the Python package behind `buildnml`, how its configuration is
validated, and how to add support for a new mesh.

## `omega_buildnml` package

`components/omega/cime_config/omega_buildnml/` implements the merging and
validation logic used by `buildnml`:

- `read_write.py` reads/writes YAML config files and packaged data.
- `config.py` merges the configuration layers, resolves mesh input files, and
  prefixes analysis filenames with the case name.
- `validate.py` validates overrides and `IOStreams` against `Default.yml`.
- `data/input_files.yaml` maps each supported mesh to its input file(s).
- `data/config_overrides.yaml` holds the coupled and mesh-specific overrides.

Validation runs whenever these files are read, so a bad edit fails fast at
`case.setup` rather than surfacing as a confusing runtime error. Stream
names are checked dynamically against `Default.yml` rather than a
hardcoded list: the `coupled` section of `config_overrides.yaml` and a
case's `user_nl_omega` are free to define brand-new `IOStreams` entries,
but a stream referenced in `input_files.yaml` must already exist in
`Default.yml`, since `input_files.yaml` only supplies a `Filename`
override for an existing `IOStreams` entry. `IOStreams` is not permitted
under a mesh entry in `config_overrides.yaml`'s `meshes` section, since
per-mesh IOStreams overrides aren't a supported use case; put IOStreams
that apply to every mesh under `coupled` instead. The required streams
(`HorzMeshIn`, `InitialVertCoord`, `InitialState`) are always enforced.

```{note}
Because CIME imports and executes `buildnml` in-process within the Python
interpreter running CIME commands (`case.setup`, `case.submit`, `create_test`),
`pyyaml` must be importable in that interpreter. On supported machines, source
the machine's CIME environment loader before invoking CIME commands (see the
{ref}`CIME environment table <omega-cime-env>`). For discussion on
shared Python environments across E3SM, see
[E3SM Discussion #8906](https://github.com/E3SM-Project/E3SM/discussions/8906).
```

## Validation and CI

`components/omega/cime_config/validate_config.py` runs the same validation
against every mesh and override entry (rather than just the ones needed for
one case), so it catches problems anywhere in the packaged configuration:

```bash
cd components/omega/cime_config
./validate_config.py
```

The `omega-buildnml` GitHub Actions workflow runs this script, along with
`omega_buildnml`'s pytest unit tests across supported Python versions, on
every pull request touching `cime_config/` or `Default.yml`. Run both
locally before opening a PR that changes either:

```bash
cd components/omega/cime_config
./validate_config.py
python -m pytest omega_buildnml/tests -v
```

## Adding a supported mesh

1. Confirm the mesh's grid alias is already defined for E3SM in
   `cime_config/config_grids.xml` at the repository root; `buildnml` looks
   up input files by the case's `OCN_GRID` value.
2. Add an entry for the mesh to `data/input_files.yaml`, listing the input
   file(s) that provide its `HorzMeshIn`, `InitialVertCoord`, and
   `InitialState` streams (split across multiple `inputs` entries if the
   initial condition is a separate file from the mesh).
3. If the mesh needs overrides beyond Omega's defaults and the coupled
   overrides, add a `meshes.<mesh name>` entry to `data/config_overrides.yaml`.
4. Run `./validate_config.py` to confirm the new entries are complete and
   consistent with `Default.yml`.
