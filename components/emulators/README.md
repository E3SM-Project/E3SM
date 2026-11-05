# E3SM Emulator Package

This package provides the infrastructure needed to couple an emulator with E3SM.
The emulators are fully configurable via `yaml` configuration files that can be generated
using the python script in `./scripts/create_emulator_config_files.py`.
See example in `./scripts/create_coupling_test_files.py` that is used to generate config files for tests.

The required files for running an emulator are:

1. `emulator_config.yaml`: 
  - describes one emulator: its model, inference backend, and input and output fields. Each emulator has its own file.
2. `coupling_plan.yaml`
  - describes which component exports each coupled field and which components import it. A coupled configuration has one plan.

These names are a convention. The code takes the paths as arguments, so the files can be named anything.

> **Status:** Emulators can be coupled to each other and to external buffers through the C++ API, as exercised by the tests.
> The connection to `CPL7` is still in progress; see [Current limitations](#current-limitations).

## Package Structure

- `common/`: implementation of the emulator class and its utilities.
- `inference/`: implementation of available inference backends (`LIBTORCH` and `PYTHON`, plus a `STUB` for testing) for invoking an emulator.
- `coupler_api/`: API layer for emulators to exchange data with other emulators.
- `interface/`: implementation of the `Harness` to allow emulators to exchange with `CPL7`.
- `scripts/`: python helpers for generating configuration files.
- `tests/`: tests for the emulator and inference layers, along with their fixtures. `coupler_api/` and `interface/` have their own `tests/` directories.

## Configuration

### Generating the files

The scripts require Python 3.10+ and PyYAML. Generating the files through them also catches mistakes such as
duplicate field names, unknown dimensions, or a coupled field that a component doesn't declare.

This example couples a land emulator to `CPL7`:

```python
from create_emulator_config_files import (
    CoupledField,
    CouplingPlan,
    EmulatorConfig,
    Field,
    Model,
)

temperature = Field(name="temperature", dimensions=["column"])
soil_moisture = Field(name="soil_moisture", dimensions=["column", "soil_levels"])
surface_flux = Field(name="surface_flux", dimensions=["column"])

land = EmulatorConfig(
    name="lnd-emulator",
    model=Model(
        backend="libtorch", path="land_emulator.pt", options={"dtype": "float64"}
    ),
    dimensions=["column", "soil_levels"],
    inputs=[temperature, soil_moisture],
    outputs=[surface_flux],
)

# Only used to name cpl7 in the plan and check its fields. No config file is written for it.
cpl7 = EmulatorConfig(
    name="cpl7",
    model=Model(backend="stub"),
    dimensions=[],
    inputs=[surface_flux],
    outputs=[temperature, soil_moisture],
)

land.write("emulator_config.yaml")

CouplingPlan(
    fields=[
        CoupledField(field=temperature, kind="state", sources=[cpl7], destinations=[land]),
        CoupledField(field=soil_moisture, kind="state", sources=[cpl7], destinations=[land]),
        CoupledField(field=surface_flux, kind="flux", sources=[land], destinations=[cpl7]),
    ]
).write("coupling_plan.yaml")
```

To regenerate the test configurations, run `python scripts/create_coupling_test_files.py`. The files are written to the current directory.

A field that an emulator both reads and writes needs a `PrognosticField`, which creates the plan entries described in [Prognostic fields](#prognostic-fields). `prognostic_fields` finds those fields in a config. The `cpl7` config must list each one in both its inputs and outputs:

```python
from create_emulator_config_files import PrognosticField, prognostic_fields

plan = CouplingPlan(
    fields=[
        *(
            PrognosticField(field=field, component=atmosphere, external=cpl7)
            for field in prognostic_fields(atmosphere)
        ),
        CoupledField(field=sst, kind="state", sources=[cpl7], destinations=[atmosphere]),
    ]
)
```

### Getting the fields from an `fme` checkpoint

For a coupled `fme` model such as SamudrACE, `scripts/extract_fme_fields.py` writes the model's inputs and outputs, in order, to a YAML field list. It also records each component's timestep and forcing fields, and the fields the components exchange. It must run in an environment with `fme` installed:

```bash
python scripts/extract_fme_fields.py checkpoints/C1A0.tar -o fields.yaml
```

`emulator_config_from_field_list` then creates a config from that list, without needing `fme` or torch:

```python
from create_emulator_config_files import Model, emulator_config_from_field_list

atmosphere = emulator_config_from_field_list(
    "fields.yaml",
    "atmosphere",
    name="atm-emulator",
    model=Model(backend="python", options={"python_module": "my_model"}),
    dimensions=["lat", "lon"],
)
atmosphere.write("atm_emulator_config.yaml")
```

Use the `field_dimensions` argument for any field whose dimensions differ from `dimensions`.

### `emulator_config.yaml`

```yaml
emulator:
  name: lnd-emulator
  type: component
  model:
    path: land_emulator.pt
    input_layout: FeatureMajor
    output_layout: FeatureMajor
    backend: libtorch
    options:
      dtype: float64
sample_space:
  dimensions:
  - column
  - soil_levels
input_order:
- temperature
- soil_moisture
output_order:
- surface_flux
inputs:
  temperature:
    name: temperature
    dimensions:
    - column
  soil_moisture:
    name: soil_moisture
    dimensions:
    - column
    - soil_levels
outputs:
  surface_flux:
    name: surface_flux
    dimensions:
    - column
```

| Key | Description |
| --- | --- |
| `emulator.name` | Name of the emulator, used to refer to it in the coupling plan. |
| `emulator.type` | `process` for a process emulator. Any other value is treated as a component. |
| `model.path` | Path to the model file. Relative paths are resolved against the directory containing the YAML file. |
| `model.input_layout`, `model.output_layout` | `FeatureMajor` or `SampleMajor`. Required, but not yet used by the backends. |
| `model.backend` | `libtorch`, `python`, or `stub` (case-insensitive). |
| `model.options` | Backend-specific settings. All values are strings. See [Inference backends](#inference-backends). |
| `sample_space.dimensions` | Names of the dimensions used by the emulator's fields. |
| `input_order`, `output_order` | Order of the fields in the emulator's buffers. `LIBTORCH` also uses this order for the model's arguments and results. |
| `inputs`, `outputs` | Field name and dimensions, keyed by field name. Each dimension must be defined in the emulator's `GridDesc`. |

Field data is stored as `double` in row-major order, so the last dimension varies fastest.
For example, `soil_moisture` above is indexed as `[column * nlevels + level]`.

### `coupling_plan.yaml`

```yaml
fields:
  temperature:
    kind: state
    name: temperature
    merge_type: direct
    sources:
    - cpl7
    destinations:
    - lnd-emulator
  surface_flux:
    kind: flux
    name: surface_flux
    merge_type: direct
    sources:
    - lnd-emulator
    destinations:
    - cpl7
```

Each key under `fields` is a route ID, which must be unique in the plan. `RouteBufferMap` also uses it.

| Key | Description |
| --- | --- |
| `name` | Name of the field the route carries. Several routes can carry the same field. |
| `kind` | Type of coupling, such as `state` or `flux`. Required, but doesn't affect the coupling yet. |
| `merge_type` | How sources are merged. Only `direct` is currently supported, which requires exactly one source. |
| `sources` | Component that exports the field. |
| `destinations` | Components that import the field. |

The rules below apply when the coupler is initialized:

- Every component listed for a route must have registered its `name` field: as an export if it's a source, or as an import if it's a destination.
- A route needs at least one destination, and a component can only be listed once in a route's `sources` or `destinations`.
- Two routes can't deliver the same field to the same component.
- Every component added to the driver must appear in the plan.
- The name of the `CPL7` component is always `cpl7`.

## Inference backends

Backends other than `STUB` are enabled at build time. See [Building and testing](#building-and-testing).

| Backend | Model | Options |
| --- | --- | --- |
| `libtorch` | TorchScript module saved with `torch.jit.save` | `device`, `dtype`, `num_threads`, `jit_optimize` |
| `python` | Python module run in an embedded interpreter | `python_module` (required), `python_factory`, `python_path` |
| `stub` | None. Outputs are left unchanged, which is useful for testing the data pipeline. | None |

### LibTorch

Inputs are passed to `forward()` positionally, in `input_order`. The module returns a tensor, or a tuple or list of tensors,
in `output_order`. Output shapes must match the configured dimensions exactly.

| Option | Default | Description |
| --- | --- | --- |
| `device` | `cpu` | `cpu`, `cuda`, or `cuda:N`. An unavailable CUDA device is an error rather than a fallback to the CPU. |
| `dtype` | `float32` | Model precision, `float32` or `float64`. Tensors are converted on the way in and out. |
| `num_threads` | `0` | Intra-op threads. `0` keeps the torch default. |
| `jit_optimize` | `false` | Lets TorchScript re-optimize the graph after the first calls. This is faster, but it changes results between early and later calls, which breaks exact restarts. |

### Python

The backend imports `python_module` and calls its factory with a dict holding the configuration (`model_path`, `input_channels`,
`output_channels`, `verbose`, and every option as a string). The factory returns an object with an `infer` method and, optionally, a `finalize` method:

```python
def create_emulator(config):
    return MyModel(config["model_path"])


class MyModel:
    def infer(self, inputs, outputs):
        outputs["dT"][:] = self.net(inputs["T"])
```

`inputs` and `outputs` are dicts of `float64` numpy arrays keyed by field name. They view the emulator's memory directly:
inputs are read-only, and the model must write its results in place.

| Option | Default | Description |
| --- | --- | --- |
| `python_module` | Required | Module to import. |
| `python_factory` | `create_emulator` | Factory function in the module. |
| `python_path` | Empty | Colon-separated directories to prepend to `sys.path`. |

numpy is needed at run time, but not to build. If the process already runs an interpreter, it is shared.

## Coupling with the `Harness`

The `Harness` owns a `CouplerDriver`, the emulators, and a `CPL7Component` that represents the external buffers.
Each buffer is bound to a field by name.

```cpp
using namespace e3sm::emulator;

Harness harness;

// x2l and l2x spans are owned by the caller.
harness.transfer_bindings(ExternalBindings{
    .to_external = {{"surface_flux", flux}},
    .from_external = {{"temperature", temperature}, {"soil_moisture", moisture}},
});

harness.add_process(ProcessDescriptor{
    .run_opts = ProcessRunOpts{.dt = 1., .coupling_interval = 1},
    .grid = grid,
    .emulator_config = read_emulator_config("emulator_config.yaml"),
});

harness.initialize("coupling_plan.yaml");

// Each coupling step
harness.export_from_external("cpl7");  // external buffers -> routes
harness.run_process("lnd-emulator");   // import, run inference, export
harness.import_to_external("cpl7");    // routes -> external buffers
```

Call these in order: `transfer_bindings` and `add_process` must both come before `initialize`, because components can't be added afterwards.

`initialize` also accepts a `RouteBufferMap`, which uses caller-owned storage for the named routes instead of allocating it.
Each buffer must be the same size as the route's source field. See `interface/tests/test_harness.cpp` for a complete example.

Destroying the `Harness` finalizes its emulators, including calling the Python model's `finalize` method. A failure there is logged to stderr, because destructors can't throw.

To couple emulators without `CPL7`, or to add another kind of component, use `CouplerDriver` directly.
Any class that satisfies the `CoupledProcess` concept in `coupler_api/include/process.hpp` can be added with `add_component`.
`coupler_api/tests/test_coupler_driver.cpp` shows this with fake components.

### Prognostic fields

A prognostic field is both an input and an output of an emulator, such as an atmosphere's `T_0`. Nothing copies an output back to its input automatically, so the field leaves through `cpl7` and returns on the next step, using two routes:

```yaml
fields:
  T_0_out:
    kind: state
    name: T_0
    merge_type: direct
    sources:
    - atm-emulator
    destinations:
    - cpl7
  T_0_in:
    kind: state
    name: T_0
    merge_type: direct
    sources:
    - cpl7
    destinations:
    - atm-emulator
```

The caller owns the field's storage. Bind the same span in both directions, and fill it with the initial state before the first step:

```cpp
std::vector<double> t_0 = initial_t_0();

harness.transfer_bindings(ExternalBindings{
    .to_external = {{"T_0", t_0}},
    .from_external = {{"T_0", t_0}},
});
```

The usual coupling step then advances the field: `export_from_external` passes the span's values to the emulator's input, `run_process` writes the new output, and `import_to_external` copies that output back to the span.

These fields don't need to exist in the real `CPL7` coupler, because the span is only storage on the caller's side.

A name that is both an input and an output is ambiguous in `Emulator::get_field(name)`, which throws. Use `get_field(name, ModelFieldRole::Input)` or `ModelFieldRole::Output`. `use_moab` doesn't support these fields yet, because a MOAB tag name can only be defined once.

## Building and testing

The package builds standalone unless it is configured from E3SM, which defines `E3SM_SOURCE_DIR`. A standalone build requires:

- CMake 3.28+
- C, C++20, and Fortran compilers
- MPI
- Network access on the first configure, to fetch EKAT and MOAB (and LibTorch if needed) into `downloads/`

### Using the test script

The `test` script configures a Debug build in `build/` with both the LibTorch and Python backends enabled, builds it, and runs the tests:

```bash
./test            # run all tests
./test moab       # run tests with names matching a regular expression
./test --help     # list the options, such as --build-only, --verbose, and -D<var>=<val>
```

### Using CMake directly

```bash
cmake -S . -B build -DEMULATOR_ENABLE_LIBTORCH=ON -DEMULATOR_ENABLE_PYTHON=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

| Option | Default | Description |
| --- | --- | --- |
| `EMULATOR_ENABLE_TESTS` | `ON` | Build the tests. |
| `EMULATOR_ENABLE_LIBTORCH` | `OFF` | Build the LibTorch backend. |
| `EMULATOR_ENABLE_PYTHON` | `OFF` | Build the Python backend. Requires the Python development files. |
| `EMULATOR_TORCH_PYTHON` | `python3` | Interpreter used to locate LibTorch and write the test fixtures. |
| `EMULATOR_LIBTORCH_VERSION` | `2.5.1` | LibTorch version to download. |
| `EMULATOR_LIBTORCH_VARIANT` | `auto` | LibTorch build to download: `auto`, `cpu`, `cu118`, `cu121`, `cu124`, or `rocm6.2`. |
| `EMULATOR_LIBTORCH_URL` | Empty | LibTorch archive to download. Overrides the variant. |

CMake looks for LibTorch in `Torch_DIR` or `CMAKE_PREFIX_PATH`, then in the torch package of `EMULATOR_TORCH_PYTHON`, and otherwise downloads it.
When using both backends, build against the torch that the Python backend will import, because loading two LibTorch versions in one process fails with symbol errors.

### Tests

| Tests | Coverage | Requires |
| --- | --- | --- |
| `inference_config_tests`, `inference_tensor_tests`, `inference_stub_backend_tests`, `fpe_guard_tests` | Inference layer | |
| `inference_libtorch_backend_tests` | LibTorch backend | `EMULATOR_ENABLE_LIBTORCH`, torch in `EMULATOR_TORCH_PYTHON` |
| `inference_python_backend_tests` | Python backend | `EMULATOR_ENABLE_PYTHON`, numpy |
| `test_emulator_config` | Parsing `emulator_config.yaml` | |
| `emulator_coupling` | Coupled emulators on each enabled backend | At least one of the backend tests above |
| `test_coupler_driver` | `CouplerDriver` with fake components | |
| `moab_adapter_1`, `moab_adapter_2`, `moab_transport_6` | MOAB adapter and transport, on 1, 2, and 6 MPI ranks | MPI |
| `test_harness` | `Harness` with external buffers | `EMULATOR_ENABLE_LIBTORCH` |

A backend's tests are skipped with a CMake warning if its Python dependency can't be imported.

## Current limitations

- The `CPL7` Fortran interface is incomplete. `interface/fortran/emu_comp_mct.F90.in` is configured but not compiled, it hardcodes the `atm` component and `ace2` emulator name, and the `Harness` has no Fortran bindings yet.
- Emulators can't be selected as components in an E3SM case yet. `cime_config/config_files.xml` points to `components/emulators/emulatoratm`, which doesn't exist, and this package doesn't define the `emulatoratm` and `emulatorocn` targets that `components/cmake/build_emulator_comps.cmake` looks for.
- Only `direct` merges are supported. `scale` is accepted when the plan is read, but initialization then fails.
- `input_layout`, `output_layout`, and `kind` are required in the configuration files but don't change behavior yet.
