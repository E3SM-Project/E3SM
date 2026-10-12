# P3's warm-rain emulator

P3 can take the rates of its four warm-rain collision processes (autoconversion,
accretion, droplet self-collection, rain self-collection) from a model: a python
module, run through the python inference backend of `components/emulators`.

## How it works

```
[emulator: state -> rates, masks] -> p3_main: part1 -> part2 (merge) -> sedimentation -> part3
```

The emulator runs before `p3_main`, outside of any kernel, so it works with
monolithic and with small kernels. Its inputs are the state P3 computes its own
warm-rain rates from: the state at the start of the step, with what
`p3_main_part1` does to it first (dry density, clipping of trace cloud and
rain, activation or prescription of the droplet number, melting of trace ice
into rain; see `run_warm_rain_emulator`, which must be kept in sync with
part1). Its outputs are inputs of `p3_main`, like the nucleation rates P3 gets
from other parameterizations: part2 computes its rates as usual, maps them to
cell averages, then merges the emulated ones (`warm_rain_emulator_merge`):

- where the emulator says it applies (`use_cloud`, `use_rain`), its rates
  replace P3's;
- elsewhere P3's rates stay, with autoconversion multiplied by
  `warm_rain_emulator_autoconv_fallback_scale`;
- emulated rates are limited so that the conservation checks that follow hold
  (sources are not negative, self-collection does not remove more number than
  there is in a step), then go through those checks as P3's do.

With `warm_rain_emulator_replaces_p3`, P3 does not compute its own warm-rain
rates and uses the emulator's in every cell, whatever its masks say. This is
for models meant to apply everywhere; it saves about 3 % of `p3_main`.

Without the emulator, nothing changes: results are bit-for-bit those of P3.

## The model

A python module with `create_emulator(config)`, which returns an object with
`infer(inputs, outputs)`. `inputs` and `outputs` are dicts of `(ncol, nlev)`
numpy arrays of doubles, by name. They view memory of P3 (nothing is copied on
a CPU; on a GPU the arrays are copied to the host and back), the inputs are
read-only, and the outputs are written in place (`outputs[name][...] = ...`).

| | names |
|---|---|
| inputs | `qc`, `nc`, `qr`, `nr` (kg/kg, #/kg; cell averages, dry), `rho` (kg/m3, dry) |
| outputs | `qc2qr_autoconv_tend`, `qc2qr_accret_tend`, `ncautr`, `nc2nr_autoconv_tend`, `nc_accret_tend`, `nc_selfcollect_tend` (≤ 0), `nr_selfcollect_tend`: cell averages, with P3's units and signs; `use_cloud`, `use_rain`: 1 where to use the first six, and the last one |

`config` is a dict of strings, with `model_path` (the parameter
`warm_rain_emulator_model_path`) and `num_threads` (the threads of the MPI
rank: the model must not use more).

`sdm_warm_rain.py` is the SDM-trained emulator (ERF super-droplet LES of
shallow cumulus): a small torch network, with its normalization, its training
envelope (the masks) and its number-rate closures, all read from the model file
`sdm_warm_rain.pt`. It runs the network only on the cells where P3 uses its
rates (inside the envelope), and returns zeros elsewhere: most cells are clear
or outside the envelope, so this makes a call several times cheaper.
`SdmWarmRain(file, compact=False)` computes every cell, as in the model's
reference vectors.

Its rates are means over 100 s, from the state at the start of that period: it
is meant for P3 time steps close to that. Its authors suggest a fallback scale
of 0.36.

## Use

Build with

```
-DEAMXX_ENABLE_INFERENCE=ON -DPython3_EXECUTABLE=/path/to/python
```

where the python has its development files, numpy, and what the module imports
(torch, for `sdm_warm_rain`). It can be the python of a virtual environment.
Then set, in `p3`:

```yaml
use_warm_rain_emulator: true
warm_rain_emulator_python_module: sdm_warm_rain
warm_rain_emulator_python_path: /path/to/eamxx/src/physics/p3/warm_rain_emulator
warm_rain_emulator_model_path: /path/to/eamxx/src/physics/p3/warm_rain_emulator/sdm_warm_rain.pt
warm_rain_emulator_autoconv_fallback_scale: 0.36
```

In a CIME case: `./xmlchange --append SCREAM_CMAKE_OPTIONS="..."` and
`./atmchange p3::use_warm_rain_emulator=true`, etc. (see the
`eamxx-p3_warm_rain_emulator` test modifier). The other parameters
(`warm_rain_emulator_cloud_self_collection`, `warm_rain_emulator_replaces_p3`)
are in `namelist_defaults_eamxx.xml`.

The executable embeds the python it was configured with, and tells it where its
standard library and its packages are (`PYTHONHOME`, if set in the environment,
is used instead of the former).

Not run yet: single precision, and GPUs other than CUDA. Nothing checks the P3
time step against the model's period.

## Reproducibility

A run with the emulator is reproducible, and restarts exactly. It is not
guaranteed to be bit-for-bit independent of the number of MPI ranks: torch does
not always compute a grid point in exactly the same way in batches of different
sizes, so the rates can differ in the last bits of the (single precision)
network when the columns are distributed differently.

## Tests

- `physics/p3/tests/p3_warm_rain_emulator_tests.cpp`: the merge.
- `tests/single-process/p3_warm_rain_emulator`: P3 with an emulator that never
  applies is bit-for-bit P3; with one of zero rates that applies everywhere,
  there is no autoconversion and no accretion; with one that converts the
  cloud water it is given to rain, there is none left (these are numpy
  modules, in `synthetic_emulators.py`); the SDM-trained emulator on different
  numbers of ranks.
- CIME: the `eamxx-p3_warm_rain_emulator` test modifier.
