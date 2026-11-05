# SamudrACE demo

This demo runs the atmosphere (ACE) from the [E3SMv3 SamudrACE emulators](https://huggingface.co/mahf708/e3smv3-lr-hist-samudrace) through the emulator package. A small C++ program plays the part of the coupled system: it owns the external buffers, supplies every input, and takes every output.

The ocean (Samudra) isn't part of the demo yet.

| File | Purpose |
| --- | --- |
| `fme_component.py` | Runs one `fme` step for each `infer` call of the `Python` backend. |
| `prepare_demo.py` | Writes the configs, inputs, and reference results from `fme`'s own rollout. |
| `samudrace_demo.cpp` | Steps the emulator through the `Harness`. |
| `compare_demo.py` | Compares the demo's outputs with the reference results. |
| `pyproject.toml`, `uv.lock` | Pin the demo's Python environment, including the `fme` commit. |

## Running the demo

From `components/emulators`:

```bash
./run_demos
```

The script needs `uv`, CMake 3.28+, MPI, and C, C++20, and Fortran compilers. On the first run it does the following, and on later runs it reuses whatever already exists:

1. Creates a Python environment from `pyproject.toml` and `uv.lock`, using `uv`'s own Python.
2. Downloads the checkpoint, initial conditions, and 1985 forcing from Hugging Face, if any are missing.
3. Writes the model's field list, using `scripts/config/extract_fme_fields.py`.
4. Builds the driver in `build-demo`.
5. Prepares the inputs and reference results, runs the demo, and compares them.

This takes about 9 GB of disk, in `demo_data/SamudrACE` by default. Run `./run_demos --help` for its options, such as `--n-steps`, which can be up to 20.

The driver is built against the environment's Python (`.venv/bin/python`), because it has to embed the same Python the model runs under. LibTorch stays off, because `fme` brings its own torch, and two LibTorch versions in one process fail with symbol errors.

## Expected results

```
step_0: 50/50 fields identical; largest difference 0.00e+00 of a standard deviation (surface_precipitation_rate)
step_1: 0/50 fields identical; largest difference 2.36e-04 of a standard deviation (STW_1)
```

The first step is bit-for-bit identical to `fme`, and `compare_demo.py` exits with an error if it isn't within a small tolerance. Later steps differ slightly because the demo holds `global_mean_co2` at its first value, while `fme` uses a new value each step.

## How inputs are handled

- **Prognostic fields** go through `cpl7` as `<name>_out` and `<name>_in` routes. See the prognostic fields section of the package README.
- **`SOLIN`** must hold its value at the end of the step, because it's in the model's `next_step_forcing_names`. The sun moves a long way in 6 hours, so this matters.
- **`sst`** is an extra input, which the model prescribes as `TS` at the end of the step. `fme` holds it constant across the atmosphere's 20 steps, as it does `OCNFRAC` and `ICEFRAC`.
- **`global_mean_co2`** is a scalar passed in the `scalar_global_mean_co2` option, and is applied at every grid point. An option can't change during a run, so longer runs need another approach.
- **`seed`** is passed as an option, and seeds the model's random state.
- **Corrector and random state** stay in the wrapper. They aren't saved, so a run can't be restarted yet.

`prepare_demo.py` and `fme_component.py` use some private `fme` attributes, such as `_step_obj` and `_get_atmosphere_forcings`, so they may need changes for other `fme` versions.
