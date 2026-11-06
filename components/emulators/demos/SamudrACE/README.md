# SamudrACE demo

This demo runs the atmosphere (ACE) and ocean (Samudra) from the [E3SMv3 SamudrACE emulators](https://huggingface.co/mahf708/e3smv3-lr-hist-samudrace) through the emulator package, at their different rates. A small C++ program plays the part of the coupled system: it owns the external buffers and supplies the forcing and initial state.

Each coupled step runs 20 atmosphere steps of 6 hours, then 1 ocean step of 5 days, about 3 minutes in total on a CPU.

| File | Purpose |
| --- | --- |
| `fme_component.py` | Runs one `fme` step for each `infer` call of the `Python` backend. |
| `prepare_demo.py` | Writes the configs, inputs, and reference results from `fme`'s own coupled rollout. |
| `samudrace_demo.cpp` | Steps both emulators through the `Harness`. |
| `compare_demo.py` | Compares the demo's outputs with the reference results. |
| `pyproject.toml`, `uv.lock` | Pin the demo's Python environment, including the `fme` commit. |

## Running the demo

From `components/emulators`:

```bash
./run_demos
```

The script needs `uv`, CMake 3.28+, MPI, and C, C++20, and Fortran compilers. On the first run it does the following, and on later runs it reuses whatever already exists:

1. Creates a Python environment from `pyproject.toml` and `uv.lock`, using `uv`'s own Python.
2. Downloads the checkpoint, initial conditions, and 1985 atmosphere and ocean forcing from Hugging Face, if any are missing.
3. Writes the model's field list, using `scripts/config/extract_fme_fields.py`.
4. Builds the driver in `build-demo`.
5. Prepares the inputs and reference results, runs the demo, and compares them.

This takes about 9 GB of disk, in `demo_data/SamudrACE` by default. Run `./run_demos --help` for its options, such as `--n-coupled-steps`, which defaults to 2. The run takes about twice as long as the model time, because `prepare_demo.py` also runs `fme` to produce the reference results.

The driver is built against the environment's Python (`.venv/bin/python`), because it has to embed the same Python the model runs under. LibTorch stays off, because `fme` brings its own torch, and two LibTorch versions in one process fail with symbol errors.

## Expected results

With 2 coupled steps:

| Steps | Result, as the largest difference relative to the field's standard deviation |
| --- | --- |
| Atmosphere steps 0 to 19 | 50/50 fields bit-for-bit identical to `fme` |
| Ocean step 0 | 1.87e-04, with an RMS of about 5e-06 at most |
| Atmosphere step 20, the first after the ocean's update | 4.66e-05 |
| Atmosphere step 39 | 2.46e-03, with an RMS of about 3e-05 at most |
| Ocean step 1 | 2.32e-04 |

The ocean can't match exactly because of rounding. `fme` averages the fluxes in `float32`, and also averages `LANDFRAC` over the window, which changes that constant field by up to 1.8e-07. This package averages in `double`, and passes `LANDFRAC` unchanged. Using `fme`'s rounding for both makes the first ocean step identical in all 80 fields. Those small differences then grow, as they do in any chaotic model.

`compare_demo.py` exits with an error if the first window isn't identical, or if the first ocean step or atmosphere step 20 differs by more than 1e-03.

## How the coupling works

The plan has these kinds of routes:

| Route | Fields | Behavior |
| --- | --- | --- |
| Atmosphere to ocean | 9 fluxes | Averaged over the 20 atmosphere steps, because the ocean runs less often |
| Ocean to atmosphere | `sst`, `ocean_sea_ice_fraction` | Held between ocean steps. They start from the ocean's initial state, because the atmosphere runs first |
| Through `cpl7` | Prognostic fields | `<name>_out` and `<name>_in` routes. The ocean's use the state's own buffer, so nothing is copied between its steps |
| From `cpl7` | Forcing | Latest value |

`fme_component.py` holds the steps that `fme`'s coupled stepper applies between components:

- **Atmosphere:** derives `OCNFRAC` and `ICEFRAC` from `sst`, `ocean_sea_ice_fraction`, and `LANDFRAC`. It also overwrites input `TS` from `sst` at the start of each 20-step window, and only then.
- **Ocean:** scales the averaged fluxes using the open-water fraction, `atmosphere_flux_fraction`, and `momentum_flux_fraction`.

## How inputs are handled

- **`SOLIN`** must hold its value at the end of the step, because it's in the atmosphere's `next_step_forcing_names`. The sun moves a long way in 6 hours, so this matters. It's read from a file so results can be compared with `fme`.
- **`global_mean_co2`** is a uniform field that `samudrace_demo.cpp` computes each step. Its formula is only valid for 1985.
- **`seed`** is passed as an option. As in `fme`, the atmosphere uses the seed itself, and the ocean uses `alternate_seed(seed)`.
- **Corrector and random state** stay in the wrapper. They aren't saved, so a run can't be restarted yet.

`prepare_demo.py` and `fme_component.py` use some private `fme` attributes, such as `_step_obj` and `_get_atmosphere_forcings`, so they may need changes for other `fme` versions.
