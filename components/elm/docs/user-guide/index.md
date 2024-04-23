# ELM User's Guide

This User's Guide describes how to set up and run ELM.

## Table of Conents

1. [Steps to build and run ELM](#steps-to-build-and-run-elm)
    1. [Scientifically supported compsets](#scientifically-supported-compsets)
    2. [Supported grids](#supported-grid)
    3. [Model spin-up for pre-industrial condition](model-spin-up-for-pre-industrial-condition)
2. [Customizing runs](customizing-runs)
    1. [Changing monthly output file](changing-monthly-output-file)
    2. [Saving additional output files](saving-additional-output-files)
3. [Running with FATES (optional)](#running-with-fates)

## Steps to build and run ELM

A step-by-step instruction on how to run fully coupled E3SM can be found [here](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/2309226536). Here we describe running ELM driven by atmospheric forcings provided via the data atmosphere (DATM) model for configurations that are used in the E3SM water cycle v3 campaign.

The water cycle campaigns of E3SM v1 and v2 used ELM's satellite phenology mode (SP-mode) in which a prescribed leaf area index is used in ELM. However, the E3SM v3 water cycle campaign uses an interactive phenology by including an active biogeochemistry (BGC) cycle in ELM. Additionally, a parameterization of sub-grid topographical effects on solar radiation is included within ELM.

### Scientifically supported compsets

The land-only compsets are referred to as "I"-compset and are supported for the following time periods: pre-industrial (1850) and historical transient (20TR). Additionally, multiple atmospheric forcing datasets can be used to drive the ELM simulations. The supported compsets are:

1. `I1850GSWCNPRDCTCBCTOP`: Climatological pre-industrial using GSWP atmospheric forcings
2. `I20TRGSWCNPRDCTCBCTOP`: Historical ELM simulation using GSWP atmospheric forcings with time varying greenhouse gas forcing and land use, land cover dataset (year 1850-2014).

### Supported grid

The `r05_r05` is the supported grid resolution for performing offline ELM simulation.

### Model spin-up for pre-industrial condition

***Add notes on how to spin-up the model.***

## Customizing runs

Few useful changes to `user_nl_elm`

### Changing monthly output file

ELM by default outputs monthly history file in `*elm.h0.**.nc` files
that include many variables (>200). At times, many of the default output
variables may not be of interest, thus one could remove all default variables
(via `hist_empty_htapes`) and only include select variables (via `hist_fincl1`)
to the monthly history files by

```fortran
   &elm_inparm
     hist_empty_htapes = .true.
     hist_fincl1 = 'TG', 'TV', 'FSA'
   /
```

#### Saving additional output files

ELM can output additional history files (such as `*elm.h1.*.nc`, `*elm.h2.*.nc`)
that have different temporal averaging (e.g. daily, hourly, every model timestep) via
`hist_nhtfrq` where

- `-24` corresponds to daily average
- `-1` corresponds to hourly average
- `0` corresponds to monthly average
- `1` corresponds to each model time step

The number of time slices in these additional files can be controlled
vai `hist_mfilt`.

```fortran
   &elm_inparm
     hist_fincl2 = 'TG'
     hist_fincl3 = 'TV'
     hist_fincl4 = 'TG', 'TV', 'FSA'
     hist_nhtfrq = 0, -24, -1, 1
     hist_mfilt  = 12, 30, 24, 48
   /
```

Using the above-mentioned settings:

- Each `*.elm.h1.*.nc` will include 30 daily average values of `TG`
- Each `*.elm.h2.*.nc` will include 24 hourly average values of `TV`
- Each `*.elm.h3.*.nc` will include 48 values of `TG`, `TV`, and `FSA` at
  each model time step, which is typically is 30 min.
  
## Running with FATES

[FATES](fates.md) can be run in various modes with ELM through the use of namelist settings.  The [FATES User's Guide section on namelist options](https://fates-users-guide.readthedocs.io/en/latest/user/Namelist-Options-and-Run-Time-Modes.html) provides guidance on enabling these different FATES run modes.
