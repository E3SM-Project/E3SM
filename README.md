# E3SM Diagnostics Package

[![Anaconda Version](https://anaconda.org/e3sm/e3sm_diags/badges/version.svg)](https://anaconda.org/e3sm/e3sm_diags)
[![Anaconda Downloads](https://anaconda.org/e3sm/e3sm_diags/badges/downloads.svg)](https://anaconda.org/e3sm/e3sm_diags)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1009157.svg)](https://doi.org/10.5281/zenodo.1009157)
[![Docker Version](https://images.microbadger.com/badges/version/e3sm/e3sm_diags.svg)](https://hub.docker.com/r/e3sm/e3sm_diags/)

[![CI/CD Build Workflow](https://github.com/E3SM-Project/e3sm_diags/actions/workflows/build_workflow.yml/badge.svg)](https://github.com/E3SM-Project/e3sm_diags/actions/workflows/build_workflow.yml)
[![CI/CD Release Workflow](https://github.com/E3SM-Project/e3sm_diags/actions/workflows/release_workflow.yml/badge.svg)](https://github.com/E3SM-Project/e3sm_diags/actions/workflows/release_workflow.yml)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![flake8](https://img.shields.io/badge/flake8-enabled-green)](https://github.com/PyCQA/flake8)
[![Checked with mypy](http://www.mypy-lang.org/static/mypy_badge.svg)](http://mypy-lang.org/)

## Table of Contents

1. [Documentation](#doc)
2. [Overview](#overview)
3. [Current State](#current-state)

## Documentation <a name="doc"></a>

- [Documentation Website](https://e3sm-project.github.io/e3sm_diags)
- [Sample Output, Model vs Observations](https://web.lcrc.anl.gov/public/e3sm/e3sm_diags_test_data/unit_test_complete_run/expected/all_sets/all_sets_v250_20210608_1878984/viewer/)
- Quick Start Guides:
  - [Quick Start Guide for NERSC Cori](https://e3sm-project.github.io/e3sm_diags/_build/html/master/quickguides/quick-guide-cori.html)
  - [Quick Start Guide for COMPY](https://e3sm-project.github.io/e3sm_diags/_build/html/master/quickguides/quick-guide-compy.html)
  - [Quick Start Guide for LCRC-Anvil](https://e3sm-project.github.io/e3sm_diags/_build/html/master/quickguides/quick-guide-anvil.html)
  - [Quick Start Guide for LCRC-Chrysalis](https://e3sm-project.github.io/e3sm_diags/_build/html/master/quickguides/quick-guide-chrysalis.html)
- Example run scripts:
  - [Model Climo vs Observation Climo Comparison](https://github.com/E3SM-Project/e3sm_diags/blob/master/examples/ex5-model-vs-obs)
  - [Model Climo vs Model Climo Comparison](https://github.com/E3SM-Project/e3sm_diags/blob/master/examples/ex4-model-vs-model)
  - [Model Time-series vs Model Time-series](https://github.com/E3SM-Project/e3sm_diags/blob/master/examples/ex1-model_ts-vs-model_ts)
  - [Model Time-series vs Model Time-series with CMIP data](https://github.com/E3SM-Project/e3sm_diags/blob/master/examples/ex2-model_ts-vs-model_ts-cmip)
  - [Model Time-series vs Observation Time-series with CMIP data](https://github.com/E3SM-Project/e3sm_diags/blob/master/examples/ex3-model_ts-vs-obs_ts-cmip)
  - [Observation vs Observation Comparison](https://github.com/E3SM-Project/e3sm_diags/tree/master/examples/ex7-obs-vs-obs)

## Overview<a name="overview"></a>

This diagnostics package is constructed for supporting the diagnostics task of DOE's [Energy Exascale Earth System Model (E3SM) project](https://climatemodeling.science.energy.gov/projects/accelerated-climate-modeling-energy). The goal of this work is to develop a comprehensive diagnostics package that:

- fully integrates the functionality of NCAR's AMWG diagnostics package.
- utilizes most updated observational datasets, including remote sensing, reanalysis and in-situ datasets.
- interfaces with diagnostics developed from different E3SM focus groups: atmosphere group, coupled simulation group, land group.
- interacts effectively with the PCMDI's metrics package and the ARM diagnostics package through a unifying framework: [Community Diagnostics Package (CDP)](https://github.com/CDAT/cdp).
- is flexible for user specified diagnostics and being configured for use by other climate models.

## Current State <a name="current-state"></a>

Algorithm and visualization codes for the AMWG Set 5, 7, 4, 3, 13, 1, 14 diagnostics, namely lat-lon contour plots (Figure 1), polar contour plots (Figure 2), zonal mean 2d plots (Figure 3), zonal mean line plots (Figure 4), 2d joint histogram for COSP cloud simulator output (Figure 5), tables (Figure 6) and Taylor Diagrams (Figure 7) summarizing metrics, for climatology seasonal means, are implemented.

The package features built-in user diagnostics, by specifying user desired diagnostics regions and pressure levels for variables with the vertical dimension.

In addition to default model versus observation comparison, the package also provide support for model versus model and obs versus obs comparisons.

<img src="misc/example_fig1.png" alt="Figure1" style="width: 280px;"/>
<h5 align="center">Figure 1: An example of the lat-lon contour plots for air temperature at 850 mb with tropical ocean region considered</h5>

<img src="misc/example_fig2.png" alt="Figure2" style="width: 280px;"/>
<h5 align="center">Figure 2: An example of the polar contour plots for precipitation rate</h5>

<img src="misc/example_fig3.png" alt="Figure3" style="width: 280px;"/>
<h5 align="center">Figure 3: An example of the pressure-lat contour plots for air temperature </h5>

<img src="misc/example_fig4.png" alt="Figure4" style="width: 280px;"/>
<h5 align="center">Figure 4: An example of the zonal mean surface air temperature line plot </h5>

<img src="misc/example_fig5.png" alt="Figure5" style="width: 280px;"/>
<h5 align="center">Figure 5: An example of 2d joint histogram plot using COSP simulator output</h5>

<img src="misc/example_fig6.png" alt="Figure6" style="width: 280px;"/>
<h5 align="center">Figure 6: An example of table summarizing metrics calculated based on lat-lon contour plots diagnostics</h5>

<img src="misc/example_fig7.png" alt="Figure7" style="width: 280px;"/>
<h5 align="center">Figure 7: An example of Taylor diagram summarizing metrics calculated based on lat-lon contour plots diagnostics of several key variables</h5>

## License

Copyright (c) 2018-2021, Energy Exascale Earth System Model Project
All rights reserved

SPDX-License-Identifier: (BSD-3-Clause)

See [LICENSE](./LICENSE) for details

Unlimited Open Source - BSD 3-clause Distribution
`LLNL-CODE-819717`
