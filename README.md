# E3SM Diagnostics Package

[![Install](https://anaconda.org/e3sm/acme_diags/badges/installer/conda.svg)](https://anaconda.org/e3sm/acme_diags)
[![Downloads](https://anaconda.org/e3sm/acme_diags/badges/downloads.svg)](https://anaconda.org/e3sm/acme_diags)
[![Version](https://anaconda.org/e3sm/acme_diags/badges/version.svg)](https://anaconda.org/e3sm/acme_diags)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1009157.svg)](https://doi.org/10.5281/zenodo.1009157)


## Table of Contents
1. [Documentation](#doc)
2. [Overview](#overview)
3. [Current State](#current-state)

## Documentation <a name="doc"></a>
* [Documentation Website](https://acme-climate.github.io/acme_diags)
* Examples:
  * [Quick Start Guide for AIMS4/ACME1](https://acme-climate.github.io/acme_diags/docs/html/quick-guide-aims4.html)
  * [Model vs Model Comparison](https://github.com/ACME-Climate/acme_diags/blob/master/examples/model-vs-model/model-vs-model.ipynb)
  * [Observation vs Observation Comparison](https://github.com/ACME-Climate/acme_diags/blob/master/examples/obs-vs-obs/obs-vs-obs.ipynb)
  * [Model vs Observation: comparing temperature at 200mb and 800mb with different variable names in the model and obs](https://github.com/ACME-Climate/acme_diags/blob/master/examples/model-vs-obs/model-vs-obs.ipynb)

## Overview<a name="overview"></a>
This diagnostics package is constructed for supporting the diagnostics task of DOE's [Energy Exascale Earth System Model (E3SM) project](https://climatemodeling.science.energy.gov/projects/accelerated-climate-modeling-energy). The goal of this work is to develop a comprehensive diagnostics package that:

* fully integrates the functionality of NCAR's AMWG diagnostics package.
* utilizes most updated observational datasets, including remote sensing, reanalysis and in-situ datasets. 
* interfaces with diagnostics developed from different ACME focus groups: atmosphere group, coupled simulation group, land group.
* interacts effectively with the PCMDI's metrics package and the ARM diagnostics package through a unifying framework: [Community Diagnostics Package (CDP)](https://github.com/UV-CDAT/cdp).
* is flexible for user specified diagnostics and being configured for use by other climate models.

## Current State <a name="current-state"></a>
Algorithm and visulization codes for the AMWG Set 5, 7, 4, 3, 13, 1, 14 diagnostics, namely lat-lon contour plots (Figure 1), polar contour plots (Figure 2), zonal mean 2d plots (Figure 3), zonal mean line plots (Figure 4), 2d joint histogram for COSP cloud simulator output (Figure 5), tables (Figure 6) and Taylor Diagrams (Figure 7) summarizing metrics, for climatology seasonal means, are implemented. 

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
