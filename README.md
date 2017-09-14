# ACME Diagnostics Package

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
This diagnostics package is constructed for supporting the diagnostics task of DOE's [Accelerated Climate Modeling for Energy (ACME) project](https://climatemodeling.science.energy.gov/projects/accelerated-climate-modeling-energy). The goal of this work is to develop a comprehensive diagnostics package that:

* fully integrates the functionality of NCAR's AMWG diagnostics package.
* utilizes most updated observational datasets, including remote sensing, reanalysis and in-situ datasets. 
* interfaces with diagnostics developed from different ACME focus groups: atmosphere group, coupled simulation group, land group.
* interacts effectively with the PCMDI's metrics package and the ARM diagnostics package through a unifying framework: [Community Diagnostics Package (CDP)](https://github.com/UV-CDAT/cdp).
* is flexible for user specified diagnostics and being configured for use by other climate models.

## Current State <a name="current-state"></a>
Algorithm and visulization codes for the AMWG Set 5, 7, 4, 3, 13 diagnostics, namely lat-lon contour plots (Figure 1), polar contour plots (Figure 2), zonal mean 2d plots (Figure 3), zonal mean line plots (Figure 4), and 2d joint histogram for COSP cloud simulator output (Figure 5), for climatology seasonal means, are done for testing. 

The package features built-in user diagnostics, by specifying user desired diagnostics regions and pressure levels for variables with the vertical dimension. 

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
