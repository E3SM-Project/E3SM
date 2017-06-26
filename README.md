# ACME Diagnostics Package

## Table of Contents
1. [Overview](#overview)
2. [Current State](#current-state)
3. [Documentation](#doc)

## Overview<a name="overview"></a>
This diagnostics package is constructed for supporting the diagnostics task of DOE's [Accelerated Climate Modeling for Energy (ACME) project](https://climatemodeling.science.energy.gov/projects/accelerated-climate-modeling-energy). The goal of this work is to develop a comprehensive diagnostics package that:

* fully integrates the functionality of NCAR's AMWG diagnostics package.
* utilizes most updated observational datasets, including remote sensing, reanalysis and in-situ datasets. 
* interfaces with diagnostics developed from different ACME focus groups: atmosphere group, coupled simulation group, land group.
* interacts effectively with the PCMDI's metrics package and the ARM diagnostics package through a unifying framework: [Community Diagnostics Package (CDP)](https://github.com/UV-CDAT/cdp).
* is flexible for user specified diagnostics and being configured for use by other climate models.

## Current State <a name="current-state"></a>
Algorithm and visulization codes for the AMWG Set 5, 7, 4, 3 diagnostics, namely lat-lon contour plots (Figure 1), polar contour plots (Figure 2), pressure-lat contour plots (Figure 3), and zonal mean line plots (Figure 4) for climatology seasonal means, are done for testing. 

The package enables the AMWG Set 5 diagnostics and features built-in user diagnostics, by specifying user desired diagnostics regions and pressure levels for variables with the vertical dimension. 

<img src="docs/example_fig1.png" alt="Figure1" style="width: 280px;"/>
<h5 align="center">Figure 1: An example of the lat-lon contour plots for air temperature at 850 mb with tropical ocean region considered</h5> 

<img src="docs/example_fig2.png" alt="Figure2" style="width: 280px;"/>
<h5 align="center">Figure 2: An example of the polar contour plots for precipitation rate</h5> 

<img src="docs/example_fig3.png" alt="Figure3" style="width: 280px;"/>
<h5 align="center">Figure 3: An example of the pressure-lat contour plots for air temperature </h5> 

<img src="docs/example_fig4.png" alt="Figure4" style="width: 280px;"/>
<h5 align="center">Figure 4: An example of the zonal mean surface air temperature line plot </h5> 

## Documentation <a name="doc"></a>

* [Quickstart guide for AIMS4](docs/quick-guide-aims4.ipynb)
* [Installation, Basic Configuration, and Running](docs/install-config-run.ipynb)
* More Configuration
  * [Defining Parameters and All Available Parameters](docs/available-parameters.ipynb)
  * [Adding User Diagnostics](docs/add-new-diagnostics.ipynb)
  * [List of Built-in Derived Variables](acme_diags/derivations/acme.py) - look for `derived_variables` around line 190
* [CDP Viewer documentation](https://github.com/UV-CDAT/cdp/blob/master/jupyter/using-the-cdp-viewer.ipynb)
