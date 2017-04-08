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
Algorithm and visulization codes for the AMWG Set 5 diagnostics, namely 2D-field contour plots for climatology seasonal means Figure 1, is done for testing. 

The package enables the AMWG Set 5 diagnostics and features built-in user diagnostics, by specifying user desired diagnostics regions and pressure levels for variables with the vertical dimension. `

<img src="docs/example_fig1.png" alt="Figure1" style="width: 500px;"/>
<h5 align="center">Figure 1: An example of the 2D-field contour plots for air temperature at 850 mb with tropical ocean region considered</h5> 

## Documentation <a name="doc"></a>

* [Quickstart guide for AIMS4](docs/quick-guide-aims4.ipynb)
* [Installation, Basic Configuration, and Running](docs/install-config-run.ipynb)
* More Configuration
  * [List of parameters for config file](docs/available-parameters.ipynb)
  * [Adding user diagnostics](docs/add-new-diagnostics.ipynb)
  * [List of built-in derived variables](jill-add-this)
* [CDP Viewer documentation](https://github.com/UV-CDAT/cdp/blob/master/jupyter/using-the-cdp-viewer.ipynb)
