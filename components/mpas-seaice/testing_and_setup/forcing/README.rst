======================================
MPAS-Seaice forcing generation scripts
======================================

Overview
========

These scripts generate standalone atmospheric and oceanic forcing files for
MPAS-Seaice.

create_atmos_forcing.py
=======================

Usage
-----

This script creates atmospheric forcing using six hourly CORE-II data and
monthly AOMIP climatologies.

.. code::

   > python create_atmos_forcing.py configFilename

where configFilename is a python config file with the following example format:

.. code::

   [forcing_generation]
   filenameMPASGrid = /location/of/MPAS/grid
   outputDir = /location/to/put/output/forcing
   startYear = 1948
   endYear = 2007
   dataDirSixHourly = /location/of/CORE-II/data
   dataDirMonthly = /location/of/AOMIP/climatologies
   scripDir = /location/of/SCRIP/executable

SCRIP
-----

This script requires the SCRIP package to be installed.
SCRIP is a software package which computes addresses and weights for remapping
and interpolating fields between grids in spherical coordinates. It can be
obtained from https://github.com/SCRIP-Project/SCRIP

CORE-II data
------------

Six-hourly air temperature, velocity and specific humidity comes from CORE-II.
Data files can be obtained from
https://data1.gfdl.noaa.gov/nomads/forms/core/COREv2/CIAF_v2.html.
To generate forcing for a given year YYYY, the following files are required:

* ${dataDirSixHourly}/t_10/t_10.YYYY.*.nc

* ${dataDirSixHourly}/q_10/q_10.YYYY.*.nc

* ${dataDirSixHourly}/u_10/u_10.YYYY.*.nc

* ${dataDirSixHourly}/v_10/v_10.YYYY.*.nc

where ${dataDirSixHourly} is the local location of the six hourly data.

AOMIP climatologies
-------------------

Monthly climatologies of cloudiness and precipitation comes from AOMIP.
The following data files are required:

* ${dataDirMonthly}/cldf.omip.nc

* ${dataDirMonthly}/prec.nmyr.nc

where ${dataDirMonthly} is the local location of the monthly data.
These files can be obtained from https://web.lcrc.anl.gov/public/e3sm/mpas_standalonedata/mpas-seaice/forcing/MPAS-Seaice_forcing_gen_data.tar.gz


create_ocean_forcing.py
=======================

Usage
-----

This script creates ocean forcing using CESM output.

.. code::

   > python create_ocean_forcing.py configFilename

where configFilename is a python config file with the following example format:

.. code::

   [forcing_generation]
   filenameMPASGrid = /location/of/MPAS/grid
   filenameGx1Grid = /location/of/gx1/grid
   filenameGx1OceanMixed = /location/of/gx1/ocean_mixed_file
   filenameMPASOceanMixed = /location/of/output/ocean_mixed_file
   scripDir = /location/of/SCRIP/executable

SCRIP
-----

This script requires the SCRIP package to be installed.
SCRIP is a software package which computes addresses and weights for remapping
and interpolating fields between grids in spherical coordinates. It can be
obtained from https://github.com/SCRIP-Project/SCRIP

gx1 input data
--------------

This script requires a gx1 grid file and ocean mixed file as input. These can be
obtained from https://web.lcrc.anl.gov/public/e3sm/mpas_standalonedata/mpas-seaice/forcing/MPAS-Seaice_forcing_gen_data.tar.gz
