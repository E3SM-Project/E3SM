
# EAM User Guide

This User Guide describes how to set up and run EAM.

## Steps to build and run EAM

EAM is not available as a standalone model.  Instead, EAM can be run in an atmosphere-only configuration.
The difference when running in atmosphere-only mode (without an interactive ocean or sea-ice) would be to
change the *compset* and *grid*.
See the
[Case Control System Basic Usage](https://esmci.github.io/cime/versions/master/html/users_guide/index.html#case-control-system-part-1-basic-usage) for general descriptions of compsets and grids.

Certain namelist paramaters, input data files, and output file specifcations can also be modified.
These are described below as ways to customize runs.

Step-by-step instructions on how to run and analyze E3SM with a script can be found at
[E3SM step-by-step guide](https://docs.e3sm.org/running-e3sm-guide/)

## Scientifically supported compsets and grids

### Compsets

All of the compsets below run with the complete set of E3SM atmospheric configuration of EAMV3. For more information on the schemes in EAMv3, see the [Technical Guide](../tech-guide/index.md)

`F2010` - Climatological present day climate (year 2010)

`F1850` - Climatological pre-industrial-day climate (year 1850)

`F20TR` - Historical EAM simulation with time varying sea-surface temperatures, aerosol emissions, and greenhouse gas forcings (year 1850-2014)

### Grids

Only one grid combination is currently supported for the above compsets:

`ne30pg2_r05_IcoswISC30E3r5` - ne30pg2 atmosphere, 0.5deg x 0.5deg land grid, and Icosahedral 30 km mesh with ice shelves cavities (wISC), E3SMv3 (E3) revision r5

## Customizing runs

### Compile-time options

Some customizations require making changes before the model is built.

### Run-time options

Run-time customization is enabled by a Fortran namelist.

Namelist parameters can be changed from default values by putting them in the `user_nl_eam` file in the case directory
with the desired new value.

This [Table of Namelist Parameters](namelist_parameters.md) includes many of the paramaters that control
physics schemes described in the [tech-guide](../tech-guide/index.md)

#### History File Namelist Parameters

By default, EAM will output a set of monthly-averaged variables. Additional output files can be specified using the following flags in the `user_nl_eam` file:

`finclX` - List of variables (in single quotes and separated by commas) that are added to tape X.

`fexclX` - List of variables (in single quotes and separated by commas) that will be excluded in tape X.

`nhtfrq` - List of write frequencies for each of the history files. A value of 0 denotes a monthly frequency. Negative values denote hourly frequencies (e.g., `-3` will write an output every 3 hours). Positive values denotes the frequency in model timesteps (e.g., `4` will write an output every 4 timesteps).

`mfilt` - List that sets the number of timesteps to write in a single file before starting a new file.

`avgflag_pertape` - List that sets the type of output to write. Choices are `'A'` for time-averaged output, `'A'` for instantaneous output, `'M'` for time-minimum output, and `'X'` for time-maximum output.

#### Example output specification

```fortran
nhtfrq = 0,-24,-6,-3
mfilt  = 1,30,120,24
avgflag_pertape = 'A','A','A','I'

fexcl1 = 'U10' # Removes U10 output from monthly files
fincl2 = 'PS','FLUT','PRECT','U200','V200','U850',
          'V850','TCO','SCO','TREFHT','QREFHT'  # Output files of daily-averaged output, which includes 30 days of output in each file
fincl3 = 'PS', 'PSL','PRECT','TUQ','TVQ','UBOT',
          'VBOT','TREFHT','FLUT','OMEGA500','TBOT',
          'U850','V850','U200','V200','T200','T500',
          'Z700'  # Output files of 6-hour-averaged output, which includes 30 days of output in each file
fincl4 = 'PRECT' # Output files of 3-hourly output with 3 days of output in every file
```

## Input datasets

Many properties of the simulated atmosphere are controlled by data sets read in by the model during initialization.
Changing the content of these files is another way to customize a run.

### Greenhouse gases (non-reacting)

Greenhouse gas concentration inputs of non-reacting species are taken from CMIP6 Forcing Datasets provided from the input4MIPs data collection. In addition to what is provided by the input4MIPS, 2015 and 2016 have been added by extrapolating from 2013 and 2014.

```fortran
inputdata/atm/cam/ggas/GHG_CMIP-1-2-0_Annual_Global_0000-2014_c20180105.nc
```

### [Aerosol physical properties](aerosol_phys_prop.md)

The aerosol properties files provide aerosol refractive index, density, and aerosol hygroscopicty information for
each aerosol species, as well as information about lognormal mode definition and lookup tables of polynomial
expression coefficients for aerosol optics calculation for each mode. These aerosol physical and chemical properties
are used by the radiation, aerosol microphysics and other related source and sink processes, and droplet
activation/ice nucleation schemes.

### [Aerosol and gas emission and oxidant files](emission_oxidant_files.md)

Details of the aerosol and gas emission and oxidant files used in various historical, present-day, and future scenarios.

### Linoz v3 input files

Linozv3 uses the ozone tendency, (net production minus loss) calculated from its climatological mean state (function of month, latitude, altitude) and a first-order Taylor series expansion about the local ozone, temperature, and overhead ozone column. For historical simulation, Linozv3 uses the linozv3 data files with monthly resolution, spanning the dates 1849-01 -- 2014-12.

#### Historical files

```fortran
 linoz_data_file                = ‘linv3_1849-2101_CMIP6_Hist_10deg_58km_c20231207.nc’
 linoz_data_path                = '/lcrc/group/e3sm/data/inputdata/atm/cam/chem/trop_mozart/ub'
 linoz_data_type                = 'INTERP_MISSING_MONTHS'
```

Refer to [this page](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/3764486280/Production+of+the+Linoz+v3+data) for more details on Linoz v3 input files.

### Topography

The global elevation on the atmosphere grid is a key input dataset.  The dataset is grid dependent. It contains the geopotential data (on both the GLL dynamics and PG2 physics grids) and two surface roughness quantities, `SGH` and `SGH30`.

EAMv3 NE30 data:

```fortran
inputdata/atm/cam/topo/USGS-gtopo30_ne30np4pg2_x6t-SGH.c20210614.nc'
```

This file is computed via a complex procedure that starts with high resolution dataset (for EAMv3, we use the GTOPO30, a 30 arc-second resolution data set on a lat-lon grid) that is then downsampled to a 3km cubed-sphere grid (cube3000) and then downsampled to the atmosphere  grid on the (GLL nodes), and then smoothed with the same viscosity operator used by the dycore. The smoothed GLL topography is then mapped to the PG2 grid. Finally, two different surface roughness fields are computed:

- `SGH30`: the variance between GTOPO30 and GTOPO30-downsampled-to-PG2 (independent of any dycore specific smoothing).     (used by CLUBB, TMS, vertical_diffusion)
- `SGH`:  the variance between the cube3000 data and the smoothed PG2 data.  (used by GWD parameterizations)

### Land-use / land cover change

Info needed on land-use land cover change / land surface data

Refer to [ELM documentation](https://docs.e3sm.org/E3SM/ELM/).

### Ocean/sea ice

The sea surface temperature and sea-ice coverage data used in F-case simulations are based on CMIP6 Forcing Datasets provided from the input4MIPs data collection. The file for the `F2010` compset is created from taking a monthly climatology of years 2005-2014. The file for the `F1850` compset is created from taking a monthly climatology of years 1870-1879.*

*[Provenenance for F1850 file](https://acme-climate.atlassian.net/wiki/spaces/ATM/pages/201525378/Provenance+for+CMIP6+DECK+SST+Sea-Ice+input+data+for+E3SM)

`F20TR`

```fortran
inputdata/ocn/docn7/SSTDATA/sst_ice_CMIP6_DECK_E3SM_1x1_c20180213.nc
```

`F2010`

```fortran
inputdata/ocn/docn7/SSTDATA/sst_ice_CMIP6_DECK_E3SM_1x1_2010_clim_c20190821.nc 
```

`F1850`

```fortran
inputdata/ocn/docn7/SSTDATA/sst_ice_CMIP6_DECK_E3SM_1x1_1850_clim_c20190125.nc
```

### Solar input

As with greenhouse gas emissions, solar input files are taken from the input4MIPs data collection that were prepared for CMIP6 Forcing Datasets.

```fortran
inputdata/atm/cam/solar/Solar_1850-2299_input4MIPS_c20181106.nc
```
