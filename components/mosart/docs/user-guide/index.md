# MOSART User's Guide

This guide describes how to set up and run MOSART.

## Steps to build and run MOSART

A step-by-step instruction on how to run fully coupled E3SM can be found [here](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/2309226536). Here we describe running MOSART driven by runoff forcings provided via the data land mode(DLND). Although the default runoff forcing is predefined by the DLND compsets, the users is able to use other runoff forcing to drive MOSART.

MOSART in the water cycle campaigns of E3SM v1, v2, and v3 was.

### Scientifically supported compsets

The mosart-only compsets are supported for multiple runoff forcing datasets that covers different domains and time periods:

1. `RMOSQIAN`: A 57-year (1948-2004, no leap year) global Qian runoff forcings. The native resolutions are 192*288 at daily time scale. The path to the stream file is `$input_data_dir/lnd/dlnd7/hcru_hcru/QIAN.daily.1948-2004.nc`
2. `RMOSGPCC`: A one-year (1979) global GPCC runoff forcing. The native resoultions are 192*288 at daily time scale. The path to the stream file is `$input_data_dir/lnd/dlnd7/hcru_hcru/GPCC.daily.nc`
3. `RMOSNLDAS`: A one-year (1975) CLM runoff forcing for [NLDAS](https://ldas.gsfc.nasa.gov/nldas) domain (North America between 25N and 53N) and spatial resolution (1/8th degree). The time scale is 3-hourly. The path to the stream file is `$input_data_dir/lnd/dlnd7/NLDAS/mosart_nldas_hist_1975.clm2.h1.1975-01-01-00000.nc`.

### Supported grid

The `r05_r05`, `NLDAS_NLDAS`, and `MOS_USRDAT` are the supported grid resolutions for performing offline MOSART simulations.

### User-defined runoff forcing

Once the case is created based on the above compsets. The users can create a `user_dlnd.streams.txt` file in the `case_scripts` directory following the format below:

```xml
<?xml version="1.0"?>
<file id="stream" version="1.0">
<dataSource>
  GENERIC
</dataSource>
<domainInfo>
  <variableNames>
    time     time
    lon      lon
    lat      lat
  </variableNames>
  <filePath>
    /path/to/forcing/file
  </filePath>
  <fileNames>
    Runoff_forcing.nc
  </fileNames>
</domainInfo>
<fieldInfo>
   <variableNames>
    QDRAI     rofsub
    QOVER     rofsur
   </variableNames>
   <filePath>
    /path/to/forcing/file
   </filePath>
   <fileNames>
    Runoff_forcing.nc
   </fileNames>
   <offset>
    0
   </offset>
</fieldInfo>
</file>
```

## Customizing runs

For default river routing simulation in MOSART (i.e. natural flow routing), one [parameter file](../tech-guide/index.md#parameters) defined by `frivinp_rtm` is required. The geographic domain and spatial resolution must match the domain defined in the simulation case. Additional parameter/input files will be needed for some extra features.

### Output variables and frequency

MOSART by default outputs monthly history file in `*mosart.h0.**.nc` files that include all key variables. User could choose to output additional history files (such as `*mosart.h1.*.nc`, `*mosart.h2.*.nc`) that have different temporal averaging (e.g. daily, hourly, every model timestep) via `rtmhist_nhtfrq` where

- `-24` corresponds to daily average
- `-1` corresponds to hourly average
- `0` corresponds to monthly average
- `1` corresponds to each model time step

The number of time slices in these additional files can be controlled
vai `rtmhist_mfilt`.

```fortran
   &mosart_inparm
     rtmhist_fincl2 = 'FLOODPLAIN_FRACTION'
     rtmhist_fincl3 = 'RIVER_DISCHARGE_OVER_LAND_LIQ'
     rtmhist_nhtfrq = 0, -24, -1
     rtmhist_mfilt  = 12, 30, 24
   /
```

Using the above-mentioned settings:

- Each `*.mosart.h1.*.nc` will include 30 daily average values of `FLOODPLAIN_FRACTION`
- Each `*.mosart.h2.*.nc` will include 24 hourly average values of `RIVER_DISCHARGE_OVER_LAND_LIQ`

### Additional options

The table below lists avaiable options defined by users through `user_nl_mosart`.

| Flag Name         | Description                                                                           |
|-------------------|---------------------------------------------------------------------------------------|
| `routingmethod`   | `1` for kenematic wave routing (default); `2` for diffusion wave routing              |
| `wrmflag`         | `.true.` for active water management; `.false.` for no water management (default)     |
| `inundflag`       | `.true.` for active flood inundation; `.false.` for no flood inundation (default)     |
| `sediflag`        | `.true.` for active sediment transport; `.false.` for no sediment transport (default) |
| `heatflag`        | `.true.` for active heat transport; `.false.` for no heat transport (default)         |

## Additional MOSART features

There are some options only made available for specific features. They can be defined through `user_nl_mosart`.

### Water management

- Parameter file: when water management is active, one additional parameter file `parafile` is required. This file defines the location and specifics for the dams/reservoirs simulated in this scheme.

- `damconstructionflag`: `0` - dam always exist; `1` - dam construction year is considered.

- `externaldemandflag`: `0` - no external water demand for the WM scheme; `1` - external water demand files are required.

  - Note if `externaldemandflag` is set to `1`, paths to monthly water demand files are requried in the `usr_nl_mosart` file through `demandpath = '/path/to/demand/files/`.

### Flood inundation

- `opt_eleprof`: `1` - use actural elevation profiles derived from DEM; `2` - use hypothetical elevation profile.

  - Note if `opt_eleprof` is set to `1`, the elevation profile data must be included in the MOSART parameter file.

### Sediment transport

- If sediment feature is activated, D50 data must be included in the MOSART parameter file. In addition, `rof_sed = .true.` has to be defined in `./user_nl_cpl` to allow sediment flux passing into the river model through the coupler.
