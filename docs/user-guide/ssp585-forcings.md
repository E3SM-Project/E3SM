# SSP585 Forcing data

These are the prescribed inputdata specifically for the SSP585 scenario, in place of the files for the historical period.

## Solar constant

`\$DIN_LOC_ROOT/atm/cam/solar/Solar_1850-2299_input4MIPS_c20181106.nc`

## Greenhouse gas concentrations

`\$DIN_LOC_ROOT/atm/cam/ggas/GHG_CMIP_SSP585-1-2-1_Annual_Global_2015-2500_c20190310.nc`

## Elevated external forcings

```fortran
no2_ext_file         \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/emissions-cmip6_ssp585_e3sm_NO2_aircraft_vertical_2015-2100_1.9x2.5_c20240304.nc 
so2_ext_file         \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/cmip6_ssp585_mam4_so2_volc_elev_2015-2100_c240331.nc 
soag0_ext_file       \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/emissions-cmip6_ssp585_e3sm_SOAG0_elev_2015-2100_1.9x2.5_c20240304.nc 
bc_a4_ext_file       \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/cmip6_ssp585_mam4_bc_a4_elev_2015-2100_c190828.nc 
mam7_num_a1_ext_file \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/cmip6_ssp585_mam4_num_a1_elev_2015-2100_c190828.nc 
num_a2_ext_file      \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/cmip6_ssp585_mam4_num_a2_elev_2015-2100_c190828.nc 
mam7_num_a3_ext_file \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/cmip6_ssp585_mam4_num_a4_elev_2015-2100_c190828.nc 
pom_a4_ext_file      \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/cmip6_ssp585_mam4_pom_a4_elev_2015-2100_c190828.nc 
so4_a1_ext_file      \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/cmip6_ssp585_mam4_so4_a1_elev_2015-2100_c190828.nc 
so4_a2_ext_file      \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/cmip6_ssp585_mam4_so4_a2_elev_2015-2100_c190828.nc 
```

## Surface emissions

```fortran
c2h4_emis_file        \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/emissions-cmip6_ssp585_e3sm_C2H4_surface_2015-2100_1.9x2.5_c20240304.nc 
c2h6_emis_file        \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/emissions-cmip6_ssp585_e3sm_C2H6_surface_2015-2100_1.9x2.5_c20240304.nc 
c3h8_emis_file        \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/emissions-cmip6_ssp585_e3sm_C3H8_surface_2015-2100_1.9x2.5_c20240304.nc 
ch2o_emis_file        \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/emissions-cmip6_ssp585_e3sm_CH2O_surface_2015-2100_1.9x2.5_c20240304.nc 
ch3cho_emis_file      \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/emissions-cmip6_ssp585_e3sm_CH3CHO_surface_2015-2100_1.9x2.5_c20240304.nc 
ch3coch3_emis_file    \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/emissions-cmip6_ssp585_e3sm_CH3COCH3_surface_2015-2100_1.9x2.5_c20240304.nc 
co_emis_file          \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/emissions-cmip6_ssp585_e3sm_CO_surface_2015-2100_1.9x2.5_c20240304.nc 
isop_emis_file        \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/emissions-cmip6_ssp585_e3sm_ISOP_surface_2015-2100_1.9x2.5_c20240304.nc 
isop_vbs_emis_file    \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/emissions-cmip6_ssp585_e3sm_ISOP_surface_2015-2100_1.9x2.5_c20240304.nc 
c10h16_emis_file      \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/emissions-cmip6_ssp585_e3sm_MTERP_surface_2015-2100_1.9x2.5_c20240304.nc 
nox_emis_file         \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/emissions-cmip6_ssp585_e3sm_NO_surface_2015-2100_1.9x2.5_c20240304.nc 
dms_emis_file         \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DMSflux.1850-2100.1deg_latlon_conserv.POPmonthlyClimFromACES4BGC_c20160727.nc
so2_emis_file         \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/cmip6_ssp585_mam4_so2_surf_2015-2100_c190828.nc 
soag0_emis_file       \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/emissions-cmip6_ssp585_e3sm_SOAG0_surf_2015-2100_1.9x2.5_c20240304.nc 
bc_a4_emis_file       \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/cmip6_ssp585_mam4_bc_a4_surf_2015-2100_c190828.nc 
mam7_num_a1_emis_file \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/cmip6_ssp585_mam4_num_a1_surf_2015-2100_c190828.nc 
num_a2_emis_file      \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/cmip6_ssp585_mam4_num_a2_surf_2015-2100_c190828.nc 
mam7_num_a3_emis_file \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/cmip6_ssp585_mam4_num_a4_surf_2015-2100_c190828.nc 
pom_a4_emis_file      \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/cmip6_ssp585_mam4_pom_a4_surf_2015-2100_c190828.nc 
so4_a1_emis_file      \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/cmip6_ssp585_mam4_so4_a1_surf_2015-2100_c190828.nc 
so4_a2_emis_file      \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP585_ne30/cmip6_ssp585_mam4_so4_a2_surf_2015-2100_c190828.nc 
e90_emis_file         \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart/ub/emissions_E90_surface_1750-2101_1.9x2.5_c20231222.nc 
```

## Prescribed oxidant for aerosol chemistry

`\$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/oxid/oxid_SSP585_1.9x2.5_L70_2014-2101_c20240228.nc`

## Stratospheric ozone (linoz) and chlorine loading data

```fortran
\$DIN_LOC_ROOT/atm/cam/chem/trop_mozart/ub/Linoz_Chlorine_Loading_CMIP6_Hist_SSP585_0003-2503_c20190414.nc
\$DIN_LOC_ROOT/atm/cam/chem/trop_mozart/ub/linv3_1849-2101_CMIP6_Hist_SSP585_10deg_58km_c20230705.nc
```

## Land use and land cover

`\$DIN_LOC_ROOT/lnd/clm2/surfdata_map/landuse.timeseries_0.5x0.5_ssp5_rcp85_simyr2015-2100_c240408.nc`
