# SSP370 Forcing data

These are the prescribed inputdata specifically for the SSP370 scenario, in place of the files for the historical period.

## Solar constant

`\$DIN_LOC_ROOT/atm/cam/solar/Solar_1850-2299_input4MIPS_c20181106.nc`

## Greenhouse gas concentrations

`\$DIN_LOC_ROOT/atm/cam/ggas/GHG_CMIP_SSP370-1-2-1_Annual_Global_2015-2500_c20210509.nc`

## Elevated external forcings

```fortran
NO2         \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/emissions-cmip6_ssp370_e3sm_NO2_aircraft_vertical_2015-2100_1.9x2.5_c20240208.nc 
SO2         \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_so2_volc_elev_2015-2100_c240331.nc 
SOAG0       \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/emissions-cmip6_ssp370_e3sm_SOAG0_elev_2015-2100_1.9x2.5_c20240208.nc 
bc_a4       \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_bc_a4_elev_2015-2100_c210216.nc 
num_a1      \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_num_a1_elev_2015-2100_c210216.nc 
num_a2      \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_num_a2_elev_2015-2100_c210216.nc 
num_a4      \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_num_a4_elev_2015-2100_c210216.nc 
pom_a4      \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_pom_a4_elev_2015-2100_c210216.nc 
so4_a1      \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_so4_a1_elev_2015-2100_c210216.nc 
so4_a2      \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_so4_a2_elev_2015-2100_c210216.nc 
```

## Surface emissions

```fortran
C2H4        \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/emissions-cmip6_ssp370_e3sm_C2H4_surface_2015-2100_1.9x2.5_c20240208.nc 
C2H6        \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/emissions-cmip6_ssp370_e3sm_C2H6_surface_2015-2100_1.9x2.5_c20240208.nc 
C3H8        \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/emissions-cmip6_ssp370_e3sm_C3H8_surface_2015-2100_1.9x2.5_c20240208.nc 
CH2O        \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/emissions-cmip6_ssp370_e3sm_CH2O_surface_2015-2100_1.9x2.5_c20240208.nc 
CH3CHO      \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/emissions-cmip6_ssp370_e3sm_CH3CHO_surface_2015-2100_1.9x2.5_c20240208.nc 
CH3COCH3    \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/emissions-cmip6_ssp370_e3sm_CH3COCH3_surface_2015-2100_1.9x2.5_c20240208.nc 
CO          \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/emissions-cmip6_ssp370_e3sm_CO_surface_2015-2100_1.9x2.5_c20240208.nc 
ISOP        \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/emissions-cmip6_ssp370_e3sm_ISOP_surface_2015-2100_1.9x2.5_c20240208.nc 
ISOP_VBS    \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/emissions-cmip6_ssp370_e3sm_ISOP_surface_2015-2100_1.9x2.5_c20240208.nc 
C10H16      \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/emissions-cmip6_ssp370_e3sm_MTERP_surface_2015-2100_1.9x2.5_c20240208.nc 
NOX         \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/emissions-cmip6_ssp370_e3sm_NO_surface_2015-2100_1.9x2.5_c20240208.nc 
DMS         \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DMSflux.1850-2100.1deg_latlon_conserv.POPmonthlyClimFromACES4BGC_c20160727.nc
SO2         \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_so2_surf_2015-2100_c210216.nc 
SOAG0       \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/emissions-cmip6_ssp370_e3sm_SOAG0_surf_2015-2100_1.9x2.5_c20240208.nc 
bc_a4       \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_bc_a4_surf_2015-2100_c210216.nc 
num_a1      \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_num_a1_surf_2015-2100_c210216.nc 
num_a2      \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_num_a2_surf_2015-2100_c210216.nc 
num_a4      \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_num_a4_surf_2015-2100_c210216.nc 
pom_a4      \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_pom_a4_surf_2015-2100_c210216.nc 
so4_a1      \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_so4_a1_surf_2015-2100_c210216.nc 
so4_a2      \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_so4_a2_surf_2015-2100_c210216.nc 
E90         \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart/ub/emissions_E90_surface_1750-2101_1.9x2.5_c20231222.nc 
```

## Prescribed oxidant for aerosol chemistry

`\$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/oxid/oxid_SSP370_1.9x2.5_L70_1849-2101_c20240228.nc`

## Stratospheric ozone (linoz) and chlorine loading data

```fortran
\$DIN_LOC_ROOT/atm/cam/chem/trop_mozart/ub/Linoz_Chlorine_Loading_CMIP6_Hist_SSP370_0003-2503_c20210202.nc
\$DIN_LOC_ROOT/atm/cam/chem/trop_mozart/ub/linv3_1849-2101_CMIP6_Hist_SSP370_10deg_58km_c20230705.nc
```

## Land use and land cover

`\$DIN_LOC_ROOT/lnd/clm2/surfdata_map/landuse.timeseries_0.5x0.5_ssp3_rcp70_simyr2015-2100_c240308.nc`
