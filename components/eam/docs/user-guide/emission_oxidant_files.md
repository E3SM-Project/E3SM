
# Emission files for EAMv3 gas and aerosol species and oxidant file for VBS SOA and stratosphere sulfate formation

## Overview

This page documents emissions data for all required aerosols and precursor gases as well as oxidants input data for running EAMv3, mostly for the MAM4 aerosol scheme, from anthropogenic (i.e., industrial, energy, transportation, domestic, and agriculture activity sectors) and natural (i.e., sea spray, vegetation, fire smoke, volcano) sources. Sulfur from agriculture, domestic, transportation, waste, and shipping sectors is emitted at the surface while sulfur from energy and industry sectors is emitted at 100-300 m above the surface, and sulfur from forest fire and grass fire is emitted at higher elevations (0-6 km). POM and BC from forest fire and grass fire are emitted at 0-6 km, while those from other sources (domestic, energy, industry, transportation, waste, and shipping) are emitted at the surface. Injection height profiles for fire emissions are derived from the corresponding AeroCom profiles (Dentener et al., 2006)[@dentener_emissions_2006], which give emissions in 6 altitude ranges (0-0.1, 0.1-0.5, 0.5-1, 1-2, 2-3, and 3-6 km). Otherwise, species emissions are assumed to be at the surface (bottom layer). Number emission fluxes each mode are calculated from mass emission fluxes based on AeroCom prescribed lognormal size distributions.

## Aerosols and gas precursors (common for EAMv1/v2/v3)

* Species: SO2, SOAG0, DMS, bc_a4, pom_a4, so4_a1, so4_a2, num_a1, num_a2, num_a4
* Data sources
  * Most of the original data are directly from input4MIPs with the following exceptions (E3SM specific treatments)
  * SO2 takes 97.5% from the input4MIPs data (all SO2-em-anthro_input4MIPs sectors)
  * SO4_a1 surf takes 2.5% from the corresponding surface sectors of input4MIPs data (SO2-em-anthro_input4MIPs sectors: AGR, SLV, WST, SHP)
  * SO4_a2 surf takes 2.5% from the corresponding surface sectors of input4MIPs data (SO2-em-anthro_input4MIPs sectors: TRA, RCO)
  * SO4_a1 elev takes 2.5% from the corresponding elevated sectors of input4MIPs data (SO2-em-anthro_input4MIPs sectors: ENE, IND)
  * SO2 (97.5%) and SO4_a1 (2.5%) also take emissions from AR5 input file for sector contvolc (constant volcanic degassing)
  * SOAG0 emissions are obtained by scaling OC (POM) emissions with a tunable factor.
  * num_a1, num_a2 and num_a4 are determined by mass concentration of aerosols species in the corresponding sectors and modes.

## Marine organic sea spray

Marine organic sea spray aerosol contributions are parameterized following the OCEANFILMS parameterization (E3SMv1; Burrows et al., 2014;[@burrows_physically_2014] 2022[@burrows_oceanfilms_2022]).  The input file for this parameterization provides a climatology of the ocean surface concentrations of several groups of organic macromolecules.  Briefly, the Parallel Ocean Program (POP; Maltrud et al., 1998)[@maltrud_global_1998] and its biogeochemical elemental cycling (BEC) routines (Moore et al., 2004)[@moore_upper_2004] were used to simulate marine biogeochemistry fields, including particulate organic matter (POC), chlorophyll, and zooplankton concentrations; these fields were used to generate maps of the estimated surface distributions of classes of macromolecules following the methods described in Burrows et al. (2014).[@burrows_physically_2014]  The scripts used to accomplish this translation are available [here](https://github.com/E3SM-Project/PreAndPostProcessingScripts/blob/devel/prepare_model_inputfiles/emis/marine_organic_aerosol/JAN_1850_MASTERLANG.jnl).

The file used as an input to E3SM is available here:
[https://web.lcrc.anl.gov/public/e3sm/inputdata/atm/cam/chem/trop_mam/marine_BGC/monthly_macromolecules_0.1deg_bilinear_latlon_year01_merge_date.nc](https://web.lcrc.anl.gov/public/e3sm/inputdata/atm/cam/chem/trop_mam/marine_BGC/monthly_macromolecules_0.1deg_bilinear_latlon_year01_merge_date.nc)

And is also published as a citeable dataset on Zenodo:

Elliott, S. M., Maltrud, M., & Burrows, S. M. (2015). Macromolecule distributions input file for the OCEANFILMS parameterization (Version v1) [Data set]. [Zenodo](https://doi.org/10.5281/zenodo.6320812).

## Oceanic dimethyl sulfide concentrations

Dimethyl sulfide (DMS) fluxes to the atmosphere are calculated in E3SM as a function of prescribed surface oceanic DMS concentrations, and an air-sea flux piston velocity that is a function of wind speed.

E3SM uses a DMS surface concentration dataset developed from a dynamic ocean biogeochemistry simulation; the methods and underlying assumptions used to produce this dataset are documented in Wang, et al. (2015).  The resolution in the DMS dataset is 1.9x2.5 degrees.

## New gas species for chemUCI in EAMv3

* Species: C2H4, C2H6, C3H8, CH2O, CH3CHO, CH3COCH3, CO, E90, ISOP, NO, NO2
* Data sources
  * anthropogenic, biomass burning, and aircraft emissions (NO2) are regridded from NCAR CESM2 emission files. They  are time-dependent during historical period and in the future scenarios.
  * biogenic emissions (C2H4, C2H6, C3H8, CH2O, CH3CHO, CH3COCH3, CO, ISOP) are from MEGAN-MACC offline data
    * 1850-1979: monthly input cycled yearly from 30-year mean (1980-2009)
    * 1980-2014: time-varying MEGAN-MACC data (historical)
    * 2015-2100: monthly input cycled yearly from 30-year mean (1980-2009)
  * natural emissions from oceans (C2H4, C2H6, C3H8, CO) and soil (NO) are regridded from NCAR CESM2 emission files. Just cycled yearly during the historical period and in the future scenarios.
  * E90 emissions?

## Oxidants file needed for VBS SOA and stratosphere sulfate formation

* Species: prsd_O3, prsd_NO3, prsd_OH
* Data sources
  * the file for historical simulation is the same as v1 and v2, inherited from CESM
  * files for SSPs are regridded from NCAR CESM2 tracer files

## Namelist setting for emissions input

### Historical (WCYCL20TR)/AMIP (F20TR)

```fortran
 ext_frc_specifier              = 'NO2         -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_NO2_aircraft_vertical_1750-2015_1.9x2.5_c20170608.nc',
         'SO2         -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so2_elev_1850-2014_c180205_kzm_1850_2014_volcano.nc',
         'SOAG0       -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/emissions-cmip6_e3sm_SOAG0_elev_1850-2014_1.9x2.5_c20230201.nc',
         'bc_a4       -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_bc_a4_elev_1850-2014_c180205.nc',
         'num_a1      -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a1_elev_1850-2014_c180205.nc',
         'num_a2      -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a2_elev_1850-2014_c180205.nc',
         'num_a4      -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a4_elev_1850-2014_c180205.nc',
         'pom_a4      -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_pom_a4_elev_1850-2014_c180205.nc',
         'so4_a1      -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so4_a1_elev_1850-2014_c180205.nc',
         'so4_a2      -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so4_a2_elev_1850-2014_c180205.nc'
         
 srf_emis_specifier             = 'C10H16 -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_MTERP_surface_1850-2014_1.9x2.5_c20230126.nc',
         'C2H4      -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_C2H4_surface_1850-2014_1.9x2.5_c20210323.nc',
         'C2H6      -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_C2H6_surface_1850-2014_1.9x2.5_c20210323.nc',
         'C3H8      -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_C3H8_surface_1850-2014_1.9x2.5_c20210323.nc',
         'CH2O   -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_CH2O_surface_1850-2014_1.9x2.5_c20210323.nc',
         'CH3CHO    -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_CH3CHO_surface_1850-2014_1.9x2.5_c20210323.nc',
         'CH3COCH3  -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_CH3COCH3_surface_1850-2014_1.9x2.5_c20210323.nc',
         'CO     -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_CO_surface_1850-2014_1.9x2.5_c20210323.nc',
         'DMS    -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DMSflux.1850-2100.1deg_latlon_conserv.POPmonthlyClimFromACES4BGC_c20160727.nc',
         'E90       -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions_E90_surface_1750-2015_1.9x2.5_c20210408.nc',
         'ISOP   -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_ISOP_surface_1850-2014_1.9x2.5_c20210323.nc',
         'ISOP_VBS -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_ISOP_surface_1850-2014_1.9x2.5_c20210323.nc',
         'NO     -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_NO_surface_1850-2014_1.9x2.5_c20220425.nc',
         'SO2    -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so2_surf_1850-2014_c180205.nc',
         'SOAG0  -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/emissions-cmip6_e3sm_SOAG0_surf_1850-2014_1.9x2.5_c20230201.nc',
         'bc_a4  -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_bc_a4_surf_1850-2014_c180205.nc',
         'num_a1 -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a1_surf_1850-2014_c180205.nc',
         'num_a2 -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a2_surf_1850-2014_c180205.nc',
         'num_a4 -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a4_surf_1850-2014_c180205.nc',
         'pom_a4 -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_pom_a4_surf_1850-2014_c180205.nc',
         'so4_a1 -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so4_a1_surf_1850-2014_c180205.nc',
         'so4_a2 -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so4_a2_surf_1850-2014_c180205.nc'
         
 tracer_cnst_datapath           = '\$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/oxid'
 tracer_cnst_file               = 'oxid_1.9x2.5_L26_1850-2015_c20181106.nc'
 tracer_cnst_filelist           = ''
 tracer_cnst_specifier          = 'prsd_O3:O3','prsd_NO3:NO3','prsd_OH:OH'
 tracer_cnst_type               = 'INTERP_MISSING_MONTHS'
```

### F2010

```fortran
 ext_frc_cycle_yr               = 2010
 ext_frc_specifier              = 'NO2         -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_NO2_aircraft_vertical_2010_clim_1.9x2.5_c20230213.nc',
         'SO2         -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so2_elev_1x1_2010_clim_c20190821.nc',
         'SOAG0       -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/emissions-cmip6_e3sm_SOAG0_elev_2010_clim_1.9x2.5_c20230213.nc',
         'bc_a4       -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_bc_a4_elev_1x1_2010_clim_c20190821.nc',
         'num_a1      -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a1_elev_1x1_2010_clim_c20190821.nc',
         'num_a2      -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a2_elev_1x1_2010_clim_c20190821.nc',
         'num_a4      -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a4_elev_1x1_2010_clim_c20190821.nc',
         'pom_a4      -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_pom_a4_elev_1x1_2010_clim_c20190821.nc',
         'so4_a1      -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so4_a1_elev_1x1_2010_clim_c20190821.nc',
         'so4_a2      -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so4_a2_elev_1x1_2010_clim_c20190821.nc'
 ext_frc_type           = 'CYCLICAL'
 
 srf_emis_cycle_yr              = 2010
 srf_emis_specifier             = 'C10H16 -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_MTERP_surface_2010_clim_1.9x2.5_c20230213.nc',
         'C2H4      -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_C2H4_surface_2010_clim_1.9x2.5_c20230213.nc',
         'C2H6      -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_C2H6_surface_2010_clim_1.9x2.5_c20230213.nc',
         'C3H8      -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_C3H8_surface_2010_clim_1.9x2.5_c20230213.nc',
         'CH2O   -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_CH2O_surface_2010_clim_1.9x2.5_c20230213.nc',
         'CH3CHO    -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_CH3CHO_surface_2010_clim_1.9x2.5_c20230213.nc',
         'CH3COCH3  -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_CH3COCH3_surface_2010_clim_1.9x2.5_c20230213.nc',
         'CO     -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_CO_surface_2010_clim_1.9x2.5_c20230213.nc',
         'DMS    -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DMSflux.2010.1deg_latlon_conserv.POPmonthlyClimFromACES4BGC_c20190220.nc',
         'E90       -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions_E90_surface_2010_clim_1.9x2.5_c20230213.nc',
         'ISOP   -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_ISOP_surface_2010_clim_1.9x2.5_c20230213.nc',
         'ISOP_VBS -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_ISOP_surface_2010_clim_1.9x2.5_c20230213.nc',
         'NO     -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/chem_gases/2degrees/emissions-cmip6_e3sm_NO_surface_2010_clim_1.9x2.5_c20230213.nc',
         'SO2    -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so2_surf_1x1_2010_clim_c20190821.nc',
         'SOAG0  -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/emissions-cmip6_e3sm_SOAG0_surf_2010_clim_1.9x2.5_c20230213.nc',
         'bc_a4  -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_bc_a4_surf_1x1_2010_clim_c20190821.nc',
         'num_a1 -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a1_surf_1x1_2010_clim_c20190821.nc',
         'num_a2 -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a2_surf_1x1_2010_clim_c20190821.nc',
         'num_a4 -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a4_surf_1x1_2010_clim_c20190821.nc',
         'pom_a4 -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_pom_a4_surf_1x1_2010_clim_c20190821.nc',
         'so4_a1 -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so4_a1_surf_1x1_2010_clim_c20190821.nc',
         'so4_a2 -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so4_a2_surf_1x1_2010_clim_c20190821.nc'
 srf_emis_type          = 'CYCLICAL'
 
 tracer_cnst_cycle_yr           = 2015
 tracer_cnst_datapath           = '\$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/oxid'
 tracer_cnst_file               = 'oxid_1.9x2.5_L26_1850-2015_c20181106.nc'
 tracer_cnst_filelist           = ''
 tracer_cnst_specifier          = 'prsd_O3:O3','prsd_NO3:NO3','prsd_OH:OH'
 tracer_cnst_type               = 'CYCLICAL'
 ```

### Future Scenarios

#### SSP370

```fortran
ext_frc_specifier              = 'NO2         -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/emissions-cmip6_ssp370_e3sm_NO2_aircraft_vertical_2015-2100_1.9x2.5_c20240208.nc',
         'SO2         -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_so2_elev_2015-2100_c210216.nc',
         'SOAG0       -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/emissions-cmip6_ssp370_e3sm_SOAG0_elev_2015-2100_1.9x2.5_c20240208.nc',
         'bc_a4       -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_bc_a4_elev_2015-2100_c210216.nc',
         'num_a1      -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_num_a1_elev_2015-2100_c210216.nc',
         'num_a2      -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_num_a2_elev_2015-2100_c210216.nc',
         'num_a4      -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_num_a4_elev_2015-2100_c210216.nc',
         'pom_a4      -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_pom_a4_elev_2015-2100_c210216.nc',
         'so4_a1      -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_so4_a1_elev_2015-2100_c210216.nc',
         'so4_a2      -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_so4_a2_elev_2015-2100_c210216.nc'
 
srf_emis_specifier             = 'C10H16 -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/emissions-cmip6_ssp370_e3sm_MTERP_surface_2015-2100_1.9x2.5_c20240208.nc',
         'C2H4      -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/emissions-cmip6_ssp370_e3sm_C2H4_surface_2015-2100_1.9x2.5_c20240208.nc',
         'C2H6      -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/emissions-cmip6_ssp370_e3sm_C2H6_surface_2015-2100_1.9x2.5_c20240208.nc',
         'C3H8      -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/emissions-cmip6_ssp370_e3sm_C3H8_surface_2015-2100_1.9x2.5_c20240208.nc',
         'CH2O      -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/emissions-cmip6_ssp370_e3sm_CH2O_surface_2015-2100_1.9x2.5_c20240208.nc',
         'CH3CHO    -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/emissions-cmip6_ssp370_e3sm_CH3CHO_surface_2015-2100_1.9x2.5_c20240208.nc',
         'CH3COCH3  -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/emissions-cmip6_ssp370_e3sm_CH3COCH3_surface_2015-2100_1.9x2.5_c20240208.nc',
         'CO        -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/emissions-cmip6_ssp370_e3sm_CO_surface_2015-2100_1.9x2.5_c20240208.nc ',
         'DMS       -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/DMSflux.1850-2100.1deg_latlon_conserv.POPmonthlyClimFromACES4BGC_c20160727.nc',
         'E90       -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart/ub/emissions_E90_surface_1750-2101_1.9x2.5_c20231222.nc',
         'ISOP      -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/emissions-cmip6_ssp370_e3sm_ISOP_surface_2015-2100_1.9x2.5_c20240208.nc',
         'ISOP_VBS  -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/emissions-cmip6_ssp370_e3sm_ISOP_surface_2015-2100_1.9x2.5_c20240208.nc',
         'NO        -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/emissions-cmip6_ssp370_e3sm_NO_surface_2015-2100_1.9x2.5_c20240208.nc',
         'SO2       -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_so2_surf_2015-2100_c210216.nc',
         'SOAG0     -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/emissions-cmip6_ssp370_e3sm_SOAG0_surf_2015-2100_1.9x2.5_c20240208.nc',
         'bc_a4     -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_bc_a4_surf_2015-2100_c210216.nc',
         'num_a1    -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_num_a1_surf_2015-2100_c210216.nc',
         'num_a2    -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_num_a2_surf_2015-2100_c210216.nc',
         'num_a4    -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_num_a4_surf_2015-2100_c210216.nc',
         'pom_a4    -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_pom_a4_surf_2015-2100_c210216.nc',
         'so4_a1    -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_so4_a1_surf_2015-2100_c210216.nc',
         'so4_a2    -> \$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/emis/CMIP6_SSP370_ne30/cmip6_ssp370_mam4_so4_a2_surf_2015-2100_c210216.nc'
         
 tracer_cnst_datapath           = '\$DIN_LOC_ROOT/atm/cam/chem/trop_mozart_aero/oxid'
 tracer_cnst_file               = 'oxid_SSP370_1.9x2.5_L70_1849-2101_c20240228.nc'
 tracer_cnst_filelist           = ''
 tracer_cnst_specifier          = 'prsd_O3:O3','prsd_NO3:NO3','prsd_OH:OH'
 tracer_cnst_type               = 'INTERP_MISSING_MONTHS'
```
