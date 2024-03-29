
# Aerosol physical properties

## Description

Key information

- Original physical properties of aerosols are documented in Liu et al. (2012). Detailed information is included in the supplement.

- Physical properties of aerosols used in E3SMv1 are documented in Wang et al. (2020).

- Marine organic aerosol properties are defined in Burrows et al. (2022).

- Dust refractive index and longwave absorption coefficients are updated by Feng et al. (2022).

- BC and POM hygroscopicity values are updated by Shan et al. (2024).

- <span style="color:red">Some description/reference about the 5th mode.</span>

## Included fields

### Mode properties

- Geometric standard deviation of each lognormal mode 

- Nominal geometric mean diameter and its lower/upper bound of values for each mode 

- Coefficients of polynomial expression (lookup-table) for extinction, absorption, and asymmetry parameter calculations. 

- Aerosol refractive index table needed for interpolation (lookup-table) calculation (for different wavelengths) 

- Crystalization and deliquesence relative humidity thresholds for aerosol wateruptake calculations 

### Species properties

- Aerosol refractive index for each species 

- Density for each species 

- Aerosol hygroscopicity for each species 

- Note that some of variables in the species property file are for bulk aerosols, so we ignore the description for them. 

## Files

```
Species properties

/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/sulfate_rrtmg_c080918.nc
/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/ocpho_rrtmg_c130709_kPOM0.04.nc
/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/ocphi_rrtmg_c100508.nc
/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/bcpho_rrtmg_c100508.nc
/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/dust_aeronet_rrtmg_c141106.nc
/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/ssam_rrtmg_c100508.nc
/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/poly_rrtmg_c130816.nc

Mode properties

/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/mam4_mode1_rrtmg_aeronetdust_c141106.nc', 
/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/mam4_mode2_rrtmg_c130628.nc',
/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/mam4_mode3_rrtmg_aeronetdust_c141106.nc', 
/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/mam4_mode4_rrtmg_c130628.nc',
/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/mam5_mode5_rrtmg_sig1.2_dgnl.40_c03072023.nc'
```

## Namelist

```
 mode_defs		= 'mam5_mode1:accum:=', 
         'A:num_a1:N:num_c1:num_mr:+',
         'A:so4_a1:N:so4_c1:sulfate:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/sulfate_rrtmg_c080918.nc:+', 
         'A:pom_a1:N:pom_c1:p-organic:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/ocpho_rrtmg_c130709_kPOM0.04.nc:+',
         'A:soa_a1:N:soa_c1:s-organic:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/ocphi_rrtmg_c100508.nc:+', 
         'A:bc_a1:N:bc_c1:black-c:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/bcpho_rrtmg_c100508.nc:+',
         'A:dst_a1:N:dst_c1:dust:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/dust_aeronet_rrtmg_c141106.nc:+', 
         'A:ncl_a1:N:ncl_c1:seasalt:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/ssam_rrtmg_c100508.nc:+',
         'A:mom_a1:N:mom_c1:m-organic:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/poly_rrtmg_c130816.nc', 
         'mam5_mode2:aitken:=',
         'A:num_a2:N:num_c2:num_mr:+', 
         'A:so4_a2:N:so4_c2:sulfate:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/sulfate_rrtmg_c080918.nc:+',
         'A:soa_a2:N:soa_c2:s-organic:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/ocphi_rrtmg_c100508.nc:+', 
         'A:ncl_a2:N:ncl_c2:seasalt:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/ssam_rrtmg_c100508.nc:+',
         'A:mom_a2:N:mom_c2:m-organic:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/poly_rrtmg_c130816.nc', 
         'mam5_mode3:coarse:=',
         'A:num_a3:N:num_c3:num_mr:+', 
         'A:dst_a3:N:dst_c3:dust:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/dust_aeronet_rrtmg_c141106.nc:+',
         'A:ncl_a3:N:ncl_c3:seasalt:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/ssam_rrtmg_c100508.nc:+', 
         'A:so4_a3:N:so4_c3:sulfate:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/sulfate_rrtmg_c080918.nc:+',
         'A:bc_a3:N:bc_c3:black-c:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/bcpho_rrtmg_c100508.nc:+', 
         'A:pom_a3:N:pom_c3:p-organic:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/ocpho_rrtmg_c130709_kPOM0.04.nc:+',
         'A:soa_a3:N:soa_c3:s-organic:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/ocphi_rrtmg_c100508.nc:+', 
         'A:mom_a3:N:mom_c3:m-organic:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/poly_rrtmg_c130816.nc',
         'mam5_mode4:primary_carbon:=', 
         'A:num_a4:N:num_c4:num_mr:+',
         'A:pom_a4:N:pom_c4:p-organic:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/ocpho_rrtmg_c130709_kPOM0.04.nc:+', 
         'A:bc_a4:N:bc_c4:black-c:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/bcpho_rrtmg_c100508.nc:+',
         'A:mom_a4:N:mom_c4:m-organic:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/poly_rrtmg_c130816.nc', 
         'mam5_mode5:strat_coarse:=',
         'A:num_a5:N:num_c5:num_mr:+', 
         'A:so4_a5:N:so4_c5:sulfate:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/sulfate_rrtmg_c080918.nc'
 rad_climate = 'A:H2OLNZ:H2O', 'N:O2:O2','N:CO2:CO2', 'A:O3:O3',
         'A:N2OLNZ:N2O', 'A:CH4LNZ:CH4','N:CFC11:CFC11', 'N:CFC12:CFC12',
         'M:mam5_mode1:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/mam4_mode1_rrtmg_aeronetdust_c141106.nc', 
         'M:mam5_mode2:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/mam4_mode2_rrtmg_c130628.nc',
         'M:mam5_mode3:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/mam4_mode3_rrtmg_aeronetdust_c141106.nc', 
         'M:mam5_mode4:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/mam4_mode4_rrtmg_c130628.nc',
         'M:mam5_mode5:/lcrc/group/e3sm/data/inputdata/atm/cam/physprops/mam5_mode5_rrtmg_sig1.2_dgnl.40_c03072023.nc'
         
```