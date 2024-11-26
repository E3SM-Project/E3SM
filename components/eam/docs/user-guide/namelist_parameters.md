
# Namelist parameters associated with atmosphere schemes

## chemUCI and Linoz v3

| Parameter                    | Description                                                              | Default value*         |
| ---------------------------- | ------------------------------------------------------------------------ | ---------------------- |
| `airpl_emis_file`            | Aviation emission                                                        |                        |
| `chlorine_loading_file`      | Chlorine loading                                                         |                        |
| `chlorine_loading_fixed_ymd` |                                                                          |                        |
| `chlorine_loading_type`      |                                                                          |                        |
| `ext_frc_specifier`          | 3-D emissions                                                            |                        |
| `ext_frc_cycle_yr`           |                                                                          |                        |
| `ext_frc_type`               |                                                                          |                        |
| `srf_emis_specifier`         | Surface emissions                                                        |                        |
| `srf_emis_cycle_yr`          |                                                                          |                        |
| `srf_emis_type`              | Upper bound of mean raindrop diameter                                    |                        |
| `linoz_data_file`            | Linoz data file                                                          |                        |
| `linoz_data_cycle_yr`        |                                                                          |                        |
| `linoz_data_path`            |                                                                          |                        |
| `linoz_data_type`            |                                                                          |                        |
| `lght_no_prd_factor`         | Lightning NOx emission factor                                            | `5.0`                  |
| `fstrat_efold_list`          | Tracer (from troposphere) list with e-folding decay in the stratosphere  |                        |

* Many of these namelist parameters specify input data files. Check the `atm_in` file for examples or refer to the [Users' Guide](../user-guide/index.md).

## Cloud Layers Unified By Binormals

| Parameter      | Description                                                                                 | Default value  |
| -------------- | ------------------------------------------------------------------------------------------- | -------------- |
| `gamma_coef`   | Width of vertical velocity within a Gaussian PDF component at low skewness                  | `0.12`         |
| `gamma_coefb`  | Width of vertical velocity within a Gaussian PDF component at high skewness                 | `0.28`         |
| `C8`           | Coefficient of damping of third moment of vertical velocity, w’3                            | `5.2`          |
| `C1`           | Coefficient of damping of second vertical moment of vertical velocity, w’2, at low skewness | `2.4`          |
| `C14`          | Coefficient of damping of second horizontal moments of vertical velocity, u’2 and v’2       | `2.0`          |
| `c_k10`        | Ratio of diffusivity of momentum to heat                                                    | `0.35`         |

## Dust aerosol

| Parameter                 | Description                                    | Default value                                     |
| ------------------------- | ---------------------------------------------- | ------------------------------------------------- |
| `dust_emis_scheme`*        | The v3 dust emission scheme (Kok et al., 2014) | `2` <!-- markdownlint-disable MD033 --><br> (set to 1 to switch to the v1/v2 scheme) |

*This parameter is set in `user_nl_drv`

## HOMME

| Parameter        | Description                                                                                 | Default value  |
| ---------------- | ------------------------------------------------------------------------------------------- | -------------- |
| `se_tstep`       | Main dycore timestep. Additional parameters control the hyper viscsosity, trancer and vertical remap timesteps, which are derived from se_tstep. <!-- markdownlint-disable MD033 --><br> units = seconds | Scales linearly with horizontal resolution. <br> NE30 default: `300` |
| `nu`             | Tensor hyperviscosity coefficient, independent of spatial resolution. <br> units = 1/s      | `3.4e-8` |
| `nu_top`         | Scalar viscosity at model top. <br> units = m^2/s                                           | Horizontal resolution dependent <br> NE30 default: `2.5e5` |
| `transport_alg`  | Select between semi-lagrangian and Eulerian based transport schemes                         | `12` = semi-lagranian method with monotinicity and mass preservation |
| `statefreq`      | print a varieity of dycore metrics to the atm.log file every “statefreq” timesteps          | `480`          |
| `vert_remap_alg` | Algorithm used to remap the vertically lagrangian levels back to the reference levels       | `10` = strict monotonicity applied on top of a 2nd order accurate PPM method  |
| `se_ftype`       | Controls how physics tendencies are applied.  0=”dribbled” in during dynamics timesteps.  1=”hard adjustment” after each physics timestep.  2=hybrid approach: hard adjustment for tracers, dribbled for remaining tendencies | `2`          |

## Modal Aerosol Module

| Parameter                | Description                                                                         | Default value               |
| ------------------------ | ----------------------------------------------------------------------------------- | --------------------------- |
| `is_output_interactive_volc`    | Switch for diagnostic output of the stratospheric aerosol optics | `.false.`  |
| `mam_amicphys_optaa`     | Recommended option of the new time-splitting treatment of H2SO4 production and loss | `1` <!-- markdownlint-disable MD033 --><br> (0 to turn it off) |
| `n_so4_monolayers_pcage` | Number of monolayers required to age primary-carbon mode particles                  | `3`                         |
| `seasalt_emis_scale`     | Tuning parameter for sea salt emission                                              | `0.55`                      |

## OCEANFILMS

| Parameter                 | Description                                                       | Default value          |
| ------------------------- | ----------------------------------------------------------------- | ---------------------- |
| `mam_mom_cycle_yr`       |                                                                    | `1`                    |
| `mam_mom_datapath` | Full pathname of the directory that contains the files specified in mam_mom_filelist  | `'atm/cam/chem/trop_mam/marine_BGC/'`                 |
| `mam_mom_filename`     | Filename of file that contains a sequence of filenames for prescribed marine organic matter ocean concentrations.  The filenames in this file are relative to the directory specified by mam_mom_datapath.| `'monthly_macromolecules_0.1deg_bilinear_latlon_year01_merge_date.nc'` |
| `mam_mom_rmfile`   | Remove the file containing prescribed aerosol deposition fluxes from local disk when no longer needed. | `FALSE`                |
| `mam_mom_specifier`     | Names of variables containing aerosol data in the prescribed aerosol datasets. | `'chla:CHL1','mpoly:TRUEPOLYC','mprot:TRUEPROTC','mlip:TRUELIPC'`                 |
| `mam_mom_datatype`       | Type of time interpolation for data in mam_mom files. Can be set to `'CYCLICAL'`, `'SERIAL'`, `'INTERP_MISSING_MONTHS'`, or `'FIXED'`. | `'CYCLICAL'`                |
| `mam_mom_cycle_yr`         | The  cycle year of the prescribed aerosol flux data if mam_mom_type is `'CYCLICAL'`. Format: YYYY   | `1`               |
| `mam_mom_fixed_ymd`        | The date at which the prescribed aerosol flux data is fixed if mam_mom_type is `'FIXED'`. Format: YYYYMMDD | `0`                |
| `mam_mom_fixed_tod`  | The time of day (seconds) corresponding to mam_mom_fixed_ymd at which the prescribed aerosol flux data is fixed if mam_mom_type is 'FIXED'. | `0`           |
| `mam_mom_bubble_thickness`   | Bubble film thickness (in m) for marine organic aerosol emission mechanism.  The physically reasonable range is approximately (0.1 - 1) x 10^ -6. | `0.1e-6`            |
| `mam_mom_mixing_state`              | Switch to select mixing state assumption in marine organic aerosol code. Currently implemented options: 0 : total external mixture, add to mass; 1 : total external mixture, replace mass; 2 : total internal mixture, add to mass; 3 : total internal mixture, replace mass. | `0` [Note: set to 3 in the atm_in namelist]        |
| `mam_mom_parameterization`     | Selection of alternate parameterizations for marine organic matter emissions.  Set fmoa=1 for Burrows et al. (2014) [@burrows_physically_2014] parameterization; fmoa=2 for Gantt et al. (2011) [@gantt_wind_2011] parameterization; fmoa=3 for simple parameterization based on Quinn et al., 2014; [@quinn_contribution_2014] fmoa=4 for Rinaldi et al. (JGR, 2013).* [@rinaldi_is_2013] | `1`                 |

*Note: non-default values have not been carefully tested and may not work as expected.

## Predicted Particle Properties

| Parameter                 | Description                                                       | Default value          |
| ------------------------- | ----------------------------------------------------------------- | ---------------------- |
| `do_prescribed_ccn`       | Turn on the prescribed CCN if true                                | `false`                |
| `micro_aerosolactivation` | Turn on aerosol activation if true                                | `true`                 |
| `micro_p3_lookup_dir`     | Directory of P3 look-up tables                                    | `inputdata/atm/cam/physprops` |
| `micro_p3_tableversion`   | P3 look-up table Version                                          | `4.1.2`                |
| `micro_subgrid_cloud`     | Sub-grid cloud properties                                         | `true`                 |
| `micro_tend_output`       | Output of P3 microphysical process rates                          | `false`                |
| `p3_accret_coeff`         | Tunable parameter for adjusting rain accretion efficiency         | `117.25`               |
| `p3_autocon_coeff`        | Tunable parameter for adjusting droplet autoconversion efficiency | `30500`                |
| `p3_embryonic_rain_size`  | Radius of embryomic raindrops from auto-conversion                | `0.000025` (m)         |
| `p3_max_mean_rain_size`   | Upper bound of mean raindrop diameter                             | `0.005` (m)            |
| `p3_mincdnc`              | Lower bound of droplet number concentration                       | `20.d6` (# m-3)        |
| `p3_nc_autocon_expon`     | Nc exponent in droplet auto-conversion                            | `-1.1`                 |
| `p3_qc_accret_expon`      | Qc exponent in rain accretion                                     | `1.15`                 |
| `p3_qc_autocon_expon`     | Qc exponeent in droplet autoconversion                            | `3.19`                 |
| `p3_wbf_coeff`            | Tunable parameter for adjusting WBF efficiency                    | `1.0`                  |
| `do_cooper_inp3`          | Turn on Cooper ice nucleation scheme if true                      | `false`                |

## Rapid Radiative Transfer Model for GCMs

| Parameter                 | Description                                                       | Default value          |
| ------------------------- | ----------------------------------------------------------------- | ---------------------- |
| `iradsw`                  | Frequency for updating shortwave fluxes and heating rate; iradsw > 0 interpreted as number of timesteps, iradsw < 0 interpreted as hours; iradsw = 0 disables shortwave radiation entirely | `-1`                 |
| `iradlw`                  | Frequency for updating longwave fluxes and heating rate; iradlw > 0 interpreted as number of timesteps, iradlw < 0 interpreted as hours; iradlw = 0 disables longwave radiation entirely   | `-1`                 |
| `irad_always`             | Length of time in timesteps (irad_always > 0) or in hours (irad_always < 0) SW/LW radiation will be run continuously from the start of an initial or restart run                           | `0`                  |
| `use_rad_dt_cosz`         | If true, use the radiation dt for all cosz calculations; calculates solar zenith angle averaged over a time step. In default model solar zenith angle is held constant over time           | `.true.`  <!-- markdownlint-disable MD033 --><br> (set by namelist_defaults_eam.xml for default physics)   |
| `spectralflux`            | Calculate fluxes (up and down) per band                                                                                                                                                    | `.false.`            |
| `liqcldoptics`            | Choice of cloud optical property parameterization for liquid clouds. Valid options are ‘slingo’ or ‘gammadist’                                                                             | `gammadist`          |
| `icecldoptics`            | Choice of cloud optical property parameterization for ice clouds. Valid options are ‘ebertcurry’ or ‘mitchell’                                                                             | `mitchell`           |

## Zhang and McFarlane deep convection scheme

| ZM Parameters             | Description                                                       | Default value          |
| ------------------------- | ----------------------------------------------------------------- | ---------------------- |
| `zmconv_ke`               | Tunable evaporation efficiency in ZM deep convection scheme       | `2.5E-6`               |
| `zmconv_tau`              | Relaxation time in ZM deep convection scheme                      | `3600`                 |
| `zmconv_dmpdz`            | Parcel fractional mass entrainment rate                           | `-0.7E-3`              |
| `zmconv_alfa`             | Initial downdraft mass flux fraction                              | `0.14D0`               |
| `zmconv_tiedke_add`       | Temperature perturbation of an air parcel                         | `0.8D0`                |
| `zmconv_cape_cin`         | Number of negative buoyancy regions that are allowed              | `1`                    |

| dCAPE-ULL Parameters      | Description                                                       | Default value          |
| ------------------------- | ----------------------------------------------------------------- | ---------------------- |
| `zmconv_trigdcape_ull`    | DCAPE trigger along with unrestricted launching level for ZM deep convection scheme        | `.true.`  |
| `zmconv_trig_dcape_only`  | DCAPE only trigger for ZM deep convection scheme                  | `.false.` <!-- markdownlint-disable MD033 --><br> If true, zmconv_trigdcape_ull must be false to use the dcape only trigger. |
| `zmconv_trig_ull_only`    | Use unrestricted launching level (ULL) only trigger for ZM deep convection scheme          | `.false.` <!-- markdownlint-disable MD033 --><br> If true, zmconv_trigdcape_ull must be false to use the ull only trigger. |

| Conv. micro. Parameters   | Description                                                       | Default value          |
| ------------------------- | ----------------------------------------------------------------- | ---------------------- |
| `zmconv_microp`           | Convective microphysics option in ZM convection scheme            | `true`                 |
| `zmconv_auto_fac`         | Cloud droplet-rain autoconversion enhancement factor in the convective microphysics scheme | `7.0`     |
| `zmconv_accr_fac`         | Cloud droplet-rain accretion enhancement factor in the convective microphysics scheme      | `1.5`     |
| `zmconv_micro_dcs`        | Autoconversion size threshold for cloud ice to snow (m)           | `150.E-6`                 |

| Mass flux adj. Parameters | Description                                                       | Default value          |
| ------------------------- | ----------------------------------------------------------------- | ---------------------- |
| `zmconv_clos_dyn_adj`     | Apply mass flux adjustment to ZM convection scheme               | `true`                 |

| MCSP Parameters              | Description                                                       | Default value          |
| ---------------------------- | ----------------------------------------------------------------- | ---------------------- |
| `zmconv_mcsp_heat_coeff`     | MCSP heating coefficient                                          | `0.3`                  |
| `zmconv_mcsp_moisture_coeff` | MCSP moisture coefficient                                         | `0.0`                  |
| `zmconv_mcsp_uwind_coeff`    | MCSP zonal wind coefficient                                       | `0.0`                  |
| `zmconv_mcsp_vwind_coeff`    | MCSP meridional wind coefficient                                  | `0.0`                  |

## Cloud Feedback Model Intercomparison Project (CFMIP) Observation Simulator Package

| Parameter                 | Description                                                       | Default value          |
| ------------------------- | ----------------------------------------------------------------- | ---------------------- |
| `cosp_lite`       | This namelist sets cosp_ncolumns=10 and cosp_nradsteps=3 (appropriate for COSP statistics derived from seasonal averages), and runs MISR, ISCCP, MODIS, and CALIPSO lidar simulators (cosp_lmisr_sim=.true.,cosp_lisccp_sim=.true., cosp_lmodis_sim=.true.,cosp_llidar_sim=.true.).  | `false`                |

## Orographic drag schemes

| Parameter                 | Description                                                       | Default value          |
| ------------------------- | ----------------------------------------------------------------- | ---------------------- |
| `use_gw_oro`       | This namelist controls the default linear orographic gravity wave drag (oGWD) for E3SM, if used, the default oGWD is turned on.                                | `true`                |
| `do_tms`           | This namelist controls the default Turbulent Mountain Stress (TMS) for E3SM, if used, the default TMS is turned on.                                 | `false`                |
| `effgw_oro`        | Efficiency associated with orographic gravity waves.                                                                     | `0.375`                |
| `tms_orocnst`      | Turbulent mountain stress parameter used when TMS calculation is turned on             | `1.0`                |
| `tms_z0fac`        | Factor determining z_0 from orographic standard deviation [ no unit ] for TMS.                                                   | `0.75`                |
| `use_od_ls`        | This namelist controls the new nonlinear oGWD, if used, the nonlinear oGWD is turned on. use_od_ls should not be used at the same time with use_gw_oro.                  | `false`                |
| `use_od_bl`       | This namelist controls the Flow-blocking drag (FBD) scheme, if used, the FBD scheme is turned on.        | `false`                |
| `use_od_ss`       | This namelist controls the small-scale GWD (sGWD) scheme, if used, the sGWD scheme is turned on.        | `false`                |
| `use_od_fd`       | This namelist controls the Turbulent orographic form drag (TOFD) scheme, if used, the TOFD scheme is turned on.        | `false`                |
| `od_ls_ncleff`    | Tuning parameter of nonlinear oGWD. Stands for effective resolution of the grid for oGWD. Scales the magnitude of nonlinear oGWD.         | `3`                |
| `od_bl_ncd`       | Tuning parameter of FBD. Stands for bulk drag coefficient. Scales the magnitude of FBD.       | `3`                |
| `od_ss_sncleff`   | Tuning parameter of sGWD. Stands for effective resolution of the grid for sGWD.Scales the magnitude of sGWD.        | `1`                |
