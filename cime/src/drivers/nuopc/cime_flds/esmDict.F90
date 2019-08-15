module esmDict

  ! Establish the NUOPC dictionary for field names
  ! The call to the dictionary initialization needs to be done on all PETS

  implicit none
  public

  public :: esmDict_Init

  character(*), parameter :: u_FILE_u = &
       __FILE__

!================================================================================
contains
!================================================================================

  subroutine esmDict_Init(rc)

    use ESMF                  , only : ESMF_SUCCESS
    use med_constants_mod     , only : CS
    use glc_elevclass_mod     , only : glc_elevclass_as_string
    use shr_nuopc_scalars_mod , only : flds_scalar_name
    use shr_nuopc_fldlist_mod , only : shr_nuopc_fldList_AddMetadata

    ! input/output variables:
    integer, intent(inout) :: rc

    ! local variables:
    integer                     :: ice_ncat                   ! number of sea ice thickness categories
    integer                     :: glc_nec                    ! number of land-ice elevation classes
    integer                     :: max_megan
    integer                     :: max_ddep   
    integer                     :: max_fire
    logical                     :: flds_i2o_per_cat
    integer                     :: n, num
    character(len= 2)           :: cnum
    character(len=CS)           :: units
    character(len=CS)           :: longname
    character(len=CS)           :: stdname
    character(len=CS)           :: name, fldname
    character(len=*), parameter :: subname='(esmDict_Init)'
    !--------------------------------------

    rc = ESMF_SUCCESS

    !---------------------------
    ! For now hardwire these
    !---------------------------

    max_megan = 20
    max_ddep  = 80
    max_fire  = 10
    glc_nec   = 10
    ice_ncat  =  5 
    flds_i2o_per_cat = .true.

    !---------------------------
    ! Create dictionary names
    !---------------------------

    longname = trim(flds_scalar_name)
    stdname  = trim(flds_scalar_name)
    units    = 'unitless'
    call shr_nuopc_fldList_AddMetadata(trim(flds_scalar_name), longname, stdname, units)

    longname = 'latitude'
    stdname  = 'latitude'
    units    = 'degrees north'
    call shr_nuopc_fldList_AddMetadata('lat', longname, stdname, units)

    longname = 'longitude'
    stdname  = 'longitude'
    units    = 'degrees east'
    call shr_nuopc_fldList_AddMetadata('lon', longname, stdname, units)

    longname = 'height'
    stdname  = 'height, depth, or levels'
    units    = 'unitless'
    call shr_nuopc_fldList_AddMetadata('hgt', longname, stdname, units)

    longname = 'cell_area_model'
    stdname  = 'cell area from model'
    units    = 'radian^2'
    call shr_nuopc_fldList_AddMetadata('area', longname, stdname, units)

    longname = 'cell_area_mapping'
    stdname  = 'cell area from mapping file'
    units    = 'radian^2'
    call shr_nuopc_fldList_AddMetadata('aream', longname, stdname, units)

    longname = 'mask'
    stdname  = 'mask'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata('mask', longname, stdname, units)

    longname = 'area_fraction'
    stdname  = 'area fraction'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata('frac', longname, stdname, units)

    !----------------------------------------------------------
    ! Masks from components
    !----------------------------------------------------------

    longname = 'Surface fraction in land'
    stdname  = 'land_fraction_from_land'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata("Sl_lfrin", longname, stdname, units)

    longname = 'Sea surface mask'
    stdname  = 'sea_surface_mask'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata("So_omask", longname, stdname, units)

    longname = 'Sea ice mask'
    stdname  = 'sea_ice_mask'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata("Si_imask", longname, stdname, units)

    !----------------------------------------------------------
    ! Fractions sent to atm
    !----------------------------------------------------------

    longname = 'Surface land fraction'
    stdname  = 'land_area_fraction'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata("Sl_lfrac", longname, stdname, units)

    longname = 'Surface ice fraction'
    stdname  = 'sea_ice_area_fraction'
    call shr_nuopc_fldList_AddMetadata("Si_ifrac", longname, stdname, units)

    longname = 'Surface ocean fraction'
    stdname  = 'sea_area_fraction'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata("So_ofrac", longname, stdname, units)

    !----------------------------------------------------------
    ! Fractional ice coverage wrt ocean sent to ocn and wav
    !----------------------------------------------------------

    longname = 'Fractional ice coverage wrt ocean'
    stdname  = 'sea_ice_area_fraction'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata("Si_ifrac", longname, stdname, units)

    !----------------------------------------------------------
    ! Fields from atm
    !----------------------------------------------------------

    longname = 'Height at the lowest model level'
    stdname  = 'height'
    units    = 'm'
    call shr_nuopc_fldList_AddMetadata('Sa_z', longname, stdname, units)

    longname = 'Surface height'
    stdname  = 'height'
    units    = 'm'
    call shr_nuopc_fldList_AddMetadata('Sa_topo', longname, stdname, units)

    longname = 'Zonal wind at the lowest model level'
    stdname  = 'eastward_wind'
    units    = 'm s-1'
    call shr_nuopc_fldList_AddMetadata('Sa_u', longname, stdname, units)

    longname = 'Meridional wind at the lowest model level'
    stdname  = 'northward_wind'
    units    = 'm s-1'
    call shr_nuopc_fldList_AddMetadata('Sa_v', longname, stdname, units)

    longname = 'Temperature at the lowest model level'
    stdname  = 'air_temperature'
    units    = 'K'
    call shr_nuopc_fldList_AddMetadata('Sa_tbot', longname, stdname, units)

    longname = 'Potential temperature at the lowest model level'
    stdname  = 'air_potential_temperature'
    units    = 'K'
    call shr_nuopc_fldList_AddMetadata('Sa_ptem', longname, stdname, units)

    longname = 'Specific humidity at the lowest model level'
    stdname  = 'specific_humidity'
    units    = 'kg kg-1'
    call shr_nuopc_fldList_AddMetadata('Sa_shum', longname, stdname, units)

    longname = 'Pressure at the lowest model level'
    stdname  = 'air_pressure'
    units    = 'Pa'
    call shr_nuopc_fldList_AddMetadata('Sa_pbot', longname, stdname, units)

    longname = 'Density at the lowest model level'
    stdname  = 'air_density'
    units    = 'kg m-3'
    call shr_nuopc_fldList_AddMetadata('Sa_dens', longname, stdname, units)

    longname = 'Convective precipitation rate'
    stdname  = 'convective_precipitation_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faxa_rainc', longname, stdname, units)

    longname = 'Large-scale (stable) precipitation rate' ! water equivalent
    stdname  = 'large_scale_precipitation_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faxa_rainl', longname, stdname, units)

    longname = 'Water flux due to rain'
    stdname  = 'rainfall_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faxa_rain', longname, stdname, units)

    longname = 'Convective snow rate'
    stdname  = 'convective_snowfall_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faxa_snowc', longname, stdname, units)

    longname = 'Large-scale (stable) snow rate'
    stdname  = 'large_scale_snowfall_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faxa_snowl', longname, stdname, units)

    longname = 'Water flux due to snow'
    stdname  = 'surface_snow_melt_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faxa_snow', longname, stdname, units)

    ! total precipitation to ocean (derived rain + snow, done AFTER mappings)
    longname = 'Water flux (rain+snow)'
    stdname  = 'precipitation_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faxa_prec', longname, stdname, units)

    longname = 'Downward longwave heat flux'
    stdname  = 'downwelling_longwave_flux'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata('Faxa_lwdn', longname, stdname, units)

    longname = 'Direct near-infrared incident solar radiation'
    stdname  = 'surface_downward_direct_shortwave_flux_due_to_near_infrared_radiation'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata("Faxa_swndr", longname, stdname, units)

    longname = 'Direct visible incident solar radiation'
    stdname  = 'surface_downward_direct_shortwave_flux_due_to_visible_radiation'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata("Faxa_swvdr", longname, stdname, units)

    longname = 'Diffuse near-infrared incident solar radiation'
    stdname  = 'surface_downward_diffuse_shortwave_flux_due_to_near_infrared_radiation'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata("Faxa_swndf", longname, stdname, units)

    longname = 'Diffuse visible incident solar radiation'
    stdname  = 'surface_downward_diffuse_shortwave_flux_due_to_visible_radiation'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata('Faxa_swvdf', longname, stdname, units)

    longname = 'Net shortwave radiation'
    stdname  = 'surface_net_shortwave_flux'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata("Faxa_swnet", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Fall_swnet", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Faii_swnet", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Foxx_swnet", longname, stdname, units)

    longname = 'Net shortwave radiation penetrating into ice and ocean'
    stdname  = 'net_downward_shortwave_flux_in_sea_ice_due_to_penetration'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata('Fioi_swpen', longname, stdname, units)

    longname ='Hydrophylic black carbon dry deposition flux'
    stdname  = 'dry_deposition_flux_of_hydrophylic_black_carbon'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faxa_bcphidry', longname, stdname, units)

    longname = 'Hydrophobic black carbon dry deposition flux'
    stdname  = 'dry_deposition_flux_of_hydrophobic_black_carbon'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Faxa_bcphodry", longname, stdname, units)

    longname = 'Hydrophylic black carbon wet deposition flux'
    stdname  = 'wet_deposition_flux_of_hydrophylic_black_carbon'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Faxa_bcphiwet", longname, stdname, units)

    longname = 'Hydrophylic organic carbon dry deposition flux'
    stdname  = 'dry_deposition_flux_of_hydrophylic_organic_carbon'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Faxa_ocphidry", longname, stdname, units)

    longname = 'Hydrophobic organic carbon dry deposition flux'
    stdname  = 'dry_deposition_flux_of_hydrophobic_organic_carbon'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Faxa_ocphodry", longname, stdname, units)

    longname = 'Hydrophylic organic carbon wet deposition flux'
    stdname  = 'wet_deposition_flux_of_hydrophylic_organic_carbon'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Faxa_ocphiwet", longname, stdname, units)

    longname = 'Dust wet deposition flux (size 1)'
    stdname  = 'wet_deposition_flux_of_dust'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Faxa_dstwet1", longname, stdname, units)

    longname = 'Dust wet deposition flux (size 2)'
    stdname  = 'wet_deposition_flux_of_dust'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Faxa_dstwet2", longname, stdname, units)

    longname = 'Dust wet deposition flux (size 3)'
    stdname  = 'wet_deposition_flux_of_dust'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Faxa_dstwet3", longname, stdname, units)

    longname = 'Dust wet deposition flux (size 4)'
    stdname  = 'wet_deposition_flux_of_dust'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Faxa_dstwet4", longname, stdname, units)

    longname = 'Dust dry deposition flux (size 1)'
    stdname  = 'dry_deposition_flux_of_dust'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Faxa_dstdry1", longname, stdname, units)

    longname = 'Dust dry deposition flux (size 2)'
    stdname  = 'dry_deposition_flux_of_dust'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Faxa_dstdry2", longname, stdname, units)

    longname = 'Dust dry deposition flux (size 3)'
    stdname  = 'dry_deposition_flux_of_dust'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Faxa_dstdry3", longname, stdname, units)

    longname = 'Dust dry deposition flux (size 4)'
    stdname  = 'dry_deposition_flux_of_dust'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Faxa_dstdry4", longname, stdname, units)

    !----------------------------------------------------------
    ! states/fluxes to atm (and ocean)
    !----------------------------------------------------------

    longname = 'Direct albedo (visible radiation)'
    stdname  = 'surface_direct_albedo_due_to_visible_radiation'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata("Si_avsdr", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Sl_avsdr", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("So_avsdr", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Sx_avsdr", longname, stdname, units)

    longname = 'Direct albedo (near-infrared radiation)'
    stdname  = 'surface_direct_albedo_due_to_near_infrared_radiation'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata("Si_anidr", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Sl_anidr", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("So_anidr", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Sx_anidr", longname, stdname, units)

    longname = 'Diffuse albedo (visible radiation)'
    stdname  = 'surface_diffuse_albedo_due_to_visible_radiation'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata("Si_avsdf", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Sl_avsdf", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("So_avsdf", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Sx_avsdf", longname, stdname, units)

    longname = 'Diffuse albedo (near-infrared radiation)'
    stdname  = 'surface_diffuse_albedo_due_to_near_infrared_radiation'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata("Si_anidf", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Sl_anidf", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("So_anidf", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Sx_anidf", longname, stdname, units)

    longname = 'Reference temperature at 2 meters'
    stdname  = 'air_temperature'
    units    = 'K'
    call shr_nuopc_fldList_AddMetadata("Si_tref", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Sl_tref", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("So_tref", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Sx_tref", longname, stdname, units)

    longname = 'Reference specific humidity at 2 meters'
    stdname  = 'specific_humidity'
    units    = 'kg kg-1'
    call shr_nuopc_fldList_AddMetadata("Si_qref", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Sl_qref", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("So_qref", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Sx_qref", longname, stdname, units)

    longname = 'Surface temperature'
    stdname  = 'surface_temperature'
    units    = 'K'
    call shr_nuopc_fldList_AddMetadata("Si_t", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Sl_t", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("So_t", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Sx_t", longname, stdname, units)

    longname = 'Surface fraction velocity in land'
    stdname  = 'fraction_velocity'
    units    = 'm s-1'
    call shr_nuopc_fldList_AddMetadata("Sl_fv", longname, stdname, units)

    longname = 'Aerodynamic resistance'
    stdname  = 'aerodynamic_resistance'
    units    = 's/m'
    call shr_nuopc_fldList_AddMetadata("Sl_ram1", longname, stdname, units)

    longname = 'Surface snow water equivalent'
    stdname  = 'surface_snow_water_equivalent'
    units    = 'm'
    call shr_nuopc_fldList_AddMetadata("Sl_snowh", longname, stdname, units)

    longname = 'Surface snow depth'
    stdname  = 'surface_snow_thickness'
    units    = 'm'
    call shr_nuopc_fldList_AddMetadata("Si_snowh", longname, stdname, units)

    longname = 'Surface saturation specific humidity in ocean'
    stdname  = 'specific_humidity_at_saturation'
    units    = 'kg kg-1'
    call shr_nuopc_fldList_AddMetadata("So_ssq", longname, stdname, units)

    longname = 'Square of exch. coeff (tracers)'
    stdname  = 'square_of_exch_coeff'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata("So_re", longname, stdname, units)

    longname = '10m wind'
    stdname  = '10m_wind'
    units    = 'm'
    call shr_nuopc_fldList_AddMetadata("Sl_u10", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Si_u10", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("So_u10", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Sx_u10", longname, stdname, units)

    longname = 'Zonal surface stress'
    stdname  = 'surface_downward_eastward_stress'
    units    = 'N m-2'
    call shr_nuopc_fldList_AddMetadata("Fall_taux", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Faox_taux", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Faii_taux", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Fioi_taux", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Faxx_taux", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Foxx_taux", longname, stdname, units)

    longname = 'Meridional surface stress'
    stdname  = 'surface_downward_northward_stress'
    units    = 'N m-2'
    call shr_nuopc_fldList_AddMetadata("Fall_tauy", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Faox_tauy", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Faii_tauy", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Fioi_tauy", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Faxx_tauy", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Foxx_tauy", longname, stdname, units)

    longname = 'Surface latent heat flux'
    stdname  = 'surface_upward_latent_heat_flux'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata("Fall_lat", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Faox_lat", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Faii_lat", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Faxx_lat", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Foxx_lat", longname, stdname, units)

    longname = 'Sensible heat flux'
    stdname  = 'surface_upward_sensible_heat_flux'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata("Fall_sen", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Faox_sen", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Faii_sen", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Faxx_sen", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Foxx_sen", longname, stdname, units)

    longname = 'Surface upward longwave heat flux'
    stdname  = 'surface_net_upward_longwave_flux'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata("Fall_lwup", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Faox_lwup", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Faii_lwup", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Faxx_lwup", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Foxx_lwup", longname, stdname, units)

    longname = 'Evaporation water flux'
    stdname  = 'water_evaporation_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Fall_evap", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Faox_evap", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Faii_evap", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Faxx_evap", longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata("Foxx_evap", longname, stdname, units)

    longname = 'Dust flux (particle bin number 1)'
    stdname  = 'dust_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Fall_flxdst1", longname, stdname, units)

    longname = 'Dust flux (particle bin number 2)'
    stdname  = 'dust_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Fall_flxdst2", longname, stdname, units)

    longname = 'Dust flux (particle bin number 3)'
    stdname  = 'dust_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Fall_flxdst3", longname, stdname, units)

    longname = 'Dust flux (particle bin number 4)'
    stdname  = 'dust_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Fall_flxdst4", longname, stdname, units)

    !-----------------------------
    ! atm<->ocn only exchange
    !-----------------------------

    longname = 'Sea level pressure'
    stdname  = 'air_pressure_at_sea_level'
    units    = 'Pa'
    call shr_nuopc_fldList_AddMetadata("Sa_pslv", longname, stdname, units)

    longname = 'Wind speed squared at 10 meters'
    stdname  = 'square_of_wind_speed'
    units    = 'm2 s-2'
    call shr_nuopc_fldList_AddMetadata("So_duu10n", longname, stdname, units)

    longname = 'Surface fraction velocity in ocean'
    stdname  = 'fraction_velocity'
    units    = 'm s-1'
    call shr_nuopc_fldList_AddMetadata("So_ustar", longname, stdname, units)

    !-----------------------------
    ! ice->ocn exchange
    !-----------------------------

    longname = 'Heat flux from melting'
    stdname  = 'surface_snow_melt_heat_flux'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata("Fioi_melth", longname, stdname, units)

    longname = 'Water flux due to melting'
    stdname  = 'surface_snow_melt_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Fioi_meltw", longname, stdname, units)

    longname = 'Salt flux'
    stdname  = 'virtual_salt_flux_into_sea_water'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Fioi_salt", longname, stdname, units)

    longname = 'Hydrophylic black carbon deposition flux'
    stdname  = 'deposition_flux_of_hydrophylic_black_carbon'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Fioi_bcphi", longname, stdname, units)

    longname = 'Hydrophobic black carbon deposition flux'
    stdname  = 'deposition_flux_of_hydrophobic_black_carbon'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Fioi_bcpho", longname, stdname, units)

    longname = 'Dust flux'
    stdname  = 'dust_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Fioi_flxdst", longname, stdname, units)

    !-----------------------------
    ! ocn -> ice exchange (some of these fields are also used in the atm/ocn flux computation)
    !-----------------------------

    longname = 'Sea surface salinity'
    stdname  = 'sea_surface_salinity'
    units    = 'g kg-1'
    call shr_nuopc_fldList_AddMetadata("So_s", longname, stdname, units)

    longname = 'Zonal sea water velocity'
    stdname  = 'eastward_sea_water_velocity'
    units    = 'm s-1'
    call shr_nuopc_fldList_AddMetadata("So_u", longname, stdname, units)

    longname = 'Fraction of sw penetrating surface layer for diurnal cycle'
    stdname  = 'Fraction_of_sw_penetrating_surface_layer'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata("So_fswpen", longname, stdname, units)

    longname = 'Meridional sea water velocity'
    stdname  = 'northward_sea_water_velocity'
    units    = 'm s-1'
    call shr_nuopc_fldList_AddMetadata("So_v", longname, stdname, units)

    longname = 'Zonal sea surface slope'
    stdname  = 'sea_surface_eastward_slope'
    units    = 'm m-1'
    call shr_nuopc_fldList_AddMetadata("So_dhdx", longname, stdname, units)

    longname = 'Meridional sea surface slope'
    stdname  = 'sea_surface_northward_slope'
    units    = 'm m-1'
    call shr_nuopc_fldList_AddMetadata("So_dhdy", longname, stdname, units)

    longname = 'Ocean Boundary Layer Depth'
    stdname  = 'ocean_boundary_layer_depth'
    units    = 'm'
    call shr_nuopc_fldList_AddMetadata("So_bldepth", longname, stdname, units)

    longname = 'Ocean melt and freeze potential'
    stdname  = 'surface_snow_and_ice_melt_heat_flux'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata("Fioo_q", longname, stdname, units)

    !-----------------------------
    ! lnd->rof exchange
    !-----------------------------

    longname = 'Water flux from land (liquid surface)'
    stdname  = 'water_flux_into_runoff_surface'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Flrl_rofsur", longname, stdname, units)

    longname = 'Water flux from land (liquid glacier, wetland, and lake)'
    stdname  = 'water_flux_into_runoff_from_gwl'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Flrl_rofgwl", longname, stdname, units)

    longname = 'Water flux from land (liquid subsurface)'
    stdname  = 'water_flux_into_runoff_subsurface'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Flrl_rofsub", longname, stdname, units)

    longname = 'Water flux from land direct to ocean'
    stdname  = 'water_flux_direct_to_ocean'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Flrl_rofdto", longname, stdname, units)

    longname = 'Water flux from land (frozen)'
    stdname  = 'frozen_water_flux_into_runoff'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Flrl_rofi", longname, stdname, units)

    ! Irrigation flux (land/rof only)
    longname = 'Irrigation flux (withdrawal from rivers)'
    stdname  = 'irrigation'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Flrl_irrig", longname, stdname, units)

    !-----------------------------
    ! rof->lnd
    !-----------------------------

    longname = 'Waterflux back to land due to flooding'
    stdname  = 'flooding_water_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Flrr_flood", longname, stdname, units)

    longname = 'River channel total water volume'
    stdname  = 'rtm_volr'
    units    = 'm'
    call shr_nuopc_fldList_AddMetadata("Flrr_volr", longname, stdname, units)

    longname = 'River channel main channel water volume'
    stdname  = 'rtm_volrmch'
    units    = 'm'
    call shr_nuopc_fldList_AddMetadata("Flrr_volrmch", longname, stdname, units)

    !-----------------------------
    ! rof->ocn (liquid and frozen) and glc->ocn
    !-----------------------------

    longname = 'glc liquid runoff flux to ocean'
    stdname  = 'glacier_liquid_runoff_flux_to_ocean'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Fogg_rofl', longname, stdname, units)

    longname = 'Water flux into sea water due to runoff (liquid)'
    stdname  = 'water_flux_into_sea_water'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Forr_rofl", longname, stdname, units)

    longname = 'Total Water flux into sea water due to runoff (liquid)'
    stdname  = 'total_water_flux_into_sea_water'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Foxx_rofl", longname, stdname, units)

    longname = 'glc frozen runoff flux to ocean'
    stdname  = 'glacier_frozen_runoff_flux_to_ocean'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Fogg_rofi', longname, stdname, units)

    longname = 'Water flux into sea water due to runoff (frozen)'
    stdname  = 'frozen_water_flux_into_sea_water'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Forr_rofi", longname, stdname, units)

    longname = 'Total Water flux into sea water due to runoff (frozen)'
    stdname  = 'total_frozen_water_flux_into_sea_water'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Foxx_rofi", longname, stdname, units)

    !-----------------------------
    ! rof(frozen)->ice and glc->ice
    !-----------------------------

    longname = 'Water flux into sea ice due to runoff (frozen)'
    stdname  = 'frozen_water_flux_into_sea_ice'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Firr_rofi", longname, stdname, units)

    longname = 'glc frozen runoff_iceberg flux to ice'
    stdname  = 'glacier_frozen_runoff_flux_to_seaice'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Figg_rofi', longname, stdname, units)

    longname = 'Total frozen water flux into sea ice '
    stdname  = 'total_frozen_water_flux_into_sea_ice'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata("Fixx_rofi", longname, stdname, units)

    !-----------------------------
    ! wav->ocn
    !-----------------------------

    longname = 'Langmuir multiplier'
    stdname  = 'wave_model_langmuir_multiplier'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata('Sw_lamult', longname, stdname, units)

    longname = 'Stokes drift u component'
    stdname  = 'wave_model_stokes_drift_eastward_velocity'
    units    = 'm/s'
    call shr_nuopc_fldList_AddMetadata('Sw_ustokes', longname, stdname, units)

    longname = 'Stokes drift v component'
    stdname  = 'wave_model_stokes_drift_northward_velocity'
    units    = 'm/s'
    call shr_nuopc_fldList_AddMetadata('Sw_vstokes', longname, stdname, units)

    longname = 'Stokes drift depth'
    stdname  = 'wave_model_stokes_drift_depth'
    units    = 'm'
    call shr_nuopc_fldList_AddMetadata('Sw_hstokes', longname, stdname, units)

    longname = 'Downward solar radiation'
    stdname  = 'surface_downward_shortwave_flux'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata("Faox_swdn", longname, stdname, units)

    longname = 'Upward solar radiation'
    stdname  = 'surface_upward_shortwave_flux'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata("Faox_swup", longname, stdname, units)

    !-----------------------------
    ! glc -> ocn
    !-----------------------------

    !-----------------------------
    ! glc -> lnd
    !-----------------------------

    longname = 'Ice sheet grid coverage on global grid'
    stdname  = 'ice_sheet_grid_mask'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata("Sg_icemask", longname, stdname, units)

    longname = 'Ice sheet mask where we are potentially sending non-zero fluxes'
    stdname  = 'icemask_coupled'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata("Sg_icemask_coupled_fluxes", longname, stdname, units)

    longname = 'Fraction of glacier area'
    stdname  = 'glacier_area_fraction'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata('Sg_ice_covered', longname, stdname, units)
    if (glc_nec > 0) then
       name = 'Sg_ice_covered'
       do num = 0, glc_nec
          cnum = glc_elevclass_as_string(num)
          call shr_nuopc_fldList_AddMetadata( 'Sg_ice_covered'//trim(cnum), &
               trim(longname)//' of elevation class '//trim(cnum), stdname , units)
       end do
    end if

    longname = 'Surface height of glacier'
    stdname  = 'height'
    units    = 'm'
    call shr_nuopc_fldList_AddMetadata('Sg_topo', longname, stdname, units)
    if (glc_nec > 0) then
       name = 'Sg_topo'
       do num = 0, glc_nec
          cnum = glc_elevclass_as_string(num)
          call shr_nuopc_fldList_AddMetadata( 'Sg_topo'//trim(cnum), &
               trim(longname)//' of elevation class '//trim(cnum), stdname , units)
       end do
    end if

    longname = 'Downward heat flux from glacier interior'
    stdname  = 'downward_heat_flux_in_glacier'
    units    = 'W m-2'
    call shr_nuopc_fldList_AddMetadata('Flgg_hflx', longname, stdname, units)
    if (glc_nec > 0) then
       name = 'Flgg_hflx'
       do num = 0, glc_nec
          cnum = glc_elevclass_as_string(num)
          call shr_nuopc_fldList_AddMetadata( 'Flgg_hflx'//trim(cnum), &
               trim(longname)//' of elevation class '//trim(cnum), stdname, units)
       end do
    end if

    !-----------------------------
    ! lnd -> glc
    !-----------------------------

    ! glc fields with multiple elevation classes: lnd->glc
    ! - fields sent from lnd->med are in multiple elevation classes
    ! - fields sent from med->glc do NOT have elevation classes
    ! - need to keep track of the l2x fields destined for glc in the
    !   additional variables, l2x_to_glc. This is needed so that can set up
    !   additional fields holding accumulated quantities of just these fields.

    ! Sets a coupling field for all glc elevation classes (1:glc_nec) plus bare land (index 0).
    ! Note that, if glc_nec = 0, then we don't create any coupling fields (not even the bare land (0) fldindex)

    longname = 'New glacier ice flux'
    stdname  = 'ice_flux_out_of_glacier'
    units    = 'kg m-2 s-1'
    if (glc_nec > 0) then
       do num = 0, glc_nec
          cnum = glc_elevclass_as_string(num)
          call shr_nuopc_fldList_AddMetadata('Flgl_qice'//trim(cnum), &
               trim(longname)//' of elevation class '//trim(cnum), stdname, units)
       end do
    end if
    call shr_nuopc_fldList_AddMetadata( 'Flgl_qice', longname, stdname, units)

    longname = 'Surface temperature of glacier'
    stdname  = 'surface_temperature'
    units    = 'deg C'
    if (glc_nec > 0) then
       do num = 0, glc_nec
          cnum = glc_elevclass_as_string(num)
          call shr_nuopc_fldList_AddMetadata('Sl_tsrf'//trim(cnum), &
               trim(longname)//' of elevation class '//trim(cnum), stdname, units)
       end do
    end if
    call shr_nuopc_fldList_AddMetadata( 'Sl_tsrf', longname, stdname, units)

    ! Sl_topo is sent from lnd -> med, but is NOT sent to glc (it is only used for the
    ! remapping in the mediator)

    longname = 'Surface height'
    stdname  = 'height'
    units    = 'm'
    if (glc_nec > 0) then
       do num = 0, glc_nec
          cnum = glc_elevclass_as_string(num)
          call shr_nuopc_fldList_AddMetadata('Sl_topo'//trim(cnum), &
               trim(longname)//' of elevation class '//trim(cnum), stdname, units)
       end do
    end if
    call shr_nuopc_fldList_AddMetadata( 'Sl_topo', longname, stdname, units)

    longname = 'Surface flux of CO2 from land'
    stdname  = 'surface_upward_flux_of_carbon_dioxide_where_land'
    units    = 'moles m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Fall_fco2_lnd', longname, stdname, units)

    longname = 'Prognostic CO2 at the lowest model level'
    stdname  = 'prognostic_CO2_lowest_level'
    units    = '1e-6 mol/mol'
    call shr_nuopc_fldList_AddMetadata('Sa_co2prog', longname, stdname, units)

    longname = 'Diagnostic CO2 at the lowest model level'
    stdname  = 'diagnostic_CO2_lowest_level'
    units    = '1e-6 mol/mol'
    call shr_nuopc_fldList_AddMetadata('Sa_co2diag', longname, stdname, units)

    longname = 'Surface flux of CO2 from land'
    stdname  = 'surface_upward_flux_of_carbon_dioxide_where_land'
    units    = 'moles m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Fall_fco2_lnd', longname, stdname, units)

    longname = 'Surface flux of CO2 from ocean'
    stdname  = 'surface_upward_flux_of_carbon_dioxide_where_open_sea'
    units    = 'moles m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faoo_fco2_ocn', longname, stdname, units)

    !-----------------------------
    ! water isotope fields
    !-----------------------------

    longname = 'Ratio of ocean surface level abund. H2_16O/H2O/Rstd'
    stdname  = 'ratio_ocean_surface_16O_abund'
    units    = '1'
    call shr_nuopc_fldList_AddMetadata('So_roce_16O', longname, stdname, units)

    longname = 'Ratio of ocean surface level abund. HDO/H2O/Rstd'
    stdname  = 'ratio_ocean_surface_HDO_abund'
    call shr_nuopc_fldList_AddMetadata('So_roce_HDO', longname, stdname, units)

    !------------------------
    ! Atmospheric specific humidty at lowest level:
    !------------------------

    longname = 'Specific humidty of H216O at the lowest model level'
    stdname  = 'H216OV'
    units    = 'kg kg-1'
    call shr_nuopc_fldList_AddMetadata('Sa_shum_16O', longname, stdname, units)
    longname = 'Specific humidty of H218O at the lowest model level'
    stdname  = 'H218OV'
    call shr_nuopc_fldList_AddMetadata('Sa_shum_18O', longname, stdname, units)
    longname = 'Specific humidty of HD16O at the lowest model level'
    stdname  = 'HD16OV'
    call shr_nuopc_fldList_AddMetadata('Sa_shum_HDO', longname, stdname, units)

    !------------------------
    ! Isotopic surface snow water equivalent (land/atm only)
    !------------------------

    longname = 'Isotopic surface snow water equivalent'
    stdname  = 'surface_snow_water_equivalent'
    units    = 'm'
    call shr_nuopc_fldList_AddMetadata('Sl_snowh_16O', longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata('Sl_snowh_HDO', longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata('Sl_snowh_18O', longname, stdname, units)

    !------------------------
    ! Isotopic Precipitation Fluxes:
    !------------------------

    longname = 'H216O Convective precipitation rate'
    stdname  = 'H2_16O_convective_precipitation_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faxa_rainc_16O', longname, stdname, units)
    longname = 'H216O Large-scale (stable) precipitation rate'
    stdname  = 'H2_16O_large_scale_precipitation_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faxa_rainl_16O', longname, stdname, units)
    longname = 'Water flux due to H216O rain' !equiv. to bulk
    stdname  = 'H2_16O_rainfall_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faxa_rain_16O', longname, stdname, units)

    longname = 'H218O Convective precipitation rate'
    stdname  = 'H2_18O_convective_precipitation_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faxa_rainc_18O', longname, stdname, units)
    longname = 'H218O Large-scale (stable) precipitation rate'
    stdname  = 'H2_18O_large_scale_precipitation_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faxa_rainl_18O', longname, stdname, units)
    longname = 'Water flux due to H218O rain'
    stdname  = 'h2_18o_rainfall_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faxa_rain_18O', longname, stdname, units)

    longname = 'HDO Convective precipitation rate'
    stdname  = 'HDO_convective_precipitation_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faxa_rainc_HDO', longname, stdname, units)
    longname = 'HDO Large-scale (stable) precipitation rate'
    stdname  = 'HDO_large_scale_precipitation_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faxa_rainl_HDO', longname, stdname, units)
    longname = 'Water flux due to HDO rain'
    stdname  = 'hdo_rainfall_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faxa_rain_HDO', longname, stdname, units)

    !------------------------
    ! Isotopic Snow Fluxes:
    !------------------------

    longname = 'H216O Convective snow rate (water equivalent)'
    stdname  = 'H2_16O_convective_snowfall_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faxa_snowc_16O', longname, stdname, units)
    longname = 'H216O Large-scale (stable) snow rate'
    stdname  = 'H2_16O_large_scale_snowfall_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faxa_snowl_16O', longname, stdname, units)
    longname = 'Water flux due to H216O snow' !equiv. to bulk
    stdname  = 'H2_16O_snowfall_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faxa_snow_16O', longname, stdname, units)

    longname = 'H218O Convective snow rate'
    stdname  = 'H2_18O_convective_snowfall_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faxa_snowc_18O', longname, stdname, units)
    longname = 'H218O Large-scale (stable) snow rate'
    stdname  = 'H2_18O_large_scale_snowfall_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faxa_snowl_18O', longname, stdname, units)
    longname = 'Water flux due to H218O snow'
    stdname  = 'h2_18o_snowfall_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faxa_snow_18O', longname, stdname, units)

    longname = 'HDO Convective snow rate'
    stdname  = 'HDO_convective_snowfall_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faxa_snowc_HDO', longname, stdname, units)
    longname = 'HDO Large-scale (stable) snow rate'
    stdname  = 'HDO_large_scale_snowfall_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faxa_snowl_HDO', longname, stdname, units)
    longname = 'Water flux due to HDO snow'
    stdname  = 'hdo_snowfall_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faxa_snow_HDO', longname, stdname, units)

    !------------------------
    ! Isotopic precipitation (rain + snow)
    !------------------------

    longname = 'Isotopic Water flux (rain+snow) for H2_16O'
    stdname  = 'h2_16o_precipitation_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faxa_prec_16O', longname, stdname, units)
    longname = 'Isotopic Water flux (rain+snow) for H2_18O'
    stdname  = 'h2_18o_precipitation_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faxa_prec_18O', longname, stdname, units)
    longname = 'Isotopic Water flux (rain+snow) for H2_HDO'
    stdname  = 'h2_HDo_precipitation_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Faxa_prec_HDO', longname, stdname, units)

    !-------------------------------------
    ! Isotopic two meter reference humidity:
    !-------------------------------------

    longname = 'Reference H216O specific humidity at 2 meters'
    stdname  = 'H216O_specific_humidity'
    units    = 'kg kg-1'
    call shr_nuopc_fldList_AddMetadata('Sl_qref_16O', longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata('Si_qref_16O', longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata('So_qref_16O', longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata('Sx_qref_16O', longname, stdname, units)

    longname = 'Reference H218O specific humidity at 2 meters'
    stdname  = 'H218O_specific_humidity'
    units    = 'kg kg-1'
    call shr_nuopc_fldList_AddMetadata('Sl_qref_18O', longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata('Si_qref_18O', longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata('So_qref_18O', longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata('Sx_qref_18O', longname, stdname, units)

    longname = 'Reference H2HDO specific humidity at 2 meters'
    stdname  = 'H2HDO_specific_humidity'
    units    = 'kg kg-1'
    call shr_nuopc_fldList_AddMetadata('Sl_qref_HDO', longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata('Si_qref_HDO', longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata('So_qref_HDO', longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata('Sx_qref_HDO', longname, stdname, units)

    !-------------------------
    ! Isotopic Evaporation flux:
    !-------------------------

    longname = 'Evaporation H216O flux'
    stdname  = 'H216O_evaporation_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Fall_evap_16O', longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata('Faii_evap_16O', longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata('Faox_evap_16O', longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata('Faxx_evap_16O', longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata('Foxx_evap_16O', longname, stdname, units)

    longname = 'Evaporation H216O flux'
    stdname  = 'H216O_evaporation_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Fall_evap_18O', longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata('Faii_evap_18O', longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata('Faox_evap_18O', longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata('Faxx_evap_18O', longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata('Foxx_evap_18O', longname, stdname, units)

    longname = 'Evaporation H2HDO flux'
    stdname  = 'H2HDO_evaporation_flux'
    units    = 'kg m-2 s-1'
    call shr_nuopc_fldList_AddMetadata('Fall_evap_HDO', longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata('Faii_evap_HDO', longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata('Faox_evap_HDO', longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata('Faxx_evap_HDO', longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata('Foxx_evap_HDO', longname, stdname, units)

    !-------------------------
    ! Isotopic sea ice melting flux
    !-------------------------

    ! 'Heat flux from melting'
    units    = 'kg m-2 s-1'
    longname = 'H2_16O heat flux due to melting'
    stdname  = 'h2_16o_surface_snow_melt_hflux'
    call shr_nuopc_fldList_AddMetadata('Fioi_melth_16O', longname, stdname, units)
    units    = 'kg m-2 s-1'
    longname = 'H2_18O heat flux due to melting'
    stdname  = 'h2_18o_surface_snow_melt_hflux'
    call shr_nuopc_fldList_AddMetadata('Fioi_melth_18O', longname, stdname, units)
    units    = 'kg m-2 s-1'
    longname = 'H2_18O heat flux due to melting'
    stdname  = 'h2_HDo_surface_snow_melt_hflux'
    call shr_nuopc_fldList_AddMetadata('Fioi_melth_HDO', longname, stdname, units)

    ! 'Water flux from melting'
    units    = 'kg m-2 s-1'
    longname = 'H2_16O water flux due to melting'
    stdname  = 'h2_16o_surface_snow_melt_wflux'
    call shr_nuopc_fldList_AddMetadata('Fioi_meltw_16O', longname, stdname, units)
    units    = 'kg m-2 s-1'
    longname = 'H2_18O water flux due to melting'
    stdname  = 'h2_18o_surface_snow_melt_wflux'
    call shr_nuopc_fldList_AddMetadata('Fioi_meltw_18O', longname, stdname, units)
    units    = 'kg m-2 s-1'
    longname = 'H2_18O water flux due to melting'
    stdname  = 'h2_HDo_surface_snow_melt_wflux'
    call shr_nuopc_fldList_AddMetadata('Fioi_meltw_HDO', longname, stdname, units)

    !-----------------------------------------------------------------------------
    ! optional per thickness category fields
    !-----------------------------------------------------------------------------

    if (flds_i2o_per_cat) then
    
       do num = 1, ice_ncat
          write(cnum,'(i2.2)') num

          ! Fractional ice coverage wrt ocean
          longname = 'fractional ice coverage wrt ocean for thickness category ' // cnum
          stdname  = 'sea_ice_area_fraction'
          units    = '1'
          name = 'Si_ifrac_' // cnum
          call shr_nuopc_fldList_AddMetadata(trim(name), longname, stdname, units)

          ! Net shortwave radiation
          longname = 'net shortwave radiation penetrating into ice and ocean times ice fraction for thickness category ' // cnum
          stdname  = 'product_of_net_downward_shortwave_flux_at_sea_water_surface_and_sea_ice_area_fraction'
          units    = 'W m-2'
          name = 'PFioi_swpen_ifrac_' // cnum
          call shr_nuopc_fldList_AddMetadata(trim(name), longname, stdname, units)
       end do

       longname = 'fractional atmosphere coverage wrt ocean'
       stdname  = 'atmosphere_area_fraction'
       units    = '1'
       call shr_nuopc_fldList_AddMetadata('Sf_afrac', longname, stdname, units)

       longname = 'fractional atmosphere coverage used in radiation computations wrt ocean'
       stdname  = 'atmosphere_area_fraction'
       units    = '1'
       call shr_nuopc_fldList_AddMetadata('Sf_afracr', longname, stdname, units)

       longname = 'net shortwave radiation times atmosphere fraction'
       stdname = 'product_of_net_downward_shortwave_flux_at_sea_water_surface_and_atmosphere_area_fraction'
       units = 'W m-2'
       call shr_nuopc_fldList_AddMetadata('Foxx_swnet_afracr', longname, stdname, units)

    end if

    !-----------------------------------------------------------------------------
    ! CARMA fields
    ! if carma_flds are specified then setup fields for CLM to CAM communication
    !-----------------------------------------------------------------------------

    ! TODO: fill this in
    ! longname = 'Volumetric soil water'
    ! stdname  = 'soil_water'
    ! units    = 'm3/m3'
    ! carma_fields = 
    ! do n = 1,shr_string_listGetNum(carma_fields)
    !    call shr_string_listGetName(carma_fields, n, fldname)
    !    call shr_nuopc_fldList_AddMetadata(trim(fldname), longname, stdname, units)
    ! endif

    !-----------------------------------------------------------------------------
    ! MEGAN emissions fluxes fields
    ! if MEGAN emission are specified then setup fields for CLM to CAM communication
    !-----------------------------------------------------------------------------

    longname = 'MEGAN emission fluxes'
    stdname  = 'megan'
    units    = 'molecules/m2/sec'
    do num = 1, max_megan
       write(cnum,'(i2.2)') num
       fldname = 'Fall_voc' // cnum
       call shr_nuopc_fldList_AddMetadata(trim(fldname), longname, stdname, units)
    end do

    !-----------------------------------------------------------------------------
    ! Fire emissions fluxes fields
    !-----------------------------------------------------------------------------

    longname = 'wild fire emission fluxes'
    stdname  = 'fire_emis'
    units    = 'kg/m2/sec'
    do num = 1, max_fire
       write(cnum,'(i2.2)') num
       fldname  = 'Fall_fire' // cnum
       call shr_nuopc_fldList_AddMetadata(trim(fldname), longname, stdname, units)
    enddo

    longname = 'wild fire plume height'
    stdname  = 'fire_plume_top'
    units    = 'm'
    call shr_nuopc_fldList_AddMetadata('Sl_fztop', longname, stdname, units)

    !-----------------------------------------------------------------------------
    ! Dry Deposition fields
    !-----------------------------------------------------------------------------

    longname = 'dry deposition velocity'
    stdname  = 'drydep_vel'
    units    = 'cm/sec'
    do num = 1, max_ddep
       write(cnum,'(i2.2)') num
       fldname  = 'Sl_dd' // cnum
       call shr_nuopc_fldList_AddMetadata(trim(fldname), longname, stdname, units)
    end do

    !-----------------------------------------------------------------------------
    ! Nitrogen Deposition fields
    !-----------------------------------------------------------------------------

    longname = 'nitrogen deposition flux'
    stdname  = 'nitrogen_deposition'
    units    = 'kg(N)/m2/sec'
    call shr_nuopc_fldList_AddMetadata('Faxa_noy', longname, stdname, units)
    call shr_nuopc_fldList_AddMetadata('Faxa_nhx', longname, stdname, units)

  end subroutine esmDict_Init

end module esmDict
