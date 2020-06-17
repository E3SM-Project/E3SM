module dynConsBiogeophysMod

#include "shr_assert.h"

  !---------------------------------------------------------------------------
  !
  ! !DESCRIPTION:
  ! Handle conservation of biogeophysical quantities (water & energy) with dynamic land
  ! cover.
  !
  ! !USES:
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use decompMod         , only : bounds_type
  use UrbanParamsType   , only : urbanparams_type
  use EnergyFluxType    , only : energyflux_type
  use LakeStateType     , only : lakestate_type
  use SoilHydrologyType , only : soilhydrology_type
  use SoilStateType     , only : soilstate_type
  use TotalWaterAndHeatMod, only : ComputeLiqIceMassNonLake, ComputeLiqIceMassLake
  use TotalWaterAndHeatMod, only : ComputeHeatNonLake, ComputeHeatLake
  use TotalWaterAndHeatMod, only : AdjustDeltaHeatForDeltaLiq
  use TotalWaterAndHeatMod, only : heat_base_temp
  use clm_varcon        , only : tfrz, cpliq
  use subgridAveMod     , only : p2c, c2g
  !use dynSubgridControlMod, only : get_for_testing_zero_dynbal_fluxes
  use clm_varcon        , only : spval
  use GridcellDataType  , only : grc_es, grc_ef, grc_ws, grc_wf
  use LandunitType      , only : lun_pp
  use ColumnType        , only : col_pp
  use VegetationType    , only : veg_pp
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  save
  public :: dyn_hwcontent_init            ! compute grid-level heat and water content, before land cover change
  public :: dyn_hwcontent_final           ! compute grid-level heat and water content, after land cover change; also compute dynbal fluxes
  !
  ! !PRIVATE MEMBER FUNCTIONS
  private :: dyn_water_content            ! compute gridcell total liquid and ice water contents
  private :: dyn_heat_content             ! compute gridcell total heat contents
  !---------------------------------------------------------------------------

contains

  !---------------------------------------------------------------------------
  subroutine dyn_hwcontent_init(bounds,                                      &
       num_nolakec, filter_nolakec,                                          &
       num_lakec, filter_lakec,                                              &
       urbanparams_vars, soilstate_vars, soilhydrology_vars, lakestate_vars, &
       energyflux_vars)
    !
    ! !DESCRIPTION:
    ! Initialize variables used for dyn_hwcontent, and compute grid cell-level heat
    ! and water content before land cover change
    !
    ! Should be called BEFORE any subgrid weight updates this time step
    !
    ! !ARGUMENTS:
      !$acc routine seq
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_nolakec
    integer                  , intent(in)    :: filter_nolakec(:)
    integer                  , intent(in)    :: num_lakec
    integer                  , intent(in)    :: filter_lakec(:)
    type(urbanparams_type)   , intent(in)    :: urbanparams_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    type(soilhydrology_type) , intent(in)    :: soilhydrology_vars
    type(lakestate_type)     , intent(in)    :: lakestate_vars
    type(energyflux_type)    , intent(inout) :: energyflux_vars
    !
    ! !LOCAL VARIABLES:
    integer :: g   ! grid cell index
    !-------------------------------------------------------------------------------
    associate( &
             liq1 => grc_ws%liq1 ,&
             ice1 => grc_ws%ice1 ,&
             heat1 => grc_es%heat1 ,&
             liquid_water_temp1 => grc_es%liquid_water_temp1 &
             )
    call dyn_water_content(bounds,                                        &
         num_nolakec, filter_nolakec,                                     &
         num_lakec, filter_lakec,                                         &
         soilhydrology_vars, lakestate_vars,             &
         liquid_mass = liq1(bounds%begg:bounds%endg), &
         ice_mass    = ice1(bounds%begg:bounds%endg))

    call dyn_heat_content( bounds,                                        &
         num_nolakec, filter_nolakec,                                     &
         num_lakec, filter_lakec,                                         &
         urbanparams_vars, soilstate_vars, soilhydrology_vars,            &
         heat_grc = heat1(bounds%begg:bounds%endg),  &
         liquid_water_temp_grc = liquid_water_temp1(bounds%begg:bounds%endg))
    end associate

  end subroutine dyn_hwcontent_init

  !---------------------------------------------------------------------------
  subroutine dyn_hwcontent_final(bounds, &
       num_nolakec, filter_nolakec, &
       num_lakec, filter_lakec, &
       urbanparams_vars, soilstate_vars, soilhydrology_vars, lakestate_vars &
       , energyflux_vars, dtime)
    !
    ! Should be called AFTER all subgrid weight updates this time step
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    !$acc routine seq
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_nolakec
    integer                  , intent(in)    :: filter_nolakec(:)
    integer                  , intent(in)    :: num_lakec
    integer                  , intent(in)    :: filter_lakec(:)
    type(urbanparams_type)   , intent(in)    :: urbanparams_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    type(soilhydrology_type) , intent(in)    :: soilhydrology_vars
    type(lakestate_type)     , intent(in)    :: lakestate_vars
    type(energyflux_type)    , intent(inout) :: energyflux_vars
    real(r8)                 , intent(in)    :: dtime ! land model time step (sec)

    !
    ! !LOCAL VARIABLES:
    integer  :: begg, endg
    integer  :: g     ! grid cell index
    real(r8) :: delta_liq(bounds%begg:bounds%endg)  ! change in gridcell h2o liq content
    real(r8) :: delta_ice(bounds%begg:bounds%endg)  ! change in gridcell h2o ice content
    real(r8) :: delta_heat(bounds%begg:bounds%endg) ! change in gridcell heat content
    !---------------------------------------------------------------------------

    associate( &
             liq1 => grc_ws%liq1 ,&
             ice1 => grc_ws%ice1 ,&
             heat1 => grc_es%heat1 ,&
             liq2 => grc_ws%liq2 ,&
             ice2 => grc_ws%ice2 ,&
             heat2 => grc_es%heat2 ,&
             liquid_water_temp2 => grc_es%liquid_water_temp2 &
             )
    begg = bounds%begg
    endg = bounds%endg

    call dyn_water_content(bounds, &
         num_nolakec, filter_nolakec, &
         num_lakec, filter_lakec, &
         soilhydrology_vars, lakestate_vars, &
         liquid_mass = liq2(bounds%begg:bounds%endg), &
         ice_mass    = ice2(bounds%begg:bounds%endg))

    call dyn_heat_content( bounds,                                &
         num_nolakec, filter_nolakec, &
         num_lakec, filter_lakec, &
         urbanparams_vars, soilstate_vars, soilhydrology_vars, &
         heat_grc = heat2(bounds%begg:bounds%endg), &
         liquid_water_temp_grc = liquid_water_temp2(bounds%begg:bounds%endg))

    !TODO: FIX
    if (.false.) then !get_for_testing_zero_dynbal_fluxes()) then
       do g = begg, endg
          delta_liq(g) = 0._r8
          delta_ice(g) = 0._r8
          delta_heat(g) = 0._r8
      end do
    else
       do g = begg, endg
          delta_liq(g)  = liq2(g) - liq1(g)
          delta_ice(g)  = ice2(g) - ice1(g)
          delta_heat(g) = heat2(g) - heat1(g)
          grc_wf%qflx_liq_dynbal (g) = delta_liq(g)/dtime
          grc_wf%qflx_ice_dynbal (g) = delta_ice(g)/dtime
          grc_ef%eflx_dynbal    (g) = delta_heat(g)/dtime
       end do
    end if
    
    end associate
  end subroutine dyn_hwcontent_final

  !-----------------------------------------------------------------------
  subroutine dyn_water_content(bounds, &
       num_nolakec, filter_nolakec, &
       num_lakec, filter_lakec, &
       soilhydrology_vars,  lakestate_vars, &
       liquid_mass, ice_mass)
    !
    ! !DESCRIPTION:
    ! Compute gridcell total liquid and ice water contents
    !
    ! !ARGUMENTS:
    !$acc routine seq
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_nolakec
    integer                  , intent(in)    :: filter_nolakec(:)
    integer                  , intent(in)    :: num_lakec
    integer                  , intent(in)    :: filter_lakec(:)
    type(soilhydrology_type) , intent(in)    :: soilhydrology_vars
    type(lakestate_type)     , intent(in)    :: lakestate_vars
    real(r8)                 , intent(out)   :: liquid_mass( bounds%begg: ) ! kg m-2
    real(r8)                 , intent(out)   :: ice_mass( bounds%begg: )    ! kg m-2
    !
    ! !LOCAL VARIABLES:
    real(r8) :: liquid_mass_col(bounds%begc:bounds%endc) ! kg m-2
    real(r8) :: ice_mass_col(bounds%begc:bounds%endc)    ! kg m-2


    !-----------------------------------------------------------------------
    call ComputeLiqIceMassNonLake(bounds, num_nolakec, filter_nolakec, &
         soilhydrology_vars,  &
         liquid_mass_col(bounds%begc:bounds%endc), &
         ice_mass_col(bounds%begc:bounds%endc))

    call ComputeLiqIceMassLake(bounds, num_lakec, filter_lakec, &
         lakestate_vars, &
         liquid_mass_col(bounds%begc:bounds%endc), &
         ice_mass_col(bounds%begc:bounds%endc))

    call c2g(bounds, &
         carr = liquid_mass_col(bounds%begc:bounds%endc), &
         garr = liquid_mass(bounds%begg:bounds%endg), &
         c2l_scale_type = 0, &
         l2g_scale_type = 0)

    call c2g(bounds, &
         carr = ice_mass_col(bounds%begc:bounds%endc), &
         garr = ice_mass(bounds%begg:bounds%endg), &
         c2l_scale_type = 0, &
         l2g_scale_type = 0)

  end subroutine dyn_water_content


  !---------------------------------------------------------------------------
  subroutine dyn_heat_content(bounds, &
       num_nolakec, filter_nolakec, &
       num_lakec, filter_lakec, &
       urbanparams_vars, soilstate_vars, soilhydrology_vars, &
       heat_grc, liquid_water_temp_grc)
    ! !DESCRIPTION:
    ! Compute grid-level heat and water content to track conservation with respect to
    ! dynamic land cover.
    !
    ! Heat content is computed relative to a baseline of 0 C. So temperatures above 0 C
    ! lead to a positive heat content, temperatures below 0 C lead to a negative heat
    ! content. For water, the baseline is considered to be ice at 0 C, so for liquid water
    ! we include the latent heat of fusion.
    !
    ! !ARGUMENTS:
      !$acc routine seq
    type(bounds_type)        , intent(in)  :: bounds
    integer                  , intent(in)  :: num_nolakec
    integer                  , intent(in)  :: filter_nolakec(:)
    integer                  , intent(in)  :: num_lakec
    integer                  , intent(in)  :: filter_lakec(:)
    type(urbanparams_type)   , intent(in)  :: urbanparams_vars
    type(soilstate_type)     , intent(in)  :: soilstate_vars
    type(soilhydrology_type) , intent(in)  :: soilhydrology_vars

    real(r8)                 , intent(out) :: heat_grc( bounds%begg: ) ! total heat content for each grid cell [J/m^2]
    real(r8)                 , intent(out) :: liquid_water_temp_grc( bounds%begg: ) ! weighted average liquid water temperature for each grid cell (K)

    !
    ! !LOCAL VARIABLES:
    integer  :: g

    real(r8) :: heat_col(bounds%begc:bounds%endc)  ! sum of heat content for all columns [J/m^2]
    real(r8) :: heat_liquid_col(bounds%begc:bounds%endc) ! sum of heat content for all columns: liquid water, excluding latent heat [J/m^2]
    real(r8) :: cv_liquid_col(bounds%begc:bounds%endc) ! sum of liquid heat capacity for all columns [J/(m^2 K)]

    real(r8) :: heat_liquid_grc(bounds%begg:bounds%endg) ! heat_liquid_col averaged to grid cell [J/m^2]
    real(r8) :: cv_liquid_grc(bounds%begg:bounds%endg) ! cv_liquid_col averaged to grid cell [J/(m^2 K)]
    !-------------------------------------------------------------------------------

    ! Enforce expected array sizes
    !!SHR_ASSERT_ALL((ubound(heat_grc) == (/bounds%endg/)), errMsg(__FILE__, __LINE__))
    !!SHR_ASSERT_ALL((ubound(liquid_water_temp_grc) == (/bounds%endg/)), errMsg(__FILE__, __LINE__))

    heat_col(bounds%begc:bounds%endc)        = spval
    heat_liquid_col(bounds%begc:bounds%endc) = spval
    cv_liquid_col(bounds%begc:bounds%endc)   = spval

    call ComputeHeatNonLake(bounds, num_nolakec, filter_nolakec, &
         urbanparams_vars, soilstate_vars, soilhydrology_vars, &
         heat = heat_col(bounds%begc:bounds%endc), &
         heat_liquid = heat_liquid_col(bounds%begc:bounds%endc), &
         cv_liquid = cv_liquid_col(bounds%begc:bounds%endc))

    call ComputeHeatLake(bounds, num_lakec, filter_lakec, &
         soilstate_vars, &
         heat = heat_col(bounds%begc:bounds%endc), &
         heat_liquid = heat_liquid_col(bounds%begc:bounds%endc), &
         cv_liquid = cv_liquid_col(bounds%begc:bounds%endc))

    call c2g(bounds, &
         carr = heat_col(bounds%begc:bounds%endc), &
         garr = heat_grc(bounds%begg:bounds%endg), &
         c2l_scale_type = 0, &
         l2g_scale_type = 0)

    call c2g(bounds, &
         carr = heat_liquid_col(bounds%begc:bounds%endc), &
         garr = heat_liquid_grc(bounds%begg:bounds%endg), &
         c2l_scale_type = 0, &
         l2g_scale_type = 0)

    call c2g(bounds, &
         carr = cv_liquid_col(bounds%begc:bounds%endc), &
         garr = cv_liquid_grc(bounds%begg:bounds%endg), &
         c2l_scale_type = 0, &
         l2g_scale_type = 0)

    do g = bounds%begg, bounds%endg
       if (cv_liquid_grc(g) > 0._r8) then
          liquid_water_temp_grc(g) = &
               (heat_liquid_grc(g) / cv_liquid_grc(g)) + heat_base_temp
       else
          ! 0 or negative water mass in this grid cell: set an arbitrary temperature
          liquid_water_temp_grc(g) = tfrz
       end if
    end do

  end subroutine dyn_heat_content

end module dynConsBiogeophysMod
