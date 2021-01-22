module VegStructUpdateMod

  !-----------------------------------------------------------------------
  ! Module for vegetation structure updates (LAI, SAI, htop, hbot)
  !
  ! !USES:
  use shr_kind_mod         , only: r8 => shr_kind_r8
  use shr_sys_mod          , only : shr_sys_flush
  use shr_const_mod        , only : SHR_CONST_PI
  use elm_varctl           , only : iulog
  use VegetationPropertiesType     , only : veg_vp
  use WaterStateType       , only : waterstate_type
  use FrictionVelocityType , only : frictionvel_type
  use CNStateType          , only : cnstate_type
  use CNCarbonStateType    , only : carbonstate_type
  use CanopyStateType      , only : canopystate_type
  use CropType             , only : crop_type
  use ColumnDataType       , only : col_ws
  use VegetationType       , only : veg_pp
  use VegetationDataType   , only : veg_cs  
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: VegStructUpdate
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine VegStructUpdate(num_soilp, filter_soilp, &
       waterstate_vars, frictionvel_vars, cnstate_vars, &
       carbonstate_vars, canopystate_vars, crop_vars)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, use C state variables and epc to diagnose
    ! vegetation structure (LAI, SAI, height)
    !
    ! !USES:
    use pftvarcon        , only : noveg, nc3crop, nc3irrig, nbrdlf_evr_shrub, nbrdlf_dcd_brl_shrub
    use pftvarcon        , only : ncorn, ncornirrig, npcropmin, ztopmx, laimx
    use clm_time_manager , only : get_rad_step_size
    use elm_varctl       , only : spinup_state, spinup_mortality_factor
    !
    ! !ARGUMENTS:
    integer                , intent(in)    :: num_soilp       ! number of column soil points in pft filter
    integer                , intent(in)    :: filter_soilp(:) ! patch filter for soil points
    type(waterstate_type)  , intent(in)    :: waterstate_vars
    type(frictionvel_type) , intent(in)    :: frictionvel_vars
    type(cnstate_type)     , intent(inout) :: cnstate_vars
    type(carbonstate_type) , intent(in)    :: carbonstate_vars
    type(canopystate_type) , intent(inout) :: canopystate_vars
    type(crop_type)        , intent(inout) :: crop_vars
    !
    ! !REVISION HISTORY:
    ! 10/28/03: Created by Peter Thornton
    ! 2/29/08, David Lawrence: revised snow burial fraction for short vegetation
    !
    ! !LOCAL VARIABLES:
    integer  :: p,c,g      ! indices
    integer  :: fp         ! lake filter indices
    real(r8) :: taper      ! ratio of height:radius_breast_height (tree allometry)
    real(r8) :: stocking   ! #stems / ha (stocking density)
    real(r8) :: ol         ! thickness of canopy layer covered by snow (m)
    real(r8) :: fb         ! fraction of canopy layer covered by snow
    real(r8) :: tlai_old   ! for use in Zeng tsai formula
    real(r8) :: tsai_old   ! for use in Zeng tsai formula
    real(r8) :: tsai_min   ! PATCH derived minimum tsai
    real(r8) :: tsai_alpha ! monthly decay rate of tsai
    real(r8) :: dt         ! radiation time step (sec)

    real(r8), parameter :: dtsmonth = 2592000._r8 ! number of seconds in a 30 day month (60x60x24x30)
    !-----------------------------------------------------------------------
    ! tsai formula from Zeng et. al. 2002, Journal of Climate, p1835
    !
    ! tsai(p) = max( tsai_alpha(ivt(p))*tsai_old + max(tlai_old-tlai(p),0_r8), tsai_min(ivt(p)) )
    ! notes:
    ! * RHS tsai & tlai are from previous timestep
    ! * should create tsai_alpha(ivt(p)) & tsai_min(ivt(p)) in pftvarcon.F90 - slevis
    ! * all non-crop patches use same values:
    !   crop    tsai_alpha,tsai_min = 0.0,0.1
    !   noncrop tsai_alpha,tsai_min = 0.5,1.0  (includes bare soil and urban)
    !-------------------------------------------------------------------------------
    
    associate(                                                            & 
         ivt                =>  veg_pp%itype                         ,       & ! Input:  [integer  (:) ] pft vegetation type                                
         woody              =>  veg_vp%woody                  ,       & ! Input:  [real(r8) (:) ] binary flag for woody lifeform (1=woody, 0=not woody)
         slatop             =>  veg_vp%slatop                 ,       & ! Input:  [real(r8) (:) ] specific leaf area at top of canopy, projected area basis [m^2/gC]
         dsladlai           =>  veg_vp%dsladlai               ,       & ! Input:  [real(r8) (:) ] dSLA/dLAI, projected area basis [m^2/gC]           
         z0mr               =>  veg_vp%z0mr                   ,       & ! Input:  [real(r8) (:) ] ratio of momentum roughness length to canopy top height (-)
         displar            =>  veg_vp%displar                ,       & ! Input:  [real(r8) (:) ] ratio of displacement height to canopy top height (-)
         dwood              =>  veg_vp%dwood                  ,       & ! Input:  [real(r8) (:) ] density of wood (gC/m^3)                          

         snow_depth         =>  col_ws%snow_depth    ,       & ! Input:  [real(r8) (:) ] snow height (m)                                   

         forc_hgt_u_patch   =>  frictionvel_vars%forc_hgt_u_patch ,       & ! Input:  [real(r8) (:) ] observational height of wind at pft-level [m]     

         leafc              =>  veg_cs%leafc      ,       & ! Input:  [real(r8) (:) ] (gC/m2) leaf C                                    
         deadstemc          =>  veg_cs%deadstemc  ,       & ! Input:  [real(r8) (:) ] (gC/m2) dead stem C                               

         farea_burned       =>  cnstate_vars%farea_burned_col     ,       & ! Input:  [real(r8) (:) ] F. Li and S. Levis                                 
         harvdate           =>  crop_vars%harvdate_patch          ,       & ! Input:  [integer  (:) ] harvest date                                       
         htmx               =>  cnstate_vars%htmx_patch           ,       & ! Output: [real(r8) (:) ] max hgt attained by a crop during yr (m)          
         peaklai            =>  cnstate_vars%peaklai_patch        ,       & ! Output: [integer  (:) ] 1: max allowed lai; 0: not at max                  

         ! *** Key Output from CN***
         tlai               =>  canopystate_vars%tlai_patch       ,       & ! Output: [real(r8) (:) ] one-sided leaf area index, no burying by snow      
         tsai               =>  canopystate_vars%tsai_patch       ,       & ! Output: [real(r8) (:) ] one-sided stem area index, no burying by snow      
         htop               =>  canopystate_vars%htop_patch       ,       & ! Output: [real(r8) (:) ] canopy top (m)                                     
         hbot               =>  canopystate_vars%hbot_patch       ,       & ! Output: [real(r8) (:) ] canopy bottom (m)                                  
         elai               =>  canopystate_vars%elai_patch       ,       & ! Output: [real(r8) (:) ] one-sided leaf area index with burying by snow    
         esai               =>  canopystate_vars%esai_patch       ,       & ! Output: [real(r8) (:) ] one-sided stem area index with burying by snow    
         frac_veg_nosno_alb =>  canopystate_vars%frac_veg_nosno_alb_patch & ! Output: [integer  (:) ] frac of vegetation not covered by snow [-]         
         )

      dt = real( get_rad_step_size(), r8 )

      ! constant allometric parameters
      taper = 200._r8
      stocking = 1000._r8

      ! convert from stems/ha -> stems/m^2
      stocking = stocking / 10000._r8

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)
         c = veg_pp%column(p)
         g = veg_pp%gridcell(p)

         if (ivt(p) /= noveg) then

            tlai_old = tlai(p) ! n-1 value
            tsai_old = tsai(p) ! n-1 value

            ! update the leaf area index based on leafC and SLA
            ! Eq 3 from Thornton and Zimmerman, 2007, J Clim, 20, 3902-3923. 
            if (dsladlai(ivt(p)) > 0._r8) then
               tlai(p) = (slatop(ivt(p))*(exp(leafc(p)*dsladlai(ivt(p))) - 1._r8))/dsladlai(ivt(p))
            else
               tlai(p) = slatop(ivt(p)) * leafc(p)
            end if
            tlai(p) = max(0._r8, tlai(p))

            ! update the stem area index and height based on LAI, stem mass, and veg type.

            ! tsai formula from Zeng et. al. 2002, Journal of Climate, p1835 (see notes)
            ! Assumes doalb time step .eq. CLM time step, SAI min and monthly decay factor
            ! alpha are set by PFT, and alpha is scaled to CLM time step by multiplying by
            ! dt and dividing by dtsmonth (seconds in average 30 day month)
            ! tsai_min scaled by 0.5 to match MODIS satellite derived values
            if (ivt(p) == nc3crop .or. ivt(p) == nc3irrig) then ! generic crops

               tsai_alpha = 1.0_r8-1.0_r8*dt/dtsmonth
               tsai_min = 0.1_r8
            else
               tsai_alpha = 1.0_r8-0.5_r8*dt/dtsmonth
               tsai_min = 1.0_r8
            end if
            tsai_min = tsai_min * 0.5_r8
            tsai(p) = max(tsai_alpha*tsai_old+max(tlai_old-tlai(p),0._r8),tsai_min)

            if (woody(ivt(p)) == 1._r8) then

               ! trees and shrubs

               ! if shrubs have a squat taper 
               if (ivt(p) >= nbrdlf_evr_shrub .and. ivt(p) <= nbrdlf_dcd_brl_shrub) then
                  taper = 10._r8
                  ! otherwise have a tall taper
               else
                  taper = 200._r8
               end if

               ! trees and shrubs for now have a very simple allometry, with hard-wired
               ! stem taper (height:radius) and hard-wired stocking density (#individuals/area)
               if (spinup_state >= 1) then 
                 htop(p) = ((3._r8 * deadstemc(p) * spinup_mortality_factor * taper * taper)/ &
                      (SHR_CONST_PI * stocking * dwood(ivt(p))))**(1._r8/3._r8)
               else
                 htop(p) = ((3._r8 * deadstemc(p) * taper * taper)/ &
                      (SHR_CONST_PI * stocking * dwood(ivt(p))))**(1._r8/3._r8)
               end if

               ! Peter Thornton, 5/3/2004
               ! Adding test to keep htop from getting too close to forcing height for windspeed
               ! Also added for grass, below, although it is not likely to ever be an issue.
               htop(p) = min(htop(p),(forc_hgt_u_patch(p)/(displar(ivt(p))+z0mr(ivt(p))))-3._r8)

               ! Peter Thornton, 8/11/2004
               ! Adding constraint to keep htop from going to 0.0.
               ! This becomes an issue when fire mortality is pushing deadstemc
               ! to 0.0.
               htop(p) = max(htop(p), 0.01_r8)

               hbot(p) = max(0._r8, min(3._r8, htop(p)-1._r8))

            else if (ivt(p) >= npcropmin) then ! prognostic crops

               if (tlai(p) >= laimx(ivt(p))) peaklai(p) = 1 ! used in CNAllocation

               if (ivt(p) == ncorn .or. ivt(p) == ncornirrig) then
                  tsai(p) = 0.1_r8 * tlai(p)
               else
                  tsai(p) = 0.2_r8 * tlai(p)
               end if

               ! "stubble" after harvest
               if (harvdate(p) < 999 .and. tlai(p) == 0._r8) then
                  tsai(p) = 0.25_r8*(1._r8-farea_burned(c)*0.90_r8)    !changed by F. Li and S. Levis
                  htmx(p) = 0._r8
                  peaklai(p) = 0
               end if
               !if (harvdate(p) < 999 .and. tlai(p) > 0._r8) write(iulog,*) 'VegStructUpdate: tlai>0 after harvest!' ! remove after initial debugging?

               ! canopy top and bottom heights
               htop(p) = ztopmx(ivt(p)) * (min(tlai(p)/(laimx(ivt(p))-1._r8),1._r8))**2
               htmx(p) = max(htmx(p), htop(p))
               htop(p) = max(0.05_r8, max(htmx(p),htop(p)))
               hbot(p) = 0.02_r8

            else ! generic crops and ...

               ! grasses

               ! height for grasses depends only on LAI
               htop(p) = max(0.25_r8, tlai(p) * 0.25_r8)

               htop(p) = min(htop(p),(forc_hgt_u_patch(p)/(displar(ivt(p))+z0mr(ivt(p))))-3._r8)

               ! Peter Thornton, 8/11/2004
               ! Adding constraint to keep htop from going to 0.0.
               htop(p) = max(htop(p), 0.01_r8)

               hbot(p) = max(0.0_r8, min(0.05_r8, htop(p)-0.20_r8))
            end if

         else

            tlai(p) = 0._r8
            tsai(p) = 0._r8
            htop(p) = 0._r8
            hbot(p) = 0._r8

         end if

         ! adjust lai and sai for burying by snow. 
         ! snow burial fraction for short vegetation (e.g. grasses) as in
         ! Wang and Zeng, 2007.
         if (ivt(p) > noveg .and. ivt(p) <= nbrdlf_dcd_brl_shrub ) then
            ol = min( max(snow_depth(c)-hbot(p), 0._r8), htop(p)-hbot(p))
            fb = 1._r8 - ol / max(1.e-06_r8, htop(p)-hbot(p))
         else
            fb = 1._r8 - max(min(snow_depth(c),0.2_r8),0._r8)/0.2_r8   ! 0.2m is assumed
            !depth of snow required for complete burial of grasses
         endif

         elai(p) = max(tlai(p)*fb, 0.0_r8)
         esai(p) = max(tsai(p)*fb, 0.0_r8)

         ! Fraction of vegetation free of snow
         if ((elai(p) + esai(p)) > 0._r8) then
            frac_veg_nosno_alb(p) = 1
         else
            frac_veg_nosno_alb(p) = 0
         end if

      end do

    end associate 

 end subroutine VegStructUpdate

end module VegStructUpdateMod
