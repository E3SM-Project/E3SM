module SFMainMod

  ! ============================================================================
  ! All subroutines realted to the SPITFIRE fire routine. 
  ! Code originally developed by Allan Spessa & Rosie Fisher as part of the NERC-QUEST project.  
  ! ============================================================================

  use shr_kind_mod          , only : r8 => shr_kind_r8;
  use spmdMod               , only : masterproc
  use clm_varctl            , only : iulog
  use atm2lndType           , only : atm2lnd_type
  use TemperatureType       , only : temperature_type
  use EcophysconType        , only : ecophyscon
  use EDEcophysconType      , only : EDecophyscon
  use EDtypesMod            , only : site, patch, cohort, AREA, DG_SF, FIRE_THRESHOLD
  use EDtypesMod            , only : LB_SF, LG_SF, NCWD, TR_SF

  implicit none
  save
  private

  public :: fire_model
  public :: fire_danger_index 
  public :: charecteristics_of_fuel
  public :: rate_of_spread
  public :: ground_fuel_consumption
  public :: fire_intensity
  public :: wind_effect
  public :: area_burnt
  public :: crown_scorching
  public :: crown_damage
  public :: cambial_damage_kill
  public :: post_fire_mortality

  integer :: write_SF = 0     ! for debugging
  logical :: DEBUG = .false.  ! for debugging

  ! ============================================================================
  ! ============================================================================

contains

  ! ============================================================================
  !        Area of site burned by fire           
  ! ============================================================================
  subroutine fire_model( site_in, atm2lnd_vars, temperature_vars)

    use clm_varctl,   only : use_ed_spit_fire

    implicit none

    type (site)             , intent(inout), pointer :: site_in
    type(atm2lnd_type)      , intent(in)             :: atm2lnd_vars
    type(temperature_type)  , intent(in)             :: temperature_vars

    type (patch), pointer :: currentPatch

    integer temporary_SF_switch
    !zero fire things
    currentPatch => site_in%youngest_patch
    temporary_SF_switch = 0
    do while(associated(currentPatch))
       currentPatch%frac_burnt = 0.0_r8
       currentPatch%AB         = 0.0_r8
       currentPatch%fire       = 0
       currentPatch => currentPatch%older
    enddo

    if(write_SF==1)then
       write(iulog,*) 'use_ed_spit_fire',use_ed_spit_fire
    endif

    if(use_ed_spit_fire.and.temporary_SF_switch==1)then
       call fire_danger_index(site_in, temperature_vars, atm2lnd_vars)
       call wind_effect(site_in, atm2lnd_vars)
       call charecteristics_of_fuel(site_in)
       call rate_of_spread(site_in)
       call ground_fuel_consumption(site_in)
       call fire_intensity(site_in)
       call area_burnt(site_in)
       call crown_scorching(site_in)
       call crown_damage(site_in)
       call cambial_damage_kill(site_in)
       call post_fire_mortality(site_in)
    end if

  end subroutine fire_model

    !*****************************************************************
    subroutine  fire_danger_index ( currentSite, temperature_vars, atm2lnd_vars)

    !*****************************************************************
   ! currentSite%acc_NI is the accumulated Nesterov fire danger index

    use clm_varcon       , only : tfrz

    use SFParamsMod, only  : SF_val_fdi_a, SF_val_fdi_b

    implicit none

    type(site)              , intent(inout), pointer :: currentSite
    type(temperature_type)  , intent(in)             :: temperature_vars
    type(atm2lnd_type)      , intent(in)             :: atm2lnd_vars
    
    real(r8) :: temp_in_C ! daily averaged temperature in celcius
    real(r8) :: rainfall  ! daily precip
    real(r8) :: rh        ! daily rh 
    
    real yipsolon; !intermediate varable for dewpoint calculation
    real dewpoint; !dewpoint in K 
    real d_NI;     !daily change in Nesterov Index. C^2 
  
    associate(                                                &
         t_veg24          => temperature_vars%t_veg24_patch , & ! Input:  [real(r8) (:)]  avg pft vegetation temperature for last 24 hrs    

         prec24           => atm2lnd_vars%prec24_patch      , & ! Input:  [real(r8) (:)]  avg pft rainfall for last 24 hrs    
         rh24             => atm2lnd_vars%rh24_patch          & ! Input:  [real(r8) (:)]  avg pft relative humidity for last 24 hrs    
         )

      ! NOTE: t_veg24(:), prec24(:) and rh24(:) are p level temperatures, precipitation and RH, 
      ! which probably won't have much inpact, unless we decide to ever calculated the NI for each patch.  

      temp_in_C  = t_veg24(currentSite%oldest_patch%clm_pno) - tfrz
      rainfall   = prec24(currentSite%oldest_patch%clm_pno) *24_r8*3600_r8
      rh         = rh24(currentSite%oldest_patch%clm_pno)

      if (rainfall > 3.0_r8) then !rezero NI if it rains... 
         d_NI = 0.0_r8
         currentSite%acc_NI = 0.0_r8
      else 
         yipsolon = (SF_val_fdi_a* temp_in_C)/(SF_val_fdi_b+ temp_in_C)+log(rh/100.0_r8) 
         dewpoint = (SF_val_fdi_b*yipsolon)/(SF_val_fdi_a-yipsolon) !Standard met. formula
         d_NI = ( temp_in_C-dewpoint)* temp_in_C !follows Nesterov 1968.  Equation 5. Thonicke et al. 2010.
         if (d_NI < 0.0_r8) then !Change in NI cannot be negative. 
            d_NI = 0.0_r8 !check 
         endif
      endif
      currentSite%acc_NI = currentSite%acc_NI + d_NI         !Accumulate Nesterov index over the fire season. 

    end associate

  end subroutine fire_danger_index


  !*****************************************************************
  subroutine  charecteristics_of_fuel ( currentSite )
  !*****************************************************************

    use SFParamsMod, only  : SF_val_alpha_FMC, SF_val_SAV, SF_val_FBD

    implicit none

    type(site), intent(in),   pointer :: currentSite

    type(patch),  pointer :: currentPatch
    type(cohort), pointer :: currentCohort

    real(r8) timeav_swc 
    real(r8) fuel_moisture(ncwd+2) ! Scaled moisture content of small litter fuels. 
    real(r8) MEF(ncwd+2)           ! Moisture extinction factor of fuels     integer n 

    fuel_moisture(:) = 0.0_r8
    
    currentPatch => currentSite%oldest_patch; 
    do while(associated(currentPatch))  
       ! How much live grass is there? 
       currentPatch%livegrass = 0.0_r8 
       currentCohort => currentPatch%tallest
       do while(associated(currentCohort))
          if(ecophyscon%woody(currentCohort%pft) == 0)then 
             currentPatch%livegrass = currentPatch%livegrass + currentCohort%bl*currentCohort%n/currentPatch%area
          endif
          currentCohort => currentCohort%shorter
       enddo
       
       ! There are SIX fuel classes
       ! 1) Leaf litter, 2:5) four CWD_AG pools (twig, s branch, l branch, trunk) and  6) live grass
       ! NCWD =4 
       ! dg_sf = 1, lb_sf, = 4, tr_sf = 5, lg_sf = 6,
     
            ! zero fire arrays. 
       currentPatch%fuel_eff_moist = 0.0_r8 
       currentPatch%fuel_bulkd     = 0.0_r8 
       currentPatch%fuel_sav       = 0.0_r8 
       currentPatch%fuel_frac(:)   = 0.0_r8 
       currentPatch%fuel_mef       = 0.0_r8
       currentPatch%sum_fuel       = 0.0_r8
       currentPatch%fuel_frac      = 0.0_r8

       if(write_sf == 1)then
          if (masterproc) write(iulog,*) ' leaf_litter1 ',currentPatch%leaf_litter
          if (masterproc) write(iulog,*) ' leaf_litter2 ',sum(currentPatch%CWD_AG)
          if (masterproc) write(iulog,*) ' leaf_litter3 ',currentPatch%livegrass
          if (masterproc) write(iulog,*) ' sum fuel', currentPatch%sum_fuel
       endif

       currentPatch%sum_fuel =  sum(currentPatch%leaf_litter) + sum(currentPatch%CWD_AG) + currentPatch%livegrass
       if(write_SF == 1)then
          if (masterproc) write(iulog,*) 'sum fuel', currentPatch%sum_fuel,currentPatch%area
       endif
       ! ===============================================
       ! Average moisture, bulk density, surface area-volume and moisture extinction of fuel
       ! ================================================   
                  
       if (currentPatch%sum_fuel > 0.0) then        
          ! Fraction of fuel in litter classes
          currentPatch%fuel_frac(dg_sf)       = sum(currentPatch%leaf_litter)/ currentPatch%sum_fuel
          currentPatch%fuel_frac(dg_sf+1:tr_sf) = currentPatch%CWD_AG          / currentPatch%sum_fuel    

          if(write_sf == 1)then
             if (masterproc) write(iulog,*) 'ff1 ',currentPatch%fuel_frac
             if (masterproc) write(iulog,*) 'ff2 ',currentPatch%fuel_frac
             if (masterproc) write(iulog,*) 'ff2a ',lg_sf,currentPatch%livegrass,currentPatch%sum_fuel
          endif

          currentPatch%fuel_frac(lg_sf)       = currentPatch%livegrass       / currentPatch%sum_fuel   
          MEF(1:ncwd+2)               = 0.524_r8 - 0.066_r8 * log10(SF_val_SAV(1:ncwd+2)) 

          !Equation 6 in Thonicke et al. 2010. 
          fuel_moisture(dg_sf+1:tr_sf)  = exp(-1.0_r8 * SF_val_alpha_FMC(dg_sf+1:tr_sf) * currentSite%acc_NI)  
          if(write_SF == 1)then
             if (masterproc) write(iulog,*) 'ff3 ',currentPatch%fuel_frac
             if (masterproc) write(iulog,*) 'fm ',fuel_moisture
             if (masterproc) write(iulog,*) 'csa ',currentSite%acc_NI
             if (masterproc) write(iulog,*) 'sfv ',SF_val_alpha_FMC
          endif
          ! FIX(RF,032414): needs refactoring. 
          ! average water content !is this the correct metric?         
          timeav_swc                  = sum(currentSite%water_memory(1:10)) / 10_r8 
          ! Equation B2 in Thonicke et al. 2010
          fuel_moisture(dg_sf)        = max(0.0_r8, 10.0_r8/9._r8 * timeav_swc - 1.0_r8/9.0_r8)           
 
          ! Average properties over the first four litter pools (dead leaves, twigs, s branches, l branches) 
          currentPatch%fuel_bulkd     = sum(currentPatch%fuel_frac(dg_sf:lb_sf) * SF_val_FBD(dg_sf:lb_sf))     
          currentPatch%fuel_sav       = sum(currentPatch%fuel_frac(dg_sf:lb_sf) * SF_val_SAV(dg_sf:lb_sf))              
          currentPatch%fuel_mef       = sum(currentPatch%fuel_frac(dg_sf:lb_sf) * MEF(dg_sf:lb_sf))              
          currentPatch%fuel_eff_moist = sum(currentPatch%fuel_frac(dg_sf:lb_sf) * fuel_moisture(dg_sf:lb_sf))         
          if(write_sf == 1)then
             if (masterproc) write(iulog,*) 'ff4 ',currentPatch%fuel_eff_moist
          endif
          ! Add on properties of live grass multiplied by grass fraction. (6)
          currentPatch%fuel_bulkd     = currentPatch%fuel_bulkd     + currentPatch%fuel_frac(lg_sf)  * SF_val_FBD(lg_sf)      
          currentPatch%fuel_sav       = currentPatch%fuel_sav       + currentPatch%fuel_frac(lg_sf)  * SF_val_SAV(lg_sf)
          currentPatch%fuel_mef       = currentPatch%fuel_mef       + currentPatch%fuel_frac(lg_sf)  * MEF(lg_sf)            
          currentPatch%fuel_eff_moist = currentPatch%fuel_eff_moist + currentPatch%fuel_frac(lg_sf)  * fuel_moisture(lg_sf)

          ! Correct averaging for the fact that we are not using the trunks pool (5)
          currentPatch%fuel_bulkd     = currentPatch%fuel_bulkd     * (1.0_r8/(1.0_r8-currentPatch%fuel_frac(tr_sf)))
          currentPatch%fuel_sav       = currentPatch%fuel_sav       * (1.0_r8/(1.0_r8-currentPatch%fuel_frac(tr_sf)))
          currentPatch%fuel_mef       = currentPatch%fuel_mef       * (1.0_r8/(1.0_r8-currentPatch%fuel_frac(tr_sf)))
          currentPatch%fuel_eff_moist = currentPatch%fuel_eff_moist * (1.0_r8/(1.0_r8-currentPatch%fuel_frac(tr_sf)))
          
          ! Convert from biomass to carbon. Which variables is this needed for?          
          currentPatch%fuel_bulkd = currentPatch%fuel_bulkd * 0.45_r8  
     
          ! Pass litter moisture into the fuel burning routine
          ! (wo/me term in Thonicke et al. 2010) 
          currentPatch%litter_moisture(dg_sf:lb_sf) = fuel_moisture(dg_sf:lb_sf)/MEF(dg_sf:lb_sf)  
          currentPatch%litter_moisture(tr_sf)       = 0.0_r8
          currentPatch%litter_moisture(lg_sf)       = fuel_moisture(lg_sf)/MEF(lg_sf)  

       else

          if(write_SF == 1)then

             if (masterproc) write(iulog,*) 'no litter fuel at all',currentPatch%patchno, &
                  currentPatch%sum_fuel,sum(currentPatch%cwd_ag),                         &
                  sum(currentPatch%cwd_bg),sum(currentPatch%leaf_litter)

          endif
          currentPatch%fuel_sav = sum(SF_val_SAV(1:ncwd+2))/(ncwd+2) ! make average sav to avoid crashing code. 

          if (masterproc) write(iulog,*) 'problem with spitfire fuel averaging'

          ! FIX(SPM,032414) refactor...should not have 0 fuel unless everything is burnt
          ! off.
          currentPatch%fuel_eff_moist = 0.0000000001_r8 
          currentPatch%fuel_bulkd     = 0.0000000001_r8 
          currentPatch%fuel_frac(:)   = 0.0000000001_r8 
          currentPatch%fuel_mef       = 0.0000000001_r8
          currentPatch%sum_fuel       = 0.0000000001_r8
          currentPatch%fuel_frac      = 0.0000000001_r8

       endif
       ! check values. 
       ! FIX(SPM,032414) refactor...
       if(write_SF == 1.and.currentPatch%fuel_sav <= 0.0_r8.or.currentPatch%fuel_bulkd <=  &
            0.0_r8.or.currentPatch%fuel_mef <= 0.0_r8.or.currentPatch%fuel_eff_moist <= 0.0_r8)then
            if (masterproc) write(iulog,*) 'problem with spitfire fuel averaging'
       endif 
       
       currentPatch => currentPatch%younger

    enddo !end patch loop
    
  end subroutine charecteristics_of_fuel


  !*****************************************************************
  subroutine  wind_effect ( currentSite, atm2lnd_vars)
  !*****************************************************************.

    ! Routine called daily from within ED within a site loop.
    ! Calculates the effective windspeed based on vegetation charecteristics. 

    implicit none

    type(site)         , intent(inout), pointer :: currentSite
    type(atm2lnd_type) , intent(in)             :: atm2lnd_vars

    type(patch) , pointer :: currentPatch
    type(cohort), pointer :: currentCohort

    ! note - this is a p level temperature, which probably won't have much inpact, 
    ! unless we decide to ever calculated the NI for each patch.  
    real(r8), pointer :: wind24(:) 

    real(r8) :: wind  ! daily wind
    real(r8) :: total_grass_area ! per patch,in m2
    real(r8) :: tree_fraction  !  site level. no units
    real(r8) :: grass_fraction !  site level. no units
    real(r8) :: bare_fraction  ! site level. no units 

    wind24  => atm2lnd_vars%wind24_patch  ! Input:  [real(r8) (:)]  avg pft windspeed (m/s)

    wind = wind24(currentSite%oldest_patch%clm_pno) * 60_r8             ! Convert to m/min for SPITFIRE units.
    if(write_SF == 1)then
       if (masterproc) write(iulog,*) 'wind24', wind24(currentSite%oldest_patch%clm_pno)
    endif
    ! --- influence of wind speed, corrected for surface roughness----
    ! --- averaged over the whole grid cell to prevent extreme divergence 
    ! average_wspeed = 0.0_r8   
    tree_fraction = 0.0_r8
    grass_fraction = 0.0_r8
    currentPatch=>currentSite%oldest_patch;  
    do while(associated(currentPatch))
       currentPatch%total_tree_area = 0.0_r8
       total_grass_area = 0.0_r8
       currentCohort => currentPatch%tallest
 
       do while(associated(currentCohort))
          write(iulog,*) 'SF currentCohort%c_area ',currentCohort%c_area
          if(ecophyscon%woody(currentCohort%pft) == 1)then
             currentPatch%total_tree_area = currentPatch%total_tree_area + currentCohort%c_area
          else
             total_grass_area = total_grass_area + currentCohort%c_area
          endif
          currentCohort => currentCohort%shorter
       enddo
       tree_fraction = tree_fraction + min(currentPatch%area,currentPatch%total_tree_area)/AREA
       grass_fraction = grass_fraction + min(currentPatch%area,total_grass_area)/AREA 
       
       if(DEBUG)then
         !write(iulog,*) 'SF  currentPatch%area ',currentPatch%area
         !write(iulog,*) 'SF  currentPatch%total_area ',currentPatch%total_tree_area
         !write(iulog,*) 'SF  total_grass_area ',tree_fraction,grass_fraction
         !write(iulog,*) 'SF  AREA ',AREA
       endif
       
       currentPatch => currentPatch%younger
    enddo !currentPatch loop

    !if there is a cover of more than one, then the grasses are under the trees
    grass_fraction = min(grass_fraction,1.0_r8-tree_fraction) 
    bare_fraction = 1.0 - tree_fraction - grass_fraction
    if(write_sf == 1)then
       if (masterproc) write(iulog,*) 'grass, trees, bare',grass_fraction, tree_fraction, bare_fraction
    endif

    currentPatch=>currentSite%oldest_patch;

    do while(associated(currentPatch))       
       currentPatch%total_tree_area = min(currentPatch%total_tree_area,currentPatch%area)      
       currentPatch%effect_wspeed = wind * (tree_fraction*0.6+grass_fraction*0.4+bare_fraction*1.0)
      
       currentPatch => currentPatch%younger
    enddo !end patch loop

  end subroutine wind_effect

  !*****************************************************************
  subroutine rate_of_spread ( currentSite ) 
    !*****************************************************************.
    !Routine called daily from within ED within a site loop.
    !Returns the updated currentPatch%ROS_front value for each patch.

    use SFParamsMod, only  : SF_val_miner_total, SF_val_part_dens, &
         SF_val_miner_damp, SF_val_fuel_energy

    implicit none

    type(site), intent(in), pointer :: currentSite

    type(patch), pointer :: currentPatch

    real(r8) dummy

    ! Rothermal fire spread model parameters. 
    real(r8) beta
    real(r8) ir !reaction intensity
    real(r8) xi,eps,q_ig,phi_wind
    real(r8) gamma_aptr,gamma_max
    real(r8) moist_damp,mw_weight
    real(r8) bet,beta_op
    real(r8) a,b,c,e

    currentPatch=>currentSite%oldest_patch;  

    do while(associated(currentPatch))
              
        ! ---initialise parameters to zero.--- 
       bet = 0.0_r8;   q_ig = 0.0_r8;   eps = 0.0_r8;   a = 0.0_r8;   b = 0.0_r8;   c = 0.0_r8;   e = 0.0_r8
       phi_wind = 0.0_r8;   xi = 0.0_r8;   gamma_max = 0.0_r8;   gamma_aptr = 0.0_r8;   mw_weight = 0.0_r8
       moist_damp = 0.0_r8;   ir = 0.0_r8;   dummy = 0.0_r8;     
       currentPatch%ROS_front = 0.0_r8
       currentPatch%sum_fuel  = currentPatch%sum_fuel * (1.0_r8 - SF_val_miner_total) !net of minerals

       ! ----start spreading---
       if (masterproc.and.DEBUG) write(iulog,*) 'SF - currentPatch%fuel_bulkd ',currentPatch%fuel_bulkd
       if (masterproc.and.DEBUG) write(iulog,*) 'SF - SF_val_part_dens ',SF_val_part_dens

       beta = (currentPatch%fuel_bulkd / 0.45_r8) / SF_val_part_dens
       
       ! Equation A6 in Thonicke et al. 2010
       beta_op = 0.200395_r8 *(currentPatch%fuel_sav**(-0.8189_r8))
       if (masterproc.and.DEBUG) write(iulog,*) 'SF - beta ',beta
       if (masterproc.and.DEBUG) write(iulog,*) 'SF - beta_op ',beta_op
       bet = beta/beta_op
       if(write_sf == 1)then
          if (masterproc) write(iulog,*) 'esf ',currentPatch%fuel_eff_moist
       endif
       ! ---heat of pre-ignition---
       !  Equation A4 in Thonicke et al. 2010 
       q_ig = 581.0_r8 +2594.0_r8 * currentPatch%fuel_eff_moist

       ! ---effective heating number---
       ! Equation A3 in Thonicke et al. 2010.  
       eps = exp(-4.528_r8 / currentPatch%fuel_sav)     
       ! Equation A7 in Thonicke et al. 2010
       b = 0.15988_r8 * (currentPatch%fuel_sav**0.54_r8)
       ! Equation A8 in Thonicke et al. 2010
       c = 7.47_r8 * (exp(-0.8711_r8 * (currentPatch%fuel_sav**0.55_r8)))
       ! Equation A9 in Thonicke et al. 2010. 
       e = 0.715_r8 * (exp(-0.01094_r8 * currentPatch%fuel_sav))
       ! Equation A5 in Thonicke et al. 2010

       if (DEBUG) then
          if (masterproc.and.DEBUG) write(iulog,*) 'SF - c ',c
          if (masterproc.and.DEBUG) write(iulog,*) 'SF - currentPatch%effect_wspeed ',currentPatch%effect_wspeed
          if (masterproc.and.DEBUG) write(iulog,*) 'SF - b ',b
          if (masterproc.and.DEBUG) write(iulog,*) 'SF - bet ',bet
          if (masterproc.and.DEBUG) write(iulog,*) 'SF - e ',e
       endif

       ! convert from m/min to ft/min for Rothermel ROS eqn
       phi_wind = c * ((3.281_r8*currentPatch%effect_wspeed)**b)*(bet**(-e)) 

       ! ---propagating flux----
       ! Equation A2 in Thonicke et al.        

       xi = (exp((0.792_r8 + 3.7597_r8 * (currentPatch%fuel_sav**0.5_r8)) * (beta+0.1_r8))) / &
            (192_r8+7.9095_r8 * currentPatch%fuel_sav)      
      
       ! ---reaction intensity----
       ! Equation in table A1 Thonicke et al. 2010. 
       a = 8.9033_r8 * (currentPatch%fuel_sav**(-0.7913_r8))
       dummy = exp(a*(1-bet))
       ! Equation in table A1 Thonicke et al. 2010. 
       gamma_max  = 1.0_r8 / (0.0591_r8 + 2.926_r8* (currentPatch%fuel_sav**(-1.5_r8)))
       gamma_aptr = gamma_max*(bet**a)*dummy

       mw_weight = currentPatch%fuel_eff_moist/currentPatch%fuel_mef
       
       ! Equation in table A1 Thonicke et al. 2010. 
       moist_damp = max(0.0_r8,(1.0_r8 - (2.59_r8 * mw_weight) + (5.11_r8 * (mw_weight**2.0_r8)) - &
            (3.52_r8*(mw_weight**3.0_r8))))

       ! FIX(SPM, 040114) ask RF if this should be an endrun
       ! if(write_SF == 1)then
       ! write(iulog,*) 'moist_damp' ,moist_damp,mw_weight,currentPatch%fuel_eff_moist,currentPatch%fuel_mef
       ! endif

       ir = gamma_aptr*(currentPatch%sum_fuel/0.45_r8)*SF_val_fuel_energy*moist_damp*SF_val_miner_damp 
       ! currentPatch%sum_fuel needs to be converted from kgC/m2 to kgBiomass/m2
       ! write(iulog,*) 'ir',gamma_aptr,moist_damp,SF_val_fuel_energy,SF_val_miner_damp
       if (((currentPatch%fuel_bulkd/0.45_r8) <= 0.0_r8).or.(eps <= 0.0_r8).or.(q_ig <= 0.0_r8)) then
          currentPatch%ROS_front = 0.0_r8
       else ! Equation 9. Thonicke et al. 2010. 
          currentPatch%ROS_front = (ir*xi*(1.0_r8+phi_wind)) / (currentPatch%fuel_bulkd/0.45_r8*eps*q_ig)
          ! write(iulog,*) 'ROS',currentPatch%ROS_front,phi_wind,currentPatch%effect_wspeed
          ! write(iulog,*) 'ros calcs',currentPatch%fuel_bulkd,ir,xi,eps,q_ig
       endif
       ! Equation 10 in Thonicke et al. 2010
       ! Can FBP System in m/min
       currentPatch%ROS_back = currentPatch%ROS_front*exp(-0.012_r8*currentPatch%effect_wspeed) 

       currentPatch => currentPatch%younger

    enddo !end patch loop

  end subroutine  rate_of_spread

  !*****************************************************************
  subroutine  ground_fuel_consumption ( cs_pnt ) 
  !*****************************************************************
    !returns the  the hypothetic fuel consumed by the fire

    use SFParamsMod, only : SF_val_miner_total, SF_val_min_moisture, &
         SF_val_mid_moisture, SF_val_low_moisture_C, SF_val_low_moisture_S, &
         SF_val_mid_moisture_C, SF_val_mid_moisture_S

    implicit none

    type(site), intent(in) :: cs_pnt

    type(patch), pointer      :: currentPatch

    real(r8) :: moist             !effective fuel moisture
    real(r8) :: tau_b(ncwd+2)     !lethal heating rates for each fuel class (min) 
    real(r8) :: fc_ground(ncwd+2) !propn of fuel consumed

    integer  :: c

    currentPatch => cs_pnt%oldest_patch;  

    do while(associated(currentPatch))
       currentPatch%burnt_frac_litter = 1.0_r8       
       ! Calculate fraction of litter is burnt for all classes. 
       ! Equation B1 in Thonicke et al. 2010---
       do c = 1, ncwd+2    !work out the burnt fraction for all pools, even if those pools dont exist.         
          moist = currentPatch%litter_moisture(c)                  
          ! 1. Very dry litter
          if (moist <= SF_val_min_moisture(c)) then
             currentPatch%burnt_frac_litter(c) = 1.0_r8  
          endif
          ! 2. Low to medium moistures
          if (moist > SF_val_min_moisture(c).and.moist <= SF_val_mid_moisture(c)) then
             currentPatch%burnt_frac_litter(c) = max(0.0_r8,min(1.0_r8,SF_val_low_moisture_C(c)- &
                  SF_val_low_moisture_S(c)*moist)) 
          else
          ! For medium to high moistures. 
             if (moist > SF_val_mid_moisture(c).and.moist <= 1.0_r8) then
                currentPatch%burnt_frac_litter(c) = max(0.0_r8,min(1.0_r8,SF_val_mid_moisture_C(c)- &
                     SF_val_mid_moisture_S(c)*moist))
             endif

          endif
          ! Very wet litter        
          if (moist >= 1.0_r8) then !this shouldn't happen? 
             currentPatch%burnt_frac_litter(c) = 0.0_r8  
          endif          
       enddo !c   

       ! we can't ever kill -all- of the grass. 
       currentPatch%burnt_frac_litter(lg_sf) = min(0.8_r8,currentPatch%burnt_frac_litter(lg_sf ))  
       ! reduce burnt amount for mineral content. 
       currentPatch%burnt_frac_litter = currentPatch%burnt_frac_litter * (1.0_r8-SF_val_miner_total) 

       !---Calculate amount of fuel burnt.---    
       FC_ground(dg_sf)   = currentPatch%burnt_frac_litter(dg_sf)   * sum(currentPatch%leaf_litter)
       FC_ground(2:tr_sf) = currentPatch%burnt_frac_litter(2:tr_sf) * currentPatch%CWD_AG
       FC_ground(lg_sf)   = currentPatch%burnt_frac_litter(lg_sf)   * currentPatch%livegrass      

       ! Following used for determination of cambial kill follows from Peterson & Ryan (1986) scheme 
       ! less empirical cf current scheme used in SPITFIRE which attempts to mesh Rothermel 
       ! and P&R, and while solving potential inconsistencies, actually results in BIG values for 
       ! fire residence time, thus lots of vegetation death!   
       ! taul is the duration of the lethal heating.  
       ! The /10 is to convert from kgC/m2 into gC/cm2, as in the Peterson and Ryan paper #Rosie,Jun 2013
        
       do c = 1,ncwd+2  
          tau_b(c)   =  39.4_r8 *(currentPatch%fuel_frac(c)*currentPatch%sum_fuel/0.45_r8/10._r8)* &
               (1.0_r8-((1.0_r8-currentPatch%burnt_frac_litter(c))**0.5_r8))  
       enddo
       tau_b(tr_sf)   =  0.0_r8
       ! Cap the residence time to 8mins, as suggested by literature survey by P&R (1986).
       currentPatch%tau_l = min(8.0_r8,sum(tau_b)) 

       !---calculate overall fuel consumed by spreading fire --- 
       ! ignore 1000hr fuels. Just interested in fuels affecting ROS   
       currentPatch%TFC_ROS = sum(FC_ground)-FC_ground(tr_sf)  

       currentPatch=>currentPatch%younger;
    enddo !end patch loop

  end subroutine ground_fuel_consumption

  !*****************************************************************
  subroutine  fire_intensity ( currentSite ) 
    !*****************************************************************
    !returns the updated currentPatch%FI value for each patch.

    !currentPatch%FI  average fire intensity of flaming front during day.  Backward ROS plays no role here. kJ/m/s or kW/m.
    !currentPatch%ROS_front  forward ROS (m/min) 
    !currentPatch%TFC_ROS total fuel consumed by flaming front (kgC/m2)

    use clm_varctl,   only : use_ed_spit_fire
    use SFParamsMod,  only : SF_val_fdi_alpha,SF_val_fuel_energy, &
         SF_val_max_durat, SF_val_durat_slope

    implicit none

    type(site), intent(in), pointer :: currentSite

    type(patch), pointer :: currentPatch

    real(r8) ROS !m/s
    real(r8) W !  kgBiomass/m2
    real(r8) :: d_fdi      !change in the NI on this day to give fire duration. 

    currentPatch => currentSite%oldest_patch;  

    do while(associated(currentPatch))
       ROS   = currentPatch%ROS_front / 60.0_r8 !m/min to m/sec 
       W     = currentPatch%TFC_ROS / 0.45_r8 !kgC/m2 to kgbiomass/m2
       currentPatch%FI = SF_val_fuel_energy * W * ROS !kj/m/s, or kW/m
       if(write_sf == 1)then
          if(masterproc) write(iulog,*) 'fire_intensity',currentPatch%fi,W,currentPatch%ROS_front
       endif
       !'decide_fire' subroutine shortened and put in here... 
       if (currentPatch%FI >= fire_threshold) then  ! 50kW/m is the threshold for a self-sustaining fire
          currentPatch%fire = 1 ! Fire...    :D
          
          ! This is like but not identical to equation 7 in Thonicke et al. 2010.  WHY? 
          d_FDI  = 1.0_r8 - exp(-SF_val_fdi_alpha*currentSite%acc_NI) !follows Venevsky et al GCB 2002 
          ! Equation 14 in Thonicke et al. 2010
          currentPatch%FD = SF_val_max_durat / (1.0_r8 + SF_val_max_durat * exp(SF_val_durat_slope*d_FDI))
          if(write_SF == 1)then
             if (masterproc) write(iulog,*) 'fire duration minutes',currentPatch%fd
          endif
          !equation 15 in Arora and Boer CTEM model.Average fire is 1 day long.
          !currentPatch%FD = 60.0_r8 * 24.0_r8 !no minutes in a day      
       else     
          currentPatch%fire = 0 ! No fire... :-/
          currentPatch%FD   = 0.0_r8      
       endif
       !  FIX(SPM,032414) needs a refactor
       !  FIX(RF,032414) : should happen outside of SF loop - doing all spitfire code is inefficient otherwise. 
       if(.not. use_ed_spit_fire)then   
          currentPatch%fire = 0 !fudge to turn fire off
       endif

       currentPatch => currentPatch%younger;
    enddo !end patch loop

  end subroutine fire_intensity


  !*****************************************************************
  subroutine  area_burnt ( currentSite ) 
    !*****************************************************************
    !currentPatch%AB  daily area burnt (m2)
    !currentPatch%NF    !Daily number of ignitions (lightning and human-caused), adjusted for size of patch. 

    use domainMod,     only : ldomain
    use EDParamsMod,   only : ED_val_nfires
    implicit none

    type(site), intent(inout), pointer :: currentSite
    type(patch), pointer :: currentPatch

    real lb !length to breadth ratio of fire ellipse
    real df  !distance fire has travelled forward
    real db !distance fire has travelled backward
    real(r8) gridarea
    real(r8) size_of_fire
    integer g

    currentSite%frac_burnt = 0.0_r8

    currentPatch => currentSite%oldest_patch;  
    do while(associated(currentPatch))
       currentPatch%AB = 0.0_r8
       currentPatch%frac_burnt = 0.0_r8
       lb = 0.0_r8; db = 0.0_r8; df = 0.0_r8

       if (currentPatch%fire == 1) then
       ! The feedback between vegetation structure and ellipse size if turned off for now, 
       ! to reduce the positive feedback in the syste,
       ! This will also be investigated by William Hoffmans proposal. 
          !      if (currentPatch%effect_wspeed < 16.67_r8) then !16.67m/min = 1km/hr 
          lb = 1.0_r8
          !     else 
          !FIX(RF,032414) FOR NO GRASS
          !        lb = currentPatch%total_canopy_area/currentPatch%area*(1.0_r8)+(8.729_r8 * &
          ! ((1.0_r8 -(exp(-0.03_r8 * 0.06_r8 * currentPatch%effect_wspeed)))**2.155_r8)) !&
          !&       +currentPatch%fpc_grass*(1.1_r8+((0.06_r8*currentPatch%effect_wspeed)**0.0464))

          !      endif

          !     if (lb > 8.0_r8)then
          !       lb = 8.0_r8  !Constraint Canadian Fire Behaviour System
          !     endif
          ! ---- calculate length of major axis---
          db = currentPatch%ROS_back  * currentPatch%FD !m
          df = currentPatch%ROS_front * currentPatch%FD !m

          ! --- calculate area burnt---
          if(lb > 0.0_r8) then
             g = currentSite%clmgcell
             gridarea = ldomain%area(g) *1000000_r8 !convert from km2 into m2
             currentPatch%NF = ldomain%area(g) * ED_val_nfires * currentPatch%area/area /365
             ! If there are 15  lightening strickes per year, per km2. (approx from NASA product) 
             ! then there are 15/365 s/km2 each day. 
     
             ! Equation 1 in Thonicke et al. 2010
             ! To Do: Connect here with the Li & Levis GDP fire suppression algorithm. 
             ! Equation 16 in arora and boer model.
             !currentPatch%ab = currentPatch%ab *3.0_r8
             size_of_fire = ((3.1416_r8/(4.0_r8*lb))*((df+db)**2.0_r8))
             currentPatch%AB = size_of_fire * currentPatch%nf 
             if (currentPatch%AB > gridarea*currentPatch%area/area) then !all of patch burnt. 

                if (masterproc) write(iulog,*) 'burnt all of patch',currentPatch%patchno, &
                     currentPatch%area/area,currentPatch%ab,currentPatch%area/area*gridarea   
                if (masterproc) write(iulog,*) 'ros',currentPatch%ROS_front,currentPatch%FD, &
                     currentPatch%NF,currentPatch%FI,size_of_fire

                if (masterproc) write(iulog,*) 'litter',currentPatch%sum_fuel,currentPatch%CWD_AG,currentPatch%leaf_litter
                ! turn km2 into m2. work out total area burnt. 
                currentPatch%AB = currentPatch%area *  gridarea/AREA 
             endif
             currentPatch%frac_burnt = currentPatch%AB / (gridarea*currentPatch%area/area)
             if(write_SF == 1)then
                if (masterproc) write(iulog,*) 'frac_burnt',currentPatch%frac_burnt
             endif
          endif
       endif! fire
       currentSite%frac_burnt = currentSite%frac_burnt + currentPatch%frac_burnt     

       currentPatch => currentPatch%younger;

    enddo !end patch loop

  end subroutine area_burnt

  !*****************************************************************
  subroutine  crown_scorching ( currentSite ) 
  !*****************************************************************
    !currentPatch%SH !average scorch height for the patch(m)
    !currentPatch%FI  average fire intensity of flaming front during day.  kW/m.

    use SFParamsMod,  only : SF_val_alpha_SH
    use EDParamsMod,  only : ED_val_ag_biomass

    implicit none

    type(site), intent(in), pointer :: currentSite

    type(patch), pointer :: currentPatch
    type(cohort), pointer :: currentCohort

    real f_ag_bmass      !fraction of a tree cohort's above-ground biomass as a proportion of total patch ag tree biomass.
    real tree_ag_biomass !total amount of above-ground tree biomass in patch. kgC/m2

    currentPatch => currentSite%oldest_patch;  
    do while(associated(currentPatch)) 

       tree_ag_biomass = 0.0_r8
       f_ag_bmass = 0.0_r8
       if (currentPatch%fire == 1) then
          currentCohort => currentPatch%tallest;
          do while(associated(currentCohort))  
             if (ecophyscon%woody(currentCohort%pft) == 1) then !trees only
                tree_ag_biomass = tree_ag_biomass+(currentCohort%bl+ED_val_ag_biomass* &
                     (currentCohort%bsw + currentCohort%bdead))*currentCohort%n
             endif !trees only

             currentCohort=>currentCohort%shorter;

          enddo !end cohort loop

          !This loop weights the scorch height for the contribution of each cohort to the overall biomass.   
          currentPatch%SH = 0.0_r8
          currentCohort => currentPatch%tallest;
          do while(associated(currentCohort))
             if (ecophyscon%woody(currentCohort%pft) == 1.and.(tree_ag_biomass > 0.0_r8)) then !trees only
                f_ag_bmass = ((currentCohort%bl+ED_val_ag_biomass*(currentCohort%bsw + &
                     currentCohort%bdead))*currentCohort%n)/tree_ag_biomass
                !equation 16 in Thonicke et al. 2010
                if(write_SF == 1)then
                   if (masterproc) write(iulog,*) 'currentPatch%SH',currentPatch%SH,f_ag_bmass
                endif
                !2/3 Byram (1959)
                currentPatch%SH = currentPatch%SH + f_ag_bmass * SF_val_alpha_SH * (currentPatch%FI**0.667_r8) 
             endif !trees only
             currentCohort=>currentCohort%shorter;
          enddo !end cohort loop
       endif !fire

       currentPatch => currentPatch%younger;  
    enddo !end patch loop

  end subroutine crown_scorching

  !*****************************************************************
  subroutine  crown_damage ( site_in )
    !*****************************************************************

    !returns the updated currentCohort%cfa value for each tree cohort within each patch.
    !currentCohort%cfa  proportion of crown affected by fire

    implicit none

    type(site), intent(in) :: site_in

    type(patch) , pointer :: currentPatch
    type(cohort), pointer :: currentCohort

    currentPatch => site_in%oldest_patch

    do while(associated(currentPatch)) 
       if (currentPatch%fire == 1) then

          currentCohort=>currentPatch%tallest

          do while(associated(currentCohort))  
             currentCohort%cfa = 0.0_r8
             if (ecophyscon%woody(currentCohort%pft) == 1) then !trees only
                ! Flames lower than bottom of canopy. 
                ! c%hite is height of cohort
                if (currentPatch%SH < (currentCohort%hite-currentCohort%hite*EDecophyscon%crown(currentCohort%pft))) then 
                   currentCohort%cfa = 0.0_r8
                else
                   ! Flames part of way up canopy. 
                   ! Equation 17 in Thonicke et al. 2010. 
                   ! flames over bottom of canopy but not over top.
                   if ((currentCohort%hite > 0.0_r8).and.(currentPatch%SH >=  &
                        (currentCohort%hite-currentCohort%hite*EDecophyscon%crown(currentCohort%pft)))) then 

                           currentCohort%cfa =  (currentPatch%SH-currentCohort%hite* &
                                EDecophyscon%crown(currentCohort%pft))/(currentCohort%hite-currentCohort%hite* &
                                EDecophyscon%crown(currentCohort%pft)) 

                   else 
                      ! Flames over top of canopy. 
                      currentCohort%cfa =  1.0_r8 
                   endif

                endif
                ! Check for strange values. 
                currentCohort%cfa = min(1.0_r8, max(0.0_r8,currentCohort%cfa))              
             endif !trees only
             !shrink canopy to account for burnt section.     
             !currentCohort%canopy_trim = min(currentCohort%canopy_trim,(1.0_r8-currentCohort%cfa)) 

             currentCohort => currentCohort%shorter;

          enddo !end cohort loop
       endif !fire?

       currentPatch => currentPatch%younger;

    enddo !end patch loop

  end subroutine crown_damage

  !*****************************************************************
  subroutine  cambial_damage_kill ( currentSite ) 
    !*****************************************************************
    ! routine description.
    ! returns the probability that trees dies due to cambial char
    ! currentPatch%tau_l = duration of lethal stem heating (min). Calculated at patch level.

    implicit none

    type(site), intent(in), pointer   :: currentSite

    type(patch), pointer  :: currentPatch
    type(cohort), pointer :: currentCohort

    real(r8) :: tau_c !critical time taken to kill cambium (minutes) 
    real(r8) :: bt    !bark thickness in cm.

    currentPatch => currentSite%oldest_patch;  

    do while(associated(currentPatch)) 

       if (currentPatch%fire == 1) then
          currentCohort => currentPatch%tallest;
          do while(associated(currentCohort))  
             if (ecophyscon%woody(currentCohort%pft) == 1) then !trees only
                ! Equation 21 in Thonicke et al 2010
                bt = EDecophyscon%bark_scaler(currentCohort%pft)*currentCohort%dbh ! bark thickness. 
                ! Equation 20 in Thonicke et al. 2010. 
                tau_c = 2.9_r8*bt**2.0_r8 !calculate time it takes to kill cambium (min)
                ! Equation 19 in Thonicke et al. 2010
                if ((currentPatch%tau_l/tau_c) >= 2.0_r8) then
                   currentCohort%cambial_mort = 1.0_r8
                else
                   if ((currentPatch%tau_l/tau_c) > 0.22_r8) then
                      currentCohort%cambial_mort = (0.563_r8*(currentPatch%tau_l/tau_c)) - 0.125_r8
                   else
                      currentCohort%cambial_mort = 0.0_r8
                   endif
                endif
             endif !trees 

             currentCohort => currentCohort%shorter;

          enddo !end cohort loop
       endif !fire?

       currentPatch=>currentPatch%younger;

    enddo !end patch loop

  end subroutine cambial_damage_kill

  !*****************************************************************
  subroutine  post_fire_mortality ( site_in )
  !*****************************************************************

    !  returns the updated currentCohort%fire_mort value for each tree cohort within each patch.
    !  currentCohort%cfa  proportion of crown affected by fire
    !  currentCohort%crownfire_mort  probability of tree post-fire mortality due to crown scorch
    !  currentCohort%cambial_mort  probability of tree post-fire mortality due to cambial char
    !  currentCohort%fire_mort  post-fire mortality from cambial and crown damage assuming two are independent.

    implicit none

    type(site), intent(in) :: site_in

    type(patch),  pointer :: currentPatch
    type(cohort), pointer :: currentCohort

    currentPatch => site_in%oldest_patch

    do while(associated(currentPatch)) 

       if (currentPatch%fire == 1) then 
          currentCohort => currentPatch%tallest
          do while(associated(currentCohort))  
             currentCohort%fire_mort = 0.0_r8
             currentCohort%crownfire_mort = 0.0_r8
             if (ecophyscon%woody(currentCohort%pft) == 1) then
                ! Equation 22 in Thonicke et al. 2010. 
                currentCohort%crownfire_mort = EDecophyscon%crown_kill(currentCohort%pft)*currentCohort%cfa**3.0_r8
                ! Equation 18 in Thonicke et al. 2010. 
                currentCohort%fire_mort = currentCohort%crownfire_mort+currentCohort%cambial_mort- &
                     (currentCohort%crownfire_mort*currentCohort%cambial_mort)  !joint prob.   
             else
                currentCohort%fire_mort = 0.0_r8 !I have changed this to zero and made the mode of death removal of leaves... 
             endif !trees

             currentCohort => currentCohort%shorter

          enddo !end cohort loop
       endif !fire?

       currentPatch => currentPatch%younger

    enddo !end patch loop

  end subroutine post_fire_mortality

  ! ============================================================================
end module SFMainMod
