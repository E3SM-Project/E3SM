module EDMainMod

  ! ===========================================================================
  ! Main ED module.    
  ! ============================================================================

  use shr_kind_mod         , only : r8 => shr_kind_r8
  use decompMod            , only : bounds_type
  use clm_varctl           , only : iulog
  use atm2lndType          , only : atm2lnd_type
  use SoilStateType        , only : soilstate_type  
  use TemperatureType      , only : temperature_type
  use WaterStateType       , only : waterstate_type
  use EDCohortDynamicsMod  , only : allocate_live_biomass, terminate_cohorts, fuse_cohorts, sort_cohorts, count_cohorts
  use EDPatchDynamicsMod   , only : disturbance_rates, fuse_patches, spawn_patches, terminate_patches
  use EDPhysiologyMod      , only : canopy_derivs, non_canopy_derivs, phenology, recruitment, trim_canopy
  use SFMainMod            , only : fire_model
  use EDtypesMod           , only : ncwd, n_sub, numpft_ed, udata
  use EDtypesMod           , only : ed_site_type, ed_patch_type, ed_cohort_type
  use EDPhenologyType      , only : ed_phenology_type
  use EDCLMLinkMod         , only : ed_clm_type

  implicit none
  private

  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: ed_driver
  public  :: ed_update_site
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: ed_ecosystem_dynamics
  private :: ed_integrate_state_variables
  private :: ed_total_balance_check
  
  logical :: DEBUG_main = .false.
  !
  ! 10/30/09: Created by Rosie Fisher
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine ed_driver( bounds, ed_allsites_inst, ed_clm_inst, ed_phenology_inst, &
       atm2lnd_inst, soilstate_inst, temperature_inst, waterstate_inst, canopystate_inst)
    !
    ! !DESCRIPTION:
    ! Main ed model routine containing gridcell loop   
    !
    ! !USES:
    use clm_time_manager     , only : get_days_per_year, get_curr_date
    use clm_time_manager     , only : get_ref_date, timemgr_datediff 
    use CanopySTateType      , only : canopystate_type
    !
    ! !ARGUMENTS:
    type(bounds_type)       , intent(in)            :: bounds           
    type(ed_site_type)      , intent(inout), target :: ed_allsites_inst( bounds%begg: )
    type(ed_clm_type)       , intent(inout)         :: ed_clm_inst
    type(ed_phenology_type) , intent(inout)         :: ed_phenology_inst
    type(atm2lnd_type)      , intent(in)            :: atm2lnd_inst
    type(soilstate_type)    , intent(in)            :: soilstate_inst
    type(temperature_type)  , intent(in)            :: temperature_inst
    type(waterstate_type)   , intent(inout)         :: waterstate_inst
    type(canopystate_type)  , intent(inout)         :: canopystate_inst
    !
    ! !LOCAL VARIABLES:
    type(ed_site_type), pointer :: currentSite
    real(r8) :: dayDiff                  ! day of run
    integer  :: dayDiffInt               ! integer of day of run
    integer  :: g                        ! gridcell  
    integer  :: yr                       ! year (0, ...)
    integer  :: mon                      ! month (1, ..., 12)
    integer  :: day                      ! day of month (1, ..., 31)
    integer  :: sec                      ! seconds of the day
    integer  :: ncdate                   ! current date
    integer  :: nbdate                   ! base date (reference date)
    !-----------------------------------------------------------------------

    call ed_clm_inst%SetValues( bounds, 0._r8 )

    ! timing statements. 
    n_sub = get_days_per_year()
    udata%deltat = 1.0_r8/n_sub !for working out age of patches in years        
    if(udata%time_period == 0)then             
       udata%time_period = n_sub
    endif
    
    call get_curr_date(yr, mon, day, sec)
    ncdate = yr*10000 + mon*100 + day
    call get_ref_date(yr, mon, day, sec)
    nbdate = yr*10000 + mon*100 + day

    call timemgr_datediff(nbdate, 0, ncdate, sec, dayDiff)

    dayDiffInt = floor(dayDiff)
    udata%time_period = mod( dayDiffInt , n_sub )

    ! where most things happen
    do g = bounds%begg,bounds%endg
       if (ed_allsites_inst(g)%istheresoil) then
          currentSite => ed_allsites_inst(g)
          call ed_ecosystem_dynamics(currentSite, &
               ed_clm_inst, ed_phenology_inst, atm2lnd_inst, &
               soilstate_inst, temperature_inst, waterstate_inst)

          call ed_update_site( ed_allsites_inst(g))
       endif
    enddo

    ! updates site & patch information

    ! link to CLM structures
    call ed_clm_inst%ed_clm_link( bounds, ed_allsites_inst(bounds%begg:bounds%endg),  &
         ed_phenology_inst, waterstate_inst, canopystate_inst)

    write(iulog,*) 'leaving ed model',bounds%begg,bounds%endg,dayDiffInt

  end subroutine ed_driver

  !-------------------------------------------------------------------------------!
  subroutine ed_ecosystem_dynamics(currentSite, &
       ed_clm_inst, ed_phenology_inst, atm2lnd_inst, &
       soilstate_inst, temperature_inst, waterstate_inst)
    !
    ! !DESCRIPTION:
    !  Core of ed model, calling all subsequent vegetation dynamics routines         
    !
    ! !ARGUMENTS:
    type(ed_site_type)      , intent(inout), pointer :: currentSite
    type(ed_phenology_type) , intent(in)             :: ed_phenology_inst
    type(ed_clm_type)       , intent(in)             :: ed_clm_inst
    type(atm2lnd_type)      , intent(in)             :: atm2lnd_inst
    type(soilstate_type)    , intent(in)             :: soilstate_inst
    type(temperature_type)  , intent(in)             :: temperature_inst
    type(waterstate_type)   , intent(in)             :: waterstate_inst
    !
    ! !LOCAL VARIABLES:
    type(ed_patch_type), pointer :: currentPatch
    !-----------------------------------------------------------------------

    !**************************************************************************
    ! Fire, growth, biogeochemistry. 
    !**************************************************************************

    !FIX(SPM,032414) take this out.  On startup these values are all zero and on restart it
    !zeros out values read in the restart file
   
    call ed_total_balance_check(currentSite, 0)
    
    call phenology(currentSite, ed_phenology_inst, temperature_inst, waterstate_inst)

    call fire_model(currentSite, atm2lnd_inst, temperature_inst)

    ! Calculate disturbance and mortality based on previous timestep vegetation.
    call disturbance_rates(currentSite)

    ! Integrate state variables from annual rates to daily timestep
    call ed_integrate_state_variables(currentSite, soilstate_inst, temperature_inst, waterstate_inst) 

    !******************************************************************************
    ! Reproduction, Recruitment and Cohort Dynamics : controls cohort organisation 
    !******************************************************************************

    currentPatch => currentSite%oldest_patch
    do while (associated(currentPatch))                 

       ! adds small cohort of each PFT
       call recruitment(0,currentPatch)                

       currentPatch => currentPatch%younger
    enddo
       
    call ed_total_balance_check(currentSite,1)

    currentPatch => currentSite%oldest_patch
    do while (associated(currentPatch))

       ! kills cohorts that are too small
       call terminate_cohorts(currentPatch)       

       ! puts cohorts in right order
       call sort_cohorts(currentPatch)            

       ! fuses similar cohorts
       call fuse_cohorts(currentPatch)            

       currentPatch => currentPatch%younger
    enddo
   
    call ed_total_balance_check(currentSite,2)

    !*********************************************************************************
    ! Patch dynamics sub-routines: fusion, new patch creation (spwaning), termination.
    !*********************************************************************************

    ! make new patches from disturbed land
    call spawn_patches(currentSite)       
   
    call ed_total_balance_check(currentSite,3)

    ! fuse on the spawned patches.
    call fuse_patches(currentSite)        
   
    call ed_total_balance_check(currentSite,4)

    ! kill patches that are too small
    call terminate_patches(currentSite)   
   
    call ed_total_balance_check(currentSite,5)

  end subroutine ed_ecosystem_dynamics

  !-------------------------------------------------------------------------------!
  subroutine ed_integrate_state_variables(currentSite, soilstate_inst, temperature_inst, waterstate_inst)
    !
    ! !DESCRIPTION:
    ! FIX(SPM,032414) refactor so everything goes through interface
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(ed_site_type)     , intent(in) :: currentSite
    type(soilstate_type)   , intent(in) :: soilstate_inst
    type(temperature_type) , intent(in) :: temperature_inst
    type(waterstate_type)  , intent(in) :: waterstate_inst
    !
    ! !LOCAL VARIABLES:
    type(ed_patch_type)  , pointer :: currentPatch
    type(ed_cohort_type) , pointer :: currentCohort

    integer  :: c                     ! Counter for litter size class 
    integer  :: p                     ! Counter for PFT
    real(r8) :: small_no              ! to circumvent numerical errors that cause negative values of things that can't be negative
    real(r8) :: cohort_biomass_store  ! remembers the biomass in the cohort for balance checking
    !-----------------------------------------------------------------------

    small_no = 0.0000000000_r8  ! Obviously, this is arbitrary.  RF - changed to zero

    currentPatch => currentSite%youngest_patch

    do while(associated(currentPatch))

       currentPatch%age = currentPatch%age + udata%deltat
       ! FIX(SPM,032414) valgrind 'Conditional jump or move depends on uninitialised value'
       if( currentPatch%age  <  0._r8 )then
          write(iulog,*) 'negative patch age?',currentSite%clmgcell, currentPatch%age, &
               currentPatch%patchno,currentPatch%area
       endif

       ! Find the derivatives of the growth and litter processes. 
       call canopy_derivs(currentPatch)
       
       ! Update Canopy Biomass Pools
       currentCohort => currentPatch%shortest
       do while(associated(currentCohort)) 

          cohort_biomass_store  = (currentCohort%balive+currentCohort%bdead+currentCohort%bstore)
          currentCohort%dbh    = max(small_no,currentCohort%dbh    + currentCohort%ddbhdt    * udata%deltat )
          currentCohort%balive = currentCohort%balive + currentCohort%dbalivedt * udata%deltat 
          currentCohort%bdead  = max(small_no,currentCohort%bdead  + currentCohort%dbdeaddt  * udata%deltat )
          currentCohort%bstore = currentCohort%bstore + currentCohort%dbstoredt * udata%deltat 

          if( (currentCohort%balive+currentCohort%bdead+currentCohort%bstore)*currentCohort%n<0._r8)then
            write(iulog,*) 'biomass is negative', currentCohort%n,currentCohort%balive, &
                 currentCohort%bdead,currentCohort%bstore
          endif

          if(abs((currentCohort%balive+currentCohort%bdead+currentCohort%bstore+udata%deltat*(currentCohort%md+ &
               currentCohort%seed_prod)-cohort_biomass_store)-currentCohort%npp_acc) > 1e-8_r8)then
             write(iulog,*) 'issue with c balance in integration', abs(currentCohort%balive+currentCohort%bdead+ &
                  currentCohort%bstore+udata%deltat* &
                 (currentCohort%md+currentCohort%seed_prod)-cohort_biomass_store-currentCohort%npp_acc)
          endif  
          !do we need these any more?        
          currentCohort%npp_acc  = 0.0_r8
          currentCohort%gpp_acc  = 0.0_r8
          currentCohort%resp_acc = 0.0_r8
          
          call allocate_live_biomass(currentCohort)
  
          currentCohort => currentCohort%taller

       enddo
      
       write(6,*)'DEBUG18: calling non_canopy_derivs with pno= ',currentPatch%clm_pno
       call non_canopy_derivs( currentPatch, temperature_inst, soilstate_inst, waterstate_inst )

       !update state variables simultaneously according to derivatives for this time period. 
       do p = 1,numpft_ed
          currentPatch%seed_bank(p) = currentPatch%seed_bank(p) + currentPatch%dseed_dt(p)*udata%deltat
       enddo

       do c = 1,ncwd
          currentPatch%cwd_ag(c) =  currentPatch%cwd_ag(c) + currentPatch%dcwd_ag_dt(c)* udata%deltat
          currentPatch%cwd_bg(c) =  currentPatch%cwd_bg(c) + currentPatch%dcwd_bg_dt(c)* udata%deltat
       enddo

       do p = 1,numpft_ed
          currentPatch%leaf_litter(p) = currentPatch%leaf_litter(p) + currentPatch%dleaf_litter_dt(p)* udata%deltat
          currentPatch%root_litter(p) = currentPatch%root_litter(p) + currentPatch%droot_litter_dt(p)* udata%deltat
       enddo

       ! Check for negative values. Write out warning to show carbon balance. 
       do p = 1,numpft_ed
          if(currentPatch%seed_bank(p)<small_no)then
            write(iulog,*) 'negative seedbank', currentPatch%seed_bank(p)
            currentPatch%seed_bank(p) = small_no
          endif
       enddo

       do c = 1,ncwd
          if(currentPatch%cwd_ag(c)<small_no)then
            write(iulog,*) 'negative CWD_AG', currentPatch%cwd_ag(c),CurrentSite%lat,currentSite%lon
            currentPatch%cwd_ag(c) = small_no
          endif
          if(currentPatch%cwd_bg(c)<small_no)then
            write(iulog,*) 'negative CWD_BG', currentPatch%cwd_bg(c),CurrentSite%lat,CurrentSite%lon
            currentPatch%cwd_bg(c) = small_no
          endif
       enddo

       do p = 1,numpft_ed
          if(currentPatch%leaf_litter(p)<small_no)then
            write(iulog,*) 'negative leaf litter', currentPatch%leaf_litter(p),CurrentSite%lat,CurrentSite%lon
            currentPatch%leaf_litter(p) = small_no
          endif
          if(currentPatch%root_litter(p)<small_no)then
               write(iulog,*) 'negative root litter', currentPatch%root_litter(p), &
               currentPatch%droot_litter_dt(p)* udata%deltat, &
               CurrentSite%lat,CurrentSite%lon
            currentPatch%root_litter(p) = small_no
          endif
       enddo

     
       ! update cohort number. This needs to happen after the CWD_input and seed_input calculations as they 
       ! assume the pre-mortality currentCohort%n. 
       currentCohort => currentPatch%shortest
       do while(associated(currentCohort)) 
         currentCohort%n = max(small_no,currentCohort%n + currentCohort%dndt * udata%deltat )  
         currentCohort => currentCohort%taller
       enddo

       currentPatch => currentPatch%older

    enddo

  end subroutine ed_integrate_state_variables

  !-------------------------------------------------------------------------------!
  subroutine ed_update_site( currentSite )
    !
    ! !DESCRIPTION:
    ! Calls routines to consolidate the ED growth process.
    ! Canopy Structure to assign canopy layers to cohorts
    ! Canopy Spread to figure out the size of tree crowns
    ! Trim_canopy to figure out the target leaf biomass. 
    ! Extra recruitment to fill empty patches.  
    !
    ! !USES:
    use EDCanopyStructureMod , only : canopy_spread, canopy_structure
    !
    ! !ARGUMENTS:
    type(ed_site_type) , intent(inout), target :: currentSite
    !
    ! !LOCAL VARIABLES:
    type (ed_patch_type) , pointer :: currentPatch   
    integer :: cohort_number ! To print out the number of cohorts.  
    integer :: g             ! Counter for sites
    !-----------------------------------------------------------------------

    call canopy_spread(currentSite)

    call ed_total_balance_check(currentSite,6)

    call canopy_structure(currentSite)

    call ed_total_balance_check(currentSite,7)

    currentPatch => currentSite%oldest_patch
    do while(associated(currentPatch))

       call terminate_cohorts(currentPatch) 

       ! FIX(SPM,040314) why is this needed for BFB restarts? Look into this at some point
       cohort_number = count_cohorts(currentPatch)  
       if (DEBUG_main) then
          write(iulog,*) 'tempCount ',cohort_number
       endif
       ! Note (RF)
       ! This breaks the balance check, but if we leave it out, then 
       ! the first new patch that isn't fused has no cohorts at the end of the spawn process
       ! and so there are radiation errors instead. 
       ! Fixing this would likely require a re-work of how seed germination works which would be tricky. 
       if(currentPatch%countcohorts < 1)then
          !write(iulog,*) 'ED: calling recruitment for no cohorts',currentPatch%siteptr%clmgcell,currentPatch%patchno
          !call recruitment(1,currentPatch)
          ! write(iulog,*) 'patch empty',currentPatch%area,currentPatch%age
       endif

       currentPatch => currentPatch%younger    

    enddo

    ! FIX(RF,032414). This needs to be monthly, not annual
    if((udata%time_period == N_SUB-1))then 
       write(iulog,*) 'calling trim canopy' 
       call trim_canopy(currentSite)  
    endif

  end subroutine ed_update_site

  !-------------------------------------------------------------------------------!
  subroutine ed_total_balance_check (currentSite, call_index )
    !
    ! !DESCRIPTION:
    ! This routine looks at the carbon in and out of the ED model and compares it to 
    ! the change in total carbon stocks. 
    ! Fluxes in are NPP. Fluxes out are decay of CWD and litter into SOM pools.  
    ! ed_allsites_inst%flux_out and ed_allsites_inst%flux_in are set where they occur 
    ! in the code. 
    !
    ! !ARGUMENTS:
    type(ed_site_type) , intent(inout) :: currentSite
    integer            , intent(in)    :: call_index
    !
    ! !LOCAL VARIABLES:
    real(r8) :: biomass_stock   ! total biomass   in KgC/site
    real(r8) :: litter_stock    ! total litter    in KgC/site
    real(r8) :: seed_stock      ! total seed mass in KgC/site
    real(r8) :: total_stock     ! total ED carbon in KgC/site
    real(r8) :: change_in_stock ! Change since last time we set ed_allsites_inst%old_stock in this routine.  KgC/site
    real(r8) :: error           ! How much carbon did we gain or lose (should be zero!) 
    real(r8) :: net_flux        ! Difference between recorded fluxes in and out. KgC/site

    ! nb. There is no time associated with these variables 
    ! because this routine can be called between any two 
    ! arbitrary points in code, even if no time has passed. 
    ! Also, the carbon pools are per site/gridcell, so that 
    ! we can account for the changing areas of patches. 

    type(ed_patch_type)  , pointer :: currentPatch
    type(ed_cohort_type) , pointer :: currentCohort
    !-----------------------------------------------------------------------

    change_in_stock = 0.0_r8
    biomass_stock   = 0.0_r8
    litter_stock    = 0.0_r8
    seed_stock      = 0.0_r8

    if (currentSite%istheresoil) then
       currentPatch => currentSite%oldest_patch 
       do while(associated(currentPatch))

          litter_stock = litter_stock + currentPatch%area * (sum(currentPatch%cwd_ag)+ &
               sum(currentPatch%cwd_bg)+sum(currentPatch%leaf_litter)+sum(currentPatch%root_litter))
          seed_stock   = seed_stock   + currentPatch%area * sum(currentPatch%seed_bank)
          currentCohort => currentPatch%tallest;

          do while(associated(currentCohort))

             biomass_stock =  biomass_stock + (currentCohort%bdead + currentCohort%balive + &
                  currentCohort%bstore) * currentCohort%n
             currentCohort => currentCohort%shorter;

          enddo !end cohort loop 

          currentPatch => currentPatch%younger

       enddo !end patch loop

    endif

    total_stock     = biomass_stock + seed_stock +litter_stock
    change_in_stock = total_stock - currentSite%old_stock  
    net_flux        = currentSite%flux_in - currentSite%flux_out
    error           = abs(net_flux - change_in_stock)   

    if ( abs(error) > 10e-6 ) then
       write(iulog,*) 'total error:in,out,net,dstock,error',call_index, currentSite%flux_in, & 
            currentSite%flux_out,net_flux,change_in_stock,error
       write(iulog,*) 'biomass,litter,seeds', biomass_stock,litter_stock,seed_stock
       write(iulog,*) 'lat lon',currentSite%lat,currentSite%lon
    endif

    currentSite%flux_in   = 0.0_r8
    currentSite%flux_out  = 0.0_r8  
    currentSite%old_stock = total_stock

  end subroutine ed_total_balance_check

end module EDMainMod
