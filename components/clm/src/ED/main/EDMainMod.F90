module EDMainMod
  ! ===========================================================================
  ! Main ED module.    
  ! ============================================================================

  use shr_kind_mod         , only : r8 => shr_kind_r8
  use clm_varctl           , only : iulog
  use atm2lndType          , only : atm2lnd_type
  use SoilStateType        , only : soilstate_type  
  use TemperatureType      , only : temperature_type
  use WaterStateType       , only : waterstate_type

  use EDCohortDynamicsMod  , only : allocate_live_biomass, terminate_cohorts, fuse_cohorts, sort_cohorts, count_cohorts
  use EDPatchDynamicsMod   , only : disturbance_rates, fuse_patches, spawn_patches, terminate_patches
  use EDCLMLinkMod         , only : clm_ed_link
  use EDPhysiologyMod      , only : canopy_derivs, non_canopy_derivs, phenology, recruitment, trim_canopy
  use EDCanopyStructureMod , only : canopy_spread, canopy_structure
  use EDtypesMod           , only : site, patch, cohort, gridcell_edstate_type, ncwd, gridCellEdState, n_sub, numpft_ed, udata
  use SFMainMod            , only : fire_model

  implicit none
  save
  private

  public :: edmodel
  public :: ecosystem_dynamics
  public :: integrate_state_variables
  public :: ed_update_sites
  public :: total_balance_check

  logical :: DEBUG_main = .false.

  ! ============================================================================
  ! 10/30/09: Created by Rosie Fisher
  ! ============================================================================

contains

  subroutine edmodel( bounds, atm2lnd_vars, soilstate_vars, temperature_vars, waterstate_vars )
    ! ============================================================================
    !        Main ed model routine containing gridcell loop   
    ! ============================================================================
    !
    ! !USES:
    use clm_time_manager, only : get_days_per_year, get_curr_date
    use clm_time_manager, only : get_ref_date, timemgr_datediff 
    use decompMod       , only : bounds_type
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type)      , intent(in) :: bounds            ! bounds
    type(atm2lnd_type)     , intent(in) :: atm2lnd_vars
    type(soilstate_type)   , intent(in) :: soilstate_vars
    type(temperature_type) , intent(in) :: temperature_vars
    type(waterstate_type)  , intent(in) :: waterstate_vars
    !
    ! !LOCAL VARIABLES:
    type(site),pointer :: currentSite
    
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

    !timing statements. 
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

    !timing statements. 
    udata%time_period = mod( dayDiffInt , n_sub )

    ! where most things happen
    do g = bounds%begg,bounds%endg
       currentSite=> gridCellEdState(g)%spnt  
       if (currentSite%istheresoil == 1) then
          call ecosystem_dynamics(currentSite, &
               atm2lnd_vars, soilstate_vars, temperature_vars, waterstate_vars) 
       endif
    enddo

    ! updates site & patch information
    call ed_update_sites( bounds, gridCellEdState ) 

    write(iulog,*) 'leaving ed model',bounds%begg,bounds%endg,dayDiffInt

    ! FIX(SPM,032414) debug npp_acc exit

  end subroutine edmodel

  !-------------------------------------------------------------------------------!
  subroutine ecosystem_dynamics(site_in, &
       atm2lnd_vars, soilstate_vars, temperature_vars, waterstate_vars)
  
    ! ============================================================================
    !  Core of ed model, calling all subsequent vegetation dynamics routines         
    ! ============================================================================
    !
    ! !ARGUMENTS:
    implicit none 
    type(site)             , intent(inout), pointer :: site_in
    type(patch)            , pointer                :: currentPatch
    type(atm2lnd_type)     , intent(in)             :: atm2lnd_vars
    type(soilstate_type)   , intent(in)             :: soilstate_vars
    type(temperature_type) , intent(in)             :: temperature_vars
    type(waterstate_type)  , intent(in)             :: waterstate_vars
    !
    !*******************************************************************************************************/
    ! Fire, growth, biogeochemistry. 
    !*******************************************************************************************************/

    !FIX(SPM,032414) take this out.  On startup these values are all zero and on restart it
    !zeros out values read in the restart file
   
    call total_balance_check(site_in,0)
    
    call phenology(site_in%clmgcell, temperature_vars, waterstate_vars)

    call fire_model(site_in, atm2lnd_vars, temperature_vars)

    ! Calculate disturbance and mortality based on previous timestep vegetation.
    call disturbance_rates(site_in)         

    ! Integrate state variables from annual rates to daily timestep
    call integrate_state_variables(site_in, soilstate_vars, temperature_vars, waterstate_vars) 

    !*******************************************************************************************************/
    ! Reproduction, Recruitment and Cohort Dynamics : controls cohort organisation */
    !*******************************************************************************************************/

    currentPatch => site_in%oldest_patch
    do while (associated(currentPatch))                 

       ! adds small cohort of each PFT
       call recruitment(0,currentPatch)                

       currentPatch => currentPatch%younger
    enddo
       
    call total_balance_check(site_in,1)
    currentPatch => site_in%oldest_patch
    do while (associated(currentPatch))

       ! kills cohorts that are too small
       call terminate_cohorts(currentPatch)       

       ! puts cohorts in right order
       call sort_cohorts(currentPatch)            

       ! fuses similar cohorts
       call fuse_cohorts(currentPatch)            

       currentPatch => currentPatch%younger
    enddo
   
    call total_balance_check(site_in,2)

    !*******************************************************************************************************/
    ! Patch dynamics sub-routines: fusion, new patch creation (spwaning), termination.  */
    !*******************************************************************************************************/

    !make new patches from disturbed land
    call spawn_patches(site_in)       
   
    call total_balance_check(site_in,3)

    !fuse on the spawned patches.
    call fuse_patches(site_in)        
   
    call total_balance_check(site_in,4)

    !kill patches that are too small
    call terminate_patches(site_in)   
   
    call total_balance_check(site_in,5)

  end subroutine ecosystem_dynamics

  !-------------------------------------------------------------------------------!
  subroutine  integrate_state_variables(site_in, soilstate_vars, temperature_vars, waterstate_vars)
    !FIX(SPM,032414) refactor so everything goes through interface
    !
    ! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8;
    !
    ! !ARGUMENTS:
    implicit none
    type(site)             , intent(in), target :: site_in
    type(soilstate_type)   , intent(in)         :: soilstate_vars
    type(temperature_type) , intent(in)         :: temperature_vars
    type(waterstate_type)  , intent(in)         :: waterstate_vars
    !
    ! !LOCAL VARIABLES:
    type(patch)  , pointer  :: currentPatch
    type(cohort) , pointer :: currentCohort

    integer  :: c                     ! Counter for litter size class 
    integer  :: p                     ! Counter for PFT
    real(r8) :: small_no              ! to circumvent numerical errors that cause negative values of things that can't be negative
    real(r8) :: cohort_biomass_store  ! remembers the biomass in the cohort for balance checking
    !-----------------------------------------------------------------------

    small_no = 0.0000000000_r8  ! Obviously, this is arbitrary.  RF - changed to zero

    currentPatch => site_in%youngest_patch
    do while(associated(currentPatch))
       currentPatch%age = currentPatch%age + udata%deltat
       ! FIX(SPM,032414) valgrind 'Conditional jump or move depends on uninitialised value'
       if( currentPatch%age  <  0._r8 )then
          write(iulog,*) 'negative patch age?',site_in%clmgcell, currentPatch%age, &
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
      
       call non_canopy_derivs( currentPatch, temperature_vars, soilstate_vars, waterstate_vars )

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
            write(iulog,*) 'negative CWD_AG', currentPatch%cwd_ag(c),Site_in%lat,site_In%lon
            currentPatch%cwd_ag(c) = small_no
          endif
          if(currentPatch%cwd_bg(c)<small_no)then
            write(iulog,*) 'negative CWD_BG', currentPatch%cwd_bg(c),Site_in%lat,Site_in%lon
            currentPatch%cwd_bg(c) = small_no
          endif
       enddo

       do p = 1,numpft_ed
          if(currentPatch%leaf_litter(p)<small_no)then
            write(iulog,*) 'negative leaf litter', currentPatch%leaf_litter(p),Site_in%lat,Site_in%lon
            currentPatch%leaf_litter(p) = small_no
          endif
          if(currentPatch%root_litter(p)<small_no)then
               write(iulog,*) 'negative root litter', currentPatch%root_litter(p), &
               currentPatch%droot_litter_dt(p)* udata%deltat, &
               Site_in%lat,Site_in%lon
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

  end subroutine integrate_state_variables

  !-------------------------------------------------------------------------------!
  subroutine ed_update_sites( bounds, geds_local )
    ! ============================================================================
    ! Calls routines to consolidate the ED growth process.
    ! Canopy Structure to assign canopy layers to cohorts
    ! Canopy Spread to figure out the size of tree crowns
    ! Trim_canopy to figure out the target leaf biomass. 
    ! Extra recruitment to fill empty patches.  
    ! ============================================================================
    !
    ! !USES:
    use clm_varctl       , only : iulog
    use clm_time_manager , only : is_restart
    use decompMod        , only : bounds_type
    !
    ! !ARGUMENTS:
    implicit none 
    type(bounds_type), intent(in) :: bounds  ! clump bounds
    type(gridcell_edstate_type), target, intent(inout) :: geds_local( bounds%begg: )
    !
    ! !LOCAL VARIABLES:
    type (site) ,  pointer :: currentSite
    type (patch) , pointer :: currentPatch   

    integer :: g             ! Counter for gridcell
    integer :: cohort_number ! To print out the number of cohorts.  
    !-----------------------------------------------------------------------

    do g = bounds%begg,bounds%endg
       currentSite => geds_local(g)%spnt
       if(currentSite%istheresoil == 1)then
        
          call canopy_spread(currentSite)
          call total_balance_check(currentSite,6)
          call canopy_structure(currentSite)
          call total_balance_check(currentSite,7)

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

       endif
    enddo

  end subroutine ed_update_sites

  !-------------------------------------------------------------------------------!
  subroutine total_balance_check (cs_pnt,call_index )
    ! ============================================================================
    ! This routine looks at the carbon in and out of the ED model and compares it to 
    ! the change in total carbon stocks. 
    ! Fluxes in are NPP. Fluxes out are decay of CWD and litter into SOM pools.  
    ! currentSite%flux_out and currentSite%flux_in are set where they occur 
    ! in the code. 
    ! ============================================================================
    !
    ! !ARGUMENTS:
    implicit none
    type(site), intent(inout), target :: cs_pnt
    integer,    intent(in)            :: call_index
    !
    ! !LOCAL VARIABLES:
    real(r8) :: biomass_stock   ! total biomass   in KgC/site
    real(r8) :: litter_stock    ! total litter    in KgC/site
    real(r8) :: seed_stock      ! total seed mass in KgC/site
    real(r8) :: total_stock     ! total ED carbon in KgC/site
    real(r8) :: change_in_stock ! Change since last time we set currentSite%old_stock in this routine.  KgC/site
    real(r8) :: error           ! How much carbon did we gain or lose (should be zero!) 
    real(r8) :: net_flux        ! Difference between recorded fluxes in and out. KgC/site

    ! nb. There is no time associated with these variables because this routine can be called between any two arbitrary points in
    ! code, even if no time has passed. 
    ! Also, the carbon pools are per site/gridcell, so that we can account for the changing areas of patches. 

    type(site)   , pointer :: currentSite  
    type(patch)  , pointer :: currentPatch
    type(cohort) , pointer :: currentCohort
    !-----------------------------------------------------------------------

    change_in_stock = 0.0_r8
    biomass_stock   = 0.0_r8
    litter_stock    = 0.0_r8
    seed_stock      = 0.0_r8

    currentSite => cs_pnt
    if (currentSite%istheresoil == 1) then
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

  end subroutine total_balance_check

  ! ============================================================================
end module EDMainMod
