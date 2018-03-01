module EDPhysiologyMod

#include "shr_assert.h"

  ! ============================================================================
  ! Miscellaneous physiology routines from ED. 
  ! ============================================================================

  use FatesGlobals, only         : fates_log
  use FatesInterfaceMod, only    : hlm_days_per_year
  use FatesInterfaceMod, only    : hlm_model_day
  use FatesInterfaceMod, only    : hlm_freq_day
  use FatesInterfaceMod, only    : hlm_day_of_year
  use FatesInterfaceMod, only    : numpft
  use FatesConstantsMod, only    : r8 => fates_r8
  use EDPftvarcon      , only    : EDPftvarcon_inst
  use FatesInterfaceMod, only    : bc_in_type
  use EDCohortDynamicsMod , only : allocate_live_biomass, zero_cohort
  use EDCohortDynamicsMod , only : create_cohort, sort_cohorts

  use EDTypesMod          , only : numWaterMem
  use EDTypesMod          , only : dl_sf, dinc_ed
  use EDTypesMod          , only : external_recruitment
  use EDTypesMod          , only : ncwd
  use EDTypesMod          , only : nlevleaf
  use EDTypesMod          , only : senes
  use EDTypesMod          , only : maxpft
  use EDTypesMod          , only : ed_site_type, ed_patch_type, ed_cohort_type

  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use FatesGlobals          , only : fates_log
  use FatesGlobals          , only : endrun => fates_endrun
  use EDParamsMod           , only : fates_mortality_disturbance_fraction
  use FatesConstantsMod        , only : itrue,ifalse

  implicit none
  private

  public :: canopy_derivs
  public :: non_canopy_derivs
  public :: trim_canopy
  public :: phenology
  private :: phenology_leafonoff
  private :: Growth_Derivatives
  public :: recruitment
  private :: cwd_input
  private :: cwd_out
  private :: fragmentation_scaler
  private :: seeds_in
  private :: seed_decay
  private :: seed_germination
  public :: flux_into_litter_pools

  logical, parameter :: DEBUG  = .false. ! local debug flag
  character(len=*), parameter, private :: sourcefile = &
        __FILE__
  ! ============================================================================

contains

  ! ============================================================================
  subroutine canopy_derivs( currentSite, currentPatch, bc_in )
    !
    ! !DESCRIPTION:
    ! spawn new cohorts of juveniles of each PFT             
    !
    ! !USES:
    !
    ! !ARGUMENTS    
    type(ed_site_type), intent(inout), target  :: currentSite
    type(ed_patch_type) , intent(inout), target :: currentPatch
    type(bc_in_type), intent(in)               :: bc_in
    !
    ! !LOCAL VARIABLES:
    type(ed_cohort_type), pointer ::currentCohort
    !----------------------------------------------------------------------

    ! call plant growth functions

    currentCohort => currentPatch%shortest

    do while(associated(currentCohort))
       call Growth_Derivatives(currentSite, currentCohort, bc_in )
       currentCohort => currentCohort%taller
    enddo

  end subroutine canopy_derivs

  ! ============================================================================
  subroutine non_canopy_derivs( currentSite, currentPatch, bc_in )
    !
    ! !DESCRIPTION:
    ! Returns time differentials of the state vector
    !
    ! !USES:
    use EDTypesMod, only : AREA
    !
    ! !ARGUMENTS    
    type(ed_site_type), intent(inout), target  :: currentSite
    type(ed_patch_type), intent(inout)         :: currentPatch
    type(bc_in_type), intent(in)               :: bc_in

    !
    ! !LOCAL VARIABLES:
    integer c,p
    !----------------------------------------------------------------------

    currentPatch%leaf_litter_in(:)   = 0.0_r8
    currentPatch%root_litter_in(:)   = 0.0_r8
    currentPatch%dleaf_litter_dt(:)  = 0.0_r8
    currentPatch%droot_litter_dt(:)  = 0.0_r8
    currentPatch%leaf_litter_out(:)  = 0.0_r8
    currentPatch%root_litter_out(:)  = 0.0_r8
    currentPatch%cwd_AG_in(:)        = 0.0_r8
    currentPatch%cwd_BG_in(:)        = 0.0_r8
    currentPatch%cwd_AG_out(:)       = 0.0_r8
    currentPatch%cwd_BG_out(:)       = 0.0_r8
    currentPatch%seeds_in(:)         = 0.0_r8  
    currentPatch%seed_decay(:)       = 0.0_r8
    currentPatch%seed_germination(:) = 0.0_r8

    ! update seed fluxes 
    call seeds_in(currentSite, currentPatch)
    call seed_decay(currentSite, currentPatch)
    call seed_germination(currentSite, currentPatch)

    ! update fragmenting pool fluxes
    call cwd_input( currentSite, currentPatch)
    call cwd_out( currentSite, currentPatch, bc_in)

    do p = 1,numpft
       currentSite%dseed_dt(p) = currentSite%dseed_dt(p) + &
            (currentPatch%seeds_in(p) - currentPatch%seed_decay(p) - &
            currentPatch%seed_germination(p)) * currentPatch%area/AREA
    enddo   
    
    do c = 1,ncwd
       currentPatch%dcwd_AG_dt(c) = currentPatch%cwd_AG_in(c) - currentPatch%cwd_AG_out(c) 
       currentPatch%dcwd_BG_dt(c) = currentPatch%cwd_BG_in(c) - currentPatch%cwd_BG_out(c) 
    enddo

    do p = 1,numpft
       currentPatch%dleaf_litter_dt(p) = currentPatch%leaf_litter_in(p) - &
             currentPatch%leaf_litter_out(p) 
       currentPatch%droot_litter_dt(p) = currentPatch%root_litter_in(p) - &
             currentPatch%root_litter_out(p) 
    enddo

  end subroutine non_canopy_derivs

  ! ============================================================================
  subroutine trim_canopy( currentSite )
    !
    ! !DESCRIPTION:
    ! Canopy trimming / leaf optimisation. Removes leaves in negative annual carbon balance. 
    !
    ! !USES:
    !
    use EDGrowthFunctionsMod, only : tree_lai
    !
    ! !ARGUMENTS    
    type (ed_site_type),intent(inout), target :: currentSite
    !
    ! !LOCAL VARIABLES:
    type (ed_cohort_type) , pointer :: currentCohort
    type (ed_patch_type)  , pointer :: currentPatch

    integer  :: z          ! leaf layer
    integer  :: trimmed    ! was this layer trimmed in this year? If not expand the canopy. 

    !----------------------------------------------------------------------

    currentPatch => currentSite%youngest_patch

    do while(associated(currentPatch))
       currentCohort => currentPatch%tallest
       do while (associated(currentCohort)) 
          trimmed = 0    
          currentCohort%treelai = tree_lai(currentCohort)    
          currentCohort%nv = ceiling((currentCohort%treelai+currentCohort%treesai)/dinc_ed)
          if (currentCohort%nv > nlevleaf)then
             write(fates_log(),*) 'nv > nlevleaf',currentCohort%nv,currentCohort%treelai,currentCohort%treesai, &
                  currentCohort%c_area,currentCohort%n,currentCohort%bl
          endif

          !Leaf cost vs netuptake for each leaf layer. 
          do z = 1,nlevleaf
             if (currentCohort%year_net_uptake(z) /= 999._r8)then !there was activity this year in this leaf layer. 
                !Leaf Cost kgC/m2/year-1
                !decidous costs. 
                if (EDPftvarcon_inst%season_decid(currentCohort%pft) == 1.or. &
                     EDPftvarcon_inst%stress_decid(currentCohort%pft) == 1)then 
                   currentCohort%leaf_cost =  1._r8/(EDPftvarcon_inst%slatop(currentCohort%pft)*1000.0_r8)
                   currentCohort%leaf_cost = currentCohort%leaf_cost + &
                        1.0_r8/(EDPftvarcon_inst%slatop(currentCohort%pft)*1000.0_r8) * &
                        EDPftvarcon_inst%allom_l2fr(currentCohort%pft) / EDPftvarcon_inst%root_long(currentCohort%pft)
                   currentCohort%leaf_cost = currentCohort%leaf_cost * (EDPftvarcon_inst%grperc(currentCohort%pft) + 1._r8)
                else !evergreen costs
                   currentCohort%leaf_cost = 1.0_r8/(EDPftvarcon_inst%slatop(currentCohort%pft)* &
                        EDPftvarcon_inst%leaf_long(currentCohort%pft)*1000.0_r8) !convert from sla in m2g-1 to m2kg-1
                   currentCohort%leaf_cost = currentCohort%leaf_cost + &
                        1.0_r8/(EDPftvarcon_inst%slatop(currentCohort%pft)*1000.0_r8) * &
                        EDPftvarcon_inst%allom_l2fr(currentCohort%pft) / EDPftvarcon_inst%root_long(currentCohort%pft)
                   currentCohort%leaf_cost = currentCohort%leaf_cost * (EDPftvarcon_inst%grperc(currentCohort%pft) + 1._r8)
                endif
                if (currentCohort%year_net_uptake(z) < currentCohort%leaf_cost)then
                   if (currentCohort%canopy_trim > EDPftvarcon_inst%trim_limit(currentCohort%pft))then

                      if ( DEBUG ) then
                         write(fates_log(),*) 'trimming leaves',currentCohort%canopy_trim,currentCohort%leaf_cost
                      endif

                      ! keep trimming until none of the canopy is in negative carbon balance.              
                      if (currentCohort%hite > EDPftvarcon_inst%hgt_min(currentCohort%pft))then
                         currentCohort%canopy_trim = currentCohort%canopy_trim - EDPftvarcon_inst%trim_inc(currentCohort%pft)
                         if (EDPftvarcon_inst%evergreen(currentCohort%pft) /= 1)then
                            currentCohort%laimemory = currentCohort%laimemory*(1.0_r8 - EDPftvarcon_inst%trim_inc(currentCohort%pft)) 
                         endif
                         trimmed = 1
                      endif
                   endif
                endif
             endif !leaf activity? 
          enddo !z
          if (currentCohort%NV.gt.2)then
             ! leaf_cost may be uninitialized, removing its diagnostic from the log
             ! to allow checking with fpe_traps (RGK)
             write(fates_log(),*) 'nv>4',currentCohort%year_net_uptake(1:6),currentCohort%canopy_trim
          endif

          currentCohort%year_net_uptake(:) = 999.0_r8
          if (trimmed == 0.and.currentCohort%canopy_trim < 1.0_r8)then
             currentCohort%canopy_trim = currentCohort%canopy_trim + EDPftvarcon_inst%trim_inc(currentCohort%pft)
          endif 

          if ( DEBUG ) then
             write(fates_log(),*) 'trimming',currentCohort%canopy_trim
          endif
         
          ! currentCohort%canopy_trim = 1.0_r8 !FIX(RF,032414) this turns off ctrim for now. 
          currentCohort => currentCohort%shorter
       enddo
       currentPatch => currentPatch%older
    enddo

  end subroutine trim_canopy

  ! ============================================================================
  subroutine phenology( currentSite, bc_in )
    !
    ! !DESCRIPTION:
    ! Phenology. 
    !
    ! !USES:
    use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm
    use EDParamsMod, only : ED_val_phen_drought_threshold, ED_val_phen_doff_time
    use EDParamsMod, only : ED_val_phen_a, ED_val_phen_b, ED_val_phen_c, ED_val_phen_chiltemp
    use EDParamsMod, only : ED_val_phen_mindayson, ED_val_phen_ncolddayslim, ED_val_phen_coldtemp
         

    !
    ! !ARGUMENTS:
    type(ed_site_type), intent(inout), target :: currentSite
    type(bc_in_type),   intent(in)            :: bc_in

    !
    ! !LOCAL VARIABLES:

    integer  :: t            ! day of year
    integer  :: ncolddays    ! no days underneath the threshold for leaf drop
    integer  :: i
    integer  :: timesincedleafon,timesincedleafoff,timesinceleafon,timesinceleafoff
    integer  :: refdate
    integer  :: curdate
    
    integer  :: yr                       ! year (0, ...)
    integer  :: mon                      ! month (1, ..., 12)
    integer  :: day                      ! day of month (1, ..., 31)
    integer  :: sec                      ! seconds of the day

    real(r8) :: gdd_threshold
    integer  :: ncdstart     ! beginning of counting period for chilling degree days.
    integer  :: gddstart     ! beginning of counting period for growing degree days.
    real(r8) :: temp_in_C    ! daily averaged temperature in celcius

    real(r8), parameter :: canopy_leaf_lifespan = 365.0_r8    ! Mean lifespan canopy leaves
                                                              ! FIX(RGK 07/10/17)
                                                              ! This is a band-aid on unusual code
                                                              

    ! Parameter of drought decid leaf loss in mm in top layer...FIX(RF,032414) 
    ! - this is arbitrary and poorly understood. Needs work. ED_

    !Parameters: defaults from Botta et al. 2000 GCB,6 709-725 
    !Parameters, default from from SDGVM model of senesence

    t  = hlm_day_of_year
    temp_in_C = bc_in%t_veg24_si - tfrz

    !-----------------Cold Phenology--------------------!              

    !Zero growing degree and chilling day counters
    if (currentSite%lat > 0)then
       ncdstart = 270  !Northern Hemisphere begining November
       gddstart = 1    !Northern Hemisphere begining January
    else
       ncdstart = 120  !Southern Hemisphere beginning May
       gddstart = 181  !Northern Hemisphere begining July
    endif
    
    ! FIX(SPM,032414) - this will only work for the first year, no?
    if (t == ncdstart)then
       currentSite%ncd = 0._r8
    endif

    !Accumulate growing/chilling days after start of counting period
    if (temp_in_C  <  ED_val_phen_chiltemp)then
       currentSite%ncd = currentSite%ncd + 1.0_r8
    endif

    !GDD accumulation function, which also depends on chilling days.
    gdd_threshold = ED_val_phen_a + ED_val_phen_b*exp(ED_val_phen_c*currentSite%ncd)

    !Accumulate temperature of last 10 days.
    currentSite%last_n_days(2:senes) =  currentSite%last_n_days(1:senes-1)
    currentSite%last_n_days(1) = temp_in_C                                      
    !count number of days for leaves off
    ncolddays = 0
    do i = 1,senes
       if (currentSite%last_n_days(i) < ED_val_phen_coldtemp)then
          ncolddays = ncolddays + 1
       endif
    enddo

    ! Here is where we do the GDD accumulation calculation
    !
    ! reset GDD on set dates
    if (t == gddstart)then
       currentSite%ED_GDD_site = 0._r8
    endif
    !
    ! accumulate the GDD using daily mean temperatures
    if (bc_in%t_veg24_si .gt. tfrz) then
       currentSite%ED_GDD_site = currentSite%ED_GDD_site + bc_in%t_veg24_si - tfrz
    endif
    

    timesinceleafoff = hlm_model_day - currentSite%leafoffdate
    !LEAF ON: COLD DECIDUOUS. Needs to
    !1) have exceeded the growing degree day threshold 
    !2) The leaves should not be on already
    !3) There should have been at least on chilling day in the counting period.  
    if (currentSite%ED_GDD_site > gdd_threshold)then
       if (currentSite%status == 1) then
          if (currentSite%ncd >= 1) then
             currentSite%status = 2     !alter status of site to 'leaves on'
             ! NOTE(bja, 2015-01) should leafondate = model_day to be consistent with leaf off?
             currentSite%leafondate = t !record leaf on date   
             if ( DEBUG ) write(fates_log(),*) 'leaves on'
          endif !ncd
       endif !status
    endif !GDD

    timesinceleafon = hlm_model_day - currentSite%leafondate


    !LEAF OFF: COLD THRESHOLD
    !Needs to:
    !1) have exceeded the number of cold days threshold
    !2) have exceeded the minimum leafon time.
    !3) The leaves should not be off already
    !4) The day of the year should be larger than the counting period. (not sure if we need this/if it will break the restarting)
    
    if (ncolddays > ED_val_phen_ncolddayslim)then
     if (timesinceleafon > ED_val_phen_mindayson)then
       if (currentSite%status == 2)then
          currentSite%status = 1        !alter status of site to 'leaves on'
          currentSite%leafoffdate = hlm_model_day   !record leaf off date   
          if ( DEBUG ) write(fates_log(),*) 'leaves off'
       endif
    endif
    endif

    !LEAF OFF: COLD LIFESPAN THRESHOLD
    if(timesinceleafoff > 400)then !remove leaves after a whole year when there is no 'off' period.  
       if(currentSite%status == 2)then
          currentSite%status = 1        !alter status of site to 'leaves on'
          currentSite%leafoffdate = hlm_model_day   !record leaf off date   
          if ( DEBUG ) write(fates_log(),*) 'leaves off'
       endif
    endif

    !-----------------Drought Phenology--------------------!
    ! Principles of drought-deciduos phenology model...
    ! The 'dstatus' flag is 2 when leaves are on, and 1 when leaves area off. 
    ! The following sets those site-level flags, which are acted on in phenology_deciduos. 
    ! A* The leaves live for either the length of time the soil moisture is over the threshold 
    ! or the lifetime of the leaves, whichever is shorter. 
    ! B*: If the soil is only wet for a very short time, then the leaves stay on for 100 days
    ! C*: The leaves are only permitted to come ON for a 60 day window around when they last came on, 
    ! to prevent 'flickering' on in response to wet season storms
    ! D*: We don't allow anything to happen in the first ten days to allow the water memory window to come into equlibirum. 
    ! E*: If the soil is always wet, the leaves come on at the beginning of the window, and then last for their lifespan. 
    ! ISSUES
    ! 1. It's not clear what water content we should track. Here we are tracking the top layer, 
    ! but we probably should track something like BTRAN,
    ! but BTRAN is defined for each PFT, and there could potentially be more than one stress-dec PFT.... ?
    ! 2. In the beginning, the window is set at an arbitrary time of the year, so the leaves might come on 
    ! in the dry season, using up stored reserves
    ! for the stress-dec plants, and potentially killing them. To get around this, we need to read in the 
    ! 'leaf on' date from some kind of start-up file
    ! but we would need that to happen for every resolution, etc. 
    ! 3. Will this methodology properly kill off the stress-dec trees where there is no water stress? 
    ! What about where the wet period coincides with the
    ! warm period? We would just get them overlapping with the cold-dec trees, even though that isn't appropriate.... 
    ! Why don't the drought deciduous trees grow
    ! in the North? Is cold decidousness maybe even the same as drought deciduosness there (and so does this 
    ! distinction actually matter??).... 

    !Accumulate surface water memory of last 10 days.

    do i = 1,numWaterMem-1 !shift memory along one
       currentSite%water_memory(numWaterMem+1-i) = currentSite%water_memory(numWaterMem-i)
    enddo
    currentSite%water_memory(1) = bc_in%h2o_liqvol_gl(1)   !waterstate_inst%h2osoi_vol_col(coli,1)

    !In drought phenology, we often need to force the leaves to stay on or off as moisture fluctuates...     
    timesincedleafoff = 0
    if (currentSite%dstatus == 1)then !the leaves are off. How long have they been off? 
       !leaves have come on, but last year, so at a later date than now.
       if (currentSite%dleafoffdate > 0.and.currentSite%dleafoffdate > t)then 
          timesincedleafoff = t + (360 - currentSite%dleafoffdate)
       else
          timesincedleafoff = t - currentSite%dleafoffdate    
       endif
    endif

    timesincedleafon = 0
    !the leaves are on. How long have they been on? 
    if (currentSite%dstatus == 2)then  
       !leaves have come on, but last year, so at a later date than now.
       if (currentSite%dleafondate > 0.and.currentSite%dleafondate > t)then 
          timesincedleafon = t + (360 - currentSite%dleafondate)
       else
          timesincedleafon = t - currentSite%dleafondate      
       endif
    endif

    !LEAF ON: DROUGHT DECIDUOUS WETNESS
    !Here, we used a window of oppurtunity to determine if we are close to the time when then leaves came on last year
    if ((t >= currentSite%dleafondate - 30.and.t <= currentSite%dleafondate + 30).or.(t > 360 - 15.and. &
         currentSite%dleafondate < 15))then ! are we in the window?
       ! TODO: CHANGE THIS MATH, MOVE THE DENOMENATOR OUTSIDE OF THE SUM (rgk 01-2017)
       if (sum(currentSite%water_memory(1:numWaterMem)/dble(numWaterMem)) &
            >= ED_val_phen_drought_threshold.and.currentSite%dstatus == 1.and.t >= 10)then 
          ! leave some minimum time between leaf off and leaf on to prevent 'flickering'.  
          if (timesincedleafoff > ED_val_phen_doff_time)then  
             currentSite%dstatus = 2     !alter status of site to 'leaves on'
             currentSite%dleafondate = t   !record leaf on date
          endif
       endif
    endif

   !we still haven't done budburst by end of window
    if (t == currentSite%dleafondate+30.and.currentSite%dstatus == 1)then 
       currentSite%dstatus = 2    ! force budburst!
       currentSite%dleafondate = t   ! record leaf on date
    endif

    !LEAF OFF: DROUGHT DECIDUOUS LIFESPAN - if the leaf gets to the end of its useful life. A*, E*
    if (currentSite%dstatus == 2.and.t >= 10)then  !D*
       !Are the leaves at the end of their lives? 
       !FIX(RF,0401014)- this is hardwiring....
       !FIX(RGK:changed from hard-coded pft 7 leaf lifespan to labeled constant (1 year)
       if ( timesincedleafon > canopy_leaf_lifespan )then 
          currentSite%dstatus = 1         !alter status of site to 'leaves on'
          currentSite%dleafoffdate = t    !record leaf on date          
       endif
    endif

    !LEAF OFF: DROUGHT DECIDUOUS DRYNESS - if the soil gets too dry, and the leaves have already been on a while... 
    if (currentSite%dstatus == 2.and.t >= 10)then  !D*
       if (sum(currentSite%water_memory(1:10)/10._r8) <= ED_val_phen_drought_threshold)then 
          if (timesincedleafon > 100)then !B* Have the leaves been on for some reasonable length of time? To prevent flickering. 
             currentSite%dstatus = 1      !alter status of site to 'leaves on'
             currentSite%dleafoffdate = t !record leaf on date           
          endif
       endif
    endif

    call phenology_leafonoff(currentSite)

  end subroutine phenology

  ! ============================================================================
  subroutine phenology_leafonoff(currentSite)
    !
    ! !DESCRIPTION:
    ! Controls the leaf on and off economics
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(ed_site_type), intent(inout), target :: currentSite
    !
    ! !LOCAL VARIABLES:
    type(ed_patch_type) , pointer :: currentPatch     
    type(ed_cohort_type), pointer :: currentCohort  

    real(r8)           :: store_output ! the amount of the store to put into leaves -
                                       ! is a barrier against negative storage and C starvation. 

    !------------------------------------------------------------------------

    currentPatch => CurrentSite%oldest_patch   

    store_output  = 0.5_r8

    do while(associated(currentPatch))    
       currentCohort => currentPatch%tallest
       do while(associated(currentCohort))        
                
          !COLD LEAF ON
          if (EDPftvarcon_inst%season_decid(currentCohort%pft) == 1)then
             if (currentSite%status == 2)then !we have just moved to leaves being on . 
                if (currentCohort%status_coh == 1)then !Are the leaves currently off?        
                   currentCohort%status_coh = 2    !Leaves are on, so change status to stop flow of carbon out of bstore. 
                   if (currentCohort%laimemory <= currentCohort%bstore)then
                      currentCohort%bl = currentCohort%laimemory !extract stored carbon to make new leaves.
                   else
                      ! we can only put on as much carbon as there is in the store...
                      ! nb. Putting all of bstore into leaves is C-starvation suicidal. 
                      ! The tendency for this could be parameterized
                      currentCohort%bl = currentCohort%bstore * store_output
                   endif

                   ! Add deployed carbon to alive biomass pool
                   currentCohort%balive = currentCohort%balive + currentCohort%bl

                   if ( DEBUG ) write(fates_log(),*) 'EDPhysMod 1 ',currentCohort%bstore

                   currentCohort%bstore = currentCohort%bstore - currentCohort%bl  ! Drain store

                   if ( DEBUG ) write(fates_log(),*) 'EDPhysMod 2 ',currentCohort%bstore

                   currentCohort%laimemory = 0.0_r8

                endif !pft phenology
             endif ! growing season 

             !COLD LEAF OFF
             currentCohort%leaf_litter = 0.0_r8 !zero leaf litter for today. 
             if (currentSite%status == 1)then !past leaf drop day? Leaves still on tree?  
                if (currentCohort%status_coh == 2)then ! leaves have not dropped
                   currentCohort%status_coh      = 1                  
                   !remember what the lai was this year to put the same amount back on in the spring... 
                   currentCohort%laimemory   = currentCohort%bl  
                   ! decrement balive for leaf litterfall        
                   currentCohort%balive      = currentCohort%balive - currentCohort%bl 
                   ! add lost carbon to litter
                   currentCohort%leaf_litter = currentCohort%bl 
                   currentCohort%bl          = 0.0_r8                            
                endif !leaf status
             endif !currentSite status
          endif  !season_decid

          !DROUGHT LEAF ON
          if (EDPftvarcon_inst%stress_decid(currentCohort%pft) == 1)then
             if (currentSite%dstatus == 2)then !we have just moved to leaves being on . 
                if (currentCohort%status_coh == 1)then !is it the leaf-on day? Are the leaves currently off?       
                   currentCohort%status_coh = 2    !Leaves are on, so change status to stop flow of carbon out of bstore. 
                   if (currentCohort%laimemory <= currentCohort%bstore)then
                      currentCohort%bl = currentCohort%laimemory !extract stored carbon to make new leaves.
                   else
                    currentCohort%bl = currentCohort%bstore * store_output    !we can only put on as much carbon as there is in the store...
                    endif
                   currentCohort%balive = currentCohort%balive + currentCohort%bl

                   if ( DEBUG ) write(fates_log(),*) 'EDPhysMod 3 ',currentCohort%bstore

                   currentCohort%bstore = currentCohort%bstore - currentCohort%bl ! empty store

                   if ( DEBUG ) write(fates_log(),*) 'EDPhysMod 4 ',currentCohort%bstore

                   currentCohort%laimemory = 0.0_r8

                endif !currentCohort status again?
             endif   !currentSite status

             !DROUGHT LEAF OFF
             if (currentSite%dstatus == 1)then        
                if (currentCohort%status_coh == 2)then ! leaves have not dropped
                   currentCohort%status_coh      = 1   
                   currentCohort%laimemory   = currentCohort%bl
                   ! decrement balive for leaf litterfall  
                   currentCohort%balive      = currentCohort%balive - currentCohort%bl   
                   ! add retranslocated carbon (very small) to store.      
                   currentCohort%bstore      = currentCohort%bstore        
                   ! add falling leaves to litter pools . convert to KgC/m2                    
                   currentCohort%leaf_litter = currentCohort%bl  
                   currentCohort%bl          = 0.0_r8                                        
                endif
             endif !status
          endif !drought dec.
          currentCohort => currentCohort%shorter
       enddo !currentCohort

       currentPatch => currentPatch%younger

    enddo !currentPatch

  end subroutine phenology_leafonoff


  ! ============================================================================
  subroutine seeds_in( currentSite, cp_pnt )
    !
    ! !DESCRIPTION:
    !  Flux from plants into seed pool. 
    !
    ! !USES:
    use EDTypesMod, only : AREA
    use EDTypesMod, only : homogenize_seed_pfts
    !
    ! !ARGUMENTS    
    type(ed_site_type), intent(inout), target  :: currentSite
    type(ed_patch_type), intent(inout), target :: cp_pnt ! seeds go to these patches.
    !
    ! !LOCAL VARIABLES:
    type(ed_patch_type),  pointer :: currentPatch
    type(ed_cohort_type), pointer :: currentCohort
    integer :: p
    logical :: pft_present(maxpft)
    real(r8) :: npfts_present
    !----------------------------------------------------------------------

    currentPatch => cp_pnt
   
    currentPatch%seeds_in(:) = 0.0_r8

    if ( homogenize_seed_pfts ) then
       ! special mode to remove intergenerational filters on PFT existence: each PFT seeds all PFTs
       ! first loop over all patches and cohorts to see what and how many PFTs are present on this site
       pft_present(:) = .false.
       npfts_present =  0._r8
       currentPatch => currentSite%oldest_patch
       do while(associated(currentPatch))
          currentCohort => currentPatch%tallest
          do while (associated(currentCohort))
             p = currentCohort%pft
             if (.not. pft_present(p)) then
                pft_present(p) = .true.
                npfts_present = npfts_present + 1._r8
             endif
             currentCohort => currentCohort%shorter
          enddo !cohort loop                        
          currentPatch => currentPatch%younger
       enddo ! patch loop
       
       ! now calculate the homogenized seed flux into each PFT pool
       currentPatch => cp_pnt
       currentCohort => currentPatch%tallest
       do while (associated(currentCohort))
          do p = 1, numpft
             if (pft_present(p)) then
                currentPatch%seeds_in(p) = currentPatch%seeds_in(p) +  currentCohort%seed_prod * currentCohort%n / &
                     (currentPatch%area * npfts_present)
             endif
          end do
          currentCohort => currentCohort%shorter
       enddo !cohort loop                  
    else

    ! normal case: each PFT seeds its own type
    currentCohort => currentPatch%tallest
    do while (associated(currentCohort))
       p = currentCohort%pft
       currentPatch%seeds_in(p) = currentPatch%seeds_in(p) +  &
             currentCohort%seed_prod * currentCohort%n/currentPatch%area
       currentCohort => currentCohort%shorter
    enddo !cohort loop

    endif

    currentPatch => currentSite%oldest_patch

    do while(associated(currentPatch))
       if (external_recruitment == 1) then !external seed rain - needed to prevent extinction  
          do p = 1,numpft
           currentPatch%seeds_in(p) = currentPatch%seeds_in(p) + &
                 EDPftvarcon_inst%seed_rain(p) !KgC/m2/year
           currentSite%seed_rain_flux(p) = currentSite%seed_rain_flux(p) + &
                 EDPftvarcon_inst%seed_rain(p) * currentPatch%area/AREA !KgC/m2/year
          enddo
       endif
       currentPatch => currentPatch%younger
    enddo

  end subroutine seeds_in
  
  ! ============================================================================
  subroutine seed_decay( currentSite, currentPatch )
    !
    ! !DESCRIPTION:
    !  Flux from seed pool into leaf litter pool    
    !
    ! !USES:
    use EDPftvarcon       , only : EDPftvarcon_inst
    !
    ! !ARGUMENTS    
    type(ed_site_type), intent(inout), target  :: currentSite
    type(ed_patch_type),intent(inout) :: currentPatch ! seeds go to these patches.
    !
    ! !LOCAL VARIABLES:
    integer  ::  p
    !----------------------------------------------------------------------

    ! default value from Liscke and Loffler 2006 ; making this a PFT-specific parameter
    ! decays the seed pool according to exponential model
    ! seed_decay_turnover is in yr-1
    do p = 1,numpft 
       currentPatch%seed_decay(p) =  currentSite%seed_bank(p) * EDPftvarcon_inst%seed_decay_turnover(p)
    enddo
 
  end subroutine seed_decay

  ! ============================================================================
  subroutine seed_germination( currentSite, currentPatch ) 
    !
    ! !DESCRIPTION:
    !  Flux from seed pool into sapling pool    
    !
    ! !USES:
    use EDPftvarcon       , only : EDPftvarcon_inst
    !
    ! !ARGUMENTS    
    type(ed_site_type), intent(inout), target  :: currentSite
    type(ed_patch_type),intent(inout) :: currentPatch ! seeds go to these patches.
    !
    ! !LOCAL VARIABLES:
    integer :: p
    real(r8) max_germination !cap on germination rates. KgC/m2/yr Lishcke et al. 2009
    !----------------------------------------------------------------------

    max_germination = 1.0_r8 !this is arbitrary

    ! germination_timescale is being pulled to PFT parameter; units are 1/yr
    ! thus the mortality rate of seed -> recruit (in units of carbon) is seed_decay_turnover(p)/germination_timescale(p)
    ! and thus the mortlaity rate (in units of individuals) is the product of that times the ratio of (hypothetical) seed mass to recruit biomass
    do p = 1,numpft
       currentPatch%seed_germination(p) =  min(currentSite%seed_bank(p) * &
             EDPftvarcon_inst%germination_timescale(p),max_germination)
    enddo

  end subroutine seed_germination

  ! ============================================================================
  subroutine Growth_Derivatives( currentSite, currentCohort, bc_in)
    !
    ! !DESCRIPTION:
    !  Main subroutine controlling growth and allocation derivatives    
    !
    ! !USES:
    use EDGrowthFunctionsMod , only : Bleaf, dDbhdBd, dhdbd, hite, mortality_rates,dDbhdBl
    use FatesInterfaceMod, only : hlm_use_ed_prescribed_phys
    use EDLoggingMortalityMod, only : LoggingMortality_frac

    !
    ! !ARGUMENTS    
    type(ed_site_type), intent(inout), target  :: currentSite
    type(ed_cohort_type),intent(inout), target :: currentCohort
    type(bc_in_type), intent(in)               :: bc_in
    !
    ! !LOCAL VARIABLES:
    real(r8) :: dbldbd   !rate of change of dead biomass per unit dbh 
    real(r8) :: dbrdbd   !rate of change of root biomass per unit dbh
    real(r8) :: dbswdbd  !rate of change of sapwood biomass per unit dbh
    real(r8) :: dhdbd_fn !rate of change of height per unit dbh
    real(r8) :: va       !fraction of growth going to alive biomass
    real(r8) :: vs       !fraction of growth going to structural biomass
    real(r8) :: u,h      !intermediates 
    real(r8) :: frac     !fraction the stored carbon is of target store amount
    real(r8) :: f_store  !fraction of NPP allocated to storage in this timestep (functionf of stored pool)
    real(r8) :: gr_fract !fraction of carbon balance that is allocated to growth (not reproduction)
    real(r8) :: target_balive  !target leaf biomass under allometric optimum.  
    real(r8) :: cmort    ! starvation mortality rate (fraction per year)
    real(r8) :: bmort    ! background mortality rate (fraction per year)
    real(r8) :: hmort    ! hydraulic failure mortality rate (fraction per year)

    real(r8) :: lmort_logging     ! Mortality fraction associated with direct logging
    real(r8) :: lmort_collateral  ! Mortality fraction associated with logging collateral damage
    real(r8) :: lmort_infra       ! Mortality fraction associated with logging infrastructure
    real(r8) :: dndt_logging      ! Mortality rate (per day) associated with the a logging event
    
    real(r8) :: balive_loss
    !----------------------------------------------------------------------

    ! Mortality for trees in the understorey. 
    !if trees are in the canopy, then their death is 'disturbance'. This probably needs a different terminology
    call mortality_rates(currentCohort,cmort,hmort,bmort)
    call LoggingMortality_frac(currentCohort%pft, currentCohort%dbh, &
                               currentCohort%lmort_logging,                       &
                               currentCohort%lmort_collateral,                    &
                               currentCohort%lmort_infra )

    if (currentCohort%canopy_layer > 1)then 
       
       ! Include understory logging mortality rates not associated with disturbance
       dndt_logging = (currentCohort%lmort_logging    + &
                       currentCohort%lmort_collateral + &
                       currentCohort%lmort_infra)/hlm_freq_day

       currentCohort%dndt = -1.0_r8 * (cmort+hmort+bmort+dndt_logging) * currentCohort%n
    else
       currentCohort%dndt = -(1.0_r8 - fates_mortality_disturbance_fraction) &
            * (cmort+hmort+bmort) * currentCohort%n
    endif

    ! Height
    currentCohort%hite = Hite(currentCohort) 
    h = currentCohort%hite
                       
    call allocate_live_biomass(currentCohort,0)

   ! calculate target size of living biomass compartment for a given dbh.   
    target_balive = Bleaf(currentCohort) * (1.0_r8 + EDPftvarcon_inst%allom_l2fr(currentCohort%pft) + &
         EDpftvarcon_inst%allom_latosa_int(currentCohort%pft)*h)
    !target balive without leaves. 
    if (currentCohort%status_coh == 1)then 
       target_balive = Bleaf(currentCohort) * (EDPftvarcon_inst%allom_l2fr(currentCohort%pft) + &
            EDpftvarcon_inst%allom_latosa_int(currentCohort%pft) * h)
    endif

    ! NPP 
    if ( DEBUG ) write(fates_log(),*) 'EDphys 716 ',currentCohort%npp_acc

    ! convert from kgC/indiv/day into kgC/indiv/year
    ! TODO: CONVERT DAYS_PER_YEAR TO DBLE (HOLDING FOR B4B COMPARISONS, RGK-01-2017)
    currentCohort%npp_acc_hold  = currentCohort%npp_acc  * hlm_days_per_year 
    currentCohort%gpp_acc_hold  = currentCohort%gpp_acc  * hlm_days_per_year
    currentCohort%resp_acc_hold = currentCohort%resp_acc * hlm_days_per_year

    if (hlm_use_ed_prescribed_phys .eq. itrue) then
       if (currentCohort%canopy_layer .eq. 1) then
          currentCohort%npp_acc_hold = EDPftvarcon_inst%prescribed_npp_canopy(currentCohort%pft) * currentCohort%c_area / currentCohort%n
       else
          currentCohort%npp_acc_hold = EDPftvarcon_inst%prescribed_npp_understory(currentCohort%pft) * currentCohort%c_area / currentCohort%n
       endif
    endif

    currentSite%flux_in = currentSite%flux_in + currentCohort%npp_acc * currentCohort%n

    ! Maintenance demands     
    if (EDPftvarcon_inst%evergreen(currentCohort%pft) == 1)then !grass and EBT
       currentCohort%leaf_md = currentCohort%bl / EDPftvarcon_inst%leaf_long(currentCohort%pft)
       currentCohort%root_md = currentCohort%br / EDPftvarcon_inst%root_long(currentCohort%pft)
       currentCohort%md      = currentCohort%root_md + currentCohort%leaf_md
    endif

    !FIX(RF,032414) - I took out the stem turnover demand as it seemed excesively high and caused odd size-reated 
    ! decline affect
    !with which I am not especially comfortable, particularly as the concept of sapwood turnover is unclear for trees that 
    !are still in an expansion phase. 

    if (EDPftvarcon_inst%season_decid(currentCohort%pft) == 1)then 
       currentCohort%root_md = currentCohort%br /EDPftvarcon_inst%root_long(currentCohort%pft)
       currentCohort%leaf_md = 0._r8
       currentCohort%md = currentCohort%root_md + currentCohort%leaf_md
    endif

    if (EDPftvarcon_inst%stress_decid(currentCohort%pft) == 1)then 
       currentCohort%root_md = currentCohort%br /EDPftvarcon_inst%root_long(currentCohort%pft)
       currentCohort%leaf_md = 0._r8
       currentCohort%md = currentCohort%root_md + currentCohort%leaf_md
    endif

    if (EDPftvarcon_inst%stress_decid(currentCohort%pft) /= 1 &
          .and.EDPftvarcon_inst%season_decid(currentCohort%pft) /= 1.and. &
         EDPftvarcon_inst%evergreen(currentCohort%pft) /= 1)then
       write(fates_log(),*) 'problem with phenology definitions',currentCohort%pft, &
            EDPftvarcon_inst%stress_decid(currentCohort%pft), &
            EDPftvarcon_inst%season_decid(currentCohort%pft),EDPftvarcon_inst%evergreen(currentCohort%pft)
    endif

    ! FIX(RF,032414) -turned off for now as it makes balive go negative....
    ! FIX(RF,032414) jan2012 0.01_r8 * currentCohort%bdead
    currentCohort%woody_turnover = 0.0_r8
    currentCohort%md = currentCohort%md + currentCohort%woody_turnover

    ! Calculate carbon balance 
    ! this is the fraction of maintenance demand we -have- to do...

    if ( DEBUG ) write(fates_log(),*) 'EDphys 760 ',currentCohort%npp_acc_hold, currentCohort%md, &
                   EDPftvarcon_inst%leaf_stor_priority(currentCohort%pft)

    currentCohort%carbon_balance = currentCohort%npp_acc_hold - &
          currentCohort%md * EDPftvarcon_inst%leaf_stor_priority(currentCohort%pft)

    ! Allowing only carbon from NPP pool to account for npp flux into the maintenance turnover pools
    ! ie this does not include any use of storage carbon or balive to make up for missing carbon balance in the transfer
    currentCohort%npp_leaf  = max(0.0_r8,min(currentCohort%npp_acc_hold*currentCohort%leaf_md/currentCohort%md, &
                                  currentCohort%leaf_md*EDPftvarcon_inst%leaf_stor_priority(currentCohort%pft)))
    currentCohort%npp_froot = max(0.0_r8,min(currentCohort%npp_acc_hold*currentCohort%root_md/currentCohort%md, &
                                  currentCohort%root_md*EDPftvarcon_inst%leaf_stor_priority(currentCohort%pft)))

    if (Bleaf(currentCohort) > 0._r8)then

       if ( DEBUG ) write(fates_log(),*) 'EDphys A ',currentCohort%carbon_balance

       if (currentCohort%carbon_balance > 0._r8)then !spend C on growing and storing

          !what fraction of the target storage do we have? 
          frac = max(0.0_r8,currentCohort%bstore/(Bleaf(currentCohort) * EDPftvarcon_inst%cushion(currentCohort%pft)))
          ! FIX(SPM,080514,fstore never used ) 
          f_store = max(exp(-1.*frac**4._r8) - exp( -1.0_r8 ),0.0_r8)  
          !what fraction of allocation do we divert to storage?
          !what is the flux into the store?
          currentCohort%storage_flux = currentCohort%carbon_balance * f_store

          currentCohort%npp_store = currentCohort%carbon_balance * f_store         
          if ( DEBUG ) write(fates_log(),*) 'EDphys B ',f_store

          !what is the tax on the carbon available for growth? 
          currentCohort%carbon_balance = currentCohort%carbon_balance * (1.0_r8 - f_store)  
       else  !cbalance is negative. Take C out of store to pay for maintenance respn.

          currentCohort%storage_flux = currentCohort%carbon_balance 

          ! Note that npp_store only tracks the flux between NPP and storage.  Storage can 
          ! also be drawn down to support some turnover demand.
          currentCohort%npp_store = min(0.0_r8,currentCohort%npp_acc_hold)

          currentCohort%carbon_balance = 0._r8 
       endif

    else

       write(fates_log(),*) 'No target leaf area in GrowthDerivs? Bleaf(cohort) <= 0?'
       call endrun(msg=errMsg(sourcefile, __LINE__))

    endif

    !Do we have enough carbon left over to make up the rest of the turnover demand? 
    balive_loss = 0._r8
    if (currentCohort%carbon_balance > currentCohort%md*(1.0_r8- EDPftvarcon_inst%leaf_stor_priority(currentCohort%pft)))then ! Yes...
       currentCohort%carbon_balance = currentCohort%carbon_balance - currentCohort%md * (1.0_r8 - &
             EDPftvarcon_inst%leaf_stor_priority(currentCohort%pft))

       currentCohort%npp_leaf  = currentCohort%npp_leaf  + &
            currentCohort%leaf_md *  (1.0_r8-EDPftvarcon_inst%leaf_stor_priority(currentCohort%pft))
       currentCohort%npp_froot = currentCohort%npp_froot + &
            currentCohort%root_md *  (1.0_r8-EDPftvarcon_inst%leaf_stor_priority(currentCohort%pft))

    else ! we can't maintain constant leaf area and root area. Balive is reduced

       currentCohort%npp_leaf  = currentCohort%npp_leaf  + &
             max(0.0_r8,currentCohort%carbon_balance*(currentCohort%leaf_md/currentCohort%md))
       currentCohort%npp_froot = currentCohort%npp_froot + &
             max(0.0_r8,currentCohort%carbon_balance*(currentCohort%root_md/currentCohort%md))

       balive_loss = currentCohort%md *(1.0_r8- EDPftvarcon_inst%leaf_stor_priority(currentCohort%pft))- currentCohort%carbon_balance
       currentCohort%carbon_balance = 0._r8
    endif

    !********************************************/
    ! Allometry & allocation of remaining carbon*/
    !********************************************/
    !Use remaining carbon to refill balive or to get larger. 

    !only if carbon balance is +ve
    if ((currentCohort%balive >= target_balive).AND.(currentCohort%carbon_balance >  0._r8))then 
       ! fraction of carbon going into active vs structural carbon        
       if (currentCohort%dbh <= EDPftvarcon_inst%allom_dbh_maxheight(currentCohort%pft))then ! cap on leaf biomass
          dbldbd = dDbhdBd(currentCohort)/dDbhdBl(currentCohort) 
          dbrdbd = EDPftvarcon_inst%allom_l2fr(currentCohort%pft) * dbldbd
          dhdbd_fn = dhdbd(currentCohort)
          dbswdbd = EDpftvarcon_inst%allom_latosa_int(currentCohort%pft) * (h*dbldbd + currentCohort%bl*dhdbd_fn)
          u  = 1.0_r8 / (dbldbd + dbrdbd + dbswdbd)     
          va = 1.0_r8 / (1.0_r8 + u)
          vs = u / (1.0_r8 + u)
          gr_fract = 1.0_r8 - EDPftvarcon_inst%seed_alloc(currentCohort%pft)
       else
          dbldbd = 0._r8; dbrdbd = 0._r8 ;dbswdbd = 0._r8      
          va = 0.0_r8
          vs = 1.0_r8
          gr_fract = 1.0_r8 - (EDPftvarcon_inst%seed_alloc(currentCohort%pft) + EDPftvarcon_inst%clone_alloc(currentCohort%pft))
       endif

       !FIX(RF,032414) - to fix high bl's. needed to prevent numerical errors without the ODEINT.  
       if (currentCohort%balive > target_balive*1.1_r8)then  
          va = 0.0_r8; vs = 1._r8
          if (DEBUG) write(fates_log(),*) 'using high bl cap',target_balive,currentCohort%balive                        
       endif

    else         
       dbldbd = 0._r8; dbrdbd = 0._r8; dbswdbd = 0._r8
       va = 1.0_r8; vs = 0._r8                        
       gr_fract = 1.0_r8
    endif

    ! calculate derivatives of living and dead carbon pools  
    currentCohort%dbalivedt = gr_fract * va * currentCohort%carbon_balance - balive_loss
    currentCohort%dbdeaddt  = gr_fract * vs * currentCohort%carbon_balance
    currentCohort%dbstoredt = currentCohort%storage_flux

    if ( DEBUG ) write(fates_log(),*) 'EDPhys dbstoredt I ',currentCohort%dbstoredt

    currentCohort%seed_prod = (1.0_r8 - gr_fract) * currentCohort%carbon_balance

    if (abs(currentCohort%npp_acc_hold-(currentCohort%dbalivedt+currentCohort%dbdeaddt+currentCohort%dbstoredt+ &
         currentCohort%seed_prod+currentCohort%md)) > 0.0000000001_r8)then
       write(fates_log(),*) 'error in carbon check growth derivs',currentCohort%npp_acc_hold- &
            (currentCohort%dbalivedt+currentCohort%dbdeaddt+currentCohort%dbstoredt+currentCohort%seed_prod+currentCohort%md)
       write(fates_log(),*) 'cohort fluxes',currentCohort%pft,currentCohort%canopy_layer,currentCohort%n, &
            currentCohort%npp_acc_hold,currentCohort%dbalivedt,balive_loss, &
            currentCohort%dbdeaddt,currentCohort%dbstoredt,currentCohort%seed_prod,currentCohort%md * &
            EDPftvarcon_inst%leaf_stor_priority(currentCohort%pft)
       write(fates_log(),*) 'proxies' ,target_balive,currentCohort%balive,currentCohort%dbh,va,vs,gr_fract
    endif

    ! prevent negative leaf pool (but not negative store pool). This is also a numerical error prevention, 
    ! but it shouldn't happen actually... 
    if (-1.0_r8*currentCohort%dbalivedt * hlm_freq_day > currentCohort%balive*0.99)then 
       write(fates_log(),*) 'using non-neg leaf mass cap',currentCohort%balive , currentCohort%dbalivedt,currentCohort%dbstoredt, &
            currentCohort%carbon_balance
       currentCohort%dbstoredt = currentCohort%dbstoredt + currentCohort%dbalivedt

       if ( DEBUG ) write(fates_log(),*) 'EDPhys dbstoredt II ',currentCohort%dbstoredt

       currentCohort%dbalivedt = 0._r8 
    endif

    currentCohort%npp_bseed = currentCohort%seed_prod

    ! calculate change in diameter and height 
    currentCohort%ddbhdt = currentCohort%dbdeaddt * dDbhdBd(currentCohort)
    currentCohort%dhdt   = currentCohort%dbdeaddt * dHdBd(currentCohort)

    ! If the cohort has grown, it is not new
    currentCohort%isnew=.false.

  end subroutine Growth_Derivatives

  ! ============================================================================
  subroutine recruitment( currentSite, currentPatch, bc_in )
    !
    ! !DESCRIPTION:
    ! spawn new cohorts of juveniles of each PFT             
    !
    ! !USES:
    use EDGrowthFunctionsMod, only : bdead,dbh, Bleaf
    use FatesInterfaceMod, only : hlm_use_ed_prescribed_phys
    !
    ! !ARGUMENTS    
    type(ed_site_type), intent(inout), target  :: currentSite
    type(ed_patch_type), intent(inout), pointer :: currentPatch
    type(bc_in_type), intent(in)                :: bc_in
    !
    ! !LOCAL VARIABLES:
    integer :: ft
    type (ed_cohort_type) , pointer :: temp_cohort
    integer :: cohortstatus
    !----------------------------------------------------------------------

    allocate(temp_cohort) ! create temporary cohort
    call zero_cohort(temp_cohort)

    do ft = 1,numpft

       temp_cohort%canopy_trim = 0.8_r8  !starting with the canopy not fully expanded 
       temp_cohort%pft         = ft
       temp_cohort%hite        = EDPftvarcon_inst%hgt_min(ft)
       temp_cohort%dbh         = Dbh(temp_cohort)
       temp_cohort%bdead       = Bdead(temp_cohort)
       temp_cohort%balive      = Bleaf(temp_cohort)*(1.0_r8 + EDPftvarcon_inst%allom_l2fr(ft) &
            + EDpftvarcon_inst%allom_latosa_int(ft)*temp_cohort%hite)
       temp_cohort%bstore      = EDPftvarcon_inst%cushion(ft)*(temp_cohort%balive/ (1.0_r8 + EDPftvarcon_inst%allom_l2fr(ft) &
            + EDpftvarcon_inst%allom_latosa_int(ft)*temp_cohort%hite))

       if (hlm_use_ed_prescribed_phys .eq. ifalse) then
          temp_cohort%n           = currentPatch%area * currentPatch%seed_germination(ft)*hlm_freq_day &
               / (temp_cohort%bdead+temp_cohort%balive+temp_cohort%bstore)
       else
          ! prescribed recruitment rates. number per sq. meter per year
          temp_cohort%n        = currentPatch%area * EDPftvarcon_inst%prescribed_recruitment(ft) * hlm_freq_day
       endif

       temp_cohort%laimemory = 0.0_r8     
       if (EDPftvarcon_inst%season_decid(temp_cohort%pft) == 1.and.currentSite%status == 1)then
         temp_cohort%laimemory = (1.0_r8/(1.0_r8 + EDPftvarcon_inst%allom_l2fr(ft) + &
              EDpftvarcon_inst%allom_latosa_int(ft)*temp_cohort%hite))*temp_cohort%balive
       endif
       if (EDPftvarcon_inst%stress_decid(temp_cohort%pft) == 1.and.currentSite%dstatus == 1)then
         temp_cohort%laimemory = (1.0_r8/(1.0_r8 + EDPftvarcon_inst%allom_l2fr(ft) + &
            EDpftvarcon_inst%allom_latosa_int(ft)*temp_cohort%hite))*temp_cohort%balive
       endif

       cohortstatus = currentSite%status
       if (EDPftvarcon_inst%stress_decid(ft) == 1)then !drought decidous, override status. 
          cohortstatus = currentSite%dstatus
       endif

       if (temp_cohort%n > 0.0_r8 )then
           if ( DEBUG ) write(fates_log(),*) 'EDPhysiologyMod.F90 call create_cohort '
           call create_cohort(currentPatch, temp_cohort%pft, temp_cohort%n, temp_cohort%hite, temp_cohort%dbh, &
                temp_cohort%balive, temp_cohort%bdead, temp_cohort%bstore,  &
                temp_cohort%laimemory, cohortstatus, temp_cohort%canopy_trim, currentPatch%NCL_p, &
                bc_in)

           ! keep track of how many individuals were recruited for passing to history
           currentSite%recruitment_rate(ft) = currentSite%recruitment_rate(ft) + temp_cohort%n

       endif

    enddo  !pft loop

    deallocate(temp_cohort) ! delete temporary cohort

  end subroutine recruitment

  ! ============================================================================
  subroutine CWD_Input( currentSite, currentPatch)
    !
    ! !DESCRIPTION:
    ! Generate litter fields from turnover.  
    !
    ! !USES:
    use SFParamsMod , only : SF_val_CWD_frac

    !
    ! !ARGUMENTS    
    type(ed_site_type), intent(inout), target :: currentSite
    type(ed_patch_type),intent(inout), target :: currentPatch
    !
    ! !LOCAL VARIABLES:
    type(ed_cohort_type), pointer :: currentCohort
    integer  :: c,p
    real(r8) :: dead_n          ! total understorey dead tree density
    real(r8) :: dead_n_dlogging ! direct logging understory dead-tree density
    real(r8) :: dead_n_ilogging ! indirect understory dead-tree density (logging)
    real(r8) :: dead_n_natural  ! understory dead density not associated
                                ! with direct logging
    real(r8) :: trunk_product   ! carbon flux into trunk products kgC/day/site
    integer  :: pft
    !----------------------------------------------------------------------

    ! ================================================        
    ! Other direct litter fluxes happen in phenology and in spawn_patches. 
    ! ================================================   

    currentCohort => currentPatch%shortest

    do while(associated(currentCohort))
      pft = currentCohort%pft        
      ! ================================================        
      ! Litter from tissue turnover. KgC/m2/year
      ! ================================================   
      currentPatch%leaf_litter_in(pft) = currentPatch%leaf_litter_in(pft) + &
               currentCohort%leaf_md * currentCohort%n/currentPatch%area !turnover

      currentPatch%root_litter_in(pft) = currentPatch%root_litter_in(pft) + &
               currentCohort%root_md * currentCohort%n/currentPatch%area !turnover
      currentPatch%leaf_litter_in(pft) = currentPatch%leaf_litter_in(pft) + &
         currentCohort%leaf_litter * currentCohort%n/currentPatch%area/hlm_freq_day

      !daily leaf loss needs to be scaled up to the annual scale here. 
      
      do c = 1,ncwd
         currentPatch%cwd_AG_in(c) = currentPatch%cwd_AG_in(c) + currentCohort%woody_turnover * &
              SF_val_CWD_frac(c) * currentCohort%n/currentPatch%area *EDPftvarcon_inst%allom_agb_frac(currentCohort%pft)
         currentPatch%cwd_BG_in(c) = currentPatch%cwd_BG_in(c) + currentCohort%woody_turnover * &
              SF_val_CWD_frac(c) * currentCohort%n/currentPatch%area *(1.0_r8-EDPftvarcon_inst%allom_agb_frac(currentCohort%pft))
      enddo

      if (currentCohort%canopy_layer > 1)then   

          ! ================================================        
          ! Litter fluxes for understorey  mortality. KgC/m2/year
          ! ================================================

          ! Total number of dead understory (n/m2)
          dead_n = -1.0_r8 * currentCohort%dndt / currentPatch%area

          ! Total number of dead understory from direct logging
          ! (it is possible that large harvestable trees are in the understory)
          dead_n_dlogging = ( currentCohort%lmort_logging) * &
                currentCohort%n/hlm_freq_day/currentPatch%area
          
          ! Total number of dead understory from indirect logging
          dead_n_ilogging = ( currentCohort%lmort_collateral + currentCohort%lmort_infra) * &
                currentCohort%n/hlm_freq_day/currentPatch%area
          
          dead_n_natural = dead_n - dead_n_dlogging - dead_n_ilogging


          currentPatch%leaf_litter_in(pft) = currentPatch%leaf_litter_in(pft) + &
               currentCohort%bl * dead_n          
          currentPatch%root_litter_in(pft) = currentPatch%root_litter_in(pft) + &
               (currentCohort%br+currentCohort%bstore)     * dead_n

          ! Update diagnostics that track resource management
          currentSite%resources_management%delta_litter_stock  = &
                currentSite%resources_management%delta_litter_stock + &
                (currentCohort%bl+currentCohort%br+currentCohort%bstore) * &
                (dead_n_ilogging+dead_n_dlogging) * & 
                hlm_freq_day * currentPatch%area
          ! Update diagnostics that track resource management
          currentSite%resources_management%delta_biomass_stock = &
                currentSite%resources_management%delta_biomass_stock + &
                (currentCohort%bl+currentCohort%br+currentCohort%bstore) * &
                (dead_n_ilogging+dead_n_dlogging) * & 
                hlm_freq_day * currentPatch%area

          do c = 1,ncwd
             
             currentPatch%cwd_BG_in(c) = currentPatch%cwd_BG_in(c) + (currentCohort%bdead+currentCohort%bsw) * &
                   SF_val_CWD_frac(c) * dead_n * (1.0_r8-EDPftvarcon_inst%allom_agb_frac(currentCohort%pft))

             ! Send AGB component of boles from non direct-logging activities to AGB litter pool
             if (c==ncwd) then
                
                currentPatch%cwd_AG_in(c) = currentPatch%cwd_AG_in(c) + (currentCohort%bdead+currentCohort%bsw) * &
                     SF_val_CWD_frac(c) * (dead_n_natural+dead_n_ilogging)  * &
                     EDPftvarcon_inst%allom_agb_frac(currentCohort%pft)
                
             else

                currentPatch%cwd_AG_in(c) = currentPatch%cwd_AG_in(c) + (currentCohort%bdead+currentCohort%bsw) * &
                     SF_val_CWD_frac(c) * dead_n  * &
                     EDPftvarcon_inst%allom_agb_frac(currentCohort%pft)

                ! Send AGB component of boles from direct-logging activities to export/harvest pool
                ! Generate trunk product (kgC/day/site)
                trunk_product = (currentCohort%bdead+currentCohort%bsw) * &
                      SF_val_CWD_frac(c) * dead_n_dlogging * EDPftvarcon_inst%allom_agb_frac(currentCohort%pft) * &
                      hlm_freq_day * currentPatch%area
                
                currentSite%flux_out = currentSite%flux_out + trunk_product

                ! Update diagnostics that track resource management
                currentSite%resources_management%trunk_product_site  = &
                      currentSite%resources_management%trunk_product_site + &
                      trunk_product
                ! Update diagnostics that track resource management
                currentSite%resources_management%trunk_product_site  = &
                      currentSite%resources_management%trunk_product_site + &
                      trunk_product
             end if

             ! Update diagnostics that track resource management
             currentSite%resources_management%delta_litter_stock  = &
                   currentSite%resources_management%delta_litter_stock + &
                   (currentCohort%bdead+currentCohort%bsw) * &
                   SF_val_CWD_frac(c) * (dead_n_natural+dead_n_ilogging) * & 
                   hlm_freq_day * currentPatch%area
             ! Update diagnostics that track resource management
             currentSite%resources_management%delta_biomass_stock = &
                   currentSite%resources_management%delta_biomass_stock + &
                   (currentCohort%bdead+currentCohort%bsw) * &
                   SF_val_CWD_frac(c) * dead_n * & 
                   hlm_freq_day * currentPatch%area
             
             if (currentPatch%cwd_AG_in(c) < 0.0_r8)then
                write(fates_log(),*) 'negative CWD in flux',currentPatch%cwd_AG_in(c), &
                      (currentCohort%bdead+currentCohort%bsw), dead_n
             endif

          end do
          ! Update diagnostics that track resource management
          currentSite%resources_management%delta_individual    = &
                currentSite%resources_management%delta_individual + &
                (dead_n_dlogging+dead_n_ilogging) * hlm_freq_day * currentPatch%area
          
       endif !canopy layer
       
       currentCohort => currentCohort%taller
    enddo  ! end loop over cohorts 

    do p = 1,numpft
       currentPatch%leaf_litter_in(p) = currentPatch%leaf_litter_in(p) + currentPatch%seed_decay(p) !KgC/m2/yr
    enddo

  end subroutine CWD_Input

  ! ============================================================================
  subroutine fragmentation_scaler( currentPatch, bc_in) 
    !
    ! !DESCRIPTION:
    ! Simple CWD fragmentation Model
    ! FIX(SPM, 091914) this should be a function as it returns a value in 
    ! currentPatch%fragmentation_scaler
    !
    ! !USES:

    use FatesSynchronizedParamsMod  , only : FatesSynchronizedParamsInst
    use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm
    use FatesConstantsMod, only : pi => pi_const
    !
    ! !ARGUMENTS    
    type(ed_patch_type), intent(inout) :: currentPatch
    type(bc_in_type),    intent(in)    :: bc_in

    !
    ! !LOCAL VARIABLES:
    logical  :: use_century_tfunc = .false.
    integer  :: j
    integer  :: ifp                   ! Index of a FATES Patch "ifp"
    real(r8) :: t_scalar
    real(r8) :: w_scalar
    real(r8) :: catanf                ! hyperbolic temperature function from CENTURY
    real(r8) :: catanf_30             ! hyperbolic temperature function from CENTURY
    real(r8) :: t1                    ! temperature argument
    real(r8) :: Q10                   ! temperature dependence
    real(r8) :: froz_q10              ! separate q10 for frozen soil respiration rates.
                                      ! default to same as above zero rates
    !----------------------------------------------------------------------

    catanf(t1) = 11.75_r8 +(29.7_r8 / pi) * atan( pi * 0.031_r8  * ( t1 - 15.4_r8 ))
    catanf_30 = catanf(30._r8)
    
    ifp = currentPatch%patchno 
    
    ! set "froz_q10" parameter
    froz_q10  = FatesSynchronizedParamsInst%froz_q10  
    Q10       = FatesSynchronizedParamsInst%Q10

    if ( .not. use_century_tfunc ) then
    !calculate rate constant scalar for soil temperature,assuming that the base rate constants 
    !are assigned for non-moisture limiting conditions at 25C. 
      if (bc_in%t_veg24_pa(ifp)  >=  tfrz) then
        t_scalar = Q10**((bc_in%t_veg24_pa(ifp)-(tfrz+25._r8))/10._r8)
                 !  Q10**((t_soisno(c,j)-(tfrz+25._r8))/10._r8)
      else
        t_scalar = (Q10**(-25._r8/10._r8))*(froz_q10**((bc_in%t_veg24_pa(ifp)-tfrz)/10._r8))
                  !Q10**(-25._r8/10._r8))*(froz_q10**((t_soisno(c,j)-tfrz)/10._r8)
      endif
    else
      ! original century uses an arctangent function to calculate the 
      ! temperature dependence of decomposition      
      t_scalar = max(catanf(bc_in%t_veg24_pa(ifp)-tfrz)/catanf_30,0.01_r8)
    endif    
   
    !Moisture Limitations   
    !BTRAN APPROACH - is quite simple, but max's out decomp at all unstressed 
    !soil moisture values, which is not realistic.  
    !litter decomp is proportional to water limitation on average... 
    w_scalar = sum(currentPatch%btran_ft(1:numpft))/numpft

    currentPatch%fragmentation_scaler =  min(1.0_r8,max(0.0_r8,t_scalar * w_scalar))
    
  end subroutine fragmentation_scaler
  
  ! ============================================================================
  subroutine cwd_out( currentSite, currentPatch, bc_in )
    !
    ! !DESCRIPTION:
    ! Simple CWD fragmentation Model
    ! spawn new cohorts of juveniles of each PFT             
    !
    ! !USES:
    use SFParamsMod, only : SF_val_max_decomp

    !
    ! !ARGUMENTS    
    type(ed_site_type), intent(inout), target  :: currentSite
    type(ed_patch_type), intent(inout), target :: currentPatch
    type(bc_in_type), intent(in)               :: bc_in
    
    !
    ! !LOCAL VARIABLES:
    integer :: c,ft
    !----------------------------------------------------------------------

    currentPatch%root_litter_out(:) = 0.0_r8
    currentPatch%leaf_litter_out(:) = 0.0_r8
    
    call fragmentation_scaler(currentPatch, bc_in)

    !Flux of coarse woody debris into decomposing litter pool. 

    currentPatch%cwd_ag_out(1:ncwd) = 0.0_r8
    currentPatch%cwd_bg_out(1:ncwd) = 0.0_r8
    currentPatch%leaf_litter_out(:) = 0.0_r8
    currentPatch%root_litter_out(:) = 0.0_r8
    
    do c = 1,ncwd  
       currentPatch%cwd_ag_out(c)      = max(0.0_r8,   currentPatch%cwd_ag(c) * &
            SF_val_max_decomp(c+1) * currentPatch%fragmentation_scaler )  
       currentPatch%cwd_bg_out(c)      = max(0.0_r8,   currentPatch%cwd_bg(c) * &
            SF_val_max_decomp(c+1) * currentPatch%fragmentation_scaler )
    enddo

    ! this is the rate at which dropped leaves stop being part of the burnable pool and begin to be part of the 
    ! decomposing pool. This should probably be highly sensitive to moisture, but also to the type of leaf 
    ! thick leaves can dry out before they are decomposed, for example. 
    ! this section needs further scientific input. 

    do ft = 1,numpft
       currentPatch%leaf_litter_out(ft) = max(0.0_r8,currentPatch%leaf_litter(ft)* SF_val_max_decomp(dl_sf) * &
            currentPatch%fragmentation_scaler )
       currentPatch%root_litter_out(ft) = max(0.0_r8,currentPatch%root_litter(ft)* SF_val_max_decomp(dl_sf) * &
            currentPatch%fragmentation_scaler )
       if ( currentPatch%leaf_litter_out(ft)<0.0_r8.or.currentPatch%root_litter_out(ft)<0.0_r8)then
         write(fates_log(),*) 'root or leaf out is negative?',SF_val_max_decomp(dl_sf),currentPatch%fragmentation_scaler
       endif
    enddo

    !add up carbon going into fragmenting pools
    currentSite%flux_out = currentSite%flux_out + sum(currentPatch%leaf_litter_out) * &
         currentPatch%area *hlm_freq_day!kgC/site/day
    currentSite%flux_out = currentSite%flux_out + sum(currentPatch%root_litter_out) * &
         currentPatch%area *hlm_freq_day!kgC/site/day
    currentSite%flux_out = currentSite%flux_out + sum(currentPatch%cwd_ag_out) * &
         currentPatch%area *hlm_freq_day!kgC/site/day
    currentSite%flux_out = currentSite%flux_out + sum(currentPatch%cwd_bg_out) * &
         currentPatch%area *hlm_freq_day!kgC/site/day

  end subroutine cwd_out



  subroutine flux_into_litter_pools(nsites, sites, bc_in, bc_out)
    ! Created by Charlie Koven and Rosie Fisher, 2014-2015
    ! take the flux out of the fragmenting litter pools and port into the decomposing litter pools. 
    ! in this implementation, decomposing pools are assumed to be humus and non-flammable, whereas fragmenting pools
    ! are assumed to be physically fragmenting but not respiring. This is a simplification, but allows us to 
    ! a) reconcile the need to track both chemical fractions (lignin, cellulose, labile) and size fractions (trunk, branch, etc.)
    ! b) to impose a realistic delay on the surge of nutrients into the litter pools when large CWD is added to the system via mortality
    
    ! because of the different subgrid structure, this subroutine includes the functionality that in the big-leaf BGC model, is calculated in SoilBiogeochemVerticalProfileMod
    
    ! The ED code is resolved at a daily timestep, but all of the CN-BGC fluxes are passed in as derivatives per second, 
    ! and then accumulated in the CNStateUpdate routines. One way of doing this is to pass back the CN fluxes per second, 
    ! and keep them constant for the whole day (making sure they are not overwritten.
    ! This means that the carbon gets passed back and forth between the photosynthesis code (fast timestepping) to the ED code (slow timestepping), back to the BGC code (fast timestepping).
    ! This means that the state update for the litter pools and for the CWD pools occurs at different timescales. 
    

    use EDTypesMod, only : AREA
    use FatesInterfaceMod, only : hlm_numlevdecomp_full
    use FatesInterfaceMod, only : hlm_numlevdecomp
    use EDPftvarcon, only : EDPftvarcon_inst
    use FatesConstantsMod, only : sec_per_day
    use FatesInterfaceMod, only : bc_in_type, bc_out_type
    use FatesInterfaceMod, only : hlm_use_vertsoilc
    use FatesConstantsMod, only : itrue
    use FatesGlobals, only : endrun => fates_endrun
    use EDParamsMod , only : ED_val_cwd_flig, ED_val_cwd_fcel


    implicit none   

    ! !ARGUMENTS    
    integer                 , intent(in)            :: nsites
    type(ed_site_type)      , intent(inout), target :: sites(nsites)
    type(bc_in_type)        , intent(in)            :: bc_in(:)
    type(bc_out_type)       , intent(inout)           :: bc_out(:)
    !
    ! !LOCAL VARIABLES:
    type (ed_patch_type)  , pointer :: currentPatch
    type (ed_cohort_type) , pointer :: currentCohort
    type(ed_site_type), pointer :: cs
    integer p,ci,j,s
    real(r8) time_convert    ! from year to seconds
    real(r8) mass_convert    ! ED uses kg, CLM uses g
    integer           :: begp,endp
    integer           :: begc,endc                                    !bounds 
    !------------------------------------------------------------------------
    real(r8) :: cinput_rootfr(1:maxpft, 1:hlm_numlevdecomp_full)      ! column by pft root fraction used for calculating inputs
    real(r8) :: croot_prof_perpatch(1:hlm_numlevdecomp_full)
    real(r8) :: surface_prof(1:hlm_numlevdecomp_full)
    integer  :: ft
    real(r8) :: rootfr_tot(1:maxpft), biomass_bg_ft(1:maxpft)
    real(r8) :: surface_prof_tot, leaf_prof_sum, stem_prof_sum, froot_prof_sum, biomass_bg_tot
    real(r8) :: delta

    ! NOTE(bja, 201608) these were removed from clm in clm4_5_10_r187
    logical, parameter :: exponential_rooting_profile = .true.
    logical, parameter :: pftspecific_rootingprofile = .true.

    ! NOTE(bja, 201608) as of clm4_5_10_r187 rootprof_exp is now a
    ! private function level parameter in RootBiophysMod.F90::exponential_rootfr()
    real(r8), parameter :: rootprof_exp  = 3.  ! how steep profile is
    ! for root C inputs (1/ e-folding depth) (1/m)

    ! NOTE(rgk, 201705) this parameter was brought over from SoilBiogeochemVerticalProfile
    ! how steep profile is for surface components (1/ e_folding depth) (1/m) 
    real(r8),  parameter :: surfprof_exp  = 10.

    ! NOTE(bja, 201608) as of clm4_5_10_r187 rootprof_beta is now a
    ! two dimensional array with the second dimension being water,1,
    ! or carbon,2,. These are currently hard coded, but may be
    ! overwritten by the namelist.

    ! Note cdk 2016/08 we actually want to use the carbon index here rather than the water index.  
    ! Doing so will be answer changing though so perhaps easiest to do this in steps.
    integer, parameter :: rooting_profile_varindex_water = 1

    real(r8) :: leaf_prof(1:nsites, 1:hlm_numlevdecomp)
    real(r8) :: froot_prof(1:nsites,  1:maxpft, 1:hlm_numlevdecomp)
    real(r8) :: croot_prof(1:nsites, 1:hlm_numlevdecomp)
    real(r8) :: stem_prof(1:nsites, 1:hlm_numlevdecomp)

    delta = 0.001_r8    
    !no of seconds in a year. 
    time_convert =  365.0_r8*sec_per_day

    ! number of grams in a kilogram
    mass_convert = 1000._r8
    
      
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! first calculate vertical profiles
    ! define two types of profiles: 
    ! (1) a surface profile, for leaves and stem inputs, which is the same for each
    ! pft but differs from one site to the next to avoid inputting any C into permafrost or bedrock
    ! (2) a fine root profile, which is indexed by both site and pft, differs for 
    ! each pft and also from one site to the next to avoid inputting any C into permafrost or bedrock
    ! (3) a coarse root profile, which is the root-biomass=weighted average of the fine root profiles
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    if (hlm_use_vertsoilc == itrue) then

       ! initialize profiles to zero
       leaf_prof(1:nsites, :)               = 0._r8
       froot_prof(1:nsites, 1:maxpft, :)    = 0._r8
       stem_prof(1:nsites, :)               = 0._r8
       
       do s = 1,nsites
          ! define a single shallow surface profile for surface additions (leaves, stems, and N deposition)
          surface_prof(:) = 0._r8
          do j = 1, hlm_numlevdecomp
             surface_prof(j) = exp(-surfprof_exp * bc_in(s)%z_sisl(j)) / bc_in(s)%dz_decomp_sisl(j)
          end do
          
          cinput_rootfr(:,:)     = 0._r8
          
          ! calculate pft-specific rooting profiles in the absence of permafrost or bedrock limitations
          if ( exponential_rooting_profile ) then
             if ( .not. pftspecific_rootingprofile ) then
                ! define rooting profile from exponential parameters
                do ft = 1, numpft
                   do j = 1, hlm_numlevdecomp
                      cinput_rootfr(ft,j) = exp(-rootprof_exp *  bc_in(s)%z_sisl(j)) / bc_in(s)%dz_decomp_sisl(j)
                   end do
                end do
             else
                ! use beta distribution parameter from Jackson et al., 1996
                do ft = 1, numpft
                   do j = 1, hlm_numlevdecomp
                      cinput_rootfr(ft,j) = &
                            ( EDPftvarcon_inst%rootprof_beta(ft, rooting_profile_varindex_water) ** & 
                            (bc_in(s)%zi_sisl(j-1)*100._r8) - &
                            EDPftvarcon_inst%rootprof_beta(ft, rooting_profile_varindex_water) ** & 
                            (bc_in(s)%zi_sisl(j)*100._r8) ) &
                            / bc_in(s)%dz_decomp_sisl(j)
                   end do
                end do
             endif
          else
             do ft = 1,numpft 
                do j = 1, hlm_numlevdecomp
                   ! use standard CLM root fraction profiles;
                   cinput_rootfr(ft,j) =  ( .5_r8*( &
                         exp(-EDPftvarcon_inst%roota_par(ft) * bc_in(s)%zi_sisl(j-1))  &
                         + exp(-EDPftvarcon_inst%rootb_par(ft) * bc_in(s)%zi_sisl(j-1))  &
                         - exp(-EDPftvarcon_inst%roota_par(ft) * bc_in(s)%zi_sisl(j))    &
                         - exp(-EDPftvarcon_inst%rootb_par(ft) * bc_in(s)%zi_sisl(j))))  &
                         / bc_in(s)%dz_decomp_sisl(j)
                end do
             end do
          endif
          
          !
          ! now add permafrost constraint: integrate rootfr over active layer of soil site,
          ! truncate below permafrost or bedrock table where present, and rescale so that integral = 1
          rootfr_tot(:) = 0._r8
          
          surface_prof_tot = 0._r8
          !
          do j = 1, min(max(bc_in(s)%max_rooting_depth_index_col, 1), hlm_numlevdecomp)
             surface_prof_tot = surface_prof_tot + surface_prof(j)  * bc_in(s)%dz_decomp_sisl(j)
          end do
          do ft = 1,numpft
             do j = 1, min(max(bc_in(s)%max_rooting_depth_index_col, 1), hlm_numlevdecomp)
                rootfr_tot(ft) = rootfr_tot(ft) + cinput_rootfr(ft,j) * bc_in(s)%dz_decomp_sisl(j)
             end do
          end do
          !
          ! rescale the fine root profile
          do ft = 1,numpft
             if ( (bc_in(s)%max_rooting_depth_index_col > 0) .and. (rootfr_tot(ft) > 0._r8) ) then
                ! where there is not permafrost extending to the surface, integrate the profiles over the active layer
                ! this is equivalent to integrating over all soil layers outside of permafrost regions
                do j = 1, min(max(bc_in(s)%max_rooting_depth_index_col, 1), hlm_numlevdecomp)
                   froot_prof(s,ft,j) = cinput_rootfr(ft,j) / rootfr_tot(ft)
                end do
             else
                ! if fully frozen, or no roots, put everything in the top layer
                froot_prof(s,ft,1) = 1._r8/bc_in(s)%dz_decomp_sisl(1)
             endif
          end do
          !
          ! rescale the shallow profiles
          if ( (bc_in(s)%max_rooting_depth_index_col > 0) .and. (surface_prof_tot > 0._r8) ) then
             ! where there is not permafrost extending to the surface, integrate the profiles over the active layer
             ! this is equivalent to integrating over all soil layers outside of permafrost regions
             do j = 1, min(max(bc_in(s)%max_rooting_depth_index_col, 1), hlm_numlevdecomp)
                ! set all surface processes to shallower profile
                leaf_prof(s,j) = surface_prof(j)/ surface_prof_tot
                stem_prof(s,j) = surface_prof(j)/ surface_prof_tot
             end do
          else
             ! if fully frozen, or no roots, put everything in the top layer
             leaf_prof(s,1) = 1._r8/bc_in(s)%dz_decomp_sisl(1)
             stem_prof(s,1) = 1._r8/bc_in(s)%dz_decomp_sisl(1)
             do j = 2, hlm_numlevdecomp
                leaf_prof(s,j) = 0._r8
                stem_prof(s,j) = 0._r8
             end do
          endif
       end do
       
    else
       
       ! for one layer decomposition model, set profiles to unity
       leaf_prof(1:nsites, :) = 1._r8
       froot_prof(1:nsites, 1:numpft, :) = 1._r8
       stem_prof(1:nsites, :) = 1._r8
       
    end if
    
    ! sanity check to ensure they integrate to 1
    do s = 1, nsites
       ! check the leaf and stem profiles
       leaf_prof_sum = 0._r8
       stem_prof_sum = 0._r8
       do j = 1, hlm_numlevdecomp
          leaf_prof_sum = leaf_prof_sum + leaf_prof(s,j) *  bc_in(s)%dz_decomp_sisl(j)
          stem_prof_sum = stem_prof_sum + stem_prof(s,j) *  bc_in(s)%dz_decomp_sisl(j)
       end do
       if ( ( abs(stem_prof_sum - 1._r8) > delta ) .or.  ( abs(leaf_prof_sum - 1._r8) > delta ) ) then
          write(fates_log(), *) 'profile sums: ',  leaf_prof_sum, stem_prof_sum
          write(fates_log(), *) 'surface_prof: ', surface_prof
          write(fates_log(), *) 'surface_prof_tot: ', surface_prof_tot
          write(fates_log(), *) 'leaf_prof: ',  leaf_prof(s,:)
          write(fates_log(), *) 'stem_prof: ',  stem_prof(s,:)
          write(fates_log(), *) 'max_rooting_depth_index_col: ', bc_in(s)%max_rooting_depth_index_col
          write(fates_log(), *) 'bc_in(s)%dz_decomp_sisl: ',  bc_in(s)%dz_decomp_sisl            
          call endrun(msg=errMsg(sourcefile, __LINE__))
       endif
       ! now check each fine root profile
       do ft = 1,numpft 
          froot_prof_sum = 0._r8
          do j = 1, hlm_numlevdecomp
             froot_prof_sum = froot_prof_sum + froot_prof(s,ft,j) *  bc_in(s)%dz_decomp_sisl(j)
          end do
          if ( ( abs(froot_prof_sum - 1._r8) > delta ) ) then
             write(fates_log(), *) 'profile sums: ', froot_prof_sum
             call endrun(msg=errMsg(sourcefile, __LINE__))
          endif
       end do
    end do
    
    ! zero the site-level C input variables
    do s = 1, nsites
       do j = 1, hlm_numlevdecomp
          bc_out(s)%FATES_c_to_litr_lab_c_col(j) = 0._r8
          bc_out(s)%FATES_c_to_litr_cel_c_col(j) = 0._r8
          bc_out(s)%FATES_c_to_litr_lig_c_col(j) = 0._r8
          croot_prof(s,j)         = 0._r8
       end do
    end do
    
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! now disaggregate the inputs vertically, using the vertical profiles
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do s = 1,nsites
         
         !      do g = bounds%begg,bounds%endg
         !         if (firstsoilpatch(g) >= 0 .and. ed_allsites_inst(g)%istheresoil) then 
         currentPatch => sites(s)%oldest_patch
         
         do while(associated(currentPatch))
            
            ! the CWD pools lose information about which PFT they came from; for the stems this doesn't matter as they all have the same profile, 
            ! however for the coarse roots they may have different profiles.  to approximately recover this information, loop over all cohorts in patch 
            ! to calculate the total root biomass in that patch of each pft, and then rescale the croot_prof as the weighted average of the froot_prof
            biomass_bg_ft(:) = 0._r8
            currentCohort => currentPatch%tallest
            do while(associated(currentCohort))      
               biomass_bg_ft(currentCohort%pft) = biomass_bg_ft(currentCohort%pft) + &
                    currentCohort%b * (currentCohort%n / currentPatch%area) * (1.0_r8-EDPftvarcon_inst%allom_agb_frac(currentCohort%pft))
               currentCohort => currentCohort%shorter
            enddo !currentCohort
            ! 
            biomass_bg_tot = 0._r8
            do ft = 1,numpft 
               biomass_bg_tot = biomass_bg_tot + biomass_bg_ft(ft)
            end do
            !         
            do j = 1, hlm_numlevdecomp
               ! zero this for each patch
               croot_prof_perpatch(j) = 0._r8
            end do
            !
            if ( biomass_bg_tot .gt. 0._r8) then
               do ft = 1,numpft 
                  do j = 1, hlm_numlevdecomp
                     croot_prof_perpatch(j) = croot_prof_perpatch(j) + froot_prof(s,ft,j) * biomass_bg_ft(ft) / biomass_bg_tot
                  end do
               end do
            else ! no biomass
               croot_prof_perpatch(1) = 1./bc_in(s)%dz_decomp_sisl(1)
            end if

            !
            ! add croot_prof as weighted average (weighted by patch area) of croot_prof_perpatch
            do j = 1, hlm_numlevdecomp
               croot_prof(s, j) = croot_prof(s, j) + croot_prof_perpatch(j) * currentPatch%area / AREA
            end do
            !
            ! now disaggregate, vertically and by decomposition substrate type, the actual fluxes from CWD and litter pools
            !
            ! do c = 1, ncwd
            !    write(fates_log(),*)'cdk CWD_AG_out', c, currentpatch%CWD_AG_out(c), ED_val_cwd_fcel, currentpatch%area/AREA
            !    write(fates_log(),*)'cdk CWD_BG_out', c, currentpatch%CWD_BG_out(c), ED_val_cwd_fcel, currentpatch%area/AREA
            ! end do
            ! do ft = 1,numpft
            !    write(fates_log(),*)'cdk leaf_litter_out', ft, currentpatch%leaf_litter_out(ft), ED_val_cwd_fcel, currentpatch%area/AREA
            !    write(fates_log(),*)'cdk root_litter_out', ft, currentpatch%root_litter_out(ft), ED_val_cwd_fcel, currentpatch%area/AREA
            ! end do
            ! !
            ! CWD pools fragmenting into decomposing litter pools. 
            do ci = 1, ncwd
               do j = 1, hlm_numlevdecomp
                  bc_out(s)%FATES_c_to_litr_cel_c_col(j) = bc_out(s)%FATES_c_to_litr_cel_c_col(j) + &
                       currentpatch%CWD_AG_out(ci) * ED_val_cwd_fcel * currentpatch%area/AREA * stem_prof(s,j)  
                  bc_out(s)%FATES_c_to_litr_lig_c_col(j) = bc_out(s)%FATES_c_to_litr_lig_c_col(j) + &
                       currentpatch%CWD_AG_out(ci) * ED_val_cwd_flig * currentpatch%area/AREA * stem_prof(s,j)
                  !
                  bc_out(s)%FATES_c_to_litr_cel_c_col(j) = bc_out(s)%FATES_c_to_litr_cel_c_col(j) + &
                       currentpatch%CWD_BG_out(ci) * ED_val_cwd_fcel * currentpatch%area/AREA * croot_prof_perpatch(j)
                  bc_out(s)%FATES_c_to_litr_lig_c_col(j) = bc_out(s)%FATES_c_to_litr_lig_c_col(j) + &
                       currentpatch%CWD_BG_out(ci) * ED_val_cwd_flig * currentpatch%area/AREA * croot_prof_perpatch(j)
               end do
            end do
            
            ! leaf and fine root pools. 
            do ft = 1,numpft
               do j = 1, hlm_numlevdecomp
                  bc_out(s)%FATES_c_to_litr_lab_c_col(j) = bc_out(s)%FATES_c_to_litr_lab_c_col(j) + &
                       currentpatch%leaf_litter_out(ft) * EDPftvarcon_inst%lf_flab(ft) * currentpatch%area/AREA * leaf_prof(s,j)
                  bc_out(s)%FATES_c_to_litr_cel_c_col(j) = bc_out(s)%FATES_c_to_litr_cel_c_col(j) + &
                       currentpatch%leaf_litter_out(ft) * EDPftvarcon_inst%lf_fcel(ft) * currentpatch%area/AREA * leaf_prof(s,j)
                  bc_out(s)%FATES_c_to_litr_lig_c_col(j) = bc_out(s)%FATES_c_to_litr_lig_c_col(j) + &
                       currentpatch%leaf_litter_out(ft) * EDPftvarcon_inst%lf_flig(ft) * currentpatch%area/AREA * leaf_prof(s,j)
                  !
                  bc_out(s)%FATES_c_to_litr_lab_c_col(j) = bc_out(s)%FATES_c_to_litr_lab_c_col(j) + &
                       currentpatch%root_litter_out(ft) * EDPftvarcon_inst%fr_flab(ft) * currentpatch%area/AREA * froot_prof(s,ft,j)
                  bc_out(s)%FATES_c_to_litr_cel_c_col(j) = bc_out(s)%FATES_c_to_litr_cel_c_col(j) + &
                       currentpatch%root_litter_out(ft) * EDPftvarcon_inst%fr_fcel(ft) * currentpatch%area/AREA * froot_prof(s,ft,j)
                  bc_out(s)%FATES_c_to_litr_lig_c_col(j) = bc_out(s)%FATES_c_to_litr_lig_c_col(j) + &
                       currentpatch%root_litter_out(ft) * EDPftvarcon_inst%fr_flig(ft) * currentpatch%area/AREA * froot_prof(s,ft,j)
                  !
                  !! and seed_decay too.  for now, use the same lability fractions as for leaf litter
                  bc_out(s)%FATES_c_to_litr_lab_c_col(j) = bc_out(s)%FATES_c_to_litr_lab_c_col(j) + &
                       currentpatch%seed_decay(ft) * EDPftvarcon_inst%lf_flab(ft) * currentpatch%area/AREA * leaf_prof(s,j)
                  bc_out(s)%FATES_c_to_litr_cel_c_col(j) = bc_out(s)%FATES_c_to_litr_cel_c_col(j) + &
                       currentpatch%seed_decay(ft) * EDPftvarcon_inst%lf_fcel(ft) * currentpatch%area/AREA * leaf_prof(s,j)
                  bc_out(s)%FATES_c_to_litr_lig_c_col(j) = bc_out(s)%FATES_c_to_litr_lig_c_col(j) + &
                       currentpatch%seed_decay(ft) * EDPftvarcon_inst%lf_flig(ft) * currentpatch%area/AREA * leaf_prof(s,j)
                  !
               enddo
            end do
              
              currentPatch => currentPatch%younger
           end do !currentPatch

        end do  ! do sites(s)
     
        do s = 1, nsites
           do j = 1, hlm_numlevdecomp                    
              ! time unit conversion
              bc_out(s)%FATES_c_to_litr_lab_c_col(j)=bc_out(s)%FATES_c_to_litr_lab_c_col(j) * mass_convert / time_convert
              bc_out(s)%FATES_c_to_litr_cel_c_col(j)=bc_out(s)%FATES_c_to_litr_cel_c_col(j) * mass_convert / time_convert
              bc_out(s)%FATES_c_to_litr_lig_c_col(j)=bc_out(s)%FATES_c_to_litr_lig_c_col(j) * mass_convert / time_convert
           end do
        end do
        
        ! write(fates_log(),*)'cdk FATES_c_to_litr_lab_c: ', FATES_c_to_litr_lab_c
        ! write_col(fates_log(),*)'cdk FATES_c_to_litr_cel_c: ', FATES_c_to_litr_cel_c    
        ! write_col(fates_log(),*)'cdk FATES_c_to_litr_lig_c: ', FATES_c_to_litr_lig_c
        ! write_col(fates_log(),*)'cdk hlm_numlevdecomp_full,  bounds%begc, bounds%endc: ', hlm_numlevdecomp_full, bounds%begc, bounds%endc
        ! write(fates_log(),*)'cdk leaf_prof: ', leaf_prof
        ! write(fates_log(),*)'cdk stem_prof: ', stem_prof    
        ! write(fates_log(),*)'cdk froot_prof: ', froot_prof
        ! write(fates_log(),*)'cdk croot_prof_perpatch: ', croot_prof_perpatch
        ! write(fates_log(),*)'cdk croot_prof: ', croot_prof

    end subroutine flux_into_litter_pools

end module EDPhysiologyMod
