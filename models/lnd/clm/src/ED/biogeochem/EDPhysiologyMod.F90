module EDPhysiologyMod

#include "shr_assert.h"

  ! ============================================================================
  ! Miscellaneous physiology routines from ED. 
  ! ============================================================================

  use shr_kind_mod        , only : r8 => shr_kind_r8
  use clm_varctl          , only : iulog 
  use TemperatureType     , only : temperature_type
  use SoilStateType       , only : soilstate_type
  use WaterstateType      , only : waterstate_type
  use pftconMod           , only : pftcon
  use EDEcophysContype    , only : EDecophyscon
  use EDCohortDynamicsMod , only : allocate_live_biomass, zero_cohort, create_cohort, fuse_cohorts, sort_cohorts
  use EDPhenologyType     , only : ed_phenology_type
  use EDTypesMod          , only : dg_sf, dinc_ed, external_recruitment
  use EDTypesMod          , only : ncwd, nlevcan_ed, n_sub, numpft_ed, senes
  use EDTypesMod          , only : ed_site_type, ed_patch_type, ed_cohort_type

  implicit none
  private

  public :: canopy_derivs
  public :: non_canopy_derivs
  public :: trim_canopy
  public :: phenology
  public :: phenology_leafonoff
  public :: Growth_Derivatives
  public :: recruitment
  public :: cwd_input
  public :: cwd_out
  public :: fragmentation_scaler
  public :: seeds_in
  public :: seed_decay
  public :: seed_germination
  ! ============================================================================

contains

  ! ============================================================================
  subroutine canopy_derivs( currentPatch )
    !
    ! !DESCRIPTION:
    ! spawn new cohorts of juveniles of each PFT             
    !
    ! !USES:
    !
    ! !ARGUMENTS    
    type(ed_patch_type) , intent(inout), target :: currentPatch
    !
    ! !LOCAL VARIABLES:
    type(ed_cohort_type), pointer ::currentCohort
    !----------------------------------------------------------------------

    ! call plant growth functions

    currentCohort => currentPatch%shortest

    do while(associated(currentCohort))
       call Growth_Derivatives(currentCohort)
       currentCohort => currentCohort%taller
    enddo

  end subroutine canopy_derivs

  ! ============================================================================
  subroutine non_canopy_derivs( currentPatch, temperature_inst, soilstate_inst, waterstate_inst)
    !
    ! !DESCRIPTION:
    ! Returns time differentials of the state vector
    !
    ! !USES:
    !
    ! !ARGUMENTS    
    type(ed_patch_type)    , intent(inout) :: currentPatch
    type(temperature_type) , intent(in)    :: temperature_inst
    type(soilstate_type)   , intent(in)    :: soilstate_inst
    type(waterstate_type)  , intent(in)    :: waterstate_inst
    !
    ! !LOCAL VARIABLES:
    integer c,p
    !----------------------------------------------------------------------

    currentPatch%leaf_litter_in(:)   = 0.0_r8
    currentPatch%root_litter_in(:)   = 0.0_r8
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
    call seeds_in(currentPatch)
    call seed_decay(currentPatch)
    call seed_germination(currentPatch)

    ! update fragmenting pool fluxes
    call cwd_input(currentPatch)
    call cwd_out( currentPatch, temperature_inst, soilstate_inst, waterstate_inst)

    do p = 1,numpft_ed
       currentPatch%dseed_dt(p) = currentPatch%seeds_in(p) - currentPatch%seed_decay(p) - currentPatch%seed_germination(p)
    enddo   
    
    do c = 1,ncwd
       currentPatch%dcwd_AG_dt(c) = currentPatch%cwd_AG_in(c) - currentPatch%cwd_AG_out(c) 
       currentPatch%dcwd_BG_dt(c) = currentPatch%cwd_BG_in(c) - currentPatch%cwd_BG_out(c) 
    enddo

    do p = 1,numpft_ed
       currentPatch%dleaf_litter_dt(p) = currentPatch%leaf_litter_in(p) - currentPatch%leaf_litter_out(p) 
       currentPatch%droot_litter_dt(p) = currentPatch%root_litter_in(p) - currentPatch%root_litter_out(p) 
    enddo

    currentPatch%leaf_litter_in(:)  = 0.0_r8
    currentPatch%root_litter_in(:)  = 0.0_r8
    currentPatch%leaf_litter_out(:) = 0.0_r8
    currentPatch%root_litter_out(:) = 0.0_r8
    currentPatch%CWD_AG_in(:)       = 0.0_r8
    currentPatch%cwd_bg_in(:)       = 0.0_r8
    currentPatch%CWD_AG_out(:)      = 0.0_r8
    currentPatch%cwd_bg_out(:)      = 0.0_r8

  end subroutine non_canopy_derivs

  ! ============================================================================
  subroutine trim_canopy( currentSite )
    !
    ! !DESCRIPTION:
    ! Canopy trimming / leaf optimisation. Removes leaves in negative annual carbon balance. 
    !
    ! !USES:
    !
    use EDParamsMod,          only : ED_val_grperc
    use EDGrowthFunctionsMod, only : tree_lai
    !
    ! !ARGUMENTS    
    type (ed_site_type),intent(inout), target :: currentSite
    !
    ! !LOCAL VARIABLES:
    type (ed_cohort_type) , pointer :: currentCohort
    type (ed_patch_type)  , pointer :: currentPatch

    real(r8) :: inc        ! rate at which canopy acclimates to uptake 
    real(r8) :: trim_limit ! this is the limit of the canopy trimming routine, so that trees 
                           ! can't just lose all their leaves and have no reproductive costs.
    integer  :: z          ! leaf layer
    integer  :: trimmed    ! was this layer trimmed in this year? If not expand the canopy. 

    trim_limit = 0.3_r8    ! Arbitrary limit to reductions in leaf area with stress. Without this nothing ever dies.  
    inc = 0.03_r8          ! Arbitrary incremental change in trimming function. Controls 
                           ! rate at which leaves are optimised to their environment. 
    !----------------------------------------------------------------------

    currentPatch => currentSite%youngest_patch

    do while(associated(currentPatch))
       currentCohort => currentPatch%tallest
       do while (associated(currentCohort)) 
          trimmed = 0    
          currentCohort%treelai = tree_lai(currentCohort)    
          currentCohort%nv = ceiling((currentCohort%treelai+currentCohort%treesai)/dinc_ed)
          if (currentCohort%nv > nlevcan_ed)then
             write(iulog,*) 'nv > nlevcan_ed',currentCohort%nv,currentCohort%treelai,currentCohort%treesai, &
                  currentCohort%c_area,currentCohort%n,currentCohort%bl
          endif

          !Leaf cost vs netuptake for each leaf layer. 
          do z = 1,nlevcan_ed
             if (currentCohort%year_net_uptake(z) /= 999._r8)then !there was activity this year in this leaf layer. 
                !Leaf Cost kgC/m2/year-1
                !decidous costs. 
                if (pftcon%season_decid(currentCohort%pft) == 1.or.pftcon%stress_decid(currentCohort%pft) == 1)then 
                   currentCohort%leaf_cost =  1._r8/(pftcon%slatop(currentCohort%pft)*1000_r8)
                   currentCohort%leaf_cost = currentCohort%leaf_cost + 1.0_r8/(pftcon%slatop(currentCohort%pft)*1000_r8) * &
                        pftcon%froot_leaf(currentCohort%pft) / EDecophyscon%root_long(currentCohort%pft)
                   currentCohort%leaf_cost = currentCohort%leaf_cost * (ED_val_grperc+1._r8)
                else !evergreen costs
                   currentCohort%leaf_cost = 1.0_r8/(pftcon%slatop(currentCohort%pft)* &
                        pftcon%leaf_long(currentCohort%pft)*1000_r8) !convert from sla in m2g-1 to m2kg-1 
                   currentCohort%leaf_cost = currentCohort%leaf_cost + 1.0_r8/(pftcon%slatop(currentCohort%pft)*1000_r8) * &
                        pftcon%froot_leaf(currentCohort%pft) / EDecophyscon%root_long(currentCohort%pft)
                   currentCohort%leaf_cost = currentCohort%leaf_cost * (ED_val_grperc+1._r8)
                endif
                if (currentCohort%year_net_uptake(z) < currentCohort%leaf_cost)then
                   if (currentCohort%canopy_trim > trim_limit)then
                      !  write(iulog,*) 'trimming leaves',currentCohort%canopy_trim,currentCohort%leaf_cost
                      ! keep trimming until none of the canopy is in negative carbon balance.              
                      if (currentCohort%hite > EDecophyscon%hgt_min(currentCohort%pft))then
                         currentCohort%canopy_trim = currentCohort%canopy_trim - inc    
                         if (pftcon%evergreen(currentCohort%pft) /= 1)then
                            currentCohort%laimemory = currentCohort%laimemory*(1.0_r8 - inc) 
                         endif
                         trimmed = 1
                      endif
                   endif
                endif
             endif !leaf activity? 
          enddo !z
          if (currentCohort%NV.gt.2)then
             write(iulog,*) 'nv>4',currentCohort%year_net_uptake(1:6),currentCohort%leaf_cost,&
             currentCohort%canopy_trim
          endif
       
          currentCohort%year_net_uptake(:) = 999.0_r8
          if (trimmed == 0.and.currentCohort%canopy_trim < 1.0_r8)then
             currentCohort%canopy_trim = currentCohort%canopy_trim + inc
          endif 
          ! write(iulog,*) 'trimming',currentCohort%canopy_trim
         
          ! currentCohort%canopy_trim = 1.0_r8 !FIX(RF,032414) this turns off ctrim for now. 
          currentCohort => currentCohort%shorter
       enddo
       currentPatch => currentPatch%older
    enddo

  end subroutine trim_canopy

  ! ============================================================================
  subroutine phenology( currentSite, ed_phenology_inst, temperature_inst, waterstate_inst)
    !
    ! !DESCRIPTION:
    ! Phenology. 
    !
    ! !USES:
    use clm_varcon, only : tfrz
    use EDTypesMod, only : udata
    !
    ! !ARGUMENTS:
    type(ed_site_type)      , intent(inout), pointer:: currentSite
    type(ed_phenology_type) , intent(in)            :: ed_phenology_inst
    type(temperature_type)  , intent(in)            :: temperature_inst
    type(waterstate_type)   , intent(in)            :: waterstate_inst
    !
    ! !LOCAL VARIABLES:
    real(r8), pointer :: t_veg24(:) 
    real(r8), pointer :: ED_GDD_patch(:)     
    integer  :: g            ! grid point  
    integer  :: t            ! day of year
    integer  :: ncolddays    ! no days underneath the threshold for leaf drop
    integer  :: ncolddayslim ! critical no days underneath the threshold for leaf drop
    integer  :: i
    integer  :: timesincedleafon,timesincedleafoff,timesinceleafon,timesinceleafoff
    real(r8) :: gdd_threshold
    real(r8) :: a,b,c        ! params of leaf-pn model from botta et al. 2000. 
    real(r8) :: cold_t       ! threshold below which cold days are counted 
    real(r8) :: coldday      ! definition of a 'chilling day' for botta model 
    real(r8) :: ncdstart     ! beginning of counting period for growing degree days.
    real(r8) :: drought_threshold
    real(r8) :: off_time     ! minimum number of days between leaf off and leaf on for drought phenology 
    real(r8) :: temp_in_C    ! daily averaged temperature in celcius
    real(r8) :: mindayson 
    !------------------------------------------------------------------------

    t_veg24       => temperature_inst%t_veg24_patch ! Input:  [real(r8) (:)]  avg pft vegetation temperature for last 24 hrs    
    ED_GDD_patch  => ed_phenology_inst%ED_GDD_patch ! Input:  [real(r8) (:)]  growing deg. days base 0 deg C (ddays)

    g = currentSite%clmgcell

    ! Parameter of drought decid leaf loss in mm in top layer...FIX(RF,032414) 
    ! - this is arbitrary and poorly understood. Needs work. ED_
    drought_threshold = 0.15 
    off_time = 100.0_r8

    !Parameters of Botta et al. 2000 GCB,6 709-725 
    a = -68.0_r8
    b = 638.0_r8
    c = -0.001_r8
    coldday = 5.0_r8

    mindayson = 30
     
    !Parameters from SDGVM model of senesence
    ncolddayslim = 5
    cold_t   = 7.5_r8

    t  = udata%time_period
    temp_in_C = t_veg24(currentSite%oldest_patch%clm_pno-1) - tfrz

    !-----------------Cold Phenology--------------------!              

    !Zero growing degree and chilling day counters
    if (currentSite%lat > 0)then
       ncdstart = 270._r8; !Northern Hemisphere begining November
    else
       ncdstart = 120._r8;  !Southern Hemisphere beginning May
    endif
    
    ! FIX(SPM,032414) - this will only work for the first year, no?
    if (t == ncdstart)then
       currentSite%ncd = 0._r8
    endif

    !Accumulate growing/chilling days after start of counting period
    if (temp_in_C  <  coldday)then
       currentSite%ncd = currentSite%ncd + 1.0_r8
    endif

    gdd_threshold = a + b*exp(c*currentSite%ncd) !GDD accumulation function, which also depends on chilling days.

    !Accumulate temperature of last 10 days.
    currentSite%last_n_days(2:senes) =  currentSite%last_n_days(1:senes-1)
    currentSite%last_n_days(1) = temp_in_C                                      
    !count number of days for leaves off
    ncolddays = 0
    do i = 1,senes
       if (currentSite%last_n_days(i) < cold_t)then
          ncolddays = ncolddays + 1
       endif
    enddo

    timesinceleafoff = t - currentSite%leafoffdate
    if (t < currentSite%leafoffdate)then
       timesinceleafoff = t +(365-currentSite%leafoffdate)
    endif

    !LEAF ON: COLD DECIDUOUS. Needs to
    !1) have exceeded the growing degree day threshold 
    !2) The leaves should not be on already
    !3) There should have been at least on chilling day in the counting period.  
    if (ED_GDD_patch(currentSite%oldest_patch%clm_pno) > gdd_threshold)then
       if (currentSite%status == 1)then
             if (currentSite%ncd >= 1)then
          currentSite%status = 2     !alter status of site to 'leaves on'
          currentSite%leafondate = t  !record leaf on date   
                write(iulog,*) 'leaves on'
             endif !ncd
       endif !status
    endif !GDD

    timesinceleafon = t - currentSite%leafondate
    if (t < currentSite%leafondate)then
       timesinceleafon = t +(365-currentSite%leafondate)
       endif

    !LEAF OFF: COLD THRESHOLD
    !Needs to:
    !1) have exceeded the number of cold days threshold
    !2) have exceeded the minimum leafon time.
    !3) The leaves should not be off already
    !4) The day of the year should be larger than the counting period. (not sure if we need this/if it will break the restarting)
    
    if (ncolddays > ncolddayslim)then
     if (timesinceleafon > mindayson)then
       if (currentSite%status == 2)then
          currentSite%status = 1        !alter status of site to 'leaves on'
          currentSite%leafoffdate = t   !record leaf off date   
          write(iulog,*) 'leaves off'
       endif
    endif
    endif

    !LEAF OFF: COLD LIFESPAN THRESHOLD
    if (timesinceleafoff > 360)then !remove leaves after a whole year when there is no 'off' period.  
       if (currentSite%status == 2)then
          currentSite%status = 1        !alter status of site to 'leaves on'
          currentSite%leafoffdate = t   !record leaf off date   
          write(iulog,*) 'leaves off'
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
    currentSite%water_memory(1) = waterstate_inst%h2osoi_vol_col(currentSite%clmcolumn,1) 
    do i = 1,9 !shift memory along one
       currentSite%water_memory(11-i) = currentSite%water_memory(10-i)
    enddo

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
       if (sum(currentSite%water_memory(1:10)/10._r8) >= drought_threshold.and.currentSite%dstatus == 1.and.t >= 10)then 
          ! leave some minimum time between leaf off and leaf on to prevent 'flickering'.  
          if (timesincedleafoff > off_time)then  
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
       !Are the leaves at the end of their lives? !FIX(RF,0401014)- this is hardwiring....
       if (timesincedleafon > 365.0*pftcon%leaf_long(7))then 
          currentSite%dstatus = 1         !alter status of site to 'leaves on'
          currentSite%dleafoffdate = t    !record leaf on date          
       endif
    endif

    !LEAF OFF: DROUGHT DECIDUOUS DRYNESS - if the soil gets too dry, and the leaves have already been on a while... 
    if (currentSite%dstatus == 2.and.t >= 10)then  !D*
       if (sum(currentSite%water_memory(1:10)/10._r8) <= drought_threshold)then 
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
    type(ed_site_type), intent(inout), pointer:: currentSite
    !
    ! !LOCAL VARIABLES:
    type(ed_patch_type) , pointer :: currentPatch     
    type(ed_cohort_type), pointer :: currentCohort  
    !------------------------------------------------------------------------

    currentPatch => CurrentSite%oldest_patch   

    do while(associated(currentPatch))    
       currentCohort => currentPatch%tallest
       do while(associated(currentCohort))        
                
          !COLD LEAF ON
          if (pftcon%season_decid(currentCohort%pft) == 1)then
             if (currentSite%status == 2)then !we have just moved to leaves being on . 
                if (currentCohort%status_coh == 1)then !Are the leaves currently off?        
                   currentCohort%status_coh = 2    !Leaves are on, so change status to stop flow of carbon out of bstore. 
                   if (currentCohort%laimemory <= currentCohort%bstore)then
                      currentCohort%bl = currentCohort%laimemory !extract stored carbon to make new leaves.
                   else
                      currentCohort%bl = currentCohort%bstore    !we can only put on as much carbon as there is in the store...
                    !nb. Putting all of bstore into leaves is C-starvation suicidal. The tendency for this could be parameterized
                   endif
                   currentCohort%balive = currentCohort%balive + currentCohort%bl  ! Add deployed carbon to alive biomass pool
                   currentCohort%bstore = currentCohort%bstore - currentCohort%bl  ! Drain store
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
          if (pftcon%stress_decid(currentCohort%pft) == 1)then
             if (currentSite%dstatus == 2)then !we have just moved to leaves being on . 
                if (currentCohort%status_coh == 1)then !is it the leaf-on day? Are the leaves currently off?       
                   currentCohort%status_coh = 2    !Leaves are on, so change status to stop flow of carbon out of bstore. 
                   if (currentCohort%laimemory <= currentCohort%bstore)then
                      currentCohort%bl = currentCohort%laimemory !extract stored carbon to make new leaves.
                   else
                      currentCohort%bl = currentCohort%bstore !we can only put on as much carbon as there is in the store...
                   endif
                   currentCohort%balive = currentCohort%balive + currentCohort%bl
                   currentCohort%bstore = currentCohort%bstore - currentCohort%bl ! empty store
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
  subroutine seeds_in( cp_pnt )
    !
    ! !DESCRIPTION:
    !  Flux from plants into seed pool. 
    !
    ! !USES:
    !
    ! !ARGUMENTS    
    type(ed_patch_type), intent(inout), target :: cp_pnt ! seeds go to these patches.
    !
    ! !LOCAL VARIABLES:
    type(ed_patch_type),  pointer :: currentPatch
    type(ed_site_type),   pointer :: currentSite
    type(ed_cohort_type), pointer :: currentCohort
    integer :: p
    !----------------------------------------------------------------------

    currentPatch => cp_pnt
    currentSite  => currentPatch%siteptr
   
    currentPatch%seeds_in(:) = 0.0_r8
    currentCohort => currentPatch%tallest
    do while (associated(currentCohort))
       p = currentCohort%pft
       currentPatch%seeds_in(p) = currentPatch%seeds_in(p) +  currentCohort%seed_prod * currentCohort%n/currentPatch%area
       currentCohort => currentCohort%shorter
    enddo !cohort loop

    currentPatch => currentSite%oldest_patch

    do while(associated(currentPatch))
       if (EXTERNAL_RECRUITMENT == 1) then !external seed rain - needed to prevent extinction  
          do p = 1,numpft_ed
           currentPatch%seeds_in(p) = currentPatch%seeds_in(p) + EDecophyscon%seed_rain(p) !KgC/m2/year
          enddo
       endif
       currentPatch => currentPatch%younger
    enddo

  end subroutine seeds_in
  
  ! ============================================================================
  subroutine seed_decay( currentPatch )
    !
    ! !DESCRIPTION:
    !  Flux from seed pool into leaf litter pool    
    !
    ! !USES:
    !
    ! !ARGUMENTS    
    type(ed_patch_type),intent(inout) :: currentPatch ! seeds go to these patches.
    !
    ! !LOCAL VARIABLES:
    integer  ::  p
    real(r8) :: seed_turnover !complete seed turnover rate in yr-1. 
    !----------------------------------------------------------------------

    seed_turnover = 0.51_r8  ! from Liscke and Loffler 2006  
    ! decays the seed pool according to exponential model
    ! sd_mort is in yr-1
    do p = 1,numpft_ed 
       currentPatch%seed_decay(p) =  currentPatch%seed_bank(p) * seed_turnover
    enddo
 
  end subroutine seed_decay

  ! ============================================================================
  subroutine seed_germination( currentPatch ) 
    !
    ! !DESCRIPTION:
    !  Flux from seed pool into sapling pool    
    !
    ! !USES:
    !
    ! !ARGUMENTS    
    type(ed_patch_type),intent(inout) :: currentPatch ! seeds go to these patches.
    !
    ! !LOCAL VARIABLES:
    integer :: p
    real(r8) max_germination !cap on germination rates. KgC/m2/yr Lishcke et al. 2009
    real(r8) germination_timescale !yr-1
    !----------------------------------------------------------------------

    germination_timescale = 0.5_r8 !this is arbitrary
    max_germination = 1.0_r8 !this is arbitrary

    do p = 1,numpft_ed
       currentPatch%seed_germination(p) =  min(currentPatch%seed_bank(p) * germination_timescale,max_germination)
    enddo

  end subroutine seed_germination

  ! ============================================================================
  subroutine Growth_Derivatives( currentCohort)
    !
    ! !DESCRIPTION:
    !  Main subroutine controlling growth and allocation derivatives    
    !
    ! !USES:
    use EDGrowthFunctionsMod , only : Bleaf, dDbhdBd, dhdbd, hite, mortality_rates,dDbhdBl
    use EDTypesMod           , only : udata
    !
    ! !ARGUMENTS    
    type(ed_cohort_type),intent(inout), target :: currentCohort
    !
    ! !LOCAL VARIABLES:
    type(ed_site_type),  pointer :: currentSite
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
    real(r8) :: balive_loss
    !----------------------------------------------------------------------

    currentSite => currentCohort%siteptr

    ! Mortality for trees in the understorey. 
    !if trees are in the canopy, then their death is 'disturbance'. This probably needs a different terminology
    if (currentCohort%canopy_layer > 1)then 
       currentCohort%dndt = -1.0_r8 * mortality_rates(currentCohort) * currentCohort%n
    else
       currentCohort%dndt = 0._r8
    endif

    ! Height
    currentCohort%hite = Hite(currentCohort) 
    h = currentCohort%hite
                       
    call allocate_live_biomass(currentCohort)

   ! calculate target size of living biomass compartment for a given dbh.   
    target_balive = Bleaf(currentCohort) * (1.0_r8 + pftcon%froot_leaf(currentCohort%pft) + &
         EDecophyscon%sapwood_ratio(currentCohort%pft)*h)
    !target balive without leaves. 
    if (currentCohort%status_coh == 1)then 
       target_balive = Bleaf(currentCohort) * (pftcon%froot_leaf(currentCohort%pft) + &
            EDecophyscon%sapwood_ratio(currentCohort%pft) * h)
    endif

    ! NPP 
    currentCohort%npp  = currentCohort%npp_acc  * N_SUB   !Link to CLM. convert from kgC/indiv/day into kgC/indiv/year
    currentCohort%gpp  = currentCohort%gpp_acc  * N_SUB   !Link to CLM. convert from kgC/indiv/day into kgC/indiv/year
    currentCohort%resp = currentCohort%resp_acc * N_SUB   !Link to CLM. convert from kgC/indiv/day into kgC/indiv/year

    currentSite%flux_in = currentSite%flux_in + currentCohort%npp_acc * currentCohort%n

    ! Maintenance demands     
    if (pftcon%evergreen(currentCohort%pft) == 1)then !grass and EBT
       currentCohort%leaf_md = currentCohort%bl / pftcon%leaf_long(currentCohort%pft)
       currentCohort%root_md = currentCohort%br / EDecophyscon%root_long(currentCohort%pft)
       currentCohort%md      = currentCohort%root_md + currentCohort%leaf_md
    endif

    !FIX(RF,032414) - I took out the stem turnover demand as it seemed excesively high and caused odd size-reated 
    ! decline affect
    !with which I am not especially comfortable, particularly as the concept of sapwood turnover is unclear for trees that 
    !are still in an expansion phase. 

    if (pftcon%season_decid(currentCohort%pft) == 1)then 
       currentCohort%root_md = currentCohort%br /EDecophyscon%root_long(currentCohort%pft)
       currentCohort%leaf_md = 0._r8
       currentCohort%md = currentCohort%root_md + currentCohort%leaf_md
    endif

    if (pftcon%stress_decid(currentCohort%pft) == 1)then 
       currentCohort%root_md = currentCohort%br /EDecophyscon%root_long(currentCohort%pft)
       currentCohort%leaf_md = 0._r8
       currentCohort%md = currentCohort%root_md + currentCohort%leaf_md
    endif

    if (pftcon%stress_decid(currentCohort%pft) /= 1.and.pftcon%season_decid(currentCohort%pft) /= 1.and. &
         pftcon%evergreen(currentCohort%pft) /= 1)then
       write(iulog,*) 'problem with phenology definitions',currentCohort%pft,pftcon%stress_decid(currentCohort%pft), &
            pftcon%season_decid(currentCohort%pft),pftcon%evergreen(currentCohort%pft)
    endif

    ! FIX(RF,032414) -turned off for now as it makes balive go negative....
    ! FIX(RF,032414) jan2012 0.01_r8 * currentCohort%bdead
    currentCohort%woody_turnover = 0.0_r8
    currentCohort%md = currentCohort%md + currentCohort%woody_turnover

    ! Calculate carbon balance 
    ! this is the fraction of maintenance demand we -have- to do...

    currentCohort%carbon_balance = currentCohort%npp - currentCohort%md *  EDecophyscon%leaf_stor_priority(currentCohort%pft)

    if (Bleaf(currentCohort) > 0._r8)then

       if (currentCohort%carbon_balance > 0._r8)then !spend C on growing and storing

          !what fraction of the target storage do we have? 
          frac = max(0.0_r8,currentCohort%bstore/(Bleaf(currentCohort) * EDecophyscon%cushion(currentCohort%pft)))
          ! FIX(SPM,080514,fstore never used ) 
          f_store = max(exp(-1.*frac**4._r8) - exp( -1.0_r8 ),0.0_r8)  
          !what fraction of allocation do we divert to storage?
          !what is the flux into the store?
          currentCohort%storage_flux = currentCohort%carbon_balance * f_store                     
          !what is the tax on the carbon available for growth? 
          currentCohort%carbon_balance = currentCohort%carbon_balance * (1.0_r8 - f_store)  
       else  !cbalance is negative. Take C out of store to pay for maintenance respn.
          currentCohort%storage_flux = currentCohort%carbon_balance 
          currentCohort%carbon_balance = 0._r8 
       endif

    else

       currentCohort%storage_flux = 0._r8
       currentCohort%carbon_balance = 0._r8
       write(iulog,*) 'ED: no leaf area in gd', currentCohort%indexnumber,currentCohort%n,currentCohort%bdead, &
             currentCohort%dbh,currentCohort%balive

    endif

    !Do we have enough carbon left over to make up the rest of the turnover demand? 
    balive_loss = 0._r8
    if (currentCohort%carbon_balance > currentCohort%md*(1.0_r8- EDecophyscon%leaf_stor_priority(currentCohort%pft)))then ! Yes...
       currentCohort%carbon_balance = currentCohort%carbon_balance - currentCohort%md * (1.0_r8 - &
             EDecophyscon%leaf_stor_priority(currentCohort%pft))
    else ! we can't maintain constant leaf area and root area. Balive is reduced
       balive_loss = currentCohort%md *(1.0_r8- EDecophyscon%leaf_stor_priority(currentCohort%pft))- currentCohort%carbon_balance
       currentCohort%carbon_balance = 0._r8
    endif

    !********************************************/
    ! Allometry & allocation of remaining carbon*/
    !********************************************/
    !Use remaining carbon to refill balive or to get larger. 

    !only if carbon balance is +ve
    if ((currentCohort%balive >= target_balive).AND.(currentCohort%carbon_balance >  0._r8))then 
       ! fraction of carbon going into active vs structural carbon        
       if (currentCohort%dbh <= EDecophyscon%max_dbh(currentCohort%pft))then ! cap on leaf biomass
          dbldbd = dDbhdBd(currentCohort)/dDbhdBl(currentCohort) 
          dbrdbd = pftcon%froot_leaf(currentCohort%pft) * dbldbd
          dhdbd_fn = dhdbd(currentCohort)
          dbswdbd = EDecophyscon%sapwood_ratio(currentCohort%pft) * (h*dbldbd + currentCohort%bl*dhdbd_fn)
          u  = 1.0_r8 / (dbldbd + dbrdbd + dbswdbd)     
          va = 1.0_r8 / (1.0_r8 + u)
          vs = u / (1.0_r8 + u)
          gr_fract = 1.0_r8 - EDecophyscon%seed_alloc(currentCohort%pft)
       else
          dbldbd = 0._r8; dbrdbd = 0._r8 ;dbswdbd = 0._r8      
          va = 0.0_r8
          vs = 1.0_r8
          gr_fract = 1.0_r8 - (EDecophyscon%seed_alloc(currentCohort%pft) + EDecophyscon%clone_alloc(currentCohort%pft))
       endif

       !FIX(RF,032414) - to fix high bl's. needed to prevent numerical errors without the ODEINT.  
       if (currentCohort%balive > target_balive*1.1_r8)then  
          va = 0.0_r8; vs = 1._r8
          write(iulog,*) 'using high bl cap',target_balive,currentCohort%balive                        
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
    currentCohort%seed_prod = (1.0_r8 - gr_fract) * currentCohort%carbon_balance
    if (abs(currentCohort%npp-(currentCohort%dbalivedt+currentCohort%dbdeaddt+currentCohort%dbstoredt+ &
         currentCohort%seed_prod+currentCohort%md)) > 0.0000000001_r8)then
       write(iulog,*) 'error in carbon check growth derivs',currentCohort%npp- &
            (currentCohort%dbalivedt+currentCohort%dbdeaddt+currentCohort%dbstoredt+currentCohort%seed_prod+currentCohort%md)
       write(iulog,*) 'cohort fluxes',currentCohort%pft,currentCohort%canopy_layer,currentCohort%n, &
            currentCohort%npp,currentCohort%dbalivedt,balive_loss, &
            currentCohort%dbdeaddt,currentCohort%dbstoredt,currentCohort%seed_prod,currentCohort%md * &
            EDecophyscon%leaf_stor_priority(currentCohort%pft)
       write(iulog,*) 'proxies' ,target_balive,currentCohort%balive,currentCohort%dbh,va,vs,gr_fract
    endif

    ! prevent negative leaf pool (but not negative store pool). This is also a numerical error prevention, 
    ! but it shouldn't happen actually... 
    if (-1.0_r8*currentCohort%dbalivedt * udata%deltat > currentCohort%balive*0.99)then 
       write(iulog,*) 'using non-neg leaf mass cap',currentCohort%balive , currentCohort%dbalivedt,currentCohort%dbstoredt, &
            currentCohort%carbon_balance
       currentCohort%dbstoredt = currentCohort%dbstoredt + currentCohort%dbalivedt
       currentCohort%dbalivedt = 0._r8 
    endif

    ! calculate change in diameter and height 
    currentCohort%ddbhdt = currentCohort%dbdeaddt * dDbhdBd(currentCohort)
    currentCohort%dhdt   = currentCohort%dbdeaddt * dHdBd(currentCohort)

  end subroutine Growth_Derivatives

  ! ============================================================================
  subroutine recruitment( t, currentPatch )
    !
    ! !DESCRIPTION:
    ! spawn new cohorts of juveniles of each PFT             
    !
    ! !USES:
    use EDGrowthFunctionsMod, only : bdead,dbh, Bleaf
    use EDTypesMod, only : udata
    !
    ! !ARGUMENTS    
    integer, intent(in) :: t
    type(ed_patch_type), intent(inout), pointer :: currentPatch
    !
    ! !LOCAL VARIABLES:
    integer :: ft
    type (ed_cohort_type) , pointer :: temp_cohort
    integer :: cohortstatus
    !----------------------------------------------------------------------

    allocate(temp_cohort) ! create temporary cohort
    call zero_cohort(temp_cohort)

    do ft = 1,numpft_ed

       temp_cohort%canopy_trim = 0.8_r8  !starting with the canopy not fully expanded 
       temp_cohort%pft         = ft
       temp_cohort%hite        = EDecophyscon%hgt_min(ft)
       temp_cohort%dbh         = Dbh(temp_cohort)
       temp_cohort%bdead       = Bdead(temp_cohort)
       temp_cohort%balive      = Bleaf(temp_cohort)*(1.0_r8 + pftcon%froot_leaf(ft) &
            + EDecophyscon%sapwood_ratio(ft)*temp_cohort%hite)
       temp_cohort%bstore      = EDecophyscon%cushion(ft)*(temp_cohort%balive/ (1.0_r8 + pftcon%froot_leaf(ft) &
            + EDecophyscon%sapwood_ratio(ft)*temp_cohort%hite))
       temp_cohort%n           = currentPatch%area * currentPatch%seed_germination(ft)*udata%deltat &
            / (temp_cohort%bdead+temp_cohort%balive+temp_cohort%bstore)
 
       if (t == 1)then
          write(iulog,*) 'filling in cohorts where there are none left; this will break carbon balance', &
               currentPatch%patchno,currentPatch%area
          temp_cohort%n = 0.1_r8*currentPatch%area
          write(iulog,*) 'cohort n',ft,temp_cohort%n
       endif

       temp_cohort%laimemory = 0.0_r8     
       if (pftcon%season_decid(temp_cohort%pft) == 1.and.currentPatch%siteptr%status == 1)then
         temp_cohort%laimemory = (1.0_r8/(1.0_r8 + pftcon%froot_leaf(ft) + &
              EDecophyscon%sapwood_ratio(ft)*temp_cohort%hite))*temp_cohort%balive
       endif
       if (pftcon%stress_decid(temp_cohort%pft) == 1.and.currentPatch%siteptr%dstatus == 1)then
         temp_cohort%laimemory = (1.0_r8/(1.0_r8 + pftcon%froot_leaf(ft) + &
            EDecophyscon%sapwood_ratio(ft)*temp_cohort%hite))*temp_cohort%balive
       endif

       cohortstatus = currentPatch%siteptr%status
       if (pftcon%stress_decid(ft) == 1)then !drought decidous, override status. 
          cohortstatus = currentPatch%siteptr%dstatus
       endif

       if (temp_cohort%n > 0.0_r8)then
          call create_cohort(currentPatch, temp_cohort%pft, temp_cohort%n, temp_cohort%hite, temp_cohort%dbh, &
               temp_cohort%balive, temp_cohort%bdead, temp_cohort%bstore,  &
               temp_cohort%laimemory, cohortstatus, temp_cohort%canopy_trim, currentPatch%NCL_p)
       endif
    enddo  !pft loop

    deallocate(temp_cohort) ! delete temporary cohort

    call fuse_cohorts(currentPatch)
    call sort_cohorts(currentPatch)

  end subroutine recruitment

  ! ============================================================================
  subroutine CWD_Input( currentPatch)
    !
    ! !DESCRIPTION:
    ! Generate litter fields from turnover.  
    !
    ! !USES:
    use SFParamsMod , only : SF_val_CWD_frac
    use EDParamsMod , only : ED_val_ag_biomass
    use EDTypesMod  , only : udata
    !
    ! !ARGUMENTS    
    type(ed_patch_type),intent(inout), target :: currentPatch
    !
    ! !LOCAL VARIABLES:
    type(ed_cohort_type), pointer :: currentCohort
    integer  :: c,p
    real(r8) :: not_dead_n !projected remaining number of trees in understorey cohort after turnover
    real(r8) :: dead_n !understorey dead tree density
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
         currentCohort%leaf_litter * currentCohort%n/currentPatch%area/udata%deltat

      !daily leaf loss needs to be scaled up to the annual scale here. 
      
      do c = 1,ncwd
         currentPatch%cwd_AG_in(c) = currentPatch%cwd_AG_in(c) + currentCohort%woody_turnover * &
              SF_val_CWD_frac(c) * currentCohort%n/currentPatch%area *ED_val_ag_biomass
         currentPatch%cwd_BG_in(c) = currentPatch%cwd_BG_in(c) + currentCohort%woody_turnover * &
              SF_val_CWD_frac(c) * currentCohort%n/currentPatch%area *(1.0_r8-ED_val_ag_biomass)
      enddo

      if (currentCohort%canopy_layer > 1)then   

          ! ================================================        
          ! Litter fluxes for understorey  mortality. KgC/m2/year
          ! ================================================
          dead_n = -1.0_r8 * currentCohort%dndt / currentPatch%area

          currentPatch%leaf_litter_in(pft) = currentPatch%leaf_litter_in(pft) + &
               (currentCohort%bl+currentCohort%leaf_litter/udata%deltat)* dead_n          
          currentPatch%root_litter_in(pft) = currentPatch%root_litter_in(pft) + &
               (currentCohort%br+currentCohort%bstore)     * dead_n

          do c = 1,ncwd
             currentPatch%cwd_AG_in(c) = currentPatch%cwd_AG_in(c) + (currentCohort%bdead+currentCohort%bsw) * &
                   SF_val_CWD_frac(c) * dead_n * ED_val_ag_biomass
             currentPatch%cwd_BG_in(c) = currentPatch%cwd_BG_in(c) + (currentCohort%bdead+currentCohort%bsw) * &
                  SF_val_CWD_frac(c) * dead_n * (1.0_r8-ED_val_ag_biomass)

             if (currentPatch%cwd_AG_in(c) < 0.0_r8)then
                write(iulog,*) 'negative CWD in flux',currentPatch%cwd_AG_in(c), &
                     (currentCohort%bdead+currentCohort%bsw), dead_n
             endif
          enddo

       endif !canopy layer

       currentCohort => currentCohort%taller

    enddo  ! end loop over cohorts 

    do p = 1,numpft_ed
       currentPatch%leaf_litter_in(p) = currentPatch%leaf_litter_in(p) + currentPatch%seed_decay(p) !KgC/m2/yr
    enddo

  end subroutine CWD_Input

  ! ============================================================================
  subroutine fragmentation_scaler( currentPatch, temperature_inst )
    !
    ! !DESCRIPTION:
    ! Simple CWD fragmentation Model
    ! FIX(SPM, 091914) this should be a function as it returns a value in currentPatch%fragmentation_scaler
    !
    ! !USES:
    use shr_const_mod      , only : SHR_CONST_PI, SHR_CONST_TKFRZ
    use EDSharedParamsMod  , only : EDParamsShareInst
    !
    ! !ARGUMENTS    
    type(ed_patch_type)    , intent(inout) :: currentPatch
    type(temperature_type) , intent(in)    :: temperature_inst
    !
    ! !LOCAL VARIABLES:
    logical  :: use_century_tfunc = .false.
    type(ed_site_type), pointer :: currentSite
    integer  :: c,p,j
    real(r8) :: t_scalar
    real(r8) :: w_scalar
    real(r8) :: catanf                ! hyperbolic temperature function from CENTURY
    real(r8) :: catanf_30             ! hyperbolic temperature function from CENTURY
    real(r8) :: t1                    ! temperature argument
    real(r8) :: Q10                   ! temperature dependence
    real(r8) :: froz_q10              ! separate q10 for frozen soil respiration rates.  default to same as above zero rates
    real(r8), pointer :: t_veg24(:)
    !----------------------------------------------------------------------

    catanf(t1) = 11.75_r8 +(29.7_r8 / SHR_CONST_PI) * atan( SHR_CONST_PI * 0.031_r8  * ( t1 - 15.4_r8 ))

    t_veg24 => temperature_inst%t_veg24_patch      ! Input:  [real(r8) (:)]  avg pft vegetation temperature for last 24 hrs

    catanf_30 = catanf(30._r8)
    
    c = currentPatch%siteptr%clmcolumn
    p = currentPatch%clm_pno
    
    ! set "froz_q10" parameter
    froz_q10  = EDParamsShareInst%froz_q10  
    Q10       = EDParamsShareInst%Q10

    if ( .not. use_century_tfunc ) then
    !calculate rate constant scalar for soil temperature,assuming that the base rate constants 
    !are assigned for non-moisture limiting conditions at 25C. 
      if (t_veg24(p)  >=  SHR_CONST_TKFRZ) then
        t_scalar = Q10**((t_veg24(p)-(SHR_CONST_TKFRZ+25._r8))/10._r8)
                 !  Q10**((t_soisno(c,j)-(SHR_CONST_TKFRZ+25._r8))/10._r8)
      else
        t_scalar = (Q10**(-25._r8/10._r8))*(froz_q10**((t_veg24(p)-SHR_CONST_TKFRZ)/10._r8))
                  !Q10**(-25._r8/10._r8))*(froz_q10**((t_soisno(c,j)-SHR_CONST_TKFRZ)/10._r8)
      endif
    else
      ! original century uses an arctangent function to calculate the temperature dependence of decomposition      
      t_scalar = max(catanf(t_veg24(p)-SHR_CONST_TKFRZ)/catanf_30,0.01_r8)
    endif    
   
    !Moisture Limitations   
    !BTRAN APPROACH - is quite simple, but max's out decomp at all unstressed soil moisture values, which is not realistic.  
    !litter decomp is proportional to water limitation on average... 
    w_scalar = sum(currentPatch%btran_ft(1:numpft_ed))/numpft_ed 

    currentPatch%fragmentation_scaler =  min(1.0_r8,max(0.0_r8,t_scalar * w_scalar))
    
  end subroutine fragmentation_scaler
  
  ! ============================================================================
  subroutine cwd_out( currentPatch, temperature_inst, soilstate_inst, waterstate_inst)
    !
    ! !DESCRIPTION:
    ! Simple CWD fragmentation Model
    ! spawn new cohorts of juveniles of each PFT             
    !
    ! !USES:
    use SFParamsMod, only : SF_val_max_decomp
    use EDTypesMod , only : udata
    !
    ! !ARGUMENTS    
    type(ed_patch_type)    , intent(inout), target :: currentPatch
    type(temperature_type) , intent(in)            :: temperature_inst
    type(soilstate_type)   , intent(in)            :: soilstate_inst
    type(waterstate_type)  , intent(in)            :: waterstate_inst
    !
    ! !LOCAL VARIABLES:
    type(ed_site_type), pointer :: currentSite
    integer :: c,ft
    !----------------------------------------------------------------------

    currentSite => currentPatch%siteptr
    currentPatch%root_litter_out = 0.0_r8
    currentPatch%leaf_litter_out = 0.0_r8

    call fragmentation_scaler(currentPatch, temperature_inst)

    !Flux of coarse woody debris into decomposing litter pool. 

    currentPatch%cwd_ag_out(1:ncwd) = 0.0_r8
    currentPatch%cwd_bg_out(1:ncwd) = 0.0_r8
    currentPatch%leaf_litter_out(1:numpft_ed) = 0.0_r8
    currentPatch%root_litter_out(1:numpft_ed) = 0.0_r8
    
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

    do ft = 1,numpft_ed
       currentPatch%leaf_litter_out(ft) = max(0.0_r8,currentPatch%leaf_litter(ft)* SF_val_max_decomp(dg_sf) * &
            currentPatch%fragmentation_scaler )
       currentPatch%root_litter_out(ft) = max(0.0_r8,currentPatch%root_litter(ft)* SF_val_max_decomp(dg_sf) * &
            currentPatch%fragmentation_scaler )
       if ( currentPatch%leaf_litter_out(ft)<0.0_r8.or.currentPatch%root_litter_out(ft)<0.0_r8)then
         write(iulog,*) 'root or leaf out is negative?',SF_val_max_decomp(dg_sf),currentPatch%fragmentation_scaler
       endif
    enddo

    !add up carbon going into fragmenting pools
    currentSite%flux_out = currentSite%flux_out + sum(currentPatch%leaf_litter_out) * &
         currentPatch%area *udata%deltat!kgC/site/day
    currentSite%flux_out = currentSite%flux_out + sum(currentPatch%root_litter_out) * &
         currentPatch%area *udata%deltat!kgC/site/day
    currentSite%flux_out = currentSite%flux_out + sum(currentPatch%cwd_ag_out) * &
         currentPatch%area *udata%deltat!kgC/site/day
    currentSite%flux_out = currentSite%flux_out + sum(currentPatch%cwd_bg_out) * &
         currentPatch%area *udata%deltat!kgC/site/day

  end subroutine cwd_out

end module EDPhysiologyMod
