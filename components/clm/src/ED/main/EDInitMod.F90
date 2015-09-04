module EDInitMod

  ! ============================================================================
  ! Contains all modules to set up the ED structure. 
  ! ============================================================================

  use spmdMod               , only : masterproc
  use decompMod             , only : bounds_type
  use clm_varpar            , only : nclmax
  use clm_varctl            , only : iulog 
  use CNCarbonFluxType      , only : carbonflux_type
  use CNCarbonStateType     , only : carbonstate_type
  use CNNitrogenStateType   , only : nitrogenstate_type
  use CanopyStateType       , only : canopystate_type
  use WaterStateType        , only : waterstate_type
  use GridcellType          , only : grc
  use EcophysConType        , only : ecophyscon
  use EDBioType             , only : EDbio_type
  use EDEcophysConType      , only : EDecophyscon
  use EDGrowthFunctionsMod  , only : bdead, bleaf, dbh
  use EDCohortDynamicsMod   , only : create_cohort, fuse_cohorts, sort_cohorts
  use EDPatchDynamicsMod    , only : create_patch
  use EDCLMLinkMod          , only : clm_ed_link, clm_indices
  use EDMainMod             , only : ed_update_sites
  use EDtypesMod            , only : site, patch, cohort, gridcell_edstate_type, area
  use EDtypesMod            , only : cohorts_per_gcell, gridCellEdState, ncwd, numpft_ed, udata

  implicit none
  save
  private

  public :: ed_init
  public :: ed_init_sites
  public :: allocate_ed_datastructure
  public :: zero_site
  public :: set_site_properties_bg
  public :: init_patches
  public :: init_cohorts

  ! ============================================================================
  ! ============================================================================

contains

  ! ============================================================================
  subroutine ed_init( bounds, waterstate_vars, canopystate_vars, EDbio_vars, & 
       carbonstate_vars, nitrogenstate_vars, carbonflux_vars) 
     !
     ! use gridCellEdState at the top level, then pass it through arg. list.  then we can
     ! actually use intents
     !
     use clm_varctl,       only : use_ed, use_ed_spit_fire
     use clm_time_manager, only : is_restart

     implicit none
     type(bounds_type)        , intent(in)    :: bounds  ! clump bounds
     type(waterstate_type)    , intent(inout) :: waterstate_vars
     type(canopystate_type)   , intent(inout) :: canopystate_vars
     type(EDbio_type)         , intent(inout) :: EDbio_vars
     type(carbonstate_type)   , intent(inout) :: carbonstate_vars
     type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
     type(carbonflux_type)    , intent(inout) :: carbonflux_vars

     if (masterproc) then
        write(iulog,*) 'ED: restart ? = ' ,is_restart()                                     ! FIX(SPM,032414) debug
        write(iulog,*) 'ED_Mod.F90 :: use_ed ',use_ed                                       ! FIX(SPM,032414) debug
        write(iulog,*) 'ED_Mod.F90 :: SPITFIRE_SWITCH (use_ed_spit_fire) ',use_ed_spit_fire ! FIX(SPM,032414) debug
        write(iulog,*) 'ED_Mod.F90 :: cohorts_per_gcell ',cohorts_per_gcell                 ! FIX(SPM,032414) debug
     end if

     if ( .not. is_restart() ) then
        ! 
        ! allocate gridcell edstate type structure at the top level for startup runs
        ! on restart this is allocated in EDRest
        !
        allocate (gridCellEdState(bounds%begg:bounds%endg))

        call ed_init_sites( bounds, gridCellEdState )

        call ed_update_sites( bounds, gridCellEdState )

        call clm_ed_link( bounds, gridCellEdState, waterstate_vars, canopystate_vars, EDbio_vars, &
             carbonstate_vars, nitrogenstate_vars, carbonflux_vars) 

     endif

  end subroutine ed_init

  ! ============================================================================
  subroutine ed_init_sites( bounds, geds_local )

     use clm_time_manager , only : is_restart

     implicit none

     type(bounds_type),                      intent(in) :: bounds
     type(gridcell_edstate_type), target, intent(inout) :: geds_local( bounds%begg: )

     call allocate_ed_datastructure( bounds, geds_local )
     !set the 'istheresoil' value before going into update...
     call clm_indices( bounds, geds_local )
     call set_site_properties_bg( bounds, geds_local )     
     !
     ! on restart, this functionality is handled in EDRestVectorMod::createPatchCohortStructure
     !
     if (.not. is_restart() ) then
        call init_patches( bounds, geds_local )
     endif

  end subroutine ed_init_sites

  ! ============================================================================
  subroutine allocate_ed_datastructure( bounds, geds_local )

    implicit none

    type(bounds_type),              intent(in) :: bounds
    type(gridcell_edstate_type), intent(inout) :: geds_local( bounds%begg: )

    type(site), pointer :: currentSite
    type(site), pointer :: store_fs
    integer :: first_site_flag
    integer :: g

    ! INITIALISE THE SITE STRUCTURES
    first_site_flag = 1
    udata%currentindex = 0 !Makes unique cohort identifiers. Needs zeroing at beginning of run. 

    do g = bounds%begg,bounds%endg
       allocate(currentSite) 
       call zero_site(currentSite)
       !create clm mapping to ED structure
       geds_local(g)%spnt => currentSite
       currentSite%clmgcell = g 
       currentSite%lat = grc%latdeg(g)  
       currentSite%lon = grc%londeg(g)

       ! Join sites together in linked list
       if(first_site_flag == 1) then !if first site do this, otherwise link it
          udata%firstsite_pnt => currentSite
          store_fs => currentSite
          first_site_flag = 0
       endif

    enddo

    udata%firstsite_pnt => store_fs

  end subroutine allocate_ed_datastructure

  ! ============================================================================
  subroutine zero_site( site_in )

    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use shr_kind_mod   , only : r8 => shr_kind_r8

    implicit none

    type(site), intent(inout) ::  site_in

    site_in%oldest_patch     => null() ! pointer to oldest patch at the site
    site_in%youngest_patch   => null() ! pointer to yngest patch at the site

    ! INDICES 
    site_in%lat              = nan
    site_in%lon              = nan
    site_in%clmgcell         = 0
    site_in%clmcolumn        = 0
    site_in%istheresoil      = 0

    ! DISTURBANCE
    site_in%disturbance_rate = 0._r8  ! site level disturbance rates from mortality and fire.
    site_in%dist_type        = 0      ! disturbance dist_type id.

    ! PHENOLOGY 
    site_in%status           = 0    ! are leaves in this pixel on or off?
    site_in%dstatus          = 0
    site_in%gdd              = nan  ! growing degree days
    site_in%ncd              = nan  ! no chilling days
    site_in%last_n_days(:)   = 999  ! record of last 10 days temperature for senescence model.
    site_in%leafondate       = 999  ! doy of leaf on
    site_in%leafoffdate      = 999  ! doy of leaf off
    site_in%dleafondate      = 999  ! doy of leaf on drought
    site_in%dleafoffdate     = 999  ! doy of leaf on drought
    site_in%water_memory(:)  = nan

    ! FIRE 
    site_in%acc_ni           = 0.0_r8     ! daily nesterov index accumulating over time. time unlimited theoretically.
    site_in%frac_burnt       = 0.0_r8     ! burn area read in from external file
 
  end subroutine zero_site

  ! ============================================================================
  subroutine set_site_properties_bg( bounds, geds_local )

     use shr_kind_mod     , only : r8 => shr_kind_r8;
     use clm_time_manager , only : is_restart

     implicit none

     type(bounds_type),                      intent(in) :: bounds  ! clump bounds
     type(gridcell_edstate_type), target, intent(inout) :: geds_local( bounds%begg: )

     type(site), pointer :: currentSite

     integer  :: i,g !beginning and end of these data clumps.

     real(r8) :: leafon(bounds%begg:bounds%endg)
     real(r8) :: leafoff(bounds%begg:bounds%endg)
     real(r8) :: stat(bounds%begg:bounds%endg)
     real(r8) :: NCD(bounds%begg:bounds%endg)
     real(r8) :: GDD(bounds%begg:bounds%endg)
     real(r8) :: dstat(bounds%begg:bounds%endg)
     real(r8) :: acc_NI(bounds%begg:bounds%endg)
     real(r8) :: watermem(bounds%begg:bounds%endg)

     integer  :: dleafoff(bounds%begg:bounds%endg)
     integer  :: dleafon(bounds%begg:bounds%endg)

     if ( .not. is_restart() ) then
        !initial guess numbers for site condition.
        do i                              = bounds%begg,bounds%endg
           NCD(i)                         = 0.0_r8
           GDD(i)                         = 30_r8
           leafon(i)                      = 100_r8
           leafoff(i)                     = 300_r8
           stat(i)                        = 2
           acc_NI(i)                      = 0.0_r8
           dstat(i)                       = 2
           dleafoff(i)                    = 300
           dleafon(i)                     = 100
           watermem(i)                    = 0.5_r8
        enddo
     else ! assignements for restarts
        do i                              = bounds%begg,bounds%endg
           NCD(i)                         = 1.0_r8 ! NCD should be 1 on restart
           !GDD(i)                        = 0.0_r8
           leafon(i)                      = 0.0_r8
           leafoff(i)                     = 0.0_r8
           stat(i)                        = 1
           acc_NI(i)                      = 0.0_r8
           dstat(i)                       = 2
           dleafoff(i)                    = 300
           dleafon(i)                     = 100
           watermem(i)                    = 0.5_r8
        enddo
     endif

     do g                                 = bounds%begg,bounds%endg
        currentSite                       => geds_local(g)%spnt

        currentSite%gdd                   = GDD(g)
        currentSite%ncd                   = NCD(g)
        currentSite%leafondate            = leafon(g)
        currentSite%leafoffdate           = leafoff(g)
        currentSite%dleafoffdate          = dleafoff(g)
        currentSite%dleafondate           = dleafon(g)

        if ( .not. is_restart() ) then
           currentSite%water_memory(1:10) = watermem(g)
        end if

        currentSite%status                = stat(g)
        !start off with leaves off to initialise
        currentSite%dstatus               = dstat(g)

        currentSite%acc_NI                = acc_NI(g)
        currentSite%frac_burnt            = 0.0_r8
        currentSite%old_stock             = 0.0_r8

     enddo

  end subroutine set_site_properties_bg

  ! ============================================================================
  subroutine init_patches( bounds, geds_local )
     !
     !initialize patches on new ground
     !
     use shr_kind_mod,  only : r8 => shr_kind_r8
     use EDParamsMod ,  only : ED_val_maxspread

     implicit none

     type(bounds_type),                      intent(in) :: bounds  ! clump bounds
     type(gridcell_edstate_type), target, intent(inout) :: geds_local( bounds%begg: )

     type(site) , pointer :: currentSite
     type(patch), pointer :: newp

     ! these are local to this routine
     real(r8) cwd_ag_local(ncwd)
     real(r8) cwd_bg_local(ncwd)
     real(r8) spread_local(nclmax)
     real(r8) leaf_litter_local(numpft_ed)
     real(r8) root_litter_local(numpft_ed)
     real(r8) seed_bank_local(numpft_ed)
     real(r8) age !notional age of this patch

     integer :: g

     cwd_ag_local      = 0.0_r8 !ED_val_init_litter -- arbitrary value for litter pools. kgC m-2
     cwd_bg_local      = 0.0_r8 !ED_val_init_litter
     leaf_litter_local = 0.0_r8
     root_litter_local = 0.0_r8
     age               = 0.0_r8
     spread_local      = ED_val_maxspread


     !FIX(SPM,032414) clean this up...inits out of this loop
     do g = bounds%begg,bounds%endg

        currentSite => geds_local(g)%spnt

        allocate(newp)

        currentSite%youngest_patch => newp

        call create_patch(currentSite,newp,age,AREA,spread_local,cwd_ag_local,cwd_bg_local,leaf_litter_local, &
             root_litter_local,seed_bank_local) !make new patch...

        newp%patchno = 1

        call init_cohorts(newp)

        newp%younger => null()
        newp%older   => null()
        currentSite%youngest_patch => newp
        currentSite%oldest_patch   => newp

     enddo !gridcells

  end subroutine init_patches

  subroutine init_cohorts( patch_in )
     ! ============================================================================
     !              initialize new cohorts on bare ground
     ! ============================================================================

     use shr_kind_mod, only : r8 => shr_kind_r8

     implicit none

     type(patch), intent(inout), pointer  :: patch_in

     type(cohort),pointer :: dc

     integer :: cstatus
     integer :: pft

     patch_in%tallest  => null()
     patch_in%shortest => null()

     do pft =  1,numpft_ed !FIX(RF,032414) - turning off veg dynamics

        allocate(dc)

        dc%pft                        = pft
        dc%n                          = EDecophyscon%initd(pft) * patch_in%area
        dc%hite                       = EDecophyscon%hgt_min(pft)
        dc%dbh                        = Dbh(dc) ! FIX(RF, 090314) - comment out addition of ' + 0.0001_r8*pft   '  - seperate out PFTs a little bit...
        dc%canopy_trim                = 1.0_r8
        dc%bdead                      = Bdead(dc)
        dc%balive                     = Bleaf(dc)*(1.0_r8 + ecophyscon%froot_leaf(pft) +EDecophyscon%sapwood_ratio(dc%pft)*dc%hite)
        dc%b                          = dc%balive + dc%bdead

        if( ecophyscon%evergreen(pft)     == 1) then
           dc%bstore                  = Bleaf(dc) * EDecophyscon%cushion(pft)
           dc%laimemory               = 0._r8
           cstatus                    = 2
        endif

        if( ecophyscon%season_decid(pft)  == 1 ) then !for dorment places
           dc%bstore                  = Bleaf(dc) * EDecophyscon%cushion(pft) !stored carbon in new seedlings.
           if(patch_in%siteptr%status==2)then 
             dc%laimemory             = 0.0_r8
           else
             dc%laimemory             = Bleaf(dc)
           endif
           ! reduce biomass according to size of store, this will be recovered when elaves com on.
           dc%balive                  = dc%balive - dc%laimemory
           cstatus                    = patch_in%siteptr%status
        endif

        if ( ecophyscon%stress_decid(pft) == 1 ) then
           dc%bstore                  = Bleaf(dc) * EDecophyscon%cushion(pft)
           dc%laimemory               = Bleaf(dc)
           dc%balive                  = dc%balive - dc%laimemory
           cstatus                    = patch_in%siteptr%dstatus
        endif

        call create_cohort(pft,dc%n,dc%hite,dc%dbh,dc%balive,dc%bdead,dc%bstore,dc%laimemory, &
             cstatus,dc%canopy_trim,1,patch_in)

        deallocate(dc) ! get rid of temporary cohort

        patch_in%tallest              => udata%storebigcohort
        patch_in%shortest             => udata%storesmallcohort

     enddo !numpft

     call fuse_cohorts(patch_in)
     call sort_cohorts(patch_in)

  end subroutine init_cohorts

  ! ============================================================================
end module EDInitMod
