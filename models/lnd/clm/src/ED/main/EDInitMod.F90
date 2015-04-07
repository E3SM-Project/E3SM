module EDInitMod

  ! ============================================================================
  ! Contains all modules to set up the ED structure. 
  ! ============================================================================

  use shr_kind_mod              , only : r8 => shr_kind_r8;
  use spmdMod                   , only : masterproc
  use decompMod                 , only : bounds_type
  use abortutils                , only : endrun
  use clm_varpar                , only : nclmax
  use clm_varctl                , only : iulog, use_ed_spit_fire 
  use clm_time_manager          , only : is_restart
  use CanopyStateType           , only : canopystate_type
  use WaterStateType            , only : waterstate_type
  use GridcellType              , only : grc
  use pftconMod                 , only : pftcon
  use EDPhenologyType           , only : ed_phenology_type
  use EDEcophysConType          , only : EDecophyscon
  use EDGrowthFunctionsMod      , only : bdead, bleaf, dbh
  use EDCohortDynamicsMod       , only : create_cohort, fuse_cohorts, sort_cohorts
  use EDPatchDynamicsMod        , only : create_patch
  use EDMainMod                 , only : ed_update_site
  use EDTypesMod                , only : ed_site_type, ed_patch_type, ed_cohort_type, area
  use EDTypesMod                , only : cohorts_per_gcell, ncwd, numpft_ed, udata
  use EDCLMLinkMod              , only : ed_clm_type

  implicit none
  private

  public  :: ed_init
  public  :: ed_init_sites
  public  :: zero_site

  private :: set_site_properties
  private :: init_patches
  private :: init_cohorts
  ! ============================================================================

contains

  ! ============================================================================
  subroutine ed_init( bounds, ed_allsites_inst, ed_clm_inst, &
       ed_phenology_inst, waterstate_inst, canopystate_inst)
    !
    ! !DESCRIPTION:
    ! use ed_allsites_inst at the top level, then pass it through arg. list.  then we can
    ! actually use intents
    !
    ! !USES: 
    !
    ! !ARGUMENTS    
    type(bounds_type)       , intent(in)            :: bounds  ! clump bounds
    type(ed_site_type)      , intent(inout), target :: ed_allsites_inst( bounds%begg: )
    type(ed_clm_type)       , intent(inout)         :: ed_clm_inst
    type(ed_phenology_type) , intent(inout)         :: ed_phenology_inst
    type(waterstate_type)   , intent(inout)         :: waterstate_inst
    type(canopystate_type)  , intent(inout)         :: canopystate_inst
    !
    ! !LOCAL VARIABLES:
    integer :: g
    !----------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'ED: restart ? = ' ,is_restart()                                     ! FIX(SPM,032414) debug
       write(iulog,*) 'ED_Mod.F90 :: SPITFIRE_SWITCH (use_ed_spit_fire) ',use_ed_spit_fire ! FIX(SPM,032414) debug
       write(iulog,*) 'ED_Mod.F90 :: cohorts_per_gcell ',cohorts_per_gcell                 ! FIX(SPM,032414) debug
    end if

    if ( .not. is_restart() ) then
       call ed_init_sites( bounds, ed_allsites_inst(bounds%begg:bounds%endg))

       do g = bounds%begg,bounds%endg
          if (ed_allsites_inst(g)%istheresoil) then
             call ed_update_site(ed_allsites_inst(g))
          end if
       end do

       call ed_clm_inst%ed_clm_link( bounds, ed_allsites_inst(bounds%begg:bounds%endg), &
            ed_phenology_inst, waterstate_inst, canopystate_inst)
    endif

  end subroutine ed_init

  ! ============================================================================
  subroutine ed_init_sites( bounds, ed_allsites_inst )
    !
    ! !DESCRIPTION:
    ! Intialize all ED sites
    !
    ! !USES: 
    use ColumnType      , only : col
    use landunit_varcon , only : istsoil
    !
    ! !ARGUMENTS    
    type(bounds_type)  , intent(in)            :: bounds
    type(ed_site_type) , intent(inout), target :: ed_allsites_inst( bounds%begg: )
    !
    ! !LOCAL VARIABLES:
    integer  :: g,l,c
    logical  :: istheresoil(bounds%begg:bounds%endg) 
    !----------------------------------------------------------------------

    ! INITIALISE THE SITE STRUCTURES
    udata%cohort_number = 0 !Makes unique cohort identifiers. Needs zeroing at beginning of run. 

    do g = bounds%begg,bounds%endg
       ! zero the site
       call zero_site(ed_allsites_inst(g))

       !create clm mapping to ED structure
       ed_allsites_inst(g)%clmgcell = g 
       ed_allsites_inst(g)%lat      = grc%latdeg(g)  
       ed_allsites_inst(g)%lon      = grc%londeg(g)
    enddo

    istheresoil(bounds%begg:bounds%endg) = .false.
    do c = bounds%begc,bounds%endc
       g = col%gridcell(c)   
       if (col%itype(c) == istsoil) then  
          istheresoil(g) = .true.
       endif
       ed_allsites_inst(g)%istheresoil = istheresoil(g)
    enddo

    call set_site_properties( bounds, ed_allsites_inst(bounds%begg:bounds%endg) )     

    ! on restart, this functionality is handled in EDRestVectorMod::createPatchCohortStructure
    if (.not. is_restart() ) then
       call init_patches( bounds, ed_allsites_inst(bounds%begg:bounds%endg) )
    endif

  end subroutine ed_init_sites

  ! ============================================================================
  subroutine zero_site( site_in )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS    
    type(ed_site_type), intent(inout) ::  site_in
    !
    ! !LOCAL VARIABLES:
    !----------------------------------------------------------------------

    site_in%oldest_patch     => null() ! pointer to oldest patch at the site
    site_in%youngest_patch   => null() ! pointer to yngest patch at the site

    ! INDICES 
    site_in%lat              = nan
    site_in%lon              = nan
    site_in%clmgcell         = 0
    site_in%clmcolumn        = 0
    site_in%istheresoil      = .false.

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
  subroutine set_site_properties( bounds, ed_allsites_inst )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    !
    ! !ARGUMENTS    
    type(bounds_type)  , intent(in)            :: bounds  
    type(ed_site_type) , intent(inout), target :: ed_allsites_inst( bounds%begg: )
    !
    ! !LOCAL VARIABLES:
    integer  :: i,g !beginning and end of these data clumps.
    real(r8) :: leafon   (bounds%begg:bounds%endg)
    real(r8) :: leafoff  (bounds%begg:bounds%endg)
    real(r8) :: stat     (bounds%begg:bounds%endg)
    real(r8) :: NCD      (bounds%begg:bounds%endg)
    real(r8) :: GDD      (bounds%begg:bounds%endg)
    real(r8) :: dstat    (bounds%begg:bounds%endg)
    real(r8) :: acc_NI   (bounds%begg:bounds%endg)
    real(r8) :: watermem (bounds%begg:bounds%endg)
    integer  :: dleafoff (bounds%begg:bounds%endg)
    integer  :: dleafon  (bounds%begg:bounds%endg)
    !----------------------------------------------------------------------

    if ( .not. is_restart() ) then
       !initial guess numbers for site condition.
       do i = bounds%begg,bounds%endg
          NCD(i)      = 0.0_r8
          GDD(i)      = 30_r8
          leafon(i)   = 100_r8
          leafoff(i)  = 300_r8
          stat(i)     = 2
          acc_NI(i)   = 0.0_r8
          dstat(i)    = 2
          dleafoff(i) = 300
          dleafon(i)  = 100
          watermem(i) = 0.5_r8
       enddo
    else ! assignements for restarts
       do i = bounds%begg,bounds%endg
          NCD(i)      = 1.0_r8 ! NCD should be 1 on restart
          !GDD(i)     = 0.0_r8
          leafon(i)   = 0.0_r8
          leafoff(i)  = 0.0_r8
          stat(i)     = 1
          acc_NI(i)   = 0.0_r8
          dstat(i)    = 2
          dleafoff(i) = 300
          dleafon(i)  = 100
          watermem(i) = 0.5_r8
       enddo
    endif

    do g = bounds%begg,bounds%endg
       ed_allsites_inst(g)%gdd          = GDD(g)
       ed_allsites_inst(g)%ncd          = NCD(g)
       ed_allsites_inst(g)%leafondate   = leafon(g)
       ed_allsites_inst(g)%leafoffdate  = leafoff(g)
       ed_allsites_inst(g)%dleafoffdate = dleafoff(g)
       ed_allsites_inst(g)%dleafondate  = dleafon(g)

       if ( .not. is_restart() ) then
          ed_allsites_inst(g)%water_memory(1:10) = watermem(g)
       end if

       ed_allsites_inst(g)%status = stat(g)
       !start off with leaves off to initialise
       ed_allsites_inst(g)%dstatus= dstat(g)

       ed_allsites_inst(g)%acc_NI     = acc_NI(g)
       ed_allsites_inst(g)%frac_burnt = 0.0_r8
       ed_allsites_inst(g)%old_stock  = 0.0_r8
    enddo

  end subroutine set_site_properties

  ! ============================================================================
  subroutine init_patches( bounds, ed_allsites_inst )
    !
    ! !DESCRIPTION:
    !initialize patches on new ground
    !
    ! !USES:
    use EDParamsMod ,  only : ED_val_maxspread
    !
    ! !ARGUMENTS    
    type(bounds_type)  , intent(in)    :: bounds 
    type(ed_site_type) , intent(inout), target :: ed_allsites_inst( bounds%begg: )
    !
    ! !LOCAL VARIABLES:
    integer  :: g
    real(r8) :: cwd_ag_local(ncwd)
    real(r8) :: cwd_bg_local(ncwd)
    real(r8) :: spread_local(nclmax)
    real(r8) :: leaf_litter_local(numpft_ed)
    real(r8) :: root_litter_local(numpft_ed)
    real(r8) :: seed_bank_local(numpft_ed)
    real(r8) :: age !notional age of this patch
    type(ed_patch_type), pointer :: newp
    !----------------------------------------------------------------------

    cwd_ag_local(:)      = 0.0_r8 !ED_val_init_litter -- arbitrary value for litter pools. kgC m-2
    cwd_bg_local(:)      = 0.0_r8 !ED_val_init_litter
    leaf_litter_local(:) = 0.0_r8
    root_litter_local(:) = 0.0_r8
    spread_local(:)      = ED_val_maxspread
    seed_bank_local(:)   = 0.0_r8 !Note (mv,11-04-2014, this is a bug fix - this line was missing)
    age                  = 0.0_r8

    !FIX(SPM,032414) clean this up...inits out of this loop
    do g = bounds%begg,bounds%endg

       allocate(newp)
!       call zero_patch(newp) !Note (mv,11-04-2014, this is a bug fix - this line was missing)

       newp%patchno = 1
       newp%younger => null()
       newp%older   => null()

       ed_allsites_inst(g)%youngest_patch => newp
       ed_allsites_inst(g)%youngest_patch => newp
       ed_allsites_inst(g)%oldest_patch   => newp

       ! make new patch...
       call create_patch(ed_allsites_inst(g), newp, age, AREA, &
            spread_local, cwd_ag_local, cwd_bg_local, leaf_litter_local,  &
            root_litter_local, seed_bank_local) 

       call init_cohorts(newp)

    enddo !gridcells

  end subroutine init_patches

  ! ============================================================================
  subroutine init_cohorts( patch_in )
    !
    ! !DESCRIPTION:
    ! initialize new cohorts on bare ground
    !
    ! !USES:
    !
    ! !ARGUMENTS    
    type(ed_patch_type), intent(inout), pointer  :: patch_in
    !
    ! !LOCAL VARIABLES:
    type(ed_cohort_type),pointer :: temp_cohort
    integer :: cstatus
    integer :: pft
    !----------------------------------------------------------------------

    patch_in%tallest  => null()
    patch_in%shortest => null()

    do pft =  1,numpft_ed !FIX(RF,032414) - turning off veg dynamics

       allocate(temp_cohort) ! temporary cohort

       temp_cohort%pft         = pft
       temp_cohort%n           = EDecophyscon%initd(pft) * patch_in%area
       temp_cohort%hite        = EDecophyscon%hgt_min(pft)
       temp_cohort%dbh         = Dbh(temp_cohort) ! FIX(RF, 090314) - comment out addition of ' + 0.0001_r8*pft   '  - seperate out PFTs a little bit...
       temp_cohort%canopy_trim = 1.0_r8
       temp_cohort%bdead       = Bdead(temp_cohort)
       temp_cohort%balive      = Bleaf(temp_cohort)*(1.0_r8 + pftcon%froot_leaf(pft) &
            + EDecophyscon%sapwood_ratio(temp_cohort%pft)*temp_cohort%hite)
       temp_cohort%b           = temp_cohort%balive + temp_cohort%bdead

       if( pftcon%evergreen(pft) == 1) then
          temp_cohort%bstore = Bleaf(temp_cohort) * EDecophyscon%cushion(pft)
          temp_cohort%laimemory = 0._r8
          cstatus = 2
       endif

       if( pftcon%season_decid(pft) == 1 ) then !for dorment places
          temp_cohort%bstore = Bleaf(temp_cohort) * EDecophyscon%cushion(pft) !stored carbon in new seedlings.
          if(patch_in%siteptr%status == 2)then 
             temp_cohort%laimemory = 0.0_r8
          else
             temp_cohort%laimemory = Bleaf(temp_cohort)
          endif
          ! reduce biomass according to size of store, this will be recovered when elaves com on.
          temp_cohort%balive = temp_cohort%balive - temp_cohort%laimemory
          cstatus = patch_in%siteptr%status
       endif

       if ( pftcon%stress_decid(pft) == 1 ) then
          temp_cohort%bstore = Bleaf(temp_cohort) * EDecophyscon%cushion(pft)
          temp_cohort%laimemory = Bleaf(temp_cohort)
          temp_cohort%balive = temp_cohort%balive - temp_cohort%laimemory
          cstatus = patch_in%siteptr%dstatus
       endif

       call create_cohort(patch_in, pft, temp_cohort%n, temp_cohort%hite, temp_cohort%dbh, &
            temp_cohort%balive, temp_cohort%bdead, temp_cohort%bstore, &
            temp_cohort%laimemory,  cstatus, temp_cohort%canopy_trim, 1)

       deallocate(temp_cohort) ! get rid of temporary cohort

    enddo !numpft

    call fuse_cohorts(patch_in)
    call sort_cohorts(patch_in)

  end subroutine init_cohorts

end module EDInitMod
