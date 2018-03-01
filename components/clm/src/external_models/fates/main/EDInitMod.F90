module EDInitMod

  ! ============================================================================
  ! Contains all modules to set up the ED structure. 
  ! ============================================================================

  use FatesConstantsMod         , only : r8 => fates_r8
  use FatesConstantsMod         , only : ifalse
  use FatesConstantsMod         , only : itrue
  use FatesGlobals              , only : endrun => fates_endrun
  use EDTypesMod                , only : nclmax
  use FatesGlobals              , only : fates_log
  use FatesInterfaceMod         , only : hlm_is_restart
  use EDPftvarcon               , only : EDPftvarcon_inst
  use EDGrowthFunctionsMod      , only : bdead, bleaf, dbh
  use EDCohortDynamicsMod       , only : create_cohort, fuse_cohorts, sort_cohorts
  use EDPatchDynamicsMod        , only : create_patch
  use EDTypesMod                , only : ed_site_type, ed_patch_type, ed_cohort_type
  use EDTypesMod                , only : ncwd
  use EDTypesMod                , only : nuMWaterMem
  use EDTypesMod                , only : maxpft
  use EDTypesMod                , only : AREA
  use FatesInterfaceMod         , only : bc_in_type
  use FatesInterfaceMod         , only : hlm_use_planthydro
  use FatesInterfaceMod         , only : hlm_use_inventory_init
  use FatesInterfaceMod         , only : numpft
  use ChecksBalancesMod         , only : SiteCarbonStock
  use FatesInterfaceMod         , only : nlevsclass

  ! CIME GLOBALS
  use shr_log_mod               , only : errMsg => shr_log_errMsg

  implicit none
  private

  logical   ::  DEBUG = .false.

  character(len=*), parameter, private :: sourcefile = &
        __FILE__

  public  :: zero_site
  public  :: init_site_vars
  public  :: init_patches
  public  :: set_site_properties
  private :: init_cohorts

  ! ============================================================================

contains

  ! ============================================================================

  subroutine init_site_vars( site_in )
    !
    ! !DESCRIPTION:
    !
    !
    ! !ARGUMENTS    
    type(ed_site_type), intent(inout) ::  site_in
    !
    ! !LOCAL VARIABLES:
    !----------------------------------------------------------------------
    !
    allocate(site_in%terminated_nindivs(1:nlevsclass,1:numpft,2))
    allocate(site_in%demotion_rate(1:nlevsclass))
    allocate(site_in%promotion_rate(1:nlevsclass))
    !
    end subroutine init_site_vars

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
    
    ! DISTURBANCE
    site_in%total_burn_flux_to_atm = 0._r8

    ! PHENOLOGY 
    site_in%status           = 0    ! are leaves in this pixel on or off?
    site_in%dstatus          = 0
    site_in%ED_GDD_site      = nan  ! growing degree days
    site_in%ncd              = nan  ! no chilling days
    site_in%last_n_days(:)   = 999  ! record of last 10 days temperature for senescence model.
    site_in%leafondate       = 999  ! doy of leaf on
    site_in%leafoffdate      = 999  ! doy of leaf off
    site_in%dleafondate      = 999  ! doy of leaf on drought
    site_in%dleafoffdate     = 999  ! doy of leaf on drought
    site_in%water_memory(:)  = nan


    ! SEED
    site_in%seed_bank(:)     = 0._r8

    ! FIRE 
    site_in%acc_ni           = 0.0_r8     ! daily nesterov index accumulating over time. time unlimited theoretically.
    site_in%frac_burnt       = 0.0_r8     ! burn area read in from external file

    ! BGC Balance Checks
    site_in%fates_to_bgc_this_ts = 0.0_r8
    site_in%fates_to_bgc_last_ts = 0.0_r8

    ! termination and recruitment info
    site_in%terminated_nindivs(:,:,:) = 0._r8
    site_in%termination_carbonflux(:) = 0._r8
    site_in%recruitment_rate(:) = 0._r8

    ! demotion/promotion info
    site_in%demotion_rate(:) = 0._r8
    site_in%demotion_carbonflux = 0._r8
    site_in%promotion_rate(:) = 0._r8
    site_in%promotion_carbonflux = 0._r8

    ! diagnostic site-level cwd and litter fluxes
    site_in%CWD_AG_diagnostic_input_carbonflux(:) = 0._r8
    site_in%CWD_BG_diagnostic_input_carbonflux(:) = 0._r8
    site_in%leaf_litter_diagnostic_input_carbonflux(:) = 0._r8
    site_in%root_litter_diagnostic_input_carbonflux(:) = 0._r8
    
    ! Resources management (logging/harvesting, etc)
    site_in%resources_management%trunk_product_site  = 0.0_r8

    ! canopy spread
    site_in%spread = 0._r8

  end subroutine zero_site

  ! ============================================================================
  subroutine set_site_properties( nsites, sites)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    !
    ! !ARGUMENTS    

    integer, intent(in)                        :: nsites
    type(ed_site_type) , intent(inout), target :: sites(nsites)
    !
    ! !LOCAL VARIABLES:
    integer  :: s
    real(r8) :: leafon
    real(r8) :: leafoff
    real(r8) :: stat
    real(r8) :: NCD
    real(r8) :: GDD
    real(r8) :: dstat
    real(r8) :: acc_NI
    real(r8) :: watermem
    integer  :: dleafoff
    integer  :: dleafon
    !----------------------------------------------------------------------

    if ( hlm_is_restart == ifalse ) then
       !initial guess numbers for site condition.
       NCD      = 0.0_r8
       GDD      = 30.0_r8
       leafon   = 100.0_r8
       leafoff  = 300.0_r8
       stat     = 2
       acc_NI   = 0.0_r8
       dstat    = 2
       dleafoff = 300
       dleafon  = 100
       watermem = 0.5_r8

    else ! assignements for restarts

       NCD      = 1.0_r8 ! NCD should be 1 on restart
       GDD      = 0.0_r8
       leafon   = 0.0_r8
       leafoff  = 0.0_r8
       stat     = 1
       acc_NI   = 0.0_r8
       dstat    = 2
       dleafoff = 300
       dleafon  = 100
       watermem = 0.5_r8

    endif

    do s = 1,nsites
       sites(s)%ncd          = NCD
       sites(s)%leafondate   = leafon
       sites(s)%leafoffdate  = leafoff
       sites(s)%dleafoffdate = dleafoff
       sites(s)%dleafondate  = dleafon
       sites(s)%ED_GDD_site  = GDD

       if ( hlm_is_restart == ifalse ) then
          sites(s)%water_memory(1:numWaterMem) = watermem
       end if

       sites(s)%status = stat
       !start off with leaves off to initialise
       sites(s)%dstatus= dstat
       
       sites(s)%acc_NI     = acc_NI
       sites(s)%frac_burnt = 0.0_r8
       sites(s)%old_stock  = 0.0_r8

       sites(s)%spread     = 1.0_r8
    end do

    return
  end subroutine set_site_properties

  ! ============================================================================
  subroutine init_patches( nsites, sites, bc_in)
     !
     ! !DESCRIPTION:
     ! initialize patches
     ! This may be call a near bare ground initialization, or it may
     ! load patches from an inventory.

     !
     

     use FatesPlantHydraulicsMod, only : updateSizeDepRhizHydProps 
     use FatesInventoryInitMod,   only : initialize_sites_by_inventory

     !
     ! !ARGUMENTS    
     integer, intent(in)                        :: nsites
     type(ed_site_type) , intent(inout), target :: sites(nsites)
     type(bc_in_type), intent(in)               :: bc_in(nsites)
     !
     ! !LOCAL VARIABLES:
     integer  :: s
     real(r8) :: cwd_ag_local(ncwd)
     real(r8) :: cwd_bg_local(ncwd)
     real(r8) :: leaf_litter_local(maxpft)
     real(r8) :: root_litter_local(maxpft)
     real(r8) :: age !notional age of this patch

     ! dummy locals
     real(r8) :: biomass_stock
     real(r8) :: litter_stock
     real(r8) :: seed_stock

     type(ed_patch_type), pointer :: newp

     ! List out some nominal patch values that are used for Near Bear Ground initializations
     ! as well as initializing inventory
     ! ---------------------------------------------------------------------------------------------
     cwd_ag_local(:)      = 0.0_r8 !ED_val_init_litter -- arbitrary value for litter pools. kgC m-2
     cwd_bg_local(:)      = 0.0_r8 !ED_val_init_litter
     leaf_litter_local(:) = 0.0_r8
     root_litter_local(:) = 0.0_r8
     age                  = 0.0_r8
     ! ---------------------------------------------------------------------------------------------

     ! ---------------------------------------------------------------------------------------------
     ! Two primary options, either a Near Bear Ground (NBG) or Inventory based cold-start
     ! ---------------------------------------------------------------------------------------------

     if ( hlm_use_inventory_init.eq.itrue ) then

        call initialize_sites_by_inventory(nsites,sites,bc_in)

        do s = 1, nsites
           if (hlm_use_planthydro.eq.itrue) then
              call updateSizeDepRhizHydProps(sites(s), bc_in(s))
           end if
           ! For carbon balance checks, we need to initialize the 
           ! total carbon stock
           call SiteCarbonStock(sites(s),sites(s)%old_stock,biomass_stock,litter_stock,seed_stock)
           
        enddo
        
     else

        !FIX(SPM,032414) clean this up...inits out of this loop
        do s = 1, nsites

           allocate(newp)

           newp%patchno = 1
           newp%younger => null()
           newp%older   => null()

           sites(s)%youngest_patch => newp
           sites(s)%youngest_patch => newp
           sites(s)%oldest_patch   => newp

           ! make new patch...
           call create_patch(sites(s), newp, age, AREA, &
                 cwd_ag_local, cwd_bg_local, leaf_litter_local,  &
                 root_litter_local) 

           call init_cohorts(newp, bc_in(s))

           ! This sets the rhizosphere shells based on the plant initialization
           ! The initialization of the plant-relevant hydraulics variables
           ! were set from a call inside of the init_cohorts()->create_cohort() subroutine
           if (hlm_use_planthydro.eq.itrue) then
              call updateSizeDepRhizHydProps(sites(s), bc_in(s))
           end if


           ! For carbon balance checks, we need to initialize the 
           ! total carbon stock
           call SiteCarbonStock(sites(s),sites(s)%old_stock,biomass_stock,litter_stock,seed_stock)

        enddo

     end if

  end subroutine init_patches

  ! ============================================================================
  subroutine init_cohorts( patch_in, bc_in)
    !
    ! !DESCRIPTION:
    ! initialize new cohorts on bare ground
    !
    ! !USES:
    !
    ! !ARGUMENTS    
    type(ed_patch_type), intent(inout), pointer  :: patch_in
    type(bc_in_type), intent(in)                 :: bc_in
    !
    ! !LOCAL VARIABLES:
    type(ed_cohort_type),pointer :: temp_cohort
    integer :: cstatus
    integer :: pft
    !----------------------------------------------------------------------

    patch_in%tallest  => null()
    patch_in%shortest => null()

    do pft =  1,numpft

       if(EDPftvarcon_inst%initd(pft)>1.0E-7) then

       allocate(temp_cohort) ! temporary cohort

       temp_cohort%pft         = pft
       temp_cohort%n           = EDPftvarcon_inst%initd(pft) * patch_in%area
       temp_cohort%hite        = EDPftvarcon_inst%hgt_min(pft)
       !temp_cohort%n           = 0.5_r8 * 0.0028_r8 * patch_in%area  ! BOC for fixed size runs EDPftvarcon_inst%initd(pft) * patch_in%area
       !temp_cohort%hite        = 28.65_r8                            ! BOC translates to DBH of 50cm. EDPftvarcon_inst%hgt_min(pft)
       temp_cohort%dbh         = Dbh(temp_cohort) ! FIX(RF, 090314) - comment out addition of ' + 0.0001_r8*pft   '  - seperate out PFTs a little bit...
       temp_cohort%canopy_trim = 1.0_r8
       temp_cohort%bdead       = Bdead(temp_cohort)
       temp_cohort%balive      = Bleaf(temp_cohort)*(1.0_r8 + EDPftvarcon_inst%allom_l2fr(pft) &
            + EDPftvarcon_inst%allom_latosa_int(temp_cohort%pft)*temp_cohort%hite)
       temp_cohort%b           = temp_cohort%balive + temp_cohort%bdead

       if( EDPftvarcon_inst%evergreen(pft) == 1) then
          temp_cohort%bstore = Bleaf(temp_cohort) * EDPftvarcon_inst%cushion(pft)
          temp_cohort%laimemory = 0._r8
          cstatus = 2
       endif

       if( EDPftvarcon_inst%season_decid(pft) == 1 ) then !for dorment places
          temp_cohort%bstore = Bleaf(temp_cohort) * EDPftvarcon_inst%cushion(pft) !stored carbon in new seedlings.
          if(patch_in%siteptr%status == 2)then 
             temp_cohort%laimemory = 0.0_r8
          else
             temp_cohort%laimemory = Bleaf(temp_cohort)
          endif
          ! reduce biomass according to size of store, this will be recovered when elaves com on.
          temp_cohort%balive = temp_cohort%balive - temp_cohort%laimemory
          cstatus = patch_in%siteptr%status
       endif

       if ( EDPftvarcon_inst%stress_decid(pft) == 1 ) then
          temp_cohort%bstore = Bleaf(temp_cohort) * EDPftvarcon_inst%cushion(pft)
          temp_cohort%laimemory = Bleaf(temp_cohort)
          temp_cohort%balive = temp_cohort%balive - temp_cohort%laimemory
          cstatus = patch_in%siteptr%dstatus
       endif

       if ( DEBUG ) write(fates_log(),*) 'EDInitMod.F90 call create_cohort '

       call create_cohort(patch_in, pft, temp_cohort%n, temp_cohort%hite, temp_cohort%dbh, &
            temp_cohort%balive, temp_cohort%bdead, temp_cohort%bstore, &
            temp_cohort%laimemory,  cstatus, temp_cohort%canopy_trim, 1, bc_in)

       deallocate(temp_cohort) ! get rid of temporary cohort

       endif

    enddo !numpft

    call fuse_cohorts(patch_in,bc_in)
    call sort_cohorts(patch_in)

  end subroutine init_cohorts

  ! ===============================================================================================


end module EDInitMod
