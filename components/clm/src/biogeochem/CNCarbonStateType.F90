module CNCarbonStateType

#include "shr_assert.h"

  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_infnan_mod         , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use clm_varpar             , only : ndecomp_cascade_transitions, ndecomp_pools, nlevcan
  use clm_varpar             , only : nlevdecomp_full, crop_prog, nlevdecomp
  use clm_varcon             , only : spval, ispval, dzsoi_decomp, zisoi, zsoi
  use landunit_varcon        , only : istcrop 
  use clm_varctl             , only : iulog, use_vertsoilc, use_cndv, spinup_state 
  use decompMod              , only : bounds_type
  use CNStateType            , only : cnstate_type
  use pftvarcon              , only : npcropmin
  use CNDecompCascadeConType , only : decomp_cascade_con
  use VegetationPropertiesType         , only : veg_vp
  use abortutils             , only : endrun
  use spmdMod                , only : masterproc 
  use subgridAveMod          , only : p2c
  use LandunitType           , only : lun_pp                
  use ColumnType             , only : col_pp                
  use clm_varctl             , only : nu_com, use_fates, use_crop
  use VegetationType         , only : veg_pp
  use CNSpeciesMod           , only : species_from_string, species_name_from_string
  use dynPatchStateUpdaterMod, only : patch_state_updater_type
  use NutrientStateType      , only : nutrientstate_type
  use NutrientStateType      , only : NutrientStateInitHistory, NutrientStateInitAllocate
  use NutrientStateType      , only : NutrientStateDynamicPatchAdjustments
  use NutrientStateType      , only : NutrientStateRestart
  use NutrientStateType      , only : NutrientStatePatchSummary
  use NutrientStateType      , only : NutrientStateColumnSummary

  ! bgc interface & pflotran
  use clm_varctl             , only : use_clm_interface, use_pflotran, pf_cmode
  
  ! 
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  type, public, extends(nutrientstate_type) :: carbonstate_type
     
     real(r8), pointer :: gresp_storage_patch      (:)     ! (gC/m2) growth respiration storage
     real(r8), pointer :: gresp_xfer_patch         (:)     ! (gC/m2) growth respiration transfer
     real(r8), pointer :: xsmrpool_patch           (:)     ! (gC/m2) abstract C pool to meet excess MR demand
     real(r8), pointer :: woodc_patch              (:)     ! (gC/m2) wood C
     real(r8), pointer :: leafcmax_patch           (:)     ! (gC/m2) ann max leaf C
     real(r8), pointer :: rootc_col                (:)     ! (gC/m2) root carbon at column level (fire)
     real(r8), pointer :: leafc_col                (:)     ! (gC/m2) column-level leafc (fire)
     real(r8), pointer :: deadstemc_col            (:)     ! (gC/m2) column-level deadstemc (fire)
     real(r8), pointer :: fuelc_col                (:)     ! fuel avalability factor for Reg.C (0-1)
     real(r8), pointer :: fuelc_crop_col           (:)     ! fuel avalability factor for Reg.A (0-1)
     real(r8), pointer :: decomp_som2c_vr_col(:,:)

   contains

     procedure , public  :: Init   
     procedure , public  :: SetValues 
     procedure , public  :: ZeroDWT
     procedure , public  :: Restart
     procedure , public  :: Summary
     procedure , public  :: DynamicPatchAdjustments
     procedure , public  :: DynamicColumnAdjustments
     procedure , private :: InitAllocate 
     procedure , private :: InitHistory  
     procedure , private :: InitCold     

  end type carbonstate_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, carbon_type, ratio, c12_carbonstate_vars)

    class(carbonstate_type)                       :: this
    type(bounds_type)      , intent(in)           :: bounds  
    character(len=3)       , intent(in)           :: carbon_type
    real(r8)               , intent(in)           :: ratio
    type(carbonstate_type) , intent(in), optional :: c12_carbonstate_vars

    this%species      = species_from_string(carbon_type)
    this%name         = species_name_from_string(carbon_type)
    this%restart_name = 'c'

    call this%InitAllocate ( bounds)
    call this%InitHistory  ( bounds, carbon_type)
    call this%InitCold     ( bounds, ratio, c12_carbonstate_vars)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !ARGUMENTS:
    class (carbonstate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begp,endp
    integer           :: begc,endc
    integer           :: begg,endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg

    call NutrientStateInitAllocate(this, bounds)

    if ( .not. use_fates ) then
       allocate(this%gresp_storage_patch (begp :endp)); this%gresp_storage_patch (:)   = nan
       allocate(this%gresp_xfer_patch    (begp :endp)); this%gresp_xfer_patch    (:)   = nan
       allocate(this%xsmrpool_patch      (begp :endp)); this%xsmrpool_patch      (:)   = nan
       allocate(this%leafcmax_patch      (begp :endp)); this%leafcmax_patch      (:)   = nan
       allocate(this%woodc_patch         (begp :endp)); this%woodc_patch         (:)   = nan
    endif
    allocate(this%rootc_col      (begc :endc)); this%rootc_col      (:)   = nan
    allocate(this%leafc_col      (begc :endc)); this%leafc_col      (:)   = nan
    allocate(this%deadstemc_col  (begc :endc)); this%deadstemc_col  (:)   = nan
    allocate(this%fuelc_col      (begc :endc)); this%fuelc_col      (:)   = nan
    allocate(this%fuelc_crop_col (begc :endc)); this%fuelc_crop_col (:)   = nan

    allocate(this%decomp_som2c_vr_col(begc:endc,1:nlevdecomp_full)); this%decomp_som2c_vr_col(:,:)= nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds, carbon_type)
    !
    ! !DESCRIPTION:
    ! add history fields for all CN variables, always set as default='inactive'
    !
    ! !USES:
    use clm_varpar , only : ndecomp_cascade_transitions, ndecomp_pools
    use clm_varpar , only : nlevdecomp, nlevdecomp_full, nlevgrnd
    use clm_varctl , only : use_c13, use_c14
    use histFileMod, only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp 
    !
    ! !ARGUMENTS:
    class (carbonstate_type) :: this
    type(bounds_type)         , intent(in) :: bounds 
    character(len=3)          , intent(in) :: carbon_type ! one of ['c12', c13','c14']
    !
    ! !LOCAL VARIABLES:
    integer           :: k,l,ii,jj 
    character(10)     :: active
    integer           :: begp,endp
    integer           :: begc,endc
    integer           :: begg,endg 
    character(24)     :: fieldname
    character(100)    :: longname
    real(r8), pointer :: data1dptr(:)   ! temp. pointer for slicing larger arrays
    real(r8), pointer :: data2dptr(:,:) ! temp. pointer for slicing larger arrays
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg

    call NutrientStateInitHistory(this, bounds)

       this%woodc_patch(begp:endp) = spval
       call hist_addfld1d (fname=trim(this%name)//'_WOOD', units='g' // trim(this%name) // '/m^2', &
             avgflag='A', long_name=trim(this%name)//' wood ', &
             ptr_patch=this%woodc_patch)

       this%gresp_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname=trim(this%name)//'_GRESP_STORAGE', units='g' // trim(this%name) // '/m^2', &
             avgflag='A', long_name=trim(this%name)//' growth respiration storage', &
             ptr_patch=this%gresp_storage_patch, default='inactive')

       this%gresp_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname=trim(this%name)//'_GRESP_XFER', units='g' // trim(this%name) // '/m^2', &
             avgflag='A', long_name=trim(this%name)//' growth respiration transfer', &
             ptr_patch=this%gresp_xfer_patch, default='inactive')

       this%xsmrpool_patch(begp:endp) = spval
       call hist_addfld1d (fname=trim(this%name)//'_XSMRPOOL', units='g' // trim(this%name) // '/m^2', &
             avgflag='A', long_name=trim(this%name)//' temporary photosynthate pool', &
             ptr_patch=this%xsmrpool_patch, default='active')

       this%fuelc_col(begc:endc) = spval
       call hist_addfld1d (fname=trim(this%name)//'_FUEL', units='g' // trim(this%name) // '/m^2', &
             avgflag='A', long_name=trim(this%name)//' fuel load', &
             ptr_col=this%fuelc_col, default='inactive')


 end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, ratio, c12_carbonstate_vars)
    !
    ! !DESCRIPTION:
    ! Initializes time varying variables used only in coupled carbon-nitrogen mode (CN):
    !
    ! !USES:
    use landunit_varcon , only: istsoil
    use pftvarcon       , only: noveg, npcropmin
    !
    ! !ARGUMENTS:
    class(carbonstate_type) :: this 
    type(bounds_type), intent(in) :: bounds  
    real(r8), intent(in) :: ratio
    type(carbonstate_type), optional, intent(in) :: c12_carbonstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: p,c,l,g,j,k
    integer :: fc                                        ! filter index
    integer :: num_special_col                           ! number of good values in special_col filter
    integer :: num_special_patch                         ! number of good values in special_patch filter
    integer :: special_col(bounds%endc-bounds%begc+1)    ! special landunit filter - columns
    integer :: special_patch(bounds%endp-bounds%begp+1)  ! special landunit filter - patches
    !-----------------------------------------------------------------------

    ! Set column filters

    num_special_col = 0
    do c = bounds%begc, bounds%endc
       l = col_pp%landunit(c)
       if (lun_pp%ifspecial(l)) then
          num_special_col = num_special_col + 1
          special_col(num_special_col) = c
       end if
    end do

    ! Set patch filters

    num_special_patch = 0
    do p = bounds%begp,bounds%endp
       l = veg_pp%landunit(p)
       if (lun_pp%ifspecial(l)) then
          num_special_patch = num_special_patch + 1
          special_patch(num_special_patch) = p
       end if
    end do


    if ( .not. use_fates ) then
       !-----------------------------------------------
       ! initialize patch-level carbon state variables
       !-----------------------------------------------

       do p = bounds%begp,bounds%endp

          this%leafcmax_patch(p) = 0._r8

          l = veg_pp%landunit(p)
          if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then

             if (veg_pp%itype(p) == noveg) then
                this%leaf_patch(p)         = 0._r8
                this%leaf_storage_patch(p) = 0._r8
             else
                if (veg_vp%evergreen(veg_pp%itype(p)) == 1._r8) then
                   this%leaf_patch(p)         = 1._r8 * ratio
                   this%leaf_storage_patch(p) = 0._r8
                else if (veg_pp%itype(p) >= npcropmin) then ! prognostic crop types
                   this%leaf_patch(p) = 0._r8
                   this%leaf_storage_patch(p) = 0._r8
                else
                   this%leaf_patch(p) = 0._r8
                   this%leaf_storage_patch(p) = 1._r8 * ratio
                end if
             end if
             this%leaf_xfer_patch(p) = 0._r8

             this%froot_patch(p)            = 0._r8 
             this%froot_storage_patch(p)    = 0._r8 
             this%froot_xfer_patch(p)       = 0._r8 

             this%livestem_patch(p)         = 0._r8 
             this%livestem_storage_patch(p) = 0._r8 
             this%livestem_xfer_patch(p)    = 0._r8 

             if (veg_vp%woody(veg_pp%itype(p)) == 1._r8) then
                this%deadstem_patch(p) = 0.1_r8 * ratio
             else
                this%deadstem_patch(p) = 0._r8 
             end if
             this%deadstem_storage_patch(p)  = 0._r8 
             this%deadstem_xfer_patch(p)     = 0._r8

             if (nu_com .ne. 'RD') then
                ! ECA competition calculate root NP uptake as a function of fine root biomass
                ! better to initialize root CNP pools with a non-zero value
                if (veg_pp%itype(p) .ne. noveg) then
                   if (veg_vp%evergreen(veg_pp%itype(p)) == 1._r8) then
                      this%leaf_patch(p) = 20._r8 * ratio
                      this%leaf_storage_patch(p) = 0._r8
                      this%froot_patch(p) = 20._r8 * ratio
                      this%froot_storage_patch(p) = 0._r8
                   else
                      this%leaf_patch(p) = 0._r8 
                      this%leaf_storage_patch(p) = 20._r8 * ratio
                      this%froot_patch(p) = 0._r8
                      this%froot_storage_patch(p) = 20._r8 * ratio
                   end if
                end if
             end if

             this%livecroot_patch(p)         = 0._r8 
             this%livecroot_storage_patch(p) = 0._r8 
             this%livecroot_xfer_patch(p)    = 0._r8 

             this%deadcroot_patch(p)         = 0._r8 
             this%deadcroot_storage_patch(p) = 0._r8 
             this%deadcroot_xfer_patch(p)    = 0._r8 

             this%gresp_storage_patch(p)      = 0._r8 
             this%gresp_xfer_patch(p)         = 0._r8 

             this%pool_patch(p)              = 0._r8 
             this%xsmrpool_patch(p)           = 0._r8 
             this%veg_trunc_patch(p)             = 0._r8 
             this%dispveg_patch(p)           = 0._r8 
             this%storveg_patch(p)           = 0._r8 
             this%totpft_patch(p)            = 0._r8 
             this%woodc_patch(p)              = 0._r8

             if ( crop_prog )then
                this%grain_patch(p)            = 0._r8 
                this%grain_storage_patch(p)    = 0._r8 
                this%grain_xfer_patch(p)       = 0._r8 
                this%cropseed_deficit_patch(p) = 0._r8
             end if

             ! calculate totvegc explicitly so that it is available for the isotope 
             ! code on the first time step.

             this%totveg_patch(p) = &
                  this%leaf_patch(p)              + &
                  this%leaf_storage_patch(p)      + &
                  this%leaf_xfer_patch(p)         + &
                  this%froot_patch(p)             + &
                  this%froot_storage_patch(p)     + &
                  this%froot_xfer_patch(p)        + &
                  this%livestem_patch(p)          + &
                  this%livestem_storage_patch(p)  + &
                  this%livestem_xfer_patch(p)     + &
                  this%deadstem_patch(p)          + &
                  this%deadstem_storage_patch(p)  + &
                  this%deadstem_xfer_patch(p)     + &
                  this%livecroot_patch(p)         + &
                  this%livecroot_storage_patch(p) + &
                  this%livecroot_xfer_patch(p)    + &
                  this%deadcroot_patch(p)         + &
                  this%deadcroot_storage_patch(p) + &
                  this%deadcroot_xfer_patch(p)    + &
                  this%gresp_storage_patch(p)      + &
                  this%gresp_xfer_patch(p)         + &
                  this%pool_patch(p)

             if ( crop_prog )then
                this%totveg_patch(p) =  this%totveg_patch(p) + &
                     this%grain_patch(p)                            + &
                     this%grain_storage_patch(p)                    + &
                     this%grain_xfer_patch(p)
             end if
          endif

       end do
    endif ! .not. use_fates
    
    ! initialize column-level variables
    do c = bounds%begc, bounds%endc
       l = col_pp%landunit(c)
       if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then

          if (.not. present(c12_carbonstate_vars)) then !c12

             do j = 1, nlevdecomp
                do k = 1, ndecomp_pools
                   if (zsoi(j) < 0.3 ) then  !! only initialize upper soil column
                      this%decomp_pools_vr_col(c,j,k) = decomp_cascade_con%initial_stock(k)
                   else
                      this%decomp_pools_vr_col(c,j,k) = 0._r8
                   endif
                end do
                this%soil_trunc_vr_col(c,j) = 0._r8
             end do
             if ( nlevdecomp > 1 ) then
                do j = nlevdecomp+1, nlevdecomp_full
                   do k = 1, ndecomp_pools
                      this%decomp_pools_vr_col(c,j,k) = 0._r8
                   end do
                   this%soil_trunc_vr_col(c,j) = 0._r8
                end do
             end if
             this%decomp_pools_col(c,1:ndecomp_pools)    = decomp_cascade_con%initial_stock(1:ndecomp_pools)
             this%decomp_pools_1m_col(c,1:ndecomp_pools) = decomp_cascade_con%initial_stock(1:ndecomp_pools)

          else

             do j = 1, nlevdecomp
                do k = 1, ndecomp_pools
                   this%decomp_pools_vr_col(c,j,k) = c12_carbonstate_vars%decomp_pools_vr_col(c,j,k) * ratio
                end do
                this%soil_trunc_vr_col(c,j) = c12_carbonstate_vars%soil_trunc_vr_col(c,j) * ratio
             end do
             if ( nlevdecomp > 1 ) then
                do j = nlevdecomp+1, nlevdecomp_full
                   do k = 1, ndecomp_pools
                      this%decomp_pools_vr_col(c,j,k) = 0._r8
                   end do
                   this%soil_trunc_vr_col(c,j) = 0._r8
                end do
             end if
             this%cwd_col(c) = c12_carbonstate_vars%cwd_col(c) * ratio
             do k = 1, ndecomp_pools
                this%decomp_pools_col(c,k)    = c12_carbonstate_vars%decomp_pools_col(c,k) * ratio
                this%decomp_pools_1m_col(c,k) = c12_carbonstate_vars%decomp_pools_1m_col(c,k) * ratio
             end do

          endif

          this%cwd_col(c)       = 0._r8
          this%veg_trunc_col(c)     = 0._r8
          this%totlit_col(c)    = 0._r8
          this%totsom_col(c)    = 0._r8
          this%totlit_1m_col(c) = 0._r8
          this%totsom_1m_col(c) = 0._r8
          this%totecosys_col(c) = 0._r8
          this%totcol_col(c)    = 0._r8

          ! dynamic landcover state variables
          this%seed_col(c)      = 0._r8
          this%prod10_col(c)    = 0._r8
          this%prod100_col(c)   = 0._r8
          this%prod1_col(c)     = 0._r8
          this%totprod_col(c)   = 0._r8

       end if

    end do

    ! now loop through special filters and explicitly set the variables that
    ! have to be in place for biogeophysics
    
    do fc = 1,num_special_col
       c = special_col(fc)

       this%seed_col(c)      = 0._r8
       this%prod10_col(c)    = 0._r8
       this%prod100_col(c)   = 0._r8
       this%prod1_col(c)     = 0._r8
       this%totprod_col(c)   = 0._r8
    end do

    do g = bounds%begg, bounds%endg
       this%seed_grc(g) = 0._r8
    end do

    ! initialize fields for special filters

    call this%SetValues (&
         num_patch=num_special_patch, filter_patch=special_patch, value_patch=0._r8, &
         num_column=num_special_col, filter_column=special_col, value_column=0._r8)

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine Restart ( this,  bounds, ncid, flag, carbon_type, c12_carbonstate_vars, cnstate_vars)
    !
    ! !DESCRIPTION: 
    ! Read/write CN restart data for carbon state
    !
    ! !USES:
    use shr_infnan_mod   , only : isnan => shr_infnan_isnan, nan => shr_infnan_nan, assignment(=)
    use shr_const_mod    , only : SHR_CONST_PDB
    use clm_time_manager , only : is_restart, get_nstep
    use clm_varcon       , only : c13ratio, c14ratio
    use clm_varctl       , only : spinup_mortality_factor, spinup_state

    use restUtilMod
    use ncdio_pio
    use tracer_varcon  , only : is_active_betr_bgc
    !
    ! !ARGUMENTS:
    class (carbonstate_type) :: this
    type(bounds_type)         , intent(in)           :: bounds 
    type(file_desc_t)         , intent(inout)        :: ncid   ! netcdf id
    character(len=*)          , intent(in)           :: flag   !'read' or 'write'
    character(len=3)          , intent(in)           :: carbon_type ! 'c12' or 'c13' or 'c14'
    type (carbonstate_type)   , intent(in), optional :: c12_carbonstate_vars 
    type (cnstate_type)       , intent(in)           :: cnstate_vars

    !
    ! !LOCAL VARIABLES:
    integer  :: i,j,k,l,c
    real(r8) :: vwc,psi             ! for calculating soilpsi
    real(r8) :: c3_del13c           ! typical del13C for C3 photosynthesis (permil, relative to PDB)
    real(r8) :: c4_del13c           ! typical del13C for C4 photosynthesis (permil, relative to PDB)
    real(r8) :: c3_r1               ! isotope ratio (13c/12c) for C3 photosynthesis
    real(r8) :: c4_r1               ! isotope ratio (13c/12c) for C4 photosynthesis
    real(r8) :: c3_r2               ! isotope ratio (13c/[12c+13c]) for C3 photosynthesis
    real(r8) :: c4_r2               ! isotope ratio (13c/[12c+13c]) for C4 photosynthesis
    real(r8) :: m, m_veg            ! multiplier for the exit_spinup code
    real(r8), pointer :: ptr2d(:,:) ! temp. pointers for slicing larger arrays
    real(r8), pointer :: ptr1d(:)   ! temp. pointers for slicing larger arrays
    character(len=128) :: varname   ! temporary
    logical  :: readvar
    integer  :: idata
    logical  :: exit_spinup  = .false.
    logical  :: enter_spinup = .false.
    ! spinup state as read from restart file, for determining whether to enter or exit spinup mode.
    integer  :: restart_file_spinup_state 
    ! flags for comparing the model and restart decomposition cascades
    integer  :: decomp_cascade_state, restart_file_decomp_cascade_state 
    !------------------------------------------------------------------------

    if (carbon_type == 'c13' .or. carbon_type == 'c14') then
       if (.not. present(c12_carbonstate_vars)) then
          call endrun(msg=' ERROR: for C14 must pass in c12_carbonstate_vars as argument' //&
               errMsg(__FILE__, __LINE__))
       end if
    end if

    if (carbon_type == 'c12') call NutrientStateRestart(this, bounds, ncid, flag)

    if ( .not. use_fates ) then

       if (carbon_type == 'c12') then
          call restartvar(ncid=ncid, flag=flag, varname='gresp_storage', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%gresp_storage_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='gresp_xfer', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%gresp_xfer_patch) 

          call restartvar(ncid=ncid, flag=flag, varname='xsmrpool', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%xsmrpool_patch) 
       end if

       !--------------------------------
       ! C13 pft carbon state variables 
       !--------------------------------

       if ( carbon_type == 'c13')  then

          if ( .not. is_restart() .and. get_nstep() == 1 ) then
             c3_del13c = -28._r8
             c4_del13c = -13._r8
             c3_r1 = SHR_CONST_PDB + ((c3_del13c*SHR_CONST_PDB)/1000._r8)
             c3_r2 = c3_r1/(1._r8 + c3_r1)
             c4_r1 = SHR_CONST_PDB + ((c4_del13c*SHR_CONST_PDB)/1000._r8)
             c4_r2 = c4_r1/(1._r8 + c4_r1)

             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%grain_patch(i)            = c12_carbonstate_vars%grain_patch(i)         * c3_r2
                   this%grain_storage_patch(i)    = c12_carbonstate_vars%grain_storage_patch(i) * c3_r2
                   this%grain_xfer_patch(i)       = c12_carbonstate_vars%grain_xfer_patch(i)    * c3_r2
                   this%dispveg_patch(i)          = c12_carbonstate_vars%dispveg_patch(i)       * c3_r2
                   this%storveg_patch(i)          = c12_carbonstate_vars%storveg_patch(i)       * c3_r2
                   this%totveg_patch(i)           = c12_carbonstate_vars%totveg_patch(i)        * c3_r2
                   this%totpft_patch(i)           = c12_carbonstate_vars%totpft_patch(i)        * c3_r2
                   this%woodc_patch(i)             = c12_carbonstate_vars%woodc_patch(i)          * c3_r2
                else
                   this%grain_patch(i)            = c12_carbonstate_vars%grain_patch(i)         * c4_r2
                   this%grain_storage_patch(i)    = c12_carbonstate_vars%grain_storage_patch(i) * c4_r2
                   this%grain_xfer_patch(i)       = c12_carbonstate_vars%grain_xfer_patch(i)    * c4_r2
                   this%dispveg_patch(i)          = c12_carbonstate_vars%dispveg_patch(i)       * c4_r2
                   this%storveg_patch(i)          = c12_carbonstate_vars%storveg_patch(i)       * c4_r2
                   this%totveg_patch(i)           = c12_carbonstate_vars%totveg_patch(i)        * c4_r2
                   this%totpft_patch(i)           = c12_carbonstate_vars%totpft_patch(i)        * c4_r2
                   this%woodc_patch(i)             = c12_carbonstate_vars%woodc_patch(i)          * c4_r2
                end if
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='leafc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leaf_patch)
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%leafc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%leaf_patch(i) = c12_carbonstate_vars%leaf_patch(i) * c3_r2
                else
                   this%leaf_patch(i) = c12_carbonstate_vars%leaf_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='leafc_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leaf_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%leafc_storage with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%leaf_storage_patch(i) = c12_carbonstate_vars%leaf_storage_patch(i) * c3_r2
                else
                   this%leaf_storage_patch(i) = c12_carbonstate_vars%leaf_storage_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='leafc_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leaf_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%leafc_xfer with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%leaf_xfer_patch(i) = c12_carbonstate_vars%leaf_xfer_patch(i) * c3_r2
                else
                   this%leaf_xfer_patch(i) = c12_carbonstate_vars%leaf_xfer_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='frootc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%froot_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%frootc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%froot_patch(i) = c12_carbonstate_vars%froot_patch(i) * c3_r2
                else
                   this%froot_patch(i) = c12_carbonstate_vars%froot_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='frootc_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%froot_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%frootc_storage with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%froot_storage_patch(i) = c12_carbonstate_vars%froot_storage_patch(i) * c3_r2
                else
                   this%froot_storage_patch(i) = c12_carbonstate_vars%froot_storage_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='frootc_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%froot_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%frootc_xfer with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%froot_xfer_patch(i) = c12_carbonstate_vars%froot_xfer_patch(i) * c3_r2
                else
                   this%froot_xfer_patch(i) = c12_carbonstate_vars%froot_xfer_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='livestemc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestem_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livestemc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%livestem_patch(i) = c12_carbonstate_vars%livestem_patch(i) * c3_r2
                else
                   this%livestem_patch(i) = c12_carbonstate_vars%livestem_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='livestemc_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestem_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livestemc_storage with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%livestem_storage_patch(i) = c12_carbonstate_vars%livestem_storage_patch(i) * c3_r2
                else
                   this%livestem_storage_patch(i) = c12_carbonstate_vars%livestem_storage_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='livestemc_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestem_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livestemc_xfer with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%livestem_xfer_patch(i) = c12_carbonstate_vars%livestem_xfer_patch(i) * c3_r2
                else
                   this%livestem_xfer_patch(i) = c12_carbonstate_vars%livestem_xfer_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstem_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadstemc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%deadstem_patch(i) = c12_carbonstate_vars%deadstem_patch(i) * c3_r2
                else
                   this%deadstem_patch(i) = c12_carbonstate_vars%deadstem_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstem_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadstemc_storage with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%deadstem_storage_patch(i) = c12_carbonstate_vars%deadstem_storage_patch(i) * c3_r2
                else
                   this%deadstem_storage_patch(i) = c12_carbonstate_vars%deadstem_storage_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstem_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadstemc_xfer with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%deadstem_xfer_patch(i) = c12_carbonstate_vars%deadstem_xfer_patch(i) * c3_r2
                else
                   this%deadstem_xfer_patch(i) = c12_carbonstate_vars%deadstem_xfer_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecroot_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livecrootc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%livecroot_patch(i) = c12_carbonstate_vars%livecroot_patch(i) * c3_r2
                else
                   this%livecroot_patch(i) = c12_carbonstate_vars%livecroot_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecroot_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livecrootc_storage with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%livecroot_storage_patch(i) = c12_carbonstate_vars%livecroot_storage_patch(i) * c3_r2
                else
                   this%livecroot_storage_patch(i) = c12_carbonstate_vars%livecroot_storage_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecroot_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livecrootc_xfer with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%livecroot_xfer_patch(i) = c12_carbonstate_vars%livecroot_xfer_patch(i) * c3_r2
                else
                   this%livecroot_xfer_patch(i) = c12_carbonstate_vars%livecroot_xfer_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcroot_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadcrootc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%deadcroot_patch(i) = c12_carbonstate_vars%deadcroot_patch(i) * c3_r2
                else
                   this%deadcroot_patch(i) = c12_carbonstate_vars%deadcroot_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcroot_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadcrootc_storage with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%deadcroot_storage_patch(i) = c12_carbonstate_vars%deadcroot_storage_patch(i) * c3_r2
                else
                   this%deadcroot_storage_patch(i) = c12_carbonstate_vars%deadcroot_storage_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_xfer_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcroot_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadcrootc_xfer with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%deadcroot_xfer_patch(i) = c12_carbonstate_vars%deadcroot_xfer_patch(i) * c3_r2
                else
                   this%deadcroot_xfer_patch(i) = c12_carbonstate_vars%deadcroot_xfer_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='gresp_storage_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%gresp_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%gresp_storage with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%gresp_storage_patch(i) = c12_carbonstate_vars%gresp_storage_patch(i) * c3_r2
                else
                   this%gresp_storage_patch(i) = c12_carbonstate_vars%gresp_storage_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='gresp_xfer_13', xtype=ncd_double,  &
               dim1name='pft', &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%gresp_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%gresp_xfer with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%gresp_xfer_patch(i) = c12_carbonstate_vars%gresp_xfer_patch(i) * c3_r2
                else
                   this%gresp_xfer_patch(i) = c12_carbonstate_vars%gresp_xfer_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='cpool_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%pool_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%cpool with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%pool_patch(i) = c12_carbonstate_vars%pool_patch(i) * c3_r2
                else
                   this%pool_patch(i) = c12_carbonstate_vars%pool_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='xsmrpool_13', xtype=ncd_double,  &
               dim1name='pft', &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%xsmrpool_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%xsmrpool with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%xsmrpool_patch(i) = c12_carbonstate_vars%xsmrpool_patch(i) * c3_r2
                else
                   this%xsmrpool_patch(i) = c12_carbonstate_vars%xsmrpool_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='pft_ctrunc_13', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%veg_trunc_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%ctrunc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%veg_trunc_patch(i) = c12_carbonstate_vars%veg_trunc_patch(i) * c3_r2
                else
                   this%veg_trunc_patch(i) = c12_carbonstate_vars%veg_trunc_patch(i) * c4_r2
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c13')  then
          call restartvar(ncid=ncid, flag=flag, varname='totvegc_13', xtype=ncd_double,  &
               dim1name='pft', &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%totveg_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing carbonstate_vars %totvegc with atmospheric c13 value'
             do i = bounds%begp,bounds%endp
                if (veg_vp%c3psn(veg_pp%itype(i)) == 1._r8) then
                   this%totveg_patch(i) = c12_carbonstate_vars%totveg_patch(i) * c3_r2
                else
                   this%totveg_patch(i) = c12_carbonstate_vars%totveg_patch(i) * c4_r2
                endif
             end do
          end if
       endif

       !--------------------------------
       ! C14 pft carbon state variables 
       !--------------------------------

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='leafc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leaf_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%leaf_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%leaf_patch(i) /= spval .and. &
                     .not. isnan(this%leaf_patch(i)) ) then
                   this%leaf_patch(i) = c12_carbonstate_vars%leaf_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='leafc_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leaf_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%leaf_storage_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%leaf_storage_patch(i) /= spval .and. &
                     .not. isnan(this%leaf_storage_patch(i)) ) then
                   this%leaf_storage_patch(i) = c12_carbonstate_vars%leaf_storage_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='leafc_xfer_14', xtype=ncd_double,  &
               dim1name='pft',    long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%leaf_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%leaf_xfer_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%leaf_xfer_patch(i) /= spval .and. .not. isnan(this%leaf_xfer_patch(i)) ) then
                   this%leaf_xfer_patch(i) = c12_carbonstate_vars%leaf_xfer_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='frootc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%froot_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%froot_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%froot_patch(i) /= spval .and. &
                     .not. isnan(this%froot_patch(i)) ) then
                   this%froot_patch(i) = c12_carbonstate_vars%froot_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='frootc_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%froot_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%froot_storage_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%froot_storage_patch(i) /= spval .and. &
                     .not. isnan(this%froot_storage_patch(i)) ) then
                   this%froot_storage_patch(i) = c12_carbonstate_vars%froot_storage_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='frootc_xfer_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%froot_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%froot_xfer_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%froot_xfer_patch(i) /= spval .and. &
                     .not. isnan(this%froot_xfer_patch(i)) ) then
                   this%froot_xfer_patch(i) = c12_carbonstate_vars%froot_xfer_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='livestemc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestem_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livestem_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%livestem_patch(i) /= spval .and. .not. isnan(this%livestem_patch(i)) ) then
                   this%livestem_patch(i) = c12_carbonstate_vars%livestem_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='livestemc_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestem_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livestem_storage_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%livestem_storage_patch(i) /= spval .and. .not. isnan(this%livestem_storage_patch(i)) ) then
                   this%livestem_storage_patch(i) = c12_carbonstate_vars%livestem_storage_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='livestemc_xfer_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livestem_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livestem_xfer_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%livestem_xfer_patch(i) /= spval .and. .not. isnan(this%livestem_xfer_patch(i)) ) then
                   this%livestem_xfer_patch(i) = c12_carbonstate_vars%livestem_xfer_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstem_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadstem_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%deadstem_patch(i) /= spval .and. .not. isnan(this%deadstem_patch(i)) ) then
                   this%deadstem_patch(i) = c12_carbonstate_vars%deadstem_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstem_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadstem_storage_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%deadstem_storage_patch(i) /= spval .and. .not. isnan(this%deadstem_storage_patch(i)) ) then
                   this%deadstem_storage_patch(i) = c12_carbonstate_vars%deadstem_storage_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadstemc_xfer_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadstem_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadstem_xfer_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%deadstem_xfer_patch(i) /= spval .and. .not. isnan(this%deadstem_xfer_patch(i)) ) then
                   this%deadstem_xfer_patch(i) = c12_carbonstate_vars%deadstem_xfer_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecroot_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livecroot_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%livecroot_patch(i) /= spval .and. .not. isnan(this%livecroot_patch(i)) ) then
                   this%livecroot_patch(i) = c12_carbonstate_vars%livecroot_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecroot_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livecroot_storage_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%livecroot_storage_patch(i) /= spval .and. .not. isnan(this%livecroot_storage_patch(i)) ) then
                   this%livecroot_storage_patch(i) = c12_carbonstate_vars%livecroot_storage_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='livecrootc_xfer_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%livecroot_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%livecroot_xfer_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%livecroot_xfer_patch(i) /= spval .and. .not. isnan(this%livecroot_xfer_patch(i)) ) then
                   this%livecroot_xfer_patch(i) = c12_carbonstate_vars%livecroot_xfer_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcroot_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadcroot_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%deadcroot_patch(i) /= spval .and. .not. isnan(this%deadcroot_patch(i)) ) then
                   this%deadcroot_patch(i) = c12_carbonstate_vars%deadcroot_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcroot_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%deadcroot_storage_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%deadcroot_storage_patch(i) /= spval .and. .not. isnan(this%deadcroot_storage_patch(i)) ) then
                   this%deadcroot_storage_patch(i) = c12_carbonstate_vars%deadcroot_storage_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_xfer_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%deadcroot_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog) 'initializing this%deadcroot_xfer_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%deadcroot_xfer_patch(i) /= spval .and. .not. isnan(this%deadcroot_xfer_patch(i)) ) then
                   this%deadcroot_xfer_patch(i) = c12_carbonstate_vars%deadcroot_xfer_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='gresp_storage_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%gresp_storage_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%gresp_storage_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%gresp_storage_patch(i) /= spval .and. .not. isnan(this%gresp_storage_patch(i)) ) then
                   this%gresp_storage_patch(i) = c12_carbonstate_vars%gresp_storage_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='gresp_xfer_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%gresp_xfer_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%gresp_xfer_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%gresp_xfer_patch(i) /= spval .and. .not. isnan(this%gresp_xfer_patch(i)) ) then
                   this%gresp_xfer_patch(i) = c12_carbonstate_vars%gresp_xfer_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='cpool_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%pool_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%pool_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%pool_patch(i) /= spval .and. .not. isnan(this%pool_patch(i)) ) then
                   this%pool_patch(i) = c12_carbonstate_vars%pool_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='xsmrpool_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%xsmrpool_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%xsmrpool_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%xsmrpool_patch(i) /= spval .and. .not. isnan(this%xsmrpool_patch(i)) ) then
                   this%xsmrpool_patch(i) = c12_carbonstate_vars%xsmrpool_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='pft_ctrunc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%veg_trunc_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%veg_trunc_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%veg_trunc_patch(i) /= spval .and. .not. isnan(this%veg_trunc_patch(i)) ) then
                   this%veg_trunc_patch(i) = c12_carbonstate_vars%veg_trunc_patch(i) * c14ratio
                endif
             end do
          end if
       end if

       if ( carbon_type == 'c14')  then
          call restartvar(ncid=ncid, flag=flag, varname='totvegc_14', xtype=ncd_double,  &
               dim1name='pft', long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=this%totveg_patch) 
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%totveg_patch with atmospheric c14 value'
             do i = bounds%begp,bounds%endp
                if (this%totveg_patch(i) /= spval .and. .not. isnan(this%totveg_patch(i)) ) then
                   this%totveg_patch(i) = c12_carbonstate_vars%totveg_patch(i) * c14ratio
                endif
             end do
          end if
       endif

    endif  ! .not. use_fates
    

    !--------------------------------
    ! C13 column carbon state variables
    !--------------------------------

    if ( carbon_type == 'c13' ) then
       do k = 1, ndecomp_pools
          varname = trim(decomp_cascade_con%decomp_pool_name_restart(k))//'c_13'
          if (use_vertsoilc) then
             ptr2d => this%decomp_pools_vr_col(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vr", xtype=ncd_double,  &
                  dim1name='column', dim2name='levgrnd', switchdim=.true., &
                  long_name='',  units='', fill_value=spval, &
                  interpinic_flag='interp', readvar=readvar, data=ptr2d)
          else
             ptr1d => this%decomp_pools_vr_col(:,1,k) ! nlevdecomp = 1; so treat as 1D variable
             call restartvar(ncid=ncid, flag=flag, varname=varname, xtype=ncd_double,  &
                  dim1name='column', long_name='',  units='', fill_value=spval, &
                  interpinic_flag='interp' , readvar=readvar, data=ptr1d)
          end if
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%decomp_pools_vr_col with atmospheric c13 value for: '//varname
             do i = bounds%begc,bounds%endc
                do j = 1, nlevdecomp
                   if (c12_carbonstate_vars%decomp_pools_vr_col(i,j,k) /= spval .and. &
                        .not. isnan(this%decomp_pools_vr_col(i,j,k)) ) then
                         this%decomp_pools_vr_col(i,j,k) = c12_carbonstate_vars%decomp_pools_vr_col(i,j,k) * c3_r2
                   endif
                end do
             end do
          end if
       end do
    end if

    if (carbon_type == 'c13') then
       call restartvar(ncid=ncid, flag=flag, varname='seedc_13', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%seed_col) 
       if (flag=='read' .and. .not. readvar) then
          if (c12_carbonstate_vars%seed_col(i) /= spval .and. &
               .not. isnan(c12_carbonstate_vars%seed_col(i)) ) then
             this%seed_col(i) = c12_carbonstate_vars%seed_col(i) * c3_r2
          end if
       end if
    end if

    if (carbon_type == 'c13') then
       if (use_vertsoilc) then
          ptr2d => this%soil_trunc_vr_col
          call restartvar(ncid=ncid, flag=flag, varname="col_ctrunc_c13_vr", xtype=ncd_double,  &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
       else
          ptr1d => this%soil_trunc_vr_col(:,1)
          call restartvar(ncid=ncid, flag=flag, varname="col_ctrunc_c13", xtype=ncd_double,  &
               dim1name='column', long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp' , readvar=readvar, data=ptr1d)
       end if
    end if

    if (carbon_type == 'c13') then
       call restartvar(ncid=ncid, flag=flag, varname='totlitc_13', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totlit_col) 
       if (flag=='read' .and. .not. readvar) then
          if (c12_carbonstate_vars%totlit_col(i) /= spval .and. &
               .not. isnan( c12_carbonstate_vars%totlit_col(i) ) ) then
             this%totlit_col(i) = c12_carbonstate_vars%totlit_col(i) * c3_r2
          end if
       end if
    end if

    if (carbon_type == 'c13') then
       call restartvar(ncid=ncid, flag=flag, varname='totcolc_13', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totcol_col) 
       if (flag=='read' .and. .not. readvar) then
          if (c12_carbonstate_vars%totcol_col(i) /= spval .and. &
               .not. isnan (c12_carbonstate_vars%totcol_col(i) ) ) then
             this%totcol_col(i) = c12_carbonstate_vars%totcol_col(i) * c3_r2
          end if
       end if
    end if

    if (carbon_type == 'c13') then
       call restartvar(ncid=ncid, flag=flag, varname='prod10c_13', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod10_col) 
       if (flag=='read' .and. .not. readvar) then
          if (c12_carbonstate_vars%prod10_col(i) /= spval .and. &
               .not. isnan( c12_carbonstate_vars%prod10_col(i) ) ) then
             this%prod10_col(i) = c12_carbonstate_vars%prod10_col(i) * c3_r2
          endif
       end if
    end if

    if (carbon_type == 'c13') then
       call restartvar(ncid=ncid, flag=flag, varname='prod100c_13', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod100_col) 
       if (flag=='read' .and. .not. readvar) then
          if (c12_carbonstate_vars%prod100_col(i) /= spval .and. &
               .not. isnan( c12_carbonstate_vars%prod100_col(i) ) ) then
             this%prod100_col(i) = c12_carbonstate_vars%prod100_col(i) * c3_r2
          endif
       end if
    endif

    if (carbon_type == 'c13') then
       call restartvar(ncid=ncid, flag=flag, varname='prod1c_13', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod1_col)
       if (flag=='read' .and. .not. readvar) then
          if (c12_carbonstate_vars%prod1_col(i) /= spval .and. &
               .not. isnan( c12_carbonstate_vars%prod1_col(i) ) ) then
             this%prod1_col(i) = c12_carbonstate_vars%prod1_col(i) * c3_r2
          endif
       end if
    end if

    !--------------------------------
    ! C14 column carbon state variables
    !--------------------------------

    if ( carbon_type == 'c14' ) then
       do k = 1, ndecomp_pools
          varname = trim(decomp_cascade_con%decomp_pool_name_restart(k))//'c_14'
          if (use_vertsoilc) then
             ptr2d => this%decomp_pools_vr_col(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vr", xtype=ncd_double,  &
                  dim1name='column', dim2name='levgrnd', switchdim=.true., &
                  long_name='',  units='', fill_value=spval, &
                  interpinic_flag='interp', readvar=readvar, data=ptr2d)
          else
             ptr1d => this%decomp_pools_vr_col(:,1,k) ! nlevdecomp = 1; so treat as 1D variable
             call restartvar(ncid=ncid, flag=flag, varname=varname, xtype=ncd_double,  &
                  dim1name='column', &
                  long_name='',  units='', fill_value=spval, &
                  interpinic_flag='interp' , readvar=readvar, data=ptr1d)
          end if
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing this%decomp_pools_vr_col with atmospheric c14 value for: '//trim(varname)
             do i = bounds%begc,bounds%endc
                do j = 1, nlevdecomp
                   if (c12_carbonstate_vars%decomp_pools_vr_col(i,j,k) /= spval .and. &
                        .not. isnan(c12_carbonstate_vars%decomp_pools_vr_col(i,j,k)) ) then
                         this%decomp_pools_vr_col(i,j,k) = c12_carbonstate_vars%decomp_pools_vr_col(i,j,k) * c3_r2
                   endif
                end do
             end do
          end if
       end do
    end if

    if ( carbon_type == 'c14' ) then
       call restartvar(ncid=ncid, flag=flag, varname='seedc_14', xtype=ncd_double,  &
            dim1name='column', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%seed_col) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%seed_col with atmospheric c14 value'
          do i = bounds%begc,bounds%endc
             if (c12_carbonstate_vars%seed_col(i) /= spval .and. &
                  .not. isnan(c12_carbonstate_vars%seed_col(i)) ) then
                this%seed_col(i) = c12_carbonstate_vars%seed_col(i) * c14ratio
             endif
          end do
       end if
    end if

    if ( carbon_type == 'c14' ) then
       if (use_vertsoilc) then 
          ptr2d => this%soil_trunc_vr_col
          call restartvar(ncid=ncid, flag=flag, varname="col_ctrunc_c14_vr", xtype=ncd_double,  &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
       else
          ptr1d => this%soil_trunc_vr_col(:,1)
          call restartvar(ncid=ncid, flag=flag, varname="col_ctrunc_c14", xtype=ncd_double,  &
               dim1name='column', long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp' , readvar=readvar, data=ptr1d)
       end if
    end if

    if ( carbon_type == 'c14' ) then
       call restartvar(ncid=ncid, flag=flag, varname='totlitc_14', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totlit_col) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%totlit_col with atmospheric c14 value'
          if (c12_carbonstate_vars%totlit_col(i) /= spval .and. &
               .not. isnan(c12_carbonstate_vars%totlit_col(i)) ) then
             this%totlit_col(i) = c12_carbonstate_vars%totlit_col(i) * c14ratio
          endif
       end if
    end if

    if ( carbon_type == 'c14' ) then
       call restartvar(ncid=ncid, flag=flag, varname='totcolc_14', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totcol_col) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%totcol_col with atmospheric c14 value'
          if (c12_carbonstate_vars%totcol_col(i) /= spval .and. &
               .not. isnan(c12_carbonstate_vars%totcol_col(i)) ) then
             this%totcol_col(i) = c12_carbonstate_vars%totcol_col(i) * c14ratio
          endif
       end if
    end if

    if ( carbon_type == 'c14' ) then
       call restartvar(ncid=ncid, flag=flag, varname='prod10c_14', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod10_col) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%prod10_col with atmospheric c14 value'
          if (c12_carbonstate_vars%prod10_col(i) /= spval .and. &
               .not. isnan(c12_carbonstate_vars%prod10_col(i)) ) then
             this%prod10_col(i) = c12_carbonstate_vars%prod10_col(i) * c14ratio
          endif
       end if
    end if

    if ( carbon_type == 'c14' ) then
       call restartvar(ncid=ncid, flag=flag, varname='prod100c_14', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod100_col) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%prod100_col with atmospheric c14 value'
          if (c12_carbonstate_vars%prod100_col(i) /= spval .and. &
               .not. isnan(c12_carbonstate_vars%prod100_col(i)) ) then
             this%prod100_col(i) = c12_carbonstate_vars%prod100_col(i) * c14ratio
          endif
       end if
    endif

    if ( carbon_type == 'c14' ) then
       call restartvar(ncid=ncid, flag=flag, varname='prod1c_14', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod1_col)
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%prod1_col with atmospheric c14 value'
          if (c12_carbonstate_vars%prod1_col(i) /= spval .and. &
               .not. isnan(c12_carbonstate_vars%prod1_col(i)) ) then
             this%prod1_col(i) = c12_carbonstate_vars%prod1_col(i) * c14ratio
          endif
       end if
    end if

    !--------------------------------
    ! Spinup state
    !--------------------------------

    if (carbon_type == 'c12'  .or. carbon_type == 'c14') then
        if (flag == 'write') then
           idata = spinup_state
        end if
        if (carbon_type == 'c12' .or. (carbon_type == 'c14' .and. flag == 'read')) then
        call restartvar(ncid=ncid, flag=flag, varname='spinup_state', xtype=ncd_int,  &
             long_name='Spinup state of the model that wrote this restart file: ' &
             // ' 0 = normal model mode, 1 = AD spinup', units='', &
             interpinic_flag='copy', readvar=readvar,  data=idata)
        end if

        if (flag == 'read') then
           if (readvar) then
              restart_file_spinup_state = idata
           else
              ! assume, for sake of backwards compatibility, that if spinup_state is not in 
              ! the restart file then current model state is the same as prior model state
              restart_file_spinup_state = spinup_state
              if ( masterproc ) then
                 write(iulog,*) ' CNRest: WARNING!  Restart file does not contain info ' &
                      // ' on spinup state used to generate the restart file. '
                 write(iulog,*) '   Assuming the same as current setting: ', spinup_state
              end if
           end if
        end if

        ! now compare the model and restart file spinup states, and either take the 
        ! model into spinup mode or out of it if they are not identical
        ! taking model out of spinup mode requires multiplying each decomposing pool 
        ! by the associated AD factor.
        ! putting model into spinup mode requires dividing each decomposing pool 
        ! by the associated AD factor.
        ! only allow this to occur on first timestep of model run.
        
        if (flag == 'read' .and. spinup_state /= restart_file_spinup_state ) then
           if (spinup_state == 0 .and. restart_file_spinup_state == 1 ) then
              if ( masterproc ) write(iulog,*) ' CNRest: taking SOM pools out of AD spinup mode'
              exit_spinup = .true.
           else if (spinup_state == 1 .and. restart_file_spinup_state == 0 ) then
              if ( masterproc ) write(iulog,*) ' CNRest: taking SOM pools into AD spinup mode'
              enter_spinup = .true.
           else
              call endrun(msg=' CNRest: error in entering/exiting spinup.  spinup_state ' &
                   // ' != restart_file_spinup_state, but do not know what to do'//&
                   errMsg(__FILE__, __LINE__))
           end if
           if (get_nstep() >= 2) then
              call endrun(msg=' CNRest: error in entering/exiting spinup - should occur only when nstep = 1'//&
                   errMsg(__FILE__, __LINE__))
           endif
           do k = 1, ndecomp_pools
              do c = bounds%begc, bounds%endc
                 do j = 1, nlevdecomp
		    if ( exit_spinup ) then
		       m = decomp_cascade_con%spinup_factor(k)
                       if (decomp_cascade_con%spinup_factor(k) > 1) m = m / cnstate_vars%scalaravg_col(c,j)
                    else if ( enter_spinup ) then 
		       m = 1. / decomp_cascade_con%spinup_factor(k)
		       if (decomp_cascade_con%spinup_factor(k) > 1) m = m * cnstate_vars%scalaravg_col(c,j)
		    end if
                    this%decomp_pools_vr_col(c,j,k) = this%decomp_pools_vr_col(c,j,k) * m
                 end do
              end do
           end do
           do i = bounds%begp, bounds%endp
              if (exit_spinup) then 
                 m_veg = spinup_mortality_factor
              else if (enter_spinup) then 
                 m_veg = 1._r8 / spinup_mortality_factor
              end if
              this%deadstem_patch(i)  = this%deadstem_patch(i) * m_veg
              this%deadcroot_patch(i) = this%deadcroot_patch(i) * m_veg
           end do
        end if
     end if

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine SetValues ( this, &
       num_patch, filter_patch, value_patch, &
       num_column, filter_column, value_column)
    !
    ! !DESCRIPTION:
    ! Set carbon state variables
    !
    ! !ARGUMENTS:
    class (carbonstate_type) :: this
    integer , intent(in) :: num_patch
    integer , intent(in) :: filter_patch(:)
    real(r8), intent(in) :: value_patch
    integer , intent(in) :: num_column
    integer , intent(in) :: filter_column(:)
    real(r8), intent(in) :: value_column
    !
    ! !LOCAL VARIABLES:
    integer :: fi,i,j,k,l     ! loop index
    !------------------------------------------------------------------------

    if ( .not. use_fates ) then
       do fi = 1,num_patch
          i = filter_patch(fi)

          this%leaf_patch(i)              = value_patch
          this%leaf_storage_patch(i)      = value_patch
          this%leaf_xfer_patch(i)         = value_patch
          this%froot_patch(i)             = value_patch
          this%froot_storage_patch(i)     = value_patch
          this%froot_xfer_patch(i)        = value_patch
          this%livestem_patch(i)          = value_patch
          this%livestem_storage_patch(i)  = value_patch
          this%livestem_xfer_patch(i)     = value_patch
          this%deadstem_patch(i)          = value_patch
          this%deadstem_storage_patch(i)  = value_patch
          this%deadstem_xfer_patch(i)     = value_patch
          this%livecroot_patch(i)         = value_patch
          this%livecroot_storage_patch(i) = value_patch
          this%livecroot_xfer_patch(i)    = value_patch
          this%deadcroot_patch(i)         = value_patch
          this%deadcroot_storage_patch(i) = value_patch
          this%deadcroot_xfer_patch(i)    = value_patch
          this%gresp_storage_patch(i)      = value_patch
          this%gresp_xfer_patch(i)         = value_patch
          this%pool_patch(i)              = value_patch
          this%xsmrpool_patch(i)           = value_patch
          this%veg_trunc_patch(i)             = value_patch
          this%dispveg_patch(i)           = value_patch
          this%storveg_patch(i)           = value_patch
          this%totveg_patch(i)            = value_patch
          this%totpft_patch(i)            = value_patch
          this%woodc_patch(i)              = value_patch
          this%totveg_abg_patch(i)        = value_patch
       end do
       
       if ( crop_prog ) then
          do fi = 1,num_patch
             i = filter_patch(fi)
             this%grain_patch(i)            = value_patch
             this%grain_storage_patch(i)    = value_patch
             this%grain_xfer_patch(i)       = value_patch
             this%cropseed_deficit_patch(i) = value_patch
          end do
       endif
    endif ! .not. use_fates
    
    do fi = 1,num_column
       i = filter_column(fi)
       this%cwd_col(i)       = value_column
       this%veg_trunc_col(i)     = value_column
       this%totlit_col(i)    = value_column
       this%totsom_col(i)    = value_column
       this%totecosys_col(i) = value_column
       this%totcol_col(i)    = value_column
       this%rootc_col(i)      = value_column
       this%totveg_col(i)    = value_column
       this%leafc_col(i)      = value_column
       this%deadstemc_col(i)  = value_column
       this%fuelc_col(i)      = value_column
       this%fuelc_crop_col(i) = value_column
       this%totlit_1m_col(i) = value_column
       this%totsom_1m_col(i) = value_column
    end do

    do j = 1,nlevdecomp_full
       do fi = 1,num_column
          i = filter_column(fi)
          this%soil_trunc_vr_col(i,j) = value_column
       end do
    end do

    do k = 1, ndecomp_pools
       do fi = 1,num_column
          i = filter_column(fi)
          this%decomp_pools_col(i,k) = value_column
          this%decomp_pools_1m_col(i,k) = value_column
       end do
    end do

    do j = 1,nlevdecomp_full
       do k = 1, ndecomp_pools
          do fi = 1,num_column
             i = filter_column(fi)
             this%decomp_pools_vr_col(i,j,k) = value_column
          end do
       end do
    end do

  end subroutine SetValues

  !-----------------------------------------------------------------------
  subroutine ZeroDwt( this, bounds )
    !
    ! !DESCRIPTION
    ! Initialize variables needed for dynamic land use.
    !
    ! !ARGUMENTS:
    class(carbonstate_type) :: this
    type(bounds_type), intent(in)  :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer  :: p          ! indices
    !-----------------------------------------------------------------------

    do p = bounds%begp,bounds%endp
       this%dispveg_patch(p)   = 0._r8
       this%storveg_patch(p)   = 0._r8
       this%totpft_patch(p)    = 0._r8
    end do

  end subroutine ZeroDwt

  !-----------------------------------------------------------------------
  subroutine Summary(this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, perform patch and column-level carbon summary calculations
    !
    ! !USES:
    use clm_varctl       , only: iulog
    use clm_time_manager , only: get_step_size
    use clm_varcon       , only: secspday
    use clm_varpar       , only: nlevdecomp, ndecomp_pools, nlevdecomp_full
    !
    ! !ARGUMENTS:
    class(carbonstate_type) :: this
    type(bounds_type)      , intent(in)    :: bounds          
    integer                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    !
    ! !LOCAL VARIABLES:
    real(r8) :: nfixlags, dtime ! temp variables for making lagged npp
    integer  :: c,p,j,k,l       ! indices
    integer  :: fp,fc           ! lake filter indices
    real(r8) :: maxdepth        ! depth to integrate soil variables
    integer  :: nlev
    real(r8) :: cropseedc_deficit_col(bounds%begc:bounds%endc)
    !-----------------------------------------------------------------------

    ! calculate patch -level summary of carbon state

    if (use_fates) return

    call NutrientStatePatchSummary( this, bounds, num_soilp, filter_soilp)
    call NutrientStateColumnSummary(this, bounds, num_soilc, filter_soilc)

    ! add carbon-specific pools to patch-level variables
    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! stored vegetation carbon, excluding cpool (STORVEGC)
       this%storveg_patch(p) =          &
            this%storveg_patch(p)       + &
            this%gresp_storage_patch(p) + &
            this%gresp_xfer_patch(p)

       ! total vegetation carbon, excluding cpool (TOTVEGC)
       this%totveg_patch(p) =       &
            this%dispveg_patch(p) + &
            this%storveg_patch(p)

       ! total pft-level carbon, including xsmrpool, ctrunc
       this%totpft_patch(p) =        &
            this%totveg_patch(p)   + &
            this%xsmrpool_patch(p) + &
            this%veg_trunc_patch(p)

       ! wood C
       this%woodc_patch(p) = &
            this%deadstem_patch(p)    + &
            this%livestem_patch(p)    + &
            this%deadcroot_patch(p)   + &
            this%livecroot_patch(p)
    end do

    ! aggregate patch-level state variables to column-level
    call p2c(bounds, num_soilc, filter_soilc, &
         this%totpft_patch(bounds%begp:bounds%endp), &
         this%totpft_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%totveg_patch(bounds%begp:bounds%endp), &
         this%totveg_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%totveg_abg_patch(bounds%begp:bounds%endp), &
         this%totveg_abg_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%cropseed_deficit_patch(bounds%begp:bounds%endp), &
         cropseedc_deficit_col(bounds%begc:bounds%endc))

    ! column level summary
    do fc = 1,num_soilc
       c = filter_soilc(fc)

       ! total product carbon
       this%totprod_col(c) = &
            this%prod10_col(c)  + &
            this%prod100_col(c) + &
            this%prod1_col(c) 

       ! total ecosystem carbon, including veg but excluding cpool (TOTECOSYSC)
       this%totecosys_col(c) = &
            this%cwd_col(c)     + &
            this%totlit_col(c)  + &
            this%totsom_col(c)  + &
            this%totprod_col(c) + &
            this%totveg_col(c)

       ! total column carbon, including veg and cpool (TOTCOLC)
       ! adding col_ctrunc, seedc
       this%totcol_col(c) = &
            this%totpft_col(c)  + &
            this%cwd_col(c)     + &
            this%totlit_col(c)  + &
            this%totsom_col(c)  + &
            this%prod1_col(c)   + &
            this%veg_trunc_col(c)   + &
            cropseedc_deficit_col(c)

       this%totabg_col(c) = &
            this%totpft_col(c)  + &
            this%totprod_col(c) + &
            this%seed_col(c)    + &
            this%veg_trunc_col(c)    
    end do

  end subroutine Summary

  !-----------------------------------------------------------------------
  subroutine DynamicPatchAdjustments( this, &
       bounds,                              &
       num_filterp_with_inactive,           &
       filterp_with_inactive,               &
       prior_weights,                       &
       patch_state_updater,                 &
       dwt_leafc_seed,                      &
       dwt_deadstemc_seed,                  &
       conv_cflux,                          &
       dwt_frootc_to_litter,                &
       dwt_livecrootc_to_litter,            &
       dwt_deadcrootc_to_litter,            &
       prod10_cflux,                        &
       prod100_cflux,                       &
       crop_product_cflux                   &
       )
    !
    ! !DESCRIPTION:
    ! Adjust state variables and compute associated fluxes when patch areas change due to
    ! dynamic landuse
    !
    ! !USES:
    use pftvarcon          , only : pconv, pprod10, pprod100
    use dynPriorWeightsMod , only : prior_weights_type
    use landunit_varcon    , only : istsoil, istcrop
    use CNComputeSeedMod   , only : ComputeSeedAmounts
    !
    ! !ARGUMENTS:
    class(carbonstate_type)        , intent(inout) :: this
    type(bounds_type)              , intent(in)    :: bounds
    integer                        , intent(in)    :: num_filterp_with_inactive ! number of points in filterp_with_inactive
    integer                        , intent(in)    :: filterp_with_inactive(:) ! patch filter that includes inactive points
    type(prior_weights_type)       , intent(in)    :: prior_weights
    type(patch_state_updater_type) , intent(in)    :: patch_state_updater
    real(r8)                       , intent(inout) :: dwt_leafc_seed           (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_deadstemc_seed       (bounds%begp:)
    real(r8)                       , intent(inout) :: conv_cflux               (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_frootc_to_litter     (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_livecrootc_to_litter (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_deadcrootc_to_litter (bounds%begp:)
    real(r8)                       , intent(inout) :: prod10_cflux             (bounds%begp:)
    real(r8)                       , intent(inout) :: prod100_cflux            (bounds%begp:)
    real(r8)                       , intent(inout) :: crop_product_cflux       (bounds%begp:)
    !
    ! !LOCAL VARIABLES:
    integer                     :: begp, endp
    integer                     :: l, c, p
    logical                     :: old_weight_was_zero(bounds%begp:bounds%endp)
    logical                     :: patch_grew(bounds%begp:bounds%endp)

    ! The following are only set for growing patches:
    real(r8)                    :: seed_leaf_patch(bounds%begp:bounds%endp)
    real(r8)                    :: seed_leaf_storage_patch(bounds%begp:bounds%endp)
    real(r8)                    :: seed_leaf_xfer_patch(bounds%begp:bounds%endp)
    real(r8)                    :: seed_deadstem_patch(bounds%begp:bounds%endp)

    real(r8)                    :: wood_product_cflux(bounds%begp:bounds%endp)
    real(r8)                    :: deadstem_patch_temp(bounds%begp:bounds%endp)

    character(len=*), parameter :: subname = 'CStateDynamicPatchAdjustments'
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    SHR_ASSERT_ALL((ubound(dwt_leafc_seed           ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_deadstemc_seed       ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(conv_cflux               ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_frootc_to_litter     ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_livecrootc_to_litter ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_deadcrootc_to_litter ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(prod10_cflux             ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(prod100_cflux            ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(crop_product_cflux       ) == (/endp/)), errMsg(__FILE__, __LINE__))
   
    call NutrientStateDynamicPatchAdjustments( &
       this,                                   &
       bounds,                                 &
       num_filterp_with_inactive,              &
       filterp_with_inactive,                  &
       prior_weights,                          &
       patch_state_updater,                    &
       this%species,                           &
       dwt_leafc_seed,                         &
       dwt_deadstemc_seed,                     &
       conv_cflux,                             &
       dwt_frootc_to_litter,                   &
       dwt_livecrootc_to_litter,               &
       dwt_deadcrootc_to_litter,               &
       prod10_cflux,                           &
       prod100_cflux,                          &
       crop_product_cflux                      &
    )

    ! 1) GRESP_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                   &
         bounds                                                  , &
         num_filterp_with_inactive                               , &
         filterp_with_inactive                                   , &
         var               = this%gresp_storage_patch(begp:endp) , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 2) GRESP_XFER_STORAGE
    call patch_state_updater%update_patch_state(                   &
         bounds                                                  , &
         num_filterp_with_inactive                               , &
         filterp_with_inactive                                   , &
         var               = this%gresp_xfer_patch(begp:endp)    , &
         flux_out_grc_area = conv_cflux(begp:endp))

    ! 3) XSMRPOOL_PATCH
    call patch_state_updater%update_patch_state(                   &
         bounds                                                  , &
         num_filterp_with_inactive                               , &
         filterp_with_inactive                                   , &
         var               = this%xsmrpool_patch(begp:endp)      , &
         flux_out_grc_area = conv_cflux(begp:endp))

  end subroutine DynamicPatchAdjustments
  
  !-----------------------------------------------------------------------
  subroutine DynamicColumnAdjustments( this, &
       bounds, clump_index, column_state_updater)
    !
    ! !DESCRIPTION:
    ! Adjust state variables and compute associated fluxes when patch areas change due to
    ! dynamic landuse
    !
    ! !USES:
    use dynPriorWeightsMod , only : prior_weights_type
    use landunit_varcon    , only : istsoil, istcrop
    use dynColumnStateUpdaterMod, only : column_state_updater_type
    !
    ! !ARGUMENTS:
    class(carbonstate_type)         , intent(inout) :: this
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: clump_index
    type(column_state_updater_type) , intent(in)    :: column_state_updater
    !
    ! !LOCAL VARIABLES:
    integer                     :: l, j
    integer                     :: begc, endc
    real(r8)                    :: adjustment_one_level(bounds%begc:bounds%endc)

    character(len=*), parameter :: subname = 'CStateDynamicColumnAdjustments'
    !-----------------------------------------------------------------------

    begc = bounds%begc
    endc = bounds%endc

    this%dyn_bal_adjustments_col(begc:endc) = 0._r8

    do l = 1, ndecomp_pools
       do j = 1, nlevdecomp
          call column_state_updater%update_column_state_no_special_handling( &
               bounds      = bounds,                                         &
               clump_index = clump_index,                                    &
               var         = this%decomp_pools_vr_col(begc:endc, j, l),     &
               adjustment  = adjustment_one_level(begc:endc))

          this%dyn_bal_adjustments_col(begc:endc) = &
               this%dyn_bal_adjustments_col(begc:endc) + &
               adjustment_one_level(begc:endc) * dzsoi_decomp(j)

       end do
    end do

    do j = 1, nlevdecomp
       call column_state_updater%update_column_state_no_special_handling( &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = this%soil_trunc_vr_col(begc:endc,j),     &
            adjustment  = adjustment_one_level(begc:endc))

       this%dyn_bal_adjustments_col(begc:endc) = &
            this%dyn_bal_adjustments_col(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)

    end do

  end subroutine DynamicColumnAdjustments

end module CNCarbonStateType
