module PhosphorusStateType

#include "shr_assert.h"

  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_infnan_mod         , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use clm_varpar             , only : ndecomp_cascade_transitions, ndecomp_pools, nlevcan
  use clm_varpar             , only : nlevdecomp_full, nlevdecomp, crop_prog
  use clm_varcon             , only : spval, ispval, dzsoi_decomp, zisoi, zsoi
  use landunit_varcon        , only : istcrop, istsoil 
  use clm_varctl             , only : use_nitrif_denitrif, use_vertsoilc, use_century_decomp
  use clm_varctl             , only : iulog, override_bgc_restart_mismatch_dump, spinup_state
  use decompMod              , only : bounds_type
  use pftvarcon              , only : npcropmin, nstor
  use CNDecompCascadeConType , only : decomp_cascade_con
  use VegetationPropertiesType         , only : veg_vp
  use abortutils             , only : endrun
  use spmdMod                , only : masterproc 
  use LandunitType           , only : lun_pp                
  use ColumnType             , only : col_pp                
  use VegetationType              , only : veg_pp
  use clm_varctl             , only : nu_com, use_crop
  ! soil phosphorus initialization Qing Z. 2017
  use pftvarcon              , only : VMAX_MINSURF_P_vr, KM_MINSURF_P_vr
  use soilorder_varcon       , only : smax, ks_sorption
  use dynPatchStateUpdaterMod      , only : patch_state_updater_type
  use CNSpeciesMod           , only : NUTRIENT_SPECIES_P
  use NutrientStateType      , only : nutrientstate_type, NutrientStateInitAllocate
  use NutrientStateType      , only : NutrientStateInitHistory
  use NutrientStateType      , only : NutrientStateDynamicPatchAdjustments
  use NutrientStateType      , only : NutrientStateRestart
  use NutrientStateType      , only : NutrientStatePatchSummary
  use NutrientStateType      , only : NutrientStateColumnSummary
  use CNSpeciesMod           , only : species_from_string, species_name_from_string
  use CNSpeciesMod           , only : species_history_name_prefix_from_string
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private

  real(r8) , parameter :: ppool_seed_param     = 0.01_r8

  type, public, extends (nutrientstate_type) :: phosphorusstate_type

     real(r8), pointer :: retransp_patch        (:)     ! patch (gP/m2) plant pool of retranslocated P
     real(r8), pointer :: plant_p_buffer_patch  (:)     ! patch (gP/m2) pft-level abstract p storage
     real(r8), pointer :: solutionp_vr_col      (:,:)   ! col (gP/m3) vertically-resolved soil solution P
     real(r8), pointer :: labilep_vr_col        (:,:)   ! col (gP/m3) vertically-resolved soil labile mineral P
     real(r8), pointer :: secondp_vr_col        (:,:)   ! col (gP/m3) vertically-resolved soil secondary mineralP
     real(r8), pointer :: occlp_vr_col          (:,:)   ! col (gP/m3) vertically-resolved soil occluded mineral P
     real(r8), pointer :: primp_vr_col          (:,:)   ! col (gP/m3) vertically-resolved soil parimary mineral P
     real(r8), pointer :: sminp_vr_col          (:,:)   ! col (gP/m3) vertically-resolved soil parimary mineral P

     ! summary (diagnostic) state variables, not involved in mass balance
     real(r8), pointer :: sminp_col             (:)     ! col (gP/m2) soil mineral P
     real(r8), pointer :: solutionp_col         (:)     ! col (gP/m2) soil solution P
     real(r8), pointer :: labilep_col           (:)     ! col (gP/m2) soil labile mineral P
     real(r8), pointer :: secondp_col           (:)     ! col (gP/m2) soil secondary mineralP
     real(r8), pointer :: occlp_col             (:)     ! col (gP/m2) soil occluded mineral P
     real(r8), pointer :: primp_col             (:)     ! col (gP/m2) soil parimary mineral P

     real(r8), pointer :: solutionp_vr_col_cur  (:,:)
     real(r8), pointer :: solutionp_vr_col_prev (:,:)
     real(r8), pointer :: labilep_vr_col_cur    (:,:)
     real(r8), pointer :: labilep_vr_col_prev   (:,:)
     real(r8), pointer :: secondp_vr_col_cur    (:,:)
     real(r8), pointer :: secondp_vr_col_prev   (:,:)
     real(r8), pointer :: occlp_vr_col_cur      (:,:)
     real(r8), pointer :: occlp_vr_col_prev     (:,:)
     real(r8), pointer :: primp_vr_col_cur      (:,:)
     real(r8), pointer :: primp_vr_col_prev     (:,:)

   contains

     procedure , public  :: Init   
     procedure , public  :: Restart
     procedure , public  :: SetValues
     procedure , public  :: ZeroDWT
     procedure , public  :: Summary
     procedure , public  :: DynamicPatchAdjustments
     procedure , public  :: DynamicColumnAdjustments
     procedure , private :: InitAllocate
     procedure , private :: InitHistory  
     procedure , private :: InitCold

  end type phosphorusstate_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds,                           &
       leafc_patch, leafc_storage_patch, frootc_patch, frootc_storage_patch, &
       deadstemc_patch, decomp_cpools_vr_col, decomp_cpools_col, decomp_pools_1m_col)

    class(phosphorusstate_type)         :: this
    type(bounds_type) , intent(in)    :: bounds  
    real(r8)          , intent(in)    :: leafc_patch          (bounds%begp:)
    real(r8)          , intent(in)    :: leafc_storage_patch  (bounds%begp:)
    real(r8)          , intent(in)    :: frootc_patch         (bounds%begp:)
    real(r8)          , intent(in)    :: frootc_storage_patch (bounds%begp:)
    real(r8)          , intent(in)    :: deadstemc_patch      (bounds%begp:)
    real(r8)          , intent(in)    :: decomp_cpools_vr_col (bounds%begc:, 1:, 1:)
    real(r8)          , intent(in)    :: decomp_cpools_col    (bounds%begc:, 1:)
    real(r8)          , intent(in)    :: decomp_pools_1m_col (bounds%begc:, 1:)

    this%species             = species_from_string('p')
    this%name                = species_name_from_string('p')
    this%history_name_prefix = species_history_name_prefix_from_string('p')
    this%restart_name        = 'p'

    call this%InitAllocate (bounds )

    call this%InitHistory (bounds)

    call this%InitCold ( bounds, leafc_patch, leafc_storage_patch, &
         frootc_patch, frootc_storage_patch, deadstemc_patch, &
         decomp_cpools_vr_col, decomp_cpools_col, decomp_pools_1m_col)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !ARGUMENTS:
    class (phosphorusstate_type) :: this
    type(bounds_type) , intent(in) :: bounds  
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

    allocate(this%retransp_patch        (begp:endp))                   ; this%retransp_patch        (:)   = nan
    allocate(this%plant_p_buffer_patch  (begp:endp))                   ; this%plant_p_buffer_patch  (:)   = nan

    allocate(this%solutionp_vr_col      (begc:endc,1:nlevdecomp_full)) ; this%solutionp_vr_col      (:,:) = nan
    allocate(this%labilep_vr_col        (begc:endc,1:nlevdecomp_full)) ; this%labilep_vr_col        (:,:) = nan
    allocate(this%secondp_vr_col        (begc:endc,1:nlevdecomp_full)) ; this%secondp_vr_col        (:,:) = nan
    allocate(this%occlp_vr_col          (begc:endc,1:nlevdecomp_full)) ; this%occlp_vr_col          (:,:) = nan
    allocate(this%primp_vr_col          (begc:endc,1:nlevdecomp_full)) ; this%primp_vr_col          (:,:) = nan
    allocate(this%sminp_vr_col          (begc:endc,1:nlevdecomp_full)) ; this%sminp_vr_col          (:,:) = nan

    allocate(this%solutionp_col         (begc:endc))                   ; this%solutionp_col         (:)   = nan
    allocate(this%labilep_col           (begc:endc))                   ; this%labilep_col           (:)   = nan
    allocate(this%secondp_col           (begc:endc))                   ; this%secondp_col           (:)   = nan
    allocate(this%occlp_col             (begc:endc))                   ; this%occlp_col             (:)   = nan
    allocate(this%primp_col             (begc:endc))                   ; this%primp_col             (:)   = nan
    allocate(this%sminp_col             (begc:endc))                   ; this%sminp_col             (:)   = nan

    allocate(this%solutionp_vr_col_cur  (begc:endc,1:nlevdecomp_full)) ; this%solutionp_vr_col_cur  (:,:) = nan
    allocate(this%solutionp_vr_col_prev (begc:endc,1:nlevdecomp_full)) ; this%solutionp_vr_col_prev (:,:) = nan
    allocate(this%labilep_vr_col_cur    (begc:endc,1:nlevdecomp_full)) ; this%labilep_vr_col_cur    (:,:) = nan
    allocate(this%labilep_vr_col_prev   (begc:endc,1:nlevdecomp_full)) ; this%labilep_vr_col_prev   (:,:) = nan
    allocate(this%secondp_vr_col_cur    (begc:endc,1:nlevdecomp_full)) ; this%secondp_vr_col_cur    (:,:) = nan
    allocate(this%secondp_vr_col_prev   (begc:endc,1:nlevdecomp_full)) ; this%secondp_vr_col_prev   (:,:) = nan
    allocate(this%occlp_vr_col_cur      (begc:endc,1:nlevdecomp_full)) ; this%occlp_vr_col_cur      (:,:) = nan
    allocate(this%occlp_vr_col_prev     (begc:endc,1:nlevdecomp_full)) ; this%occlp_vr_col_prev     (:,:) = nan
    allocate(this%primp_vr_col_cur      (begc:endc,1:nlevdecomp_full)) ; this%primp_vr_col_cur      (:,:) = nan
    allocate(this%primp_vr_col_prev     (begc:endc,1:nlevdecomp_full)) ; this%primp_vr_col_prev     (:,:) = nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! add history fields for all CN variables, always set as default='inactive'
    !
    ! !USES:
    use clm_varpar , only : ndecomp_cascade_transitions, ndecomp_pools
    use clm_varpar , only : nlevdecomp, nlevdecomp_full,crop_prog, nlevgrnd
    use histFileMod, only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp 
    use decompMod  , only : bounds_type
    !
    ! !ARGUMENTS:
    class(phosphorusstate_type) :: this
    type(bounds_type)         , intent(in) :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer           :: k,l,ii,jj 
    character(10)     :: active
    character(8)      :: vr_suffix
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

    this%retransp_patch(begp:endp) = spval
    call hist_addfld1d (fname='RETRANSP', units='gP/m^2', &
         avgflag='A', long_name='plant pool of retranslocated P', &
         ptr_patch=this%retransp_patch)

    this%plant_p_buffer_patch(begp:endp) = spval
    call hist_addfld1d (fname='PLANTP_BUFFER', units='gP/m^2', &
            avgflag='A', long_name='plant phosphorus stored as buffer', &
            ptr_col=this%plant_p_buffer_patch,default='inactive')

    if ( nlevdecomp_full > 1 ) then

       this%sminp_col(begc:endc) = spval
       call hist_addfld1d (fname='SMINP', units='gP/m^2', &
            avgflag='A', long_name='soil mineral P', &
            ptr_col=this%sminp_col)

    endif

    ! add suffix if number of soil decomposition depths is greater than 1
    if (nlevdecomp > 1) then
       vr_suffix = "_vr"
    else 
       vr_suffix = ""
    endif

    this%solutionp_vr_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='SOLUTIONP'//trim(vr_suffix), units='gP/m^3',  type2d='levdcmp', &
         avgflag='A', long_name='soil solution P (vert. res.)', &
         ptr_col=this%solutionp_vr_col)

    this%labilep_vr_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='LABILEP'//trim(vr_suffix), units='gP/m^3',  type2d='levdcmp', &
         avgflag='A', long_name='soil labile P (vert. res.)', &
         ptr_col=this%labilep_vr_col)

    this%secondp_vr_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='SECONDP'//trim(vr_suffix), units='gP/m^3',  type2d='levdcmp', &
         avgflag='A', long_name='soil secondary P (vert. res.)', &
         ptr_col=this%secondp_vr_col)

    this%occlp_vr_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='OCCLP'//trim(vr_suffix), units='gP/m^3',  type2d='levdcmp', &
         avgflag='A', long_name='soil occluded P (vert. res.)', &
         ptr_col=this%occlp_vr_col)

    this%primp_vr_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='PRIMP'//trim(vr_suffix), units='gP/m^3',  type2d='levdcmp', &
         avgflag='A', long_name='soil primary P (vert. res.)', &
         ptr_col=this%primp_vr_col)

    this%sminp_vr_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='SMINP'//trim(vr_suffix), units='gP/m^3',  type2d='levdcmp', &
         avgflag='A', long_name='soil mineral P (vert. res.)', &
         ptr_col=this%sminp_vr_col)

    if ( nlevdecomp_full > 1 ) then

       this%solutionp_col(begc:endc) = spval
       call hist_addfld1d (fname='SOLUTIONP', units='gP/m^2', &
            avgflag='A', long_name='soil solution P', &
            ptr_col=this%solutionp_col)

       this%labilep_col(begc:endc) = spval
       call hist_addfld1d (fname='LABILEP', units='gP/m^2', &
            avgflag='A', long_name='soil Labile P', &
            ptr_col=this%labilep_col)

       this%secondp_col(begc:endc) = spval
       call hist_addfld1d (fname='SECONDP', units='gP/m^2', &
            avgflag='A', long_name='soil secondary P', &
            ptr_col=this%secondp_col)

       this%occlp_col(begc:endc) = spval
       call hist_addfld1d (fname='OCCLP', units='gP/m^2', &
            avgflag='A', long_name='soil occluded P', &
            ptr_col=this%occlp_col)

       this%primp_col(begc:endc) = spval
       call hist_addfld1d (fname='PRIMP', units='gP/m^2', &
            avgflag='A', long_name='soil primary P', &
            ptr_col=this%primp_col)
    endif

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, &
       leafc_patch, leafc_storage_patch, frootc_patch, frootc_storage_patch, &
       deadstemc_patch, decomp_cpools_vr_col, decomp_cpools_col, decomp_pools_1m_col)
    !
    ! !DESCRIPTION:
    ! Initializes time varying variables used only in coupled carbon-phosphorus mode (CN):
    !
    ! !USES:
    use clm_varpar     , only : crop_prog
    use decompMod      , only : bounds_type
    use pftvarcon      , only : noveg, npcropmin
    !
    ! !ARGUMENTS:
    class(phosphorusstate_type)      :: this
    type(bounds_type) , intent(in) :: bounds  
    real(r8)          , intent(in) :: leafc_patch(bounds%begp:)
    real(r8)          , intent(in) :: leafc_storage_patch(bounds%begp:)
    real(r8)          , intent(in) :: frootc_patch(bounds%begp:)
    real(r8)          , intent(in) :: frootc_storage_patch(bounds%begp:)
    real(r8)          , intent(in) :: deadstemc_patch(bounds%begp:)
    real(r8)          , intent(in) :: decomp_cpools_vr_col(bounds%begc:,:,:)
    real(r8)          , intent(in) :: decomp_cpools_col(bounds%begc:,:)
    real(r8)          , intent(in) :: decomp_pools_1m_col(bounds%begc:,:)
    !
    ! !LOCAL VARIABLES:
    integer :: fc,fp,g,l,c,p,j,k                       ! indices
    integer :: num_special_col                         ! number of good values in special_col filter
    integer :: num_special_patch                         ! number of good values in special_patch filter
    integer :: special_col   (bounds%endc-bounds%begc+1) ! special landunit filter - columns
    integer :: special_patch (bounds%endp-bounds%begp+1) ! special landunit filter - patches
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(leafc_patch)          == (/bounds%endp/)),                               errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(leafc_storage_patch)  == (/bounds%endp/)),                               errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(frootc_patch)         == (/bounds%endp/)),                               errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(frootc_storage_patch) == (/bounds%endp/)),                               errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(deadstemc_patch)      == (/bounds%endp/)),                               errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(decomp_cpools_col)    == (/bounds%endc,ndecomp_pools/)),                 errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(decomp_pools_1m_col) == (/bounds%endc,ndecomp_pools/)),                 errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(decomp_cpools_vr_col) == (/bounds%endc,nlevdecomp_full,ndecomp_pools/)), errMsg(__FILE__, __LINE__))

    ! Set column filters

    num_special_patch = 0
    do p = bounds%begp,bounds%endp
       l = veg_pp%landunit(p)
       if (lun_pp%ifspecial(l)) then
          num_special_patch = num_special_patch + 1
          special_patch(num_special_patch) = p
       end if
    end do

    ! Set patch filters

    num_special_col = 0
    do c = bounds%begc, bounds%endc
       l = col_pp%landunit(c)
       if (lun_pp%ifspecial(l)) then
          num_special_col = num_special_col + 1
          special_col(num_special_col) = c
       end if
    end do

    !-------------------------------------------
    ! initialize pft-level variables
    !-------------------------------------------

    do p = bounds%begp,bounds%endp

       l = veg_pp%landunit(p)
       if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then

          if (veg_pp%itype(p) == noveg) then
             this%leaf_patch(p) = 0._r8
             this%leaf_storage_patch(p) = 0._r8
          else
             this%leaf_patch(p)         = leafc_patch(p)         / veg_vp%leafcp(veg_pp%itype(p))
             this%leaf_storage_patch(p) = leafc_storage_patch(p) / veg_vp%leafcp(veg_pp%itype(p))
          end if

          this%leaf_xfer_patch(p)        = 0._r8
          if ( crop_prog )then
             this%grain_patch(p)            = 0._r8
             this%grain_storage_patch(p)    = 0._r8
             this%grain_xfer_patch(p)       = 0._r8
             this%cropseed_deficit_patch(p) = 0._r8
          end if
          this%froot_patch(p)            = 0._r8
          this%froot_storage_patch(p)    = 0._r8
          this%froot_xfer_patch(p)       = 0._r8
          this%livestem_patch(p)         = 0._r8
          this%livestem_storage_patch(p) = 0._r8
          this%livestem_xfer_patch(p)    = 0._r8

          ! tree types need to be initialized with some stem mass so that
          ! roughness length is not zero in canopy flux calculation

          if (veg_vp%woody(veg_pp%itype(p)) == 1._r8) then
             this%deadstem_patch(p) = deadstemc_patch(p) / veg_vp%deadwdcp(veg_pp%itype(p))
          else
             this%deadstem_patch(p) = 0._r8
          end if
          
          if (nu_com .ne. 'RD') then
              ! ECA competition calculate root NP uptake as a function of fine root biomass
              ! better to initialize root CNP pools with a non-zero value
              if (veg_pp%itype(p) .ne. noveg) then
                 this%froot_patch(p) = frootc_patch(p) / veg_vp%frootcp(veg_pp%itype(p))
                 this%froot_storage_patch(p) = frootc_storage_patch(p) / veg_vp%frootcp(veg_pp%itype(p))
              end if
          end if
           
          this%deadstem_storage_patch(p)  = 0._r8
          this%deadstem_xfer_patch(p)     = 0._r8
          this%livecroot_patch(p)         = 0._r8
          this%livecroot_storage_patch(p) = 0._r8
          this%livecroot_xfer_patch(p)    = 0._r8
          this%deadcroot_patch(p)         = 0._r8
          this%deadcroot_storage_patch(p) = 0._r8
          this%deadcroot_xfer_patch(p)    = 0._r8
          this%retransp_patch(p)           = 0._r8
          this%pool_patch(p)              = 0._r8
          if (nstor(veg_pp%itype(p)) .gt. 1e-6_r8) then
              this%pool_patch(p)          = 1.0_r8
          end if
          this%veg_trunc_patch(p)             = 0._r8
          this%dispveg_patch(p)           = 0._r8
          this%storveg_patch(p)           = 0._r8
          this%totveg_patch(p)            = 0._r8
          this%totpft_patch(p)            = 0._r8
       end if
       this%plant_p_buffer_patch(p)= 1.e-4_r8 
    end do

    !-------------------------------------------
    ! initialize column-level variables
    !-------------------------------------------

    do c = bounds%begc, bounds%endc
       l = col_pp%landunit(c)
       if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then

          ! column phosphorus state variables
          this%veg_trunc_col(c) = 0._r8
          this%sminp_col(c) = 0._r8
          do j = 1, nlevdecomp
             do k = 1, ndecomp_pools
                this%decomp_pools_vr_col(c,j,k) = decomp_cpools_vr_col(c,j,k) / decomp_cascade_con%initial_cp_ratio(k)
             end do
             this%sminp_vr_col(c,j) = 0._r8
             this%soil_trunc_vr_col(c,j) = 0._r8
          end do
          if ( nlevdecomp > 1 ) then
             do j = nlevdecomp+1, nlevdecomp_full
                do k = 1, ndecomp_pools
                   this%decomp_pools_vr_col(c,j,k) = 0._r8
                end do
                this%sminp_vr_col(c,j) = 0._r8
                this%soil_trunc_vr_col(c,j) = 0._r8
             end do
          end if
          do k = 1, ndecomp_pools
             this%decomp_pools_col(c,k)    = decomp_cpools_col(c,k)    / decomp_cascade_con%initial_cp_ratio(k)
             this%decomp_pools_1m_col(c,k) = decomp_pools_1m_col(c,k) / decomp_cascade_con%initial_cp_ratio(k)
          end do

          do j = 1, nlevdecomp_full
             this%solutionp_vr_col(c,j) = 0._r8
             this%labilep_vr_col(c,j)   = 0._r8
             this%secondp_vr_col(c,j)   = 0._r8
             this%occlp_vr_col(c,j)     = 0._r8
             this%primp_vr_col(c,j)     = 0._r8
             this%sminp_vr_col(c,j) = 0._r8
          end do
          this%solutionp_col(c) = 0._r8
          this%labilep_col(c)   = 0._r8
          this%secondp_col(c)   = 0._r8
          this%occlp_col(c)     = 0._r8
          this%primp_col(c)     = 0._r8

          this%totlit_col(c)    = 0._r8
          this%totsom_col(c)    = 0._r8
          this%totlit_1m_col(c) = 0._r8
          this%totsom_1m_col(c) = 0._r8
          this%totecosys_col(c) = 0._r8
          this%totcol_col(c)    = 0._r8
          this%cwd_col(c)       = 0._r8

          ! dynamic landcover state variables
          this%seed_col(c)         = 0._r8
          this%prod1_col(c)        = 0._r8
          this%prod10_col(c)       = 0._r8
          this%prod100_col(c)      = 0._r8
          this%totprod_col(c)      = 0._r8
       end if
    end do

    ! now loop through special filters and explicitly set the variables that
    ! have to be in place for biogeophysics

    do fc = 1,num_special_col
       c = special_col(fc)

       this%seed_col(c)    = 0._r8
       this%prod1_col(c)   = 0._r8
       this%prod10_col(c)  = 0._r8	  
       this%prod100_col(c) = 0._r8	  
       this%totprod_col(c) = 0._r8	  
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
  subroutine Restart ( this,  bounds, ncid, flag, cnstate_vars)

    !
    ! !DESCRIPTION: 
    ! Read/write CN restart data for carbon state
    !
    ! !USES:
    use shr_infnan_mod      , only : isnan => shr_infnan_isnan, nan => shr_infnan_nan, assignment(=)
    use clm_time_manager    , only : is_restart, get_nstep
    use clm_varctl          , only : spinup_mortality_factor
    use CNStateType         , only : cnstate_type
    use restUtilMod
    use ncdio_pio
    !
    ! !ARGUMENTS:
    class (phosphorusstate_type) :: this
    type(bounds_type)          , intent(in)    :: bounds 
    type(file_desc_t)          , intent(inout) :: ncid   
    type(cnstate_type)         , intent(in)    :: cnstate_vars
    character(len=*)           , intent(in)    :: flag   !'read' or 'write' or 'define'
    !
    ! !LOCAL VARIABLES:
    integer            :: i,j,k,l,c,a,b,d
    logical            :: readvar
    integer            :: idata
    logical            :: exit_spinup = .false.
    logical            :: enter_spinup = .false.
    real(r8)           :: m, m_veg         ! multiplier for the exit_spinup code
    real(r8), pointer  :: ptr2d(:,:) ! temp. pointers for slicing larger arrays
    real(r8), pointer  :: ptr1d(:)   ! temp. pointers for slicing larger arrays
    character(len=128) :: varname    ! temporary
    integer            :: itemp      ! temporary 
    integer , pointer  :: iptemp(:)  ! pointer to memory to be allocated
    ! spinup state as read from restart file, for determining whether to enter or exit spinup mode.
    integer            :: restart_file_spinup_state 
    ! flags for comparing the model and restart decomposition cascades
    integer            :: decomp_cascade_state, restart_file_decomp_cascade_state 
    real(r8)           :: smax_c, ks_sorption_c
    real(r8)           :: rootfr(1:nlevdecomp)
    real(r8)           :: pinit_prof(1:nlevdecomp)
    real(r8)           :: rootfr_tot

    !------------------------------------------------------------------------

    !--------------------------------
    ! patch phosphorus state variables
    !--------------------------------
    associate(&
         isoilorder     => cnstate_vars%isoilorder  &
         )

    call NutrientStateRestart(this, bounds, ncid, flag)

    call restartvar(ncid=ncid, flag=flag, varname='retransp', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%retransp_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='plant_p_buffer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%plant_p_buffer_patch)

    !--------------------------------
    ! column phosphorus state variables
    !--------------------------------

    ! sminn
    if (use_vertsoilc) then
       ptr2d => this%solutionp_vr_col
       call restartvar(ncid=ncid, flag=flag, varname="solutionp_vr", xtype=ncd_double,  &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)
       ptr2d => this%labilep_vr_col
       call restartvar(ncid=ncid, flag=flag, varname="labilep_vr", xtype=ncd_double,  &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)
       ptr2d => this%secondp_vr_col
       call restartvar(ncid=ncid, flag=flag, varname="secondp_vr", xtype=ncd_double,  &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)
       ptr2d => this%occlp_vr_col
       call restartvar(ncid=ncid, flag=flag, varname="occlp_vr", xtype=ncd_double,  &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)
       ptr2d => this%primp_vr_col
       call restartvar(ncid=ncid, flag=flag, varname="primp_vr", xtype=ncd_double,  &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)

    else

       ptr1d => this%solutionp_vr_col(:,1)
       call restartvar(ncid=ncid, flag=flag, varname="solutionp", xtype=ncd_double,&
            dim1name='column', &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp' , readvar=readvar, data=ptr1d)

       ptr1d => this%labilep_vr_col(:,1)
       call restartvar(ncid=ncid, flag=flag, varname="labilep", xtype=ncd_double,&
            dim1name='column', &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp' , readvar=readvar, data=ptr1d)

       ptr1d => this%secondp_vr_col(:,1)
       call restartvar(ncid=ncid, flag=flag, varname="secondp", xtype=ncd_double,&
            dim1name='column', &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp' , readvar=readvar, data=ptr1d)

       ptr1d => this%occlp_vr_col(:,1)
       call restartvar(ncid=ncid, flag=flag, varname="occlp", xtype=ncd_double,&
            dim1name='column', &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp' , readvar=readvar, data=ptr1d)

       ptr1d => this%primp_vr_col(:,1)
       call restartvar(ncid=ncid, flag=flag, varname="primp", xtype=ncd_double,&
            dim1name='column', &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp' , readvar=readvar, data=ptr1d)

    end if
    if (flag=='read' .and. .not. readvar) then
       call endrun(msg='ERROR::'//trim(varname)//' is required on an initialization dataset'//&
            errMsg(__FILE__, __LINE__))
    end if

    !--------------------------------
    ! Spinup state
    !--------------------------------

    ! Do nothing for write
    ! Note that the call to write spinup_state out was done in CNCarbonStateType and
    ! cannot be called again because it will try to define the variable twice
    ! when the flag below is set to define
    if (flag == 'read') then
       call restartvar(ncid=ncid, flag=flag, varname='spinup_state', xtype=ncd_int,  &
            long_name='Spinup state of the model that wrote this restart file: ' &
            // ' 0 = normal model mode, 1 = AD spinup', units='', &
            interpinic_flag='copy', readvar=readvar,  data=idata)
       if (readvar) then
          restart_file_spinup_state = idata
       else
          ! assume, for sake of backwards compatibility, that if spinup_state is not in 
          ! the restart file then current model state is the same as prior model state
          restart_file_spinup_state = spinup_state
          if ( masterproc ) then
             write(iulog,*) ' WARNING!  Restart file does not contain info ' &
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
          if ( masterproc ) write(iulog,*) ' NitrogenStateType Restart: taking SOM pools out of AD spinup mode'
          exit_spinup = .true.
       else if (spinup_state == 1 .and. restart_file_spinup_state == 0 ) then
          if ( masterproc ) write(iulog,*) ' NitrogenStateType Restart: taking SOM pools into AD spinup mode'
          enter_spinup = .true.
       else
          call endrun(msg=' Error in entering/exiting spinup.  spinup_state ' &
               // ' != restart_file_spinup_state, but do not know what to do'//&
               errMsg(__FILE__, __LINE__))
       end if
       if (get_nstep() >= 2) then
          call endrun(msg=' Error in entering/exiting spinup - should occur only when nstep = 1'//&
               errMsg(__FILE__, __LINE__))
       endif
       do k = 1, ndecomp_pools
          do c = bounds%begc, bounds%endc
            do j = 1, nlevdecomp
	       if ( exit_spinup ) then
		 m = decomp_cascade_con%spinup_factor(k)
                 if (decomp_cascade_con%spinup_factor(k) > 1) m = m  / cnstate_vars%scalaravg_col(c,j)
               else if ( enter_spinup ) then 
                 m = 1. / decomp_cascade_con%spinup_factor(k)
		 if (decomp_cascade_con%spinup_factor(k) > 1) m = m  * cnstate_vars%scalaravg_col(c,j)
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
          if (nu_com == 'RD' .and. exit_spinup) then
              !Initialize plant P storage pool when exiting spinup from CN only mode 
              if (this%pool_patch(i) .lt. this%leaf_patch(i)) then 
                  this%pool_patch(i) = this%leaf_patch(i)
              else if (this%pool_patch(i) .lt. this%leaf_storage_patch(i)) then 
                  this%pool_patch(i) = this%leaf_storage_patch(i)
              end if
          end if
       end do

       ! soil phosphorus initialization when exit AD spinup Qing Z. 2017
       if ( exit_spinup) then ! AD spinup -> RG spinup
          if (.not. cnstate_vars%pdatasets_present) then
              call endrun(msg='ERROR:: P pools are required on surface dataset'//&
              errMsg(__FILE__, __LINE__))
          end if

          do j = 1, nlevdecomp
             rootfr(j) = exp(-3.0* zsoi(j))
          end do
          rootfr_tot = 0._r8
          do j = 1, nlevdecomp
             rootfr_tot = rootfr_tot + rootfr(j)
          end do
          do j = 1, nlevdecomp
             pinit_prof(j) = rootfr(j) / rootfr_tot / dzsoi_decomp(j) ! 1/m
          end do

          do c = bounds%begc, bounds%endc
             if (use_vertsoilc) then
                do j = 1, nlevdecomp
                   ! solve equilibrium between loosely adsorbed and solution
                   ! phosphorus
                   ! the P maps used in the initialization are generated for the top 50cm soils
                   ! assume soil below 50 cm has the same p pool concentration
                   ! divide 0.5m when convert p pools from g/m2 to g/m3
                   ! assume p pools evenly distributed at dif layers
                   ! Prescribe P initial profile based on exponential rooting profile [need to improve]
                   if ((nu_com .eq. 'ECA') .or. (nu_com .eq. 'MIC')) then
                      a = 1.0_r8
                      b = VMAX_MINSURF_P_vr(j,cnstate_vars%isoilorder(c)) + &
                          KM_MINSURF_P_vr(j,cnstate_vars%isoilorder(c)) - cnstate_vars%labp_col(c)*pinit_prof(j)
                      d = -1.0_r8* cnstate_vars%labp_col(c)*pinit_prof(j) * KM_MINSURF_P_vr(j,cnstate_vars%isoilorder(c))

                      this%solutionp_vr_col(c,j) = (-b+(b**2.0_r8-4.0_r8*a*d)**0.5_r8)/(2.0_r8*a)
                      this%labilep_vr_col(c,j) = cnstate_vars%labp_col(c)*pinit_prof(j) - this%solutionp_vr_col(c,j)
                      this%secondp_vr_col(c,j) = cnstate_vars%secp_col(c)*pinit_prof(j)
                      this%occlp_vr_col(c,j) = cnstate_vars%occp_col(c)*pinit_prof(j)
                      this%primp_vr_col(c,j) = cnstate_vars%prip_col(c)*pinit_prof(j)
                   else if (nu_com .eq. 'RD') then
                      a = 1.0_r8
                      b = smax(cnstate_vars%isoilorder(c)) + &
                          ks_sorption(cnstate_vars%isoilorder(c)) - cnstate_vars%labp_col(c)/0.5_r8
                      d = -1.0_r8* cnstate_vars%labp_col(c)/0.5_r8 * ks_sorption(cnstate_vars%isoilorder(c))

                      this%solutionp_vr_col(c,j) = (-b+(b**2.0_r8-4.0_r8*a*d)**0.5_r8)/(2.0_r8*a)
                      this%labilep_vr_col(c,j) = cnstate_vars%labp_col(c)/0.5_r8 - this%solutionp_vr_col(c,j)
                      this%secondp_vr_col(c,j) = cnstate_vars%secp_col(c)/0.5_r8
                      this%occlp_vr_col(c,j) = cnstate_vars%occp_col(c)/0.5_r8
                      this%primp_vr_col(c,j) = cnstate_vars%prip_col(c)/0.5_r8
                   end if

                   if (nu_com .eq. 'RD') then 
                       smax_c = smax(isoilorder(c))
                       ks_sorption_c = ks_sorption(isoilorder(c))
                       this%solutionp_vr_col(c,j) = (cnstate_vars%labp_col(c)/0.5_r8*ks_sorption_c)/&
                                    (smax_c-cnstate_vars%labp_col(c)/0.5_r8)
                       this%labilep_vr_col(c,j) = cnstate_vars%labp_col(c)/0.5_r8
                       this%secondp_vr_col(c,j) = cnstate_vars%secp_col(c)/0.5_r8
                       this%occlp_vr_col(c,j) = cnstate_vars%occp_col(c)/0.5_r8
                       this%primp_vr_col(c,j) = cnstate_vars%prip_col(c)/0.5_r8
                   end if

                end do
             else
                if ((nu_com .eq. 'ECA') .or. (nu_com .eq. 'MIC')) then
                   a = 1.0_r8
                   b = VMAX_MINSURF_P_vr(1,cnstate_vars%isoilorder(c)) + &
                       KM_MINSURF_P_vr(1,cnstate_vars%isoilorder(c)) - cnstate_vars%labp_col(c)/zisoi(nlevdecomp)
                   d = -1.0_r8* cnstate_vars%labp_col(c)/zisoi(nlevdecomp) * KM_MINSURF_P_vr(j,cnstate_vars%isoilorder(c))

                   this%solutionp_vr_col(c,1) = (-b+(b**2.0_r8-4.0_r8*a*d)**0.5_r8)/(2.0_r8*a) * zisoi(nlevdecomp) ! convert to g/m2
                   this%labilep_vr_col(c,1) = cnstate_vars%labp_col(c) - this%solutionp_vr_col(c,1)
                   this%secondp_vr_col(c,1) = cnstate_vars%secp_col(c)
                   this%occlp_vr_col(c,1) = cnstate_vars%occp_col(c)
                   this%primp_vr_col(c,1) = cnstate_vars%prip_col(c)
                else if (nu_com .eq. 'RD') then
                   a = 1.0_r8
                   b = smax(cnstate_vars%isoilorder(c)) + &
                       ks_sorption(cnstate_vars%isoilorder(c)) - cnstate_vars%labp_col(c)/0.5_r8
                   d = -1.0_r8* cnstate_vars%labp_col(c)/0.5_r8 * ks_sorption(cnstate_vars%isoilorder(c))

                   this%solutionp_vr_col(c,1) = (-b+(b**2.0_r8-4.0_r8*a*d)**0.5_r8)/(2.0_r8*a) * 0.5_r8 ! convert to g/m2
                   this%labilep_vr_col(c,1) = cnstate_vars%labp_col(c) - this%solutionp_vr_col(c,1)
                   this%secondp_vr_col(c,1) = cnstate_vars%secp_col(c)
                   this%occlp_vr_col(c,1) = cnstate_vars%occp_col(c)
                   this%primp_vr_col(c,1) = cnstate_vars%prip_col(c)
                end if
             end if
          end do
       end if
       
    end if
    end associate

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine SetValues ( this, &
       num_patch, filter_patch, value_patch, &
       num_column, filter_column, value_column)
    !
    ! !DESCRIPTION:
    ! Set phosphorus state variables
    !
    ! !ARGUMENTS:
    class (phosphorusstate_type) :: this
    integer , intent(in) :: num_patch
    integer , intent(in) :: filter_patch(:)
    real(r8), intent(in) :: value_patch
    integer , intent(in) :: num_column
    integer , intent(in) :: filter_column(:)
    real(r8), intent(in) :: value_column
    !
    ! !LOCAL VARIABLES:
    integer :: fi,i     ! loop index
    integer :: j,k      ! indices
    !------------------------------------------------------------------------

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
       this%retransp_patch(i)           = value_patch
       this%pool_patch(i)              = value_patch
       this%veg_trunc_patch(i)             = value_patch
       this%dispveg_patch(i)           = value_patch
       this%storveg_patch(i)           = value_patch
       this%totveg_patch(i)            = value_patch
       this%totpft_patch(i)            = value_patch
    end do

    if ( crop_prog )then
       do fi = 1,num_patch
          i = filter_patch(fi)
          this%grain_patch(i)            = value_patch
          this%grain_storage_patch(i)    = value_patch
          this%grain_xfer_patch(i)       = value_patch
          this%cropseed_deficit_patch(i) = value_patch
       end do
    end if

    do fi = 1,num_column
       i = filter_column(fi)

       this%sminp_col(i)       = value_column
       this%solutionp_col(i)   = value_column
       this%labilep_col(i)     = value_column
       this%secondp_col(i)     = value_column
       this%occlp_col(i)       = value_column
       this%primp_col(i)       = value_column
       this%veg_trunc_col(i)      = value_column
       this%cwd_col(i)        = value_column
       this%totlit_col(i)     = value_column
       this%totsom_col(i)     = value_column
       this%totecosys_col(i)  = value_column
       this%totcol_col(i)     = value_column
       this%totsom_1m_col(i)  = value_column
       this%totlit_1m_col(i)  = value_column
    end do

    do j = 1,nlevdecomp_full
       do fi = 1,num_column
          i = filter_column(fi)
          this%sminp_vr_col(i,j)       = value_column
          this%solutionp_vr_col(i,j)   = value_column
          this%labilep_vr_col(i,j)     = value_column
          this%secondp_vr_col(i,j)     = value_column
          this%occlp_vr_col(i,j)       = value_column
          this%primp_vr_col(i,j)       = value_column
          this%soil_trunc_vr_col(i,j)      = value_column
       end do
    end do

    ! column and decomp_pools
    do k = 1, ndecomp_pools
       do fi = 1,num_column
          i = filter_column(fi)
          this%decomp_pools_col(i,k)    = value_column
          this%decomp_pools_1m_col(i,k) = value_column
       end do
    end do

    ! column levdecomp, and decomp_pools
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
    class(phosphorusstate_type) :: this
    type(bounds_type), intent(in)  :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer  :: p          ! indices
    !-----------------------------------------------------------------------

    do p = bounds%begp,bounds%endp
       this%dispveg_patch(p) = 0._r8
       this%storveg_patch(p) = 0._r8
       this%totveg_patch(p)  = 0._r8
       this%totpft_patch(p)  = 0._r8
    end do

  end subroutine ZeroDwt

  !-----------------------------------------------------------------------
  subroutine Summary(this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
    !
    ! !USES:
    use clm_varpar    , only: nlevdecomp,ndecomp_cascade_transitions,ndecomp_pools
    use clm_varctl    , only: use_nitrif_denitrif
    use subgridAveMod , only: p2c
    !
    ! !ARGUMENTS:
    class (phosphorusstate_type) :: this
    type(bounds_type) , intent(in) :: bounds  
    integer           , intent(in) :: num_soilc       ! number of soil columns in filter
    integer           , intent(in) :: filter_soilc(:) ! filter for soil columns
    integer           , intent(in) :: num_soilp       ! number of soil patches in filter
    integer           , intent(in) :: filter_soilp(:) ! filter for soil patches
    !
    ! !LOCAL VARIABLES:
    integer  :: c,p,j,k,l   ! indices
    integer  :: fp,fc       ! lake filter indices
    real(r8) :: maxdepth    ! depth to integrate soil variables
    real(r8) :: cropseedp_deficit_col(bounds%begc:bounds%endc)
    !-----------------------------------------------------------------------

    call NutrientStatePatchSummary( this, bounds, num_soilp, filter_soilp)
    call NutrientStateColumnSummary(this, bounds, num_soilc, filter_soilc)

    ! add phosphorus-specific pools to patch-level variables
    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! stored vegetation phosphorus, including retranslocated N pool (STORVEGN)
       this%storveg_patch(p) = &
            this%storveg_patch(p) + &
            this%retransp_patch(p)

       ! total vegetation phosphorus (TOTVEGN)
       this%totveg_patch(p) = &
            this%dispveg_patch(p) + &
            this%storveg_patch(p)

       ! total pft-level carbon (add pft_ntrunc)
       this%totpft_patch(p) = &
            this%totveg_patch(p) + &
            this%veg_trunc_patch(p)

    end do

    ! aggregate patch-level state variables to column-level
    call p2c(bounds, num_soilc, filter_soilc, &
         this%totveg_patch(bounds%begp:bounds%endp), &
         this%totveg_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%totpft_patch(bounds%begp:bounds%endp), &
         this%totpft_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%cropseed_deficit_patch(bounds%begp:bounds%endp), &
         cropseedp_deficit_col(bounds%begc:bounds%endc))

    ! vertically integrate soil mineral P pools
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%solutionp_col(c) = 0._r8
       this%labilep_col(c)   = 0._r8
       this%secondp_col(c)   = 0._r8
       this%occlp_col(c)     = 0._r8
       this%primp_col(c)     = 0._r8
    end do

    do j = 1, nlevdecomp
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%solutionp_col(c) = &
               this%solutionp_col(c) + &
               this%solutionp_vr_col(c,j) * dzsoi_decomp(j)
          this%labilep_col(c) = &
               this%labilep_col(c) + &
               this%labilep_vr_col(c,j) * dzsoi_decomp(j)
          this%secondp_col(c) = &
               this%secondp_col(c) + &
               this%secondp_vr_col(c,j) * dzsoi_decomp(j)
          this%occlp_col(c) = &
               this%occlp_col(c) + &
               this%occlp_vr_col(c,j) * dzsoi_decomp(j)
          this%primp_col(c) = &
               this%primp_col(c) + &
               this%primp_vr_col(c,j) * dzsoi_decomp(j)
       end do
    end do

    ! total sminp
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%sminp_col(c)      = 0._r8
    end do
    do j = 1, nlevdecomp
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%sminp_vr_col(c,j) = this%solutionp_vr_col(c,j)+ &
               this%labilep_vr_col(c,j)+ &
               this%secondp_vr_col(c,j)
       end do
    end do
    do j = 1, nlevdecomp
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%sminp_col(c) =  this%sminp_col(c) + &
               this%sminp_vr_col(c,j) * dzsoi_decomp(j)
       end do
    end do

    ! column-level summary
    do fc = 1,num_soilc
       c = filter_soilc(fc)

       ! total wood product phosphorus
       this%totprod_col(c) = &
            this%prod1_col(c)  + &
            this%prod10_col(c) + &
            this%prod100_col(c)	 

       ! total ecosystem phosphorus, including veg (TOTECOSYSP)
       this%totecosys_col(c) = &
            this%cwd_col(c)       + &
            this%totlit_col(c)    + &
            this%totsom_col(c)    + &
            this%solutionp_col(c) + &
            this%labilep_col(c)   + &
            this%secondp_col(c)   + &
            this%primp_col(c)     + &
            this%occlp_col(c)     + &
            this%totprod_col(c)   + &
            this%totveg_col(c)

       ! total column phosphorus, including pft (TOTCOLP)
       this%totcol_col(c) = &
            this%totpft_col(c)    + &
            this%cwd_col(c)       + &
            this%totlit_col(c)    + &
            this%totsom_col(c)    + &
            this%prod1_col(c)     + &
            this%solutionp_col(c) + &
            this%labilep_col(c)   + &
            this%secondp_col(c)   + &
            this%veg_trunc_col(c) + &
            cropseedp_deficit_col(c)
    end do

 end subroutine Summary

   !-----------------------------------------------------------------------
  subroutine DynamicPatchAdjustments( this, &
       bounds,                              &
       num_filterp_with_inactive,           &
       filterp_with_inactive,               &
       prior_weights,                       &
       patch_state_updater,                 &
       dwt_leafp_seed,                      &
       dwt_deadstemp_seed,                  &
       dwt_ppool_seed,                      &
       conv_pflux,                          &
       dwt_frootp_to_litter,                &
       dwt_livecrootp_to_litter,            &
       dwt_deadcrootp_to_litter,            &
       prod10_pflux,                        &
       prod100_pflux,                       &
       crop_product_pflux                   &
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
    class(phosphorusstate_type)    , intent(inout) :: this
    type(bounds_type)              , intent(in)    :: bounds
    integer                        , intent(in)    :: num_filterp_with_inactive
    integer                        , intent(in)    :: filterp_with_inactive(:)
    type(prior_weights_type)       , intent(in)    :: prior_weights
    type(patch_state_updater_type) , intent(in)    :: patch_state_updater
    real(r8)                       , intent(inout) :: dwt_leafp_seed           (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_deadstemp_seed       (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_ppool_seed           (bounds%begp:)
    real(r8)                       , intent(inout) :: conv_pflux               (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_frootp_to_litter     (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_livecrootp_to_litter (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_deadcrootp_to_litter (bounds%begp:)
    real(r8)                       , intent(inout) :: prod10_pflux             (bounds%begp:)
    real(r8)                       , intent(inout) :: prod100_pflux            (bounds%begp:)
    real(r8)                       , intent(inout) :: crop_product_pflux       (bounds%begp:)
    !
    ! !LOCAL VARIABLES:
    integer                     :: begp, endp
    integer                     :: l, c, p
    logical                     :: old_weight_was_zero      (bounds%begp:bounds%endp)
    logical                     :: patch_grew               (bounds%begp:bounds%endp)

    ! The following are only set for growing patches:
    real(r8)                    :: seed_leafp_patch         (bounds%begp:bounds%endp)
    real(r8)                    :: seed_leafp_storage_patch (bounds%begp:bounds%endp)
    real(r8)                    :: seed_leafp_xfer_patch    (bounds%begp:bounds%endp)
    real(r8)                    :: seed_deadstemp_patch     (bounds%begp:bounds%endp)
    real(r8)                    :: seed_ppool_patch         (bounds%begp:bounds%endp)

    real(r8)                    :: wood_product_pflux       (bounds%begp:bounds%endp)
    real(r8)                    :: deadstemp_patch_temp     (bounds%begp:bounds%endp)

    character(len=*), parameter :: subname = 'PStateDynamicPatchAdjustments'
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    SHR_ASSERT_ALL((ubound(dwt_leafp_seed           ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_deadstemp_seed       ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_ppool_seed           ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(conv_pflux               ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_frootp_to_litter     ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_livecrootp_to_litter ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_deadcrootp_to_litter ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(prod10_pflux             ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(prod100_pflux            ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(crop_product_pflux       ) == (/endp/)), errMsg(__FILE__, __LINE__))

    call NutrientStateDynamicPatchAdjustments( &
       this,                                   &
       bounds,                                 &
       num_filterp_with_inactive,              &
       filterp_with_inactive,                  &
       prior_weights,                          &
       patch_state_updater,                    &
       NUTRIENT_SPECIES_P,                     &
       dwt_leafp_seed,                         &
       dwt_deadstemp_seed,                     &
       conv_pflux,                             &
       dwt_frootp_to_litter,                   &
       dwt_livecrootp_to_litter,               &
       dwt_deadcrootp_to_litter,               &
       prod10_pflux,                           &
       prod100_pflux,                          &
       crop_product_pflux,                     &
       dwt_ppool_seed,                         &
       ppool_seed_param,                       &
       seed_ppool_patch                        &
    )

    ! 1) RETRANSP_PATCH
    call patch_state_updater%update_patch_state(              &
         bounds                                             , &
         num_filterp_with_inactive                          , &
         filterp_with_inactive                              , &
         var               = this%retransp_patch(begp:endp) , &
         flux_out_grc_area = conv_pflux(begp:endp))

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
    use dynPriorWeightsMod       , only : prior_weights_type
    use landunit_varcon          , only : istsoil, istcrop
    use dynColumnStateUpdaterMod , only : column_state_updater_type
    !
    ! !ARGUMENTS:
    class(phosphorusstate_type)     , intent(inout) :: this
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: clump_index
    type(column_state_updater_type) , intent(in)    :: column_state_updater
    !
    ! !LOCAL VARIABLES:
    integer                     :: l, j
    integer                     :: begc, endc
    real(r8)                    :: adjustment_one_level(bounds%begc:bounds%endc)

    character(len=*), parameter :: subname = 'NStateDynamicColumnAdjustments'
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

          this%dyn_bal_adjustments_col(begc:endc) =      &
               this%dyn_bal_adjustments_col(begc:endc) + &
               adjustment_one_level(begc:endc) * dzsoi_decomp(j)

       end do
    end do

    do j = 1, nlevdecomp
       call column_state_updater%update_column_state_no_special_handling( &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = this%soil_trunc_vr_col(begc:endc,j),                &
            adjustment  = adjustment_one_level(begc:endc))

       this%dyn_bal_adjustments_col(begc:endc) =      &
            this%dyn_bal_adjustments_col(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)

       call column_state_updater%update_column_state_no_special_handling( &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = this%solutionp_vr_col(begc:endc,j),             &
            adjustment  = adjustment_one_level(begc:endc))

       this%dyn_bal_adjustments_col(begc:endc) =      &
            this%dyn_bal_adjustments_col(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)

       call column_state_updater%update_column_state_no_special_handling( &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = this%labilep_vr_col(begc:endc,j),               &
            adjustment  = adjustment_one_level(begc:endc))

       this%dyn_bal_adjustments_col(begc:endc) =      &
            this%dyn_bal_adjustments_col(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)

       call column_state_updater%update_column_state_no_special_handling( &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = this%secondp_vr_col(begc:endc,j),               &
            adjustment  = adjustment_one_level(begc:endc))

       this%dyn_bal_adjustments_col(begc:endc) =      &
            this%dyn_bal_adjustments_col(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)
    end do

  end subroutine DynamicColumnAdjustments

end module PhosphorusStateType
