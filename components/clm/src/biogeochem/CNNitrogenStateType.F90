module CNNitrogenStateType

#include "shr_assert.h"

  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_infnan_mod         , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use clm_varpar             , only : ndecomp_cascade_transitions, ndecomp_pools, nlevcan
  use clm_varpar             , only : nlevdecomp_full, nlevdecomp, crop_prog
  use clm_varcon             , only : spval, ispval, dzsoi_decomp, zisoi
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
  use VegetationType         , only : veg_pp
  use clm_varctl             , only : use_pflotran, pf_cmode
  use clm_varctl             , only : nu_com, use_crop
  use dynPatchStateUpdaterMod, only : patch_state_updater_type               
  use CNSpeciesMod           , only : NUTRIENT_SPECIES_N
  use NutrientStateType      , only : nutrientstate_type, NutrientStateInitAllocate
  use NutrientStateType      , only : NutrientStateInitHistory
  use CNSpeciesMod           , only : species_from_string, species_name_from_string
  ! 
  ! !PUBLIC TYPES:
  implicit none
  save
  private

  real(r8) , parameter :: npool_seed_param     = 0.1_r8

  type, public, extends (nutrientstate_type) :: nitrogenstate_type

     real(r8), pointer :: retransn_patch            (:)     ! patch (gN/m2) plant pool of retranslocated N
     real(r8), pointer :: plant_n_buffer_patch      (:)   ! patch (gN/m2) pft-level abstract N storage
     real(r8), pointer :: plant_n_buffer_col        (:)   ! patch (gN/m2) col-level abstract N storage
     real(r8), pointer :: sminn_vr_col              (:,:) ! col (gN/m3) vertically-resolved soil mineral N

     ! NITRIF_DENITRIF
     real(r8), pointer :: smin_no3_vr_col           (:,:) ! col (gN/m3) vertically-resolved soil mineral NO3
     real(r8), pointer :: smin_no3_col              (:)   ! col (gN/m2) soil mineral NO3 pool
     real(r8), pointer :: smin_nh4_vr_col           (:,:) ! col (gN/m3) vertically-resolved soil mineral NH4
     real(r8), pointer :: smin_nh4_col              (:)   ! col (gN/m2) soil mineral NH4 pool

     ! summary (diagnostic) state variables, not involved in mass balance
     real(r8), pointer :: sminn_col                 (:)   ! col (gN/m2) soil mineral N

     ! for newly-added coupled codes with pflotran (it should be included in total 'sminn' defined above when doing summation)
     real(r8), pointer :: smin_nh4sorb_vr_col       (:,:) ! col (gN/m3) vertically-resolved soil mineral NH4 absorbed
     real(r8), pointer :: smin_nh4sorb_col          (:)   ! col (gN/m2) soil mineral NH4 pool absorbed

     real(r8), pointer :: plant_nbuffer_col         (:)   ! col plant nitrogen buffer, (gN/m2), used to exchange info with betr 

     ! for dynamic C/N/P allocation cost-benefit analysis
     real(r8), pointer :: npimbalance_patch         (:)
     real(r8), pointer :: pnup_pfrootc_patch        (:)
     real(r8), pointer :: ppup_pfrootc_patch        (:)
     real(r8), pointer :: ptlai_pleafc_patch        (:)
     
     real(r8), pointer :: ppsnsun_ptlai_patch       (:)
     real(r8), pointer :: ppsnsun_pleafn_patch      (:)
     real(r8), pointer :: ppsnsun_pleafp_patch      (:)
     
     real(r8), pointer :: plmrsun_ptlai_patch       (:)
     real(r8), pointer :: plmrsun_pleafn_patch      (:)
     real(r8), pointer :: plaisun_ptlai_patch       (:)
     
     real(r8), pointer :: ppsnsha_ptlai_patch       (:)
     real(r8), pointer :: ppsnsha_pleafn_patch      (:)
     real(r8), pointer :: ppsnsha_pleafp_patch      (:)
     
     real(r8), pointer :: plmrsha_ptlai_patch       (:)
     real(r8), pointer :: plmrsha_pleafn_patch      (:)
     real(r8), pointer :: plaisha_ptlai_patch       (:)
     
     real(r8), pointer :: benefit_pgpp_pleafc_patch (:)   ! partial gpp / partial leaf carbon (used by symbiotic n2 fixation and dynamic allocation)
     real(r8), pointer :: benefit_pgpp_pleafn_patch (:)   ! partial gpp / partial leaf nitrogen (used by phosphatase activity and dynamic allocation)
     real(r8), pointer :: benefit_pgpp_pleafp_patch (:)   ! partial gpp / partial leaf phosphorus (used by phosphatase activity and dynamic allocation)
     real(r8), pointer :: cost_pgpp_pfrootc_patch   (:)   ! partial gpp /  partial fine root carbon (used by dynamic allocation)
     real(r8), pointer :: cost_plmr_pleafc_patch    (:)   ! partial maintenance respiration /  partial leaf carbon (used by dynamic allocation)
     real(r8), pointer :: cost_plmr_pleafn_patch    (:)   ! partial maintenance respiration /  partial leaf nitrogen (used by dynamic allocation)
     
     real(r8), pointer :: ppsn_ptlai_z              (:,:)
     real(r8), pointer :: ppsn_pleafn_z             (:,:)
     real(r8), pointer :: ppsn_pleafp_z             (:,:)
     
     real(r8), pointer :: ppsn_ptlai_z_vcmax        (:,:)
     real(r8), pointer :: ppsn_pleafn_z_vcmax       (:,:)
     real(r8), pointer :: ppsn_pleafp_z_vcmax       (:,:)
     
     real(r8), pointer :: ppsn_ptlai_z_jmax         (:,:)
     real(r8), pointer :: ppsn_pleafn_z_jmax        (:,:)
     real(r8), pointer :: ppsn_pleafp_z_jmax        (:,:)
     
     real(r8), pointer :: ppsn_ptlai_z_tpu          (:,:)
     real(r8), pointer :: ppsn_pleafn_z_tpu         (:,:)
     real(r8), pointer :: ppsn_pleafp_z_tpu         (:,:)
    
     real(r8), pointer :: plmr_ptlai_z              (:,:)
     real(r8), pointer :: plmr_pleafn_z             (:,:)

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

  end type nitrogenstate_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds,                           &
       leafc_patch, leafc_storage_patch, frootc_patch, frootc_storage_patch, &
       deadstemc_patch, decomp_cpools_vr_col, decomp_cpools_col, decomp_pools_1m_col)

    class(nitrogenstate_type)         :: this
    type(bounds_type) , intent(in)    :: bounds  
    real(r8)          , intent(in)    :: leafc_patch          (bounds%begp:)
    real(r8)          , intent(in)    :: leafc_storage_patch  (bounds%begp:)
    real(r8)          , intent(in)    :: frootc_patch         (bounds%begp:)
    real(r8)          , intent(in)    :: frootc_storage_patch (bounds%begp:)
    real(r8)          , intent(in)    :: deadstemc_patch      (bounds%begp:)
    real(r8)          , intent(in)    :: decomp_cpools_vr_col (bounds%begc:, 1:, 1:)
    real(r8)          , intent(in)    :: decomp_cpools_col    (bounds%begc:, 1:)
    real(r8)          , intent(in)    :: decomp_pools_1m_col (bounds%begc:, 1:)

    this%species = species_from_string('n')
    this%name    = species_name_from_string('n')

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
    class (nitrogenstate_type) :: this
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

    allocate(this%retransn_patch            (begp:endp))                   ;  this%retransn_patch            (:)   = nan
    allocate(this%plant_n_buffer_patch      (begp:endp))                   ;  this%plant_n_buffer_patch      (:)   = nan
    allocate(this%plant_n_buffer_col        (begc:endc))                   ;  this%plant_n_buffer_col        (:)   = nan
    allocate(this%sminn_vr_col              (begc:endc,1:nlevdecomp_full)) ;  this%sminn_vr_col              (:,:) = nan
    allocate(this%smin_no3_vr_col           (begc:endc,1:nlevdecomp_full)) ;  this%smin_no3_vr_col           (:,:) = nan
    allocate(this%smin_nh4_vr_col           (begc:endc,1:nlevdecomp_full)) ;  this%smin_nh4_vr_col           (:,:) = nan
    allocate(this%smin_no3_col              (begc:endc))                   ;  this%smin_no3_col              (:)   = nan
    allocate(this%smin_nh4_col              (begc:endc))                   ;  this%smin_nh4_col              (:)   = nan
    allocate(this%sminn_col                 (begc:endc))                   ;  this%sminn_col                 (:)   = nan

    ! for dynamic C/N/P allocation
    allocate(this%npimbalance_patch         (begp:endp))                   ;  this%npimbalance_patch         (:)   = nan
    allocate(this%pnup_pfrootc_patch        (begp:endp))                   ;  this%pnup_pfrootc_patch        (:)   = nan
    allocate(this%ppup_pfrootc_patch        (begp:endp))                   ;  this%ppup_pfrootc_patch        (:)   = nan
    allocate(this%ptlai_pleafc_patch        (begp:endp))                   ;  this%ptlai_pleafc_patch        (:)   = nan
    allocate(this%ppsnsun_ptlai_patch       (begp:endp))                   ;  this%ppsnsun_ptlai_patch       (:)   = nan
    allocate(this%ppsnsun_pleafn_patch      (begp:endp))                   ;  this%ppsnsun_pleafn_patch      (:)   = nan
    allocate(this%ppsnsun_pleafp_patch      (begp:endp))                   ;  this%ppsnsun_pleafp_patch      (:)   = nan
    allocate(this%plmrsun_ptlai_patch       (begp:endp))                   ;  this%plmrsun_ptlai_patch       (:)   = nan
    allocate(this%plmrsun_pleafn_patch      (begp:endp))                   ;  this%plmrsun_pleafn_patch      (:)   = nan
    allocate(this%plaisun_ptlai_patch       (begp:endp))                   ;  this%plaisun_ptlai_patch       (:)   = nan
    allocate(this%ppsnsha_ptlai_patch       (begp:endp))                   ;  this%ppsnsha_ptlai_patch       (:)   = nan
    allocate(this%ppsnsha_pleafn_patch      (begp:endp))                   ;  this%ppsnsha_pleafn_patch      (:)   = nan
    allocate(this%ppsnsha_pleafp_patch      (begp:endp))                   ;  this%ppsnsha_pleafp_patch      (:)   = nan
    allocate(this%plmrsha_ptlai_patch       (begp:endp))                   ;  this%plmrsha_ptlai_patch       (:)   = nan
    allocate(this%plmrsha_pleafn_patch      (begp:endp))                   ;  this%plmrsha_pleafn_patch      (:)   = nan
    allocate(this%plaisha_ptlai_patch       (begp:endp))                   ;  this%plaisha_ptlai_patch       (:)   = nan
    allocate(this%benefit_pgpp_pleafc_patch (begp:endp))                   ;  this%benefit_pgpp_pleafc_patch (:)   = nan
    allocate(this%benefit_pgpp_pleafn_patch (begp:endp))                   ;  this%benefit_pgpp_pleafn_patch (:)   = nan
    allocate(this%benefit_pgpp_pleafp_patch (begp:endp))                   ;  this%benefit_pgpp_pleafp_patch (:)   = nan
    allocate(this%cost_pgpp_pfrootc_patch   (begp:endp))                   ;  this%cost_pgpp_pfrootc_patch   (:)   = nan
    allocate(this%cost_plmr_pleafc_patch    (begp:endp))                   ;  this%cost_plmr_pleafc_patch    (:)   = nan
    allocate(this%cost_plmr_pleafn_patch    (begp:endp))                   ;  this%cost_plmr_pleafn_patch    (:)   = nan
    allocate(this%ppsn_ptlai_z              (begp:endp,1:nlevcan))         ;  this%ppsn_ptlai_z              (:,:) = nan
    allocate(this%ppsn_pleafn_z             (begp:endp,1:nlevcan))         ;  this%ppsn_pleafn_z             (:,:) = nan
    allocate(this%ppsn_pleafp_z             (begp:endp,1:nlevcan))         ;  this%ppsn_pleafp_z             (:,:) = nan
    allocate(this%ppsn_ptlai_z_vcmax        (begp:endp,1:nlevcan))         ;  this%ppsn_ptlai_z_vcmax        (:,:) = nan
    allocate(this%ppsn_pleafn_z_vcmax       (begp:endp,1:nlevcan))         ;  this%ppsn_pleafn_z_vcmax       (:,:) = nan
    allocate(this%ppsn_pleafp_z_vcmax       (begp:endp,1:nlevcan))         ;  this%ppsn_pleafp_z_vcmax       (:,:) = nan
    allocate(this%ppsn_ptlai_z_jmax         (begp:endp,1:nlevcan))         ;  this%ppsn_ptlai_z_jmax         (:,:) = nan
    allocate(this%ppsn_pleafn_z_jmax        (begp:endp,1:nlevcan))         ;  this%ppsn_pleafn_z_jmax        (:,:) = nan
    allocate(this%ppsn_pleafp_z_jmax        (begp:endp,1:nlevcan))         ;  this%ppsn_pleafp_z_jmax        (:,:) = nan
    allocate(this%ppsn_ptlai_z_tpu          (begp:endp,1:nlevcan))         ;  this%ppsn_ptlai_z_tpu          (:,:) = nan
    allocate(this%ppsn_pleafn_z_tpu         (begp:endp,1:nlevcan))         ;  this%ppsn_pleafn_z_tpu         (:,:) = nan
    allocate(this%ppsn_pleafp_z_tpu         (begp:endp,1:nlevcan))         ;  this%ppsn_pleafp_z_tpu         (:,:) = nan
    allocate(this%plmr_ptlai_z              (begp:endp,1:nlevcan))         ;  this%plmr_ptlai_z              (:,:) = nan
    allocate(this%plmr_pleafn_z             (begp:endp,1:nlevcan))         ;  this%plmr_pleafn_z             (:,:) = nan

    allocate(this%smin_nh4sorb_vr_col       (begc:endc,1:nlevdecomp_full)) ;  this%smin_nh4sorb_vr_col       (:,:) = nan
    allocate(this%smin_nh4sorb_col          (begc:endc))                   ;  this%smin_nh4sorb_col          (:)   = nan

    allocate(this%plant_nbuffer_col         (begc:endc))                   ;  this%plant_nbuffer_col         (:)   = nan

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
    class(nitrogenstate_type) :: this
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

    !-------------------------------
    ! N state variables - native to PFT
    !-------------------------------
    
    this%retransn_patch(begp:endp) = spval
    call hist_addfld1d (fname='N_RETRANS', units='gN/m^2', &
         avgflag='A', long_name='plant pool of retranslocated N', &
         ptr_patch=this%retransn_patch)

    this%npimbalance_patch(begp:endp) = spval
    call hist_addfld1d (fname='LEAF_NPIMBALANCE', units='gN/gP', &
         avgflag='A', long_name='leaf np imbalance partial C partial P/partial C partial N', &
         ptr_patch=this%npimbalance_patch)
     
    if ( nlevdecomp_full > 1 ) then

       this%sminn_col(begc:endc) = spval
       call hist_addfld1d (fname='N_SMIN', units='gN/m^2', &
            avgflag='A', long_name='soil mineral N', &
            ptr_col=this%sminn_col)

    endif

    this%plant_n_buffer_patch(begp:endp) = spval
    call hist_addfld1d (fname='N_PLANT_BUFFER', units='gN/m^2', &
            avgflag='A', long_name='plant nitrogen stored as buffer', &
            ptr_col=this%plant_n_buffer_patch,default='inactive')
    
    ! add suffix if number of soil decomposition depths is greater than 1
    if (nlevdecomp > 1) then
       vr_suffix = "_vr"
    else 
       vr_suffix = ""
    endif

    if (use_nitrif_denitrif .or. (use_pflotran .and. pf_cmode)) then
       this%smin_no3_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SMIN_NO3'//trim(vr_suffix), units='gN/m^3',  type2d='levdcmp', &
            avgflag='A', long_name='soil mineral NO3 (vert. res.)', &
            ptr_col=this%smin_no3_vr_col)

       this%smin_nh4_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SMIN_NH4'//trim(vr_suffix), units='gN/m^3',  type2d='levdcmp', &
            avgflag='A', long_name='soil mineral NH4 (vert. res.)', &
            ptr_col=this%smin_nh4_vr_col)

       ! pflotran
       if(use_pflotran .and. pf_cmode) then
          this%smin_nh4sorb_vr_col(begc:endc,:) = spval
          call hist_addfld_decomp (fname='SMIN_NH4SORB'//trim(vr_suffix), units='gN/m^3',  type2d='levdcmp', &
            avgflag='A', long_name='soil mineral NH4 absorbed (vert. res.)', &
            ptr_col=this%smin_nh4sorb_vr_col)
       end if

       if ( nlevdecomp_full > 1 ) then
          this%smin_no3_col(begc:endc) = spval
          call hist_addfld1d (fname='SMIN_NO3', units='gN/m^2', &
               avgflag='A', long_name='soil mineral NO3', &
               ptr_col=this%smin_no3_col)

          this%smin_nh4_col(begc:endc) = spval
          call hist_addfld1d (fname='SMIN_NH4', units='gN/m^2', &
               avgflag='A', long_name='soil mineral NH4', &
               ptr_col=this%smin_nh4_col)

          ! pflotran
          if(use_pflotran .and. pf_cmode) then
            this%smin_nh4sorb_col(begc:endc) = spval
            call hist_addfld1d (fname='SMIN_NH4SORB', units='gN/m^2', &
               avgflag='A', long_name='soil mineral NH4 absorbed', &
               ptr_col=this%smin_nh4sorb_col)
          end if
       end if

       this%sminn_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='N_SMIN'//trim(vr_suffix), units='gN/m^3',  type2d='levdcmp', &
            avgflag='A', long_name='soil mineral N', &
            ptr_col=this%sminn_vr_col, default = 'inactive')
    else
       this%sminn_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='N_SMIN'//trim(vr_suffix), units='gN/m^3',  type2d='levdcmp', &
            avgflag='A', long_name='soil mineral N', &
            ptr_col=this%sminn_vr_col)
    end if

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, &
       leafc_patch, leafc_storage_patch, frootc_patch, frootc_storage_patch, &
       deadstemc_patch, decomp_cpools_vr_col, decomp_cpools_col, decomp_pools_1m_col)
    !
    ! !DESCRIPTION:
    ! Initializes time varying variables used only in coupled carbon-nitrogen mode (CN):
    !
    ! !USES:
    use clm_varpar     , only : crop_prog
    use decompMod      , only : bounds_type
    use pftvarcon      , only : noveg, npcropmin
    !
    ! !ARGUMENTS:
    class(nitrogenstate_type)      :: this
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
             this%leaf_patch(p)         = leafc_patch(p)         / veg_vp%leafcn(veg_pp%itype(p))
             this%leaf_storage_patch(p) = leafc_storage_patch(p) / veg_vp%leafcn(veg_pp%itype(p))
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
             this%deadstem_patch(p) = deadstemc_patch(p) / veg_vp%deadwdcn(veg_pp%itype(p))
          else
             this%deadstem_patch(p) = 0._r8
          end if
          
          if (nu_com .ne. 'RD') then
              ! ECA competition calculate root NP uptake as a function of fine root biomass
              ! better to initialize root CNP pools with a non-zero value
              if (veg_pp%itype(p) .ne. noveg) then
                 this%froot_patch(p) = frootc_patch(p) / veg_vp%frootcn(veg_pp%itype(p))
                 this%froot_storage_patch(p) = frootc_storage_patch(p) / veg_vp%frootcn(veg_pp%itype(p))
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
          this%retransn_patch(p)           = 0._r8
          this%pool_patch(p)              = 0._r8
          if (nstor(veg_pp%itype(p)) .gt. 1e-6_r8) then 
              this%pool_patch(p)          = 10.0_r8
          end if
          this%veg_trunc_patch(p)             = 0._r8
          this%dispveg_patch(p)           = 0._r8
          this%storveg_patch(p)           = 0._r8
          this%totveg_patch(p)            = 0._r8
          this%totpft_patch(p)            = 0._r8          
          this%plant_n_buffer_patch(p)     = 1._r8
       end if

       this%npimbalance_patch(p) = 0.0_r8
       this%pnup_pfrootc_patch(p) = 0.0_r8 
       this%benefit_pgpp_pleafc_patch(p) = 0.0_r8   
    end do

    !-------------------------------------------
    ! initialize column-level variables
    !-------------------------------------------

    do c = bounds%begc, bounds%endc
       l = col_pp%landunit(c)
       if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then

          ! column nitrogen state variables
          this%veg_trunc_col(c) = 0._r8
          this%sminn_col(c) = 0._r8
          do j = 1, nlevdecomp
             do k = 1, ndecomp_pools
                this%decomp_pools_vr_col(c,j,k) = decomp_cpools_vr_col(c,j,k) / decomp_cascade_con%initial_cn_ratio(k)
             end do
             this%sminn_vr_col(c,j) = 0._r8
             this%soil_trunc_vr_col(c,j) = 0._r8
          end do
          if ( nlevdecomp > 1 ) then
             do j = nlevdecomp+1, nlevdecomp_full
                do k = 1, ndecomp_pools
                   this%decomp_pools_vr_col(c,j,k) = 0._r8
                end do
                this%sminn_vr_col(c,j) = 0._r8
                this%soil_trunc_vr_col(c,j) = 0._r8
             end do
          end if
          do k = 1, ndecomp_pools
             this%decomp_pools_col(c,k)    = decomp_cpools_col(c,k)    / decomp_cascade_con%initial_cn_ratio(k)
             this%decomp_pools_1m_col(c,k) = decomp_pools_1m_col(c,k) / decomp_cascade_con%initial_cn_ratio(k)
          end do
          if (use_nitrif_denitrif .or. (use_pflotran .and. pf_cmode)) then
             do j = 1, nlevdecomp_full
                this%smin_nh4_vr_col(c,j) = 0._r8
                this%smin_no3_vr_col(c,j) = 0._r8
                if(use_pflotran .and. pf_cmode) then
                    this%smin_nh4sorb_vr_col(c,j) = 0._r8
                end if
             end do
             this%smin_nh4_col(c) = 0._r8
             this%smin_no3_col(c) = 0._r8
             if(use_pflotran .and. pf_cmode) then
                this%smin_nh4sorb_col(c) = 0._r8
             end if
          end if
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

  subroutine Restart ( this,  bounds, ncid, flag, cnstate_vars )
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
    class (nitrogenstate_type) :: this
    type(bounds_type)          , intent(in)    :: bounds 
    type(file_desc_t)          , intent(inout) :: ncid   
    type(cnstate_type)         , intent(in)    :: cnstate_vars
    character(len=*)           , intent(in)    :: flag   !'read' or 'write' or 'define'
    !
    ! !LOCAL VARIABLES:
    integer            :: i,j,k,l,c
    logical            :: readvar
    integer            :: idata
    logical            :: exit_spinup = .false.
    logical            :: enter_spinup = .false.
    real(r8)           :: m, m_veg          ! multiplier for the exit_spinup code
    real(r8), pointer  :: ptr2d(:,:) ! temp. pointers for slicing larger arrays
    real(r8), pointer  :: ptr1d(:)   ! temp. pointers for slicing larger arrays
    character(len=128) :: varname    ! temporary
    integer            :: itemp      ! temporary 
    integer , pointer  :: iptemp(:)  ! pointer to memory to be allocated
    ! spinup state as read from restart file, for determining whether to enter or exit spinup mode.
    integer            :: restart_file_spinup_state 
    ! flags for comparing the model and restart decomposition cascades
    integer            :: decomp_cascade_state, restart_file_decomp_cascade_state 
    !------------------------------------------------------------------------

    !--------------------------------
    ! patch nitrogen state variables
    !--------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='leafn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%leaf_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='leafn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%leaf_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='leafn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%leaf_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='frootn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%froot_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='frootn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%froot_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='frootn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%froot_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livestemn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livestem_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livestemn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livestem_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livestemn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livestem_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadstemn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadstem_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadstemn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadstem_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadstemn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadstem_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livecrootn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livecroot_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livecrootn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livecroot_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livecrootn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livecroot_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadcrootn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadcroot_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadcrootn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadcroot_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadcrootn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadcroot_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='retransn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%retransn_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='npool', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%pool_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='pft_ntrunc', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%veg_trunc_patch) 

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag,  varname='grainn', xtype=ncd_double,  &
            dim1name='pft',    long_name='grain N', units='gN/m2', &
            interpinic_flag='interp', readvar=readvar, data=this%grain_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='grainn_storage', xtype=ncd_double,  &
            dim1name='pft',    long_name='grain N storage', units='gN/m2', &
            interpinic_flag='interp', readvar=readvar, data=this%grain_storage_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='grainn_xfer', xtype=ncd_double,  &
            dim1name='pft',    long_name='grain N transfer', units='gN/m2', &
            interpinic_flag='interp', readvar=readvar, data=this%grain_xfer_patch)
    end if
    
    call restartvar(ncid=ncid, flag=flag,  varname='npimbalance_patch', xtype=ncd_double,  &
        dim1name='pft',    long_name='npimbalance_patch', units='-', &
        interpinic_flag='interp', readvar=readvar, data=this%npimbalance_patch)
     call restartvar(ncid=ncid, flag=flag,  varname='pnup_pfrootc_patch', xtype=ncd_double,  &
        dim1name='pft',    long_name='pnup_pfrootc_patch', units='-', &
        interpinic_flag='interp', readvar=readvar, data=this%pnup_pfrootc_patch)
    call restartvar(ncid=ncid, flag=flag,  varname='benefit_pgpp_pleafc_patch', xtype=ncd_double,  &
        dim1name='pft',    long_name='benefit_pgpp_pleafc_patch', units='-', &
        interpinic_flag='interp', readvar=readvar, data=this%benefit_pgpp_pleafc_patch)
 
    !--------------------------------
    ! column nitrogen state variables
    !--------------------------------

    ! sminn
    if (use_vertsoilc) then
       ptr2d => this%sminn_vr_col
       call restartvar(ncid=ncid, flag=flag, varname="sminn_vr", xtype=ncd_double,  &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)
    else
       ptr1d => this%sminn_vr_col(:,1)
       call restartvar(ncid=ncid, flag=flag, varname="sminn", xtype=ncd_double,  &
            dim1name='column', &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp' , readvar=readvar, data=ptr1d)
    end if
    if (flag=='read' .and. .not. readvar) then
       call endrun(msg='ERROR::'//trim(varname)//' is required on an initialization dataset'//&
            errMsg(__FILE__, __LINE__))
    end if

    ! decomposing N pools
    do k = 1, ndecomp_pools
       varname=trim(decomp_cascade_con%decomp_pool_name_restart(k))//'n'
       if (use_vertsoilc) then
          ptr2d => this%decomp_pools_vr_col(:,:,k)
          call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vr", xtype=ncd_double, &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=ptr2d) 
       else
          ptr1d => this%decomp_pools_vr_col(:,1,k)
          call restartvar(ncid=ncid, flag=flag, varname=varname, xtype=ncd_double,  &
               dim1name='column', &
               long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp' , readvar=readvar, data=ptr1d)
       end if
       if (flag=='read' .and. .not. readvar) then
          call endrun(msg='ERROR:: '//trim(varname)//' is required on an initialization dataset'//&
               errMsg(__FILE__, __LINE__))
       end if
    end do

    if (use_vertsoilc) then
       ptr2d => this%soil_trunc_vr_col
       call restartvar(ncid=ncid, flag=flag, varname="col_ntrunc_vr", xtype=ncd_double,  &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)
    else
       ptr1d => this%soil_trunc_vr_col(:,1)
       call restartvar(ncid=ncid, flag=flag, varname="col_ntrunc", xtype=ncd_double,  &
            dim1name='column', &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp' , readvar=readvar, data=ptr1d)
    end if

    if (use_nitrif_denitrif .or. (use_pflotran .and. pf_cmode)) then
       ! smin_no3_vr
       if (use_vertsoilc) then
          ptr2d => this%smin_no3_vr_col(:,:)
          call restartvar(ncid=ncid, flag=flag, varname='smin_no3_vr', xtype=ncd_double, &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
       else
          ptr1d => this%smin_no3_vr_col(:,1)
          call restartvar(ncid=ncid, flag=flag, varname='smin_no3', xtype=ncd_double, &
               dim1name='column', &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=ptr1d)
       end if
       if (flag=='read' .and. .not. readvar) then
          call endrun(msg= 'ERROR:: smin_no3_vr'//' is required on an initialization dataset' )
       end if
    end if

    if (use_nitrif_denitrif .or. (use_pflotran .and. pf_cmode)) then
       ! smin_nh4
       if (use_vertsoilc) then
          ptr2d => this%smin_nh4_vr_col(:,:)
          call restartvar(ncid=ncid, flag=flag, varname='smin_nh4_vr', xtype=ncd_double, &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=ptr2d) 
       else
          ptr1d => this%smin_nh4_vr_col(:,1)
          call restartvar(ncid=ncid, flag=flag, varname='smin_nh4', xtype=ncd_double, &
               dim1name='column', &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=ptr1d)
       end if
       if (flag=='read' .and. .not. readvar) then
          call endrun(msg= 'ERROR:: smin_nh4_vr'//' is required on an initialization dataset' )
       end if
    end if

    ! pflotran: smin_nh4sorb
    if (use_pflotran .and. pf_cmode) then
       if (use_vertsoilc) then
          ptr2d => this%smin_nh4sorb_vr_col(:,:)
          call restartvar(ncid=ncid, flag=flag, varname='smin_nh4sorb_vr', xtype=ncd_double, &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
        else
          ptr1d => this%smin_nh4sorb_vr_col(:,1)
          call restartvar(ncid=ncid, flag=flag, varname='smin_nh4sorb', xtype=ncd_double, &
               dim1name='column', &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=ptr1d)
       end if
       if (flag=='read' .and. .not. readvar) then
          call endrun(msg= 'ERROR:: smin_nh4sorb_vr'//' is required on an initialization dataset' )
       end if
    end if

    ! Set the integrated sminn based on sminn_vr, as is done in CNSummaryMod (this may
    ! not be the most appropriate method or place to do this)

    this%sminn_col(bounds%begc:bounds%endc) = 0._r8
    do j = 1, nlevdecomp
       do c = bounds%begc, bounds%endc
          this%sminn_col(c) = &
               this%sminn_col(c) + &
               this%sminn_vr_col(c,j) * dzsoi_decomp(j)
       end do
    end do

    call restartvar(ncid=ncid, flag=flag, varname='totcoln', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%totcol_col) 

    call restartvar(ncid=ncid, flag=flag, varname='seedn', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%seed_col) 

    call restartvar(ncid=ncid, flag=flag, varname='prod10n', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%prod10_col) 

    call restartvar(ncid=ncid, flag=flag, varname='prod100n', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%prod100_col) 

    call restartvar(ncid=ncid, flag=flag, varname='prod1n', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%prod1_col)

    ! decomp_cascade_state - the purpose of this is to check to make sure the bgc used 
    ! matches what the restart file was generated with.  
    ! add info about the SOM decomposition cascade

    if (use_century_decomp) then
       decomp_cascade_state = 1
    else
       decomp_cascade_state = 0
    end if
    ! add info about the nitrification / denitrification state
    if (use_nitrif_denitrif .or. (use_pflotran .and. pf_cmode)) then
       decomp_cascade_state = decomp_cascade_state + 10
    end if
    if (flag == 'write') itemp = decomp_cascade_state    
    call restartvar(ncid=ncid, flag=flag, varname='decomp_cascade_state', xtype=ncd_int,  &
         long_name='BGC of the model that wrote this restart file:' &
         // '  1s column: 0 = CLM-CN cascade, 1 = Century cascade;' &
         // ' 10s column: 0 = CLM-CN denitrification, 10 = Century denitrification', units='', &
         interpinic_flag='skip', readvar=readvar, data=itemp)
    if (flag=='read') then
       if (.not. readvar) then
          ! assume, for sake of backwards compatibility, that if decomp_cascade_state 
          ! is not in the restart file, then the current model state is the same as 
          ! the prior model state
          restart_file_decomp_cascade_state = decomp_cascade_state
          if ( masterproc ) write(iulog,*) ' CNRest: WARNING!  Restart file does not ' &
               // ' contain info on decomp_cascade_state used to generate the restart file.  '
          if ( masterproc ) write(iulog,*) '   Assuming the same as current setting: ', decomp_cascade_state
       else
          restart_file_decomp_cascade_state = itemp  
          if (decomp_cascade_state /= restart_file_decomp_cascade_state ) then
             if ( masterproc ) then
                write(iulog,*) 'CNRest: ERROR--the decomposition cascade differs between the current ' &
                     // ' model state and the model that wrote the restart file. '
                write(iulog,*) 'The model will be horribly out of equilibrium until after a lengthy spinup. '
                write(iulog,*) 'Stopping here since this is probably an error in configuring the run. '
                write(iulog,*) 'If you really wish to proceed, then override by setting '
                write(iulog,*) 'override_bgc_restart_mismatch_dump to .true. in the namelist'
                if ( .not. override_bgc_restart_mismatch_dump ) then
                   call endrun(msg= ' CNRest: Stopping. Decomposition cascade mismatch error.'//&
                        errMsg(__FILE__, __LINE__))
                endif
             endif
          endif
       end if
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

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine SetValues ( this, &
       num_patch, filter_patch, value_patch, &
       num_column, filter_column, value_column)
    !
    ! !DESCRIPTION:
    ! Set nitrogen state variables
    !
    ! !ARGUMENTS:
    class (nitrogenstate_type) :: this
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
       this%retransn_patch(i)           = value_patch
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

       this%sminn_col(i)       = value_column
       this%veg_trunc_col(i)      = value_column
       this%cwd_col(i)        = value_column
       if (use_nitrif_denitrif .or. (use_pflotran .and. pf_cmode)) then
          this%smin_no3_col(i) = value_column
          this%smin_nh4_col(i) = value_column
          if(use_pflotran .and. pf_cmode) then
             this%smin_nh4sorb_col(i) = value_column
          end if
       end if
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
          this%sminn_vr_col(i,j)       = value_column
          this%soil_trunc_vr_col(i,j)      = value_column
          if (use_nitrif_denitrif  .or. (use_pflotran .and. pf_cmode)) then
             this%smin_no3_vr_col(i,j) = value_column
             this%smin_nh4_vr_col(i,j) = value_column
             if(use_pflotran .and. pf_cmode) then
               this%smin_nh4sorb_vr_col(i,j) = value_column
             end if
          end if
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
    class(nitrogenstate_type) :: this
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
    use clm_varpar    , only: nlevdecomp_full
    !
    ! !ARGUMENTS:
    class (nitrogenstate_type) :: this
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
    integer  :: nlev
    real(r8) :: cropseedn_deficit_col(bounds%begc:bounds%endc)
    !-----------------------------------------------------------------------

    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! displayed vegetation nitrogen, excluding storage (DISPVEGN)
       this%dispveg_patch(p) = &
            this%leaf_patch(p)      + &
            this%froot_patch(p)     + &
            this%livestem_patch(p)  + &
            this%deadstem_patch(p)  + &
            this%livecroot_patch(p) + &
            this%deadcroot_patch(p)
       
      ! stored vegetation nitrogen, including retranslocated N pool (STORVEGN)
      this%storveg_patch(p) = &
           this%leaf_storage_patch(p)      + &
           this%froot_storage_patch(p)     + &
           this%livestem_storage_patch(p)  + &
           this%deadstem_storage_patch(p)  + &
           this%livecroot_storage_patch(p) + &
           this%deadcroot_storage_patch(p) + &
           this%leaf_xfer_patch(p)         + &
           this%froot_xfer_patch(p)        + &
           this%livestem_xfer_patch(p)     + &
           this%deadstem_xfer_patch(p)     + &
           this%livecroot_xfer_patch(p)    + &
           this%deadcroot_xfer_patch(p)    + &
           this%pool_patch(p)              + &
           this%retransn_patch(p)

      if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
         this%dispveg_patch(p) = &
              this%dispveg_patch(p) + &
              this%grain_patch(p)

         this%storveg_patch(p) = &
              this%storveg_patch(p) + &
              this%grain_storage_patch(p)     + &
              this%grain_xfer_patch(p)
      end if

      ! total vegetation nitrogen (TOTVEGN)
      this%totveg_patch(p) = &
           this%dispveg_patch(p) + &
           this%storveg_patch(p)

      ! total pft-level carbon (add pft_ntrunc)
      this%totpft_patch(p) = &
           this%totveg_patch(p) + &
           this%veg_trunc_patch(p)

   end do

   call p2c(bounds, num_soilc, filter_soilc, &
        this%plant_n_buffer_patch(bounds%begp:bounds%endp), &
        this%plant_n_buffer_col(bounds%begc:bounds%endc))

   call p2c(bounds, num_soilc, filter_soilc, &
        this%totveg_patch(bounds%begp:bounds%endp), &
        this%totveg_col(bounds%begc:bounds%endc))

   call p2c(bounds, num_soilc, filter_soilc, &
        this%totpft_patch(bounds%begp:bounds%endp), &
        this%totpft_col(bounds%begc:bounds%endc))

   call p2c(bounds, num_soilc, filter_soilc, &
        this%cropseed_deficit_patch(bounds%begp:bounds%endp), &
        cropseedn_deficit_col(bounds%begc:bounds%endc))

   ! vertically integrate NO3 NH4 N2O pools
   nlev = nlevdecomp
   if (use_pflotran .and. pf_cmode) nlev = nlevdecomp_full

   if (use_nitrif_denitrif .or. (use_pflotran .and. pf_cmode)) then
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%smin_no3_col(c) = 0._r8
         this%smin_nh4_col(c) = 0._r8
         if(use_pflotran .and. pf_cmode) then
            this%smin_nh4sorb_col(c) = 0._r8
         end if
      end do
      do j = 1, nlev
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            this%smin_no3_col(c) = &
                 this%smin_no3_col(c) + &
                 this%smin_no3_vr_col(c,j) * dzsoi_decomp(j)
            
            this%smin_nh4_col(c) = &
                 this%smin_nh4_col(c) + &
                 this%smin_nh4_vr_col(c,j) * dzsoi_decomp(j)
            if(use_pflotran .and. pf_cmode) then
               this%smin_nh4sorb_col(c) = &
                 this%smin_nh4sorb_col(c) + &
                 this%smin_nh4sorb_vr_col(c,j) * dzsoi_decomp(j)
            end if
          end do 
       end do

    end if

   ! vertically integrate each of the decomposing N pools
   do l = 1, ndecomp_pools
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%decomp_pools_col(c,l) = 0._r8
      end do
      do j = 1, nlev
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            this%decomp_pools_col(c,l) = &
                 this%decomp_pools_col(c,l) + &
                 this%decomp_pools_vr_col(c,j,l) * dzsoi_decomp(j)
         end do
      end do
   end do

   ! for vertically-resolved soil biogeochemistry, calculate some diagnostics of carbon pools to a given depth
   if ( nlevdecomp > 1) then

      do l = 1, ndecomp_pools
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            this%decomp_pools_1m_col(c,l) = 0._r8
         end do
      end do

      ! vertically integrate each of the decomposing n pools to 1 meter
      maxdepth = 1._r8
      do l = 1, ndecomp_pools
         do j = 1, nlevdecomp
            if ( zisoi(j) <= maxdepth ) then
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  this%decomp_pools_1m_col(c,l) = &
                       this%decomp_pools_1m_col(c,l) + &
                       this%decomp_pools_vr_col(c,j,l) * dzsoi_decomp(j)
               end do
            elseif ( zisoi(j-1) < maxdepth ) then
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  this%decomp_pools_1m_col(c,l) = &
                       this%decomp_pools_1m_col(c,l) + &
                       this%decomp_pools_vr_col(c,j,l) * (maxdepth - zisoi(j-1))
               end do
            endif
         end do
      end do
      
      ! total litter nitrogen to 1 meter (TOTLITN_1m)
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%totlit_1m_col(c) = 0._r8
      end do
      do l = 1, ndecomp_pools
         if ( decomp_cascade_con%is_litter(l) ) then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               this%totlit_1m_col(c) = &
                    this%totlit_1m_col(c) + &
                    this%decomp_pools_1m_col(c,l)
            end do
         end if
      end do
      
      ! total soil organic matter nitrogen to 1 meter (TOTSOMN_1m)
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%totsom_1m_col(c) = 0._r8
      end do
      do l = 1, ndecomp_pools
         if ( decomp_cascade_con%is_soil(l) ) then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               this%totsom_1m_col(c) = &
                    this%totsom_1m_col(c) + &
                    this%decomp_pools_1m_col(c,l)
            end do
         end if
      end do
      
   endif
   
   ! total litter nitrogen (TOTLITN)
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      this%totlit_col(c)    = 0._r8
   end do
   do l = 1, ndecomp_pools
      if ( decomp_cascade_con%is_litter(l) ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            this%totlit_col(c) = &
                 this%totlit_col(c) + &
                 this%decomp_pools_col(c,l)
         end do
      end if
   end do
   
   ! total soil organic matter nitrogen (TOTSOMN)
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      this%totsom_col(c)    = 0._r8
   end do
   do l = 1, ndecomp_pools
      if ( decomp_cascade_con%is_soil(l) ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            this%totsom_col(c) = &
                 this%totsom_col(c) + &
                 this%decomp_pools_col(c,l)
         end do
      end if
   end do
   
   ! total cwdn
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      this%cwd_col(c) = 0._r8
   end do
   do l = 1, ndecomp_pools
      if ( decomp_cascade_con%is_cwd(l) ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            this%cwd_col(c) = &
                 this%cwd_col(c) + &
                 this%decomp_pools_col(c,l)
         end do
      end if
   end do

   ! total sminn
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      this%sminn_col(c)      = 0._r8
   end do
   do j = 1, nlev
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%sminn_col(c) = &
              this%sminn_col(c) + &
              this%sminn_vr_col(c,j) * dzsoi_decomp(j)
      end do
   end do

   ! total col_ntrunc
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      this%veg_trunc_col(c) = 0._r8
   end do
   do j = 1, nlev
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%veg_trunc_col(c) = &
              this%veg_trunc_col(c) + &
              this%soil_trunc_vr_col(c,j) * dzsoi_decomp(j)
      end do
   end do

   do fc = 1,num_soilc
      c = filter_soilc(fc)

      ! total wood product nitrogen
      this%totprod_col(c) = &
           this%prod1_col(c) + &
           this%prod10_col(c) + &
           this%prod100_col(c)	 

      ! total ecosystem nitrogen, including veg (TOTECOSYSN)
      this%totecosys_col(c) = &
           this%cwd_col(c) + &
           this%totlit_col(c) + &
           this%totsom_col(c) + &
           this%sminn_col(c) + &
           this%totprod_col(c) + &
           this%totveg_col(c)

      ! total column nitrogen, including pft (TOTCOLN)
      this%totcol_col(c) = &
           this%totpft_col(c) + &
           this%cwd_col(c) + &
           this%totlit_col(c) + &
           this%totsom_col(c) + &
           this%sminn_col(c) + &
           this%prod1_col(c) + &
           this%veg_trunc_col(c)+ &
           this%plant_n_buffer_col(c) + &
           cropseedn_deficit_col(c)
           
      this%totabg_col (c) =  &
           this%totpft_col(c) + &
           this%totprod_col(c) + &
           this%seed_col(c) + &
           this%veg_trunc_col(c)+ &
           this%plant_n_buffer_col(c) 

      this%totblg_col(c) = &
           this%cwd_col(c) + &
           this%totlit_col(c) + &
           this%totsom_col(c) + &
           this%sminn_col(c) 
           
   end do

 end subroutine Summary
 
  !-----------------------------------------------------------------------
  subroutine DynamicPatchAdjustments( this, &
       bounds,                              &
       num_filterp_with_inactive,           &
       filterp_with_inactive,               &
       prior_weights,                       &
       patch_state_updater,                 &
       dwt_leafn_seed,                      &
       dwt_deadstemn_seed,                  &
       dwt_npool_seed,                      &
       conv_nflux,                          &
       dwt_frootn_to_litter,                &
       dwt_livecrootn_to_litter,            &
       dwt_deadcrootn_to_litter,            &
       prod10_nflux,                        &
       prod100_nflux,                       &
       crop_product_nflux                   &
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
    class(nitrogenstate_type)      , intent(inout) :: this
    type(bounds_type)              , intent(in)    :: bounds
    integer                        , intent(in)    :: num_filterp_with_inactive
    integer                        , intent(in)    :: filterp_with_inactive(:)
    type(prior_weights_type)       , intent(in)    :: prior_weights
    type(patch_state_updater_type) , intent(in)    :: patch_state_updater
    real(r8)                       , intent(inout) :: dwt_leafn_seed           (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_deadstemn_seed       (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_npool_seed           (bounds%begp:)
    real(r8)                       , intent(inout) :: conv_nflux               (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_frootn_to_litter     (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_livecrootn_to_litter (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_deadcrootn_to_litter (bounds%begp:)
    real(r8)                       , intent(inout) :: prod10_nflux             (bounds%begp:)
    real(r8)                       , intent(inout) :: prod100_nflux            (bounds%begp:)
    real(r8)                       , intent(inout) :: crop_product_nflux       (bounds%begp:)
    !
    ! !LOCAL VARIABLES:
    integer                     :: begp, endp
    integer                     :: l, c, p
    logical                     :: old_weight_was_zero      (bounds%begp:bounds%endp)
    logical                     :: patch_grew               (bounds%begp:bounds%endp)

    ! The following are only set for growing patches:
    real(r8)                    :: seed_leafn_patch         (bounds%begp:bounds%endp)
    real(r8)                    :: seed_leafn_storage_patch (bounds%begp:bounds%endp)
    real(r8)                    :: seed_leafn_xfer_patch    (bounds%begp:bounds%endp)
    real(r8)                    :: seed_deadstemn_patch     (bounds%begp:bounds%endp)
    real(r8)                    :: seed_npool_patch         (bounds%begp:bounds%endp)

    real(r8)                    :: wood_product_nflux       (bounds%begp:bounds%endp)
    real(r8)                    :: deadstemn_patch_temp     (bounds%begp:bounds%endp)

    character(len=*), parameter :: subname = 'NStateDynamicPatchAdjustments'
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    SHR_ASSERT_ALL((ubound(dwt_leafn_seed           ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_deadstemn_seed       ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_npool_seed           ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(conv_nflux               ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_frootn_to_litter     ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_livecrootn_to_litter ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_deadcrootn_to_litter ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(prod10_nflux             ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(prod100_nflux            ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(crop_product_nflux       ) == (/endp/)), errMsg(__FILE__, __LINE__))
   
    old_weight_was_zero = patch_state_updater%old_weight_was_zero(bounds)
    patch_grew          = patch_state_updater%patch_grew(bounds)

    call ComputeSeedAmounts(bounds                                        , &
         species                    = NUTRIENT_SPECIES_N                  , &
         leaf_patch                 = this%leaf_patch(begp:endp)          , &
         leaf_storage_patch         = this%leaf_storage_patch(begp:endp)  , &
         leaf_xfer_patch            = this%leaf_xfer_patch(begp:endp)     , &

         ! Calculations only needed for patches that grew:
         compute_here_patch         = patch_grew(begp:endp)               , &

         ! For patches that previously had zero area, ignore the current state for the
         ! sake of computing leaf proportions:
         ignore_current_state_patch = old_weight_was_zero(begp:endp)      , &

         seed_leaf_patch            = seed_leafn_patch(begp:endp)         , &
         seed_leaf_storage_patch    = seed_leafn_storage_patch(begp:endp) , &
         seed_leaf_xfer_patch       = seed_leafn_xfer_patch(begp:endp)    , &
         seed_deadstem_patch        = seed_deadstemn_patch(begp:endp)     , &
         pool_seed_param            = npool_seed_param                    , &
         pool_seed_patch            = seed_npool_patch(begp:endp))

    ! 1) LEAFN_PATCH
    call patch_state_updater%update_patch_state(            &
         bounds                                           , &
         num_filterp_with_inactive                        , &
         filterp_with_inactive                            , &
         var               = this%leaf_patch   (begp:endp) , &
         flux_out_grc_area = conv_nflux       (begp:endp) , &
         seed              = seed_leafn_patch (begp:endp) , &
         seed_addition     = dwt_leafn_seed   (begp:endp))

    ! 2) LEAFN_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                    &
         bounds                                                   , &
         num_filterp_with_inactive                                , &
         filterp_with_inactive                                    , &
         var               = this%leaf_storage_patch   (begp:endp) , &
         flux_out_grc_area = conv_nflux               (begp:endp) , &
         seed              = seed_leafn_storage_patch (begp:endp) , &
         seed_addition     = dwt_leafn_seed           (begp:endp))

    ! 3) LEAFN_XFER_PATCH
    call patch_state_updater%update_patch_state( &
         bounds                                                        , &
         num_filterp_with_inactive                                     , &
         filterp_with_inactive                                         , &
         var               = this%leaf_xfer_patch   (begp:endp), &
         flux_out_grc_area = conv_nflux            (begp:endp), &
         seed              = seed_leafn_xfer_patch (begp:endp), &
         seed_addition     = dwt_leafn_seed        (begp:endp))

    ! 4) FROOTN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%froot_patch(begp:endp)             , &
         flux_out_col_area = dwt_frootn_to_litter(begp:endp))

    ! 5) FROOTN_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%froot_storage_patch(begp:endp)     , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 6) FROOTN_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%froot_xfer_patch(begp:endp)        , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 7) LIVESTEMN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livestem_patch(begp:endp)          , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 8) LIVESTEMN_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livestem_storage_patch(begp:endp)  , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 9) LIVESTEMN_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livestem_xfer_patch(begp:endp)     , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 10) PROD10_NFLUX
    wood_product_nflux(begp:endp)      = 0._r8
    deadstemn_patch_temp(begp:endp)    = this%deadstem_patch(begp:endp)
    call patch_state_updater%update_patch_state_partition_flux_by_type(     &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pprod10                             , &
         var                        = deadstemn_patch_temp    (begp:endp) , &
         flux1_out                  = prod10_nflux            (begp:endp) , &
         flux2_out                  = wood_product_nflux      (begp:endp) , &
         seed                       = seed_deadstemn_patch    (begp:endp) )

    ! 11) PROD100_NFLUX
    wood_product_nflux(begp:endp)      = 0._r8
    deadstemn_patch_temp(begp:endp)    = this%deadstem_patch(begp:endp)
    call patch_state_updater%update_patch_state_partition_flux_by_type(     &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pprod100                            , &
         var                        = deadstemn_patch_temp    (begp:endp) , &
         flux1_out                  = prod100_nflux           (begp:endp) , &
         flux2_out                  = wood_product_nflux      (begp:endp) , &
         seed                       = seed_deadstemn_patch    (begp:endp))

    ! 12) DEADSTEMN_PATCH
    wood_product_nflux(begp:endp)      = 0._r8
    call patch_state_updater%update_patch_state_partition_flux_by_type(     &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pconv                               , &
         var                        = this%deadstem_patch   (begp:endp)    , &
         flux1_out                  = conv_nflux           (begp:endp)    , &
         flux2_out                  = wood_product_nflux   (begp:endp)    , &
         seed                       = seed_deadstemn_patch (begp:endp)    , &
         seed_addition              = dwt_deadstemn_seed   (begp:endp))

    ! 13) DEADSTEMN_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadstem_storage_patch(begp:endp)  , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 14) DEADSTEMN_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadstem_xfer_patch(begp:endp)     , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 15) LIVECROOTN_PATCH 
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livecroot_patch(begp:endp)         , &
         flux_out_col_area = dwt_livecrootn_to_litter(begp:endp))

    ! 16) LIVECROOTN_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livecroot_storage_patch(begp:endp) , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 17) LIVECROOTN_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livecroot_xfer_patch(begp:endp)    , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 18) DEADCROOTN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadcroot_patch(begp:endp)         , &
         flux_out_col_area = dwt_deadcrootn_to_litter(begp:endp))

    ! 19) DEADCROOTN_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadcroot_storage_patch(begp:endp) , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 20) DEADCROOT_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadcroot_xfer_patch(begp:endp)    , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 21) RETRANSN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%retransn_patch(begp:endp)           , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 22) NTRUNC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%veg_trunc_patch(begp:endp)             , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 23) CPOOL_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%pool_patch(begp:endp)              , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 24) DISPVEGN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%dispveg_patch(begp:endp))

    ! 25) STORVEGN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%storveg_patch(begp:endp))

    ! 26) TOTVEGN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%totveg_patch(begp:endp))

    ! 27) TOTPFTN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%totpft_patch(begp:endp))

    if (use_crop) then
       ! This is a negative pool. So any deficit that we haven't repaid gets sucked out
       ! of the atmosphere.
       call patch_state_updater%update_patch_state(         &
            bounds                                        , &
            num_filterp_with_inactive                     , &
            filterp_with_inactive                         , &
            var = this%cropseed_deficit_patch(begp:endp) , &
            flux_out_grc_area = conv_nflux(begp:endp))
    end if

    ! These fluxes are computed as negative quantities, but are expected to be positive,
    ! so flip the signs
    do p = begp,endp
       dwt_frootn_to_litter(p)     = -1._r8 * dwt_frootn_to_litter(p)
       dwt_livecrootn_to_litter(p) = -1._r8 * dwt_livecrootn_to_litter(p)
       dwt_deadcrootn_to_litter(p) = -1._r8 * dwt_deadcrootn_to_litter(p)
    end do

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
    class(nitrogenstate_type)       , intent(inout) :: this
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

    !this%dyn_bal_adjustments_col(begc:endc) = 0._r8

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

       call column_state_updater%update_column_state_no_special_handling( &
           bounds      = bounds                          , &
           clump_index = clump_index                     , &
           var         = this%sminn_vr_col(begc:endc, j) , &
           adjustment  = adjustment_one_level(begc:endc))

       this%dyn_bal_adjustments_col(begc:endc) = &
           this%dyn_bal_adjustments_col(begc:endc) + &
           adjustment_one_level(begc:endc) * dzsoi_decomp(j)
    end do

  end subroutine DynamicColumnAdjustments

end module CNNitrogenStateType
