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
  use CNSpeciesMod           , only : species_from_string, species_name_from_string
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

    this%species = species_from_string('p')
    this%name    = species_name_from_string('p')

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

    !-------------------------------
    ! P state variables - native to PFT
    !-------------------------------
    
    if (crop_prog) then
       this%grain_patch(begp:endp) = spval
       call hist_addfld1d (fname='GRAINP', units='gP/m^2', &
            avgflag='A', long_name='grain P', &
            ptr_patch=this%grain_patch, default='inactive')

       this%cropseed_deficit_patch(begp:endp) = spval
       call hist_addfld1d (fname='CROPSEEDP_DEFICIT', units='gP/m^2', &
            avgflag='A', long_name='P used for crop seed that needs to be repaid', &
            ptr_patch=this%cropseed_deficit_patch)
    end if

    this%leaf_patch(begp:endp) = spval
    call hist_addfld1d (fname='LEAFP', units='gP/m^2', &
         avgflag='A', long_name='leaf P', &
         ptr_patch=this%leaf_patch)

    this%leaf_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='LEAFP_STORAGE', units='gP/m^2', &
         avgflag='A', long_name='leaf P storage', &
         ptr_patch=this%leaf_storage_patch, default='inactive')

    this%leaf_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='LEAFP_XFER', units='gP/m^2', &
         avgflag='A', long_name='leaf P transfer', &
         ptr_patch=this%leaf_xfer_patch, default='inactive')

    this%froot_patch(begp:endp) = spval
    call hist_addfld1d (fname='FROOTP', units='gP/m^2', &
         avgflag='A', long_name='fine root P', &
         ptr_patch=this%froot_patch)

    this%froot_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='FROOTP_STORAGE', units='gP/m^2', &
         avgflag='A', long_name='fine root P storage', &
         ptr_patch=this%froot_storage_patch, default='inactive')

    this%froot_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='FROOTP_XFER', units='gP/m^2', &
         avgflag='A', long_name='fine root P transfer', &
         ptr_patch=this%froot_xfer_patch, default='inactive')

    this%livestem_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMP', units='gP/m^2', &
         avgflag='A', long_name='live stem P', &
         ptr_patch=this%livestem_patch)

    this%livestem_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMP_STORAGE', units='gP/m^2', &
         avgflag='A', long_name='live stem P storage', &
         ptr_patch=this%livestem_storage_patch, default='inactive')

    this%livestem_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMP_XFER', units='gP/m^2', &
         avgflag='A', long_name='live stem P transfer', &
         ptr_patch=this%livestem_xfer_patch, default='inactive')

    this%deadstem_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMP', units='gP/m^2', &
         avgflag='A', long_name='dead stem P', &
         ptr_patch=this%deadstem_patch)

    this%deadstem_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMP_STORAGE', units='gP/m^2', &
         avgflag='A', long_name='dead stem P storage', &
         ptr_patch=this%deadstem_storage_patch, default='inactive')

    this%deadstem_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMP_XFER', units='gP/m^2', &
         avgflag='A', long_name='dead stem P transfer', &
         ptr_patch=this%deadstem_xfer_patch, default='inactive')

    this%livecroot_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTP', units='gP/m^2', &
         avgflag='A', long_name='live coarse root P', &
         ptr_patch=this%livecroot_patch)

    this%livecroot_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTP_STORAGE', units='gP/m^2', &
         avgflag='A', long_name='live coarse root P storage', &
         ptr_patch=this%livecroot_storage_patch, default='inactive')

    this%livecroot_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTP_XFER', units='gP/m^2', &
         avgflag='A', long_name='live coarse root P transfer', &
         ptr_patch=this%livecroot_xfer_patch, default='inactive')

    this%deadcroot_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTP', units='gP/m^2', &
         avgflag='A', long_name='dead coarse root P', &
         ptr_patch=this%deadcroot_patch)

    this%deadcroot_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTP_STORAGE', units='gP/m^2', &
         avgflag='A', long_name='dead coarse root P storage', &
         ptr_patch=this%deadcroot_storage_patch, default='inactive')

    this%deadcroot_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTP_XFER', units='gP/m^2', &
         avgflag='A', long_name='dead coarse root P transfer', &
         ptr_patch=this%deadcroot_xfer_patch, default='inactive')

    this%retransp_patch(begp:endp) = spval
    call hist_addfld1d (fname='RETRANSP', units='gP/m^2', &
         avgflag='A', long_name='plant pool of retranslocated P', &
         ptr_patch=this%retransp_patch)

    this%pool_patch(begp:endp) = spval
    call hist_addfld1d (fname='PPOOL', units='gP/m^2', &
         avgflag='A', long_name='temporary plant P pool', &
         ptr_patch=this%pool_patch, default='inactive')

    this%veg_trunc_patch(begp:endp) = spval
    call hist_addfld1d (fname='PFT_PTRUNC', units='gP/m^2', &
         avgflag='A', long_name='pft-level sink for P truncation', &
         ptr_patch=this%veg_trunc_patch, default='inactive')

    this%dispveg_patch(begp:endp) = spval
    call hist_addfld1d (fname='DISPVEGP', units='gP/m^2', &
         avgflag='A', long_name='displayed vegetation phosphorus', &
         ptr_patch=this%dispveg_patch)

    this%storveg_patch(begp:endp) = spval
    call hist_addfld1d (fname='STORVEGP', units='gP/m^2', &
         avgflag='A', long_name='stored vegetation phosphorus', &
         ptr_patch=this%storveg_patch)

    this%totveg_patch(begp:endp) = spval
    call hist_addfld1d (fname='TOTVEGP', units='gP/m^2', &
         avgflag='A', long_name='total vegetation phosphorus', &
         ptr_patch=this%totveg_patch)

    this%totpft_patch(begp:endp) = spval
    call hist_addfld1d (fname='TOTPFTP', units='gP/m^2', &
         avgflag='A', long_name='total PFT-level phosphorus', &
         ptr_patch=this%totpft_patch)

    this%plant_p_buffer_patch(begp:endp) = spval
    call hist_addfld1d (fname='PLANTP_BUFFER', units='gP/m^2', &
            avgflag='A', long_name='plant phosphorus stored as buffer', &
            ptr_col=this%plant_p_buffer_patch,default='inactive')
    !-------------------------------
    ! P state variables - native to column
    !-------------------------------

    if ( nlevdecomp_full > 1 ) then
       this%decomp_pools_vr_col(begc:endc,:,:) = spval
       this%decomp_pools_1m_col(begc:endc,:) = spval
    end if
    this%decomp_pools_col(begc:endc,:) = spval
    do l  = 1, ndecomp_pools
       if ( nlevdecomp_full > 1 ) then
          data2dptr => this%decomp_pools_vr_col(:,:,l)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'P_vr'
          longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' P (vertically resolved)'
          call hist_addfld2d (fname=fieldname, units='gP/m^3',  type2d='levdcmp', &
               avgflag='A', long_name=longname, &
               ptr_col=data2dptr)
       endif

       data1dptr => this%decomp_pools_col(:,l)
       fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'P'
       longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' P'
       call hist_addfld1d (fname=fieldname, units='gP/m^2', &
            avgflag='A', long_name=longname, &
            ptr_col=data1dptr)

       if ( nlevdecomp_full > 1 ) then
          data1dptr => this%decomp_pools_1m_col(:,l)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'P_1m'
          longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' P to 1 meter'
          call hist_addfld1d (fname=fieldname, units='gP/m^2', &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr, default = 'inactive')
       endif
    end do


    if ( nlevdecomp_full > 1 ) then

       this%sminp_col(begc:endc) = spval
       call hist_addfld1d (fname='SMINP', units='gP/m^2', &
            avgflag='A', long_name='soil mineral P', &
            ptr_col=this%sminp_col)

       this%totlit_1m_col(begc:endc) = spval
       call hist_addfld1d (fname='TOTLITP_1m', units='gP/m^2', &
            avgflag='A', long_name='total litter P to 1 meter', &
            ptr_col=this%totlit_1m_col)

       this%totsom_1m_col(begc:endc) = spval
       call hist_addfld1d (fname='TOTSOMP_1m', units='gP/m^2', &
            avgflag='A', long_name='total soil organic matter P to 1 meter', &
            ptr_col=this%totsom_1m_col)
    endif

    this%veg_trunc_col(begc:endc) = spval
    call hist_addfld1d (fname='COL_PTRUNC', units='gP/m^2',  &
         avgflag='A', long_name='column-level sink for P truncation', &
         ptr_col=this%veg_trunc_col)

    ! add suffix if number of soil decomposition depths is greater than 1
    if (nlevdecomp > 1) then
       vr_suffix = "_vr"
    else 
       vr_suffix = ""
    endif

    this%solutionp_vr_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='SOLUTIONP'//trim(vr_suffix), units='gp/m^3',  type2d='levdcmp', &
         avgflag='A', long_name='soil solution P (vert. res.)', &
         ptr_col=this%solutionp_vr_col)

    this%labilep_vr_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='LABILEP'//trim(vr_suffix), units='gp/m^3',  type2d='levdcmp', &
         avgflag='A', long_name='soil labile P (vert. res.)', &
         ptr_col=this%labilep_vr_col)

    this%secondp_vr_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='SECONDP'//trim(vr_suffix), units='gp/m^3',  type2d='levdcmp', &
         avgflag='A', long_name='soil secondary P (vert. res.)', &
         ptr_col=this%secondp_vr_col)

    this%occlp_vr_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='OCCLP'//trim(vr_suffix), units='gp/m^3',  type2d='levdcmp', &
         avgflag='A', long_name='soil occluded P (vert. res.)', &
         ptr_col=this%occlp_vr_col)

    this%primp_vr_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='PRIMP'//trim(vr_suffix), units='gp/m^3',  type2d='levdcmp', &
         avgflag='A', long_name='soil primary P (vert. res.)', &
         ptr_col=this%primp_vr_col)

    this%sminp_vr_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='SMINP'//trim(vr_suffix), units='gp/m^3',  type2d='levdcmp', &
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

    this%totlit_col(begc:endc) = spval
    call hist_addfld1d (fname='TOTLITP', units='gP/m^2', &
         avgflag='A', long_name='total litter P', &
         ptr_col=this%totlit_col)

    this%totsom_col(begc:endc) = spval
    call hist_addfld1d (fname='TOTSOMP', units='gP/m^2', &
         avgflag='A', long_name='total soil organic matter P', &
         ptr_col=this%totsom_col)

    this%totecosys_col(begc:endc) = spval
    call hist_addfld1d (fname='TOTECOSYSP', units='gP/m^2', &
         avgflag='A', long_name='total ecosystem P but excl product pools', &
         ptr_col=this%totecosys_col)

    this%totcol_col(begc:endc) = spval
    call hist_addfld1d (fname='TOTCOLP', units='gP/m^2', &
         avgflag='A', long_name='total column-level P but excl product pools', &
         ptr_col=this%totcol_col)

    this%seed_grc(begg:endg) = spval
    call hist_addfld1d (fname='SEEDP_GRC', units='gP/m^2', &
         avgflag='A', long_name='P pool for seeding new PFTs ', &
         ptr_gcell=this%seed_grc, default='inactive')

    this%seed_col(begc:endc) = spval
    call hist_addfld1d (fname='SEEDP', units='gP/m^2', &
         avgflag='A', long_name='P pool for seeding new PFTs ', &
         ptr_col=this%seed_col, default='inactive')

    this%prod10_col(begc:endc) = spval
    call hist_addfld1d (fname='PROD10P', units='gP/m^2', &
         avgflag='A', long_name='10-yr wood product P', &
         ptr_col=this%prod10_col, default='inactive')

    this%prod100_col(begc:endc) = spval
    call hist_addfld1d (fname='PROD100P', units='gP/m^2', &
         avgflag='A', long_name='100-yr wood product P', &
         ptr_col=this%prod100_col, default='inactive')

    this%prod1_col(begc:endc) = spval
    call hist_addfld1d (fname='PROD1P', units='gP/m^2', &
         avgflag='A', long_name='1-yr crop product P', &
         ptr_col=this%prod1_col, default='inactive')

    this%totprod_col(begc:endc) = spval
    call hist_addfld1d (fname='TOTPRODP', units='gP/m^2', &
         avgflag='A', long_name='total wood product P', &
         ptr_col=this%totprod_col, default='inactive')

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


    call restartvar(ncid=ncid, flag=flag, varname='leafp', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%leaf_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='leafp_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%leaf_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='leafp_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%leaf_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='frootp', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%froot_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='frootp_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%froot_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='frootp_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%froot_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livestemp', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livestem_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livestemp_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livestem_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livestemp_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livestem_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadstemp', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadstem_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadstemp_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadstem_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadstemp_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadstem_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livecrootp', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livecroot_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livecrootp_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livecroot_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livecrootp_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livecroot_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadcrootp', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadcroot_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadcrootp_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadcroot_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadcrootp_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadcroot_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='retransp', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%retransp_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='ppool', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%pool_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='pft_ptrunc', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%veg_trunc_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='plant_p_buffer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%plant_p_buffer_patch)

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag,  varname='grainp', xtype=ncd_double,  &
            dim1name='pft',    long_name='grain P', units='gP/m2', &
            interpinic_flag='interp', readvar=readvar, data=this%grain_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='grainp_storage', xtype=ncd_double,  &
            dim1name='pft',    long_name='grain P storage', units='gP/m2', &
            interpinic_flag='interp', readvar=readvar, data=this%grain_storage_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='grainp_xfer', xtype=ncd_double,  &
            dim1name='pft',    long_name='grain P transfer', units='gP/m2', &
            interpinic_flag='interp', readvar=readvar, data=this%grain_xfer_patch)
    end if

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

    ! decomposing P pools
    do k = 1, ndecomp_pools
       varname=trim(decomp_cascade_con%decomp_pool_name_restart(k))//'p'
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
       call restartvar(ncid=ncid, flag=flag, varname="col_ptrunc_vr", xtype=ncd_double,  &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)
    else
       ptr1d => this%soil_trunc_vr_col(:,1)
       call restartvar(ncid=ncid, flag=flag, varname="col_ptrunc", xtype=ncd_double,  &
            dim1name='column', &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp' , readvar=readvar, data=ptr1d)
    end if

    !!! Debug balance 
    call restartvar(ncid=ncid, flag=flag, varname='totsomp', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%totsom_col) 
    call restartvar(ncid=ncid, flag=flag, varname='cwdp', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%cwd_col) 
    call restartvar(ncid=ncid, flag=flag, varname='totlitp', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%totlit_col) 

    call restartvar(ncid=ncid, flag=flag, varname='totcolp', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%totcol_col) 

    call restartvar(ncid=ncid, flag=flag, varname='seedp', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%seed_col) 

    call restartvar(ncid=ncid, flag=flag, varname='prod10p', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%prod10_col) 

    call restartvar(ncid=ncid, flag=flag, varname='prod100p', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%prod100_col) 

    call restartvar(ncid=ncid, flag=flag, varname='prod1p', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%prod1_col)

    ! decomp_cascade_state - the purpose of this is to check to make sure the bgc used 
    ! matches what the restart file was generated with.  
    ! add info about the SOM decomposition cascade

!    if (use_century_decomp) then
!       decomp_cascade_state = 1
!    else
!       decomp_cascade_state = 0
!    end if
    ! add info about the nitrification / denitrification state
!    if (use_nitrif_denitrif) then
!       decomp_cascade_state = decomp_cascade_state + 10
!    end if
!    if (flag == 'write') itemp = decomp_cascade_state    
!    call restartvar(ncid=ncid, flag=flag, varname='decomp_cascade_state', xtype=ncd_int,  &
!         long_name='BGC of the model that wrote this restart file:' &
!         // '  1s column: 0 = CLM-CN cascade, 1 = Century cascade;' &
!         // ' 10s column: 0 = CLM-CN denitrification, 10 = Century denitrification', units='', &
!         interpinic_flag='skip', readvar=readvar, data=itemp)
!    if (flag=='read') then
!       if (.not. readvar) then
!          ! assume, for sake of backwards compatibility, that if decomp_cascade_state 
!          ! is not in the restart file, then the current model state is the same as 
!          ! the prior model state
!          restart_file_decomp_cascade_state = decomp_cascade_state
!          if ( masterproc ) write(iulog,*) ' CNRest: WARNING!  Restart file does not ' &
!               // ' contain info on decomp_cascade_state used to generate the restart file.  '
!          if ( masterproc ) write(iulog,*) '   Assuming the same as current setting: ', decomp_cascade_state
!       else
!          restart_file_decomp_cascade_state = itemp  
!          if (decomp_cascade_state /= restart_file_decomp_cascade_state ) then
!             if ( masterproc ) then
!                write(iulog,*) 'CNRest: ERROR--the decomposition cascade differs between the current ' &
!                     // ' model state and the model that wrote the restart file. '
!                write(iulog,*) 'The model will be horribly out of equilibrium until after a lengthy spinup. '
!                write(iulog,*) 'Stopping here since this is probably an error in configuring the run. '
!                write(iulog,*) 'If you really wish to proceed, then override by setting '
!                write(iulog,*) 'override_bgc_restart_mismatch_dump to .true. in the namelist'
!                if ( .not. override_bgc_restart_mismatch_dump ) then
!                   call endrun(msg= ' CNRest: Stopping. Decomposition cascade mismatch error.'//&
!                        errMsg(__FILE__, __LINE__))
!                endif
!             endif
!          endif
!       end if
!    end if

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

    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! displayed vegetation phosphorus, excluding storage (DISPVEGN)
       this%dispveg_patch(p) = &
            this%leaf_patch(p)      + &
            this%froot_patch(p)     + &
            this%livestem_patch(p)  + &
            this%deadstem_patch(p)  + &
            this%livecroot_patch(p) + &
            this%deadcroot_patch(p)
       
      ! stored vegetation phosphorus, including retranslocated N pool (STORVEGN)
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
           this%retransp_patch(p)

      if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
         this%dispveg_patch(p) = &
              this%dispveg_patch(p) + &
              this%grain_patch(p)

         this%storveg_patch(p) = &
              this%storveg_patch(p) + &
              this%grain_storage_patch(p)     + &
              this%grain_xfer_patch(p)
      end if

      ! total vegetation phosphorus (TOTVEGN)
      this%totveg_patch(p) = &
           this%dispveg_patch(p) + &
           this%storveg_patch(p)

      ! total pft-level carbon (add pft_ntrunc)
      this%totpft_patch(p) = &
           this%totveg_patch(p) + &
           this%veg_trunc_patch(p)

   end do

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
      this%solutionp_col(c) =0._r8
      this%labilep_col(c)   =0._r8
      this%secondp_col(c)   =0._r8
      this%occlp_col(c)     =0._r8
      this%primp_col(c)     =0._r8
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


   ! vertically integrate each of the decomposing P pools
   do l = 1, ndecomp_pools
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%decomp_pools_col(c,l) = 0._r8
      end do
      do j = 1, nlevdecomp
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
      
      ! total litter phosphorus to 1 meter (TOTLITN_1m)
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
      
      ! total soil organic matter phosphorus to 1 meter (TOTSOMN_1m)
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
   
   ! total litter phosphorus (TOTLITN)
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
   
   ! total soil organic matter phosphorus (TOTSOMN)
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
!                                  this%occlp_vr_col(c,j)+ &
!                                  this%primp_vr_col(c,j)
      end do
   end do
   do j = 1, nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%sminp_col(c) =  this%sminp_col(c) + &
         this%sminp_vr_col(c,j) * dzsoi_decomp(j)
      end do
   end do


   ! total col_ntrunc
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      this%veg_trunc_col(c) = 0._r8
   end do
   do j = 1, nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%veg_trunc_col(c) = &
              this%veg_trunc_col(c) + &
              this%soil_trunc_vr_col(c,j) * dzsoi_decomp(j)
      end do
   end do

   do fc = 1,num_soilc
      c = filter_soilc(fc)

      ! total wood product phosphorus
      this%totprod_col(c) = &
           this%prod1_col(c) + &
           this%prod10_col(c) + &
           this%prod100_col(c)	 

      ! total ecosystem phosphorus, including veg (TOTECOSYSP)
      this%totecosys_col(c) = &
           this%cwd_col(c) + &
           this%totlit_col(c) + &
           this%totsom_col(c) + &
           this%solutionp_col(c) + &
           this%labilep_col(c) + &
           this%secondp_col(c) + &
           this%primp_col(c) + &
           this%occlp_col(c) + &
           this%totprod_col(c) + &
           this%totveg_col(c)

      ! total column phosphorus, including pft (TOTCOLP)
      this%totcol_col(c) = &
           this%totpft_col(c) + &
           this%cwd_col(c) + &
           this%totlit_col(c) + &
           this%totsom_col(c) + &
           this%prod1_col(c) + &
           this%solutionp_col(c) + &
           this%labilep_col(c) + &
           this%secondp_col(c) + &
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
   
    old_weight_was_zero = patch_state_updater%old_weight_was_zero(bounds)
    patch_grew          = patch_state_updater%patch_grew(bounds)

    call ComputeSeedAmounts(bounds                                        , &
         species                    = NUTRIENT_SPECIES_P                  , &
         leaf_patch                 = this%leaf_patch(begp:endp)          , &
         leaf_storage_patch         = this%leaf_storage_patch(begp:endp)  , &
         leaf_xfer_patch            = this%leaf_xfer_patch(begp:endp)     , &

         ! Calculations only needed for patches that grew:
         compute_here_patch         = patch_grew(begp:endp)               , &

         ! For patches that previously had zero area, ignore the current state for the
         ! sake of computing leaf proportions:
         ignore_current_state_patch = old_weight_was_zero(begp:endp)      , &

         seed_leaf_patch            = seed_leafp_patch(begp:endp)         , &
         seed_leaf_storage_patch    = seed_leafp_storage_patch(begp:endp) , &
         seed_leaf_xfer_patch       = seed_leafp_xfer_patch(begp:endp)    , &
         seed_deadstem_patch        = seed_deadstemp_patch(begp:endp)     , &
         pool_seed_param            = ppool_seed_param                    , &
         pool_seed_patch            = seed_ppool_patch(begp:endp))

    ! 1) LEAFP_PATCH
    call patch_state_updater%update_patch_state(            &
         bounds                                           , &
         num_filterp_with_inactive                        , &
         filterp_with_inactive                            , &
         var               = this%leaf_patch   (begp:endp) , &
         flux_out_grc_area = conv_pflux       (begp:endp) , &
         seed              = seed_leafp_patch (begp:endp) , &
         seed_addition     = dwt_leafp_seed   (begp:endp))

    ! 2) LEAFP_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                    &
         bounds                                                   , &
         num_filterp_with_inactive                                , &
         filterp_with_inactive                                    , &
         var               = this%leaf_storage_patch   (begp:endp) , &
         flux_out_grc_area = conv_pflux               (begp:endp) , &
         seed              = seed_leafp_storage_patch (begp:endp) , &
         seed_addition     = dwt_leafp_seed           (begp:endp))

    ! 3) LEAFP_XFER_PATCH
    call patch_state_updater%update_patch_state( &
         bounds                                                        , &
         num_filterp_with_inactive                                     , &
         filterp_with_inactive                                         , &
         var               = this%leaf_xfer_patch   (begp:endp), &
         flux_out_grc_area = conv_pflux            (begp:endp), &
         seed              = seed_leafp_xfer_patch (begp:endp), &
         seed_addition     = dwt_leafp_seed        (begp:endp))

    ! 4) FROOTP_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%froot_patch(begp:endp)             , &
         flux_out_col_area = dwt_frootp_to_litter(begp:endp))

    ! 5) FROOTP_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%froot_storage_patch(begp:endp)     , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 6) FROOTP_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%froot_xfer_patch(begp:endp)        , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 7) PROD10_PFLUX
    wood_product_pflux(begp:endp)      = 0._r8
    deadstemp_patch_temp(begp:endp)    = this%deadstem_patch(begp:endp)
    call patch_state_updater%update_patch_state_partition_flux_by_type(     &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pprod10                             , &
         var                        = deadstemp_patch_temp    (begp:endp) , &
         flux1_out                  = prod10_pflux            (begp:endp) , &
         flux2_out                  = wood_product_pflux      (begp:endp) , &
         seed                       = seed_deadstemp_patch    (begp:endp) )

    ! 8) PROD100_PFLUX
    wood_product_pflux(begp:endp)      = 0._r8
    deadstemp_patch_temp(begp:endp)    = this%deadstem_patch(begp:endp)
    call patch_state_updater%update_patch_state_partition_flux_by_type(     &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pprod100                            , &
         var                        = deadstemp_patch_temp    (begp:endp) , &
         flux1_out                  = prod100_pflux           (begp:endp) , &
         flux2_out                  = wood_product_pflux      (begp:endp) , &
         seed                       = seed_deadstemp_patch    (begp:endp))

    ! 9) DEADSTEMP_PATCH
    wood_product_pflux(begp:endp)      = 0._r8
    call patch_state_updater%update_patch_state_partition_flux_by_type(     &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pconv                               , &
         var                        = this%deadstem_patch   (begp:endp)    , &
         flux1_out                  = conv_pflux           (begp:endp)    , &
         flux2_out                  = wood_product_pflux   (begp:endp)    , &
         seed                       = seed_deadstemp_patch (begp:endp)    , &
         seed_addition              = dwt_deadstemp_seed   (begp:endp))

    ! 10) DEADSTEMP_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadstem_storage_patch(begp:endp)  , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 11) DEADSTEMP_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadstem_xfer_patch(begp:endp)     , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 12) LIVESTEMP_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livestem_patch(begp:endp)          , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 13) LIVESTEMP_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livestem_storage_patch(begp:endp)  , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 14) LIVESTEMP_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livestem_xfer_patch(begp:endp)     , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 15) LIVECROOTP_PATCH 
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livecroot_patch(begp:endp)         , &
         flux_out_col_area = dwt_livecrootp_to_litter(begp:endp))

    ! 16) LIVECROOTP_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livecroot_storage_patch(begp:endp) , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 17) LIVECROOTP_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livecroot_xfer_patch(begp:endp)    , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 18) DEADCROOTP_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadcroot_patch(begp:endp)         , &
         flux_out_col_area = dwt_deadcrootp_to_litter(begp:endp))

    ! 19) DEADCROOTP_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadcroot_storage_patch(begp:endp) , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 20) DEADCROOT_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadcroot_xfer_patch(begp:endp)    , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 21) RETRANSP_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%retransp_patch(begp:endp)           , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 22) PTRUNC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%veg_trunc_patch(begp:endp)             , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 23) PPOOL_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%pool_patch(begp:endp)              , &
         flux_out_grc_area = conv_pflux(begp:endp))

    ! 24) DISPVEGP_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%dispveg_patch(begp:endp))

    ! 25) STORVEGP_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%storveg_patch(begp:endp))

    ! 26) TOTVEGP_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%totveg_patch(begp:endp))

    ! 27) TOTPFTP_PATCH
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
            flux_out_grc_area = conv_pflux(begp:endp))
    end if

    ! These fluxes are computed as negative quantities, but are expected to be positive,
    ! so flip the signs
    do p = begp,endp
       dwt_frootp_to_litter(p)     = -1._r8 * dwt_frootp_to_litter(p)
       dwt_livecrootp_to_litter(p) = -1._r8 * dwt_livecrootp_to_litter(p)
       dwt_deadcrootp_to_litter(p) = -1._r8 * dwt_deadcrootp_to_litter(p)
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

    !this%dyn_bal_adjustments_col(begc:endc) = 0._r8

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
