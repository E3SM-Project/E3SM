module PhosphorusStateType

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
  use pftvarcon              , only : npcropmin
  use CNDecompCascadeConType , only : decomp_cascade_con
  use EcophysConType         , only : ecophyscon
  use abortutils             , only : endrun
  use spmdMod                , only : masterproc 
  use LandunitType           , only : lun                
  use ColumnType             , only : col                
  use PatchType              , only : pft
  use clm_varctl             , only : nu_com
              
  ! 
  ! !PUBLIC TYPES:
  implicit none
  save
  private

  type, public :: phosphorusstate_type

     real(r8), pointer :: grainp_patch                 (:)     ! patch (gP/m2) grain P (crop)
     real(r8), pointer :: grainp_storage_patch         (:)     ! patch (gP/m2) grain P storage (crop)
     real(r8), pointer :: grainp_xfer_patch            (:)     ! patch (gP/m2) grain P transfer (crop)
     real(r8), pointer :: leafp_patch                  (:)     ! patch (gP/m2) leaf P 
     real(r8), pointer :: leafp_storage_patch          (:)     ! patch (gP/m2) leaf P storage
     real(r8), pointer :: leafp_xfer_patch             (:)     ! patch (gP/m2) leaf P transfer
     real(r8), pointer :: frootp_patch                 (:)     ! patch (gP/m2) fine root P
     real(r8), pointer :: frootp_storage_patch         (:)     ! patch (gP/m2) fine root P storage
     real(r8), pointer :: frootp_xfer_patch            (:)     ! patch (gP/m2) fine root P transfer
     real(r8), pointer :: livestemp_patch              (:)     ! patch (gP/m2) live stem P
     real(r8), pointer :: livestemp_storage_patch      (:)     ! patch (gP/m2) live stem P storage
     real(r8), pointer :: livestemp_xfer_patch         (:)     ! patch (gP/m2) live stem P transfer
     real(r8), pointer :: deadstemp_patch              (:)     ! patch (gP/m2) dead stem P
     real(r8), pointer :: deadstemp_storage_patch      (:)     ! patch (gP/m2) dead stem P storage
     real(r8), pointer :: deadstemp_xfer_patch         (:)     ! patch (gP/m2) dead stem P transfer
     real(r8), pointer :: livecrootp_patch             (:)     ! patch (gP/m2) live coarse root P
     real(r8), pointer :: livecrootp_storage_patch     (:)     ! patch (gP/m2) live coarse root P storage
     real(r8), pointer :: livecrootp_xfer_patch        (:)     ! patch (gP/m2) live coarse root P transfer
     real(r8), pointer :: deadcrootp_patch             (:)     ! patch (gP/m2) dead coarse root P
     real(r8), pointer :: deadcrootp_storage_patch     (:)     ! patch (gP/m2) dead coarse root P storage
     real(r8), pointer :: deadcrootp_xfer_patch        (:)     ! patch (gP/m2) dead coarse root P transfer
     real(r8), pointer :: retransp_patch               (:)     ! patch (gP/m2) plant pool of retranslocated P
     real(r8), pointer :: ppool_patch                  (:)     ! patch (gP/m2) temporary plant P pool
     real(r8), pointer :: ptrunc_patch                 (:)     ! patch (gP/m2) pft-level sink for P truncation

     real(r8), pointer :: decomp_ppools_vr_col         (:,:,:)     ! col (gP/m3) vertically-resolved decomposing (litter, cwd, soil) P pools 
     real(r8), pointer :: solutionp_vr_col             (:,:)       ! col (gP/m3) vertically-resolved soil solution P
     real(r8), pointer :: labilep_vr_col               (:,:)       ! col (gP/m3) vertically-resolved soil labile mineral P
     real(r8), pointer :: secondp_vr_col               (:,:)       ! col (gP/m3) vertically-resolved soil secondary mineralP
     real(r8), pointer :: occlp_vr_col                 (:,:)       ! col (gP/m3) vertically-resolved soil occluded mineral P
     real(r8), pointer :: primp_vr_col                 (:,:)       ! col (gP/m3) vertically-resolved soil parimary mineral P
     real(r8), pointer :: sminp_vr_col                 (:,:)       ! col (gP/m3) vertically-resolved soil parimary mineral P
     real(r8), pointer :: ptrunc_vr_col                (:,:)       ! col (gP/m3) vertically-resolved column-level sink for P truncation

     ! wood product pools, for dynamic landcover
     real(r8), pointer :: seedp_col                    (:)     ! col (gP/m2) column-level pool for seeding new Patches
     real(r8), pointer :: prod1p_col                   (:)     ! col (gN/m2) crop product N pool, 1-year lifespan
     real(r8), pointer :: prod10p_col                  (:)     ! col (gP/m2) wood product P pool, 10-year lifespan
     real(r8), pointer :: prod100p_col                 (:)     ! col (gP/m2) wood product P pool, 100-year lifespan
     real(r8), pointer :: totprodp_col                 (:)     ! col (gP/m2) total wood product P

     ! summary (diagnostic) state variables, not involved in mass balance
     real(r8), pointer :: dispvegp_patch               (:)     ! patch (gP/m2) displayed veg phosphorus, excluding storage
     real(r8), pointer :: storvegp_patch               (:)     ! patch (gP/m2) stored vegetation phosphorus
     real(r8), pointer :: totvegp_patch                (:)     ! patch (gP/m2) total vegetation phosphorus
     real(r8), pointer :: totpftp_patch                (:)     ! patch (gP/m2) total pft-level phosphorus
     real(r8), pointer :: decomp_ppools_col            (:,:)   ! col (gP/m2)  decomposing (litter, cwd, soil) P pools
     real(r8), pointer :: decomp_ppools_1m_col         (:,:)   ! col (gP/m2)  diagnostic: decomposing (litter, cwd, soil) P pools to 1 meter
     real(r8), pointer :: sminp_col                    (:)     ! col (gP/m2) soil mineral P
     real(r8), pointer :: solutionp_col                (:)         ! col (gP/m2) soil solution P
     real(r8), pointer :: labilep_col                  (:)         ! col (gP/m2) soil labile mineral P
     real(r8), pointer :: secondp_col                  (:)         ! col (gP/m2) soil secondary mineralP
     real(r8), pointer :: occlp_col                    (:)         ! col (gP/m2) soil occluded mineral P
     real(r8), pointer :: primp_col                    (:)         ! col (gP/m2) soil parimary mineral P
     real(r8), pointer :: ptrunc_col                   (:)     ! col (gP/m2) column-level sink for P truncation
     real(r8), pointer :: cwdp_col                     (:)     ! col (gP/m2) Diagnostic: coarse woody debris P
     real(r8), pointer :: totlitp_col                  (:)     ! col (gP/m2) total litter phosphorus
     real(r8), pointer :: totsomp_col                  (:)     ! col (gP/m2) total soil organic matter phosphorus
     real(r8), pointer :: totlitp_1m_col               (:)     ! col (gP/m2) total litter phosphorus to 1 meter
     real(r8), pointer :: totsomp_1m_col               (:)     ! col (gP/m2) total soil organic matter phosphorus to 1 meter
     real(r8), pointer :: totecosysp_col               (:)     ! col (gP/m2) total ecosystem phosphorus, incl veg 
     real(r8), pointer :: totcolp_col                  (:)     ! col (gP/m2) total column phosphorus, incl veg

     ! patch averaged to column variables 
     real(r8), pointer :: totvegp_col                  (:)     ! col (gP/m2) total vegetation phosphorus (p2c)
     real(r8), pointer :: totpftp_col                  (:)     ! col (gP/m2) total pft-level phosphorus (p2c)

     ! col balance checks
     real(r8), pointer :: begpb_patch                  (:)     ! patch phosphorus mass, beginning of time step (gP/m**2)
     real(r8), pointer :: endpb_patch                  (:)     ! patch phosphorus mass, end of time step (gP/m**2)
     real(r8), pointer :: errpb_patch                  (:)     ! patch phosphorus balance error for the timestep (gP/m**2)
     real(r8), pointer :: begpb_col                    (:)     ! col phosphorus mass, beginning of time step (gP/m**2)
     real(r8), pointer :: endpb_col                    (:)     ! col phosphorus mass, end of time step (gP/m**2)
     real(r8), pointer :: errpb_col                    (:)     ! colphosphorus balance error for the timestep (gP/m**2)

     real(r8), pointer :: actual_leafcp                (:)     ! dynamic leaf cp ratio
     real(r8), pointer :: actual_frootcp               (:)     ! dynamic fine root cp ratio
     real(r8), pointer :: actual_livewdcp              (:)     ! dynamic live wood cp ratio
     real(r8), pointer :: actual_deadwdcp              (:)     ! dynamic dead wood cp ratio
     real(r8), pointer :: actual_graincp               (:)     ! dynamic grain cp ratio
     
     ! debug
     real(r8), pointer :: totpftp_beg_col              (:)
     real(r8), pointer :: solutionp_beg_col            (:)
     real(r8), pointer :: labilep_beg_col              (:)
     real(r8), pointer :: secondp_beg_col              (:)
     real(r8), pointer :: totlitp_beg_col              (:)
     real(r8), pointer :: cwdp_beg_col                 (:)
     real(r8), pointer :: totsomp_beg_col              (:)
         
     real(r8), pointer :: totlitp_end_col              (:)
     real(r8), pointer :: totpftp_end_col              (:)
     real(r8), pointer :: solutionp_end_col            (:)
     real(r8), pointer :: labilep_end_col              (:)
     real(r8), pointer :: secondp_end_col              (:)
     real(r8), pointer :: cwdp_end_col                 (:)
     real(r8), pointer :: totsomp_end_col              (:)
     
   contains

     procedure , public  :: Init   
     procedure , public  :: Restart
     procedure , public  :: SetValues
     procedure , public  :: ZeroDWT
     procedure , public  :: Summary
     procedure , private :: InitAllocate
     procedure , private :: InitHistory  
     procedure , private :: InitCold     

  end type phosphorusstate_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds,                           &
       leafc_patch, leafc_storage_patch, frootc_patch, frootc_storage_patch, &
       deadstemc_patch, decomp_cpools_vr_col, decomp_cpools_col, decomp_cpools_1m_col)

    class(phosphorusstate_type)         :: this
    type(bounds_type) , intent(in)    :: bounds  
    real(r8)          , intent(in)    :: leafc_patch          (bounds%begp:)
    real(r8)          , intent(in)    :: leafc_storage_patch  (bounds%begp:)
    real(r8)          , intent(in)    :: frootc_patch         (bounds%begp:)
    real(r8)          , intent(in)    :: frootc_storage_patch (bounds%begp:)
    real(r8)          , intent(in)    :: deadstemc_patch      (bounds%begp:)
    real(r8)          , intent(in)    :: decomp_cpools_vr_col (bounds%begc:, 1:, 1:)
    real(r8)          , intent(in)    :: decomp_cpools_col    (bounds%begc:, 1:)
    real(r8)          , intent(in)    :: decomp_cpools_1m_col (bounds%begc:, 1:)

    call this%InitAllocate (bounds )

    call this%InitHistory (bounds)

    call this%InitCold ( bounds, leafc_patch, leafc_storage_patch, &
         frootc_patch, frootc_storage_patch, deadstemc_patch, &
         decomp_cpools_vr_col, decomp_cpools_col, decomp_cpools_1m_col)

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
    !------------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc

    allocate(this%grainp_patch             (begp:endp))                   ; this%grainp_patch             (:)   = nan     
    allocate(this%grainp_storage_patch     (begp:endp))                   ; this%grainp_storage_patch     (:)   = nan
    allocate(this%grainp_xfer_patch        (begp:endp))                   ; this%grainp_xfer_patch        (:)   = nan     
    allocate(this%leafp_patch              (begp:endp))                   ; this%leafp_patch              (:)   = nan
    allocate(this%leafp_storage_patch      (begp:endp))                   ; this%leafp_storage_patch      (:)   = nan     
    allocate(this%leafp_xfer_patch         (begp:endp))                   ; this%leafp_xfer_patch         (:)   = nan     
    allocate(this%frootp_patch             (begp:endp))                   ; this%frootp_patch             (:)   = nan
    allocate(this%frootp_storage_patch     (begp:endp))                   ; this%frootp_storage_patch     (:)   = nan     
    allocate(this%frootp_xfer_patch        (begp:endp))                   ; this%frootp_xfer_patch        (:)   = nan     
    allocate(this%livestemp_patch          (begp:endp))                   ; this%livestemp_patch          (:)   = nan
    allocate(this%livestemp_storage_patch  (begp:endp))                   ; this%livestemp_storage_patch  (:)   = nan
    allocate(this%livestemp_xfer_patch     (begp:endp))                   ; this%livestemp_xfer_patch     (:)   = nan
    allocate(this%deadstemp_patch          (begp:endp))                   ; this%deadstemp_patch          (:)   = nan
    allocate(this%deadstemp_storage_patch  (begp:endp))                   ; this%deadstemp_storage_patch  (:)   = nan
    allocate(this%deadstemp_xfer_patch     (begp:endp))                   ; this%deadstemp_xfer_patch     (:)   = nan
    allocate(this%livecrootp_patch         (begp:endp))                   ; this%livecrootp_patch         (:)   = nan
    allocate(this%livecrootp_storage_patch (begp:endp))                   ; this%livecrootp_storage_patch (:)   = nan
    allocate(this%livecrootp_xfer_patch    (begp:endp))                   ; this%livecrootp_xfer_patch    (:)   = nan
    allocate(this%deadcrootp_patch         (begp:endp))                   ; this%deadcrootp_patch         (:)   = nan
    allocate(this%deadcrootp_storage_patch (begp:endp))                   ; this%deadcrootp_storage_patch (:)   = nan
    allocate(this%deadcrootp_xfer_patch    (begp:endp))                   ; this%deadcrootp_xfer_patch    (:)   = nan
    allocate(this%retransp_patch           (begp:endp))                   ; this%retransp_patch           (:)   = nan
    allocate(this%ppool_patch              (begp:endp))                   ; this%ppool_patch              (:)   = nan
    allocate(this%ptrunc_patch             (begp:endp))                   ; this%ptrunc_patch             (:)   = nan
    allocate(this%dispvegp_patch           (begp:endp))                   ; this%dispvegp_patch           (:)   = nan
    allocate(this%storvegp_patch           (begp:endp))                   ; this%storvegp_patch           (:)   = nan
    allocate(this%totvegp_patch            (begp:endp))                   ; this%totvegp_patch            (:)   = nan
    allocate(this%totpftp_patch            (begp:endp))                   ; this%totpftp_patch            (:)   = nan

    allocate(this%ptrunc_vr_col            (begc:endc,1:nlevdecomp_full)) ; this%ptrunc_vr_col            (:,:) = nan
    
    allocate(this%solutionp_vr_col         (begc:endc,1:nlevdecomp_full)) ; this%solutionp_vr_col         (:,:) = nan  
    allocate(this%labilep_vr_col           (begc:endc,1:nlevdecomp_full)) ; this%labilep_vr_col           (:,:) = nan  
    allocate(this%secondp_vr_col           (begc:endc,1:nlevdecomp_full)) ; this%secondp_vr_col           (:,:) = nan  
    allocate(this%occlp_vr_col             (begc:endc,1:nlevdecomp_full)) ; this%occlp_vr_col             (:,:) = nan  
    allocate(this%primp_vr_col             (begc:endc,1:nlevdecomp_full)) ; this%primp_vr_col             (:,:) = nan  
    allocate(this%sminp_vr_col             (begc:endc,1:nlevdecomp_full)) ; this%sminp_vr_col             (:,:) = nan  

    allocate(this%solutionp_col            (begc:endc))                   ; this%solutionp_col            (:)   = nan
    allocate(this%labilep_col              (begc:endc))                   ; this%labilep_col              (:)   = nan
    allocate(this%secondp_col              (begc:endc))                   ; this%secondp_col              (:)   = nan
    allocate(this%occlp_col                (begc:endc))                   ; this%occlp_col                (:)   = nan
    allocate(this%primp_col                (begc:endc))                   ; this%primp_col                (:)   = nan
    allocate(this%cwdp_col                 (begc:endc))                   ; this%cwdp_col                 (:)   = nan
    allocate(this%sminp_col                (begc:endc))                   ; this%sminp_col                (:)   = nan
    allocate(this%ptrunc_col               (begc:endc))                   ; this%ptrunc_col               (:)   = nan
    allocate(this%seedp_col                (begc:endc))                   ; this%seedp_col                (:)   = nan
    allocate(this%prod1p_col               (begc:endc))                   ; this%prod1p_col               (:)   = nan
    allocate(this%prod10p_col              (begc:endc))                   ; this%prod10p_col              (:)   = nan
    allocate(this%prod100p_col             (begc:endc))                   ; this%prod100p_col             (:)   = nan
    allocate(this%totprodp_col             (begc:endc))                   ; this%totprodp_col             (:)   = nan
    allocate(this%totlitp_col              (begc:endc))                   ; this%totlitp_col              (:)   = nan
    allocate(this%totsomp_col              (begc:endc))                   ; this%totsomp_col              (:)   = nan
    allocate(this%totlitp_1m_col           (begc:endc))                   ; this%totlitp_1m_col           (:)   = nan
    allocate(this%totsomp_1m_col           (begc:endc))                   ; this%totsomp_1m_col           (:)   = nan
    allocate(this%totecosysp_col           (begc:endc))                   ; this%totecosysp_col           (:)   = nan
    allocate(this%totcolp_col              (begc:endc))                   ; this%totcolp_col              (:)   = nan
    allocate(this%decomp_ppools_col        (begc:endc,1:ndecomp_pools))   ; this%decomp_ppools_col        (:,:) = nan
    allocate(this%decomp_ppools_1m_col     (begc:endc,1:ndecomp_pools))   ; this%decomp_ppools_1m_col     (:,:) = nan
    allocate(this%totpftp_col              (begc:endc))                   ; this%totpftp_col              (:)   = nan
    allocate(this%totvegp_col              (begc:endc))                   ; this%totvegp_col              (:)   = nan

    allocate(this%decomp_ppools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools));
    this%decomp_ppools_vr_col(:,:,:)= nan

    allocate(this%begpb_patch (begp:endp));     this%begpb_patch (:) =nan
    allocate(this%begpb_col   (begc:endc));     this%begpb_col   (:) =nan
    allocate(this%endpb_patch (begp:endp));     this%endpb_patch (:) =nan
    allocate(this%endpb_col   (begc:endc));     this%endpb_col   (:) =nan
    allocate(this%errpb_patch (begp:endp));     this%errpb_patch (:) =nan
    allocate(this%errpb_col   (begc:endc));     this%errpb_col   (:) =nan 

    allocate(this%actual_leafcp       (begp:endp)); this%actual_leafcp       (:) = nan
    allocate(this%actual_frootcp      (begp:endp)); this%actual_frootcp      (:) = nan
    allocate(this%actual_livewdcp     (begp:endp)); this%actual_livewdcp     (:) = nan
    allocate(this%actual_deadwdcp     (begp:endp)); this%actual_deadwdcp     (:) = nan
    allocate(this%actual_graincp      (begp:endp)); this%actual_graincp      (:) = nan
    
    ! debug
    allocate(this%totpftp_beg_col    (begc:endc)); this%totpftp_beg_col      (:) = nan
    allocate(this%solutionp_beg_col  (begc:endc)); this%solutionp_beg_col    (:) = nan
    allocate(this%labilep_beg_col    (begc:endc)); this%labilep_beg_col      (:) = nan
    allocate(this%secondp_beg_col    (begc:endc)); this%secondp_beg_col      (:) = nan
    allocate(this%totlitp_beg_col    (begc:endc)); this%totlitp_beg_col      (:) = nan
    allocate(this%cwdp_beg_col       (begc:endc)); this%cwdp_beg_col         (:) = nan
    allocate(this%totsomp_beg_col    (begc:endc)); this%totsomp_beg_col      (:) = nan
         
    allocate(this%totlitp_end_col    (begc:endc)); this%totlitp_end_col      (:) = nan
    allocate(this%totpftp_end_col    (begc:endc)); this%totpftp_end_col      (:) = nan
    allocate(this%labilep_end_col    (begc:endc)); this%labilep_end_col      (:) = nan
    allocate(this%secondp_end_col    (begc:endc)); this%secondp_end_col      (:) = nan
    allocate(this%solutionp_end_col  (begc:endc)); this%solutionp_end_col    (:) = nan
    allocate(this%cwdp_end_col       (begc:endc)); this%cwdp_end_col         (:) = nan
    allocate(this%totsomp_end_col    (begc:endc)); this%totsomp_end_col      (:) = nan
    
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
       this%grainp_patch(begp:endp) = spval
       call hist_addfld1d (fname='GRAINP', units='gP/m^2', &
            avgflag='A', long_name='grain P', &
            ptr_patch=this%grainp_patch, default='inactive')
    end if

    this%leafp_patch(begp:endp) = spval
    call hist_addfld1d (fname='LEAFP', units='gP/m^2', &
         avgflag='A', long_name='leaf P', &
         ptr_patch=this%leafp_patch)

    this%leafp_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='LEAFP_STORAGE', units='gP/m^2', &
         avgflag='A', long_name='leaf P storage', &
         ptr_patch=this%leafp_storage_patch, default='inactive')

    this%leafp_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='LEAFP_XFER', units='gP/m^2', &
         avgflag='A', long_name='leaf P transfer', &
         ptr_patch=this%leafp_xfer_patch, default='inactive')

    this%frootp_patch(begp:endp) = spval
    call hist_addfld1d (fname='FROOTP', units='gP/m^2', &
         avgflag='A', long_name='fine root P', &
         ptr_patch=this%frootp_patch)

    this%frootp_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='FROOTP_STORAGE', units='gP/m^2', &
         avgflag='A', long_name='fine root P storage', &
         ptr_patch=this%frootp_storage_patch, default='inactive')

    this%frootp_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='FROOTP_XFER', units='gP/m^2', &
         avgflag='A', long_name='fine root P transfer', &
         ptr_patch=this%frootp_xfer_patch, default='inactive')

    this%livestemp_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMP', units='gP/m^2', &
         avgflag='A', long_name='live stem P', &
         ptr_patch=this%livestemp_patch)

    this%livestemp_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMP_STORAGE', units='gP/m^2', &
         avgflag='A', long_name='live stem P storage', &
         ptr_patch=this%livestemp_storage_patch, default='inactive')

    this%livestemp_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMP_XFER', units='gP/m^2', &
         avgflag='A', long_name='live stem P transfer', &
         ptr_patch=this%livestemp_xfer_patch, default='inactive')

    this%deadstemp_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMP', units='gP/m^2', &
         avgflag='A', long_name='dead stem P', &
         ptr_patch=this%deadstemp_patch)

    this%deadstemp_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMP_STORAGE', units='gP/m^2', &
         avgflag='A', long_name='dead stem P storage', &
         ptr_patch=this%deadstemp_storage_patch, default='inactive')

    this%deadstemp_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMP_XFER', units='gP/m^2', &
         avgflag='A', long_name='dead stem P transfer', &
         ptr_patch=this%deadstemp_xfer_patch, default='inactive')

    this%livecrootp_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTP', units='gP/m^2', &
         avgflag='A', long_name='live coarse root P', &
         ptr_patch=this%livecrootp_patch)

    this%livecrootp_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTP_STORAGE', units='gP/m^2', &
         avgflag='A', long_name='live coarse root P storage', &
         ptr_patch=this%livecrootp_storage_patch, default='inactive')

    this%livecrootp_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTP_XFER', units='gP/m^2', &
         avgflag='A', long_name='live coarse root P transfer', &
         ptr_patch=this%livecrootp_xfer_patch, default='inactive')

    this%deadcrootp_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTP', units='gP/m^2', &
         avgflag='A', long_name='dead coarse root P', &
         ptr_patch=this%deadcrootp_patch)

    this%deadcrootp_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTP_STORAGE', units='gP/m^2', &
         avgflag='A', long_name='dead coarse root P storage', &
         ptr_patch=this%deadcrootp_storage_patch, default='inactive')

    this%deadcrootp_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTP_XFER', units='gP/m^2', &
         avgflag='A', long_name='dead coarse root P transfer', &
         ptr_patch=this%deadcrootp_xfer_patch, default='inactive')

    this%retransp_patch(begp:endp) = spval
    call hist_addfld1d (fname='RETRANSP', units='gP/m^2', &
         avgflag='A', long_name='plant pool of retranslocated P', &
         ptr_patch=this%retransp_patch)

    this%ppool_patch(begp:endp) = spval
    call hist_addfld1d (fname='PPOOL', units='gP/m^2', &
         avgflag='A', long_name='temporary plant P pool', &
         ptr_patch=this%ppool_patch, default='inactive')

    this%ptrunc_patch(begp:endp) = spval
    call hist_addfld1d (fname='PFT_PTRUNC', units='gP/m^2', &
         avgflag='A', long_name='pft-level sink for P truncation', &
         ptr_patch=this%ptrunc_patch, default='inactive')

    this%dispvegp_patch(begp:endp) = spval
    call hist_addfld1d (fname='DISPVEGP', units='gP/m^2', &
         avgflag='A', long_name='displayed vegetation phosphorus', &
         ptr_patch=this%dispvegp_patch)

    this%storvegp_patch(begp:endp) = spval
    call hist_addfld1d (fname='STORVEGP', units='gP/m^2', &
         avgflag='A', long_name='stored vegetation phosphorus', &
         ptr_patch=this%storvegp_patch)

    this%totvegp_patch(begp:endp) = spval
    call hist_addfld1d (fname='TOTVEGP', units='gP/m^2', &
         avgflag='A', long_name='total vegetation phosphorus', &
         ptr_patch=this%totvegp_patch)

    this%totpftp_patch(begp:endp) = spval
    call hist_addfld1d (fname='TOTPFTP', units='gP/m^2', &
         avgflag='A', long_name='total PFT-level phosphorus', &
         ptr_patch=this%totpftp_patch)

    call hist_addfld1d (fname='actual_leafcp', units='gC/gP', &
         avgflag='A', long_name='flexible leafCP', &
         ptr_patch=this%actual_leafcp)
    call hist_addfld1d (fname='actual_frootcp', units='gC/gP', &
         avgflag='A', long_name='flexible frootCP', &
         ptr_patch=this%actual_frootcp)
    call hist_addfld1d (fname='actual_livewdcp', units='gC/gP', &
         avgflag='A', long_name='flexible livewdCP', &
         ptr_patch=this%actual_livewdcp)
    call hist_addfld1d (fname='actual_deadwdcp', units='gC/gP', &
         avgflag='A', long_name='flexible deadwdCP', &
         ptr_patch=this%actual_deadwdcp)
    call hist_addfld1d (fname='actual_graincp', units='gC/gP', &
         avgflag='A', long_name='flexible grainCP', &
         ptr_patch=this%actual_graincp)
         
    !-------------------------------
    ! P state variables - native to column
    !-------------------------------

    if ( nlevdecomp_full > 1 ) then
       this%decomp_ppools_vr_col(begc:endc,:,:) = spval
       this%decomp_ppools_1m_col(begc:endc,:) = spval
    end if
    this%decomp_ppools_col(begc:endc,:) = spval
    do l  = 1, ndecomp_pools
       if ( nlevdecomp_full > 1 ) then
          data2dptr => this%decomp_ppools_vr_col(:,:,l)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'P_vr'
          longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' P (vertically resolved)'
          call hist_addfld2d (fname=fieldname, units='gP/m^3',  type2d='levdcmp', &
               avgflag='A', long_name=longname, &
               ptr_col=data2dptr)
       endif

       data1dptr => this%decomp_ppools_col(:,l)
       fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'P'
       longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' P'
       call hist_addfld1d (fname=fieldname, units='gP/m^2', &
            avgflag='A', long_name=longname, &
            ptr_col=data1dptr)

       if ( nlevdecomp_full > 1 ) then
          data1dptr => this%decomp_ppools_1m_col(:,l)
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

       this%totlitp_1m_col(begc:endc) = spval
       call hist_addfld1d (fname='TOTLITP_1m', units='gP/m^2', &
            avgflag='A', long_name='total litter P to 1 meter', &
            ptr_col=this%totlitp_1m_col)

       this%totsomp_1m_col(begc:endc) = spval
       call hist_addfld1d (fname='TOTSOMP_1m', units='gP/m^2', &
            avgflag='A', long_name='total soil organic matter P to 1 meter', &
            ptr_col=this%totsomp_1m_col)
    endif

    this%ptrunc_col(begc:endc) = spval
    call hist_addfld1d (fname='COL_PTRUNC', units='gP/m^2',  &
         avgflag='A', long_name='column-level sink for P truncation', &
         ptr_col=this%ptrunc_col)

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

    this%totlitp_col(begc:endc) = spval
    call hist_addfld1d (fname='TOTLITP', units='gP/m^2', &
         avgflag='A', long_name='total litter P', &
         ptr_col=this%totlitp_col)

    this%totsomp_col(begc:endc) = spval
    call hist_addfld1d (fname='TOTSOMP', units='gP/m^2', &
         avgflag='A', long_name='total soil organic matter P', &
         ptr_col=this%totsomp_col)

    this%totecosysp_col(begc:endc) = spval
    call hist_addfld1d (fname='TOTECOSYSP', units='gP/m^2', &
         avgflag='A', long_name='total ecosystem P', &
         ptr_col=this%totecosysp_col)

    this%totcolp_col(begc:endc) = spval
    call hist_addfld1d (fname='TOTCOLP', units='gP/m^2', &
         avgflag='A', long_name='total column-level P', &
         ptr_col=this%totcolp_col)

    this%seedp_col(begc:endc) = spval
    call hist_addfld1d (fname='SEEDP', units='gP/m^2', &
         avgflag='A', long_name='P pool for seeding new PFTs ', &
         ptr_col=this%seedp_col, default='inactive')

    this%prod10p_col(begc:endc) = spval
    call hist_addfld1d (fname='PROD10P', units='gP/m^2', &
         avgflag='A', long_name='10-yr wood product P', &
         ptr_col=this%prod10p_col, default='inactive')

    this%prod100p_col(begc:endc) = spval
    call hist_addfld1d (fname='PROD100P', units='gP/m^2', &
         avgflag='A', long_name='100-yr wood product P', &
         ptr_col=this%prod100p_col, default='inactive')

    this%prod1p_col(begc:endc) = spval
    call hist_addfld1d (fname='PROD1P', units='gP/m^2', &
         avgflag='A', long_name='1-yr crop product P', &
         ptr_col=this%prod1p_col, default='inactive')

    this%totprodp_col(begc:endc) = spval
    call hist_addfld1d (fname='TOTPRODP', units='gP/m^2', &
         avgflag='A', long_name='total wood product P', &
         ptr_col=this%totprodp_col, default='inactive')

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, &
       leafc_patch, leafc_storage_patch, frootc_patch, frootc_storage_patch, &
       deadstemc_patch, decomp_cpools_vr_col, decomp_cpools_col, decomp_cpools_1m_col)
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
    real(r8)          , intent(in) :: decomp_cpools_1m_col(bounds%begc:,:)
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
    SHR_ASSERT_ALL((ubound(decomp_cpools_1m_col) == (/bounds%endc,ndecomp_pools/)),                 errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(decomp_cpools_vr_col) == (/bounds%endc,nlevdecomp_full,ndecomp_pools/)), errMsg(__FILE__, __LINE__))

    ! Set column filters

    num_special_patch = 0
    do p = bounds%begp,bounds%endp
       l = pft%landunit(p)
       if (lun%ifspecial(l)) then
          num_special_patch = num_special_patch + 1
          special_patch(num_special_patch) = p
       end if
    end do

    ! Set patch filters

    num_special_col = 0
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%ifspecial(l)) then
          num_special_col = num_special_col + 1
          special_col(num_special_col) = c
       end if
    end do

    !-------------------------------------------
    ! initialize pft-level variables
    !-------------------------------------------

    do p = bounds%begp,bounds%endp

       l = pft%landunit(p)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then

          if (pft%itype(p) == noveg) then
             this%leafp_patch(p) = 0._r8
             this%leafp_storage_patch(p) = 0._r8
          else
             this%leafp_patch(p)         = leafc_patch(p)         / ecophyscon%leafcp(pft%itype(p))
             this%leafp_storage_patch(p) = leafc_storage_patch(p) / ecophyscon%leafcp(pft%itype(p))
          end if

          this%leafp_xfer_patch(p)        = 0._r8
          if ( crop_prog )then
             this%grainp_patch(p)         = 0._r8
             this%grainp_storage_patch(p) = 0._r8
             this%grainp_xfer_patch(p)    = 0._r8
          end if
          this%frootp_patch(p)            = 0._r8
          this%frootp_storage_patch(p)    = 0._r8
          this%frootp_xfer_patch(p)       = 0._r8
          this%livestemp_patch(p)         = 0._r8
          this%livestemp_storage_patch(p) = 0._r8
          this%livestemp_xfer_patch(p)    = 0._r8

          ! tree types need to be initialized with some stem mass so that
          ! roughness length is not zero in canopy flux calculation

          if (ecophyscon%woody(pft%itype(p)) == 1._r8) then
             this%deadstemp_patch(p) = deadstemc_patch(p) / ecophyscon%deadwdcp(pft%itype(p))
          else
             this%deadstemp_patch(p) = 0._r8
          end if
          
          if (nu_com .ne. 'RD') then
              ! ECA competition calculate root NP uptake as a function of fine root biomass
              ! better to initialize root CNP pools with a non-zero value
              if (pft%itype(p) .ne. noveg) then
                 this%frootp_patch(p) = frootc_patch(p) / ecophyscon%frootcp(pft%itype(p))
                 this%frootp_storage_patch(p) = frootc_storage_patch(p) / ecophyscon%frootcp(pft%itype(p))
              end if
          end if
           
          this%deadstemp_storage_patch(p)  = 0._r8
          this%deadstemp_xfer_patch(p)     = 0._r8
          this%livecrootp_patch(p)         = 0._r8
          this%livecrootp_storage_patch(p) = 0._r8
          this%livecrootp_xfer_patch(p)    = 0._r8
          this%deadcrootp_patch(p)         = 0._r8
          this%deadcrootp_storage_patch(p) = 0._r8
          this%deadcrootp_xfer_patch(p)    = 0._r8
          this%retransp_patch(p)           = 0._r8
          this%ppool_patch(p)              = 0._r8
          this%ptrunc_patch(p)             = 0._r8
          this%dispvegp_patch(p)           = 0._r8
          this%storvegp_patch(p)           = 0._r8
          this%totvegp_patch(p)            = 0._r8
          this%totpftp_patch(p)            = 0._r8
       end if
       
       this%actual_leafcp(p)       = ecophyscon%leafcp(pft%itype(p))
       this%actual_frootcp(p)      = ecophyscon%frootcp(pft%itype(p))
       this%actual_livewdcp(p)     = ecophyscon%livewdcp(pft%itype(p))
       this%actual_deadwdcp(p)     = ecophyscon%deadwdcp(pft%itype(p))
       this%actual_graincp(p)      = ecophyscon%graincp(pft%itype(p))
       
    end do

    !-------------------------------------------
    ! initialize column-level variables
    !-------------------------------------------

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then

          ! column phosphorus state variables
          this%ptrunc_col(c) = 0._r8
          this%sminp_col(c) = 0._r8
          do j = 1, nlevdecomp
             do k = 1, ndecomp_pools
                this%decomp_ppools_vr_col(c,j,k) = decomp_cpools_vr_col(c,j,k) / decomp_cascade_con%initial_cp_ratio(k)
             end do
             this%sminp_vr_col(c,j) = 0._r8
             this%ptrunc_vr_col(c,j) = 0._r8
          end do
          if ( nlevdecomp > 1 ) then
             do j = nlevdecomp+1, nlevdecomp_full
                do k = 1, ndecomp_pools
                   this%decomp_ppools_vr_col(c,j,k) = 0._r8
                end do
                this%sminp_vr_col(c,j) = 0._r8
                this%ptrunc_vr_col(c,j) = 0._r8
             end do
          end if
          do k = 1, ndecomp_pools
             this%decomp_ppools_col(c,k)    = decomp_cpools_col(c,k)    / decomp_cascade_con%initial_cp_ratio(k)
             this%decomp_ppools_1m_col(c,k) = decomp_cpools_1m_col(c,k) / decomp_cascade_con%initial_cp_ratio(k)
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

          this%totlitp_col(c)    = 0._r8
          this%totsomp_col(c)    = 0._r8
          this%totlitp_1m_col(c) = 0._r8
          this%totsomp_1m_col(c) = 0._r8
          this%totecosysp_col(c) = 0._r8
          this%totcolp_col(c)    = 0._r8
          this%cwdp_col(c)       = 0._r8

          ! dynamic landcover state variables
          this%seedp_col(c)         = 0._r8
          this%prod1p_col(c)        = 0._r8
          this%prod10p_col(c)       = 0._r8
          this%prod100p_col(c)      = 0._r8
          this%totprodp_col(c)      = 0._r8
       end if
    end do

    ! now loop through special filters and explicitly set the variables that
    ! have to be in place for biogeophysics

    do fc = 1,num_special_col
       c = special_col(fc)

       this%seedp_col(c)    = 0._r8
       this%prod1p_col(c)   = 0._r8
       this%prod10p_col(c)  = 0._r8	  
       this%prod100p_col(c) = 0._r8	  
       this%totprodp_col(c) = 0._r8	  
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
    integer            :: i,j,k,l,c
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
    !------------------------------------------------------------------------

    !--------------------------------
    ! patch phosphorus state variables
    !--------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='leafp', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%leafp_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='leafp_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%leafp_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='leafp_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%leafp_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='frootp', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%frootp_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='frootp_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%frootp_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='frootp_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%frootp_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livestemp', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livestemp_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livestemp_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livestemp_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livestemp_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livestemp_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadstemp', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadstemp_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadstemp_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadstemp_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadstemp_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadstemp_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livecrootp', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livecrootp_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livecrootp_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livecrootp_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livecrootp_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livecrootp_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadcrootp', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadcrootp_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadcrootp_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadcrootp_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadcrootp_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadcrootp_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='retransp', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%retransp_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='ppool', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%ppool_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='pft_ptrunc', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%ptrunc_patch) 

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag,  varname='grainp', xtype=ncd_double,  &
            dim1name='pft',    long_name='grain P', units='gP/m2', &
            interpinic_flag='interp', readvar=readvar, data=this%grainp_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='grainp_storage', xtype=ncd_double,  &
            dim1name='pft',    long_name='grain P storage', units='gP/m2', &
            interpinic_flag='interp', readvar=readvar, data=this%grainp_storage_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='grainp_xfer', xtype=ncd_double,  &
            dim1name='pft',    long_name='grain P transfer', units='gP/m2', &
            interpinic_flag='interp', readvar=readvar, data=this%grainp_xfer_patch)
    end if

    call restartvar(ncid=ncid, flag=flag,  varname='actual_leafcp', xtype=ncd_double,  &
        dim1name='pft',    long_name='flexible leafCP', units='gC/gP', &
        interpinic_flag='interp', readvar=readvar, data=this%actual_leafcp)
    call restartvar(ncid=ncid, flag=flag,  varname='actual_frootcp', xtype=ncd_double,  &
        dim1name='pft',    long_name='flexible frootCP', units='gC/gP', &
        interpinic_flag='interp', readvar=readvar, data=this%actual_frootcp)
    call restartvar(ncid=ncid, flag=flag,  varname='actual_livewdcp', xtype=ncd_double,  &
        dim1name='pft',    long_name='flexible livewdCP', units='gC/gP', &
        interpinic_flag='interp', readvar=readvar, data=this%actual_livewdcp)
    call restartvar(ncid=ncid, flag=flag,  varname='actual_deadwdcp', xtype=ncd_double,  &
        dim1name='pft',    long_name='flexible deadwdCP', units='gC/gP', &
        interpinic_flag='interp', readvar=readvar, data=this%actual_deadwdcp)
    call restartvar(ncid=ncid, flag=flag,  varname='actual_graincp', xtype=ncd_double,  &
        dim1name='pft',    long_name='flexible grainCP', units='gC/gP', &
        interpinic_flag='interp', readvar=readvar, data=this%actual_graincp)
    
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
          ptr2d => this%decomp_ppools_vr_col(:,:,k)
          call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vr", xtype=ncd_double, &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=ptr2d) 
       else
          ptr1d => this%decomp_ppools_vr_col(:,1,k)
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
       ptr2d => this%ptrunc_vr_col
       call restartvar(ncid=ncid, flag=flag, varname="col_ptrunc_vr", xtype=ncd_double,  &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)
    else
       ptr1d => this%ptrunc_vr_col(:,1)
       call restartvar(ncid=ncid, flag=flag, varname="col_ptrunc", xtype=ncd_double,  &
            dim1name='column', &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp' , readvar=readvar, data=ptr1d)
    end if

    !!! Debug balance 
    call restartvar(ncid=ncid, flag=flag, varname='totsomp', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%totsomp_col) 
    call restartvar(ncid=ncid, flag=flag, varname='cwdp', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%cwdp_col) 
    call restartvar(ncid=ncid, flag=flag, varname='totlitp', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%totlitp_col) 

    call restartvar(ncid=ncid, flag=flag, varname='totcolp', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%totcolp_col) 

    call restartvar(ncid=ncid, flag=flag, varname='seedp', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%seedp_col) 

    call restartvar(ncid=ncid, flag=flag, varname='prod10p', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%prod10p_col) 

    call restartvar(ncid=ncid, flag=flag, varname='prod100p', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%prod100p_col) 

    call restartvar(ncid=ncid, flag=flag, varname='prod1p', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%prod1p_col)

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
                 if (decomp_cascade_con%spinup_factor(k) > 1) m = m  / cnstate_vars%scalaravg_col(c)
               else if ( enter_spinup ) then 
                 m = 1. / decomp_cascade_con%spinup_factor(k)
		 if (decomp_cascade_con%spinup_factor(k) > 1) m = m  * cnstate_vars%scalaravg_col(c)
               end if 
               this%decomp_ppools_vr_col(c,j,k) = this%decomp_ppools_vr_col(c,j,k) * m
             end do
          end do
       end do

       do i = bounds%begp, bounds%endp
          if (exit_spinup) then
             m_veg = spinup_mortality_factor
          else if (enter_spinup) then
             m_veg = 1._r8 / spinup_mortality_factor
          end if
          this%deadstemp_patch(i)  = this%deadstemp_patch(i) * m_veg
          this%deadcrootp_patch(i) = this%deadcrootp_patch(i) * m_veg
       end do

    end if

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

       this%leafp_patch(i)              = value_patch
       this%leafp_storage_patch(i)      = value_patch
       this%leafp_xfer_patch(i)         = value_patch
       this%frootp_patch(i)             = value_patch
       this%frootp_storage_patch(i)     = value_patch
       this%frootp_xfer_patch(i)        = value_patch
       this%livestemp_patch(i)          = value_patch
       this%livestemp_storage_patch(i)  = value_patch
       this%livestemp_xfer_patch(i)     = value_patch
       this%deadstemp_patch(i)          = value_patch
       this%deadstemp_storage_patch(i)  = value_patch
       this%deadstemp_xfer_patch(i)     = value_patch
       this%livecrootp_patch(i)         = value_patch
       this%livecrootp_storage_patch(i) = value_patch
       this%livecrootp_xfer_patch(i)    = value_patch
       this%deadcrootp_patch(i)         = value_patch
       this%deadcrootp_storage_patch(i) = value_patch
       this%deadcrootp_xfer_patch(i)    = value_patch
       this%retransp_patch(i)           = value_patch
       this%ppool_patch(i)              = value_patch
       this%ptrunc_patch(i)             = value_patch
       this%dispvegp_patch(i)           = value_patch
       this%storvegp_patch(i)           = value_patch
       this%totvegp_patch(i)            = value_patch
       this%totpftp_patch(i)            = value_patch
       
       this%actual_leafcp(i)            = value_patch  
       this%actual_frootcp(i)           = value_patch  
       this%actual_livewdcp(i)          = value_patch  
       this%actual_deadwdcp(i)          = value_patch  
    end do

    if ( crop_prog )then
       do fi = 1,num_patch
          i = filter_patch(fi)
          this%grainp_patch(i)          = value_patch
          this%grainp_storage_patch(i)  = value_patch
          this%grainp_xfer_patch(i)     = value_patch   
          this%actual_graincp(i)        = value_patch   
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
       this%ptrunc_col(i)      = value_column
       this%cwdp_col(i)        = value_column
       this%totlitp_col(i)     = value_column
       this%totsomp_col(i)     = value_column
       this%totecosysp_col(i)  = value_column
       this%totcolp_col(i)     = value_column
       this%totsomp_1m_col(i)  = value_column
       this%totlitp_1m_col(i)  = value_column
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
          this%ptrunc_vr_col(i,j)      = value_column
       end do
    end do

    ! column and decomp_pools
    do k = 1, ndecomp_pools
       do fi = 1,num_column
          i = filter_column(fi)
          this%decomp_ppools_col(i,k)    = value_column
          this%decomp_ppools_1m_col(i,k) = value_column
       end do
    end do

    ! column levdecomp, and decomp_pools
    do j = 1,nlevdecomp_full
       do k = 1, ndecomp_pools
          do fi = 1,num_column
             i = filter_column(fi)
             this%decomp_ppools_vr_col(i,j,k) = value_column
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
       this%dispvegp_patch(p) = 0._r8
       this%storvegp_patch(p) = 0._r8
       this%totvegp_patch(p)  = 0._r8
       this%totpftp_patch(p)  = 0._r8
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
    !-----------------------------------------------------------------------

    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! displayed vegetation phosphorus, excluding storage (DISPVEGN)
       this%dispvegp_patch(p) = &
            this%leafp_patch(p)      + &
            this%frootp_patch(p)     + &
            this%livestemp_patch(p)  + &
            this%deadstemp_patch(p)  + &
            this%livecrootp_patch(p) + &
            this%deadcrootp_patch(p)
       
      ! stored vegetation phosphorus, including retranslocated N pool (STORVEGN)
      this%storvegp_patch(p) = &
           this%leafp_storage_patch(p)      + &
           this%frootp_storage_patch(p)     + &
           this%livestemp_storage_patch(p)  + &
           this%deadstemp_storage_patch(p)  + &
           this%livecrootp_storage_patch(p) + &
           this%deadcrootp_storage_patch(p) + &
           this%leafp_xfer_patch(p)         + &
           this%frootp_xfer_patch(p)        + &
           this%livestemp_xfer_patch(p)     + &
           this%deadstemp_xfer_patch(p)     + &
           this%livecrootp_xfer_patch(p)    + &
           this%deadcrootp_xfer_patch(p)    + &
           this%ppool_patch(p)              + &
           this%retransp_patch(p)

      if ( crop_prog .and. pft%itype(p) >= npcropmin )then
         this%dispvegp_patch(p) = &
              this%dispvegp_patch(p) + &
              this%grainp_patch(p)

         this%storvegp_patch(p) = &
              this%storvegp_patch(p) + &
              this%grainp_storage_patch(p)     + &
              this%grainp_xfer_patch(p)
      end if

      ! total vegetation phosphorus (TOTVEGN)
      this%totvegp_patch(p) = &
           this%dispvegp_patch(p) + &
           this%storvegp_patch(p)

      ! total pft-level carbon (add pft_ntrunc)
      this%totpftp_patch(p) = &
           this%totvegp_patch(p) + &
           this%ptrunc_patch(p)

   end do

   call p2c(bounds, num_soilc, filter_soilc, &
        this%totvegp_patch(bounds%begp:bounds%endp), &
        this%totvegp_col(bounds%begc:bounds%endc))

   call p2c(bounds, num_soilc, filter_soilc, &
        this%totpftp_patch(bounds%begp:bounds%endp), &
        this%totpftp_col(bounds%begc:bounds%endc))

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
         this%decomp_ppools_col(c,l) = 0._r8
      end do
      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            this%decomp_ppools_col(c,l) = &
                 this%decomp_ppools_col(c,l) + &
                 this%decomp_ppools_vr_col(c,j,l) * dzsoi_decomp(j)
         end do
      end do
   end do

   ! for vertically-resolved soil biogeochemistry, calculate some diagnostics of carbon pools to a given depth
   if ( nlevdecomp > 1) then

      do l = 1, ndecomp_pools
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            this%decomp_ppools_1m_col(c,l) = 0._r8
         end do
      end do

      ! vertically integrate each of the decomposing n pools to 1 meter
      maxdepth = 1._r8
      do l = 1, ndecomp_pools
         do j = 1, nlevdecomp
            if ( zisoi(j) <= maxdepth ) then
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  this%decomp_ppools_1m_col(c,l) = &
                       this%decomp_ppools_1m_col(c,l) + &
                       this%decomp_ppools_vr_col(c,j,l) * dzsoi_decomp(j)
               end do
            elseif ( zisoi(j-1) < maxdepth ) then
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  this%decomp_ppools_1m_col(c,l) = &
                       this%decomp_ppools_1m_col(c,l) + &
                       this%decomp_ppools_vr_col(c,j,l) * (maxdepth - zisoi(j-1))
               end do
            endif
         end do
      end do
      
      ! total litter phosphorus to 1 meter (TOTLITN_1m)
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%totlitp_1m_col(c) = 0._r8
      end do
      do l = 1, ndecomp_pools
         if ( decomp_cascade_con%is_litter(l) ) then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               this%totlitp_1m_col(c) = &
                    this%totlitp_1m_col(c) + &
                    this%decomp_ppools_1m_col(c,l)
            end do
         end if
      end do
      
      ! total soil organic matter phosphorus to 1 meter (TOTSOMN_1m)
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%totsomp_1m_col(c) = 0._r8
      end do
      do l = 1, ndecomp_pools
         if ( decomp_cascade_con%is_soil(l) ) then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               this%totsomp_1m_col(c) = &
                    this%totsomp_1m_col(c) + &
                    this%decomp_ppools_1m_col(c,l)
            end do
         end if
      end do
      
   endif
   
   ! total litter phosphorus (TOTLITN)
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      this%totlitp_col(c)    = 0._r8
   end do
   do l = 1, ndecomp_pools
      if ( decomp_cascade_con%is_litter(l) ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            this%totlitp_col(c) = &
                 this%totlitp_col(c) + &
                 this%decomp_ppools_col(c,l)
         end do
      end if
   end do
   
   ! total soil organic matter phosphorus (TOTSOMN)
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      this%totsomp_col(c)    = 0._r8
   end do
   do l = 1, ndecomp_pools
      if ( decomp_cascade_con%is_soil(l) ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            this%totsomp_col(c) = &
                 this%totsomp_col(c) + &
                 this%decomp_ppools_col(c,l)
         end do
      end if
   end do
   
   ! total cwdn
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      this%cwdp_col(c) = 0._r8
   end do
   do l = 1, ndecomp_pools
      if ( decomp_cascade_con%is_cwd(l) ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            this%cwdp_col(c) = &
                 this%cwdp_col(c) + &
                 this%decomp_ppools_col(c,l)
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
      this%ptrunc_col(c) = 0._r8
   end do
   do j = 1, nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%ptrunc_col(c) = &
              this%ptrunc_col(c) + &
              this%ptrunc_vr_col(c,j) * dzsoi_decomp(j)
      end do
   end do

   do fc = 1,num_soilc
      c = filter_soilc(fc)

      ! total wood product phosphorus
      this%totprodp_col(c) = &
           this%prod1p_col(c) + &
           this%prod10p_col(c) + &
           this%prod100p_col(c)	 

      ! total ecosystem phosphorus, including veg (TOTECOSYSP)
      this%totecosysp_col(c) = &
           this%cwdp_col(c) + &
           this%totlitp_col(c) + &
           this%totsomp_col(c) + &
           this%solutionp_col(c) + &
           this%labilep_col(c) + &
           this%secondp_col(c) + &
           this%primp_col(c) + &
           this%occlp_col(c) + &
           this%totprodp_col(c) + &
           this%totvegp_col(c)

      ! total column phosphorus, including pft (TOTCOLP)
      this%totcolp_col(c) = &
           this%totpftp_col(c) + &
           this%cwdp_col(c) + &
           this%totlitp_col(c) + &
           this%totsomp_col(c) + &
           this%solutionp_col(c) + &
           this%labilep_col(c) + &
           this%secondp_col(c) + &
           this%totprodp_col(c) + &
           this%seedp_col(c)    + &
           this%ptrunc_col(c)
   end do

 end subroutine Summary

end module PhosphorusStateType
