module PhosphorusStateType

#include "shr_assert.h"

  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_infnan_mod         , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use elm_varpar             , only : ndecomp_cascade_transitions, ndecomp_pools, nlevcan
  use elm_varpar             , only : nlevdecomp_full, nlevdecomp, crop_prog
  use elm_varcon             , only : spval, ispval, dzsoi_decomp, zisoi, zsoi
  use landunit_varcon        , only : istcrop, istsoil 
  use elm_varctl             , only : use_nitrif_denitrif, use_vertsoilc, use_century_decomp
  use elm_varctl             , only : iulog, override_bgc_restart_mismatch_dump, spinup_state
  use decompMod              , only : bounds_type
  use pftvarcon              , only : npcropmin, nstor
  use CNDecompCascadeConType , only : decomp_cascade_con
  use VegetationPropertiesType         , only : veg_vp
  use abortutils             , only : endrun
  use spmdMod                , only : masterproc 
  use LandunitType           , only : lun_pp                
  use ColumnType             , only : col_pp                
  use VegetationType              , only : veg_pp
  use elm_varctl             , only : nu_com, use_crop
  ! soil phosphorus initialization Qing Z. 2017
  use pftvarcon              , only : VMAX_MINSURF_P_vr, KM_MINSURF_P_vr
  use soilorder_varcon       , only : smax, ks_sorption
  use dynPatchStateUpdaterMod      , only : patch_state_updater_type
  use SpeciesMod           , only : CN_SPECIES_P
  ! 
  ! !PUBLIC TYPES:
  implicit none
  save
  private

  real(r8) , parameter :: ppool_seed_param     = 0.01_r8

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
     real(r8), pointer :: plant_p_buffer_patch        (:)     ! patch (gP/m2) pft-level abstract p storage
     real(r8), pointer :: decomp_ppools_vr_col         (:,:,:)     ! col (gP/m3) vertically-resolved decomposing (litter, cwd, soil) P pools 
     real(r8), pointer :: solutionp_vr_col             (:,:)       ! col (gP/m3) vertically-resolved soil solution P
     real(r8), pointer :: labilep_vr_col               (:,:)       ! col (gP/m3) vertically-resolved soil labile mineral P
     real(r8), pointer :: secondp_vr_col               (:,:)       ! col (gP/m3) vertically-resolved soil secondary mineralP
     real(r8), pointer :: occlp_vr_col                 (:,:)       ! col (gP/m3) vertically-resolved soil occluded mineral P
     real(r8), pointer :: primp_vr_col                 (:,:)       ! col (gP/m3) vertically-resolved soil parimary mineral P
     real(r8), pointer :: sminp_vr_col                 (:,:)       ! col (gP/m3) vertically-resolved soil parimary mineral P
     real(r8), pointer :: ptrunc_vr_col                (:,:)       ! col (gP/m3) vertically-resolved column-level sink for P truncation

     ! wood product pools, for dynamic landcover
     real(r8), pointer :: cropseedp_deficit_patch      (:)     ! (gP/m2) pool for seeding new crop growth; this is a NEGATIVE term, indicating the amount of seed usage that needs to be repaid
     real(r8), pointer :: seedp_grc                    (:)     ! (gP/m2) gridcell-level pool for seeding new PFTs via dynamic landcover
     real(r8), pointer :: seedp_col                    (:)     ! col (gP/m2) column-level pool for seeding new Patches
     real(r8), pointer :: prod1p_col                   (:)     ! col (gN/m2) crop product N pool, 1-year lifespan
     real(r8), pointer :: prod10p_col                  (:)     ! col (gP/m2) wood product P pool, 10-year lifespan
     real(r8), pointer :: prod100p_col                 (:)     ! col (gP/m2) wood product P pool, 100-year lifespan
     real(r8), pointer :: totprodp_col                 (:)     ! col (gP/m2) total wood product P
     real(r8), pointer :: dyn_pbal_adjustments_col     (:)     ! (gP/m2) adjustments to each column made in this timestep via dynamic column area adjustments

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
     real(r8), pointer :: begpb_grc                    (:)     ! grid cell phosphorus mass, beginning of time step (gP/m**2)
     real(r8), pointer :: endpb_grc                    (:)     ! grid cell phosphorus mass, end of time step (gP/m**2)
     real(r8), pointer :: errpb_grc                    (:)     ! grid cell phosphorus balance error for the timestep (gP/m**2)

     real(r8), pointer :: solutionp_vr_col_cur         (:,:)
     real(r8), pointer :: solutionp_vr_col_prev        (:,:)
     real(r8), pointer :: labilep_vr_col_cur           (:,:)
     real(r8), pointer :: labilep_vr_col_prev          (:,:)
     real(r8), pointer :: secondp_vr_col_cur           (:,:)
     real(r8), pointer :: secondp_vr_col_prev          (:,:)
     real(r8), pointer :: occlp_vr_col_cur             (:,:)
     real(r8), pointer :: occlp_vr_col_prev            (:,:)
     real(r8), pointer :: primp_vr_col_cur             (:,:)
     real(r8), pointer :: primp_vr_col_prev            (:,:)

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
    integer           :: begg,endg
    !------------------------------------------------------------------------

    
  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! add history fields for all CN variables, always set as default='inactive'
    !
    ! !USES:
    use elm_varpar , only : ndecomp_cascade_transitions, ndecomp_pools
    use elm_varpar , only : nlevdecomp, nlevdecomp_full,crop_prog, nlevgrnd
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
    integer           :: begp,endp
    integer           :: begc,endc
    integer           :: begg,endg 
    character(24)     :: fieldname
    character(100)    :: longname
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg

    !-------------------------------
    ! P state variables - native to PFT
    !-------------------------------
    
    !-------------------------------
    ! P state variables - native to column
    !-------------------------------

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
    use elm_varpar     , only : crop_prog
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

    ! Set patch filters

    !-------------------------------------------
    ! initialize pft-level variables
    !-------------------------------------------


    !-------------------------------------------
    ! initialize column-level variables
    !-------------------------------------------


    ! initialize fields for special filters


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
    use elm_varctl          , only : spinup_mortality_factor
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
    use elm_varpar    , only: nlevdecomp,ndecomp_cascade_transitions,ndecomp_pools
    use elm_varctl    , only: use_nitrif_denitrif
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

   ! vertically integrate soil mineral P pools

 end subroutine Summary

end module PhosphorusStateType
