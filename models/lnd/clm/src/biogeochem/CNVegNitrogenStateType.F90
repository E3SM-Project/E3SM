module CNVegNitrogenStateType

#include "shr_assert.h"

  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use shr_infnan_mod                     , only : isnan => shr_infnan_isnan, nan => shr_infnan_nan, assignment(=)
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  use clm_varpar                         , only : ndecomp_cascade_transitions, ndecomp_pools, nlevcan
  use clm_varpar                         , only : nlevdecomp_full, nlevdecomp, crop_prog
  use clm_varcon                         , only : spval, ispval, dzsoi_decomp, zisoi
  use landunit_varcon                    , only : istcrop, istsoil 
  use clm_varctl                         , only : use_nitrif_denitrif, use_vertsoilc, use_century_decomp
  use clm_varctl                         , only : iulog, override_bgc_restart_mismatch_dump
  use decompMod                          , only : bounds_type
  use pftconMod                          , only : npcropmin, noveg, pftcon
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con
  use abortutils                         , only : endrun
  use spmdMod                            , only : masterproc 
  use LandunitType                       , only : lun                
  use ColumnType                         , only : col                
  use PatchType                          , only : patch                
  ! 
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  type, public :: cnveg_nitrogenstate_type

     real(r8), pointer :: grainn_patch             (:) ! (gN/m2) grain N (crop)
     real(r8), pointer :: grainn_storage_patch     (:) ! (gN/m2) grain N storage (crop)
     real(r8), pointer :: grainn_xfer_patch        (:) ! (gN/m2) grain N transfer (crop)
     real(r8), pointer :: leafn_patch              (:) ! (gN/m2) leaf N 
     real(r8), pointer :: leafn_storage_patch      (:) ! (gN/m2) leaf N storage
     real(r8), pointer :: leafn_xfer_patch         (:) ! (gN/m2) leaf N transfer
     real(r8), pointer :: frootn_patch             (:) ! (gN/m2) fine root N
     real(r8), pointer :: frootn_storage_patch     (:) ! (gN/m2) fine root N storage
     real(r8), pointer :: frootn_xfer_patch        (:) ! (gN/m2) fine root N transfer
     real(r8), pointer :: livestemn_patch          (:) ! (gN/m2) live stem N
     real(r8), pointer :: livestemn_storage_patch  (:) ! (gN/m2) live stem N storage
     real(r8), pointer :: livestemn_xfer_patch     (:) ! (gN/m2) live stem N transfer
     real(r8), pointer :: deadstemn_patch          (:) ! (gN/m2) dead stem N
     real(r8), pointer :: deadstemn_storage_patch  (:) ! (gN/m2) dead stem N storage
     real(r8), pointer :: deadstemn_xfer_patch     (:) ! (gN/m2) dead stem N transfer
     real(r8), pointer :: livecrootn_patch         (:) ! (gN/m2) live coarse root N
     real(r8), pointer :: livecrootn_storage_patch (:) ! (gN/m2) live coarse root N storage
     real(r8), pointer :: livecrootn_xfer_patch    (:) ! (gN/m2) live coarse root N transfer
     real(r8), pointer :: deadcrootn_patch         (:) ! (gN/m2) dead coarse root N
     real(r8), pointer :: deadcrootn_storage_patch (:) ! (gN/m2) dead coarse root N storage
     real(r8), pointer :: deadcrootn_xfer_patch    (:) ! (gN/m2) dead coarse root N transfer
     real(r8), pointer :: retransn_patch           (:) ! (gN/m2) plant pool of retranslocated N
     real(r8), pointer :: npool_patch              (:) ! (gN/m2) temporary plant N pool
     real(r8), pointer :: ntrunc_patch             (:) ! (gN/m2) patch-level sink for N truncation

     ! wood product pools, for dynamic landcover
     real(r8), pointer :: seedn_col                (:) ! (gN/m2) column-level pool for seeding new Patches
     real(r8), pointer :: prod10n_col              (:) ! (gN/m2) wood product N pool, 10-year lifespan
     real(r8), pointer :: prod100n_col             (:) ! (gN/m2) wood product N pool, 100-year lifespan
     real(r8), pointer :: totprodn_col             (:) ! (gN/m2) total wood product N

     ! summary (diagnostic) state variables, not involved in mass balance
     real(r8), pointer :: dispvegn_patch           (:) ! (gN/m2) displayed veg nitrogen, excluding storage
     real(r8), pointer :: storvegn_patch           (:) ! (gN/m2) stored vegetation nitrogen
     real(r8), pointer :: totvegn_patch            (:) ! (gN/m2) total vegetation nitrogen
     real(r8), pointer :: totvegn_col              (:) ! (gN/m2) total vegetation nitrogen (p2c)
     real(r8), pointer :: totn_patch               (:) ! (gN/m2) total patch-level nitrogen
     real(r8), pointer :: totn_col                 (:) ! (gN/m2) total column nitrogen, incl veg   
     real(r8), pointer :: totecosysn_col           (:) ! (gN/m2) total ecosystem nitrogen, incl veg  
     real(r8), pointer :: begnb_col                (:) ! (gN/m2) nitrogen mass, beginning of time step 
     real(r8), pointer :: endnb_col                (:) ! (gN/m2) nitrogen mass, end of time step 

   contains

     procedure , public  :: Init   
     procedure , public  :: Restart
     procedure , public  :: SetValues
     procedure , public  :: ZeroDWT
     procedure , public  :: Summary => Summary_nitrogenstate
     procedure , private :: InitAllocate 
     procedure , private :: InitHistory  
     procedure , private :: InitCold     

  end type cnveg_nitrogenstate_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds,                           &
       leafc_patch, leafc_storage_patch, deadstemc_patch)

    class(cnveg_nitrogenstate_type)   :: this
    type(bounds_type) , intent(in)    :: bounds  
    real(r8)          , intent(in)    :: leafc_patch          (bounds%begp:)
    real(r8)          , intent(in)    :: leafc_storage_patch  (bounds%begp:)
    real(r8)          , intent(in)    :: deadstemc_patch      (bounds%begp:)

    call this%InitAllocate (bounds )
    call this%InitHistory (bounds)
    call this%InitCold ( bounds, &
         leafc_patch, leafc_storage_patch, deadstemc_patch)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !ARGUMENTS:
    class (cnveg_nitrogenstate_type) :: this
    type(bounds_type) , intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begp,endp
    integer           :: begc,endc
    !------------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc

    allocate(this%grainn_patch             (begp:endp)) ; this%grainn_patch             (:) = nan     
    allocate(this%grainn_storage_patch     (begp:endp)) ; this%grainn_storage_patch     (:) = nan
    allocate(this%grainn_xfer_patch        (begp:endp)) ; this%grainn_xfer_patch        (:) = nan     
    allocate(this%leafn_patch              (begp:endp)) ; this%leafn_patch              (:) = nan
    allocate(this%leafn_storage_patch      (begp:endp)) ; this%leafn_storage_patch      (:) = nan     
    allocate(this%leafn_xfer_patch         (begp:endp)) ; this%leafn_xfer_patch         (:) = nan     
    allocate(this%frootn_patch             (begp:endp)) ; this%frootn_patch             (:) = nan
    allocate(this%frootn_storage_patch     (begp:endp)) ; this%frootn_storage_patch     (:) = nan     
    allocate(this%frootn_xfer_patch        (begp:endp)) ; this%frootn_xfer_patch        (:) = nan     
    allocate(this%livestemn_patch          (begp:endp)) ; this%livestemn_patch          (:) = nan
    allocate(this%livestemn_storage_patch  (begp:endp)) ; this%livestemn_storage_patch  (:) = nan
    allocate(this%livestemn_xfer_patch     (begp:endp)) ; this%livestemn_xfer_patch     (:) = nan
    allocate(this%deadstemn_patch          (begp:endp)) ; this%deadstemn_patch          (:) = nan
    allocate(this%deadstemn_storage_patch  (begp:endp)) ; this%deadstemn_storage_patch  (:) = nan
    allocate(this%deadstemn_xfer_patch     (begp:endp)) ; this%deadstemn_xfer_patch     (:) = nan
    allocate(this%livecrootn_patch         (begp:endp)) ; this%livecrootn_patch         (:) = nan
    allocate(this%livecrootn_storage_patch (begp:endp)) ; this%livecrootn_storage_patch (:) = nan
    allocate(this%livecrootn_xfer_patch    (begp:endp)) ; this%livecrootn_xfer_patch    (:) = nan
    allocate(this%deadcrootn_patch         (begp:endp)) ; this%deadcrootn_patch         (:) = nan
    allocate(this%deadcrootn_storage_patch (begp:endp)) ; this%deadcrootn_storage_patch (:) = nan
    allocate(this%deadcrootn_xfer_patch    (begp:endp)) ; this%deadcrootn_xfer_patch    (:) = nan
    allocate(this%retransn_patch           (begp:endp)) ; this%retransn_patch           (:) = nan
    allocate(this%npool_patch              (begp:endp)) ; this%npool_patch              (:) = nan
    allocate(this%ntrunc_patch             (begp:endp)) ; this%ntrunc_patch             (:) = nan
    allocate(this%dispvegn_patch           (begp:endp)) ; this%dispvegn_patch           (:) = nan
    allocate(this%storvegn_patch           (begp:endp)) ; this%storvegn_patch           (:) = nan
    allocate(this%totvegn_patch            (begp:endp)) ; this%totvegn_patch            (:) = nan
    allocate(this%totn_patch               (begp:endp)) ; this%totn_patch               (:) = nan

    allocate(this%seedn_col                (begc:endc)) ; this%seedn_col                (:) = nan
    allocate(this%prod10n_col              (begc:endc)) ; this%prod10n_col              (:) = nan
    allocate(this%prod100n_col             (begc:endc)) ; this%prod100n_col             (:) = nan
    allocate(this%totprodn_col             (begc:endc)) ; this%totprodn_col             (:) = nan
    allocate(this%totvegn_col              (begc:endc)) ; this%totvegn_col              (:) = nan
    allocate(this%totn_col                 (begc:endc)) ; this%totn_col                 (:) = nan
    allocate(this%totecosysn_col           (begc:endc)) ; this%totecosysn_col           (:) = nan
    allocate(this%begnb_col                (begc:endc)) ; this%begnb_col                (:) = nan
    allocate(this%endnb_col                (begc:endc)) ; this%endnb_col                (:) = nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! add history fields for all CN variables, always set as default='inactive'
    !
    ! !USES:
    use histFileMod, only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(cnveg_nitrogenstate_type) :: this
    type(bounds_type)         , intent(in) :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer           :: k,l,ii,jj 
    integer           :: begp,endp
    integer           :: begc,endc
    character(24)     :: fieldname
    character(100)    :: longname
    real(r8), pointer :: data1dptr(:)   ! temp. pointer for slicing larger arrays
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc

    !-------------------------------
    ! patch state variables 
    !-------------------------------
    
    if (crop_prog) then
       this%grainn_patch(begp:endp) = spval
       call hist_addfld1d (fname='GRAINN', units='gN/m^2', &
            avgflag='A', long_name='grain N', &
            ptr_patch=this%grainn_patch)
    end if

    this%leafn_patch(begp:endp) = spval
    call hist_addfld1d (fname='LEAFN', units='gN/m^2', &
         avgflag='A', long_name='leaf N', &
         ptr_patch=this%leafn_patch)

    this%leafn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='LEAFN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='leaf N storage', &
         ptr_patch=this%leafn_storage_patch, default='inactive')

    this%leafn_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='LEAFN_XFER', units='gN/m^2', &
         avgflag='A', long_name='leaf N transfer', &
         ptr_patch=this%leafn_xfer_patch, default='inactive')

    this%frootn_patch(begp:endp) = spval
    call hist_addfld1d (fname='FROOTN', units='gN/m^2', &
         avgflag='A', long_name='fine root N', &
         ptr_patch=this%frootn_patch)

    this%frootn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='FROOTN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='fine root N storage', &
         ptr_patch=this%frootn_storage_patch, default='inactive')

    this%frootn_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='FROOTN_XFER', units='gN/m^2', &
         avgflag='A', long_name='fine root N transfer', &
         ptr_patch=this%frootn_xfer_patch, default='inactive')

    this%livestemn_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMN', units='gN/m^2', &
         avgflag='A', long_name='live stem N', &
         ptr_patch=this%livestemn_patch)

    this%livestemn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='live stem N storage', &
         ptr_patch=this%livestemn_storage_patch, default='inactive')

    this%livestemn_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMN_XFER', units='gN/m^2', &
         avgflag='A', long_name='live stem N transfer', &
         ptr_patch=this%livestemn_xfer_patch, default='inactive')

    this%deadstemn_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMN', units='gN/m^2', &
         avgflag='A', long_name='dead stem N', &
         ptr_patch=this%deadstemn_patch)

    this%deadstemn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='dead stem N storage', &
         ptr_patch=this%deadstemn_storage_patch, default='inactive')

    this%deadstemn_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMN_XFER', units='gN/m^2', &
         avgflag='A', long_name='dead stem N transfer', &
         ptr_patch=this%deadstemn_xfer_patch, default='inactive')

    this%livecrootn_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTN', units='gN/m^2', &
         avgflag='A', long_name='live coarse root N', &
         ptr_patch=this%livecrootn_patch)

    this%livecrootn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='live coarse root N storage', &
         ptr_patch=this%livecrootn_storage_patch, default='inactive')

    this%livecrootn_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTN_XFER', units='gN/m^2', &
         avgflag='A', long_name='live coarse root N transfer', &
         ptr_patch=this%livecrootn_xfer_patch, default='inactive')

    this%deadcrootn_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTN', units='gN/m^2', &
         avgflag='A', long_name='dead coarse root N', &
         ptr_patch=this%deadcrootn_patch)

    this%deadcrootn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='dead coarse root N storage', &
         ptr_patch=this%deadcrootn_storage_patch, default='inactive')

    this%deadcrootn_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTN_XFER', units='gN/m^2', &
         avgflag='A', long_name='dead coarse root N transfer', &
         ptr_patch=this%deadcrootn_xfer_patch, default='inactive')

    this%retransn_patch(begp:endp) = spval
    call hist_addfld1d (fname='RETRANSN', units='gN/m^2', &
         avgflag='A', long_name='plant pool of retranslocated N', &
         ptr_patch=this%retransn_patch)

    this%npool_patch(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL', units='gN/m^2', &
         avgflag='A', long_name='temporary plant N pool', &
         ptr_patch=this%npool_patch, default='inactive')

    this%ntrunc_patch(begp:endp) = spval
    call hist_addfld1d (fname='PFT_NTRUNC', units='gN/m^2', &
         avgflag='A', long_name='patch-level sink for N truncation', &
         ptr_patch=this%ntrunc_patch)

    this%dispvegn_patch(begp:endp) = spval
    call hist_addfld1d (fname='DISPVEGN', units='gN/m^2', &
         avgflag='A', long_name='displayed vegetation nitrogen', &
         ptr_patch=this%dispvegn_patch)

    this%storvegn_patch(begp:endp) = spval
    call hist_addfld1d (fname='STORVEGN', units='gN/m^2', &
         avgflag='A', long_name='stored vegetation nitrogen', &
         ptr_patch=this%storvegn_patch)

    this%totvegn_patch(begp:endp) = spval
    call hist_addfld1d (fname='TOTVEGN', units='gN/m^2', &
         avgflag='A', long_name='total vegetation nitrogen', &
         ptr_patch=this%totvegn_patch)

    this%totn_patch(begp:endp) = spval
    call hist_addfld1d (fname='TOTPFTN', units='gN/m^2', &
         avgflag='A', long_name='total patch-level nitrogen', &
         ptr_patch=this%totn_patch)

    !-------------------------------
    ! column state variables 
    !-------------------------------

    this%seedn_col(begc:endc) = spval
    call hist_addfld1d (fname='SEEDN', units='gN/m^2', &
         avgflag='A', long_name='pool for seeding new patches', &
         ptr_col=this%seedn_col)

    this%prod10n_col(begc:endc) = spval
    call hist_addfld1d (fname='PROD10N', units='gN/m^2', &
         avgflag='A', long_name='10-yr wood product N', &
         ptr_col=this%prod10n_col)

    this%prod100n_col(begc:endc) = spval
    call hist_addfld1d (fname='PROD100N', units='gN/m^2', &
         avgflag='A', long_name='100-yr wood product N', &
         ptr_col=this%prod100n_col)

    this%totprodn_col(begc:endc) = spval
    call hist_addfld1d (fname='TOTPRODN', units='gN/m^2', &
         avgflag='A', long_name='total wood product N', &
         ptr_col=this%totprodn_col)

    this%totecosysn_col(begc:endc) = spval
    call hist_addfld1d (fname='TOTECOSYSN', units='gN/m^2', &
         avgflag='A', long_name='total ecosystem N', &
         ptr_col=this%totecosysn_col)

    this%totn_col(begc:endc) = spval
    call hist_addfld1d (fname='TOTCOLN', units='gN/m^2', &
         avgflag='A', long_name='total column-level N', &
         ptr_col=this%totn_col)

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, &
       leafc_patch, leafc_storage_patch, deadstemc_patch)
    !
    ! !DESCRIPTION:
    ! Initializes time varying variables used only in coupled carbon-nitrogen mode (CN):
    !
    ! !ARGUMENTS:
    class(cnveg_nitrogenstate_type) :: this
    type(bounds_type) , intent(in) :: bounds  
    real(r8)          , intent(in) :: leafc_patch(bounds%begp:)
    real(r8)          , intent(in) :: leafc_storage_patch(bounds%begp:)
    real(r8)          , intent(in) :: deadstemc_patch(bounds%begp:)
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
    SHR_ASSERT_ALL((ubound(deadstemc_patch)      == (/bounds%endp/)),                               errMsg(__FILE__, __LINE__))

    ! Set column filters

    num_special_patch = 0
    do p = bounds%begp,bounds%endp
       l = patch%landunit(p)
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
    ! initialize patch-level variables
    !-------------------------------------------

    do p = bounds%begp,bounds%endp

       l = patch%landunit(p)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then

          if (patch%itype(p) == noveg) then
             this%leafn_patch(p) = 0._r8
             this%leafn_storage_patch(p) = 0._r8
          else
             this%leafn_patch(p)         = leafc_patch(p)         / pftcon%leafcn(patch%itype(p))
             this%leafn_storage_patch(p) = leafc_storage_patch(p) / pftcon%leafcn(patch%itype(p))
          end if

          this%leafn_xfer_patch(p)        = 0._r8
          if ( crop_prog )then
             this%grainn_patch(p)         = 0._r8
             this%grainn_storage_patch(p) = 0._r8
             this%grainn_xfer_patch(p)    = 0._r8
          end if
          this%frootn_patch(p)            = 0._r8
          this%frootn_storage_patch(p)    = 0._r8
          this%frootn_xfer_patch(p)       = 0._r8
          this%livestemn_patch(p)         = 0._r8
          this%livestemn_storage_patch(p) = 0._r8
          this%livestemn_xfer_patch(p)    = 0._r8

          ! tree types need to be initialized with some stem mass so that
          ! roughness length is not zero in canopy flux calculation

          if (pftcon%woody(patch%itype(p)) == 1._r8) then
             this%deadstemn_patch(p) = deadstemc_patch(p) / pftcon%deadwdcn(patch%itype(p))
          else
             this%deadstemn_patch(p) = 0._r8
          end if

          this%deadstemn_storage_patch(p)  = 0._r8
          this%deadstemn_xfer_patch(p)     = 0._r8
          this%livecrootn_patch(p)         = 0._r8
          this%livecrootn_storage_patch(p) = 0._r8
          this%livecrootn_xfer_patch(p)    = 0._r8
          this%deadcrootn_patch(p)         = 0._r8
          this%deadcrootn_storage_patch(p) = 0._r8
          this%deadcrootn_xfer_patch(p)    = 0._r8
          this%retransn_patch(p)           = 0._r8
          this%npool_patch(p)              = 0._r8
          this%ntrunc_patch(p)             = 0._r8
          this%dispvegn_patch(p)           = 0._r8
          this%storvegn_patch(p)           = 0._r8
          this%totvegn_patch(p)            = 0._r8
          this%totn_patch(p)               = 0._r8
       end if
    end do

    !-------------------------------------------
    ! initialize column-level variables
    !-------------------------------------------

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then

          ! dynamic landcover state variables
          this%seedn_col(c)      = 0._r8
          this%prod10n_col(c)    = 0._r8
          this%prod100n_col(c)   = 0._r8
          this%totprodn_col(c)   = 0._r8

          ! total nitrogen pools
          this%totecosysn_col(c) = 0._r8
          this%totn_col(c)       = 0._r8
       end if
    end do

    ! now loop through special filters and explicitly set the variables that
    ! have to be in place for biogeophysics

    do fc = 1,num_special_col
       c = special_col(fc)
       this%seedn_col(c)    = 0._r8
       this%prod10n_col(c)  = 0._r8	  
       this%prod100n_col(c) = 0._r8	  
       this%totprodn_col(c) = 0._r8	  
    end do

    ! initialize fields for special filters

    call this%SetValues (&
         num_patch=num_special_patch, filter_patch=special_patch, value_patch=0._r8, &
         num_column=num_special_col, filter_column=special_col, value_column=0._r8)

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine Restart ( this,  bounds, ncid, flag )
    !
    ! !DESCRIPTION: 
    ! Read/write restart data 
    !
    ! !USES:
    use restUtilMod
    use ncdio_pio
    !
    ! !ARGUMENTS:
    class (cnveg_nitrogenstate_type) :: this
    type(bounds_type)          , intent(in)    :: bounds 
    type(file_desc_t)          , intent(inout) :: ncid   
    character(len=*)           , intent(in)    :: flag   !'read' or 'write' or 'define'
    !
    ! !LOCAL VARIABLES:
    integer            :: i,j,k,l,c
    logical            :: readvar
    real(r8), pointer  :: ptr1d(:)   ! temp. pointers for slicing larger arrays
    character(len=128) :: varname    ! temporary
    !------------------------------------------------------------------------

    !--------------------------------
    ! patch nitrogen state variables
    !--------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='leafn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%leafn_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='leafn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%leafn_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='leafn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%leafn_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='frootn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%frootn_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='frootn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%frootn_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='frootn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%frootn_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livestemn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livestemn_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livestemn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livestemn_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livestemn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livestemn_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadstemn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadstemn_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadstemn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadstemn_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadstemn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadstemn_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livecrootn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livecrootn_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livecrootn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livecrootn_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livecrootn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livecrootn_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadcrootn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadcrootn_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadcrootn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadcrootn_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadcrootn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadcrootn_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='retransn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%retransn_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='npool', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%npool_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='pft_ntrunc', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%ntrunc_patch) 

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag,  varname='grainn', xtype=ncd_double,  &
            dim1name='pft',    long_name='grain N', units='gN/m2', &
            interpinic_flag='interp', readvar=readvar, data=this%grainn_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='grainn_storage', xtype=ncd_double,  &
            dim1name='pft',    long_name='grain N storage', units='gN/m2', &
            interpinic_flag='interp', readvar=readvar, data=this%grainn_storage_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='grainn_xfer', xtype=ncd_double,  &
            dim1name='pft',    long_name='grain N transfer', units='gN/m2', &
            interpinic_flag='interp', readvar=readvar, data=this%grainn_xfer_patch)
    end if

    !--------------------------------
    ! column nitrogen state variables
    !--------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='seedn', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%seedn_col) 

    call restartvar(ncid=ncid, flag=flag, varname='prod10n', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%prod10n_col) 

    call restartvar(ncid=ncid, flag=flag, varname='prod100n', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%prod100n_col) 

    call restartvar(ncid=ncid, flag=flag, varname='totcoln', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%totn_col) 

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
    class (cnveg_nitrogenstate_type) :: this
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

       this%leafn_patch(i)              = value_patch
       this%leafn_storage_patch(i)      = value_patch
       this%leafn_xfer_patch(i)         = value_patch
       this%frootn_patch(i)             = value_patch
       this%frootn_storage_patch(i)     = value_patch
       this%frootn_xfer_patch(i)        = value_patch
       this%livestemn_patch(i)          = value_patch
       this%livestemn_storage_patch(i)  = value_patch
       this%livestemn_xfer_patch(i)     = value_patch
       this%deadstemn_patch(i)          = value_patch
       this%deadstemn_storage_patch(i)  = value_patch
       this%deadstemn_xfer_patch(i)     = value_patch
       this%livecrootn_patch(i)         = value_patch
       this%livecrootn_storage_patch(i) = value_patch
       this%livecrootn_xfer_patch(i)    = value_patch
       this%deadcrootn_patch(i)         = value_patch
       this%deadcrootn_storage_patch(i) = value_patch
       this%deadcrootn_xfer_patch(i)    = value_patch
       this%retransn_patch(i)           = value_patch
       this%npool_patch(i)              = value_patch
       this%ntrunc_patch(i)             = value_patch
       this%dispvegn_patch(i)           = value_patch
       this%storvegn_patch(i)           = value_patch
       this%totvegn_patch(i)            = value_patch
       this%totn_patch(i)               = value_patch
    end do

    if ( crop_prog )then
       do fi = 1,num_patch
          i = filter_patch(fi)
          this%grainn_patch(i)          = value_patch
          this%grainn_storage_patch(i)  = value_patch
          this%grainn_xfer_patch(i)     = value_patch   
       end do
    end if

    do fi = 1,num_column
       i = filter_column(fi)

       this%totecosysn_col(i) = value_column
       this%totn_col(i)       = value_column
    end do

  end subroutine SetValues

  !-----------------------------------------------------------------------
  subroutine ZeroDwt( this, bounds )
    !
    ! !DESCRIPTION
    ! Initialize variables needed for dynamic land use.
    !
    ! !ARGUMENTS:
    class(cnveg_nitrogenstate_type) :: this
    type(bounds_type), intent(in)  :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer  :: p          ! indices
    !-----------------------------------------------------------------------

    do p = bounds%begp,bounds%endp
       this%dispvegn_patch(p) = 0._r8
       this%storvegn_patch(p) = 0._r8
       this%totvegn_patch(p)  = 0._r8
       this%totn_patch(p)     = 0._r8
    end do

  end subroutine ZeroDwt

  !-----------------------------------------------------------------------
  subroutine Summary_nitrogenstate(this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp,&
       soilbiogeochem_nitrogenstate_inst)
    !
    ! !USES:
    use subgridAveMod, only : p2c
    use SoilBiogeochemNitrogenStateType, only : soilbiogeochem_nitrogenstate_type
    !
    ! !ARGUMENTS:
    class(cnveg_nitrogenstate_type)                      :: this
    type(bounds_type)                       , intent(in) :: bounds  
    integer                                 , intent(in) :: num_soilc       ! number of soil columns in filter
    integer                                 , intent(in) :: filter_soilc(:) ! filter for soil columns
    integer                                 , intent(in) :: num_soilp       ! number of soil patches in filter
    integer                                 , intent(in) :: filter_soilp(:) ! filter for soil patches
    type(soilbiogeochem_nitrogenstate_type) , intent(in) :: soilbiogeochem_nitrogenstate_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: c,p,j,k,l   ! indices
    integer  :: fp,fc       ! lake filter indices
    real(r8) :: maxdepth    ! depth to integrate soil variables
    !-----------------------------------------------------------------------

    ! --------------------------------------------
    ! patch level summary
    ! --------------------------------------------
    
    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! displayed vegetation nitrogen, excluding storage (DISPVEGN)
       this%dispvegn_patch(p) = &
            this%leafn_patch(p)      + &
            this%frootn_patch(p)     + &
            this%livestemn_patch(p)  + &
            this%deadstemn_patch(p)  + &
            this%livecrootn_patch(p) + &
            this%deadcrootn_patch(p)

       ! stored vegetation nitrogen, including retranslocated N pool (STORVEGN)
       this%storvegn_patch(p) = &
            this%leafn_storage_patch(p)      + &
            this%frootn_storage_patch(p)     + &
            this%livestemn_storage_patch(p)  + &
            this%deadstemn_storage_patch(p)  + &
            this%livecrootn_storage_patch(p) + &
            this%deadcrootn_storage_patch(p) + &
            this%leafn_xfer_patch(p)         + &
            this%frootn_xfer_patch(p)        + &
            this%livestemn_xfer_patch(p)     + &
            this%deadstemn_xfer_patch(p)     + &
            this%livecrootn_xfer_patch(p)    + &
            this%deadcrootn_xfer_patch(p)    + &
            this%npool_patch(p)              + &
            this%retransn_patch(p)

       if ( crop_prog .and. patch%itype(p) >= npcropmin )then
          this%dispvegn_patch(p) = &
               this%dispvegn_patch(p) + &
               this%grainn_patch(p)

          this%storvegn_patch(p) = &
               this%storvegn_patch(p) + &
               this%grainn_storage_patch(p)     + &
               this%grainn_xfer_patch(p)
       end if

       ! total vegetation nitrogen (TOTVEGN)
       this%totvegn_patch(p) = &
            this%dispvegn_patch(p) + &
            this%storvegn_patch(p)

       ! total patch-level carbon (add ntrunc)
       this%totn_patch(p) = &
            this%totvegn_patch(p) + &
            this%ntrunc_patch(p)

    end do

    ! --------------------------------------------
    ! column level summary
    ! --------------------------------------------

    call p2c(bounds, num_soilc, filter_soilc, &
         this%totvegn_patch(bounds%begp:bounds%endp), &
         this%totvegn_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%totn_patch(bounds%begp:bounds%endp), &
         this%totn_col(bounds%begc:bounds%endc))

    ! total wood product nitrogen
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%totprodn_col(c) = &
            this%prod10n_col(c) + &
            this%prod100n_col(c)	 

       ! total ecosystem nitrogen, including veg (TOTECOSYSN)
       this%totecosysn_col(c) =    &
            soilbiogeochem_nitrogenstate_inst%cwdn_col(c)    + &
            soilbiogeochem_nitrogenstate_inst%totlitn_col(c) + &
            soilbiogeochem_nitrogenstate_inst%totsomn_col(c) + &
            soilbiogeochem_nitrogenstate_inst%sminn_col(c)   + &
            this%totprodn_col(c)                             + &
            this%totvegn_col(c)                              

       ! total column nitrogen, including patch (TOTCOLN)

       this%totn_col(c) = this%totn_col(c)                   + &
            soilbiogeochem_nitrogenstate_inst%cwdn_col(c)    + &
            soilbiogeochem_nitrogenstate_inst%totlitn_col(c) + &
            soilbiogeochem_nitrogenstate_inst%totsomn_col(c) + &
            soilbiogeochem_nitrogenstate_inst%sminn_col(c)   + &
            this%totprodn_col(c)                             + &
            this%seedn_col(c)                                + &
            soilbiogeochem_nitrogenstate_inst%ntrunc_col(c)
    end do

  end subroutine Summary_nitrogenstate

end module CNVegNitrogenStateType
