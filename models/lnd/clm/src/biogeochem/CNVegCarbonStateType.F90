module CNVegCarbonStateType

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use shr_const_mod  , only : SHR_CONST_PDB
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use pftconMod	     , only : noveg, npcropmin, pftcon
  use clm_varpar     , only : crop_prog
  use clm_varcon     , only : spval
  use clm_varctl     , only : iulog, use_cndv
  use decompMod      , only : bounds_type
  use abortutils     , only : endrun
  use spmdMod        , only : masterproc 
  use LandunitType   , only : lun                
  use ColumnType     , only : col                
  use PatchType      , only : patch
  ! 
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  type, public :: cnveg_carbonstate_type
     
     real(r8), pointer :: grainc_patch             (:) ! (gC/m2) grain C (crop model)
     real(r8), pointer :: grainc_storage_patch     (:) ! (gC/m2) grain C storage (crop model)
     real(r8), pointer :: grainc_xfer_patch        (:) ! (gC/m2) grain C transfer (crop model)
     real(r8), pointer :: leafc_patch              (:) ! (gC/m2) leaf C
     real(r8), pointer :: leafc_storage_patch      (:) ! (gC/m2) leaf C storage
     real(r8), pointer :: leafc_xfer_patch         (:) ! (gC/m2) leaf C transfer
     real(r8), pointer :: frootc_patch             (:) ! (gC/m2) fine root C
     real(r8), pointer :: frootc_storage_patch     (:) ! (gC/m2) fine root C storage
     real(r8), pointer :: frootc_xfer_patch        (:) ! (gC/m2) fine root C transfer
     real(r8), pointer :: livestemc_patch          (:) ! (gC/m2) live stem C
     real(r8), pointer :: livestemc_storage_patch  (:) ! (gC/m2) live stem C storage
     real(r8), pointer :: livestemc_xfer_patch     (:) ! (gC/m2) live stem C transfer
     real(r8), pointer :: deadstemc_patch          (:) ! (gC/m2) dead stem C
     real(r8), pointer :: deadstemc_storage_patch  (:) ! (gC/m2) dead stem C storage
     real(r8), pointer :: deadstemc_xfer_patch     (:) ! (gC/m2) dead stem C transfer
     real(r8), pointer :: livecrootc_patch         (:) ! (gC/m2) live coarse root C
     real(r8), pointer :: livecrootc_storage_patch (:) ! (gC/m2) live coarse root C storage
     real(r8), pointer :: livecrootc_xfer_patch    (:) ! (gC/m2) live coarse root C transfer
     real(r8), pointer :: deadcrootc_patch         (:) ! (gC/m2) dead coarse root C
     real(r8), pointer :: deadcrootc_storage_patch (:) ! (gC/m2) dead coarse root C storage
     real(r8), pointer :: deadcrootc_xfer_patch    (:) ! (gC/m2) dead coarse root C transfer
     real(r8), pointer :: gresp_storage_patch      (:) ! (gC/m2) growth respiration storage
     real(r8), pointer :: gresp_xfer_patch         (:) ! (gC/m2) growth respiration transfer
     real(r8), pointer :: cpool_patch              (:) ! (gC/m2) temporary photosynthate C pool
     real(r8), pointer :: xsmrpool_patch           (:) ! (gC/m2) abstract C pool to meet excess MR demand
     real(r8), pointer :: ctrunc_patch             (:) ! (gC/m2) patch-level sink for C truncation
     real(r8), pointer :: woodc_patch              (:) ! (gC/m2) wood C
     real(r8), pointer :: leafcmax_patch           (:) ! (gC/m2) ann max leaf C
     real(r8), pointer :: totc_patch               (:) ! (gC/m2) total patch-level carbon, including cpool
     real(r8), pointer :: rootc_col                (:) ! (gC/m2) root carbon at column level (fire)
     real(r8), pointer :: leafc_col                (:) ! (gC/m2) column-level leafc (fire)
     real(r8), pointer :: fuelc_col                (:) ! (0-1) fuel avalability factor for Reg.C 
     real(r8), pointer :: fuelc_crop_col           (:) ! (0-1) fuel avalability factor for Reg.A 

     ! pools for dynamic landcover
     real(r8), pointer :: seedc_col                (:) ! (gC/m2) column-level pool for seeding new Patches
     real(r8), pointer :: prod10c_col              (:) ! (gC/m2) wood product C pool, 10-year lifespan
     real(r8), pointer :: prod100c_col             (:) ! (gC/m2) wood product C pool, 100-year lifespan
     real(r8), pointer :: totprodc_col             (:) ! (gC/m2) total wood product C

     ! summary (diagnostic) state variables, not involved in mass balance
     real(r8), pointer :: dispvegc_patch           (:) ! (gC/m2) displayed veg carbon, excluding storage and cpool
     real(r8), pointer :: storvegc_patch           (:) ! (gC/m2) stored vegetation carbon, excluding cpool
     real(r8), pointer :: totvegc_patch            (:) ! (gC/m2) total vegetation carbon, excluding cpool
     real(r8), pointer :: totvegc_col              (:) ! (gC/m2) total vegetation carbon, excluding cpool averaged to column (p2c)

     ! Total C pools       
     real(r8), pointer :: totc_col                 (:) ! (gC/m2) total column carbon, incl veg and cpool 
     real(r8), pointer :: totecosysc_col           (:) ! (gC/m2) total ecosystem carbon, incl veg but excl cpool 
     real(r8), pointer :: begcb_col                (:) ! (gC/m2) carbon mass, beginning of time step 
     real(r8), pointer :: endcb_col                (:) ! (gC/m2) carbon mass, end of time step

   contains

     procedure , public  :: Init   
     procedure , public  :: SetValues 
     procedure , public  :: ZeroDWT
     procedure , public  :: Restart
     procedure , public  :: Summary => Summary_carbonstate
     procedure , private :: InitAllocate 
     procedure , private :: InitHistory  
     procedure , private :: InitCold     

  end type cnveg_carbonstate_type

  real(r8) :: c3_r2 ! isotope ratio (13c/[12c+13c]) for C3 photosynthesis
  real(r8) :: c4_r2 ! isotope ratio (13c/[12c+13c]) for C4 photosynthesis
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, carbon_type, ratio, c12_cnveg_carbonstate_inst)

    class(cnveg_carbonstate_type)                       :: this
    type(bounds_type)            , intent(in)           :: bounds  
    real(r8)                     , intent(in)           :: ratio
    character(len=3)             , intent(in)           :: carbon_type
    type(cnveg_carbonstate_type) , intent(in), optional :: c12_cnveg_carbonstate_inst
    !-----------------------------------------------------------------------

    call this%InitAllocate ( bounds)
    call this%InitHistory ( bounds, carbon_type)
    if (present(c12_cnveg_carbonstate_inst)) then
       call this%InitCold  ( bounds, ratio, carbon_type, c12_cnveg_carbonstate_inst)
    else
       call this%InitCold  ( bounds, ratio, carbon_type)
    end if

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !ARGUMENTS:
    class (cnveg_carbonstate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begp,endp
    integer           :: begc,endc
    !------------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc

    allocate(this%leafc_patch              (begp:endp)) ; this%leafc_patch              (:) = nan
    allocate(this%leafc_storage_patch      (begp:endp)) ; this%leafc_storage_patch      (:) = nan
    allocate(this%leafc_xfer_patch         (begp:endp)) ; this%leafc_xfer_patch         (:) = nan
    allocate(this%frootc_patch             (begp:endp)) ; this%frootc_patch             (:) = nan
    allocate(this%frootc_storage_patch     (begp:endp)) ; this%frootc_storage_patch     (:) = nan
    allocate(this%frootc_xfer_patch        (begp:endp)) ; this%frootc_xfer_patch        (:) = nan
    allocate(this%livestemc_patch          (begp:endp)) ; this%livestemc_patch          (:) = nan
    allocate(this%livestemc_storage_patch  (begp:endp)) ; this%livestemc_storage_patch  (:) = nan
    allocate(this%livestemc_xfer_patch     (begp:endp)) ; this%livestemc_xfer_patch     (:) = nan
    allocate(this%deadstemc_patch          (begp:endp)) ; this%deadstemc_patch          (:) = nan
    allocate(this%deadstemc_storage_patch  (begp:endp)) ; this%deadstemc_storage_patch  (:) = nan
    allocate(this%deadstemc_xfer_patch     (begp:endp)) ; this%deadstemc_xfer_patch     (:) = nan
    allocate(this%livecrootc_patch         (begp:endp)) ; this%livecrootc_patch         (:) = nan
    allocate(this%livecrootc_storage_patch (begp:endp)) ; this%livecrootc_storage_patch (:) = nan
    allocate(this%livecrootc_xfer_patch    (begp:endp)) ; this%livecrootc_xfer_patch    (:) = nan
    allocate(this%deadcrootc_patch         (begp:endp)) ; this%deadcrootc_patch         (:) = nan
    allocate(this%deadcrootc_storage_patch (begp:endp)) ; this%deadcrootc_storage_patch (:) = nan
    allocate(this%deadcrootc_xfer_patch    (begp:endp)) ; this%deadcrootc_xfer_patch    (:) = nan
    allocate(this%gresp_storage_patch      (begp:endp)) ; this%gresp_storage_patch      (:) = nan
    allocate(this%gresp_xfer_patch         (begp:endp)) ; this%gresp_xfer_patch         (:) = nan
    allocate(this%cpool_patch              (begp:endp)) ; this%cpool_patch              (:) = nan
    allocate(this%xsmrpool_patch           (begp:endp)) ; this%xsmrpool_patch           (:) = nan
    allocate(this%ctrunc_patch             (begp:endp)) ; this%ctrunc_patch             (:) = nan
    allocate(this%dispvegc_patch           (begp:endp)) ; this%dispvegc_patch           (:) = nan
    allocate(this%storvegc_patch           (begp:endp)) ; this%storvegc_patch           (:) = nan
    allocate(this%leafcmax_patch           (begp:endp)) ; this%leafcmax_patch           (:) = nan
    allocate(this%totc_patch               (begp:endp))  ; this%totc_patch               (:) = nan
    allocate(this%grainc_patch             (begp:endp)) ; this%grainc_patch             (:) = nan
    allocate(this%grainc_storage_patch     (begp:endp)) ; this%grainc_storage_patch     (:) = nan
    allocate(this%grainc_xfer_patch        (begp:endp)) ; this%grainc_xfer_patch        (:) = nan
    allocate(this%woodc_patch              (begp:endp)) ; this%woodc_patch              (:) = nan     

    allocate(this%seedc_col                (begc:endc)) ; this%seedc_col                (:) = nan
    allocate(this%prod10c_col              (begc:endc)) ; this%prod10c_col              (:) = nan
    allocate(this%prod100c_col             (begc:endc)) ; this%prod100c_col             (:) = nan
    allocate(this%totprodc_col             (begc:endc)) ; this%totprodc_col             (:) = nan
    allocate(this%rootc_col                (begc:endc)) ; this%rootc_col                (:) = nan
    allocate(this%leafc_col                (begc:endc)) ; this%leafc_col                (:) = nan
    allocate(this%fuelc_col                (begc:endc)) ; this%fuelc_col                (:) = nan
    allocate(this%fuelc_crop_col           (begc:endc)) ; this%fuelc_crop_col           (:) = nan

    allocate(this%totvegc_patch            (begp:endp)) ; this%totvegc_patch            (:) = nan
    allocate(this%totvegc_col              (begc:endc)) ; this%totvegc_col              (:) = nan

    allocate(this%totc_col                 (begc:endc)) ; this%totc_col                 (:) = nan
    allocate(this%totecosysc_col           (begc:endc)) ; this%totecosysc_col           (:) = nan
    allocate(this%begcb_col                (begc:endc)) ; this%begcb_col                (:) = nan
    allocate(this%endcb_col                (begc:endc)) ; this%endcb_col                (:) = nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds, carbon_type)
    !
    ! !DESCRIPTION:
    ! add history fields for all CN variables, always set as default='inactive'
    !
    ! !USES:
    use clm_varctl , only : use_c13, use_c14
    use histFileMod, only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp 
    !
    ! !ARGUMENTS:
    class (cnveg_carbonstate_type) :: this
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

    !-------------------------------
    ! C12 state variables
    !-------------------------------

    if (carbon_type == 'c12') then

       if (crop_prog) then
          this%grainc_patch(begp:endp) = spval
          call hist_addfld1d (fname='GRAINC', units='gC/m^2', &
               avgflag='A', long_name='grain C', &
               ptr_patch=this%grainc_patch)
       end if
       
       this%woodc_patch(begp:endp) = spval
       call hist_addfld1d (fname='WOODC', units='gC/m^2', &
            avgflag='A', long_name='wood C', &
            ptr_patch=this%woodc_patch)

       this%leafc_patch(begp:endp) = spval
       call hist_addfld1d (fname='LEAFC', units='gC/m^2', &
            avgflag='A', long_name='leaf C', &
            ptr_patch=this%leafc_patch)

       this%leafc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='LEAFC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='leaf C storage', &
            ptr_patch=this%leafc_storage_patch, default='inactive')

       this%leafc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='LEAFC_XFER', units='gC/m^2', &
            avgflag='A', long_name='leaf C transfer', &
            ptr_patch=this%leafc_xfer_patch, default='inactive')

       this%frootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='FROOTC', units='gC/m^2', &
            avgflag='A', long_name='fine root C', &
            ptr_patch=this%frootc_patch)

       this%frootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='FROOTC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='fine root C storage', &
            ptr_patch=this%frootc_storage_patch, default='inactive')

       this%frootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='FROOTC_XFER', units='gC/m^2', &
            avgflag='A', long_name='fine root C transfer', &
            ptr_patch=this%frootc_xfer_patch, default='inactive')

       this%livestemc_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVESTEMC', units='gC/m^2', &
            avgflag='A', long_name='live stem C', &
            ptr_patch=this%livestemc_patch)

       this%livestemc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVESTEMC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='live stem C storage', &
            ptr_patch=this%livestemc_storage_patch, default='inactive')

       this%livestemc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVESTEMC_XFER', units='gC/m^2', &
            avgflag='A', long_name='live stem C transfer', &
            ptr_patch=this%livestemc_xfer_patch, default='inactive')

       this%deadstemc_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADSTEMC', units='gC/m^2', &
            avgflag='A', long_name='dead stem C', &
            ptr_patch=this%deadstemc_patch)

       this%deadstemc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADSTEMC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='dead stem C storage', &
            ptr_patch=this%deadstemc_storage_patch, default='inactive')

       this%deadstemc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADSTEMC_XFER', units='gC/m^2', &
            avgflag='A', long_name='dead stem C transfer', &
            ptr_patch=this%deadstemc_xfer_patch, default='inactive')

       this%livecrootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVECROOTC', units='gC/m^2', &
            avgflag='A', long_name='live coarse root C', &
            ptr_patch=this%livecrootc_patch)

       this%livecrootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVECROOTC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='live coarse root C storage', &
            ptr_patch=this%livecrootc_storage_patch, default='inactive')

       this%livecrootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='LIVECROOTC_XFER', units='gC/m^2', &
            avgflag='A', long_name='live coarse root C transfer', &
            ptr_patch=this%livecrootc_xfer_patch, default='inactive')

       this%deadcrootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADCROOTC', units='gC/m^2', &
            avgflag='A', long_name='dead coarse root C', &
            ptr_patch=this%deadcrootc_patch)

       this%deadcrootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADCROOTC_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='dead coarse root C storage', &
            ptr_patch=this%deadcrootc_storage_patch,  default='inactive')

       this%deadcrootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='DEADCROOTC_XFER', units='gC/m^2', &
            avgflag='A', long_name='dead coarse root C transfer', &
            ptr_patch=this%deadcrootc_xfer_patch, default='inactive')

       this%gresp_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='GRESP_STORAGE', units='gC/m^2', &
            avgflag='A', long_name='growth respiration storage', &
            ptr_patch=this%gresp_storage_patch, default='inactive')

       this%gresp_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='GRESP_XFER', units='gC/m^2', &
            avgflag='A', long_name='growth respiration transfer', &
            ptr_patch=this%gresp_xfer_patch, default='inactive')

       this%cpool_patch(begp:endp) = spval
       call hist_addfld1d (fname='CPOOL', units='gC/m^2', &
            avgflag='A', long_name='temporary photosynthate C pool', &
            ptr_patch=this%cpool_patch)

       this%xsmrpool_patch(begp:endp) = spval
       call hist_addfld1d (fname='XSMRPOOL', units='gC/m^2', &
            avgflag='A', long_name='temporary photosynthate C pool', &
            ptr_patch=this%xsmrpool_patch)

       this%ctrunc_patch(begp:endp) = spval
       call hist_addfld1d (fname='PFT_CTRUNC', units='gC/m^2', &
            avgflag='A', long_name='patch-level sink for C truncation', &
            ptr_patch=this%ctrunc_patch)

       this%dispvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='DISPVEGC', units='gC/m^2', &
            avgflag='A', long_name='displayed veg carbon, excluding storage and cpool', &
            ptr_patch=this%dispvegc_patch)

       this%storvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='STORVEGC', units='gC/m^2', &
            avgflag='A', long_name='stored vegetation carbon, excluding cpool', &
            ptr_patch=this%storvegc_patch)

       this%totvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='TOTVEGC', units='gC/m^2', &
            avgflag='A', long_name='total vegetation carbon, excluding cpool', &
            ptr_patch=this%totvegc_patch)

       this%totc_patch(begp:endp) = spval
       call hist_addfld1d (fname='TOTPFTC', units='gC/m^2', &
            avgflag='A', long_name='total patch-level carbon, including cpool', &
            ptr_patch=this%totc_patch)

       this%seedc_col(begc:endc) = spval
       call hist_addfld1d (fname='SEEDC', units='gC/m^2', &
            avgflag='A', long_name='pool for seeding new Patches', &
            ptr_col=this%seedc_col)

       this%prod10c_col(begc:endc) = spval
       call hist_addfld1d (fname='PROD10C', units='gC/m^2', &
            avgflag='A', long_name='10-yr wood product C', &
            ptr_col=this%prod10c_col)

       this%prod100c_col(begc:endc) = spval
       call hist_addfld1d (fname='PROD100C', units='gC/m^2', &
            avgflag='A', long_name='100-yr wood product C', &
            ptr_col=this%prod100c_col)

       this%totprodc_col(begc:endc) = spval
       call hist_addfld1d (fname='TOTPRODC', units='gC/m^2', &
            avgflag='A', long_name='total wood product C', &
            ptr_col=this%totprodc_col)

       this%fuelc_col(begc:endc) = spval
       call hist_addfld1d (fname='FUELC', units='gC/m^2', &
            avgflag='A', long_name='fuel load', &
            ptr_col=this%fuelc_col)

       this%totc_col(begc:endc) = spval
       call hist_addfld1d (fname='TOTCOLC', units='gC/m^2', &
            avgflag='A', long_name='total column carbon, incl veg and cpool', &
            ptr_col=this%totc_col)

       this%totecosysc_col(begc:endc) = spval
       call hist_addfld1d (fname='TOTECOSYSC', units='gC/m^2', &
            avgflag='A', long_name='total ecosystem carbon, incl veg but excl cpool', &
            ptr_col=this%totecosysc_col)

    end if

    !-------------------------------
    ! C13 state variables 
    !-------------------------------

    if ( carbon_type == 'c13' ) then

       this%leafc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LEAFC', units='gC13/m^2', &
            avgflag='A', long_name='C13 leaf C', &
            ptr_patch=this%leafc_patch)

       this%leafc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LEAFC_STORAGE', units='gC13/m^2', &
            avgflag='A', long_name='C13 leaf C storage', &
            ptr_patch=this%leafc_storage_patch, default='inactive')

       this%leafc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LEAFC_XFER', units='gC13/m^2', &
            avgflag='A', long_name='C13 leaf C transfer', &
            ptr_patch=this%leafc_xfer_patch, default='inactive')

       this%frootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_FROOTC', units='gC13/m^2', &
            avgflag='A', long_name='C13 fine root C', &
            ptr_patch=this%frootc_patch)

       this%frootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_FROOTC_STORAGE', units='gC13/m^2', &
            avgflag='A', long_name='C13 fine root C storage', &
            ptr_patch=this%frootc_storage_patch, default='inactive')

       this%frootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_FROOTC_XFER', units='gC13/m^2', &
            avgflag='A', long_name='C13 fine root C transfer', &
            ptr_patch=this%frootc_xfer_patch, default='inactive')

       this%livestemc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVESTEMC', units='gC13/m^2', &
            avgflag='A', long_name='C13 live stem C', &
            ptr_patch=this%livestemc_patch)

       this%livestemc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVESTEMC_STORAGE', units='gC13/m^2', &
            avgflag='A', long_name='C13 live stem C storage', &
            ptr_patch=this%livestemc_storage_patch, default='inactive')

       this%livestemc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVESTEMC_XFER', units='gC13/m^2', &
            avgflag='A', long_name='C13 live stem C transfer', &
            ptr_patch=this%livestemc_xfer_patch, default='inactive')

       this%deadstemc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADSTEMC', units='gC13/m^2', &
            avgflag='A', long_name='C13 dead stem C', &
            ptr_patch=this%deadstemc_patch)

       this%deadstemc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADSTEMC_STORAGE', units='gC13/m^2', &
            avgflag='A', long_name='C13 dead stem C storage', &
            ptr_patch=this%deadstemc_storage_patch, default='inactive')

       this%deadstemc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADSTEMC_XFER', units='gC13/m^2', &
            avgflag='A', long_name='C13 dead stem C transfer', &
            ptr_patch=this%deadstemc_xfer_patch, default='inactive')

       this%livecrootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVECROOTC', units='gC13/m^2', &
            avgflag='A', long_name='C13 live coarse root C', &
            ptr_patch=this%livecrootc_patch)

       this%livecrootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVECROOTC_STORAGE', units='gC13/m^2', &
            avgflag='A', long_name='C13 live coarse root C storage', &
            ptr_patch=this%livecrootc_storage_patch, default='inactive')

       this%livecrootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_LIVECROOTC_XFER', units='gC13/m^2', &
            avgflag='A', long_name='C13 live coarse root C transfer', &
            ptr_patch=this%livecrootc_xfer_patch, default='inactive')

       this%deadcrootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADCROOTC', units='gC13/m^2', &
            avgflag='A', long_name='C13 dead coarse root C', &
            ptr_patch=this%deadcrootc_patch)

       this%deadcrootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADCROOTC_STORAGE', units='gC13/m^2', &
            avgflag='A', long_name='C13 dead coarse root C storage', &
            ptr_patch=this%deadcrootc_storage_patch,  default='inactive')

       this%deadcrootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DEADCROOTC_XFER', units='gC13/m^2', &
            avgflag='A', long_name='C13 dead coarse root C transfer', &
            ptr_patch=this%deadcrootc_xfer_patch, default='inactive')

       this%gresp_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_GRESP_STORAGE', units='gC13/m^2', &
            avgflag='A', long_name='C13 growth respiration storage', &
            ptr_patch=this%gresp_storage_patch, default='inactive')

       this%gresp_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_GRESP_XFER', units='gC13/m^2', &
            avgflag='A', long_name='C13 growth respiration transfer', &
            ptr_patch=this%gresp_xfer_patch, default='inactive')

       this%cpool_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_CPOOL', units='gC13/m^2', &
            avgflag='A', long_name='C13 temporary photosynthate C pool', &
            ptr_patch=this%cpool_patch)

       this%xsmrpool_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_XSMRPOOL', units='gC13/m^2', &
            avgflag='A', long_name='C13 temporary photosynthate C pool', &
            ptr_patch=this%xsmrpool_patch)

       this%ctrunc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_PFT_CTRUNC', units='gC13/m^2', &
            avgflag='A', long_name='C13 patch-level sink for C truncation', &
            ptr_patch=this%ctrunc_patch)

       this%dispvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_DISPVEGC', units='gC13/m^2', &
            avgflag='A', long_name='C13 displayed veg carbon, excluding storage and cpool', &
            ptr_patch=this%dispvegc_patch)

       this%storvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_STORVEGC', units='gC13/m^2', &
            avgflag='A', long_name='C13 stored vegetation carbon, excluding cpool', &
            ptr_patch=this%storvegc_patch)

       this%totvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_TOTVEGC', units='gC13/m^2', &
            avgflag='A', long_name='C13 total vegetation carbon, excluding cpool', &
            ptr_patch=this%totvegc_patch)

       this%totc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C13_TOTPFTC', units='gC13/m^2', &
            avgflag='A', long_name='C13 total patch-level carbon, including cpool', &
            ptr_patch=this%totc_patch)

       this%seedc_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_SEEDC', units='gC13/m^2', &
            avgflag='A', long_name='C13 pool for seeding new Patches', &
            ptr_col=this%seedc_col)

       this%prod10c_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_PROD10C', units='gC13/m^2', &
            avgflag='A', long_name='C13 10-yr wood product C', &
            ptr_col=this%prod10c_col)

       this%prod100c_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_PROD100C', units='gC13/m^2', &
            avgflag='A', long_name='C13 100-yr wood product C', &
            ptr_col=this%prod100c_col)

       this%totprodc_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_TOTPRODC', units='gC13/m^2', &
            avgflag='A', long_name='C13 total wood product C', &
            ptr_col=this%totprodc_col)

       this%totc_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_TOTCOLC', units='gC13/m^2', &
            avgflag='A', long_name='C13 total column carbon, incl veg and cpool', &
            ptr_col=this%totc_col)

       this%totecosysc_col(begc:endc) = spval
       call hist_addfld1d (fname='C13_TOTECOSYSC', units='gC13/m^2', &
            avgflag='A', long_name='C13 total ecosystem carbon, incl veg but excl cpool', &
            ptr_col=this%totecosysc_col)

    endif

    !-------------------------------
    ! C14 state variables 
    !-------------------------------

    if ( carbon_type == 'c14') then

       this%leafc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LEAFC', units='gC14/m^2', &
            avgflag='A', long_name='C14 leaf C', &
            ptr_patch=this%leafc_patch)

       this%leafc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LEAFC_STORAGE', units='gC14/m^2', &
            avgflag='A', long_name='C14 leaf C storage', &
            ptr_patch=this%leafc_storage_patch, default='inactive')

       this%leafc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LEAFC_XFER', units='gC14/m^2', &
            avgflag='A', long_name='C14 leaf C transfer', &
            ptr_patch=this%leafc_xfer_patch, default='inactive')

       this%frootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_FROOTC', units='gC14/m^2', &
            avgflag='A', long_name='C14 fine root C', &
            ptr_patch=this%frootc_patch)

       this%frootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_FROOTC_STORAGE', units='gC14/m^2', &
            avgflag='A', long_name='C14 fine root C storage', &
            ptr_patch=this%frootc_storage_patch, default='inactive')

       this%frootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_FROOTC_XFER', units='gC14/m^2', &
            avgflag='A', long_name='C14 fine root C transfer', &
            ptr_patch=this%frootc_xfer_patch, default='inactive')

       this%livestemc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVESTEMC', units='gC14/m^2', &
            avgflag='A', long_name='C14 live stem C', &
            ptr_patch=this%livestemc_patch)

       this%livestemc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVESTEMC_STORAGE', units='gC14/m^2', &
            avgflag='A', long_name='C14 live stem C storage', &
            ptr_patch=this%livestemc_storage_patch, default='inactive')

       this%livestemc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVESTEMC_XFER', units='gC14/m^2', &
            avgflag='A', long_name='C14 live stem C transfer', &
            ptr_patch=this%livestemc_xfer_patch, default='inactive')

       this%deadstemc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADSTEMC', units='gC14/m^2', &
            avgflag='A', long_name='C14 dead stem C', &
            ptr_patch=this%deadstemc_patch)

       this%deadstemc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADSTEMC_STORAGE', units='gC14/m^2', &
            avgflag='A', long_name='C14 dead stem C storage', &
            ptr_patch=this%deadstemc_storage_patch, default='inactive')

       this%deadstemc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADSTEMC_XFER', units='gC14/m^2', &
            avgflag='A', long_name='C14 dead stem C transfer', &
            ptr_patch=this%deadstemc_xfer_patch, default='inactive')

       this%livecrootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVECROOTC', units='gC14/m^2', &
            avgflag='A', long_name='C14 live coarse root C', &
            ptr_patch=this%livecrootc_patch)

       this%livecrootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVECROOTC_STORAGE', units='gC14/m^2', &
            avgflag='A', long_name='C14 live coarse root C storage', &
            ptr_patch=this%livecrootc_storage_patch, default='inactive')

       this%livecrootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_LIVECROOTC_XFER', units='gC14/m^2', &
            avgflag='A', long_name='C14 live coarse root C transfer', &
            ptr_patch=this%livecrootc_xfer_patch, default='inactive')

       this%deadcrootc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADCROOTC', units='gC14/m^2', &
            avgflag='A', long_name='C14 dead coarse root C', &
            ptr_patch=this%deadcrootc_patch)

       this%deadcrootc_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADCROOTC_STORAGE', units='gC14/m^2', &
            avgflag='A', long_name='C14 dead coarse root C storage', &
            ptr_patch=this%deadcrootc_storage_patch,  default='inactive')

       this%deadcrootc_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DEADCROOTC_XFER', units='gC14/m^2', &
            avgflag='A', long_name='C14 dead coarse root C transfer', &
            ptr_patch=this%deadcrootc_xfer_patch, default='inactive')

       this%gresp_storage_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_GRESP_STORAGE', units='gC14/m^2', &
            avgflag='A', long_name='C14 growth respiration storage', &
            ptr_patch=this%gresp_storage_patch, default='inactive')

       this%gresp_xfer_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_GRESP_XFER', units='gC14/m^2', &
            avgflag='A', long_name='C14 growth respiration transfer', &
            ptr_patch=this%gresp_xfer_patch, default='inactive')

       this%cpool_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_CPOOL', units='gC14/m^2', &
            avgflag='A', long_name='C14 temporary photosynthate C pool', &
            ptr_patch=this%cpool_patch)

       this%xsmrpool_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_XSMRPOOL', units='gC14/m^2', &
            avgflag='A', long_name='C14 temporary photosynthate C pool', &
            ptr_patch=this%xsmrpool_patch)

       this%ctrunc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_PFT_CTRUNC', units='gC14/m^2', &
            avgflag='A', long_name='C14 patch-level sink for C truncation', &
            ptr_patch=this%ctrunc_patch)

       this%dispvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_DISPVEGC', units='gC14/m^2', &
            avgflag='A', long_name='C14 displayed veg carbon, excluding storage and cpool', &
            ptr_patch=this%dispvegc_patch)

       this%storvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_STORVEGC', units='gC14/m^2', &
            avgflag='A', long_name='C14 stored vegetation carbon, excluding cpool', &
            ptr_patch=this%storvegc_patch)

       this%totvegc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_TOTVEGC', units='gC14/m^2', &
            avgflag='A', long_name='C14 total vegetation carbon, excluding cpool', &
            ptr_patch=this%totvegc_patch)

       this%totc_patch(begp:endp) = spval
       call hist_addfld1d (fname='C14_TOTPFTC', units='gC14/m^2', &
            avgflag='A', long_name='C14 total patch-level carbon, including cpool', &
            ptr_patch=this%totc_patch)

       this%seedc_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_SEEDC', units='gC14/m^2', &
            avgflag='A', long_name='C14 pool for seeding new Patches', &
            ptr_col=this%seedc_col)

       this%prod10c_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_PROD10C', units='gC14/m^2', &
            avgflag='A', long_name='C14 10-yr wood product C', &
            ptr_col=this%prod10c_col)

       this%prod100c_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_PROD100C', units='gC14/m^2', &
            avgflag='A', long_name='C14 100-yr wood product C', &
            ptr_col=this%prod100c_col)

       this%totprodc_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_TOTPRODC', units='gC14/m^2', &
            avgflag='A', long_name='C14 total wood product C', &
            ptr_col=this%totprodc_col)

       this%totc_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_TOTCOLC', units='gC14/m^2', &
            avgflag='A', long_name='C14 total column carbon, incl veg and cpool', &
            ptr_col=this%totc_col)

       this%totecosysc_col(begc:endc) = spval
       call hist_addfld1d (fname='C14_TOTECOSYSC', units='gC14/m^2', &
            avgflag='A', long_name='C14 total ecosystem carbon, incl veg but excl cpool', &
            ptr_col=this%totecosysc_col)

    endif

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, ratio, carbon_type, c12_cnveg_carbonstate_inst)
    !
    ! !DESCRIPTION:
    ! Initializes time varying variables used only in coupled carbon-nitrogen mode (CN):
    !
    ! !USES:
    use landunit_varcon	 , only : istsoil, istcrop 
    use clm_time_manager , only : is_restart, get_nstep
    !
    ! !ARGUMENTS:
    class(cnveg_carbonstate_type)                       :: this 
    type(bounds_type)            , intent(in)           :: bounds  
    real(r8)                     , intent(in)           :: ratio
    character(len=3)             , intent(in)           :: carbon_type ! 'c12' or 'c13' or 'c14'
    type(cnveg_carbonstate_type) , optional, intent(in) :: c12_cnveg_carbonstate_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: p,c,l,j,k,i
    integer  :: fc                                       ! filter index
    real(r8) :: c3_r1                                    ! isotope ratio (13c/12c) for C3 photosynthesis
    real(r8) :: c4_r1                                    ! isotope ratio (13c/12c) for C4 photosynthesis
    real(r8) :: c3_del13c                                ! typical del13C for C3 photosynthesis (permil, relative to PDB)
    real(r8) :: c4_del13c                                ! typical del13C for C4 photosynthesis (permil, relative to PDB)
    integer  :: num_special_col                          ! number of good values in special_col filter
    integer  :: num_special_patch                        ! number of good values in special_patch filter
    integer  :: special_col(bounds%endc-bounds%begc+1)   ! special landunit filter - columns
    integer  :: special_patch(bounds%endp-bounds%begp+1) ! special landunit filter - patches
    !-----------------------------------------------------------------------

    if (carbon_type == 'c13' .or. carbon_type == 'c14') then
       if (.not. present(c12_cnveg_carbonstate_inst)) then
          call endrun(msg=' ERROR: for C13 or C14 must pass in c12_cnveg_carbonstate_inst as argument' //&
               errMsg(__FILE__, __LINE__))
       end if
    end if

    ! Set column filters

    num_special_col = 0
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%ifspecial(l)) then
          num_special_col = num_special_col + 1
          special_col(num_special_col) = c
       end if
    end do

    ! Set patch filters

    num_special_patch = 0
    do p = bounds%begp,bounds%endp
       l = patch%landunit(p)
       if (lun%ifspecial(l)) then
          num_special_patch = num_special_patch + 1
          special_patch(num_special_patch) = p
       end if
    end do

    !-----------------------------------------------
    ! initialize patch-level carbon state variables
    !-----------------------------------------------

    do p = bounds%begp,bounds%endp

       this%leafcmax_patch(p) = 0._r8

       l = patch%landunit(p)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then

          if (patch%itype(p) == noveg) then
             this%leafc_patch(p)         = 0._r8
             this%leafc_storage_patch(p) = 0._r8
          else
             if (pftcon%evergreen(patch%itype(p)) == 1._r8) then
                this%leafc_patch(p)         = 1._r8 * ratio
                this%leafc_storage_patch(p) = 0._r8
             else if (patch%itype(p) >= npcropmin) then ! prognostic crop types
                this%leafc_patch(p) = 0._r8
                this%leafc_storage_patch(p) = 0._r8
             else
                this%leafc_patch(p) = 0._r8
                this%leafc_storage_patch(p) = 1._r8 * ratio
             end if
          end if
          this%leafc_xfer_patch(p) = 0._r8

          this%frootc_patch(p)            = 0._r8 
          this%frootc_storage_patch(p)    = 0._r8 
          this%frootc_xfer_patch(p)       = 0._r8 

          this%livestemc_patch(p)         = 0._r8 
          this%livestemc_storage_patch(p) = 0._r8 
          this%livestemc_xfer_patch(p)    = 0._r8 

          if (pftcon%woody(patch%itype(p)) == 1._r8) then
             this%deadstemc_patch(p) = 0.1_r8 * ratio
          else
             this%deadstemc_patch(p) = 0._r8 
          end if
          this%deadstemc_storage_patch(p)  = 0._r8 
          this%deadstemc_xfer_patch(p)     = 0._r8 

          this%livecrootc_patch(p)         = 0._r8 
          this%livecrootc_storage_patch(p) = 0._r8 
          this%livecrootc_xfer_patch(p)    = 0._r8 

          this%deadcrootc_patch(p)         = 0._r8 
          this%deadcrootc_storage_patch(p) = 0._r8 
          this%deadcrootc_xfer_patch(p)    = 0._r8 

          this%gresp_storage_patch(p)      = 0._r8 
          this%gresp_xfer_patch(p)         = 0._r8 

          this%cpool_patch(p)              = 0._r8 
          this%xsmrpool_patch(p)           = 0._r8 
          this%ctrunc_patch(p)             = 0._r8 
          this%dispvegc_patch(p)           = 0._r8 
          this%storvegc_patch(p)           = 0._r8 
          this%woodc_patch(p)              = 0._r8
          this%totc_patch(p)               = 0._r8 

          if ( crop_prog )then
             this%grainc_patch(p)         = 0._r8 
             this%grainc_storage_patch(p) = 0._r8 
             this%grainc_xfer_patch(p)    = 0._r8 
          end if

          ! calculate totvegc explicitly so that it is available for the isotope 
          ! code on the first time step.

          this%totvegc_patch(p) = &
               this%leafc_patch(p)              + &
               this%leafc_storage_patch(p)      + &
               this%leafc_xfer_patch(p)         + &
               this%frootc_patch(p)             + &
               this%frootc_storage_patch(p)     + &
               this%frootc_xfer_patch(p)        + &
               this%livestemc_patch(p)          + &
               this%livestemc_storage_patch(p)  + &
               this%livestemc_xfer_patch(p)     + &
               this%deadstemc_patch(p)          + &
               this%deadstemc_storage_patch(p)  + &
               this%deadstemc_xfer_patch(p)     + &
               this%livecrootc_patch(p)         + &
               this%livecrootc_storage_patch(p) + &
               this%livecrootc_xfer_patch(p)    + &
               this%deadcrootc_patch(p)         + &
               this%deadcrootc_storage_patch(p) + &
               this%deadcrootc_xfer_patch(p)    + &
               this%gresp_storage_patch(p)      + &
               this%gresp_xfer_patch(p)         + &
               this%cpool_patch(p)

          if ( crop_prog )then
             this%totvegc_patch(p) =         &
                  this%totvegc_patch(p)    + &
                  this%grainc_patch(p)         + &
                  this%grainc_storage_patch(p) + &
                  this%grainc_xfer_patch(p)
          end if

       endif

    end do

    ! -----------------------------------------------
    ! initialize column-level variables
    ! -----------------------------------------------

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          ! dynamic landcover state variables
          this%seedc_col(c)      = 0._r8
          this%prod10c_col(c)    = 0._r8
          this%prod100c_col(c)   = 0._r8
          this%totprodc_col(c)   = 0._r8

          ! total carbon pools
          this%totecosysc_col(c) = 0._r8
          this%totc_col(c)       = 0._r8
       end if
    end do

    ! now loop through special filters and explicitly set the variables that
    ! have to be in place for biogeophysics
    
    do fc = 1,num_special_col
       c = special_col(fc)
       this%seedc_col(c)      = 0._r8
       this%prod10c_col(c)    = 0._r8
       this%prod100c_col(c)   = 0._r8
       this%totprodc_col(c)   = 0._r8
    end do

    if ( .not. is_restart() .and. get_nstep() == 1 ) then
       c3_del13c = -28._r8
       c4_del13c = -13._r8
       c3_r1 = SHR_CONST_PDB + ((c3_del13c*SHR_CONST_PDB)/1000._r8)
       c3_r2 = c3_r1/(1._r8 + c3_r1)
       c4_r1 = SHR_CONST_PDB + ((c4_del13c*SHR_CONST_PDB)/1000._r8)
       c4_r2 = c4_r1/(1._r8 + c4_r1)

       do p = bounds%begp,bounds%endp
          if (pftcon%c3psn(patch%itype(p)) == 1._r8) then
             this%grainc_patch(p)            = c12_cnveg_carbonstate_inst%grainc_patch(p)         * c3_r2
             this%grainc_storage_patch(p)    = c12_cnveg_carbonstate_inst%grainc_storage_patch(p) * c3_r2
             this%grainc_xfer_patch(p)       = c12_cnveg_carbonstate_inst%grainc_xfer_patch(p)    * c3_r2
             this%dispvegc_patch(p)          = c12_cnveg_carbonstate_inst%dispvegc_patch(p)       * c3_r2
             this%storvegc_patch(p)          = c12_cnveg_carbonstate_inst%storvegc_patch(p)       * c3_r2
             this%totvegc_patch(p)           = c12_cnveg_carbonstate_inst%totvegc_patch(p)        * c3_r2
             this%totc_patch(p)              = c12_cnveg_carbonstate_inst%totc_patch(p)           * c3_r2
             this%woodc_patch(p)             = c12_cnveg_carbonstate_inst%woodc_patch(p)          * c3_r2
          else
             this%grainc_patch(p)            = c12_cnveg_carbonstate_inst%grainc_patch(p)         * c4_r2
             this%grainc_storage_patch(p)    = c12_cnveg_carbonstate_inst%grainc_storage_patch(p) * c4_r2
             this%grainc_xfer_patch(p)       = c12_cnveg_carbonstate_inst%grainc_xfer_patch(p)    * c4_r2
             this%dispvegc_patch(p)          = c12_cnveg_carbonstate_inst%dispvegc_patch(p)       * c4_r2
             this%storvegc_patch(p)          = c12_cnveg_carbonstate_inst%storvegc_patch(p)       * c4_r2
             this%totvegc_patch(p)           = c12_cnveg_carbonstate_inst%totvegc_patch(p)        * c4_r2
             this%totc_patch(p)              = c12_cnveg_carbonstate_inst%totc_patch(p)           * c4_r2
             this%woodc_patch(p)             = c12_cnveg_carbonstate_inst%woodc_patch(p)          * c4_r2
          end if
       end do
    end if

    ! initialize fields for special filters

    call this%SetValues (&
         num_patch=num_special_patch, filter_patch=special_patch, value_patch=0._r8, &
         num_column=num_special_col, filter_column=special_col, value_column=0._r8)

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine Restart ( this,  bounds, ncid, flag, carbon_type, c12_cnveg_carbonstate_inst)
    !
    ! !DESCRIPTION: 
    ! Read/write CN restart data for carbon state
    !
    ! !USES:
    use shr_infnan_mod   , only : isnan => shr_infnan_isnan, nan => shr_infnan_nan, assignment(=)
    use clm_varcon       , only : c14ratio
    use restUtilMod
    use ncdio_pio
    !
    ! !ARGUMENTS:
    class (cnveg_carbonstate_type)                               :: this
    type(bounds_type)                     , intent(in)           :: bounds 
    type(file_desc_t)                     , intent(inout)        :: ncid   ! netcdf id
    character(len=*)                      , intent(in)           :: flag   !'read' or 'write'
    character(len=3)                      , intent(in)           :: carbon_type ! 'c12' or 'c13' or 'c14'
    type (cnveg_carbonstate_type)         , intent(in), optional :: c12_cnveg_carbonstate_inst 
    !
    ! !LOCAL VARIABLES:
    integer            :: i,j,k,l,c
    real(r8)           :: m         ! multiplier for the exit_spinup code
    character(len=128) :: varname   ! temporary
    logical            :: readvar
    integer            :: idata
    logical            :: exit_spinup  = .false.
    logical            :: enter_spinup = .false.
    real(r8)           :: c3_del13c ! typical del13C for C3 photosynthesis (permil, relative to PDB)
    real(r8)           :: c4_del13c ! typical del13C for C4 photosynthesis (permil, relative to PDB)
    real(r8)           :: c3_r1     ! isotope ratio (13c/12c) for C3 photosynthesis
    real(r8)           :: c4_r1     ! isotope ratio (13c/12c) for C4 photosynthesis
    real(r8)           :: c3_r2     ! isotope ratio (13c/[12c+13c]) for C3 photosynthesis
    real(r8)           :: c4_r2     ! isotope ratio (13c/[12c+13c]) for C4 photosynthesis
    ! flags for comparing the model and restart decomposition cascades
    integer            :: decomp_cascade_state, restart_file_decomp_cascade_state 
    !------------------------------------------------------------------------

    if (carbon_type == 'c13' .or. carbon_type == 'c14') then
       if (.not. present(c12_cnveg_carbonstate_inst)) then
          call endrun(msg=' ERROR: for C14 must pass in c12_cnveg_carbonstate_inst as argument' //&
               errMsg(__FILE__, __LINE__))
       end if
    end if

    !--------------------------------
    ! patch carbon state variables (c12)
    !--------------------------------

    if (carbon_type == 'c12') then
       call restartvar(ncid=ncid, flag=flag, varname='leafc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='leafc_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_storage_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='leafc_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_xfer_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='frootc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='frootc_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_storage_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='frootc_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_xfer_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='livestemc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='livestemc_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_storage_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='livestemc_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_xfer_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_storage_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_xfer_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_storage_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_xfer_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_storage_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_xfer_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='gresp_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%gresp_storage_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='gresp_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%gresp_xfer_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='cpool', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%cpool_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='xsmrpool', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%xsmrpool_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='pft_ctrunc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%ctrunc_patch) 

       call restartvar(ncid=ncid, flag=flag, varname='leafcmax', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafcmax_patch)

       call restartvar(ncid=ncid, flag=flag, varname='totcolc', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totc_col) 
    end if

    !--------------------------------
    ! C13 patch carbon state variables 
    !--------------------------------

    if ( carbon_type == 'c13')  then
       call restartvar(ncid=ncid, flag=flag, varname='leafc_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_patch)
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%leafc with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%leafc_patch(i) = c12_cnveg_carbonstate_inst%leafc_patch(i) * c3_r2
             else
                this%leafc_patch(i) = c12_cnveg_carbonstate_inst%leafc_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='leafc_storage_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%leafc_storage with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%leafc_storage_patch(i) = c12_cnveg_carbonstate_inst%leafc_storage_patch(i) * c3_r2
             else
                this%leafc_storage_patch(i) = c12_cnveg_carbonstate_inst%leafc_storage_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='leafc_xfer_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%leafc_xfer with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%leafc_xfer_patch(i) = c12_cnveg_carbonstate_inst%leafc_xfer_patch(i) * c3_r2
             else
                this%leafc_xfer_patch(i) = c12_cnveg_carbonstate_inst%leafc_xfer_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='frootc_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%frootc with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%frootc_patch(i) = c12_cnveg_carbonstate_inst%frootc_patch(i) * c3_r2
             else
                this%frootc_patch(i) = c12_cnveg_carbonstate_inst%frootc_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='frootc_storage_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%frootc_storage with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%frootc_storage_patch(i) = c12_cnveg_carbonstate_inst%frootc_storage_patch(i) * c3_r2
             else
                this%frootc_storage_patch(i) = c12_cnveg_carbonstate_inst%frootc_storage_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='frootc_xfer_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%frootc_xfer with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%frootc_xfer_patch(i) = c12_cnveg_carbonstate_inst%frootc_xfer_patch(i) * c3_r2
             else
                this%frootc_xfer_patch(i) = c12_cnveg_carbonstate_inst%frootc_xfer_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livestemc_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%livestemc with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%livestemc_patch(i) = c12_cnveg_carbonstate_inst%livestemc_patch(i) * c3_r2
             else
                this%livestemc_patch(i) = c12_cnveg_carbonstate_inst%livestemc_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livestemc_storage_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%livestemc_storage with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%livestemc_storage_patch(i) = c12_cnveg_carbonstate_inst%livestemc_storage_patch(i) * c3_r2
             else
                this%livestemc_storage_patch(i) = c12_cnveg_carbonstate_inst%livestemc_storage_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livestemc_xfer_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%livestemc_xfer with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%livestemc_xfer_patch(i) = c12_cnveg_carbonstate_inst%livestemc_xfer_patch(i) * c3_r2
             else
                this%livestemc_xfer_patch(i) = c12_cnveg_carbonstate_inst%livestemc_xfer_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%deadstemc with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%deadstemc_patch(i) = c12_cnveg_carbonstate_inst%deadstemc_patch(i) * c3_r2
             else
                this%deadstemc_patch(i) = c12_cnveg_carbonstate_inst%deadstemc_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_storage_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%deadstemc_storage with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%deadstemc_storage_patch(i) = c12_cnveg_carbonstate_inst%deadstemc_storage_patch(i) * c3_r2
             else
                this%deadstemc_storage_patch(i) = c12_cnveg_carbonstate_inst%deadstemc_storage_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_xfer_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%deadstemc_xfer with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%deadstemc_xfer_patch(i) = c12_cnveg_carbonstate_inst%deadstemc_xfer_patch(i) * c3_r2
             else
                this%deadstemc_xfer_patch(i) = c12_cnveg_carbonstate_inst%deadstemc_xfer_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%livecrootc with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%livecrootc_patch(i) = c12_cnveg_carbonstate_inst%livecrootc_patch(i) * c3_r2
             else
                this%livecrootc_patch(i) = c12_cnveg_carbonstate_inst%livecrootc_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_storage_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%livecrootc_storage with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%livecrootc_storage_patch(i) = c12_cnveg_carbonstate_inst%livecrootc_storage_patch(i) * c3_r2
             else
                this%livecrootc_storage_patch(i) = c12_cnveg_carbonstate_inst%livecrootc_storage_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_xfer_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%livecrootc_xfer with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%livecrootc_xfer_patch(i) = c12_cnveg_carbonstate_inst%livecrootc_xfer_patch(i) * c3_r2
             else
                this%livecrootc_xfer_patch(i) = c12_cnveg_carbonstate_inst%livecrootc_xfer_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%deadcrootc with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%deadcrootc_patch(i) = c12_cnveg_carbonstate_inst%deadcrootc_patch(i) * c3_r2
             else
                this%deadcrootc_patch(i) = c12_cnveg_carbonstate_inst%deadcrootc_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_storage_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%deadcrootc_storage with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%deadcrootc_storage_patch(i) = c12_cnveg_carbonstate_inst%deadcrootc_storage_patch(i) * c3_r2
             else
                this%deadcrootc_storage_patch(i) = c12_cnveg_carbonstate_inst%deadcrootc_storage_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_xfer_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%deadcrootc_xfer with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%deadcrootc_xfer_patch(i) = c12_cnveg_carbonstate_inst%deadcrootc_xfer_patch(i) * c3_r2
             else
                this%deadcrootc_xfer_patch(i) = c12_cnveg_carbonstate_inst%deadcrootc_xfer_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='gresp_storage_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%gresp_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%gresp_storage with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%gresp_storage_patch(i) = c12_cnveg_carbonstate_inst%gresp_storage_patch(i) * c3_r2
             else
                this%gresp_storage_patch(i) = c12_cnveg_carbonstate_inst%gresp_storage_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='gresp_xfer_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%gresp_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%gresp_xfer with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%gresp_xfer_patch(i) = c12_cnveg_carbonstate_inst%gresp_xfer_patch(i) * c3_r2
             else
                this%gresp_xfer_patch(i) = c12_cnveg_carbonstate_inst%gresp_xfer_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='cpool_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%cpool_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%cpool with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%cpool_patch(i) = c12_cnveg_carbonstate_inst%cpool_patch(i) * c3_r2
             else
                this%cpool_patch(i) = c12_cnveg_carbonstate_inst%cpool_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='xsmrpool_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%xsmrpool_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%xsmrpool with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%xsmrpool_patch(i) = c12_cnveg_carbonstate_inst%xsmrpool_patch(i) * c3_r2
             else
                this%xsmrpool_patch(i) = c12_cnveg_carbonstate_inst%xsmrpool_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='pft_ctrunc_13', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%ctrunc_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%ctrunc with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%ctrunc_patch(i) = c12_cnveg_carbonstate_inst%ctrunc_patch(i) * c3_r2
             else
                this%ctrunc_patch(i) = c12_cnveg_carbonstate_inst%ctrunc_patch(i) * c4_r2
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='totcolc_13', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totc_col) 
    end if

    !--------------------------------
    ! C14 patch carbon state variables 
    !--------------------------------

    if ( carbon_type == 'c14')  then
       call restartvar(ncid=ncid, flag=flag, varname='leafc_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%leafc_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%leafc_patch(i) /= spval .and. &
                  .not. isnan(this%leafc_patch(i)) ) then
                this%leafc_patch(i) = c12_cnveg_carbonstate_inst%leafc_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='leafc_storage_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%leafc_storage_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%leafc_storage_patch(i) /= spval .and. &
                  .not. isnan(this%leafc_storage_patch(i)) ) then
                this%leafc_storage_patch(i) = c12_cnveg_carbonstate_inst%leafc_storage_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='leafc_xfer_14', xtype=ncd_double,  &
            dim1name='pft',    long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%leafc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%leafc_xfer_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%leafc_xfer_patch(i) /= spval .and. .not. isnan(this%leafc_xfer_patch(i)) ) then
                this%leafc_xfer_patch(i) = c12_cnveg_carbonstate_inst%leafc_xfer_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='frootc_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%frootc_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%frootc_patch(i) /= spval .and. &
                  .not. isnan(this%frootc_patch(i)) ) then
                this%frootc_patch(i) = c12_cnveg_carbonstate_inst%frootc_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='frootc_storage_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%frootc_storage_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%frootc_storage_patch(i) /= spval .and. &
                  .not. isnan(this%frootc_storage_patch(i)) ) then
                this%frootc_storage_patch(i) = c12_cnveg_carbonstate_inst%frootc_storage_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='frootc_xfer_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%frootc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%frootc_xfer_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%frootc_xfer_patch(i) /= spval .and. &
                  .not. isnan(this%frootc_xfer_patch(i)) ) then
                this%frootc_xfer_patch(i) = c12_cnveg_carbonstate_inst%frootc_xfer_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livestemc_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%livestemc_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%livestemc_patch(i) /= spval .and. .not. isnan(this%livestemc_patch(i)) ) then
                this%livestemc_patch(i) = c12_cnveg_carbonstate_inst%livestemc_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livestemc_storage_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%livestemc_storage_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%livestemc_storage_patch(i) /= spval .and. .not. isnan(this%livestemc_storage_patch(i)) ) then
                this%livestemc_storage_patch(i) = c12_cnveg_carbonstate_inst%livestemc_storage_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livestemc_xfer_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%livestemc_xfer_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%livestemc_xfer_patch(i) /= spval .and. .not. isnan(this%livestemc_xfer_patch(i)) ) then
                this%livestemc_xfer_patch(i) = c12_cnveg_carbonstate_inst%livestemc_xfer_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%deadstemc_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%deadstemc_patch(i) /= spval .and. .not. isnan(this%deadstemc_patch(i)) ) then
                this%deadstemc_patch(i) = c12_cnveg_carbonstate_inst%deadstemc_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_storage_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%deadstemc_storage_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%deadstemc_storage_patch(i) /= spval .and. .not. isnan(this%deadstemc_storage_patch(i)) ) then
                this%deadstemc_storage_patch(i) = c12_cnveg_carbonstate_inst%deadstemc_storage_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_xfer_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadstemc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%deadstemc_xfer_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%deadstemc_xfer_patch(i) /= spval .and. .not. isnan(this%deadstemc_xfer_patch(i)) ) then
                this%deadstemc_xfer_patch(i) = c12_cnveg_carbonstate_inst%deadstemc_xfer_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%livecrootc_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%livecrootc_patch(i) /= spval .and. .not. isnan(this%livecrootc_patch(i)) ) then
                this%livecrootc_patch(i) = c12_cnveg_carbonstate_inst%livecrootc_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_storage_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%livecrootc_storage_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%livecrootc_storage_patch(i) /= spval .and. .not. isnan(this%livecrootc_storage_patch(i)) ) then
                this%livecrootc_storage_patch(i) = c12_cnveg_carbonstate_inst%livecrootc_storage_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_xfer_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%livecrootc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%livecrootc_xfer_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%livecrootc_xfer_patch(i) /= spval .and. .not. isnan(this%livecrootc_xfer_patch(i)) ) then
                this%livecrootc_xfer_patch(i) = c12_cnveg_carbonstate_inst%livecrootc_xfer_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%deadcrootc_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%deadcrootc_patch(i) /= spval .and. .not. isnan(this%deadcrootc_patch(i)) ) then
                this%deadcrootc_patch(i) = c12_cnveg_carbonstate_inst%deadcrootc_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_storage_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%deadcrootc_storage_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%deadcrootc_storage_patch(i) /= spval .and. .not. isnan(this%deadcrootc_storage_patch(i)) ) then
                this%deadcrootc_storage_patch(i) = c12_cnveg_carbonstate_inst%deadcrootc_storage_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_xfer_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%deadcrootc_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog) 'initializing this%deadcrootc_xfer_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%deadcrootc_xfer_patch(i) /= spval .and. .not. isnan(this%deadcrootc_xfer_patch(i)) ) then
                this%deadcrootc_xfer_patch(i) = c12_cnveg_carbonstate_inst%deadcrootc_xfer_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='gresp_storage_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%gresp_storage_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%gresp_storage_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%gresp_storage_patch(i) /= spval .and. .not. isnan(this%gresp_storage_patch(i)) ) then
                this%gresp_storage_patch(i) = c12_cnveg_carbonstate_inst%gresp_storage_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='gresp_xfer_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%gresp_xfer_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%gresp_xfer_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%gresp_xfer_patch(i) /= spval .and. .not. isnan(this%gresp_xfer_patch(i)) ) then
                this%gresp_xfer_patch(i) = c12_cnveg_carbonstate_inst%gresp_xfer_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='cpool_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%cpool_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%cpool_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%cpool_patch(i) /= spval .and. .not. isnan(this%cpool_patch(i)) ) then
                this%cpool_patch(i) = c12_cnveg_carbonstate_inst%cpool_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='xsmrpool_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%xsmrpool_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%xsmrpool_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%xsmrpool_patch(i) /= spval .and. .not. isnan(this%xsmrpool_patch(i)) ) then
                this%xsmrpool_patch(i) = c12_cnveg_carbonstate_inst%xsmrpool_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='pft_ctrunc_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%ctrunc_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%ctrunc_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%ctrunc_patch(i) /= spval .and. .not. isnan(this%ctrunc_patch(i)) ) then
                this%ctrunc_patch(i) = c12_cnveg_carbonstate_inst%ctrunc_patch(i) * c14ratio
             endif
          end do
       end if

       call restartvar(ncid=ncid, flag=flag, varname='totcolc_14', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totc_col) 
    end if

    !--------------------------------
    ! patch prognostic crop variables
    !--------------------------------

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag,  varname='grainc', xtype=ncd_double,  &
            dim1name='pft', long_name='grain C', units='gC/m2', &
            interpinic_flag='interp', readvar=readvar, data=this%grainc_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='grainc_storage', xtype=ncd_double,  &
            dim1name='pft', long_name='grain C storage', units='gC/m2', &
            interpinic_flag='interp', readvar=readvar, data=this%grainc_storage_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='grainc_xfer', xtype=ncd_double,  &
            dim1name='pft', long_name='grain C transfer', units='gC/m2', &
            interpinic_flag='interp', readvar=readvar, data=this%grainc_xfer_patch)
    end if

    !--------------------------------
    ! column carbon state variables
    !--------------------------------

    if (carbon_type == 'c12') then
       call restartvar(ncid=ncid, flag=flag, varname='seedc', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%seedc_col) 
    end if

    if (carbon_type == 'c12') then
       call restartvar(ncid=ncid, flag=flag, varname='prod10c', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod10c_col) 
    end if

    if (carbon_type == 'c12') then
       call restartvar(ncid=ncid, flag=flag, varname='prod100c', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod100c_col) 
    end if

    !--------------------------------
    ! C13 column carbon state variables
    !--------------------------------

    if (carbon_type == 'c13') then
       call restartvar(ncid=ncid, flag=flag, varname='seedc_13', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%seedc_col) 
       if (flag=='read' .and. .not. readvar) then
          if (this%seedc_col(i) /= spval .and. &
               .not. isnan(this%seedc_col(i)) ) then
             this%seedc_col(i) = c12_cnveg_carbonstate_inst%seedc_col(i) * c3_r2
          end if
       end if
    end if

    if (carbon_type == 'c13') then
       call restartvar(ncid=ncid, flag=flag, varname='prod10c_13', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod10c_col) 
       if (flag=='read' .and. .not. readvar) then
          if (this%prod10c_col(i) /= spval .and. &
               .not. isnan( this%prod10c_col(i) ) ) then
             this%prod10c_col(i) = c12_cnveg_carbonstate_inst%prod10c_col(i) * c3_r2
          endif
       end if
    end if

    if (carbon_type == 'c13') then
       call restartvar(ncid=ncid, flag=flag, varname='prod100c_13', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod100c_col) 
       if (flag=='read' .and. .not. readvar) then
          if (this%prod100c_col(i) /= spval .and. &
               .not. isnan( this%prod100c_col(i) ) ) then
             this%prod100c_col(i) = c12_cnveg_carbonstate_inst%prod100c_col(i) * c3_r2
          endif
       end if
    endif

    !--------------------------------
    ! C14 column carbon state variables
    !--------------------------------

    if ( carbon_type == 'c14' ) then
       call restartvar(ncid=ncid, flag=flag, varname='seedc_14', xtype=ncd_double,  &
            dim1name='column', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%seedc_col) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%seedc_col with atmospheric c14 value'
          do i = bounds%begc,bounds%endc
             if (this%seedc_col(i) /= spval .and. &
                  .not. isnan(this%seedc_col(i)) ) then
                this%seedc_col(i) = c12_cnveg_carbonstate_inst%seedc_col(i) * c14ratio
             endif
          end do
       end if
    end if

    if ( carbon_type == 'c14' ) then
       call restartvar(ncid=ncid, flag=flag, varname='prod10c_14', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod10c_col) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%prod10c_col with atmospheric c14 value'
          if (this%prod10c_col(i) /= spval .and. &
               .not. isnan(this%prod10c_col(i)) ) then
             this%prod10c_col(i) = c12_cnveg_carbonstate_inst%prod10c_col(i) * c14ratio
          endif
       end if
    end if

    if ( carbon_type == 'c14' ) then
       call restartvar(ncid=ncid, flag=flag, varname='prod100c_14', xtype=ncd_double,  &
            dim1name='column', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%prod100c_col) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%prod100c_col with atmospheric c14 value'
          if (this%prod100c_col(i) /= spval .and. &
               .not. isnan(this%prod100c_col(i)) ) then
             this%prod100c_col(i) = c12_cnveg_carbonstate_inst%prod100c_col(i) * c14ratio
          endif
       end if
    endif

    if (carbon_type == 'c13' .or. carbon_type == 'c14') then
       if (.not. present(c12_cnveg_carbonstate_inst)) then
          call endrun(msg=' ERROR: for C13 or C14 must pass in c12_cnveg_carbonstate_inst as argument' //&
               errMsg(__FILE__, __LINE__))
       end if
    end if

    c3_del13c = -28._r8
    c3_r1 = SHR_CONST_PDB + ((c3_del13c*SHR_CONST_PDB)/1000._r8)
    c3_r2 = c3_r1/(1._r8 + c3_r1)

    c4_del13c = -13._r8
    c4_r1 = SHR_CONST_PDB + ((c4_del13c*SHR_CONST_PDB)/1000._r8)
    c4_r2 = c4_r1/(1._r8 + c4_r1)

    !--------------------------------
    ! C12 carbon state variables
    !--------------------------------

    if (carbon_type == 'c12') then
       call restartvar(ncid=ncid, flag=flag, varname='totvegc', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totvegc_patch) 
    end if

    !--------------------------------
    ! C13 carbon state variables 
    !--------------------------------

    if ( carbon_type == 'c13')  then
       call restartvar(ncid=ncid, flag=flag, varname='totvegc_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totvegc_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing cnveg_carbonstate_inst%totvegc with atmospheric c13 value'
          do i = bounds%begp,bounds%endp
             if (pftcon%c3psn(patch%itype(i)) == 1._r8) then
                this%totvegc_patch(i) = c12_cnveg_carbonstate_inst%totvegc_patch(i) * c3_r2
             else
                this%totvegc_patch(i) = c12_cnveg_carbonstate_inst%totvegc_patch(i) * c4_r2
             endif
          end do
       end if
    endif

    !--------------------------------
    ! C14 patch carbon state variables 
    !--------------------------------

    if ( carbon_type == 'c14')  then
       call restartvar(ncid=ncid, flag=flag, varname='totvegc_14', xtype=ncd_double,  &
            dim1name='pft', long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%totvegc_patch) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing this%totvegc_patch with atmospheric c14 value'
          do i = bounds%begp,bounds%endp
             if (this%totvegc_patch(i) /= spval .and. &
                  .not. isnan(this%totvegc_patch(i)) ) then
                this%totvegc_patch(i) = c12_cnveg_carbonstate_inst%totvegc_patch(i) * c14ratio
             endif
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
    class (cnveg_carbonstate_type) :: this
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

    do fi = 1,num_patch
       i  = filter_patch(fi)
       this%leafc_patch(i)              = value_patch
       this%leafc_storage_patch(i)      = value_patch
       this%leafc_xfer_patch(i)         = value_patch
       this%frootc_patch(i)             = value_patch
       this%frootc_storage_patch(i)     = value_patch
       this%frootc_xfer_patch(i)        = value_patch
       this%livestemc_patch(i)          = value_patch
       this%livestemc_storage_patch(i)  = value_patch
       this%livestemc_xfer_patch(i)     = value_patch
       this%deadstemc_patch(i)          = value_patch
       this%deadstemc_storage_patch(i)  = value_patch
       this%deadstemc_xfer_patch(i)     = value_patch
       this%livecrootc_patch(i)         = value_patch
       this%livecrootc_storage_patch(i) = value_patch
       this%livecrootc_xfer_patch(i)    = value_patch
       this%deadcrootc_patch(i)         = value_patch
       this%deadcrootc_storage_patch(i) = value_patch
       this%deadcrootc_xfer_patch(i)    = value_patch
       this%gresp_storage_patch(i)      = value_patch
       this%gresp_xfer_patch(i)         = value_patch
       this%cpool_patch(i)              = value_patch
       this%xsmrpool_patch(i)           = value_patch
       this%ctrunc_patch(i)             = value_patch
       this%dispvegc_patch(i)           = value_patch
       this%storvegc_patch(i)           = value_patch
       this%woodc_patch(i)              = value_patch
       this%totvegc_patch(i)            = value_patch
       this%totc_patch(i)               = value_patch
       if ( crop_prog ) then
          this%grainc_patch(i)          = value_patch
          this%grainc_storage_patch(i)  = value_patch
          this%grainc_xfer_patch(i)     = value_patch
       end if
    end do

    do fi = 1,num_column
       i  = filter_column(fi)
       this%rootc_col(i)                = value_column
       this%leafc_col(i)                = value_column
       this%fuelc_col(i)                = value_column
       this%fuelc_crop_col(i)           = value_column
       this%totvegc_col(i)              = value_column
       this%totc_col(i)                 = value_column
       this%totecosysc_col(i)           = value_column
    end do

  end subroutine SetValues

  !-----------------------------------------------------------------------
  subroutine ZeroDwt( this, bounds )
    !
    ! !DESCRIPTION
    ! Initialize variables needed for dynamic land use.
    !
    ! !ARGUMENTS:
    class(cnveg_carbonstate_type) :: this
    type(bounds_type), intent(in)  :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer  :: p          ! indices
    !-----------------------------------------------------------------------

    do p = bounds%begp,bounds%endp
       this%dispvegc_patch(p)   = 0._r8
       this%storvegc_patch(p)   = 0._r8
       this%totc_patch(p)       = 0._r8
    end do

  end subroutine ZeroDwt

  !-----------------------------------------------------------------------
  subroutine Summary_carbonstate(this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       soilbiogeochem_cwdc_col, soilbiogeochem_totlitc_col, soilbiogeochem_totsomc_col, &
       soilbiogeochem_ctrunc_col)
    !
    ! !USES:
    use subgridAveMod, only : p2c
    !
    ! !DESCRIPTION:
    ! On the radiation time step, perform patch and column-level carbon summary calculations
    !
    ! !ARGUMENTS:
    class(cnveg_carbonstate_type)  :: this
    type(bounds_type) , intent(in) :: bounds          
    integer           , intent(in) :: num_soilc       ! number of soil columns in filter
    integer           , intent(in) :: filter_soilc(:) ! filter for soil columns
    integer           , intent(in) :: num_soilp       ! number of soil patches in filter
    integer           , intent(in) :: filter_soilp(:) ! filter for soil patches
    real(r8)          , intent(in) :: soilbiogeochem_cwdc_col(bounds%begc:)   
    real(r8)          , intent(in) :: soilbiogeochem_totlitc_col(bounds%begc:)
    real(r8)          , intent(in) :: soilbiogeochem_totsomc_col(bounds%begc:)
    real(r8)          , intent(in) :: soilbiogeochem_ctrunc_col(bounds%begc:)
    !
    ! !LOCAL VARIABLES:
    integer  :: c,p,j,k,l       ! indices
    integer  :: fp,fc           ! lake filter indices
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(soilbiogeochem_cwdc_col)    == (/bounds%endc/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(soilbiogeochem_totlitc_col) == (/bounds%endc/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(soilbiogeochem_totsomc_col) == (/bounds%endc/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(soilbiogeochem_ctrunc_col)  == (/bounds%endc/)), errMsg(__FILE__, __LINE__))

    ! calculate patch -level summary of carbon state

    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! displayed vegetation carbon, excluding storage and cpool (DISPVEGC)
       this%dispvegc_patch(p) =        &
            this%leafc_patch(p)      + &
            this%frootc_patch(p)     + &
            this%livestemc_patch(p)  + &
            this%deadstemc_patch(p)  + &
            this%livecrootc_patch(p) + &
            this%deadcrootc_patch(p)

       ! stored vegetation carbon, excluding cpool (STORVEGC)
       this%storvegc_patch(p) =                &
            this%cpool_patch(p)              + &
            this%leafc_storage_patch(p)      + &
            this%frootc_storage_patch(p)     + &
            this%livestemc_storage_patch(p)  + &
            this%deadstemc_storage_patch(p)  + &
            this%livecrootc_storage_patch(p) + &
            this%deadcrootc_storage_patch(p) + &
            this%leafc_xfer_patch(p)         + &
            this%frootc_xfer_patch(p)        + &
            this%livestemc_xfer_patch(p)     + &
            this%deadstemc_xfer_patch(p)     + &
            this%livecrootc_xfer_patch(p)    + &
            this%deadcrootc_xfer_patch(p)    + &
            this%gresp_storage_patch(p)      + &
            this%gresp_xfer_patch(p)

       if ( crop_prog .and. patch%itype(p) >= npcropmin )then
          this%storvegc_patch(p) =            &
               this%storvegc_patch(p)       + &
               this%grainc_storage_patch(p) + &
               this%grainc_xfer_patch(p)

          this%dispvegc_patch(p) =            &
               this%dispvegc_patch(p)       + &
               this%grainc_patch(p)
       end if

       ! total vegetation carbon, excluding cpool (TOTVEGC)
       this%totvegc_patch(p) = &
            this%dispvegc_patch(p) + &
            this%storvegc_patch(p)

       ! total patch-level carbon, including xsmrpool, ctrunc
       this%totc_patch(p) = &
            this%totvegc_patch(p) + &
            this%xsmrpool_patch(p) + &
            this%ctrunc_patch(p)

       ! (WOODC) - wood C
       this%woodc_patch(p) = &
            this%deadstemc_patch(p)    + &
            this%livestemc_patch(p)    + &
            this%deadcrootc_patch(p)   + &
            this%livecrootc_patch(p)

    end do

    ! --------------------------------------------
    ! column level summary
    ! --------------------------------------------

    call p2c(bounds, num_soilc, filter_soilc, &
         this%totvegc_patch(bounds%begp:bounds%endp), &
         this%totvegc_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%totc_patch(bounds%begp:bounds%endp), &
         this%totc_col(bounds%begc:bounds%endc))

    do fc = 1,num_soilc
       c = filter_soilc(fc)

       ! total wood product carbon
       this%totprodc_col(c) =      &
            this%prod10c_col(c)  + &
            this%prod100c_col(c)	 

       ! total ecosystem carbon, including veg but excluding cpool (TOTECOSYSC)
       this%totecosysc_col(c) =    &
            soilbiogeochem_cwdc_col(c)    + &
            soilbiogeochem_totlitc_col(c) + &
            soilbiogeochem_totsomc_col(c) + &
            this%totprodc_col(c)          + &
            this%totvegc_col(c)

       ! total column carbon, including veg and cpool (TOTCOLC)
       this%totc_col(c) =  this%totc_col(c) + &
            soilbiogeochem_cwdc_col(c)      + &
            soilbiogeochem_totlitc_col(c)   + &
            soilbiogeochem_totsomc_col(c)   + &
            this%totprodc_col(c)            + &
            this%seedc_col(c)               + &
            soilbiogeochem_ctrunc_col(c)

    end do

  end subroutine Summary_carbonstate

end module CNVegCarbonStateType
