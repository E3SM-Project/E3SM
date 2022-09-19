module ocn2lndType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle ocn2lnd, lnd2ocn mapping
  !
  ! !USES:
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod   , only : errMsg => shr_log_errMsg
  use shr_megan_mod , only : shr_megan_mechcomps_n
  use elm_varpar    , only : numrad, ndst, nlevgrnd !ndst = number of dust bins.
  use elm_varcon    , only : rair, grav, cpair, hfus, tfrz, spval
  use elm_varctl    , only : iulog
  use seq_drydep_mod, only : n_drydep, drydep_method, DD_XLND
  use decompMod     , only : bounds_type
  use abortutils    , only : endrun
  use seq_flds_mod  , only : ocn_lnd_one_way
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  ! !PUBLIC DATA TYPES
  !----------------------------------------------------
  ! ocean -> land variables structure
  !
  type, public :: ocn2lnd_type

      real(r8), pointer :: frac_h2oocn_grc(:)   => null() ! sea surface height [m]

    contains

      procedure, public  :: Init
      procedure, private :: InitAllocate 
      procedure, private :: InitHistory  
      procedure, public  :: Restart

   end type ocn2lnd_type

  contains

    !------------------------------------------------------------------------
    subroutine Init(this, bounds)

      class(ocn2lnd_type) :: this
      type(bounds_type), intent(in) :: bounds  

      call this%InitAllocate(bounds)
      call this%InitHistory(bounds)

    end subroutine Init

    !------------------------------------------------------------------------
    subroutine InitAllocate(this, bounds)
      !
      ! !DESCRIPTION:
      ! Initialize ocn2lnd derived type
      !
      ! !ARGUMENTS:
      class(ocn2lnd_type) :: this
      type(bounds_type), intent(in) :: bounds  
      !
      ! !LOCAL VARIABLES:
      real(r8) :: ival  = 0.0_r8  ! initial value
      real     :: ival_float = 0.0
      integer  :: ival_int = 0
      integer*2 :: ival_short = 0
      integer  :: begg, endg
      integer  :: begc, endc
      integer  :: begp, endp
      !------------------------------------------------------------------------

      begg = bounds%begg; endg= bounds%endg
      begc = bounds%begc; endc= bounds%endc
      begp = bounds%begp; endp= bounds%endp

      allocate(this%frac_h2oocn_grc(begg:endg)); this%frac_h2oocn_grc(:)   = ival

    end subroutine InitAllocate

    !------------------------------------------------------------------------
    subroutine InitHistory(this, bounds)
      !
      ! !USES:
      use histFileMod, only : hist_addfld1d
      !
      ! !ARGUMENTS:
      class(ocn2lnd_type) :: this
      type(bounds_type), intent(in) :: bounds  
      !
      ! !LOCAL VARIABLES:
      integer  :: begg, endg
      integer  :: begp, endp
      !---------------------------------------------------------------------

      begg = bounds%begg; endg= bounds%endg
      begp = bounds%begp; endp= bounds%endp

      if (ocn_lnd_one_way) then
         this%frac_h2oocn_grc(begg:endg) = spval
         call hist_addfld1d (fname='FRAC_H2OOCN',  units='-'      , &
              avgflag='A', long_name='coastal inundation fraction', &
              ptr_lnd=this%frac_h2oocn_grc)
      endif

    end subroutine InitHistory

    !------------------------------------------------------------------------
    subroutine Restart(this, bounds, ncid, flag)
      ! 
      ! !USES:
      use restUtilMod
      use ncdio_pio
      !
      ! !ARGUMENTS:
      class(ocn2lnd_type) :: this
      type(bounds_type), intent(in) :: bounds  
      type(file_desc_t), intent(inout) :: ncid   
      character(len=*) , intent(in)    :: flag   
      !
      ! !LOCAL VARIABLES:
      logical            :: readvar 
      !------------------------------------------------------------------------

      ! TODO: do we need restart for ocean-land one way coupling?

    end subroutine Restart

end module ocn2lndType