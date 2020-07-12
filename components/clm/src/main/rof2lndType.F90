
!50==================================================
 ! Author: Chang Liao( changliao at pnnl.gov )
 ! Module: H2SC (hillslope based soil column drainage function)
 ! rof->lnd exchange
 ! First edit: 20180530
 ! Revison: 20200626 Chang Liao
 !50==================================================
module rof2lndType

    !-----------------------------------------------------------------------
    
    !
    ! !USES:
    use shr_kind_mod  , only : r8 => shr_kind_r8
    use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
    use shr_log_mod   , only : errMsg => shr_log_errMsg  
    use clm_varctl    , only : iulog  
    use clm_varcon    , only : spval
    use decompMod     , only : bounds_type
    use abortutils    , only : endrun

    !
    ! !PUBLIC TYPES:
    implicit none
    private
    save
    !
  ! !PUBLIC DATA TYPES:
  !----------------------------------------------------
  ! mosart -> land variables structure
  !
 
  !----------------------------------------------------
  type, public :: rof2lnd_type
  !DMR additions for CPL_BYPASS option
#ifdef CPL_BYPASS
  
#endif
real(r8), pointer :: channel_depth  (:)   => null() 
real(r8), pointer :: gage_height  (:)   => null() !rof stream gage height
real(r8), pointer :: hillslope_slope  (:)   => null() 
real(r8), pointer :: hillslope_length  (:)   => null()
real(r8), pointer :: elevation_profile1  (:)   => null() 
real(r8), pointer :: elevation_profile2  (:)   => null() 
real(r8), pointer :: elevation_profile3  (:)   => null() 
real(r8), pointer :: elevation_profile4  (:)   => null() 
real(r8), pointer :: elevation_profile5  (:)   => null() 
real(r8), pointer :: elevation_profile6  (:)   => null() 
real(r8), pointer :: elevation_profile7  (:)   => null() 
real(r8), pointer :: elevation_profile8  (:)   => null() 
real(r8), pointer :: elevation_profile9  (:)   => null() 
real(r8), pointer :: elevation_profile10  (:)   => null()
real(r8), pointer :: elevation_profile11  (:)   => null()
contains

 procedure, public  :: Init
 procedure, private :: InitAllocate 
 procedure, private :: InitHistory  
 procedure, public  :: InitAccBuffer
 procedure, public  :: InitAccVars
 procedure, public  :: UpdateAccVars
 procedure, public  :: Restart

end type rof2lnd_type


!----------------------------------------------------

contains

!------------------------------------------------------------------------
subroutine Init(this, bounds)

class(rof2lnd_type) :: this
type(bounds_type), intent(in) :: bounds  

call this%InitAllocate(bounds)
call this%InitHistory(bounds)

end subroutine Init

!------------------------------------------------------------------------
subroutine InitAllocate(this, bounds)
!
! !DESCRIPTION:
! Initialize rof2lnd derived type
!
! !ARGUMENTS:
class(rof2lnd_type) :: this
type(bounds_type), intent(in) :: bounds  
!
! !LOCAL VARIABLES:
real(r8) :: ival  = 0.0_r8  ! initial value
integer  :: begg, endg
integer  :: begc, endc
integer  :: begp, endp
!------------------------------------------------------------------------

begg = bounds%begg; endg= bounds%endg
begc = bounds%begc; endc= bounds%endc
begp = bounds%begp; endp= bounds%endp
#ifdef CPL_BYPASS
#endif

allocate(this%channel_depth   (begg:endg) ); this%channel_depth  (:)   = ival
allocate(this%gage_height   (begg:endg) ); this%gage_height  (:)   = ival
allocate(this%hillslope_slope   (begg:endg) ); this%hillslope_slope  (:)   = ival
allocate(this%hillslope_length   (begg:endg) ); this%hillslope_length (:)   = ival
allocate(this%elevation_profile1   (begg:endg) ); this%elevation_profile1 (:)   = ival
allocate(this%elevation_profile2   (begg:endg) ); this%elevation_profile2 (:)   = ival
allocate(this%elevation_profile3   (begg:endg) ); this%elevation_profile3 (:)   = ival
allocate(this%elevation_profile4   (begg:endg) ); this%elevation_profile4 (:)   = ival
allocate(this%elevation_profile5   (begg:endg) ); this%elevation_profile5 (:)   = ival
allocate(this%elevation_profile6   (begg:endg) ); this%elevation_profile6 (:)   = ival
allocate(this%elevation_profile7   (begg:endg) ); this%elevation_profile7 (:)   = ival
allocate(this%elevation_profile8   (begg:endg) ); this%elevation_profile8 (:)   = ival
allocate(this%elevation_profile9   (begg:endg) ); this%elevation_profile9 (:)   = ival
allocate(this%elevation_profile10   (begg:endg) ); this%elevation_profile10 (:)   = ival
allocate(this%elevation_profile11   (begg:endg) ); this%elevation_profile11 (:)   = ival

end subroutine InitAllocate

!------------------------------------------------------------------------
subroutine InitHistory(this, bounds)
!
! !USES:
use histFileMod, only : hist_addfld1d
!
! !ARGUMENTS:
class(rof2lnd_type) :: this
type(bounds_type), intent(in) :: bounds  
!
! !LOCAL VARIABLES:
integer  :: begg, endg
integer  :: begp, endp
!---------------------------------------------------------------------

begg = bounds%begg; endg= bounds%endg
begp = bounds%begp; endp= bounds%endp
this%gage_height(begg:endg) = spval
call hist_addfld1d (fname='GAGEH',  units='meter',  &
avgflag='A', long_name='gage height', &
ptr_lnd=this%gage_height)
#ifdef CPL_BYPASS
#endif
end subroutine InitHistory

!-----------------------------------------------------------------------
subroutine InitAccBuffer (this, bounds)

use clm_varcon  , only : spval
use accumulMod  , only : init_accum_field
!
! !ARGUMENTS:
class(rof2lnd_type) :: this
type(bounds_type), intent(in) :: bounds  
!---------------------------------------------------------------------


end subroutine InitAccBuffer

!-----------------------------------------------------------------------
subroutine InitAccVars(this, bounds)
!
! !DESCRIPTION:
! Initialize module variables that are associated with
! time accumulated fields. This routine is called for both an initial run
! and a restart run (and must therefore must be called after the restart file 
! is read in and the accumulation buffer is obtained)
!
! !USES 
use accumulMod       , only : extract_accum_field
use clm_time_manager , only : get_nstep
!
! !ARGUMENTS:
class(rof2lnd_type) :: this
type(bounds_type), intent(in) :: bounds  
!
! !LOCAL VARIABLES:
integer  :: begp, endp
integer  :: nstep
integer  :: ier
real(r8), pointer :: rbufslp(:)  ! temporary
!---------------------------------------------------------------------

begp = bounds%begp; endp = bounds%endp

! Determine time step
nstep = get_nstep()



end subroutine InitAccVars

!-----------------------------------------------------------------------
subroutine UpdateAccVars (this, bounds)
!
! USES
use clm_time_manager, only : get_nstep
use accumulMod      , only : update_accum_field, extract_accum_field
!
! !ARGUMENTS:
class(rof2lnd_type)                 :: this
type(bounds_type)      , intent(in) :: bounds  
!
! !LOCAL VARIABLES:
integer :: g,c,p                     ! indices
integer :: dtime                     ! timestep size [seconds]
integer :: nstep                     ! timestep number
integer :: ier                       ! error status
integer :: begp, endp
real(r8), pointer :: rbufslp(:)      ! temporary single level - pft level
!---------------------------------------------------------------------

begp = bounds%begp; endp = bounds%endp

nstep = get_nstep()

! Allocate needed dynamic memory for single level pft field


end subroutine UpdateAccVars

!------------------------------------------------------------------------
subroutine Restart(this, bounds, ncid, flag)
! 
! !USES:
use restUtilMod
use ncdio_pio
!
! !ARGUMENTS:
class(rof2lnd_type) :: this
type(bounds_type), intent(in) :: bounds  
type(file_desc_t), intent(inout) :: ncid   
character(len=*) , intent(in)    :: flag   
!
! !LOCAL VARIABLES:
logical            :: readvar 
!------------------------------------------------------------------------
call restartvar(ncid=ncid, flag=flag, varname='gage_height', xtype=ncd_double, &
         dim1name='gridcell', &
         long_name='gage height', units='meter', &
         interpinic_flag='skip', readvar=readvar, data=this%gage_height)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, readvar=readvar, not restart: initialize flood to zero
       this%gage_height = 0._r8
    endif


end subroutine Restart

end module rof2lndType