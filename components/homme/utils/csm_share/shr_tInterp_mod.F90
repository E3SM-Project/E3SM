!===============================================================================
! SVN $Id: shr_tInterp_mod.F90 239 2006-02-08 19:02:33Z kauff $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070423/shr/shr_tInterp_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: shr_tInterp_mod -- time interpolation routines
!
! !DESCRIPTION:
!     shared time interpolation factor routines
!
! !REVISION HISTORY:
!     2004-Dec-10 - J. Schramm - first version
!     2005-Apr-10 - T. Craig - updated for shr bundles
!
! !INTERFACE: ------------------------------------------------------------------

module shr_tInterp_mod
 
! !USES:

   use shr_sys_mod   ! shared system calls
   use shr_cal_mod   ! shared calendar type and methods
   use shr_kind_mod  ! kinds for strong typing
   use shr_const_mod ! shared constants

   implicit none

   private ! except

! !PUBLIC TYPES:

   ! no public types

! !PUBLIC MEMBER FUNCTIONS:

   public :: shr_tInterp_getFactors  ! get time-interp factors
   public :: shr_tInterp_setAbort    ! set abort on error
   public :: shr_tInterp_setDebug    ! set debug level
   public :: shr_tInterp_getDebug    ! get debug level

! !PUBLIC DATA MEMBERS:

   ! no public data

!EOP

   real(SHR_KIND_R8),parameter :: c0 = 0.0_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: c1 = 1.0_SHR_KIND_R8
   real(SHR_KIND_R8),parameter :: eps = 1.0E-12_SHR_KIND_R8
   logical ,save               :: doabort = .true.
   integer ,save               :: debug = 0

!===============================================================================
contains
!===============================================================================

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_tInterp_getFactors -- calculate time interpolation factors
!
! !DESCRIPTION:
!     Returns two interpolation factors
!     Legal algorithms are (algo):
!       lower   - sets factors to data 1 (lower-bound), f1=1, f2=0
!       upper   - sets factors to data 2 (upper-bound), f1=0, f2=1
!       nearest - sets factors to nearest data in time
!       linear  - sets factors to linear interpolation between lb and ub
!     \newline
!     call shr\_tInterp\_getFactors(D1,s1,D2,s2,D,s,f1,f2,'linear',rc)
!     \newline
!     time of 2 >= time of 1  for all algos
!     time of 2 >= time of model data >= time of 1  for linear
!
! !REVISION HISTORY:
!     2005-Apr-10 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_tInterp_getFactors(D1,S1,D2,S2,Din,Sin,f1,f2,algo,rc)

   implicit none

! !USES:

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in)           :: D1,S1   ! LB date & sec (20010115,3600)
   integer(SHR_KIND_IN),intent(in)           :: D2,S2   ! UB date & sec
   integer(SHR_KIND_IN),intent(in)           :: Din,Sin ! desired/model date & sec
   real(SHR_KIND_R8)   ,intent(out)          :: f1      ! wgt for 1
   real(SHR_KIND_R8)   ,intent(out)          :: f2      ! wgt for 2
   character(*)        ,intent(in) ,optional :: algo    ! algorithm
   integer(SHR_KIND_IN),intent(out),optional :: rc      ! return code

!EOP

   !----- local  ------
   real(SHR_KIND_R8)      :: etime1      ! elapsed days for 1 (lower-bound)
   real(SHR_KIND_R8)      :: etime2      ! elapsed days for 2 (upper-bound)
   real(SHR_KIND_R8)      :: etimeIn     ! elapsed days for model date
   integer(SHR_KIND_IN)   :: eday        ! elapsed days
   character(SHR_KIND_CS) :: lalgo       ! local algo variable
   integer(SHR_KIND_IN)   :: lrc         ! local rc

   !----- formats -----
   character(*),parameter :: subName = "(shr_tInterp_getFactors)"
   character(*),parameter :: F00   = "('(shr_tInterp_getFactors) ',8a)" 
   character(*),parameter :: F01   = "('(shr_tInterp_getFactors) ',a,2f17.8)" 
   character(*),parameter :: F02   = "('(shr_tInterp_getFactors) ',a,3f17.8)" 
   character(*),parameter :: F03   = "('(shr_tInterp_getFactors) ',2a,3(i9.8,i6))" 

!-------------------------------------------------------------------------------
! Computes time interpolation factors
!-------------------------------------------------------------------------------

   lrc = 0

   if (present(algo)) then
     lalgo = algo
   else
     lalgo = 'linear'
   endif

   !--- compute elapsed time ---
   call shr_cal_date2eday(D1, eday)
   etime1 = real(eday,SHR_KIND_R8) + real(s1,SHR_KIND_R8)/shr_const_cDay

   call shr_cal_date2eday(D2, eday)
   etime2 = real(eday,SHR_KIND_R8) + real(s2,SHR_KIND_R8)/shr_const_cDay

   call shr_cal_date2eday(Din, eday)
   etimein = real(eday,SHR_KIND_R8) + real(sin,SHR_KIND_R8)/shr_const_cDay

   ! --- always check that 1 <= 2, although we could relax this requirement ---
   if (etime2 < etime1) then
     write(6,F01) ' ERROR: etime2 < etime1 ',etime2,etime1
     lrc = 1
     call shr_tInterp_abort(subName//' etime2 < etime1 ')
   endif

   f1 = c0
   f2 = c0
   ! --- set interpolation factors ---
   if (trim(lalgo) == 'lower') then
     if (etime1 < etime2) then
       f1 = c1
     else
       f2 = c1
     endif
   elseif (trim(lalgo) == 'upper') then
     if (etime1 < etime2) then
       f2 = c1
     else
       f1 = c1
     endif
   elseif (trim(lalgo) == 'nearest') then
     if (abs(etimein-etime1) <= abs(etime2-etimein)) then
       f1 = c1
       f2 = c0
     else
       f1 = c0
       f2 = c1
     endif
   elseif (trim(lalgo) == 'linear') then
     !--- check that etimein is between etime1 and etime2 ---
     if (etime2 < etimein .or. etime1 > etimein) then
       write(6,F02) ' ERROR illegal linear times: ',etime1,etimein,etime2
       lrc = 1
       call shr_tInterp_abort(subName//' illegal algo option '//trim(algo))
     endif
     if (etime2 == etime1) then
       f1 = 0.5_SHR_KIND_R8
       f2 = 0.5_SHR_KIND_R8
     else
       f1 = (etime2-etimein)/(etime2-etime1)
       f2 = c1 - f1
     endif
   else
     write(6,F00) 'ERROR: illegal lalgo option: ',trim(lalgo)
     lrc = 1
     call shr_tInterp_abort(subName//' illegal algo option '//trim(lalgo))
   endif

   !--- check that f1 and f2 are OK, each between 0 and 1 and they sum to 1 ---
   if (f1 < c0-eps .or. f1 > c1+eps .or. &
       f2 < c0-eps .or. f2 > c1+eps .or. &
       abs(f1+f2-c1) > eps) then
     write(6,F01) 'ERROR: illegal tInterp values ',f1,f2
     lrc = 1
     call shr_tInterp_abort(subName//' illegal tInterp values ')
   endif

   if (debug > 0) then
     write(6,F03) 'DEBUG: algo,D1,S1,Din,Sin,D2,S2=',trim(lAlgo),D1,S1,Din,Sin,D2,S2
     write(6,F01) 'DEBUG: algo,f1,f2= '//trim(lAlgo),f1,f2
   endif

   if (present(rc)) rc = lrc

end subroutine shr_tInterp_getFactors

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_tInterp_setAbort -- Set local shr_tInterp abort flag
!
! !DESCRIPTION:
!     Set local shr_tInterp abort flag, true = abort, false = print and continue
!     \newline
!     call shr\_tInterp\_setAbort(.false.)
!
! !REVISION HISTORY:
!     2005-Apr-10  - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_tInterp_setAbort(flag)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  logical,intent(in) :: flag

!EOP

  !--- formats ---
  character(*),parameter :: subName =   "(shr_tInterp_setAbort) "
  character(*),parameter :: F00     = "('(shr_tInterp_setAbort) ',a) "

!-------------------------------------------------------------------------------

  doabort = flag

end subroutine shr_tInterp_setAbort

!===============================================================================
!XXBOP ===========================================================================
!
! !IROUTINE: shr_tInterp_abort -- local interface for abort
!
! !DESCRIPTION:
!     Local interface for shr\_tInterp abort calls
!     \newline
!     call shr\_tInterp\_abort(subName//' ERROR illegal option')
!
! !REVISION HISTORY:
!     2005-Apr-10  - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_tInterp_abort(string)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  character(*),optional,intent(in) :: string

!XXEOP

  !--- formats ---
  character(SHR_KIND_CL) :: lstring
  character(*),parameter :: subName =   "(shr_tInterp_abort) "
  character(*),parameter :: F00     = "('(shr_tInterp_abort) ',a) "

!-------------------------------------------------------------------------------

  lstring = ''
  if (present(string)) lstring = string

  if (doabort) then
    call shr_sys_abort(lstring)
  else
    write(6,F00) trim(lstring)
  endif

end subroutine shr_tInterp_abort

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_tInterp_getDebug -- Get local shr_tInterp debug level
!
! !DESCRIPTION:
!     Get local shr_tInterp debug level, 0 = production
!     \newline
!     call shr\_tInterp\_getDebug(level)
!
! !REVISION HISTORY:
!     2005-Jun-14  - B. Kauffman
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_tInterp_getDebug(level)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  integer,intent(out) :: level

!EOP

  !--- formats ---
  character(*),parameter :: subName =   "(shr_tInterp_getDebug) "
  character(*),parameter :: F00     = "('(shr_tInterp_getDebug) ',a) "

!-------------------------------------------------------------------------------

  level = debug

end subroutine shr_tInterp_getDebug

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_tInterp_setDebug -- Set local shr_tInterp debug level
!
! !DESCRIPTION:
!     Set local shr_tInterp debug level, 0 = production
!     \newline
!     call shr\_tInterp\_setDebug(2)
!
! !REVISION HISTORY:
!     2005-Apr-10  - T. Craig - first prototype
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_tInterp_setDebug(iflag)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  integer,intent(in) :: iflag

!EOP

  !--- formats ---
  character(*),parameter :: subName =   "(shr_tInterp_setDebug) "
  character(*),parameter :: F01     = "('(shr_tInterp_setDebug) ',a,i3) "

!-------------------------------------------------------------------------------

  debug = iflag
  if (debug>0) write(6,F01) "DEBUG: level changed to ",debug

end subroutine shr_tInterp_setDebug

!===============================================================================
!===============================================================================

end module shr_tInterp_mod
