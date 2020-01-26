module shr_tInterp_mod

  ! shared time interpolation factor routines

  use ESMF
  use shr_sys_mod   ! shared system calls
  use shr_cal_mod   ! shared calendar type and methods
  use shr_kind_mod  ! kinds for strong typing
  use shr_const_mod ! shared constants
  use shr_kind_mod , only : i8=>shr_kind_i8, r8=>shr_kind_r8, cs=>shr_kind_cs, cl=>shr_kind_cl, cxx=>shr_kind_cxx
  use shr_log_mod  , only : s_loglev  => shr_log_Level
  use shr_log_mod  , only : s_logunit => shr_log_Unit
  use shr_orb_mod  , only : shr_orb_cosz, shr_orb_decl

  implicit none
  private ! except

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: shr_tInterp_getFactors  ! get time-interp factors
  public :: shr_tInterp_getAvgCosz  ! get cosz, time avg of
  public :: shr_tInterp_getCosz     ! get cosz
  public :: shr_tInterp_setAbort    ! set abort on error
  public :: shr_tInterp_setDebug    ! set debug level
  public :: shr_tInterp_getDebug    ! get debug level

  real(R8),parameter :: c0 = 0.0_R8
  real(R8),parameter :: c1 = 1.0_R8
  real(R8),parameter :: eps = 1.0E-12_R8

  logical :: doabort = .true.
  integer :: debug = 0

!===============================================================================
contains
!===============================================================================

  subroutine shr_tInterp_getFactors(D1,S1,D2,S2,Din,Sin,f1,f2,calendar,algo,rc)

    ! calculate time interpolation factors
    ! Returns two interpolation factors
    ! Legal algorithms are (algo):
    !   lower   - sets factors to data 1 (lower-bound), f1=1, f2=0
    !   upper   - sets factors to data 2 (upper-bound), f1=0, f2=1
    !   nearest - sets factors to nearest data in time
    !   linear  - sets factors to linear interpolation between lb and ub
    !   time of 2 >= time of 1  for all algos
    !   time of 2 >= time of model data >= time of 1  for linear

    ! !input/output parameters:
    integer      ,intent(in)           :: D1,S1   ! LB date & sec (20010115,3600)
    integer      ,intent(in)           :: D2,S2   ! UB date & sec
    integer      ,intent(in)           :: Din,Sin ! desired/model date & sec
    real(R8)     ,intent(out)          :: f1      ! wgt for 1
    real(R8)     ,intent(out)          :: f2      ! wgt for 2
    character(*) ,intent(in)           :: calendar!calendar type
    character(*) ,intent(in) ,optional :: algo    ! algorithm
    integer      ,intent(out),optional :: rc      ! return code

    !----- local  ------
    type(ESMF_Time)        :: itime1      ! time for 1 (lower-bound)
    type(ESMF_Time)        :: itime2      ! time for 2 (upper-bound)
    type(ESMF_Time)        :: itimeIn     ! time for model date
    type(ESMF_TimeInterval):: timeint     ! timeinterval
    integer(I8)            :: snum, sden  ! delta times in seconds
    integer(I8)            :: sint1,sint2 ! delta times in seconds
    character(CS)          :: lalgo       ! local algo variable
    integer                :: lrc         ! local rc
    character(*),parameter :: subName = "(shr_tInterp_getFactors) "
    character(*),parameter :: F00   = "('(shr_tInterp_getFactors) ',8a)"
    character(*),parameter :: F01   = "('(shr_tInterp_getFactors) ',a,2f17.8)"
    character(*),parameter :: F02   = "('(shr_tInterp_getFactors) ',a,3i9)"
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

    call shr_cal_timeSet(itimein,Din,Sin,calendar)
    call shr_cal_timeSet(itime1 ,D1 ,S1 ,calendar)
    call shr_cal_timeSet(itime2 ,D2 ,S2 ,calendar)

    !DML
    !  write(s_logunit,*) subName,' DTIME ',Din,D1,D2
    !  write(s_logunit,*) subName,' STIME ',sin,s1,s2
    !  write(s_logunit,*) subName,' ETIME ',etime1,etime2,etimein
    !  write(s_logunit,*) subName,' ITIME ',itime1,itime2,itimein

    ! --- always check that 1 <= 2, although we could relax this requirement ---
    if (itime2 < itime1) then
       if (s_loglev > 0) write(s_logunit,F01) ' ERROR: itime2 < itime1 D=',D1,S1,D2,S2
       lrc = 1
       call shr_tInterp_abort(subName//' itime2 < itime1 ')
    endif

    f1 = -1.0
    ! --- set interpolation factors ---
    if (trim(lalgo) == 'lower') then
       if (itime1 < itime2) then
          f1 = c1
       else
          f1 = c0
       endif
    elseif (trim(lalgo) == 'upper') then
       if (itime1 < itime2) then
          f1 = c0
       else
          f1 = c1
       endif
    elseif (trim(lalgo) == 'nearest') then
       timeint = itime1-itimein
       call ESMF_TimeIntervalGet(timeint,StartTimeIn=itimein,s_i8=sint1)
       timeint = itime2-itimein
       call ESMF_TimeIntervalGet(timeint,StartTimeIn=itimein,s_i8=sint2)
       if (abs(sint1) <= abs(sint2)) then
          f1 = c1
       else
          f1 = c0
       endif
    elseif (trim(lalgo) == 'linear') then
       !--- check that itimein is between itime1 and itime2 ---
       if (itime2 < itimein .or. itime1 > itimein) then
          write(s_logunit,F02) ' ERROR illegal linear times: ',D1,S1,Din,Sin,D2,S2
          lrc = 1
          call shr_tInterp_abort(subName//' illegal itimes ')
       endif
       if (itime2 == itime1) then
          f1 = 0.5_R8
       else
          timeint = itime2 - itimein
          call ESMF_TimeIntervalGet(timeint,StartTimeIn=itimein,s_i8=snum)
          timeint = itime2 - itime1
          call ESMF_TimeIntervalGet(timeint,StartTimeIn=itime1,s_i8=sden)
          f1 = real(snum,R8)/real(sden,R8)
       endif
    else
       if (s_loglev > 0) write(s_logunit,F00) 'ERROR: illegal lalgo option: ',trim(lalgo)
       lrc = 1
       call shr_tInterp_abort(subName//' illegal algo option '//trim(lalgo))
    endif

    f2 = c1 - f1

    !--- check that f1 and f2 are OK, each between 0 and 1 and they sum to 1 ---
    if (f1 < c0-eps .or. f1 > c1+eps .or. &
         f2 < c0-eps .or. f2 > c1+eps .or. &
         abs(f1+f2-c1) > eps) then
       if (s_loglev > 0) write(s_logunit,F01) 'ERROR: illegal tInterp values ',f1,f2
       lrc = 1
       call shr_tInterp_abort(subName//' illegal tInterp values ')
    endif

    if (debug > 0 .and. s_loglev > 0) then
       write(s_logunit,F03) 'DEBUG: algo,D1,S1,Din,Sin,D2,S2=',trim(lAlgo),D1,S1,Din,Sin,D2,S2
       write(s_logunit,F01) 'DEBUG: algo,f1,f2= '//trim(lAlgo),f1,f2
    endif

    if (present(rc)) rc = lrc

  end subroutine shr_tInterp_getFactors

  !===============================================================================
  subroutine shr_tInterp_getAvgCosz(tavCosz, lonr, latr,&
       ymd1, tod1, ymd2, tod2, eccen, mvelpp, lambm0, obliqr, dt, calendar)

    !   Returns a time-average of cos(z) over the time interval [LB,UB].
    !   The time-avg is calculated via a right-sum Riemann sum where the partitian
    !   width is the model dt, and the left-most partitian starts at the LB.
    !
    ! NOTE: For cosine of solar zenith angle forcing the time-stamps MUST be for
    !       the beginning of the interval.

    ! !input/output parameters:
    real(R8)     ,intent(out)  :: tavCosz(:) ! t-avg of cosz over [LB,UB]
    real(R8)     ,intent(in)   :: latr(:)    ! latitude
    real(R8)     ,intent(in)   :: lonr(:)    ! longitude
    integer      ,intent(in)   :: ymd1,tod1  ! date of lb
    integer      ,intent(in)   :: ymd2,tod2  ! date of ub
    real(R8)     ,intent(in)   :: eccen      ! orb param
    real(R8)     ,intent(in)   :: mvelpp     ! orb param
    real(R8)     ,intent(in)   :: lambm0     ! orb param
    real(R8)     ,intent(in)   :: obliqr     ! orb param
    integer      ,intent(in)   :: dt         ! model time step in secs
    character(*) ,intent(in)   :: calendar   ! calendar type

    !----- local  ------
    integer                 :: lsize             ! size of local data
    type(ESMF_Time)         :: reday1, reday2    ! LB, UB time
    type(ESMF_TimeInterval) :: timeint           ! time interval
    integer(I8)             :: n                 ! number of time-samples in t-avg
    type(ESMF_Time)         :: reday,reday0      ! time sample's elapsed seconds
    integer                 :: ymd,tod,ymd0,tod0 ! used to compute time of time-sample
    integer                 :: ldt               ! local dt as needed
    integer(I8)             :: ldt8              ! local dt as needed in i8
    integer(I8)             :: dtsec             ! delta time from timeint
    real(R8),pointer        :: cosz(:)           ! cos(zenith angle)
    character(*),parameter :: subName = "(shr_tInterp_getAvgCosz) "
    character(*),parameter :: F00   = "('(shr_tInterp_getAvgCosz) ',8a)"

    !-------------------------------------------------------------------------------
    ! Computes time avg cosz over interval [LB,UB]
    !-------------------------------------------------------------------------------

    lsize = size(lonr)
    allocate(cosz(lsize))
    if (lsize < 1 .or. size(latr) /= lsize .or. size(tavCosz) /= lsize) then
       call shr_sys_abort(subname//' ERROR: lon lat tavCosz sizes disagree')
    endif

    ldt = dt

    !--- get LB & UB dates ---
    call shr_cal_timeSet(reday1,ymd1,tod1,calendar)
    call shr_cal_timeSet(reday2,ymd2,tod2,calendar)
    if (reday1 > reday2) call shr_sys_abort(subname//'ERROR: lower-bound > upper-bound')

    timeint = reday2-reday1
    call ESMF_TimeIntervalGet(timeint,s_i8=dtsec)
    ldt8 = ldt
    if (mod(dtsec,ldt8) /= 0) then
       ldt8 = (dtsec)/((dtsec)/ldt8+1)
       ldt = ldt8
    endif

    !--- compute time average ---
    tavCosz = 0.0_R8 ! initialize partial sum
    n       = 0               ! dt weighted average in t-avg
    reday   = reday1          ! mid [LB,UB] interval t-step starts at LB
    ymd     = ymd1
    tod     = tod1
    do while( reday < reday2) ! mid-interval t-steps thru interval [LB,UB]

       !--- advance to next time in [LB,UB] ---
       ymd0 = ymd
       tod0 = tod
       reday0 = reday
       call shr_cal_advDateInt(ldt,'seconds',ymd0,tod0,ymd,tod,calendar)
       call shr_cal_timeSet(reday,ymd,tod,calendar)

       if (reday > reday2) then
          ymd = ymd2
          tod = tod2
          timeint = reday2-reday0
          call ESMF_TimeIntervalGet(timeint,s_i8=dtsec)
          ldt = dtsec
       endif

       !--- get next cosz value for t-avg ---
       call shr_tInterp_getCosz(cosz,lonr,latr,ymd,tod,eccen,mvelpp,lambm0,obliqr,calendar)
       n = n + ldt
       tavCosz = tavCosz + cosz*real(ldt,R8)  ! add to partial sum

    end do
    tavCosz = tavCosz/real(n,R8) ! form t-avg

    deallocate( cosz )

  end subroutine shr_tInterp_getAvgCosz

  !===============================================================================
  subroutine shr_tInterp_getCosz(cosz,lonr,latr,ymd,tod,eccen,mvelpp,lambm0,obliqr,calendar)

    ! Calculate the cos(solar-zenith angle).

    ! !input/output parameters:
    real(R8)     , intent(out)   :: cosz(:)    ! cos(zenith angle)
    real(R8)     , intent(in)    :: latr(:)    ! latitude
    real(R8)     , intent(in)    :: lonr(:)    ! longitude
    integer      , intent(in)    :: ymd,tod    ! date of interest
    real(R8)     , intent(in)    :: eccen      ! orb param
    real(R8)     , intent(in)    :: mvelpp     ! orb param
    real(R8)     , intent(in)    :: lambm0     ! orb param
    real(R8)     , intent(in)    :: obliqr     ! orb param
    character(*) , intent(in)    :: calendar   ! calendar type

    !----- local  ------
    integer                 :: n
    integer                 :: lsize
    real(R8)                :: calday      ! julian days
    real(R8)                :: declin,eccf ! orb params
    real(R8)     ,parameter :: solZenMin = 0.001_R8 ! min solar zenith angle
    character(*) ,parameter :: subName = "(shr_tInterp_getCosz) "
    character(*) ,parameter :: F00   = "('(shr_tInterp_getCosz) ',8a)"

    !-------------------------------------------------------------------------------
    ! Returns cos(zenith angle)
    !-------------------------------------------------------------------------------

    lsize = size(lonr)
    if (lsize < 1 .or. size(latr) /= lsize .or. size(cosz) /= lsize) then
       call shr_sys_abort(subname//' ERROR: lon lat cosz sizes disagree')
    endif

    call shr_cal_date2julian(ymd,tod,calday,calendar)
    call shr_orb_decl( calday, eccen, mvelpp, lambm0, obliqr, declin, eccf)
    do n = 1,lsize
       cosz(n) = max(solZenMin, shr_orb_cosz( calday, latr(n), lonr(n), declin))
    end do

  end subroutine shr_tInterp_getCosz

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
    !
    !-------------------------------------------------------------------------------

    doabort = flag

  end subroutine shr_tInterp_setAbort

  !===============================================================================
  subroutine shr_tInterp_abort(string)

    ! local interface for abort

    ! !input/output parameters:
    character(*),optional,intent(in) :: string

    !--- formats ---
    character(CL) :: lstring
    character(*),parameter :: subName =   "(shr_tInterp_abort) "
    character(*),parameter :: F00     = "('(shr_tInterp_abort) ',a) "
    !-------------------------------------------------------------------------------

    lstring = ''
    if (present(string)) lstring = string
    if (doabort) then
       call shr_sys_abort(lstring)
    else
       write(s_logunit,F00) trim(lstring)
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
    !
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
    !
    !-------------------------------------------------------------------------------

    debug = iflag
    if (debug>0 .and. s_loglev > 0) write(s_logunit,F01) "DEBUG: level changed to ",debug

  end subroutine shr_tInterp_setDebug

  !===============================================================================
  !===============================================================================

end module shr_tInterp_mod
