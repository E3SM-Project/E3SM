!-------------------------------------------------------------------------------
! $Id: clubb_precision.F90 6849 2014-04-22 21:52:30Z charlass@uwm.edu $
!===============================================================================
module clubb_precision

  implicit none

  public :: stat_nknd, stat_rknd, time_precision, dp, core_rknd

  private ! Default scope

  ! This definition of double precision must use a real type that is 64 bits
  ! wide, because (at least) the LAPACK routines depend on this definition being
  ! accurate. Otherwise, LAPACK must be recompiled, or some other trickery must
  ! be done.
  integer, parameter :: &
    dp = selected_real_kind( p=12 )    ! double precision

  ! The precisions below are arbitrary, and could be adjusted as
  ! needed for long simulations or time averaging.  Note that on
  ! most machines 12 digits of precision will use a data type
  ! which is 8 bytes long.
  integer, parameter ::  & 
    stat_nknd = selected_int_kind( 8 ), & 
    stat_rknd = selected_real_kind( p=12 ), & 
    time_precision = selected_real_kind( p=12 ), &
    core_rknd = CLUBB_REAL_TYPE ! Value from the preprocessor directive

end module clubb_precision
!-------------------------------------------------------------------------------
