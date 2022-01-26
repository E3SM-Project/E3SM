module SimpleMathMod

#include "shr_assert.h"
 
use shr_kind_mod , only : r8 => shr_kind_r8
!------------------------------------------------------------------------------
  !
  ! DESCRIPTIONS:
  ! module contains simple mathematical functions for arrays
  ! Created by Jinyun Tang, Feb., 2014

implicit none

  interface array_normalization
    module procedure array_normalization_2d, array_normalization_2d_filter
  end interface array_normalization

  interface array_div_vector
    module procedure array_div_vector_filter, array_div_vector_nofilter
  end interface array_div_vector

  public :: shr_flux_update_stress_elm
contains
!--------------------------------------------------------------------------------
  subroutine array_normalization_2d(which_dim, arr2d_inout)
  !$acc routine seq
  !DESCRIPTIONS
  !do normalization for the input array along dimension which_dim
  !
  !USES
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none

  integer,  intent(in) :: which_dim     !do normalization along which dimension?
  real(r8), intent(inout) :: arr2d_inout(:,:)   !input 2d array


  !local variables
  integer  :: sz1, sz2     !array size
  integer  :: j1, j2       !indices
  real(r8) :: arr_sum

  sz1 = size(arr2d_inout,1)
  sz2 = size(arr2d_inout,2)

  if(which_dim==1)then
    !normalize along dimension 1, so loop along dimension 2
    do j2 = 1, sz2
      !obtain the total
      arr_sum=0._r8
      do j1 = 1, sz1
        arr_sum=arr_sum+arr2d_inout(j1,j2)
      enddo
      !normalize with the total if arr_sum is non-zero
      if(arr_sum/=0._r8)then
        do j1 = 1, sz1
          arr2d_inout(j1,j2) = arr2d_inout(j1,j2)/arr_sum
        enddo
      endif
    enddo
  elseif(which_dim==2)then
    !normalize along dimension 2, so loop along dimension 1
    do j1 = 1, sz1
      !obtain the total
      arr_sum=0._r8
      do j2 = 1, sz2
        arr_sum=arr_sum+arr2d_inout(j1,j2)
      enddo
      !normalize with the total if arr_sum is non-zero
      !I think there should be a safer mask for this to screen off spval values
      !Jinyun Tang, May 30, 2014
      if(arr_sum>0._r8 .or. arr_sum < 0._r8)then
        do j2 = 1, sz2
          arr2d_inout(j1,j2) = arr2d_inout(j1,j2)/arr_sum
        enddo
      endif
    enddo
  endif
  return
  end subroutine array_normalization_2d

!--------------------------------------------------------------------------------
subroutine array_normalization_2d_filter(lbj1, ubj1, lbj2, ubj2, numf, filter, arr2d_inout)
  !$acc routine seq
  !DESCRIPTIONS
  !do normalization with filter for the input array along dimension 2

  !
  !USES
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  integer,  intent(in) :: lbj1         !left bound of dim 1
  integer,  intent(in) :: lbj2         !left bound of dim 2
  integer,  intent(in) :: ubj1         !right bound of dim 1
  integer,  intent(in) :: ubj2         !right bound of dim 2
  integer,  intent(in) :: numf         !filter size
  integer,  intent(in) :: filter(:)    !filter
  real(r8), intent(inout) :: arr2d_inout(lbj1: , lbj2: )   !input 2d array


  !local variables
  integer  :: sz1, sz2     !array size
  integer  :: j2           !indices
  integer  :: f, p         !indices
  real(r8) :: arr_sum(lbj1:ubj1)

  ! Enforce expected array sizes


  arr_sum(:) = 0._r8
  do j2 = lbj2, ubj2
    do f = 1, numf
      p = filter(f)
      !obtain the total
      arr_sum(p)=arr_sum(p)+arr2d_inout(p,j2)
    enddo
  enddo

    !normalize with the total if arr_sum is non-zero
  do j2 = lbj2, ubj2
    do f = 1, numf
      p = filter(f)
      !I found I have to ensure >0._r8 because of some unknown reason, jyt May 23, 2014
      !I will test this later with arr_sum(p)/=0._r8
      if(arr_sum(p)>0._r8 .or. arr_sum(p)<0._r8)then
        arr2d_inout(p,j2) = arr2d_inout(p,j2)/arr_sum(p)
      endif
    enddo
  enddo
  return
  end subroutine array_normalization_2d_filter
!--------------------------------------------------------------------------------

  subroutine array_div_vector_filter(lbj1, ubj1, lbj2, ubj2, &
       arr1d_in, fn, filter,  arr2d_inout)
  !$acc routine seq
  !DESCRIPTIONS
  !array divided by a vector, arr2d_in is divided by one
  !element in arr1d_in
  !It always assumes the filter is along with dimenion 1
  !
  ! USES
  !
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  integer,  intent(in) :: lbj1         !left bound of dim 1
  integer,  intent(in) :: lbj2         !left bound of dim 2
  integer,  intent(in) :: ubj1         !right bound of dim 1
  integer,  intent(in) :: ubj2         !right bound of dim 2
  real(r8), intent(in) :: arr1d_in(lbj1: )   !1d scaling factor
  integer , intent(in) :: fn
  integer , intent(in) :: filter(:)     !filter
  real(r8), intent(inout) :: arr2d_inout(lbj1: ,lbj2: ) !2d array to be scaled

  integer :: sz
  integer :: j, f, p

  ! Enforce expected array sizes


  do j = lbj2, ubj2
     do f = 1, fn
        p = filter(f)
        if (arr1d_in(p) > 0._r8 .or. arr1d_in(p) < 0._r8) then
           arr2d_inout(p,j) = arr2d_inout(p,j)/arr1d_in(p)
        else
           arr2d_inout(p,j) = 0._r8
        end if
     end do
  end do
  return
  end subroutine array_div_vector_filter

!--------------------------------------------------------------------------------

  subroutine array_div_vector_nofilter(arr1d_in, which_dim, arr2d_inout)
  !$acc routine seq
  !DESCRIPTIONS
  !array divided by a vector, each row in arr2d_in is divided by one
  !element in arr1d_in
  !
  !USES
  !
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  real(r8), intent(in) :: arr1d_in(:)     !scaling factor
  integer,  intent(in) :: which_dim        !which dimension is scaled
  real(r8), intent(inout) :: arr2d_inout(:,:)   !2d array to be scaled

  integer :: sz1, sz2
  integer :: j1, j2

  sz1=size(arr2d_inout,1)
  sz2=size(arr2d_inout,2)

  if(which_dim==1)then
    ! Enforce expected array sizes
    !#py !#py call shr_assert(sz1    == size(arr1d_in), errMsg(__FILE__, __LINE__))

    do j2 = 1, sz2
      do j1 = 1, sz1
        if(arr1d_in(j1)>0._r8)then
          arr2d_inout(j1,j2) = arr2d_inout(j1,j2)/arr1d_in(j1)
        endif
      enddo
    enddo
  else
    ! Enforce expected array sizes
    !#py !#py call shr_assert(sz2    == size(arr1d_in), errMsg(__FILE__, __LINE__))

    do j2 = 1, sz2
      do j1 = 1, sz1
        if(arr1d_in(j2)>0._r8 .or. arr1d_in(j2)<0._r8)then
          arr2d_inout(j1,j2) = arr2d_inout(j1,j2)/arr1d_in(j2)
        endif
      enddo
    enddo

  endif
  return
  end subroutine array_div_vector_nofilter

  PURE INTEGER FUNCTION MAXLOC_(ARR)
        !$ACC ROUTINE SEQ
        use shr_kind_mod, only: r8 => shr_kind_r8
        REAL(R8),DIMENSION(:),INTENT(IN) :: ARR
        REAL(R8),DIMENSION(SIZE(ARR)) :: LIS
        REAL(R8) :: VAL
        INTEGER :: I
        VAL=MAXVAL(ARR)
        DO I=1,SIZE(ARR),1
            IF(ARR(I)==VAL)THEN
                MAXLOC_=I
                RETURN
            END IF
        END DO
    END FUNCTION MAXLOC_

    subroutine matvec_acc(START,END_,RES,A,X)
      !$acc routine seq
      !As of Cuda 10.1 calling cuBlas functions from device code
      !is not supported.  So must create any blas routines with
      !acc routine seq manually.  This is for square matrices Matrix Vector Multiplication
      !used only in PhotosynthesisMod::calcstressroot so far

      use shr_kind_mod , only : r8 => shr_kind_r8
      INTEGER, INTENT(IN) :: START, END_ !section of matrices/vector to multiply
      REAL(R8) , INTENT(INOUT) :: RES(START:END_)
      REAL(R8) , INTENT(IN)  :: A(START:END_,START:END_), X(START:END_)
      REAL(R8)  :: transA(START:END_,START:END_)
      REAL(R8) :: SUM
      INTEGER :: COL, ROW

      transA = transpose(A)

      DO COL = START, END_

        RES(COL) = DOT_PRODUCT(transA(start:end_,col),X(start:end_))

      END DO

    end subroutine matvec_acc



!===============================================================================

! After the stress has been recalculated in iteration loops, call this routine
! to adjust the value used in the next iteration for improved convergence.
!
! We use the sign convention where all arguments are assumed to be positive,
! except for tau_diff and prev_tau_diff, which can be of either sign.
subroutine shr_flux_update_stress_elm(wind0, wsresp, tau_est, tau, prev_tau, &
     tau_diff, prev_tau_diff, wind_adj)
  !$acc routine seq 
  ! Wind speed from atmosphere (not updated by iteration) [m/s]
  real(r8), intent(in) :: wind0
  ! Response of boundary layer wind to stress changes in a time step [m/s/Pa]
  real(r8), intent(in) :: wsresp
  ! Estimated tau that would be in equilibrium with boundary layer wind [Pa]
  real(r8), intent(in) :: tau_est
  ! Stress that has just been calculated in this loop [Pa]
  real(r8), intent(inout) :: tau
  ! Stress assumed when calculating this tau (use tau_est for first iter.) [Pa]
  real(r8), intent(inout) :: prev_tau
  ! Difference between last two values of tau used for iterations [Pa]
  ! (Use a very large value for first iteration)
  real(r8), intent(inout) :: tau_diff
  ! Copy of the input tau_diff (for diagnostic purposes) [Pa]
  real(r8), intent(out) :: prev_tau_diff
  ! Wind updated by iteration [m/s]
  real(r8), intent(out) :: wind_adj

  real(r8) :: wsresp_applied

  ! maximum ratio between abs(tau_diff) and abs(prev_tau_diff)
  real(r8) :: tau_diff_fac

  ! Using wsresp_applied = 0.5 * wsresp improves accuracy somewhat, similar to
  ! using the trapezoidal method rather than backward Euler.
  wsresp_applied = 0.5 * wsresp

  ! Forbid removing more than 99% of wind speed.
  ! This is applied before anything else to ensure that the applied stress is
  ! not so strong that it reverses the direction of the winds.
  if ( (tau - tau_est) * wsresp_applied > 0.99_r8 * wind0 ) then
     tau = tau_est + 0.99_r8 * wind0 / wsresp_applied
  end if

  prev_tau_diff = tau_diff
  tau_diff = tau - prev_tau

  ! damp large changes each iteration for convergence
  if (tau_diff * prev_tau_diff < 0._r8) then
     ! if oscillating about the solution, use stronger damping
     ! this should not be set below 0.5
     tau_diff_fac = 0.6_r8
  else
     tau_diff_fac = 0.95_r8
  end if

  if (abs(tau_diff) > abs(tau_diff_fac*prev_tau_diff)) then
     tau_diff = sign(tau_diff_fac*prev_tau_diff, tau_diff)
     tau = prev_tau + tau_diff
  end if
  prev_tau = tau

  wind_adj = wind0 - (tau - tau_est) * wsresp_applied

end subroutine shr_flux_update_stress_elm

end module SimpleMathMod
