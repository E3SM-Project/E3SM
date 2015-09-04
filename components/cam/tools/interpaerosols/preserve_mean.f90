module preserve_mean

   use shr_kind_mod, only : r8=>shr_kind_r8
!   use globals, only: nxo, nyo, ntime

   implicit none
   private
   public ::  monthly_to_midmonth
contains

subroutine monthly_to_midmonth (arr, nlev, nlonout,nxo,nyo,ntime)
!-------------------------------------------------------------------
!
!  Input/Output:
!     AEROSOL, M_ps_cam from monthly averages to mid-month values
!
!--------------------------------------------------------------------
!
! Arguments
!
   integer, intent(in) :: nlev
   integer, intent(in) :: nxo,nyo,ntime
   integer, intent(in) :: nlonout(nyo)
   real(r8), intent(inout) :: arr(nxo,nyo,nlev,ntime)  ! input/output array
!
! Local workspace
!
   real(r8) Matrix(1:12,1:12)   ! map from mid-month to monthly average
   integer pivot(12)            ! pivot from LU factor
   real(r8) vec(12)             ! vector to be used for LU solve
   integer i,j,k                ! spatial indices
   integer :: m                 ! constituent index
!
! ntime must be 12 
!
   if (ntime /= 12) then
      write(6,*) 'monthly_to_midmonth: ntime must be 12 got ', ntime
      stop 999
   end if
!
! construct and factor linear map from 
! mid month values to monthly averages
!
   call construct_matrix (Matrix, pivot)
!
! convert monthly averages to mid-month values
! could be sped up by solving many rhs at once
!  rather than one at a time.
!
   do k=1,nlev
      do j=1,nyo
         do i=1,nlonout(j)
            vec(:) = arr(i,j,k,:)
            call LUSolve (vec, Matrix, pivot)
            arr(i,j,k,:) = vec(:)
         end do
      end do
   end do

   return
end subroutine monthly_to_midmonth

subroutine construct_matrix (Matrix, pivot)
!--------------------------------------------------------------------
!  Input: None
!    
!  Output:
!    Matrix (LU Factored form)
!    pivot (pivots corresponding to LU Factored form of matrix)
!
!  Method:
!    construct matrix representing averaging 
!      from mid-monthly values
!      to monthly average values
!    factor matrix (matrix => LUFactoredMatrix, pivots) so that the 
!      backsolves represent map from averages to mid-month values
!--------------------------------------------------------------------

  real(r8), intent(out) :: Matrix(12,12) !map from mid-month values to monthly average
  integer, intent(out) :: pivot(12) ! pivots from LU factorization of Matrix

  integer row, col, colp, colm  ! indexes ino matrix
  real(r8) :: Len(12)   ! number of days in month
  data Len /31.0, 28.0, 31.0, 30.0, 31.0, 30.0, 31.0, 31.0, 30.0, 31.0, 30.0, 31.0 /

  Matrix(:,:) = 0.0  
! 
! construct tridiagonal cyclic matrix
!
  do row = 1,12
    col  = row
    colm = row - 1
    colp = row + 1
    if (colm ==  0) then 
      colm = 12
    endif
    if (colp == 13) then 
      colp =  1
    endif
    Matrix(row,colm) = Len(row) / (4*(Len(row )+Len(colm)))
    Matrix(row,colp) = Len(row) / (4*(Len(colp)+Len(row )))
    Matrix(row,col ) = 1.0 - Matrix(row,colm) - Matrix(row,colp)
  enddo

  call LUFactor (Matrix, pivot);

  return 
end subroutine construct_matrix


subroutine LUFactor (Matrix, pivot)
!-----------------------------------------------------------------
! interface to matrix factorization routine
!  here, LAPACK
!-----------------------------------------------------------------

  real(r8) Matrix(12,12) ! matrix to be factored
  integer pivot(12)      ! pivots from factorization
  integer :: info = 0    ! did factorization fail?

  call DGETRF (12, 12, Matrix, 12, pivot, info)
  if (info < 0) then
    write(6,*)"AEROSOL_INITIALIZE:DGETRF (LAPACK) factor argument ",-info," has illegal value"
    stop 999
  elseif(info > 0) then
    write(6,*)"AEROSOL_INITIALIZE:DGETRF (LAPACK) factor matrix upper element ",info," was zero"
    stop 999
  endif

  return
end subroutine LUFactor


subroutine LUSolve (vec, Matrix, pivot)
!-----------------------------------------------------------------
! interface to matrix back-solve routine
!  here, LAPACK
!-----------------------------------------------------------------

  real(r8)  vec(12)        ! vector to be solved.  vector then solution
  real(r8)  Matrix(12,12)  ! factored matrix
  integer  pivot(12)       ! pivots from factorization
  integer :: info = 0      ! did solve fail?

  call DGETRS ('N', 12, 1, Matrix, 12, pivot, vec, 12, info)
  if (info < 0) then
    write(6,*)"AEROSOL_INITIALIZE:DGETRS (LAPACK) solve argument ",-info," has illegal value"
    stop 999
  endif

  return
end subroutine LUSolve

end module preserve_mean
