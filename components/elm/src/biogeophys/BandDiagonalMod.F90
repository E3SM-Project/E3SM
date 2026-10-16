module BandDiagonalMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Band Diagonal matrix solution
  !
  ! !USES:
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use decompMod      , only : bounds_type
  use abortutils     , only : endrun
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use elm_varctl     , only : iulog
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: BandDiagonal
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine BandDiagonal(bounds, lbj, ubj, jtop, jbot, numf, filter, b, r, u)
    !
    ! !DESCRIPTION:
    ! Tridiagonal matrix solution
    use lapack_acc_seq
    implicit none
    type(bounds_type), intent(in) :: bounds
    integer , intent(in)    :: lbj, ubj                        ! lbinning and ubing level indices
    integer , intent(in)    :: jtop(bounds%begc:)  
    integer , intent(in)    :: jbot( bounds%begc: )            ! bottom level for each column [col]
    integer , intent(in)    :: numf                            ! filter dimension
    integer , intent(in)    :: filter(:)                       ! filter
    real(r8), intent(in)    :: b( bounds%begc: , 1:   , lbj: ) ! compact band matrix [col, nband, j]
    real(r8), intent(in)    :: r( bounds%begc: , lbj: )        ! "r" rhs of linear system [col, j]
    real(r8), intent(inout) :: u( bounds%begc: , lbj: )        ! solution [col, j]
    !
    ! ! LOCAL VARIABLES:
    integer , parameter    :: nband = 5       ! band width
    integer  :: j,ci,fc,info,n,nmax, k,jstart,jstop        !indices
    integer, parameter  :: kl=(nband-1)/2,ku=kl              !number of sub/super diagonals
    integer, parameter  :: m=2*kl+ku+1
    real(r8),allocatable :: ab(:,:,:),temp(:,:,:) !compact storage array
    real(r8),allocatable :: result(:,:)

    !-----------------------------------------------------------------------

!Set up input matrix AB
!An m-by-n band matrix with kl subdiagonals and ku superdiagonals
!may be stored compactly in a two-dimensional array with
!kl+ku+1 rows and n columns
!AB(KL+KU+1+i-j,j) = A(i,j)

   ! calculate the maximum N-dimension needed (jbot-jtop +1)
   nmax = -100
   !$acc parallel loop gang vector default(present) reduction(max:nmax) copy(nmax)
   do fc = 1, numf
     ci = filter(fc) 
     n = jbot(ci)-jtop(ci)+1
     nmax = max(nmax,n)
   end do 

   ! m is the number of rows required for storage space by dgbsv
   ! n is the number of levels (snow/soil)
   !scs: replace ubj with jbot
   
   allocate(ab(1:m,1:nmax,1:numf), temp(1:m,1:nmax,1:numf), result(1:nmax,1:numf) ) 
   !$acc enter data create(ab(:,:,:) ,temp(:,:,:), result(:,:) )
    
    ! initialize the matrices
   !$acc parallel loop independent gang vector default(present) collapse(3)
   do fc = 1, numf  
      do j = 1, nmax 
         do k = 1, m
            ab(k,j,fc) = 0._r8 
            temp(k,j,fc) = 0._r8
         end do 
      end do 
   end do 

   !$acc parallel loop independent gang vector default(present)
   do fc = 1, numf
      ci = filter(fc)
      n = jbot(ci)-jtop(ci)+1
      result(1:n,fc) = r(ci, jtop(ci):jbot(ci))
   end do

   !$acc parallel loop independent gang vector default(present)  
   do fc = 1, numf
      ci = filter(fc) 
      n = jbot(ci)-jtop(ci)+1
      !band 1 :
      !$acc loop seq 
      do j = 0, n-3
         ! ab(kl+ku-1,3:n)=b(ci,1,jtop(ci):jbot(ci)-2)   ! 2nd superdiagonal
         ab(kl+ku-1,3+j,fc) = b(ci,1,jtop(ci)+j)   ! 2nd superdiagonal
      end do 
      !$acc loop seq 
      do j = 0, n-2
         !ab(kl+ku+0,2:n,fc) = b(ci,2,jtop(ci):jbot(ci)-1)   ! 1st superdiagonal
         ab(kl+ku+0,2+j,fc) = b(ci,2,jtop(ci)+j) 
      end do 
      !$acc loop seq 
      do j = 0, n-1 
         !ab(kl+ku+1,1:n,fc) = b(ci,3,jtop(ci):jbot(ci)  )   ! diagonal
         ab(kl+ku+1,1+j,fc) = b(ci,3,jtop(ci)+j)   ! diagonal
      end do 
      !$acc loop seq 
      do j = 0 , n-2
         !ab(kl+ku+2,1:n-1,fc) = b(ci,4,jtop(ci)+1:jbot(ci)) ! 1st subdiagonal
         ab(kl+ku+2,1+j,fc) = b(ci,4,jtop(ci)+1+j ) ! 1st subdiagonal
      end do
      !$acc loop seq 
      do j = 0, n-3 
         !ab(kl+ku+3,1:n-2,fc) = b(ci,5,jtop(ci)+2:jbot(ci)) ! 2nd subdiagonal 
         ab(kl+ku+3,1+j,fc) = b(ci,5,jtop(ci)+2+j) ! 2nd subdiagonal 
      end do 
   end do

   !$acc parallel loop independent gang vector default(present) collapse(3)
   do fc = 1, numf  
      do j = 1, nmax 
         do k = 1, m
            temp(k,j,fc) = ab(k,j,fc)
         end do 
      end do 
   end do 
 
   !$acc parallel loop independent gang vector default(present)  
   do fc = 1,numf
      ci = filter(fc)
      n = jbot(ci)-jtop(ci)+1
     call dgbsv_oacc(n, kl, ku, 1, ab(1:m,1:n,fc), m , result(1:n,fc),n,info)
     ! ! ! DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
     ! !  call dgbsv( n, kl, ku, 1, ab, m, ipiv, result, n, info )

    end do
    
    !$acc parallel loop independent gang vector default(present) 
    do fc = 1, numf
      ci = filter(fc) 
      n = jbot(ci)-jtop(ci)+1
      !$acc loop seq 
      do j =1, n
        u(ci, jtop(ci)+j-1) = result(j,fc)
      end do 
    end do 

    !$acc exit data delete(ab(:,:,:) ,temp(:,:,:), result(:,:) )  !!!!!!ipiv(:nmax,:numf)
    deallocate(temp)
    deallocate(ab)
    deallocate(result)
  end subroutine BandDiagonal

end module BandDiagonalMod
