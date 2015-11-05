module TridiagonalMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Tridiagonal matrix solution
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: Tridiagonal
  public :: trisim
  interface Tridiagonal
    module procedure Tridiagonal_sr
    module procedure Tridiagonal_mr
  end interface Tridiagonal

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Tridiagonal_sr (bounds, lbj, ubj, jtop, numf, filter, a, b, c, r, u, is_col_active)
    !
    ! !DESCRIPTION:
    ! Tridiagonal matrix solution
    ! A x = r
    ! where x and r are vectors
    ! !USES:
    use shr_kind_mod   , only: r8 => shr_kind_r8
    use clm_varctl     , only : iulog
    use decompMod      , only : bounds_type
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type) , intent(in)    :: bounds                                   ! bounds
    integer           , intent(in)    :: lbj, ubj                                 ! lbinning and ubing level indices
    integer           , intent(in)    :: jtop( bounds%begc: bounds%endc)          ! top level for each column [col]
    integer           , intent(in)    :: numf                                     ! filter dimension
    integer           , intent(in)    :: filter(:)                                ! filter
    real(r8)          , intent(in)    :: a( bounds%begc:bounds%endc , lbj:ubj)    ! "a" left off diagonal of tridiagonal matrix [col , j]
    real(r8)          , intent(in)    :: b( bounds%begc:bounds%endc , lbj:ubj)    ! "b" diagonal column for tridiagonal matrix [col  , j]
    real(r8)          , intent(in)    :: c( bounds%begc:bounds%endc , lbj:ubj)    ! "c" right off diagonal tridiagonal matrix [col   , j]
    real(r8)          , intent(in)    :: r( bounds%begc:bounds%endc , lbj:ubj)    ! "r" forcing term of tridiagonal matrix [col      , j]
    real(r8)          , intent(inout) :: u( bounds%begc:bounds%endc , lbj:ubj)    ! solution [col                                    , j]
                                                                                  !
    integer                           :: j,ci,fc                                  ! indices
    logical, optional, intent(in)     :: is_col_active(bounds%begc:bounds%endc)   !
    logical                           :: l_is_col_active(bounds%begc:bounds%endc) !
    real(r8)                          :: gam(bounds%begc:bounds%endc,lbj:ubj)     ! temporary
    real(r8)                          :: bet(bounds%begc:bounds%endc)             ! temporary

    character(len=255)                :: subname ='Tridiagonal_sr'
    !-----------------------------------------------------------------------


    ! Solve the matrix
    if(present(is_col_active))then
       l_is_col_active(:) = is_col_active(:)
    else
       l_is_col_active(:) = .true.
    endif

    do fc = 1,numf
        ci = filter(fc)
        if(l_is_col_active(ci))then
            bet(ci) = b(ci,jtop(ci))
        endif
    end do

    do j = lbj, ubj
       do fc = 1,numf
           ci = filter(fc)
           if(l_is_col_active(ci))then
             if (j >= jtop(ci)) then
               if (j == jtop(ci)) then
                 u(ci,j) = r(ci,j) / bet(ci)
               else
                 gam(ci,j) = c(ci,j-1) / bet(ci)
                 bet(ci) = b(ci,j) - a(ci,j) * gam(ci,j)
                 u(ci,j) = (r(ci,j) - a(ci,j)*u(ci,j-1)) / bet(ci)
               end if
             end if
           endif
        end do
    end do

    do j = ubj-1,lbj,-1
        do fc = 1,numf
           ci = filter(fc)
           if(l_is_col_active(ci))then
             if (j >= jtop(ci)) then
               u(ci,j) = u(ci,j) - gam(ci,j+1) * u(ci,j+1)
             end if
           endif
        end do
    end do


  end subroutine Tridiagonal_sr
  !-----------------------------------------------------------------------
  subroutine Tridiagonal_mr (bounds, lbj, ubj, jtop, numf, filter, ntrcs, a, b, c, r, u, is_col_active)
    !
    ! !DESCRIPTION:
    ! Tridiagonal matrix solution
    ! A X = R
    ! where A, X and R are all matrices.
    ! !USES:
    use shr_kind_mod   , only: r8 => shr_kind_r8
    use clm_varctl     , only : iulog
    use decompMod      , only : bounds_type
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type) , intent(in)    :: bounds                                         ! bounds
    integer           , intent(in)    :: lbj, ubj                                       ! lbinning and ubing level indices
    integer           , intent(in)    :: jtop( bounds%begc: bounds%endc)                ! top level for each column [col]
    integer           , intent(in)    :: numf                                           ! filter dimension
    integer           , intent(in)    :: ntrcs                                          !
    integer           , intent(in)    :: filter(:)                                      ! filter
    real(r8)          , intent(in)    :: a( bounds%begc:bounds%endc , lbj:ubj)          ! "a" left off diagonal of tridiagonal matrix [col , j]
    real(r8)          , intent(in)    :: b( bounds%begc:bounds%endc , lbj:ubj)          ! "b" diagonal column for tridiagonal matrix [col  , j]
    real(r8)          , intent(in)    :: c( bounds%begc:bounds%endc , lbj:ubj)          ! "c" right off diagonal tridiagonal matrix [col   , j]
    real(r8)          , intent(in)    :: r( bounds%begc:bounds%endc , lbj:ubj, 1:ntrcs) ! "r" forcing term of tridiagonal matrix [col , j]
    real(r8)          , intent(inout) :: u( bounds%begc:bounds%endc , lbj:ubj, 1:ntrcs) ! solution [col, j]
                                                                                        !
    integer                           :: j,ci,fc,k                                      ! indices
    logical, optional, intent(in)     :: is_col_active(bounds%begc:bounds%endc)         !
    logical                           :: l_is_col_active(bounds%begc:bounds%endc)       !
    real(r8)                          :: gam(bounds%begc:bounds%endc,lbj:ubj)           ! temporary
    real(r8)                          :: bet(bounds%begc:bounds%endc)                   ! temporary

    character(len=255) :: subname ='Tridiagonal_sr'
    !-----------------------------------------------------------------------

    ! Solve the matrix
    if (present(is_col_active)) then
       l_is_col_active(:) = is_col_active(:)
    else
       l_is_col_active(:) = .true.
    endif

    do fc = 1,numf
       ci = filter(fc)
       if (l_is_col_active(ci))then
          bet(ci) = b(ci,jtop(ci))
       endif
    end do

    do j = lbj, ubj
       do fc = 1,numf
          ci = filter(fc)
          if (l_is_col_active(ci))then
             if (j >= jtop(ci)) then
                if (j == jtop(ci))then
                   do k = 1, ntrcs
                     u(ci,j,k) = r(ci,j,k)/bet(ci)
                   enddo
                else
                   gam(ci,j) = c(ci,j-1) / bet(ci)
                   bet(ci) = b(ci,j) - a(ci,j) * gam(ci,j)
                   do k = 1, ntrcs
                      u(ci,j,k) = (r(ci,j, k) - a(ci,j)*u(ci,j-1, k)) / bet(ci)
                    end do
                 end if
             end if
          end if
       end do
    end do


    do j = ubj-1,lbj,-1
       do fc = 1,numf
          ci = filter(fc)
          if (l_is_col_active(ci)) then
             if (j >= jtop(ci)) then
                do k = 1, ntrcs
                  u(ci,j, k) = u(ci,j, k) - gam(ci,j+1) * u(ci,j+1, k)
                end do
             end if
          end if
       end do
    end do


  end subroutine Tridiagonal_mr

!----------------
  subroutine Trisim(bounds, lbj, ubj, numf, filter, a1,b1,c1,d1,e1,a2,b2,c2,d2,e2,w1, w2)
     !
     !DESCRIPTIONS
     ! This subroutine solves two coupled tridiagonal equations
     ! A1*W1(J-1)+B1*W1(j)+C1*W1(J+1) = D1*W2(j) + E1 AND
     ! A2*W2(J-1)+B2*W2(j)+C2*W2(J+1) = D2*W1(j) + E2
     ! BOUNDARY CONDITIONS ARE APPLIED AT J=1 AND J=JFLT
     ! for the input coefficients array
     ! A1(1)=A2(1)=B1(1)=B2(1)=C1(1)=C2(1)=D1(1)=D2(1)=E1(1)=E2(1)=0
     ! A1(2)=A2(2)=C1(M)=C2(M)=0
     ! the solution is solved at 2 ... M
     ! but I have changed solution location to 1.. M-1 in the code implemented below
     ! Reference:
     ! Deshpande, M. D. and Giddens, D. P. (1977), Direct solution of two linear systems of equations
     ! forming coupled tridiagonal-type matrices. Int. J. Numer. Meth. Engng.,
     ! 11: 1049�1052. doi: 10.1002/nme.1620110612
     !Created by Jinyun Tang
     !Attention: Now the code is specifically written for the soil water coupling with hydraulic
     !redistribution. Idealy, the code should fit well with the purpose of solving coupled heat and
     !water transport equation, but I did not make any attempt here.
     use shr_kind_mod , only : r8 => shr_kind_r8
     use clm_varctl   , only : iulog
     use decompMod    , only : bounds_type
     implicit none
     type(bounds_type) , intent(in)  :: bounds                                 ! bounds
     integer           , intent(in)  :: lbj, ubj                               ! lbinning and ubing level indices
     integer           , intent(in)  :: numf                                   ! filter dimension
     integer           , intent(in)  :: filter(:)                              ! filter
     real(r8)          , intent(in)  :: a1(bounds%begc:bounds%endc, lbj:ubj)   !
     real(r8)          , intent(in)  :: b1(bounds%begc:bounds%endc, lbj:ubj)   !
     real(r8)          , intent(in)  :: c1(bounds%begc:bounds%endc, lbj:ubj)   !
     real(r8)          , intent(in)  :: d1(bounds%begc:bounds%endc, lbj:ubj)   !
     real(r8)          , intent(in)  :: e1(bounds%begc:bounds%endc, lbj:ubj)   !
     real(r8)          , intent(in)  :: a2(bounds%begc:bounds%endc, lbj:ubj)   !
     real(r8)          , intent(in)  :: b2(bounds%begc:bounds%endc, lbj:ubj)   !
     real(r8)          , intent(in)  :: c2(bounds%begc:bounds%endc, lbj:ubj)   !
     real(r8)          , intent(in)  :: d2(bounds%begc:bounds%endc, lbj:ubj)   !
     real(r8)          , intent(in)  :: e2(bounds%begc:bounds%endc, lbj:ubj)   !
     real(r8)          , intent(out) :: w1(bounds%begc:bounds%endc, lbj:ubj-1) !
     real(r8)          , intent(out) :: w2(bounds%begc:bounds%endc, lbj:ubj-1) !
     !local variables
     real(r8)                        :: gam1(bounds%begc:bounds%endc, lbj:ubj) !
     real(r8)                        :: gam2(bounds%begc:bounds%endc, lbj:ubj) !
     real(r8)                        :: xi1(bounds%begc:bounds%endc, lbj:ubj)  !
     real(r8)                        :: xi2(bounds%begc:bounds%endc, lbj:ubj)  !
     real(r8)                        :: phi1(bounds%begc:bounds%endc, lbj:ubj) !
     real(r8)                        :: phi2(bounds%begc:bounds%endc, lbj:ubj) !
     real(r8)                        :: bp1, bp2, dg1, dg2, ee1, ee2           !
     real(r8)                        :: q                                      !
     integer                         :: m, l1, j, j2, k, fc, ci                !
     character(len=255)              :: subname='Trisim'

   m = ubj
     l1= ubj - 1
     j = lbj + 1
     do fc = 1,numf
        ci = filter(fc)

        bp1 = b1(ci,j)
        dg1 = d1(ci,j)
        bp2 = b2(ci,j)
        dg2 = d2(ci,j)

        q = 1._r8 / (bp1 * bp2 - dg1 * dg2)

        phi1(ci,j) = q * bp2 * c1(ci,j)
        phi2(ci,j) = q * bp1 * c2(ci,j)
        gam1(ci,j) =-q * dg1 * c2(ci,j)
        gam2(ci,j) =-q * dg2 * c1(ci,j)

        ee1 = e1(ci,j)
        ee2 = e2(ci,j)

        xi1(ci,j) = q* (bp2 * ee1+dg1 * ee2)
        xi2(ci,j) = q* (dg2 * ee1+bp1 * ee2)
     enddo

     do j= lbj+2 , l1
        do fc = 1, numf
           ci = filter(fc)

           bp1 =b1(ci,j) -a1(ci,j) * phi1(ci,j - 1)
           dg1 =d1(ci,j) -a1(ci,j) * gam1(ci,j - 1)
           bp2 =b2(ci,j) -a2(ci,j) * phi2(ci,j - 1)
           dg2 =d2(ci,j) -a2(ci,j) * gam2(ci,j - 1)

           q = 1._r8/(bp1 * bp2 - dg1 * dg2)

           ee1 = e1(ci,j)-a1(ci,j) * xi1(ci,j - 1)
           ee2 = e2(ci,j)-a2(ci,j) * xi2(ci,j - 1)

           phi1(ci,j) =q * bp2 * c1(ci,j)
           phi2(ci,j) =q * bp1 * c2(ci,j)

           gam1(ci,j) = - q * dg1 * c2(ci,j)
           gam2(ci,j) = - q * dg2 * c1(ci,j)

           xi1(ci,j)  = q * (bp2 * ee1 + dg1 * ee2)
           xi2(ci,j)  = q * (dg2 * ee1 + bp1 * ee2)
        enddo
     enddo
     j=m
     do fc = 1, numf

        ci = filter(fc)

        bp1 =b1(ci,j) -a1(ci,j) * phi1(ci, j - 1)
        dg1 =d1(ci,j) -a1(ci,j) * gam1(ci, j - 1)
        bp2 =b2(ci,j) -a2(ci,j) * phi2(ci, j - 1)
        dg2= d2(ci,j) -a2(ci,j) * gam2(ci, j - 1)

        q = 1.0/(bp1*bp2-dg1*dg2)

        ee1 = e1(ci,j) - a1(ci,j) * xi1(ci,j - 1)
        ee2 = e2(ci,j) - a2(ci,j) * xi2(ci,j - 1)

        xi1(ci,j) = q* ( bp2* ee1+dg1* ee2)
        xi2(ci,j) = q* ( dg2* ee1+bp1* ee2)

        w1(ci,m-1) = xi1(ci,m)
        w2(ci,m-1) = xi2(ci,m)
     enddo

     do j2 = lbj+1, l1
        j=m+ 1-j2
        k =j -1
        do fc = 1, numf
           ci = filter(fc)
           w1 (ci, k )=-phi1(ci, j) *w1(ci, k+ 1) + gam1(ci, j)*w2(ci, k+1)+xi1(ci, j)
           w2 (ci, k )=-phi2(ci, j) *w2(ci, k+ 1) + gam2(ci, j)*w1(ci, k+1)+xi2(ci, j)
        enddo
     enddo

   end subroutine trisim
end module TridiagonalMod
