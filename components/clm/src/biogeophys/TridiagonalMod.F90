module TridiagonalMod

#include "shr_assert.h"

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
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Tridiagonal (bounds, lbj, ubj, jtop, numf, filter, a, b, c, r, u)
    !
    ! !DESCRIPTION:
    ! Tridiagonal matrix solution
    !
    ! !USES:
    use shr_kind_mod   , only : r8 => shr_kind_r8
    use shr_log_mod    , only : errMsg => shr_log_errMsg
    use clm_varpar     , only : nlevurb
    use column_varcon  , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varctl     , only : iulog
    use decompMod      , only : bounds_type
    use ColumnType     , only : col                
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds             
    integer , intent(in)    :: lbj, ubj                 ! lbinning and ubing level indices
    integer , intent(in)    :: jtop( bounds%begc: )     ! top level for each column [col]
    integer , intent(in)    :: numf                     ! filter dimension
    integer , intent(in)    :: filter(:)                ! filter
    real(r8), intent(in)    :: a( bounds%begc: , lbj: ) ! "a" left off diagonal of tridiagonal matrix [col, j]
    real(r8), intent(in)    :: b( bounds%begc: , lbj: ) ! "b" diagonal column for tridiagonal matrix [col, j]
    real(r8), intent(in)    :: c( bounds%begc: , lbj: ) ! "c" right off diagonal tridiagonal matrix [col, j]
    real(r8), intent(in)    :: r( bounds%begc: , lbj: ) ! "r" forcing term of tridiagonal matrix [col, j]
    real(r8), intent(inout) :: u( bounds%begc: , lbj: ) ! solution [col, j]
    !
    integer  :: j,ci,fc                   !indices
    real(r8) :: gam(bounds%begc:bounds%endc,lbj:ubj)      !temporary
    real(r8) :: bet(bounds%begc:bounds%endc)              !temporary
    !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(jtop) == (/bounds%endc/)),      errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(a)    == (/bounds%endc, ubj/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(b)    == (/bounds%endc, ubj/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(c)    == (/bounds%endc, ubj/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(r)    == (/bounds%endc, ubj/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(u)    == (/bounds%endc, ubj/)), errMsg(__FILE__, __LINE__))

    ! Solve the matrix

    do fc = 1,numf
       ci = filter(fc)
       bet(ci) = b(ci,jtop(ci))
    end do

    do j = lbj, ubj
       do fc = 1,numf
          ci = filter(fc)
          if ((col%itype(ci) == icol_sunwall .or. col%itype(ci) == icol_shadewall &
              .or. col%itype(ci) == icol_roof) .and. j <= nlevurb) then
             if (j >= jtop(ci)) then
                if (j == jtop(ci)) then
                   u(ci,j) = r(ci,j) / bet(ci)
                else
                   gam(ci,j) = c(ci,j-1) / bet(ci)
                   bet(ci) = b(ci,j) - a(ci,j) * gam(ci,j)
                   u(ci,j) = (r(ci,j) - a(ci,j)*u(ci,j-1)) / bet(ci)
                end if
             end if
          else if (col%itype(ci) /= icol_sunwall .and. col%itype(ci) /= icol_shadewall &
                   .and. col%itype(ci) /= icol_roof) then
             if (j >= jtop(ci)) then
                if (j == jtop(ci)) then
                   u(ci,j) = r(ci,j) / bet(ci)
                else
                   gam(ci,j) = c(ci,j-1) / bet(ci)
                   bet(ci) = b(ci,j) - a(ci,j) * gam(ci,j)
                   u(ci,j) = (r(ci,j) - a(ci,j)*u(ci,j-1)) / bet(ci)
                end if
             end if
          end if
       end do
    end do

    do j = ubj-1,lbj,-1
       do fc = 1,numf
          ci = filter(fc)
          if ((col%itype(ci) == icol_sunwall .or. col%itype(ci) == icol_shadewall &
              .or. col%itype(ci) == icol_roof) .and. j <= nlevurb-1) then
             if (j >= jtop(ci)) then
                u(ci,j) = u(ci,j) - gam(ci,j+1) * u(ci,j+1)
             end if
          else if (col%itype(ci) /= icol_sunwall .and. col%itype(ci) /= icol_shadewall &
                   .and. col%itype(ci) /= icol_roof) then
             if (j >= jtop(ci)) then
                u(ci,j) = u(ci,j) - gam(ci,j+1) * u(ci,j+1)
             end if
          end if
       end do
    end do

  end subroutine Tridiagonal

end module TridiagonalMod
