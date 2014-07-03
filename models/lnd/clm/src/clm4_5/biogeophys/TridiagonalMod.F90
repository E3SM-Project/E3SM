module TridiagonalMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: TridiagonalMod
!
! !DESCRIPTION:
! Tridiagonal matrix solution
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: Tridiagonal
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Tridiagonal
!
! !INTERFACE:
  subroutine Tridiagonal (lbc, ubc, lbj, ubj, jtop, numf, filter, &
                          a, b, c, r, u)
!
! !DESCRIPTION:
! Tridiagonal matrix solution
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use clm_varpar    , only : nlevurb
    use clm_varcon    , only : icol_roof, icol_sunwall, icol_shadewall
    use clm_varctl    , only : iulog
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)    :: lbc, ubc               ! lbinning and ubing column indices
    integer , intent(in)    :: lbj, ubj               ! lbinning and ubing level indices
    integer , intent(in)    :: jtop(lbc:ubc)          ! top level for each column
    integer , intent(in)    :: numf                   ! filter dimension
    integer , intent(in)    :: filter(1:numf)         ! filter
    real(r8), intent(in)    :: a(lbc:ubc, lbj:ubj)    ! "a" left off diagonal of tridiagonal matrix
    real(r8), intent(in)    :: b(lbc:ubc, lbj:ubj)    ! "b" diagonal column for tridiagonal matrix
    real(r8), intent(in)    :: c(lbc:ubc, lbj:ubj)    ! "c" right off diagonal tridiagonal matrix
    real(r8), intent(in)    :: r(lbc:ubc, lbj:ubj)    ! "r" forcing term of tridiagonal matrix
    real(r8), intent(inout) :: u(lbc:ubc, lbj:ubj)    ! solution
! local pointers to original implicit in arguments
!
    integer , pointer :: ctype(:)           ! column type

!
! !CALLED FROM:
! subroutine BiogeophysicsLake in module BiogeophysicsLakeMod
! subroutine SoilTemperature in module SoilTemperatureMod
! subroutine SoilWater in module HydrologyMod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
!  1 July 2003: Mariana Vertenstein; modified for vectorization
!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
    integer  :: j,ci,fc                   !indices
    real(r8) :: gam(lbc:ubc,lbj:ubj)      !temporary
    real(r8) :: bet(lbc:ubc)              !temporary
!-----------------------------------------------------------------------

    ! Assign local pointers to derived subtypes components (column-level)

    ctype          => clm3%g%l%c%itype

    ! Solve the matrix

    do fc = 1,numf
       ci = filter(fc)
       bet(ci) = b(ci,jtop(ci))
    end do

    do j = lbj, ubj
       do fc = 1,numf
          ci = filter(fc)
          if ((ctype(ci) == icol_sunwall .or. ctype(ci) == icol_shadewall &
              .or. ctype(ci) == icol_roof) .and. j <= nlevurb) then
             if (j >= jtop(ci)) then
                if (j == jtop(ci)) then
                   u(ci,j) = r(ci,j) / bet(ci)
                else
                   gam(ci,j) = c(ci,j-1) / bet(ci)
                   bet(ci) = b(ci,j) - a(ci,j) * gam(ci,j)
                   u(ci,j) = (r(ci,j) - a(ci,j)*u(ci,j-1)) / bet(ci)
                end if
             end if
          else if (ctype(ci) /= icol_sunwall .and. ctype(ci) /= icol_shadewall &
                   .and. ctype(ci) /= icol_roof) then
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
          if ((ctype(ci) == icol_sunwall .or. ctype(ci) == icol_shadewall &
              .or. ctype(ci) == icol_roof) .and. j <= nlevurb-1) then
             if (j >= jtop(ci)) then
                u(ci,j) = u(ci,j) - gam(ci,j+1) * u(ci,j+1)
             end if
          else if (ctype(ci) /= icol_sunwall .and. ctype(ci) /= icol_shadewall &
                   .and. ctype(ci) /= icol_roof) then
             if (j >= jtop(ci)) then
                u(ci,j) = u(ci,j) - gam(ci,j+1) * u(ci,j+1)
             end if
          end if
       end do
    end do

  end subroutine Tridiagonal

end module TridiagonalMod
