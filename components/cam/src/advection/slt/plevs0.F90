
subroutine plevs0 (ncol    , ncold   ,nver    ,ps      ,pint    , &
                   pmid    ,pdel)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Define the pressures of the interfaces and midpoints from the
! coordinate definitions and the surface pressure.
! 
! Method: 
! 
! Author: B. Boville
! 
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid,       only: plev, plevp
  use hycoef,       only: hyai, hybi, ps0, hyam, hybm
  implicit none


!-----------------------------------------------------------------------
  integer , intent(in)  :: ncol               ! Longitude dimension
  integer , intent(in)  :: ncold              ! Declared longitude dimension
  integer , intent(in)  :: nver               ! vertical dimension
  real(r8), intent(in)  :: ps(ncold)          ! Surface pressure (pascals)
  real(r8), intent(out) :: pint(ncold,nver+1) ! Pressure at model interfaces
  real(r8), intent(out) :: pmid(ncold,nver)   ! Pressure at model levels
  real(r8), intent(out) :: pdel(ncold,nver)   ! Layer thickness (pint(k+1) - pint(k))
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
  integer i,k             ! Longitude, level indices
!-----------------------------------------------------------------------
!
! Set interface pressures
!
!$OMP PARALLEL DO PRIVATE (K, I)
  do k=1,nver+1
     do i=1,ncol
        pint(i,k) = hyai(k)*ps0 + hybi(k)*ps(i)
     end do
  end do
!
! Set midpoint pressures and layer thicknesses
!
!$OMP PARALLEL DO PRIVATE (K, I)
  do k=1,nver
     do i=1,ncol
        pmid(i,k) = hyam(k)*ps0 + hybm(k)*ps(i)
        pdel(i,k) = pint(i,k+1) - pint(i,k)
     end do
  end do

  return
end subroutine plevs0

