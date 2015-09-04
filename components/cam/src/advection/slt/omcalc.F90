
subroutine omcalc(rcoslat ,d       ,u       ,v       ,dpsl    , &
                  dpsm    ,pmid    ,pdel    ,rpmid   ,pbot    , &
                  omga    ,nlon    )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Calculate vertical pressure velocity (omga = dp/dt)
! 
! Method: 
! First evaluate the expressions for omega/p, then rescale to omega at
! the end.
! 
! Author: CCM1
! 
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid,       only: plev, plon, plevp
  use pspect
  use hycoef,       only: hybm, hybd, nprlev
  implicit none


!------------------------------Arguments--------------------------------
  integer , intent(in) :: nlon                 ! lonitude dimension
  real(r8), intent(in) :: rcoslat(nlon)        ! 1 / cos(lat)
  real(r8), intent(in) :: d(plon,plev)         ! divergence
  real(r8), intent(in) :: u(plon,plev)         ! zonal wind * cos(lat)
  real(r8), intent(in) :: v(plon,plev)         ! meridional wind * cos(lat)
  real(r8), intent(in) :: dpsl(plon)           ! longitudinal component of grad ln(ps)
  real(r8), intent(in) :: dpsm(plon)           ! latitudinal component of grad ln(ps)
  real(r8), intent(in) :: pmid(plon,plev)      ! mid-level pressures
  real(r8), intent(in) :: pdel(plon,plev)      ! layer thicknesses (pressure)
  real(r8), intent(in) :: rpmid(plon,plev)     ! 1./pmid
  real(r8), intent(in) :: pbot(plon)           ! bottom interface pressure
  real(r8), intent(out):: omga(plon,plev)      ! vertical pressure velocity
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
  integer i,k               ! longitude, level indices
  real(r8) d_i(plev)        ! divergence (single colummn)
  real(r8) u_i(plev)        ! zonal wind * cos(lat) (single colummn)
  real(r8) v_i(plev)        ! meridional wind * cos(lat) (single colummn)
  real(r8) pmid_i(plev)     ! mid-level pressures (single colummn)
  real(r8) pdel_i(plev)     ! layer thicknesses (pressure) (single colummn)
  real(r8) rpmid_i(plev)    ! 1./pmid (single colummn)
  real(r8) omga_i(plev)     ! vertical pressure velocity (single colummn)
  real(r8) hkk              ! diagonal element of hydrostatic matrix
  real(r8) hlk              ! super diagonal element
  real(r8) suml             ! partial sum over l = (1, k-1)
  real(r8) vgpk             ! v dot grad ps
  real(r8) tmp              ! vector temporary
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE (I, SUML, D_I, U_I, V_I, PMID_I, PDEL_I, RPMID_I, &
!$OMP                      OMGA_I, HKK, VGPK, TMP, HLK)
  do i=1,nlon
!
! Zero partial sum
!
     suml = 0._r8
!
! Collect column data
!
     d_i = d(i,:)
     u_i = u(i,:)
     v_i = v(i,:)
     pmid_i = pmid(i,:)
     pdel_i = pdel(i,:)
     rpmid_i = rpmid(i,:)
!
! Pure pressure part: top level
!
     hkk = 0.5_r8*rpmid_i(1)
     omga_i(1) = -hkk*d_i(1)*pdel_i(1)
     suml = suml + d_i(1)*pdel_i(1)
!
! sum(k)(v(j)*ps*grad(lnps)*db(j)) part. Not normally invoked since 
! the top layer is normally a pure pressure layer.
!
     if (1>=nprlev) then
        vgpk = rcoslat(i)*(u_i(1)*dpsl(i) + v_i(1)*dpsm(i))*pbot(i)
        tmp = vgpk*hybd(1)
        omga_i(1) = omga_i(1) + hybm(1)*rpmid_i(1)*vgpk - hkk*tmp
        suml = suml + tmp
     end if
!
! Integrals to level above bottom
!
     do k=2,plev-1
!
! Pure pressure part
!
        hkk = 0.5_r8*rpmid_i(k)
        hlk = rpmid_i(k)
        omga_i(k) = -hkk*d_i(k)*pdel_i(k) - hlk*suml
        suml = suml + d_i(k)*pdel_i(k)
!
! v(j)*grad(lnps) part
!
        if (k>=nprlev) then
           vgpk = rcoslat(i)*(u_i(k)*dpsl(i) + v_i(k)*dpsm(i))*pbot(i)
           tmp = vgpk*hybd(k)
           omga_i(k) = omga_i(k) + hybm(k)*rpmid_i(k)*vgpk - hkk*tmp
           suml = suml + tmp
        end if
     end do
!
! Pure pressure part: bottom level
!
     hkk = 0.5_r8*rpmid_i(plev)
     hlk =     rpmid_i(plev)
     omga_i(plev) = -hkk*d_i(plev)*pdel_i(plev) - hlk*suml
!
! v(j)*grad(lnps) part. Normally invoked, but omitted if the model is
! running in pure pressure coordinates throughout (e.g. stratospheric 
! mechanistic model).
!
     if (plev>=nprlev) then
        vgpk = rcoslat(i)*(u_i(plev)*dpsl(i) + v_i(plev)*dpsm(i))* pbot(i)
        omga_i(plev) = omga_i(plev) + hybm(plev)*rpmid_i(plev)*vgpk - &
             hkk*vgpk*hybd(plev)
     end if
!
! The above expressions give omega/p. Rescale to omega.
!
     do k=1,plev
        omga_i(k) = omga_i(k)*pmid_i(k)
     end do
!
! Save results
!
     omga(i,:) = omga_i(:)
!
   end do
!
  return
end subroutine omcalc

