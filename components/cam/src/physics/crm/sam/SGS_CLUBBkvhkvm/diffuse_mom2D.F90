
subroutine diffuse_mom2D
	
!        momentum tendency due to SGS diffusion

use vars
use sgs, only: tk, grdf_x, grdf_z
use params, only: docolumn
#ifdef CLUBB_CRM
use params, only: doclubb, doclubb_sfc_fluxes, docam_sfc_fluxes
#endif
implicit none

real rdx2,rdz2,rdz,rdx25,rdz25,rdx21,rdx251
real dxz,dzx

integer i,j,k,ic,ib,kc,kcu
real tkx, tkz, rhoi, iadzw, iadz
real fu(0:nx,1,nz),fv(0:nx,1,nz),fw(0:nx,1,nz)

rdx2=1./dx/dx
rdx25=0.25*rdx2

dxz=dx/dz

j=1

if(.not.docolumn) then


do k=1,nzm

 kc=k+1
 kcu=min(kc,nzm)
 dxz=dx/(dz*adzw(kc))
 rdx21=rdx2 * grdf_x(k)
 rdx251=rdx25 * grdf_x(k)
 
   do i=0,nx
    ic=i+1
    tkx=rdx21*tk(i,j,k)
    fu(i,j,k)=-2.*tkx*(u(ic,j,k)-u(i,j,k))
    fv(i,j,k)=-tkx*(v(ic,j,k)-v(i,j,k))
    tkx=rdx251*(tk(i,j,k)+tk(ic,j,k)+tk(i,j,kcu)+tk(ic,j,kcu)) 	
    fw(i,j,k)=-tkx*(w(ic,j,kc)-w(i,j,kc)+(u(ic,j,kcu)-u(ic,j,k))*dxz)
   end do 
   do i=1,nx
    ib=i-1
    dudt(i,j,k,na)=dudt(i,j,k,na)-(fu(i,j,k)-fu(ib,j,k))
    dvdt(i,j,k,na)=dvdt(i,j,k,na)-(fv(i,j,k)-fv(ib,j,k))
    dwdt(i,j,kc,na)=dwdt(i,j,kc,na)-(fw(i,j,k)-fw(ib,j,k))
   end do  

end do 

end if 

!-------------------------
rdz=1./dz
dzx=dz/dx

do k=1,nzm-1
 kc=k+1
 uwsb(kc)=0.
 vwsb(kc)=0.
 iadz = 1./adz(k)
 iadzw= 1./adzw(kc)
 rdz2=rdz*rdz *grdf_z(k)
 rdz25=0.25*rdz2
   do i=1,nx
    ib=i-1
    tkz=rdz2*tk(i,j,k)
    fw(i,j,kc)=-2.*tkz*(w(i,j,kc)-w(i,j,k))*rho(k)*iadz
    tkz=rdz25*(tk(i,j,k)+tk(ib,j,k)+tk(i,j,kc)+tk(ib,j,kc))
    fu(i,j,kc)=-tkz*( (u(i,j,kc)-u(i,j,k))*iadzw + &
                       (w(i,j,kc)-w(ib,j,kc))*dzx)*rhow(kc) 	
    fv(i,j,kc)=-tkz*(v(i,j,kc)-v(i,j,k))*iadzw*rhow(kc)
    uwsb(kc)=uwsb(kc)+fu(i,j,kc)
    vwsb(kc)=vwsb(kc)+fv(i,j,kc)
  end do 
end do

uwsb(1) = 0.
vwsb(1) = 0.
	
do i=1,nx
 tkz=rdz2*grdf_z(nzm)*tk(i,j,nzm)
 fw(i,j,nz)=-2.*tkz*(w(i,j,nz)-w(i,j,nzm))/adz(nzm)*rho(nzm)
#ifdef CLUBB_CRM
 ! Add in the surface flux later -dschanen UWM 27 Aug 2008
 if ( doclubb .and. (doclubb_sfc_fluxes .or. docam_sfc_fluxes) ) then
   fu(i,j,1) = 0.0  ! This is handled by CLUBB
   fv(i,j,1) = 0.0  ! " "
 else
   fu(i,j,1)=fluxbu(i,j) * rdz * rhow(1)
   fv(i,j,1)=fluxbv(i,j) * rdz * rhow(1)
 end if
#else
 fu(i,j,1)=fluxbu(i,j) * rdz * rhow(1)
 fv(i,j,1)=fluxbv(i,j) * rdz * rhow(1)
#endif /*CLUBB_CRM*/
 fu(i,j,nz)=fluxtu(i,j) * rdz * rhow(nz)
 fv(i,j,nz)=fluxtv(i,j) * rdz * rhow(nz)
 uwsb(1) = uwsb(1) + fu(i,j,1)
 vwsb(1) = vwsb(1) + fv(i,j,1)
end do
	 

do k=1,nzm
  kc=k+1
  rhoi = 1./(rho(k)*adz(k))
  do i=1,nx
    dudt(i,j,k,na)=dudt(i,j,k,na)-(fu(i,j,kc)-fu(i,j,k))*rhoi
    dvdt(i,j,k,na)=dvdt(i,j,k,na)-(fv(i,j,kc)-fv(i,j,k))*rhoi
  end do
end do ! k

do k=2,nzm
  rhoi = 1./(rhow(k)*adzw(k))
  do i=1,nx	 
    dwdt(i,j,k,na)=dwdt(i,j,k,na)-(fw(i,j,k+1)-fw(i,j,k))*rhoi
  end do
end do ! k


end subroutine diffuse_mom2D


