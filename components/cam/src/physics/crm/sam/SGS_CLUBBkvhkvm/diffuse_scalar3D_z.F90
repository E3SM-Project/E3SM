subroutine diffuse_scalar3D_z (field,fluxb,fluxt,tkh,rho,rhow,flux)

use grid
use params, only: docolumn,dowallx,dowally,dosgs
use sgs, only: grdf_x,grdf_y,grdf_z
implicit none
! input	
real field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)	! scalar
real tkh(0:nxp1,1-YES3D:nyp1,nzm)	! eddy conductivity
real fluxb(nx,ny)		! bottom flux
real fluxt(nx,ny)		! top flux
real rho(nzm)
real rhow(nz)
real flux(nz)
! local        
real flx(0:nx,0:ny,0:nzm)
real dfdt(nx,ny,nz)
real rdx2,rdy2,rdz2,rdz,rdx5,rdy5,rdz5,tmp
real dxy,dxz,dyx,dyz,dzx,dzy,tkx,tky,tkz,rhoi
integer i,j,k,ib,ic,jb,jc,kc,kb


if(.not.dosgs) return

rdx2=1./(dx*dx)
rdy2=1./(dy*dy)
rdz2=1./(dz*dz)
rdz=1./dz
dxy=dx/dy
dxz=dx/dz
dyx=dy/dx
dyz=dy/dz
dzx=dz/dx
dzy=dz/dy

dfdt(:,:,:)=0.

!  Vertical diffusion:

flux(1) = 0.
tmp=1./adzw(nz)
do j=1,ny
 do i=1,nx	
   flx(i,j,0)=fluxb(i,j)*rdz*rhow(1)
   flx(i,j,nzm)=fluxt(i,j)*rdz*tmp*rhow(nz)
   flux(1) = flux(1) + flx(i,j,0)
 end do
end do


do k=1,nzm-1
 kc=k+1	
 flux(kc)=0. 
 rhoi = rhow(kc)/adzw(kc)
 rdz5=0.5*rdz2 * grdf_z(k)
 do j=1,ny
  do i=1,nx
    tkz=rdz5*(tkh(i,j,k)+tkh(i,j,kc))
    flx(i,j,k)=-tkz*(field(i,j,kc)-field(i,j,k))*rhoi
    flux(kc) = flux(kc) + flx(i,j,k)
  end do 
 end do
end do

do k=1,nzm
 kb=k-1
 rhoi = 1./(adz(k)*rho(k))
 do j=1,ny
  do i=1,nx		 
   dfdt(i,j,k)=dtn*(dfdt(i,j,k)-(flx(i,j,k)-flx(i,j,kb))*rhoi)
   field(i,j,k)=field(i,j,k)+dfdt(i,j,k)
  end do 
 end do 	 
end do 

end subroutine diffuse_scalar3D_z
