
subroutine diffuse_mom2D_xy
	
!        momentum tendency due to SGS diffusion

use vars
use sgs, only: tk, grdf_x, grdf_z
use params, only: docolumn
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

end subroutine diffuse_mom2D_xy


