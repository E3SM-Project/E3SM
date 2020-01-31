subroutine diffuse_scalar3D_xy (field,fluxb,fluxt,tkh,rho,rhow,flux)

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

!-----------------------------------------
if(dowallx) then

  if(mod(rank,nsubdomains_x).eq.0) then
    do k=1,nzm
     do j=1,ny
         field(0,j,k) = field(1,j,k)
     end do
    end do
  end if
  if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
    do k=1,nzm
     do j=1,ny
         field(nx+1,j,k) = field(nx,j,k)
     end do
    end do
  end if

end if

if(dowally) then

  if(rank.lt.nsubdomains_x) then
    do k=1,nzm
       do i=1,nx
         field(i,1-YES3D,k) = field(i,1,k)
       end do
    end do
  end if
  if(rank.gt.nsubdomains-nsubdomains_x-1) then
    do k=1,nzm
       do i=1,ny
         field(i,ny+YES3D,k) = field(i,ny,k)
       end do
    end do
  end if

end if



if(dowally) then

 call task_rank_to_index(rank, ib, jb)
 if(jb.eq.0) then
   do k=1,nzm
     do i=1,nx
       field(i,1-YES3D,k) = field(i,1,k)
     end do
   end do
 end if
 if(jb.eq.nsubdomains_y-1) then
   do k=1,nzm
     do i=1,nx
       field(i,ny+YES3D,k) = field(i,ny,k)
     end do
   end do
 end if

end if

!-----------------------------------------


!  Horizontal diffusion:


do k=1,nzm
	
 rdx5=0.5*rdx2  * grdf_x(k)
 rdy5=0.5*rdy2  * grdf_y(k)

 do j=1,ny
  do i=0,nx
    ic=i+1
    tkx=rdx5*(tkh(i,j,k)+tkh(ic,j,k)) 	
    flx(i,j,k)=-tkx*(field(ic,j,k)-field(i,j,k))
  end do 
  do i=1,nx
    ib=i-1
    dfdt(i,j,k)=dfdt(i,j,k)-(flx(i,j,k)-flx(ib,j,k))
  end do 
 end do 

 do j=0,ny
  jc=j+1
  do i=1,nx
   tky=rdy5*(tkh(i,j,k)+tkh(i,jc,k)) 	
   flx(i,j,k)=-tky*(field(i,jc,k)-field(i,j,k))
  end do 
 end do
 do j=1,ny
  jb=j-1
  do i=1,nx	    
    dfdt(i,j,k)=dfdt(i,j,k)-(flx(i,j,k)-flx(i,jb,k))
  end do 
 end do 

 do j=1, ny
  do i=1, nx
   field(i,j,k) = field(i,j,k) + dfdt(i,j,k) * dtn
  end do
 end do
 
end do ! k

flux = 0.0

end subroutine diffuse_scalar3D_xy
