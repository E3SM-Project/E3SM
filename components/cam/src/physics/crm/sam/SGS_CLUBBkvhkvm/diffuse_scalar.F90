subroutine diffuse_scalar (f,fluxb,fluxt, &
                          fdiff,flux,f2lediff,f2lediss,fwlediff,doit)

use grid
use vars, only: rho, rhow
use sgs, only: tkh
implicit none

! input:	
real f(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)	! scalar
real fluxb(nx,ny)		! bottom flux
real fluxt(nx,ny)		! top flux
real flux(nz)
real f2lediff(nz),f2lediss(nz),fwlediff(nz)
real fdiff(nz)
logical doit
! Local
real df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)	! scalar
integer i,j,k

!call t_startf ('diffuse_scalars')

df(:,:,:) = f(:,:,:)

if(RUN3D) then
!  call diffuse_scalar3D (f,fluxb,fluxt,tkh,rho,rhow,flux)
  call diffuse_scalar3D_xy (f,fluxb,fluxt,tkh,rho,rhow,flux)
  call diffuse_scalar3D_z (f,fluxb,fluxt,tkh,rho,rhow,flux)
else  
!  call diffuse_scalar2D (f,fluxb,fluxt,tkh,rho,rhow,flux)
  call diffuse_scalar2D_xy (f,fluxb,fluxt,tkh,rho,rhow,flux)
  call diffuse_scalar2D_z (f,fluxb,fluxt,tkh,rho,rhow,flux)
endif

do k=1,nzm
   fdiff(k)=0.
   do j=1,ny
    do i=1,nx
     fdiff(k)=fdiff(k)+f(i,j,k)-df(i,j,k)
    end do
   end do
end do

!call t_stopf ('diffuse_scalars')

end subroutine diffuse_scalar 
