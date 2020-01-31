subroutine diffuse_scalar_z (f,fluxb,fluxt, &
                          fdiff,flux,f2lediff,f2lediss,fwlediff,doit)

use grid
use vars, only: rho, rhow
use sgs, only: tkh
#ifdef CLUBB_CRM
use sgs, only: tkh_clubb 
use params, only: doclubb
#endif
implicit none

! input:	
real f(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)	! scalar
real fluxb(nx,ny)		! bottom flux
real fluxt(nx,ny)		! top flux
real flux(nz)
real fdiff(nz)
real f2lediff(nz)
real f2lediss(nz)
real fwlediff(nz)
logical doit
! Local
real df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)	! scalar
real tkhtemp(0:nxp1, 1-YES3D:nyp1, nzm) ! eddy diffusivity
real f0(nzm),df0(nzm),factor_xy
real f2lediss_z(nzm)
real r2dx,r2dy,r2dx0,r2dy0,r2dz
integer i,j,k,kb,kc,jb,jc

!call t_startf ('diffuse_scalars_z')

tkhtemp = 0.0
#ifndef  CLUBB_CRM
tkhtemp(:, :, :) = tkh(:, :, :)
#else
if(doclubb) then
  tkhtemp(:, :, :) = tkh_clubb(:, :, :)
else
  tkhtemp(:, :, :) = tkh(:, :, :)
endif
#endif
	
  do k=1,nzm
    do j=dimy1_s,dimy2_s
     do i=dimx1_s,dimx2_s
      df(i,j,k) = f(i,j,k)
     end do
    end do
  end do


if(RUN3D) then
  call diffuse_scalar3D_z (f,fluxb,fluxt,tkhtemp,rho,rhow,flux)
else  
  call diffuse_scalar2D_z (f,fluxb,fluxt,tkhtemp,rho,rhow,flux)
endif
	
  do k=1,nzm
    fdiff(k)=0.
    do j=1,ny
     do i=1,nx
      fdiff(k)=fdiff(k)+f(i,j,k)-df(i,j,k)
     end do
    end do
  end do

!call t_stopf ('diffuse_scalars_z')

end subroutine diffuse_scalar_z
