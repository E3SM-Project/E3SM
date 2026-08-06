#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module viscosity_mod
  use viscosity_base, only: compute_zeta_C0, compute_div_C0, &
    compute_zeta_C0_contra, compute_eta_C0_contra, compute_div_C0_contra, make_c0, neighbor_minmax, &
    make_c0_vector

use dimensions_mod, only : np, nlev,qsize,nelemd
use hybrid_mod, only : hybrid_t
use element_mod, only : element_t
use kinds, only : real_kind, iulog
use edgetype_mod, only : EdgeBuffer_t
use edge_mod, only : edgevpack, edgevunpack
use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk, vorticity_sphere, derivinit
use bndry_mod, only : bndry_exchangev
use control_mod, only : hypervis_scaling, nu, nu_div
use parallel_mod, only : parallel_t
use physical_constants, only : g


implicit none

public :: compute_pv_C0_contra

contains


subroutine biharmonic_wk(elem,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute weak biharmonic operator
!    input:  h,v (stored in elem()%, in lat-lon coordinates
!    output: ptens,vtens  overwritten with weak biharmonic of h,v (output in lat-lon coordinates)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(inout), target :: elem(:)
integer :: nt,nets,nete
real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)  :: vtens
real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: ptens
type (EdgeBuffer_t)  , intent(inout) :: edge3
type (derivative_t)  , intent(in) :: deriv

! local
integer :: k,kptr,i,j,ie,ic
real (kind=real_kind), dimension(:,:), pointer :: rspheremv
real (kind=real_kind), dimension(np,np) :: lap_ps
real (kind=real_kind), dimension(np,np,nlev) :: T
real (kind=real_kind), dimension(np,np,2) :: v
real (kind=real_kind) ::  nu_ratio1,nu_ratio2
logical var_coef1

   !if tensor hyperviscosity with tensor V is used, then biharmonic operator is (\grad\cdot V\grad) (\grad \cdot \grad) 
   !so tensor is only used on second call to laplace_sphere_wk
   var_coef1 = .true.
   if(hypervis_scaling > 0)  var_coef1= .false.

   ! note: there is a scaling bug in the treatment of nu_div
   ! nu_ratio is applied twice, once in each laplace operator
   ! so in reality:   nu_div_actual = (nu_div/nu)**2 nu
   ! We should fix this, but it requires adjusting all CAM defaults
   nu_ratio1=1
   nu_ratio2=1
   if (nu_div/=nu) then
      if(hypervis_scaling /= 0) then
         ! we have a problem with the tensor in that we cant seperate
         ! div and curl components.  So we do, with tensor V:
         ! nu * (del V del ) * ( nu_ratio * grad(div) - curl(curl))
         nu_ratio1=(nu_div/nu)**2   ! preserve buggy scaling
         nu_ratio2=1
      else
         nu_ratio1=nu_div/nu
         nu_ratio2=nu_div/nu
      endif
   endif


   do ie=nets,nete
      do k=1,nlev
         do j=1,np
            do i=1,np
               ! filter surface height, not thickness
               T(i,j,k)=elem(ie)%state%p(i,j,k,nt) + elem(ie)%state%ps(i,j)
            enddo
         enddo
        
         ptens(:,:,k,ie)=laplace_sphere_wk(T(:,:,k),deriv,elem(ie),var_coef=var_coef1)
         vtens(:,:,:,k,ie)=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,&
              elem(ie),var_coef=var_coef1,nu_ratio=nu_ratio1)

      enddo
      kptr=0
      call edgeVpack(edge3, ptens(1,1,1,ie),nlev,kptr,ie)
      kptr=nlev
      call edgeVpack(edge3, vtens(1,1,1,1,ie),2*nlev,kptr,ie)

   enddo
   
   call bndry_exchangeV(hybrid,edge3)
   
   do ie=nets,nete
      rspheremv     => elem(ie)%rspheremp(:,:)
      
      kptr=0
      call edgeVunpack(edge3, ptens(1,1,1,ie), nlev, kptr, ie)
      kptr=nlev
      call edgeVunpack(edge3, vtens(1,1,1,1,ie), 2*nlev, kptr, ie)
      
      ! apply inverse mass matrix, then apply laplace again
      do k=1,nlev
         do j=1,np
            do i=1,np
               T(i,j,k)=rspheremv(i,j)*ptens(i,j,k,ie)
               v(i,j,1)=rspheremv(i,j)*vtens(i,j,1,k,ie)
               v(i,j,2)=rspheremv(i,j)*vtens(i,j,2,k,ie)
            enddo
         enddo
         ptens(:,:,k,ie)=laplace_sphere_wk(T(:,:,k),deriv,elem(ie),var_coef=.true.)
         vtens(:,:,:,k,ie)=vlaplace_sphere_wk(v(:,:,:),deriv,elem(ie),var_coef=.true.,&
              nu_ratio=nu_ratio2)
      enddo
   enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine


subroutine compute_pv_C0_contra(pv,elem,phimean,par,nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 potential vorticity for swe.  That is, solve:
!     < PHI, pv > = <PHI, (curl(elem%state%v) + coriolis)/h >
!
!    input:  v (stored in elem()%, in contra-variant coordinates)
!    output: pv(:,:,:,:)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (parallel_t)      , intent(in) :: par
type (element_t)     , intent(in), target :: elem(:)
integer :: nt
real (kind=real_kind) :: phimean
real (kind=real_kind), dimension(np,np,nlev,nelemd) :: pv
real (kind=real_kind), dimension(np,np,2) :: ulatlon
real (kind=real_kind), dimension(np,np) :: v1,v2

! local
integer :: k,ie
type (derivative_t)          :: deriv

call derivinit(deriv)

do k=1,nlev
do ie=1,nelemd
    v1 = elem(ie)%state%v(:,:,1,k,nt)
    v2 = elem(ie)%state%v(:,:,2,k,nt)
    ulatlon(:,:,1) = elem(ie)%D(:,:,1,1)*v1 + elem(ie)%D(:,:,1,2)*v2
    ulatlon(:,:,2) = elem(ie)%D(:,:,2,1)*v1 + elem(ie)%D(:,:,2,2)*v2
   pv(:,:,k,ie)=vorticity_sphere(ulatlon,deriv,elem(ie))
   pv(:,:,k,ie)=pv(:,:,k,ie) + elem(ie)%fcor(:,:)
   pv(:,:,k,ie)=pv(:,:,k,ie)/(elem(ie)%state%p(:,:,k,nt) + elem(ie)%state%ps(:,:) + phimean)*g
enddo
enddo

call make_C0(pv,elem,par)

end subroutine

end module viscosity_mod
