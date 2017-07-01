#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module viscosity_preqx_base
!
!  This module should be renamed "global_deriv_mod.F90"
! 
!  It is a collection of derivative operators that must be applied to the field 
!  over the sphere (as opposed to derivative operators that can be applied element 
!  by element)
!
!
use kinds, only : real_kind, iulog
use dimensions_mod, only : np, nlev
use hybrid_mod, only : hybrid_t
use parallel_mod, only : parallel_t
use element_mod, only : element_t
use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk
use edgetype_mod, only : EdgeBuffer_t, EdgeDescriptor_t
use edge_mod, only : edgevpack, edgevunpack

use bndry_mod, only : bndry_exchangev
use control_mod, only : hypervis_scaling, nu, nu_div
use perf_mod, only: t_startf, t_stopf

implicit none
private
save


public :: biharmonic_wk_dp3d


contains


subroutine biharmonic_wk_dp3d(elem,dptens,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute weak biharmonic operator
!    input:  h,v (stored in elem()%, in lat-lon coordinates
!    output: ptens,vtens  overwritten with weak biharmonic of h,v (output in lat-lon coordinates)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(inout), target :: elem(:)
integer              , intent(in)  :: nt,nets,nete
real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)  :: vtens
real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: ptens,dptens
type (EdgeBuffer_t)  , intent(inout) :: edge3
type (derivative_t)  , intent(in) :: deriv

! local
integer :: i,j,k,kptr,ie
real (kind=real_kind), dimension(:,:), pointer :: rspheremv
real (kind=real_kind), dimension(np,np) :: tmp
real (kind=real_kind), dimension(np,np) :: tmp2
real (kind=real_kind), dimension(np,np,2) :: v
real (kind=real_kind) :: nu_ratio1, nu_ratio2
logical var_coef1

   !if tensor hyperviscosity with tensor V is used, then biharmonic operator is (\grad\cdot V\grad) (\grad \cdot \grad) 
   !so tensor is only used on second call to laplace_sphere_wk
   var_coef1 = .true.
   if(hypervis_scaling > 0)    var_coef1 = .false.

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

#if (defined COLUMN_OPENMP)
!$omp parallel do default(shared), private(k,tmp)
#endif
      do k=1,nlev
         tmp=elem(ie)%state%T(:,:,k,nt) 
         ptens(:,:,k,ie)=laplace_sphere_wk(tmp,deriv,elem(ie),var_coef=var_coef1)
         tmp=elem(ie)%state%dp3d(:,:,k,nt) 
         dptens(:,:,k,ie)=laplace_sphere_wk(tmp,deriv,elem(ie),var_coef=var_coef1)
         vtens(:,:,:,k,ie)=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),&
              var_coef=var_coef1,nu_ratio=nu_ratio1)
      enddo
      kptr=0
      call edgeVpack(edge3, ptens(1,1,1,ie),nlev,kptr,ie)
      kptr=nlev
      call edgeVpack(edge3, vtens(1,1,1,1,ie),2*nlev,kptr,ie)
      kptr=3*nlev
      call edgeVpack(edge3, dptens(1,1,1,ie),nlev,kptr,ie)

   enddo
   
   call t_startf('biwkdp3d_bexchV')
   call bndry_exchangeV(hybrid,edge3)
   call t_stopf('biwkdp3d_bexchV')
   
   do ie=nets,nete
      rspheremv     => elem(ie)%rspheremp(:,:)
      
      kptr=0
      call edgeVunpack(edge3, ptens(1,1,1,ie), nlev, kptr, ie)
      kptr=nlev
      call edgeVunpack(edge3, vtens(1,1,1,1,ie), 2*nlev, kptr, ie)
      kptr=3*nlev
      call edgeVunpack(edge3, dptens(1,1,1,ie), nlev, kptr, ie)
      
      ! apply inverse mass matrix, then apply laplace again
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,v,tmp,tmp2)
#endif
      do k=1,nlev
         tmp(:,:)=rspheremv(:,:)*ptens(:,:,k,ie)
         ptens(:,:,k,ie)=laplace_sphere_wk(tmp,deriv,elem(ie),var_coef=.true.)
         tmp2(:,:)=rspheremv(:,:)*dptens(:,:,k,ie)
         dptens(:,:,k,ie)=laplace_sphere_wk(tmp2,deriv,elem(ie),var_coef=.true.)
         v(:,:,1)=rspheremv(:,:)*vtens(:,:,1,k,ie)
         v(:,:,2)=rspheremv(:,:)*vtens(:,:,2,k,ie)
         vtens(:,:,:,k,ie)=vlaplace_sphere_wk(v(:,:,:),deriv,elem(ie),&
              var_coef=.true.,nu_ratio=nu_ratio2)

      enddo
   enddo
#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine


end module 
