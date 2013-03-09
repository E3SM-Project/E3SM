#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module viscosity_mod
!
!  This module should be renamed "global_deriv_mod.F90"
! 
!  It is a collection of derivative operators that must be applied to the field 
!  over the sphere (as opposed to derivative operators that can be applied element 
!  by element)
!
!
use kinds, only : real_kind, iulog
use dimensions_mod, only : np, nlev,qsize
use hybrid_mod, only : hybrid_t
use element_mod, only : element_t
use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk, vorticity_sphere, derivinit, divergence_sphere
use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack, edgevunpackmin, &
    edgevunpackmax, initEdgeBuffer, FreeEdgeBuffer
use bndry_mod, only : bndry_exchangev

implicit none
public :: biharmonic_wk
#ifdef _PRIM
public :: biharmonic_wk_scalar
public :: biharmonic_wk_scalar_minmax
#endif
public :: compute_zeta_C0
public :: compute_div_C0
public :: compute_zeta_C0_2d
public :: compute_div_C0_2d
public :: test_ibyp

interface compute_zeta_C0_2d
    module procedure compute_zeta_C0_2d_sphere
    module procedure compute_zeta_C0_2d_contra
end interface

interface compute_div_C0_2d
    module procedure compute_div_C0_2d_sphere
    module procedure compute_div_C0_2d_contra
end interface


type (EdgeBuffer_t)          :: edge1

contains

#ifdef _PRIM
subroutine biharmonic_wk(elem,pstens,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete,nu_ratio)
#else
subroutine biharmonic_wk(elem,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete,nu_ratio)
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute weak biharmonic operator
!    input:  h,v (stored in elem()%, in lat-lon coordinates
!    output: ptens,vtens  overwritten with weak biharmonic of h,v (output in lat-lon coordinates)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(inout), target :: elem(:)
real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)  :: vtens
real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: ptens
type (EdgeBuffer_t)  , intent(inout) :: edge3
type (derivative_t)  , intent(in) :: deriv
integer :: nt,nets,nete
real (kind=real_kind) ::  nu_ratio
#ifdef _PRIM
real (kind=real_kind), dimension(np,np,nets:nete) :: pstens
#endif

! local
integer :: k,kptr,i,j,ie,ic
real (kind=real_kind), dimension(:,:), pointer :: rspheremv, viscosity
real (kind=real_kind), dimension(np,np) :: lap_ps
real (kind=real_kind), dimension(np,np,nlev) :: T
real (kind=real_kind), dimension(np,np,2) :: v

   do ie=nets,nete
      
      viscosity    => elem(ie)%variable_hyperviscosity

#ifdef _PRIM
      ! should filter lnps + PHI_s/RT?
      pstens(:,:,ie)=laplace_sphere_wk(elem(ie)%state%ps_v(:,:,nt),deriv,elem(ie),viscosity)
#endif
      
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k, j, i)
#endif
      do k=1,nlev
         do j=1,np
            do i=1,np
#ifdef _PRIM
               T(i,j,k)=elem(ie)%state%T(i,j,k,nt) 
#elif defined _PRIMDG
            T(i,j,k)=elem(ie)%state%p(i,j,k,nt) + elem(ie)%state%phis(i,j)
#else            
               ! filter surface height, not thickness
               T(i,j,k)=elem(ie)%state%p(i,j,k,nt) + elem(ie)%state%ps(i,j)
#endif
            enddo
         enddo
         ptens(:,:,k,ie)=laplace_sphere_wk(T(:,:,k),deriv,elem(ie),viscosity)
         vtens(:,:,:,k,ie)=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),viscosity,nu_ratio)
      enddo
      kptr=0
      call edgeVpack(edge3, ptens(1,1,1,ie),nlev,kptr,elem(ie)%desc)
      kptr=nlev
      call edgeVpack(edge3, vtens(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)

#ifdef _PRIM
      kptr=3*nlev
      call edgeVpack(edge3, pstens(1,1,ie),1,kptr,elem(ie)%desc)
#endif
   enddo
   
   call bndry_exchangeV(hybrid,edge3)
   
   do ie=nets,nete
      rspheremv     => elem(ie)%rspheremp(:,:)
      viscosity     => elem(ie)%variable_hyperviscosity
      
      kptr=0
      call edgeVunpack(edge3, ptens(1,1,1,ie), nlev, kptr, elem(ie)%desc)
      kptr=nlev
      call edgeVunpack(edge3, vtens(1,1,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
      
      ! apply inverse mass matrix, then apply laplace again
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k, j, i, v)
#endif
      do k=1,nlev
         do j=1,np
            do i=1,np
               T(i,j,k)=rspheremv(i,j)*ptens(i,j,k,ie)
               v(i,j,1)=rspheremv(i,j)*vtens(i,j,1,k,ie)
               v(i,j,2)=rspheremv(i,j)*vtens(i,j,2,k,ie)
            enddo
         enddo
         ptens(:,:,k,ie)=laplace_sphere_wk(T(:,:,k),deriv,elem(ie),viscosity)
         vtens(:,:,:,k,ie)=vlaplace_sphere_wk(v(:,:,:),deriv,elem(ie),viscosity,nu_ratio)
      enddo
         
#ifdef _PRIM
      kptr=3*nlev
      call edgeVunpack(edge3, pstens(1,1,ie), 1, kptr, elem(ie)%desc)
      ! apply inverse mass matrix, then apply laplace again
      lap_ps(:,:)=rspheremv(:,:)*pstens(:,:,ie)
      pstens(:,:,ie)=laplace_sphere_wk(lap_ps,deriv,elem(ie),viscosity)
#endif

   enddo
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine


#ifdef _PRIM
subroutine biharmonic_wk_dp3d(elem,dptens,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete,nu_ratio)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute weak biharmonic operator
!    input:  h,v (stored in elem()%, in lat-lon coordinates
!    output: ptens,vtens  overwritten with weak biharmonic of h,v (output in lat-lon coordinates)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(inout), target :: elem(:)
real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)  :: vtens
real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: ptens,dptens
type (EdgeBuffer_t)  , intent(inout) :: edge3
type (derivative_t)  , intent(in) :: deriv
integer :: nt,nets,nete
real (kind=real_kind) ::  nu_ratio

! local
integer :: k,kptr,ie
real (kind=real_kind), dimension(:,:), pointer :: rspheremv, viscosity
real (kind=real_kind), dimension(np,np) :: tmp
real (kind=real_kind), dimension(np,np,2) :: v

   do ie=nets,nete
      viscosity    => elem(ie)%variable_hyperviscosity

#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
      do k=1,nlev
         tmp=elem(ie)%state%T(:,:,k,nt) 
         ptens(:,:,k,ie)=laplace_sphere_wk(tmp,deriv,elem(ie),viscosity)
         tmp=elem(ie)%state%dp3d(:,:,k,nt) 
         dptens(:,:,k,ie)=laplace_sphere_wk(tmp,deriv,elem(ie),viscosity)
         vtens(:,:,:,k,ie)=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),viscosity,nu_ratio)
      enddo
      kptr=0
      call edgeVpack(edge3, ptens(1,1,1,ie),nlev,kptr,elem(ie)%desc)
      kptr=nlev
      call edgeVpack(edge3, vtens(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)
      kptr=3*nlev
      call edgeVpack(edge3, dptens(1,1,1,ie),nlev,kptr,elem(ie)%desc)

   enddo
   
   call bndry_exchangeV(hybrid,edge3)
   
   do ie=nets,nete
      rspheremv     => elem(ie)%rspheremp(:,:)
      viscosity     => elem(ie)%variable_hyperviscosity
      
      kptr=0
      call edgeVunpack(edge3, ptens(1,1,1,ie), nlev, kptr, elem(ie)%desc)
      kptr=nlev
      call edgeVunpack(edge3, vtens(1,1,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
      kptr=3*nlev
      call edgeVunpack(edge3, dptens(1,1,1,ie), nlev, kptr, elem(ie)%desc)
      
      ! apply inverse mass matrix, then apply laplace again
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k,  v)
#endif
      do k=1,nlev
         tmp(:,:)=rspheremv(:,:)*ptens(:,:,k,ie)
         ptens(:,:,k,ie)=laplace_sphere_wk(tmp,deriv,elem(ie),viscosity)
         tmp(:,:)=rspheremv(:,:)*dptens(:,:,k,ie)
         dptens(:,:,k,ie)=laplace_sphere_wk(tmp,deriv,elem(ie),viscosity)

         v(:,:,1)=rspheremv(:,:)*vtens(:,:,1,k,ie)
         v(:,:,2)=rspheremv(:,:)*vtens(:,:,2,k,ie)
         vtens(:,:,:,k,ie)=vlaplace_sphere_wk(v(:,:,:),deriv,elem(ie),viscosity,nu_ratio)

      enddo
   enddo
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine


subroutine biharmonic_wk_scalar(elem,qtens,deriv,edgeq,hybrid,nets,nete)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute weak biharmonic operator
!    input:  qtens = Q
!    output: qtens = weak biharmonic of Q
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(inout), target :: elem(:)
real (kind=real_kind), dimension(np,np,nlev,qsize,nets:nete) :: qtens
type (EdgeBuffer_t)  , intent(inout) :: edgeq
type (derivative_t)  , intent(in) :: deriv
integer :: nets,nete

! local
integer :: k,kptr,i,j,ie,ic,q
real (kind=real_kind), dimension(np,np) :: lap_p
real (kind=real_kind), dimension(:,:), pointer :: viscosity

   do ie=nets,nete
      viscosity    => elem(ie)%variable_hyperviscosity
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k, q, lap_p)
#endif
      do q=1,qsize      
         do k=1,nlev    !  Potential loop inversion (AAM)
            lap_p(:,:)=qtens(:,:,k,q,ie)
! Original use of qtens on left and right hand sides caused OpenMP errors (AAM)
           qtens(:,:,k,q,ie)=laplace_sphere_wk(lap_p,deriv,elem(ie),viscosity)
         enddo
      enddo
      call edgeVpack(edgeq, qtens(:,:,:,:,ie),qsize*nlev,0,elem(ie)%desc)
   enddo

   call bndry_exchangeV(hybrid,edgeq)
   
   do ie=nets,nete
      call edgeVunpack(edgeq, qtens(:,:,:,:,ie),qsize*nlev,0,elem(ie)%desc)

      ! apply inverse mass matrix, then apply laplace again
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k, q, lap_p)
#endif
      do q=1,qsize      
      do k=1,nlev    !  Potential loop inversion (AAM)
         lap_p(:,:)=elem(ie)%rspheremp(:,:)*qtens(:,:,k,q,ie)
         qtens(:,:,k,q,ie)=laplace_sphere_wk(lap_p,deriv,elem(ie),viscosity)
      enddo
      enddo
   enddo
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine

subroutine biharmonic_wk_scalar_minmax(elem,qtens,deriv,edgeq,hybrid,nets,nete,emin,emax)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute weak biharmonic operator
!    input:  qtens = Q
!    output: qtens = weak biharmonic of Q and Q element min/max
!
!    note: emin/emax must be initialized with Q element min/max.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(inout), target :: elem(:)
real (kind=real_kind), dimension(np,np,nlev,qsize,nets:nete) :: qtens
type (EdgeBuffer_t)  , intent(inout) :: edgeq
type (derivative_t)  , intent(in) :: deriv
real (kind=real_kind), intent(out), dimension(nlev,qsize,nets:nete) :: emin,emax
integer :: nets,nete

! local
integer :: k,kptr,i,j,ie,ic,q
real (kind=real_kind), dimension(np,np) :: lap_p
real (kind=real_kind) :: Qmin(np,np,nlev,qsize)
real (kind=real_kind) :: Qmax(np,np,nlev,qsize)
real (kind=real_kind), dimension(:,:), pointer :: viscosity 


   do ie=nets,nete
      viscosity    => elem(ie)%variable_hyperviscosity
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k, q, lap_p)
#endif
      do q=1,qsize      
      do k=1,nlev    !  Potential loop inversion (AAM)
         Qmin(:,:,k,q)=emin(k,q,ie)  ! need to set all values in element for
         Qmax(:,:,k,q)=emax(k,q,ie)  ! edgeVpack routine below
         lap_p(:,:) = qtens(:,:,k,q,ie)
! Original use of qtens on left and right hand sides caused OpenMP errors (AAM)
         qtens(:,:,k,q,ie)=laplace_sphere_wk(lap_p,deriv,elem(ie),viscosity)
      enddo
      enddo
      call edgeVpack(edgeq, qtens(:,:,:,:,ie),qsize*nlev,0,elem(ie)%desc)
      call edgeVpack(edgeq,Qmin,nlev*qsize,nlev*qsize,elem(ie)%desc)
      call edgeVpack(edgeq,Qmax,nlev*qsize,2*nlev*qsize,elem(ie)%desc)
   enddo
   
   call bndry_exchangeV(hybrid,edgeq)
   
   do ie=nets,nete
      do q=1,qsize      
      do k=1,nlev
         Qmin(:,:,k,q)=emin(k,q,ie)  ! restore element data.  we could avoid
         Qmax(:,:,k,q)=emax(k,q,ie)  ! this by adding a "ie" index to Qmin/max
      enddo
      enddo
! WARNING - edgeVunpackMin/Max take second argument as input/ouput
      call edgeVunpack(edgeq, qtens(:,:,:,:,ie),qsize*nlev,0,elem(ie)%desc)
      call edgeVunpackMin(edgeq, Qmin,qsize*nlev,qsize*nlev,elem(ie)%desc)
      call edgeVunpackMax(edgeq, Qmax,qsize*nlev,2*qsize*nlev,elem(ie)%desc)

      ! apply inverse mass matrix, then apply laplace again
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k, q, lap_p)
#endif
      do q=1,qsize      
      do k=1,nlev    !  Potential loop inversion (AAM)
         lap_p(:,:)=elem(ie)%rspheremp(:,:)*qtens(:,:,k,q,ie)
         qtens(:,:,k,q,ie)=laplace_sphere_wk(lap_p,deriv,elem(ie),viscosity)
         ! note: only need to consider the corners, since the data we packed was
         ! constant within each element
         emin(k,q,ie)=min(qmin(1,1,k,q),qmin(1,np,k,q),qmin(np,1,k,q),qmin(np,np,k,q))
         emin(k,q,ie)=max(emin(k,q,ie),0d0)
         emax(k,q,ie)=max(qmax(1,1,k,q),qmax(1,np,k,q),qmax(np,1,k,q),qmax(np,np,k,q))
      enddo
      enddo
   enddo
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine

#endif





subroutine make_C0_2d(zeta,elem,hybrid,nets,nete)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! apply DSS (aka assembly procedure) to zeta.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
real (kind=real_kind), dimension(np,np,nets:nete) :: zeta
integer :: nets,nete

! local
integer :: k,i,j,ie,ic,kptr


call initEdgeBuffer(edge1,1)

do ie=nets,nete
   zeta(:,:,ie)=zeta(:,:,ie)*elem(ie)%spheremp(:,:)
   kptr=0
   call edgeVpack(edge1, zeta(1,1,ie),1,kptr,elem(ie)%desc)
enddo
call bndry_exchangeV(hybrid,edge1)
do ie=nets,nete
   kptr=0
   call edgeVunpack(edge1, zeta(1,1,ie),1,kptr,elem(ie)%desc)
   zeta(:,:,ie)=zeta(:,:,ie)*elem(ie)%rspheremp(:,:)
enddo

call FreeEdgeBuffer(edge1) 
end subroutine


subroutine make_C0(zeta,elem,hybrid,nets,nete)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! apply DSS (aka assembly procedure) to zeta.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: zeta
integer :: nets,nete

! local
integer :: k,i,j,ie,ic,kptr


call initEdgeBuffer(edge1,nlev)

do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
   do k=1,nlev
      zeta(:,:,k,ie)=zeta(:,:,k,ie)*elem(ie)%spheremp(:,:)
   enddo
   kptr=0
   call edgeVpack(edge1, zeta(1,1,1,ie),nlev,kptr,elem(ie)%desc)
enddo
call bndry_exchangeV(hybrid,edge1)
do ie=nets,nete
   kptr=0
   call edgeVunpack(edge1, zeta(1,1,1,ie),nlev,kptr,elem(ie)%desc)
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
   do k=1,nlev
      zeta(:,:,k,ie)=zeta(:,:,k,ie)*elem(ie)%rspheremp(:,:)
   enddo
enddo
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif

call FreeEdgeBuffer(edge1) 
end subroutine


subroutine make_C0_vector(v,elem,hybrid,nets,nete)
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
real (kind=real_kind), dimension(np,np,2,nlev,nets:nete) :: v
integer :: nets,nete

! local
integer :: k,i,j,ie,ic,kptr
type (EdgeBuffer_t)          :: edge2


call initEdgeBuffer(edge2,2*nlev)

do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
   do k=1,nlev
      v(:,:,1,k,ie)=v(:,:,1,k,ie)*elem(ie)%spheremp(:,:)
      v(:,:,2,k,ie)=v(:,:,2,k,ie)*elem(ie)%spheremp(:,:)
   enddo
   kptr=0
   call edgeVpack(edge2, v(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)
enddo
call bndry_exchangeV(hybrid,edge2)
do ie=nets,nete
   kptr=0
   call edgeVunpack(edge2, v(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
   do k=1,nlev
      v(:,:,1,k,ie)=v(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)
      v(:,:,2,k,ie)=v(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
   enddo
enddo
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif

call FreeEdgeBuffer(edge2) 
end subroutine





subroutine compute_zeta_C0_2d_sphere(zeta,elem,hybrid,nets,nete,nt,k)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 vorticity.  That is, solve:  
!     < PHI, zeta > = <PHI, curl(elem%state%v >
!
!    input:  v (stored in elem()%, in lat-lon coordinates)
!    output: zeta(:,:,:,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
real (kind=real_kind), dimension(np,np,nets:nete) :: zeta
integer :: nt,nets,nete,k

! local
integer :: i,j,ie,ic
type (derivative_t)          :: deriv

call derivinit(deriv)

do ie=nets,nete
   !    zeta(:,:,k,ie)=elem(ie)%state%zeta(:,:,k)
   zeta(:,:,ie)=vorticity_sphere(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie))
enddo

call make_C0_2d(zeta,elem,hybrid,nets,nete)

end subroutine

subroutine compute_zeta_C0_2d_contra(zeta,elem,hybrid,nets,nete,nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 vorticity.  That is, solve:  
!     < PHI, zeta > = <PHI, curl(elem%state%v >
!
!    input:  v (stored in elem()%, in contra-variant coordinates)
!    output: zeta(:,:,:,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: zeta
real (kind=real_kind), dimension(np,np,2) :: ulatlon
real (kind=real_kind), dimension(np,np) :: v1,v2
integer :: nt,nets,nete

! local
integer :: k,ie
type (derivative_t)          :: deriv

call derivinit(deriv)

do k=1,nlev
do ie=nets,nete
    v1 = elem(ie)%state%v(:,:,1,k,nt)
    v2 = elem(ie)%state%v(:,:,2,k,nt)
    ulatlon(:,:,1) = elem(ie)%D(1,1,:,:)*v1 + elem(ie)%D(1,2,:,:)*v2
    ulatlon(:,:,2) = elem(ie)%D(2,1,:,:)*v1 + elem(ie)%D(2,2,:,:)*v2
   !    zeta(:,:,k,ie)=elem(ie)%state%zeta(:,:,k)
   zeta(:,:,k,ie)=vorticity_sphere(ulatlon,deriv,elem(ie))
enddo
enddo

call make_C0(zeta,elem,hybrid,nets,nete)

end subroutine


subroutine compute_div_C0_2d_sphere(zeta,elem,hybrid,nets,nete,nt,k)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 divergence. That is, solve:  
!     < PHI, zeta > = <PHI, div(elem%state%v >
!
!    input:  v (stored in elem()%, in lat-lon coordinates)
!    output: zeta(:,:,:,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
real (kind=real_kind), dimension(np,np,nets:nete) :: zeta
integer :: nt,nets,nete,k

! local
integer :: i,j,ie,ic
type (derivative_t)          :: deriv

call derivinit(deriv)

do ie=nets,nete
   !    zeta(:,:,k,ie)=elem(ie)%state%zeta(:,:,k)
   zeta(:,:,ie)=divergence_sphere(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie))
enddo

call make_C0_2d(zeta,elem,hybrid,nets,nete)

end subroutine

subroutine compute_div_C0_2d_contra(zeta,elem,hybrid,nets,nete,nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 divergence. That is, solve:  
!     < PHI, zeta > = <PHI, div(elem%state%v >
!
!    input:  v (stored in elem()%, in contra-variant coordinates)
!    output: zeta(:,:,:,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: zeta
real (kind=real_kind), dimension(np,np,2) :: ulatlon
real (kind=real_kind), dimension(np,np) :: v1,v2
integer :: nt,nets,nete

! local
integer :: k,ie
type (derivative_t)          :: deriv

call derivinit(deriv)

do k=1,nlev
do ie=nets,nete
    v1 = elem(ie)%state%v(:,:,1,k,nt)
    v2 = elem(ie)%state%v(:,:,2,k,nt)
    ulatlon(:,:,1) = elem(ie)%D(1,1,:,:)*v1 + elem(ie)%D(1,2,:,:)*v2
    ulatlon(:,:,2) = elem(ie)%D(2,1,:,:)*v1 + elem(ie)%D(2,2,:,:)*v2
   !    zeta(:,:,k,ie)=elem(ie)%state%zeta(:,:,k)
   zeta(:,:,k,ie)=divergence_sphere(ulatlon,deriv,elem(ie))
enddo
enddo

call make_C0(zeta,elem,hybrid,nets,nete)

end subroutine

subroutine compute_zeta_C0(zeta,elem,hybrid,nets,nete,nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 vorticity.  That is, solve:  
!     < PHI, zeta > = <PHI, curl(elem%state%v >
!
!    input:  v (stored in elem()%, in lat-lon coordinates)
!    output: zeta(:,:,:,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: zeta
integer :: nt,nets,nete

! local
integer :: k,i,j,ie,ic
type (derivative_t)          :: deriv

call derivinit(deriv)

do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
do k=1,nlev
   !    zeta(:,:,k,ie)=elem(ie)%state%zeta(:,:,k)
   zeta(:,:,k,ie)=vorticity_sphere(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie))
enddo
enddo

call make_C0(zeta,elem,hybrid,nets,nete)

end subroutine


subroutine compute_div_C0(zeta,elem,hybrid,nets,nete,nt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute C0 divergence. That is, solve:  
!     < PHI, zeta > = <PHI, div(elem%state%v >
!
!    input:  v (stored in elem()%, in lat-lon coordinates)
!    output: zeta(:,:,:,:)   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(in), target :: elem(:)
real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: zeta
integer :: nt,nets,nete

! local
integer :: k,i,j,ie,ic
type (derivative_t)          :: deriv

call derivinit(deriv)

do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
do k=1,nlev
   !    zeta(:,:,k,ie)=elem(ie)%state%zeta(:,:,k)
   zeta(:,:,k,ie)=divergence_sphere(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie))
enddo
enddo

call make_C0(zeta,elem,hybrid,nets,nete)

end subroutine








#ifdef _PRIM
subroutine neighbor_minmax(elem,hybrid,edgeMinMax,nets,nete,min_neigh,max_neigh)
!
! compute Q min&max over the element and all its neighbors
!
!
integer :: nets,nete
type (element_t)     , intent(in) :: elem(:)
type (hybrid_t)      , intent(in) :: hybrid
type (EdgeBuffer_t)  , intent(inout) :: edgeMinMax
real (kind=real_kind) :: min_neigh(nlev,qsize,nets:nete)
real (kind=real_kind) :: max_neigh(nlev,qsize,nets:nete)

! local
integer :: ie,k,q
real (kind=real_kind) :: Qmin(np,np,nlev,qsize)
real (kind=real_kind) :: Qmax(np,np,nlev,qsize)


    ! compute Qmin, Qmax
    do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k, q)
#endif
       do q=1,qsize	
          do k=1,nlev
             Qmin(:,:,k,q)=min_neigh(k,q,ie)
             Qmax(:,:,k,q)=max_neigh(k,q,ie)
          enddo
       end do
       call edgeVpack(edgeMinMax,Qmin,nlev*qsize,0,elem(ie)%desc)
       call edgeVpack(edgeMinMax,Qmax,nlev*qsize,nlev*qsize,elem(ie)%desc)
    enddo

    call bndry_exchangeV(hybrid,edgeMinMax)
       
    do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k, q)
#endif
       do q=1,qsize	
          do k=1,nlev         
             Qmin(:,:,k,q)=min_neigh(k,q,ie) ! restore element data.  we could avoid
             Qmax(:,:,k,q)=max_neigh(k,q,ie) ! this by adding a "ie" index to Qmin/max
          enddo
       end do
! WARNING - edgeVunpackMin/Max take second argument as input/ouput
       call edgeVunpackMin(edgeMinMax,Qmin,nlev*qsize,0,elem(ie)%desc)
       call edgeVunpackMax(edgeMinMax,Qmax,nlev*qsize,nlev*qsize,elem(ie)%desc)
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k, q)
#endif
       do q=1,qsize
          do k=1,nlev
             ! note: only need to consider the corners, since the data we packed was
             ! constant within each element
             min_neigh(k,q,ie)=min(qmin(1,1,k,q),qmin(1,np,k,q),qmin(np,1,k,q),qmin(np,np,k,q))
             min_neigh(k,q,ie)=max(min_neigh(k,q,ie),0d0)
             max_neigh(k,q,ie)=max(qmax(1,1,k,q),qmax(1,np,k,q),qmax(np,1,k,q),qmax(np,np,k,q))
          enddo
       end do
    end do
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif

end subroutine


#else


subroutine neighbor_minmax(elem,hybrid,edgeMinMax,nets,nete,nt,min_neigh,max_neigh,min_var,max_var,kmass)
!
! compute Q min&max over the element and all its neighbors
!
!
integer :: nets,nete,nt
type (hybrid_t)      , intent(in) :: hybrid
type (element_t)     , intent(inout) :: elem(:)
type (EdgeBuffer_t)  , intent(in) :: edgeMinMax
real (kind=real_kind) :: min_neigh(nlev,nets:nete)
real (kind=real_kind) :: max_neigh(nlev,nets:nete)
real (kind=real_kind),optional :: min_var(nlev,nets:nete)
real (kind=real_kind),optional :: max_var(nlev,nets:nete)
real (kind=real_kind) :: Qmin(np,np,nlev)
real (kind=real_kind) :: Qmax(np,np,nlev)
real (kind=real_kind) :: Qvar(np,np,nlev)
type (EdgeBuffer_t)          :: edgebuf
integer, optional :: kmass

! local
integer :: ie,k,q

  if(present(kmass))then
!the check if kmass is a valid number is done in sweq_mod
    do k=1,nlev
      if(k.ne.kmass)then
         do ie=nets,nete
            elem(ie)%state%p(:,:,k,nt)=elem(ie)%state%p(:,:,k,nt)/&
            elem(ie)%state%p(:,:,kmass,nt)
         enddo
      endif
    enddo
  endif



    ! create edge buffer for 3 fields
    call initEdgeBuffer(edgebuf,3*nlev)


    ! compute p min, max
    do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
       do k=1,nlev
          Qmin(:,:,k)=minval(elem(ie)%state%p(:,:,k,nt))
          Qmax(:,:,k)=maxval(elem(ie)%state%p(:,:,k,nt))
          ! max - min - crude approximation to TV within the element:
          Qvar(:,:,k)=Qmax(1,1,k)-Qmin(1,1,k)
       enddo
       call edgeVpack(edgebuf,Qmax,nlev,0,elem(ie)%desc)
       call edgeVpack(edgebuf,Qmin,nlev,nlev,elem(ie)%desc)
       call edgeVpack(edgebuf,Qvar,nlev,2*nlev,elem(ie)%desc)
    enddo
    
    call bndry_exchangeV(hybrid,edgebuf)
       
    do ie=nets,nete
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
       do k=1,nlev
          Qmin(:,:,k)=minval(elem(ie)%state%p(:,:,k,nt))
          Qmax(:,:,k)=maxval(elem(ie)%state%p(:,:,k,nt))
       enddo

       ! now unpack the min
       if (present(min_var)) then
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
          do k=1,nlev
             Qvar(:,:,k)=Qmax(1,1,k)-Qmin(1,1,k)
          enddo
! WARNING - edgeVunpackMin/Max take second argument as input/ouput
          call edgeVunpackMin(edgebuf,Qvar,nlev,2*nlev,elem(ie)%desc)
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
          do k=1,nlev
             min_var(k,ie)=minval(Qvar(:,:,k))
          enddo
       endif

       ! now unpack the max
       if (present(max_var)) then
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
          do k=1,nlev
             Qvar(:,:,k)=Qmax(1,1,k)-Qmin(1,1,k)
          enddo
! WARNING - edgeVunpackMin/Max take second argument as input/ouput
          call edgeVunpackMax(edgebuf,Qvar,nlev,2*nlev,elem(ie)%desc)
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
          do k=1,nlev
             max_var(k,ie)=maxval(Qvar(:,:,k))
          enddo
       endif


! WARNING - edgeVunpackMin/Max take second argument as input/ouput
       call edgeVunpackMax(edgebuf,Qmax,nlev,0,elem(ie)%desc)
       call edgeVunpackMin(edgebuf,Qmin,nlev,nlev,elem(ie)%desc)
#if (defined ELEMENT_OPENMP)
!$omp parallel do private(k)
#endif
       do k=1,nlev
          max_neigh(k,ie)=maxval(Qmax(:,:,k))
          min_neigh(k,ie)=minval(Qmin(:,:,k))
       enddo
       
    end do

    call FreeEdgeBuffer(edgebuf) 
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif

  if(present(kmass))then
    do k=1,nlev
       if(k.ne.kmass)then
          do ie=nets,nete
             elem(ie)%state%p(:,:,k,nt)=elem(ie)%state%p(:,:,k,nt)*&
             elem(ie)%state%p(:,:,kmass,nt)
          enddo
       endif
    enddo
  endif
end subroutine
#endif





  subroutine test_ibyp(elem, hybrid,  nets,   nete)
!
! Note: vector test functions should be co-variant since u is contra-variant
!  PHIvec = PHIcov  (test function)
!  PHIcon = DtD PHIcov
!
! weak grad:
!  < PHIcov du/dt > = < PHIcon grad(p) >    (output of grad is covariant)
!  < PHIcov du/dt > = -< div(PHIcon) p >    (input of div is contra)
!  verify:
!    gradient_sphere_wk(p) = - <div(PHIcon) p >
!    gradient_sphere_wk(p) + MASS*grad(p) = b.c. (b.c. are covariant)
!
! weak div:
!   < PHI div(u) > = -< grad(PHI) dot u >     u=contra, output of grad is covariant
! verify:
!   divergence_sphere_wk(u) = -<grad(PHI) dot u>
!   divergence_sphere_wk(u) + MASS*div(u) = b.c.  (b.c. are scalars)
!
! weak curl:
!  < PHIcov du/dt > = < PHIcov curl( a ) >    (output of curl is contra)
!  < PHIcov du/dt > = < vor(PHIcov) a >       (input to vor is covariant)
! verify:
!    curl_sphere_wk(a) = < vor(PHIcov) a >
!    curl_sphere_wk(a) - MASS*curl(a) = b.c. (b.c. are contra)
!
    ! ---------------------
    use kinds, only : real_kind
    ! ---------------------
    use physical_constants, only : rearth 
    ! ---------------------
    use dimensions_mod, only : np, nlev
    ! ---------------------
    use element_mod, only : element_t
    ! ---------------------
    use hybrid_mod, only : hybrid_t
    ! ---------------------
    use derivative_mod, only : derivative_t, gradient_sphere, divergence_sphere,vorticity_sphere,&
                               divergence_sphere_wk, curl_sphere
    use global_norms_mod

    implicit none

    type (element_t)     , intent(inout), target :: elem(:)

    type (hybrid_t)      , intent(in) :: hybrid

    integer              , intent(in) :: nets
    integer              , intent(in) :: nete

#if 0
#undef CURLGRAD_TEST
#define IBYP_TEST
    ! =================
    ! Local
    ! =================
    ! pointer ...
    real (kind=real_kind), dimension(:,:), pointer :: rspheremv,spheremv

    ! Thread private working set ...

    real (kind=real_kind), dimension(np,np,nets:nete) :: ptens
    real (kind=real_kind), dimension(np,np,nets:nete) :: ptens2
    real (kind=real_kind), dimension(np,np,nets:nete) :: ptens3

    real (kind=real_kind), dimension(np,np,2,nets:nete)    :: pv      ! p*v lat-lon
    real (kind=real_kind), dimension(np,np,nets:nete)            :: E          ! kinetic energy term
    real (kind=real_kind), dimension(np,np,nets:nete)  :: divbig
    real (kind=real_kind), dimension(np,np,nets:nete)  :: wdivbig
    real (kind=real_kind), dimension(np,np,2,nets:nete)    :: gradbig

    real (kind=real_kind), dimension(np,np,2)      :: ulatlon
    real (kind=real_kind), dimension(np,np,2)      :: grade
    real (kind=real_kind), dimension(np,np,2)      :: grade2
    real (kind=real_kind), dimension(np,np)      :: vor
    real (kind=real_kind), dimension(np,np)      :: div  
    real (kind=real_kind), dimension(np,np)      :: wdiv  

    real (kind=real_kind) ::  v1,v2,v3

    real*8                :: st,et, time_adv
    integer    :: i,j,k,ie,iie,ii,jj
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep

    type (derivative_t)          :: deriv
    call derivinit(deriv)



    ! ===================================
    ! construct test function
    ! ===================================

    do iie=nets,nete
       do jj=1,np
          do ii=1,np
             ! test for cardinal function at iie,jj,ii

    write(iulog,'(a,3i4,2e20.10)') 'carinal function:  ',iie,ii,jj

    ! construct two global cardinal functions  pv and E    
    do ie=nets,nete
       do j=1,np
          do i=1,np


             E(i,j,ie)=0
#ifdef CURLGRAD_TEST
             if (ie==iie .and. i==ii .and. j==jj) E(i,j,ie)=1
#else
             if (ie==1 .and. i==1 .and. j==1) E(i,j,ie)=1
             if (ie==1 .and. i==1 .and. j==2) E(i,j,ie)=1
             if (ie==1 .and. i==2 .and. j==1) E(i,j,ie)=1
             if (ie==1 .and. i==3 .and. j==3) E(i,j,ie)=1
#endif
             
             ! delta function in contra coordinates
             v1     = 0
             v2     = 0
             if (ie==iie .and. i==ii .and. j==jj) then
                !v1=rearth
                v2=rearth
             endif
             
             ulatlon(i,j,1)=elem(ie)%D(1,1,i,j)*v1 + elem(ie)%D(1,2,i,j)*v2   ! contra->latlon
             ulatlon(i,j,2)=elem(ie)%D(2,1,i,j)*v1 + elem(ie)%D(2,2,i,j)*v2   ! contra->latlon
             pv(i,j,1,ie) = ulatlon(i,j,1)
             pv(i,j,2,ie) = ulatlon(i,j,2)
             
          end do
       end do
    enddo
    call make_C0(E,elem,hybrid,nets,nete)
    call make_C0_vector(pv,elem,hybrid,nets,nete) 


#ifdef CURLGRAD_TEST
    ! check curl(grad(E)) 
    do ie=nets,nete
       if ( maxval(abs(E(:,:,ie))) > 0 ) then
       !write(iulog,'(a,i4,2e20.10)') 'maxval: E =',ie,maxval(E(:,:,ie))
       grade=curl_sphere(E(:,:,ie),deriv,elem(ie))
       div=divergence_sphere(grade,deriv,elem(ie))
       vor=vorticity_sphere(grade,deriv,elem(ie))
       if (maxval(abs(div))*rearth**2 > .2e-11) then
          write(iulog,'(a,i4,2e20.10)') 'maxval: div(curl),  vor(curl)=',ie,maxval(abs(div))*rearth**2,maxval(abs(vor))*rearth**2
       endif

       grade=gradient_sphere(E(:,:,ie),deriv,elem(ie)%Dinv)
       vor=vorticity_sphere(grade,deriv,elem(ie))
       div=divergence_sphere(grade,deriv,elem(ie))
       if (maxval(abs(vor))*rearth**2 > .2e-11) then
          write(iulog,'(a,i4,2e20.10)') 'maxval: curl(grad), div(grad)=',ie,maxval(abs(vor))*rearth**2,maxval(abs(div))*rearth**2
       endif
       endif
    enddo

    ! check div(curl(E)) with DSS 
    do ie=nets,nete
       gradbig(:,:,:,ie)=curl_sphere(E(:,:,ie),deriv,elem(ie))
    enddo
    call make_C0_vector(gradbig,elem,hybrid,nets,nete)
    do ie=nets,nete
       divbig(:,:,ie)=divergence_sphere(gradbig(:,:,:,ie),deriv,elem(ie))
    enddo
    call make_C0(divbig,elem,hybrid,nets,nete)
    do ie=nets,nete
       if (maxval(abs(divbig(:,:,ie)))*rearth**2 > .8e-12) then
          write(iulog,'(a,i4,2e20.10)') 'maxval: [div([curl])]=',ie,maxval(abs(divbig(:,:,ie)))*rearth**2
       endif
    enddo


    ! check curl(grad(E)) with DSS 
    do ie=nets,nete
       gradbig(:,:,:,ie)=gradient_sphere(E(:,:,ie),deriv,elem(ie)%Dinv)
    enddo
    call make_C0_vector(gradbig,elem,hybrid,nets,nete)
    do ie=nets,nete
       divbig(:,:,ie)=vorticity_sphere(gradbig(:,:,:,ie),deriv,elem(ie))
    enddo
    call make_C0(divbig,elem,hybrid,nets,nete)
    
    do ie=nets,nete
       if (maxval(abs(divbig(:,:,ie)))*rearth**2 > .8e-12) then
          write(iulog,'(a,i4,2e20.10)') 'maxval: [curl([gradl])]=',ie,maxval(abs(divbig(:,:,ie)))*rearth**2
       endif
    enddo
#endif


#ifdef IBYP_TEST
    ! compare <grad(E) dot pv> and <E div(pv)>  < E weak_div(pv) >
    v2=0
    do ie=nets,nete
       spheremv     => elem(ie)%spheremp(:,:)

       div = divergence_sphere(pv(1,1,1,ie),deriv,elem(ie))      ! latlon vector -> scalar 
       grade = gradient_sphere(E(1,1,ie),deriv,elem(ie)%Dinv)
       wdiv = divergence_sphere_wk(pv(1,1,1,ie),deriv,elem(ie)) 



       do j=1,np
          do i=1,np
!             write(iulog,'(3i3,3e22.14)') ie,i,j,pv(i,j,1,ie),pv(i,j,2,ie),div(i,j)

             ! (grad(E) dot pv )
             ptens3(i,j,ie) = grade(i,j,1)*pv(i,j,1,ie) + grade(i,j,2)*pv(i,j,2,ie)
             v2=v2+wdiv(i,j)*E(i,j,ie)
          end do
       end do
       ptens(:,:,ie)=div(:,:)*E(:,:,ie)   ! < E div(pv) >
       ptens2(:,:,ie)=wdiv(:,:)*E(:,:,ie)/spheremv(:,:)   ! < wdiv E >
       ! ===================================================
       ! Pack cube edges of tendencies, rotate velocities
       ! ===================================================
       divbig(:,:,ie)=div(:,:)*spheremv(:,:)
       wdivbig(:,:,ie)=wdiv(:,:)
    end do
    call make_C0(divbig,elem,hybrid,nets,nete)
    call make_C0(wdivbig,elem,hybrid,nets,nete)

    
!!$    v1=global_integral(elem,ptens,hybrid,np,nets,nete)
!!$    v3=global_integral(elem,ptens3,hybrid,np,nets,nete)
!!$    print *,'< E div(pv) >   =',v1
!!$    print *,'< E div_wk(pv) >=',v2/(4*4*atan(1d0))
!!$    v2=global_integral(elem,ptens2,hybrid,np,nets,nete)
!!$    print *,'< E div_wk(pv) >=',v2
!!$    print *,'-<grad(E),pv >  =',-v3
!!$!    print *,'sum1-sum2/max',(v1-v2)/max(abs(v1),abs(v2))
!!$


    do ie=nets,nete
       div(:,:)=divbig(:,:,ie)
       wdiv(:,:)=wdivbig(:,:,ie)
       ! ===========================================================
       ! Compute velocity and pressure tendencies for all levels
       ! ===========================================================
          do j=1,np
             do i=1,np
                ! < div(pv) >   vs < div_wk(pv) >
                if ( abs(div(i,j)-wdiv(i,j)) > .15e-17) then
!                   write(iulog,'(3i3,4e22.14)') ie,i,j,div(i,j),wdiv(i,j),div(i,j)-wdiv(i,j),E(i,j,ie)
                endif
                if ( E(i,j,ie)/=0 ) then
!                   write(iulog,'(3i3,4e22.14)') ie,i,j,div(i,j),wdiv(i,j),div(i,j)-wdiv(i,j),E(i,j,ie)
                endif
                
             end do
          end do


    end do
#ifdef DEBUGOMP
#if (! defined ELEMENT_OPENMP)
!$OMP BARRIER
#endif
#endif
    write(iulog,'(a,3i4,2e20.10)') 'max diff div-wdiv: ',iie,ii,jj,maxval(abs(divbig(:,:,:)-wdivbig(:,:,:))),maxval(divbig(:,:,:))
#endif


    enddo
    enddo
    enddo
    stop
#endif
  end subroutine test_ibyp









  subroutine check_edge_flux(elem,deriv,nets,nete)
!
!  check two identities:
!  1. div and weak div are adjoints:
!     integral[  p div(u) ] + integral[ grad(p) dot u ] = boundary_integral[ u p]
!  1. grad and weak grad are adjoints:
!     integral[  p div(u) ] + integral[ grad(p) dot u ] = boundary_integral[ u p]
!
!
  use dimensions_mod, only : np, np, nlev
  use element_mod, only    : element_t
  use derivative_mod, only  : derivative_t, divergence_sphere, divergence_sphere_wk, &
                             element_boundary_integral, gradient_sphere, gradient_sphere_wk
  use physical_constants, only : rrearth

  implicit none
  
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv
  integer :: nets,nete
  ! local 
  real (kind=real_kind), dimension(np,np,2) :: ucontra,ulatlon,gradp,gradp_wk,ucov
  real (kind=real_kind), dimension(np,np) :: phidivu,ugradphi,rhs,lhs,p
  real (kind=real_kind), dimension(np,np) :: rhs2,lhs2
  integer :: i,j,ie


  print *,'integration by parts identity: check div/weak div:'
  ! test integration by parts identity for each Cardinal function PHI:
  ! div(u)*spheremp - div_wk(u) = boundary integral phi u dot n
  do ie=nets,nete
     call random_number(ucontra)
  ! contra->latlon
     ulatlon(:,:,1)=(elem(ie)%D(1,1,:,:)*ucontra(:,:,1) + elem(ie)%D(1,2,:,:)*ucontra(:,:,2))
     ulatlon(:,:,2)=(elem(ie)%D(2,1,:,:)*ucontra(:,:,1) + elem(ie)%D(2,2,:,:)*ucontra(:,:,2))
     phidivu = elem(ie)%spheremp(:,:)*divergence_sphere(ulatlon,deriv,elem(ie))
     ugradphi = divergence_sphere_wk(ulatlon,deriv,elem(ie))
     lhs = phidivu - ugradphi
     
     rhs = element_boundary_integral(ulatlon,deriv,elem(ie))
     
     
     do j=1,np
        do i=1,np
           if ( abs(lhs(i,j)-rhs(i,j)) .gt. 1d-20) then
              write(*,'(a)') 'ERROR: div/div_wk integration by parts failure!'
              write(*,'(a,2i3,a,3e12.5)') 'for test function (i,j)=',i,j,' lhs,rhs=',lhs(i,j),rhs(i,j),lhs(i,j)-rhs(i,j)
           endif
        enddo
     enddo
  enddo

  print *,'integration by parts identity: check grad/weak grad:'
  ! PHIVEC = contra cardinal function 
  !          check each contra component seperately

  do ie=nets,nete
     call random_number(p)
     
     ! grad(p)  (lat/lon vector)
     gradp = gradient_sphere(p,deriv,elem(ie)%Dinv)
     gradp(:,:,1)=gradp(:,:,1)*elem(ie)%spheremp(:,:)  
     gradp(:,:,2)=gradp(:,:,2)*elem(ie)%spheremp(:,:)
     gradp_wk = gradient_sphere_wk(p,deriv,elem(ie))
     
     ucontra(:,:,1)=p
     ucontra(:,:,2)=0
     ! contra->latlon
     ulatlon(:,:,1)=(elem(ie)%D(1,1,:,:)*ucontra(:,:,1) + elem(ie)%D(1,2,:,:)*ucontra(:,:,2))
     ulatlon(:,:,2)=(elem(ie)%D(2,1,:,:)*ucontra(:,:,1) + elem(ie)%D(2,2,:,:)*ucontra(:,:,2))

     rhs = element_boundary_integral(ulatlon,deriv,elem(ie))
     lhs = gradp(:,:,1)-gradp_wk(:,:,1)

     ucontra(:,:,1)=0
     ucontra(:,:,2)=p
     ! contra->latlon
     ulatlon(:,:,1)=(elem(ie)%D(1,1,:,:)*ucontra(:,:,1) + elem(ie)%D(1,2,:,:)*ucontra(:,:,2))
     ulatlon(:,:,2)=(elem(ie)%D(2,1,:,:)*ucontra(:,:,1) + elem(ie)%D(2,2,:,:)*ucontra(:,:,2))
     rhs2 = element_boundary_integral(ulatlon,deriv,elem(ie))
     lhs2 = gradp(:,:,2)-gradp_wk(:,:,2)  


     ! gradient_sphere() and gradient_sphere_wk() compute covariant vectors
     ! and then convert to latlon.  Thus to get the identity to work, the
     ! same transformation has to be applied to the vector of boundary integrals:
     ! tread (rhs,rhs2) as covariant vector and convert to latlon:
     gradp(:,:,1)=rhs
     gradp(:,:,2)=rhs2
     rhs(:,:)=elem(ie)%Dinv(1,1,:,:)*gradp(:,:,1) + elem(ie)%Dinv(2,1,:,:)*gradp(:,:,2)
     rhs2(:,:)=elem(ie)%Dinv(1,2,:,:)*gradp(:,:,1) + elem(ie)%Dinv(2,2,:,:)*gradp(:,:,2)


     do j=1,np
        do i=1,np
           if ( abs(lhs(i,j)-rhs(i,j)) .gt. 1d-20) then
              write(*,'(a)') 'ERROR: grad/grad_wk (1) integration by parts failure!'
              write(*,'(a,2i3,a,4e13.5)') 'for test function (i,j)=',i,j,&
                   ' lhs,rhs=',lhs(i,j),rhs(i,j),lhs(i,j)-rhs(i,j),lhs(i,j)/rhs(i,j)
           endif
        enddo
     enddo


     do j=1,np
        do i=1,np
           if ( abs(lhs2(i,j)-rhs2(i,j)) .gt. 1d-20) then
              write(*,'(a)') 'ERROR: grad/grad_wk (2) integration by parts failure!'
              write(*,'(a,2i3,a,3e12.5)') 'for test function (i,j)=',i,j,' lhs,rhs=',lhs(i,j),rhs(i,j),lhs(i,j)-rhs(i,j)
           endif
        enddo
     enddo

  enddo
  print *,'done. integration by parts identity check:'
  stop
  end subroutine




end module
