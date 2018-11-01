#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


module sl_advection
!
!  classic semi-lagrange advection
!  with optional global iteration for pseudo-local mass conservation
!
!  Author:  James Overfelt   3/2015
!
  use kinds, only              : real_kind
  use dimensions_mod, only     : nlev, nlevp, np, qsize
  use physical_constants, only : rgas, Rwater_vapor, kappa, g, rearth, rrearth, cp
  use derivative_mod, only     : derivative_t, gradient_sphere, divergence_sphere
  use element_mod, only        : element_t
  use hybvcoord_mod, only      : hvcoord_t
  use time_mod, only           : TimeLevel_t, TimeLevel_Qdp
  use control_mod, only        : integration, test_case, hypervis_order, use_semi_lagrange_transport
  use edge_mod, only           : edgevpack_nlyr, edgevunpack_nlyr, edge_g, &
                                 initghostbuffer3D, ghostVpack_unoriented, ghostVunpack_unoriented
  use edgetype_mod, only       : EdgeDescriptor_t, EdgeBuffer_t, ghostbuffer3D_t
  use hybrid_mod, only         : hybrid_t
  use bndry_mod, only          : bndry_exchangev
  use viscosity_mod, only      : biharmonic_wk_scalar, neighbor_minmax, &
                                 neighbor_minmax_start, neighbor_minmax_finish
  use perf_mod, only           : t_startf, t_stopf, t_barrierf ! _EXTERNAL
  use parallel_mod, only   : abortmp, parallel_t

  implicit none

  private
  save

  type (ghostBuffer3D_t)   :: ghostbuf_tr

  public :: Prim_Advec_Tracers_remap_ALE, sl_init1

contains


!=================================================================================================!
  subroutine sl_init1(par, elem)
    use interpolate_mod,        only : interpolate_tracers_init
    type(parallel_t) :: par
    type (element_t) :: elem(:)

    if  (use_semi_lagrange_transport) then
       call initghostbuffer3D(ghostbuf_tr,nlev*qsize,np)
       call interpolate_tracers_init()
    endif

  end subroutine 


subroutine  Prim_Advec_Tracers_remap_ALE( elem , deriv , hybrid , dt , tl , nets , nete )
  use coordinate_systems_mod, only : cartesian3D_t, cartesian2D_t
  use dimensions_mod,         only : max_neigh_edges
  use bndry_mod,              only : ghost_exchangevfull
  use interpolate_mod,        only : interpolate_tracers, minmax_tracers
  use control_mod   ,         only : qsplit
  use global_norms_mod,       only : wrap_repro_sum
  use parallel_mod,           only: global_shared_buf, global_shared_sum, syncmp
  use control_mod,            only : use_semi_lagrange_transport_local_conservation



  implicit none
  type (element_t)     , intent(inout) :: elem(:)
  type (derivative_t)  , intent(in   ) :: deriv
  type (hybrid_t)      , intent(in   ) :: hybrid
  real(kind=real_kind) , intent(in   ) :: dt
  type (TimeLevel_t)   , intent(in   ) :: tl
  integer              , intent(in   ) :: nets
  integer              , intent(in   ) :: nete


  type(cartesian3D_t)                           :: dep_points  (np,np)
  integer                                       :: elem_indexes(np,np)
  type(cartesian2D_t)                           :: para_coords (np,np)
  real(kind=real_kind)                          :: Que         (np,np,nlev,qsize,nets:nete)
  real(kind=real_kind)                          :: Que_t       (np,np,nlev,qsize,nets:nete)
  real(kind=real_kind)                          :: minq        (np,np,nlev,qsize,nets:nete)
  real(kind=real_kind)                          :: maxq        (np,np,nlev,qsize,nets:nete)
  real(kind=real_kind)                          :: f                      (qsize)
  real(kind=real_kind)                          :: g                      (qsize)
  real(kind=real_kind)                          :: mass              (nlev,qsize)
  real(kind=real_kind)                          :: elem_mass         (nlev,qsize,nets:nete)
  real(kind=real_kind)                          :: rho         (np,np,nlev,      nets:nete)

  real(kind=real_kind)                          :: neigh_q     (np,np,qsize,max_neigh_edges+1)
  real(kind=real_kind)                          :: u           (np,np,qsize)

  integer                                       :: i,j,k,l,n,q,ie,kptr, n0_qdp, np1_qdp
  integer                                       :: num_neighbors

  call t_barrierf('Prim_Advec_Tracers_remap_ALE', hybrid%par%comm)
  call t_startf('Prim_Advec_Tracers_remap_ALE')

  call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp)

  ! compute displacements for departure grid
  ! store in elem%derived%vstar
  call ALE_RKdss (elem, nets, nete, hybrid, deriv, dt, tl)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  run ghost exchange to get global ID of all neighbors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  call t_startf('Prim_Advec_Tracers_remap_ALE_ghost_exchange')
  do ie=nets,nete
     kptr=0
     do k=1,nlev
     do q=1,qsize
        ! note: pack so that tracers per level are contiguous so we can unpack into
        ! array neigh_q()
        elem(ie)%state%Q(:,:,k,q) = elem(ie)%state%Qdp(:,:,k,q,n0_qdp) / elem(ie)%derived%dp(:,:,k)
        call ghostVpack_unoriented(ghostbuf_tr, elem(ie)%state%Q(:,:,k,q),np,1,kptr,elem(ie)%desc)
        kptr=kptr+1
     enddo
     enddo
  end do

  call t_startf('pat_remap_ale_gexchV')
  call ghost_exchangeVfull(hybrid%par,hybrid%ithr,ghostbuf_tr)
  call t_stopf('pat_remap_ale_gexchV')

  do ie=nets,nete
     num_neighbors = elem(ie)%desc%actual_neigh_edges+1
     do k=1,nlev

        ! find departure points
        call ALE_departure_from_gll     (dep_points, elem(ie)%derived%vstar(:,:,:,k), elem(ie), dt)

        ! find element containing departure point
        call ALE_elems_with_dep_points  (elem_indexes, dep_points, num_neighbors, elem(ie)%desc%neigh_corners)

        ! compute the parametric points
        call ALE_parametric_coords      (para_coords, elem_indexes, dep_points, num_neighbors, elem(ie)%desc%neigh_corners)

        ! for each level k, unpack all tracer neighbor data on that level
        kptr=(k-1)*qsize
        neigh_q=0
        u(:,:,:) = elem(ie)%state%Q(:,:,k,1:qsize)
        call ghostVunpack_unoriented (ghostbuf_tr, neigh_q, np, qsize, kptr, elem(ie)%desc, elem(ie)%GlobalId, u)

        do i=1,np
        do j=1,np
          ! interpolate tracers to deperature grid
          call interpolate_tracers     (para_coords(i,j), neigh_q(:,:,:,elem_indexes(i,j)),f)
          elem(ie)%state%Q(i,j,k,:) = f 
!         call minmax_tracers          (para_coords(i,j), neigh_q(:,:,:,elem_indexes(i,j)),f,g)
          do q=1,qsize
            f(q) = MINVAL(neigh_q(:,:,q,elem_indexes(i,j)))
            g(q) = MAXVAL(neigh_q(:,:,q,elem_indexes(i,j)))
          end do 
          minq(i,j,k,:,ie) = f 
          maxq(i,j,k,:,ie) = g 
        enddo
        enddo

     end do
  end do


  call t_stopf('Prim_Advec_Tracers_remap_ALE_ghost_exchange')
  ! compute original mass, at tl_1%n0
  elem_mass = 0
  do ie=nets,nete
    n=0
    do k=1,nlev
    do q=1,qsize
      n=n+1
      global_shared_buf(ie,n) = 0
      do j=1,np
        global_shared_buf(ie,n) = global_shared_buf(ie,n) + DOT_PRODUCT(elem(ie)%state%Qdp(:,j,k,q,n0_qdp),elem(ie)%spheremp(:,j))
        elem_mass(k,q,ie)       = elem_mass(k,q,ie)       + DOT_PRODUCT(elem(ie)%state%Qdp(:,j,k,q,n0_qdp),elem(ie)%spheremp(:,j))
      end do
    end do
    end do
  end do
  call wrap_repro_sum(nvars=n, comm=hybrid%par%comm)
  n=0
  do k=1,nlev
  do q=1,qsize
    n=n+1
    mass(k,q) = global_shared_sum(n)
  enddo
  enddo


  do ie=nets,nete
  do k=1,nlev
    rho(:,:,k,ie) = elem(ie)%spheremp(:,:)*elem(ie)%state%dp3d(:,:,k,tl%np1)
  end do
  end do

  do ie=nets,nete
    Que_t(:,:,:,1:qsize,ie) = elem(ie)%state%Q(:,:,:,1:qsize)
  end do
  call t_startf('Prim_Advec_Tracers_remap_ALE_Cobra')
  if (use_semi_lagrange_transport_local_conservation) then
    call Cobra_Elem (Que, Que_t, rho, minq, maxq, elem_mass, hybrid, nets, nete)
  else
    call Cobra_SLBQP(Que, Que_t, rho, minq, maxq, mass, hybrid, nets, nete)
  end if
  call t_stopf('Prim_Advec_Tracers_remap_ALE_Cobra')

  do ie=nets,nete
    elem(ie)%state%Q(:,:,:,1:qsize) =  Que(:,:,:,1:qsize,ie)
  end do


  do ie=nets,nete
  do k=1,nlev
  do q=1,qsize
     ! note: pack so that tracers per level are contiguous so we can unpack into
     ! array neigh_q()
     elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%state%Q(:,:,k,q) * elem(ie)%state%dp3d(:,:,k,tl%np1)
  enddo
  enddo
  end do


! do ie=nets,nete
!    do q = 1 , qsize
!       do k = 1 , nlev    !  Potential loop inversion (AAM)
!          elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%spheremp(:,:)* elem(ie)%state%Qdp(:,:,k,q,np1_qdp)
!       enddo
!    enddo
!     call edgeVpack(edgeAdv    , elem(ie)%state%Qdp(:,:,:,:,np1_qdp) , nlev*qsize , 0 , elem(ie)%desc )
! enddo
! call bndry_exchangeV( hybrid , edgeAdv    )
! do ie = nets , nete
!    call edgeVunpack( edgeAdv    , elem(ie)%state%Qdp(:,:,:,:,np1_qdp) , nlev*qsize , 0 , elem(ie)%desc )
!    do q = 1 , qsize
!       do k = 1 , nlev    !  Potential loop inversion (AAM)
!          elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%rspheremp(:,:) * elem(ie)%state%Qdp(:,:,k,q,np1_qdp)
!       enddo
!    enddo
! enddo
  call t_stopf('Prim_Advec_Tracers_remap_ALE')


end subroutine Prim_Advec_Tracers_remap_ALE

subroutine VDOT(rp,Que,rho,mass,hybrid,nets,nete)
  use parallel_mod,        only: global_shared_buf, global_shared_sum
  use global_norms_mod,    only: wrap_repro_sum

  implicit none
  integer             , intent(in)              :: nets
  integer             , intent(in)              :: nete
  real(kind=real_kind), intent(out)             :: rp                (nlev,qsize)
  real(kind=real_kind), intent(in)              :: Que         (np*np,nlev,qsize,nets:nete)
  real(kind=real_kind), intent(in)              :: rho         (np*np,nlev,      nets:nete)
  real(kind=real_kind), intent(in)              :: mass              (nlev,qsize)
  type (hybrid_t)     , intent(in)              :: hybrid

  integer                                       :: k,n,q,ie

  global_shared_buf = 0 
  do ie=nets,nete
    n=0
    do q=1,qsize
    do k=1,nlev
      n=n+1
      global_shared_buf(ie,n) = global_shared_buf(ie,n) + DOT_PRODUCT(Que(:,k,q,ie), rho(:,k,ie))
    end do
    end do
  end do

  call wrap_repro_sum(nvars=n, comm=hybrid%par%comm)

  n=0
  do q=1,qsize
  do k=1,nlev
    n=n+1
    rp(k,q) = global_shared_sum(n) - mass(k,q)
  enddo
  enddo
  
end subroutine VDOT

subroutine Cobra_SLBQP(Que, Que_t, rho, minq, maxq, mass, hybrid, nets, nete) 

  use parallel_mod,        only: global_shared_buf, global_shared_sum
  use global_norms_mod,    only: wrap_repro_sum

  implicit none
  integer             , intent(in)              :: nets
  integer             , intent(in)              :: nete
  real(kind=real_kind), intent(out)             :: Que         (np*np,nlev,qsize,nets:nete)
  real(kind=real_kind), intent(in)              :: Que_t       (np*np,nlev,qsize,nets:nete)
  real(kind=real_kind), intent(in)              :: rho         (np*np,nlev,      nets:nete)
  real(kind=real_kind), intent(in)              :: minq        (np*np,nlev,qsize,nets:nete)
  real(kind=real_kind), intent(in)              :: maxq        (np*np,nlev,qsize,nets:nete)
  real(kind=real_kind), intent(in)              :: mass              (nlev,qsize)
  type (hybrid_t)     , intent(in)              :: hybrid

  integer,                            parameter :: max_clip = 50
  real(kind=real_kind),               parameter :: eta = 1D-08           
  real(kind=real_kind),               parameter :: hfd = 1D-10             
  real(kind=real_kind)                          :: lambda_p          (nlev,qsize)
  real(kind=real_kind)                          :: lambda_c          (nlev,qsize)
  real(kind=real_kind)                          :: rp                (nlev,qsize)
  real(kind=real_kind)                          :: rc                (nlev,qsize)
  real(kind=real_kind)                          :: rd                (nlev,qsize)
  real(kind=real_kind)                          :: alpha             (nlev,qsize)
  integer                                       :: j,k,n,q,ie
  integer                                       :: nclip

  nclip = 0

  Que(:,:,:,:) = Que_t(:,:,:,:)

  Que = MIN(MAX(Que,minq),maxq)

  call VDOT(rp,Que,rho,mass,hybrid,nets,nete)
  nclip = nclip + 1

  if (MAXVAL(ABS(rp)).lt.eta) return

  do ie=nets,nete
  do q=1,qsize
  do k=1,nlev
     Que(:,k,q,ie) = hfd * rho(:,k,ie) + Que_t(:,k,q,ie)
  enddo
  enddo
  enddo

  Que = MIN(MAX(Que,minq),maxq)

  call VDOT(rc,Que,rho,mass,hybrid,nets,nete)

  rd = rc-rp
  if (MAXVAL(ABS(rd)).eq.0) return 
  
  alpha = 0
  WHERE (rd.ne.0) alpha = hfd / rd 

  lambda_p = 0
  lambda_c =  -alpha*rp

  do while (MAXVAL(ABS(rc)).gt.eta .and. nclip.lt.max_clip)

    do ie=nets,nete
    do q=1,qsize
    do k=1,nlev
       Que(:,k,q,ie) = (lambda_c(k,q) + hfd) * rho(:,k,ie) + Que_t(:,k,q,ie)
    enddo
    enddo
    enddo
    Que = MIN(MAX(Que,minq),maxq)

    call VDOT(rc,Que,rho,mass,hybrid,nets,nete)
    nclip = nclip + 1

    rd = rp-rc

    if (MAXVAL(ABS(rd)).eq.0) exit

    alpha = 0
    WHERE (rd.ne.0) alpha = (lambda_p - lambda_c) / rd 

    rp       = rc
    lambda_p = lambda_c

    lambda_c = lambda_c -  alpha * rc

  enddo
end subroutine Cobra_SLBQP


subroutine Cobra_Elem(Que, Que_t, rho, minq, maxq, mass, hybrid, nets, nete) 

  use parallel_mod,        only: global_shared_buf, global_shared_sum
  use global_norms_mod,    only: wrap_repro_sum

  implicit none
  integer             , intent(in)              :: nets
  integer             , intent(in)              :: nete
  real(kind=real_kind), intent(out)             :: Que         (np*np,nlev,qsize,nets:nete)
  real(kind=real_kind), intent(in)              :: Que_t       (np*np,nlev,qsize,nets:nete)
  real(kind=real_kind), intent(in)              :: rho         (np*np,nlev,      nets:nete)
  real(kind=real_kind), intent(in)              :: minq        (np*np,nlev,qsize,nets:nete)
  real(kind=real_kind), intent(in)              :: maxq        (np*np,nlev,qsize,nets:nete)
  real(kind=real_kind), intent(in)              :: mass              (nlev,qsize,nets:nete)
  type (hybrid_t)     , intent(in)              :: hybrid

  integer,                            parameter :: max_clip = 50
  real(kind=real_kind),               parameter :: eta = 1D-10           
  real(kind=real_kind),               parameter :: hfd = 1D-08             
  real(kind=real_kind)                          :: lambda_p          (nlev,qsize,nets:nete)
  real(kind=real_kind)                          :: lambda_c          (nlev,qsize,nets:nete)
  real(kind=real_kind)                          :: rp                (nlev,qsize,nets:nete)
  real(kind=real_kind)                          :: rc                (nlev,qsize,nets:nete)
  real(kind=real_kind)                          :: rd                (nlev,qsize,nets:nete)
  real(kind=real_kind)                          :: alpha             (nlev,qsize,nets:nete)
  integer                                       :: j,k,n,q,ie
  integer                                       :: nclip
  integer                                       :: mloc(3)

  nclip = 1

  Que(:,:,:,:) = Que_t(:,:,:,:)

  Que = MIN(MAX(Que,minq),maxq)

  do ie=nets,nete
  do q=1,qsize
  do k=1,nlev
    rp(k,q,ie) = DOT_PRODUCT(Que(:,k,q,ie), rho(:,k,ie)) - mass(k,q,ie)
  end do
  end do
  end do

  if (MAXVAL(ABS(rp)).lt.eta) return

  do ie=nets,nete
  do q=1,qsize
  do k=1,nlev
     Que(:,k,q,ie) = hfd * rho(:,k,ie) + Que_t(:,k,q,ie)
  enddo
  enddo
  enddo

  Que = MIN(MAX(Que,minq),maxq)

  do ie=nets,nete
  do q=1,qsize
  do k=1,nlev
    rc(k,q,ie) = DOT_PRODUCT(Que(:,k,q,ie), rho(:,k,ie)) - mass(k,q,ie)
  end do
  end do
  end do

  rd = rc-rp
  if (MAXVAL(ABS(rd)).eq.0) return 
  
  alpha = 0
  WHERE (rd.ne.0) alpha = hfd / rd 

  lambda_p = 0
  lambda_c =  -alpha*rp

! if (hybrid%par%masterproc) print *,__FILE__,__LINE__," mass(20,1,4):",mass(20,1,4)
! do k=1,np*np
!   if (hybrid%par%masterproc) print *,__FILE__,__LINE__," maxq(k,20,1,4):", &
!     maxq(k,20,1,4) ,minq(k,20,1,4),maxq(k,20,1,4)-minq(k,20,1,4)
! enddo
! do k=1,np*np
!   if (hybrid%par%masterproc) print *,__FILE__,__LINE__," Que(k,20,1,4):",Que(k,20,1,4) ,rho(k,20,4)
! enddo
  do while (MAXVAL(ABS(rc)).gt.eta .and. nclip.lt.max_clip)
    nclip = nclip + 1

    do ie=nets,nete
    do q=1,qsize
    do k=1,nlev
!      Que(:,k,q,ie) = (lambda_c(k,q,ie) + hfd) * rho(:,k,ie) + Que_t(:,k,q,ie)
       Que(:,k,q,ie) = lambda_c(k,q,ie) * rho(:,k,ie) + Que_t(:,k,q,ie)
    enddo
    enddo
    enddo

!   do ie=nets,nete
!   do q=1,qsize
!   do k=1,nlev
!     rc(k,q,ie) = DOT_PRODUCT(Que(:,k,q,ie), rho(:,k,ie)) - mass(k,q,ie)
!   end do
!   end do
!   end do
!   if (hybrid%par%masterproc) print *,__FILE__,__LINE__," rc(20,1,4):",rc(20,1,4), DOT_PRODUCT(Que(:,20,1,4), rho(:,20,4))

    Que = MIN(MAX(Que,minq),maxq)

    do ie=nets,nete
    do q=1,qsize
    do k=1,nlev
      rc(k,q,ie) = DOT_PRODUCT(Que(:,k,q,ie), rho(:,k,ie)) - mass(k,q,ie)
    end do
    end do
    end do

    
    mloc = MAXLOC(ABS(rc))
!   if (hybrid%par%masterproc) print *,__FILE__,__LINE__," MAXVAL(ABS(rc)):",MAXVAL(ABS(rc)), mloc, nclip

    rd = rp-rc

!   if (MAXVAL(ABS(rd)).eq.0) exit

    alpha = 0
    WHERE (rd.ne.0) alpha = (lambda_p - lambda_c) / rd 
!   WHERE (alpha.eq.0.and.MAXVAL(ABS(rc)).gt.eta) alpha=10;

    rp       = rc
    lambda_p = lambda_c

    lambda_c = lambda_c -  alpha * rc

!   if (hybrid%par%masterproc) print *,__FILE__,__LINE__," rc(20,1,4):",rc(20,1,4)
!   if (hybrid%par%masterproc) print *,__FILE__,__LINE__," rd(20,1,4):",rd(20,1,4)
!   if (hybrid%par%masterproc) print *,__FILE__,__LINE__," lambda_p(20,1,4):",lambda_p(20,1,4)
!   if (hybrid%par%masterproc) print *,__FILE__,__LINE__," lambda_c(20,1,4):",lambda_c(20,1,4)
!   if (hybrid%par%masterproc) print *,__FILE__,__LINE__," alpha(20,1,4):",alpha(20,1,4)
!   if (hybrid%par%masterproc) print *
  enddo
! if (hybrid%par%masterproc) print *,__FILE__,__LINE__," MAXVAL(ABS(rc)):",MAXVAL(ABS(rc)),eta," nclip:",nclip
end subroutine Cobra_Elem 

! ----------------------------------------------------------------------------------!
!SUBROUTINE ALE_RKDSS
! AUTHOR: CHRISTOPH ERATH, MARK TAYLOR, 06. December 2012
!
! DESCRIPTION: ! create a runge kutta taylor serios mixture to calculate the departure grid
!
! CALLS:
! INPUT:
!
! OUTPUT:
!-----------------------------------------------------------------------------------!

! this will calculate the velocity at time t+1/2  along the trajectory s(t) given the velocities
! at the GLL points at time t and t+1 using a second order time accurate formulation.


subroutine ALE_RKdss(elem, nets, nete, hy, deriv, dt, tl)

  use derivative_mod,  only : derivative_t, ugradv_sphere
  use edgetype_mod,    only : EdgeBuffer_t
  use bndry_mod,       only : bndry_exchangev
  use kinds,           only : real_kind
  use hybrid_mod,      only : hybrid_t
  use element_mod,     only : element_t
  use dimensions_mod,   only : np, nlev

  implicit none

  type (element_t)     , intent(inout)             :: elem(:)
  integer              , intent(in   )             :: nets
  integer              , intent(in   )             :: nete
  type (hybrid_t)      , intent(in)                :: hy ! distributed parallel structure (shared)
  type (derivative_t)  , intent(in)                :: deriv ! derivative struct
  real (kind=real_kind), intent(in)                :: dt ! timestep
  type (TimeLevel_t)   , intent(in)                :: tl

  integer                                          :: ie, k
  real (kind=real_kind), dimension(np,np,2)        :: vtmp
  integer :: np1

  np1 = tl%np1


! RK-SSP 2 stage 2nd order:
!     x*(t+1) = x(t) + U(x(t),t) dt
!     x(t+1) = x(t) +  1/2 ( U(x*(t+1),t+1) + U(x(t),t) ) dt
! apply taylor series:
!  U(x*(t+1),t+1) = U(x(t),t+1) + (x*(t+1)-x(t)) gradU(x(t),t+1)
!
! x(t+1) = x(t) +  1/2 ( U(x(t),t+1) + (x*(t+1)-x(t)) gradU(x(t),t+1) + U(x(t),t) ) dt
! (x(t+1) - x(t)) / dt =  1/2 ( U(x(t),t+1) + (x*(t+1)-x(t)) gradU(x(t),t+1) + U(x(t),t) )
! (x(t+1) - x(t)) / dt =  1/2 ( U(x(t),t+1) + U(x(t),t) + (x*(t+1)-x(t)) gradU(x(t),t+1) )
! (x(t+1) - x(t)) / dt =  1/2 ( U(x(t),t+1) + U(x(t),t) + U(x(t),t) dt  gradU(x(t),t+1) )
!
!
!  (x(t+1)-x(t))/dt =  1/2(U(x(t),t+1) + U(x(t),t) + dt U(x(t),t) gradU(x(t),t+1))
!
! suppose dt = -ts (we go backward)
!  (x(t-ts)-x(t))/-ts =  1/2( U(x(t),t-ts)+U(x(t),t)) - ts 1/2 U(x(t),t) gradU(x(t),t-ts)
!
!  x(t-ts) = x(t)) -ts * [ 1/2( U(x(t),t-ts)+U(x(t),t)) - ts 1/2 U(x(t),t) gradU(x(t),t-ts) ]
!
!    !------------------------------------------------------------------------------------

  do ie=nets,nete
     ! vstarn0 = U(x,t)
     ! vstar   = U(x,t+1)
    do k=1,nlev
      vtmp(:,:,:)=ugradv_sphere(elem(ie)%state%v(:,:,:,k,np1), elem(ie)%derived%vstar(:,:,:,k),deriv,elem(ie))

      elem(ie)%derived%vstar(:,:,:,k) = &
           (elem(ie)%state%v(:,:,:,k,np1) + elem(ie)%derived%vstar(:,:,:,k))/2 - dt*vtmp(:,:,:)/2

      elem(ie)%derived%vstar(:,:,1,k) = elem(ie)%derived%vstar(:,:,1,k)*elem(ie)%spheremp(:,:)
      elem(ie)%derived%vstar(:,:,2,k) = elem(ie)%derived%vstar(:,:,2,k)*elem(ie)%spheremp(:,:)
    enddo
    call edgeVpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%derived%vstar,2*nlev,0,2*nlev)
  enddo

  call t_startf('ALE_RKdss_bexchV')
  call bndry_exchangeV(hy,edge_g)
  call t_stopf('ALE_RKdss_bexchV')

  do ie=nets,nete
    call edgeVunpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%derived%vstar,2*nlev,0,2*nlev)
    do k=1, nlev
      elem(ie)%derived%vstar(:,:,1,k) = elem(ie)%derived%vstar(:,:,1,k)*elem(ie)%rspheremp(:,:)
      elem(ie)%derived%vstar(:,:,2,k) = elem(ie)%derived%vstar(:,:,2,k)*elem(ie)%rspheremp(:,:)
    end do
  end do
end subroutine ALE_RKdss

! ----------------------------------------------------------------------------------!
!SUBROUTINE FVM_DEP_FROM_GLL----------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, MARK TAYLOR 14. December 2011                            !
! DESCRIPTION: calculates the deparute grid for fvm coming from the gll points      !
!                                                                                   !
! CALLS:
! INPUT:
!
! OUTPUT:
!-----------------------------------------------------------------------------------!
subroutine ALE_departure_from_gll(acart, vstar, elem, dt)
  use physical_constants,     only : rearth
  use coordinate_systems_mod, only : spherical_polar_t, cartesian3D_t, change_coordinates
  use time_mod,               only : timelevel_t
  use element_mod,            only : element_t
  use kinds,                  only : real_kind
  use dimensions_mod,         only : np

  implicit none

  type(cartesian3D_t)     ,intent(out)  :: acart(np,np)
  real (kind=real_kind)   ,intent(in)   :: vstar(np,np,2)
  type (element_t)        ,intent(in)   :: elem
  real (kind=real_kind)   ,intent(in)   :: dt

  integer                               :: i,j

  real (kind=real_kind)                 :: uxyz (np,np,3)

   ! convert velocity from lat/lon to cartesian 3D

  do i=1,3
     ! Summing along the third dimension is a sum over components for each point.
     ! (This is just a faster way of doing a dot product for each grid point,
     ! since reindexing the inputs to use the intrinsic effectively would be
     ! just asking for trouble.)
     uxyz(:,:,i)=sum( elem%vec_sphere2cart(:,:,i,:)*vstar(:,:,:) ,3)
  end do
  ! interpolate velocity to fvm nodes
  ! compute departure point
  ! crude, 1st order accurate approximation.  to be improved
  do i=1,np
     do j=1,np
        acart(i,j) = change_coordinates(elem%spherep(i,j)) 
        acart(i,j)%x = acart(i,j)%x - dt*uxyz(i,j,1)/rearth
        acart(i,j)%y = acart(i,j)%y - dt*uxyz(i,j,2)/rearth
        acart(i,j)%z = acart(i,j)%z - dt*uxyz(i,j,3)/rearth
     enddo
  enddo

end subroutine ALE_departure_from_gll




subroutine ALE_elems_with_dep_points (elem_indexes, dep_points, num_neighbors, ngh_corners)

  use element_mod,            only : element_t
  use dimensions_mod,         only : np
  use coordinate_systems_mod, only : cartesian3D_t, change_coordinates
  use interpolate_mod,        only : point_inside_quad

  implicit none

  ! The ngh_corners array is a list of corners of both elem and all of it's
  ! neighor elements all sorted by global id.
  integer              , intent(in)                :: num_neighbors
  type(cartesian3D_t),intent(in)                   :: ngh_corners(4,num_neighbors)
  integer              , intent(out)               :: elem_indexes(np,np)
  type(cartesian3D_t)  , intent(in)                :: dep_points(np,np)

  integer                                          :: i,j,n
  logical                                          :: inside

  elem_indexes = -1
  do i=1,np
    do j=1,np
! Just itererate the neighbors in global id order to get the same result on every processor.
      do n = 1, num_neighbors
! Mark Taylor's handy dandy point_inside_gc check.
       inside = point_inside_quad (ngh_corners(:,n), dep_points(i,j))
       if (inside) then
         elem_indexes(i,j) = n
         exit
       end if
     end do
    end do
  end do

  if (MINVAL(elem_indexes(:,:))==-1) then
    write (*,*) __FILE__,__LINE__,"Aborting because point not found in neighbor list. Info:"
    do i=1,np
      do j=1,np
        if (elem_indexes(i,j)==-1) then
          write (*,*)   " departure point ",dep_points(i,j)
          do n = 1, num_neighbors
            write (*,*) " quad checked    ",ngh_corners(1,n)
            write (*,*) "                 ",ngh_corners(2,n)
            write (*,*) "                 ",ngh_corners(3,n)
            write (*,*) "                 ",ngh_corners(4,n)
            write (*,*)
          end do
          exit
        end if
      end do
    end do
    call abortmp("ERROR elems_with_dep_points: Can't find departure grid. Time step too long?")
  end if
end subroutine ALE_elems_with_dep_points

function  shape_fcn_deriv(pc) result(dNds)
  real (kind=real_kind), intent(in)  ::  pc(2)
  real (kind=real_kind)              :: dNds(4,2)
 
  dNds(1, 1) = - 0.25 * (1.0 - pc(2))
  dNds(1, 2) = - 0.25 * (1.0 - pc(1))

  dNds(2, 1) =   0.25 * (1.0 - pc(2))
  dNds(2, 2) = - 0.25 * (1.0 + pc(1))

  dNds(3, 1) =   0.25 * (1.0 + pc(2))
  dNds(3, 2) =   0.25 * (1.0 + pc(1))

  dNds(4, 1) = - 0.25 * (1.0 + pc(2))
  dNds(4, 2) =   0.25 * (1.0 - pc(1))
end function   

function inv_2x2(A) result(A_inv)
  real (kind=real_kind), intent(in)  :: A    (2,2)
  real (kind=real_kind)              :: A_inv(2,2)
  real (kind=real_kind) :: det, denom

  det = A(1,1) * A(2,2) - A(2,1) * A(1,2)
  denom = 1/det
  ! inverse:
  A_inv(1,1) =  denom * A(2,2)  !  dxidx
  A_inv(2,1) = -denom * A(2,1)  !  detadx
  A_inv(1,2) = -denom * A(1,2)  !  dxidy
  A_inv(2,2) =  denom * A(1,1)  !  detady
end function

function INV(dxds) result(dsdx)

  real (kind=real_kind), intent(in)  :: dxds(3,2)

  real (kind=real_kind)  ::     dsdx(2,3)
  real (kind=real_kind)  ::      ata(2,2)
  real (kind=real_kind)  ::  ata_inv(2,2)


  !     dxds = | dxdxi   dxdeta |
  !            | dydxi   dydeta |
  !            | dzdxi   dzdeta |
  ata  = MATMUL(TRANSPOSE(dxds), dxds)
  ata_inv = inv_2x2(ata)
  dsdx = MATMUL(ata_inv, TRANSPOSE(dxds))
  !     dsdx = |  dxidx   dxidy   dxidz |
  !            | detadx  detady  detadz |

end function

subroutine shape_fcn(N, pc)
  real (kind=real_kind), intent(out) :: N(4)
  real (kind=real_kind), intent(in)  :: pc(2)

  ! shape function for each node evaluated at param_coords
  N(1) = 0.25 * (1.0 - pc(1)) * (1.0 - pc(2)) 
  N(2) = 0.25 * (1.0 + pc(1)) * (1.0 - pc(2)) 
  N(3) = 0.25 * (1.0 + pc(1)) * (1.0 + pc(2)) 
  N(4) = 0.25 * (1.0 - pc(1)) * (1.0 + pc(2)) 
end subroutine


function F(coords, pc) result(x)
  real (kind=real_kind), intent(in) :: pc(2), coords(4,3)

  real (kind=real_kind)            :: N(4), x(3)
  call shape_fcn(N,pc)
  x = MATMUL(TRANSPOSE(coords), N)
  x = x/SQRT(DOT_PRODUCT(x,x))
end function

function  DF(coords, pc) result(dxds)
  real (kind=real_kind), intent(in)  :: coords(4,3)
  real (kind=real_kind), intent(in)  :: pc(2)
 
  real (kind=real_kind)              :: dxds(3,2)
  real (kind=real_kind)              :: dNds(4,2)
  real (kind=real_kind)              ::  dds(3,2)
  real (kind=real_kind)              ::    c(2)
  real (kind=real_kind)              ::    x(3)
  real (kind=real_kind)              ::   xc(3,2)
  real (kind=real_kind)              :: nx, nx2 
  integer                            :: i,j

  dNds = shape_fcn_deriv  (pc)
  dds  = MATMUL(TRANSPOSE(coords), dNds)

  x    = F(coords, pc)
  nx2  = DOT_PRODUCT(x,x)
  nx   = SQRT(nx2)
  c    = MATMUL(TRANSPOSE(dds), x)
  do j=1,2
    do i=1,3
      xc(i,j) = x(i)*c(j)
    end do
  end do
  dxds = nx2*dds - xc
  dxds = dxds/(nx*nx2)
end function   


function cartesian_parametric_coordinates(sphere, corners3D) result (ref)
  use coordinate_systems_mod, only : cartesian2d_t, cartesian3D_t, spherical_polar_t, spherical_to_cart
  implicit none
  type (spherical_polar_t), intent(in) :: sphere
  type (cartesian3D_t)    , intent(in) :: corners3D(4)  !x,y,z coords of element corners

  type (cartesian2D_t)                 :: ref

  integer,               parameter :: MAXIT = 20
  real (kind=real_kind), parameter :: TOL   = 1.0E-13
  integer,               parameter :: n     = 3

  type (cartesian3D_t)             :: cart
  real (kind=real_kind)            :: coords(4,3), dxds(3,2), dsdx(2,3)
  real (kind=real_kind)            :: p(3), pc(2), dx(3), x(3), ds(2)
  real (kind=real_kind)            :: dist, step                          
  
  integer                          :: i,j,k,iter
  do i=1,4                               
    coords(i,1) = corners3D(i)%x 
    coords(i,2) = corners3D(i)%y 
    coords(i,3) = corners3D(i)%z 
  end do

  pc = 0
  p  = 0
  cart = spherical_to_cart(sphere)

  p(1) = cart%x
  p(2) = cart%y
  p(3) = cart%z 

  dx   = 0
  ds   = 0
  dsdx = 0
  dxds = 0

  !*-------------------------------------------------------------------------*!

  ! Initial guess, center of element
  dist = 9999999.                         
  step = 9999999.                         
  iter = 0 

  do while  (TOL*TOL.lt.dist .and. iter.lt.MAXIT .and. TOL*TOL.lt.step)
    iter = iter + 1

    dxds =  DF (coords, pc)
    x    =   F (coords, pc)
    dsdx = INV (dxds)

    dx   = x - p
    dist = DOT_PRODUCT(dx,dx)
    ds   = MATMUL(dsdx, dx)
    pc   = pc - ds
    step = DOT_PRODUCT(ds,ds)
  enddo

  ref%x = pc(1)
  ref%y = pc(2)
end function


subroutine  ALE_parametric_coords (parametric_coord, elem_indexes, dep_points, num_neighbors, ngh_corners)
  use coordinate_systems_mod, only : cartesian2d_t, cartesian3D_t, spherical_polar_t, change_coordinates, distance
  use interpolate_mod,        only : parametric_coordinates
  use dimensions_mod,         only : np

  implicit none

  type(cartesian2D_t)       , intent(out)       :: parametric_coord(np,np)
  type(cartesian3D_t)       , intent(in)        :: dep_points(np,np)
  integer                   , intent(in)        :: elem_indexes(np,np)
  integer                   , intent(in)        :: num_neighbors
  type(cartesian3D_t)       , intent(in)        :: ngh_corners(4,num_neighbors)

  type (spherical_polar_t)                      :: sphere(np,np)
  integer                                       :: i,j,n
  type(cartesian2D_t)                           :: parametric_test(np,np)
  real(kind=real_kind)                          :: d

  do j=1,np
    sphere(:,j) = change_coordinates(dep_points(:,j))
  end do

  do i=1,np
    do j=1,np
      n = elem_indexes(i,j)
      parametric_coord(i,j)= parametric_coordinates(sphere(i,j),ngh_corners(:,n))
    end do
  end do
end subroutine ALE_parametric_coords


end module 
