#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module advance_mod
  use kinds, only : real_kind
  use dimensions_mod, only : np,npdg
  use physical_constants, only : rearth, rrearth 
  implicit none
  ! semi-implicit needs to be re-initialized each time dt changes 
  real (kind=real_kind) :: initialized_for_dt   = 0
contains


  subroutine advance_nonstag( elem, edge2,  edge3,  deriv,  flt,   hybrid,  &
       dt,  pmean,     tl,   nets,   nete)

    ! ---------------------
    use kinds, only : real_kind
    ! ---------------------
    use dimensions_mod, only : np, nlev
    ! ---------------------
    use element_mod, only : element_t
    ! ---------------------
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack
    ! ---------------------
    use filter_mod, only : filter_t, filter_P
    ! ---------------------
    use hybrid_mod, only : hybrid_t
    ! ---------------------
    use reduction_mod, only : reductionbuffer_ordered_1d_t
    ! ---------------------
    use derivative_mod, only : derivative_t, gradient_sphere, divergence_sphere,vorticity_sphere,&
         divergence_sphere_wk
    ! ---------------------
    use time_mod, only : timelevel_t, smooth
    ! ---------------------
    use control_mod, only :  filter_freq, filter_counter, topology, test_case, LFTfreq
    ! ---------------------
    use shallow_water_mod, only : tc1_velocity, vortex_velocity, swirl_velocity
    ! ---------------------
    use cg_mod, only : cg_t
    ! ---------------------
    use bndry_mod, only : bndry_exchangev
    use viscosity_mod, only : neighbor_minmax
    ! ---------------------
    !  FOR DEBUGING use only 
    ! ---------------------
    !    use schedule_mod
    use global_norms_mod
    ! ---------------------
    use perf_mod, only : t_startf, t_stopf ! _EXTERNAL


    implicit none

    type (element_t)     , intent(inout), target :: elem(:)
    type (EdgeBuffer_t)  , intent(in) :: edge2
    type (EdgeBuffer_t)  , intent(inout) :: edge3

    type (derivative_t)  , intent(in) :: deriv

    type (filter_t)                   :: flt
    type (hybrid_t)      , intent(in) :: hybrid

    real (kind=real_kind), intent(in) :: dt
    real (kind=real_kind), intent(in) :: pmean
    type (TimeLevel_t)   , intent(in) :: tl
    integer              , intent(in) :: nets
    integer              , intent(in) :: nete
    integer               :: ig

    ! =================
    ! Local
    ! =================
    ! Thread private working set ...
    real (kind=real_kind), dimension(:,:), pointer     :: rspheremp,spheremp
    real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)  :: vtens
    real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: ptens
    logical :: do_leapfrog


    real (kind=real_kind) :: dt2, real_time

    real*8                :: st,et
    integer    :: i,j,k,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep, steptype

    nm1   = tl%nm1
    n0    = tl%n0
    np1   = tl%np1
    nstep = tl%nstep

    dt2 = 2.0d0*dt

    real_time = dt*real(nstep,kind=real_kind)

    call t_startf('advance_nonstag')



    ! LFTfreq=0   pure Leapfrog (default)
    ! LFTfreq=1   pure Leapfrog-trapazoidal
    ! LFTfreq=n   RK2, then n-1 leapfrogs

    ! steptype = 0  leapfrog
    ! steptype = 1  leapfrog-trap
    ! steptype = N  RK2, followed by N-1 leapfrogs
    if (LFTfreq==0) then
       steptype=0
    elseif (LFTfreq==1) then
       steptype=1
    else
       steptype=2  
       if (mod(nstep,LFTfreq).ne.0) steptype=0
    endif

    ! in all cases, use LF during during bootstrap (nstep=0) phase
    ! For RK methods we should remove bootstrap procedure
    if (nstep==0) steptype=0

    if (steptype==0) then   

       ! Leapfrog timestep: u(np1) = u(nm1) + dt2*DSS [ RHS(u(n0)) ]
       call compute_and_apply_rhs(np1,nm1,n0,dt2,real_time,edge3,elem,pmean,hybrid,deriv,vtens,ptens,nets,nete)

       ! ====================================================
       ! apply viscosity  
       ! ====================================================
       call advance_hypervis(edge3,elem,hybrid,deriv,vtens,ptens,np1,nets,nete,dt2)

    else if (steptype==1)  then
       if (smooth/=0) stop 'ERROR: smooth>0 only allowed for leapfrog'
       ! leapfrog+trapazoidal
       ! 2x as expensive as LF, but 2nd order, no Robert filter needed, 
       ! dt sqrt(2) larger than LF
       ! u(*) = u(n-1) + 2dt F(u(n))   u(*) is at time level n+1
       ! u(n+1) = u(n) + dt [ F(u(n)) + F(u(*)) ] /2
       
       ! u(n+1) = u(n) + dt/2 F(u(n))   
       call compute_and_apply_rhs(np1,n0,n0,dt/2,real_time,edge3,elem,pmean,hybrid,deriv,vtens,ptens,nets,nete)

       ! u(n-1) = u(n-1) + 4( u(n+1)-u(n))     u(*) above
       do ie=nets,nete
          elem(ie)%state%v(:,:,:,:,nm1)  = elem(ie)%state%v(:,:,:,:,nm1) + &
               4*(elem(ie)%state%v(:,:,:,:,np1)-elem(ie)%state%v(:,:,:,:,n0)  )
          elem(ie)%state%p(:,:,:,nm1)  = elem(ie)%state%p(:,:,:,nm1) + &
               4*(elem(ie)%state%p(:,:,:,np1)-elem(ie)%state%p(:,:,:,n0)  )
       enddo

       ! u(n+1) = u(n+1) + dt/2 F(u(*))        
       call compute_and_apply_rhs(np1,np1,nm1,dt/2,real_time+dt,edge3,elem,pmean,hybrid,deriv,vtens,ptens,nets,nete)
       ! ====================================================
       ! apply viscosity  Note: use dt, not dt/2
       ! ====================================================
       call advance_hypervis(edge3,elem,hybrid,deriv,vtens,ptens,np1,nets,nete,dt)

    else if (steptype==2) then
       if (smooth/=0) stop 'ERROR: smooth>0 only allowed for leapfrog'
       ! RK2 (which is forward euler at dt/2 followed by LF with dt)
       ! Foward Euler  u(n0) -> u(np1) at t+.5
       call compute_and_apply_rhs(np1,n0,n0,dt/2,real_time,edge3,elem,pmean,hybrid,deriv,vtens,ptens,nets,nete)
       ! leapfrog:  u(dt) = u(n0) + dt RHS(dt/2)     (store in u(np1))
       call compute_and_apply_rhs(np1,n0,np1,dt,real_time+dt,edge3,elem,pmean,hybrid,deriv,vtens,ptens,nets,nete)

       call advance_hypervis(edge3,elem,hybrid,deriv,vtens,ptens,np1,nets,nete,dt)
    endif




    if (smooth/=0) then
    do ie=nets,nete
       ! ====================================================
       ! apply Robert filter
       ! ====================================================
       do k=1,nlev
          do j=1,np
             do i=1,np
                elem(ie)%state%v(i,j,1,k,n0)  = elem(ie)%state%v(i,j,1,k,n0) + smooth*(elem(ie)%state%v(i,j,1,k,nm1) &
                     - 2.0D0*elem(ie)%state%v(i,j,1,k,n0) + elem(ie)%state%v(i,j,1,k,np1))
                elem(ie)%state%v(i,j,2,k,n0)  = elem(ie)%state%v(i,j,2,k,n0) + smooth*(elem(ie)%state%v(i,j,2,k,nm1) &
                     - 2.0D0*elem(ie)%state%v(i,j,2,k,n0) + elem(ie)%state%v(i,j,2,k,np1))
                elem(ie)%state%p(i,j,k,n0)  = elem(ie)%state%p(i,j,k,n0) + smooth*(elem(ie)%state%p(i,j,k,nm1) &
                     - 2.0D0*elem(ie)%state%p(i,j,k,n0) + elem(ie)%state%p(i,j,k,np1))
             end do
          end do
       end do
    end do
    endif

    call t_stopf('advance_nonstag')

#if (! defined VERT_OPENMP)
    !$OMP BARRIER
#endif
  end subroutine advance_nonstag

!--------------------------------------------------------------------------------------------

  subroutine advance_nonstag_rk( MyRk, elem, edge2,  edge3,  deriv,  flt,   hybrid,  &
       dt,  pmean,     tl,   nets,   nete)

    ! ---------------------
    use kinds, only : real_kind
    ! ---------------------
    use dimensions_mod, only : np, nlev
    ! ---------------------
    use element_mod, only : element_t
    ! ---------------------
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack, edgedgvunpack
    ! ---------------------
    use filter_mod, only : filter_t, filter_P
    ! ---------------------
    use hybrid_mod, only : hybrid_t
    ! ---------------------
    use reduction_mod, only : reductionbuffer_ordered_1d_t
    ! ---------------------
    use derivative_mod, only : derivative_t, gradient_sphere, divergence_sphere,vorticity_sphere,&
         divergence_sphere_wk, edge_flux_u_cg, dgmodal_to_gll, gll_to_dgmodal
    ! ---------------------
    use time_mod, only : timelevel_t, smooth
    ! ---------------------
    use control_mod, only :  filter_freq, filter_counter, topology, test_case, sub_case, &
          limiter_option, nu, nu_div, nu_s, tracer_advection_formulation, TRACERADV_UGRADQ, &
	  kmass
    ! ---------------------
    use shallow_water_mod, only : tc1_velocity, vortex_velocity, vortex_exact, swirl_velocity
    ! ---------------------
    use cg_mod, only : cg_t
    ! ---------------------
    use bndry_mod, only : bndry_exchangev
    use viscosity_mod, only : neighbor_minmax, biharmonic_wk
    ! ---------------------
    !  FOR DEBUGING use only 
    ! ---------------------
    !    use schedule_mod
    use global_norms_mod
    ! ---------------------
    use types_mod, only : rk_t


    implicit none

    type (element_t)     , intent(inout), target :: elem(:)
    type (Rk_t)          , intent(in) :: MyRk
    type (EdgeBuffer_t)  , intent(in) :: edge2
    type (EdgeBuffer_t)  , intent(inout) :: edge3

    type (derivative_t)  , intent(in) :: deriv

    type (filter_t)                   :: flt
    type (hybrid_t)      , intent(in) :: hybrid

    real (kind=real_kind), intent(in) :: dt
    real (kind=real_kind), intent(in) :: pmean
    type (TimeLevel_t)   , intent(inout) :: tl
    integer              , intent(in) :: nets
    integer              , intent(in) :: nete
    integer                           :: ig

    ! =================
    ! Local
    ! =================

    ! pointer ...

    real (kind=real_kind), dimension(:,:), pointer     :: fcor,rspheremp,spheremp,metdet,rmetdet
    real (kind=real_kind), dimension(:,:,:,:), pointer :: met,metinv

    ! Thread private working set ...

    real (kind=real_kind), dimension(np,np,2,nlev,nets:nete) :: vtens
    real (kind=real_kind), dimension(np,np,nlev,nets:nete)   :: ptens
    real (kind=real_kind), dimension(np,np,nlev,nets:nete)   :: ptens_dg
    real (kind=real_kind), dimension(npdg,npdg)   :: phat

    real (kind=real_kind), dimension(0:np+1,0:np+1,nlev)   :: pedges
    real (kind=real_kind), dimension(np,np,nlev)   :: p_storage, ptens_tvd

    real (kind=real_kind), dimension(np,np,2)    :: grade   ! kinetic energy gradient
    real (kind=real_kind), dimension(np,np,2)    :: gradh   ! grad(h)

    real (kind=real_kind), dimension(np,np,2)    :: pv      ! p*v lat-lon
    real (kind=real_kind), dimension(np,np)      :: E          ! kinetic energy term
    real (kind=real_kind), dimension(np,np)      :: zeta       ! relative vorticity
    real (kind=real_kind), dimension(np,np)      :: div, flux, plocal
    real (kind=real_kind), dimension(np,np,2)    :: ulatlon


    real (kind=real_kind) :: v1,v2,pstar,gmn,fjmax(4)
    real (kind=real_kind) :: vtens1,vtens2
    real (kind=real_kind) :: pmin(nlev,nets:nete),pmax(nlev,nets:nete)
    real (kind=real_kind) :: plmin(np,np,nlev,nets:nete),plmax(np,np,nlev,nets:nete)

    real (kind=real_kind) :: delta,real_time

    real*8     :: st,et,dtstage
    integer    :: i,j,k,s,ie,m,n
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: ntmp

    logical :: Debug = .FALSE.
    logical :: use_advective_flux 

    real (kind=real_kind) ::  notreliable


    if (test_case=="swtc1" .or. test_case=="vortex" .or. test_case=="swirl") then
       ! advection test cases support conservation form or advective form
       use_advective_flux=.true.
    else
       ! use shallow water flux
       use_advective_flux=.false.
       ! shallow water test cases require conservation form of h equation
       if (tracer_advection_formulation==TRACERADV_UGRADQ) then
          print *,'ERROR: shallow water tests require conservation formulation:'
          stop '(tracer_advection_formulation=1)'
       endif   
    endif

    nm1   = tl%nm1
    n0    = tl%n0
    np1   = tl%np1
    nstep = tl%nstep

    ! hyperviscosity can be applied in 3 different ways:
    ! 1. added to the RHS within each RK stage ( 1 extra DSS per stage)
    ! 2. time-split after each RK stage (2 extra DSS per stage)
    ! 3. time-split after complete RK step ( 2 extra DSS per step)
#define HYPERVIS_T1
#undef HYPERVIS_T2
#undef HYPERVIS_T3


    ! We want to make this leap-frog compliant
    ! Copy u^n to u^n+1
    do ie=nets,nete
       do k=1,nlev
          do j=1,np
             do i=1,np
                elem(ie)%state%v(i,j,1,k,np1)  = elem(ie)%state%v(i,j,1,k,n0)
                elem(ie)%state%v(i,j,2,k,np1)  = elem(ie)%state%v(i,j,2,k,n0)
                elem(ie)%state%p(i,j,k, np1)   = elem(ie)%state%p(i,j,k,n0)
             end do
          end do
       end do
    enddo
    real_time = dt*real(nstep,kind=real_kind)

    do s=1,MyRk%Stages
       dtstage = dt*(MyRk%beta(s)/MyRk%alpha(s))

       if(debug)then
          if(hybrid%par%masterproc) print *,"Performing stage ", s
       endif
       ! ===================================
       ! find min/max of p over all neighbors
       ! ===================================

!group optimal based on iteration
       if (( limiter_option == 8 ).or.(limiter_option == 81 )) then
	  if(kmass.ne.-1)then
	    call neighbor_minmax(elem,hybrid,edge3,nets,nete,n0,pmin,pmax,kmass=kmass)
	  else
	    call neighbor_minmax(elem,hybrid,edge3,nets,nete,n0,pmin,pmax)
	  endif
       endif


       if(Debug) print *,'homme: adv.._rk 1'


       ! ===================================
       ! construct v tendencies and v.grad(p)
       ! on the velocity grid...
       ! ===================================

#ifdef HYPERVIS_T1
       ! compute weak biharmonic operator.  has mass matrix built in,
       ! but to fit into structure below, remove it for now:

!applying viscosity to q field

	if(kmass.ne.-1)then
	do ie=nets,nete
	    do k=1,nlev
	      if(k.ne.kmass)then
		elem(ie)%state%p(:,:,k,n0)=elem(ie)%state%p(:,:,k,n0)/&
					    elem(ie)%state%p(:,:,kmass,n0)
	      endif
	    enddo
	enddo
	endif

        do ie=nets,nete
           do k=1,nlev
              ! contra -> latlon
              do j=1,np
                 do i=1,np
                    v1     = elem(ie)%state%v(i,j,1,k,n0)   ! contra
                    v2     = elem(ie)%state%v(i,j,2,k,n0)   ! contra 
                    elem(ie)%state%v(i,j,1,k,n0)=elem(ie)%D(1,1,i,j)*v1 + elem(ie)%D(1,2,i,j)*v2   ! contra->latlon
                    elem(ie)%state%v(i,j,2,k,n0)=elem(ie)%D(2,1,i,j)*v1 + elem(ie)%D(2,2,i,j)*v2   ! contra->latlon
                 enddo
              enddo
           enddo
        enddo
	call biharmonic_wk(elem,ptens,vtens,deriv,edge3,hybrid,n0,nets,nete)
        ! convert lat-lon -> contra variant
        do ie=nets,nete
           do k=1,nlev
              do j=1,np
                 do i=1,np
                    v1=elem(ie)%state%v(i,j,1,k,n0)
                    v2=elem(ie)%state%v(i,j,2,k,n0)
                    elem(ie)%state%v(i,j,1,k,n0) = elem(ie)%Dinv(1,1,i,j)*v1 + elem(ie)%Dinv(1,2,i,j)*v2
                    elem(ie)%state%v(i,j,2,k,n0) = elem(ie)%Dinv(2,1,i,j)*v1 + elem(ie)%Dinv(2,2,i,j)*v2
                 enddo
              enddo
           enddo
        enddo
        

	do ie=nets,nete
	    spheremp     => elem(ie)%spheremp
	    do k=1,nlev
	      ptens(:,:,k,ie) =  -nu_s*ptens(:,:,k,ie)/spheremp(:,:)
	      vtens(:,:,1,k,ie) = -nu*vtens(:,:,1,k,ie)/spheremp(:,:)
	      vtens(:,:,2,k,ie) = -nu*vtens(:,:,2,k,ie)/spheremp(:,:)
	    enddo
	enddo

	if(kmass.ne.-1)then
	!we do not apply viscosity to mass field
	ptens(:,:,kmass,:)=0.0d0
	do ie=nets,nete
	    do k=1,nlev
	      if(k.ne.kmass)then
		elem(ie)%state%p(:,:,k,n0)=elem(ie)%state%p(:,:,k,n0)*&
					    elem(ie)%state%p(:,:,kmass,n0)
	      endif
	    enddo
	enddo
	endif
#else
       ptens=0
       vtens=0
! ---
#endif
! --- endif  HYPERVIS_T1

       do ie=nets,nete
          met    => elem(ie)%met
          metinv => elem(ie)%metinv
          metdet => elem(ie)%metdet
          rmetdet => elem(ie)%rmetdet
          fcor   => elem(ie)%fcor
          spheremp     => elem(ie)%spheremp

          call set_prescribed_velocity(elem(ie),n0,real_time)

          do k=1,nlev
             ! ==============================================
             ! Compute kinetic energy term
             ! ==============================================
             do j=1,np
                do i=1,np

                   v1     = elem(ie)%state%v(i,j,1,k,n0)   ! contra
                   v2     = elem(ie)%state%v(i,j,2,k,n0)   ! contra 
                   ulatlon(i,j,1)=elem(ie)%D(1,1,i,j)*v1 + elem(ie)%D(1,2,i,j)*v2   ! contra->latlon
                   ulatlon(i,j,2)=elem(ie)%D(2,1,i,j)*v1 + elem(ie)%D(2,2,i,j)*v2   ! contra->latlon

                   E(i,j) = 0.5D0*(ulatlon(i,j,1)**2 + ulatlon(i,j,2)**2)  +&
                        elem(ie)%state%p(i,j,k,n0) + elem(ie)%state%ps(i,j)

                   pv(i,j,1) = ulatlon(i,j,1)*(pmean+elem(ie)%state%p(i,j,k,n0))
                   pv(i,j,2) = ulatlon(i,j,2)*(pmean+elem(ie)%state%p(i,j,k,n0))
                end do
             end do
             grade = gradient_sphere(E,deriv,elem(ie)%Dinv)       ! scalar -> latlon vector
             !grade = gradient_sphere_wk(E,deriv,elem(ie)%Dinv)       ! scalar -> latlon vector
             zeta = vorticity_sphere(ulatlon,deriv,elem(ie)) ! latlon vector -> scalar 
             if (tracer_advection_formulation==TRACERADV_UGRADQ) then
                gradh = gradient_sphere(elem(ie)%state%p(:,:,k,n0),deriv,elem(ie)%Dinv)
                div = ulatlon(:,:,1)*gradh(:,:,1)+ulatlon(:,:,2)*gradh(:,:,2)
             else
                div = divergence_sphere(pv,deriv,elem(ie))      ! latlon vector -> scalar
             endif
             if (npdg>0) then
                ptens_dg(:,:,k,ie) = divergence_sphere_wk(pv,deriv,elem(ie))      ! latlon vector -> scalar
             endif

             ! ==============================================
             ! Compute velocity tendency terms
             ! ==============================================
             ! accumulate all RHS terms
             vtens(:,:,1,k,ie)=vtens(:,:,1,k,ie) + (ulatlon(:,:,2)*(fcor(:,:) + zeta(:,:))  - grade(:,:,1))
             vtens(:,:,2,k,ie)=vtens(:,:,2,k,ie) + (-ulatlon(:,:,1)*(fcor(:,:) + zeta(:,:)) - grade(:,:,2))
             ptens(:,:,k,ie) = ptens(:,:,k,ie) - div(:,:)

             ! take the local element timestep
             vtens(:,:,:,k,ie)=ulatlon(:,:,:) + dtstage*vtens(:,:,:,k,ie)
             ptens(:,:,k,ie) = elem(ie)%state%p(:,:,k,n0) + dtstage*ptens(:,:,k,ie)
          end do!end of loop over levels

          
          if ((limiter_option == 8))then
             call limiter_optim_iter_full4(ptens(:,:,:,ie),elem(ie)%spheremp(:,:),&
                  pmin(:,ie),pmax(:,ie),kmass)
!             call limiter_optim_iter_full3(ptens(:,:,:,ie),elem(ie)%spheremp(:,:),&
!                  pmin(:,ie),pmax(:,ie),kmass,.true.,.true.)
          endif

          if ((limiter_option == 81))then
             call limiter_optim_iter_full(ptens(:,:,:,ie),elem(ie)%spheremp(:,:),&
                  pmin(:,ie),pmax(:,ie),kmass,.true.,.false.,&
                  notreliable)
          endif
          
          if ((limiter_option == 84))then
	     pmin(:,ie)=0.0d0
             if (test_case=='swirl') then
                 pmin(1,ie)=.1d0
                 if (nlev>=3) then
                    k=3; pmin(k,ie)=.1d0
                 endif
             endif
 	     call limiter_optim_iter_full(ptens(:,:,:,ie),elem(ie)%spheremp(:,:),&
               pmin(:,ie),pmax(:,ie),kmass,.true.,.false.,&
		notreliable)
          endif

          if ( (limiter_option == 4) ) then
             call limiter2d_zero(ptens(:,:,:,ie),elem(ie)%spheremp, kmass)
          endif
          
          do k=1,nlev
             ptens(:,:,k,ie) = ptens(:,:,k,ie)*elem(ie)%spheremp(:,:)
             vtens(:,:,1,k,ie) = vtens(:,:,1,k,ie)*elem(ie)%spheremp(:,:)
             vtens(:,:,2,k,ie) = vtens(:,:,2,k,ie)*elem(ie)%spheremp(:,:)
          enddo
          
          ! ===================================================
          ! Pack cube edges of tendencies, rotate velocities
          ! ===================================================
          if (npdg>0) then
             ! for hybrid cg/dg element, pack p for flux calculation below
             do k=1,nlev
                ptens(:,:,k,ie)=elem(ie)%state%p(:,:,k,n0)
             enddo
          endif
          kptr=0
          call edgeVpack(edge3, ptens(1,1,1,ie),nlev,kptr,elem(ie)%desc)
          kptr=nlev
          call edgeVpack(edge3,vtens(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)

       end do

       if(Debug) print *,'homme: adv.._rk 2'

#if (! defined VERT_OPENMP)
       !$OMP BARRIER
#endif

       call bndry_exchangeV(hybrid,edge3)

#if (! defined VERT_OPENMP)
       !$OMP BARRIER
#endif
       
       do ie=nets,nete

          rspheremp     => elem(ie)%rspheremp

          ! ===========================================================
          ! Unpack the edges for vgradp and vtens
          ! ===========================================================
          if (npdg==0) then
             kptr=0
             call edgeVunpack(edge3, ptens(1,1,1,ie), nlev, kptr, elem(ie)%desc)
          endif
          kptr=nlev
          call edgeVunpack(edge3, vtens(1,1,1,1,ie), 2*nlev, kptr, elem(ie)%desc)

          ! ===========================================================
          ! Compute velocity and pressure tendencies for all levels
          ! ===========================================================
          do k=1,nlev
             do j=1,np
                do i=1,np
                   ptens(i,j,k,ie) = rspheremp(i,j)*ptens(i,j,k,ie)
                   vtens1=rspheremp(i,j)*vtens(i,j,1,k,ie)
                   vtens2=rspheremp(i,j)*vtens(i,j,2,k,ie)

                   ! lat-lon -> contra
                   vtens(i,j,1,k,ie) = elem(ie)%Dinv(1,1,i,j)*vtens1 + elem(ie)%Dinv(1,2,i,j)*vtens2
                   vtens(i,j,2,k,ie) = elem(ie)%Dinv(2,1,i,j)*vtens1 + elem(ie)%Dinv(2,2,i,j)*vtens2
                end do
             end do
          end do

          if (npdg>0) then
             kptr=0
             call edgeDGVunpack(edge3, pedges, nlev, kptr, elem(ie)%desc)
             pedges=pedges+pmean  ! add in mean value, to get edge flux correct
             do k=1,nlev
                plocal(:,:)=elem(ie)%state%p(:,:,k,n0)+pmean  ! add in mean value
                if (use_advective_flux) then
                   ! simple upwind flux
                   !flux=edge_flux_u_cg( elem(ie)%state%v(:,:,:,k,n0), elem(ie)%state%p(:,:,k,n0),&
                   !     pedges(:,:,k), deriv, elem(ie), u_is_contra=.true.)
                   flux=adv_flux_term(elem(ie),deriv,elem(ie)%state%v(:,:,:,k,n0),plocal,pedges(:,:,k))
                else
                   ! shallow water flux:
                   fjmax=sw_fjmax(elem(ie)%state%v(:,:,:,k,n0),plocal,pedges(:,:,k),elem(ie))
                   call swsys_flux(1,elem(ie),deriv,fjmax,plocal,pedges(:,:,k),elem(ie)%state%v(:,:,:,k,n0),flux)
                endif

                ! combine weak gradient and edge flux:
                ptens(:,:,k,ie) = ptens_dg(:,:,k,ie) + flux(:,:)

                ! advance in time. GLL quadrature, cardinal function basis, under-integrated.  
                ! local mass matrix is diagonal, with entries elem(ie)%spheremp(),
                ! so we divide through by elem(ie)%spheremp().
                ptens(:,:,k,ie) = elem(ie)%state%p(:,:,k,n0) - dtstage*ptens(:,:,k,ie)/elem(ie)%spheremp(:,:)
                if (npdg<np) then
                   ! modal timestep, with exact integration.  using prognostic variable: p*metdet
                   ! local mass matrix is diagonal assuming npdg<np so that GLL quadrature is exact)
                   ! (note: GLL/modal conversion comutes with time-stepping)

                   ! compute modal coefficients of p*metdet
                   ! (spherical inner-product of Legendre polynomial and p)
                   phat = gll_to_dgmodal(ptens(:,:,k,ie)*elem(ie)%metdet(:,:),deriv)   

                   ! modal based limiter goes here

                   ! evalute modal expanion of p*metdet on GLL points
                   ptens(:,:,k,ie)=dgmodal_to_gll(phat,deriv) 
                   
                   ! convert back to p
                   ptens(:,:,k,ie)=ptens(:,:,k,ie)*elem(ie)%rmetdet(:,:)
                endif
             enddo
             ! truncation + mass weighted redistribution:
             if (limiter_option==4) call limiter2d_zero(ptens(:,:,:,ie),elem(ie)%spheremp, kmass)
             ! Find optimal (l2 norm) solution which is closest to unlimited solution:
             if (limiter_option==8) call limiter_optim_iter_full4(ptens(:,:,:,ie),elem(ie)%spheremp(:,:),&
                  pmin(:,ie),pmax(:,ie),kmass)
             !if (limiter_option==8) call limiter_optim_iter_full3(ptens(:,:,:,ie),elem(ie)%spheremp(:,:),&
              !    pmin(:,ie),pmax(:,ie),kmass,.true.,.true.)
          endif
          do k=1,nlev
             ! ====================================================
             ! average different timelevels for RK-SSP
             ! ====================================================
             do j=1,np
                do i=1,np
                   elem(ie)%state%v(i,j,1,k,n0) = MyRk%alpha0(s)*elem(ie)%state%v(i,j,1,k,np1) &
                        + MyRk%alpha(s)*vtens(i,j,1,k,ie)
                   elem(ie)%state%v(i,j,2,k,n0) = MyRk%alpha0(s)*elem(ie)%state%v(i,j,2,k,np1) &
                        + MyRk%alpha(s)*vtens(i,j,2,k,ie)
                   elem(ie)%state%p(i,j,k,n0) = MyRk%alpha0(s)*elem(ie)%state%p(i,j,k,np1)+ &
                        MyRk%alpha(s)*ptens(i,j,k,ie)
                end do
             end do
          end do
       end do



#ifdef HYPERVIS_T2
       call advance_hypervis(edge3,elem,hybrid,deriv,vtens,ptens,n0,&
            nets,nete,dtstage) 
! ---
#endif
! --- endif  HYPERVIS_T2

       real_time =real_time + dtstage

!this is only for output reasons, if velocities are prescribed
      do ie=nets,nete
	  call set_prescribed_velocity(elem(ie),n0,real_time)
      enddo


    enddo ! stage loop

    ! ====================================================
    ! apply viscosity  
    ! ====================================================
#ifdef HYPERVIS_T3
    call advance_hypervis(edge3,elem,hybrid,deriv,vtens,ptens,n0,&
         nets,nete,dt) 
! ---
#endif
! --- endif HYPERVIS_T3

    ! make leap-frog compliant
    ! after RK loop u^n containe u^{n+1} and u^{n+1} contains u^n
    ! thus we swap them
    ntmp    = tl%np1
    tl%np1  = tl%n0
    tl%n0   = ntmp

#if (! defined VERT_OPENMP)
    !$OMP BARRIER
#endif
  end subroutine advance_nonstag_rk

!-----------------------------------------------------------------------------------------------

  subroutine limiter_optim_iter_full4(ptens,sphweights,minp,maxp,kmass)
!THIS IS A NEW VERSION OF LIM8, from 3D code

    use kinds, only : real_kind
    use dimensions_mod, only : np, np, nlev

    real (kind=real_kind), dimension(nlev), intent(inout)   :: minp, maxp
    real (kind=real_kind), dimension(np,np,nlev), intent(inout)   :: ptens
    real (kind=real_kind), dimension(np,np), intent(in)   :: sphweights
    integer, intent(in) :: kmass

    real (kind=real_kind), dimension(np,np) :: weights
    real (kind=real_kind), dimension(np,np) :: ptens_mass
    integer  k1, k, i, j, iter, i1, i2
    integer :: whois_neg(np*np), whois_pos(np*np), neg_counter, pos_counter

    real (kind=real_kind) :: addmass, weightssum, mass
    real (kind=real_kind) :: x(np*np),c(np*np)

    real (kind=real_kind) :: al_neg(np*np), al_pos(np*np), howmuch
    real (kind=real_kind) :: tol_limiter=1e-15

    integer, parameter :: maxiter = 5

    weights=sphweights

!if kmass is valid number, we first get (\rho Q)/(\rho) fields 
    if(kmass.ne.-1)then   
	ptens_mass(:,:)=ptens(:,:,kmass)
	weights=weights*ptens_mass
	do k=1,nlev
	  if(k.ne.kmass)then
	    ptens(:,:,k)=ptens(:,:,k)/ptens_mass
	  endif
	enddo
    endif

    do k=1,nlev

      if(k.ne.kmass)then

      k1=1
      do i=1,np
	do j=1,np
	  c(k1)=weights(i,j)
	  x(k1)=ptens(i,j,k)
	  k1=k1+1
	enddo
      enddo

      mass=sum(c*x)

      ! relax constraints to ensure limiter has a solution:
      ! This is only needed if runnign with the SSP CFL>1 or 
      ! due to roundoff errors
      if( (mass/sum(c))<minp(k) ) then
         minp(k)=mass/sum(c)
      endif
      if( (mass/sum(c))>maxp(k) ) then
         maxp(k)=mass/sum(c)
      endif

      addmass=0.0d0
      pos_counter=0;
      neg_counter=0;
      
      ! apply constraints, compute change in mass caused by constraints 
      do k1=1,np*np
         if((x(k1)>=maxp(k)) ) then
	    addmass=addmass+(x(k1)-maxp(k))*c(k1)
	    x(k1)=maxp(k)
	    whois_pos(k1)=-1
         else
            pos_counter=pos_counter+1;
	    whois_pos(pos_counter)=k1;
         endif
         if((x(k1)<=minp(k)) ) then
	    addmass=addmass-(minp(k)-x(k1))*c(k1)
	    x(k1)=minp(k)
	    whois_neg(k1)=-1
         else
	    neg_counter=neg_counter+1;
	    whois_neg(neg_counter)=k1;
         endif
      enddo
      
      ! iterate to find field that satifies constraints and is l2-norm closest to original 
      weightssum=0.0d0
      if(addmass>0)then
         do i2=1,maxIter
            weightssum=0.0
            do k1=1,pos_counter
               i1=whois_pos(k1)
               weightssum=weightssum+c(i1)
               al_pos(i1)=maxp(k)-x(i1)
            enddo
            
            if((pos_counter>0).and.(addmass>tol_limiter*abs(mass)))then
               do k1=1,pos_counter
                  i1=whois_pos(k1)
                  howmuch=addmass/weightssum
                  if(howmuch>al_pos(i1))then
                     howmuch=al_pos(i1)
                     whois_pos(k1)=-1
                  endif
                  addmass=addmass-howmuch*c(i1)
                  weightssum=weightssum-c(i1)
                  x(i1)=x(i1)+howmuch
               enddo
               !now sort whois_pos and get a new number for pos_counter
               !here neg_counter and whois_neg serve as temp vars
               neg_counter=pos_counter
               whois_neg=whois_pos
               whois_pos=-1
               pos_counter=0
               do k1=1,neg_counter
                  if(whois_neg(k1).ne.-1)then
                     pos_counter=pos_counter+1
                     whois_pos(pos_counter)=whois_neg(k1)
                  endif
               enddo
            else
               exit
            endif
         enddo
      else
         do i2=1,maxIter
            weightssum=0.0
            do k1=1,neg_counter
               i1=whois_neg(k1)
               weightssum=weightssum+c(i1)
               al_neg(i1)=x(i1)-minp(k)
            enddo
            
            if((neg_counter>0).and.((-addmass)>tol_limiter*abs(mass)))then
               
               do k1=1,neg_counter
                  i1=whois_neg(k1)
                  howmuch=-addmass/weightssum
                  if(howmuch>al_neg(i1))then
                     howmuch=al_neg(i1)
                     whois_neg(k1)=-1
                  endif
                  addmass=addmass+howmuch*c(i1)
                  weightssum=weightssum-c(i1)
                  x(i1)=x(i1)-howmuch
               enddo
               !now sort whois_pos and get a new number for pos_counter
               !here pos_counter and whois_pos serve as temp vars
               pos_counter=neg_counter
               whois_pos=whois_neg
               whois_neg=-1
               neg_counter=0
               do k1=1,pos_counter
                  if(whois_pos(k1).ne.-1)then
                     neg_counter=neg_counter+1
                     whois_neg(neg_counter)=whois_pos(k1)
                  endif
               enddo
            else
               exit
            endif
            
         enddo
         
      endif
      
      k1=1
      do i=1,np
         do j=1,np
            ptens(i,j,k)=x(k1)
            k1=k1+1
         enddo
      enddo
      
      endif !if k is not kmass
   enddo
   
    if(kmass.ne.-1)then   
	do k=1,nlev
	  if(k.ne.kmass)then
	    ptens(:,:,k)=ptens(:,:,k)*ptens_mass
	  endif
	enddo
    endif

  end subroutine limiter_optim_iter_full4


!----------------------------------------------------------------------------
  subroutine limiter_optim_iter_full3(ptens,sphweights,minp,maxp,kmass,checkmin,checkmax)
!The idea here is the following: We need to find a grid field which is closest
!to the initial field (in terms of a weighted sum), but satisfies the constraints.
!So, first we find values which do not satisfy constraints and bring these values
!to a closest constraint. This way we introduce some mass change (addmass),
!so, we redistribute addmass in the way that error is smallest. This redistribution might !violate constraints (though I think the solution is given by one iteration only if the !problem is well-posed) due to round off, for example; thus, we do a few iterations. 

    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev
    use control_mod, only : tol_limiter

    logical, intent(in) :: checkmin,checkmax
    real (kind=real_kind), dimension(nlev), intent(inout)   :: minp
    real (kind=real_kind), dimension(nlev), intent(inout)   :: maxp
    real (kind=real_kind), dimension(np,np,nlev), intent(inout)   :: ptens
    real (kind=real_kind), dimension(np,np), intent(in)   :: sphweights
    integer, intent(in) :: kmass
 
    real (kind=real_kind), dimension(np,np) :: weights
    real (kind=real_kind), dimension(np,np) :: ptens_mass
    integer  k1, k, i, j, iter, i1, i2
    integer :: pos_counter, neg_counter, whois_neg(np*np), whois_pos(np*np)
    real (kind=real_kind) :: addmass, weightssum, mass
    real (kind=real_kind) :: x(np*np),c(np*np)
    real (kind=real_kind) :: al_neg(np*np), al_pos(np*np), howmuch

    integer, parameter :: maxIter=10

    weights=sphweights

!!!!this is a fix for huge tracers, otherwise we'd expect undershoots do not accumulate badly
    if(minp(1)<0) minp(1)=0.0d0

!if kmass is valid number, we first get (\rho Q)/(\rho) fields 
    if(kmass.ne.-1)then   
	ptens_mass(:,:)=ptens(:,:,kmass)
	weights=weights*ptens_mass
	do k=1,nlev
	  if(k.ne.kmass)then
	    ptens(:,:,k)=ptens(:,:,k)/ptens_mass
	  endif
	enddo
    endif

    do k=1,nlev

      if(k.ne.kmass)then
	k1=1
	do i=1,np
	  do j=1,np
	    c(k1)=weights(i,j)
	    x(k1)=ptens(i,j,k)
	    k1=k1+1
	  enddo
	enddo

        mass=sum(c*x)

        if((mass-minp(k)*sum(c)<0.0).and.checkmin)then
	   minp(k)=mass/sum(c)
        endif
        if((mass-maxp(k)*sum(c)>0.0).and.checkmax)then
	   maxp(k)=mass/sum(c)
        endif

	addmass=0.0d0
	pos_counter=0;
        neg_counter=0;

	do k1=1,np*np
	  if((x(k1)>=maxp(k)).AND.(checkmax))then
	    addmass=addmass+(x(k1)-maxp(k))*c(k1)
	    x(k1)=maxp(k)
	    whois_pos(k1)=-1
          else
            pos_counter=pos_counter+1;
	    whois_pos(pos_counter)=k1;
	  endif
	  if((x(k1)<=minp(k)).AND.(checkmin))then
	    addmass=addmass-(minp(k)-x(k1))*c(k1)
	    x(k1)=minp(k)
	    whois_neg(k1)=-1
	  else
	    neg_counter=neg_counter+1;
	    whois_neg(neg_counter)=k1;
	  endif
	enddo

	weightssum=0.0d0

	if(addmass>0)then

	    do i2=1,maxIter

		weightssum=0.0
		do k1=1,pos_counter
		    i1=whois_pos(k1)
		    weightssum=weightssum+c(i1)
		    al_pos(i1)=maxp(k)-x(i1)
		enddo

!		if((pos_counter>0).and.(addmass>1e-20))then
		if((pos_counter>0).and.(addmass>tol_limiter*abs(mass)))then

		  do k1=1,pos_counter
		    i1=whois_pos(k1)
		    howmuch=addmass/weightssum
		    if(howmuch>al_pos(i1))then
			howmuch=al_pos(i1)
			whois_pos(k1)=-1
		    endif
		    addmass=addmass-howmuch*c(i1)
		    weightssum=weightssum-c(i1)
		    x(i1)=x(i1)+howmuch
		  enddo
		  !now sort whois_pos and get a new number for pos_counter
		  !here neg_counter and whois_neg serve as temp vars
		  neg_counter=pos_counter
		  whois_neg=whois_pos
		  whois_pos=-1
		  pos_counter=0
		  do k1=1,neg_counter
		    if(whois_neg(k1).ne.-1)then
		      pos_counter=pos_counter+1
		      whois_pos(pos_counter)=whois_neg(k1)
		    endif
		  enddo
		else
!			  if (abs(sum(c*x))>0) x=x*(mass+addmass)/sum(c*x)
!			  x=x+addmass/sum(c)
		  exit
		endif
	    enddo
	else
	    do i2=1,maxIter

		weightssum=0.0
		do k1=1,neg_counter
		    i1=whois_neg(k1)
		    weightssum=weightssum+c(i1)
		    al_neg(i1)=x(i1)-minp(k)
		enddo

!		if((neg_counter>0).and.((-addmass)>1e-20))then
		if((neg_counter>0).and.((-addmass)>tol_limiter*abs(mass)))then

		  do k1=1,neg_counter
		    i1=whois_neg(k1)
		    howmuch=-addmass/weightssum
		    if(howmuch>al_neg(i1))then
			howmuch=al_neg(i1)
			whois_neg(k1)=-1
		    endif
		    addmass=addmass+howmuch*c(i1)
		    weightssum=weightssum-c(i1)
		    x(i1)=x(i1)-howmuch
		  enddo
		  !now sort whois_pos and get a new number for pos_counter
		  !here pos_counter and whois_pos serve as temp vars
		  pos_counter=neg_counter
		  whois_pos=whois_neg
		  whois_neg=-1
		  neg_counter=0
		  do k1=1,pos_counter
		    if(whois_pos(k1).ne.-1)then
		      neg_counter=neg_counter+1
		      whois_neg(neg_counter)=whois_pos(k1)
		    endif
		  enddo
		else
!			 if (abs(sum(c*x))>0) x=x*(mass+addmass)/sum(c*x)
!			 x=x+addmass/sum(c)
		  exit
		endif
	    enddo
	endif


	k1=1
	do i=1,np
	  do j=1,np
	    ptens(i,j,k)=x(k1)
	    k1=k1+1
	  enddo
	enddo

      endif!if k is not kmass
    enddo

    if(kmass.ne.-1)then   
	do k=1,nlev
	  if(k.ne.kmass)then
	    ptens(:,:,k)=ptens(:,:,k)*ptens_mass
	  endif
	enddo
    endif

  end subroutine limiter_optim_iter_full3



!--------------------------------------------------------------------------


  subroutine limiter_optim_iter_full(ptens,sphweights,minp,maxp,kmass,checkmin,checkmax, notreliable)
!The idea here is the following: We need to find a grid field which is closest
!to the initial field (in terms of a weighted sum), but satisfies the constraints.
!So, first we find values which do not satisfy constraints and bring these values
!to a closest constraint. This way we introduce some mass change (addmass),
!so, we redistribute addmass in the way that error is smallest. This redistribution might !violate constraints (though I think the solution is given by one iteration only if the !problem is well-posed) due to round off, for example; thus, we do a few iterations. 

    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev
    use control_mod, only : tol_limiter

    logical, intent(in) :: checkmin,checkmax
    real (kind=real_kind), dimension(nlev), intent(inout)   :: minp
    real (kind=real_kind), dimension(nlev), intent(inout)   :: maxp
    real (kind=real_kind), dimension(np,np,nlev), intent(inout)   :: ptens
    real (kind=real_kind), dimension(np,np), intent(in)   :: sphweights
    integer, intent(in) :: kmass
    real (kind=real_kind),  intent(out) :: notreliable
 
    real (kind=real_kind), dimension(np,np) :: weights
    real (kind=real_kind), dimension(np,np) :: ptens_mass
    integer  k1, k, i, j, iter
    logical :: whois_neg(np*np), whois_pos(np*np)
    real (kind=real_kind) :: addmass, weightssum, mass
    real (kind=real_kind) :: x(np*np),c(np*np)

    integer, parameter :: maxiter = 10

    weights=sphweights
    notreliable=0

!if kmass is valid number, we first get (\rho Q)/(\rho) fields 
    if(kmass.ne.-1)then   
	ptens_mass(:,:)=ptens(:,:,kmass)
	weights=weights*ptens_mass
	do k=1,nlev
	  if(k.ne.kmass)then
	    ptens(:,:,k)=ptens(:,:,k)/ptens_mass
	  endif
	enddo
    endif

    do k=1,nlev

      if(k.ne.kmass)then
	k1=1
	do i=1,np
	  do j=1,np
	    c(k1)=weights(i,j)
	    x(k1)=ptens(i,j,k)
	    k1=k1+1
	  enddo
	enddo

        mass=sum(c*x)

!first, check if the problem has a solution
!if it seens that it does not, we initialize a flag, 
!but the limiter still might do a good job
        if((mass-minp(k)*sum(c)<-tol_limiter).and.checkmin)then
           notreliable=1
	   minp(k)=mass/sum(c)
        endif
        if((mass-maxp(k)*sum(c)>tol_limiter).and.checkmax)then
           notreliable=1
	   maxp(k)=mass/sum(c)
        endif

        iter=0
        do 
          iter=iter+1

          if(iter>maxiter)then
	    notreliable=1.0
            exit
          endif

	  addmass=0.0d0

	  do k1=1,np*np
	    whois_neg(k1)=.true.
	    whois_pos(k1)=.true.
	    if((x(k1)>=maxp(k)).AND.(checkmax))then
	      addmass=addmass+(x(k1)-maxp(k))*c(k1)
	      x(k1)=maxp(k)
	      whois_pos(k1)=.false.
	    endif
	    if((x(k1)<=minp(k)).AND.(checkmin))then
	      addmass=addmass-(minp(k)-x(k1))*c(k1)
	      x(k1)=minp(k)
	      whois_neg(k1)=.false.
	    endif
	  enddo

	  if(abs(addmass)>tol_limiter*abs(mass))then

	    weightssum=0.0d0
	    if(addmass>0)then
	      do k1=1,np*np
		if(whois_pos(k1))then
		  weightssum=weightssum+c(k1)
		endif
	      enddo
              if(weightssum>0.0)then
		do k1=1,np*np
		  if(whois_pos(k1))then
		    x(k1)=x(k1)+addmass/weightssum
		  endif
		enddo
              else
	        x=x+addmass/sum(c)
		notreliable=1
		exit
	      endif
            else
	      do k1=1,np*np
		if(whois_neg(k1))then
		  weightssum=weightssum+c(k1)
		endif
	      enddo
              if(weightssum>0.0)then
		do k1=1,np*np
		  if(whois_neg(k1))then
		    x(k1)=x(k1)+addmass/weightssum
		  endif
		enddo
              else
	        x=x+addmass/sum(c)
		notreliable=1
		exit
	      endif
            endif

	  else
            x=x+addmass/sum(c)
	    exit
	  endif

	enddo!end of iteration

	k1=1
	do i=1,np
	  do j=1,np
	    ptens(i,j,k)=x(k1)
	    k1=k1+1
	  enddo
	enddo


      endif
    enddo


    if(kmass.ne.-1)then   
	do k=1,nlev
	  if(k.ne.kmass)then
	    ptens(:,:,k)=ptens(:,:,k)*ptens_mass
	  endif
	enddo
    endif


  end subroutine limiter_optim_iter_full
!--------------------------------------------------------------------------


  subroutine limiter2d_zero(Q,spheremp,kmass)
    !
    ! mass conserving sign-preserving limiter (2D only). 
    ! uses specified global minimum
    ! 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev
    use control_mod, only :  test_case

    implicit none
    real (kind=real_kind), intent(inout)  :: Q(np,np,nlev)
    real (kind=real_kind), intent(in)  :: spheremp(np,np)
    integer, intent(in) :: kmass

    ! local
    integer k
    real (kind=real_kind) :: qmin

    if(kmass.ne.-1)then
      do k=nlev,1,-1
	qmin=0
	if (test_case=='swirl') then
	    if (k.eq.1) qmin=0.0d0   ! lifted cosine bell
	    if (k.eq.3) qmin=0.1d0   ! lifted slotted cylinder
        endif
	if(k.ne.kmass)then
	  Q(:,:,k)=Q(:,:,k)/Q(:,:,kmass)
	  call limiter2d_min_onelevel(Q(:,:,k),spheremp(:,:)*Q(:,:,kmass),qmin)
	  Q(:,:,k)=Q(:,:,k)*Q(:,:,kmass)
	endif
      enddo
    else
      do k=nlev,1,-1
	qmin=0
	if (test_case=='swirl') then
	    if (k.eq.1) qmin=0.1d0   ! lifted cosine bell
	    if (k.eq.3) qmin=0.1d0   ! lifted slotted cylinder
	endif
	call limiter2d_min_onelevel(Q(:,:,k),spheremp(:,:),qmin)
      enddo
    endif
  end subroutine limiter2d_zero


!-------------------------------------------------------------------------------
  subroutine limiter2d_min_onelevel(Q,spheremp,qmin)
    use kinds, only : real_kind
    use dimensions_mod, only : np

    implicit none
    real (kind=real_kind), intent(inout)  :: Q(np,np)
    real (kind=real_kind), intent(in)  :: spheremp(np,np)
    real (kind=real_kind), intent(in)  :: qmin

    ! local
    integer i,j
    real (kind=real_kind) :: mass,mass_added,area,mass2

    if ( minval(Q(:,:)) < qmin ) then
      mass=sum( Q(:,:)*spheremp(:,:) )
      area=sum( spheremp(:,:) )
      mass2 = mass - area*qmin
      Q(:,:)=Q(:,:)-qmin

      ! negative mass.  so reduce all postive values to zero 
      ! then increase negative values as much as possible
      if (mass2 < 0) Q(:,:)=-Q(:,:) 
      mass_added=0
      do j=1,np	
	  do i=1,np
	    if (Q(i,j)<0) then
		Q(i,j)=0
	    else
		mass_added = mass_added + Q(i,j)*spheremp(i,j)
	    endif
	  enddo
      enddo
      ! now scale the all positive values to restore mass
      if (mass_added>0) Q(:,:) = Q(:,:)*abs(mass2)/mass_added
      if (mass2 < 0) Q(:,:)=-Q(:,:)         
      Q(:,:)=Q(:,:)+qmin
    endif

  end subroutine limiter2d_min_onelevel

!-------------------------------------------------------------------------------
  subroutine limiter2d_max_onelevel(Q,spheremp,qmax)
  
    use kinds, only : real_kind
    use dimensions_mod, only : np

    implicit none
    real (kind=real_kind), intent(inout)  :: Q(np,np)
    real (kind=real_kind), intent(in)  :: spheremp(np,np)
    real (kind=real_kind), intent(in)  :: qmax

    ! local
    integer i,j
    real (kind=real_kind) :: mass,mass_added,area,mass2


    ! max limiter
    if ( maxval(Q(:,:)) > qmax ) then
       mass=sum( Q(:,:)*spheremp(:,:) )
       area=sum( spheremp(:,:) )
       mass2 = area*qmax - mass
       
      Q(:,:)=qmax-Q(:,:)

      if (mass2 < 0) Q(:,:)=-Q(:,:) 
      mass_added=0
      do j=1,np	
	  do i=1,np
	    if (Q(i,j)<0) then
		Q(i,j)=0
	    else
		mass_added = mass_added + Q(i,j)*spheremp(i,j)
	    endif
	  enddo
      enddo
      ! now scale the all positive values to restore mass
      if (mass_added>0) Q(:,:) = Q(:,:)*abs(mass2)/mass_added
      if (mass2 < 0) Q(:,:)=-Q(:,:) 
      Q(:,:)=qmax-Q(:,:)
    endif

  end subroutine limiter2d_max_onelevel

!-------------------------------------------------------------------------------------




  subroutine advance_hypervis(edge3,elem,hybrid,deriv,vtens,ptens,nt,nets,nete,dt2)
    !
    !  take one timestep of:  
    !          u(:,:,:,np) = u(:,:,:,np) +  dt2*nu*laplacian**order ( u )
    !          h(:,:,:,np) = h(:,:,:,np) +  dt2*nu_s*laplacian**order ( h )
    !
    !   For height equation, we want to dissipate P.E. = gH^2, where H = surface height = h + h_surface
    !   since dH/dt = dh/dt, we can solve:
    !          h(:,:,:,np) = h(:,:,:,np) +  dt2*nu_s*laplacian**order ( H )
    !
    !   (to understand this, think of a flow at rest with topography.  H=constant,
    !    laplace(H)=0.  But h= -h_surface and laplacian(h) <> 0
    !
    !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
    !
    !
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev
    use control_mod, only : nu, nu_div, nu_s, hypervis_order, hypervis_subcycle, limiter_option,&
	  test_case, kmass
    use hybrid_mod, only : hybrid_t
    use element_mod, only : element_t
    use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack
    use bndry_mod, only : bndry_exchangev
    use viscosity_mod, only : biharmonic_wk, neighbor_minmax
    ! ---------------------
    implicit none

    type (hybrid_t)      , intent(in) :: hybrid
    type (element_t)     , intent(inout), target :: elem(:)
    real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)  :: vtens
    real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: ptens
    type (EdgeBuffer_t)  , intent(inout) :: edge3
    type (derivative_t)  , intent(in) :: deriv
    real (kind=real_kind) :: dt2
    integer :: nt,nets,nete


    ! local
    integer :: k,kptr,i,j,ie,ic
    real (kind=real_kind), dimension(:,:), pointer :: spheremp,rspheremp
    real (kind=real_kind), dimension(np,np) :: ptot
    real (kind=real_kind), dimension(np,np) :: lap_p
    real (kind=real_kind), dimension(np,np,2) :: lap_v
    real (kind=real_kind) :: v1,v2,dt
    real (kind=real_kind) :: pmin(nlev,nets:nete),pmax(nlev,nets:nete)
    real (kind=real_kind) :: notreliable

    if (hypervis_order == 0 ) return;
    if (nu_s == 0 .and. nu == 0 ) return;

!group of lim3 limiters, redistribution
    if ( (limiter_option == 8 ).or.( limiter_option == 81 )) then
      if(kmass.ne.-1)then
	call neighbor_minmax(elem,hybrid,edge3,nets,nete,nt,pmin,pmax,kmass=kmass)
      else
	call neighbor_minmax(elem,hybrid,edge3,nets,nete,nt,pmin,pmax)
      endif
    endif



    do ie=nets,nete
       do k=1,nlev

          ! contra -> latlon
          do j=1,np
             do i=1,np
                v1     = elem(ie)%state%v(i,j,1,k,nt)   ! contra
                v2     = elem(ie)%state%v(i,j,2,k,nt)   ! contra 
                elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%D(1,1,i,j)*v1 + elem(ie)%D(1,2,i,j)*v2   ! contra->latlon
                elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%D(2,1,i,j)*v1 + elem(ie)%D(2,2,i,j)*v2   ! contra->latlon
             enddo
          enddo
       enddo
    enddo


    dt=dt2/hypervis_subcycle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  regular viscosity  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (hypervis_order == 1) then

       do ic=1,hypervis_subcycle
          do ie=nets,nete

             spheremp     => elem(ie)%spheremp

             do k=1,nlev
                ! filter surface height, not thickness
                do j=1,np
                   do i=1,np             
                      ptot(i,j)=elem(ie)%state%p(i,j,k,nt) + elem(ie)%state%ps(i,j)
                   enddo
                enddo

                lap_p=laplace_sphere_wk(ptot,deriv,elem(ie),var_coef=.false.)
                lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)

                ! advace in time.  (note: DSS commutes with time stepping, so we
                ! can time advance and then DSS.  this has the advantage of
                ! not letting any discontinuties accumulate in p,v via tol
                do j=1,np
                   do i=1,np             
                      elem(ie)%state%p(i,j,k,nt)=elem(ie)%state%p(i,j,k,nt)*spheremp(i,j)  +  dt*nu_s*lap_p(i,j) 
                      elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%state%v(i,j,1,k,nt)*spheremp(i,j) + dt*nu*lap_v(i,j,1)
                      elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%state%v(i,j,2,k,nt)*spheremp(i,j) + dt*nu*lap_v(i,j,2)
                   enddo
                enddo
             enddo

             kptr=0
             call edgeVpack(edge3, elem(ie)%state%p(:,:,:,nt),nlev,kptr,elem(ie)%desc)
             kptr=nlev
             call edgeVpack(edge3,elem(ie)%state%v(:,:,:,:,nt),2*nlev,kptr,elem(ie)%desc)
          enddo

          call bndry_exchangeV(hybrid,edge3)

          do ie=nets,nete
             rspheremp     => elem(ie)%rspheremp

             kptr=0
             call edgeVunpack(edge3, elem(ie)%state%p(:,:,:,nt), nlev, kptr, elem(ie)%desc)
             kptr=nlev
             call edgeVunpack(edge3, elem(ie)%state%v(:,:,:,:,nt), 2*nlev, kptr, elem(ie)%desc)

             ! apply inverse mass matrix
             do k=1,nlev
                do j=1,np
                   do i=1,np             
                      elem(ie)%state%p(i,j,k,nt)=rspheremp(i,j)*elem(ie)%state%p(i,j,k,nt)
                      elem(ie)%state%v(i,j,1,k,nt)=rspheremp(i,j)*elem(ie)%state%v(i,j,1,k,nt)
                      elem(ie)%state%v(i,j,2,k,nt)=rspheremp(i,j)*elem(ie)%state%v(i,j,2,k,nt)
                   enddo
                enddo
             enddo
          enddo
       enddo  ! subcycle
    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  hyper viscosity  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (hypervis_order == 2) then
       do ic=1,hypervis_subcycle
          call biharmonic_wk(elem,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete)

          do ie=nets,nete
             spheremp     => elem(ie)%spheremp
             do k=1,nlev
                ! advace in time.  
                ! note: DSS commutes with time stepping, so we can time advance and then DSS.
                ! note: weak operators alreayd have mass matrix "included"

                do j=1,np
                   do i=1,np
                      elem(ie)%state%p(i,j,k,nt)  =  elem(ie)%state%p(i,j,k,nt)*spheremp(i,j)  -dt*nu_s*ptens(i,j,k,ie)
                      elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%state%v(i,j,1,k,nt)*spheremp(i,j) - dt*nu*vtens(i,j,1,k,ie)
                      elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%state%v(i,j,2,k,nt)*spheremp(i,j) - dt*nu*vtens(i,j,2,k,ie)
                   enddo
                enddo
             enddo

             ! smooth some of the negativities introduced by diffusion:
             if (test_case=="swtc1" .or. test_case=="vortex" .or. test_case=="swirl") then

                   ! this routine takes p  (not p*spheremp)
                   do k=1,nlev
                      elem(ie)%state%p(:,:,k,nt)  =  elem(ie)%state%p(:,:,k,nt)/spheremp(:,:)
                   enddo

		   if ((limiter_option == 8))then
		      call limiter_optim_iter_full(elem(ie)%state%p(:,:,:,nt),elem(ie)%spheremp(:,:),&
		      pmin(:,ie),pmax(:,ie),kmass,.true.,.true.,&
			notreliable)
		   elseif ((limiter_option == 81))then
		      call limiter_optim_iter_full(elem(ie)%state%p(:,:,:,nt),elem(ie)%spheremp(:,:),&
		      pmin(:,ie),pmax(:,ie),kmass,.true.,.false.,notreliable)
		   elseif ((limiter_option == 84))then
		      pmin(:,ie)=0.0d0
		      call limiter_optim_iter_full(ptens(:,:,:,ie),elem(ie)%spheremp(:,:),&
			  pmin(:,ie),pmax(:,ie),kmass,.true.,.false.,notreliable)
          	   else
		      call limiter2d_zero(elem(ie)%state%p(:,:,:,nt),elem(ie)%spheremp,kmass)
		   endif

                   do k=1,nlev
                      elem(ie)%state%p(:,:,k,nt)  =  elem(ie)%state%p(:,:,k,nt)*spheremp(:,:)
                   enddo

             endif


             kptr=0
             call edgeVpack(edge3, elem(ie)%state%p(:,:,:,nt),nlev,kptr,elem(ie)%desc)
             kptr=nlev
             call edgeVpack(edge3,elem(ie)%state%v(:,:,:,:,nt),2*nlev,kptr,elem(ie)%desc)
          enddo

          call bndry_exchangeV(hybrid,edge3)

          do ie=nets,nete
             rspheremp     => elem(ie)%rspheremp

             kptr=0
             call edgeVunpack(edge3, elem(ie)%state%p(:,:,:,nt), nlev, kptr, elem(ie)%desc)
             kptr=nlev
             call edgeVunpack(edge3, elem(ie)%state%v(:,:,:,:,nt), 2*nlev, kptr, elem(ie)%desc)

             ! apply inverse mass matrix
             do k=1,nlev
                do j=1,np
                   do i=1,np
                      elem(ie)%state%p(i,j,k,nt)=rspheremp(i,j)*elem(ie)%state%p(i,j,k,nt)
                      elem(ie)%state%v(i,j,1,k,nt)=rspheremp(i,j)*elem(ie)%state%v(i,j,1,k,nt)
                      elem(ie)%state%v(i,j,2,k,nt)=rspheremp(i,j)*elem(ie)%state%v(i,j,2,k,nt)
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif


    ! convert lat-lon -> contra variant
    do ie=nets,nete
       do k=1,nlev
          do j=1,np
             do i=1,np
                v1=elem(ie)%state%v(i,j,1,k,nt)
                v2=elem(ie)%state%v(i,j,2,k,nt)
                elem(ie)%state%v(i,j,1,k,nt) = elem(ie)%Dinv(1,1,i,j)*v1 + elem(ie)%Dinv(1,2,i,j)*v2
                elem(ie)%state%v(i,j,2,k,nt) = elem(ie)%Dinv(2,1,i,j)*v1 + elem(ie)%Dinv(2,2,i,j)*v2
             enddo
          enddo
       enddo
    enddo


    ! lat-lon conversion 1740  20%
    ! opss+pack+mass     5937  68% 
    ! unpack+mass_inv    1015  12%     
    ! total 8692
  end subroutine advance_hypervis

!-----------------------------------------------------------------------------------

  subroutine advance_si_nonstag(elem, edge1,edge2,   edge3 ,  red     ,            &
       deriv,                    &
       flt  ,                   &
       cg   ,   blkjac    ,  lambdasq, &
       dt   ,   pmean ,  tl      ,            &
       nets ,    nete)

    ! ---------------------
    use kinds, only : real_kind
    ! ---------------------
    use dimensions_mod, only : np, nlev
    ! ---------------------
    use element_mod, only : element_t
    ! ---------------------
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack
    ! ---------------------
    use filter_mod, only : filter_t, filter_P
    ! ---------------------
    use hybrid_mod, only : hybrid_t
    ! ---------------------
    use reduction_mod, only : reductionbuffer_ordered_1d_t
    ! ---------------------
    !    use parallel_mod
    ! ---------------------
    use derivative_mod, only : derivative_t,  gradient_wk, divergence, &
         vorticity, gradient
    ! ---------------------
    use time_mod, only : timelevel_t, smooth
    ! ---------------------
    use control_mod, only : filter_freq, filter_counter, topology, test_case, precon_method
    ! ---------------------
    use shallow_water_mod, only : tc1_velocity
    ! ---------------------
    use cg_mod, only : cg_t
    ! ---------------------
    use solver_mod, only : pcg_solver, blkjac_t, blkjac_init
    ! ---------------------
    use bndry_mod, only : bndry_exchangev
    ! ---------------------
    use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
    ! ---------------------
    !    use schedule_mod
    ! ---------------------
    implicit none
    type (element_t), intent(inout), target :: elem(:)
    type (EdgeBuffer_t)               :: edge1
    type (EdgeBuffer_t)               :: edge2
    type (EdgeBuffer_t)               :: edge3
    type (ReductionBuffer_ordered_1d_t)     :: red
    type (derivative_t)               :: deriv
    type (filter_t)                   :: flt
    type (cg_t)                       :: cg

    integer              , intent(in) :: nets
    integer              , intent(in) :: nete
    type (blkjac_t)                   :: blkjac(nets:nete)

    real (kind=real_kind), intent(in) :: dt
    real (kind=real_kind), intent(in) :: pmean
    real (kind=real_kind), intent(inout) :: lambdasq(nlev)
    type (TimeLevel_t)   , intent(in) :: tl
    ! =================
    ! Local Variables
    ! =================

    ! pointers ...

    real (kind=real_kind), dimension(:,:), pointer     :: mp
    real (kind=real_kind), dimension(:,:), pointer     :: metdetp

    real (kind=real_kind), dimension(:,:), pointer     :: fcor
    real (kind=real_kind), dimension(:,:), pointer     :: rmp
    real (kind=real_kind), dimension(:,:), pointer     :: metdet
    real (kind=real_kind), dimension(:,:), pointer     :: rmetdet
    real (kind=real_kind), dimension(:,:,:,:), pointer :: met
    real (kind=real_kind), dimension(:,:,:,:), pointer :: metinv, Dinv

    ! Thread private working set ...

    real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)  :: Ru
    real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)  :: grad_dp
    real (kind=real_kind), dimension(np,np,nlev,nets:nete)    :: vgradp     ! v.grad(p) on velocity grid
    real (kind=real_kind), dimension(np,np,nlev,nets:nete)    :: Rs   
    real (kind=real_kind), dimension(np,np,nlev,nets:nete)    :: dp ! solution to Helmholtz equation

    real (kind=real_kind), dimension(np,np,2)    :: gradp      ! weak pressure gradient, time level (n  )
    real (kind=real_kind), dimension(np,np,2)    :: gradpm1    ! weak pressure gradient, time level (n-1)
    real (kind=real_kind), dimension(np,np,2)    :: grade      ! strong kinetic energy gradient
    real (kind=real_kind), dimension(np,np,2)    :: gvm1       ! metdet*v(n-1), v contravariant
    real (kind=real_kind), dimension(np,np,2)    :: gv         ! metdet*v(n-1), v contravariant
    real (kind=real_kind), dimension(np,np,2)    :: vco        ! covariant velocity
    real (kind=real_kind), dimension(np,np,2)    :: ulatlon    ! lat-lon velocity

    real (kind=real_kind), dimension(np,np)      :: E          ! kinetic energy term
    real (kind=real_kind), dimension(np,np)      :: zeta       ! relative vorticity

    real (kind=real_kind), dimension(np,np)      :: div        ! timestep n   velocity divergence (p-grid)
    real (kind=real_kind), dimension(np,np)      :: divm1      ! timestep n-1 velocity divergence (p-grid)

    real (kind=real_kind) ::  v1,v2
    real (kind=real_kind) ::  gradp1,gradp2
    real (kind=real_kind) ::  grad_dp1,grad_dp2
    real (kind=real_kind) ::  Ru1,Ru2,Rp

    real (kind=real_kind) ::  dt2
    real (kind=real_kind) ::  time_adv

    real*8  :: et,st
    integer i,j,k,ie
    integer kptr
    integer point
    integer iptr
    integer nm1,n0,np1
    integer nstep

    real (kind=real_kind) :: tmp1
    real (kind=real_kind) :: p0sum,v0sum

    call t_startf('advance_si_nonstag')

    if ( dt /= initialized_for_dt ) then
       if(cg%hybrid%par%masterproc) print *,'Initializing semi-implicit matricies for dt=',dt

       lambdasq(:) = pmean*dt*dt
       if (precon_method == "block_jacobi") then
          call blkjac_init(elem, deriv,lambdasq,nets,nete,blkjac)
       end if
       initialized_for_dt = dt
    endif



    nm1 = tl%nm1
    n0  = tl%n0
    np1 = tl%np1
    nstep = tl%nstep

    dt2 = 2*dt

    !DBG print *,'advance_si: point #1'

#if (! defined VERT_OPENMP)
    !$OMP BARRIER
#endif

    ! ====================
    ! Call Filter...
    ! ====================
    !DBG print *,'advance_si: point #3'

    if (nstep > 0 .and. filter_freq > 0 .and. MODULO(nstep,filter_freq) == 0 ) then

       !DBG print *,'advance_si: point #4'
       do ie=nets,nete
          do k=1,nlev
             call filter_P(elem(ie)%state%p(:,:,k,n0),flt)

             ulatlon(:,:,1)=elem(ie)%D(1,1,:,:)*elem(ie)%state%v(:,:,1,k,n0)+&
                            elem(ie)%D(1,2,:,:)*elem(ie)%state%v(:,:,2,k,n0)
             ulatlon(:,:,2)=elem(ie)%D(2,1,:,:)*elem(ie)%state%v(:,:,1,k,n0)+&
                            elem(ie)%D(2,2,:,:)*elem(ie)%state%v(:,:,2,k,n0)

             do j=1,np
                do i=1,np
                   elem(ie)%state%v(i,j,1,k,n0) = elem(ie)%mp(i,j)*ulatlon(i,j,1)
                   elem(ie)%state%v(i,j,2,k,n0) = elem(ie)%mp(i,j)*ulatlon(i,j,2)
                   elem(ie)%state%p(i,j,k,n0)   = elem(ie)%mp(i,j)*elem(ie)%state%p(i,j,k,n0)
                end do
             end do

          end do
          kptr=0
          call edgeVpack(edge3, elem(ie)%state%v(:,:,:,:,n0),2*nlev,kptr,elem(ie)%desc)
          kptr=2*nlev
          call edgeVpack(edge3, elem(ie)%state%p(:,:,:,n0),nlev,kptr,elem(ie)%desc)
          kptr=0
          !DBG print *,'advance_si: point #6'
       end do

       !DBG print *,'advance_si: point #8'
#if (! defined VERT_OPENMP)
       !$OMP BARRIER
#endif

       call bndry_exchangeV(cg%hybrid,edge3)
#if (! defined VERT_OPENMP)
       !$OMP BARRIER
#endif

       do ie=nets,nete

          kptr=0
          call edgeVunpack(edge3, elem(ie)%state%v(:,:,:,:,n0), 2*nlev, kptr, elem(ie)%desc)
          kptr=2*nlev
          call edgeVunpack(edge3, elem(ie)%state%p(:,:,:,n0), nlev, kptr, elem(ie)%desc)

          do k=1,nlev
             vco(:,:,1) = elem(ie)%Dinv(1,1,:,:)*elem(ie)%state%v(:,:,1,k,n0)+&
                          elem(ie)%Dinv(1,2,:,:)*elem(ie)%state%v(:,:,2,k,n0)
             vco(:,:,2) = elem(ie)%Dinv(2,1,:,:)*elem(ie)%state%v(:,:,1,k,n0)+&
                          elem(ie)%Dinv(2,2,:,:)*elem(ie)%state%v(:,:,2,k,n0)
             do j=1,np
                do i=1,np
                   elem(ie)%state%v(i,j,1,k,n0) = elem(ie)%rmp(i,j)*vco(i,j,1)
                   elem(ie)%state%v(i,j,2,k,n0) = elem(ie)%rmp(i,j)*vco(i,j,2)
                   elem(ie)%state%p(i,j,k,n0)   = elem(ie)%rmp(i,j)*elem(ie)%state%p(i,j,k,n0)
                end do
             end do
          end do

       end do

#if (! defined VERT_OPENMP)
       !$OMP BARRIER
#endif

    end if
    !DBG print *,'advance_si: point #10'

    do ie=nets,nete

       mp => elem(ie)%mp
       fcor => elem(ie)%fcor
       metdet => elem(ie)%metdet
       metinv => elem(ie)%metinv
       Dinv => elem(ie)%Dinv

       !JMD       TIMER_DETAIL_START(timer,2,st)
       !JMD metdet => elem(ie)%metdet
       !JMD if(TIMER_DETAIL(2,timer)) then
       !JMD	 TIMER_START(et)
       !JMD	 timer%pointers = timer%pointers + (et - st)
       !JMD       endif

       !DBG print *,'advance_si: point #11'
       do k=1,nlev

          ! ==============================================
          !
          ! Compute gradient of pressure field at time 
          ! level n-1 and n
          ! 
          !   2 x {2.0*(np+np)*np*(2.0*np-1.0) + 2*np*np} Flops
          !
          ! ==============================================
! this would become gradient_sphere

       ! scale  by rearth to get scaling right between phi and v.
#ifdef _WK_GRAD
          gradpm1(:,:,:)=gradient_wk(elem(ie)%state%p(:,:,k,nm1),deriv)*rrearth
          gradp(:,:,:)  =gradient_wk(elem(ie)%state%p(:,:,k,n0), deriv)*rrearth
#else
          gradpm1(:,:,:)=gradient(elem(ie)%state%p(:,:,k,nm1),deriv)*rrearth
          gradp(:,:,:)  =gradient(elem(ie)%state%p(:,:,k,n0), deriv)*rrearth
! ---
#endif
! --- endif _WK_GRAD

          ! ==============================================
          !
          ! Compute kinetic energy term: 10*np*np Flops
          !
          ! ==============================================
          do j=1,np
             do i=1,np
                v1     = elem(ie)%state%v(i,j,1,k,n0)
                v2     = elem(ie)%state%v(i,j,2,k,n0)

                vco(i,j,1) = elem(ie)%met(1,1,i,j)*v1 + elem(ie)%met(1,2,i,j)*v2
                vco(i,j,2) = elem(ie)%met(2,1,i,j)*v1 + elem(ie)%met(2,2,i,j)*v2

                E(i,j) = 0.5D0*( vco(i,j,1)*v1 + vco(i,j,2)*v2 )

             end do
          end do

          ! =========================================
          !
          ! Compute metdet * relative vorticity (zeta)
          !
          !   2.0*np*np*(2.0*np-1.0) + 3 np*np Flops
          !
          ! =========================================

          !DBG print *,'advance_si: point #12'
          zeta(:,:)  = vorticity(vco,deriv)*rrearth

          ! ==============================================
          !
          ! Compute vgradient of kinetic energy field
          !
          !   2.0*np*np*(2.0*np-1.0) + 2*np*np Flops
          !
          ! ==============================================
          grade(:,:,:)=gradient(E,deriv)*rrearth


          ! ==============================================
          !
          ! Compute Ru^i and v.grad(p) term
          !
          !    23*np*np Flops
          !
          ! ==============================================
          !DBG print *,'advance_si: point #13 ie:=',ie
          do j=1,np
             do i=1,np

!make sure its consistent lat lon NOT not contra or co-variant
                Ru1 =  mp(i,j)*(   elem(ie)%state%v(i,j,2,k,n0)*(metdet(i,j)*fcor(i,j) + zeta(i,j)) &
                     - grade(i,j,1))                                              &
                     + gradpm1(i,j,1) + elem(ie)%state%gradps(i,j,1)

                Ru2 =  mp(i,j)*( - elem(ie)%state%v(i,j,1,k,n0)*(metdet(i,j)*fcor(i,j) + zeta(i,j)) &
                     - grade(i,j,2))                                               &
                     + gradpm1(i,j,2) + elem(ie)%state%gradps(i,j,2)
                Ru(i,j,1,k,ie)   = dt2*(Dinv(1,1,i,j)*Ru1 + Dinv(2,1,i,j)*Ru2)
                Ru(i,j,2,k,ie)   = dt2*(Dinv(1,2,i,j)*Ru1 + Dinv(2,2,i,j)*Ru2)

                vgradp(i,j,k,ie)  =  elem(ie)%state%v(i,j,1,k,n0)*gradp(i,j,1) + &
                     elem(ie)%state%v(i,j,2,k,n0)*gradp(i,j,2)

             end do
          end do
       end do


       ! ===================================================
       !
       ! Pack cube edges of grad(p) and V.grad(V) 
       ! into edge buffer
       !
       ! ===================================================

       !DBG print *,'advance_si: point #14'
       kptr=0
       call edgeVpack(edge3, vgradp(1,1,1,ie),nlev,kptr,elem(ie)%desc)

       kptr=nlev
       call edgeVpack(edge3,Ru(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)
       !DBG print *,'advance_si: point #15'

       ! =============================================================
       !
       ! Rotate edges (if necessary, e.g. if we're on the cube)
       !
       ! ============================================================= 

       !DBG print *,'advance_si: point #15.1'
    end do

    ! =============================================================
    ! Insert communications here: for shared memory, just a single
    ! thread barrier is required
    ! =============================================================

    call bndry_exchangeV(cg%hybrid,edge3)
#if (! defined VERT_OPENMP)
    !$OMP BARRIER
#endif

    do ie=nets,nete

       rmp     => elem(ie)%rmp
       mp      => elem(ie)%mp
       Dinv    => elem(ie)%Dinv
       metdet  => elem(ie)%metdet
       rmetdet=> elem(ie)%rmetdet

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       !   3*4*(np+1)
       ! ===========================================================

       kptr=0
       call edgeVunpack(edge3, vgradp(1,1,1,ie), nlev, kptr, elem(ie)%desc)

       kptr=nlev
       call edgeVunpack(edge3, Ru(1,1,1,1,ie), 2*nlev, kptr, elem(ie)%desc)

       ! ===========================================================
       ! Compute velocity and pressure tendencies for all levels
       ! ===========================================================

       do k=1,nlev

          ! =========================================================
          !
          ! Scale velocity tendency and v.grad(p) term by inverse mass
          ! matrix, scale velocity by metric g factor.
          !
          !     11 np*np Flops
          !
          ! =========================================================

          do j=1,np
             do i=1,np
                vgradp(i,j,k,ie) = rmp(i,j)*vgradp(i,j,k,ie)

                Ru1 = Ru(i,j,1,k,ie)
                Ru2 = Ru(i,j,2,k,ie)
                Ru(i,j,1,k,ie) = rmp(i,j)*(Dinv(1,1,i,j)*Ru1+Dinv(1,2,i,j)*Ru2)
                Ru(i,j,2,k,ie) = rmp(i,j)*(Dinv(2,1,i,j)*Ru1+Dinv(2,2,i,j)*Ru2)

                gv(i,j,1)   = metdet(i,j)*elem(ie)%state%v(i,j,1,k,n0)
                gv(i,j,2)   = metdet(i,j)*elem(ie)%state%v(i,j,2,k,n0)

                gvm1(i,j,1) = metdet(i,j)*(elem(ie)%state%v(i,j,1,k,nm1) + 0.5D0*Ru(i,j,1,k,ie))
                gvm1(i,j,2) = metdet(i,j)*(elem(ie)%state%v(i,j,2,k,nm1) + 0.5D0*Ru(i,j,2,k,ie))

             end do
          end do
          ! ==========================================================
          !
          ! Compute divergence of metdet*v(n  ), v(n  ) contravariant
          ! Compute divergence of metdet*v(n-1), v(n-1) contravariant
          !
          !  2 x [ 2*(np+np)*np*(2*np-1) + 3*np*np ] Flops
          ! ==========================================================

          div    = divergence(gv  ,deriv)*rrearth
          divm1  = divergence(gvm1,deriv)*rrearth

          ! ====================================================
          ! Compute the Right Hand Side of the Helmholtz eq
          !
          !    7*np*np Flops
          !
          ! ====================================================

          iptr=1
          do j=1,np
             do i=1,np
                Rs(i,j,k,ie) = dt2*mp(i,j)*(vgradp(i,j,k,ie)*metdetp(i,j)   &
                     - elem(ie)%state%p(i,j,k,n0)*div(i,j)       &
                     -                 pmean*divm1(i,j))
             end do
          end do

       end do
       !DBG print *,'advance_si: point #14'
       kptr=0
       call edgeVpack(edge1, Rs(1,1,1,ie),nlev,kptr,elem(ie)%desc)

    end do

    call bndry_exchangeV(cg%hybrid,edge1)
#if (! defined VERT_OPENMP)
    !$OMP BARRIER
#endif

    do ie=nets,nete

       rmp     => elem(ie)%rmp
       kptr=0
       call edgeVunpack(edge1, Rs(1,1,1,ie), nlev, kptr, elem(ie)%desc)
       do k=1,nlev
          do j=1,np
             do i=1,np
                Rs(i,j,k,ie) = rmp(i,j)*Rs(i,j,k,ie)
             enddo
          enddo
       enddo
    enddo

    ! ======================================================
    ! Invoke the solver!
    ! ======================================================

    !DBG print *,'advance_si: before call to pcg_solver'
    point = 3
#ifdef _HTRACE
    !JMD    call EVENT_POINT(point)
! ---
#endif
! --- endif _HTRACE

#if (! defined VERT_OPENMP)
    !$OMP BARRIER
#endif


    dp(:,:,:,nets:nete) = pcg_solver(elem, &
	 Rs(:,:,:,nets:nete),  &     ! rhs of Helmholtz problem
         cg,              &     ! cg struct
         red,             &     ! reduction buffer
         edge1 ,          &     ! single vector edge exchange buffer
         edge2 ,          &     ! single vector edge exchange buffer
         lambdasq,        &     ! Helmholtz length scale squared
         deriv,      &     ! staggered derivative struct
         nets,            &     ! starting element number
         nete,            &     ! ending   element number
         blkjac)
    point = 4


#ifdef _HTRACE
    !JMD   call EVENT_POINT(point)
! ---
#endif
! --- endif _HTRACE

    do ie=nets,nete

       metinv => elem(ie)%metinv
       Dinv => elem(ie)%Dinv

       do k=1,nlev

          ! compute grad dp needed to back substitute for du

#ifdef _WK_GRAD
          grad_dp(:,:,:,k,ie)=gradient_wk(dp(:,:,k,ie),deriv)*rrearth
#else
          grad_dp(:,:,:,k,ie)=gradient(dp(:,:,k,ie),deriv)*rrearth
! ---
#endif
! --- endif _WK_GRAD


          ! ==================================================
          ! Rotate grad_dp to form contravariant object:
          !        6 np*np Flops
          ! ==================================================

          do j=1,np
             do i=1,np
                grad_dp1 = grad_dp(i,j,1,k,ie)
                grad_dp2 = grad_dp(i,j,2,k,ie)
                grad_dp(i,j,1,k,ie) = Dinv(1,1,i,j)*grad_dp1 + Dinv(2,1,i,j)*grad_dp2
                grad_dp(i,j,2,k,ie) = Dinv(1,2,i,j)*grad_dp1 + Dinv(2,2,i,j)*grad_dp2

             end do
          end do

       enddo

       kptr=0
       call edgeVpack(edge2, grad_dp(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)
    end do

#if (! defined VERT_OPENMP)
    !$OMP BARRIER
#endif
    call bndry_exchangeV(cg%hybrid,edge2)
#if (! defined VERT_OPENMP)
    !$OMP BARRIER
#endif
    do ie=nets,nete
       rmp => elem(ie)%rmp
       Dinv => elem(ie)%Dinv


       kptr=0      

       call edgeVunpack(edge2, grad_dp(1,1,1,1,ie), 2*nlev, kptr, elem(ie)%desc)

       do k=1,nlev

          ! ==============================
          ! Update geopotential
          !    6 np*np Flops
          ! ==============================


          iptr=1
          do j=1,np
             do i=1,np
                elem(ie)%state%p(i,j,k,np1) = elem(ie)%state%p(i,j,k,nm1) + dp(i,j,k,ie)
                elem(ie)%state%p(i,j,k,n0)  = elem(ie)%state%p(i,j,k,n0) + smooth*(elem(ie)%state%p(i,j,k,nm1) &
                     - 2.0D0*elem(ie)%state%p(i,j,k,n0) + elem(ie)%state%p(i,j,k,np1))
                iptr=iptr+1
             end do
          end do

          ! ==============================
          ! Update velocity
          !    16 np*np Flops
          ! ==============================

          do j=1,np
             do i=1,np
                grad_dp1 = rmp(i,j)* &
                          (Dinv(1,1,i,j)*grad_dp(i,j,1,k,ie)+Dinv(1,2,i,j)*grad_dp(i,j,2,k,ie))
                grad_dp2 = rmp(i,j)* &
                          (Dinv(2,1,i,j)*grad_dp(i,j,1,k,ie)+Dinv(2,2,i,j)*grad_dp(i,j,2,k,ie))

                elem(ie)%state%v(i,j,1,k,np1) = elem(ie)%state%v(i,j,1,k,nm1) + Ru(i,j,1,k,ie) + dt*grad_dp1
                elem(ie)%state%v(i,j,2,k,np1) = elem(ie)%state%v(i,j,2,k,nm1) + Ru(i,j,2,k,ie) + dt*grad_dp2
                elem(ie)%state%v(i,j,1,k,n0)  = elem(ie)%state%v(i,j,1,k,n0) + smooth*(elem(ie)%state%v(i,j,1,k,nm1) &
                     - 2.0D0*elem(ie)%state%v(i,j,1,k,n0) + elem(ie)%state%v(i,j,1,k,np1))
                elem(ie)%state%v(i,j,2,k,n0)  = elem(ie)%state%v(i,j,2,k,n0) + smooth*(elem(ie)%state%v(i,j,2,k,nm1) &
                     - 2.0D0*elem(ie)%state%v(i,j,2,k,n0) + elem(ie)%state%v(i,j,2,k,np1))

             end do
          end do
       end do
    end do

    call t_stopf('advance_si_nonstag')

#if (! defined VERT_OPENMP)
    !$OMP BARRIER
#endif

  end subroutine advance_si_nonstag


!----------------------------------------------------------------------------------------


  subroutine set_prescribed_velocity(elem,n0,time)
  use control_mod, only :  topology, test_case
  use element_mod, only : element_t
  use dimensions_mod, only : np, nlev
  use shallow_water_mod, only : tc1_velocity, vortex_velocity, swirl_velocity
  implicit none

  type (element_t)     , intent(inout) :: elem

  ! local
  integer :: n0,k
  real (kind=real_kind) :: time

  if (topology == "cube" .and. test_case=="swtc1") then
     do k=1,nlev
        elem%state%v(:,:,:,k,n0)=tc1_velocity(elem%spherep,elem%Dinv)
     end do
  else if (topology == "cube" .and. test_case=="vortex") then                
     do k=1,nlev
        elem%state%v(:,:,:,k,n0)=vortex_velocity(time,elem%spherep,elem%Dinv)
     end do
  else if (topology == "cube" .and. test_case=="swirl") then                
     do k=1,nlev
        elem%state%v(:,:,:,k,n0)=swirl_velocity(time,elem%spherep,elem%Dinv)
     end do

  end if
  end subroutine set_prescribed_velocity

!----------------------------------------------------------------------------------------


  subroutine compute_and_apply_rhs(np1,nm1,n0,dt2,real_time,edge3,elem,pmean,hybrid,deriv,vtens,ptens,nets,nete)
  ! ===================================
  ! compute the RHS, accumulate into u(np1) and apply DSS
  !
  !           u(np1) = u(nm1) + dt2*DSS[ RHS(u(n0)) ]
  !
  ! This subroutine is normally called to compute a leapfrog timestep
  ! but by adjusting np1,nm1,n0 and dt2, many other timesteps can be
  ! accomodated.  For exaple, setting nm1=np1=n0 this routine will
  ! take a forward euler step, overwriting the input with the output.
  !
  ! if  dt2=0, then the DSS'd RHS is returned in vtens,ptens
  ! and u(np1) is not changed.  
  !
  ! Combining the RHS and DSS pack operation in one routine 
  ! allows us to fuse these two loops for more cache reuse
  !
  ! Combining the dt advance and DSS unpack operation in one routine 
  ! allows us to fuse these two loops for more cache reuse
  !
  ! note: for prescribed velocity case, velocity will be computed at
  ! "real_time", which should be the time of timelevel n0.  
  ! ===================================
  use kinds, only : real_kind
  use dimensions_mod, only : np, nlev
  use hybrid_mod, only : hybrid_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t, divergence_sphere, gradient_sphere, vorticity_sphere
  use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack
  use bndry_mod, only : bndry_exchangev
  implicit none

  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)  :: vtens
  real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: ptens
  type (EdgeBuffer_t)  , intent(inout) :: edge3
  type (derivative_t)  , intent(in) :: deriv
  real (kind=real_kind) :: dt2,pmean,real_time
  integer :: nm1,np1,n0,nets,nete

  ! local
  ! pointer ...
  real (kind=real_kind), dimension(:,:), pointer :: fcor,rspheremp,spheremp,metdet,rmetdet
  real (kind=real_kind), dimension(:,:,:,:), pointer :: met,metinv
  real (kind=real_kind), dimension(np,np,2)    :: grade   ! strong kinetic energy gradient
  real (kind=real_kind), dimension(np,np,2)    :: pv      ! p*v lat-lon
  real (kind=real_kind), dimension(np,np)                     :: E          ! kinetic energy term
  real (kind=real_kind), dimension(np,np)                     :: zeta       ! relative vorticity
  real (kind=real_kind), dimension(np,np)      :: div  
  real (kind=real_kind), dimension(np,np,2)      :: ulatlon

  integer i,j,k,kptr,ie
  real (kind=real_kind) ::  v1,v2
  real (kind=real_kind) ::  vtens1,vtens2




  ! ===================================
  ! construct v tendencies and v.grad(p)
  ! on the velocity grid...
  ! ===================================
  do ie=nets,nete
     met    => elem(ie)%met
     metinv => elem(ie)%metinv
     metdet => elem(ie)%metdet
     rmetdet => elem(ie)%rmetdet
     fcor   => elem(ie)%fcor
     spheremp     => elem(ie)%spheremp

     call set_prescribed_velocity(elem(ie),n0,real_time)

     do k=1,nlev
        ! ==============================================
        ! Compute kinetic energy term
        ! ==============================================
        do j=1,np
           do i=1,np
              
              v1     = elem(ie)%state%v(i,j,1,k,n0)   ! contra
              v2     = elem(ie)%state%v(i,j,2,k,n0)   ! contra 
              ulatlon(i,j,1)=elem(ie)%D(1,1,i,j)*v1 + elem(ie)%D(1,2,i,j)*v2   ! contra->latlon
              ulatlon(i,j,2)=elem(ie)%D(2,1,i,j)*v1 + elem(ie)%D(2,2,i,j)*v2   ! contra->latlon
              
              E(i,j) = 0.5D0*(ulatlon(i,j,1)**2 + ulatlon(i,j,2)**2)  +&
                   elem(ie)%state%p(i,j,k,n0) + elem(ie)%state%ps(i,j)
              
              pv(i,j,1) = ulatlon(i,j,1)*(pmean+elem(ie)%state%p(i,j,k,n0))
              pv(i,j,2) = ulatlon(i,j,2)*(pmean+elem(ie)%state%p(i,j,k,n0))
           end do
        end do
        grade = gradient_sphere(E,deriv,elem(ie)%Dinv)       ! scalar -> latlon vector
        zeta = vorticity_sphere(ulatlon,deriv,elem(ie)) ! latlon vector -> scalar 
        div = divergence_sphere(pv,deriv,elem(ie))      ! latlon vector -> scalar 
        
        ! ==============================================
        ! Compute velocity tendency terms
        ! ==============================================
        do j=1,np
           do i=1,np
              ! accumulate strong form terms, apply mass matrix
              vtens(i,j,1,k,ie)=spheremp(i,j)*(ulatlon(i,j,2)*(fcor(i,j) + zeta(i,j))  - grade(i,j,1))
              vtens(i,j,2,k,ie)=spheremp(i,j)*(-ulatlon(i,j,1)*(fcor(i,j) + zeta(i,j)) - grade(i,j,2))
              ptens(i,j,k,ie) =  -spheremp(i,j)*div(i,j)
           end do
        end do
     end do
     
     ! ===================================================
     ! Pack cube edges of tendencies, rotate velocities
     ! ===================================================
     kptr=0
     call edgeVpack(edge3, ptens(1,1,1,ie),nlev,kptr,elem(ie)%desc)
     kptr=nlev
     call edgeVpack(edge3,vtens(1,1,1,1,ie),2*nlev,kptr,elem(ie)%desc)
  end do
  
  
#if (! defined VERT_OPENMP)
  !$OMP BARRIER
#endif
  call bndry_exchangeV(hybrid,edge3)
#if (! defined VERT_OPENMP)
  !$OMP BARRIER
#endif
  
  do ie=nets,nete
     rspheremp     => elem(ie)%rspheremp
     
     ! ===========================================================
     ! Unpack the edges for vgradp and vtens
     ! ===========================================================
     kptr=0
     call edgeVunpack(edge3, ptens(1,1,1,ie), nlev, kptr, elem(ie)%desc)
     
     kptr=nlev
     call edgeVunpack(edge3, vtens(1,1,1,1,ie), 2*nlev, kptr, elem(ie)%desc)
     
     ! ===========================================================
     ! Compute velocity and pressure tendencies for all levels
     ! ===========================================================
     do k=1,nlev
        do j=1,np
           do i=1,np
              ptens(i,j,k,ie) = rspheremp(i,j)*ptens(i,j,k,ie)
              vtens1=rspheremp(i,j)*vtens(i,j,1,k,ie)
              vtens2=rspheremp(i,j)*vtens(i,j,2,k,ie)
              
              ! lat-lon -> contra
              vtens(i,j,1,k,ie) = elem(ie)%Dinv(1,1,i,j)*vtens1 + elem(ie)%Dinv(1,2,i,j)*vtens2
              vtens(i,j,2,k,ie) = elem(ie)%Dinv(2,1,i,j)*vtens1 + elem(ie)%Dinv(2,2,i,j)*vtens2
           end do
        end do
     end do
     
     if (dt2/=0) then
     do k=1,nlev
        ! ====================================================
        ! Update
        ! ====================================================
        do j=1,np
           do i=1,np
              elem(ie)%state%v(i,j,1,k,np1) = elem(ie)%state%v(i,j,1,k,nm1) + dt2*vtens(i,j,1,k,ie)
              elem(ie)%state%v(i,j,2,k,np1) = elem(ie)%state%v(i,j,2,k,nm1) + dt2*vtens(i,j,2,k,ie)
              elem(ie)%state%p(i,j,k,np1) = elem(ie)%state%p(i,j,k,nm1) + dt2*ptens(i,j,k,ie)
           end do
        end do
     end do
     endif
  end do
  end subroutine compute_and_apply_rhs
  



!=======================================================================================!
!  Advection Flux Term							   		!
!  from DG code - modified for CG velocity
!=======================================================================================!
function adv_flux_term(elem,deriv,contrauv,si,si_neighbor) result(numflux)
!=======================================================================================!
    use derivative_mod, only : derivative_t
    use element_mod, only : element_t
    integer, parameter :: south=1, east=2, north=3, west=4
    type (derivative_t)         :: deriv
    real (kind=real_kind), dimension(np,np,2),intent(in) :: contrauv
    real (kind=real_kind), dimension(np,np),  intent(in) :: si
    real (kind=real_kind), dimension(0:np+1,0:np+1),   intent(in) :: si_neighbor
    type (element_t) :: elem

    real (kind=real_kind), dimension(np,4) :: si_senw
    real (kind=real_kind), dimension(np,np) :: numflux
    real (kind=real_kind), dimension(np,np) :: mij
    real (kind=real_kind), dimension(np)   :: lf_south,lf_north,lf_east,lf_west
    real(kind=real_kind) ::  alfa1, alfa2, ul,ur , left, right, f_left, f_right, s1,s2
    integer i,j, k
!=======================================================================================!
    mij(:,:)  = 0
    mij(1,1)  = 1
    mij(np,np)= 1

    ! convert from edgeDGVunpack variable to Ram's variables:
    si_senw(:,south) = si_neighbor(1:np,0)
    si_senw(:,north) = si_neighbor(1:np,np+1)
    si_senw(:,east) = si_neighbor(np+1,1:np)
    si_senw(:,west) = si_neighbor(0,1:np)

    ! South & North    Max flux Jacobians
!   ul = abs(contrauv(i,1,2))
    alfa1 = maxval(abs(contrauv(1:np,1,2)))
!   ur = abs(contrauv(i,np,2))
    alfa2 = maxval(abs(contrauv(1:np,np,2)))

    do i = 1, np
       ! South wall
       left = si_senw(i,south)
       right  = si(i,1)
!       f_left = fxy_halo(i,south,2)
!       f_right = fxy(i,1,2)
       f_left = left*contrauv(i,1,2)
       f_right= right*contrauv(i,1,2)
       lf_south(i) =  0.5D0 *(f_left + f_right - alfa1*(right - left))

       ! North wall
       left  = si(i,np)
       right = si_senw(i,north)
!       f_left  = fxy(i,np,2)
!       f_right = fxy_halo(i,north,2)
       f_left  = left*contrauv(i,np,2)
       f_right = right*contrauv(i,np,2)
       lf_north(i) =  0.5D0 *(f_left + f_right - alfa2*(right - left))

    enddo


    ! East & West   max of Flux Jacobians

!    ur = abs(contrauv(1,j,1))
    alfa1 = maxval(abs(contrauv(1,1:np,1)))
!    ul = abs(contrauv(np,j,1))
    alfa2 = maxval(abs(contrauv(np,1:np,1)))

    do j = 1, np
       !West wall
       left  =  si_senw(j,west)
       right  = si(1,j)
!       f_left  = fxy_halo(j,west,1)
!       f_right = fxy(1,j,1)
       f_left = left*contrauv(1,j,1)
       f_right= right*contrauv(1,j,1)
       lf_west(j) =  0.5D0 *(f_left + f_right - alfa1*(right - left))

       !East wall
       left  = si(np,j)
       right  = si_senw(j,east)
!       f_left  = fxy(np,j,1)
!       f_right = fxy_halo(j,east,1)
       f_left = left*contrauv(np,j,1)
       f_right= right*contrauv(np,j,1)
       lf_east(j) =  0.5D0 *(f_left + f_right - alfa2*(right - left))

    enddo

    !Flux integral along the element boundary
    ! note: added metdet() below which is not used in DG code.  
    ! in CG code, all integrals used in weak formulation match physical integrals used
    ! to define mass and energy:
    ! ds (arc length)   =  metdet(i,j)*w(i)
    ! dA (area measure) =  metdet(i,j)*w(i)*w(j)

    do j = 1, np
       do i = 1, np

          s1 =  (lf_east(j) *mij(i,np) - lf_west(j) *mij(i,1) )* deriv%Mvv_twt(j,j)*elem%metdet(i,j)
          s2 =  (lf_north(i)*mij(j,np) - lf_south(i)*mij(j,1) )* deriv%Mvv_twt(i,i)*elem%metdet(i,j)

          numflux(i,j) = (s1 + s2) * rrearth

       enddo
    enddo
!=======================================================================================!
end function adv_flux_term


!==================================================================================!
! Element-wise Max flux Jacobian for SW system 
! from DG code, modified for CG velocity
!----------------------------------------------------------------------------------
 Function sw_fjmax(contuv,si,si_neighbor,elem) result(fjmax)
!----------------------------------------------------------------------------------
 use element_mod, only : element_t
 Implicit None
 real (kind=real_kind),dimension(np,np,2),intent(in):: contuv
 real (kind=real_kind), dimension(np,np),  intent(in) :: si
 real (kind=real_kind), dimension(0:np+1,0:np+1),   intent(in) :: si_neighbor
 type (element_t) :: elem

 real (kind=real_kind),dimension(np,np) :: gh11,gh22,g11,g22
 real (kind=real_kind),dimension(np,4)  :: gh11_halo, gh22_halo 
 real (kind=real_kind),dimension(4)   :: fjmax
 real (kind=real_kind):: alfa1,alfa2, ul,ur 
 integer, parameter:: south=1,east=2,north=3,west=4
 integer:: i,j,wall
!========================================================
#if 0
    ! for debugging: with this, we should duplicate adv_flux_term()
    g11=0
    g22=0
#else
    g11=(elem%metinv(1,1,:,:))   ! sqrt(g11)=contra component of nhat on east/west edges
    g22=(elem%metinv(2,2,:,:))   ! sgrt(g22)=contra component of nhat on north/south edges
#endif

    gh11(:,:) = (si(:,:))*g11(:,:)
    gh22(:,:) = (si(:,:))*g22(:,:)

    ! convert from edgeDGVunpack variable to Ram's variables:
    gh11_halo(:,south) = (si_neighbor(1:np,0))*g11(1:np,1)
    gh11_halo(:,north) = (si_neighbor(1:np,np+1))*g11(1:np,np)
    gh11_halo(:,east) = (si_neighbor(np+1,1:np))*g11(np,1:np)
    gh11_halo(:,west) = (si_neighbor(0,1:np))*g11(1,1:np)

    gh22_halo(:,south) = (si_neighbor(1:np,0))*g22(1:np,1)
    gh22_halo(:,north) = (si_neighbor(1:np,np+1))*g22(1:np,np)
    gh22_halo(:,east) = (si_neighbor(np+1,1:np))*g22(np,1:np)
    gh22_halo(:,west) = (si_neighbor(0,1:np))*g22(1,1:np)




    alfa1 = 0.0D0
    alfa2 = 0.0D0
    do i = 1, np
       ul = abs(contuv(i,1,2))      + sqrt(gh22(i,1))
       ur = abs(contuv(i,1,2))      + sqrt(gh22_halo(i,south))
       alfa1= max(alfa1,ul,ur)
       
       ul = abs(contuv(i,np,2)) + sqrt(gh22_halo(i,north))
       ur = abs(contuv(i,np,2))     + sqrt(gh22(i,np))
       alfa2= max(alfa2,ul,ur)
    enddo

    fjmax(south) = alfa1
    fjmax(north) = alfa2

    alfa1 = 0.0D0
    alfa2 = 0.0D0
    do j = 1, np
       ul = abs(contuv(1,j,1))  + sqrt(gh11_halo(j,west))
       ur = abs(contuv(1,j,1))  + sqrt(gh11(1,j))
     alfa1= max(alfa1,ul,ur)
       ul = abs(contuv(np,j,1))  + sqrt(gh11(np,j))
       ur = abs(contuv(np,j,1))  + sqrt(gh11_halo(j,east))
     alfa2= max(alfa2,ul,ur)
    enddo

    fjmax(west) = alfa1
    fjmax(east) = alfa2 

   End Function sw_fjmax  
!=======================================================================================!

!=======================================================================================! 
subroutine swsys_flux(numeqn,elem,deriv,fjmax,si,si_neighbor,uvcontra,fluxout)
   use element_mod, only : element_t
   use derivative_mod, only : derivative_t
   integer, parameter :: south=1, east=2, north=3, west=4
   integer, intent(in):: numeqn
   type (derivative_t)                                  :: deriv
   type (element_t) :: elem
   real (kind=real_kind), dimension(4),    intent(in)   :: fjmax
   real (kind=real_kind), dimension(np,np,2,numeqn),intent(in) :: uvcontra
   real (kind=real_kind), dimension(np,np,numeqn),  intent(in) :: si
   real (kind=real_kind), dimension(0:np+1,0:np+1,numeqn),   intent(in) :: si_neighbor
   real (kind=real_kind), dimension(np,np,numeqn), intent(out) :: fluxout

   real (kind=real_kind), dimension(np,4,numeqn) :: si_senw
   real (kind=real_kind), dimension(np,np) :: mij
   real (kind=real_kind), dimension(np)   :: lf_south,lf_north,lf_east,lf_west
   real(kind=real_kind) ::  ul,ur , left, right, f_left, f_right, s1,s2

   integer i,j,k,eqn 

      mij(:,:) = 0.0D0
      mij(1,1) = 1.0D0
    mij(np,np) = 1.0D0

    !For SW-system  (u1,u2,dp,pt) order   (4 equations)

    fluxout(:,:,:) = 0.0D0 

      do eqn = 1, numeqn
         ! convert from edgeDGVunpack variable to Ram's variables:
         si_senw(:,south,eqn) = si_neighbor(1:np,0,eqn)
         si_senw(:,north,eqn) = si_neighbor(1:np,np+1,eqn)
         si_senw(:,east,eqn) = si_neighbor(np+1,1:np,eqn)
         si_senw(:,west,eqn) = si_neighbor(0,1:np,eqn)
         

           ! East & West   LF flux  (fjmax <- max of flux Jacobian)

         do j = 1, np

                    left  = si_senw(j,west,eqn) 
                   right  = si(1,j,eqn) 
                  f_left  = uvcontra(1,j,1,eqn)*left
                  f_right = uvcontra(1,j,1,eqn)*right
               lf_west(j) = 0.5D0 *(f_left + f_right - fjmax(west)*(right - left))

                    left  = si(np,j,eqn) 
                   right  = si_senw(j,east,eqn) 
                  f_left  = uvcontra(np,j,1,eqn)*left
                  f_right = uvcontra(np,j,1,eqn)*right
               lf_east(j) = 0.5D0 *(f_left + f_right - fjmax(east)*(right - left))

          end do

           ! North& South  LF flux  (fjmax <- max of flux Jacobian)

         do i = 1, np

                     left = si_senw(i,south,eqn) 
                   right  = si(i,1,eqn) 
                   f_left = uvcontra(i,1,2,eqn)*left
                  f_right = uvcontra(i,1,2,eqn)*right
              lf_south(i) = 0.5D0 *(f_left + f_right - fjmax(south)*(right - left))

                    left  = si(i,np,eqn) 
                    right = si_senw(i,north,eqn) 
                  f_left  = uvcontra(i,np,2,eqn)*left
                  f_right = uvcontra(i,np,2,eqn)*right
              lf_north(i) = 0.5D0 *(f_left + f_right - fjmax(north)*(right - left))

         end do

        !Flux integral along the element boundary
        ! note: added metdet() below which is not used in DG code.  
        ! in CG code, all integrals used in weak formulation match physical integrals used
        ! to define mass and energy:
        ! v dot nhat ds (arc length)   =  ucontra * metdet(i,j)*w(i)
        ! dA (area measure)            =  metdet(i,j)*w(i)*w(j)
        do j = 1, np
        do i = 1, np
            s1 = (lf_east(j) *mij(i,np) - lf_west(j) *mij(i,1) )* deriv%Mvv_twt(j,j) * elem%metdet(i,j) 
            s2 = (lf_north(i)*mij(j,np) - lf_south(i)*mij(j,1) )* deriv%Mvv_twt(i,i) * elem%metdet(i,j)
            fluxout(i,j,eqn) = (s1 + s2) * rrearth
        end do
        end do

    end do 

 end subroutine  swsys_flux 
!=======================================================================================================! 

end module advance_mod
