#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module advance_mod
  use kinds,              only: real_kind
  use dimensions_mod,     only: np
  use physical_constants, only: rearth, rrearth

  implicit none

  ! semi-implicit needs to be re-initialized each time dt changes
  real (kind=real_kind) :: initialized_for_dt   = 0

contains

  subroutine advance_nonstag( elem, edge2,  edge3,  deriv,  hybrid,  &
       dt,  pmean,     tl,   nets,   nete)

    use kinds,          only: real_kind
    use dimensions_mod, only: np, nlev
    use element_mod,    only: element_t
    use edgetype_mod,   only: EdgeBuffer_t
    use hybrid_mod,     only: hybrid_t
    use derivative_mod, only: derivative_t
    use time_mod,       only: timelevel_t, smooth
    use control_mod,    only: LFTfreq
    use perf_mod,       only : t_startf, t_stopf ! _EXTERNAL
    use global_norms_mod

    implicit none

    type (element_t)     , intent(inout), target :: elem(:)
    type (EdgeBuffer_t)  , intent(in) :: edge2
    type (EdgeBuffer_t)  , intent(inout) :: edge3

    type (derivative_t)  , intent(in) :: deriv

    type (hybrid_t)      , intent(in) :: hybrid

    real (kind=real_kind), intent(in) :: dt
    real (kind=real_kind), intent(in) :: pmean
    type (TimeLevel_t)   , intent(in) :: tl
    integer              , intent(in) :: nets
    integer              , intent(in) :: nete

    ! =================
    ! Local
    ! =================
    ! Thread private working set ...
    real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)  :: vtens
    real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: ptens


    real (kind=real_kind) :: dt2, real_time

    integer    :: i,j,k,ie
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

#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif
  end subroutine advance_nonstag

!--------------------------------------------------------------------------------------------

  subroutine advance_nonstag_rk( MyRk, elem, edge2,  edge3,  deriv,  hybrid,  &
       dt,  pmean,     tl,   nets,   nete)

    use kinds,          only: real_kind
    use dimensions_mod, only: np, nlev
    use element_mod,    only: element_t
    use edge_mod,       only: edgevpack, edgevunpack
    use edgetype_mod,   only: EdgeBuffer_t
    use hybrid_mod,     only: hybrid_t
    use derivative_mod, only: derivative_t, gradient_sphere, divergence_sphere,vorticity_sphere,&
         divergence_sphere_wk
    use time_mod,       only: timelevel_t
    use control_mod,    only:  test_case, limiter_option, nu, nu_s, kmass
    use bndry_mod,      only: bndry_exchangev
    use viscosity_mod,  only: neighbor_minmax, biharmonic_wk
    use global_norms_mod
    use types_mod,      only: rk_t


    implicit none

    type (element_t)     , intent(inout), target :: elem(:)
    type (Rk_t)          , intent(in) :: MyRk
    type (EdgeBuffer_t)  , intent(in) :: edge2
    type (EdgeBuffer_t)  , intent(inout) :: edge3

    type (derivative_t)  , intent(in) :: deriv

    type (hybrid_t)      , intent(in) :: hybrid

    real (kind=real_kind), intent(in) :: dt
    real (kind=real_kind), intent(in) :: pmean
    type (TimeLevel_t)   , intent(inout) :: tl
    integer              , intent(in) :: nets
    integer              , intent(in) :: nete

    ! =================
    ! Local
    ! =================

    ! pointer ...

    real (kind=real_kind), dimension(:,:), pointer     :: fcor,rspheremp,spheremp

    ! Thread private working set ...

    real (kind=real_kind), dimension(np,np,2,nlev,nets:nete) :: vtens
    real (kind=real_kind), dimension(np,np,nlev,nets:nete)   :: ptens
    real (kind=real_kind), dimension(np,np,nlev,nets:nete)   :: ptens_dg
    real (kind=real_kind), dimension(0:np+1,0:np+1,nlev)     :: pedges

    real (kind=real_kind), dimension(np,np,2)    :: grade   ! kinetic energy gradient
    real (kind=real_kind), dimension(np,np,2)    :: gradh   ! grad(h)

    real (kind=real_kind), dimension(np,np,2)    :: pv      ! p*v lat-lon
    real (kind=real_kind), dimension(np,np)      :: E          ! kinetic energy term
    real (kind=real_kind), dimension(np,np)      :: zeta       ! relative vorticity
    real (kind=real_kind), dimension(np,np)      :: div, plocal
    real (kind=real_kind), dimension(np,np,2)    :: ulatlon

    real (kind=real_kind) :: v1,v2,fjmax(4)
    real (kind=real_kind) :: vtens1,vtens2
    real (kind=real_kind) :: pmin(nlev,nets:nete),pmax(nlev,nets:nete)
    real (kind=real_kind) :: real_time
    real (kind=real_kind) :: dtstage

    integer :: i,j,k,s,ie
    integer :: kptr
    integer :: n0,np1
    integer :: nstep
    integer :: ntmp

    logical :: Debug = .FALSE.

    if (kmass==-1) then
       if (limiter_option ==0 .or. limiter_option==4) then
          ! when we dont advect a seperate density, we only allow limiter=0 or 4
       else
          print *,'Error: limiter can only be applied for advection tests with kmass>0'
          stop 
       endif
    endif

    n0    = tl%n0
    np1   = tl%np1
    nstep = tl%nstep

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
       if (( limiter_option == 8 ).or.(limiter_option == 81 ).or.(limiter_option == 9)) then
          call neighbor_minmax(elem,hybrid,edge3,nets,nete,n0,pmin,pmax,kmass=kmass)
       endif


       if(Debug) print *,'homme: adv.._rk 1'


       ! ===================================
       ! construct v tendencies and v.grad(p)
       ! on the velocity grid...
       ! ===================================

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
                    elem(ie)%state%v(i,j,1,k,n0)=elem(ie)%D(i,j,1,1)*v1 + elem(ie)%D(i,j,1,2)*v2   ! contra->latlon
                    elem(ie)%state%v(i,j,2,k,n0)=elem(ie)%D(i,j,2,1)*v1 + elem(ie)%D(i,j,2,2)*v2   ! contra->latlon
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
                    elem(ie)%state%v(i,j,1,k,n0) = elem(ie)%Dinv(i,j,1,1)*v1 + elem(ie)%Dinv(i,j,1,2)*v2
                    elem(ie)%state%v(i,j,2,k,n0) = elem(ie)%Dinv(i,j,2,1)*v1 + elem(ie)%Dinv(i,j,2,2)*v2
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

       do ie=nets,nete
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
                   ulatlon(i,j,1)=elem(ie)%D(i,j,1,1)*v1 + elem(ie)%D(i,j,1,2)*v2   ! contra->latlon
                   ulatlon(i,j,2)=elem(ie)%D(i,j,2,1)*v1 + elem(ie)%D(i,j,2,2)*v2   ! contra->latlon

                   E(i,j) = 0.5D0*(ulatlon(i,j,1)**2 + ulatlon(i,j,2)**2)  +&
                        elem(ie)%state%p(i,j,k,n0) + elem(ie)%state%ps(i,j)

                   pv(i,j,1) = ulatlon(i,j,1)*(pmean+elem(ie)%state%p(i,j,k,n0))
                   pv(i,j,2) = ulatlon(i,j,2)*(pmean+elem(ie)%state%p(i,j,k,n0))
                end do
             end do
             grade = gradient_sphere(E,deriv,elem(ie)%Dinv)       ! scalar -> latlon vector
             !grade = gradient_sphere_wk(E,deriv,elem(ie)%Dinv)       ! scalar -> latlon vector
             zeta = vorticity_sphere(ulatlon,deriv,elem(ie)) ! latlon vector -> scalar 
             div = divergence_sphere(pv,deriv,elem(ie))      ! latlon vector -> scalar

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
             call limiter_optim_wrap(ptens(:,:,:,ie),elem(ie)%spheremp(:,:),&
                  pmin(:,ie),pmax(:,ie),kmass)
          endif

          if ((limiter_option == 9))then
             call limiter_optim_wrap_lim9(ptens(:,:,:,ie),elem(ie)%spheremp(:,:),&
                  pmin(:,ie),pmax(:,ie),kmass)
          endif

          if ((limiter_option == 81))then
             pmax(:,ie)=pmax(:,ie)+9e19  ! disable max constraint
             call limiter_optim_wrap(ptens(:,:,:,ie),elem(ie)%spheremp(:,:),&
                  pmin(:,ie),pmax(:,ie),kmass)
          endif
          
          if ((limiter_option == 84))then
	     pmin(:,ie)=0.0d0
             if (test_case=='swirl') then
                 pmin(1,ie)=.1d0
                 if (nlev>=3) then
                    k=3; pmin(k,ie)=.1d0
                 endif
             endif
             pmax(:,ie)=pmax(:,ie)+9e19  ! disable max constraint
             call limiter_optim_wrap(ptens(:,:,:,ie),elem(ie)%spheremp(:,:),&
                  pmin(:,ie),pmax(:,ie),kmass)
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
          kptr=0
          call edgeVpack(edge3, ptens(1,1,1,ie),nlev,kptr,ie)
          kptr=nlev
          call edgeVpack(edge3,vtens(1,1,1,1,ie),2*nlev,kptr,ie)

       end do

       if(Debug) print *,'homme: adv.._rk 2'

#if (defined HORIZ_OPENMP)
       !$OMP BARRIER
#endif

       call bndry_exchangeV(hybrid,edge3)

#if (defined HORIZ_OPENMP)
       !$OMP BARRIER
#endif
       
       do ie=nets,nete

          rspheremp     => elem(ie)%rspheremp

          ! ===========================================================
          ! Unpack the edges for vgradp and vtens
          ! ===========================================================
          kptr=0
          call edgeVunpack(edge3, ptens(1,1,1,ie), nlev, kptr, ie)

          kptr=nlev
          call edgeVunpack(edge3, vtens(1,1,1,1,ie), 2*nlev, kptr, ie)

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
                   vtens(i,j,1,k,ie) = elem(ie)%Dinv(i,j,1,1)*vtens1 + elem(ie)%Dinv(i,j,1,2)*vtens2
                   vtens(i,j,2,k,ie) = elem(ie)%Dinv(i,j,2,1)*vtens1 + elem(ie)%Dinv(i,j,2,2)*vtens2
                end do
             end do
          end do

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

       real_time =real_time + dtstage
!this is only for output reasons, if velocities are prescribed
      do ie=nets,nete
	  call set_prescribed_velocity(elem(ie),n0,real_time)
      enddo
    enddo ! stage loop


    ! make leap-frog compliant
    ! after RK loop u^n containe u^{n+1} and u^{n+1} contains u^n
    ! thus we swap them
    ntmp    = tl%np1
    tl%np1  = tl%n0
    tl%n0   = ntmp

#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif
  end subroutine advance_nonstag_rk

!-----------------------------------------------------------------------------------------------

  subroutine limiter_optim_wrap(ptens,sphweights,minp,maxp,kmass)
!THIS IS A NEW VERSION OF LIM8, from 3D code

    use kinds, only : real_kind
    use dimensions_mod, only : np, np, nlev
    use derivative_mod, only : limiter_optim_iter_full

    real (kind=real_kind), dimension(nlev), intent(inout)   :: minp, maxp
    real (kind=real_kind), dimension(np,np,nlev), intent(inout)   :: ptens
    real (kind=real_kind), dimension(np,np), intent(in)   :: sphweights
    integer, intent(in) :: kmass
    ! local variables
    real (kind=real_kind), dimension(np,np,nlev) :: ptens_mass
    integer :: k

    do k=1,nlev
       ptens_mass(:,:,k)=ptens(:,:,kmass) ! save mass, since 3D limiter will apply limiter to it
    enddo
    call limiter_optim_iter_full(ptens,sphweights,minp,maxp,ptens_mass)
    ptens(:,:,kmass) = ptens_mass(:,:,1)
  end subroutine


  subroutine limiter_optim_wrap_lim9(ptens,sphweights,minp,maxp,kmass)
!THIS IS A NEW VERSION OF LIM9, from 3D code

    use kinds, only : real_kind
    use dimensions_mod, only : np, np, nlev
    use derivative_mod, only : limiter_clip_and_sum

    real (kind=real_kind), dimension(nlev), intent(inout)   :: minp, maxp
    real (kind=real_kind), dimension(np,np,nlev), intent(inout)   :: ptens
    real (kind=real_kind), dimension(np,np), intent(in)   :: sphweights
    integer, intent(in) :: kmass
    ! local variables
    real (kind=real_kind), dimension(np,np,nlev) :: ptens_mass
    integer :: k

    do k=1,nlev
       ptens_mass(:,:,k)=ptens(:,:,kmass)
    enddo
    call limiter_clip_and_sum(ptens,sphweights,minp,maxp,ptens_mass)
    ptens(:,:,kmass) = ptens_mass(:,:,1)
  end subroutine limiter_optim_wrap_lim9




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

    do k=nlev,1,-1
       qmin=0
       if (test_case=='swirl') then
          if (k.eq.1) qmin=0.0d0   ! lifted cosine bell
          if (k.eq.3) qmin=0.1d0   ! lifted slotted cylinder
       endif
       if (kmass==-1) then
          call limiter2d_min_onelevel(Q(:,:,k),spheremp(:,:),qmin)
       else
          if(k.ne.kmass)then
             Q(:,:,k)=Q(:,:,k)/Q(:,:,kmass)
             call limiter2d_min_onelevel(Q(:,:,k),spheremp(:,:)*Q(:,:,kmass),qmin)
             Q(:,:,k)=Q(:,:,k)*Q(:,:,kmass)
          endif
       endif
    enddo
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
    use control_mod, only : nu, nu_s, hypervis_order, hypervis_subcycle, limiter_option,&
	  test_case, kmass
    use hybrid_mod, only : hybrid_t
    use element_mod, only : element_t
    use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk
    use edge_mod, only : edgevpack_nlyr, edgevunpack_nlyr
    use edgetype_mod, only : EdgeBuffer_t
    use bndry_mod, only : bndry_exchangev
    use viscosity_mod, only : biharmonic_wk, neighbor_minmax
    implicit none

    type (hybrid_t)      , intent(in) :: hybrid
    type (element_t)     , intent(inout), target :: elem(:)
    integer :: nt,nets,nete
    real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)  :: vtens
    real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: ptens
    type (EdgeBuffer_t)  , intent(inout) :: edge3
    type (derivative_t)  , intent(in) :: deriv
    real (kind=real_kind) :: dt2


    ! local
    integer :: k,kptr,i,j,ie,ic
    real (kind=real_kind), dimension(:,:), pointer :: spheremp,rspheremp
    real (kind=real_kind), dimension(np,np) :: ptot
    real (kind=real_kind), dimension(np,np) :: lap_p
    real (kind=real_kind), dimension(np,np,2) :: lap_v
    real (kind=real_kind) :: v1,v2,dt

    if (hypervis_order == 0 ) return;
    if (nu_s == 0 .and. nu == 0 ) return;


    do ie=nets,nete
       do k=1,nlev

          ! contra -> latlon
          do j=1,np
             do i=1,np
                v1     = elem(ie)%state%v(i,j,1,k,nt)   ! contra
                v2     = elem(ie)%state%v(i,j,2,k,nt)   ! contra 
                elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%D(i,j,1,1)*v1 + elem(ie)%D(i,j,1,2)*v2   ! contra->latlon
                elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%D(i,j,2,1)*v1 + elem(ie)%D(i,j,2,2)*v2   ! contra->latlon
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
             call edgeVpack_nlyr(edge3,elem(ie)%desc,elem(ie)%state%p(:,:,:,nt),nlev,kptr,3*nlev)
             kptr=nlev
             call edgeVpack_nlyr(edge3,elem(ie)%desc,elem(ie)%state%v(:,:,:,:,nt),2*nlev,kptr,3*nlev)
          enddo

          call bndry_exchangeV(hybrid,edge3)

          do ie=nets,nete
             rspheremp     => elem(ie)%rspheremp

             kptr=0
             call edgeVunpack_nlyr(edge3,elem(ie)%desc,elem(ie)%state%p(:,:,:,nt), nlev, kptr,3*nlev)
             kptr=nlev
             call edgeVunpack_nlyr(edge3,elem(ie)%desc,elem(ie)%state%v(:,:,:,:,nt), 2*nlev, kptr,3*nlev)

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

             kptr=0
             call edgeVpack_nlyr(edge3,elem(ie)%desc,elem(ie)%state%p(:,:,:,nt),nlev,kptr,3*nlev)
             kptr=nlev
             call edgeVpack_nlyr(edge3,elem(ie)%desc,elem(ie)%state%v(:,:,:,:,nt),2*nlev,kptr,3*nlev)
          enddo

          call bndry_exchangeV(hybrid,edge3)

          do ie=nets,nete
             rspheremp     => elem(ie)%rspheremp

             kptr=0
             call edgeVunpack_nlyr(edge3,elem(ie)%desc,elem(ie)%state%p(:,:,:,nt), nlev, kptr,3*nlev)
             kptr=nlev
             call edgeVunpack_nlyr(edge3,elem(ie)%desc,elem(ie)%state%v(:,:,:,:,nt), 2*nlev, kptr,3*nlev)

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
                elem(ie)%state%v(i,j,1,k,nt) = elem(ie)%Dinv(i,j,1,1)*v1 + elem(ie)%Dinv(i,j,1,2)*v2
                elem(ie)%state%v(i,j,2,k,nt) = elem(ie)%Dinv(i,j,2,1)*v1 + elem(ie)%Dinv(i,j,2,2)*v2
             enddo
          enddo
       enddo
    enddo


    ! lat-lon conversion 1740  20%
    ! opss+pack+mass     5937  68% 
    ! unpack+mass_inv    1015  12%     
    ! total 8692
  end subroutine advance_hypervis


!----------------------------------------------------------------------------------------


  subroutine set_prescribed_velocity(elem,n0,time)
  use control_mod, only :  test_case
  use element_mod, only : element_t
  use dimensions_mod, only : nlev
  use shallow_water_mod, only : tc1_velocity, vortex_velocity, swirl_velocity
  implicit none

  type (element_t)     , intent(inout) :: elem

  ! local
  integer :: n0,k
  real (kind=real_kind) :: time

  if (test_case=="swtc1") then
     do k=1,nlev
        elem%state%v(:,:,:,k,n0)=tc1_velocity(elem%spherep,elem%Dinv)
     end do
  else if (test_case=="vortex") then
     do k=1,nlev
        elem%state%v(:,:,:,k,n0)=vortex_velocity(time,elem%spherep,elem%Dinv)
     end do
  else if (test_case=="swirl") then
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
  use edge_mod, only : edgevpack_nlyr, edgevunpack_nlyr
  use edgetype_mod, only : EdgeBuffer_t
  use bndry_mod, only : bndry_exchangev
  use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
  implicit none

  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  integer :: nm1,np1,n0,nets,nete
  real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)  :: vtens
  real (kind=real_kind), dimension(np,np,nlev,nets:nete) :: ptens
  type (EdgeBuffer_t)  , intent(inout) :: edge3
  type (derivative_t)  , intent(in) :: deriv
  real (kind=real_kind) :: dt2,pmean,real_time

  ! local
  ! pointer ...
  real (kind=real_kind), dimension(:,:), pointer :: fcor,rspheremp,spheremp
  real (kind=real_kind), dimension(np,np,2)    :: grade   ! strong kinetic energy gradient
  real (kind=real_kind), dimension(np,np,2)    :: pv      ! p*v lat-lon
  real (kind=real_kind), dimension(np,np)                     :: E          ! kinetic energy term
  real (kind=real_kind), dimension(np,np)                     :: zeta       ! relative vorticity
  real (kind=real_kind), dimension(np,np)      :: div  
  real (kind=real_kind), dimension(np,np,2)      :: ulatlon

  integer i,j,k,kptr,ie
  real (kind=real_kind) ::  v1,v2
  real (kind=real_kind) ::  vtens1,vtens2

  call t_startf('compute_and_apply_rhs')


  ! ===================================
  ! construct v tendencies and v.grad(p)
  ! on the velocity grid...
  ! ===================================
  do ie=nets,nete
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
              ulatlon(i,j,1)=elem(ie)%D(i,j,1,1)*v1 + elem(ie)%D(i,j,1,2)*v2   ! contra->latlon
              ulatlon(i,j,2)=elem(ie)%D(i,j,2,1)*v1 + elem(ie)%D(i,j,2,2)*v2   ! contra->latlon
              
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
     call edgeVpack_nlyr(edge3,elem(ie)%desc,ptens(1,1,1,ie),nlev,kptr,3*nlev)
     kptr=nlev
     call edgeVpack_nlyr(edge3,elem(ie)%desc,vtens(1,1,1,1,ie),2*nlev,kptr,3*nlev)
  end do
  
  
#if (defined HORIZ_OPENMP)
  !$OMP BARRIER
#endif
  call bndry_exchangeV(hybrid,edge3)
#if (defined HORIZ_OPENMP)
  !$OMP BARRIER
#endif
  
  do ie=nets,nete
     rspheremp     => elem(ie)%rspheremp
     
     ! ===========================================================
     ! Unpack the edges for vgradp and vtens
     ! ===========================================================
     kptr=0
     call edgeVunpack_nlyr(edge3, elem(ie)%desc, ptens(1,1,1,ie), nlev, kptr,3*nlev)
     
     kptr=nlev
     call edgeVunpack_nlyr(edge3, elem(ie)%desc, vtens(1,1,1,1,ie), 2*nlev, kptr, 3*nlev)
     
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
              vtens(i,j,1,k,ie) = elem(ie)%Dinv(i,j,1,1)*vtens1 + elem(ie)%Dinv(i,j,1,2)*vtens2
              vtens(i,j,2,k,ie) = elem(ie)%Dinv(i,j,2,1)*vtens1 + elem(ie)%Dinv(i,j,2,2)*vtens2
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
  call t_stopf('compute_and_apply_rhs')
  end subroutine compute_and_apply_rhs
  



end module advance_mod
