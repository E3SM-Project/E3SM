#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef TRILINOS

module prim_implicit_mod
  use kinds, only : real_kind
  use edgetype_mod, only : EdgeBuffer_t
  use parallel_mod, only : abortmp, parallel_t, iam
  use perf_mod, only: t_startf, t_stopf, t_barrierf, t_adj_detailf ! _EXTERNAL
  implicit none
  private
  save

  public :: prim_implicit_init

  type(EdgeBuffer_t) :: edge3p1

contains

  subroutine prim_implicit_init(par, elem)
    use edge_mod, only : initEdgeBuffer
    use element_mod, only : element_t
    use dimensions_mod, only : nlev, nelemd
    use control_mod, only : qsplit, rsplit
    implicit none

    type (parallel_t), intent(in) :: par
    type (element_t),  intent(in), target :: elem(:)

    if (rsplit==0) then
       call initEdgeBuffer(par,edge3p1,elem,3*nlev+1)
    else
       ! need extra buffer space for dp3d
       call initEdgeBuffer(par,edge3p1,elem,4*nlev+1)
    endif

  end subroutine prim_implicit_init



  subroutine residual(xstate, fx, nelemd, c_ptr_to_object) bind(C,name='calc_f')
    ! ===================================
    ! compute the RHS, accumulate into a residual for each dependent variable and apply DSS
    !
    ! For the Backward Euler first order option 
    !
    !           F(u)  = (u(np1) - u(n0))/dt  -  DSS[ RHS(u(np1)) ]
    !
    ! For the BDF2 second order option 
    !
    !           F(u)  = (u(np1) - u(n0))/dt  -  DSS[ RHS(u(np1)) ]
    !
    ! This subroutine is normally called to compute a leapfrog timestep
    ! but by adjusting np1,nm1,n0 and dt2, many other timesteps can be
    ! accomodated.  For example, setting nm1=np1=n0 this routine will
    ! take a forward euler step, overwriting the input with the output.
    !
    !    qn0 = timelevel used to access Qdp() in order to compute virtual Temperature
    !          qn0=-1 for the dry case
    !
    ! if  dt2<0, then the DSS'd RHS is returned in timelevel np1
    !
    ! Combining the RHS and DSS pack operation in one routine
    ! allows us to fuse these two loops for more cache reuse
    !
    ! Combining the dt advance and DSS unpack operation in one routine
    ! allows us to fuse these two loops for more cache reuse
    !
    ! note: for prescribed velocity case, velocity will be computed at
    ! "real_time", which should be the time of timelevel n0.
    !
    ! ===================================
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use hybrid_mod, only : hybrid_t
    use element_mod, only : element_t
    use derivative_mod, only : derivative_t, divergence_sphere, gradient_sphere, vorticity_sphere
    use edge_mod, only : edgevpack, edgevunpack
    use bndry_mod, only : bndry_exchangev
    use control_mod, only : moisture, qsplit, use_cpstar, rsplit
    use hybvcoord_mod, only : hvcoord_t
    use time_mod, only : TimeLevel_t

    use physical_constants, only : cp, cpwater_vapor, Rgas, kappa
    use physics_mod, only : virtual_specific_heat, virtual_temperature
    use prim_si_mod, only : preq_vertadv, preq_omega_ps, preq_hydrostatic

    use prim_derived_type_mod ,only : derived_type, initialize
    use, intrinsic :: iso_c_binding

    implicit none
    integer  :: np1,nm1,n0,nm,qn0,nets,nete,lx,n
    real*8   :: dt2, dti
    logical  :: compute_diagnostics

    real (kind=real_kind) :: eta_ave_w  ! weighting for eta_dot_dpdn mean flux
    real (kind=real_kind), pointer, dimension(:,:)   :: ps         ! surface pressure for current tiime level
    real (kind=real_kind), pointer, dimension(:,:,:) :: phi        ! geopotential 

    real (kind=real_kind), dimension(np,np,nlev)   :: omega_p
    real (kind=real_kind), dimension(np,np,nlev)   :: T_v
    real (kind=real_kind), dimension(np,np,nlev)   :: divdp
    real (kind=real_kind), dimension(np,np,nlev+1) :: eta_dot_dpdn ! half level vertical velocity on p-grid
    real (kind=real_kind), dimension(np,np)        :: sdot_sum     ! temporary field
    real (kind=real_kind), dimension(np,np,2)      :: vtemp        ! generic gradient storage
    real (kind=real_kind), dimension(np,np)        :: vgrad_T      ! v.grad(T)
    real (kind=real_kind), dimension(np,np)        :: Ephi         ! kinetic energy + PHI term
    real (kind=real_kind), dimension(np,np,2)      :: grad_ps      ! lat-lon coord version
    real (kind=real_kind), dimension(np,np,2,nlev) :: grad_p
    real (kind=real_kind), dimension(np,np,nlev)   :: vort         ! vorticity
    real (kind=real_kind), dimension(np,np,nlev)   :: p            ! pressure
    real (kind=real_kind), dimension(np,np,nlev)   :: dp           ! delta pressure
    real (kind=real_kind), dimension(np,np,nlev)   :: rdp          ! inverse of delta pressure
    real (kind=real_kind), dimension(np,np,nlev)   :: T_vadv       ! temperature vertical advection
    real (kind=real_kind), dimension(np,np,nlev)   :: vgrad_p      ! v.grad(p)
    real (kind=real_kind), dimension(np,np,nlev+1) :: ph           ! half level pressures on p-grid
    real (kind=real_kind), dimension(np,np,2,nlev) :: v_vadv       ! velocity vertical advection
    real (kind=real_kind) ::  kappa_star(np,np,nlev)
    real (kind=real_kind) ::  vtens1(np,np,nlev)
    real (kind=real_kind) ::  vtens2(np,np,nlev)
    real (kind=real_kind) ::  ttens(np,np,nlev)
    real (kind=real_kind) ::  fvtens(np,np,2,nlev,nelem)
    real (kind=real_kind) ::  fttens(np,np,nlev,nelem)
    real (kind=real_kind) ::  fdptens(np,np,nlev,nelem)
    real (kind=real_kind) ::  fpstens(np,np,nelem)

    real (kind=real_kind) ::  cp2,cp_ratio,E,de,Qt,v1,v2
    real (kind=real_kind) ::  glnps1,glnps2,gpterm
    real (kind=real_kind) ::  gam
    integer :: i,j,k,kptr,ie, n_Q

    type (element_t), target     :: elem(nelem)
    type (hvcoord_t)     :: hvcoord
    type (hybrid_t)      :: hybrid
    type (derivative_t)  :: deriv
    type (TimeLevel_t)   :: tl

    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(in)        :: xstate(nelemd)
    real (c_double) ,intent(out)       :: fx(nelemd)
    integer                            :: method, nstep

    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    ! set these to match the compute_and_apply_rhs for evenutal reincorporation, so
    ! BE  0= x(np1)-x(nm1) - (RSS*x(n0))
    ! BDF2  0= (1.5) (x(np1)-x(nm1)) - (0.5) (x(nm1)-x(nm)) - (RSS*x(n0))

    np1        = fptr%tl%np1  
    nm1        = fptr%tl%n0    
    n0         = fptr%tl%np1    
    nm         = fptr%tl%nm1    
    nstep      = fptr%tl%nstep 

    method     = fptr%method
    dt2        = fptr%dt
    dti        = 1/dt2
    hvcoord    = fptr%hvcoord
    hybrid     = fptr%hybrid
    deriv      = fptr%deriv
    nets       = fptr%nets
    nete       = fptr%nete
    compute_diagnostics = fptr%compute_diagnostics
    qn0        = fptr%n_Q
    eta_ave_w  = fptr%eta_ave_w

    if ((method==11).or.(nstep==0)) then
       gam=0.0  ! use BE for BDF as well for the first time step
    else if (method==12) then
       gam=0.5  ! BDF2 
    else ! for adding more methods, right now default to BE
       call abortmp('ERROR: method > 12 not yet coded for implicit time integration')
    end if

    !JMD  call t_barrierf('sync_residual', hybrid%par%comm)
    call t_adj_detailf(+1)
    call t_startf('residual')

    call t_startf('residual_int')
    !fvtens = 0.0d0
    !fttens = 0.0d0
    !fpstens = 0.0d0
    !fdptens = 0.0d0

    lx = 1
    do ie=nets,nete
       do k=1,nlev
          do j=1,np
             do i=1,np
                fptr%base(ie)%state%v(i,j,1,k,np1) = xstate(lx)
                lx = lx+1
             end do
          end do
       end do
    end do
    do ie=nets,nete
       do k=1,nlev
          do j=1,np
             do i=1,np
                fptr%base(ie)%state%v(i,j,2,k,np1) = xstate(lx)
                lx = lx+1
             end do
          end do
       end do
    end do
    do ie=nets,nete
       do k=1,nlev
          do j=1,np
             do i=1,np
                fptr%base(ie)%state%T(i,j,k,np1) = xstate(lx)
                lx = lx+1
             end do
          end do
       end do
    end do
    do ie=nets,nete
       do j=1,np
          do i=1,np
             fptr%base(ie)%state%ps_v(i,j,np1) = xstate(lx)
             lx = lx+1
          end do
       end do
    end do
    call t_stopf('residual_int')

    call t_startf('residual_cal')

    do ie=nets,nete
       !     ps => fptr%base(ie)%state%ps_v(:,:,n0)
       phi => fptr%base(ie)%derived%phi(:,:,:)

       ! ==================================================
       ! compute pressure (p) on half levels from ps
       ! using the hybrid coordinates relationship, i.e.
       ! e.g. equation (3.a.92) of the CCM-2 description,
       ! (NCAR/TN-382+STR), June 1993, p. 24.
       ! ==================================================
       ! vertically eulerian only needs grad(ps)
       if (rsplit==0) &
            grad_ps = gradient_sphere(fptr%base(ie)%state%ps_v(:,:,n0),deriv,fptr%base(ie)%Dinv)

       ! ============================
       ! compute p and delta p
       ! ============================
       !#if (defined COLUMN_OPENMP)
!!$omp parallel do private(k,i,j)
       !#endif
       !     do k=1,nlev+1
       !       ph(:,:,k)   = hvcoord%hyai(k)*hvcoord%ps0 + hvcoord%hybi(k)*fptr%base(ie)%state%ps_v(:,:,n0)
       !     end do

#if (defined COLUMN_OPENMP)
       !$omp parallel do private(k,i,j,v1,v2,vtemp,dp,p,grad_p,rdp,divdp,vort)
#endif
       do k=1,nlev
          if (rsplit==0) then
             dp(:,:,k) = (hvcoord%hyai(k+1)*hvcoord%ps0 + hvcoord%hybi(k+1)*fptr%base(ie)%state%ps_v(:,:,n0)) &
                  - (hvcoord%hyai(k)*hvcoord%ps0 + hvcoord%hybi(k)*fptr%base(ie)%state%ps_v(:,:,n0))
             p(:,:,k)   = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*fptr%base(ie)%state%ps_v(:,:,n0)
             grad_p(:,:,:,k) = hvcoord%hybm(k)*grad_ps(:,:,:)

          else
             ! vertically lagrangian code: we advect dp3d instead of ps_v
             ! we also need grad(p) at all levels (not just grad(ps))
             !p(k)= hyam(k)*ps0 + hybm(k)*ps
             !    = .5*(hyai(k+1)+hyai(k))*ps0 + .5*(hybi(k+1)+hybi(k))*ps
             !    = .5*(ph(k+1) + ph(k) )  = ph(k) + dp(k)/2
             !
             ! p(k+1)-p(k) = ph(k+1)-ph(k) + (dp(k+1)-dp(k))/2
             !             = dp(k) + (dp(k+1)-dp(k))/2 = (dp(k+1)+dp(k))/2
             dp(:,:,k) = fptr%base(ie)%state%dp3d(:,:,k,n0)
             if (k==1) then
                p(:,:,k)=hvcoord%hyai(k)*hvcoord%ps0 + dp(:,:,k)/2
             else
                p(:,:,k)=p(:,:,k-1) + dp(:,:,k-1)/2 + dp(:,:,k)/2
             endif
             grad_p(:,:,:,k) = gradient_sphere(p(:,:,k),deriv,fptr%base(ie)%Dinv)
          endif

          rdp(:,:,k) = 1.0D0/dp(:,:,k)

          ! ============================
          ! compute vgrad_lnps
          ! ============================
          do j=1,np
             do i=1,np
                v1 = fptr%base(ie)%state%v(i,j,1,k,n0)
                v2 = fptr%base(ie)%state%v(i,j,2,k,n0)
                !              vgrad_p(i,j,k) = &
                !                   hvcoord%hybm(k)*(v1*grad_ps(i,j,1) + v2*grad_ps(i,j,2))
                vgrad_p(i,j,k) = (v1*grad_p(i,j,1,k) + v2*grad_p(i,j,2,k))
                vtemp(i,j,1) = v1*dp(i,j,k)
                vtemp(i,j,2) = v2*dp(i,j,k)
             end do
          end do

          ! ================================
          ! Accumulate mean Vel_rho flux in vn0
          ! ================================
          fptr%base(ie)%derived%vn0(:,:,:,k)=fptr%base(ie)%derived%vn0(:,:,:,k)+eta_ave_w*vtemp(:,:,:)


          ! =========================================
          !
          ! Compute relative vorticity and divergence
          !
          ! =========================================
          divdp(:,:,k)=divergence_sphere(vtemp,deriv,fptr%base(ie))
          vort(:,:,k)=vorticity_sphere(fptr%base(ie)%state%v(:,:,:,k,n0),deriv,fptr%base(ie))

       enddo

       ! compute T_v for timelevel n0
       !if ( moisture /= "dry") then
       if (qn0 == -1 ) then
#if (defined COLUMN_OPENMP)
          !$omp parallel do private(k,i,j)
#endif
          do k=1,nlev
             do j=1,np
                do i=1,np
                   T_v(i,j,k) = fptr%base(ie)%state%T(i,j,k,n0)
                   kappa_star(i,j,k) = kappa
                end do
             end do
          end do
       else
#if (defined COLUMN_OPENMP)
          !$omp parallel do private(k,i,j,Qt)
#endif
          do k=1,nlev
             do j=1,np
                do i=1,np
                   ! Qt = fptr%base(ie)%state%Q(i,j,k,1)
                   Qt = fptr%base(ie)%state%Qdp(i,j,k,1,qn0)/dp(i,j,k)
                   T_v(i,j,k) = Virtual_Temperature(fptr%base(ie)%state%T(i,j,k,n0),Qt)
                   if (use_cpstar==1) then
                      kappa_star(i,j,k) =  Rgas/Virtual_Specific_Heat(Qt)
                   else
                      kappa_star(i,j,k) = kappa
                   endif
                end do
             end do
          end do
       end if

       ! ====================================================
       ! Compute Hydrostatic equation, modeld after CCM-3
       ! ====================================================
       !call geopotential_t(p,dp,T_v,Rgas,fptr%base(ie)%derived%phi)
       call preq_hydrostatic(fptr%base(ie)%derived%phi,fptr%base(ie)%state%phis,T_v,p,dp)

       ! ====================================================
       ! Compute omega_p according to CCM-3
       ! ====================================================
       call preq_omega_ps(omega_p,hvcoord,p,vgrad_p,divdp)

       ! ==================================================
       ! zero partial sum for accumulating sum
       !    (div(v_k) + v_k.grad(lnps))*dsigma_k = div( v dp )
       ! used by eta_dot_dpdn and lnps tendency
       ! ==================================================
       sdot_sum=0

       ! ==================================================
       ! Compute eta_dot_dpdn
       ! save sdot_sum as this is the -RHS of ps_v equation
       ! ==================================================
       if (rsplit>0) then
          ! VERTICALLY LAGRANGIAN:   no vertical motion
          eta_dot_dpdn=0
          T_vadv=0
          v_vadv=0
       else
          do k=1,nlev
             ! ==================================================
             ! add this term to PS equation so we exactly conserve dry mass
             ! ==================================================
             sdot_sum(:,:) = sdot_sum(:,:) + divdp(:,:,k)
             eta_dot_dpdn(:,:,k+1) = sdot_sum(:,:)
          end do

          ! ===========================================================
          ! at this point, eta_dot_dpdn contains integral_etatop^eta[ divdp ]
          ! compute at interfaces:
          !    eta_dot_dpdn = -dp/dt - integral_etatop^eta[ divdp ]
          ! for reference: at mid layers we have:
          !    omega = v grad p  - integral_etatop^eta[ divdp ]
          ! ===========================================================
#if (defined COLUMN_OPENMP)
          !$omp parallel do private(k)
#endif
          do k=1,nlev-1
             eta_dot_dpdn(:,:,k+1) = hvcoord%hybi(k+1)*sdot_sum(:,:) - eta_dot_dpdn(:,:,k+1)
          end do

          eta_dot_dpdn(:,:,1     ) = 0.0D0
          eta_dot_dpdn(:,:,nlev+1) = 0.0D0

          ! ===========================================================
          ! Compute vertical advection of T and v from eq. CCM2 (3.b.1)
          ! ==============================================
          call preq_vertadv(fptr%base(ie)%state%T(:,:,:,n0),fptr%base(ie)%state%v(:,:,:,:,n0), &
               eta_dot_dpdn,rdp,T_vadv,v_vadv)
       endif

       ! ================================
       ! accumulate mean vertical flux:
       ! ================================
#if (defined COLUMN_OPENMP)
       !$omp parallel do private(k)
#endif
       do k=1,nlev  !  Loop index added (AAM)
          fptr%base(ie)%derived%eta_dot_dpdn(:,:,k) = &
               fptr%base(ie)%derived%eta_dot_dpdn(:,:,k) + eta_ave_w*eta_dot_dpdn(:,:,k)
          fptr%base(ie)%derived%omega_p(:,:,k) = &
               fptr%base(ie)%derived%omega_p(:,:,k) + eta_ave_w*omega_p(:,:,k)
       enddo
       fptr%base(ie)%derived%eta_dot_dpdn(:,:,nlev+1) = &
            fptr%base(ie)%derived%eta_dot_dpdn(:,:,nlev+1) + eta_ave_w*eta_dot_dpdn(:,:,nlev+1)


       ! ==============================================
       ! Compute phi + kinetic energy term: 10*np*np Flops
       ! ==============================================
#if (defined COLUMN_OPENMP)
       !$omp parallel do private(k,i,j,v1,v2,E,Ephi,vtemp,vgrad_T,gpterm,glnps1,glnps2)
#endif
       do k=1,nlev
          do j=1,np
             do i=1,np
                v1     = fptr%base(ie)%state%v(i,j,1,k,n0)
                v2     = fptr%base(ie)%state%v(i,j,2,k,n0)
                E = 0.5D0*( v1*v1 + v2*v2 )
                Ephi(i,j)=E+fptr%base(ie)%derived%phi(i,j,k)+fptr%base(ie)%derived%pecnd(i,j,k)
             end do
          end do
          ! ================================================
          ! compute gradp term (ps/p)*(dp/dps)*T
          ! ================================================
          vtemp(:,:,:)   = gradient_sphere(fptr%base(ie)%state%T(:,:,k,n0),deriv,fptr%base(ie)%Dinv)
          do j=1,np
             do i=1,np
                v1     = fptr%base(ie)%state%v(i,j,1,k,n0)
                v2     = fptr%base(ie)%state%v(i,j,2,k,n0)
                vgrad_T(i,j) =  v1*vtemp(i,j,1) + v2*vtemp(i,j,2)
             end do
          end do


          ! vtemp = grad ( E + PHI )
          vtemp = gradient_sphere(Ephi(:,:),deriv,fptr%base(ie)%Dinv)

          do j=1,np
             do i=1,np
                !              gpterm = hvcoord%hybm(k)*T_v(i,j,k)/p(i,j,k)
                !              glnps1 = Rgas*gpterm*grad_ps(i,j,1)
                !              glnps2 = Rgas*gpterm*grad_ps(i,j,2)
                gpterm = T_v(i,j,k)/p(i,j,k)
                glnps1 = Rgas*gpterm*grad_p(i,j,1,k)
                glnps2 = Rgas*gpterm*grad_p(i,j,2,k)

                v1     = fptr%base(ie)%state%v(i,j,1,k,n0)
                v2     = fptr%base(ie)%state%v(i,j,2,k,n0)

                vtens1(i,j,k) =   - v_vadv(i,j,1,k)                           &
                     + v2*(fptr%base(ie)%fcor(i,j) + vort(i,j,k))        &
                     - vtemp(i,j,1) - glnps1

                vtens2(i,j,k) =   - v_vadv(i,j,2,k)                            &
                     - v1*(fptr%base(ie)%fcor(i,j) + vort(i,j,k))        &
                     - vtemp(i,j,2) - glnps2

                ttens(i,j,k)  = - T_vadv(i,j,k) - vgrad_T(i,j) + kappa_star(i,j,k)*T_v(i,j,k)*omega_p(i,j,k)

             end do
          end do

       end do

#ifdef ENERGY_DIAGNOSTICS
       ! =========================================================
       !
       ! diagnostics
       ! recomputes some gradients that were not saved above
       ! uses:  sdot_sum(), eta_dot_dpdn(), grad_ps()
       ! grad_phi(), dp(), p(), T_vadv(), v_vadv(), divdp()
       ! =========================================================

       ! =========================================================
       ! (AAM) - This section has accumulations over vertical levels.
       !   Be careful if implementing OpenMP
       ! =========================================================

       if (compute_diagnostics) then
          fptr%base(ie)%accum%KEhorz1=0
          fptr%base(ie)%accum%KEhorz2=0
          fptr%base(ie)%accum%IEhorz1=0
          fptr%base(ie)%accum%IEhorz2=0
          fptr%base(ie)%accum%IEhorz1_wet=0
          fptr%base(ie)%accum%IEhorz2_wet=0
          fptr%base(ie)%accum%KEvert1=0
          fptr%base(ie)%accum%KEvert2=0
          fptr%base(ie)%accum%IEvert1=0
          fptr%base(ie)%accum%IEvert2=0
          fptr%base(ie)%accum%IEvert1_wet=0
          fptr%base(ie)%accum%IEvert2_wet=0
          fptr%base(ie)%accum%T1=0
          fptr%base(ie)%accum%T2=0
          fptr%base(ie)%accum%T2_s=0
          fptr%base(ie)%accum%S1=0
          fptr%base(ie)%accum%S1_wet=0
          fptr%base(ie)%accum%S2=0

          do j=1,np
             do i=1,np
                fptr%base(ie)%accum%S2(i,j) = fptr%base(ie)%accum%S2(i,j) - &
                     sdot_sum(i,j)*fptr%base(ie)%state%phis(i,j)
             enddo
          enddo

          do k=1,nlev
             ! vtemp = grad_E(:,:,k)
             do j=1,np
                do i=1,np
                   v1     = fptr%base(ie)%state%v(i,j,1,k,n0)
                   v2     = fptr%base(ie)%state%v(i,j,2,k,n0)
                   Ephi(i,j)=0.5D0*( v1*v1 + v2*v2 )
                enddo
             enddo
             vtemp = gradient_sphere(Ephi,deriv,fptr%base(ie)%Dinv)
             do j=1,np
                do i=1,np
                   ! dp/dn u dot grad(E)
                   v1     = fptr%base(ie)%state%v(i,j,1,k,n0)
                   v2     = fptr%base(ie)%state%v(i,j,2,k,n0)
                   fptr%base(ie)%accum%KEhorz2(i,j) = fptr%base(ie)%accum%KEhorz2(i,j) + &
                        (v1*vtemp(i,j,1)  + v2*vtemp(i,j,2))*dp(i,j,k)
                   ! E div( u dp/dn )
                   fptr%base(ie)%accum%KEhorz1(i,j) = fptr%base(ie)%accum%KEhorz1(i,j) + Ephi(i,j)*divdp(i,j,k)

                   ! Cp T div( u dp/dn)   ! dry horizontal advection component
                   fptr%base(ie)%accum%IEhorz1(i,j) = fptr%base(ie)%accum%IEhorz1(i,j) + Cp*fptr%base(ie)%state%T(i,j,k,n0)*divdp(i,j,k)


                enddo
             enddo


             ! vtemp = grad_phi(:,:,k)
             vtemp = gradient_sphere(phi(:,:,k),deriv,fptr%base(ie)%Dinv)
             do j=1,np
                do i=1,np
                   v1     = fptr%base(ie)%state%v(i,j,1,k,n0)
                   v2     = fptr%base(ie)%state%v(i,j,2,k,n0)
                   E = 0.5D0*( v1*v1 + v2*v2 )
                   ! NOTE:  Cp_star = Cp + (Cpv-Cp)*q
                   ! advection terms can thus be broken into two components: dry and wet
                   ! dry components cancel exactly
                   ! wet components should cancel exactly
                   !
                   ! some diagnostics
                   ! e = eta_dot_dpdn()
                   de =  eta_dot_dpdn(i,j,k+1)-eta_dot_dpdn(i,j,k)
                   ! Cp T de/dn, integral dn:
                   fptr%base(ie)%accum%IEvert1(i,j)=fptr%base(ie)%accum%IEvert1(i,j) + Cp*fptr%base(ie)%state%T(i,j,k,n0)*de
                   ! E de/dn
                   fptr%base(ie)%accum%KEvert1(i,j)=fptr%base(ie)%accum%KEvert1(i,j) + E*de
                   ! Cp T_vadv dp/dn
                   fptr%base(ie)%accum%IEvert2(i,j)=fptr%base(ie)%accum%IEvert2(i,j) + Cp*T_vadv(i,j,k)*dp(i,j,k)
                   ! dp/dn V dot V_vadv
                   fptr%base(ie)%accum%KEvert2(i,j)=fptr%base(ie)%accum%KEvert2(i,j) + (v1*v_vadv(i,j,1,k) + v2*v_vadv(i,j,2,k)) *dp(i,j,k)

                   ! IEvert1_wet():  (Cpv-Cp) T Qdp_vadv  (Q equation)
                   ! IEvert2_wet():  (Cpv-Cp) Qdp T_vadv   T equation
                   if (use_cpstar==1) then
                      fptr%base(ie)%accum%IEvert2_wet(i,j)=fptr%base(ie)%accum%IEvert2_wet(i,j) +&
                           (Cpwater_vapor-Cp)*fptr%base(ie)%state%Q(i,j,k,1)*T_vadv(i,j,k)*dp(i,j,k)
                   endif

                   gpterm = T_v(i,j,k)/p(i,j,k)
                   fptr%base(ie)%accum%T1(i,j) = fptr%base(ie)%accum%T1(i,j) - &
                        Rgas*gpterm*(grad_p(i,j,1,k)*v1 + grad_p(i,j,2,k)*v2)*dp(i,j,k)

                   fptr%base(ie)%accum%T2(i,j) = fptr%base(ie)%accum%T2(i,j) - &
                        (vtemp(i,j,1)*v1 + vtemp(i,j,2)*v2)*dp(i,j,k)

                   ! S1 = < Cp_star dp/dn , RT omega_p/cp_star >
                   fptr%base(ie)%accum%S1(i,j) = fptr%base(ie)%accum%S1(i,j) + &
                        Rgas*T_v(i,j,k)*omega_p(i,j,k)*dp(i,j,k)

                   ! cp_star = cp + cp2
                   if (use_cpstar==1) then
                      cp2 = (Cpwater_vapor-Cp)*fptr%base(ie)%state%Q(i,j,k,1)
                      cp_ratio = cp2/(cp+cp2)
                      fptr%base(ie)%accum%S1_wet(i,j) = fptr%base(ie)%accum%S1_wet(i,j) + &
                           cp_ratio*(Rgas*T_v(i,j,k)*omega_p(i,j,k)*dp(i,j,k))
                   endif

                   fptr%base(ie)%accum%CONV(i,j,:,k)=-Rgas*gpterm*grad_p(i,j,:,k)-vtemp(i,j,:)
                enddo
             enddo

             vtemp(:,:,:) = gradient_sphere(fptr%base(ie)%state%phis(:,:),deriv,fptr%base(ie)%Dinv)
             do j=1,np
                do i=1,np
                   v1     = fptr%base(ie)%state%v(i,j,1,k,n0)
                   v2     = fptr%base(ie)%state%v(i,j,2,k,n0)
                   fptr%base(ie)%accum%T2_s(i,j) = fptr%base(ie)%accum%T2_s(i,j) - &
                        (vtemp(i,j,1)*v1 + vtemp(i,j,2)*v2)*dp(i,j,k)
                enddo
             enddo

             vtemp(:,:,:)   = gradient_sphere(fptr%base(ie)%state%T(:,:,k,n0),deriv,fptr%base(ie)%Dinv)
             do j=1,np
                do i=1,np
                   v1     = fptr%base(ie)%state%v(i,j,1,k,n0)
                   v2     = fptr%base(ie)%state%v(i,j,2,k,n0)

                   ! Cp dp/dn u dot gradT
                   fptr%base(ie)%accum%IEhorz2(i,j) = fptr%base(ie)%accum%IEhorz2(i,j) + &
                        Cp*(v1*vtemp(i,j,1) + v2*vtemp(i,j,2))*dp(i,j,k)

                   if (use_cpstar==1) then
                      fptr%base(ie)%accum%IEhorz2_wet(i,j) = fptr%base(ie)%accum%IEhorz2_wet(i,j) + &
                           (Cpwater_vapor-Cp)*fptr%base(ie)%state%Q(i,j,k,1)*&
                           (v1*vtemp(i,j,1) + v2*vtemp(i,j,2))*dp(i,j,k)
                   endif

                enddo
             enddo

          enddo
       endif
#endif

       ! =========================================================
       ! local element residual calc, store in np1.
       ! apply mass matrix
       ! =========================================================

#if (defined COLUMN_OPENMP)
       !$omp parallel do private(k)
#endif

       ! for BE, gam = 0. For BDF2, gam = 0.5
       do k=1,nlev
          fttens(:,:,k,ie) = fptr%base(ie)%spheremp(:,:)* &
               ( (1+gam)*(fptr%base(ie)%state%T(:,:,k,np1) - fptr%base(ie)%state%T(:,:,k,nm1))*dti - &
               gam*(fptr%base(ie)%state%T(:,:,k,nm1) - fptr%base(ie)%state%T(:,:,k,nm))*dti - &
               ttens(:,:,k) )

          fvtens(:,:,1,k,ie) = fptr%base(ie)%spheremp(:,:)* &
               ( (1+gam)*(fptr%base(ie)%state%v(:,:,1,k,np1) - fptr%base(ie)%state%v(:,:,1,k,nm1))*dti - &
               gam*(fptr%base(ie)%state%v(:,:,1,k,nm1) - fptr%base(ie)%state%v(:,:,1,k,nm))*dti - &
               vtens1(:,:,k) )

          fvtens(:,:,2,k,ie) = fptr%base(ie)%spheremp(:,:)* &
               ( (1+gam)*(fptr%base(ie)%state%v(:,:,2,k,np1) - fptr%base(ie)%state%v(:,:,2,k,nm1))*dti - &
               gam*(fptr%base(ie)%state%v(:,:,2,k,nm1) - fptr%base(ie)%state%v(:,:,2,k,nm))*dti - &
               vtens2(:,:,k) )
       enddo
       fpstens(:,:,ie) = fptr%base(ie)%spheremp(:,:)* &
            ( (1+gam)*(fptr%base(ie)%state%ps_v(:,:,np1) - fptr%base(ie)%state%ps_v(:,:,nm1))*dti - &
            gam*(fptr%base(ie)%state%ps_v(:,:,nm1) - fptr%base(ie)%state%ps_v(:,:,nm))*dti + &
            sdot_sum(:,:) )

       ! add for fdptens later for vert lagrangian

       ! =========================================================
       !
       ! Pack ps(np1), T, and v tendencies into comm buffer
       !
       ! =========================================================
       kptr=0
       call edgeVpack(edge3p1, fpstens(:,:,ie),1,kptr,ie)

       kptr=1
       call edgeVpack(edge3p1, fttens(:,:,:,ie),nlev,kptr,ie)

       kptr=nlev+1
       call edgeVpack(edge3p1, fvtens(:,:,:,:,ie),2*nlev,kptr,ie)

       if (rsplit>0) then
          kptr=kptr+2*nlev
          call edgeVpack(edge3p1, fdptens(:,:,:,ie),nlev,kptr,ie)
       endif
    end do

    ! =============================================================
    ! Insert communications here: for shared memory, just a single
    ! sync is required
    ! =============================================================

    call bndry_exchangeV(hybrid,edge3p1)

    do ie=nets,nete
       ! ===========================================================
       ! Unpack the edges for vgrad_T and v tendencies...
       ! ===========================================================
       kptr=0
       call edgeVunpack(edge3p1, fpstens(:,:,ie), 1, kptr, ie)

       kptr=1
       call edgeVunpack(edge3p1, fttens(:,:,:,ie), nlev, kptr, ie)

       kptr=nlev+1
       call edgeVunpack(edge3p1, fvtens(:,:,:,:,ie), 2*nlev, kptr, ie)

       if (rsplit>0) then
          kptr=kptr+2*nlev
          call edgeVunpack(edge3p1, fdptens(:,:,:,ie),nlev,kptr, ie)
       endif

       ! ====================================================
       ! Scale tendencies by inverse mass matrix
       ! ====================================================

#if (defined COLUMN_OPENMP)
       !$omp parallel do private(k)
#endif
       do k=1,nlev
          fttens(:,:,k,ie)   = fptr%base(ie)%rspheremp(:,:)*fttens(:,:,k,ie)
          fvtens(:,:,1,k,ie) = fptr%base(ie)%rspheremp(:,:)*fvtens(:,:,1,k,ie)
          fvtens(:,:,2,k,ie) = fptr%base(ie)%rspheremp(:,:)*fvtens(:,:,2,k,ie)
       end do

       if (rsplit>0) then
          ! vertically lagrangian: complete dp3d timestep:
          do k=1,nlev
             fdptens(:,:,k,ie)= fptr%base(ie)%rspheremp(:,:)*fdptens(:,:,k,ie)
          enddo
          ! when debugging: also update ps_v
          fpstens(:,:,ie) = fptr%base(ie)%rspheremp(:,:)*fpstens(:,:,ie)
       else
          ! vertically eulerian: complete ps_v timestep:
          fpstens(:,:,ie) = fptr%base(ie)%rspheremp(:,:)*fpstens(:,:,ie)
       endif

    end do

    if (rsplit==0) then
       ! forward-in-time, maybe hypervis applied to PS
       call rhs_advance_hypervis(edge3p1,fvtens,fttens,fpstens,fptr%base,hvcoord,hybrid,deriv,np1,nets,nete,dt2,eta_ave_w)
    else
       ! forward-in-time, hypervis applied to dp3d
       call rhs_advance_hypervis_dp(edge3p1,fvtens,fttens,fdptens,fptr%base,hvcoord,hybrid,deriv,np1,nets,nete,dt2,eta_ave_w)
    endif


    call t_stopf('residual_cal')
    call t_startf('residual_fin')
    lx = 1
    do ie=nets,nete
       do k=1,nlev
          do j=1,np
             do i=1,np
                fx(lx) = fvtens(i,j,1,k,ie)
                lx = lx+1
             end do
          end do
       end do
    end do
    do ie=nets,nete
       do k=1,nlev
          do j=1,np
             do i=1,np
                fx(lx) = fvtens(i,j,2,k,ie)
                lx = lx+1
             end do
          end do
       end do
    end do
    do ie=nets,nete
       do k=1,nlev
          do j=1,np
             do i=1,np
                fx(lx) = fttens(i,j,k,ie)
                lx = lx+1
             end do
          end do
       end do
    end do
    do ie=nets,nete
       do j=1,np
          do i=1,np
             fx(lx) = fpstens(i,j,ie)
             lx = lx+1
          end do
       end do
    end do
    call t_stopf('residual_fin')

    !       if (hybrid%masterthread) print*, "F(u,v,t,p)",NORM2(fvtens),NORM2(fttens),NORM2(fpstens)

#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif
#endif
    call t_stopf('residual')
    call t_adj_detailf(-1)

  end subroutine residual



  ! placeholder for precon routines
  subroutine update_prec_state(xs, nelemd, c_ptr_to_object) bind(C,name='update_prec_state')

    use, intrinsic :: iso_c_binding
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use prim_derived_type_mod ,only : derived_type, initialize
    use perf_mod, only : t_startf, t_stopf
    use hybrid_mod, only : hybrid_t

    implicit none

    real (c_double) ,intent(in)        :: xs(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer              :: ns
    integer              :: ne

    integer    :: i,j,k,n,ie
    integer    :: nm1,n0,np1,lx
    type (hybrid_t)      :: hybrid

    !    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr
    !
    !    np1        = fptr%ntl1
    !    n0         = fptr%ntl1
    !    nm1        = fptr%ntl2
    !    ns         = fptr%nets
    !    ne         = fptr%nete
    !
    !       lx = 1
    !       do ie=ns,ne
    !        if (n.le.3) then
    !        do k=1,nlev
    !          do j=1,np
    !            do i=1,np
    !             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
    !             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx)
    !             if (n==3) fptr%base(ie)%state%T(i,j,k,np1)  = xs(lx)
    !             lx = lx+1
    !            end do  !np
    !          end do  !np
    !        end do  !nlev
    !        else
    !          do j=1,np
    !            do i=1,np
    !              if (n==4) fptr%base(ie)%state%ps_v(i,j,np1)  = xs(lx)
    !              lx = lx+1
    !            end do  !np
    !          end do  !np
    !        end if ! nvar
    !       end do !ie

  end subroutine update_prec_state



  subroutine test_id(xs, nelemd, fx, c_ptr_to_object) bind(C,name='test_id')
    use ,intrinsic :: iso_c_binding
    use kinds, only : real_kind
    use prim_derived_type_mod ,only : derived_type, initialize
    use perf_mod, only : t_startf, t_stopf
    use hybrid_mod, only : hybrid_t

    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(in)        :: xs(nelemd)
    real (c_double) ,intent(out)       :: fx(nelemd)
    !type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object

    !integer    :: i,j,k,n,ie
    !type (hybrid_t)      :: hybrid

    !call flush(6)
    !write(6,*)'test id'
    !call flush(6)

    !    call t_startf('test id')
    !call flush(6)
    !write(6,*)'startf'
    !call flush(6)


    !    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr
    !call flush(6)
    !write(6,*)'cf pointer'
    !call flush(6)

    !write(6,*)'setting'
    !call flush(6)
    fx=xs

    !    call t_stopf('test id')

    !write(6,*)'returning'
    !call flush(6)
  end subroutine test_id



  subroutine rhs_advance_hypervis(edge3,fv,ft,fp,elem,hvcoord,hybrid,deriv,nt,nets,nete,dt2,eta_ave_w)
    !
    !  take one timestep of:
    !          fv(:,:,2,:,:) = fv(:,:,2,:,:) +  nu*laplacian**order ( v )
    !          ft(:,:,:,:) = ft(:,:,:,:) +  nu_s*laplacian**order ( T )
    !          fp(:,:,:,:) = fp(:,:,:,:) +  nu_s*laplacian**order ( T )
    !
    use dimensions_mod, only : np, np, nlev
    use control_mod, only : nu, nu_div, nu_s, hypervis_order, hypervis_subcycle, nu_p, nu_top, psurf_vis
    use hybrid_mod, only : hybrid_t
    use hybvcoord_mod, only : hvcoord_t
    use element_mod, only : element_t
    use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk
    use edge_mod, only : edgevpack, edgevunpack
    use edgetype_mod, only : EdgeBuffer_t
    use bndry_mod, only : bndry_exchangev
    use viscosity_mod, only : biharmonic_wk
    use physical_constants, only: Cp
    !  use time_mod, only : TimeLevel_t
    implicit none

    type (hybrid_t)      , intent(in) :: hybrid
    type (element_t)     , intent(inout), target :: elem(:)
    type (EdgeBuffer_t)  , intent(inout) :: edge3
    type (derivative_t)  , intent(in) :: deriv
    type (hvcoord_t), intent(in)      :: hvcoord
    !  type (TimeLevel_t)   , intent(in) :: tl
    real (kind=real_kind), intent(inout) :: fv(:,:,:,:,:)
    real (kind=real_kind), intent(inout) :: ft(:,:,:,:)
    real (kind=real_kind), intent(inout) :: fp(:,:,:)

    real (kind=real_kind) :: dt2
    integer :: nets,nete

    ! local
    real (kind=real_kind) :: eta_ave_w  ! weighting for mean flux terms
    real (kind=real_kind) :: nu_scale, dpdn,dpdn0, nu_scale_top
    integer :: k,kptr,i,j,ie,ic,nt
    real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)      :: vtens
    real (kind=real_kind), dimension(np,np,nlev,nets:nete)        :: ptens
    real (kind=real_kind), dimension(np,np,nets:nete) :: pstens
    real (kind=real_kind), dimension(np,np,nlev) :: p
    real (kind=real_kind), dimension(np,np) :: dptemp1,dptemp2


    ! NOTE: PGI compiler bug: when using spheremp, rspheremp and ps as pointers to elem(ie)% members,
    !       data is incorrect (offset by a few numbers actually)
    !       removed for now.
    !       real (kind=real_kind), dimension(:,:), pointer :: spheremp,rspheremp
    !       real (kind=real_kind), dimension(:,:,:), pointer   :: ps

    real (kind=real_kind), dimension(np,np) :: lap_p
    real (kind=real_kind), dimension(np,np,2) :: lap_v
    real (kind=real_kind) :: v1,v2,dt,heating,utens_tmp,vtens_tmp,ptens_tmp


    if (nu_s == 0 .and. nu == 0 .and. nu_p==0 ) return;
    !JMD  call t_barrierf('sync_advance_hypervis', hybrid%par%comm)
    call t_adj_detailf(+1)
    call t_startf('rhs_advance_hypervis')

    dt=dt2

    !took out regluar viscosity option

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  hyper viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! nu_p=0:
    !   scale T dissipaton by dp  (conserve IE, dissipate T^2)
    ! nu_p>0
    !   dont scale:  T equation IE dissipation matches (to truncation error)
    !                IE dissipation from continuity equation
    !                (1 deg: to about 0.1 W/m^2)
    !

    if (hypervis_order == 2) then
       do ic=1,hypervis_subcycle
          call biharmonic_wk(elem,pstens,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete)
          do ie=nets,nete

             ! comptue mean flux
             if (nu_p>0) then
#if 0
                elem(ie)%derived%psdiss_ave(:,:)=&
                     elem(ie)%derived%psdiss_ave(:,:)+eta_ave_w*elem(ie)%state%ps_v(:,:,nt)/hypervis_subcycle
                elem(ie)%derived%psdiss_biharmonic(:,:)=&
                     elem(ie)%derived%psdiss_biharmonic(:,:)+eta_ave_w*pstens(:,:,ie)/hypervis_subcycle
#else
                do k=1,nlev
                   dptemp1(:,:) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                        ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,nt)
                   elem(ie)%derived%dpdiss_ave(:,:,k)=elem(ie)%derived%dpdiss_ave(:,:,k)+eta_ave_w*dptemp1(:,:)/hypervis_subcycle

                   dptemp2(:,:) = (hvcoord%hybi(k+1)-hvcoord%hybi(k))*pstens(:,:,ie)
                   elem(ie)%derived%dpdiss_biharmonic(:,:,k)=&
                        elem(ie)%derived%dpdiss_biharmonic(:,:,k)+eta_ave_w*dptemp2(:,:)/hypervis_subcycle
                enddo
#endif
             endif
             nu_scale=1
#if (defined COLUMN_OPENMP)
             !$omp parallel do private(k,i,j,lap_p,lap_v,nu_scale_top,dpdn,dpdn0,nu_scale,utens_tmp,vtens_tmp,ptens_tmp)
#endif
             do k=1,nlev
                ! advance in time.
                ! note: DSS commutes with time stepping, so we can time advance and then DSS.
                ! note: weak operators alreayd have mass matrix "included"

                ! add regular diffusion in top 3 layers:
                if (nu_top>0 .and. k<=3) then
                   lap_p=laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                   lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                endif
                nu_scale_top = 1
                if (k==1) nu_scale_top=4
                if (k==2) nu_scale_top=2

                do j=1,np
                   do i=1,np
                      if (nu_p==0) then
                         ! normalize so as to conserve IE
                         ! scale by 1/rho (normalized to be O(1))
                         ! dp/dn = O(ps0)*O(delta_eta) = O(ps0)/O(nlev)
                         dpdn = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                              ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,nt)
                         dpdn0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                              ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
                         nu_scale = dpdn0/dpdn
                      endif

                      ! biharmonic terms need a negative sign:
                      if (nu_top>0 .and. k<=3) then
                         utens_tmp=(-nu*vtens(i,j,1,k,ie) + nu_scale_top*nu_top*lap_v(i,j,1))
                         vtens_tmp=(-nu*vtens(i,j,2,k,ie) + nu_scale_top*nu_top*lap_v(i,j,2))
                         ptens_tmp=nu_scale*(-nu_s*ptens(i,j,k,ie) + nu_scale_top*nu_top*lap_p(i,j) )
                      else
                         utens_tmp=-nu*vtens(i,j,1,k,ie)
                         vtens_tmp=-nu*vtens(i,j,2,k,ie)
                         ptens_tmp=-nu_scale*nu_s*ptens(i,j,k,ie)
                      endif

                      ptens(i,j,k,ie) = ptens_tmp
                      vtens(i,j,1,k,ie)=utens_tmp
                      vtens(i,j,2,k,ie)=vtens_tmp
                   enddo
                enddo
             enddo

             pstens(:,:,ie)  =  -nu_p*pstens(:,:,ie)

             kptr=0
             call edgeVpack(edge3, ptens(:,:,:,ie),nlev,kptr,ie)
             kptr=nlev
             call edgeVpack(edge3,vtens(:,:,:,:,ie),2*nlev,kptr,ie)
             kptr=3*nlev
             call edgeVpack(edge3,pstens(:,:,ie),1,kptr,ie)
          enddo


          call bndry_exchangeV(hybrid,edge3)

          do ie=nets,nete

             kptr=0
             call edgeVunpack(edge3, ptens(:,:,:,ie), nlev, kptr, ie)
             kptr=nlev
             call edgeVunpack(edge3, vtens(:,:,:,:,ie), 2*nlev, kptr, ie)


             ! apply inverse mass matrix, accumulate tendencies
#if (defined COLUMN_OPENMP)
             !$omp parallel do private(k)
#endif
             ! dt removed here since its going in the residual
             do k=1,nlev
                vtens(:,:,1,k,ie)=vtens(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)
                vtens(:,:,2,k,ie)=vtens(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
                ptens(:,:,k,ie)=ptens(:,:,k,ie)*elem(ie)%rspheremp(:,:)
             enddo

             ! apply hypervis to u -> u+utens:
             ! E0 = dpdn * .5*u dot u + dpdn * T  + dpdn*PHIS
             ! E1 = dpdn * .5*(u+utens) dot (u+utens) + dpdn * (T-X) + dpdn*PHIS
             ! E1-E0:   dpdn (u dot utens) + dpdn .5 utens dot utens   - dpdn X
             !      X = (u dot utens) + .5 utens dot utens
             !  alt:  (u+utens) dot utens
             do k=1,nlev
                do j=1,np
                   do i=1,np
                      ! update v first (gives better results than updating v after heating)

                      ! subtract because we are updating the residual
                      ! verify: is heating also of opposite sign (added not subtracted since its being accumulated in residual?)
                      fv(i,j,:,k,ie)=fv(i,j,:,k,ie) - vtens(i,j,:,k,ie)
                      !                    fv(i,j,:,k,ie)=fv(i,j,:,k,ie) ! debug

                      v1=elem(ie)%state%v(i,j,1,k,nt)
                      v2=elem(ie)%state%v(i,j,2,k,nt)
                      heating = (vtens(i,j,1,k,ie)*v1  + vtens(i,j,2,k,ie)*v2 )
                      ft(i,j,k,ie)=ft(i,j,k,ie)  - ptens(i,j,k,ie) + heating/cp
                      !                    ft(i,j,k,ie)=ft(i,j,k,ie)  ! debug

                   enddo
                enddo
             enddo

             if (nu_p>0) then
                kptr=3*nlev
                call edgeVunpack(edge3, pstens(:,:,ie), 1, kptr, ie)
                ! take out dt and substract since being accumulated in residual
                pstens(:,:,ie)=pstens(:,:,ie)*elem(ie)%rspheremp(:,:)
                fp(:,:,ie)=fp(:,:,ie) - pstens(:,:,ie)
                !              fp(:,:,ie)=fp(:,:,ie) ! debug
             endif

          enddo
#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
          !$OMP BARRIER
#endif
#endif
       enddo
    endif

    call t_stopf('rhs_advance_hypervis')
    call t_adj_detailf(-1)

  end subroutine rhs_advance_hypervis



  subroutine rhs_advance_hypervis_dp(edge3,fv,ft,fp3,elem,hvcoord,hybrid,deriv,nt,nets,nete,dt2,eta_ave_w)
    !
    !  take one timestep of:
    !          fv(:,:,:,np) = fv(:,:,:,np) +  nu*laplacian**order ( v )
    !          ft(:,:,:,np) = ft(:,:,:,np) +  nu_s*laplacian**order ( T )
    !
    use dimensions_mod, only : np, np, nlev, nc, ntrac, max_corner_elem
    use control_mod, only : nu, nu_div, nu_s, hypervis_order, hypervis_subcycle, nu_p, nu_top, psurf_vis, swest
    use hybrid_mod, only : hybrid_t
    use hybvcoord_mod, only : hvcoord_t
    use element_mod, only : element_t
    use derivative_mod, only : derivative_t, laplace_sphere_wk, vlaplace_sphere_wk
    use derivative_mod, only : subcell_Laplace_fluxes, subcell_dss_fluxes
    use edge_mod, only : edgevpack, edgevunpack, edgeDGVunpack
    use edgetype_mod, only : EdgeBuffer_t, EdgeDescriptor_t
    use bndry_mod, only : bndry_exchangev
    use viscosity_mod, only : biharmonic_wk_dp3d
    use physical_constants, only: Cp
    use derivative_mod, only : subcell_Laplace_fluxes
    !  use time_mod, only : TimeLevel_t
    use prim_advance_mod, only : distribute_flux_at_corners
    implicit none

    type (hybrid_t)      , intent(in) :: hybrid
    type (element_t)     , intent(inout), target :: elem(:)
    type (EdgeBuffer_t)  , intent(inout) :: edge3
    type (derivative_t)  , intent(in) :: deriv
    type (hvcoord_t), intent(in)      :: hvcoord
    !  type (TimeLevel_t)   , intent(in) :: tl
    real (kind=real_kind), intent(inout) :: fv(:,:,:,:,:)
    real (kind=real_kind), intent(inout) :: ft(:,:,:,:)
    real (kind=real_kind), intent(inout) :: fp3(:,:,:,:)

    real (kind=real_kind) :: dt2
    integer :: nets,nete

    ! local
    real (kind=real_kind) :: eta_ave_w  ! weighting for mean flux terms
    real (kind=real_kind) :: dpdn,dpdn0, nu_scale_top
    integer :: k,kptr,i,j,ie,ic,nt
    real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)      :: vtens
    real (kind=real_kind), dimension(np,np,nlev,nets:nete)        :: ttens
    real (kind=real_kind), dimension(np,np,nlev,nets:nete)        :: dptens
    real (kind=real_kind), dimension(0:np+1,0:np+1,nlev)          :: corners
    real (kind=real_kind), dimension(2,2,2)                       :: cflux
    real (kind=real_kind), dimension(nc,nc,4,nlev,nets:nete)      :: dpflux
    real (kind=real_kind), dimension(np,np,nlev) :: p
    real (kind=real_kind), dimension(np,np) :: dptemp1,dptemp2
    type (EdgeDescriptor_t)                                       :: desc


    ! NOTE: PGI compiler bug: when using spheremp, rspheremp and ps as pointers to elem(ie)% members,
    !       data is incorrect (offset by a few numbers actually)
    !       removed for now.
    !       real (kind=real_kind), dimension(:,:), pointer :: spheremp,rspheremp
    !       real (kind=real_kind), dimension(:,:,:), pointer   :: ps

    real (kind=real_kind), dimension(np,np) :: lap_t,lap_dp
    real (kind=real_kind), dimension(np,np,2) :: lap_v
    real (kind=real_kind) :: v1,v2,dt,heating,utens_tmp,vtens_tmp,ttens_tmp,dptens_tmp

    real (kind=real_kind)                     :: temp      (np,np,nlev)
    real (kind=real_kind)                     :: laplace_fluxes(nc,nc,4)



    if (nu_s == 0 .and. nu == 0 .and. nu_p==0 ) return;
    !JMD  call t_barrierf('sync_advance_hypervis', hybrid%par%comm)
    call t_adj_detailf(+1)
    call t_startf('rhs_advance_hypervis_dp')
    dt=dt2
    !took out regluar viscosity option

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  hyper viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! nu_p=0:
    !   scale T dissipaton by dp  (conserve IE, dissipate T^2)
    ! nu_p>0
    !   dont scale:  T equation IE dissipation matches (to truncation error)
    !                IE dissipation from continuity equation
    !                (1 deg: to about 0.1 W/m^2)
    !
    if (hypervis_order == 2) then
       do ic=1,hypervis_subcycle
          call biharmonic_wk_dp3d(elem,dptens,dpflux,ttens,vtens,deriv,edge3,hybrid,nt,nets,nete)

          do ie=nets,nete

             ! comptue mean flux
             if (nu_p>0) then
                elem(ie)%derived%dpdiss_ave(:,:,:)=elem(ie)%derived%dpdiss_ave(:,:,:)+&
                     eta_ave_w*elem(ie)%state%dp3d(:,:,:,nt)/hypervis_subcycle
                elem(ie)%derived%dpdiss_biharmonic(:,:,:)=elem(ie)%derived%dpdiss_biharmonic(:,:,:)+&
                     eta_ave_w*dptens(:,:,:,ie)/hypervis_subcycle
             endif
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,i,j,lap_t,lap_dp,lap_v,nu_scale_top,utens_tmp,vtens_tmp,ttens_tmp,dptens_tmp,laplace_fluxes)
#endif
             do k=1,nlev
                ! advance in time.
                ! note: DSS commutes with time stepping, so we can time advance and then DSS.
                ! note: weak operators alreayd have mass matrix "included"

                ! add regular diffusion in top 3 layers:
                if (nu_top>0 .and. k<=3) then
                   lap_t=laplace_sphere_wk(elem(ie)%state%T(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                   lap_dp=laplace_sphere_wk(elem(ie)%state%dp3d(:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                   lap_v=vlaplace_sphere_wk(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.)
                endif
                nu_scale_top = 1
                if (k==1) nu_scale_top=4
                if (k==2) nu_scale_top=2

                do j=1,np
                   do i=1,np
                      ! biharmonic terms need a negative sign:
                      if (nu_top>0 .and. k<=3) then
                         utens_tmp=(-nu*vtens(i,j,1,k,ie) + nu_scale_top*nu_top*lap_v(i,j,1))
                         vtens_tmp=(-nu*vtens(i,j,2,k,ie) + nu_scale_top*nu_top*lap_v(i,j,2))
                         ttens_tmp=(-nu_s*ttens(i,j,k,ie) + nu_scale_top*nu_top*lap_t(i,j) )
                         dptens_tmp=(-nu_p*dptens(i,j,k,ie) + nu_scale_top*nu_top*lap_dp(i,j) )
                      else
                         utens_tmp=-nu*vtens(i,j,1,k,ie)
                         vtens_tmp=-nu*vtens(i,j,2,k,ie)
                         ttens_tmp=-nu_s*ttens(i,j,k,ie)
                         dptens_tmp=-nu_p*dptens(i,j,k,ie)
                      endif
                      ttens(i,j,k,ie) = ttens_tmp
                      dptens(i,j,k,ie) =dptens_tmp
                      vtens(i,j,1,k,ie)=utens_tmp
                      vtens(i,j,2,k,ie)=vtens_tmp
                   enddo
                enddo
                if (0<ntrac) then
                   elem(ie)%sub_elem_mass_flux(:,:,:,k) = elem(ie)%sub_elem_mass_flux(:,:,:,k) - &
                        eta_ave_w*nu_p*dpflux(:,:,:,k,ie)/hypervis_subcycle
                   if (nu_top>0 .and. k<=3) then
                      laplace_fluxes=subcell_Laplace_fluxes(elem(ie)%state%dp3d(:,:,k,nt),deriv,elem(ie),np,nc)
                      elem(ie)%sub_elem_mass_flux(:,:,:,k) = elem(ie)%sub_elem_mass_flux(:,:,:,k) + &
                           eta_ave_w*nu_scale_top*nu_top*laplace_fluxes/hypervis_subcycle
                   endif
                endif
             enddo


             kptr=0
             call edgeVpack(edge3, ttens(:,:,:,ie),nlev,kptr,ie)
             kptr=nlev
             call edgeVpack(edge3,vtens(:,:,:,:,ie),2*nlev,kptr,ie)
             kptr=3*nlev
             call edgeVpack(edge3,dptens(:,:,:,ie),nlev,kptr,ie)
          enddo


          call bndry_exchangeV(hybrid,edge3)

          do ie=nets,nete

             kptr=0
             call edgeVunpack(edge3, ttens(:,:,:,ie), nlev, kptr, ie)
             kptr=nlev
             call edgeVunpack(edge3, vtens(:,:,:,:,ie), 2*nlev, kptr, ie)
             kptr=3*nlev
             if (0<ntrac) then
                do k=1,nlev
                   temp(:,:,k) = dptens(:,:,k,ie) / elem(ie)%spheremp
                enddo
                corners = 0.0d0
                corners(1:np,1:np,:) = dptens(:,:,:,ie)
             endif
             call edgeVunpack(edge3, dptens(:,:,:,ie), nlev, kptr, ie)



             if (0<ntrac) then
                kptr=3*nlev
                desc = elem(ie)%desc

                call edgeDGVunpack(edge3, corners, nlev, kptr, ie)
                corners = corners/dt

                do k=1,nlev
                   temp(:,:,k) =  elem(ie)%rspheremp(:,:)*elem(ie)%state%dp3d(:,:,k,nt) - temp(:,:,k)
                   temp(:,:,k) =  temp(:,:,k)/dt

                   call distribute_flux_at_corners(cflux, corners(:,:,k), desc%getmapP)

                   cflux(1,1,:)   = elem(ie)%rspheremp(1,  1) * cflux(1,1,:)
                   cflux(2,1,:)   = elem(ie)%rspheremp(np, 1) * cflux(2,1,:)
                   cflux(1,2,:)   = elem(ie)%rspheremp(1, np) * cflux(1,2,:)
                   cflux(2,2,:)   = elem(ie)%rspheremp(np,np) * cflux(2,2,:)

                   elem(ie)%sub_elem_mass_flux(:,:,:,k) = elem(ie)%sub_elem_mass_flux(:,:,:,k) + &
                        eta_ave_w*subcell_dss_fluxes(temp(:,:,k), np, nc, elem(ie)%metdet,cflux)/hypervis_subcycle
                end do
             endif



             ! apply inverse mass matrix, accumulate tendencies
#if (defined COLUMN_OPENMP)
             !$omp parallel do private(k)
#endif
             ! remove dt since its going in the residual
             do k=1,nlev
                vtens(:,:,1,k,ie)=vtens(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)
                vtens(:,:,2,k,ie)=vtens(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
                ttens(:,:,k,ie)=ttens(:,:,k,ie)*elem(ie)%rspheremp(:,:)
                dptens(:,:,k,ie)=dptens(:,:,k,ie)*elem(ie)%rspheremp(:,:)
             enddo

             ! apply hypervis to u -> u+utens:
             ! E0 = dpdn * .5*u dot u + dpdn * T  + dpdn*PHIS
             ! E1 = dpdn * .5*(u+utens) dot (u+utens) + dpdn * (T-X) + dpdn*PHIS
             ! E1-E0:   dpdn (u dot utens) + dpdn .5 utens dot utens   - dpdn X
             !      X = (u dot utens) + .5 utens dot utens
             !  alt:  (u+utens) dot utens
#if (defined COLUMN_OPENMP)
             !$omp parallel do private(k,i,j,v1,v2,heating)
#endif
             do k=1,nlev
                do j=1,np
                   do i=1,np
                      ! update v first (gives better results than updating v after heating)

                      ! subtract not add since going in the residual
                      fv(i,j,:,k,ie)=fv(i,j,:,k,nt) - vtens(i,j,:,k,ie)

                      v1=elem(ie)%state%v(i,j,1,k,nt)
                      v2=elem(ie)%state%v(i,j,2,k,nt)
                      heating = (vtens(i,j,1,k,ie)*v1  + vtens(i,j,2,k,ie)*v2 )
                      ft(i,j,k,ie)=ft(i,j,k,ie)  - ttens(i,j,k,ie) + heating/cp

                      fp3(i,j,k,nt)=fp3(i,j,k,nt) - dptens(i,j,k,ie)

                   enddo
                enddo
             enddo
          enddo
#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
          !$OMP BARRIER
#endif
#endif
       enddo
    endif

    call t_stopf('rhs_advance_hypervis_dp')
    call t_adj_detailf(-1)

  end subroutine rhs_advance_hypervis_dp

end module prim_implicit_mod

#endif ! TRILINOS
