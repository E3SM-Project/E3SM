#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module dg3d_dynamics_mod
!=======================================================================================================!
 use kinds,           only: real_kind
 use dimensions_mod,  only: ne, np, nlev
 use element_mod,     only: element_t   
 use derivative_mod,  only: derivative_t, gradient_wk
 use quadrature_mod,  only: quadrature_t, gausslobatto, jacobi
 use physical_constants, only : rearth , rrearth, g, dd_pi 
 use dg_core_mod, only : sphere2contra, sphere2cov 
 use dg3d_core_mod, only: dg3d_source_term, dg3d_gradient_mass, dg3d_gradient_mom,  &
     dg3d_laxfred_flux, dg3d_fluxjacobian, check_fj, real_vorticity 
 use control_mod, only : nu

!=======================================================================================================!
!=======================================================================================================!
implicit none
private
 public :: dg3d_rhs_terms         !includes source (forcing)  terms 
 public :: dg3d_uvform_rhs
!public :: dg3d_rhs_economy
!=======================================================================================================!
 private:: real_vorticity 
 private:: source_term 
 private:: gradient_mass   
 private:: gradient_mom
 public :: pres_grad_term
 public :: gradient_p3d
 private:: cube_laplacian
 private:: general_grad
!=======================================================================================================!    
 public :: diverge 
 public :: divergence_cov 
 public :: cov_vorticity
 public :: gradient_rhs
 public :: diffusion_uv
 public :: diffusion_temp
 public :: diffusion_hypr
 public :: diffusion_theta
 public :: biharmonic_diff
 public :: horizontal_diff
 public :: implicit_diff
!============================================= LDG Specific ===========================================!
 public :: dg3d_cent_flux
 public :: dg3d_diff_grads
 public :: dg3d_diff_grads_uv
 public :: dg3d_diff_flux
 private :: jump_fluxint
 private :: central_fluxint
 private :: ldg_grads
!=======================================================================================================!
real (kind=real_kind), dimension(:,:,:,:), pointer :: met,metinv, D, Dinv
real (kind=real_kind), dimension(:,:),     pointer :: fcor,rmv,mv,metdet,rmetdet
real (kind=real_kind), dimension(np,np) :: mmx, mmxi   , mmi
real (kind=real_kind), parameter, private:: zero = 0.0D0
real (kind=real_kind), parameter, private:: one  = 1.0D0   
real (kind=real_kind), parameter, private:: dcof = 1.5D5
real (kind=real_kind), parameter, private:: dcof2= 0.2D15
real (kind=real_kind), parameter, private:: tcof = 0.0D0
!=======================================================================================================!
 contains
!============== Horizontal Aspects of Dynamics 	(2d Slices of SW-system)================================!
!=======================================================================================================!
subroutine dg3d_uvform_rhs(elem,klev,neq,deriv,uvbuf,htbuf,dpbuf,ptbuf,qtbuf,  &
                           uv,cori,ht,dp,pt,qt,htop,pgrad,rhsf,sw3d_rhs)
!=======================================================================================================!
    Implicit None
    type(element_t), target, intent(in) :: elem
    integer,intent(in) :: klev, neq
    type (derivative_t),intent(in):: deriv 
    real (kind=real_kind), dimension(0:np+1,0:np+1),   intent(in)   :: htbuf,dpbuf, ptbuf, qtbuf 
    real (kind=real_kind), dimension(0:np+1,0:np+1,2), intent(in)   :: uvbuf
    real (kind=real_kind), dimension(np,np,2),         intent(in)   :: uv, pgrad
    real (kind=real_kind), dimension(np,np,neq),       intent(in)   :: rhsf
    real (kind=real_kind), dimension(np,np),           intent(in)   :: ht,pt,qt, dp, cori, htop
    real (kind=real_kind), dimension(np,np,neq),       intent(out)  :: sw3d_rhs
!=======================================================================================================!
    integer, parameter  :: south=1, east=2, north=3, west=4, neqn=5 !including Q 
    real (kind=real_kind), dimension(4)       :: fjmax
    real (kind=real_kind), dimension(4,np)    :: SL,SR,lambda                           
    real (kind=real_kind), dimension(np,np)   :: phi,gcori,dpg,ptg,qtg,  tdif , udif,vdif
    real (kind=real_kind), dimension(np,np)   :: energy,gh11,gh22, cdiv, ulap,vlap
    real (kind=real_kind), dimension(np,np)   :: sg, u1_rhs, u2_rhs, vor ,vorg 
    real (kind=real_kind), dimension(np,np,2) :: contrauv,couv, ghij,uvflx, sp_grad
    real (kind=real_kind), dimension(np,4,2)  :: uvflx_senw, uv_senw, contuv_senw, couv_senw 
    real (kind=real_kind), dimension(np,4)    :: vflx_senw, uflx_senw
    real (kind=real_kind), dimension(np,4)    :: gh11_senw, gh22_senw, ht_senw, energy_senw
    real (kind=real_kind), dimension(np,np,2) :: guv
    real (kind=real_kind), dimension(np,4)    :: dp_senw, pt_senw, dpg_senw, ptg_senw, phi_senw 
    real (kind=real_kind), dimension(np,4)    :: qt_senw, qtg_senw 
!=======================================================================================================!    
    real (kind=real_kind), dimension(np,np,  neqn)  :: sw_grad, flux_sw, sw_source, rhs_force
    real (kind=real_kind), dimension(np,np,  neqn)  :: sw_vec 
    real (kind=real_kind), dimension(np,4,   neqn)  :: sw_vec_halo
    real (kind=real_kind), dimension(np,np,2,neqn)  :: swsys_fluxvec 
    real (kind=real_kind), dimension(np,4,2, neqn)  :: swsys_flux_halo 
    real (kind=real_kind):: rtmp,alfa1,alfa2,u1,u2,fact,grv,damp, term 
    Integer:: i,j,k, wall, eqn 
!=======================================================================================================!
!	From HOMME
!=======================================================================================================!
 met    => elem%met
 metinv => elem%metinv
 metdet => elem%metdet
 mv     => elem%mp   
 fcor   => elem%fcor
 Dinv   => elem%Dinv
 D      => elem%D
 grv = g 
!=======================================================================================================!
!  Mass and Inverse-Mass Matrix										!
!=======================================================================================================!
!=======================================================================================================!
#if 1
     mmx(:,:) = mv(:,:)
#if defined(_USE_VECTOR)
     call vrec(mmxi(1,1),mmx(1,1),nv*nv)
#else
     do j = 1,np
     do i = 1,np
        mmxi(i,j) = one / mmx(i,j) 
     enddo
     enddo
#endif
#endif
!=======================================================================================================!
 do j= 1,np
 do i= 1,np
  mmx(i,j) =  elem%mp(i,j)      !mass matrix 
  mmi(i,j) = 1.0D0 / mmx(i,j)
  sg(i,j)  = elem%metdet(i,j)   !metric {sqrt(g)} 
 enddo
 enddo
!=======================================================================================================!
! Flux computations for the equations
!=======================================================================================================!
 gcori= zero
 dpg  = zero
 ptg  = zero
 phi  = zero
 energy= zero

 do j = 1,np
 do i = 1,np
    dpg(i,j)  = dp(i,j) * sg(i,j)    !Delta_P
    ptg(i,j)  = pt(i,j) * dpg(i,j)   !Delta_Th
    qtg(i,j)  = qt(i,j) * dpg(i,j)   !Delta_Q 
 enddo
 enddo
!=======================================================================================================!
!	Boundary values for velocity & flux terms from the neighbors
!=======================================================================================================!
 do j = 1, 2 
 do i = 1,np
    uv_senw(i,south,j) = uvbuf(i,   0,j)
    uv_senw(i,east,j)  = uvbuf(np+1,i,j)
    uv_senw(i,north,j) = uvbuf(i,np+1,j)
    uv_senw(i,west,j)  = uvbuf(0,i,   j)
 enddo
 enddo
!=======================================================================================================!
!  Contravariant components for the halo region (needed for continuity eqns) 
!=======================================================================================================!
    do k= 1, np
       contuv_senw(k,1,1) = Dinv(1,1,k,1)  * uv_senw(k,1,1) + Dinv(1,2,k,1)  * uv_senw(k,1,2)
       contuv_senw(k,1,2) = Dinv(2,1,k,1)  * uv_senw(k,1,1) + Dinv(2,2,k,1)  * uv_senw(k,1,2)
       contuv_senw(k,2,1) = Dinv(1,1,np,k) * uv_senw(k,2,1) + Dinv(1,2,np,k) * uv_senw(k,2,2)
       contuv_senw(k,2,2) = Dinv(2,1,np,k) * uv_senw(k,2,1) + Dinv(2,2,np,k) * uv_senw(k,2,2)

       contuv_senw(k,3,1) = Dinv(1,1,k,np) * uv_senw(k,3,1) + Dinv(1,2,k,np) * uv_senw(k,3,2)
       contuv_senw(k,3,2) = Dinv(2,1,k,np) * uv_senw(k,3,1) + Dinv(2,2,k,np) * uv_senw(k,3,2)
       contuv_senw(k,4,1) = Dinv(1,1,1,k)  * uv_senw(k,4,1) + Dinv(1,2,1,k)  * uv_senw(k,4,2)
       contuv_senw(k,4,2) = Dinv(2,1,1,k)  * uv_senw(k,4,1) + Dinv(2,2,1,k)  * uv_senw(k,4,2)
    enddo

!=======================================================================================================!
! Contravariant & Covariant vectors from known Sph (u,v) 
!=======================================================================================================!

    contrauv(:,:,:) = sphere2contra(uv,Dinv)
        couv(:,:,:) = sphere2cov(uv,D)

!=======================================================================================================!
!  Covariant components for the halo region
!=======================================================================================================!
    do k= 1, np
       couv_senw(k,1,1) = met(1,1,k,1)  * contuv_senw(k,1,1) + met(1,2,k,1)  * contuv_senw(k,1,2)
       couv_senw(k,1,2) = met(2,1,k,1)  * contuv_senw(k,1,1) + met(2,2,k,1)  * contuv_senw(k,1,2)
       couv_senw(k,2,1) = met(1,1,np,k) * contuv_senw(k,2,1) + met(1,2,np,k) * contuv_senw(k,2,2)
       couv_senw(k,2,2) = met(2,1,np,k) * contuv_senw(k,2,1) + met(2,2,np,k) * contuv_senw(k,2,2)

       couv_senw(k,3,1) = met(1,1,k,np) * contuv_senw(k,3,1) + met(1,2,k,np) * contuv_senw(k,3,2)
       couv_senw(k,3,2) = met(2,1,k,np) * contuv_senw(k,3,1) + met(2,2,k,np) * contuv_senw(k,3,2)
       couv_senw(k,4,1) = met(1,1,1,k)  * contuv_senw(k,4,1) + met(1,2,1,k)  * contuv_senw(k,4,2)
       couv_senw(k,4,2) = met(2,1,1,k)  * contuv_senw(k,4,1) + met(2,2,1,k)  * contuv_senw(k,4,2)
    enddo
!=======================================================================================================!
 do i = 1,np
    ht_senw(i,south) = htbuf(i,   0)
    ht_senw(i,east)  = htbuf(np+1,i)
    ht_senw(i,north) = htbuf(i,np+1)
    ht_senw(i,west)  = htbuf(0,i   )
 
    dp_senw(i,south) = dpbuf(i,   0) 
    dp_senw(i,east)  = dpbuf(np+1,i) 
    dp_senw(i,north) = dpbuf(i,np+1) 
    dp_senw(i,west)  = dpbuf(0,i   ) 

    pt_senw(i,south) = ptbuf(i,   0)
    pt_senw(i,east)  = ptbuf(np+1,i)
    pt_senw(i,north) = ptbuf(i,np+1)
    pt_senw(i,west)  = ptbuf(0,i   )

    qt_senw(i,south) = qtbuf(i,   0)
    qt_senw(i,east)  = qtbuf(np+1,i)
    qt_senw(i,north) = qtbuf(i,np+1)
    qt_senw(i,west)  = qtbuf(0,i   )
 end do
!=======================================================================================================!
 do k = 1,np
    dpg_senw(k,south) = dp_senw(k,south) *  sg(k,1)
    dpg_senw(k,east)  = dp_senw(k,east)  *  sg(np,k)
    dpg_senw(k,north) = dp_senw(k,north) *  sg(k,np)
    dpg_senw(k,west)  = dp_senw(k,west)  *  sg(1,k)

    ptg_senw(k,south) = pt_senw(k,south) * dpg_senw(k,south)
    ptg_senw(k,east)  = pt_senw(k,east)  * dpg_senw(k,east)
    ptg_senw(k,north) = pt_senw(k,north) * dpg_senw(k,north)
    ptg_senw(k,west)  = pt_senw(k,west)  * dpg_senw(k,west) 

    qtg_senw(k,south) = qt_senw(k,south) * dpg_senw(k,south)
    qtg_senw(k,east)  = qt_senw(k,east)  * dpg_senw(k,east)
    qtg_senw(k,north) = qt_senw(k,north) * dpg_senw(k,north)
    qtg_senw(k,west)  = qt_senw(k,west)  * dpg_senw(k,west) 

 end do
!=======================================================================================================!
!       Energy flux   (KE + PE) 
!=======================================================================================================!

 do j = 1,np
 do i = 1,np
    energy(i,j)= grv*ht(i,j) + 0.5D0* (uv(i,j,1)*uv(i,j,1) + uv(i,j,2)*uv(i,j,2) )
 end do
 end do
!=======================================================================================================!
!	Energy for the Halo region from the neighbours
!=======================================================================================================!
 do wall = 1, 4
 do k = 1,np
    energy_senw(k,wall) = grv* ht_senw(k,wall) + 0.5D0*(uv_senw(k,wall,1)*uv_senw(k,wall,1)   &
                                                    + uv_senw(k,wall,2)*uv_senw(k,wall,2) )
 end do
 end do  
!=======================================================================================================!
!	Flux Jacobian of the SW system {fj --> Max [|u| + sqrt(gh*G^ii)] }
!=======================================================================================================!
 do j = 1,np
 do i = 1,np
    ghij(i,j,1) = (grv * htop(i,j) )* metinv(1,1,i,j)
    ghij(i,j,2) = (grv * htop(i,j) )* metinv(2,2,i,j)
 end do
 end do

 do k = 1,np
    gh11_senw(k,1) = (grv * abs(ht_senw(k,1)) )* metinv(1,1,k,1)
    gh11_senw(k,2) = (grv * abs(ht_senw(k,2)) )* metinv(1,1,np,k)
    gh11_senw(k,3) = (grv * abs(ht_senw(k,3)) )* metinv(1,1,k,np)
    gh11_senw(k,4) = (grv * abs(ht_senw(k,4)) )* metinv(1,1,1,k)

    gh22_senw(k,1) = (grv * abs(ht_senw(k,1)) )* metinv(2,2,k,1)
    gh22_senw(k,2) = (grv * abs(ht_senw(k,2)) )* metinv(2,2,np,k)
    gh22_senw(k,3) = (grv * abs(ht_senw(k,3)) )* metinv(2,2,k,np)
    gh22_senw(k,4) = (grv * abs(ht_senw(k,4)) )* metinv(2,2,1,k)
 end do
!=======================================================================================================!
!  fjmax(:) = dg3d_fluxjacobian(uv,uv_senw,ghij(:,:,1),ghij(:,:,2),gh11_senw,gh22_senw)
   fjmax(:) = check_fj(contrauv,contuv_senw,ghij(:,:,1),ghij(:,:,2),gh11_senw,gh22_senw)
!=======================================================================================================!
!=======================================================================================================!
 sw_vec        = 0.0D0
 swsys_fluxvec = 0.0D0
 sw_vec_halo     = 0.0D0
 swsys_flux_halo = 0.0D0

 do j = 1,np
 do i = 1,np
   ! State vectors for the 3D system 
    sw_vec(i,j,1) = couv(i,j,1)
    sw_vec(i,j,2) = couv(i,j,2)
    sw_vec(i,j,3) = dpg(i,j)
    sw_vec(i,j,4) = ptg(i,j)
    sw_vec(i,j,5) = qtg(i,j)

   ! Flux vectors for the 3D system 
    swsys_fluxvec(i,j,1,1) = energy(i,j)
    swsys_fluxvec(i,j,2,1) = 0.0D0      
    swsys_fluxvec(i,j,1,2) = 0.0D0      
    swsys_fluxvec(i,j,2,2) = energy(i,j) 

    swsys_fluxvec(i,j,1,3) = contrauv(i,j,1) * dpg(i,j)
    swsys_fluxvec(i,j,2,3) = contrauv(i,j,2) * dpg(i,j)

    swsys_fluxvec(i,j,1,4) = contrauv(i,j,1) * ptg(i,j)
    swsys_fluxvec(i,j,2,4) = contrauv(i,j,2) * ptg(i,j)

    swsys_fluxvec(i,j,1,5) = contrauv(i,j,1) * qtg(i,j)
    swsys_fluxvec(i,j,2,5) = contrauv(i,j,2) * qtg(i,j)

 enddo
 enddo

!=======================================================================================================!
! Corresponding "halo" form the neighbouring elements 
!=======================================================================================================!

 do wall = 1, 4 
 do k  = 1,np 
    sw_vec_halo(k,wall,1) = couv_senw(k,wall,1)
    sw_vec_halo(k,wall,2) = couv_senw(k,wall,2)
    sw_vec_halo(k,wall,3) = dpg_senw(k,wall) 
    sw_vec_halo(k,wall,4) = ptg_senw(k,wall) 
    sw_vec_halo(k,wall,5) = qtg_senw(k,wall) 

    swsys_flux_halo(k,wall,1,1) = energy_senw(k,wall)
    swsys_flux_halo(k,wall,2,1) = 0.0D0          
    swsys_flux_halo(k,wall,1,2) = 0.0D0          
    swsys_flux_halo(k,wall,2,2) = energy_senw(k,wall)
 enddo 
 enddo

 do j = 1, 2
 do wall= 1, 4 
 do k= 1,np
    swsys_flux_halo(k,wall,j,3) = contuv_senw(k,wall,j) * dpg_senw(k,wall)
    swsys_flux_halo(k,wall,j,4) = contuv_senw(k,wall,j) * ptg_senw(k,wall)
    swsys_flux_halo(k,wall,j,5) = contuv_senw(k,wall,j) * qtg_senw(k,wall)
 end do
 end do
 end do


!=======================================================================================================!
! Divergence Damping
! do j=1,np
! do i=1,np
!     guv(i,j,1)= couv(i,j,1) * elem%metdet(i,j)
!     guv(i,j,2)= couv(i,j,2) * elem%metdet(i,j)
!     udif(i,j) = couv(i,j,1)
!     vdif(i,j) = couv(i,j,2)
! end do
! end do
!=======================================================================================================!
! Optional diffusion or forcings
!=======================================================================================================!
 cdiv(:,:) = zero
    damp = dcof2

  do j = 1,np
  do i = 1,np
       rhs_force(i,j,1) = rhsf(i,j,1)  !+ damp*udif(i,j)                  !u forcing (mometum) 
       rhs_force(i,j,2) = rhsf(i,j,2)  !+ damp*vdif(i,j)                  !v forcing    '' 
       rhs_force(i,j,3) = rhsf(i,j,3)                                    !"No" forcing (mass) 
       rhs_force(i,j,4) = rhsf(i,j,4)
       rhs_force(i,j,5) = rhsf(i,j,5)
  end do
  end do

!rhs_force(:,:,:) = 0.0D0 

!=======================================================================================================!
 flux_sw(:,:,:)= 0.0D0
 sw_grad(:,:,:)= 0.0D0
!=======================================================================================================!  
!    Lax-Friedriechs flux for the boundary points 
!
  call dg3d_laxfred_flux(neqn,deriv,fjmax,sw_vec,sw_vec_halo,swsys_fluxvec,swsys_flux_halo,flux_sw)
!=======================================================================================================!
!    Gradient  for the momentum equations
!
 call dg3d_gradient_mom(deriv,energy,sw_grad(:,:,1),sw_grad(:,:,2))
!=======================================================================================================!    
!    Gradient  for the continuity equation     (pressure)
!
        sw_grad(:,:,3)= dg3d_gradient_mass(deriv,swsys_fluxvec(:,:,:,3))
! 
!    Gradient  for the continuity equation     (Pot-temp)
!
        sw_grad(:,:,4)= dg3d_gradient_mass(deriv,swsys_fluxvec(:,:,:,4))
!
!    Gradient  for the continuity equation     (Q-moisture)
!
        sw_grad(:,:,5)= dg3d_gradient_mass(deriv,swsys_fluxvec(:,:,:,5))
!=======================================================================================================!

 sw3d_rhs(:,:,:)= 0.0D0
 
!sw_source(:,:,:)= 0.0D0

!=======================================================================================================!
! Source term Computation  (for the mometum equations) 
!=======================================================================================================!
!  vor(:,:) = real_vorticity(uv,metdet,D,deriv)
   vorg(:,:) = cov_vorticity(deriv,couv)

! converting into Sph-(u,v) grads { (u,v) = A^-T *(u_1,u_2) }

   do j=1,np
   do i=1,np
      sp_grad(i,j,1) = pgrad(i,j,1) *Dinv(1,1,i,j) + pgrad(i,j,2) * Dinv(2,1,i,j)
      sp_grad(i,j,2) = pgrad(i,j,1) *Dinv(1,2,i,j) + pgrad(i,j,2) * Dinv(2,2,i,j)
   enddo
   enddo

! (u,v)-form RHS momentum forcing

   do j = 1, np
       do i = 1, np
          term =  (vorg(i,j)/sg(i,j) +  cori(i,j))
          sw_source(i,j,1) = ( uv(i,j,2) * term - sp_grad(i,j,1) + rhs_force(i,j,1)) 
          sw_source(i,j,2) = (-uv(i,j,1) * term - sp_grad(i,j,2) + rhs_force(i,j,2)) 
          sw_source(i,j,3) =  0.0D0
          sw_source(i,j,4) =  rhs_force(i,j,4)     !th-forcing 
          sw_source(i,j,5) =  0.0D0
       enddo
   enddo

!=======================================================================================================!
!	Compute RHS of the ODE system corresponding  to the 3D  system        
!=======================================================================================================!

!   do j = 1, np
!   do i = 1, np
!       u1_rhs(i,j)= (sw_grad(i,j,1) - flux_sw(i,j,1)) * mmi(i,j) 
!       u2_rhs(i,j)= (sw_grad(i,j,2) - flux_sw(i,j,2)) * mmi(i,j) 
!    enddo
!    enddo
!
!   do j=1,np
!   do i=1,np
!     sw3d_rhs(i,j,1) = u1_rhs(i,j) *Dinv(1,1,i,j) + u2_rhs(i,j) * Dinv(2,1,i,j) + sw_source(i,j,1)
!     sw3d_rhs(i,j,2) = u1_rhs(i,j) *Dinv(1,2,i,j) + u2_rhs(i,j) * Dinv(2,2,i,j) + sw_source(i,j,2)
!   enddo
!   enddo

! (u,v)-form RHS momentum 

    do j = 1, np
    do i = 1, np
        u1 = (sw_grad(i,j,1) - flux_sw(i,j,1)) 
        u2 = (sw_grad(i,j,2) - flux_sw(i,j,2)) 
      sw3d_rhs(i,j,1) = (u1 *Dinv(1,1,i,j) + u2 * Dinv(2,1,i,j))*mmi(i,j) + sw_source(i,j,1)
      sw3d_rhs(i,j,2) = (u1 *Dinv(1,2,i,j) + u2 * Dinv(2,2,i,j))*mmi(i,j) + sw_source(i,j,2)
     enddo
     enddo
 
  do eqn=3,5
   do j= 1,np
   do i= 1,np       
    sw3d_rhs(i,j,eqn) = (sw_source(i,j,eqn) + sw_grad(i,j,eqn) - flux_sw(i,j,eqn))*mmi(i,j) 
   enddo
   enddo
  enddo

end subroutine  dg3d_uvform_rhs
!=======================================================================================================!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!          contravarint (u1,u2) - formulation 
!=======================================================================================================!
subroutine dg3d_rhs_terms(elem,klev,neq,deriv,uvbuf,htbuf,dpbuf,ptbuf,qtbuf,  &
                          uvcomp,couv,cori,ht,dp,pt,qt,htop,pgrad,rhsf,sw_rhs)
!=======================================================================================================!
    Implicit None
    type(element_t), target, intent(inout) :: elem
    integer,intent(in) :: klev, neq
    type (derivative_t),intent(in):: deriv 
    real (kind=real_kind), dimension(0:np+1,0:np+1),   intent(in)   :: htbuf,dpbuf, ptbuf, qtbuf 
    real (kind=real_kind), dimension(0:np+1,0:np+1,2), intent(in)   :: uvbuf
    real (kind=real_kind), dimension(np,np,2),         intent(in)   :: uvcomp, couv, pgrad
    real (kind=real_kind), dimension(np,np,neq),       intent(in)   :: rhsf
    real (kind=real_kind), dimension(np,np),           intent(in)   :: ht,pt,qt, dp, cori, htop
    real (kind=real_kind), dimension(np,np,neq),       intent(out)  :: sw_rhs
!=======================================================================================================!
    integer, parameter :: south=1, east=2, north=3, west=4, neqn=5 !including Q 
    real (kind=real_kind), dimension(4)       :: fjmax
    real (kind=real_kind), dimension(4,np)    :: SL,SR,lambda                           
    real (kind=real_kind), dimension(np,np)   :: phi,gcori,dpg,ptg,qtg,  tdif , udif,vdif
    real (kind=real_kind), dimension(np,np)   :: energy,gh11,gh22, cdiv, ulap,vlap
    real (kind=real_kind), dimension(np,np,2) :: ghij,uvflx
    real (kind=real_kind), dimension(np,4,2)  :: uvflx_senw, uv_senw, couv_senw
    real (kind=real_kind), dimension(np,4)    :: vflx_senw, uflx_senw
    real (kind=real_kind), dimension(np,4)    :: gh11_senw, gh22_senw, ht_senw, energy_senw
    real (kind=real_kind), dimension(np,np,2) :: guv
    real (kind=real_kind), dimension(np,4)    :: dp_senw, pt_senw, dpg_senw, ptg_senw, phi_senw 
    real (kind=real_kind), dimension(np,4)    :: qt_senw, qtg_senw 
!=======================================================================================================!    
    real (kind=real_kind), dimension(np,np,  neqn)  :: sw_grad, flux_sw, sw_source, rhs_force
    real (kind=real_kind), dimension(np,np,  neqn)  :: sw_vec 
    real (kind=real_kind), dimension(np,4,   neqn)  :: sw_vec_halo
    real (kind=real_kind), dimension(np,np,2,neqn)  :: swsys_fluxvec 
    real (kind=real_kind), dimension(np,4,2, neqn)  :: swsys_flux_halo 
    real (kind=real_kind):: rtmp,alfa1,alfa2,ul,ur,fact,grv,damp
    Integer:: i,j,k, wall, eqn 
!=======================================================================================================!
!	From HOMME
!=======================================================================================================!
 met    => elem%met
 metinv => elem%metinv
 metdet => elem%metdet
 mv     => elem%mp   
 fcor   => elem%fcor
!=======================================================================================================!
 grv = g 
!=======================================================================================================!
!  Mass and Inverse-Mass Matrix										!
!=======================================================================================================!
#if 1
    mmx(:,:) = mv(:,:)
#if defined(_USE_VECTOR)
    call vrec(mmxi(1,1),mmx(1,1),nv*nv)
#else
    do j = 1,np
    do i = 1,np
       mmxi(i,j) = one / mmx(i,j) 
    enddo
    enddo
#endif
#endif
!=======================================================================================================!
!=======================================================================================================!
! Flux computations for the equations
!=======================================================================================================!
 gcori= zero
 dpg  = zero
 ptg  = zero
 phi  = zero
 energy= zero

 do j = 1,np
 do i = 1,np
    dpg(i,j)  = dp(i,j) * elem%metdet(i,j) 
    ptg(i,j)  = pt(i,j) * dpg(i,j)
    qtg(i,j)  = qt(i,j) * dpg(i,j)
    gcori(i,j)= cori(i,j) * elem%metdet(i,j) 
 enddo
 enddo
 do j = 1,np
 do i = 1,np
    energy(i,j)= grv*ht(i,j) + 0.5D0* ( couv(i,j,1)*uvcomp(i,j,1) + couv(i,j,2)*uvcomp(i,j,2) )
 end do
 end do
!=======================================================================================================!
! Divergence Damping
  do j=1,np
  do i=1,np
      guv(i,j,1)= couv(i,j,1) * elem%metdet(i,j)
      guv(i,j,2)= couv(i,j,2) * elem%metdet(i,j)
      udif(i,j) = couv(i,j,1)
      vdif(i,j) = couv(i,j,2)
  end do
  end do

!=======================================================================================================!
! Optional diffusion or forcings
!=======================================================================================================!
 cdiv(:,:) = zero 
!cdiv(:,:) = divergence_cov(deriv,elem%metdet,guv(:,:,:)) 
! udif(:,:) = cube_laplacian(deriv,udif(:,:))
! vdif(:,:) = cube_laplacian(deriv,vdif(:,:))

! if (klev < 4 )  then
!    damp = 5.0D05
  !  udif(:,:) = cube_laplacian(deriv,udif(:,:))
  !  vdif(:,:) = cube_laplacian(deriv,vdif(:,:))
! else
    damp = dcof2
   ! udif(:,:) = cube_laplacian(deriv,udif(:,:))
   ! vdif(:,:) = cube_laplacian(deriv,vdif(:,:))
! endif 
     
 do j = 1,np
 do i = 1,np
      rhs_force(i,j,1) = rhsf(i,j,1)  !+ damp*udif(i,j)                  !u forcing (mometum) 
      rhs_force(i,j,2) = rhsf(i,j,2)  !+ damp*vdif(i,j)                  !v forcing    '' 
      rhs_force(i,j,3) = rhsf(i,j,3)                                    !"No" forcing (mass) 
      rhs_force(i,j,4) = rhsf(i,j,4)   
  !   energy(i,j)      = energy(i,j)  - dcof*cdiv(i,j)                  ! divergence damping 
 end do
 end do

 sw_source(:,:,:)= zero 
 call dg3d_source_term(mv,rmv,deriv,gcori,uvcomp,couv,pgrad,rhs_force,sw_source)
!=======================================================================================================!
!	Boundary values for velocity & flux terms from the neighbors
!=======================================================================================================!
 do j = 1, 2 
 do i = 1,np
    uv_senw(i,south,j) = uvbuf(i,   0,j)
    uv_senw(i,east,j)  = uvbuf(np+1,i,j)
    uv_senw(i,north,j) = uvbuf(i,np+1,j)
    uv_senw(i,west,j)  = uvbuf(0,i,   j)
 enddo
 enddo
!=======================================================================================================!
 do i = 1,np
    ht_senw(i,south) = htbuf(i,   0)
    ht_senw(i,east)  = htbuf(np+1,i)
    ht_senw(i,north) = htbuf(i,np+1)
    ht_senw(i,west)  = htbuf(0,i   )
 
    dp_senw(i,south) = dpbuf(i,   0) 
    dp_senw(i,east)  = dpbuf(np+1,i) 
    dp_senw(i,north) = dpbuf(i,np+1) 
    dp_senw(i,west)  = dpbuf(0,i   ) 

    pt_senw(i,south) = ptbuf(i,   0)
    pt_senw(i,east)  = ptbuf(np+1,i)
    pt_senw(i,north) = ptbuf(i,np+1)
    pt_senw(i,west)  = ptbuf(0,i   )

    qt_senw(i,south) = qtbuf(i,   0)
    qt_senw(i,east)  = qtbuf(np+1,i)
    qt_senw(i,north) = qtbuf(i,np+1)
    qt_senw(i,west)  = qtbuf(0,i   )
 end do
!=======================================================================================================!
 do k = 1,np
    dpg_senw(k,south) = dp_senw(k,south) *  elem%metdet(k,1)
    dpg_senw(k,east)  = dp_senw(k,east)  *  elem%metdet(np,k)
    dpg_senw(k,north) = dp_senw(k,north) *  elem%metdet(k,np)
    dpg_senw(k,west)  = dp_senw(k,west)  *  elem%metdet(1,k)

    ptg_senw(k,south) = pt_senw(k,south) * dpg_senw(k,south)
    ptg_senw(k,east)  = pt_senw(k,east)  * dpg_senw(k,east)
    ptg_senw(k,north) = pt_senw(k,north) * dpg_senw(k,north)
    ptg_senw(k,west)  = pt_senw(k,west)  * dpg_senw(k,west) 

    qtg_senw(k,south) = qt_senw(k,south) * dpg_senw(k,south)
    qtg_senw(k,east)  = qt_senw(k,east)  * dpg_senw(k,east)
    qtg_senw(k,north) = qt_senw(k,north) * dpg_senw(k,north)
    qtg_senw(k,west)  = qt_senw(k,west)  * dpg_senw(k,west) 

 end do
!=======================================================================================================!
!	Covariant components for the halo region 
!=======================================================================================================!
 do k = 1,np
    couv_senw(k,1,1) = elem%met(1,1,k,1)  * uv_senw(k,1,1) + elem%met(1,2,k,1)  * uv_senw(k,1,2)
    couv_senw(k,1,2) = elem%met(2,1,k,1)  * uv_senw(k,1,1) + elem%met(2,2,k,1)  * uv_senw(k,1,2)
    couv_senw(k,2,1) = elem%met(1,1,np,k) * uv_senw(k,2,1) + elem%met(1,2,np,k) * uv_senw(k,2,2)
    couv_senw(k,2,2) = elem%met(2,1,np,k) * uv_senw(k,2,1) + elem%met(2,2,np,k) * uv_senw(k,2,2)

    couv_senw(k,3,1) = elem%met(1,1,k,np) * uv_senw(k,3,1) + elem%met(1,2,k,np) * uv_senw(k,3,2)
    couv_senw(k,3,2) = elem%met(2,1,k,np) * uv_senw(k,3,1) + elem%met(2,2,k,np) * uv_senw(k,3,2)
    couv_senw(k,4,1) = elem%met(1,1,1,k)  * uv_senw(k,4,1) + elem%met(1,2,1,k)  * uv_senw(k,4,2)
    couv_senw(k,4,2) = elem%met(2,1,1,k)  * uv_senw(k,4,1) + elem%met(2,2,1,k)  * uv_senw(k,4,2)
 enddo     
!=======================================================================================================!    
!	Energy for the Halo region from the neighbours
!=======================================================================================================!
 do wall = 1, 4
 do k = 1,np
    energy_senw(k,wall) = grv * ht_senw(k,wall) + 							&
                   	  0.5D0 * (uv_senw(k,wall,1) * couv_senw(k,wall,1) + 				&
                          	   uv_senw(k,wall,2) * couv_senw(k,wall,2) )
 end do
 end do  
!=======================================================================================================!
!	Flux Jacobian of the SW system
!=======================================================================================================!
 do j = 1,np
 do i = 1,np
   !ghij(i,j,1) = (grv * abs(ht(i,j)) )* elem%metinv(1,1,i,j)
   !ghij(i,j,2) = (grv * abs(ht(i,j)) )* elem%metinv(2,2,i,j)
    ghij(i,j,1) = (grv * htop(i,j) )* elem%metinv(1,1,i,j)
    ghij(i,j,2) = (grv * htop(i,j) )* elem%metinv(2,2,i,j)
 end do
 end do

 do k = 1,np
    gh11_senw(k,1) = (grv * abs(ht_senw(k,1)) )* elem%metinv(1,1,k,1)
    gh11_senw(k,2) = (grv * abs(ht_senw(k,2)) )* elem%metinv(1,1,np,k)
    gh11_senw(k,3) = (grv * abs(ht_senw(k,3)) )* elem%metinv(1,1,k,np)
    gh11_senw(k,4) = (grv * abs(ht_senw(k,4)) )* elem%metinv(1,1,1,k)

    gh22_senw(k,1) = (grv * abs(ht_senw(k,1)) )* elem%metinv(2,2,k,1)
    gh22_senw(k,2) = (grv * abs(ht_senw(k,2)) )* elem%metinv(2,2,np,k)
    gh22_senw(k,3) = (grv * abs(ht_senw(k,3)) )* elem%metinv(2,2,k,np)
    gh22_senw(k,4) = (grv * abs(ht_senw(k,4)) )* elem%metinv(2,2,1,k)
 end do
!=======================================================================================================!
!  fjmax(:) = dg3d_fluxjacobian(uvcomp,uv_senw,ghij(:,:,1),ghij(:,:,2),gh11_senw,gh22_senw)
   fjmax(:) = check_fj(uvcomp,uv_senw,ghij(:,:,1),ghij(:,:,2),gh11_senw,gh22_senw)
!=======================================================================================================!
!=======================================================================================================!
 sw_vec       = zero
 swsys_fluxvec= zero
 sw_vec_halo    = zero
 swsys_flux_halo= zero
 do j = 1,np
 do i = 1,np
    sw_vec(i,j,1) = couv(i,j,1)
    sw_vec(i,j,2) = couv(i,j,2)
    sw_vec(i,j,3) = dpg(i,j)
    sw_vec(i,j,4) = ptg(i,j)
    sw_vec(i,j,5) = qtg(i,j)

    swsys_fluxvec(i,j,1,1) = energy(i,j)
    swsys_fluxvec(i,j,2,1) = zero      
    swsys_fluxvec(i,j,1,2) = zero      
    swsys_fluxvec(i,j,2,2) = energy(i,j) 

    swsys_fluxvec(i,j,1,3) = uvcomp(i,j,1) * dpg(i,j)
    swsys_fluxvec(i,j,2,3) = uvcomp(i,j,2) * dpg(i,j)

    swsys_fluxvec(i,j,1,4) = uvcomp(i,j,1) * ptg(i,j)
    swsys_fluxvec(i,j,2,4) = uvcomp(i,j,2) * ptg(i,j)

    swsys_fluxvec(i,j,1,5) = uvcomp(i,j,1) * qtg(i,j)
    swsys_fluxvec(i,j,2,5) = uvcomp(i,j,2) * qtg(i,j)

 enddo
 enddo

 do wall = 1, 4 
 do k  = 1,np 
    sw_vec_halo(k,wall,1) = couv_senw(k,wall,1)
    sw_vec_halo(k,wall,2) = couv_senw(k,wall,2)
    sw_vec_halo(k,wall,3) = dpg_senw(k,wall) 
    sw_vec_halo(k,wall,4) = ptg_senw(k,wall) 
    sw_vec_halo(k,wall,5) = qtg_senw(k,wall) 

    swsys_flux_halo(k,wall,1,1) = energy_senw(k,wall)
    swsys_flux_halo(k,wall,2,1) = zero          
    swsys_flux_halo(k,wall,1,2) = zero          
    swsys_flux_halo(k,wall,2,2) = energy_senw(k,wall)
 enddo 
 enddo

 do j = 1, 2
 do wall= 1, 4 
 do k= 1,np
    swsys_flux_halo(k,wall,j,3) = uv_senw(k,wall,j) * dpg_senw(k,wall)
    swsys_flux_halo(k,wall,j,4) = uv_senw(k,wall,j) * ptg_senw(k,wall)
    swsys_flux_halo(k,wall,j,5) = uv_senw(k,wall,j) * qtg_senw(k,wall)
 end do
 end do
 end do
!=======================================================================================================!
 flux_sw(:,:,:)= zero
!=======================================================================================================!
  call dg3d_laxfred_flux(neqn,deriv,fjmax,sw_vec,sw_vec_halo,swsys_fluxvec,swsys_flux_halo,flux_sw)
!=======================================================================================================!  
!=======================================================================================================!    
 sw_grad(:,:,:)= zero
!=======================================================================================================!    
!	Gradient  for the momentum equations
!=======================================================================================================!
 call dg3d_gradient_mom(deriv,energy,sw_grad(:,:,1),sw_grad(:,:,2))
!=======================================================================================================!    
!	Gradient  for the continuity equation     (pressure)
!=======================================================================================================!
 sw_grad(:,:,3)= dg3d_gradient_mass(deriv,swsys_fluxvec(:,:,:,3))
!=======================================================================================================!    
!	Gradient  for the continuity equation     (Pot-temp)
!=======================================================================================================!
 sw_grad(:,:,4)= dg3d_gradient_mass(deriv,swsys_fluxvec(:,:,:,4))
!=======================================================================================================!
!	Gradient  for the continuity equation     (Q-moisture)
!=======================================================================================================!
 sw_grad(:,:,5)= dg3d_gradient_mass(deriv,swsys_fluxvec(:,:,:,5))

!	Compute RHS of the ODE system corresponding  to the SW model         
!=======================================================================================================!
 sw_rhs(:,:,:)= zero
!
 do j= 1,np
 do i= 1,np       
  mmi(i,j) = 1.0D0 / elem%mp(i,j)
 enddo
 enddo

 do eqn=1,neqn
 do j= 1,np
 do i= 1,np       
    rtmp = mmi(i,j)
    sw_rhs(i,j,eqn) = (sw_source(i,j,eqn) + sw_grad(i,j,eqn) - flux_sw(i,j,eqn))*rtmp     
 enddo
 enddo
 enddo
!=======================================================================================================!
end subroutine  dg3d_rhs_terms
!=======================================================================================================!
!=======================================================================================================!
!	Divergence  (covariant input)  computation
!=======================================================================================================!
function  divergence_cov(deriv,mterm,couv) result(div)
   Implicit None

    type (derivative_t)         :: deriv
    real (kind=real_kind), dimension(np,np), intent(in) :: mterm
    real (kind=real_kind), dimension(np,np,2), intent(in) ::  couv  

    real (kind=real_kind), dimension(np,np) :: cdiv, div

    real(kind=real_kind) ::  dvdx00,dvdx01
    real(kind=real_kind) ::  dvdx10,dvdx11
    real(kind=real_kind) ::  dudy00,dudy01
    real(kind=real_kind) ::  dudy10,dudy11

    real (kind=real_kind) :: term,delm
    integer i,j,k,l    
    logical, parameter :: UseUnroll = .TRUE.
!=======================================================================================================!
   !  Covariant vorticity
   !
   !  div = [(u1)_x + (u2)_y] / G, computed in physical space by
   !                               collocation differentiation.
   !
if(MODULO(np,2) == 0 .and. UseUnroll) then 
    do j=1,np,2
       do l=1,np,2

          dudy00=zero
          dudy01=zero
          dudy10=zero
          dudy11=zero

          dvdx00=zero
          dvdx01=zero
          dvdx10=zero
          dvdx11=zero

          do i=1,np

             dvdx00 = dvdx00 + deriv%Dvv(i,l  )* couv(i  ,j,1)
             dvdx01 = dvdx01 + deriv%Dvv(i,l+1)* couv(i  ,j,1)
             dvdx10 = dvdx10 + deriv%Dvv(i,l  )* couv(i,j+1,1)
             dvdx11 = dvdx11 + deriv%Dvv(i,l+1)* couv(i,j+1,1)

             dudy00 = dudy00 + deriv%Dvv(i,l  )* couv(j  ,i,2)
             dudy01 = dudy01 + deriv%Dvv(i,l+1)* couv(j  ,i,2)
             dudy10 = dudy10 + deriv%Dvv(i,l  )* couv(j+1,i,2)
             dudy11 = dudy11 + deriv%Dvv(i,l+1)* couv(j+1,i,2)

          end do

          cdiv(l  ,j  ) = dvdx00
          cdiv(l+1,j  ) = dvdx01
          cdiv(l  ,j+1) = dvdx10
          cdiv(l+1,j+1) = dvdx11

          deriv%vvtemp(j  ,l  ) = dudy00
          deriv%vvtemp(j  ,l+1) = dudy01
          deriv%vvtemp(j+1,l  ) = dudy10
          deriv%vvtemp(j+1,l+1) = dudy11

        end do
    end do
else
    do j=1,np
       do l=1,np
          dudy00=zero
          dvdx00=zero
          do i=1,np
             dvdx00 = dvdx00 + deriv%Dvv(i,l  )* couv(i  ,j,1)
             dudy00 = dudy00 + deriv%Dvv(i,l  )* couv(j  ,i,2)
          end do
          cdiv(l  ,j  ) = dvdx00
          deriv%vvtemp(j  ,l  ) = dudy00
        end do
    end do
endif
   do j=1,np
       do i=1,np
           div(i,j)= rrearth*(cdiv(i,j) + deriv%vvtemp(i,j))  / mterm(i,j)
       end do
   end do


end function  divergence_cov
!=======================================================================================================!
!  Selective diffusion at various vertical levels 
!=======================================================================================================!
 function horizontal_diff(klev,deriv,sg,ginv,uv)  result(difuv)

    type (derivative_t)              :: deriv

    integer             , intent(in) :: klev 
    real(kind=real_kind), intent(in) :: sg(np,np)
    real(kind=real_kind), intent(in) :: ginv(2,2,np,np)
    real(kind=real_kind), intent(in) :: uv(np,np,2)
    real(kind=real_kind)             :: difuv(np,np,2)
    real(kind=real_kind)             :: damp              

    integer i,j 

    damp = 0.0D0 

   if (klev < 4) then
     damp = 5.0D05 
     difuv(:,:,:) = diffusion_uv(deriv,sg,ginv,uv) 
   else
      damp = dcof2  
      difuv(:,:,:) = biharmonic_diff(deriv,sg,ginv,uv) 
    ! damp = 2.5D05 
    ! difuv(:,:,:) = diffusion_uv(deriv,sg,ginv,uv) 
   endif

        do j=1,np
            do i=1,np
               difuv(i,j,1) = damp*difuv(i,j,1)
               difuv(i,j,2) = damp*difuv(i,j,2)
            enddo 
        enddo 

end function  horizontal_diff
!=======================================================================================================!
!=======================================================================================================!
!	Laplacian business  del^2 (u,v)
!=======================================================================================================!
 function diffusion_uv(deriv,sg,ginv,uv)  result(difuv)

    type (derivative_t)              :: deriv

    real(kind=real_kind), intent(in) :: sg(np,np)
    real(kind=real_kind), intent(in) :: ginv(2,2,np,np)
    real(kind=real_kind), intent(in) :: uv(np,np,2)
    real(kind=real_kind)             :: difuv(np,np,2)

    real(kind=real_kind),dimension(np,np) :: udiv, vdiv 
    real(kind=real_kind), dimension(np,np,2) :: gradu, gradv, grad_u, grad_v

    integer i
    integer j

    real(kind=real_kind) ::  v1, v2
!=======================================================================================================!
    !  Laplacian operatons for (u_1, u_2)

         gradu(:,:,:) = gradient_rhs(deriv,uv(:,:,1))
         gradv(:,:,:) = gradient_rhs(deriv,uv(:,:,2))

        do j=1,np
            do i=1,np
                   v1 = gradu(i,j,1)
                   v2 = gradu(i,j,2)
              grad_u(i,j,1) =  sg(i,j) *(ginv(1,1,i,j)*v1 + ginv(1,2,i,j)*v2)
              grad_u(i,j,2) =  sg(i,j) *(ginv(2,1,i,j)*v1 + ginv(2,2,i,j)*v2)
                   v1 = gradv(i,j,1)
                   v2 = gradv(i,j,2)
              grad_v(i,j,1) =  sg(i,j) *(ginv(1,1,i,j)*v1 + ginv(1,2,i,j)*v2)
              grad_v(i,j,2) =  sg(i,j) *(ginv(2,1,i,j)*v1 + ginv(2,2,i,j)*v2)
            end do
        end do

        udiv(:,:)  = diverge(deriv,grad_u(:,:,:))
        vdiv(:,:)  = diverge(deriv,grad_v(:,:,:))

       do j=1,np
           do i=1,np
              difuv(i,j,1) = udiv(i,j) / sg(i,j)
              difuv(i,j,2) = vdiv(i,j) / sg(i,j)
           end do
       end do

!=======================================================================================================!
end function diffusion_uv
!=======================================================================================================!
!=======================================================================================================!
 function implicit_diff(dtime,deriv,sg,ginv,tmp)  result(dift)

    type (derivative_t)              :: deriv

    real(kind=real_kind), intent(in) :: dtime 
    real(kind=real_kind), intent(in) :: sg(np,np)
    real(kind=real_kind), intent(in) :: ginv(2,2,np,np)
    real(kind=real_kind), intent(in) :: tmp(np,np)
    real(kind=real_kind)             :: dift(np,np)
    real(kind=real_kind)             :: dcoft        
    integer i
    integer j

     dcoft = 1.0D05 
     dift(:,:) = diffusion_temp(deriv,sg,ginv,tmp)

       do j=1,np
           do i=1,np
              dift(i,j) = tmp(i,j) - dtime*dcoft * dift(i,j)
           end do
       end do
!=======================================================================================================!
end function implicit_diff 
!=======================================================================================================!
!	Localized Laplacian del^2 (temp) for pot-temp 
!=======================================================================================================!
 function diffusion_theta(elem,deriv,dp,tmp)  result(dift)

    type (derivative_t)              :: deriv
    type(element_t), intent(in) :: elem

    real(kind=real_kind), intent(in) :: tmp(np,np), dp(np,np) 
    real(kind=real_kind)             :: dift(np,np)

    real(kind=real_kind) :: sg(np,np)
    real(kind=real_kind) :: ginv(2,2,np,np)

    real(kind=real_kind),dimension(np,np) :: udiv, vdiv
    real(kind=real_kind), dimension(np,np,2) :: gradt, gradv, grad_t, grad_v

    integer i
    integer j

    real(kind=real_kind) ::  v1, v2, d_tnu, tmp2
!=======================================================================================================!

    sg(:,:) = elem%metdet(:,:) 
    ginv(:,:,:,:) = elem%metinv(:,:,:,:) 

    !  Laplacian operatons for (u_1, u_2)

       ! MNL: adjustment because we're on ref elem, not cube
       if (ne.ne.0) then
           tmp2 = acos(0.0d0)/dble(2*ne)
           ! tmp2 = 40.0d0/acos(0.0d0)
           d_tnu = nu*tmp2*tmp2
       else
            d_tnu = nu
       end if

         gradt(:,:,:) = gradient_rhs(deriv,tmp(:,:))

        do j=1,np
            do i=1,np
                   v1 = gradt(i,j,1)
                   v2 = gradt(i,j,2)
              grad_t(i,j,1) = d_tnu* dp(i,j)* sg(i,j) *(ginv(1,1,i,j)*v1 + ginv(1,2,i,j)*v2)
              grad_t(i,j,2) = d_tnu* dp(i,j)* sg(i,j) *(ginv(2,1,i,j)*v1 + ginv(2,2,i,j)*v2)
            end do
        end do

        dift(:,:)  = diverge(deriv,grad_t(:,:,:))

       do j=1,np
           do i=1,np
              dift(i,j) = dift(i,j) / (sg(i,j))
           end do
       end do
!=======================================================================================================!
end function diffusion_theta
!=======================================================================================================!
!=======================================================================================================!
 function diffusion_hypr(elem,deriv,uv)  result(difuv)

    type (derivative_t)              :: deriv
    type(element_t), intent(in) :: elem

    real(kind=real_kind), intent(in) :: uv(np,np,2)
    real(kind=real_kind)             :: difuv(np,np,2)

    real(kind=real_kind) :: sg(np,np)
    real(kind=real_kind) :: ginv(2,2,np,np)

    real(kind=real_kind),dimension(np,np) :: udiv, vdiv 
    real(kind=real_kind), dimension(np,np,2) :: gradu, gradv, grad_u, grad_v

    integer i
    integer j

    real(kind=real_kind) ::  v1, v2, hy_nu, tmp
!=======================================================================================================!

    sg(:,:) = elem%metdet(:,:) 
    ginv(:,:,:,:) = elem%metinv(:,:,:,:) 

        ! MNL: adjustment because we're on ref elem, not cube
       if (ne.ne.0) then
           tmp = acos(0.0d0)/dble(2*ne)
           ! tmp = 40.0d0/acos(0.0d0)
           hy_nu =  2.0D08*tmp*tmp
        else
           hy_nu = 2.0D08
        end if

    !  Laplacian operatons for (u_1, u_2)

         gradu(:,:,:) = gradient_rhs(deriv,uv(:,:,1))
         gradv(:,:,:) = gradient_rhs(deriv,uv(:,:,2))

        do j=1,np
            do i=1,np
                   v1 = gradu(i,j,1)
                   v2 = gradu(i,j,2)
              grad_u(i,j,1) =  sg(i,j) *(ginv(1,1,i,j)*v1 + ginv(1,2,i,j)*v2)
              grad_u(i,j,2) =  sg(i,j) *(ginv(2,1,i,j)*v1 + ginv(2,2,i,j)*v2)
                   v1 = gradv(i,j,1)
                   v2 = gradv(i,j,2)
              grad_v(i,j,1) =  sg(i,j) *(ginv(1,1,i,j)*v1 + ginv(1,2,i,j)*v2)
              grad_v(i,j,2) =  sg(i,j) *(ginv(2,1,i,j)*v1 + ginv(2,2,i,j)*v2)
            end do
        end do

        udiv(:,:)  = diverge(deriv,grad_u(:,:,:))
        vdiv(:,:)  = diverge(deriv,grad_v(:,:,:))

       do j=1,np
           do i=1,np
              difuv(i,j,1) = hy_nu * udiv(i,j) / sg(i,j)
              difuv(i,j,2) = hy_nu * vdiv(i,j) / sg(i,j)
           end do
       end do

!=======================================================================================================!
end function diffusion_hypr
!=======================================================================================================!
 function diffusion_temp(deriv,sg,ginv,tmp)  result(dift)

    type (derivative_t)              :: deriv

    real(kind=real_kind), intent(in) :: sg(np,np)
    real(kind=real_kind), intent(in) :: ginv(2,2,np,np)
    real(kind=real_kind), intent(in) :: tmp(np,np)
    real(kind=real_kind)             :: dift(np,np)

    real(kind=real_kind),dimension(np,np) :: udiv, vdiv
    real(kind=real_kind), dimension(np,np,2) :: gradt, gradv, grad_t, grad_v

    integer i
    integer j

    real(kind=real_kind) ::  v1, v2
!=======================================================================================================!
    !  Laplacian operatons for (u_1, u_2)

         gradt(:,:,:) = gradient_rhs(deriv,tmp(:,:))

        do j=1,np
            do i=1,np
                   v1 = gradt(i,j,1)
                   v2 = gradt(i,j,2)
              grad_t(i,j,1) =  sg(i,j) *(ginv(1,1,i,j)*v1 + ginv(1,2,i,j)*v2)
              grad_t(i,j,2) =  sg(i,j) *(ginv(2,1,i,j)*v1 + ginv(2,2,i,j)*v2)
            end do
        end do

        dift(:,:)  = diverge(deriv,grad_t(:,:,:))

       do j=1,np
           do i=1,np
              dift(i,j) = dift(i,j) / sg(i,j)
           end do
       end do
!=======================================================================================================!
end function diffusion_temp
!=======================================================================================================!
!=======================================================================================================!
!	Laplacian of the Laplacian  del^4 (u,v)
!=======================================================================================================!
 function biharmonic_diff(deriv,sg,ginv,uv)  result(difuv)

    type (derivative_t)              :: deriv

    real(kind=real_kind), intent(in) :: sg(np,np)
    real(kind=real_kind), intent(in) :: ginv(2,2,np,np)
    real(kind=real_kind), intent(in) :: uv(np,np,2)
    real(kind=real_kind)             :: difuv(np,np,2), lapuv(np,np,2)

    real(kind=real_kind),dimension(np,np) :: udiv, vdiv 
    real(kind=real_kind), dimension(np,np,2) :: gradu, gradv, grad_u, grad_v

    integer i
    integer j

    real(kind=real_kind) ::  v1, v2
!=======================================================================================================!
    !  Laplacian operatons for (u_1, u_2)

         gradu(:,:,:) = gradient_rhs(deriv,uv(:,:,1))
         gradv(:,:,:) = gradient_rhs(deriv,uv(:,:,2))

        do j=1,np
            do i=1,np
                   v1 = gradu(i,j,1)
                   v2 = gradu(i,j,2)
              grad_u(i,j,1) =  sg(i,j) *(ginv(1,1,i,j)*v1 + ginv(1,2,i,j)*v2)
              grad_u(i,j,2) =  sg(i,j) *(ginv(2,1,i,j)*v1 + ginv(2,2,i,j)*v2)
                   v1 = gradv(i,j,1)
                   v2 = gradv(i,j,2)
              grad_v(i,j,1) =  sg(i,j) *(ginv(1,1,i,j)*v1 + ginv(1,2,i,j)*v2)
              grad_v(i,j,2) =  sg(i,j) *(ginv(2,1,i,j)*v1 + ginv(2,2,i,j)*v2)
            end do
        end do

        udiv(:,:)  = diverge(deriv,grad_u(:,:,:))
        vdiv(:,:)  = diverge(deriv,grad_v(:,:,:))

       do j=1,np
           do i=1,np
              lapuv(i,j,1) = udiv(i,j) / sg(i,j)
              lapuv(i,j,2) = vdiv(i,j) / sg(i,j)
           end do
       end do

!    begin Biharmonic step 

         gradu(:,:,:) = gradient_rhs(deriv,lapuv(:,:,1))
         gradv(:,:,:) = gradient_rhs(deriv,lapuv(:,:,2))

        do j=1,np
            do i=1,np
                   v1 = gradu(i,j,1)
                   v2 = gradu(i,j,2)
              grad_u(i,j,1) =  sg(i,j) *(ginv(1,1,i,j)*v1 + ginv(1,2,i,j)*v2)
              grad_u(i,j,2) =  sg(i,j) *(ginv(2,1,i,j)*v1 + ginv(2,2,i,j)*v2)
                   v1 = gradv(i,j,1)
                   v2 = gradv(i,j,2)
              grad_v(i,j,1) =  sg(i,j) *(ginv(1,1,i,j)*v1 + ginv(1,2,i,j)*v2)
              grad_v(i,j,2) =  sg(i,j) *(ginv(2,1,i,j)*v1 + ginv(2,2,i,j)*v2)
            end do
        end do

        udiv(:,:)  = diverge(deriv,grad_u(:,:,:))
        vdiv(:,:)  = diverge(deriv,grad_v(:,:,:))

       do j=1,np
           do i=1,np
              difuv(i,j,1) = udiv(i,j) / sg(i,j)
              difuv(i,j,2) = vdiv(i,j) / sg(i,j)
           end do
       end do
 end function biharmonic_diff
!=======================================================================================================!
 function cube_laplacian (deriv,fld)  result(diff)

    type (derivative_t)              :: deriv

   !real(kind=real_kind), intent(in) :: sg(np,np)
   !real(kind=real_kind), intent(in) :: ginv(2,2,np,np)
    real(kind=real_kind), intent(in) :: fld(np,np)
    real(kind=real_kind)             :: diff(np,np)

    real(kind=real_kind),dimension(np,np) :: udiv, vdiv
    real(kind=real_kind), dimension(np,np,2) :: gradt, gradv, grad_t, grad_v

    integer i
    integer j

    real(kind=real_kind) ::  v1, v2
!=======================================================================================================!
    !  Laplacian operatons for (u_1, u_2) on the "cube"

         gradt(:,:,:) = gradient_rhs(deriv,fld(:,:))


          diff(:,:)  = diverge(deriv,gradt(:,:,:))

!=======================================================================================================!
end function cube_laplacian 
!=======================================================================================================!
!=======================================================================================================!
 function diverge(deriv,fuv)  result(div)

    type (derivative_t)              :: deriv

    real(kind=real_kind), intent(in) :: fuv(np,np,2)

    real(kind=real_kind)             :: div(np,np)

    integer i
    integer j
    integer l

    real(kind=real_kind) ::  dvdx00,dvdx01
    real(kind=real_kind) ::  dvdx10,dvdx11
    real(kind=real_kind) ::  dudy00,dudy01
    real(kind=real_kind) ::  dudy10,dudy11    
    logical, parameter :: UseUnroll = .TRUE.
!=======================================================================================================!
    !  divergence
    !
    !  div = [(u1)_x + (u2)_y] / G, computed in physical space by
    !                               collocation differentiation.
    !
if(MODULO(np,2) == 0 .and. UseUnroll) then 
   do j=1,np,2
      do l=1,np,2

   !do j=1,np
   !   do l=1,np
          dudy00=zero
          dudy01=zero
          dudy10=zero
          dudy11=zero

          dvdx00=zero
          dvdx01=zero
          dvdx10=zero
          dvdx11=zero

          do i=1,np

             dvdx00 = dvdx00 + deriv%Dvv(i,l  )* fuv(i  ,j,1)
             dvdx01 = dvdx01 + deriv%Dvv(i,l+1)* fuv(i  ,j,1)
             dvdx10 = dvdx10 + deriv%Dvv(i,l  )* fuv(i,j+1,1)
             dvdx11 = dvdx11 + deriv%Dvv(i,l+1)* fuv(i,j+1,1)

             dudy00 = dudy00 + deriv%Dvv(i,l  )* fuv(j  ,i,2)
             dudy01 = dudy01 + deriv%Dvv(i,l+1)* fuv(j  ,i,2)
             dudy10 = dudy10 + deriv%Dvv(i,l  )* fuv(j+1,i,2)
             dudy11 = dudy11 + deriv%Dvv(i,l+1)* fuv(j+1,i,2)

          end do

          div(l  ,j  ) = dvdx00
          div(l+1,j  ) = dvdx01
          div(l  ,j+1) = dvdx10
          div(l+1,j+1) = dvdx11

          deriv%vvtemp(j  ,l  ) = dudy00
          deriv%vvtemp(j  ,l+1) = dudy01
          deriv%vvtemp(j+1,l  ) = dudy10
          deriv%vvtemp(j+1,l+1) = dudy11

       end do
    end do
else
   do j=1,np
      do l=1,np
          dudy00=zero
          dvdx00=zero
          do i=1,np
             dvdx00 = dvdx00 + deriv%Dvv(i,l  )* fuv(i  ,j,1)
             dudy00 = dudy00 + deriv%Dvv(i,l  )* fuv(j  ,i,2)
          end do
          div(l  ,j  ) = dvdx00
          deriv%vvtemp(j  ,l  ) = dudy00
       end do
    end do
endif

    do j=1,np
       do i=1,np
          div(i,j)=(div(i,j) + deriv%vvtemp(i,j))*rrearth
       end do
    end do
!=======================================================================================================!
  end function diverge
!=======================================================================================================!
!=======================================================================================================!
  function cov_vorticity(deriv,couv)  result(vor)

    type (derivative_t)              :: deriv

    real(kind=real_kind), intent(in) :: couv(np,np,2)

    real(kind=real_kind)             :: vor(np,np)

    integer i
    integer j
    integer l
    real(kind=real_kind) ::  dvdx00,dvdx01
    real(kind=real_kind) ::  dvdx10,dvdx11
    real(kind=real_kind) ::  dudy00,dudy01
    real(kind=real_kind) ::  dudy10,dudy11    
    logical, parameter :: UseUnroll = .TRUE.
!=======================================================================================================!
    !  Covariant vorticity
    !
    !  vor = [(u2)_x - (u1)_y] / G
    !
    !
if(MODULO(np,2) == 0 .and. UseUnroll) then 
    do j=1,np,2
       do l=1,np,2

          dudy00=zero
          dudy01=zero
          dudy10=zero
          dudy11=zero

          dvdx00=zero
          dvdx01=zero
          dvdx10=zero
          dvdx11=zero

          do i=1,np

             dvdx00 = dvdx00 + deriv%Dvv(i,l  )* couv(i  ,j,2)
             dvdx01 = dvdx01 + deriv%Dvv(i,l+1)* couv(i  ,j,2)
             dvdx10 = dvdx10 + deriv%Dvv(i,l  )* couv(i,j+1,2)
             dvdx11 = dvdx11 + deriv%Dvv(i,l+1)* couv(i,j+1,2)

             dudy00 = dudy00 + deriv%Dvv(i,l  )* couv(j  ,i,1)
             dudy01 = dudy01 + deriv%Dvv(i,l+1)* couv(j  ,i,1)
             dudy10 = dudy10 + deriv%Dvv(i,l  )* couv(j+1,i,1)
             dudy11 = dudy11 + deriv%Dvv(i,l+1)* couv(j+1,i,1)

          end do

          vor(l  ,j  ) = dvdx00
          vor(l+1,j  ) = dvdx01
          vor(l  ,j+1) = dvdx10
          vor(l+1,j+1) = dvdx11

          deriv%vvtemp(j  ,l  ) = dudy00
          deriv%vvtemp(j  ,l+1) = dudy01
          deriv%vvtemp(j+1,l  ) = dudy10
          deriv%vvtemp(j+1,l+1) = dudy11

       end do
    end do
else
    do j=1,np
       do l=1,np
          dudy00=zero
	  dvdx00=zero
	  do i=1,np
             dvdx00 = dvdx00 + deriv%Dvv(i,l  )* couv(i  ,j,2)
	     dudy00 = dudy00 + deriv%Dvv(i,l  )* couv(j  ,i,1)
 	  end do
          vor(l  ,j  ) = dvdx00
	  deriv%vvtemp(j  ,l  ) = dudy00
       end do
    end do
endif

    do j=1,np
       do i=1,np
          vor(i,j)=(vor(i,j)-deriv%vvtemp(i,j))*rrearth
       end do
    end do
!=======================================================================================================!
 end function cov_vorticity
!=======================================================================================================!
!=======================================================================================================!
function gradient_rhs(deriv,ff)  result(gradf)

    type (derivative_t)              :: deriv

    real(kind=real_kind), intent(in) :: ff(np,np)

    real(kind=real_kind)             :: gradf(np,np,2)
    real (kind=real_kind), dimension(np,np)   :: phi1, phi2

    integer i
    integer j
    integer l

    real(kind=real_kind) ::  dvdx00,dvdx01
    real(kind=real_kind) ::  dvdx10,dvdx11
    real(kind=real_kind) ::  dudy00,dudy01
    real(kind=real_kind) ::  dudy10,dudy11    
    logical, parameter :: UseUnroll = .TRUE.
!=======================================================================================================!
  !Geopotential ht gradients
if(MODULO(np,2) == 0 .and. UseUnroll) then 
    do j=1,np,2
       do l=1,np,2

          dudy00=zero
          dudy01=zero
          dudy10=zero
          dudy11=zero


          dvdx00=zero
          dvdx01=zero
          dvdx10=zero
          dvdx11=zero

          do i=1,np

             dvdx00 = dvdx00 + deriv%Dvv(i,l  )* ff(i  ,j)
             dvdx01 = dvdx01 + deriv%Dvv(i,l+1)* ff(i  ,j)
             dvdx10 = dvdx10 + deriv%Dvv(i,l  )* ff(i,j+1)
             dvdx11 = dvdx11 + deriv%Dvv(i,l+1)* ff(i,j+1)

             dudy00 = dudy00 + deriv%Dvv(i,l  )* ff(j  ,i)
             dudy01 = dudy01 + deriv%Dvv(i,l+1)* ff(j  ,i)
             dudy10 = dudy10 + deriv%Dvv(i,l  )* ff(j+1,i)
             dudy11 = dudy11 + deriv%Dvv(i,l+1)* ff(j+1,i)

          end do

          phi1(l  ,j  ) = dvdx00
          phi1(l+1,j  ) = dvdx01
          phi1(l  ,j+1) = dvdx10
          phi1(l+1,j+1) = dvdx11
          phi2(j  ,l  ) = dudy00
          phi2(j  ,l+1) = dudy01
          phi2(j+1,l  ) = dudy10
          phi2(j+1,l+1) = dudy11

        end do
    end do
else
    do j=1,np
       do l=1,np

          dudy00=zero
	  dvdx00=zero

          do i=1,np

             dvdx00 = dvdx00 + deriv%Dvv(i,l  )* ff(i  ,j)
	     dudy00 = dudy00 + deriv%Dvv(i,l  )* ff(j  ,i)

          end do

          phi1(l  ,j  ) = dvdx00
	  phi2(j  ,l  ) = dudy00

        end do
    end do
endif

    do j=1,np
       do i=1,np
           gradf(i,j,1)= phi1(i,j)*rrearth
           gradf(i,j,2)= phi2(i,j)*rrearth
       end do
    end do
!=======================================================================================================!
 end function gradient_rhs
!=======================================================================================================!
!=======================================================================================================!
! Source term  computation
!=======================================================================================================!
subroutine source_term(elem,deriv,dx,dy,gcori,contrauv,couv,prg,source)
   Implicit None

    type(element_t), intent(in)         :: elem
    type (derivative_t)         :: deriv
    real (kind=real_kind), intent(in) :: dx, dy
    real (kind=real_kind), dimension(np,np), intent(in) :: gcori
    real (kind=real_kind), dimension(np,np,2), intent(in) :: contrauv, couv, prg
    real (kind=real_kind), dimension(np,np,4), intent(out) :: source     

    real (kind=real_kind), dimension(np,np,2) :: gradp
    real (kind=real_kind), dimension(np,np) :: vor, vort, phi1,phi2

    real(kind=real_kind) ::  dvdx00,dvdx01
    real(kind=real_kind) ::  dvdx10,dvdx11
    real(kind=real_kind) ::  dudy00,dudy01
    real(kind=real_kind) ::  dudy10,dudy11

    real (kind=real_kind) :: term,delm
    integer i,j,k,l    
    logical, parameter :: UseUnroll = .TRUE.
!=======================================================================================================!
if(MODULO(np,2) == 0 .and. UseUnroll) then 
    do j=1,np,2
       do l=1,np,2

          dudy00=zero
          dudy01=zero
          dudy10=zero
          dudy11=zero

          dvdx00=zero
          dvdx01=zero
          dvdx10=zero
          dvdx11=zero

          do i=1,np

             dvdx00 = dvdx00 + deriv%Dvv(i,l  )* couv(i  ,j,2)
             dvdx01 = dvdx01 + deriv%Dvv(i,l+1)* couv(i  ,j,2)
             dvdx10 = dvdx10 + deriv%Dvv(i,l  )* couv(i,j+1,2)
             dvdx11 = dvdx11 + deriv%Dvv(i,l+1)* couv(i,j+1,2)

             dudy00 = dudy00 + deriv%Dvv(i,l  )* couv(j  ,i,1)
             dudy01 = dudy01 + deriv%Dvv(i,l+1)* couv(j  ,i,1)
             dudy10 = dudy10 + deriv%Dvv(i,l  )* couv(j+1,i,1)
             dudy11 = dudy11 + deriv%Dvv(i,l+1)* couv(j+1,i,1)

          end do

          vor(l  ,j  ) = dvdx00
          vor(l+1,j  ) = dvdx01
          vor(l  ,j+1) = dvdx10
          vor(l+1,j+1) = dvdx11

          deriv%vvtemp(j  ,l  ) = dudy00
          deriv%vvtemp(j  ,l+1) = dudy01
          deriv%vvtemp(j+1,l  ) = dudy10
          deriv%vvtemp(j+1,l+1) = dudy11

        end do
    end do
else
    do j=1,np
       do l=1,np
          dudy00=zero
	  dvdx00=zero
          do i=1,np
             dvdx00 = dvdx00 + deriv%Dvv(i,l  )* couv(i  ,j,2)
	     dudy00 = dudy00 + deriv%Dvv(i,l  )* couv(j  ,i,1)
          end do
          vor(l  ,j  ) = dvdx00
	  deriv%vvtemp(j  ,l  ) = dudy00
        end do
    end do
endif

    do j=1,np
       do i=1,np
           vor(i,j)=(vor(i,j)*dy - deriv%vvtemp(i,j)*dx) 
       end do
    end do


    source(:,:,:) = zero 


    do j = 1,np
      do i = 1,np

    !  term =  (vor(i,j) +  gcori(i,j)) * mmx(i,j) 
    !  source(i,j,1) =  contrauv(i,j,2) * term
    !  source(i,j,2) = -contrauv(i,j,1) * term  

       term =  (vor(i,j) +  gcori(i,j)) 
    ! source(i,j,1) =  (contrauv(i,j,2) * term - gradp(i,j,1)  ) * mmx(i,j)
    ! source(i,j,2) = -(contrauv(i,j,1) * term + gradp(i,j,2)  ) * mmx(i,j)

      source(i,j,1) =  (contrauv(i,j,2) * term - prg(i,j,1) ) * mmx(i,j)
      source(i,j,2) = -(contrauv(i,j,1) * term + prg(i,j,2) ) * mmx(i,j)

      enddo
    enddo

  end subroutine  source_term
!=======================================================================================================! 
function gradient_mass(deriv,dx,dy,uvflx)  result(gradf)
   
    type (derivative_t)              :: deriv

    real(kind=real_kind), intent(in) :: dx
    real(kind=real_kind), intent(in) :: dy
    real(kind=real_kind), intent(in) :: uvflx(np,np,2)

    real(kind=real_kind)             :: gradf(np,np)

    ! Dvv_twt ->  transpose of derivax *G-weights
    ! Mvv_twt ->  Diadgonal matrix with G-weights

    integer i
    integer j
    integer l

    real(kind=real_kind)  sumx00,sumx01
    real(kind=real_kind)  sumy00,sumy01
    real(kind=real_kind)  sumx10,sumx11
    real(kind=real_kind)  sumy10,sumy11
!=======================================================================================================!

    do j=1,np
       do l=1,np
          sumx00=zero
          sumy00=zero

          do i=1,np
             sumx00 = sumx00 + deriv%Dvv_twt(i,l) * uvflx(i,j,1)
             sumy00 = sumy00 + deriv%mvv_twt(i,l) * uvflx(i,j,2)
          end do

          deriv%vvtempt(j  ,l  ,1) = sumx00
          deriv%vvtempt(j  ,l  ,2) = sumy00

       end do
    end do

    do j=1,np
       do i=1,np
          sumx00=zero
          sumy00=zero

          do l=1,np
             sumx00 = sumx00 +  deriv%mvv_twt(l,j  )*deriv%vvtempt(l,i,1)
             sumy00 = sumy00 +  deriv%Dvv_twt(l,j  )*deriv%vvtempt(l,i,2)
          end do

           gradf(i,j) = dy*sumx00 + dx*sumy00
       end do
    end do


end function  gradient_mass
!=======================================================================================================!
!=======================================================================================================!
subroutine gradient_mom(deriv,dx,dy,energy,gradu1,gradu2)
   
    type (derivative_t)              :: deriv

    real(kind=real_kind), intent(in) :: dx
    real(kind=real_kind), intent(in) :: dy
    real(kind=real_kind), intent(in) :: energy(np,np)

    real(kind=real_kind), intent(out):: gradu1(np,np), gradu2(np,np)

    ! Dvv_twt ->  transpose of derivax *G-weights
    ! Mvv_twt ->  Diadgonal matrix with G-weights

    integer i
    integer j
    integer l

    real(kind=real_kind)  sumx00,sumx01
    real(kind=real_kind)  sumy00,sumy01
    real(kind=real_kind)  sumx10,sumx11
    real(kind=real_kind)  sumy10,sumy11
!=======================================================================================================!

   !Grad-u 

    do j=1,np
       do l=1,np
          sumx00=zero

          do i=1,np
             sumx00 = sumx00 + deriv%Dvv_twt(i,l) * energy(i,j)
          end do

          deriv%vvtempt(j  ,l  ,1) = sumx00

       end do
    end do

    do j=1,np
       do i=1,np
          sumx00=zero
          do l=1,np
             sumx00 = sumx00 +  deriv%mvv_twt(l,j  )*deriv%vvtempt(l,i,1)
          end do
           gradu1(i,j) = dy*sumx00 
       end do
    end do

   !Grad-v 

    do j=1,np
       do l=1,np
          sumy00=zero
          do i=1,np
             sumy00 = sumy00 + deriv%mvv_twt(i,l) * energy(i,j)
          end do
          deriv%vvtempt(j  ,l  ,2) = sumy00
       end do
    end do

    do j=1,np
       do i=1,np
          sumy00=zero
          do l=1,np
             sumy00 = sumy00 +  deriv%Dvv_twt(l,j  )*deriv%vvtempt(l,i,2)
          end do
           gradu2(i,j) = dx*sumy00
       end do
    end do

end subroutine  gradient_mom
!=======================================================================================================!
!=======================================================================================================!
!=======================================================================================================!
subroutine pres_grad_term(elem,klev,dtt,deriv,lnpr,temp,force,pgrad)
!=======================================================================================================! 
 implicit none
 
 type(element_t), intent(in) :: elem
 integer,              intent(in):: klev
 type (derivative_t),  intent(in):: deriv
 real (kind=real_kind), intent(in)  :: dtt  
 real (kind=real_kind), dimension(np,np),intent(in)  :: lnpr,temp, force
 real (kind=real_kind), dimension(np,np,2),intent(out):: pgrad 
!=======================================================================================================!
 real (kind=real_kind), dimension(np,np,2):: lpgrad
 real (kind=real_kind), dimension(np,np)  :: rdtm, tdif, dtmp
 real (kind=real_kind), parameter  :: p_0 = 100000.0D0   ! Initial Surface pressure
 real (kind=real_kind), parameter  :: r_d = 287.04D0      ! Gas const (dry)
 real (kind=real_kind), parameter  :: c_p = 1004.64D0    ! Cp
 real (kind=real_kind), parameter  :: kapa= r_d/c_p
 real (kind=real_kind)             :: damp           
 integer:: i,j
!=======================================================================================================!
 tdif(:,:) = 0.0D0 
!=======================================================================================================!
!tdif(:,:)= diffusion_temp(deriv,elem%metdet,elem%metinv,temp(:,:)) 

 if (klev < 3) then
    damp = 5.0D05
     tdif(:,:) = cube_laplacian(deriv,temp(:,:))
 else
    damp = 0.0D0  
  ! tdif(:,:)= diffusion_temp(deriv,elem%metdet,elem%metinv,temp(:,:)) 
 endif 

  do j=1,np
  do i=1,np
       dtmp(i,j)=  temp(i,j) - damp*dtt*tdif(i,j)        
  end do
  end do

 do j=1,np
 do i=1,np
       rdtm(i,j)= r_d*(dtmp(i,j) - force(i,j)) 
     ! rdtm(i,j)= r_d* temp(i,j) / (2.0D0 * lnpr(i,j))
     ! dtmp(i,j) = lnpr(i,j) * lnpr(i,j) 
 end do
 end do
!=======================================================================================================!
 lpgrad(:,:,:)= gradient_rhs(deriv,lnpr(:,:))

 pgrad= zero

 do j=1,np
 do i=1,np
    pgrad(i,j,1) = lpgrad(i,j,1) * rdtm(i,j)  
    pgrad(i,j,2) = lpgrad(i,j,2) * rdtm(i,j) 
 end do
 end do
!=======================================================================================================!
end subroutine pres_grad_term  
!=======================================================================================================!
!	A general gradient in (a^1,a^2) 
!=======================================================================================================!
 function general_grad(grad,g,ginv)  result(ggrad)

    real(kind=real_kind), intent(in) :: g(2,2,np,np), ginv(2,2,np,np)  
    real(kind=real_kind), intent(in) :: grad(np,np,2)
    real(kind=real_kind)             :: ggrad(np,np,2)

    integer  :: i, j 
    real(kind=real_kind) ::  v1, v2

        do j=1,np
            do i=1,np
                   v1 = grad(i,j,1)
                   v2 = grad(i,j,2)
             ! ggrad(i,j,1) =  g(1,1,i,j) *(ginv(1,1,i,j)*v1 + ginv(1,2,i,j)*v2)
             ! ggrad(i,j,2) =  g(2,2,i,j) *(ginv(2,1,i,j)*v1 + ginv(2,2,i,j)*v2)
              ggrad(i,j,1) =  (ginv(1,1,i,j)*v1 + ginv(1,2,i,j)*v2)
              ggrad(i,j,2) =  (ginv(2,1,i,j)*v1 + ginv(2,2,i,j)*v2)
            end do
        end do

 end function general_grad 
!=======================================================================================================!
subroutine gradient_p3d(elem,deriv,p3d,t3d,pgrad3)
!=======================================================================================================! 
 implicit none
 
 type(element_t), intent(in) :: elem
 type (derivative_t),  intent(in):: deriv
 real (kind=real_kind), dimension(np,np,nlev+1),intent(in)  :: p3d 
 real (kind=real_kind), dimension(np,np,nlev),intent(in)    :: t3d 
 real (kind=real_kind), dimension(np,np,2,nlev),intent(out) :: pgrad3
!=======================================================================================================!
 real (kind=real_kind), dimension(np,np,nlev+1):: lnp3d 
 real (kind=real_kind), dimension(np,np)  :: lnpr, lsum1, lsum2, deno
 real (kind=real_kind), dimension(np,np,2):: lpgrad
 real (kind=real_kind), dimension(np,np,2,nlev+1):: lpg3
 real (kind=real_kind), dimension(np,np,nlev+1):: ln2p 
 real (kind=real_kind):: trm
 real (kind=real_kind), parameter  :: p_0 = 100000.0D0   ! Initial Surface pressure
 real (kind=real_kind), parameter  :: r_d = 287.04D0      ! Gas const (dry)
 real (kind=real_kind), parameter  :: c_p = 1004.64D0    ! Cp
 real (kind=real_kind), parameter  :: kapa= r_d/c_p
 integer:: i,j,k ,l 
!=======================================================================================================!
  do k=1,nlev+1
  do j=1,np
  do i=1,np
   lnp3d(i,j,k) = log(p3d(i,j,k))
   ln2p(i,j,k) = lnp3d(i,j,k) * lnp3d(i,j,k)
  end do
  end do
    lnpr(:,:) = ln2p(:,:,k)
    lpgrad(:,:,:)= gradient_rhs(deriv,lnpr(:,:))
    lpg3(:,:,:,k) = lpgrad(:,:,:)
  end do

  do k=1,nlev

  do j=1,np
  do i=1,np
   !lnpr(i,j) = (lnp3d(i,j,k) + lnp3d(i,j,k+1))*0.5D0 
    lsum1(i,j) = (lpg3(i,j,1,k) + lpg3(i,j,1,k+1))*0.5D0 
    lsum2(i,j) = (lpg3(i,j,2,k) + lpg3(i,j,2,k+1))*0.5D0 
    deno(i,j) = (lnp3d(i,j,k) + lnp3d(i,j,k+1))*0.5D0 
  end do
  end do

  ! lpgrad(:,:,:)= gradient_rhs(deriv,lnpr(:,:))

 !do l = 1, 2
  do j=1,np
   do i=1,np
    !pgrad3(i,j,l,k) = lpgrad(i,j,l) *r_d * t3d(i,j,k)
         trm = 0.5D0 * r_d * t3d(i,j,k) / deno(i,j)
     pgrad3(i,j,1,k) = lsum1(i,j) * trm
     pgrad3(i,j,2,k) = lsum2(i,j) * trm
   end do
  end do
 !end do

  end do

!=======================================================================================================!
end subroutine gradient_p3d  
!=======================================================================================================!
!==================== LDG ZONE =========================================================================
!=======================================================================================!
subroutine dg3d_cent_flux(elem,gradbuf,grad,cgrad)
!=======================================================================================!
    Implicit None
    integer, parameter :: south=1, east=2, north=3, west=4
    type(element_t), intent(in) :: elem
    real (kind=real_kind), dimension(0:np+1,0:np+1,2), intent(in)   :: gradbuf
    real (kind=real_kind), dimension(np,np,2),         intent(in)   :: grad
    real (kind=real_kind), dimension(np,np,2),         intent(out)  :: cgrad

    real (kind=real_kind), dimension(np,4,2) :: grad_senw

    real (kind=real_kind), dimension(np)     :: cflx_south,cflx_north,cflx_east,cflx_west
    real(kind=real_kind) ::   f_left, f_right, s1,s2, x0y0
    integer i,j,k, west_edge, south_edge

!=======================================================================================================!
!  Boundary values for velocity & flux terms from the neighbors                                         !
!  H-field (scalar) on S, E, N, W sides of an element                                                   !
!  U,V-component on  (S,E,N,W) =  (1,2,3,4) sides of an element                                         !
!=======================================================================================================!
    do i = 1,np
    do j = 1, 2
       grad_senw(i,south,j)= gradbuf(i,   0,j)
       grad_senw(i,east,j) = gradbuf(np+1,i,j)
       grad_senw(i,north,j)= gradbuf(i,np+1,j)
       grad_senw(i,west,j) = gradbuf(0,i,   j)
    enddo
    enddo

 ! Copying initial grads

    do k = 1, 2
    do j = 1,np
    do i = 1,np
      cgrad(i,j,k) = grad(i,j,k)
    enddo
    enddo
    enddo

         x0y0 = -0.25D0 * dd_pi
         west_edge = 0
         south_edge = 0


    !! West/South edge of the Cubed-Sphere detection

          if (elem%cartp(1,1)%x == x0y0) west_edge=1
          if (elem%cartp(1,1)%y == x0y0) south_edge=1


    !For "Centered/sided" Flux for the grad  vectors

           ! Boundary  Central-flux    for x^1 terms

         do j = 1,np
                  f_left  = grad_senw(j,west,1)
                  f_right = grad(1,j,1)
              if (west_edge == 1) f_left = f_right
                cflx_west(j) = 0.5D0 *(f_right + f_left)
               !cflx_west(j) = f_left
           !        cgrad(1,j,1) =  cflx_west(j)

                  f_left  = grad(np,j,1)
                  f_right = grad_senw(j,east,1)
                cflx_east(j) = 0.5D0 *(f_right + f_left)
               !cflx_east(j) = f_left
           !        cgrad(np,j,1) =  cflx_east(j)

                   f_left = grad_senw(j,south,1)
                  f_right = grad(j,1,1)
              if (south_edge == 1) f_left = f_right
               cflx_south(j) = 0.5D0 *(f_right + f_left)
              !cflx_south(j) = f_left
           !        cgrad(j,1,1) =  cflx_south(j)

                  f_left  = grad(j,np,1)
                  f_right = grad_senw(j,north,1)
               cflx_north(j) = 0.5D0 *(f_right + f_left)
              !cflx_north(j) = f_left
           !        cgrad(j,np,1) =  cflx_north(j)
          end do

           ! 4-Boundary   Central-flux for x^2 term

         do i = 1,np
                  f_left  = grad_senw(i,west,2)
                  f_right = grad(1,i,2)
              if (west_edge == 1) f_left = f_right
                cflx_west(i) = 0.5D0 *(f_right + f_left)
               !cflx_west(i) = f_left
           !        cgrad(1,i,2) =  cflx_west(i)

                  f_left  = grad(np,i,2)
                  f_right = grad_senw(i,east,2)
                cflx_east(i) = 0.5D0 *(f_right + f_left)
               !cflx_east(i) = f_left
           !        cgrad(np,i,2) =  cflx_east(i)

                   f_left = grad_senw(i,south,2)
                  f_right = grad(i,1,2)
              if (south_edge == 1) f_left = f_right
               cflx_south(i) = 0.5D0 *(f_right + f_left)
              !cflx_south(i) = f_left
           !     cgrad(i,1,2) =  cflx_south(i)
                  f_left  = grad(i,np,2)
                  f_right = grad_senw(i,north,2)
               cflx_north(i) = 0.5D0 *(f_right + f_left)
              !cflx_north(i) = f_left
           !   cgrad(i,np,2) =  cflx_north(i)
         end do

!=======================================================================================================!
end subroutine dg3d_cent_flux
!=======================================================================================================!
!=++++++++++++++++======================================================================================!
! LDG & Diffusion Zone
!=++++++++++++++++======================================================================================!
subroutine dg3d_diff_flux(elem,deriv,gradbuf,grads,visflx)
!=======================================================================================!
    Implicit None
    integer, parameter :: south=1, east=2, north=3, west=4
    real (kind=real_kind), parameter :: old_d_nu = 5.0D4 !1.0D5 !2.0D5 
    type(element_t), intent(in) :: elem
    type (derivative_t),                               intent(in)   :: deriv
    real (kind=real_kind), dimension(0:np+1,0:np+1,2), intent(in)   :: gradbuf
    real (kind=real_kind), dimension(np,np,2),         intent(in)   :: grads
    real (kind=real_kind), dimension(np,np)  ,         intent(out)  :: visflx

    real (kind=real_kind), dimension(np,4,2)   :: grad_halo
    real (kind=real_kind), dimension(np,np)    :: dgrad, cflxint
    real (kind=real_kind) :: d_nu, tmp
    Integer :: i,j,k

!=======================================================================================================!
!=======================================================================================================!
!=======================================================================================================!
!=======================================================================================================!

    if (ne.ne.0) then
        tmp = dble(2*ne)/acos(0.0d0)
        d_nu = nu*tmp*tmp
    else
        d_nu = nu
    end if
    
    !print *, 'Running dg3d_diff_flux with nu = ', d_nu

!=======================================================================================================!
!  Boundary values for velocity & flux terms from the neighbors                                         !
!  H-field (scalar) on S, E, N, W sides of an element                                                   !
!  U,V-component on  (S,E,N,W) =  (1,2,3,4) sides of an element                                         !
!=======================================================================================================!
    do i = 1,np
    do j = 1, 2
       grad_halo(i,south,j)= gradbuf(i,   0,j)
       grad_halo(i,east,j) = gradbuf(np+1,i,j)
       grad_halo(i,north,j)= gradbuf(i,np+1,j)
       grad_halo(i,west,j) = gradbuf(0,i,   j)
    enddo
    enddo

    Call  central_fluxint(elem,deriv,grads,grad_halo,cflxint)

    dgrad(:,:) = dg3d_gradient_mass(deriv,grads)

    do j = 1,np
    do i = 1,np
      visflx(i,j) = d_nu * (cflxint(i,j) - dgrad(i,j)) / elem%mp(i,j)
    enddo
    enddo
!=======================================================================================================!
end subroutine dg3d_diff_flux
!=======================================================================================!
subroutine dg3d_diff_grads(elem,deriv,contrauvbuf,contrauv,couv,dif_gradu,dif_gradv)
!=======================================================================================!
    Implicit None
    integer, parameter :: south=1, east=2, north=3, west=4
    type(element_t), intent(in) :: elem
    type (derivative_t),                               intent(in)   :: deriv
    real (kind=real_kind), dimension(0:np+1,0:np+1,2), intent(in)   :: contrauvbuf
    real (kind=real_kind), dimension(np,np,2),         intent(in)   :: contrauv,couv
    real (kind=real_kind), dimension(np,np,2),         intent(out)  :: dif_gradu, dif_gradv

    real (kind=real_kind), dimension(np,4,2)   :: couv_halo, contrauv_halo
    real (kind=real_kind), dimension(np,np,2)  :: jflx
    real (kind=real_kind), dimension(np,np)    :: jfu,jfv,cu,cv
    Integer :: i,j,k, wall, eqn
!=======================================================================================================!
!  Boundary values for velocity & flux terms from the neighbors                                         !
!  H-field (scalar) on S, E, N, W sides of an element                                                   !
!  U,V-component on  (S,E,N,W) =  (1,2,3,4) sides of an element                                         !
!=======================================================================================================!
    do i = 1,np
    do j = 1, 2
       contrauv_halo(i,south,j)= contrauvbuf(i,   0,j)
       contrauv_halo(i,east,j) = contrauvbuf(np+1,i,j)
       contrauv_halo(i,north,j)= contrauvbuf(i,np+1,j)
       contrauv_halo(i,west,j) = contrauvbuf(0,i,   j)
    enddo
    enddo

    do k= 1,np
       couv_halo(k,1,1) = elem%met(1,1,k,1)  * contrauv_halo(k,1,1) + elem%met(1,2,k,1)  * contrauv_halo(k,1,2)
       couv_halo(k,1,2) = elem%met(2,1,k,1)  * contrauv_halo(k,1,1) + elem%met(2,2,k,1)  * contrauv_halo(k,1,2)
       couv_halo(k,2,1) = elem%met(1,1,np,k) * contrauv_halo(k,2,1) + elem%met(1,2,np,k) * contrauv_halo(k,2,2)
       couv_halo(k,2,2) = elem%met(2,1,np,k) * contrauv_halo(k,2,1) + elem%met(2,2,np,k) * contrauv_halo(k,2,2)

       couv_halo(k,3,1) = elem%met(1,1,k,np) * contrauv_halo(k,3,1) + elem%met(1,2,k,np) * contrauv_halo(k,3,2)
       couv_halo(k,3,2) = elem%met(2,1,k,np) * contrauv_halo(k,3,1) + elem%met(2,2,k,np) * contrauv_halo(k,3,2)
       couv_halo(k,4,1) = elem%met(1,1,1,k)  * contrauv_halo(k,4,1) + elem%met(1,2,1,k)  * contrauv_halo(k,4,2)
       couv_halo(k,4,2) = elem%met(2,1,1,k)  * contrauv_halo(k,4,1) + elem%met(2,2,1,k)  * contrauv_halo(k,4,2)
    enddo

    Call  jump_fluxint(deriv,couv,couv_halo,jflx)
        cu(:,:) = couv(:,:,1)
       jfu(:,:) = jflx(:,:,1)
    Call  ldg_grads(elem,deriv,cu,jfu,dif_gradu)

        cv(:,:) = couv(:,:,2)
       jfv(:,:) = jflx(:,:,2)
    Call  ldg_grads(elem,deriv,cv,jfv,dif_gradv)

!=======================================================================================================!
end subroutine dg3d_diff_grads
!=======================================================================================!

 subroutine dg3d_diff_grads_uv(elem,deriv,uvbuf,uv,dif_gradu,dif_gradv)
!=======================================================================================!
    Implicit None
    integer, parameter :: south=1, east=2, north=3, west=4
    type(element_t), intent(in) :: elem
    type (derivative_t),                               intent(in)   :: deriv
    real (kind=real_kind), dimension(0:np+1,0:np+1,2), intent(in)   :: uvbuf
    real (kind=real_kind), dimension(np,np,2),         intent(in)   :: uv
    real (kind=real_kind), dimension(np,np,2),         intent(out)  :: dif_gradu, dif_gradv

    real (kind=real_kind), dimension(np,np,2)  :: couv, gradu1, gradu2
    real (kind=real_kind), dimension(np,4,2)   :: couv_halo, uv_halo
    real (kind=real_kind), dimension(np,np,2)  :: jflx
    real (kind=real_kind), dimension(np,np)    :: jfu,jfv,cu,cv
    Integer :: i,j,k, wall, eqn
!=======================================================================================================!
!  Boundary values for velocity & flux terms from the neighbors
!  H-field (scalar) on S, E, N, W sides of an element
!  U,V-component on  (S,E,N,W) =  (1,2,3,4) sides of an element
!=======================================================================================================!
    do i = 1,np
    do j = 1, 2
       uv_halo(i,south,j)= uvbuf(i,   0,j)
       uv_halo(i,east,j) = uvbuf(np+1,i,j)
       uv_halo(i,north,j)= uvbuf(i,np+1,j)
       uv_halo(i,west,j) = uvbuf(0,i,   j)
    enddo
    enddo

   do k= 1, np
       couv_halo(k,1,1) = elem%D(1,1,k,1)  * uv_halo(k,1,1) + elem%D(2,1,k,1)  * uv_halo(k,1,2)
       couv_halo(k,1,2) = elem%D(1,2,k,1)  * uv_halo(k,1,1) + elem%D(2,2,k,1)  * uv_halo(k,1,2)
       couv_halo(k,2,1) = elem%D(1,1,np,k) * uv_halo(k,2,1) + elem%D(2,1,np,k) * uv_halo(k,2,2)
       couv_halo(k,2,2) = elem%D(1,2,np,k) * uv_halo(k,2,1) + elem%D(2,2,np,k) * uv_halo(k,2,2)

       couv_halo(k,3,1) = elem%D(1,1,k,np) * uv_halo(k,3,1) + elem%D(2,1,k,np) * uv_halo(k,3,2)
       couv_halo(k,3,2) = elem%D(1,2,k,np) * uv_halo(k,3,1) + elem%D(2,2,k,np) * uv_halo(k,3,2)
       couv_halo(k,4,1) = elem%D(1,1,1,k)  * uv_halo(k,4,1) + elem%D(2,1,1,k)  * uv_halo(k,4,2)
       couv_halo(k,4,2) = elem%D(1,2,1,k)  * uv_halo(k,4,1) + elem%D(2,2,1,k)  * uv_halo(k,4,2)
    enddo

        couv(:,:,:) = sphere2cov(uv,elem%D)

    Call  jump_fluxint(deriv,couv,couv_halo,jflx)
        cu(:,:) = couv(:,:,1)
       jfu(:,:) = jflx(:,:,1)
    Call  ldg_grads(elem,deriv,cu,jfu,dif_gradu)

        cv(:,:) = couv(:,:,2)
       jfv(:,:) = jflx(:,:,2)
    Call  ldg_grads(elem,deriv,cv,jfv,dif_gradv)

!=======================================================================================================!
end subroutine dg3d_diff_grads_uv 
!=======================================================================================================!

!=======================================================================================================!
subroutine central_fluxint(elem,deriv,grad,grad_senw,cfluxint)

   integer, parameter :: south=1, east=2, north=3, west=4
   type(element_t), intent(in) :: elem
   type (derivative_t)                                  :: deriv
   real (kind=real_kind), dimension(np,np,2),intent(in) :: grad
   real (kind=real_kind), dimension(np,4,2), intent(in) :: grad_senw

   real (kind=real_kind), dimension(np,np), intent(out) :: cfluxint

   real (kind=real_kind), dimension(np,np) :: mij
   real (kind=real_kind), dimension(np)   :: cflx_south,cflx_north,cflx_east,cflx_west
   real(kind=real_kind) ::   f_left, f_right, s1,s2, x0y0

   integer i,j,k, west_edge, south_edge

      mij(:,:) = 0.0D0
      mij(1,1) = 1.0D0
    mij(np,np) = 1.0D0

         x0y0 = -0.25D0 * dd_pi
         west_edge = 0
         south_edge = 0

        !do j = 1,np
        ! if (elem%cartp(1,j)%x == x0y0) then
        !   west_edge = 1
        ! endif
        !enddo
        !do i = 1,np
        ! if (elem%cartp(i,1)%y  == x0y0) then
        !   south_edge = 1
        ! endif
        !enddo

    !! West/South edge of the Cubed-Sphere detection

          if (elem%cartp(1,1)%x == x0y0) west_edge=1
          if (elem%cartp(1,1)%y == x0y0) south_edge=1

    !For "Centered/sided" Flux for the grad  vectors

           ! East & West   Central flux

         do j = 1,np
                  f_left  = grad_senw(j,west,1)
                  f_right = grad(1,j,1)
              if (west_edge == 1) f_left = f_right
               !cflx_west(j) = 0.5D0 *(f_right + f_left)
               cflx_west(j) = f_left

                  f_left  = grad(np,j,1)
                  f_right = grad_senw(j,east,1)
               !cflx_east(j) = 0.5D0 *(f_right + f_left)
              cflx_east(j) = f_left
          end do

           ! North& South  Central flux

         do i = 1,np
                   f_left = grad_senw(i,south,2)
                  f_right = grad(i,1,2)
              if (south_edge == 1) f_left = f_right
              !cflx_south(i) = 0.5D0 *(f_right + f_left)
              cflx_south(i) = f_left

                  f_left  = grad(i,np,2)
                  f_right = grad_senw(i,north,2)
              !cflx_north(i) = 0.5D0 *(f_right + f_left)
               cflx_north(i) = f_left
         end do

        !Flux integral along the element boundary

        do j = 1,np
        do i = 1,np
            s1 = (cflx_east(j) *mij(i,np) - cflx_west(j) *mij(i,1) )* deriv%mvv_twt(j,j)
            s2 = (cflx_north(i)*mij(j,np) - cflx_south(i)*mij(j,1) )* deriv%mvv_twt(i,i)
            cfluxint(i,j) = (s1 + s2) * rrearth
        end do
        end do

 end subroutine  central_fluxint
!=======================================================================================================!
!=======================================================================================================!
subroutine jump_fluxint(deriv,uv,uv_senw,jfluxint)

   integer, parameter :: south=1, east=2, north=3, west=4
   type (derivative_t)                                  :: deriv
   real (kind=real_kind), dimension(np,np,2),intent(in) :: uv
   real (kind=real_kind), dimension(np,4,2), intent(in) :: uv_senw

   real (kind=real_kind), dimension(np,np,2), intent(out) :: jfluxint

   real (kind=real_kind), dimension(np,np) :: mij
   real (kind=real_kind), dimension(np)   :: jump_south,jump_north,jump_east,jump_west
   real(kind=real_kind) ::   f_left, f_right, s1,s2

   integer i,j,k,eqn

      mij(:,:) = 0.0D0
      mij(1,1) = 1.0D0
    mij(np,np) = 1.0D0


    !For Jump Flux for the covariant vectors

      do eqn = 1, 2

           ! East & West   LF flux  (fjmax <- max of flux Jacobian)

         do j = 1,np
                  f_left  = uv_senw(j,west,eqn)
                  f_right = uv(1,j,eqn)
               jump_west(j) = 1.0D0 *(f_right - f_left)

                  f_left  = uv(np,j,eqn)
                  f_right = uv_senw(j,east,eqn)
               jump_east(j) = 1.0D0 *(f_right - f_left)
          end do

           ! North& South  LF flux  (fjmax <- max of flux Jacobian)

         do i = 1,np
                   f_left = uv_senw(i,south,eqn)
                  f_right = uv(i,1,eqn)
              jump_south(i) = 1.0D0 *(f_right - f_left)

                  f_left  = uv(i,np,eqn)
                  f_right = uv_senw(i,north,eqn)
              jump_north(i) = 1.0D0 *(f_right - f_left)
         end do

        !Flux integral along the element boundary
        do j = 1,np
        do i = 1,np
            s1 = (jump_east(j) *mij(i,np) - jump_west(j) *mij(i,1) )* deriv%mvv_twt(j,j)
            s2 = (jump_north(i)*mij(j,np) - jump_south(i)*mij(j,1) )* deriv%mvv_twt(i,i)
            jfluxint(i,j,eqn) = (s1 + s2) * rrearth
        end do
        end do

    end do

 end subroutine  jump_fluxint
!=======================================================================================================!
  subroutine ldg_grads(elem,deriv,uv,jfint,grad_u)

    type (derivative_t)   :: deriv
    type(element_t), intent(in) :: elem
    real(kind=real_kind), dimension(np,np), intent(in) :: uv, jfint
    real(kind=real_kind), dimension(np,np,2), intent(out) :: grad_u

    real(kind=real_kind), dimension(np,np) :: gr1 ,gr2, sg
    real(kind=real_kind), dimension(2,2,np,np) :: ginv

    real(kind=real_kind) ::  s1,s2, weight , v1,v2
    integer :: i,j,l
!=======================================================================================================!

    sg(:,:) = elem%metdet(:,:)
    ginv(:,:,:,:) = elem%metinv(:,:,:,:)

!!  DoubleInt[ grad(U) = (U_x1, U_x2)]

    do j=1,np
       do l=1,np
          s2=0.0D0
          s1=0.0D0
          do i=1,np
             s1 = s1 + deriv%Dvv(i,l)* uv(i,j)
             s2 = s2 + deriv%Dvv(i,l)* uv(j,i)
          end do
          gr1(l,j) = s1
          gr2(j,l) = s2
       end do
    end do

    do j=1,np
       do i=1,np
          weight = elem%mp(i,j)*rrearth      !(double integral effect)
          gr1(i,j) = gr1(i,j)*weight
          gr2(i,j) = gr2(i,j)*weight
       end do
    end do

!! Grad recovery by LDG way

    do j=1,np
       do i=1,np
          weight = 1.0D0/elem%mp(i,j)
          gr1(i,j) =  (jfint(i,j) + gr1(i,j))* weight
          gr2(i,j) =  (jfint(i,j) + gr2(i,j))* weight
       end do
    end do

 !! Tensor Gradients (for internal terms in the general Laplacian)

       do j=1,np
            do i=1,np
                   v1 = gr1(i,j)
                   v2 = gr2(i,j)
             grad_u(i,j,1) =  sg(i,j) *(ginv(1,1,i,j)*v1 + ginv(1,2,i,j)*v2)
             grad_u(i,j,2) =  sg(i,j) *(ginv(2,1,i,j)*v1 + ginv(2,2,i,j)*v2)
            end do
        end do

end subroutine ldg_grads
!=======================================================================================================!
!=======================================================================================================!
end module dg3d_dynamics_mod    


