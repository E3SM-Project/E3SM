#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!=======================================================================================!
module dg_core_mod
!=======================================================================================!
 use kinds,           only : real_kind
 use dimensions_mod,  only : np, nlev ,ne 
 use element_mod, only : element_t 
 use quadrature_mod, only: quadrature_t, gauss, gausslobatto
 use derivative_mod,  only : derivative_t, gradient_wk
 use physical_constants, only : rrearth, rearth, g, dd_pi
 use control_mod, only : nu  
!=======================================================================================!
#ifdef _PRIMDG
!=======================================================================================!
implicit none
 private
!=======================================================================================!
 public :: swsys_flux
 public :: co2contra
 public :: contra2co
 public :: sphere2contra
! public :: uv_vorticity
 public :: sphere2cov    
 public :: contra2sphere
 public :: height2phi
 public :: psi2height    
 public :: phi2height    
 public :: dp2pt
 public :: pt2dp 
 public :: mono_filter
 public :: indexer
!=======================================================================================!
#endif
!=======================================================================================!
#ifdef _SWDG
!=======================================================================================! 
 use dg_flux_mod, only : riemanntype,fnum_flux,signalwave
!=======================================================================================!
 implicit none
 private
!=======================================================================================! 
 public :: swsys_flux
 public :: dg_adv_model 
 public :: dg_sw_model     
 public :: gradient_mass   
 public :: gradient_mom    
 public :: source_term
 public :: adv_flux_term
 public :: sw_fjmax 
 public :: dgsw_uvh_rhs 
!=======================================================================================!
 public :: co2contra
 public :: contra2co
 public :: sphere2contra
 public :: uv_vorticity
 public :: sphere2cov    
 public :: cov2sphere    
 public :: contra2sphere
 public :: height2phi
 public :: phi2height 
 public :: mono_filter
 public :: indexer 
!============================================= LDG Specific ===========================================!
 public :: dg_diff_grads_uv 
 public :: dg_diff_grads
 public :: dg_diff_flux
 private :: jump_fluxint
 private :: central_fluxint
 private :: ldg_grads
!=======================================================================================!    
#endif
!=======================================================================================!
 real (kind=real_kind), dimension(:,:,:,:), pointer :: met,metinv,D,Dinv
 real (kind=real_kind), dimension(:,:),     pointer :: fcor,rmv,mv,mvi,metdet,rmetdet
 real (kind=real_kind), public, dimension(np,np) :: mmx,mmxi
!=======================================================================================!
 real (kind=real_kind), parameter, private:: scale= 1.0D0
 real (kind=real_kind), private:: zero= 0.0D0
 real (kind=real_kind), private:: one = 1.0D0        
!=======================================================================================!
 contains
!=======================================================================================!
#ifdef _SWDG
!=======================================================================================!
!  DG Shallow-Water Model (u,v) formulation RHS					   			!

    subroutine dgsw_uvh_rhs(elem,deriv,uvbuf,htbuf,uv,ht,hill,rhs)

!=======================================================================================!
  use element_mod,     only : element_t
    Implicit None
    integer, parameter :: south=1, east=2, north=3, west=4
    type(element_t) ,                        intent(inout), target :: elem
    type (derivative_t),                               intent(in)  :: deriv 
    real (kind=real_kind), dimension(0:np+1,0:np+1),   intent(in)  :: htbuf
    real (kind=real_kind), dimension(0:np+1,0:np+1,2), intent(in)  :: uvbuf
    real (kind=real_kind), dimension(np,np,2),         intent(in)  :: uv 
    real (kind=real_kind), dimension(np,np),           intent(in)  :: ht,hill      
    real (kind=real_kind), dimension(np,np,3),         intent(out) :: rhs

    real (kind=real_kind), dimension(4,np)    :: SL,SR,lambda

    real (kind=real_kind), dimension(np,np,2) :: contrauv,couv
    real (kind=real_kind), dimension(np,np)   :: phi,gcori, mmtx, gradu1,gradu2, u1_rhs, u2_rhs 
    real (kind=real_kind), dimension(np,np)   :: energy, gh, vort, gh11, gh22 
    real (kind=real_kind), dimension(np,np,2) :: ghij, fxy
    real (kind=real_kind), dimension(np,4,2)  :: fxy_halo,uv_halo,contuv_halo, couv_halo
    real (kind=real_kind), dimension(np,4)    :: fx_halo,fy_halo,phi_halo
    real (kind=real_kind), dimension(np,4)    :: gh11_halo, gh22_halo, sg_halo       
    real (kind=real_kind), dimension(np,4)    :: gh_halo,ht_halo,energy_halo
    real (kind=real_kind), dimension(np,np,3) :: grad,numflux,source
    real (kind=real_kind), dimension(np,np,3) :: vec 
    real (kind=real_kind), dimension(np,4,3)  :: vec_halo
    real (kind=real_kind), dimension(np,np,2,3):: fluxvec 
    real (kind=real_kind), dimension(np,4,2,3) :: flux_halo 
    real (kind=real_kind), dimension(4)       :: fjmax

    real (kind=real_kind) :: alfa1,alfa2,ul,ur,fact,grv
    real (kind=real_kind) :: rtmp,tmp1,tmp2,tmp3,tmp4
    Integer :: i,j,k, wall, eqn 
    real*8  :: st,et,time_dgsw,time_swflux
!=======================================================================================================!
!=======================================================================================================!
!  From HOMME to Local Pointer										!
!=======================================================================================================!
    met    => elem%met
    metinv => elem%metinv
    metdet => elem%metdet
    mv     => elem%mp
    fcor   => elem%fcor
    Dinv   => elem%Dinv
    D      => elem%D
!=======================================================================================================!
!  Length-scale redefined, gravity									!
!=======================================================================================================!
    grv = g 
!=======================================================================================================!
!  Mass and Inverse-Mass Matrix										!
!=======================================================================================================!
#if 1
    mmx(:,:) = mv(:,:)
#if defined(_USE_VECTOR)
    call vrec(mmxi(1,1),mmx(1,1),np*np)
#else
    do j = 1, np
    do i = 1, np
       mmxi(i,j) = 1.0D0 / mmx(i,j) 
    enddo
    enddo
#endif
#endif

!=======================================================================================================!
!  Boundary values for velocity & flux terms from the neighbors						!
!  H-field (scalar) on S, E, N, W sides of an element							!
!  U,V-component on  (S,E,N,W) =  (1,2,3,4) sides of an element						!
!=======================================================================================================!
    do i = 1, np
       ht_halo(i,south) = htbuf(i,   0)
       ht_halo(i,east)  = htbuf(np+1,i)
       ht_halo(i,north) = htbuf(i,np+1)
       ht_halo(i,west)  = htbuf(0,i   )
    do j = 1, 2 
       uv_halo(i,south,j)= uvbuf(i,   0,j)
       uv_halo(i,east,j) = uvbuf(np+1,i,j)
       uv_halo(i,north,j)= uvbuf(i,np+1,j)
       uv_halo(i,west,j) = uvbuf(0,i,   j)
    enddo       
    enddo

!=======================================================================================================!
!  Covariant components for the halo region [u^1,u^2] = A^-1 * [u,v] 
!=======================================================================================================!
    do k= 1, np
       contuv_halo(k,1,1) = Dinv(1,1,k,1)  * uv_halo(k,1,1) + Dinv(1,2,k,1)  * uv_halo(k,1,2)
       contuv_halo(k,1,2) = Dinv(2,1,k,1)  * uv_halo(k,1,1) + Dinv(2,2,k,1)  * uv_halo(k,1,2)
       contuv_halo(k,2,1) = Dinv(1,1,np,k) * uv_halo(k,2,1) + Dinv(1,2,np,k) * uv_halo(k,2,2)
       contuv_halo(k,2,2) = Dinv(2,1,np,k) * uv_halo(k,2,1) + Dinv(2,2,np,k) * uv_halo(k,2,2)

       contuv_halo(k,3,1) = Dinv(1,1,k,np) * uv_halo(k,3,1) + Dinv(1,2,k,np) * uv_halo(k,3,2)
       contuv_halo(k,3,2) = Dinv(2,1,k,np) * uv_halo(k,3,1) + Dinv(2,2,k,np) * uv_halo(k,3,2)
       contuv_halo(k,4,1) = Dinv(1,1,1,k)  * uv_halo(k,4,1) + Dinv(1,2,1,k)  * uv_halo(k,4,2)
       contuv_halo(k,4,2) = Dinv(2,1,1,k)  * uv_halo(k,4,1) + Dinv(2,2,1,k)  * uv_halo(k,4,2)
    enddo

    contrauv(:,:,:) = sphere2contra(uv,Dinv) 
        couv(:,:,:) = sphere2cov(uv,D) 

    do j= 1, np
    do i= 1, np
       gh(i,j) = grv*ht(i,j)
    enddo
    enddo

    do wall= 1, 4
    do k= 1, np
       gh_halo(k,wall) = grv*ht_halo(k,wall)
    enddo
    enddo


    do j= 1, np
    do i= 1, np
       gh11(i,j) = gh(i,j)*metinv(1,1,i,j)
       gh22(i,j) = gh(i,j)*metinv(2,2,i,j)
    enddo
    enddo

    do k= 1, np
       tmp1 = grv*ht_halo(k,1)
       tmp2 = grv*ht_halo(k,2)
       tmp3 = grv*ht_halo(k,3)
       tmp4 = grv*ht_halo(k,4)

       gh11_halo(k,1) = tmp1*metinv(1,1,k,1)
       gh11_halo(k,2) = tmp2*metinv(1,1,np,k)
       gh11_halo(k,3) = tmp3*metinv(1,1,k,np)
       gh11_halo(k,4) = tmp4*metinv(1,1,1,k)

       gh22_halo(k,1) = tmp1*metinv(2,2,k,1)
       gh22_halo(k,2) = tmp2*metinv(2,2,np,k)
       gh22_halo(k,3) = tmp3*metinv(2,2,k,np)
       gh22_halo(k,4) = tmp4*metinv(2,2,1,k)
    enddo
!=======================================================================================================!
!  Covariant components for the halo region
!=======================================================================================================!
    do k= 1, np
       couv_halo(k,1,1) = met(1,1,k,1)  * contuv_halo(k,1,1) + met(1,2,k,1)  * contuv_halo(k,1,2)
       couv_halo(k,1,2) = met(2,1,k,1)  * contuv_halo(k,1,1) + met(2,2,k,1)  * contuv_halo(k,1,2)
       couv_halo(k,2,1) = met(1,1,np,k) * contuv_halo(k,2,1) + met(1,2,np,k) * contuv_halo(k,2,2)
       couv_halo(k,2,2) = met(2,1,np,k) * contuv_halo(k,2,1) + met(2,2,np,k) * contuv_halo(k,2,2)

       couv_halo(k,3,1) = met(1,1,k,np) * contuv_halo(k,3,1) + met(1,2,k,np) * contuv_halo(k,3,2)
       couv_halo(k,3,2) = met(2,1,k,np) * contuv_halo(k,3,1) + met(2,2,k,np) * contuv_halo(k,3,2)
       couv_halo(k,4,1) = met(1,1,1,k)  * contuv_halo(k,4,1) + met(1,2,1,k)  * contuv_halo(k,4,2)
       couv_halo(k,4,2) = met(2,1,1,k)  * contuv_halo(k,4,1) + met(2,2,1,k)  * contuv_halo(k,4,2)
    enddo

!=======================================================================================================!
!  Flux computations for the Continuity equation 							!
!   (u,v) Fluxes for the continuity equation								!
!   Note   "ht - hill => " Depth of the fluid								!					
!=======================================================================================================!
    do j = 1, np
    do i = 1, np
       phi(i,j)= (ht(i,j) - hill(i,j)) * metdet(i,j)
       fxy(i,j,1)= contrauv(i,j,1) * phi(i,j)
       fxy(i,j,2)= contrauv(i,j,2) * phi(i,j)
    enddo
    enddo
!=======================================================================================================!
!  (u1, u2) Fluxes along the Halo  (S,E,N,W) from the neighbours					!
!=======================================================================================================!
    do k = 1, np
       phi_halo(k,south) = (ht_halo(k,south) - hill(k, 1))* metdet(k,1)
       phi_halo(k,east)  = (ht_halo(k,east) - hill(np,k)) * metdet(np,k)
       phi_halo(k,north) = (ht_halo(k,north) - hill(k,np))* metdet(k,np)
       phi_halo(k,west)  = (ht_halo(k,west) - hill(1, k)) * metdet(1,k)
    do j = 1,2
    do wall=1,4 
       fxy_halo(k,wall,j)= contuv_halo(k,wall,j) * phi_halo(k,wall)
    enddo
    enddo
    enddo
!=======================================================================================================!
! Energy for the internal elements => energy(i,j)= g*h + 0.5*(u_1*u^1 + u_2*u^2)			!
!=======================================================================================================!
    do j= 1, np
    do i= 1, np
       energy(i,j)= gh(i,j) + 0.5D0*( uv(i,j,1)*uv(i,j,1) + uv(i,j,2)*uv(i,j,2) )
    enddo
    enddo
!=======================================================================================================!
! Energy for the Halo region from the neighbours							!
!=======================================================================================================!
    do wall= 1, 4
    do k= 1, np
       energy_halo(k,wall) = gh_halo(k,wall) + 0.5D0*(uv_halo(k,wall,1)*uv_halo(k,wall,1)   & 
                                                    + uv_halo(k,wall,2)*uv_halo(k,wall,2) )
    enddo
    enddo

 ! Source terms 

    vort(:,:) = uv_vorticity(uv,metdet,D,deriv)

    do j = 1, np
    do i = 1, np
     source(i,j,1) = 0.0D0
     source(i,j,2) = (vort(i,j) + fcor(i,j) )* uv(i,j,2) !* mmx(i,j) 
     source(i,j,3) =-(vort(i,j) + fcor(i,j) )* uv(i,j,1) !* mmx(i,j) 
    enddo
    enddo

  ! Max flux-Jacobian 

    fjmax(:) =  sw_fjmax(contrauv,contuv_halo,gh11,gh22,gh11_halo,gh22_halo)

!=======================================================================================================!
    do j = 1, np
    do i = 1, np
       vec(i,j,1) = phi(i,j)
       vec(i,j,2) = couv(i,j,1)
       vec(i,j,3) = couv(i,j,2)

       fluxvec(i,j,1,1) = fxy(i,j,1)
       fluxvec(i,j,2,1) = fxy(i,j,2)

       fluxvec(i,j,1,2) = energy(i,j)
       fluxvec(i,j,2,2) = 0.0D0     

       fluxvec(i,j,1,3) = 0.0D0     
       fluxvec(i,j,2,3) = energy(i,j)
    enddo
    enddo

    vec_halo(:,:,1) = phi_halo(:,:)
    vec_halo(:,:,2) = couv_halo(:,:,1)
    vec_halo(:,:,3) = couv_halo(:,:,2)

    flux_halo(:,:,:,1) = fxy_halo(:,:,:)
    flux_halo(:,:,1,2) = energy_halo(:,:)
    flux_halo(:,:,2,2) = 0.0D0         
    flux_halo(:,:,1,3) = 0.0D0         
    flux_halo(:,:,2,3) = energy_halo(:,:)

!==
    Call swsys_flux(3,deriv,fjmax,vec,vec_halo,fluxvec,flux_halo,numflux)
!==
!--------------------------------------------------
!  Gradient  for the continuity equation     	
!----------------------------------------------
    grad(:,:,1) = gradient_mass(deriv,fxy)

!-------------------------------------------------
!  Gradient  for the momentum equations									!	
!------------------------------------------------
#if 0
   grad(:,:,2:3)=gradient_wk(energy,deriv)
#else
   Call gradient_mom(deriv,energy,gradu1,gradu2)

   !! converting to lat/lon spherical gradients [u,v] = A^-T * [u_1,u_2]  

!  do j=1,np
!  do i=1,np
!     grad(i,j,2) = gradu1(i,j) *Dinv(1,1,i,j) + gradu2(i,j) * Dinv(2,1,i,j) 
!     grad(i,j,3) = gradu1(i,j) *Dinv(1,2,i,j) + gradu2(i,j) * Dinv(2,2,i,j) 
!  enddo
!  enddo
#endif
!==
!  Compute RHS of the ODE system corresponding  to the SW model         				!
!  mvi(i,j)= 1 / mv(i,j)										!
!==
    do j = 1, np
    do i = 1, np
       rtmp = 1.0D0/elem%mp(i,j)
       rhs(i,j,1)= (source(i,j,1) + grad(i,j,1) - numflux(i,j,1)) * rtmp
       u1_rhs(i,j)= (gradu1(i,j) - numflux(i,j,2)) * rtmp
       u2_rhs(i,j)= (gradu2(i,j) - numflux(i,j,3)) * rtmp
    enddo
    enddo

   do j=1,np
   do i=1,np
      rhs(i,j,2) = u1_rhs(i,j) *Dinv(1,1,i,j) + u2_rhs(i,j) * Dinv(2,1,i,j) + source(i,j,2) 
      rhs(i,j,3) = u1_rhs(i,j) *Dinv(1,2,i,j) + u2_rhs(i,j) * Dinv(2,2,i,j) + source(i,j,3)  
   enddo
   enddo
 end subroutine  dgsw_uvh_rhs 
!==================================================================================!
! Element-wise Max flux Jacobian for SW system 
!----------------------------------------------------------------------------------
 Function sw_fjmax(contuv,contuv_halo,gh11,gh22,gh11_halo,gh22_halo) result(fjmax)
!----------------------------------------------------------------------------------
 Implicit None
 real (kind=real_kind),dimension(np,np,2),intent(in):: contuv
 real (kind=real_kind), dimension(np,4,2),intent(in):: contuv_halo
 real (kind=real_kind),dimension(np,np),  intent(in):: gh11,gh22 
 real (kind=real_kind),dimension(np,4),   intent(in):: gh11_halo, gh22_halo 

 real (kind=real_kind),dimension(4)   :: fjmax
 real (kind=real_kind):: alfa1,alfa2, ul,ur 
 integer, parameter:: south=1,east=2,north=3,west=4
 integer:: i,j,wall
!========================================================
    alfa1 = 0.0D0
    alfa2 = 0.0D0
    do i = 1, np
       ul = abs(contuv(i,1,2))      + sqrt(gh22(i,1))
       ur = abs(contuv_halo(i,south,2)) + sqrt(gh22_halo(i,south))
       alfa1= max(alfa1,ul,ur)
       
       ul = abs(contuv_halo(i,north,2)) + sqrt(gh22_halo(i,north))
       ur = abs(contuv(i,np,2))     + sqrt(gh22(i,np))
       alfa2= max(alfa2,ul,ur)
    enddo

    fjmax(south) = alfa1
    fjmax(north) = alfa2

    alfa1 = 0.0D0
    alfa2 = 0.0D0
    do j = 1, np
       ul = abs(contuv_halo(j,west,1))  + sqrt(gh11_halo(j,west))
       ur = abs(contuv(1,j,1))      + sqrt(gh11(1,j))
     alfa1= max(alfa1,ul,ur)
       ul = abs(contuv(np,j,1))     + sqrt(gh11(np,j))
       ur = abs(contuv_halo(j,east,1))  + sqrt(gh11_halo(j,east))
     alfa2= max(alfa2,ul,ur)
    enddo

    fjmax(west) = alfa1
    fjmax(east) = alfa2 

   End Function sw_fjmax  
!=======================================================================================!

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!  DG Shallow-Water Model RHS (Classic covariant formulation)					   			!
!=======================================================================================!
subroutine dg_sw_model(elem,deriv,contrauvbuf,htbuf,contrauv,couv,ht,hill,rhs)
!=======================================================================================!
  use element_mod,     only : element_t
    Implicit None
    integer, parameter :: south=1, east=2, north=3, west=4
    type(element_t) ,                          intent(inout), target :: elem
    type (derivative_t),                               intent(in)   :: deriv 
    real (kind=real_kind), dimension(0:np+1,0:np+1),   intent(in)   :: htbuf
    real (kind=real_kind), dimension(0:np+1,0:np+1,2), intent(in)   :: contrauvbuf
    real (kind=real_kind), dimension(np,np,2),         intent(in)   :: contrauv,couv
    real (kind=real_kind), dimension(np,np),           intent(in)   :: ht,hill      
    real (kind=real_kind), dimension(np,np,3),         intent(out)  :: rhs

    real (kind=real_kind), dimension(4)       :: fjmax
    real (kind=real_kind), dimension(4,np)    :: SL,SR,lambda                
    real (kind=real_kind), dimension(np,np)   :: phi,gcori, mmtx, gradu1,gradu2
    real (kind=real_kind), dimension(np,np)   :: energy, gh11,gh22
    real (kind=real_kind), dimension(np,np,2) :: ghij, fxy
    real (kind=real_kind), dimension(np,4,2)  :: fxy_halo,contrauv_halo,couv_halo
    real (kind=real_kind), dimension(np,4)    :: fx_halo,fy_halo,phi_halo
    real (kind=real_kind), dimension(np,4)    :: gh11_halo,gh22_halo,ht_halo,energy_halo
    real (kind=real_kind), dimension(np,np,3)  :: grad,numflux,source
    real (kind=real_kind), dimension(np,np,3)  :: vec 
    real (kind=real_kind), dimension(np,4,3)   :: vec_halo
    real (kind=real_kind), dimension(np,np,2,3):: fluxvec 
    real (kind=real_kind), dimension(np,4,2,3) :: flux_halo 
    real (kind=real_kind) :: alfa1,alfa2,ul,ur,fact,grv
    real (kind=real_kind) :: rtmp,tmp1,tmp2,tmp3,tmp4
    Integer :: i,j,k, wall, eqn 
    real*8  :: st,et,time_dgsw,time_swflux
!=======================================================================================================!
!=======================================================================================================!
!  From HOMME to Local Pointer										!
!=======================================================================================================!
    met    => elem%met
    metinv => elem%metinv
    metdet => elem%metdet
    mv     => elem%mp
    fcor   => elem%fcor
!=======================================================================================================!
!  Length-scale redefined, gravity									!
!=======================================================================================================!
    grv = g 
!=======================================================================================================!
!  Mass and Inverse-Mass Matrix										!
!=======================================================================================================!
#if 1
    mmx(:,:) = mv(:,:)
#if defined(_USE_VECTOR)
    call vrec(mmxi(1,1),mmx(1,1),np*np)
#else
    do j = 1, np
    do i = 1, np
       mmxi(i,j) = one / mmx(i,j) 
    enddo
    enddo
#endif
#endif
!=======================================================================================================!
! Source term	      =>  source(i,j,1)= 0								!
!			  source(i,j,2)= u^2(i,j)*(vor(i,j)+gcori(i,j))*M(i,j)				!
!			  source(i,j,3)=-u^1(i,j)*(vor(i,j)+gcori(i,j))*M(i,j)				!
!=======================================================================================================! 
! Coriolis term	=>  gcori_{ij}= sqrt(det(G_ij))*fcor_{ij}
    do j = 1, np
    do i = 1, np
       gcori(i,j)= fcor(i,j) * metdet(i,j)
    enddo
    enddo 

    Call source_term(deriv,elem%mp,gcori,contrauv,couv,source)
!=======================================================================================================!
!  Boundary values for velocity & flux terms from the neighbors						!
!  H-field (scalar) on S, E, N, W sides of an element							!
!  U,V-component on  (S,E,N,W) =  (1,2,3,4) sides of an element						!
!=======================================================================================================!
    do i = 1, np
       ht_halo(i,south) = htbuf(i,   0)
       ht_halo(i,east)  = htbuf(np+1,i)
       ht_halo(i,north) = htbuf(i,np+1)
       ht_halo(i,west)  = htbuf(0,i   )
    do j = 1, 2 
       contrauv_halo(i,south,j)= contrauvbuf(i,   0,j)
       contrauv_halo(i,east,j) = contrauvbuf(np+1,i,j)
       contrauv_halo(i,north,j)= contrauvbuf(i,np+1,j)
       contrauv_halo(i,west,j) = contrauvbuf(0,i,   j)
    enddo       
    enddo
!=======================================================================================================!
!  Flux computations for the Continuity equation 							!
!   (u,v) Fluxes for the continuity equation								!
!   Note   "ht - hill => " Depth of the fluid								!					
!=======================================================================================================!
    do j = 1, np
    do i = 1, np
       phi(i,j)= (ht(i,j) - hill(i,j)) * metdet(i,j)
       fxy(i,j,1)= contrauv(i,j,1) * phi(i,j)
       fxy(i,j,2)= contrauv(i,j,2) * phi(i,j)
    enddo
    enddo
!=======================================================================================================!
!  (u1, u2) Fluxes along the Halo  (S,E,N,W) from the neighbours					!
!=======================================================================================================!
    do k = 1, np
       phi_halo(k,south) = (ht_halo(k,south) - hill(k, 1))* metdet(k,1)
       phi_halo(k,east)  = (ht_halo(k,east) - hill(np,k)) * metdet(np,k)
       phi_halo(k,north) = (ht_halo(k,north) - hill(k,np))* metdet(k,np)
       phi_halo(k,west)  = (ht_halo(k,west) - hill(1, k)) * metdet(1,k)
    do j = 1,2
    do wall=1,4 
       fxy_halo(k,wall,j)= contrauv_halo(k,wall,j) * phi_halo(k,wall)
    enddo
    enddo
    enddo
!=======================================================================================================!
!  Signal Speed: Max of Flux Jacobians on each wall							!
!=======================================================================================================!
    do j= 1, np
    do i= 1, np
       gh11(i,j) = grv*ht(i,j)*metinv(1,1,i,j)
       gh22(i,j) = grv*ht(i,j)*metinv(2,2,i,j)
    enddo
    enddo
   
    do k= 1, np
       tmp1 = grv*ht_halo(k,1)
       tmp2 = grv*ht_halo(k,2)
       tmp3 = grv*ht_halo(k,3)
       tmp4 = grv*ht_halo(k,4)

       gh11_halo(k,1) = tmp1*metinv(1,1,k,1)
       gh11_halo(k,2) = tmp2*metinv(1,1,np,k)
       gh11_halo(k,3) = tmp3*metinv(1,1,k,np)
       gh11_halo(k,4) = tmp4*metinv(1,1,1,k)

       gh22_halo(k,1) = tmp1*metinv(2,2,k,1)
       gh22_halo(k,2) = tmp2*metinv(2,2,np,k)
       gh22_halo(k,3) = tmp3*metinv(2,2,k,np)
       gh22_halo(k,4) = tmp4*metinv(2,2,1,k)
    enddo
!=======================================================================================================!
!  Covariant components for the halo region 								!
!=======================================================================================================!
    do k= 1, np
       couv_halo(k,1,1) = met(1,1,k,1)  * contrauv_halo(k,1,1) + met(1,2,k,1)  * contrauv_halo(k,1,2)
       couv_halo(k,1,2) = met(2,1,k,1)  * contrauv_halo(k,1,1) + met(2,2,k,1)  * contrauv_halo(k,1,2)
       couv_halo(k,2,1) = met(1,1,np,k) * contrauv_halo(k,2,1) + met(1,2,np,k) * contrauv_halo(k,2,2)
       couv_halo(k,2,2) = met(2,1,np,k) * contrauv_halo(k,2,1) + met(2,2,np,k) * contrauv_halo(k,2,2)

       couv_halo(k,3,1) = met(1,1,k,np) * contrauv_halo(k,3,1) + met(1,2,k,np) * contrauv_halo(k,3,2)
       couv_halo(k,3,2) = met(2,1,k,np) * contrauv_halo(k,3,1) + met(2,2,k,np) * contrauv_halo(k,3,2)
       couv_halo(k,4,1) = met(1,1,1,k)  * contrauv_halo(k,4,1) + met(1,2,1,k)  * contrauv_halo(k,4,2)
       couv_halo(k,4,2) = met(2,1,1,k)  * contrauv_halo(k,4,1) + met(2,2,1,k)  * contrauv_halo(k,4,2)
    enddo
!=======================================================================================================!
! Energy for the internal elements => energy(i,j)= g*h + 0.5*(u_1*u^1 + u_2*u^2)			!
!=======================================================================================================!
    do j= 1, np
    do i= 1, np
       energy(i,j)= grv*ht(i,j) + 0.5D0*( couv(i,j,1)*contrauv(i,j,1) + couv(i,j,2)*contrauv(i,j,2) )
    enddo
    enddo
!=======================================================================================================!
! Energy for the Halo region from the neighbours							!
!=======================================================================================================!
    do wall= 1, 4
    do k= 1, np
       energy_halo(k,wall)=  grv*ht_halo(k,wall) 							&
       			   + 0.5D0*(contrauv_halo(k,wall,1)*couv_halo(k,wall,1)				&
			   + contrauv_halo(k,wall,2)*couv_halo(k,wall,2) )
    enddo
    enddo
!=======================================================================================================!  
!  Primitive Variable Terms: U=[sqrt(det(G(i,j)))*h, u_1, u_2]						!
!		vec(i,j,1)= phi(i,j)									!
!		vec(i,j,2)= couv(i,j,1)									!
!		vec(i,j,3)= couv(i,j,2)									!
!=======================================================================================================!
!  Flux1 term=> fluxvec(i,j,1,1)= sqrt(det(G(i,j)))*(ht-hill)*u^1(contra)				!
!		fluxvec(i,j,1,2)= energy(i,j)								!
!		fluxvec(i,j,1,3)= 0									!
!=======================================================================================================!
!  Flux2 term=>	fluxvec(i,j,2,1)= sqrt(det(G(i,j)))*(ht-hill)*u^2(contra)				!
!		fluxvec(i,j,2,2)= 0									!
!		fluxvec(i,j,2,3)= energy(i,j)								!
!=======================================================================================================!
!		energy(i,j)= g*ht + 0.5*(u_1*u^1 + u_2*u^2)						!
!		h(i,j)	   = ht(i,j)-hill(i,j)								!
!=======================================================================================================! 
!=======================================================================================================!
!    gh11(:,:)= ghij(:,:,1)
!    gh22(:,:)= ghij(:,:,2)
!=======================================================================================================!
!   Call signalwave(contrauv,contrauv_halo,ghij(:,:,1),ghij(:,:,2),gh11_halo,gh22_halo,fjmax,lambda,SL,SR)

    fjmax(:) =  sw_fjmax(contrauv,contrauv_halo,gh11,gh22,gh11_halo,gh22_halo)
!=======================================================================================================!
    do j = 1, np
    do i = 1, np
       vec(i,j,1) = phi(i,j)
       vec(i,j,2) = couv(i,j,1)
       vec(i,j,3) = couv(i,j,2)

       fluxvec(i,j,1,1) = fxy(i,j,1)
       fluxvec(i,j,2,1) = fxy(i,j,2)

       fluxvec(i,j,1,2) = energy(i,j)
       fluxvec(i,j,2,2) = zero      

       fluxvec(i,j,1,3) = zero      
       fluxvec(i,j,2,3) = energy(i,j)
    enddo
    enddo

    vec_halo(:,:,1) = phi_halo(:,:)
    vec_halo(:,:,2) = couv_halo(:,:,1)
    vec_halo(:,:,3) = couv_halo(:,:,2)

    flux_halo(:,:,:,1) = fxy_halo(:,:,:)
    flux_halo(:,:,1,2) = energy_halo(:,:)
    flux_halo(:,:,2,2) = zero          
    flux_halo(:,:,1,3) = zero          
    flux_halo(:,:,2,3) = energy_halo(:,:)
!=======================================================================================================!
    Call swsys_flux(3,deriv,fjmax,vec,vec_halo,fluxvec,flux_halo,numflux)
!=======================================================================================================!
!													!
! Gradient term	=>  grad(i,j,1)= grad^1(sqrt(det(G(i,j)))*h*u^1) + grad^2(sqrt(det(G(i,j))*h*u^2)	!
!		    grad(i,j,2)= grad^1(energy(i,j))							!
!		    grad(i,j,3)= grad^2(energy(i,j))							!
!													!
!=======================================================================================================!
!-------------------------------------------------------------------------------------------------------!
!  Gradient  for the continuity equation     								!
!-------------------------------------------------------------------------------------------------------!
    grad(:,:,1) = gradient_mass(deriv,fxy)

!-------------------------------------------------------------------------------------------------------!
!  Gradient  for the momentum equations									!	
!-------------------------------------------------------------------------------------------------------!
#if 0
   grad(:,:,2:3)=gradient_wk(energy,deriv)
#else
   Call gradient_mom(deriv,energy,gradu1,gradu2)

   do j=1,np
   do i=1,np
      grad(i,j,2) = gradu1(i,j)
      grad(i,j,3) = gradu2(i,j)
   enddo
   enddo
#endif
!=======================================================================================================!
!  Compute RHS of the ODE system corresponding  to the SW model         				!
!  mvi(i,j)= 1 / mv(i,j)										!
!=======================================================================================================!
    do j = 1, np
    do i = 1, np
       rtmp = one/elem%mp(i,j)
       rhs(i,j,1)= (source(i,j,1) + grad(i,j,1) - numflux(i,j,1)) * rtmp
       rhs(i,j,2)= (source(i,j,2) + grad(i,j,2) - numflux(i,j,2)) * rtmp
       rhs(i,j,3)= (source(i,j,3) + grad(i,j,3) - numflux(i,j,3)) * rtmp
    enddo
    enddo

end subroutine  dg_sw_model
!=======================================================================================!
!  Gradient Mass					   				!    
!  Dvv_twt ->  transpose of derivax *G-weights						!
!  Mvv_twt ->  Diadgonal matrix with G-weights						!
!=======================================================================================!
function gradient_mass(deriv,uvflx)  result(gradf)
!=======================================================================================!
    type (derivative_t)              :: deriv
    real(kind=real_kind), intent(in) :: uvflx(np,np,2)
    real(kind=real_kind)             :: gradf(np,np)

    ! Dvv_twt ->  transpose of derivax *G-weights
    ! Mvv_twt ->  Diadgonal matrix with G-weights

    integer i
    integer j
    integer l
    logical, parameter :: UseUnroll = .TRUE.

    real(kind=real_kind)  sumx00,sumx01
    real(kind=real_kind)  sumy00,sumy01
    real(kind=real_kind)  sumx10,sumx11
    real(kind=real_kind)  sumy10,sumy11
!=======================================================================================!
if(MODULO(np,2) == 0 .and. UseUnroll) then 
    do j=1,np,2
       do l=1,np,2
          sumx00=zero
          sumx01=zero
          sumx10=zero
          sumx11=zero

          sumy00=zero
          sumy01=zero
          sumy10=zero
          sumy11=zero

          do i=1,np
             sumx00  = sumx00  + deriv%Dvv_twt(i,l  ) * uvflx(i,  j,1)
             sumx01  = sumx01  + deriv%Dvv_twt(i,l+1) * uvflx(i,  j,1)
             sumx10  = sumx10  + deriv%Dvv_twt(i,l  ) * uvflx(i,j+1,1)
             sumx11  = sumx11  + deriv%Dvv_twt(i,l+1) * uvflx(i,j+1,1)

             sumy00  = sumy00  + deriv%Mvv_twt(i,  l) * uvflx(i,  j,2)
             sumy01  = sumy01  + deriv%Mvv_twt(i,l+1) * uvflx(i,  j,2)
             sumy10  = sumy10  + deriv%Mvv_twt(i,  l) * uvflx(i,j+1,2)
             sumy11  = sumy11  + deriv%Mvv_twt(i,l+1) * uvflx(i,j+1,2)
          end do

          deriv%vvtempt(j  ,l  ,1) = sumx00
          deriv%vvtempt(j  ,l+1,1) = sumx01
          deriv%vvtempt(j+1,l  ,1) = sumx10
          deriv%vvtempt(j+1,l+1,1) = sumx11

          deriv%vvtempt(j  ,l  ,2) = sumy00
          deriv%vvtempt(j  ,l+1,2) = sumy01
          deriv%vvtempt(j+1,l  ,2) = sumy10
          deriv%vvtempt(j+1,l+1,2) = sumy11

       end do
    end do


    do j=1,np,2
       do i=1,np,2
          sumx00=zero
          sumx01=zero
          sumx10=zero
          sumx11=zero

          sumy00=zero
          sumy01=zero
          sumy10=zero
          sumy11=zero

          do l=1,np
             sumx00 = sumx00 +  deriv%Mvv_twt(l,j  )*deriv%vvtempt(l,i  ,1)
             sumx01 = sumx01 +  deriv%Mvv_twt(l,j+1)*deriv%vvtempt(l,i  ,1)
             sumx10 = sumx10 +  deriv%Mvv_twt(l,j  )*deriv%vvtempt(l,i+1,1)
             sumx11 = sumx11 +  deriv%Mvv_twt(l,j+1)*deriv%vvtempt(l,i+1,1)

             sumy00 = sumy00 +  deriv%Dvv_twt(l,j  )*deriv%vvtempt(l,i  ,2)
             sumy01 = sumy01 +  deriv%Dvv_twt(l,j+1)*deriv%vvtempt(l,i  ,2)
             sumy10 = sumy10 +  deriv%Dvv_twt(l,j  )*deriv%vvtempt(l,i+1,2)
             sumy11 = sumy11 +  deriv%Dvv_twt(l,j+1)*deriv%vvtempt(l,i+1,2)
          end do

          gradf(i,j)    = (sumx00 + sumy00) * rrearth
          gradf(i,j+1)  = (sumx01 + sumy01) * rrearth
          gradf(i+1,j)  = (sumx10 + sumy10) * rrearth
          gradf(i+1,j+1)= (sumx11 + sumy11) * rrearth
       end do
    end do
else
    do j=1,np
       do l=1,np
          sumx00=zero
	  sumy00=zero
          do i=1,np
             sumx00  = sumx00  + deriv%Dvv_twt(i,l) * uvflx(i,j,1)
             sumy00  = sumy00  + deriv%Mvv_twt(i,l) * uvflx(i,j,2)
	  enddo
          deriv%vvtempt(j,l,1) = sumx00
          deriv%vvtempt(j,l,2) = sumy00
      enddo
   enddo
    do j=1,np
       do i=1,np
          sumx00=zero
	  sumy00=zero
          do l=1,np
             sumx00 = sumx00 +  deriv%Mvv_twt(l,j)*deriv%vvtempt(l,i,1)
             sumy00 = sumy00 +  deriv%Dvv_twt(l,j)*deriv%vvtempt(l,i,2)
	  enddo
          gradf(i,j) = (sumx00 + sumy00) * rrearth
       enddo
    enddo
endif

!=======================================================================================!
end function  gradient_mass
!=======================================================================================!
!  Gradient Momentum					   				!    
!  Dvv_twt ->  transpose of derivax *G-weights						!
!  Mvv_twt ->  Diadgonal matrix with G-weights						!
!=======================================================================================!
subroutine gradient_mom(deriv,energy,gradu1,gradu2)
!=======================================================================================!
    type (derivative_t)              :: deriv
    real(kind=real_kind), intent(in) :: energy(np,np)
    real(kind=real_kind), intent(out):: gradu1(np,np), gradu2(np,np)
    integer:: i,j,l    
    logical, parameter :: UseUnroll = .TRUE.
    real(kind=real_kind):: sumx00,sumx01
    real(kind=real_kind):: sumy00,sumy01
    real(kind=real_kind):: sumx10,sumx11
    real(kind=real_kind):: sumy10,sumy11
!=======================================================================================!
    !Grad-u 
    !Grad-v 

if(MODULO(np,2) == 0 .and. UseUnroll) then 
    do j=1,np,2
       do l=1,np,2
          sumx00=zero
          sumx01=zero
          sumx10=zero
          sumx11=zero

          sumy00=zero
          sumy01=zero
          sumy10=zero
          sumy11=zero

          do i=1,np
             sumx00 = sumx00 + deriv%Dvv_twt(i,l)  * energy(i,j)
             sumx01 = sumx01 + deriv%Dvv_twt(i,l+1)* energy(i,j)
             sumx10 = sumx10 + deriv%Dvv_twt(i,l)  * energy(i,j+1)
             sumx11 = sumx11 + deriv%Dvv_twt(i,l+1)* energy(i,j+1)

             sumy00 = sumy00 + deriv%Mvv_twt(i,l)  * energy(i,j)
             sumy01 = sumy01 + deriv%Mvv_twt(i,l+1)* energy(i,j)
             sumy10 = sumy10 + deriv%Mvv_twt(i,l)  * energy(i,j+1)
             sumy11 = sumy11 + deriv%Mvv_twt(i,l+1)* energy(i,j+1)
          end do

          deriv%vvtempt(j  ,l  ,1) = sumx00
          deriv%vvtempt(j  ,l+1,1) = sumx01
          deriv%vvtempt(j+1,l  ,1) = sumx10
          deriv%vvtempt(j+1,l+1,1) = sumx11

          deriv%vvtempt(j  ,l  ,2) = sumy00
          deriv%vvtempt(j  ,l+1,2) = sumy01
          deriv%vvtempt(j+1,l  ,2) = sumy10
          deriv%vvtempt(j+1,l+1,2) = sumy11

       end do
    end do




    do j=1,np,2
       do i=1,np,2
          sumx00=zero
          sumx01=zero
          sumx10=zero
          sumx11=zero

          sumy00=zero
          sumy01=zero
          sumy10=zero
          sumy11=zero

          do l=1,np
             sumx00 = sumx00 +  deriv%Mvv_twt(l,j)*deriv%vvtempt(l,i,1)
             sumx01 = sumx01 +  deriv%Mvv_twt(l,j+1)*deriv%vvtempt(l,i,1)
             sumx10 = sumx10 +  deriv%Mvv_twt(l,j)*deriv%vvtempt(l,i+1,1)
             sumx11 = sumx11 +  deriv%Mvv_twt(l,j+1)*deriv%vvtempt(l,i+1,1)

             sumy00 = sumy00 +  deriv%Dvv_twt(l,j)*deriv%vvtempt(l,i,2)
             sumy01 = sumy01 +  deriv%Dvv_twt(l,j+1)*deriv%vvtempt(l,i,2)
             sumy10 = sumy10 +  deriv%Dvv_twt(l,j)*deriv%vvtempt(l,i+1,2)
             sumy11 = sumy11 +  deriv%Dvv_twt(l,j+1)*deriv%vvtempt(l,i+1,2)
          end do
          gradu1(i,j)    = sumx00*rrearth
          gradu1(i,j+1)  = sumx01*rrearth
          gradu1(i+1,j)  = sumx10*rrearth
          gradu1(i+1,j+1)= sumx11*rrearth

          gradu2(i,j)    = sumy00*rrearth
          gradu2(i,j+1)  = sumy01*rrearth
          gradu2(i+1,j)  = sumy10*rrearth
          gradu2(i+1,j+1)= sumy11*rrearth
       end do
    end do
else
    do j=1,np
       do l=1,np
          sumx00=zero
          sumy00=zero
          do i=1,np
             sumx00 = sumx00 + deriv%Dvv_twt(i,l) * energy(i,j)
             sumy00 = sumy00 + deriv%Mvv_twt(i,l) * energy(i,j)
	  enddo
          deriv%vvtempt(j,l,1) = sumx00
          deriv%vvtempt(j,l,2) = sumy00
      enddo
    enddo
    do j=1,np
       do i=1,np
          sumx00=zero
          sumy00=zero
          do l=1,np
             sumx00 = sumx00 +  deriv%Mvv_twt(l,j)*deriv%vvtempt(l,i,1)
             sumy00 = sumy00 +  deriv%Dvv_twt(l,j)*deriv%vvtempt(l,i,2)
          enddo
          gradu1(i,j) = sumx00*rrearth
          gradu2(i,j) = sumy00*rrearth
       enddo
    enddo
endif
end subroutine  gradient_mom
!=======================================================================================!

function  uv_vorticity(uv,sg,D,deriv) result(vor) 
!==========================================================================================
!   Implicit None
    type (derivative_t)         :: deriv
    real (kind=real_kind), intent(in) :: D(2,2,np,np)
    real (kind=real_kind), dimension(np,np), intent(in) :: sg  
    real (kind=real_kind), dimension(np,np,2), intent(in) :: uv
    real (kind=real_kind), dimension(np,np,2) :: couv 
    real (kind=real_kind), dimension(np,np) :: vor
    real(kind=real_kind) ::  dvdx00, dudy00
    real (kind=real_kind):: term,delm
    integer:: i,j,k,l    

    couv(:,:,:) = sphere2cov(uv,D) 

    do j=1,np
       do l=1,np
          dudy00=0.0D0 
          dvdx00=0.0D0 
          do i=1,np
             dvdx00 = dvdx00 + deriv%Dvv(i,l)* couv(i,j,2)
             dudy00 = dudy00 + deriv%Dvv(i,l)* couv(j,i,1)
          enddo
          vor(l,j) = dvdx00
          deriv%vvtemp(j,l) = dudy00
       enddo
    enddo

    do j=1,np
       do i=1,np
          vor(i,j)=(vor(i,j)-deriv%vvtemp(i,j))*rrearth /sg(i,j) 
       end do
    end do


!==========================================================================================
end function  uv_vorticity 
!==========================================================================================
subroutine  source_term(deriv,mmx,gcori,contrauv,couv,source)
!==========================================================================================
    Implicit None
    type (derivative_t)         :: deriv
    real (kind=real_kind), intent(in), dimension(np,np) :: mmx
    real (kind=real_kind), dimension(np,np), intent(in) :: gcori
    real (kind=real_kind), dimension(np,np,2), intent(in) :: contrauv, couv  
    real (kind=real_kind), dimension(np,np,3), intent(out):: source     
    real (kind=real_kind), dimension(np,np) :: vor, vort
    real(kind=real_kind) ::  dvdx00,dvdx01
    real(kind=real_kind) ::  dvdx10,dvdx11
    real(kind=real_kind) ::  dudy00,dudy01
    real(kind=real_kind) ::  dudy10,dudy11
    real (kind=real_kind):: term,delm
    integer:: i,j,k,l    
    logical, parameter :: UseUnroll = .TRUE.
!==========================================================================================
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
             dvdx00 = dvdx00 + deriv%Dvv(i,l)* couv(i,j,2)
             dudy00 = dudy00 + deriv%Dvv(i,l)* couv(j,i,1)
	  enddo
          vor(l,j) = dvdx00
          deriv%vvtemp(j,l) = dudy00
       enddo
    enddo
endif


    do j=1,np
       do i=1,np
          vor(i,j)=(vor(i,j)-deriv%vvtemp(i,j))*rrearth 
       end do
    end do


    do j = 1, np
       do i = 1, np

          term =  (vor(i,j) +  gcori(i,j)) * mmx(i,j) 

          source(i,j,1) =  zero
          source(i,j,2) =  contrauv(i,j,2) * term
          source(i,j,3) = -contrauv(i,j,1) * term

       enddo
    enddo
!==========================================================================================
end subroutine  source_term
!=++++++++++++++++======================================================================!
! LDG & Diffusion Zone   [nu*Del^2(U)r]
!=++++++++++++++++======================================================================!
subroutine dg_diff_flux(elem,deriv,gradbuf,grads,visflx)
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
!  1/(4*dx*dy) scaling 

     if (ne.ne.0) then
         tmp = dble(2*ne)/acos(0.0d0)
         d_nu = nu*tmp*tmp
     else
        d_nu = nu   
     end if

  !     d_nu = nu  
  ! print *, 'Running dg3d_diff_flux with nu = ', d_nu

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

    dgrad(:,:) = gradient_mass(deriv,grads)

    do j = 1,np
    do i = 1,np
      visflx(i,j) = d_nu * (cflxint(i,j) - dgrad(i,j)) / elem%mp(i,j)
    enddo
    enddo
!=======================================================================================================!
end subroutine dg_diff_flux
!=======================================================================================================!
subroutine dg_diff_grads(elem,deriv,contrauvbuf,contrauv,couv,dif_gradu,dif_gradv)
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
end subroutine dg_diff_grads
!=======================================================================================================!
!=======================================================================================================!
subroutine dg_diff_grads_uv(elem,deriv,uvbuf,uv,dif_gradu,dif_gradv)
!=======================================================================================!
    Implicit None
    integer, parameter :: south=1, east=2, north=3, west=4
    type(element_t), intent(in) :: elem
    type (derivative_t),                               intent(in)   :: deriv
    real (kind=real_kind), dimension(0:np+1,0:np+1,2), intent(in)   :: uvbuf
    real (kind=real_kind), dimension(np,np,2),         intent(in)   :: uv
    real (kind=real_kind), dimension(np,np,2),         intent(out)  :: dif_gradu, dif_gradv
    real (kind=real_kind), dimension(np,np,2)  :: couv, gradu1, gradu2

    real (kind=real_kind), dimension(np,4,2)   :: uv_halo, couv_halo 
    real (kind=real_kind), dimension(np,np,2)  :: jflx
    real (kind=real_kind), dimension(np,np)    :: jfu,jfv,cu,cv
    Integer :: i,j,k, wall, eqn
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
!   Call  jump_fluxint(deriv,uv,uv_halo,jflx)
         !cu(:,:) = uv(:,:,1)
          cu(:,:) = couv(:,:,1) 
         jfu(:,:) = jflx(:,:,1)
    Call  ldg_grads(elem,deriv,cu,jfu,gradu1)

         !cv(:,:) = uv(:,:,2)
          cv(:,:) = couv(:,:,2) 
         jfv(:,:) = jflx(:,:,2)
    Call  ldg_grads(elem,deriv,cv,jfv,gradu2)

      dif_gradu(:,:,:) = gradu1(:,:,:) 
      dif_gradv(:,:,:) = gradu2(:,:,:) 

!=======================================================================================================!
end subroutine dg_diff_grads_uv
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

! do j=1,np
!  do i=1,np
!     dif_gradu(i,j,1) = gradu1(i,j,1) *elem%Dinv(1,1,i,j) + gradu1(i,j,2) * elem%Dinv(2,1,i,j) 
!     dif_gradu(i,j,2) = gradu1(i,j,1) *elem%Dinv(1,2,i,j) + gradu1(i,j,2) * elem%Dinv(2,2,i,j) 
!
!     dif_gradv(i,j,1) = gradu2(i,j,1) *elem%Dinv(1,1,i,j) + gradu2(i,j,2) * elem%Dinv(2,1,i,j) 
!     dif_gradv(i,j,2) = gradu2(i,j,1) *elem%Dinv(1,2,i,j) + gradu2(i,j,2) * elem%Dinv(2,2,i,j) 
!  enddo
!  enddo


end subroutine ldg_grads
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!=======================================================================================!
!  Advection Only Zone							   		!
!=======================================================================================!
subroutine dg_adv_model(elem,deriv,contrauvbuf,htbuf,contrauv,ht,h_rhs)
!=======================================================================================!
    use element_mod, only : element_t
    Implicit None
    type (element_t) , intent(in), target :: elem
    type (derivative_t)  , intent(in) :: deriv
    real (kind=real_kind), dimension(0:np+1,0:np+1),   intent(in)   :: htbuf
    real (kind=real_kind), dimension(0:np+1,0:np+1,2), intent(in)   :: contrauvbuf
    real (kind=real_kind), dimension(np,np,2),         intent(in)   :: contrauv
    real (kind=real_kind), dimension(np,np),           intent(in)   :: ht
    real (kind=real_kind), dimension(np,np),           intent(out)  :: h_rhs
    real (kind=real_kind), dimension(np,np,2) :: fxy
    real (kind=real_kind), dimension(np,np)   :: phi, flux,uflx,vflx, hill,grad,source
    real (kind=real_kind), dimension(np,4,2)  :: fxy_halo, contrauv_halo
    real (kind=real_kind), dimension(np,4)    :: fx_halo, fy_halo, phi_halo, ht_halo
    real (kind=real_kind), dimension(np,np)   :: gradn ,  fluxn
    real (kind=real_kind)  :: s1,s2,trm , fact
    Integer  :: i,j, wall, k, l
!=======================================================================================!
    !  From HOMME

    metdet => elem%metdet
    mv => elem%mp

    !Inverted mass-matrix 
    do j = 1, np
       do i = 1, np
          mmxi(i,j) = one / mv(i,j)
       enddo
    enddo

    hill  = zero
    source= zero

    !--------------------------------------------------------------
    ! Boundary values for velocity & flux terms from the neighbors
    !--------------------------------------------------------------


    do i = 1, np

       ! U-component on  S,E,N,W  (4) sides of an element

       contrauv_halo(i,1,1) = contrauvbuf(i,   0,1)
       contrauv_halo(i,2,1) = contrauvbuf(np+1,i,1)
       contrauv_halo(i,3,1) = contrauvbuf(i,np+1,1)
       contrauv_halo(i,4,1) = contrauvbuf(0,i,   1)

       ! V-component on  S,E,N,W  (4) sides of an element

       contrauv_halo(i,1,2) = contrauvbuf(i,   0,2)
       contrauv_halo(i,2,2) = contrauvbuf(np+1,i,2)
       contrauv_halo(i,3,2) = contrauvbuf(i,np+1,2)
       contrauv_halo(i,4,2) = contrauvbuf(0,i,   2)

       ! H-field (scalar) on S, E, N, W sides of an element

       ht_halo(i,1) = htbuf(i,   0)
       ht_halo(i,2) = htbuf(np+1,i)
       ht_halo(i,3) = htbuf(i,np+1)
       ht_halo(i,4) = htbuf(0,i   )

    enddo


    !--------------------------------------------------------------
    ! Flux computations for the Continuity equation
    !--------------------------------------------------------------

    !(u,v) Fluxes for the continuity equation
    !Note   "ht - hill => " Depth of the fluid

    do j = 1, np
       do i = 1, np
          phi(i,j) = (ht(i,j) - hill(i,j)) * metdet(i,j)
          fxy(i,j,1) = contrauv(i,j,1) * phi(i,j)
          fxy(i,j,2) = contrauv(i,j,2) * phi(i,j)
       enddo
    enddo

    !(u,v) Fluxes along the Halo  (S,E,N,W) from the neighbours

    do wall = 1, 4
       do k = 1, np

          phi_halo(k,1) = (ht_halo(k,1) - hill(k,1)) * metdet(k,1)
          phi_halo(k,2) = (ht_halo(k,2) - hill(np,k))* metdet(np,k)
          phi_halo(k,3) = (ht_halo(k,3) - hill(k,np))* metdet(k,np)
          phi_halo(k,4) = (ht_halo(k,4) - hill(1,k)) * metdet(1,k)

          fxy_halo(k,wall,1) = contrauv_halo(k,wall,1)* phi_halo(k,wall)
          fxy_halo(k,wall,2) = contrauv_halo(k,wall,2)* phi_halo(k,wall)
       enddo
    enddo


    !----------------------------------------
    ! Fluxes for the continuity equation
    !----------------------------------------

    fluxn(:,:) =  adv_flux_term(deriv,contrauv,contrauv_halo,phi,phi_halo,fxy,fxy_halo)

    !----------------------------------------
    ! Gardient term [ F.Del(xi) ]
    !----------------------------------------

    gradn(:,:) = gradient_mass(deriv,fxy)

    !----------------------------------------
    ! Rhs of the corresponding ODE
    !----------------------------------------

    do j = 1, np
       do i = 1, np
          h_rhs(i,j) = (source(i,j) + gradn(i,j) - fluxn(i,j)) * mmxi(i,j)
       enddo
    enddo

!=======================================================================================!
end subroutine  dg_adv_model
!=======================================================================================!
!  Advection Flux Term							   		!
!=======================================================================================!
function adv_flux_term(deriv,contrauv,contrauv_halo,si,si_senw,fxy,fxy_halo) result(numflux)
!=======================================================================================!
    Implicit None
    integer, parameter :: south=1, east=2, north=3, west=4
    type (derivative_t)         :: deriv
    real (kind=real_kind), dimension(np,np,2),intent(in) :: contrauv, fxy
    real (kind=real_kind), dimension(np,np),  intent(in) :: si
    real (kind=real_kind), dimension(np,4,2), intent(in) :: contrauv_halo,fxy_halo
    real (kind=real_kind), dimension(np,4),   intent(in) :: si_senw
    real (kind=real_kind), dimension(np,np) :: numflux
    real (kind=real_kind), dimension(np,np) :: mij
    real (kind=real_kind), dimension(np)   :: lf_south,lf_north,lf_east,lf_west
    real(kind=real_kind) ::  alfa1, alfa2, ul,ur , left, right, f_left, f_right, s1,s2
    integer i,j, k
!=======================================================================================!
    mij(:,:)  = zero
    mij(1,1)  = one
    mij(np,np)= one 

    ! South & North    Max flux Jacobians

    alfa1 = zero
    alfa2 = zero

    do i = 1, np
       ul = abs(contrauv(i,1,2))
       ur = abs(contrauv_halo(i,south,2))
       alfa1 = max(alfa1,ul,ur)
       ul = abs(contrauv_halo(i,north,2))
       ur = abs(contrauv(i,np,2))
       alfa2 = max(alfa2,ul,ur)
    enddo

    do i = 1, np
       ! South wall
       f_left = fxy_halo(i,south,2)
       f_right = fxy(i,1,2)
       left = si_senw(i,south)
       right  = si(i,1)
       lf_south(i) =  0.5D0 *(f_left + f_right - alfa1*(right - left))
       ! North wall
       f_left  = fxy(i,np,2)
       f_right = fxy_halo(i,north,2)
       left  = si(i,np)
       right = si_senw(i,north)
       lf_north(i) =  0.5D0 *(f_left + f_right - alfa2*(right - left))
    enddo


    ! East & West   max of Flux Jacobians

    alfa1 = zero
    alfa2 = zero

    do j = 1, np
       ul = abs(contrauv_halo(j,west,1))
       ur = abs(contrauv(1,j,1))
       alfa1 = max(alfa1,ul,ur)
       ul = abs(contrauv(np,j,1))
       ur = abs(contrauv_halo(j,east,1))
       alfa2 = max(alfa2,ul,ur)
    enddo

    do j = 1, np
       !West wall
       f_left  = fxy_halo(j,west,1)
       f_right = fxy(1,j,1)
       left  =  si_senw(j,west)
       right  = si(1,j)
       lf_west(j) =  0.5D0 *(f_left + f_right - alfa1*(right - left))

       !East wall
       f_left  = fxy(np,j,1)
       f_right = fxy_halo(j,east,1)
       left  = si(np,j)
       right  = si_senw(j,east)
       lf_east(j) =  0.5D0 *(f_left + f_right - alfa2*(right - left))

    enddo

    !Flux integral along the element boundary


    do j = 1, np
       do i = 1, np

          s1 =  (lf_east(j) *mij(i,np) - lf_west(j) *mij(i,1) )* deriv%Mvv_twt(j,j)
          s2 =  (lf_north(i)*mij(j,np) - lf_south(i)*mij(j,1) )* deriv%Mvv_twt(i,i)

          numflux(i,j) = (s1 + s2) * rrearth

       enddo
    enddo
!=======================================================================================!
end function adv_flux_term
!=======================================================================================!
#endif
!=======================================================================================!
function contra2co(vin,met) result(vout)
!=======================================================================================!
real (kind=real_kind),intent(in) :: met(2,2,np,np)
real (kind=real_kind),intent(in) :: vin(np,np,2)
real (kind=real_kind)            :: vout(np,np,2)
integer :: i,j
real (kind=real_kind):: v1,v2
!=======================================================================================!
do j=1,np
do i=1,np
   v1= vin(i,j,1)
   v2= vin(i,j,2)
   vout(i,j,1)= met(1,1,i,j)*v1 + met(1,2,i,j)*v2
   vout(i,j,2)= met(2,1,i,j)*v1 + met(2,2,i,j)*v2
enddo
enddo
end function contra2co
!=======================================================================================!
function co2contra(vin,metinv) result(vout)
!=======================================================================================!
real (kind=real_kind),intent(in) :: metinv(2,2,np,np)
real (kind=real_kind),intent(in) :: vin(np,np,2)
real (kind=real_kind)            :: vout(np,np,2)
integer :: i,j
real (kind=real_kind):: v1,v2
!=======================================================================================!
do j=1,np
do i=1,np
   v1= vin(i,j,1)
   v2= vin(i,j,2)
   vout(i,j,1)= metinv(1,1,i,j)*v1 + metinv(1,2,i,j)*v2
   vout(i,j,2)= metinv(2,1,i,j)*v1 + metinv(2,2,i,j)*v2
enddo
enddo
end function co2contra
!=======================================================================================!
function sphere2cov(vin,D) result(vout)
!=======================================================================================!
real (kind=real_kind),intent(in) :: D(2,2,np,np)
real (kind=real_kind),intent(in) :: vin(np,np,2)
real (kind=real_kind)            :: vout(np,np,2)
integer:: i,j
real (kind=real_kind):: v1,v2
!
!  [u_1,u_2] = A^T [u, v] 
!
  do j=1,np
  do i=1,np
    v1= vin(i,j,1)
    v2= vin(i,j,2)
    vout(i,j,1) = D(1,1,i,j)*v1 + D(2,1,i,j)*v2
    vout(i,j,2) = D(1,2,i,j)*v1 + D(2,2,i,j)*v2
  enddo
  enddo
 return
end function sphere2cov   
!=======================================================================================!
function cov2sphere(vin,Dinv) result(vout)
!=======================================================================================!
real (kind=real_kind),intent(in) :: Dinv(2,2,np,np)
real (kind=real_kind),intent(in) :: vin(np,np,2)
real (kind=real_kind)            :: vout(np,np,2)
integer:: i,j
real (kind=real_kind):: v1,v2
!
!  [u,v] = A^-T [u_1,u_2] 
!
  do j=1,np
  do i=1,np
    v1= vin(i,j,1)
    v2= vin(i,j,2)
    vout(i,j,1) = Dinv(1,1,i,j)*v1 + Dinv(2,1,i,j)*v2
    vout(i,j,2) = Dinv(1,2,i,j)*v1 + Dinv(2,2,i,j)*v2
  enddo
  enddo
 return
end function cov2sphere
!=======================================================================================!
function sphere2contra(vin,Dinv) result(vout)
!=======================================================================================!
real (kind=real_kind),intent(in) :: Dinv(2,2,np,np)
real (kind=real_kind),intent(in) :: vin(np,np,2)
real (kind=real_kind)            :: vout(np,np,2)
integer:: i,j
real (kind=real_kind):: v1,v2
!
!  [u^1,u^2] = A^-1 [u, v] 
!
  do j=1,np
  do i=1,np
    v1= vin(i,j,1)
    v2= vin(i,j,2)
    vout(i,j,1) = Dinv(1,1,i,j)*v1 + Dinv(1,2,i,j)*v2
    vout(i,j,2) = Dinv(2,1,i,j)*v1 + Dinv(2,2,i,j)*v2
  enddo
  enddo
 return
end function sphere2contra
!=======================================================================================!
function contra2sphere(vin,D) result(vout)
!=======================================================================================!
real (kind=real_kind),intent(in) :: D(2,2,np,np)
real (kind=real_kind),intent(in) :: vin(np,np,2)
real (kind=real_kind)            :: vout(np,np,2)
integer:: i,j
real (kind=real_kind):: v1,v2
!
!  [u,v] = A [u^1, u^2] 
!
  do j=1,np
  do i=1,np
    v1= vin(i,j,1)
    v2= vin(i,j,2)
    vout(i,j,1) = D(1,1,i,j)*v1 + D(1,2,i,j)*v2
    vout(i,j,2) = D(2,1,i,j)*v1 + D(2,2,i,j)*v2
  enddo
  enddo
!=======================================================================================!
 return
!=======================================================================================!
end function contra2sphere
!=======================================================================================!
!=======================================================================================!
function height2phi(vin,metdet) result(vout)
!=======================================================================================!
real (kind=real_kind),intent(in) :: metdet(np,np)
real (kind=real_kind),intent(in) :: vin(np,np)
real (kind=real_kind)            :: vout(np,np)
integer:: i,j
!=======================================================================================!  
  do j=1,np
  do i=1,np
    vout(i,j)= vin(i,j)*metdet(i,j)
  enddo
  enddo
!=======================================================================================!
 return
!=======================================================================================!
end function height2phi
!=======================================================================================!
!=======================================================================================!
function phi2height(vin,metdet) result(vout)
!=======================================================================================!
real (kind=real_kind),intent(in) :: metdet(np,np)
real (kind=real_kind),intent(in) :: vin(np,np)
real (kind=real_kind)            :: vout(np,np)
integer:: i,j
!=======================================================================================!  
  do j=1,np
  do i=1,np
    vout(i,j) = vin(i,j)/metdet(i,j)
  enddo
  enddo
!=======================================================================================!
 return
!=======================================================================================!
end function phi2height
!=======================================================================================!
#ifdef _PRIMDG
!=======================================================================================!
function psi2height(vin,grv) result(vout)
!=======================================================================================!
real (kind=real_kind),intent(in) :: grv
real (kind=real_kind),intent(in) :: vin(np,np)
real (kind=real_kind)            :: vout(np,np)
integer:: i,j
!=======================================================================================!  
  do j=1,np
  do i=1,np
    vout(i,j) = vin(i,j)/grv
  enddo
  enddo
!=======================================================================================!
 return
!=======================================================================================!
end function psi2height
!=======================================================================================!
function dp2pt(vin,vp) result(vout)
!=======================================================================================!
real (kind=real_kind),intent(in) :: vp(np,np)
real (kind=real_kind),intent(in) :: vin(np,np)
real (kind=real_kind)            :: vout(np,np)
integer:: i,j
!=======================================================================================!  
  do j=1,np
  do i=1,np
    vout(i,j)= vin(i,j)*vp(i,j)
  enddo
  enddo
!=======================================================================================!
 return
!=======================================================================================!
end function dp2pt
!=======================================================================================!
function pt2dp(vin,vp) result(vout)
!=======================================================================================!
real (kind=real_kind),intent(in) :: vp(np,np)
real (kind=real_kind),intent(in) :: vin(np,np)
real (kind=real_kind)            :: vout(np,np)
integer:: i,j
!=======================================================================================!  
  do j=1,np
  do i=1,np
    vout(i,j)= vin(i,j)/vp(i,j)
  enddo
  enddo
!=======================================================================================!
 return
!=======================================================================================!
end function pt2dp
!=======================================================================================!    
#endif
!=======================================================================================! 
subroutine swsys_flux(numeqn,deriv,fjmax,si,si_senw,uvflx,uvflx_senw,fluxout)

   integer, parameter :: south=1, east=2, north=3, west=4
   integer, intent(in):: numeqn
   type (derivative_t)                                  :: deriv
   real (kind=real_kind), dimension(4),    intent(in)   :: fjmax
   real (kind=real_kind), dimension(np,np,2,numeqn),intent(in) :: uvflx
   real (kind=real_kind), dimension(np,np,numeqn),  intent(in) :: si
   real (kind=real_kind), dimension(np,4,2,numeqn), intent(in) :: uvflx_senw
   real (kind=real_kind), dimension(np,4,numeqn),   intent(in) :: si_senw

   real (kind=real_kind), dimension(np,np,numeqn), intent(out) :: fluxout

   real (kind=real_kind), dimension(np,np) :: mij
   real (kind=real_kind), dimension(np)   :: lf_south,lf_north,lf_east,lf_west
   real(kind=real_kind) ::  fj(4), ul,ur , left, right, f_left, f_right, s1,s2

   integer i,j,k,eqn 

      mij(:,:) = 0.0D0
      mij(1,1) = 1.0D0
    mij(np,np) = 1.0D0

       fj(:) = scale*fjmax(:)

    !For SW-system  (u1,u2,dp,pt) order   (4 equations)

    fluxout(:,:,:) = 0.0D0 

      do eqn = 1, numeqn


           ! East & West   LF flux  (fjmax <- max of flux Jacobian)

         do j = 1, np

                  f_left  = uvflx_senw(j,west,1,eqn)
                  f_right = uvflx(1,j,1,eqn)
                    left  = si_senw(j,west,eqn)
                   right  = si(1,j,eqn)
               lf_west(j) = 0.5D0 *(f_left + f_right - fj(west)*(right - left))

                  f_left  = uvflx(np,j,1,eqn)
                  f_right = uvflx_senw(j,east,1,eqn)
                    left  = si(np,j,eqn)
                   right  = si_senw(j,east,eqn)
               lf_east(j) = 0.5D0 *(f_left + f_right - fj(east)*(right - left))

          end do

           ! North& South  LF flux  (fjmax <- max of flux Jacobian)

         do i = 1, np

                   f_left = uvflx_senw(i,south,2,eqn)
                  f_right = uvflx(i,1,2,eqn)
                     left = si_senw(i,south,eqn)
                   right  = si(i,1,eqn)
              lf_south(i) = 0.5D0 *(f_left + f_right - fj(south)*(right - left))

                  f_left  = uvflx(i,np,2,eqn)
                  f_right = uvflx_senw(i,north,2,eqn)
                    left  = si(i,np,eqn)
                    right = si_senw(i,north,eqn)
              lf_north(i) = 0.5D0 *(f_left + f_right - fj(north)*(right - left))

         end do

        !Flux integral along the element boundary

        do j = 1, np
        do i = 1, np
            s1 = (lf_east(j) *mij(i,np) - lf_west(j) *mij(i,1) )* deriv%Mvv_twt(j,j)
            s2 = (lf_north(i)*mij(j,np) - lf_south(i)*mij(j,1) )* deriv%Mvv_twt(i,i)
            fluxout(i,j,eqn) = (s1 + s2) * rrearth
        end do
        end do

    end do 

 end subroutine  swsys_flux 
!=======================================================================================================! 
!=======================================================================================================!
 subroutine mono_filter(elem,gll,dt,htbuf,contrauv,ht,ht_new)
!=======================================================================================!
    Implicit None
    type (element_t) , intent(in), target :: elem
    real (kind=real_kind), intent(in) :: dt
    real (kind=real_kind), dimension(0:np+1,0:np+1),   intent(in)   :: htbuf
    real (kind=real_kind), dimension(np,np,2),         intent(in)   :: contrauv
    real (kind=real_kind), dimension(np,np),           intent(in)   :: ht
    real (kind=real_kind), dimension(np,np),           intent(inout)  :: ht_new
    type (quadrature_t):: gll

    real (kind=real_kind), dimension(np,np,2) :: fxy
    real (kind=real_kind), dimension(0:np+1,0:np+1)  :: fld
    real (kind=real_kind), dimension(np,np)   :: phi, flux,uflx,vflx, hill,grad,source
    real (kind=real_kind), dimension(np,np)   :: uu,vv
    real (kind=real_kind), dimension(np)      :: glp
    real (kind=real_kind), dimension(0:np+1)  :: egl
    real (kind=real_kind)  :: s1,s2,s3,s4,fact,cox,coy
    real (kind=real_kind)  :: fval,fmin,fmax
    Integer  :: i,j, wall, k, l, ix,iy
!=======================================================================================!
    !  From HOMME
    metdet => elem%metdet
    mv => elem%mp
      
    do j = 1, np
       do i = 1, np
          fld(i,j) = ht(i,j)
          uu(i,j) = contrauv(i,j,1)
          vv(i,j) = contrauv(i,j,2)
       enddo
    enddo


  ! extended gll points with halo regions
       do  k = 1, np
        glp(k) = gll%points(k)
        egl(k) = glp(k)
       enddo

       egl(0) = egl(1) - (glp(2) - glp(1))
       egl(np+1) = egl(np) + (glp(np) - glp(np-1))

    !--------------------------------------------------------------
    ! Boundary values for given fields (padding)
    !--------------------------------------------------------------

    do i = 1, np
       fld(i,0) = htbuf(i,   0)
       fld(np+1,i) = htbuf(np+1,i)
       fld(i,np+1) = htbuf(i,np+1)
       fld(0,i) = htbuf(0,i   )
    enddo

      fld(0,0) = (fld(0,1) + fld(1,1) + fld(1,0))/3.0D0
      fld(np+1,np+1) = (fld(np+1,np) + fld(np,np) + fld(np,np+1))/3.0D0
      fld(np+1,0) = (fld(np,0) + fld(np,1) + fld(np+1,1))/3.0D0
      fld(0,np+1) = (fld(0,np) + fld(1,np) + fld(1,np+1))/3.0D0

    !--------------------------------------------------------------
    ! Flux computations for the Continuity equation
    !--------------------------------------------------------------


    do j = 1, np
       do i = 1, np

              cox = glp(i) -  uu(i,j) * dt
              coy = glp(j) -  vv(i,j) * dt

                ix =indexer(egl,cox)
                iy =indexer(egl,coy)

    !           idx(i,j) = ix
    !           idy(i,j) = iy

                 s1 = fld(ix,iy)
                 s2 = fld(ix+1,iy)
                 s3 = fld(ix+1,iy+1)
                 s4 = fld(ix,iy+1)

                fmin = min(s1,s2,s3,s4)
                fmax = max(s1,s2,s3,s4)
                fval = ht_new(i,j)

              if (fval < fmin) ht_new(i,j) = fmin
              if (fval > fmax) ht_new(i,j) = fmax
       enddo
    enddo

end subroutine mono_filter

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Searching the position of (xp) on "egl" grid (Bisection)
function indexer(grid,xp) result(ii)
        implicit None
          real (kind=real_kind), intent(in) :: xp
          real (kind=real_kind), intent(in), dimension(0:np+1):: grid
          integer  :: ii, nm,na,nb
            na = 0
            nb = np+1
             do
               if  ((nb-na) <=  1)  exit
               nm = (nb + na)/2
                if (xp  >  grid(nm)) then
                 na = nm
                else
                 nb = nm
                endif
             enddo

              ii = na

end function indexer
!----------------------------------------------  
end module dg_core_mod


