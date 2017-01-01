#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
!
!  getter and setter functions that must be provided by each model
!  
!  IMPORTANT NOTE:  For vertically lagrangian models, these
!  routines should ONLY be used outside the timestepping loop
!  on reference levels.  To compute these fields on floating levels
!  the model should do that directly (or we need to modify the interface)
!
!  see src/preqx/element_ops.F90 for documentation
!     
!  helper fuctions in here not required to be provided by all models:
!     get_pfull_and_exner
!     get_temperature
!     get_pottemp
!    
!
module element_ops

  use dimensions_mod, only: np, nlev, nlevp, nelemd, qsize, max_corner_elem
  use element_mod,    only: element_t
  use element_state,  only: timelevels
  use hybrid_mod,     only: hybrid_t
  use hybvcoord_mod,  only: hvcoord_t
  use kinds,          only: real_kind, iulog
  use perf_mod,       only: t_startf, t_stopf, t_barrierf, t_adj_detailf ! _EXTERNAL
  use parallel_mod,   only: abortmp
  use physical_constants, only : p0, Cp, Rgas, Rwater_vapor, Cpwater_vapor, kappa
  use control_mod, only:    use_moisture, use_cpstar
  use prim_si_mod, only: preq_hydrostatic_v2
  implicit none

contains

  subroutine get_field(elem,name,field,hvcoord,nt,ntQ)
  implicit none
  type (element_t), intent(in) :: elem
  character(len=*), intent(in) :: name
  real (kind=real_kind), intent(out)  :: field(np,np,nlev)
  type (hvcoord_t),     intent(in)    :: hvcoord          
  integer, intent(in) :: nt
  integer, intent(in) :: ntQ

  select case(name)
  case ('temperature')
     call get_temperature(elem,field,hvcoord,nt,ntQ)
  case ('pottemp')
     call get_pottemp(elem,field,hvcoord,nt,ntQ)
  case ('phi')
     field = elem%state%phi(:,:,:,nt)
     !field = elem%derived%phi(:,:,:)
  case default
     print *,'name = ',trim(name)
     call abortmp('ERROR: get_field name not supported in this model')
  end select

  end subroutine



  subroutine get_pottemp(elem,pottemp,hvcoord,nt,ntQ)
  implicit none
    
  type (element_t), intent(in)        :: elem
  real (kind=real_kind), intent(out)  :: pottemp(np,np,nlev)
  type (hvcoord_t),     intent(in)    :: hvcoord                      ! hybrid vertical coordinate struct
  integer, intent(in) :: nt
  integer, intent(in) :: ntQ
  
  !   local
  real (kind=real_kind) :: p(np,np,nlev)
  real (kind=real_kind) :: dp(np,np,nlev)
  real (kind=real_kind) :: kappa_star(np,np,nlev)
  real (kind=real_kind) :: Qt(np,np,nlev)
  integer :: k
  
  pottemp(:,:,:) = elem%state%theta(:,:,:,nt)
  
  end subroutine get_pottemp
  


  subroutine get_temperature(elem,temperature,hvcoord,nt,ntQ)
  implicit none
  
  type (element_t), intent(in)        :: elem
  real (kind=real_kind), intent(out)  :: temperature(np,np,nlev)
  type (hvcoord_t),     intent(in)    :: hvcoord                      ! hybrid vertical coordinate struct
  integer, intent(in) :: nt
  integer, intent(in) :: ntQ
  
  !   local
  real (kind=real_kind) :: dp(np,np,nlev)
  real (kind=real_kind) :: exner(np,np,nlev)
  real (kind=real_kind) :: pfull(np,np,nlev)
  real (kind=real_kind) :: pfull_i(np,np,nlevp)
  integer :: k
  
  
#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
  do k=1,nlev
     dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
          ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem%state%ps_v(:,:,nt)
  enddo

  call get_pfull_and_exner(pfull,pfull_i,exner,hvcoord,elem%state%theta(:,:,:,nt),&
       dp,elem%state%phi(:,:,:,nt),elem%state%phis,elem%state%Qdp(:,:,:,1,ntQ))


#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
  do k=1,nlev
     temperature(:,:,k)= elem%state%theta(:,:,k,nt)*exner(:,:,k)
  enddo

  end subroutine get_temperature




  subroutine get_pfull_and_exner(pfull,pfull_i,exner,hvcoord,theta,dp3d,phi,phis,Qdp)
  implicit none
!
! compute exner pressure, full presure
!
!
! hydrostatic case:
!    input:  dp3d, Qdp (if use_moisture)
!    output:  pfull = compute from dp3d
!             exner = pfull**kappa_star
!

! nonhydro formula used in hydrostatic model (used to debug this routine)
!         call this routine with phi computed from hydrostatic exner
!         
! nonhydro case:
! input:  dp3d, Qdp (if use_moisture), phi, phis, theta
! output:  
!       for k=2..nlev, use the equation of state:  pfull/e = rho*Rstar*theta
!                                                    rho   = -dp3d/dphi 
! for k=1, we cant compute rho because we dont know phi at the model top
! (our boundary condition specifies ptop) so we take p(k=1) = 2/3 ptop + 1/3 p(k=2)
!
!  
  real (kind=real_kind), intent(out) :: exner(np,np,nlev)  ! exner full pressure
  real (kind=real_kind), intent(out) :: pfull(np,np,nlev)   ! full nonhyrdo pressure
  real (kind=real_kind), intent(out) :: pfull_i(np,np,nlevp) ! full nonhyrdo pressure interfaces
  type (hvcoord_t),     intent(in)  :: hvcoord             ! hybrid vertical coordinate struct
  real (kind=real_kind), intent(in) :: theta(np,np,nlev)   
  real (kind=real_kind), intent(in) :: dp3d(np,np,nlev)   
  real (kind=real_kind), intent(in) :: phi(np,np,nlev)   
  real (kind=real_kind), intent(in) :: phis(np,np)
  real (kind=real_kind), intent(in) :: Qdp(np,np,nlev)   

  !   local
  real (kind=real_kind) :: kappa_star(np,np,nlev)
  real (kind=real_kind) :: R_star(np,np,nlev)
  real (kind=real_kind) :: Qt(np,np,nlev)
  real (kind=real_kind) :: p(np,np,nlev)
  real (kind=real_kind) :: rho(np,np,nlev)
  real (kind=real_kind) :: temp_i(np,np,nlevp)  ! temp variable on level interfaces
  real (kind=real_kind) :: dphi(np,np,nlev)
  real (kind=real_kind) :: ptop
  integer :: k


  ptop = hvcoord%hyai(1)*hvcoord%ps0

#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
  do k=1,nlev
     if (use_moisture) then
        Qt(:,:,k) = Qdp(:,:,k)/dp3d(:,:,k)
        R_star(:,:,k)=(Rgas + (Rwater_vapor - Rgas)*Qt(:,:,k))
        
        if (use_cpstar==1) then
           kappa_star(:,:,k) = (Rgas + (Rwater_vapor - Rgas)*Qt(:,:,k)) / &
                (Cp + (Cpwater_vapor-Cp)*Qt(:,:,k) )
        else
           kappa_star(:,:,k) = (Rgas + (Rwater_vapor - Rgas)*Qt(:,:,k)) / Cp
        endif
     else
        R_star(:,:,k)=Rgas
        kappa_star(:,:,k)=Rgas/Cp
     endif
  enddo

! hydrostatic model:
  pfull_i(:,:,1)=ptop
  do k=1,nlev
     pfull_i(:,:,k+1)=pfull_i(:,:,k) + dp3d(:,:,k)
  enddo
  
  do k=1,nlev
     pfull(:,:,k)=pfull_i(:,:,k) + dp3d(:,:,k)/2
  enddo
  exner  = (pfull/p0)**kappa_star
  return


! compute dphi at mid points.  
#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
  do k=2,nlev
     temp_i(:,:,k) = (phi(:,:,k-1)+phi(:,:,k))/2
  enddo
  temp_i(:,:,nlev+1) = phis(:,:)

#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
  do k=2,nlev
     dphi(:,:,k) = temp_i(:,:,k+1) - temp_i(:,:,k)
  enddo


! we dont yet know phi_i(:,:1) (model top) and so cant compute dphi(:,:,1)  
!
#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
  do k=2,nlev
     rho(:,:,k) = -dp3d(:,:,k)/dphi(:,:,k)
     ! exner = (p/p0)**kappa
     ! p/exner = rho* Rstar * theta
     ! exner**( 1/kappa - 1)  = rho* Rstar * theta / p0
     ! exner = (rho* Rstar * theta / p0) ** ( kappa / 1-kappa)
     exner(:,:,k)=(rho(:,:,k)*R_star(:,:,k)*theta(:,:,k)/p0) ** &
          (kappa_star(:,:,k) / ( 1 - kappa_star(:,:,k)))
     pfull(:,:,k) = exner(:,:,k)*rho(:,:,k)*R_star(:,:,k)*theta(:,:,k)
  enddo
!
! use top of model boundary condition to compute
! pfull (and exner) at model top
!
  pfull(:,:,1) = (2*ptop + pfull(:,:,2)) / 3
! cant use p/exner = rho*R*theta since rho(1) not yet known
  exner(:,:,1) = (pfull(:,:,1)/p0)**kappa_star(:,:,1)


#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
  do k=2,nlev
     pfull_i(:,:,k) = (pfull(:,:,k-1)+pfull(:,:,k))/2
  enddo
  pfull_i(:,:,1) = ptop  
  pfull_i(:,:,nlev+1) = ptop + sum(dp3d(:,:,:),3)  


  end subroutine 








  subroutine copy_state(elem,nin,nout)
  implicit none
  
  type (element_t), intent(inout)   :: elem
  integer :: nin,nout

  elem%state%v(:,:,:,:,nout)=elem%state%v(:,:,:,:,nin)
  elem%state%w(:,:,:,nout)   =elem%state%w(:,:,:,nin)
  elem%state%theta(:,:,:,nout)   =elem%state%theta(:,:,:,nin)
  elem%state%phi(:,:,:,nout)   =elem%state%phi(:,:,:,nin)
  elem%state%dp3d(:,:,:,nout)=elem%state%dp3d(:,:,:,nin)
  elem%state%ps_v(:,:,nout)  =elem%state%ps_v(:,:,nin)
  end subroutine copy_state



  subroutine set_thermostate(elem,temperature,hvcoord,nt,ntQ)
!
! Assuming a hydrostatic intital state and given surface pressure,
! and no moisture, compute theta and phi 
!
! input:  ps_v, temperature
! ouput:  state variables:   theta, phi
!
  implicit none
  
  type (element_t), intent(inout)   :: elem
  real (kind=real_kind), intent(in) :: temperature(np,np,nlev)
  type (hvcoord_t),     intent(in)  :: hvcoord                      ! hybrid vertical coordinate struct
  
  !   local
  real (kind=real_kind) :: p(np,np,nlev)
  real (kind=real_kind) :: dp(np,np,nlev)
  real (kind=real_kind) :: dexner(np,np,nlev)
  real (kind=real_kind) :: integrand(np,np,nlev)
  real (kind=real_kind) :: pi(np,np,nlev+1)
  integer :: k,nt,ntQ

#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
  do k=1,nlev
     p(:,:,k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem%state%ps_v(:,:,nt)
     dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
          ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem%state%ps_v(:,:,nt)
  enddo

#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
  do k=1,nlev
     elem%state%theta(:,:,k,nt)=temperature(:,:,k)*(p(:,:,k)/p0)**(-kappa)
  enddo


! use dry formula for exner to initialize model:
  pi(:,:,1)=hvcoord%hyai(1)*hvcoord%ps0
  do k=1,nlev
     pi(:,:,k+1)=pi(:,:,k) + dp(:,:,k)
     dexner(:,:,k) = (pi(:,:,k+1)/p0)**kappa - (pi(:,:,k)/p0)**kappa
  enddo
  integrand(:,:,:) = Cp*elem%state%theta(:,:,:,nt)*dexner(:,:,:)
  call preq_hydrostatic_v2(elem%state%phi(:,:,:,nt),elem%state%phis,integrand)


  end subroutine set_thermostate




  subroutine set_state(u,v,T,ps,phis,p,dp,zm, g,  i,j,k,elem,n0,n1)
!
! set state variables at node(i,j,k) at layer midpoints
!
  real(real_kind),  intent(in)    :: u,v,T,ps,phis,p,dp,zm,g
  integer,          intent(in)    :: i,j,k,n0,n1
  type(element_t),  intent(inout) :: elem

  ! set prognostic state variables at level midpoints
  elem%state%v   (i,j,1,k,n0:n1) = u
  elem%state%v   (i,j,2,k,n0:n1) = v
  elem%state%dp3d(i,j,k,n0:n1)   = dp
  elem%state%ps_v(i,j,n0:n1)     = ps
  elem%state%phis(i,j)           = phis

  if (use_moisture) then
     call abortmp('ERROR: thetah set_state not yet coded for moisture')
  else
     elem%state%theta(i,j,k,n0:n1)=T/((p/p0)**kappa)
  endif
  ! set some diagnostic variables
  elem%derived%dp(i,j,k)         = dp
  elem%derived%phi(i,j,k)        = g*zm
  end subroutine set_state




  
end module

