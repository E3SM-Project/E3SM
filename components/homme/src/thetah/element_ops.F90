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
  real (kind=real_kind) :: p(np,np,nlev)
  real (kind=real_kind) :: dp(np,np,nlev)
  real (kind=real_kind) :: kappa_star(np,np,nlev)
  real (kind=real_kind) :: Qt(np,np,nlev)
  integer :: k
  
  !   pressure on floating levels 
  !    dp  => elem%state%dp3d(:,:,:,n0)
  !    p(:,:,1)=hvcoord%hyai(1)*hvcoord%ps0 + dp(:,:,1)/2
  !    do k=2,nlev
  !       p(:,:,k)=p(:,:,k-1) + dp(:,:,k-1)/2 + dp(:,:,k)/2
  !    enddo
  
#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
  do k=1,nlev
     p(:,:,k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem%state%ps_v(:,:,nt)
     dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
          ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem%state%ps_v(:,:,nt)
  enddo

  if (use_moisture==1) then
#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
     do k=1,nlev
        Qt(:,:,k) = elem%state%Qdp(:,:,k,1,ntQ)/dp(:,:,k)
        
        if (use_cpstar==1) then
           kappa_star(:,:,k) = (Rgas + (Rwater_vapor - Rgas)*Qt(:,:,k)) / &
                (Cp + (Cpwater_vapor-Cp)*Qt(:,:,k) )
        else
           kappa_star(:,:,k) = (Rgas + (Rwater_vapor - Rgas)*Qt(:,:,k)) / Cp
        endif
        
        temperature(:,:,k)= elem%state%theta(:,:,k,nt)*(p(:,:,k)/p0)**kappa_star(:,:,k)
     enddo
  else
#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
     do k=1,nlev
        temperature(:,:,k)= elem%state%theta(:,:,k,nt)*(p(:,:,k)/p0)**kappa
     enddo
  endif

  end subroutine get_temperature


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



  subroutine set_thermostate(elem,temperature,hvcoord,n0,ntQ)
  implicit none
  
  type (element_t), intent(inout)   :: elem
  real (kind=real_kind), intent(in) :: temperature(np,np,nlev)
  type (hvcoord_t),     intent(in)  :: hvcoord                      ! hybrid vertical coordinate struct
  
  !   local
  real (kind=real_kind) :: p(np,np,nlev)
  real (kind=real_kind) :: dp(np,np,nlev)
  real (kind=real_kind) :: kappa_star(np,np,nlev)
  real (kind=real_kind) :: Qt(np,np,nlev)
  integer :: k,n0,ntQ

#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
  do k=1,nlev
     p(:,:,k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem%state%ps_v(:,:,n0)
     dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
          ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem%state%ps_v(:,:,n0)
  enddo

  if (use_moisture==1) then
#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
     do k=1,nlev
        Qt(:,:,k) = elem%state%Qdp(:,:,k,1,ntQ)/dp(:,:,k)
        
        if (use_cpstar==1) then
           kappa_star(:,:,k) = (Rgas + (Rwater_vapor - Rgas)*Qt(:,:,k)) / &
                (Cp + (Cpwater_vapor-Cp)*Qt(:,:,k) )
        else
           kappa_star(:,:,k) = (Rgas + (Rwater_vapor - Rgas)*Qt(:,:,k)) / Cp
        endif
        elem%state%theta(:,:,k,n0) = temperature(:,:,k)/((p(:,:,k)/p0)**kappa_star(:,:,k))
     enddo
  else
#if (defined COLUMN_OPENMP)
  !$omp parallel do default(shared), private(k)
#endif
     do k=1,nlev
        !temperature(:,:,k)= elem%state%theta(:,:,k,nt)*(p(:,:,k)/p0)**kappa
        elem%state%theta(:,:,k,n0)=temperature(:,:,k)/((p(:,:,k)/p0)**kappa)
     enddo
  endif

  end subroutine set_thermostate


subroutine set_state(u,v,T,ps,phis,p,dp,zm, g,  i,j,k,elem,n0,n1)

  ! set state variables at node(i,j,k) at layer midpoints

  real(real_kind),  intent(in)    :: u,v,T,ps,phis,p,dp,zm,g
  integer,          intent(in)    :: i,j,k,n0,n1
  type(element_t),  intent(inout) :: elem

  ! set prognostic state variables at level midpoints
  elem%state%v   (i,j,1,k,n0:n1) = u
  elem%state%v   (i,j,2,k,n0:n1) = v
  elem%state%dp3d(i,j,k,n0:n1)   = dp
  elem%state%ps_v(i,j,n0:n1)     = ps
  elem%state%phis(i,j)           = phis


  if (use_moisture==1) then
     call abortmp('ERROR: thetah set_state not yet coded for moisture')
  else
     elem%state%theta(i,j,k,n0:n1)=T/((p/p0)**kappa)
  endif



  ! set some diagnostic variables
  elem%derived%dp(i,j,k)         = dp
  elem%derived%phi(i,j,k)        = g*zm

  

end subroutine
  
end module

