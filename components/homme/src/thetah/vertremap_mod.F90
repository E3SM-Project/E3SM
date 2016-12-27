#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module vertremap_mod
  use vertremap_mod_base, only: remap1

  use kinds, only                  : real_kind,int_kind
  use dimensions_mod, only         : np,nlev,qsize,nlevp,npsq
  use hybvcoord_mod, only          : hvcoord_t
  use element_mod, only            : element_t
  use perf_mod, only               : t_startf, t_stopf  ! _EXTERNAL
  use parallel_mod, only           : abortmp, parallel_t
  use control_mod, only : vert_remap_q_alg

  implicit none
  private
  public :: vertical_remap

contains


  subroutine vertical_remap(hybrid,elem,hvcoord,dt,np1,np1_qdp,nets,nete)

  ! This routine is called at the end of the vertically Lagrangian
  ! dynamics step to compute the vertical flux needed to get back
  ! to reference eta levels
  !
  ! input:
  !     derived%dp()  delta p on levels at beginning of timestep
  !     state%dp3d(np1)  delta p on levels at end of timestep
  ! output:
  !     state%ps_v(np1)          surface pressure at time np1
  !     derived%eta_dot_dpdn()   vertical flux from final Lagrangian
  !                              levels to reference eta levels
  !
  use kinds,          only: real_kind
  use hybvcoord_mod,  only: hvcoord_t
  use control_mod,    only: rsplit
  use hybrid_mod,     only: hybrid_t


  type (hybrid_t),  intent(in)    :: hybrid  ! distributed parallel structure (shared)
  type (element_t), intent(inout) :: elem(:)
  type (hvcoord_t)                :: hvcoord
  real (kind=real_kind)           :: dt

  integer :: ie,i,j,k,np1,nets,nete,np1_qdp
  integer :: q

  real (kind=real_kind), dimension(np,np,nlev)  :: dp,dp_star
  real (kind=real_kind), dimension(np,np,nlev,5)  :: ttmp

  call t_startf('vertical_remap')

  ! reference levels:
  !   dp(k) = (hyai(k+1)-hyai(k))*ps0 + (hybi(k+1)-hybi(k))*ps_v(i,j)
  !   hybi(1)=0          pure pressure at top of atmosphere
  !   hyai(1)=ptop
  !   hyai(nlev+1) = 0   pure sigma at bottom
  !   hybi(nlev+1) = 1
  !
  ! sum over k=1,nlev
  !  sum(dp(k)) = (hyai(nlev+1)-hyai(1))*ps0 + (hybi(nlev+1)-hybi(1))*ps_v
  !             = -ps0 + ps_v
  !  ps_v =  ps0+sum(dp(k))
  !
  ! reference levels:
  !    dp(k) = (hyai(k+1)-hyai(k))*ps0 + (hybi(k+1)-hybi(k))*ps_v
  ! floating levels:
  !    dp_star(k) = dp(k) + dt_q*(eta_dot_dpdn(i,j,k+1) - eta_dot_dpdn(i,j,k) )
  ! hence:
  !    (dp_star(k)-dp(k))/dt_q = (eta_dot_dpdn(i,j,k+1) - eta_dot_dpdn(i,j,k) )
  !
   do ie=nets,nete
     ! update final ps_v
     elem(ie)%state%ps_v(:,:,np1) = hvcoord%hyai(1)*hvcoord%ps0 + &
          sum(elem(ie)%state%dp3d(:,:,:,np1),3)
     do k=1,nlev
        dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
             ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,np1)
        if (rsplit==0) then
           dp_star(:,:,k) = dp(:,:,k) + dt*(elem(ie)%derived%eta_dot_dpdn(:,:,k+1) -&
                elem(ie)%derived%eta_dot_dpdn(:,:,k))
        else
           dp_star(:,:,k) = elem(ie)%state%dp3d(:,:,k,np1)
        endif
     enddo
     if (minval(dp_star)<0) then
        do k=1,nlev
        do i=1,np
        do j=1,np
           if (dp_star(i,j,k ) < 0) then
              print *,'level = ',k
              print *,"column location lat,lon (radians):",elem(ie)%spherep(i,j)%lat,elem(ie)%spherep(i,j)%lon
           endif
        enddo
        enddo
        enddo
        call abortmp('negative layer thickness.  timestep or remap time too large')
     endif

     if (rsplit>0) then
        !  REMAP u,v,T from levels in dp3d() to REF levels
        ttmp(:,:,:,1)=elem(ie)%state%v(:,:,1,:,np1)*dp_star
        ttmp(:,:,:,2)=elem(ie)%state%v(:,:,2,:,np1)*dp_star
        ttmp(:,:,:,3)=elem(ie)%state%theta(:,:,:,np1)*dp_star
        ttmp(:,:,:,4)=elem(ie)%state%phi(:,:,:,np1)*dp_star
        ttmp(:,:,:,5)=elem(ie)%state%w(:,:,:,np1)*dp_star

        call t_startf('vertical_remap1_1')
        call remap1(ttmp,np,5,dp_star,dp)
        call t_stopf('vertical_remap1_1')

        elem(ie)%state%v(:,:,1,:,np1)=ttmp(:,:,:,1)/dp
        elem(ie)%state%v(:,:,2,:,np1)=ttmp(:,:,:,2)/dp
        elem(ie)%state%theta(:,:,:,np1)=ttmp(:,:,:,3)/dp
        elem(ie)%state%phi(:,:,:,np1)=ttmp(:,:,:,4)/dp
        elem(ie)%state%w(:,:,:,np1)=ttmp(:,:,:,5)/dp
     endif

     ! remap the gll tracers from lagrangian levels (dp_star)  to REF levels dp
     if (qsize>0) then

       call t_startf('vertical_remap1_3')
       call remap1(elem(ie)%state%Qdp(:,:,:,:,np1_qdp),np,qsize,dp_star,dp)
       call t_stopf('vertical_remap1_3')

     endif
  enddo
  call t_stopf('vertical_remap')
  end subroutine vertical_remap


end module 




