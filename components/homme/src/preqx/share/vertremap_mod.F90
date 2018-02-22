#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module vertremap_mod
  use vertremap_base, only: remap1

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
  real (kind=real_kind), dimension(np,np,nlev,2)  :: ttmp

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
#undef REMAP_TE
#ifdef REMAP_TE
        ! remap u,v and cp*T + .5 u^2
        ttmp(:,:,:,1)=(elem(ie)%state%v(:,:,1,:,np1)**2 + &
             elem(ie)%state%v(:,:,2,:,np1)**2)/2 + &
             elem(ie)%state%t(:,:,:,np1)*cp
#else
        ttmp(:,:,:,1)=elem(ie)%state%t(:,:,:,np1)
#endif
        ttmp(:,:,:,1)=ttmp(:,:,:,1)*dp_star

        call t_startf('vertical_remap1_1')
        call remap1(ttmp,np,1,dp_star,dp)
        call t_stopf('vertical_remap1_1')

        elem(ie)%state%t(:,:,:,np1)=ttmp(:,:,:,1)/dp

        ttmp(:,:,:,1)=elem(ie)%state%v(:,:,1,:,np1)*dp_star
        ttmp(:,:,:,2)=elem(ie)%state%v(:,:,2,:,np1)*dp_star

        call t_startf('vertical_remap1_2')
        call remap1(ttmp,np,2,dp_star,dp)
        call t_stopf('vertical_remap1_2')

        elem(ie)%state%v(:,:,1,:,np1)=ttmp(:,:,:,1)/dp
        elem(ie)%state%v(:,:,2,:,np1)=ttmp(:,:,:,2)/dp

#ifdef REMAP_TE
        ! back out T from TE
        elem(ie)%state%t(:,:,:,np1) = &
             ( elem(ie)%state%t(:,:,:,np1) - ( (elem(ie)%state%v(:,:,1,:,np1)**2 + &
             elem(ie)%state%v(:,:,2,:,np1)**2)/2))/cp
#endif
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

!routine to apply external forcing to homme on Lagrangian levels
!that is, unrelated to whether remap stage was done or not
  subroutine remap_vsplit_dyn(hybrid,elem,hvcoord,dt,np1,nets,nete)
  ! input:
  !     ps_forcing(:,:,nets,nete) is ps for levels at which forcing is obtained
  !     originally
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

  real (kind=real_kind), dimension(np,np,nlev)  :: dp_forcing,dp_star
  real (kind=real_kind), dimension(np,np,nlev,2)  :: ttmp
  real (kind=real_kind) :: ps_forcing(np,np)

  call t_startf('remap_vsplit_dyn')
   do ie=nets,nete
     ps_forcing(:,:) = hvcoord%hyai(1)*hvcoord%ps0 + sum(elem(ie)%state%dp3d(:,:,:,np1),3)
     do k=1,nlev
        !obtain dp for forcing
        dp_forcing(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
             ( hvcoord%hybi(k+1) - hvcoord%hybi(k))*ps_forcing(:,:)
        dp_star(:,:,k) = elem(ie)%state%dp3d(:,:,k,np1)
     enddo
     !since forcing was received from physics when pressure levels were
     !different from the current ones, there is a question how to treat the
     !forcing. for now we adopt the same approach as in the original code,
     !namely, assume that the forcing is for pressure levels that are based on
     !eta coordinate and *current* surface pressure. this way we do not need to
     !stretch the old pressure grid (to match sum(dp_old)=sum(dp_new) or develop
     !a special monotone and conservative mapping procedure.

     !scaling idea, abandoned
     !dp_forcing = dp_forcing*(sum(dp_star(:,:,:),3)/sum(dp_forcing(:,:,3)))

     !not sure whether to check both dp_star and dp_forcing?

    if (minval(dp_star)<0) then
        do k=1,nlev
        do i=1,np
        do j=1,np
           if (dp_star(i,j,k ) < 0) then
              print *,"In remap_vsplit_dyn: level = ",k
              print *,"In remap_vsplit_dyn: column location lat,lon(radians):",&
                      elem(ie)%spherep(i,j)%lat,elem(ie)%spherep(i,j)%lon
           endif
        enddo
        enddo
        enddo
        call abortmp('In remap_vsplit_dyn: negative layer thickness. timestep or remap time too large')
     endif

     !we don't support REMAP_TE since it is off by default

     ! remap forcing T
     ttmp(:,:,:,1)=elem(ie)%derived%FT(:,:,:)*dp_forcing

     call t_startf('vertical_remap1_1')
     !remap from dp_forcing to dp_star
     call remap1(ttmp,np,1,dp_forcing,dp_star)
     call t_stopf('vertical_remap1_1')
     !add forcing to temperature
     elem(ie)%state%t(:,:,:,np1)=elem(ie)%state%t(:,:,:,np1) + ttmp(:,:,:,1)/dp_star(:,:,:)

     !remap forcing V
     ttmp(:,:,:,1)=elem(ie)%derived%v(:,:,1,:,np1)*dp_forcing
     ttmp(:,:,:,2)=elem(ie)%derived%v(:,:,2,:,np1)*dp_forcing

     call t_startf('vertical_remap1_2')
     !remap from dp_forcing to dp_star
     call remap1(ttmp,np,2,dp_forcing,dp_star)
     call t_stopf('vertical_remap1_2')
     !add forcing back
     elem(ie)%state%v(:,:,1,:,np1)=elem(ie)%state%v(:,:,1,:,np1) + ttmp(:,:,:,1)/dp_star(:,:,:)
     elem(ie)%state%v(:,:,2,:,np1)=elem(ie)%state%v(:,:,2,:,np1) + ttmp(:,:,:,2)/dp_star(:,:,:)
  enddo
  call t_stopf('remap_vsplit_dyn')
  end subroutine remap_vsplit_dyn

end module 




