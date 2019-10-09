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
  call abortmp('ERROR: vertical remap not supported by PESE model')
  call t_stopf('vertical_remap')
  end subroutine vertical_remap


end module 




