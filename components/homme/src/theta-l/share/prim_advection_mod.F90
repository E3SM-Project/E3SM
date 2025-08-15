#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_advection_mod

  use dimensions_mod, only  : nlev, qsize, nelemd
  use kinds, only           : real_kind
  use parallel_mod, only    : parallel_t
  use derivative_mod, only  : derivative_t
  use element_mod, only     : element_t
  use hybvcoord_mod, only   : hvcoord_t
  use time_mod, only        : TimeLevel_t
  use hybrid_mod, only      : hybrid_t
  use control_mod, only     : transport_alg
  use sl_advection, only    : prim_advec_tracers_observe_velocity_ale, &
       prim_advec_tracers_remap_ALE, sl_init1
  use prim_advection_base, only: prim_advec_init2, prim_advec_init1_rk2, &
       prim_advec_tracers_remap_rk2

  implicit none
  
contains

  !  replace prim_advec_init1 and prim_advec_tracers_remap
  !
  !  prim_advec_init1:  call original eulerian init1 and the semi-lagrange init1
  !
  !  prim_advec_tracers_remap: replace with a wrapper function that can call the 
  !     eulerian code or the semi-Lagrangian code
  !
  !  TODO Is this sufficiently general that we can move this and the preqx
  !  duplicate into a base module?

  subroutine Prim_Advec_Init1(par, elem)
    type(parallel_t) :: par
    type (element_t) :: elem(:)

    call prim_advec_init1_rk2(par, elem)
    call sl_init1(par,elem)

  end subroutine Prim_Advec_Init1

  subroutine Prim_Advec_Tracers_observe_velocity(elem, tl, n, nets, nete)
    type (element_t)     , intent(inout) :: elem(:)
    type (TimeLevel_t)   , intent(in   ) :: tl
    integer              , intent(in   ) :: n       ! step in 1:dt_tracer_factor
    integer              , intent(in   ) :: nets
    integer              , intent(in   ) :: nete

    if (transport_alg /= 0) call Prim_Advec_Tracers_observe_velocity_ALE(elem, tl, n, nets, nete)
  end subroutine Prim_Advec_Tracers_observe_velocity

  subroutine Prim_Advec_Tracers_remap( elem , deriv , hvcoord ,  hybrid , dt , tl , nets , nete )
    implicit none
    type (element_t)     , intent(inout) :: elem(:)
    type (derivative_t)  , intent(in   ) :: deriv
    type (hvcoord_t)     , intent(in   ) :: hvcoord
    type (hybrid_t)      , intent(in   ) :: hybrid
    real(kind=real_kind) , intent(in   ) :: dt
    type (TimeLevel_t)   , intent(inout) :: tl
    integer              , intent(in   ) :: nets
    integer              , intent(in   ) :: nete

    if (transport_alg == 0) then
       call Prim_Advec_Tracers_remap_rk2( elem , deriv , hvcoord , hybrid , dt , tl , nets , nete )
    else
       call Prim_Advec_Tracers_remap_ALE( elem , deriv , hvcoord , hybrid , dt , tl , nets , nete )
    end if
  end subroutine Prim_Advec_Tracers_remap

end module prim_advection_mod
