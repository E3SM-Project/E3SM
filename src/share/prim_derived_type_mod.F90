#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_derived_type_mod

   use kinds, only : real_kind
   use element_mod, only : element_t
   use hybvcoord_mod, only : hvcoord_t
   use hybrid_mod, only : hybrid_t
   use derivative_mod, only : derivative_t
   use time_mod, only : TimeLevel_t

   implicit none

  type ,public :: derived_type

    integer                 :: method
    type (element_t)      ,allocatable ,dimension(:) :: base
    type (hvcoord_t)        :: hvcoord
    type (hybrid_t)         :: hybrid
    logical                 :: compute_diagnostics
    integer                 :: n_Q
    real (kind=real_kind)   :: eta_ave_w
    type (derivative_t)     :: deriv
    real (kind=real_kind)   :: dt
    type (TimeLevel_t)      :: tl
    integer                 :: nets
    integer                 :: nete

  end type derived_type

!  type ,public :: precon_type
!     sequence
!
!    integer                 :: n
!    type (element_t)      ,allocatable ,dimension(:) :: base
!    type (blkjac_t)       ,allocatable ,dimension(:) :: bc
!    type (reductionbuffer_ordered_1d_t)              :: red
!    type (ref_state_t)      :: refstate
!    type (hvcoord_t)        :: hvcoord
!    type (hybrid_t)         :: hybrid
!    type (derivative_t)     :: deriv
!    type (filter_t)         :: flt
!    type (cg_t)             :: cg
!    real (kind=real_kind)   :: dt
!    real (kind=real_kind)   :: pmean
!    integer                 :: ntl1
!    integer                 :: ntl2
!    integer                 :: ntl3
!    integer                 :: nets
!    integer                 :: nete
!
!  end type precon_type

 contains
  subroutine initialize(object, method, elem, hvcoord, compute_diagnostics, &
           n_Q, eta_ave_w, hybrid, deriv, dt, tl, nets, nete)

   use kinds, only : real_kind
   use element_mod, only : element_t
   use hybvcoord_mod, only : hvcoord_t
   use hybrid_mod, only : hybrid_t
   use derivative_mod, only : derivative_t
   use time_mod, only : TimeLevel_t

   implicit none 

    integer :: ie

    integer                ,intent(in)  :: method
    integer                ,intent(in)  :: nets
    integer                ,intent(in)  :: nete
    type (element_t)       ,intent(in)  :: elem(nets:nete)
    type (hvcoord_t)       ,intent(in)  :: hvcoord
    type (hybrid_t)        ,intent(in)  :: hybrid
    logical                ,intent(in)  :: compute_diagnostics
    integer                ,intent(in)  :: n_Q
    real (kind=real_kind)  ,intent(in)  :: eta_ave_w
    type (derivative_t)    ,intent(in)  :: deriv
    real (kind=real_kind)  ,intent(in)  :: dt
    type (TimeLevel_t)     ,intent(in)  :: tl
    type(derived_type)     ,intent(out) :: object

   allocate(object%base(nets:nete))
   object%method = method
    do ie=nets,nete
   object%base(ie)= elem(ie)
    end do
   object%hvcoord = hvcoord
   object%hybrid = hybrid
   object%compute_diagnostics = compute_diagnostics
   object%n_Q = n_Q
   object%eta_ave_w = eta_ave_w
   object%deriv = deriv
   object%dt = dt
   object%tl = tl
   object%nets = nets
   object%nete = nete

  end subroutine initialize

!  subroutine init_precon(combo, lenx, elem, blkjac, refstate, hvcoord, flt, hybrid, &
!         red, deriv, cg, dt, ntl1, ntl2, ntl3, nets, nete)
!
!   use dimensions_mod, only : nelemd, nlev
!   use kinds, only : real_kind
!   use element_mod, only : element_t
!   use derivative_mod, only : derivative_t
!   use reduction_mod, only : reductionbuffer_ordered_1d_t
!   use cg_mod, only : cg_t
!   use solver_mod, only : blkjac_t
!   use filter_mod, only : filter_t
!   use hybvcoord_mod, only : hvcoord_t
!   use hybrid_mod, only : hybrid_t
!   use prim_si_ref_mod, only : ref_state_t

!   implicit none 

!    integer :: ie, k

!    integer                ,intent(in)              :: lenx
!    type (element_t)       ,intent(inout)           :: elem(nets:nete)
!    type (blkjac_t)        ,intent(inout)           :: blkjac(nets:nete)
!    type (ref_state_t)     ,intent(inout)           :: refstate
!    type (hvcoord_t)       ,intent(inout)           :: hvcoord
!    type (filter_t)        ,intent(inout)           :: flt
!    type (hybrid_t)        ,intent(inout)           :: hybrid
!    type (reductionbuffer_ordered_1d_t) ,intent(in) :: red
!    type (derivative_t)    ,intent(in)              :: deriv
!    type (cg_t)            ,intent(in)              :: cg
!    real (kind=real_kind)  ,intent(in)              :: dt
!    integer                ,intent(in)              :: ntl1
!    integer                ,intent(in)              :: ntl2
!    integer                ,intent(in)              :: ntl3
!    integer                ,intent(in)              :: nets
!    integer                ,intent(in)              :: nete
!    type(precon_type)      ,intent(out)             :: combo
!
!    allocate(combo%base(nets:nete))
!    allocate(combo%bc(nets:nete))
!
!    combo%n = lenx
!     do ie=nets,nete
!    combo%base(ie)   = elem(ie)
!    combo%bc(ie)     = blkjac(ie)
!     end do
!    combo%refstate   = refstate
!    combo%hvcoord    = hvcoord
!    combo%flt        = flt
!    combo%hybrid     = hybrid
!    combo%red        = red
!    combo%deriv      = deriv
!    combo%cg         = cg
!    combo%dt         = dt
!    combo%ntl1       = ntl1
!    combo%ntl2       = ntl2
!    combo%ntl3       = ntl3
!    combo%nets       = nets
!    combo%nete       = nete

!  end subroutine init_precon

end module prim_derived_type_mod
