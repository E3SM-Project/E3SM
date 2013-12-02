#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


module derived_type_mod

   use dimensions_mod, only : nelemd, nlev
   use kinds, only : real_kind
   use element_mod, only : element_t
   use edge_mod, only : EdgeBuffer_t
   use hybrid_mod, only : hybrid_t
   use derivative_mod, only : derivative_t
   use time_mod, only : TimeLevel_t
   use reduction_mod, only : reductionbuffer_ordered_1d_t
   use cg_mod, only : cg_t

   implicit none

! precon needs to be re-initialized each time dt changes 
   real (kind=real_kind) :: initialized_for_dt   = 0

  type ,public :: derived_type

    integer                 :: n
    type (element_t)      ,allocatable ,dimension(:) :: base
    type (EdgeBuffer_t)     :: edge1
    type (EdgeBuffer_t)     :: edge2
    type (EdgeBuffer_t)     :: edge3
    type (hybrid_t)         :: hybrid
    type (derivative_t)     :: deriv
    real (kind=real_kind)   :: dt
    real (kind=real_kind)   :: pmean
    type (TimeLevel_t)      :: tl
    integer                 :: nets
    integer                 :: nete

  end type derived_type


 contains


  subroutine initialize(object, lenx, elem, pmean, edge1,edge2,edge3, &
         hybrid, deriv, dt, tl, nets, nete)

   use dimensions_mod, only : nelemd
   use kinds, only : real_kind
   use element_mod, only : element_t
   use edge_mod, only : EdgeBuffer_t
   use hybrid_mod, only : hybrid_t
   use derivative_mod, only : derivative_t
   use time_mod, only : TimeLevel_t

   implicit none 

    integer :: ie

    integer                ,intent(in)  :: lenx
    type (element_t)       ,intent(in)  :: elem(nets:nete)
    real (kind=real_kind)  ,intent(in)  :: pmean
    type (EdgeBuffer_t)    ,intent(in)  :: edge1
    type (EdgeBuffer_t)    ,intent(in)  :: edge2
    type (EdgeBuffer_t)    ,intent(in)  :: edge3
    type (hybrid_t)        ,intent(in)  :: hybrid
    type (derivative_t)    ,intent(in)  :: deriv
    real (kind=real_kind)  ,intent(in)  :: dt
    type (TimeLevel_t)     ,intent(in)  :: tl
    integer                ,intent(in)  :: nets
    integer                ,intent(in)  :: nete
    type(derived_type)     ,intent(out) :: object

    allocate(object%base(nets:nete))
    object%n = lenx
     do ie=nets,nete
    object%base(ie)= elem(ie)
     end do
    object%pmean = pmean
    object%edge1 = edge1
    object%edge2 = edge2
    object%edge3 = edge3
    object%hybrid = hybrid
    object%deriv = deriv
    object%dt = dt
    object%tl = tl
    object%nets = nets
    object%nete = nete

  end subroutine initialize


end module derived_type_mod
