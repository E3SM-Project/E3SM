#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


module precon_type_mod

   use dimensions_mod, only : nelemd, nlev
   use kinds, only : real_kind
   use element_mod, only : element_t
   use edgetype_mod, only : EdgeBuffer_t
   use hybrid_mod, only : hybrid_t
   use derivative_mod, only : derivative_t
   use time_mod, only : TimeLevel_t
   use reduction_mod, only : reductionbuffer_ordered_1d_t
   use cg_mod, only : cg_t
   use solver_mod, only : blkjac_t

   implicit none

! precon needs to be re-initialized each time dt changes 
   real (kind=real_kind) :: initialized_for_dt   = 0

  type ,public :: precon_type

    integer                 :: n
    type (element_t)      ,allocatable ,dimension(:) :: base
    type (blkjac_t)       ,allocatable ,dimension(:) :: bc
    type (EdgeBuffer_t)     :: edge1
    type (EdgeBuffer_t)     :: edge2
    type (EdgeBuffer_t)     :: edge3
    type (reductionbuffer_ordered_1d_t)              :: red
    type (derivative_t)     :: deriv
    type (cg_t)             :: cg
    real (kind=real_kind) ,allocatable ,dimension(:) :: lambdasq
    real (kind=real_kind)   :: dt
    real (kind=real_kind)   :: pmean
    type (TimeLevel_t)      :: tl
    integer                 :: nets
    integer                 :: nete

  end type precon_type

 contains


  subroutine init_precon(combo, lenx, elem, blkjac, edge1, edge2, edge3, &
         red, deriv, cg, lambdasq, dt, pmean, tl, nets, nete)

   use dimensions_mod, only : nelemd, nlev
   use kinds, only : real_kind
   use element_mod, only : element_t
   use edgetype_mod, only : EdgeBuffer_t
   use derivative_mod, only : derivative_t
   use reduction_mod, only : reductionbuffer_ordered_1d_t
   use time_mod, only : TimeLevel_t
   use cg_mod, only : cg_t
   use solver_mod, only : blkjac_t

   implicit none 

    integer :: ie, k

    integer                ,intent(in)              :: lenx
    integer                ,intent(in)              :: nets
    integer                ,intent(in)              :: nete
    type (element_t)       ,intent(inout)           :: elem(nets:nete)
    type (blkjac_t)        ,intent(inout)           :: blkjac(nets:nete)
    type (EdgeBuffer_t)    ,intent(in)              :: edge1
    type (EdgeBuffer_t)    ,intent(in)              :: edge2
    type (EdgeBuffer_t)    ,intent(in)              :: edge3
    type (reductionbuffer_ordered_1d_t) ,intent(in) :: red
    type (derivative_t)    ,intent(in)              :: deriv
    type (cg_t)            ,intent(in)              :: cg
    real (kind=real_kind)  ,intent(inout)           :: lambdasq(nlev)
    real (kind=real_kind)  ,intent(in)              :: dt
    real (kind=real_kind)  ,intent(in)              :: pmean
    type (TimeLevel_t)     ,intent(in)              :: tl
    type(precon_type)      ,intent(out)             :: combo

    allocate(combo%base(nets:nete))
    allocate(combo%bc(nets:nete))
    allocate(combo%lambdasq(nlev))

    if ( dt /= initialized_for_dt ) then
     if(cg%hybrid%par%masterproc) print *,'Initializing SI matricies for dt=',dt
       lambdasq(:) = pmean*dt*dt
       initialized_for_dt = dt
    endif

    combo%n = lenx
     do ie=nets,nete
    combo%base(ie)   = elem(ie)
    combo%bc(ie)     = blkjac(ie)
     end do
    combo%edge1      = edge1
    combo%edge2      = edge2
    combo%edge3      = edge3
    combo%red        = red
    combo%deriv      = deriv
    combo%cg         = cg
     do k=1,nlev
    combo%lambdasq(k)= lambdasq(k)
     end do
    combo%dt         = dt
    combo%pmean      = pmean
    combo%tl         = tl
    combo%nets       = nets
    combo%nete       = nete

  end subroutine init_precon

end module precon_type_mod
