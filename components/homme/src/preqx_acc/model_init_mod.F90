#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!
!  model specific initialization code, called from prim_init2
!
module model_init_mod
  use element_mod,    only: element_t
  use hybrid_mod,         only: hybrid_t
  use derivative_mod, only: derivative_t
  use time_mod,       only: timelevel_t
  use hybvcoord_mod,  only: hvcoord_t

  implicit none

  public :: model_init2

contains


  !_____________________________________________________________________
  subroutine vertical_mesh_init2(elem, nets, nete, hybrid, hvcoord)

    ! additional solver specific initializations (called from prim_init2)

    type (element_t),			intent(inout), target :: elem(:)! array of element_t structures
    integer,				intent(in) :: nets,nete		! start and end element indices
    type (hybrid_t),			intent(in) :: hybrid		! mpi/omp data struct
    type (hvcoord_t),			intent(inout)	:: hvcoord	! hybrid vertical coord data struct

  end subroutine vertical_mesh_init2

    
  subroutine model_init2( elem , hybrid, deriv ,hvcoord,tl,nets,nete)
    use element_state, only: state_qdp, derived_vn0, derived_divdp, derived_divdp_proj
    use dimensions_mod, only: nelemd

    implicit none
    type(element_t)   , intent(in) :: elem(:)
    type(hybrid_t)    , intent(in) :: hybrid
    type(derivative_t), intent(in) :: deriv
    type (hvcoord_t)  , intent(in) :: hvcoord
    type (TimeLevel_t), intent(in) :: tl
    integer                        :: nets,nete

    integer :: ie

    !$omp barrier
    !$omp master

  
    !$acc enter data pcreate(state_Qdp,derived_vn0,derived_divdp,derived_divdp_proj)
    !$acc enter data pcopyin(elem(1:nelemd),deriv)
    do ie = 1 , nelemd
      !$acc enter data pcopyin(elem(ie)%desc)
      !$acc enter data pcopyin(elem(ie)%desc%putmapP,elem(ie)%desc%getmapP,elem(ie)%desc%reverse)
    enddo

    !$omp end master
    !$omp barrier
  end subroutine model_init2


end module model_init_mod
