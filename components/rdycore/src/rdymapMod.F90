module RDyMapMod

#include <petsc/finclude/petsc.h>
#include <petsc/finclude/petscvec.h>
  use petsc
  use RtmSpmd      , only : mpicom_rof, masterproc
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_sys_mod  , only : shr_sys_flush


  type, public :: rdy_map_type
     character(len=PETSC_MAX_PATH_LEN) :: filename

     PetscInt          :: s_ncells_loc      ! number of local source grid cells
     PetscInt, pointer :: s_ids_loc_nidx(:) ! natural IDs of source grid cells

     PetscInt          :: d_ncells_loc      ! number of destination source grid cells
     PetscInt, pointer :: d_ids_loc_nidx(:) ! natural IDs of destination grid cells

     Vec               :: s_vec             ! Vec on source grid
     Vec               :: d_vec             ! Vec on destrination grid
     VecScatter        :: s2d_scatter       ! VecScatter to transferring data from source grid to destination grid

  end type rdy_map_type

  public :: RDyMapCreate

contains

  ! ************************************************************************** !
  function RDyMapCreate()
    !
    ! This function creates a rdymapping_type
    !
    implicit none

    type(rdy_map_type), pointer :: RDyMapCreate
    type(rdy_map_type), pointer :: map

    allocate(map)

    map%s_vec       = PETSC_NULL_VEC
    map%d_vec       = PETSC_NULL_VEC
    map%s2d_scatter = PETSC_NULL_VEC_SCATTER

    RDyMapCreate => map

  end function RDyMapCreate

end module RDyMapMod
