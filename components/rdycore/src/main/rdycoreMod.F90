module rdycoreMod

#include <petsc/finclude/petsc.h>

  use petsc
  !use rdycore
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_sys_mod  , only : shr_sys_flush

  implicit none

  private

 !type(RDy)                               :: rdy_                      ! RDycore data structure

  PetscInt              , public          :: num_cells_owned           ! number of cells that locally owned
  PetscInt              , public          :: num_cells_global          ! total number of cells in the mesh
  PetscInt              , pointer         :: natural_id_cells_owned(:) ! natural IDs of cells that are locally owned

  integer               , public, pointer :: rdycore_pocn(:)           ! PE rank for each grid cell

  PetscReal             , pointer         :: total_runoff_data(:)      ! the water source to RDycore's SWE

  integer               , public          :: iulog = 6

end module rdycoreMod
