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
  PetscInt              , public, pointer :: natural_id_cells_owned(:) ! natural IDs of cells that are locally owned

  integer               , public, pointer :: rdycore_pocn(:)           ! PE rank for each grid cell

  PetscReal             , pointer         :: total_runoff_data(:)      ! the water source to RDycore's SWE

  integer               , public          :: iulog = 6
  character(len=16)     , public          :: inst_name
  character(len=16)     , public          :: inst_suffix         ! char string associated with instance (ie. "_0001" or "")
  integer               , public          :: inst_index          ! number of current instance (ie. 1)


end module rdycoreMod
