module rdy_import_export

#include <petsc/finclude/petsc.h>

  use petsc
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use rdydecompMod    , only : rdy_bounds_type
  use lnd2rdyType     , only : lnd2rdy_type
  use rdycoreMod      , only : rtm2rdy_map, rtm2rdy_nvars
  use rdy_cpl_indices

  implicit none
  private

  public :: rdy_import_mct

contains

  !------------------------------------------------------------------------
  subroutine rdy_import_mct(rdy_bounds, x2rd, ncells_rtm, lnd2rdy_vars)
    !
    ! !DESCRIPTION:
    ! Imports data from the coupler for the RDycore
    !
    implicit none
    !
    ! !ARGUMENTS:
    type(rdy_bounds_type) , intent (in)    :: rdy_bounds
    real(r8)              , intent (in)    :: x2rd(:,:)
    integer               , intent (in)    :: ncells_rtm
    type(lnd2rdy_type)    , intent (inout) :: lnd2rdy_vars
    !
    ! !LOCAL VARIABLES:
    integer              :: ii, g, idx                     ! indices
    real(r8)             :: runoff_unit_conversion = 1.d-3 ! [mm/s] to [m/s]
    PetscScalar, pointer :: v_loc(:)
    character(len=1024)  :: string
    integer              :: myrank
    PetscViewer          :: viewer
    PetscErrorCode       :: ierr

    ! pack the data
    PetscCallA(VecGetArrayF90(rtm2rdy_map%s_vec, v_loc, ierr))
    do ii = 1, ncells_rtm
       v_loc((ii-1)*rtm2rdy_nvars + 1) = x2rd(index_x2rdy_Flrl_qsur, ii) * runoff_unit_conversion
       v_loc((ii-1)*rtm2rdy_nvars + 2) = x2rd(index_x2rdy_Flrl_qsub, ii) * runoff_unit_conversion
    end do
    PetscCallA(VecRestoreArrayF90(rtm2rdy_map%s_vec, v_loc, ierr))

    ! scatter the data
    PetscCallA(VecScatterBegin(rtm2rdy_map%s2d_scatter, rtm2rdy_map%s_vec, rtm2rdy_map%d_vec, INSERT_VALUES, SCATTER_FORWARD, ierr))
    PetscCallA(VecScatterEnd(rtm2rdy_map%s2d_scatter, rtm2rdy_map%s_vec, rtm2rdy_map%d_vec, INSERT_VALUES, SCATTER_FORWARD, ierr))

#if 0
    ! find the rank for the processor
    PetscCallA(mpi_comm_rank(PETSC_COMM_WORLD, myrank, ierr))

    string = 'src_vec.txt'
    PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD, trim(string), viewer, ierr))
    PetscCallA(VecView(rtm2rdy_map%s_vec, viewer, ierr))
    PetscCallA(PetscViewerDestroy(viewer, ierr))

    write(string,*)myrank
    string = 'dst_vec_' // trim(adjustl(string)) // '.txt'
    PetscCallA(PetscViewerASCIIOpen(PETSC_COMM_WORLD, trim(string), viewer, ierr))
    PetscCallA(VecView(rtm2rdy_map%d_vec, viewer, ierr))
    PetscCallA(PetscViewerDestroy(viewer, ierr))
#endif

    ! unpack the data
    PetscCallA(VecGetArrayF90(rtm2rdy_map%d_vec, v_loc, ierr))

    do g = rdy_bounds%begg, rdy_bounds%endg
       idx = (g - rdy_bounds%begg)*rtm2rdy_nvars

       lnd2rdy_vars%forc_qsur(g) = v_loc(idx + 1)
       lnd2rdy_vars%forc_qsub(g) = v_loc(idx + 2)
    end do

    PetscCallA(VecRestoreArrayF90(rtm2rdy_map%d_vec, v_loc, ierr))

  end subroutine rdy_import_mct

end module rdy_import_export
