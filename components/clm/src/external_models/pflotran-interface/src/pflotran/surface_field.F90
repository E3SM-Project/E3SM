module Surface_Field_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
#include "petsc/finclude/petscvec.h"
  use petscvec
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: surface_field_type

    Vec :: mannings0, mannings_loc

    Vec :: work, work_loc

    Vec :: area
    
    Vec :: press_subsurf         ! MPI
    Vec :: temp_subsurf          ! MPI

    Vec :: subsurf_temp_vec_1dof ! MPI
    Vec :: subsurf_temp_vec_ndof ! MPI
    Vec :: subsurf_avg_vdarcy    ! MPI +ve value => Flow from surface to subsurface

    ! residual vectors
    Vec :: flow_r

    ! Solution vectors (yy = previous solution, xx = current iterate)
    Vec :: flow_xx, flow_xx_loc, flow_dxx, flow_yy, flow_accum

    ! vectors to save temporally average quantities
    Vec, pointer :: avg_vars_vec(:)
    PetscInt :: nvars

    ! vectors to save temporally average flowrates
    Vec :: flowrate_inst
    Vec :: flowrate_aveg

  end type surface_field_type

  public :: SurfaceFieldCreate, &
            SurfaceFieldDestroy

contains

! ************************************************************************** !

function SurfaceFieldCreate()
  ! 
  ! Allocates and initializes a new surface Field object
  ! 
  ! Author: Gautam Bisht
  ! Date: 01/17/2012
  ! 

  implicit none
  
  type(surface_field_type), pointer :: SurfaceFieldCreate
  
  type(surface_field_type), pointer :: surface_field
  
  allocate(surface_field)

  ! nullify PetscVecs
  surface_field%mannings0 = PETSC_NULL_VEC
  surface_field%mannings_loc = PETSC_NULL_VEC

  surface_field%work = PETSC_NULL_VEC
  surface_field%work_loc = PETSC_NULL_VEC

  surface_field%area = PETSC_NULL_VEC
  
  surface_field%flow_r = PETSC_NULL_VEC
  surface_field%flow_xx = PETSC_NULL_VEC
  surface_field%flow_xx_loc = PETSC_NULL_VEC
  surface_field%flow_dxx = PETSC_NULL_VEC
  surface_field%flow_yy = PETSC_NULL_VEC
  surface_field%flow_accum = PETSC_NULL_VEC
  
  surface_field%press_subsurf = PETSC_NULL_VEC

  surface_field%subsurf_temp_vec_1dof = PETSC_NULL_VEC
  surface_field%subsurf_temp_vec_ndof = PETSC_NULL_VEC

  nullify(surface_field%avg_vars_vec)
  surface_field%nvars = 0

  surface_field%flowrate_inst = PETSC_NULL_VEC
  surface_field%flowrate_aveg = PETSC_NULL_VEC

  surface_field%temp_subsurf = PETSC_NULL_VEC

  SurfaceFieldCreate => surface_field

end function SurfaceFieldCreate

! ************************************************************************** !

subroutine SurfaceFieldDestroy(surface_field)
  ! 
  ! Deallocates a field object
  ! 
  ! Author: Gautam Bisht
  ! Date: 01/17/2012
  ! 

  implicit none
  
  type(surface_field_type), pointer :: surface_field
  
  PetscErrorCode :: ierr
  PetscInt :: ivar

  ! Destroy PetscVecs
  if (surface_field%mannings0 /= PETSC_NULL_VEC) then
    call VecDestroy(surface_field%mannings0,ierr);CHKERRQ(ierr)
  endif
  if (surface_field%mannings_loc /= PETSC_NULL_VEC) then
    call VecDestroy(surface_field%mannings_loc,ierr);CHKERRQ(ierr)
  endif

  if (surface_field%work /= PETSC_NULL_VEC) then
    call VecDestroy(surface_field%work,ierr);CHKERRQ(ierr)
  endif
  if (surface_field%work_loc  /= PETSC_NULL_VEC) then
    call VecDestroy(surface_field%work_loc,ierr);CHKERRQ(ierr)
  endif

  if (surface_field%area  /= PETSC_NULL_VEC) then
    call VecDestroy(surface_field%area,ierr);CHKERRQ(ierr)
  endif
  
  if (surface_field%press_subsurf /= PETSC_NULL_VEC) then
    call VecDestroy(surface_field%press_subsurf,ierr);CHKERRQ(ierr)
  endif

  if (surface_field%flow_r /= PETSC_NULL_VEC) then
    call VecDestroy(surface_field%flow_r,ierr);CHKERRQ(ierr)
  endif
  if (surface_field%flow_xx /= PETSC_NULL_VEC) then
    call VecDestroy(surface_field%flow_xx,ierr);CHKERRQ(ierr)
  endif
  if (surface_field%flow_xx_loc /= PETSC_NULL_VEC) then
    call VecDestroy(surface_field%flow_xx_loc,ierr);CHKERRQ(ierr)
  endif
  if (surface_field%flow_dxx /= PETSC_NULL_VEC) then
    call VecDestroy(surface_field%flow_dxx,ierr);CHKERRQ(ierr)
  endif
  if (surface_field%flow_yy /= PETSC_NULL_VEC) then
    call VecDestroy(surface_field%flow_yy,ierr);CHKERRQ(ierr)
  endif
  if (surface_field%flow_accum /= PETSC_NULL_VEC) then
    call VecDestroy(surface_field%flow_accum,ierr);CHKERRQ(ierr)
  endif
  
  if (surface_field%subsurf_temp_vec_1dof/= PETSC_NULL_VEC) then
    call VecDestroy(surface_field%subsurf_temp_vec_1dof,ierr);CHKERRQ(ierr)
  endif
  if (surface_field%subsurf_temp_vec_ndof/= PETSC_NULL_VEC) then
    call VecDestroy(surface_field%subsurf_temp_vec_ndof,ierr);CHKERRQ(ierr)
  endif

  do ivar = 1,surface_field%nvars
    call VecDestroy(surface_field%avg_vars_vec(ivar),ierr);CHKERRQ(ierr)
  enddo

  if (surface_field%flowrate_inst/= PETSC_NULL_VEC) then
    call VecDestroy(surface_field%flowrate_inst,ierr);CHKERRQ(ierr)
  endif
  if (surface_field%flowrate_aveg/= PETSC_NULL_VEC) then
    call VecDestroy(surface_field%flowrate_aveg,ierr);CHKERRQ(ierr)
  endif

  if (surface_field%temp_subsurf /= PETSC_NULL_VEC ) then
    call VecDestroy(surface_field%temp_subsurf,ierr);CHKERRQ(ierr)
  endif

  if (associated(surface_field)) deallocate(surface_field)
  nullify(surface_field)

end subroutine SurfaceFieldDestroy

end module Surface_Field_module
