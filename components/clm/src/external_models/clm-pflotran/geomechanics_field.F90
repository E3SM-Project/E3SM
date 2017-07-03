module Geomechanics_Field_module

  use PFLOTRAN_Constants_module
! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
  implicit none

  private

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  type, public :: geomech_field_type
    Vec :: work
    Vec :: work_loc
    ! residual vectors
    Vec :: disp_r
    Vec :: press          ! store pressure from subsurf
    Vec :: press_loc
    Vec :: press_init_loc ! store initial pressure
    Vec :: temp           ! store temperature from subsurf
    Vec :: temp_loc
    Vec :: temp_init_loc  ! store initial temperature
    Vec :: subsurf_vec_1dof ! MPI
    Vec :: imech_loc
    Vec :: strain
    Vec :: strain_loc
    Vec :: stress
    Vec :: stress_loc
    Vec :: strain_subsurf  ! Stores strains after scattering from geomech to subsurf
    Vec :: stress_subsurf  ! Stores stresses after scattering from geomech to subsurf
    Vec :: strain_subsurf_loc
    Vec :: stress_subsurf_loc
    
    Vec :: porosity_init_loc
    
    ! Solution vectors (xx = current iterate)
    Vec :: disp_xx, disp_xx_loc
    Vec :: disp_xx_init_loc
  end type geomech_field_type

  public :: GeomechFieldCreate, &
            GeomechFieldDestroy

contains

! ************************************************************************** !

function GeomechFieldCreate()
  ! 
  ! Allocates and initializes a new geomechanics
  ! Field object
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/05/13
  ! 

  implicit none
  
  type(geomech_field_type), pointer :: GeomechFieldCreate
  type(geomech_field_type), pointer :: geomech_field
  
  allocate(geomech_field)

  ! nullify PetscVecs
  geomech_field%work = 0
  geomech_field%work_loc = 0
  
  geomech_field%disp_r = 0
  geomech_field%disp_xx = 0
  geomech_field%disp_xx_loc = 0
  geomech_field%disp_xx_init_loc = 0
  
  geomech_field%press = 0
  geomech_field%press_loc = 0
  geomech_field%press_init_loc = 0
  geomech_field%temp = 0
  geomech_field%temp_loc = 0
  geomech_field%temp_init_loc = 0
  geomech_field%subsurf_vec_1dof = 0
  geomech_field%imech_loc = 0

  geomech_field%strain = 0
  geomech_field%strain_loc = 0
  geomech_field%stress = 0
  geomech_field%stress_loc = 0 

  geomech_field%strain_subsurf = 0
  geomech_field%stress_subsurf = 0
  geomech_field%strain_subsurf_loc = 0
  geomech_field%stress_subsurf_loc = 0
  
  geomech_field%porosity_init_loc = 0

  GeomechFieldCreate => geomech_field

end function GeomechFieldCreate

! ************************************************************************** !

subroutine GeomechFieldDestroy(geomech_field)
  ! 
  ! Deallocates a geomechanics field object
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/05/13
  ! 

  implicit none
  
  type(geomech_field_type), pointer :: geomech_field
  PetscErrorCode :: ierr
  
  if (.not.associated(geomech_field)) return
  
  ! Destroy PetscVecs
  if (geomech_field%work /= 0) then
    call VecDestroy(geomech_field%work,ierr);CHKERRQ(ierr)
  endif
  if (geomech_field%work_loc  /= 0) then
    call VecDestroy(geomech_field%work_loc,ierr);CHKERRQ(ierr)
  endif

  if (geomech_field%disp_r /= 0) then
    call VecDestroy(geomech_field%disp_r,ierr);CHKERRQ(ierr)
  endif
  if (geomech_field%disp_xx /= 0) then
    call VecDestroy(geomech_field%disp_xx,ierr);CHKERRQ(ierr)
  endif
  if (geomech_field%disp_xx_loc /= 0) then
    call VecDestroy(geomech_field%disp_xx_loc,ierr);CHKERRQ(ierr)
  endif
  if (geomech_field%disp_xx_init_loc /= 0) then
    call VecDestroy(geomech_field%disp_xx_init_loc,ierr);CHKERRQ(ierr)
  endif
  
  if (geomech_field%press /= 0) then
    call VecDestroy(geomech_field%press,ierr);CHKERRQ(ierr)
  endif
  if (geomech_field%press_loc /= 0) then
    call VecDestroy(geomech_field%press_loc,ierr);CHKERRQ(ierr)
  endif
  if (geomech_field%press_init_loc /= 0) then
    call VecDestroy(geomech_field%press_init_loc,ierr);CHKERRQ(ierr)
  endif
  if (geomech_field%temp /= 0) then
    call VecDestroy(geomech_field%temp,ierr);CHKERRQ(ierr)
  endif
  if (geomech_field%temp_loc /= 0) then
    call VecDestroy(geomech_field%temp_loc,ierr);CHKERRQ(ierr)
  endif
  if (geomech_field%temp_init_loc /= 0) then
    call VecDestroy(geomech_field%temp_init_loc,ierr);CHKERRQ(ierr)
  endif

  if (geomech_field%subsurf_vec_1dof /= 0 ) then
    call VecDestroy(geomech_field%subsurf_vec_1dof,ierr);CHKERRQ(ierr)
  endif
  if (geomech_field%imech_loc /= 0) then
    call VecDestroy(geomech_field%imech_loc,ierr);CHKERRQ(ierr)
  endif

  if (geomech_field%strain /= 0) then
    call VecDestroy(geomech_field%strain,ierr);CHKERRQ(ierr)
  endif
  if (geomech_field%strain_loc /= 0) then
    call VecDestroy(geomech_field%strain_loc,ierr);CHKERRQ(ierr)
  endif
  if (geomech_field%stress /= 0) then
    call VecDestroy(geomech_field%stress,ierr);CHKERRQ(ierr)
  endif
  if (geomech_field%stress_loc /= 0) then
    call VecDestroy(geomech_field%stress_loc,ierr);CHKERRQ(ierr)
  endif

  if (geomech_field%strain_subsurf /= 0) then
    call VecDestroy(geomech_field%strain_subsurf,ierr);CHKERRQ(ierr)
  endif
  if (geomech_field%stress_subsurf /= 0) then
    call VecDestroy(geomech_field%stress_subsurf,ierr);CHKERRQ(ierr)
  endif
  if (geomech_field%strain_subsurf_loc /= 0) then
    call VecDestroy(geomech_field%strain_subsurf_loc,ierr);CHKERRQ(ierr)
  endif
  if (geomech_field%stress_subsurf_loc /= 0) then
    call VecDestroy(geomech_field%stress_subsurf_loc,ierr);CHKERRQ(ierr)
  endif

  if (geomech_field%porosity_init_loc /= 0) then
    call VecDestroy(geomech_field%porosity_init_loc,ierr);CHKERRQ(ierr)
  endif

  if (associated(geomech_field)) deallocate(geomech_field)
  nullify(geomech_field)

end subroutine GeomechFieldDestroy

end module Geomechanics_Field_module
