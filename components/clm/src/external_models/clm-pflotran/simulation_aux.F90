module Simulation_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  type,public :: simulation_aux_type

    ! Note: These are GLOBAL vectors (i.e. they do not contain ghost control
    !       volumes)

    ! Size of entire subsurface domain
    Vec :: subsurf_pres
    Vec :: subsurf_temp
    Vec :: subsurf_sat
    Vec :: subsurf_den

    Vec :: subsurf_por0
    Vec :: subsurf_por
    Vec :: subsurf_strain
    Vec :: subsurf_stress

    ! Size of surface cells of subsurface domain
    Vec :: subsurf_pres_top_bc
    Vec :: subsurf_temp_top_bc
    Vec :: subsurf_mflux_exchange_with_surf
    Vec :: subsurf_hflux_exchange_with_surf


    ! Size of entire surface domain
    Vec :: surf_head
    Vec :: surf_temp
    Vec :: surf_mflux_exchange_with_subsurf
    Vec :: surf_hflux_exchange_with_subsurf

    VecScatter :: surf_to_subsurf
    VecScatter :: subsurf_to_surf
    VecScatter :: subsurf_to_hydrogeophyics
    VecScatter :: geomechanics_to_subsurf
    VecScatter :: subsurf_to_geomechanics

  end type simulation_aux_type

  public :: SimAuxCreate, &
            SimAuxCopyVecScatter, &
            SimAuxCopySubsurfVec, &
            SimAuxCopySubsurfTopBCVec, &
            SimAuxCopySurfVec, &
            SimAuxCopySubsurfGeomechVec, &
            SimAuxDestroy

contains

! ************************************************************************** !

function SimAuxCreate()
  ! 
  ! This routine allocates auxillary object.
  ! 
  ! Author: Gautam Bisht,LBNL
  ! Date: 08/20/13
  ! 

  use Option_module

  implicit none

  type (simulation_aux_type),pointer :: SimAuxCreate

  type (simulation_aux_type),pointer :: aux

  allocate(aux)
  aux%subsurf_pres = 0
  aux%subsurf_temp = 0
  aux%subsurf_sat = 0
  aux%subsurf_den = 0
  aux%subsurf_por0 = 0
  aux%subsurf_por = 0
  aux%subsurf_strain = 0
  aux%subsurf_stress = 0
  aux%subsurf_pres_top_bc = 0
  aux%subsurf_temp_top_bc = 0
  aux%subsurf_mflux_exchange_with_surf = 0
  aux%subsurf_hflux_exchange_with_surf = 0

  aux%surf_head = 0
  aux%surf_temp = 0
  aux%surf_mflux_exchange_with_subsurf = 0
  aux%surf_hflux_exchange_with_subsurf = 0

  aux%surf_to_subsurf = 0
  aux%subsurf_to_surf = 0
  aux%subsurf_to_hydrogeophyics = 0
  aux%subsurf_to_geomechanics = 0
  aux%geomechanics_to_subsurf = 0

  SimAuxCreate => aux

end function SimAuxCreate

! ************************************************************************** !

subroutine SimAuxCopyVecScatter(aux, vscat, vscat_index)
  ! 
  ! This routine copies VectorScatter to an appropriate context.
  ! 
  ! Author: Gautam Bisht,LBNL
  ! Date: 10/02/13
  ! 

  implicit none

  type (simulation_aux_type),pointer :: aux
  VecScatter :: vscat
  PetscInt :: vscat_index

  PetscErrorCode :: ierr

  select case (vscat_index)
    case(SURF_TO_SUBSURF)
      call VecScatterCopy(vscat, aux%surf_to_subsurf, ierr);CHKERRQ(ierr)
    case(SUBSURF_TO_SURF)
      call VecScatterCopy(vscat, aux%subsurf_to_surf, ierr);CHKERRQ(ierr)
    case(SUBSURF_TO_HYDROGEOPHY)
      call VecScatterCopy(vscat, aux%subsurf_to_hydrogeophyics,  &
                          ierr);CHKERRQ(ierr)
    case(SUBSURF_TO_GEOMECHANICS)
      call VecScatterCopy(vscat, aux%subsurf_to_geomechanics,  &
                          ierr);CHKERRQ(ierr)
    case(GEOMECHANICS_TO_SUBSURF)
      call VecScatterCopy(vscat, aux%geomechanics_to_subsurf,  &
                          ierr);CHKERRQ(ierr)
  end select  

end subroutine SimAuxCopyVecScatter

! ************************************************************************** !

subroutine SimAuxCopySubsurfVec(aux, subsurf_vec)
  ! 
  ! This routine creates 3D vectors related with subsurface-flow.
  ! 
  ! Author: Gautam Bisht,LBNL
  ! Date: 10/02/13
  ! 

  implicit none

  type (simulation_aux_type),pointer :: aux
  Vec :: subsurf_vec

  PetscErrorCode :: ierr

  call VecDuplicate(subsurf_vec,aux%subsurf_pres,ierr);CHKERRQ(ierr)
  call VecDuplicate(subsurf_vec,aux%subsurf_temp,ierr);CHKERRQ(ierr)
  call VecDuplicate(subsurf_vec,aux%subsurf_sat,ierr);CHKERRQ(ierr)
  call VecDuplicate(subsurf_vec,aux%subsurf_den,ierr);CHKERRQ(ierr)
  call VecDuplicate(subsurf_vec,aux%subsurf_por0,ierr);CHKERRQ(ierr)
  call VecDuplicate(subsurf_vec,aux%subsurf_por,ierr);CHKERRQ(ierr)

end subroutine SimAuxCopySubsurfVec

! ************************************************************************** !

subroutine SimAuxCopySubsurfTopBCVec(aux, subsurf_top_bc_vec)
  ! 
  ! This routine creates vectors associated with surface of subsurface domain
  ! related with subsurface-flow.
  ! 
  ! Author: Gautam Bisht,LBNL
  ! Date: 10/02/13
  ! 

  implicit none

  type (simulation_aux_type),pointer :: aux
  Vec :: subsurf_top_bc_vec

  PetscErrorCode :: ierr

  call VecDuplicate(subsurf_top_bc_vec,aux%subsurf_pres_top_bc, &
                    ierr);CHKERRQ(ierr)
  call VecDuplicate(subsurf_top_bc_vec,aux%subsurf_temp_top_bc, &
                    ierr);CHKERRQ(ierr)
  call VecDuplicate(subsurf_top_bc_vec,aux%subsurf_mflux_exchange_with_surf, &
                    ierr);CHKERRQ(ierr)
  call VecDuplicate(subsurf_top_bc_vec,aux%subsurf_hflux_exchange_with_surf, &
                    ierr);CHKERRQ(ierr)

end subroutine SimAuxCopySubsurfTopBCVec

! ************************************************************************** !

subroutine SimAuxCopySurfVec(aux, surf_head_vec)
  ! 
  ! This routine creates vectors associated with surface-flow.
  ! 
  ! Author: Gautam Bisht,LBNL
  ! Date: 10/02/13
  ! 

  implicit none

  type (simulation_aux_type),pointer :: aux
  Vec :: surf_head_vec

  PetscErrorCode :: ierr

  call VecDuplicate(surf_head_vec,aux%surf_head,ierr);CHKERRQ(ierr)
  call VecDuplicate(surf_head_vec,aux%surf_temp,ierr);CHKERRQ(ierr)
  call VecDuplicate(surf_head_vec,aux%surf_mflux_exchange_with_subsurf, &
                    ierr);CHKERRQ(ierr)
  call VecDuplicate(surf_head_vec,aux%surf_hflux_exchange_with_subsurf, &
                    ierr);CHKERRQ(ierr)

end subroutine SimAuxCopySurfVec

! ************************************************************************** !

subroutine SimAuxCopySubsurfGeomechVec(aux, subsurf_geomech_vec)
  ! 
  ! This routine creates vectors associated with geomechanics.
  ! 
  ! Author: Gautam Bisht,LBNL
  ! Date: 10/02/13
  ! 
  implicit none

  type (simulation_aux_type),pointer :: aux
  Vec :: subsurf_geomech_vec

  PetscErrorCode :: ierr

  call VecDuplicate(subsurf_geomech_vec, aux%subsurf_stress,  &
                    ierr);CHKERRQ(ierr)
  call VecDuplicate(subsurf_geomech_vec, aux%subsurf_strain,  &
                    ierr);CHKERRQ(ierr)

end subroutine SimAuxCopySubsurfGeomechVec

! ************************************************************************** !

subroutine SimAuxDestroy(aux)
  ! 
  ! This routine deallocates auxillary object.
  ! 
  ! Author: Gautam Bisht,LBNL
  ! Date: 08/20/13
  ! 

  implicit none

  type(simulation_aux_type), pointer :: aux

  PetscErrorCode :: ierr

  if (.not.associated(aux)) return

  if (aux%subsurf_pres /= 0) then
    call VecDestroy(aux%subsurf_pres,ierr);CHKERRQ(ierr)
  endif
  if (aux%subsurf_temp /= 0) then
    call VecDestroy(aux%subsurf_temp,ierr);CHKERRQ(ierr)
  endif
  if (aux%subsurf_sat /= 0) then
    call VecDestroy(aux%subsurf_sat,ierr);CHKERRQ(ierr)
  endif
  if (aux%subsurf_den /= 0) then
    call VecDestroy(aux%subsurf_den,ierr);CHKERRQ(ierr)
  endif
  if (aux%subsurf_por0 /= 0) then
    call VecDestroy(aux%subsurf_por0,ierr);CHKERRQ(ierr)
  endif
  if (aux%subsurf_por /= 0) then
    call VecDestroy(aux%subsurf_por,ierr);CHKERRQ(ierr)
  endif
  if (aux%subsurf_stress /= 0) then
    call VecDestroy(aux%subsurf_stress,ierr);CHKERRQ(ierr)
  endif
  if (aux%subsurf_strain /= 0) then
    call VecDestroy(aux%subsurf_strain,ierr);CHKERRQ(ierr)
  endif

  if (aux%subsurf_pres_top_bc /= 0) then
    call VecDestroy(aux%subsurf_pres_top_bc,ierr);CHKERRQ(ierr)
  endif
  if (aux%subsurf_temp_top_bc /= 0) then
    call VecDestroy(aux%subsurf_temp_top_bc,ierr);CHKERRQ(ierr)
  endif
  if (aux%subsurf_mflux_exchange_with_surf /= 0) then
    call VecDestroy(aux%subsurf_mflux_exchange_with_surf,ierr);CHKERRQ(ierr)
  endif
  if (aux%subsurf_hflux_exchange_with_surf /= 0) then
    call VecDestroy(aux%subsurf_hflux_exchange_with_surf,ierr);CHKERRQ(ierr)
  endif

  if (aux%surf_head /= 0) then
    call VecDestroy(aux%surf_head,ierr);CHKERRQ(ierr)
  endif
  if (aux%surf_temp /= 0) then
    call VecDestroy(aux%surf_temp,ierr);CHKERRQ(ierr)
  endif
  if (aux%surf_mflux_exchange_with_subsurf /= 0) then
    call VecDestroy(aux%surf_mflux_exchange_with_subsurf,ierr);CHKERRQ(ierr)
  endif
  if (aux%surf_hflux_exchange_with_subsurf /= 0) then
    call VecDestroy(aux%surf_hflux_exchange_with_subsurf,ierr);CHKERRQ(ierr)
  endif

  if (aux%surf_to_subsurf /= 0) then
    call VecScatterDestroy(aux%surf_to_subsurf,ierr);CHKERRQ(ierr)
  endif
  if (aux%subsurf_to_surf /= 0) then
    call VecScatterDestroy(aux%subsurf_to_surf,ierr);CHKERRQ(ierr)
  endif
  if (aux%subsurf_to_hydrogeophyics /= 0) then
    call VecScatterDestroy(aux%subsurf_to_hydrogeophyics,ierr);CHKERRQ(ierr)
  endif
  if (aux%subsurf_to_geomechanics /= 0) then
    call VecScatterDestroy(aux%subsurf_to_geomechanics, ierr);CHKERRQ(ierr)
  endif
  if (aux%geomechanics_to_subsurf /= 0) then
    call VecScatterDestroy(aux%geomechanics_to_subsurf, ierr);CHKERRQ(ierr)
  endif

  deallocate(aux)
  nullify(aux)

end subroutine SimAuxDestroy

end module Simulation_Aux_module
