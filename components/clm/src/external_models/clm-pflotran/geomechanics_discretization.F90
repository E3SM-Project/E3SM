module Geomechanics_Discretization_module

  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use PFLOTRAN_Constants_module
  
  implicit none

  private
 
#include "petsc/finclude/petscsys.h"

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscmat.h90"
#include "petsc/finclude/petscdm.h"
#include "petsc/finclude/petscdm.h90"
#include "petsc/finclude/petscdmda.h"
#include "petsc/finclude/petscdmshell.h90"

  type, public :: gmdm_ptr_type
    DM :: dm  ! PETSc DM
    type(gmdm_type), pointer :: gmdm
  end type gmdm_ptr_type

  type, public :: geomech_discretization_type
    PetscInt :: itype                          ! type of discretization (e.g. structured, unstructured, etc.)
    character(len=MAXWORDLENGTH) :: ctype      ! name of discretization
    PetscReal :: origin(3)                     ! origin of global domain
    type(geomech_grid_type), pointer :: grid   ! pointer to a grid object
    character(len=MAXSTRINGLENGTH) :: filename
    PetscInt :: dm_index_to_ndof(3)            ! mapping between a dm_ptr to the number of degrees of freedom
    type(gmdm_ptr_type), pointer :: dm_1dof
    type(gmdm_ptr_type), pointer :: dm_ngeodof 
    type(gmdm_ptr_type), pointer :: dm_n_stress_strain_dof    ! For stress and strain
  end type geomech_discretization_type

  public :: GeomechDiscretizationCreate, &
            GeomechDiscretizationDestroy, &
            GeomechDiscretizationCreateVector, &
            GeomechDiscretizationDuplicateVector, &         
            GeomechDiscretizationCreateJacobian, &
            GeomechDiscretizationGlobalToLocal, &
            GeomechDiscretizationLocalToGlobal, &
            GeomechDiscretizationLocalToGlobalAdd, &
            GeomechDiscretizationLocalToLocal, &
            GeomechDiscretizationGlobalToNatural, &
            GeomechDiscretizationNaturalToGlobal, &
            GeomechDiscretizationGlobalToLocalBegin, &
            GeomechDiscretizationGlobalToLocalEnd, &
            GeomechDiscretizationLocalToLocalBegin, &
            GeomechDiscretizationLocalToLocalEnd, &
            GeomechDiscretizGlobalToNaturalBegin, &
            GeomechDiscretizGlobalToNaturalEnd, &
            GeomechDiscretizNaturalToGlobalBegin, &
            GeomechDiscretizNaturalToGlobalEnd, &
            GeomechDiscretizationCreateDMs,&
            GeomechDiscretizationGetDMPtrFromIndex, &
            GeomechDiscretAOApplicationToPetsc
            
contains

! ************************************************************************** !

function GeomechDiscretizationCreate()
  ! 
  ! Creates a geomechanics discretization
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 05/23/2013
  ! 

  implicit none
  
  type(geomech_discretization_type), pointer :: GeomechDiscretizationCreate
  type(geomech_discretization_type), pointer :: geomech_discretization
  
  allocate(geomech_discretization)
  geomech_discretization%ctype = ''
  geomech_discretization%itype = 0
  geomech_discretization%origin = 0.d0
  geomech_discretization%filename = ''

  ! nullify DM pointers
  allocate(geomech_discretization%dm_1dof)
  allocate(geomech_discretization%dm_ngeodof)
  allocate(geomech_discretization%dm_n_stress_strain_dof)
  geomech_discretization%dm_1dof%dm = 0
  geomech_discretization%dm_ngeodof%dm = 0
  geomech_discretization%dm_n_stress_strain_dof%dm = 0
  nullify(geomech_discretization%dm_1dof%gmdm)
  nullify(geomech_discretization%dm_ngeodof%gmdm)
  nullify(geomech_discretization%dm_n_stress_strain_dof%gmdm)
  nullify(geomech_discretization%grid)
  
  GeomechDiscretizationCreate => geomech_discretization

end function GeomechDiscretizationCreate

! ************************************************************************** !

subroutine GeomechDiscretizationCreateDMs(geomech_discretization,option)
  ! 
  ! creates distributed, parallel meshes/grids
  ! If there are multiple degrees of freedom per grid cell, this will call
  ! GeomechDiscretizationCreateDM() multiple times to create the DMs corresponding
  ! to one degree of freedom grid cell and those corresponding to multiple
  ! degrees of freedom per cell for geomechanics.
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/02/13
  ! 
      
  use Option_module    
      
  implicit none
  
  type(geomech_discretization_type) :: geomech_discretization
  type(option_type) :: option
      
  PetscInt :: ndof
  PetscErrorCode :: ierr
  type(geomech_grid_type), pointer :: geomech_grid

  !-----------------------------------------------------------------------
  ! Generate the DM objects that will manage communication.
  !-----------------------------------------------------------------------
  ndof = 1
  call GeomechDiscretizationCreateDM(geomech_discretization, &
                                     geomech_discretization%dm_1dof, &
                                     ndof,option)
  
  if (option%ngeomechdof > 0) then
    ndof = option%ngeomechdof
    call GeomechDiscretizationCreateDM(geomech_discretization, &
                                       geomech_discretization%dm_ngeodof, &
                                       ndof,option)

    call GeomechDiscretizationCreateDM(geomech_discretization, &
                              geomech_discretization%dm_n_stress_strain_dof, &
                              option%n_stress_strain_dof,option)
  endif


end subroutine GeomechDiscretizationCreateDMs

! ************************************************************************** !

subroutine GeomechDiscretizationCreateDM(geomech_discretization,dm_ptr, &
                                         ndof,option)
  ! 
  ! creates a distributed, parallel mesh/grid
  ! for geomechanics
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/02/13
  ! 

  use Option_module
  
  implicit none
  
  type(geomech_discretization_type) :: geomech_discretization
  type(gmdm_ptr_type), pointer :: dm_ptr
  type(option_type) :: option
  PetscInt :: ndof
  PetscErrorCode :: ierr

  select case(geomech_discretization%itype)
    case(STRUCTURED_GRID)
      option%io_buffer = &
        'Geomechanics currently works only with unstructured grid.'
      call printErrMsg(option)
    case(UNSTRUCTURED_GRID)
#if !defined(PETSC_HAVE_PARMETIS)
            option%io_buffer = &
             'Must compile with Parmetis in order to use Geomechanics ' // &
             'unstructured grids.'
            call printErrMsg(option)
#endif
      call GMCreateGMDM(geomech_discretization%grid, &
                        dm_ptr%gmdm,ndof,option)
      call DMShellCreate(option%mycomm,dm_ptr%dm,ierr);CHKERRQ(ierr)
      call DMShellSetGlobalToLocalVecScatter(dm_ptr%dm, &
                                             dm_ptr%gmdm%scatter_gtol, &
                                             ierr);CHKERRQ(ierr)
  end select

end subroutine GeomechDiscretizationCreateDM

! ************************************************************************** !

subroutine GeomechDiscretizationCreateVector(geomech_discretization, &
                                             dm_index,vector, &
                                             vector_type,option)
  ! 
  ! Creates a PETSc vector for the nodes
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/02/13
  ! 
  use Option_module                                      

  implicit none
  
  type(geomech_discretization_type) :: geomech_discretization
  type(option_type) :: option
  type(gmdm_ptr_type), pointer :: dm_ptr
  PetscInt :: dm_index
  Vec :: vector
  PetscInt :: vector_type
  PetscInt :: ndof
  PetscErrorCode :: ierr
  
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(geomech_discretization, &
                                                   dm_index)

  call GMGridDMCreateVector(geomech_discretization%grid,dm_ptr%gmdm,vector, &
                            vector_type,option)
                            
  call VecSet(vector,0.d0,ierr);CHKERRQ(ierr)
  
end subroutine GeomechDiscretizationCreateVector

! ************************************************************************** !

subroutine GeomechDiscretizationDuplicateVector(geomech_discretization, &
                                                vector1,vector2)
  ! 
  ! Duplicates a Petsc vector
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/02/13
  ! 

  implicit none
  
  type(geomech_discretization_type) :: geomech_discretization
  Vec :: vector1
  Vec :: vector2
  PetscErrorCode :: ierr
  
  call VecDuplicate(vector1,vector2,ierr);CHKERRQ(ierr)
  call VecCopy(vector1,vector2,ierr);CHKERRQ(ierr)
  
end subroutine GeomechDiscretizationDuplicateVector

! ************************************************************************** !

function GeomechDiscretizationGetDMPtrFromIndex(geomech_discretization,dm_index)
  ! 
  ! Returns the integer pointer for
  ! the Geomech DM referenced
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/02/13
  ! 

  implicit none
  
  type(geomech_discretization_type) :: geomech_discretization
  type(gmdm_ptr_type), pointer :: GeomechDiscretizationGetDMPtrFromIndex
  PetscInt :: dm_index
  
  select case (dm_index)
    case(ONEDOF)
      GeomechDiscretizationGetDMPtrFromIndex => &
        geomech_discretization%dm_1dof
    case(NGEODOF)
      GeomechDiscretizationGetDMPtrFromIndex => &
        geomech_discretization%dm_ngeodof
    case(SIX_INTEGER)
      GeomechDiscretizationGetDMPtrFromIndex => &
        geomech_discretization%dm_n_stress_strain_dof
  end select  
  
end function GeomechDiscretizationGetDMPtrFromIndex

! ************************************************************************** !

subroutine GeomechDiscretizationCreateJacobian(geomech_discretization, &
                                               dm_index, &
                                               mat_type,Jacobian,option)
  ! 
  ! Creates Jacobian matrix associated
  ! with geomechanics discretization
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/05/13
  ! 

  use Option_module
  
  implicit none

  type(geomech_discretization_type) :: geomech_discretization
  type(option_type) :: option
  type(gmdm_ptr_type), pointer :: dm_ptr
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  MatType :: mat_type
  Mat :: Jacobian

  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(geomech_discretization, &
                                                    dm_index)


  call GMGridDMCreateJacobian(geomech_discretization%grid,dm_ptr%gmdm, &
                              mat_type,Jacobian,option)
  call MatSetOption(Jacobian,MAT_KEEP_NONZERO_PATTERN,PETSC_FALSE, &
                    ierr);CHKERRQ(ierr)
  call MatSetOption(Jacobian,MAT_ROW_ORIENTED,PETSC_FALSE,ierr);CHKERRQ(ierr)

end subroutine GeomechDiscretizationCreateJacobian

! ************************************************************************** !

subroutine GeomechDiscretizationGlobalToLocal(geomech_discretization, &
                                              global_vec, &
                                              local_vec,dm_index)
  ! 
  ! Performs global to local communication
  ! with geomech DM
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/02/13
  ! 

  implicit none

  type(geomech_discretization_type) :: geomech_discretization
  type(gmdm_ptr_type), pointer :: dm_ptr
  Vec :: global_vec
  Vec :: local_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(geomech_discretization, &
                                                   dm_index)
    
  call DMGlobalToLocalBegin(dm_ptr%dm,global_vec,INSERT_VALUES,local_vec, &
                            ierr);CHKERRQ(ierr)
  call DMGlobalToLocalEnd(dm_ptr%dm,global_vec,INSERT_VALUES,local_vec, &
                          ierr);CHKERRQ(ierr)
  
end subroutine GeomechDiscretizationGlobalToLocal

! ************************************************************************** !

subroutine GeomechDiscretizationLocalToGlobal(geomech_discretization, & 
                                              local_vec, &
                                              global_vec,dm_index)
  ! 
  ! Performs local to global communication
  ! with DM
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/02/13
  ! 

  implicit none
  
  type(geomech_discretization_type) :: geomech_discretization
  type(gmdm_ptr_type), pointer :: dm_ptr
  Vec :: local_vec
  Vec :: global_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(geomech_discretization, &
                                                   dm_index)
  
  call VecScatterBegin(dm_ptr%gmdm%scatter_ltog,local_vec,global_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(dm_ptr%gmdm%scatter_ltog,local_vec,global_vec, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  
end subroutine GeomechDiscretizationLocalToGlobal

! ************************************************************************** !

subroutine GeomechDiscretizationLocalToGlobalAdd(geomech_discretization, &
                                                 local_vec, &
                                                 global_vec,dm_index)
  ! 
  ! Performs local to global communication
  ! with DM and adds
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 09/17/13
  ! 

  implicit none
  
  type(geomech_discretization_type) :: geomech_discretization
  type(gmdm_ptr_type), pointer :: dm_ptr
  Vec :: local_vec
  Vec :: global_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(geomech_discretization, &
                                                   dm_index)
  
  call VecScatterBegin(dm_ptr%gmdm%scatter_ltog,local_vec,global_vec, &
                       ADD_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(dm_ptr%gmdm%scatter_ltog,local_vec,global_vec, &
                     ADD_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  
end subroutine GeomechDiscretizationLocalToGlobalAdd

! ************************************************************************** !

subroutine GeomechDiscretizationLocalToLocal(geomech_discretization, &
                                             local_vec1, &
                                             local_vec2,dm_index)
  ! 
  ! Performs local to local communication
  ! with geomech DM
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/02/13
  ! 

  implicit none
  
  type(geomech_discretization_type) :: geomech_discretization
  type(gmdm_ptr_type), pointer :: dm_ptr
  Vec :: local_vec1
  Vec :: local_vec2
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(geomech_discretization, &
                                                   dm_index)
  
  call VecScatterBegin(dm_ptr%gmdm%scatter_ltol,local_vec1,local_vec2, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(dm_ptr%gmdm%scatter_ltol,local_vec1,local_vec2, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  
end subroutine GeomechDiscretizationLocalToLocal

! ************************************************************************** !

subroutine GeomechDiscretizationGlobalToNatural(geomech_discretization, &
                                                global_vec, &
                                                natural_vec,dm_index)
  ! 
  ! Performs global to natural
  ! communication with geomech DM
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/02/13
  ! 

  implicit none
  
  type(geomech_discretization_type) :: geomech_discretization
  type(gmdm_ptr_type), pointer :: dm_ptr
  Vec :: global_vec
  Vec :: natural_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(geomech_discretization, &
                                                   dm_index)

  call VecScatterBegin(dm_ptr%gmdm%scatter_gton,global_vec,natural_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(dm_ptr%gmdm%scatter_gton,global_vec,natural_vec, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  
end subroutine GeomechDiscretizationGlobalToNatural

! ************************************************************************** !

subroutine GeomechDiscretizationNaturalToGlobal(geomech_discretization, &
                                                natural_vec, &
                                                global_vec,dm_index)
  ! 
  ! Performs natural to global
  ! communication with DM
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/02/13
  ! 

  implicit none
  
  type(geomech_discretization_type) :: geomech_discretization
  type(gmdm_ptr_type), pointer :: dm_ptr
  Vec :: global_vec
  Vec :: natural_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(geomech_discretization, &
                                                   dm_index)
  
  call VecScatterBegin(dm_ptr%gmdm%scatter_gton,natural_vec,global_vec, &
                       INSERT_VALUES,SCATTER_REVERSE,ierr);CHKERRQ(ierr)
  call VecScatterEnd(dm_ptr%gmdm%scatter_gton,natural_vec,global_vec, &
                     INSERT_VALUES,SCATTER_REVERSE,ierr);CHKERRQ(ierr)
  
end subroutine GeomechDiscretizationNaturalToGlobal

! ************************************************************************** !

subroutine GeomechDiscretizationGlobalToLocalBegin(geomech_discretization, &
                                                   global_vec, &
                                                   local_vec,dm_index)
  ! 
  ! Begins global to local
  ! communication with geomech DM
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/02/13
  ! 

  implicit none
  
  type(geomech_discretization_type) :: geomech_discretization
  type(gmdm_ptr_type), pointer :: dm_ptr
  Vec :: local_vec
  Vec :: global_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(geomech_discretization, &
                                                   dm_index)
  
  call VecScatterBegin(dm_ptr%gmdm%scatter_gtol,global_vec,local_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  
end subroutine GeomechDiscretizationGlobalToLocalBegin

! ************************************************************************** !

subroutine GeomechDiscretizationGlobalToLocalEnd(geomech_discretization, &
                                                 global_vec, &
                                                 local_vec,dm_index)
  ! 
  ! Ends global to local communication
  ! with geomech DM
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/02/13
  ! 

 implicit none
  
  type(geomech_discretization_type) :: geomech_discretization
  type(gmdm_ptr_type), pointer :: dm_ptr
  Vec :: local_vec
  Vec :: global_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(geomech_discretization, &
                                                   dm_index)
  
  call VecScatterEnd(dm_ptr%gmdm%scatter_gtol,global_vec,local_vec, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  
end subroutine GeomechDiscretizationGlobalToLocalEnd

! ************************************************************************** !

subroutine GeomechDiscretizationLocalToLocalBegin(geomech_discretization, &
                                                  local_vec1, &
                                                  local_vec2,dm_index)
  ! 
  ! Begins local to local communication
  ! with geomech DM
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/02/13
  ! 

  implicit none
  
  type(geomech_discretization_type) :: geomech_discretization
  type(gmdm_ptr_type), pointer :: dm_ptr
  Vec :: local_vec1
  Vec :: local_vec2
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(geomech_discretization, &
                                                   dm_index)
  
  call VecScatterBegin(dm_ptr%gmdm%scatter_ltol,local_vec1,local_vec2, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)

end subroutine GeomechDiscretizationLocalToLocalBegin

! ************************************************************************** !

subroutine GeomechDiscretizationLocalToLocalEnd(geomech_discretization, &
                                                local_vec1, &
                                                local_vec2,dm_index)
  ! 
  ! Ends local to local communication
  ! with geomech DM
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/02/13
  ! 

  implicit none
  
  type(geomech_discretization_type) :: geomech_discretization
  type(gmdm_ptr_type), pointer :: dm_ptr
  Vec :: local_vec1
  Vec :: local_vec2
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(geomech_discretization, &
                                                   dm_index)
  
  call VecScatterEnd(dm_ptr%gmdm%scatter_ltol,local_vec1,local_vec2, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)

end subroutine GeomechDiscretizationLocalToLocalEnd

! ************************************************************************** !

subroutine GeomechDiscretizGlobalToNaturalBegin(geomech_discretization, &
                                                global_vec, &
                                                natural_vec,dm_index)
  ! 
  ! Begins global to natural communication
  ! with geomech DM
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/02/13
  ! 

  implicit none
  
  type(geomech_discretization_type) :: geomech_discretization
  type(gmdm_ptr_type), pointer :: dm_ptr
  Vec :: global_vec
  Vec :: natural_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(geomech_discretization, &
                                                   dm_index)
  
  call VecScatterBegin(dm_ptr%gmdm%scatter_gton,global_vec,natural_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  
end subroutine GeomechDiscretizGlobalToNaturalBegin

! ************************************************************************** !

subroutine GeomechDiscretizGlobalToNaturalEnd(geomech_discretization, &
                                              global_vec, &
                                              natural_vec,dm_index)
  ! 
  ! Ends global to natural communication
  ! with geomech DM
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/02/13
  ! 

  implicit none
  
  type(geomech_discretization_type) :: geomech_discretization
  type(gmdm_ptr_type), pointer :: dm_ptr
  Vec :: global_vec
  Vec :: natural_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(geomech_discretization, &
                                                   dm_index)
  
  call VecScatterEnd(dm_ptr%gmdm%scatter_gton,global_vec,natural_vec, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  
end subroutine GeomechDiscretizGlobalToNaturalEnd

! ************************************************************************** !

subroutine GeomechDiscretizNaturalToGlobalBegin(geomech_discretization, &
                                                natural_vec, &
                                                global_vec,dm_index)
  ! 
  ! Begins natural to global communication
  ! with geomech DM
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/02/13
  ! 

  implicit none
  
  type(geomech_discretization_type) :: geomech_discretization
  type(gmdm_ptr_type), pointer :: dm_ptr
  Vec :: global_vec
  Vec :: natural_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(geomech_discretization, &
                                                   dm_index)
    
end subroutine GeomechDiscretizNaturalToGlobalBegin

! ************************************************************************** !

subroutine GeomechDiscretizNaturalToGlobalEnd(geomech_discretization, &
                                              natural_vec, &
                                              global_vec,dm_index)
  ! 
  ! Ends natural to global communication
  ! with geomech DM
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/02/13
  ! 

  implicit none
  
  type(geomech_discretization_type) :: geomech_discretization
  type(gmdm_ptr_type), pointer :: dm_ptr
  Vec :: global_vec
  Vec :: natural_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex(geomech_discretization, &
                                                   dm_index)
  
end subroutine GeomechDiscretizNaturalToGlobalEnd

! ************************************************************************** !

subroutine GeomechDiscretAOApplicationToPetsc(geomech_discretization,int_array)
  ! 
  ! Maps application ordering to petsc
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/02/13
  ! 

  implicit none
  
#include "petsc/finclude/petscao.h"  
  
  type(geomech_discretization_type) :: geomech_discretization
  PetscInt :: int_array(:)
  PetscErrorCode :: ierr
  AO :: ao
  
  ao = geomech_discretization%grid%ao_natural_to_petsc_nodes
  
  call AOApplicationToPetsc(ao,size(int_array),int_array,ierr);CHKERRQ(ierr)
  
end subroutine GeomechDiscretAOApplicationToPetsc

! ************************************************************************** !

subroutine GeomechDiscretizationDestroy(geomech_discretization)
  ! 
  ! Deallocates a geomechanics discretization
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 05/23/2013
  ! 

  implicit none
  
  type(geomech_discretization_type), pointer :: geomech_discretization
  
  PetscErrorCode :: ierr
  PetscInt :: i

  if (.not.associated(geomech_discretization)) return
      
  if (associated(geomech_discretization%dm_1dof%gmdm)) &
    call GMDMDestroy(geomech_discretization%dm_1dof%gmdm)
  if (associated(geomech_discretization%dm_ngeodof%gmdm)) &
    call GMDMDestroy(geomech_discretization%dm_ngeodof%gmdm)
  if (associated(geomech_discretization%dm_n_stress_strain_dof%gmdm)) &
    call GMDMDestroy(geomech_discretization%dm_n_stress_strain_dof%gmdm)

  if (associated(geomech_discretization%dm_1dof)) &
    deallocate(geomech_discretization%dm_1dof)
  nullify(geomech_discretization%dm_1dof)
  if (associated(geomech_discretization%dm_ngeodof)) &
    deallocate(geomech_discretization%dm_ngeodof)
  nullify(geomech_discretization%dm_ngeodof)
  if (associated(geomech_discretization%dm_n_stress_strain_dof)) &
    deallocate(geomech_discretization%dm_n_stress_strain_dof)
  nullify(geomech_discretization%dm_n_stress_strain_dof)
  
  call GMGridDestroy(geomech_discretization%grid)

  if (associated(geomech_discretization)) &
    deallocate(geomech_discretization)
  nullify(geomech_discretization)


end subroutine GeomechDiscretizationDestroy

end module Geomechanics_Discretization_module
