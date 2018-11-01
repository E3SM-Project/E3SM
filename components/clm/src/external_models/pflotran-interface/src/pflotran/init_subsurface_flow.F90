module Init_Subsurface_Flow_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: InitSubsurfFlowSetupRealization
  
contains

! ************************************************************************** !

subroutine InitSubsurfFlowSetupRealization(realization)
  ! 
  ! Initializes material property data structres and assign them to the domain.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/04/14
  ! 
  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Init_Common_module
  use Material_module
  
  use TH_module
  use Condition_Control_module
  
  implicit none
  
  class(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  PetscErrorCode :: ierr
  
  option => realization%option
  patch => realization%patch
  
  ! initialize FLOW
  ! set up auxillary variable arrays
  if (option%nflowdof > 0) then
    select case(option%iflowmode)
      case(TH_MODE)
        call MaterialSetup(realization%patch%aux%Material%material_parameter, &
                           patch%material_property_array, &
                           patch%characteristic_curves_array, &
                           realization%option)
    end select
    select case(option%iflowmode)
      case(TH_MODE)
        call THSetup(realization)
    end select
  
    ! assign initial conditionsRealizAssignFlowInitCond
    call CondControlAssignFlowInitCond(realization)

    ! override initial conditions if they are to be read from a file
    if (len_trim(option%initialize_flow_filename) > 1) then
      call InitSubsurfFlowReadInitCond(realization, &
                                       option%initialize_flow_filename)
    endif
  
    select case(option%iflowmode)
      case(TH_MODE)
        call THUpdateAuxVars(realization)
    end select
  else ! no flow mode specified
    if (len_trim(realization%nonuniform_velocity_filename) > 0) then
#if defined(PETSC_HAVE_HDF5)
      call InitCommonReadVelocityField(realization)
#else
      option%io_buffer = 'PFLOTRAN must be compiled with HDF5 to read HDF5 &
        &formatted fluxes in for transport with no flow.'
      call printErrMsg(option)
#endif
    endif
  endif  
  
end subroutine InitSubsurfFlowSetupRealization

! ************************************************************************** !

subroutine InitSubsurfFlowReadInitCond(realization,filename)
  ! 
  ! Assigns flow initial condition from HDF5 file
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/05/10, 12/04/14
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Grid_module
  use Patch_module
  use Discretization_module
  use HDF5_module
  
  implicit none
  
  class(realization_subsurface_type) :: realization
  character(len=MAXSTRINGLENGTH) :: filename
  
  PetscInt :: local_id, idx, offset
  PetscReal, pointer :: xx_p(:)
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscReal, pointer :: vec_p(:)  
  PetscErrorCode :: ierr
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(patch_type), pointer :: cur_patch

  option => realization%option
  discretization => realization%discretization
  field => realization%field
  patch => realization%patch

  if (option%iflowmode /= TH_MODE) then
    option%io_buffer = 'Reading of flow initial conditions from HDF5 ' // &
                       'file (' // trim(filename) // &
                       'not currently not supported for mode: ' // &
                       trim(option%flowmode)
  endif      

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit

    grid => cur_patch%grid

      ! assign initial conditions values to domain
    call VecGetArrayF90(field%flow_xx, xx_p, ierr);CHKERRQ(ierr)

    ! Pressure for all modes 
    if (option%iflowmode == TH_MODE) then
      offset = TH_PRESSURE_DOF
    endif
    !offset = 1
    group_name = ''
    dataset_name = 'Pressure'
    call HDF5ReadCellIndexedRealArray(realization,field%work, &
                                      filename,group_name, &
                                      dataset_name,option%id>0)
    call VecGetArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)
    do local_id=1, grid%nlmax
      if (cur_patch%imat(grid%nL2G(local_id)) <= 0) cycle
      if (dabs(vec_p(local_id)) < 1.d-40) then
        print *,  option%myrank, grid%nG2A(grid%nL2G(local_id)), &
              ': Potential error - zero pressure in Initial Condition read from file.'
      endif
      idx = (local_id-1)*option%nflowdof + offset
      xx_p(idx) = vec_p(local_id)
    enddo
    call VecRestoreArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)

    ! Temperature for TH mode
    if (option%iflowmode == TH_MODE) then
      offset = TH_TEMPERATURE_DOF
      group_name = ''
      dataset_name = 'Temperature'
      call HDF5ReadCellIndexedRealArray(realization,field%work, &
                                      filename,group_name, &
                                      dataset_name,option%id>0)
      call VecGetArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)
      do local_id=1, grid%nlmax
        if (cur_patch%imat(grid%nL2G(local_id)) <= 0) cycle
        if (dabs(vec_p(local_id)) < 1.d-40) then
          print *,  option%myrank, grid%nG2A(grid%nL2G(local_id)), &
              ': Potential error - zero pressure in Initial Condition read from file.'
        endif
        idx = (local_id-1)*option%nflowdof + offset
        xx_p(idx) = vec_p(local_id)
      enddo
      call VecRestoreArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)
    endif

    call VecRestoreArrayF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)
        
    cur_patch => cur_patch%next
  enddo
   
  ! update dependent vectors
  call DiscretizationGlobalToLocal(discretization,field%flow_xx, &
                                   field%flow_xx_loc,NFLOWDOF)  
  call VecCopy(field%flow_xx, field%flow_yy, ierr);CHKERRQ(ierr)

end subroutine InitSubsurfFlowReadInitCond

end module Init_Subsurface_Flow_module
