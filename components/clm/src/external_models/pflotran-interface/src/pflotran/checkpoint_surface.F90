module Checkpoint_Surface_Header_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  implicit none

  private

! We must manually specify the number of bytes required for the 
! checkpoint header ('surface_bagsize'), since sizeof() is not supported by 
! some Fortran compilers.  To be on the safe side, we assume an integer is 8 
! bytes.
! Currently:
!   PetscReal: 5 
!   PetscInt: 5
!   Total: 10 * 8 = 80
! IMPORTANT: If you change the contents of the header, you MUST update 
! 'surface_bagsize' or risk corrupting memory.
#ifdef PetscSizeT
  PetscSizeT, parameter :: surface_bagsize = 80
#else
  ! PETSC_SIZEOF_SIZE_T isn't defined, so we just have to assume that it 
  ! is 8 bytes.  This is dangerous, but what can we do?
  integer*8, parameter :: surface_bagsize = 80
#endif

  public :: surface_bagsize
  type, public :: surface_checkpoint_header_type

    integer*8 :: revision_number  ! increment this every time there is a change

    integer*8 :: grid_discretization_type

    integer*8 :: nsurfflowdof
    integer*8 :: surface_flow_formulation
    real*8 :: surf_flow_time
    real*8 :: surf_flow_dt
    real*8 :: surf_flow_prev_dt

    integer*8 :: subsurf_surf_coupling
    real*8 :: surf_subsurf_coupling_time
    real*8 :: surf_subsurf_coupling_flow_dt

  end type surface_checkpoint_header_type
end module Checkpoint_Surface_Header_module

module Checkpoint_Surface_module
#include "petsc/finclude/petscdm.h"
  use petscdm
  use Checkpoint_Surface_Header_module

  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: SurfaceCheckpointBinary, SurfaceRestartBinary
  public :: SurfaceCheckpointProcessModelBinary, &
            SurfaceRestartProcessModelBinary

Interface PetscBagGetData
Subroutine PetscBagGetData(bag,ctx,ierr)
      use Checkpoint_Surface_Header_module
      PetscBag bag
      type(surface_checkpoint_header_type), pointer :: ctx
      PetscErrorCode ierr
End Subroutine
End Interface PetscBagGetData

contains

! ************************************************************************** !

subroutine SurfaceCheckpointBinary(surf_realization, &
                                   surf_flow_prev_dt, &
                                   id)
  ! 
  ! This subroutine writes a checkpoint file for surface realization.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/11/13
  ! 

  use Realization_Surface_class
  use Surface_Field_module
  use Grid_module
  use Discretization_module
  use Output_Aux_module
  use Option_module

  implicit none

  class(realization_surface_type) :: surf_realization
  PetscReal :: surf_flow_prev_dt
  PetscInt, intent(in) :: id ! id should not be altered within this subroutine

  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXWORDLENGTH) :: id_string
  PetscViewer :: viewer
  PetscBag :: bag
  type(surface_checkpoint_header_type), pointer :: surf_header
  PetscErrorCode :: ierr
  PetscLogDouble :: tstart, tend  
  
  Vec :: global_vec
  PetscInt :: int_flag
  
  type(surface_field_type), pointer :: surf_field
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(output_option_type), pointer :: output_option
  PetscInt :: i, j, k

  surf_field => surf_realization%surf_field
  option => surf_realization%option
  discretization => surf_realization%discretization
  grid => discretization%grid 

  option%io_buffer = 'Checkpointing of surface flow must be updated.'
  call printErrMsg(option)

  ! Open the checkpoint file.
  call PetscTime(tstart,ierr);CHKERRQ(ierr)
  if (id < 0) then
    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
               '-restart-surf.chk'
  else 
    write(id_string,'(i8)') id
    filename = trim(option%global_prefix) // trim(option%group_prefix) // &
               '-surf.chk' // trim(adjustl(id_string))
  endif

  call PetscViewerCreate(option%mycomm,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerSetType(viewer,PETSCVIEWERBINARY,ierr);CHKERRQ(ierr)
  call PetscViewerFileSetMode(viewer,FILE_MODE_WRITE,ierr);CHKERRQ(ierr)
  call PetscViewerBinarySkipInfo(viewer,ierr);CHKERRQ(ierr)
  call PetscViewerFileSetName(viewer,filename,ierr);CHKERRQ(ierr)

  call PetscBagCreate(option%mycomm,surface_bagsize, bag, ierr);CHKERRQ(ierr)
  call PetscBagGetData(bag, surf_header, ierr);CHKERRQ(ierr)

  call SurfCheckpointRegisterBagHeader(bag,surf_header)

  ! Register
  surf_header%grid_discretization_type = grid%itype

  surf_header%nsurfflowdof = option%nsurfflowdof
  surf_header%surface_flow_formulation = option%surface_flow_formulation
  surf_header%surf_flow_time = option%surf_flow_time
  surf_header%surf_flow_dt = option%surf_flow_dt
  surf_header%surf_flow_prev_dt = surf_flow_prev_dt

  surf_header%subsurf_surf_coupling = option%subsurf_surf_coupling
  surf_header%surf_subsurf_coupling_time = option%surf_subsurf_coupling_time
  surf_header%surf_subsurf_coupling_flow_dt = option%surf_subsurf_coupling_flow_dt

  ! Actually write the components of the PetscBag and then free it.
  call PetscBagView(bag, viewer, ierr);CHKERRQ(ierr)
  call PetscBagDestroy(bag, ierr);CHKERRQ(ierr)

  !--------------------------------------------------------------------
  ! Dump all the relevant vectors.
  !--------------------------------------------------------------------

  if (option%nflowdof > 0) then
    call DiscretizationCreateVector(discretization,ONEDOF, &
                                    global_vec,GLOBAL,option)
    ! grid%flow_xx is the vector into which all of the primary variables are 
    ! packed for the TSSolve().
    call VecView(surf_field%flow_xx, viewer, ierr);CHKERRQ(ierr)

    ! Mannings coefficient.
    call DiscretizationLocalToGlobal(discretization,surf_field%mannings_loc, &
                                     global_vec,ONEDOF)
    call VecView(global_vec,viewer,ierr);CHKERRQ(ierr)
  endif

  if (global_vec /= PETSC_NULL_VEC) then
    call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  endif

  ! We are finished, so clean up.
  call PetscViewerDestroy(viewer, ierr);CHKERRQ(ierr)

  write(option%io_buffer,'(" --> Dump checkpoint file: ", a32)') trim(filename)
  call printMsg(option)

  call PetscTime(tend,ierr);CHKERRQ(ierr)
  write(option%io_buffer, &
        '("      Seconds to write to checkpoint file: ", f10.2)') tend-tstart
  call printMsg(option)

end subroutine SurfaceCheckpointBinary

! ************************************************************************** !

subroutine SurfaceRestartBinary(surf_realization, surf_flow_prev_dt, surf_flow_read)
  ! 
  ! This subroutine restarts surface-realization simulation by reading a
  ! checkpoint file.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/11/13
  ! 


  use Realization_Surface_class
  use Surface_Field_module
  use Grid_module
  use Discretization_module
  use Output_Aux_module
  use Option_module

  implicit none

  class(realization_surface_type) :: surf_realization
  PetscReal :: surf_flow_prev_dt
  PetscBool :: surf_flow_read
  
  type(surface_field_type), pointer :: surf_field
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization

  PetscViewer :: viewer
  PetscBag :: bag
  type(surface_checkpoint_header_type), pointer :: surf_header
  PetscLogDouble :: tstart, tend
  PetscErrorCode :: ierr

  Vec :: global_vec
  character(len=MAXSTRINGLENGTH) :: string

  surf_field => surf_realization%surf_field
  option => surf_realization%option
  discretization => surf_realization%discretization
  grid => discretization%grid 

  call PetscTime(tstart,ierr);CHKERRQ(ierr)
  option%io_buffer = '--> Open checkpoint file: ' // &
                     trim(option%surf_restart_filename)
  call printMsg(option)
  call PetscViewerBinaryOpen(option%mycomm,option%surf_restart_filename, &
                             FILE_MODE_READ,viewer,ierr);CHKERRQ(ierr)
  ! skip reading info file when loading, but not working
  call PetscViewerBinarySetSkipOptions(viewer,PETSC_TRUE,ierr);CHKERRQ(ierr)

  ! Get the header data.
  call PetscBagCreate(option%mycomm, surface_bagsize, bag, ierr);CHKERRQ(ierr)
  call PetscBagGetData(bag, surf_header, ierr);CHKERRQ(ierr)
  call SurfCheckpointRegisterBagHeader(bag,surf_header)
  call PetscBagLoad(viewer, bag, ierr);CHKERRQ(ierr)
  
  if (surf_header%revision_number /= CHECKPOINT_REVISION_NUMBER) then
    write(string,*) surf_header%revision_number
    option%io_buffer = 'The revision number # of checkpoint file (' // &
                       trim(option%surf_restart_filename) // ', rev=' // &
                       trim(adjustl(string)) // &
                       ') does not match the current revision number' // &
                       ' of PFLOTRAN checkpoint files ('
    write(string,*) CHECKPOINT_REVISION_NUMBER
    option%io_buffer = trim(option%io_buffer) // trim(adjustl(string)) // ').'
    call printErrMsg(option)
  endif

   if (surf_header%grid_discretization_type /= grid%itype) then
     write(string,*) surf_header%grid_discretization_type
     option%io_buffer = 'The discretization of checkpoint file (' // &
                       trim(option%restart_filename) // ', grid_type=' // &
                       trim(adjustl(string)) // &
                       ') does not match the discretization of the current problem' // &
                       ' grid_type= ('
    write(string,*) grid%itype
    option%io_buffer = trim(option%io_buffer) // trim(adjustl(string)) // ').'
    call printErrMsg(option)
  endif

  ! Check DOFs in surface-flow
  if (option%nsurfflowdof /= surf_header%nsurfflowdof) then
    write(string,*) surf_header%nsurfflowdof
    option%io_buffer = 'Number of surface-flow dofs in restart file (' // &
                       trim(adjustl(string)) // &
           ') does not match the number of surface-flow dofs in the input file ('
    write(string,*) option%nsurfflowdof
    option%io_buffer = trim(option%io_buffer) // string // ')'
    call printErrMsg(option)
  endif

  ! Check surface-flow formulation
  if (option%surface_flow_formulation /= surf_header%surface_flow_formulation) then
    option%io_buffer = 'Surface flow formulation in restart file ' // &
      ' does not match the surface flow formulation the input file '
    option%io_buffer = trim(option%io_buffer)
    call printErrMsg(option)
  endif

  ! Save values from header
  if (option%surf_flow_on .and. &
      option%nsurfflowdof == surf_header%nsurfflowdof) then

    option%surf_flow_time = surf_header%surf_flow_time
    option%surf_flow_dt = surf_header%surf_flow_dt
    surf_flow_prev_dt = surf_header%surf_flow_prev_dt

    option%subsurf_surf_coupling = int(surf_header%subsurf_surf_coupling)
    option%surf_subsurf_coupling_time = surf_header%surf_subsurf_coupling_time
    option%surf_subsurf_coupling_flow_dt = surf_header%surf_subsurf_coupling_flow_dt

    surf_flow_read = PETSC_TRUE
  endif

  if (surf_flow_read) then
    call DiscretizationCreateVector(discretization,ONEDOF, &
                                    global_vec,GLOBAL,option)
    ! Load the PETSc vectors.
    call VecLoad(surf_field%flow_xx,viewer,ierr);CHKERRQ(ierr)
    call DiscretizationGlobalToLocal(discretization,surf_field%flow_xx, &
                                     surf_field%flow_xx_loc,NFLOWDOF)
    call VecCopy(surf_field%flow_xx,surf_field%flow_yy,ierr);CHKERRQ(ierr)

    call VecLoad(global_vec,viewer,ierr);CHKERRQ(ierr)
    call DiscretizationGlobalToLocal(discretization,global_vec, &
                                     surf_field%mannings_loc,ONEDOF)
  endif

  ! We are finished, so clean up.
  if (global_vec /= PETSC_NULL_VEC) then
    call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  endif

  call PetscViewerDestroy(viewer, ierr);CHKERRQ(ierr)
  call PetscTime(tend,ierr);CHKERRQ(ierr)

  call PetscBagDestroy(bag, ierr);CHKERRQ(ierr)

  write(option%io_buffer, &
        '("      Seconds to read to checkpoint file: ", f6.2)') tend-tstart
  call printMsg(option)


end subroutine SurfaceRestartBinary

! ************************************************************************** !

subroutine SurfCheckpointRegisterBagHeader(bag,header)
  ! 
  ! This subroutine registers entities within the PETSc bag to header for
  ! surface-realization.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/11/13
  ! 

  implicit none
  
  PetscBag :: bag
  type(surface_checkpoint_header_type), pointer :: header

  PetscInt :: i
  PetscErrorCode :: ierr
  
  i = CHECKPOINT_REVISION_NUMBER
  call PetscBagRegisterInt(bag,header%revision_number,i, &
                           "revision_number", &
                           "revision_number", &
                           ierr);CHKERRQ(ierr)
  ! Register variables that are passed into timestepper().
  call PetscBagRegisterInt(bag,header%grid_discretization_type, 0, &
                           "grid_discretization_type", &
                           "grid_discretization_type", &
                           ierr);CHKERRQ(ierr)
  ! Surface-flow
  call PetscBagRegisterInt(bag,header%nsurfflowdof,0, &
                           "nsurfflowdof", &
                           "Number of surface flow degrees of freedom", &
                           ierr);CHKERRQ(ierr)
  call PetscBagRegisterInt(bag,header%surface_flow_formulation,0, &
                           "surface_flow_formulation", &
                           "Type of surface-flow formulation", &
                           ierr);CHKERRQ(ierr)
  
  call PetscBagRegisterReal(bag,header%surf_flow_time,0.d0, &
                            "surf_flow_time", &
                            "Surface Flow Simulation time (seconds)", &
                            ierr);CHKERRQ(ierr)
  call PetscBagRegisterReal(bag,header%surf_flow_dt,0.d0, &
                            "surf_flow_dt", &
                            "Current size of surface flow timestep (seconds)", &
                            ierr);CHKERRQ(ierr)
  call PetscBagRegisterReal(bag,header%surf_flow_prev_dt,0.d0, &
                            "surf_flow_prev_dt", &
                            "Previous size of surfae flow timestep (seconds)", &
                            ierr);CHKERRQ(ierr)

  ! Surface-Subsurface coupling
  call PetscBagRegisterReal(bag,header%subsurf_surf_coupling,0, &
                            "subsurf_surf_coupling", &
                            "Type of surface-subsurface coupling", &
                            ierr);CHKERRQ(ierr)
  call PetscBagRegisterReal(bag,header%surf_subsurf_coupling_time,0.d0, &
                            "surf_subsurf_coupling_time", &
                            "Surface-Subsurface coupling time (seconds)", &
                            ierr);CHKERRQ(ierr)
  call PetscBagRegisterReal(bag,header%surf_subsurf_coupling_flow_dt,0.d0, &
                            "surf_subsurf_coupling_flow_dt", &
                            "Surface-Subsurface coupling timestep (seconds)", &
                            ierr);CHKERRQ(ierr)

end subroutine SurfCheckpointRegisterBagHeader

! ************************************************************************** !

subroutine SurfaceCheckpointProcessModelBinary(viewer, surf_realization)
  ! 
  ! This subroutine writes a checkpoint file for surface realization using
  ! process model approach.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 09/19/13
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Realization_Surface_class
  use Surface_Field_module
  use Grid_module
  use Discretization_module
  use Output_Aux_module
  use Option_module

  implicit none

  class(realization_surface_type) :: surf_realization
  PetscViewer :: viewer

  type(surface_field_type), pointer :: surf_field
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  PetscErrorCode :: ierr
  Vec :: global_vec

  surf_field => surf_realization%surf_field
  option => surf_realization%option
  discretization => surf_realization%discretization
  grid => discretization%grid

  global_vec = PETSC_NULL_VEC
  !--------------------------------------------------------------------
  ! Dump all the relevant vectors.
  !--------------------------------------------------------------------

  if (option%nflowdof > 0) then
    call DiscretizationCreateVector(discretization,ONEDOF, &
                                    global_vec,GLOBAL,option)
    ! grid%flow_xx is the vector into which all of the primary variables are
    ! packed for the TSSolve().
    call VecView(surf_field%flow_xx, viewer, ierr);CHKERRQ(ierr)

    ! Mannings coefficient.
    call DiscretizationLocalToGlobal(discretization,surf_field%mannings_loc, &
                                     global_vec,ONEDOF)
    call VecView(global_vec,viewer,ierr);CHKERRQ(ierr)
  endif

  if (global_vec /= PETSC_NULL_VEC) then
    call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  endif

end subroutine SurfaceCheckpointProcessModelBinary

! ************************************************************************** !

subroutine SurfaceRestartProcessModelBinary(viewer,surf_realization)
  ! 
  ! This subroutine reads a checkpoint file for surface realization using
  ! process model approach.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 09/19/13
  ! 

  use Realization_Surface_class
  use Surface_Field_module
  use Grid_module
  use Discretization_module
  use Output_Aux_module
  use Option_module

  implicit none

  class(realization_surface_type) :: surf_realization
  PetscViewer :: viewer

  type(surface_field_type), pointer :: surf_field
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  PetscErrorCode :: ierr
  Vec :: global_vec

  surf_field => surf_realization%surf_field
  option => surf_realization%option
  discretization => surf_realization%discretization
  grid => discretization%grid

  global_vec = PETSC_NULL_VEC

  if (option%surf_flow_on) then
    call DiscretizationCreateVector(discretization,ONEDOF, &
                                    global_vec,GLOBAL,option)
    ! Load the PETSc vectors.
    call VecLoad(surf_field%flow_xx,viewer,ierr);CHKERRQ(ierr)
    call DiscretizationGlobalToLocal(discretization,surf_field%flow_xx, &
                                     surf_field%flow_xx_loc,NFLOWDOF)
    call VecCopy(surf_field%flow_xx,surf_field%flow_yy,ierr);CHKERRQ(ierr)

    call VecLoad(global_vec,viewer,ierr);CHKERRQ(ierr)
    call DiscretizationGlobalToLocal(discretization,global_vec, &
                                     surf_field%mannings_loc,ONEDOF)
  endif

  ! We are finished, so clean up.
  if (global_vec /= PETSC_NULL_VEC) then
    call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  endif

end subroutine SurfaceRestartProcessModelBinary

end module Checkpoint_Surface_module
