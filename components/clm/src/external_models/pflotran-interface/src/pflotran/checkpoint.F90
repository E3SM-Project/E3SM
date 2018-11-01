module Checkpoint_module
#include "petsc/finclude/petscdm.h"
  use petscdm
  use PFLOTRAN_Constants_module

  implicit none
  
  private

  type :: checkpoint_header_type
    PetscInt :: version
    PetscInt :: test_header_size
  end type checkpoint_header_type

  type :: base_test_header_type
    PetscInt :: int1
    PetscReal :: real1
    PetscInt :: int2
    PetscReal :: real2
    PetscInt :: int3
    PetscReal :: real3
    PetscInt :: int4
  end type base_test_header_type

  type, extends(base_test_header_type) :: extended_test_header_type
    PetscReal :: real4
    PetscInt :: int5
    PetscReal :: real5
  end type extended_test_header_type

  interface PetscBagGetData
    subroutine PetscBagGetData(bag,header,ierr)
#include "petsc/finclude/petscsys.h"
      use petscsys
      import :: checkpoint_header_type
      implicit none
      PetscBag :: bag
      type(checkpoint_header_type), pointer :: header
      PetscErrorCode :: ierr
    end subroutine
  end interface PetscBagGetData

  public :: CheckpointFilename, &
            CheckpointAppendNameAtTime, &
            CheckpointAppendNameAtTimestep, &
            CheckpointOpenFileForWriteBinary, &
            CheckPointWriteCompatibilityBinary, &
            CheckPointReadCompatibilityBinary, &
            CheckpointFlowProcessModelBinary, &
            RestartFlowProcessModelBinary, &
#if defined(PETSC_HAVE_HDF5)
            RestartFlowProcessModelHDF5, &
            CheckpointOpenFileForWriteHDF5, &
            CheckPointWriteCompatibilityHDF5, &
            CheckpointFlowProcessModelHDF5, &
            CheckPointWriteIntDatasetHDF5, &
            CheckPointReadRealDatasetHDF5, &
            CheckPointWriteRealDatasetHDF5, &
            CheckPointReadIntDatasetHDF5, &
            CheckpointOpenFileForReadHDF5, &
            CheckPointReadCompatibilityHDF5, &
#endif
            CheckpointPeriodicTimeWaypoints, &
            CheckpointInputRecord, &
            CheckpointRead

contains

! ************************************************************************** !

function CheckpointFilename(append_name, option)
  !
  ! This subroutine creates the filename of a checkpoint file without a suffix
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/30/15
  ! 

  use Option_module
  use String_module, only : StringNull

  character(len=MAXSTRINGLENGTH) :: append_name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: CheckpointFilename

  CheckpointFilename = trim(option%global_prefix) // &
                       trim(option%group_prefix) // &
                       trim(adjustl(append_name))

  CheckpointFilename = adjustl(CheckpointFilename)

end function CheckpointFilename

! ************************************************************************** !

function CheckpointAppendNameAtTime(checkpoint_option,time,option)
  !
  ! This subroutine forms the appendage to the checkpoint filename.
  !
  ! Author: Jenn Frederick
  ! Date: 1/29/2016
  ! 

  use Output_Aux_module
  use Units_module
  use Option_module

  implicit none

  type(checkpoint_option_type) :: checkpoint_option
  PetscReal :: time
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: CheckpointAppendNameAtTime
  character(len=MAXWORDLENGTH) :: word
  PetscReal :: temp_time

  ! time is actually option%time. do not overwrite it.
  temp_time = time * checkpoint_option%tconv
  !write(time_string,'(1pe12.4)') time
  write(word,'(f15.4)') temp_time
  CheckpointAppendNameAtTime = '-' // trim(adjustl(word)) // &
                             trim(adjustl(checkpoint_option%tunit))
    
end function CheckpointAppendNameAtTime

! ************************************************************************** !

function CheckpointAppendNameAtTimestep(checkpoint_option,timestep,option)
  !
  ! This subroutine forms the appendage to the checkpoint filename.
  !
  ! Author: Jenn Frederick
  ! Date: 1/29/2016
  ! 

  use Output_Aux_module
  use Units_module
  use Option_module

  implicit none

  type(checkpoint_option_type) :: checkpoint_option
  PetscInt :: timestep
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: CheckpointAppendNameAtTimestep
  character(len=MAXWORDLENGTH) :: word

  write(word,'(i9)') timestep
  CheckpointAppendNameAtTimestep = '-' // 'ts' // trim(adjustl(word))

end function CheckpointAppendNameAtTimestep

! ************************************************************************** !

subroutine CheckpointOpenFileForWriteBinary(viewer,append_name,option)
  ! 
  ! Opens checkpoint file; sets format
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/26/13
  ! 

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module

  implicit none

  PetscViewer :: viewer
  character(len=MAXSTRINGLENGTH) :: append_name
  type(option_type) :: option

  PetscErrorCode :: ierr
  character(len=MAXSTRINGLENGTH) :: filename

  filename = CheckpointFilename(append_name,option)
  filename = trim(filename) // '.chk'

  !geh: To skip .info file, need to split PetscViewerBinaryOpen() 
  !     into the routines it calls so that PetscViewerBinarySkipInfo()
  !     can be called after PetscViewerSetType(), but before
  !     PetscViewerFileSetName().  See note in PETSc docs.
  !call PetscViewerBinaryOpen(option%mycomm, filename, FILE_MODE_WRITE, &
  !                           viewer, ierr)
  call PetscViewerCreate(option%mycomm,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerSetType(viewer,PETSCVIEWERBINARY,ierr);CHKERRQ(ierr)
  call PetscViewerFileSetMode(viewer,FILE_MODE_WRITE,ierr);CHKERRQ(ierr)
  call PetscViewerBinarySkipInfo(viewer,ierr);CHKERRQ(ierr)
  call PetscViewerFileSetName(viewer,filename,ierr);CHKERRQ(ierr)
  
  write(option%io_buffer,'(" --> Dump checkpoint file: ", a64)') &
    trim(adjustl(filename))
  call printMsg(option)

end subroutine CheckpointOpenFileForWriteBinary

! ************************************************************************** !

subroutine CheckPointWriteCompatibilityBinary(viewer,option)
  ! 
  ! Writes a PetscBag holding the version number and the size of a
  ! complex extended class to ensure that the size of the class matches.
  ! The purpose of this test is to catch incompatibility.  
  !
  ! Technically, the BagSize should be 8 * the number of objects (int, real,
  ! etc.).  If we use 4 for PetscInt, the size is incorrect (due to padding
  ! in the OS???).  Anyway, using the following test sets a size sufficiently
  ! large:
  !
  ! see PETSC_DIR/src/sys/examples/tutorials/ex5f90.F90
  !
  ! class(whatever_type), pointer :: header
  ! type(whatever_type) :: dummy_header
  ! character(len=1),pointer :: dummy_char(:)
  ! PetscSizeT :: bagsize = size(transfer(dummy_header,dummy_char)) 
  ! 
  ! Author: Glenn Hammond
  ! Date: 003/26/15
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  
  implicit none

  PetscViewer :: viewer
  type(option_type) :: option

  type(checkpoint_header_type), pointer :: header
  type(checkpoint_header_type) :: dummy_header
  character(len=1),pointer :: dummy_char(:)
  PetscBag :: bag
  PetscSizeT :: bagsize
  PetscErrorCode :: ierr

  ! solely for test purposes here
  type(extended_test_header_type) :: test_header

  bagsize = size(transfer(dummy_header,dummy_char))

  call PetscBagCreate(option%mycomm,bagsize,bag,ierr);CHKERRQ(ierr)
  call PetscBagGetData(bag,header,ierr);CHKERRQ(ierr)
  call PetscBagRegisterInt(bag,header%version,0, &
                           "checkpoint_version","",ierr);CHKERRQ(ierr)
  call PetscBagRegisterInt(bag,header%test_header_size,0, &
                           "test_header_size","",ierr);CHKERRQ(ierr)
  header%version = CHECKPOINT_REVISION_NUMBER
  header%test_header_size = size(transfer(test_header,dummy_char))
  call PetscBagView(bag,viewer,ierr);CHKERRQ(ierr)
  call PetscBagDestroy(bag,ierr);CHKERRQ(ierr)

end subroutine CheckPointWriteCompatibilityBinary

! ************************************************************************** !

subroutine CheckPointReadCompatibilityBinary(viewer,option)
  ! 
  ! Reads in a PetscBag holding the version number and the size of a
  ! complex extended class to ensure that the size of the class matches.
  ! The purpose of this test is to catch incompatibility.  
  !
  ! Technically, the BagSize should be 8 * the number of objects (int, real,
  ! etc.).  If we use 4 for PetscInt, the size is incorrect (due to padding
  ! in the OS???).  Anyway, using the following test sets a size sufficiently
  ! large:
  !
  ! class(whatever_type), pointer :: header
  ! type(whatever_type) :: dummy_header
  ! character(len=1),pointer :: dummy_char(:)
  ! PetscSizeT :: bagsize = size(transfer(dummy_header,dummy_char)) 
  ! 
  ! Author: Glenn Hammond
  ! Date: 003/26/15
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  
  implicit none

  PetscViewer :: viewer
  type(option_type) :: option

  type(checkpoint_header_type), pointer :: header
  type(checkpoint_header_type) :: dummy_header
  character(len=1),pointer :: dummy_char(:)
  PetscBag :: bag
  PetscSizeT :: bagsize
  PetscErrorCode :: ierr
  character(len=MAXWORDLENGTH) :: word, word2
  PetscInt :: temp_int

  ! solely for test purposes here
  type(extended_test_header_type) :: test_header

  bagsize = size(transfer(dummy_header,dummy_char))

  call PetscBagCreate(option%mycomm,bagsize,bag,ierr);CHKERRQ(ierr)
  call PetscBagGetData(bag,header,ierr);CHKERRQ(ierr)
  call PetscBagRegisterInt(bag,header%version,0, &
                           "checkpoint_version","",ierr);CHKERRQ(ierr)
  call PetscBagRegisterInt(bag,header%test_header_size,0, &
                           "test_header_size","",ierr);CHKERRQ(ierr)
  call PetscBagLoad(viewer,bag,ierr);CHKERRQ(ierr)

  ! check compatibility
  if (header%version /= CHECKPOINT_REVISION_NUMBER) then
    write(word,*) header%version
    write(word2,*) CHECKPOINT_REVISION_NUMBER
    option%io_buffer = 'Incorrect checkpoint file format (' // &
      trim(adjustl(word)) // ' vs ' // &
      trim(adjustl(word2)) // ').'
    call printErrMsg(option)
  endif
  
  temp_int = size(transfer(test_header,dummy_char))
  if (header%test_header_size /= temp_int) then
    write(word,*) header%test_header_size
    write(word2,*) temp_int
    option%io_buffer = 'Inconsistent PetscBagSize (' // &
      trim(adjustl(word)) // ' vs ' // &
      trim(adjustl(word2)) // ').'
    call printErrMsg(option)
  endif

  call PetscBagDestroy(bag,ierr);CHKERRQ(ierr)

end subroutine CheckPointReadCompatibilityBinary

! ************************************************************************** !

subroutine CheckpointFlowProcessModelBinary(viewer,realization)
  ! 
  ! Checkpoints flow process model vectors
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/26/13
  ! 
#include "petsc/finclude/petscvec.h"
  use petscvec
  use Option_module
  use Realization_Subsurface_class
  use Field_module
  use Discretization_module
  use Grid_module
  use Global_module
  use Material_module
  use Variables_module, only : POROSITY, PERMEABILITY_X, PERMEABILITY_Y, &
                               PERMEABILITY_Z, STATE
  
  implicit none

  PetscViewer :: viewer
  class(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  Vec :: global_vec
  
  option => realization%option
  field => realization%field
  discretization => realization%discretization
  grid => realization%patch%grid
  
  global_vec = PETSC_NULL_VEC
  
  if (option%nflowdof > 0) then
    call DiscretizationCreateVector(realization%discretization,ONEDOF, &
                                    global_vec,GLOBAL,option)
    ! grid%flow_xx is the vector into which all of the primary variables are 
    ! packed for the SNESSolve().
    call VecView(field%flow_xx, viewer, ierr);CHKERRQ(ierr)


    ! If we are running with multiple phases, we need to dump the vector 
    ! that indicates what phases are present, as well as the 'var' vector 
    ! that holds variables derived from the primary ones via the translator.
    select case(option%iflowmode)
      case(TH_MODE)
        call DiscretizationLocalToGlobal(realization%discretization, &
                                         field%iphas_loc,global_vec,ONEDOF)
        call VecView(global_vec, viewer, ierr);CHKERRQ(ierr)
      case default
    end select 

    ! Porosity and permeability.
    ! (We only write diagonal terms of the permeability tensor for now, 
    ! since we have yet to add the full-tensor formulation.)
    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material, &
                                  field%work_loc,POROSITY,ZERO_INTEGER)
    call DiscretizationLocalToGlobal(discretization,field%work_loc, &
                                      global_vec,ONEDOF)
    call VecView(global_vec,viewer,ierr);CHKERRQ(ierr)
    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material, &
                                  field%work_loc,PERMEABILITY_X,ZERO_INTEGER)
    call DiscretizationLocalToGlobal(discretization,field%work_loc, &
                                      global_vec,ONEDOF)
    call VecView(global_vec,viewer,ierr);CHKERRQ(ierr)
    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material, &
                                  field%work_loc,PERMEABILITY_Y,ZERO_INTEGER)
    call DiscretizationLocalToGlobal(discretization,field%work_loc, &
                                      global_vec,ONEDOF)
    call VecView(global_vec,viewer,ierr);CHKERRQ(ierr)
    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material, &
                                  field%work_loc,PERMEABILITY_Z,ZERO_INTEGER)
    call DiscretizationLocalToGlobal(discretization,field%work_loc, &
                                      global_vec,ONEDOF)
    call VecView(global_vec,viewer,ierr);CHKERRQ(ierr)
  
  endif
  
  if (global_vec /= PETSC_NULL_VEC) then
    call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  endif  
  
end subroutine CheckpointFlowProcessModelBinary

! ************************************************************************** !

subroutine RestartFlowProcessModelBinary(viewer,realization)
  ! 
  ! Restarts flow process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/26/13
  ! 
#include "petsc/finclude/petscvec.h"
  use petscvec
      
  use Option_module
  use Realization_Subsurface_class
  use Field_module
  use Discretization_module
  use Grid_module
  use Global_module
  use Material_module
  use Variables_module, only : POROSITY, PERMEABILITY_X, PERMEABILITY_Y, &
                               PERMEABILITY_Z, STATE
  
  implicit none

  PetscViewer :: viewer
  class(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  Vec :: global_vec
  
  option => realization%option
  field => realization%field
  discretization => realization%discretization
  grid => realization%patch%grid
  
  global_vec = PETSC_NULL_VEC
  
  if (option%nflowdof > 0) then
    call DiscretizationCreateVector(realization%discretization,ONEDOF, &
                                    global_vec,GLOBAL,option)
  ! Load the PETSc vectors.
    call VecLoad(field%flow_xx,viewer,ierr);CHKERRQ(ierr)
    call DiscretizationGlobalToLocal(discretization,field%flow_xx, &
                                     field%flow_xx_loc,NFLOWDOF)
    call VecCopy(field%flow_xx,field%flow_yy,ierr);CHKERRQ(ierr)

    select case(option%iflowmode)
      case(TH_MODE)
        call VecLoad(global_vec,viewer,ierr);CHKERRQ(ierr)
        call DiscretizationGlobalToLocal(discretization,global_vec, &
                                         field%iphas_loc,ONEDOF)
        call VecCopy(field%iphas_loc,field%iphas_old_loc,ierr);CHKERRQ(ierr)
        call DiscretizationLocalToLocal(discretization,field%iphas_loc, &
                                        field%iphas_old_loc,ONEDOF)
      case default
    end select
    
    call VecLoad(global_vec,viewer,ierr);CHKERRQ(ierr)
    call DiscretizationGlobalToLocal(discretization,global_vec, &
                                      field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(realization%patch%aux%Material, &
                                  field%work_loc,POROSITY,ZERO_INTEGER)
    call VecLoad(global_vec,viewer,ierr);CHKERRQ(ierr)
    call DiscretizationGlobalToLocal(discretization,global_vec, &
                                      field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(realization%patch%aux%Material, &
                                  field%work_loc,PERMEABILITY_X,ZERO_INTEGER)
    call VecLoad(global_vec,viewer,ierr);CHKERRQ(ierr)
    call DiscretizationGlobalToLocal(discretization,global_vec, &
                                      field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(realization%patch%aux%Material, &
                                  field%work_loc,PERMEABILITY_Y,ZERO_INTEGER)
    call VecLoad(global_vec,viewer,ierr);CHKERRQ(ierr)
    call DiscretizationGlobalToLocal(discretization,global_vec, &
                                      field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(realization%patch%aux%Material, &
                                  field%work_loc,PERMEABILITY_Z,ZERO_INTEGER)
  endif
  
  if (global_vec /= PETSC_NULL_VEC) then
    call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  endif  
  
end subroutine RestartFlowProcessModelBinary

! ************************************************************************** !

#if defined(PETSC_HAVE_HDF5)
subroutine CheckpointOpenFileForWriteHDF5(file_id,grp_id,append_name,option, &
                                          id_stamp)
  !
  ! Opens checkpoint file; sets format
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/30/15
  !
  use Option_module
  use hdf5

  implicit none

  type(option_type) :: option
  character(len=MAXWORDLENGTH), optional, intent(in) :: id_stamp
  character(len=MAXSTRINGLENGTH) :: append_name
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: filename
  PetscErrorCode :: ierr
  PetscMPIInt :: hdf5_err

#if defined(SCORPIO_WRITE)
  integer, intent(out) :: file_id
  integer :: prop_id
  integer,intent(out) :: grp_id
#else
  integer(HID_T), intent(out) :: file_id
  integer(HID_T) :: prop_id
  integer(HID_T), intent(out) :: grp_id
#endif

  filename = CheckpointFilename(append_name, option)
  filename = trim(filename) // '.h5'

#if defined(SCORPIO_WRITE)
    filename = trim(filename) // CHAR(0)
    call scorpio_open_file(filename, option%iowrite_group_id, &
                              SCORPIO_FILE_CREATE, file_id, ierr)
#else

    ! initialize fortran interface
  call h5open_f(hdf5_err)

  call h5pcreate_f(H5P_FILE_ACCESS_F, prop_id, hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id, option%mycomm, MPI_INFO_NULL, hdf5_err)
#endif
  call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, hdf5_err, &
                   H5P_DEFAULT_F, prop_id)
  call h5pclose_f(prop_id, hdf5_err)

#endif

  string = "Checkpoint"
  call h5gcreate_f(file_id, string, grp_id, hdf5_err, OBJECT_NAMELEN_DEFAULT_F)

  write(option%io_buffer,'(" --> Dump checkpoint file: ", a64)') &
    trim(adjustl(filename))
  call printMsg(option)

end subroutine CheckpointOpenFileForWriteHDF5

! ************************************************************************** !

subroutine CheckpointOpenFileForReadHDF5(filename, file_id, grp_id, option)
  !
  ! Opens HDF5 checkpoint file for reading
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 08/09/15
  !
  use Option_module
  use hdf5

  implicit none

  character(len=MAXSTRINGLENGTH),intent(in) :: filename
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr
  PetscMPIInt :: hdf5_err

#if defined(SCORPIO)
  integer, intent(out) :: file_id
  integer :: prop_id
  integer,intent(out) :: grp_id
#else
  integer(HID_T), intent(out) :: file_id
  integer(HID_T) :: prop_id
  integer(HID_T), intent(out) :: grp_id
#endif

#if defined(SCORPIO)
  write(option%io_buffer, &
        '("Checkpoint from HDF5 not supported for SCORPIO. Darn.")')
  call printErrMsg(option)
#else

  ! initialize fortran interface
  call h5open_f(hdf5_err)

  call h5pcreate_f(H5P_FILE_ACCESS_F, prop_id, hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id, option%mycomm, MPI_INFO_NULL, hdf5_err)
#endif
  call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, hdf5_err, prop_id)
  call h5pclose_f(prop_id, hdf5_err)

  string = "Checkpoint"
  call h5gopen_f(file_id, string, grp_id, hdf5_err)
#endif

end subroutine CheckpointOpenFileForReadHDF5

! ************************************************************************** !

subroutine CheckPointWriteIntDatasetHDF5(chk_grp_id, dataset_name, dataset_rank, &
     dims, start, length, stride, data_int_array, option)
  !
  ! Within a HDF5 group (chk_grp_id), creates a new dataset (named dataset_name)
  ! and writes integer data type.
  !
  ! Author: Gautam Bisht
  ! Date: 07/30/15
  ! 
  use Option_module
  use hdf5
  use HDF5_module, only : trick_hdf5
  
  implicit none

#if defined(SCORPIO_WRITE)
  integer :: chk_grp_id
  PetscMPIInt :: dataset_rank
  integer, pointer :: dims(:)
  integer, pointer :: start(:)
  integer, pointer :: stride(:)
  integer, pointer :: length(:)
#else
  integer(HID_T) :: chk_grp_id
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscMPIInt :: dataset_rank
  integer(HSIZE_T), pointer :: dims(:)
  integer(HSIZE_T), pointer :: start(:)
  integer(HSIZE_T), pointer :: stride(:)
  integer(HSIZE_T), pointer :: length(:)
#endif
  type(option_type) :: option

  integer(HID_T) :: data_set_id
  integer(HID_T) :: grp_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: prop_id
  PetscErrorCode :: hdf5_err
  PetscErrorCode :: hdf5_flag
  PetscMPIInt, parameter :: ON=1, OFF=0

  ! must be 'integer' so that ibuffer does not switch to 64-bit integers 
  ! when PETSc is configured with --with-64-bit-indices=yes.
  integer, pointer :: data_int_array(:)

  call h5screate_simple_f(dataset_rank, dims, memory_space_id, hdf5_err, dims)

  dataset_name = trim(adjustl(dataset_name)) // CHAR(0)

  call h5eset_auto_f(OFF, hdf5_err)
  call h5dopen_f(chk_grp_id, dataset_name, data_set_id, hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON, hdf5_err)

  if (hdf5_flag < 0) then
    call h5pcreate_f(H5P_DATASET_CREATE_F, prop_id, hdf5_err)
    call h5screate_simple_f(dataset_rank, dims, grp_space_id, hdf5_err, dims)
    call h5dcreate_f(chk_grp_id, dataset_name, H5T_NATIVE_INTEGER, grp_space_id, &
                     data_set_id, hdf5_err, prop_id)
    call h5pclose_f(prop_id, hdf5_err)
  else
    call h5dget_space_f(data_set_id, grp_space_id, hdf5_err)
  endif

  call h5sselect_hyperslab_f(grp_space_id, H5S_SELECT_SET_F, start, length, &
                             hdf5_err, stride, stride)

  ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F, prop_id, hdf5_err)
#ifndef SERIAL_HDF5
  if (trick_hdf5) then
    call h5pset_dxpl_mpio_f(prop_id, H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
  else
    call h5pset_dxpl_mpio_f(prop_id, H5FD_MPIO_COLLECTIVE_F, &
                            hdf5_err)
  endif
#endif

  call h5dwrite_f(data_set_id, H5T_NATIVE_INTEGER, data_int_array, dims, &
                  hdf5_err, memory_space_id, grp_space_id, prop_id)

  call h5sclose_f(memory_space_id, hdf5_err)
  call h5sclose_f(grp_space_id, hdf5_err)
  call h5pclose_f(prop_id, hdf5_err)
  call h5dclose_f(data_set_id, hdf5_err)

end subroutine CheckPointWriteIntDatasetHDF5

! ************************************************************************** !

subroutine CheckPointWriteRealDatasetHDF5(chk_grp_id, dataset_name, dataset_rank, &
     dims, start, length, stride, data_real_array, option)
  !
  ! Within a HDF5 group (chk_grp_id), creates a new dataset (named dataset_name)
  ! and writes integer data type.
  !
  ! Author: Gautam Bisht
  ! Date: 07/30/15
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use hdf5
  use HDF5_module, only : trick_hdf5
  
  implicit none

#if defined(SCORPIO_WRITE)
  integer :: chk_grp_id
  PetscMPIInt :: dataset_rank
  integer, pointer :: dims(:)
  integer, pointer :: start(:)
  integer, pointer :: stride(:)
  integer, pointer :: length(:)
#else
  integer(HID_T) :: chk_grp_id
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscMPIInt :: dataset_rank
  integer(HSIZE_T), pointer :: dims(:)
  integer(HSIZE_T), pointer :: start(:)
  integer(HSIZE_T), pointer :: stride(:)
  integer(HSIZE_T), pointer :: length(:)
#endif
  type(option_type) :: option

  integer(HID_T) :: data_set_id
  integer(HID_T) :: grp_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: prop_id
  PetscErrorCode :: hdf5_err
  PetscErrorCode :: hdf5_flag
  PetscMPIInt, parameter :: ON=1, OFF=0

  PetscReal, pointer :: data_real_array(:)

  call h5screate_simple_f(dataset_rank, dims, memory_space_id, hdf5_err, dims)

  dataset_name = trim(adjustl(dataset_name)) // CHAR(0)

  call h5eset_auto_f(OFF, hdf5_err)
  call h5dopen_f(chk_grp_id, dataset_name, data_set_id, hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON, hdf5_err)

  if (hdf5_flag < 0) then
    call h5pcreate_f(H5P_DATASET_CREATE_F, prop_id, hdf5_err)
    call h5screate_simple_f(dataset_rank, dims, grp_space_id, hdf5_err, dims)
    call h5dcreate_f(chk_grp_id, dataset_name, H5T_NATIVE_DOUBLE, grp_space_id, &
                     data_set_id, hdf5_err, prop_id)
    call h5pclose_f(prop_id, hdf5_err)
  else
    call h5dget_space_f(data_set_id, grp_space_id, hdf5_err)
  endif

  call h5sselect_hyperslab_f(grp_space_id, H5S_SELECT_SET_F, start, length, &
                             hdf5_err, stride, stride)

  ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F, prop_id, hdf5_err)
#ifndef SERIAL_HDF5
  if (trick_hdf5) then
    call h5pset_dxpl_mpio_f(prop_id, H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
  else
    call h5pset_dxpl_mpio_f(prop_id, H5FD_MPIO_COLLECTIVE_F, &
                            hdf5_err)
  endif
#endif

  call h5dwrite_f(data_set_id, H5T_NATIVE_DOUBLE, data_real_array, dims, &
                  hdf5_err, memory_space_id, grp_space_id, prop_id)

  call h5sclose_f(memory_space_id, hdf5_err)
  call h5sclose_f(grp_space_id, hdf5_err)
  call h5pclose_f(prop_id, hdf5_err)
  call h5dclose_f(data_set_id, hdf5_err)

end subroutine CheckPointWriteRealDatasetHDF5

! ************************************************************************** !

subroutine CheckPointReadIntDatasetHDF5(chk_grp_id, dataset_name, dataset_rank, &
     dims, start, length, stride, data_int_array, option)
  !
  ! Within a HDF5 group (chk_grp_id), reads data from a dataset (named dataset_name)
  !
  ! Author: Gautam Bisht
  ! Date: 08/16/15
  ! 
  use Option_module
  use hdf5
  use HDF5_module, only : trick_hdf5
  
  implicit none


#if defined(SCORPIO_WRITE)
  integer :: chk_grp_id
  PetscMPIInt :: dataset_rank
  integer, pointer :: dims(:)
  integer, pointer :: start(:)
  integer, pointer :: stride(:)
  integer, pointer :: length(:)
#else
  integer(HID_T) :: chk_grp_id
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscMPIInt :: dataset_rank
  integer(HSIZE_T), pointer :: dims(:)
  integer(HSIZE_T), pointer :: start(:)
  integer(HSIZE_T), pointer :: stride(:)
  integer(HSIZE_T), pointer :: length(:)
#endif
  type(option_type) :: option

  integer(HID_T) :: data_set_id
  integer(HID_T) :: grp_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: prop_id
  PetscErrorCode :: hdf5_err
  PetscErrorCode :: hdf5_flag
  PetscMPIInt, parameter :: ON=1, OFF=0

  ! must be 'integer' so that ibuffer does not switch to 64-bit integers 
  ! when PETSc is configured with --with-64-bit-indices=yes.
  integer, pointer :: data_int_array(:)

  call h5screate_simple_f(dataset_rank, dims, memory_space_id, hdf5_err, dims)

  dataset_name = trim(adjustl(dataset_name)) // CHAR(0)

  call h5eset_auto_f(OFF, hdf5_err)
  call h5dopen_f(chk_grp_id, dataset_name, data_set_id, hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON, hdf5_err)

  call h5dget_space_f(data_set_id, grp_space_id, hdf5_err)

  call h5sselect_hyperslab_f(grp_space_id, H5S_SELECT_SET_F, start, length, &
                             hdf5_err, stride, stride)

  ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F, prop_id, hdf5_err)
#ifndef SERIAL_HDF5
  if (trick_hdf5) then
    call h5pset_dxpl_mpio_f(prop_id, H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
  else
    call h5pset_dxpl_mpio_f(prop_id, H5FD_MPIO_COLLECTIVE_F, &
                            hdf5_err)
  endif
#endif

  call h5dread_f(data_set_id, H5T_NATIVE_INTEGER, data_int_array, dims, &
                  hdf5_err, memory_space_id, grp_space_id, prop_id)

  call h5sclose_f(memory_space_id, hdf5_err)
  call h5sclose_f(grp_space_id, hdf5_err)
  call h5pclose_f(prop_id, hdf5_err)
  call h5dclose_f(data_set_id, hdf5_err)

end subroutine CheckPointReadIntDatasetHDF5

! ************************************************************************** !

subroutine CheckPointReadRealDatasetHDF5(chk_grp_id, dataset_name, dataset_rank, &
     dims, start, length, stride, data_real_array, option)
  !
  ! Within a HDF5 group (chk_grp_id), reads data from a dataset (named dataset_name)
  !
  ! Author: Gautam Bisht
  ! Date: 08/16/15
  ! 
  use Option_module
  use hdf5
  use HDF5_module, only : trick_hdf5
  
  implicit none


#if defined(SCORPIO_WRITE)
  integer :: chk_grp_id
  PetscMPIInt :: dataset_rank
  integer, pointer :: dims(:)
  integer, pointer :: start(:)
  integer, pointer :: stride(:)
  integer, pointer :: length(:)
#else
  integer(HID_T) :: chk_grp_id
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscMPIInt :: dataset_rank
  integer(HSIZE_T), pointer :: dims(:)
  integer(HSIZE_T), pointer :: start(:)
  integer(HSIZE_T), pointer :: stride(:)
  integer(HSIZE_T), pointer :: length(:)
#endif
  type(option_type) :: option

  integer(HID_T) :: data_set_id
  integer(HID_T) :: grp_space_id
  integer(HID_T) :: memory_space_id
  integer(HID_T) :: prop_id
  PetscErrorCode :: hdf5_err
  PetscErrorCode :: hdf5_flag
  PetscMPIInt, parameter :: ON=1, OFF=0

  PetscReal, pointer :: data_real_array(:)

  call h5screate_simple_f(dataset_rank, dims, memory_space_id, hdf5_err, dims)

  dataset_name = trim(adjustl(dataset_name)) // CHAR(0)

  call h5eset_auto_f(OFF, hdf5_err)
  call h5dopen_f(chk_grp_id, dataset_name, data_set_id, hdf5_err)
  hdf5_flag = hdf5_err
  call h5eset_auto_f(ON, hdf5_err)

  call h5dget_space_f(data_set_id, grp_space_id, hdf5_err)

  call h5sselect_hyperslab_f(grp_space_id, H5S_SELECT_SET_F, start, length, &
                             hdf5_err, stride, stride)

  ! write the data
  call h5pcreate_f(H5P_DATASET_XFER_F, prop_id, hdf5_err)
#ifndef SERIAL_HDF5
  if (trick_hdf5) then
    call h5pset_dxpl_mpio_f(prop_id, H5FD_MPIO_INDEPENDENT_F, &
                            hdf5_err)
  else
    call h5pset_dxpl_mpio_f(prop_id, H5FD_MPIO_COLLECTIVE_F, &
                            hdf5_err)
  endif
#endif

  call h5dread_f(data_set_id, H5T_NATIVE_DOUBLE, data_real_array, dims, &
                  hdf5_err, memory_space_id, grp_space_id, prop_id)

  call h5sclose_f(memory_space_id, hdf5_err)
  call h5sclose_f(grp_space_id, hdf5_err)
  call h5pclose_f(prop_id, hdf5_err)
  call h5dclose_f(data_set_id, hdf5_err)

end subroutine CheckPointReadRealDatasetHDF5

! ************************************************************************** !

subroutine CheckPointWriteCompatibilityHDF5(chk_grp_id, option)
  !
  ! Write the PFLOTRAN checkpoint version number. The purpose of this is to
  ! catch incompatibility.
  !
  ! Author: Gautam Bisht
  ! Date: 08/30/15
  !
  use Option_module
  use hdf5

  implicit none

#if defined(SCORPIO_WRITE)
  integer :: chk_grp_id
  integer, pointer :: dims(:)
  integer, pointer :: start(:)
  integer, pointer :: stride(:)
  integer, pointer :: length(:)
#else
  integer(HID_T) :: chk_grp_id
  integer(HSIZE_T), pointer :: dims(:)
  integer(HSIZE_T), pointer :: start(:)
  integer(HSIZE_T), pointer :: stride(:)
  integer(HSIZE_T), pointer :: length(:)
#endif
  type(option_type) :: option


  PetscMPIInt :: dataset_rank
  character(len=MAXSTRINGLENGTH) :: dataset_name
  ! must be 'integer' so that ibuffer does not switch to 64-bit integers 
  ! when PETSc is configured with --with-64-bit-indices=yes.
  integer, pointer :: int_array(:)

  dataset_name = "Revision Number" // CHAR(0)

  allocate(start(1))
  allocate(dims(1))
  allocate(length(1))
  allocate(stride(1))
  allocate(int_array(1))

  dataset_rank = 1
  dims(1) = ONE_INTEGER
  start(1) = 0
  length(1) = ONE_INTEGER
  stride(1) = ONE_INTEGER

  int_array(1) = CHECKPOINT_REVISION_NUMBER

  call CheckPointWriteIntDatasetHDF5(chk_grp_id, dataset_name, dataset_rank, &
     dims, start, length, stride, int_array, option)

  deallocate(start)
  deallocate(dims)
  deallocate(length)
  deallocate(stride)
  deallocate(int_array)

end subroutine CheckPointWriteCompatibilityHDF5

! ************************************************************************** !

subroutine CheckPointReadCompatibilityHDF5(chk_grp_id, option)
  !
  ! Reads the PFLOTRAN checkpoint version number. The purpose of this is to
  ! catch incompatibility.
  !
  ! Author: Gautam Bisht
  ! Date: 08/16/15
  !
  use Option_module
  use hdf5

  implicit none

#if defined(SCORPIO_WRITE)
  integer :: chk_grp_id
  integer, pointer :: dims(:)
  integer, pointer :: start(:)
  integer, pointer :: stride(:)
  integer, pointer :: length(:)
#else
  integer(HID_T) :: chk_grp_id
  integer(HSIZE_T), pointer :: dims(:)
  integer(HSIZE_T), pointer :: start(:)
  integer(HSIZE_T), pointer :: stride(:)
  integer(HSIZE_T), pointer :: length(:)
#endif
  type(option_type) :: option


  PetscMPIInt :: dataset_rank
  character(len=MAXSTRINGLENGTH) :: dataset_name
  ! must be 'integer' so that ibuffer does not switch to 64-bit integers 
  ! when PETSc is configured with --with-64-bit-indices=yes.
  integer, pointer :: int_array(:)
  character(len=MAXWORDLENGTH) :: word, word2

  dataset_name = "Revision Number" // CHAR(0)

  allocate(start(1))
  allocate(dims(1))
  allocate(length(1))
  allocate(stride(1))
  allocate(int_array(1))

  dataset_rank = 1
  dims(1) = ONE_INTEGER
  start(1) = 0
  length(1) = ONE_INTEGER
  stride(1) = ONE_INTEGER

  call CheckPointReadIntDatasetHDF5(chk_grp_id, dataset_name, dataset_rank, &
       dims, start, length, stride, int_array, option)
  
  if (int_array(1) /= CHECKPOINT_REVISION_NUMBER) then
    write(word,*) int_array(1)
    write(word2,*) CHECKPOINT_REVISION_NUMBER
    option%io_buffer = 'Incorrect checkpoint file format (' // &
      trim(adjustl(word)) // ' vs ' // &
      trim(adjustl(word2)) // ').'
    call printErrMsg(option)
  endif

  deallocate(start)
  deallocate(dims)
  deallocate(length)
  deallocate(stride)
  deallocate(int_array)

end subroutine CheckPointReadCompatibilityHDF5

! ************************************************************************** !

subroutine CheckpointFlowProcessModelHDF5(pm_grp_id, realization)
  !
  ! Checkpoints flow process model vectors
  !
  ! Author: Glenn Hammond
  ! Date: 07/26/13
  !
#include "petsc/finclude/petscvec.h"
  use petscvec
  use Option_module
  use Realization_Subsurface_class
  use Field_module
  use Discretization_module
  use Grid_module
  use Global_module
  use Material_module
  use Variables_module, only : POROSITY, PERMEABILITY_X, PERMEABILITY_Y, &
                               PERMEABILITY_Z, STATE
  use hdf5
  use HDF5_module, only : HDF5WriteDataSetFromVec
  implicit none

#if defined(SCORPIO_WRITE)
  integer :: pm_grp_id
#else
  integer(HID_T) :: pm_grp_id
#endif
  class(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  Vec :: global_vec
  Vec :: natural_vec
  character(len=MAXSTRINGLENGTH) :: dataset_name

  option => realization%option
  field => realization%field
  discretization => realization%discretization
  grid => realization%patch%grid

  global_vec = PETSC_NULL_VEC

  if (option%nflowdof > 0) then
     call DiscretizationCreateVector(realization%discretization, NFLOWDOF, &
                                    natural_vec, NATURAL, option)

    call DiscretizationGlobalToNatural(discretization, field%flow_xx, &
                                       natural_vec, NFLOWDOF)

    dataset_name = "Primary_Variables" // CHAR(0)
    call HDF5WriteDataSetFromVec(dataset_name, option, natural_vec, &
         pm_grp_id, H5T_NATIVE_DOUBLE)
    call VecDestroy(natural_vec, ierr);CHKERRQ(ierr)

    call DiscretizationCreateVector(realization%discretization, ONEDOF, &
                                    global_vec, GLOBAL,option)
     call DiscretizationCreateVector(realization%discretization, ONEDOF, &
                                    natural_vec, NATURAL, option)

    ! If we are running with multiple phases, we need to dump the vector
    ! that indicates what phases are present, as well as the 'var' vector
    ! that holds variables derived from the primary ones via the translator.
    select case(option%iflowmode)
      case(TH_MODE)
        call DiscretizationLocalToGlobal(realization%discretization, &
                                         field%iphas_loc,global_vec,ONEDOF)

        call DiscretizationGlobalToNatural(discretization, global_vec, &
                                           natural_vec, ONEDOF)
        dataset_name = "State" // CHAR(0)
        call HDF5WriteDataSetFromVec(dataset_name, option, natural_vec, &
            pm_grp_id, H5T_NATIVE_DOUBLE)
       case default
    end select

    ! Porosity and permeability.
    ! (We only write diagonal terms of the permeability tensor for now,
    ! since we have yet to add the full-tensor formulation.)
    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material, &
                                 field%work_loc,POROSITY,ZERO_INTEGER)
    call DiscretizationLocalToGlobal(discretization,field%work_loc, &
                                     global_vec,ONEDOF)
    call DiscretizationGlobalToNatural(discretization, global_vec, &
                                       natural_vec, ONEDOF)
    dataset_name = "Porosity" // CHAR(0)
    call HDF5WriteDataSetFromVec(dataset_name, option, natural_vec, &
                                             pm_grp_id, H5T_NATIVE_DOUBLE)

    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material, &
                                  field%work_loc,PERMEABILITY_X,ZERO_INTEGER)
    call DiscretizationLocalToGlobal(discretization,field%work_loc, &
                                      global_vec,ONEDOF)
    call DiscretizationGlobalToNatural(discretization, global_vec, &
                                       natural_vec, ONEDOF)
    dataset_name = "Permeability_X" // CHAR(0)
    call HDF5WriteDataSetFromVec(dataset_name, option, natural_vec, &
                                             pm_grp_id, H5T_NATIVE_DOUBLE)

    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material, &
                                  field%work_loc,PERMEABILITY_Y,ZERO_INTEGER)
    call DiscretizationLocalToGlobal(discretization,field%work_loc, &
                                      global_vec,ONEDOF)
    call DiscretizationGlobalToNatural(discretization, global_vec, &
                                       natural_vec, ONEDOF)
    dataset_name = "Permeability_Y" // CHAR(0)
    call HDF5WriteDataSetFromVec(dataset_name, option, natural_vec, &
                                             pm_grp_id, H5T_NATIVE_DOUBLE)

    call MaterialGetAuxVarVecLoc(realization%patch%aux%Material, &
                                  field%work_loc,PERMEABILITY_Z,ZERO_INTEGER)
    call DiscretizationLocalToGlobal(discretization,field%work_loc, &
                                      global_vec,ONEDOF)
    call DiscretizationGlobalToNatural(discretization, global_vec, &
                                       natural_vec, ONEDOF)
    dataset_name = "Permeability_Z" // CHAR(0)
    call HDF5WriteDataSetFromVec(dataset_name, option, natural_vec, &
                                             pm_grp_id, H5T_NATIVE_DOUBLE)

    call VecDestroy(global_vec, ierr);CHKERRQ(ierr)
    call VecDestroy(natural_vec, ierr);CHKERRQ(ierr)
  endif

end subroutine CheckpointFlowProcessModelHDF5

! ************************************************************************** !

subroutine RestartFlowProcessModelHDF5(pm_grp_id, realization)
  !
  ! Restarts flow process model vectors
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 08/16/2015
  !
#include "petsc/finclude/petscvec.h"
  use petscvec
  use Option_module
  use Realization_Subsurface_class
  use Field_module
  use Discretization_module
  use Grid_module
  use Global_module
  use Material_module
  use Variables_module, only : POROSITY, PERMEABILITY_X, PERMEABILITY_Y, &
                               PERMEABILITY_Z, STATE
  use hdf5
  use HDF5_module, only : HDF5ReadDataSetInVec
  implicit none

#if defined(SCORPIO_WRITE)
  integer :: pm_grp_id
#else
  integer(HID_T) :: pm_grp_id
#endif
  class(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  Vec :: global_vec
  Vec :: natural_vec
  character(len=MAXSTRINGLENGTH) :: dataset_name

  option => realization%option
  field => realization%field
  discretization => realization%discretization
  grid => realization%patch%grid

  global_vec = PETSC_NULL_VEC

  if (option%nflowdof > 0) then
    call DiscretizationCreateVector(realization%discretization, NFLOWDOF, &
                                    natural_vec, NATURAL, option)

    dataset_name = "Primary_Variables" // CHAR(0)
    call HDF5ReadDataSetInVec(dataset_name, option, natural_vec, &
         pm_grp_id, H5T_NATIVE_DOUBLE)

    call DiscretizationNaturalToGlobal(discretization, natural_vec, field%flow_xx, &
                                       NFLOWDOF)
    call DiscretizationGlobalToLocal(discretization,field%flow_xx, &
                                     field%flow_xx_loc,NFLOWDOF)
    call VecCopy(field%flow_xx,field%flow_yy,ierr);CHKERRQ(ierr)
    
    call VecDestroy(natural_vec, ierr);CHKERRQ(ierr)

    call DiscretizationCreateVector(realization%discretization, ONEDOF, &
                                    global_vec, GLOBAL,option)
    call DiscretizationCreateVector(realization%discretization, ONEDOF, &
                                    natural_vec, NATURAL, option)

    ! If we are running with multiple phases, we need to dump the vector
    ! that indicates what phases are present, as well as the 'var' vector
    ! that holds variables derived from the primary ones via the translator.
    dataset_name = "State" // CHAR(0)
    select case(option%iflowmode)
      case(TH_MODE)
        call HDF5ReadDataSetInVec(dataset_name, option, natural_vec, &
             pm_grp_id, H5T_NATIVE_DOUBLE)
        call DiscretizationNaturalToGlobal(discretization, natural_vec, &
                                           global_vec, ONEDOF)
        call DiscretizationGlobalToLocal(realization%discretization, &
                                         global_vec, field%iphas_loc, ONEDOF)
        call VecCopy(field%iphas_loc,field%iphas_old_loc,ierr);CHKERRQ(ierr)
        call DiscretizationLocalToLocal(discretization,field%iphas_loc, &
                                        field%iphas_old_loc,ONEDOF)
     case default
    end select

    ! Porosity and permeability.
    ! (We only write diagonal terms of the permeability tensor for now,
    ! since we have yet to add the full-tensor formulation.)
    dataset_name = "Porosity" // CHAR(0)
    call HDF5ReadDataSetInVec(dataset_name, option, natural_vec, &
                              pm_grp_id, H5T_NATIVE_DOUBLE)
    call DiscretizationNaturalToGlobal(discretization, natural_vec, global_vec, &
                                       ONEDOF)
    call DiscretizationGlobalToLocal(discretization, global_vec, field%work_loc, &
                                     ONEDOF)
    call MaterialSetAuxVarVecLoc(realization%patch%aux%Material, &
                                 field%work_loc,POROSITY,ZERO_INTEGER)

    dataset_name = "Permeability_X" // CHAR(0)
    call HDF5ReadDataSetInVec(dataset_name, option, natural_vec, &
                              pm_grp_id, H5T_NATIVE_DOUBLE)
    call DiscretizationNaturalToGlobal(discretization, natural_vec, global_vec, &
                                       ONEDOF)
    call DiscretizationGlobalToLocal(discretization, global_vec, field%work_loc, &
                                     ONEDOF)
    call MaterialSetAuxVarVecLoc(realization%patch%aux%Material, &
                                 field%work_loc,PERMEABILITY_X,ZERO_INTEGER)

    dataset_name = "Permeability_Y" // CHAR(0)
    call HDF5ReadDataSetInVec(dataset_name, option, natural_vec, &
                              pm_grp_id, H5T_NATIVE_DOUBLE)
    call DiscretizationNaturalToGlobal(discretization, natural_vec, global_vec, &
                                       ONEDOF)
    call DiscretizationGlobalToLocal(discretization, global_vec, field%work_loc, &
                                     ONEDOF)
    call MaterialSetAuxVarVecLoc(realization%patch%aux%Material, &
                                 field%work_loc,PERMEABILITY_Y,ZERO_INTEGER)

    dataset_name = "Permeability_Z" // CHAR(0)
    call HDF5ReadDataSetInVec(dataset_name, option, natural_vec, &
                              pm_grp_id, H5T_NATIVE_DOUBLE)
    call DiscretizationNaturalToGlobal(discretization, natural_vec, global_vec, &
                                       ONEDOF)
    call DiscretizationGlobalToLocal(discretization, global_vec, field%work_loc, &
                                     ONEDOF)
    call MaterialSetAuxVarVecLoc(realization%patch%aux%Material, &
                                 field%work_loc,PERMEABILITY_Z,ZERO_INTEGER)

    call VecDestroy(global_vec, ierr);CHKERRQ(ierr)
    call VecDestroy(natural_vec, ierr);CHKERRQ(ierr)
  endif

end subroutine RestartFlowProcessModelHDF5
#endif

! ************************************************************************** !

subroutine CheckpointRead(input,option,checkpoint_option,waypoint_list)
  ! 
  ! Reads the CHECKPOINT card in an input file.
  ! 
  ! Author: Jenn Frederick
  ! Date: 01/29/2016
  !  

  use Option_module
  use Input_Aux_module
  use Output_Aux_module
  use Waypoint_module
  use String_module
  use Units_module

  implicit none

  type(input_type),pointer :: input
  type(option_type) :: option
  type(checkpoint_option_type), pointer :: checkpoint_option
  type(waypoint_list_type) :: waypoint_list
  
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: card
  character(len=MAXSTRINGLENGTH) :: temp_string
  character(len=MAXWORDLENGTH) :: internal_units
  character(len=MAXWORDLENGTH) :: default_time_units
  type(waypoint_type), pointer :: waypoint
  PetscReal :: units_conversion
  PetscReal :: temp_real
  PetscReal, pointer :: temp_real_array(:)
  PetscInt :: i
  PetscBool :: format_binary
  PetscBool :: format_hdf5

  if (.not.associated(checkpoint_option)) then
    checkpoint_option => CheckpointOptionCreate()
  endif
  
  format_binary = PETSC_FALSE
  format_hdf5 = PETSC_FALSE
  default_time_units = ''
  do
    call InputReadPflotranString(input,option)
    call InputReadStringErrorMsg(input,option,'CHECKPOINT')
    if (InputCheckExit(input,option)) exit
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'checkpoint option or value', &
                        'CHECKPOINT')
    call StringToUpper(word)
    select case(trim(word))
      case ('PERIODIC')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'time increment', &
                            'CHECKPOINT,PERIODIC')
        select case(trim(word))
          case('TIME')
            call InputReadDouble(input,option,temp_real)
            call InputErrorMsg(input,option,'time increment', &
                                'CHECKPOINT,PERIODIC,TIME')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'time increment units', &
                                'CHECKPOINT,PERIODIC,TIME')
            internal_units = 'sec'
            units_conversion = UnitsConvertToInternal(word, &
                                internal_units,option)
            checkpoint_option%tconv = 1.d0/units_conversion
            checkpoint_option%tunit = trim(word)
            checkpoint_option%periodic_time_incr = temp_real*units_conversion
          case('TIMESTEP')
            call InputReadInt(input,option,checkpoint_option%periodic_ts_incr)
            call InputErrorMsg(input,option,'timestep increment', &
                                'CHECKPOINT,PERIODIC,TIMESTEP')
          case default
            call InputKeywordUnrecognized(word,'CHECKPOINT,PERIODIC', &
                                          option)
        end select
      case ('TIMES')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'time units', &
                            'CHECKPOINT,TIMES')
        internal_units = 'sec'
        units_conversion = UnitsConvertToInternal(word,internal_units, &
                                                  option)
        checkpoint_option%tconv = 1.d0/units_conversion
        checkpoint_option%tunit = trim(word)
!geh: this needs to be tested.
#if 0
        temp_string = 'CHECKPOINT,TIMES'
        nullify(temp_real_array)
        call UtilityReadArray(temp_real_array,NEG_ONE_INTEGER, &
                              temp_string,input,option)
        do i = 1, size(temp_real_array)
          waypoint => WaypointCreate()
          waypoint%time = temp_real_array(i)*units_conversion
          waypoint%print_checkpoint = PETSC_TRUE
          call WaypointInsertInList(waypoint,waypoint_list)
        enddo
        call DeallocateArray(temp_real_array)
#else
        do
          call InputReadDouble(input,option,temp_real)
          if (input%ierr /= 0) exit
          call InputErrorMsg(input,option,'checkpoint time', &
                              'CHECKPOINT,TIMES') 
          waypoint => WaypointCreate()
          waypoint%time = temp_real * units_conversion
          waypoint%print_checkpoint = PETSC_TRUE
          call WaypointInsertInList(waypoint,waypoint_list)     
        enddo
#endif
      case ('FORMAT')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'format type', &
                            'CHECKPOINT,FORMAT')
        call StringToUpper(word)
        select case(trim(word))
          case('BINARY')
            format_binary = PETSC_TRUE
          case('HDF5')
            format_hdf5 = PETSC_TRUE
          case default
            call InputKeywordUnrecognized(word,'CHECKPOINT,FORMAT', &
                                          option)
        end select
      case ('TIME_UNITS')
        call InputReadWord(input,option,default_time_units,PETSC_TRUE)
        call InputErrorMsg(input,option,'time units','CHECKPOINT')
      case default
        temp_string = 'Must specify PERIODIC TIME, PERIODIC TIMESTEP, &
                      &TIMES, or FORMAT'
        call InputKeywordUnrecognized(word,'CHECKPOINT',temp_string,option)
    end select
  enddo
  if (len_trim(default_time_units) > 0) then
    internal_units = 'sec'
    units_conversion = UnitsConvertToInternal(default_time_units, &
                                              internal_units,option)
    checkpoint_option%tconv = 1.d0/units_conversion
    checkpoint_option%tunit = trim(default_time_units)
  endif
  if (format_binary .and. format_hdf5) then
    checkpoint_option%format = CHECKPOINT_BOTH
  else if (format_hdf5) then
    checkpoint_option%format = CHECKPOINT_HDF5
  else ! default
    checkpoint_option%format = CHECKPOINT_BINARY
  endif
  
end subroutine CheckpointRead

! ************************************************************************** !

subroutine CheckpointPeriodicTimeWaypoints(checkpoint_option,waypoint_list)
  ! 
  ! Inserts periodic time waypoints into list
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/03/16
  !  

  use Option_module
  use Waypoint_module
  use Output_Aux_module
  use Utility_module

  implicit none

  type(option_type) :: option
  type(checkpoint_option_type), pointer :: checkpoint_option
  type(waypoint_list_type) :: waypoint_list
  type(waypoint_type), pointer :: waypoint
  character(len=MAXWORDLENGTH) :: word
  PetscReal :: final_time
  PetscReal :: temp_real
  PetscReal :: num_waypoints, warning_num_waypoints
  PetscInt :: k
  
  final_time = WaypointListGetFinalTime(waypoint_list)
  warning_num_waypoints = 15000.0

  if (final_time < 1.d-40) then
    option%io_buffer = 'No final time specified in waypoint list. &
      &Send your input deck to pflotran-dev.'
    call printMsg(option)
  endif
  
  ! add waypoints for periodic checkpoint
  if (associated(checkpoint_option)) then
    if (Initialized(checkpoint_option%periodic_time_incr)) then
      temp_real = 0.d0
      num_waypoints = final_time / checkpoint_option%periodic_time_incr
      if ((num_waypoints > warning_num_waypoints) .and. &
          OptionPrintToScreen(option)) then
        write(word,*) floor(num_waypoints)
        write(*,*) 'WARNING: Large number (' // trim(adjustl(word)) // &
                   ') of periodic checkpoints requested.'
        write(*,'(a68)',advance='no') '         Creating periodic checkpoint &
                                      &waypoints . . . Progress: 0%-'
      endif
      k = 0
      do
        k = k + 1
        temp_real = temp_real + checkpoint_option%periodic_time_incr
        if (temp_real > final_time) exit
        waypoint => WaypointCreate()
        waypoint%time = temp_real
        waypoint%print_checkpoint = PETSC_TRUE
        call WaypointInsertInList(waypoint,waypoint_list)
        if ((num_waypoints > warning_num_waypoints) .and. &
            OptionPrintToScreen(option)) then
          call PrintProgressBarInt(num_waypoints,TEN_INTEGER,k)
        endif
      enddo
    endif
  endif

end subroutine CheckpointPeriodicTimeWaypoints
  
! ************************************************************************** !

subroutine CheckpointInputRecord(checkpoint_option,waypoint_list)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/17/2016
  !  
  use Output_Aux_module
  use Waypoint_module

  implicit none

  type(checkpoint_option_type), pointer :: checkpoint_option
  type(waypoint_list_type), pointer :: waypoint_list
  
  type(waypoint_type), pointer :: cur_waypoint
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: checkpoints_found
  PetscInt :: id = INPUT_RECORD_UNIT

  write(id,'(a)') ' '
    write(id,'(a)') '---------------------------------------------------------&
                    &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'CHECKPOINTS'

  if (associated(checkpoint_option)) then
    write(id,'(a29)',advance='no') 'periodic timestep: '
    if (checkpoint_option%periodic_ts_incr == 0) then
      write(id,'(a)') 'OFF'
    else
      write(id,'(a)') 'ON'
      write(id,'(a29)',advance='no') 'timestep increment: '
      write(word,*) checkpoint_option%periodic_ts_incr
      write(id,'(a)') adjustl(trim(word))
    endif

    write(id,'(a29)',advance='no') 'periodic time: '
    if (checkpoint_option%periodic_time_incr <= 0) then
      write(id,'(a)') 'OFF'
    else
      write(id,'(a)') 'ON'
      write(id,'(a29)',advance='no') 'time increment: '
      write(word,*) checkpoint_option%periodic_time_incr * &
                    checkpoint_option%tconv
      write(id,'(a)') adjustl(trim(word)) // &
                      adjustl(trim(checkpoint_option%tunit))
    endif
  endif

  string = ''
  checkpoints_found = PETSC_FALSE
  write(id,'(a29)',advance='no') 'specific times: '
  cur_waypoint => waypoint_list%first
  do
    if (.not.associated(cur_waypoint)) exit
    if (cur_waypoint%print_checkpoint) then
      checkpoints_found = PETSC_TRUE
      write(word,*) cur_waypoint%time*checkpoint_option%tconv
      string = trim(string) // adjustl(trim(word)) // ','
    endif
    cur_waypoint => cur_waypoint%next
  enddo
  if (checkpoints_found) then
    write(id,'(a)') 'ON'
    write(id,'(a29)',advance='no') 'times (' // &
                                    trim(checkpoint_option%tunit) // '): '
    write(id,'(a)') trim(string)
  else
    write(id,'(a)') 'OFF'
  endif
  
end subroutine CheckpointInputRecord

end module Checkpoint_module
