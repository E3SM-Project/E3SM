module Wrapper_Hydrogeophysics_module
 
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  public :: HydrogeophysicsWrapperInit, &
            HydrogeophysicsWrapperStart, &
            HydrogeophysicsWrapperStep, &
            HydrogeophysicsWrapperStop, &
            HydrogeophysicsWrapperDestroy

contains

! ************************************************************************** !

subroutine HydrogeophysicsWrapperInit(option, &
                                      pflotran_tracer_vec_mpi_, &
                                      pflotran_tracer_vec_seq_, &
                                      pflotran_saturation_vec_mpi_, &
                                      pflotran_saturation_vec_seq_, &
                                      pflotran_temperature_vec_mpi_, &
                                      pflotran_temperature_vec_seq_, &
                                      pflotran_scatter_, &
                                      pf_e4d_master_comm)
  ! 
  ! Initializes the hydrogeophysics module
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/02/13
  ! 
  
  use Option_module

  use e4d_vars, only : E4D_COMM, my_rank, n_rank, PFE4D_MASTER_COMM, &
               pflotran_tracer_vec_mpi, pflotran_tracer_vec_seq, &
               pflotran_saturation_vec_mpi, pflotran_saturation_vec_seq, &
               pflotran_temperature_vec_mpi, pflotran_temperature_vec_seq, &
               pflotran_scatter, pflotran_vec_size, &
               pflotran_group_prefix
  use e4d_setup, only : setup_e4d, destroy_e4d
  use e4d_run, only: run_e4D
  
  implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  type(option_type) :: option
  Vec :: pflotran_tracer_vec_mpi_
  Vec :: pflotran_tracer_vec_seq_
  Vec :: pflotran_saturation_vec_mpi_
  Vec :: pflotran_saturation_vec_seq_
  Vec :: pflotran_temperature_vec_mpi_
  Vec :: pflotran_temperature_vec_seq_
  VecScatter :: pflotran_scatter_
  PetscMPIInt :: pf_e4d_master_comm
  PetscErrorCode :: ierr
  
!  call printMsg(option,'HydrogeophysicsWrapperInit()')
  
  ! from e4d_vars.F90
  E4D_COMM = option%mycomm
  PFE4D_MASTER_COMM = pf_e4d_master_comm
  my_rank = option%myrank
  n_rank = option%mycommsize
  pflotran_tracer_vec_mpi = pflotran_tracer_vec_mpi_
  pflotran_tracer_vec_seq = pflotran_tracer_vec_seq_
  pflotran_saturation_vec_mpi = pflotran_saturation_vec_mpi_
  pflotran_saturation_vec_seq = pflotran_saturation_vec_seq_
  pflotran_temperature_vec_mpi = pflotran_temperature_vec_mpi_
  pflotran_temperature_vec_seq = pflotran_temperature_vec_seq_
  pflotran_scatter = pflotran_scatter_
  pflotran_group_prefix = option%group_prefix
  ! pflotran_tracer_vec_seq only defined on master E4D process
  if (pflotran_tracer_vec_seq > 0) then
    call VecGetSize(pflotran_tracer_vec_seq,pflotran_vec_size, &
                    ierr);CHKERRQ(ierr)
    ! don't need saturation size as it is identical   
  endif

  call setup_e4d
  call run_e4d
  call destroy_e4d
  
end subroutine HydrogeophysicsWrapperInit

! ************************************************************************** !

subroutine HydrogeophysicsWrapperStart(option)
  ! 
  ! Starts the hydrogeophysics forward simulation
  ! loop
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/02/13
  ! 
  
  use Option_module

  implicit none

  type(option_type) :: option
  
!  call printMsg(option,'HydrogeophysicsWrapperStart()')

end subroutine HydrogeophysicsWrapperStart

! ************************************************************************** !

subroutine HydrogeophysicsWrapperStep(time, &
                                      tracer_mpi, &
                                      tracer_seq, &
                                      saturation_mpi, &
                                      saturation_seq, &
                                      temperature_mpi, &
                                      temperature_seq, &
                                      scatter,comm,option)
  ! 
  ! Performs a forward simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/02/13
  ! 
  use Option_module

  implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  PetscReal :: time
  Vec :: tracer_mpi
  Vec :: tracer_seq
  Vec :: saturation_mpi
  Vec :: saturation_seq
  Vec :: temperature_mpi
  Vec :: temperature_seq
  VecScatter :: scatter
  PetscMPIInt :: comm
  type(option_type) :: option
  
  PetscReal :: temp_time
  PetscErrorCode :: ierr
  
!  call printMsg(option,'HydrogeophysicsWrapperStep()')
  
  ! Bcasting a 1 to E4D tells it to run a forward simulation 
  if (comm /= MPI_COMM_NULL) then
    call MPI_Bcast(ONE_INTEGER_MPI,ONE_INTEGER_MPI,MPI_INTEGER, &
                   ZERO_INTEGER_MPI,comm,ierr)
    temp_time = time
    call MPI_Bcast(temp_time,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                   ZERO_INTEGER_MPI,comm,ierr)
  endif
  call VecScatterBegin(scatter,tracer_mpi,tracer_seq, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(scatter,tracer_mpi,tracer_seq, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterBegin(scatter,saturation_mpi,saturation_seq, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(scatter,saturation_mpi,saturation_seq, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  if (temperature_mpi /= 0) then
    call VecScatterBegin(scatter,temperature_mpi,temperature_seq, &
                         INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
    call VecScatterEnd(scatter,temperature_mpi,temperature_seq, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  endif
  
end subroutine HydrogeophysicsWrapperStep

! ************************************************************************** !

subroutine HydrogeophysicsWrapperStop(option,comm)
  ! 
  ! Stops the hydrogeophysics forward simulation
  ! loop
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/02/13
  ! 
  
  use Option_module

  implicit none

  type(option_type) :: option
  PetscMPIInt :: comm
  PetscErrorCode :: ierr
  
!  call printMsg(option,'HydrogeophysicsWrapperStop()')
  
  ! Bcasting a 0 to E4D tells it to stop
  call MPI_Bcast(ZERO_INTEGER_MPI,ONE_INTEGER_MPI,MPI_INTEGER, &
                 ZERO_INTEGER_MPI,comm,ierr)

end subroutine HydrogeophysicsWrapperStop

! ************************************************************************** !

recursive subroutine HydrogeophysicsWrapperDestroy(option)
  ! 
  ! Destroys the contents of the hydrogeophysics
  ! module
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/02/13
  ! 

  use Option_module

  implicit none
  
  type(option_type) :: option
  
!  call printMsg(option,'HydrogeophysicsWrapperDestroy()')

end subroutine HydrogeophysicsWrapperDestroy

end module Wrapper_Hydrogeophysics_module
