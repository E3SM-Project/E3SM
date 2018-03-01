
module spmdMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: spmdMod
!
! !DESCRIPTION:
! SPMD initialization
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varctl  , only: iulog
  implicit none

  private

#include "mpif.h"


  ! Default settings valid even if there is no spmd

  logical, public :: masterproc      ! proc 0 logical for printing msgs
  integer, public :: mpicom
  !
  ! Public methods
  !
  public :: spmd_init                ! Initialization

  public :: MPI_INTEGER
  public :: MPI_REAL8
  public :: MPI_LOGICAL
  public :: MPI_SUM
  public :: MPI_MIN
  public :: MPI_MAX
  public :: MPI_LOR
  public :: MPI_STATUS_SIZE
  public :: MPI_ANY_SOURCE
  public :: MPI_CHARACTER
  public :: MPI_COMM_WORLD
  public :: MPI_MAX_PROCESSOR_NAME

  !


contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: spmd_init( clm_mpicom )
!
! !INTERFACE:
  subroutine spmd_init(  )
!
! !DESCRIPTION:
! MPI initialization (number of cpus, processes, tids, etc)
!
! !USES
!
! !ARGUMENTS:
    implicit none

!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------


  masterproc = .true.
  end subroutine spmd_init

end module spmdMod
