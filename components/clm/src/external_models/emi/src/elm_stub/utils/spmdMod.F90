
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

  ! Default settings valid even if there is no spmd

  logical, public :: masterproc      ! proc 0 logical for printing msgs
  integer, public :: mpicom
  integer, public :: iam
  !
  ! Public methods
  !
  public :: spmd_init                ! Initialization

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
  iam = 0
  end subroutine spmd_init

end module spmdMod
