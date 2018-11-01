
module PFLOTRAN_Provenance_module

!
!-_-! write_warning_comment !-_-!
! 
! IMPORTANT NOTE: This module should have no dependencies on other modules!!!
!
! Purpose: This file is a template / dummy file for
! pflotran_provenance.F90. pflotran-provenance.py automatically
! generantes pflotran_provenance.F90 from pflotran_no_provenance.F90.
! If you can not run pflotran-provenance.py to generate the provenance
! information, then modify your build system to compile and link this
! file instead of pflotran_provenance.F90
!

#include "petsc/finclude/petscsys.h"
  use petscsys

  implicit none

  public

  ! Size of provenance information
  PetscInt, parameter :: provenance_max_str_len = 7


  ! PFLOTRAN provenance information
  character(len=*), parameter :: pflotran_compile_date_time = "unknown"
  character(len=*), parameter :: pflotran_compile_user = "unknown"
  character(len=*), parameter :: pflotran_compile_hostname = "unknown"
  character(len=*), parameter :: pflotran_changeset = "unknown"
  character(len=*), parameter :: pflotran_status = "unknown"

  PetscInt, parameter :: detail_pflotran_fflags_len = 1
  character(len=provenance_max_str_len), parameter :: &
       detail_pflotran_fflags(detail_pflotran_fflags_len) = (/ "unknown" /)

  PetscInt, parameter :: detail_pflotran_status_len = 1
  character(len=provenance_max_str_len), parameter :: &
       detail_pflotran_status(detail_pflotran_status_len) = (/ "unknown" /)

  PetscInt, parameter :: detail_pflotran_parent_len = 1
  character(len=provenance_max_str_len), parameter :: &
       detail_pflotran_parent(detail_pflotran_parent_len) = (/ "unknown" /)

  ! FIXME(bja, 2013-11-25): break gcc when diffs are present
  PetscInt, parameter :: detail_pflotran_diff_len = 1
  character(len=provenance_max_str_len), parameter :: &
       detail_pflotran_diff(detail_pflotran_diff_len) = (/ "unknown" /)

  ! PETSc provenance information
  character(len=*), parameter :: petsc_status = "unknown"
  character(len=*), parameter :: petsc_changeset = "unknown"

  PetscInt, parameter :: detail_petsc_status_len = 1
  character(len=provenance_max_str_len), parameter :: &
       detail_petsc_status(detail_petsc_status_len) = (/ "unknown" /)

  PetscInt, parameter :: detail_petsc_parent_len = 1
  character(len=provenance_max_str_len), parameter :: &
       detail_petsc_parent(detail_petsc_parent_len) = (/ "unknown" /)

  PetscInt, parameter :: detail_petsc_config_len = 1
  character(len=provenance_max_str_len), parameter :: &
       detail_petsc_config(detail_petsc_config_len) = (/ "unknown" /)

  public :: PrintProvenanceToScreen

contains

subroutine PrintProvenanceToScreen()

  implicit none

  PetscInt i

  write(*, '(''------------------------------ Provenance --------------------------------------'')')

  write(*, '(''pflotran_compile_date_time = '', a)') pflotran_compile_date_time
  write(*, '(''pflotran_compile_user = '', a)') pflotran_compile_user
  write(*, '(''pflotran_compile_hostname = '', a)') pflotran_compile_hostname
  write(*, '(''pflotran_changeset = '', a)') pflotran_changeset
  write(*, '(''pflotran_status = '', a)') pflotran_status
  write(*, '(''petsc_changeset = '', a)') petsc_changeset
  write(*, '(''petsc_status = '', a)') petsc_status

!-_-! write_provenance_details !-_-!

  write(*, '(''--------------------------------------------------------------------------------'')')

end subroutine PrintProvenanceToScreen

! ************************************************************************** !

end module PFLOTRAN_Provenance_module
