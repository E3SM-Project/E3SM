! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=============================================================================================
 module mpas_atmphys_utilities

 use mpas_kind_types
 implicit none
 private
 public:: physics_error_fatal, &
          physics_message

 character(len=StrKIND),public:: mpas_err_message

 contains

!=============================================================================================
 subroutine physics_message(str)
!=============================================================================================

!input arguments:
 character(len=*),intent(in):: str

!---------------------------------------------------------------------------------------------

!#ifdef _MPI
! write(0,*) trim(str)
!#endif
!print*,trim(str)
write(0,*) trim(str)

 end subroutine physics_message

!=============================================================================================
 subroutine physics_error_fatal(str)
!=============================================================================================

!input arguments:
 character(len=*),intent(in):: str

!---------------------------------------------------------------------------------------------

 write(0,*)
 write(0,*) ( '------------------------------ FATAL CALLED ------------------------------')
 write(0,*) trim(str)
 stop ' MPAS core_physics abort'
 
 end subroutine physics_error_fatal

!=============================================================================================
 end module mpas_atmphys_utilities
!=============================================================================================
