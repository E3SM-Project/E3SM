! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_utilities

 use mpas_kind_types

 implicit none
 private
 public:: physics_error_fatal, &
          physics_message

 character(len=StrKIND),public:: mpas_err_message

!Simple utility subroutines to broadcast messages or abort physics parameterizations.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.


 contains


!=================================================================================================================
 subroutine physics_message(str)
!=================================================================================================================

 use mpas_log, only : mpas_log_write

!input arguments:
 character(len=*),intent(in):: str

!-----------------------------------------------------------------------------------------------------------------

 call mpas_log_write(trim(str))

 end subroutine physics_message

!=================================================================================================================
 subroutine physics_error_fatal(str)
!=================================================================================================================

 use mpas_derived_types, only : MPAS_LOG_ERR, MPAS_LOG_CRIT
 use mpas_log, only : mpas_log_write

!input arguments:
 character(len=*),intent(in):: str

!--------------------------------------------------------------------------------------------------

 call mpas_log_write(' ',                                                                          messageType=MPAS_LOG_ERR)
 call mpas_log_write('------------------------------ FATAL CALLED ------------------------------', messageType=MPAS_LOG_ERR)
 call mpas_log_write(trim(str),                                                                    messageType=MPAS_LOG_ERR)
 call mpas_log_write('MPAS core_physics abort', messageType=MPAS_LOG_CRIT)
 
 end subroutine physics_error_fatal

!=================================================================================================================
 end module mpas_atmphys_utilities
!=================================================================================================================
