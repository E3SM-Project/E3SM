!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glimmer_log.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

!> module providing file logging and error/message handling
!!
!! Six levels of message/error are defined:
!! - Diagnostic messages
!! - Timestep enumeration and related information
!! - Information messages
!! - Warning messages
!! - Error messages
!! - Fatal error messages
!!
!! These are numbered 1--6, with increasing severity, and the level of
!! message output may be set to output all messages, only those above a particular 
!! severity, or none at all. It should be noted that even if all messages are
!! turned off, the model will still halt if it encounters a fatal
!! error!
!! 
!! The other point to note is that when calling the messaging routines,
!! the numerical identifier of a message level should be replaced by the
!! appropriate parameter:
!! - GM_DIAGNOSTIC
!! - GM_TIMESTEP
!! - GM_INFO
!! - GM_WARNING
!! - GM_ERROR
!! - GM_FATAL
module glimmer_log

  use glimmer_global, only : fname_length,dirsep

  implicit none

  integer,parameter :: GM_DIAGNOSTIC = 1 !< Numerical identifier for diagnostic messages.
  integer,parameter :: GM_TIMESTEP   = 2 !< Numerical identifier for timestep messages.
  integer,parameter :: GM_INFO       = 3 !< Numerical identifier for information messages.
  integer,parameter :: GM_WARNING    = 4 !< Numerical identifier for warning messages.
  integer,parameter :: GM_ERROR      = 5 !< Numerical identifier for (non-fatal) error messages.
  integer,parameter :: GM_FATAL      = 6 !< Numerical identifier for fatal error messages.

  integer, parameter            :: GM_levels = 6 !< the number of logging levels
  logical, private, dimension(GM_levels) :: gm_show = .false.

  character(len=*), parameter, dimension(0:GM_levels), private :: msg_prefix = (/ &
       '* UNKNOWN      ', &
       '*              ', &
       '*              ', &
       '               ', &
       '* WARNING:     ', &
       '* ERROR:       ', &
       '* FATAL ERROR :' /) !< array containing log level names


  character(len=fname_length),private :: glimmer_logname !< name of log file
  integer,private :: glimmer_unit = 6                    !< log unit

contains
  !> derives name of log file from file name by stripping directories and appending .log
  function logname(fname)
    implicit none
    character(len=*), intent(in) :: fname !< the file name
    character(len=fname_length) :: logname
    
    character(len=*), parameter :: suffix='.log'
    integer i
    i = scan(fname,dirsep,.True.)
    if (i /= 0) then
       logname = trim(fname(i+1:))//suffix
    else
       logname = trim(fname)//suffix
    end if
  end function logname
    
  !> opens log file
  subroutine open_log(unit,fname)
    use parallel
    implicit none
    integer, optional          :: unit   !< file unit to use
    character(len=*), optional :: fname  !< name of log file

    ! local variables
    character(len=8) :: date
    character(len=10) :: time

    if (present(unit)) then
       glimmer_unit = unit
    end if
    if (present(fname)) then
       glimmer_logname = adjustl(trim(fname))
    else
       glimmer_logname = 'glide.log'
    end if

    if ((main_task).and.(glimmer_unit /= 6)) then
       open(unit=glimmer_unit,file=glimmer_logname,status='unknown')
    end if

    call date_and_time(date,time)
    call write_log_div
    if (main_task) write(unit=glimmer_unit,fmt="(a,a4,'-',a2,'-',a2,' ',a2,':',a2,':',a6)") ' Started logging at ',&
         date(1:4),date(5:6),date(7:8),time(1:2),time(3:4),time(5:10)
    call write_log_div
  end subroutine open_log

  !> write to log
  subroutine write_log(message,type,file,line)
    use glimmer_global, only : msg_length
	use parallel
    implicit none
    integer,intent(in),optional          :: type    !< Type of error to be generated (see list above).
    character(len=*),intent(in)          :: message !< message to be written
    character(len=*),intent(in),optional :: file    !< the name of the file which triggered the message
    integer,intent(in),optional          :: line    !< the line number at the which the message was triggered

    ! local variables
    character(len=msg_length) :: msg
    integer :: local_type
    character(len=6) :: line_num

    local_type = 0
    if (present(type)) then
       if (type >= 1 .or. type <= GM_levels) then
          local_type = type
       end if
    else
       local_type = GM_INFO
    end if

    ! constructing message
    if (present(file) .and. present(line)) then
       if (main_task) write(*,*)"Logged at",file,line
       write(line_num,'(I6)')line
       write(msg,*) trim(msg_prefix(local_type))//' (',trim(file),':',trim(adjustl(line_num)),') '//trim(message)
    else
       write(msg,*) trim(msg_prefix(local_type))//' '//trim(message)
    end if

    ! messages are always written to file log
    if (main_task) write(glimmer_unit,*) trim(msg)

    ! and maybe to std out
    if (local_type /= 0) then
       if ((main_task).and.(gm_show(local_type))) write(*,*) trim(msg)
    end if

    ! stop logging if we encountered a fatal error
    if (local_type == GM_FATAL) then
       if (main_task) write(*,*) "Fatal error encountered, exiting..."
       call close_log
       stop
    end if
  end subroutine write_log

  !> start a new section
  subroutine write_log_div
    use parallel
    implicit none
    if (main_task) write(glimmer_unit,*) '*******************************************************************************'
  end subroutine write_log_div

  !> close log file
  subroutine close_log
    use parallel
    implicit none
    ! local variables
    character(len=8) :: date
    character(len=10) :: time

    call date_and_time(date,time)
    call write_log_div
    if (main_task) write(unit=glimmer_unit,fmt="(a,a4,'-',a2,'-',a2,' ',a2,':',a2,':',a6)") ' Finished logging at ',&
         date(1:4),date(5:6),date(7:8),time(1:2),time(3:4),time(5:10)
    call write_log_div
    
    if (main_task) close(glimmer_unit)
  end subroutine close_log

  !> synchronise log to disk
  subroutine sync_log
    implicit none
    close(glimmer_unit)
    open(unit=glimmer_unit,file=glimmer_logname, position="append", status='old')
  end subroutine sync_log

  !> Sets the output message level.
  subroutine glimmer_set_msg_level(level)
    integer, intent(in) :: level !< The message level (6 is all messages; 0 is no messages). 
    integer :: i

    do i=1,GM_levels
       if (i>(GM_levels-level)) then
          gm_show(i)=.true.
       else
          gm_show(i)=.false.
       endif
    enddo

  end subroutine glimmer_set_msg_level

  !> return glimmer log unit
  function glimmer_get_logunit()
    implicit none
    integer glimmer_get_logunit

    glimmer_get_logunit = glimmer_unit
  end function glimmer_get_logunit
  
  subroutine set_glimmer_unit(unit)

  ! This subroutine should be called when the log file is already open, but glimmer_unit
  ! needs to be set to a desired value (e.g. for CESM coupled runs).
	use parallel
    implicit none
    integer, optional          :: unit   !*FD file unit to use

    ! local variables
    character(len=8) :: date
    character(len=10) :: time

    if (present(unit)) then
       glimmer_unit = unit
    end if

    call date_and_time(date,time)
    call write_log_div
    if (main_task) write(unit=glimmer_unit,fmt="(a,a4,'-',a2,'-',a2,' ',a2,':',a2,':',a6)") &
         ' Started logging at ',&
         date(1:4),date(5:6),date(7:8),time(1:2),time(3:4),time(5:10)
    call write_log_div
  end subroutine set_glimmer_unit

end module glimmer_log
