!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glimmer_global.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2014
!   CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of CISM.
!
!   CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glimmer_global

  !> Module holding global variables for Glimmer. Holds real-type
  !> kind values, and other global code parameters.

  implicit none

  integer,parameter :: sp = kind(1.0) 

  !> Single precision --- Fortran single-precision real-type kind 
  !> value. Used internally.
  !>
  !> Note that if the code is being compiled with forced typing (e.g. with 
  !> the -r8 flag), then this parameter may need to be set in agreement with 
  !> that.

  integer,parameter :: dp = kind(1.0d0) 
  
  !> Double precision --- Fortran double-precision real-type kind 
  !> value. Used internally.
  !>
  !> Note that if the code is being compiled with forced typing (e.g. with
  !> the -r8 flag), then this parameter may need to be set in agreement
  !> with that

!WHL - Removed rk from the code, so commenting out these declarations
!!#ifdef GLIMMER_SP
!!  integer,parameter :: rk=sp !< Precision of glimmer module --- the general Fortran real-type kind value for the Glimmer module and its interfaces.
!!#else
!!  integer,parameter :: rk=dp !< Precision of glimmer module --- the general Fortran real-type kind value for the Glimmer module and its interfaces.
!!#endif

  integer,parameter :: size_t = kind(1)

  !> Precision of glimmer module --- the general Fortran real-type kind value 
  !> for the Glimmer module and its interfaces.
  !>
  !> Note that if the code is being compiled with forced typing (e.g. with 
  !> the -r8 flag), then this parameter must be set in agreement with that. 

  integer,parameter :: fname_length=200 !< Specifies the length of character string variables used to hold filenames.
  integer,parameter :: msg_length=500  !< lenght of message buffers

  !> Specifies the length of character string variables used to
  !> hold filenames.

  character, parameter :: dirsep = '/'
  !> directory separator

  character, parameter :: linefeed = achar(10)          !< ASCII linefeed
  character, parameter :: char_ret = achar(13)          !< ASCII carriage-return
  character(2), parameter :: cr_lf = char_ret//linefeed !< default newline appropriate for UNIX-type systems
  character, parameter :: endline = linefeed
  !> ASCII linefeed and carriage-return characters,
  !> and set up default newline appropriate for UNIX-type systems

  real(kind=dp) :: wall_start_time, wall_stop_time

end module glimmer_global
