!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glimmer_horiz_bcs_serial.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

!This is a skeleton dummy module for the serial case. It allows serial compilation.

module glimmer_horiz_bcs
  implicit none !If a variable isn't defined, throw a compiler error
  private       !Must declare all public routines in the header.

  !Enforce boundary conditions for a variety of variable types
  public :: horiz_bcs_unstag_scalar   !Unstaggered scalar variables
  public :: horiz_bcs_stag_vector_ew  !Staggered ew-direction vector
  public :: horiz_bcs_stag_vector_ns  !Staggered ns-direction vector
  public :: horiz_bcs_stag_scalar   !Unstaggered scalar variables

  interface horiz_bcs_unstag_scalar
    module procedure horiz_bcs_unstag_scalar_logical_2d
    module procedure horiz_bcs_unstag_scalar_integer_2d
    module procedure horiz_bcs_unstag_scalar_real4_2d
    module procedure horiz_bcs_unstag_scalar_real8_2d
    module procedure horiz_bcs_unstag_scalar_real8_3d
  end interface

  interface horiz_bcs_stag_vector_ew
    module procedure horiz_bcs_stag_vector_ew_real8_2d
    module procedure horiz_bcs_stag_vector_ew_real8_3d
  end interface

  interface horiz_bcs_stag_vector_ns
    module procedure horiz_bcs_stag_vector_ns_real8_2d
    module procedure horiz_bcs_stag_vector_ns_real8_3d
  end interface

  interface horiz_bcs_stag_scalar
    module procedure horiz_bcs_stag_scalar_integer_2d
    module procedure horiz_bcs_stag_scalar_real8_2d
    module procedure horiz_bcs_stag_scalar_real8_3d
  end interface

  integer, parameter, public :: ghost_shift = 0

  logical :: not_printed = .true.

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_unstag_scalar_real8_2d( a )
    implicit none
    real(8),dimension(:,:), intent(inout) :: a
    if (not_printed) write(*,*) 'Using a dummy file for serial building. Horizontal BCs are not being enforced!'
    if (not_printed) write(*,'(A,I5)') 'Printed from line: ', __LINE__
    not_printed = .false.
  end subroutine horiz_bcs_unstag_scalar_real8_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_unstag_scalar_integer_2d( a )
    implicit none
    integer,dimension(:,:), intent(inout) :: a
    if (not_printed) write(*,*) 'Using a dummy file for serial building. Horizontal BCs are not being enforced!'
    if (not_printed) write(*,'(A,I5)') 'Printed from line: ', __LINE__
    not_printed = .false.
  end subroutine horiz_bcs_unstag_scalar_integer_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_unstag_scalar_logical_2d( a )
    implicit none
    logical,dimension(:,:), intent(inout) :: a
    if (not_printed) write(*,*) 'Using a dummy file for serial building. Horizontal BCs are not being enforced!'
    if (not_printed) write(*,'(A,I5)') 'Printed from line: ', __LINE__
    not_printed = .false.
  end subroutine horiz_bcs_unstag_scalar_logical_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_unstag_scalar_real8_3d( a )
    implicit none
    real(8),dimension(:,:,:), intent(inout) :: a
    if (not_printed) write(*,*) 'Using a dummy file for serial building. Horizontal BCs are not being enforced!'
    if (not_printed) write(*,'(A,I5)') 'Printed from line: ', __LINE__
    not_printed = .false.
  end subroutine horiz_bcs_unstag_scalar_real8_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_unstag_scalar_real4_2d( a )
    implicit none
    real(4),dimension(:,:), intent(inout) :: a
    if (not_printed) write(*,*) 'Using a dummy file for serial building. Horizontal BCs are not being enforced!'
    if (not_printed) write(*,'(A,I5)') 'Printed from line: ', __LINE__
    not_printed = .false.
  end subroutine horiz_bcs_unstag_scalar_real4_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_stag_vector_ew_real8_3d( a )
    implicit none
    real(8),dimension(:,:,:), intent(inout) :: a
    if (not_printed) write(*,*) 'Using a dummy file for serial building. Horizontal BCs are not being enforced!'
    if (not_printed) write(*,'(A,I5)') 'Printed from line: ', __LINE__
    not_printed = .false.
  end subroutine horiz_bcs_stag_vector_ew_real8_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_stag_vector_ew_real8_2d( a )
    implicit none
    real(8),dimension(:,:), intent(inout) :: a
    if (not_printed) write(*,*) 'Using a dummy file for serial building. Horizontal BCs are not being enforced!'
    if (not_printed) write(*,'(A,I5)') 'Printed from line: ', __LINE__
    not_printed = .false.
  end subroutine horiz_bcs_stag_vector_ew_real8_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_stag_vector_ns_real8_3d( a )
    implicit none
    real(8),dimension(:,:,:), intent(inout) :: a
    if (not_printed) write(*,*) 'Using a dummy file for serial building. Horizontal BCs are not being enforced!'
    if (not_printed) write(*,'(A,I5)') 'Printed from line: ', __LINE__
    not_printed = .false.
  end subroutine horiz_bcs_stag_vector_ns_real8_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_stag_vector_ns_real8_2d( a )
    implicit none
    real(8),dimension(:,:), intent(inout) :: a
    if (not_printed) write(*,*) 'Using a dummy file for serial building. Horizontal BCs are not being enforced!'
    if (not_printed) write(*,'(A,I5)') 'Printed from line: ', __LINE__
    not_printed = .false.
  end subroutine horiz_bcs_stag_vector_ns_real8_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_stag_scalar_integer_2d( a )
    implicit none
    integer,dimension(:,:), intent(inout) :: a
    if (not_printed) write(*,*) 'Using a dummy file for serial building. Horizontal BCs are not being enforced!'
    if (not_printed) write(*,'(A,I5)') 'Printed from line: ', __LINE__
    not_printed = .false.
  end subroutine horiz_bcs_stag_scalar_integer_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_stag_scalar_real8_2d( a )
    implicit none
    real(8),dimension(:,:), intent(inout) :: a
    if (not_printed) write(*,*) 'Using a dummy file for serial building. Horizontal BCs are not being enforced!'
    if (not_printed) write(*,'(A,I5)') 'Printed from line: ', __LINE__
    not_printed = .false.
  end subroutine horiz_bcs_stag_scalar_real8_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_stag_scalar_real8_3d( a )
    implicit none
    real(8),dimension(:,:,:), intent(inout) :: a
    if (not_printed) write(*,*) 'Using a dummy file for serial building. Horizontal BCs are not being enforced!'
    if (not_printed) write(*,'(A,I5)') 'Printed from line: ', __LINE__
    not_printed = .false.
  end subroutine horiz_bcs_stag_scalar_real8_3d

end module glimmer_horiz_bcs

