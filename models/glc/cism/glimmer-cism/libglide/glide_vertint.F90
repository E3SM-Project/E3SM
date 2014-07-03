!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glide_vertint.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

!TODO - Remove this module?  Not currently used.

module glide_vertint

    !*FD This module contains routines to vertically integrate fields
    !*FD All 3d fields are assumed to use the (z,x,y) coordinate system,
    !*FD where the top is the minimum z and the bottom is the maximum z.

    use glimmer_global , only: dp
    implicit none

contains

    !*FD Performs vertical integration, places the result on a 3d field
    !*FD where each level in the 3d field is the integral of all levels
    !*FD above it

    subroutine vertint_output3d(infield, outfield, levels, topdown, initial_value)
        real(dp), dimension(:,:,:), intent(in) :: infield
        real(dp), dimension(:,:,:), intent(out) :: outfield
        real(dp), dimension(:), intent(in) :: levels
        logical :: topdown !*FD Controls the direction of integration.  If true,
                           !*FD outfield(1,:,:) contains zeros and each level
                           !*FD below it accumulates another part of the
                           !*FD integral.  If false, outfield(upn,:,:) contains
                           !*FD zeros and each level above it accumulates
                           !*FD another part of the integral
        real(dp), dimension(:,:), intent(in), optional :: initial_value

        integer :: upn
        integer :: i
        integer :: lower, upper, step !Loop control, parameterized based on
                                      !value of topdown
        
        real(dp) :: deltax
        
        upn = size(levels)
        if (topdown) then
            lower = 2
            upper = upn
            step  = 1
        else
            lower = upn-1
            upper = 1
            step = -1
        end if

        if (present(initial_value)) then
            outfield(lower - step,:,:) = initial_value
        else
            outfield(lower - step,:,:) = 0
        end if
        

        do i = lower, upper, step
            deltax = step*(levels(i) - levels(i - step))
            !Apply trapezoid rule
            outfield(i,:,:) = outfield(i - step,:,:) + .5 * deltax*(infield(i - step,:,:) + infield(i,:,:))
        end do
    end subroutine vertint_output3d

    subroutine vertint_output2d(infield, outfield, levels, initial_value)
        !*FD Vertically integrates the 3D field and places the result of the
        !*FD integral on a 2D field
        real(dp), dimension(:,:,:), intent(in)  :: infield
        real(dp), dimension(:,:),   intent(out) :: outfield
        real(dp), dimension(:),     intent(in)  :: levels

        real(dp), dimension(:,:), intent(in), optional :: initial_value

        integer :: upn
        integer :: i
        real(dp) :: deltax

        upn = size(levels)
        
        if (present(initial_value)) then
            outfield = initial_value
        else
            outfield = 0
        end if

        do i = 2, upn
            deltax = levels(i) - levels(i - 1)
            outfield = outfield + .5 * deltax*(infield(i-1,:,:) + infield(i,:,:))
        end do
    end subroutine


    !Contained unit test cases
    !Based around evaluation of the integral of x^2dx from 0 to 1.
    subroutine test_vertint()
        real(dp), dimension(11) :: levels
        real(dp), dimension(11,1,1) :: values
        real(dp), dimension(1,1) :: answer
        real(dp), dimension(1,1) :: ival

        integer :: i
        real(dp) :: val


        !Test case where we have evenly spaced levels
        val = 0
        do i = 1,11
            levels(i) = val
            values(i,1,1) = val ** 2
            val = val + .1
            write(*,*) levels(i),values(i,1,1)
        end do

        ival = 0

        call vertint_output2d(values, answer, levels, ival)
        write(*,*) answer(1,1)

        !Test case where we do not have evenly spaced levels
        levels(1) = 0
        levels(2) = .2
        levels(3) = .4
        levels(4) = .5
        levels(5) = .6
        levels(6) = .7
        levels(7) = .8
        levels(8) = .85
        levels(9) = .9
        levels(10) = .95
        levels(11) = 1
        do i = 1,11
            values(i,1,1) = levels(i) ** 2
            write(*,*) levels(i),values(i,1,1)
        end do
        ival = 0

        call vertint_output2d(values, answer, levels, ival)
        write(*,*) answer(1,1)

    end subroutine
end module glide_vertint
