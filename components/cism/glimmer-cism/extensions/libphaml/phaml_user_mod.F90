!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   phaml_user_mod.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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


!------------------------------------------------------------------------------
! This module contains user global data, and functions for formatting
!------------------------------------------------------------------------------
module phaml_user_mod
    use global
    use phaml
    implicit none

    integer, save :: gnsn, gewn, gdns, gdew, num_arrays, modnum
    real(my_real), save, allocatable, dimension(:,:) :: uphaml
contains        
    !------------------------------------------------------------------------------
    !SUBROUTINE: user_init
    !ARGUMENTS: ewn,nsn,dew,dns,modn
    !DESCRIPTION:
    ! This subroutine sets all of the mandatory fields in the user module.  The 
    ! num_arrays is set to 1 by default but can be changed if more than one array
    ! needs to be stored in the usermod.  This will require altering the array_init,
    ! as well as the update_usermod function specific to the module.
    !------------------------------------------------------------------------------       
    subroutine user_init(ewn,nsn,dew,dns,modn)
        implicit none
        integer, intent(in) :: nsn,ewn,modn
        real(my_real), intent(in) :: dns,dew
        gnsn = nsn  !ny
        gewn = ewn  !nx
        gdns = dns  !dy
        gdew = dew  !dx
        modnum = modn
        !this defines how many arrays are passed in update_usermod
        num_arrays = 1
        !allocate(uphaml(gewn, gnsn))
        !allocate(uphaml_one(gewn*gnsn))
    end subroutine user_init
    !------------------------------------------------------------------------------
    !SUBROUTINE: user_close
    !ARGUMENTS: 
    !DESCRIPTION:
    ! This simply calls any functions that might be needed to free memory or end 
    ! the use of the user module.  Currently, only the arrays need to be deallocated.
    !------------------------------------------------------------------------------
    subroutine user_close()
        call array_close()
    end subroutine user_close
    !------------------------------------------------------------------------------
    !SUBROUTINE: array_init
    !ARGUMENTS: 
    !DESCRIPTION:
    ! A subroutine to allocate the arrays.
    !------------------------------------------------------------------------------
    subroutine array_init()
        !if (size(uphaml) .le. 0) then
            allocate(uphaml(gewn, gnsn))
        !end if
    end subroutine array_init
    !------------------------------------------------------------------------------
    !SUBROUTINE: array_close
    !ARGUMENTS: 
    !DESCRIPTION:
    ! This simply deallocates any dynamic arrays within the usermodule.
    !------------------------------------------------------------------------------       
    subroutine array_close()
        if (size(uphaml) > 0) then
            deallocate(uphaml)
        end if
    end subroutine array_close
    !------------------------------------------------------------------------------
    !SUBROUTINE: combine_arrays
    !ARGUMENTS: array1, array2, retarray
    !DESCRIPTION:
    ! Given two arrays(array1, array2) this function simply concatenates the two 
    ! together as one array.
    !------------------------------------------------------------------------------       
    subroutine combine_arrays(array1, array2, retarray)
        implicit none
        real(my_real), intent(in),dimension(:) :: array1,array2
        real(my_real), intent(out),dimension(:) :: retarray
        integer :: i,maxind
        
        maxind = size(array1) + size(array2)
        do i=1, maxind
            retarray(i) = array1(i)
            retarray(maxind + i) = array2(i)
        end do
    end subroutine combine_arrays
    !------------------------------------------------------------------------------
    !SUBROUTINE: split_arrays
    !ARGUMENTS: array1, array2, retarray, maxind
    !DESCRIPTION:
    ! Given a single array(retarray) this functions splits it based off of the maxind
    ! variable and puts the data back into two individual arrays.
    !------------------------------------------------------------------------------       
    subroutine split_arrays(array1, array2, retarray, maxind)
        implicit none
        real(my_real), intent(out),dimension(:) :: array1,array2
        real(my_real), intent(in),dimension(:) :: retarray
        integer, intent(in) :: maxind
        integer :: i
        
        do i=1, maxind
            array1(i) = retarray(i)
            array2(i) = retarray(maxind + i)
        end do
    end subroutine split_arrays
    !------------------------------------------------------------------------------
    !SUBROUTINE: reshape_array_to_one
    !ARGUMENTS: twod, oned
    !DESCRIPTION:
    ! This function takes a two dimensional vector and reshapes it to a single one 
    ! dimensional vector based off of the CISM grid properties. 
    !------------------------------------------------------------------------------
    subroutine reshape_array_to_one(twod,oned)
        implicit none
        real(my_real), intent(in),dimension(:,:) :: twod
        real(my_real), intent(out),dimension(:) :: oned
        integer :: r,c
        
         do r=0,gnsn-1
            do c=0,gewn-1
               oned(gewn*r+c+1) = twod(c+1,r+1)
            end do
        end do
    end subroutine reshape_array_to_one
    !------------------------------------------------------------------------------
    !SUBROUTINE: reshape_array_to_two
    !ARGUMENTS: twod, oned
    !DESCRIPTION:
    ! This function takes a one dimensional vector and reshapes it to a two 
    ! dimensional vector based off of the CISM grid properties. 
    !------------------------------------------------------------------------------  
    subroutine reshape_array_to_two(twod,oned)
        implicit none
        real(my_real), intent(out),dimension(:,:) :: twod
        real(my_real), intent(in),dimension(:) :: oned
        integer :: count,r,c
        
        do count=0,size(oned) - 1
            r = count/gewn + 1
            c = int(Mod(count,gewn)) + 1
            twod(c,r) = real(oned(count+1))        
        end do
    end subroutine reshape_array_to_two
    
    !------------------------------------------------------------------------------
    !SUBROUTINE: get_xyarrays
    !ARGUMENTS: x, y 
    !DESCRIPTION:
    ! This subroutine returns x and y vectors of all the nodes in the CISM grid.  
    ! This is needed in order to get the solution from phaml.
    !------------------------------------------------------------------------------
    subroutine get_xyarrays(x,y)
        implicit none
        real(my_real), intent(inout) :: x(:),y(:)
        integer :: count,r,c

        count = 1
        !populate the arrays with the xy coordinates
        do r=0,gnsn-1
            do c=0,gewn-1
               x(count) = REAL(c*gdew)
               y(count) = REAL(r*gdns)
               count = count + 1
            end do
        end do
    end subroutine get_xyarrays
    
    
    
    !------------------------------------------------------------------------------
    !SUBROUTINE: getew
    !ARGUMENTS: x
    !DESCRIPTION:
    ! This subroutine simply truncates an x by the east/west distance to the nearest
    ! corresponding grid point in CISM.  NOTE: The points in the mesh are not
    ! quadrature and you can not use this function there without losing accuracy.
    !------------------------------------------------------------------------------
    integer function getew(x)
        implicit none
        real(my_real), intent(in) :: x
        getew = int((x/gdew)) + 1
    end function getew
    !------------------------------------------------------------------------------
    !SUBROUTINE: getns
    !ARGUMENTS: y
    !DESCRIPTION:
    ! This subroutine simply truncates a y by the north/south distance to the nearest
    ! corresponding grid point in CISM.  NOTE: The points in the mesh are not
    ! quadrature and you can not use this function there without losing accuracy.
    !------------------------------------------------------------------------------
    integer function getns(y)
        implicit none
        real(my_real), intent(in) :: y
        getns = int((y/gdns)) + 1
    end function getns
    
end module phaml_user_mod
