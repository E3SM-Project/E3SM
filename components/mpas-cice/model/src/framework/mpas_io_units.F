! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!***********************************************************************
!
!  mpas_io_units
!
!> \brief   MPAS Logging module
!> \author  Doug Jacobsen
!> \date    07/18/2013
!> \details 
!> This module contains routines for handling log and error files.
!
!-----------------------------------------------------------------------

module mpas_io_units

   use mpas_kind_types

   integer, parameter, private :: maxUnits = 99
   logical, dimension(0:maxUnits), private, save :: unitsInUse

   contains

!***********************************************************************
!
!  routine mpas_new_unit
!
!> \brief   MPAS New unit routine
!> \author  Doug Jacobsen
!> \date    07/18/2013
!> \details 
!> This routine determines a new unit that is not yet in use, and returns
!> the unit number
!
!-----------------------------------------------------------------------
    subroutine mpas_new_unit(newUnit)!{{{
        integer, intent(inout) :: newUnit

        integer :: i

        logical :: opened

        do i = 1, maxUnits
            if (.not. unitsInUse(i)) then
                inquire(i, opened=opened)
                if (opened) then
                    unitsInUse(i) = .true.
                else
                    newUnit = i
                    unitsInUse(newUnit) = .true.
                    return
                endif
            end if
        end do

    end subroutine mpas_new_unit!}}}

!***********************************************************************
!
!  routine mpas_release_unit
!
!> \brief   MPAS Release unit routine
!> \author  Doug Jacobsen
!> \date    07/18/2013
!> \details 
!> This routine releases a unit that is in use.
!
!-----------------------------------------------------------------------
    subroutine mpas_release_unit(releasedUnit)!{{{
        integer, intent(in) :: releasedUnit

        unitsInUse(releasedUnit) = .false.

    end subroutine mpas_release_unit!}}}

end module mpas_io_units
