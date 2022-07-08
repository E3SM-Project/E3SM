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

   implicit none

   private

   integer, parameter :: maxUnits = 200
   logical, dimension(0:maxUnits), save :: unitsInUse

   ! Units reserved for unformatted I/O
   integer, parameter :: unformatted_min = 101
   integer, parameter :: unformatted_max = maxUnits

   public :: mpas_new_unit, &
             mpas_release_unit

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
    subroutine mpas_new_unit(newUnit, unformatted)!{{{

        integer, intent(inout) :: newUnit
        logical, optional, intent(in) :: unformatted

        integer :: i, minsearch, maxsearch

        logical :: opened

        newUnit = -1

        !
        ! Determine the range over which to search for an unused unit
        !
        minsearch = 1
        maxsearch = unformatted_min - 1
        if ( present(unformatted) ) then
           if ( unformatted ) then
              minsearch = unformatted_min
              maxsearch = unformatted_max
           end if
        end if

        do i = minsearch, maxsearch
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

        if (0 <= releasedUnit .and. releasedUnit <= maxUnits) then
           unitsInUse(releasedUnit) = .false.
        end if

    end subroutine mpas_release_unit!}}}

end module mpas_io_units
