! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
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

   integer, save :: stdoutUnit, stderrUnit
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

        do i = 1, maxUnits
            if(.not. unitsInUse(i) ) then
                newUnit = i
                unitsInUse(newUnit) = .true.
                return
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

!***********************************************************************
!
!  routine mpas_io_units_set_stdout
!
!> \brief   MPAS Set stdout unit number
!> \author  Doug Jacobsen
!> \date    07/18/2013
!> \details 
!> This routine sets the unit number for stdout
!
!-----------------------------------------------------------------------
    subroutine mpas_io_units_set_stdout(stdoutUnit_in)!{{{
        integer, intent(in) :: stdoutUnit_in


        if(stdoutUnit /= stdoutUnit_in) then
            unitsInUse(stdoutUnit) = .false.
            stdoutUnit = stdoutUnit_in
            unitsInUse(stdoutUnit) = .true.
        end if

    end subroutine mpas_io_units_set_stdout!}}}

!***********************************************************************
!
!  routine mpas_io_units_set_stderr
!
!> \brief   MPAS Set stderr unit number
!> \author  Doug Jacobsen
!> \date    07/18/2013
!> \details 
!> This routine sets the unit number for stderr
!
!-----------------------------------------------------------------------
    subroutine mpas_io_units_set_stderr(stderrUnit_in)!{{{
        integer, intent(in) :: stderrUnit_in


        if(stderrUnit /= stderrUnit_in) then
            unitsInUse(stderrUnit) = .false.
            stderrUnit = stderrUnit_in
            unitsInUse(stderrUnit) = .true.
        end if

    end subroutine mpas_io_units_set_stderr!}}}

!***********************************************************************
!
!  routine mpas_init_logger
!
!> \brief   MPAS Logger initialization
!> \author  Doug Jacobsen
!> \date    07/18/2013
!> \details 
!> This routine initializes the MPAS logger
!
!-----------------------------------------------------------------------
    subroutine mpas_io_units_init(stdoutUnit_in, stderrUnit_in)!{{{
        integer, intent(in), optional :: stdoutUnit_in !< Input - Optional: Unit to set stdoutUnit to
        integer, intent(in), optional :: stderrUnit_in !< Input - Optional: Unit to set stderrUnit to

        if(present(stdoutUnit_in)) then
            stdoutUnit = stdoutUnit_in
        else
            stdoutUnit = 6
        end if

        if(present(stderrUnit_in)) then
            stderrUnit = stderrUnit_in
        else
            stderrUnit = 0
        end if
        
        unitsInUse(stdoutUnit) = .true.
        unitsInUse(stderrUnit) = .true.

    end subroutine mpas_io_units_init!}}}

end module mpas_io_units
