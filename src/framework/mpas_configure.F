! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!-----------------------------------------------------------------------
!  mpas_configure
!
!> \brief MPAS Configuration routines.
!> \author Michael Duda, Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This module will contain all namelist parameter definitions, as well as the routine which reads them from the namelist file.
!
!-----------------------------------------------------------------------
module mpas_configure

   use mpas_dmpar
   use mpas_io_units

#include "config_defs.inc"

   contains

!-----------------------------------------------------------------------
!  routine mpas_read_namelist
!
!> \brief MPAS read namelist routine
!> \author Michael Duda, Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine reads and broadcasts the namelist file. 
!
!-----------------------------------------------------------------------
   subroutine mpas_read_namelist(dminfo, nml_filename)

      implicit none

      type (dm_info), intent(in) :: dminfo !< Input: Domain information
      character (len=*), optional :: nml_filename !< Input - Optional: Namelist filename. Defaults to namelist.input

      integer :: funit, ierr

#include "config_namelist_defs.inc"

      funit = 21

      ! Set default values for namelist options
#include "config_set_defaults.inc"

      if (dminfo % my_proc_id == IO_NODE) then
         if (present(nml_filename)) then
            write(stderrUnit,*) 'Reading ', trim(nml_filename)
            open(funit,file=trim(nml_filename),status='old',form='formatted')
         else
            write(stderrUnit,*) 'Reading namelist.input'
            open(funit,file='namelist.input',status='old',form='formatted')
         end if

#include "config_namelist_reads.inc"
         close(funit)
         write(stderrUnit,*) ' '
      end if

#include "config_bcast_namelist.inc"

   end subroutine mpas_read_namelist

end module mpas_configure
