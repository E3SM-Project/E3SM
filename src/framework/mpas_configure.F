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

#include "namelist_call.inc"
#include "namelist_defines.inc"

end module mpas_configure
