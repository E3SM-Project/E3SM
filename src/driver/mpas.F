! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
program mpas

   use mpas_subdriver

   implicit none

   call mpas_init()

   call mpas_run() 

   call mpas_finalize()

   stop

end program mpas
