! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
program mpas

   use mpas_subdriver
   use mpas_derived_types, only : core_type, domain_type

   implicit none

   type (core_type), pointer :: corelist => null()
   type (domain_type), pointer :: domain => null()

   call mpas_init(corelist, domain)

   call mpas_run(domain) 

   call mpas_finalize(corelist, domain)

   stop

end program mpas
