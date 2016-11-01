! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 program build_tables
 use atmphys_build_tables_thompson

 implicit none

 write(0,*) '--- Building tables for Thompson cloud microphysics scheme:'
 call build_tables_thompson

 stop

 end program build_tables
!=================================================================================================================
