! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!***********************************************************************
!
!  mpas_kind_types
!
!> \brief   MPAS Kind definition module
!> \author  Michael Duda, Doug Jacobsen
!> \date    03/27/13
!> \details 
!> This module defines the kind types for basic fortran data types within MPAS.
!
!-----------------------------------------------------------------------

module mpas_kind_types

   integer, parameter :: R4KIND = selected_real_kind(6)
#ifdef SINGLE_PRECISION
   integer, parameter :: RKIND  = selected_real_kind(6)
#else
   integer, parameter :: RKIND  = selected_real_kind(12)
#endif

   integer, parameter :: I8KIND = selected_int_kind(18)

   integer, parameter :: StrKIND = 512
   integer, parameter :: ShortStrKIND = 32

   contains

!***********************************************************************
!
!  routine dummy_kinds
!
!> \brief   MPAS Dummy kind routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This is a dummy routine that doesn't do anything.
!
!-----------------------------------------------------------------------
   subroutine dummy_kinds()

   end subroutine dummy_kinds

end module mpas_kind_types
