! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!***********************************************************************
!
!  mpas_constants
!
!> \brief   MPAS Constant Module
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This module provides various constants that can be used in different parts of MPAS.
!> They may or may not be a physical quantity.
!
!-----------------------------------------------------------------------

module mpas_constants

   use mpas_kind_types

   real (kind=RKIND), parameter :: pii     = 3.141592653589793   !< Constant: Pi
   real (kind=RKIND), parameter :: a       = 6371229.0           !< Constant: Spherical Earth radius [m]
   real (kind=RKIND), parameter :: omega   = 7.29212e-5          !< Constant: Angular rotation rate of the Earth [s-1]
   real (kind=RKIND), parameter :: gravity = 9.80616             !< Constant: Acceleration due to gravity [m s-2]
   real (kind=RKIND), parameter :: rgas    = 287.0               !< Constant: Gas constant for dry air [J kg-1 K-1]
   real (kind=RKIND), parameter :: rv      = 461.6               !< Constant: Gas constant for water vapor [J kg-1 K-1]
   real (kind=RKIND), parameter :: rvord   = rv/rgas             !
   real (kind=RKIND), parameter :: cp      = 1003.0              !< Constant: Specific heat of dry air at constant pressure [J kg-1 K-1]
   real (kind=RKIND), parameter :: cv      = cp - rgas           !< Constant: Specific heat of dry air at constant volume [J kg-1 K-1]
   real (kind=RKIND), parameter :: cvpm    = -cv/cp              !
   real (kind=RKIND), parameter :: prandtl = 1.0                 !< Constant: Prandtl number


   contains


!***********************************************************************
!
!  routine dummy
!
!> \brief   MPAS Dummy Routine
!> \author  Michael Duda
!> \date    03/27/13
!> \details 
!> This is a dummy routine that doesn't do anything.
!
!-----------------------------------------------------------------------
   subroutine dummy()

   end subroutine dummy

end module mpas_constants
