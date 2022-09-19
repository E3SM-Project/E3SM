! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vertical_advection
!
!> \brief  MPAS ocean vertical advection module
!> \author Carolyn Begeman
!> \date   February 2022
!> \details
!>  This module contains the initialization for vertical advection
!>  schemes.
!
!-----------------------------------------------------------------------

module ocn_vertical_advection

   use ocn_config
   use mpas_log

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------
   integer, public ::            &
      configVertAdvMethod         ! choice of vertical advection method
   integer, parameter, public :: &! supported vertical advection methods
      vertAdvFluxForm = 1,       &! flux form
      vertAdvRemap    = 2         ! remapping
   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------
   public :: ocn_vertical_advection_init

contains

!***********************************************************************

!***********************************************************************
!
!  routine ocn_vertical_advection_init
!
!> \brief   Initializes vertical advection scheme properties
!> \author  Carolyn Begeman
!> \date    February 2022
!> \details
!>  This routine initializes quantities related to
!>  the vertical advection scheme.
!
!----------------------------------------------------------------------

   subroutine ocn_vertical_advection_init(err)!{{{

   !--------------------------------------------------------------------

      integer, intent(out) :: err


      err = 0

      select case (trim(config_vert_advection_method))

      case ('flux-form')

         configVertAdvMethod  = vertAdvFluxForm

      case ('remap')

         configVertAdvMethod  = vertAdvRemap

      case default

         call mpas_log_write( &
            'Invalid choice for config_vert_advection_method. Choices are: flux-form, remap')
         err = 1

      end select

   end subroutine ocn_vertical_advection_init!}}}

!***********************************************************************

end module ocn_vertical_advection

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
