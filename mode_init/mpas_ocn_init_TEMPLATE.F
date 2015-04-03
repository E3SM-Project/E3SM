! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_TEMPLATE
!
!> \brief MPAS ocean initialize case -- TEMPLATE
!> \author Doug Jacobsen
!> \date   03/23/2015
!> \details
!>  This module contains the routines for initializing the 
!>  the TEMPLATE test case
!
!-----------------------------------------------------------------------

module ocn_init_TEMPLATE

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants

   use ocn_constants
   use ocn_init_vertical_grids
   use ocn_init_cell_markers

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_init_setup_TEMPLATE, &
             ocn_init_validate_TEMPLATE

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_TEMPLATE
!
!> \brief   Setup for baroclinic channel test case
!> \author  Doug Jacobsen
!> \date    03/23/2015
!> \details 
!>  This routine sets up the initial conditions for the baroclinic channel test case.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_setup_TEMPLATE(domain, err)!{{{

   !--------------------------------------------------------------------

      type (domain_type), intent(inout) :: domain
      integer, intent(out) :: err

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_configuration', config_configuration)

      if(config_configuration .ne. trim('TEMPLATE')) return

      ! Setup configuration

   !--------------------------------------------------------------------

   end subroutine ocn_init_setup_TEMPLATE!}}}

!***********************************************************************
!
!  routine ocn_init_validate_TEMPLATE
!
!> \brief   Validation for baroclinic channel test case
!> \author  Doug Jacobsen
!> \date    03/23/2015
!> \details 
!>  This routine validates the configuration options for the baroclinic channel test case.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_TEMPLATE(configPool, packagePool, err)!{{{

   !--------------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: configPool, packagePool

      integer, intent(out) :: err

      character (len=StrKIND), pointer :: config_configuration
      integer, pointer :: config_vert_levels, config_TEMPLATE_vert_levels

      err = 0

      call mpas_pool_get_config(configPool, 'config_configuration', config_configuration)

      if(config_configuration .ne. trim('TEMPLATE')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_TEMPLATE_vert_levels', config_TEMPLATE_vert_levels)

      if(config_vert_levels <= 0 .and. config_TEMPLATE_vert_levels > 0) then
         config_vert_levels = config_TEMPLATE_vert_levels
      else if (config_vert_levels <= 0) then
         write(stderrUnit,*) 'ERROR: Validation failed for TEMPLATE. Not given a usable value for vertical levels.'
         err = 1
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_TEMPLATE!}}}


!***********************************************************************

end module ocn_init_TEMPLATE

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
