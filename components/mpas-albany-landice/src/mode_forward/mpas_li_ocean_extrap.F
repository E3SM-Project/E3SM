! Copyright (c) 2013-2018,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_ocean_extrap
!
!> \MPAS land-ice ocean-data extrapolation driver
!> \author Holly Han
!> \date   January 2022
!> \details
!>  This module contains the routines for extrapolating
!>  ocean data (e.g., temperature, salinity, thermal forcing)
!>  into ice draft
!
!-----------------------------------------------------------------------

module li_ocean_extrap

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_log
   use li_mask
   use li_setup

   implicit none
   private

   !--------------------------------------------------------------------
   ! Public parameters
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   ! Public member functions
   !--------------------------------------------------------------------

   public :: li_ocean_extrap_solve

   !--------------------------------------------------------------------
   ! Private module variables
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine li_ocean_extrap_solve
!
!> \brief   Initializes ocean extrapolation scheme
!> \author  Holly Han
!> \date    12 Jan 2023
!> \details
!>  This routine performs horizontal and vertical extrapolation
!>  of ocean data (e.g., temperature, salinity, thermal forcing)
!
!-----------------------------------------------------------------------

   subroutine li_ocean_extrap_solve(domain, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      logical, pointer :: config_ocean_data_extrapolation

      ! No init is needed.
      err = 0

      call mpas_pool_get_config(liConfigs, 'config_ocean_data_extrapolation', config_ocean_data_extrapolation)

      if ( config_ocean_data_extrapolation ) then
         ! call the extrapolation scheme
         call mpas_log_write('ocean data will be extrapolated into the MALI ice draft')
         
         ! process temperature field
             ! get geometry pool of data e.g., 'ismip6shelfMelt_3dThermalForcing'
             ! get a valid mask (preprocessed valid ocean data mask)
             ! get an updated seedMask
             ! call horizontal extrapolation
             ! call vertical extrapolation 
                 ! switch between horizontal and vertical extrap ??
         ! process salinity field
             ! get a valid mask (proprocessed valid ocean data mask)
             ! get an updated seedMask 
             ! call horizontal extrapolation
             ! call vertical extrapolation 
             
         ! calculate thermal forcing
      else
         ! do nothing
         call mpas_log_write('ocean data will NOT be extrapolated into the MALI ice draft')
      endif
   !--------------------------------------------------------------------

   end subroutine li_ocean_extrap_solve

end module li_ocean_extrap

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||