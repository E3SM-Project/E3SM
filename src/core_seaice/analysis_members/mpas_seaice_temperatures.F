! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_temperatures
!
!> \brief MPAS-Seaice analysis mode member: temperatures
!> \author Adrian K. Turner
!> \date   6th September 2015
!> \details
!>  MPAS-Seaice analysis mode member: temperatures
!>  Calculates ice and snow temperatures using the CICE column
!>  package routines.
!>
!-----------------------------------------------------------------------

module seaice_temperatures

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

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

   public :: seaice_bootstrap_temperatures, &
             seaice_init_temperatures, &
             seaice_precompute_temperatures, &
             seaice_compute_temperatures, &
             seaice_restart_temperatures, &
             seaice_finalize_temperatures

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine seaice_bootstrap_temperatures
!
!> \brief   Bootstrap MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    10th November 2015
!> \details
!>  This routine conducts all bootstraps required for the
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_bootstrap_temperatures(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine seaice_bootstrap_temperatures!}}}

!***********************************************************************
!
!  routine seaice_init_temperatures
!
!> \brief   Initialize MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    6th September 2015
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_init_temperatures(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine seaice_init_temperatures!}}}

!***********************************************************************
!
!  routine seaice_precompute_temperatures
!
!> \brief   Precompute MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    9th September 2015
!> \details
!>
!
!-----------------------------------------------------------------------

   subroutine seaice_precompute_temperatures(domain, instance, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine seaice_precompute_temperatures!}}}

!***********************************************************************
!
!  routine seaice_compute_temperatures
!
!> \brief   Compute MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    6th September 2015
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Seaice analysis member. It computes ice and snow layer
!>  temperatures from enthalpies and salinities using CICE column
!>  package routines.
!
!-----------------------------------------------------------------------

   subroutine seaice_compute_temperatures(domain, instance, timeLevel, err)!{{{

      use ice_colpkg, only: &
           colpkg_ice_temperature, &
           colpkg_snow_temperature

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      type (block_type), pointer :: &
           block

      type (mpas_pool_type), pointer :: &
           tracersPool, &
           temperaturesAMPool

      real(kind=RKIND), dimension(:,:,:), pointer :: &
           iceAreaCategory, &
           snowVolumeCategory, &
           iceTemperature, &
           snowTemperature, &
           iceEnthalpy, &
           snowEnthalpy, &
           iceSalinity

      integer, pointer :: &
           nCellsSolve, &
           nCategories, &
           nIceLayers, &
           nSnowLayers

      integer :: &
           iCell, &
           iCategory, &
           iIceLayer, &
           iSnowLayer

      err = 0

      block => domain % blocklist
      do while (associated(block))

         call MPAS_pool_get_subpool(block % structs, 'temperaturesAM', temperaturesAMPool)
         call MPAS_pool_get_subpool(block % structs, 'tracers', tracersPool)

         call MPAS_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call MPAS_pool_get_dimension(block % dimensions, 'nCategories', nCategories)
         call MPAS_pool_get_dimension(block % dimensions, 'nIceLayers', nIceLayers)
         call MPAS_pool_get_dimension(block % dimensions, 'nSnowLayers', nSnowLayers)

         call MPAS_pool_get_array(temperaturesAMPool, 'iceTemperature', iceTemperature)
         call MPAS_pool_get_array(temperaturesAMPool, 'snowTemperature', snowTemperature)

         call MPAS_pool_get_array(tracersPool, 'iceAreaCategory', iceAreaCategory, 1)
         call MPAS_pool_get_array(tracersPool, 'snowVolumeCategory', snowVolumeCategory, 1)
         call MPAS_pool_get_array(tracersPool, 'iceEnthalpy', iceEnthalpy, 1)
         call MPAS_pool_get_array(tracersPool, 'snowEnthalpy', snowEnthalpy, 1)
         call MPAS_pool_get_array(tracersPool, 'iceSalinity', iceSalinity, 1)

         ! initialize
         iceTemperature = 0.0_RKIND
         snowTemperature = 0.0_RKIND

         ! compute temperatures
         do iCell = 1, nCellsSolve
            do iCategory = 1, nCategories

               ! check if ice present
               if (iceAreaCategory(1,iCategory,iCell) > 1e-11_RKIND) then

                  ! ice layers
                  do iIceLayer = 1, nIceLayers
                     iceTemperature(iIceLayer, iCategory, iCell) = &
                          colpkg_ice_temperature(iceEnthalpy(iIceLayer,iCategory,iCell), &
                                                 iceSalinity(iIceLayer,iCategory,iCell))
                  enddo ! iIceLayer

                  ! snow layers
                  if (snowVolumeCategory(1,iCategory,iCell) > 1e-11_RKIND) then

                     do iSnowLayer = 1, nSnowLayers
                        snowTemperature(iSnowLayer, iCategory, iCell) = &
                             colpkg_snow_temperature(snowEnthalpy(iSnowLayer,iCategory,iCell))
                     enddo ! iIceLayer

                  endif ! snowVolumeCategory

               endif ! iceAreaCategory

            enddo ! iCategory
         enddo ! iCell

         block => block % next
      enddo

   end subroutine seaice_compute_temperatures!}}}

!***********************************************************************
!
!  routine seaice_restart_temperatures
!
!> \brief   Save restart for MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    6th September 2015
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_restart_temperatures(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine seaice_restart_temperatures!}}}

!***********************************************************************
!
!  routine seaice_finalize_temperatures
!
!> \brief   Finalize MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    6th September 2015
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_finalize_temperatures(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine seaice_finalize_temperatures!}}}

!-----------------------------------------------------------------------

end module seaice_temperatures

! vim: foldmethod=marker
