!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_short_wave_absorption
!
!> \brief MPAS ocean tracer short wave
!> \author Doug Jacobsen
!> \date   12/17/12
!> \details
!>  This module contains the routine for computing
!>  short wave tendencies
!
!-----------------------------------------------------------------------

module ocn_tracer_short_wave_absorption

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines
   use ocn_constants
   use ocn_tracer_short_wave_absorption_jerlov
   use ocn_tracer_short_wave_absorption_variable

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

   public :: ocn_tracer_short_wave_absorption_tend, &
             ocn_tracer_short_wave_absorption_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: useJerlov

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_short_wave_absorption_tend
!
!> \brief   Computes tendency term for surface fluxes
!> \author  Doug Jacobsen
!> \date    12/17/12
!> \details
!>  This routine computes the tendency for tracers based on surface fluxes.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_short_wave_absorption_tend(meshPool, swForcingPool, forcingPool, index_temperature, & !{{{
                    layerThickness, penetrativeTemperatureFlux, penetrativeTemperatureFluxOBL, tend, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool, swForcingPool, forcingPool          !< Input: mesh information

      real (kind=RKIND), dimension(:), intent(in) :: &
         penetrativeTemperatureFlux !< Input: short wave heat flux

      real (kind=RKIND), dimension(:), intent(inout) :: &
         penetrativeTemperatureFluxOBL

      real (kind=RKIND), dimension(:,:), intent(in) :: layerThickness !< Input: Layer thicknesses

      integer, intent(in) :: index_temperature

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tend          !< Input/Output: velocity tendency

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

      character(len=strKIND), pointer :: config_sw_absorption_type

      call MPAS_pool_get_config(ocnConfigs, 'config_sw_absorption_type', config_sw_absorption_type)

      if (trim(config_sw_absorption_type)=='none')  return

      call mpas_timer_start("short wave")

      err = 0
      if(useJerlov) then
         call ocn_tracer_short_wave_absorption_jerlov_tend(meshPool, forcingPool, index_temperature, layerThickness, &
                                        penetrativeTemperatureFlux, penetrativeTemperatureFluxOBL, tend, err)
      else
         call ocn_tracer_short_wave_absorption_variable_tend(meshPool,swForcingPool, forcingPool, index_temperature, &
                                        layerThickness, penetrativeTemperatureFlux, penetrativeTemperatureFluxOBL,tend,err)
      endif

      call mpas_timer_stop("short wave")
   !--------------------------------------------------------------------

   end subroutine ocn_tracer_short_wave_absorption_tend!}}}

!***********************************************************************
!
!  routine ocn_tracer_short_wave_absorption_init
!
!> \brief   Initializes ocean tracer surface flux quantities
!> \author  Doug Jacobsen
!> \date    12/17/12
!> \details
!>  This routine initializes quantities related to surface fluxes in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_short_wave_absorption_init(domain,err)!{{{

   !--------------------------------------------------------------------

      type (domain_type) :: domain

      integer, intent(out) :: err !< Output: error flag

      logical, pointer :: config_use_activeTracers_surface_bulk_forcing
      character (len=StrKind), pointer :: config_sw_absorption_type

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_sw_absorption_type', config_sw_absorption_type)
      call MPAS_pool_get_config(ocnConfigs, 'config_use_activeTracers_surface_bulk_forcing', &
                                config_use_activeTracers_surface_bulk_forcing)

      useJerlov=.false.

      if(.not. config_use_activeTracers_surface_bulk_forcing) then
         if(trim(config_sw_absorption_type) .ne. 'none') then
            call mpas_log_write( &
               'You have specified bulk_forcing off with shortwave absorption on ' // &
               'either set config_sw_absorption_type to none or enable activeTracers_surface_bulk_forcing', &
                MPAS_LOG_CRIT)
             err = 1
         endif
         return
      endif

      if ( trim( config_sw_absorption_type ) == 'jerlov') then
         useJerlov=.true.
!        call ocn_tracer_short_wave_absorption_jerlov_init(err)
      else if ( trim( config_sw_absorption_type ) == 'ohlmann00' ) then
         call ocn_tracer_short_wave_absorption_variable_init(domain, err)
      else if ( trim( config_sw_absorption_type ) .ne. 'none') then
         call mpas_log_write( &
            'Incorrect option for config_sw_absorption_type. Options are: jerlov or ohlmann00 or none', &
            MPAS_LOG_CRIT)
         err = 1
         return
      end if


   end subroutine ocn_tracer_short_wave_absorption_init!}}}

!***********************************************************************


end module ocn_tracer_short_wave_absorption

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
