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

   use mpas_grid_types
   use ocn_constants
   use ocn_tracer_short_wave_absorption_jerlov

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

   subroutine ocn_tracer_short_wave_absorption_tend(meshPool, index_temperature, layerThickness, penetrativeTemperatureFlux, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      real (kind=RKIND), dimension(:), intent(in) :: &
         penetrativeTemperatureFlux !< Input: short wave heat flux

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

      err = 0

      if ( useJerlov ) then
         call ocn_tracer_short_wave_absorption_jerlov_tend(meshPool, index_temperature, layerThickness, penetrativeTemperatureFlux, tend, err)
      end if

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

   subroutine ocn_tracer_short_wave_absorption_init(err)!{{{

   !--------------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      character (len=StrKind), pointer :: config_sw_absorption_type

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_sw_absorption_type', config_sw_absorption_type)

      useJerlov = .false.

      if ( trim( config_sw_absorption_type ) .ne. 'jerlov') then
        write(0,*) 'Incorrect option for config_sw_absorption_type. Options are: jerlov'
        err = 1
        return
      else if ( trim( config_sw_absorption_type ) == 'jerlov') then
        useJerlov = .true.
        call ocn_tracer_short_wave_absorption_jerlov_init(err)
      end if


   end subroutine ocn_tracer_short_wave_absorption_init!}}}

!***********************************************************************

end module ocn_tracer_short_wave_absorption

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
