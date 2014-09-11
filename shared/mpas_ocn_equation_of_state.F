! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_equation_of_state
!
!> \brief MPAS ocean equation of state driver
!> \author Mark Petersen
!> \date   September 2011
!> \details
!>  This module contains the main driver routine for calling
!>  the equation of state.
!
!-----------------------------------------------------------------------

module ocn_equation_of_state

   use mpas_kind_types
   use mpas_grid_types
   use ocn_equation_of_state_linear
   use ocn_equation_of_state_jm
   use ocn_constants
   use mpas_io_units

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

   public :: ocn_equation_of_state_density, &
             ocn_equation_of_state_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: linearEos, jmEos


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_equation_of_state
!
!> \brief   Calls equation of state
!> \author  Mark Petersen
!> \date    September 2011
!> \details 
!>  This routine calls the equation of state to update the density
!
!-----------------------------------------------------------------------

   subroutine ocn_equation_of_state_density(statePool, diagnosticsPool, meshPool, k_displaced, displacement_type, density, err, &
      thermalExpansionCoeff, salineContractionCoeff, timeLevelIn)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  This module contains routines necessary for computing the density
   !  from model temperature and salinity using an equation of state.
   !
   ! Input: mesh - mesh metadata
   !        s - state: tracers
   !        k_displaced 
   !
   !  If k_displaced==0, density is returned with no displacement 
   !
   !  If k_displaced~=0, density is returned, and is for
   !  a parcel adiabatically displaced from its original level to level 
   !  k_displaced.  When using the linear EOS, state % displacedDensity is 
   !  still filled, but depth (i.e. pressure) does not modify the output.
   !
   ! Output: s - state: computed density
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none

      type (mpas_pool_type), intent(in) :: statePool
      type (mpas_pool_type), intent(inout) :: diagnosticsPool
      type (mpas_pool_type), intent(in) :: meshPool
      integer, intent(in), optional :: timeLevelIn
      integer :: k_displaced
      character(len=*), intent(in) :: displacement_type
      real (kind=RKIND), dimension(:,:), intent(out) :: density
      integer, intent(out) :: err
      real (kind=RKIND), dimension(:,:), intent(out), optional :: &
         thermalExpansionCoeff,  &! Thermal expansion coefficient (alpha), defined as $-1/\rho d\rho/dT$ (note negative sign)
         salineContractionCoeff   ! Saline contraction coefficient (beta), defined as $1/\rho d\rho/dS$

      integer, dimension(:), pointer :: maxLevelCell
      real (kind=RKIND), dimension(:,:), pointer :: tracersSurfaceValue
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers
      integer :: iCell, k
      integer, pointer :: indexT, indexS
      type (dm_info) :: dminfo
      integer :: timeLevel

      err = 0

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_array(diagnosticsPool, 'tracersSurfaceValue', tracersSurfaceValue)
      call mpas_pool_get_array(statePool, 'tracers', tracers, timeLevel)
      call mpas_pool_get_dimension(statePool, 'index_temperature', indexT)
      call mpas_pool_get_dimension(statePool, 'index_salinity', indexS)

      if (linearEos) then

         call ocn_equation_of_state_linear_density(meshPool, indexT, indexS, tracers, density, err, &
            thermalExpansionCoeff, salineContractionCoeff)

      elseif (jmEos) then

         call ocn_equation_of_state_jm_density(meshPool, k_displaced, displacement_type, indexT, indexS, tracers, density, err, &
            tracersSurfaceValue, thermalExpansionCoeff, salineContractionCoeff)

      endif

   end subroutine ocn_equation_of_state_density!}}}

!***********************************************************************
!
!  routine ocn_equation_of_stateInit
!
!> \brief   Initializes ocean momentum horizontal mixing quantities
!> \author  Mark Petersen
!> \date    September 2011
!> \details 
!>  This routine initializes a variety of quantities related to 
!>  horizontal velocity mixing in the ocean. Since a variety of 
!>  parameterizations are available, this routine primarily calls the
!>  individual init routines for each parameterization. 
!
!----------------------------------------------------------------------

   subroutine ocn_equation_of_state_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err

      character (len=StrKIND), pointer :: config_eos_type

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_eos_type', config_eos_type)

      linearEos = .false.
      jmEos = .false.

      if (config_eos_type.eq.'linear') then
         linearEos = .true.
         call ocn_equation_of_state_linear_init(err)
      elseif (config_eos_type.eq.'jm') then
         jmEos = .true.
         call ocn_equation_of_state_jm_init(err)
      else
         write (stderrUnit,*) 'Invalid choice for config_eos_type.'
         write (stderrUnit,*) '  Choices are: linear, jm'
         err = 1
      endif

   !--------------------------------------------------------------------

   end subroutine ocn_equation_of_state_init!}}}

!***********************************************************************

end module ocn_equation_of_state

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
