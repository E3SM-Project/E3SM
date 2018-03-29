! Copyright (c) 2013-2018,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_diagnostic_vars
!
!> \brief MPAS land ice module for calculating diagnostic variables
!> \author Matt Hoffman
!> \date   17 April 2011
!> \details
!>  This module contains various subroutines for
!>  calculating diagnostic (time-independent) variables
!>  for the land ice core.
!>  These calculations should be general so as to be independent
!>  of time integration scheme.
!
!-----------------------------------------------------------------------

module li_diagnostic_vars

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use mpas_log

   use li_mask
   use li_constants

   implicit none
   private

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
   public :: li_calculate_diagnostic_vars, &
             li_calculate_apparent_diffusivity, &
             li_calculate_flowParamA

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine li_calculate_diagnostic_vars
!
!> \brief   Calculates diagnostic variables
!> \author  Matthew Hoffman
!> \date    12 September 2013
!> \details
!>  This routine calculates diagnostic variables using the current prognostic
!>  variables.  These should only be variables that are output for convenience
!>  but otherwise not needed by the model.
!>
!
!-----------------------------------------------------------------------

   subroutine li_calculate_diagnostic_vars(domain, err)

      use li_thermal, only: li_compute_pressure_melting_point_fields

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object
      ! Note: domain is passed in because halo updates are needed in this routine
      ! and halo updates have to happen outside block loops, which requires domain.

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
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: thermalPool
      integer, pointer :: nCells
      real (kind=RKIND), dimension(:), pointer :: layerCenterSigma
      real (kind=RKIND), dimension(:), pointer :: thickness
      real (kind=RKIND), dimension(:,:), pointer :: pmpTemperature
      real (kind=RKIND), dimension(:), pointer :: basalPmpTemperature
      real (kind=RKIND), dimension(:), pointer :: thicknessOld
      real (kind=RKIND), dimension(:), pointer :: dHdt
      real (kind=RKIND), pointer :: deltat

      integer :: err_tmp


      err = 0

      block => domain % blocklist
      do while (associated(block))
          call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
          call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
          call mpas_pool_get_subpool(block % structs, 'thermal', thermalPool)

          call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

          call mpas_pool_get_array(meshPool, 'layerCenterSigma', layerCenterSigma)
          call mpas_pool_get_array(meshPool, 'deltat', deltat)
          call mpas_pool_get_array(geometryPool, 'thickness', thickness)
          call mpas_pool_get_array(geometryPool, 'thicknessOld', thicknessOld)
          call mpas_pool_get_array(geometryPool, 'dHdt', dHdt)
          call mpas_pool_get_array(thermalPool, 'pmpTemperature', pmpTemperature)
          call mpas_pool_get_array(thermalPool, 'basalPmpTemperature', basalPmpTemperature)

         ! Calculate diagnostic PMP temperature
         ! Note: this is only to allow these fields to be output;
         !       the thermal module performs these calculations internally as needed.
         ! Note: The PMP temperature calculated here will differ from the PMP values
         !       in the thermal solve, because the thermal solver occurs before
         !       advection occurs.
         call li_compute_pressure_melting_point_fields(nCells, thickness, layerCenterSigma, &
               pmpTemperature, basalPmpTemperature)

         ! Calculate dHdt.  Do it here after all thickness adjustments are complete
         dHdt = (thickness - thicknessOld) / deltat * scyr  ! Units are m/yr
         thicknessOld = thickness  ! Reset thicknessOld for next time step

         block => block % next
      end do


      ! === error check and exit
      if (err == 1) then
          print *, "An error has occurred in li_calculate_diagnostic_vars. Aborting..."
          !call mpas_dmpar_global_abort(dminfo)
      endif

   !--------------------------------------------------------------------
   end subroutine li_calculate_diagnostic_vars


!***********************************************************************
!
!  subroutine li_calculate_apparent_diffusivity
!
!> \brief   Computes apparent diffusivity
!> \author  Matt Hoffman
!> \date    19 April 2012
!> \details
!> This routine computes the apparent diffusivity.
!> Estimate diffusivity using the relation that the 2-d flux Q=-D grad h and Q=UH,
!> where h is surface elevation, D is diffusivity, U is 2-d velocity vector, and H is thickness
!> Solving for D = UH/-grad h
!> DCFL: dt = 0.5 * dx**2 / D = 0.5 * dx**2 * slopemag / flux_downslope
!
!-----------------------------------------------------------------------
   subroutine li_calculate_apparent_diffusivity(meshPool, velocityPool, scratchPool, geometryPool, allowableDiffDt)
      use mpas_vector_reconstruction

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      type (mpas_pool_type), intent(in) :: &
         velocityPool          !< Input: velocity information

      type (mpas_pool_type), intent(in) :: &
         scratchPool          !< Input: scratch information

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: &
         geometryPool          !< Input: geometry information

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      real(kind=RKIND), intent(out) :: allowableDiffDt !< Output: allowable timestep based on diffusive CFL

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      logical, pointer :: config_print_thickness_advection_info
      real (kind=RKIND), dimension(:), pointer :: normalSlopeEdge
      type (field1dReal), pointer :: slopeReconstructXField, slopeReconstructYField, slopeReconstructZField
      !< Only needed for calling mpas_reconstruct, but not actually used here
      type (field1dReal), pointer :: slopeCellAxis1Field
      type (field1dReal), pointer :: slopeCellAxis2Field
      real (kind=RKIND), dimension(:), pointer :: slopeCellAxis1, slopeCellAxis2
      real (kind=RKIND), dimension(:), pointer :: thickness
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness
      real (kind=RKIND), dimension(:,:), pointer :: uReconstructAxis1, uReconstructAxis2
      real (kind=RKIND), dimension(:), pointer :: apparentDiffusivity
      real (kind=RKIND), dimension(:), pointer :: dcEdge
      integer, dimension(:), pointer :: cellMask
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell
      integer, pointer :: nCells, nVertLevels
      real (kind=RKIND) :: allowableDtHere
      real (kind=RKIND) :: fluxVeloAxis1, fluxVeloAxis2
      real (kind=RKIND) :: fluxDownslope
      real (kind=RKIND) :: slopeCellMagnitude
      real (kind=RKIND) :: dCell
      integer :: iCell, iEdge, iLevel
      real (kind=RKIND), parameter :: bigNumber = 1.0e16_RKIND
      !<- This is ~300 million years in seconds, but it is small enough not too overflow
      real (kind=RKIND), parameter :: smallNumber = 1.0e-36_RKIND
      logical :: divideSingularityFound

      ! Note: This routine could be broken into 2: one to calculate diffusivity
      ! and another to get the diffusive CFL timestep.  In that case, the first (and possibly the second)
      ! could be moved to diagnostic_variable_solve_after_velocity.  However, since
      ! diffusivity is only used for this check, I don't think it makes sense to separate these
      ! calculations for now.

      ! get needed variables
      call mpas_pool_get_config(liConfigs, 'config_print_thickness_advection_info', config_print_thickness_advection_info)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(geometryPool, 'normalSlopeEdge', normalSlopeEdge)
      call mpas_pool_get_array(geometryPool, 'thickness', thickness, timeLevel=1)
      call mpas_pool_get_array(geometryPool, 'layerThickness', layerThickness, timeLevel=1)
      call mpas_pool_get_array(geometryPool, 'apparentDiffusivity', apparentDiffusivity)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_array(velocityPool, 'uReconstructZonal', uReconstructAxis1)
      call mpas_pool_get_array(velocityPool, 'uReconstructMeridional', uReconstructAxis2)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_field(scratchPool, 'workCell', slopeReconstructXField)
      call mpas_allocate_scratch_field(slopeReconstructXField, .true.)
      call mpas_pool_get_field(scratchPool, 'workCell2', slopeReconstructYField)
      call mpas_allocate_scratch_field(slopeReconstructYField, .true.)
      call mpas_pool_get_field(scratchPool, 'workCell3', slopeReconstructZField)
      call mpas_allocate_scratch_field(slopeReconstructZField, .true.)
      call mpas_pool_get_field(scratchPool, 'slopeCellX', slopeCellAxis1Field)
      call mpas_allocate_scratch_field(slopeCellAxis1Field, .true.)
      slopeCellAxis1 => slopeCellAxis1Field % array
      call mpas_pool_get_field(scratchPool, 'slopeCellY', slopeCellAxis2Field)
      call mpas_allocate_scratch_field(slopeCellAxis2Field, .true.)
      slopeCellAxis2 => slopeCellAxis2Field % array

      ! given thickness, compute layerThickness
      call li_calculate_layerThickness(meshPool, thickness, layerThickness)

      ! Initialize output
      allowableDiffDt = bigNumber

      ! Approximate slope at cell centers
      ! reconstruct routines set uReconstructZonal = uReconstructX; uReconstructMeridional = uReconstructY
      ! for planar meshes, so those variables can be used as orthogonal components of the vector
      ! in either the plane or sphere.  This avoids needing to add logic for if we are on a sphere or not.
      call mpas_reconstruct(meshPool, normalSlopeEdge,        &
                            slopeReconstructXField % array, slopeReconstructYField % array, slopeReconstructZField % array, &
                            slopeCellAxis1, slopeCellAxis2)


      ! Approximate flux at cell centers
      divideSingularityFound = .false.
      do iCell = 1, nCells
         slopeCellMagnitude = sqrt(slopeCellAxis1(iCell)**2 + slopeCellAxis2(iCell)**2) + smallNumber

         if ( (slopeCellMagnitude < 1.0e-4_RKIND) .and. &
              (max(maxval(uReconstructAxis1(:,iCell)), maxval(uReconstructAxis2(:,iCell))) < 3.18e-8_RKIND) ) then
              ! 3.18e-8=1 m/yr in m/s
            ! Ignore diffusivity near 'divide-singularities'
            apparentDiffusivity(iCell) = 0.0_RKIND
            divideSingularityFound = .true.
         else
            fluxDownslope = 0.0_RKIND
            do iLevel = 1, nVertLevels
               fluxVeloAxis1 = (uReconstructAxis1(iLevel, iCell) + uReconstructAxis1(iLevel+1, iCell)) * 0.5_RKIND
               fluxVeloAxis2 = (uReconstructAxis2(iLevel, iCell) + uReconstructAxis2(iLevel+1, iCell)) * 0.5_RKIND
               fluxDownslope = fluxDownslope + &
                               (-1.0_RKIND * slopeCellAxis1(iCell) * fluxVeloAxis1 - slopeCellAxis2(iCell) * fluxVeloAxis2) &
                               * layerThickness(iLevel, iCell) / slopeCellMagnitude
            enddo
            apparentDiffusivity(iCell) = abs(fluxDownslope) / slopeCellMagnitude
         endif

         ! Calculate allowable timestep based on DCFL
         if ( li_mask_is_grounded_ice(cellMask(iCell)) .and. li_mask_is_dynamic_ice(cellMask(iCell)) ) then
            ! Find shortest distance to a neighboring cell center, dCell
            dCell = minval(dcEdge(1:nEdgesOnCell(iCell)))
            allowableDtHere = 0.5_RKIND * dCell**2 / (apparentDiffusivity(iCell) + smallNumber)
         else
            allowableDtHere = bigNumber
         endif
         allowableDiffDt = min(allowableDiffDt, allowableDtHere)
      enddo

      if (divideSingularityFound .and. config_print_thickness_advection_info) then
         call mpas_log_write('Notice: In calculating apparentDiffusivity, one or more cells have been ignored ' // &
            'due to flat slope and low velocity (assumed to be a divide where diffusivity is undefined).')
      endif

      call mpas_deallocate_scratch_field(slopeReconstructXField, .true.)
      call mpas_deallocate_scratch_field(slopeReconstructYField, .true.)
      call mpas_deallocate_scratch_field(slopeReconstructZField, .true.)
      call mpas_deallocate_scratch_field(slopeCellAxis1Field, .true.)
      call mpas_deallocate_scratch_field(slopeCellAxis2Field, .true.)
   !--------------------------------------------------------------------
   end subroutine li_calculate_apparent_diffusivity



!***********************************************************************
!
!  subroutine li_calculate_flowParamA
!
!> \brief   Calculates the flow law parameter A based on temperature
!> \author  Matt Hoffman
!> \date    23 Jan 2014
!> \details
!>  This routine calculates the flow law parameter A based on temperature
!>  depending on what option is chosen.
!>  The default option is a constant A assigned from config_default_flowParamA.
!>  The PB1982 option uses this equation from \emph{Paterson and Budd} [1982]
!>  and \emph{Paterson} [1994] (copied from CISM):
!>   \[
!>   A(T^{*})=A0 \exp \left(\frac{-Q}{RT^{*}}\right)
!>   \]
!>   This is equation 9 in {\em Payne and Dongelmans}. $A0$ is a constant of proportionality,
!>   $Q$ is the activation energy for for ice creep, and $R$ is the universal gas constant.
!>   The pressure-corrected temperature, $T^{*}$ is given by:
!>   \[
!>   T^{*} = T - T_{pmp} + T_0
!>   \]
!>   \[
!>   T_{pmp} = T_0 - \sigma \rho g H \Phi
!>   \]
!>   $T$ is the ice temperature, $T_0$ is the triple point of water,
!>   $\rho$ is the ice density, and $\Phi$ is the (constant) rate of change of
!>   melting point temperature with pressure.
!>
!>  The CP2010 option uses this equation from the 4th Edition of Physics of Glaciers (Eq. 3.35):
!>   \[
!>   A(T^{*})=A0 \exp \left(\frac{-Q}{R} ( \frac{1}{T^{*}} - \frac{1}{T_t})\right)
!>   \]
!>   where the variables are the same as above and $T_t$ is the pressure corrected
!>   transition temperature (-10 deg C at 0 pressure).
!>   Values for $A0, Q, \Phi$ differ from PB1982.
!>
!>  All options are adjusted by the enhancement factor (which defaults to 1.0).
!-----------------------------------------------------------------------
   subroutine li_calculate_flowParamA(meshPool, temperature, thickness, flowParamA, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         temperature    !< Input: temperature
      real (kind=RKIND), dimension(:), intent(in) :: &
         thickness    !< Input: thickness

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      integer, intent(out) :: err

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         flowParamA    !< Input: flowParamA

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      integer, pointer :: nCells, nVertLevels
      character (len=StrKIND), pointer :: config_flowParamA_calculation
      real (kind=RKIND), pointer :: config_default_flowParamA, &
                                    config_enhancementFactor,  &
                                    config_dynamic_thickness,  &
                                    config_ice_density
      integer :: iCell, iLevel, err_tmp
      real (kind=RKIND), dimension(:), pointer :: layerCenterSigma
      real (kind=RKIND) :: A0, Q, pressureMeltPointSlope
      real (kind=RKIND) :: temperatureCorrected, transitionTemperatureCorrected

      err = 0
      err_tmp = 0

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'layerCenterSigma', layerCenterSigma)

      call mpas_pool_get_config(liConfigs, 'config_flowParamA_calculation', config_flowParamA_calculation)
      call mpas_pool_get_config(liConfigs, 'config_enhancementFactor', config_enhancementFactor)
      call mpas_pool_get_config(liConfigs, 'config_default_flowParamA', config_default_flowParamA)
      call mpas_pool_get_config(liConfigs, 'config_dynamic_thickness', config_dynamic_thickness)
      call mpas_pool_get_config(liConfigs, 'config_ice_density', config_ice_density)


      select case(config_flowParamA_calculation)
      !-----------------------------------------------------------------
      case('constant')
          flowParamA = config_default_flowParamA
      !-----------------------------------------------------------------
      case('PB1982')
          pressureMeltPointSlope = 9.7456e-8_RKIND
          do iCell = 1, nCells
            if (thickness(iCell) > 0.0_RKIND) then
              do iLevel = 1, nVertLevels
                ! Calculate the pressure-corrected temperature
                temperatureCorrected = min(273.15_RKIND, temperature(iLevel,iCell) + pressureMeltPointSlope * &
                         thickness(iCell) * config_ice_density * gravity * layerCenterSigma(iLevel) )
                temperatureCorrected = max(223.15_RKIND, temperatureCorrected)
                ! Calculate flow A
                if (temperatureCorrected > 263.15_RKIND) then
                  A0 = 1.733e3_RKIND
                  Q  = 139.0e3_RKIND
                else
                  A0 = 3.613e-13_RKIND
                  Q  = 60.0e3_RKIND
                endif
                flowParamA(iLevel,iCell) = A0 * exp(-1.0_RKIND * Q / (idealGasConstant * temperatureCorrected))
              enddo ! levels
            else
              flowParamA(:,iCell) = 0.0_RKIND  ! non-ice cells get 0
            endif ! if dynamic ice
          enddo ! cells
      !-----------------------------------------------------------------
      case('CP2010')
          pressureMeltPointSlope = 7.0e-8_RKIND
          do iCell = 1, nCells
            if (thickness(iCell) > 0.0_RKIND) then  ! SIA solver could make use of A on thin ice
              ! if doing 2nd order averaging of flwa onto edges (otherwise this could be the dynamic thickness limit)
              do iLevel = 1, nVertLevels
                ! Calculate the pressure-corrected temperature
                temperatureCorrected = min(273.15_RKIND, temperature(iLevel,iCell) + pressureMeltPointSlope * &
                         thickness(iCell) * config_ice_density * gravity * layerCenterSigma(iLevel) )
                temperatureCorrected = max(223.15_RKIND, temperatureCorrected)
                transitionTemperatureCorrected = 263.15_RKIND + pressureMeltPointSlope * &
                         thickness(iCell) * config_ice_density * gravity * layerCenterSigma(iLevel)
                ! Calculate flow A
                A0 = 3.5e-25_RKIND
                if (temperatureCorrected > 263.15_RKIND) then
                  Q  = 115.0e3_RKIND
                else
                  Q  = 6.0e4_RKIND
                endif
                flowParamA(iLevel,iCell) = A0 * exp(-1.0_RKIND * Q / idealGasConstant * &
                   (1.0_RKIND/temperatureCorrected - 1.0_RKIND/transitionTemperatureCorrected))
              enddo ! levels
            else
              flowParamA(:,iCell) = 0.0_RKIND  ! non-ice cells get 0
            endif ! if dynamic ice
          enddo ! cells
      !-----------------------------------------------------------------
      end select

      !print *,'max flwa', maxval(flowParamA)
      !print *,'config_enhancementFactor', config_enhancementFactor

      ! Include enhancement factor
      flowParamA = flowParamA * config_enhancementFactor

      err = ior(err, err_tmp)

   end subroutine li_calculate_flowParamA


!***********************************************************************
!***********************************************************************
! Private subroutines:
!***********************************************************************
!***********************************************************************



end module li_diagnostic_vars

