! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_advection
!
!> \brief MPAS land ice advection solver
!> \author William Lipscomb
!> \date   December 2015
!> \details
!>  This module contains the routines for advecting ice thickness
!>  and tracers for land ice.
!-----------------------------------------------------------------------

module li_advection

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar
   use li_setup
   use li_mask

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
   public :: li_advection_thickness_tracers

   !TODO - can make the following private once the old advection scheme is removed
   public :: li_vertical_remap, tracer_setup, tracer_finish  
 
   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   !WHL - debug                                                                                                                                                    
   logical, parameter :: new_thickness_advection = .true.
!   logical, parameter :: new_thickness_advection = .false. 
   public :: new_thickness_advection

 contains


!***********************************************************************
!
!  subroutine li_advection_thickness_tracers
!
!> \brief   Advection for ice thickness and tracers
!> \author  William Lipscomb
!> \date    November 2015
!> \details
!>  This routine (1) computes new values of ice thickness and tracers under
!>  horizontal advection, (2) applies surface and basal mass balance
!>  terms, and (3) remaps tracers back onto the standard vertical grid.
!>  Based on subroutine li_tendency_thickness by Matthew Hoffman
!>  Note: The layerThickness and tracer fields must have had halo updates
!>        before calling this subroutine.
!-----------------------------------------------------------------------

   subroutine li_advection_thickness_tracers(&
        dt,                     &
        meshPool,               &
        velocityPool,           &
        geometryPool,           &
        thermalPool,            &
        scratchPool,            &
        layerThicknessOld,      &
        layerThicknessNew,      &
        allowableDt,            &
        err,                    &
        advectTracersIn)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(in) :: &
           dt                     !< Input: time step (s)

      type (mpas_pool_type), intent(in) :: &
           meshPool               !< Input: mesh information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
           layerThicknessOld      !< Input: old layer thickness

      logical, intent(in), optional :: &
           advectTracersIn        !< Input: if true, then advect tracers as well as thickness
                                  !         (default is true)

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: &
           velocityPool           !< Input/output: velocity information
                                  !                (needs to be inout for li_calculate_mask call

      type (mpas_pool_type), intent(inout) :: &
           geometryPool           !< Input/output: geometry information to be updated

      type (mpas_pool_type), intent(inout) :: &
           thermalPool            !< Input/output: temperature/enthalpy information to be updated

      type (mpas_pool_type), intent(inout) :: &
           scratchPool            !< Input/output: work arrays for tracer advection

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
           layerThicknessNew      !< Input: new layer thickness

      real (kind=RKIND), intent(out) :: &
           allowableDt            !< Output: maximum allowable dt based on CFL condition

      integer, intent(out) :: &
           err                    !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer, pointer :: nVertLevels

      real (kind=RKIND), dimension(:), pointer :: &
           thickness,             & ! ice thickness
           sfcMassBal,            & ! surface mass balance
           basalMassBal,          & ! basal mass balance
           groundedBasalMassBal,  & ! basal mass balance for grounded ice
           floatingBasalMassBal     ! basal mass balance for floating ice

      real (kind=RKIND), dimension(:,:), pointer :: &
           temperature,           & ! interior ice temperature
           waterfrac,             & ! interior water fraction
           enthalpy                 ! interior ice enthalpy

      real (kind=RKIND), dimension(:,:), pointer :: &
           normalVelocity,        & ! normal component of velocity on cell edges
           layerThicknessEdge       ! layer thickness on upstream edge of cell

      real (kind=RKIND), dimension(:,:,:), pointer :: &
           advectedTracersOld,    & ! old values of advected tracers
           advectedTracersNew       ! new values of advected tracers

      type (field3DReal), pointer :: &
           advectedTracersOldField, & ! scratch field containing old values of advected tracers
           advectedTracersNewField    ! scratch field containing new values of advected tracers

      integer, dimension(:), pointer :: &
           cellMask                 ! integer bitmask for cells

      character (len=StrKIND), pointer :: &
           config_thickness_advection   ! method for advecting thickness and tracers

      logical, pointer :: &
           config_print_thickness_advection_info   !TODO - change to config_print_advection_info?

      real (kind=RKIND), pointer :: &
           config_ice_density

      logical :: advectTracers     ! if true, then advect tracers as well as thickness

      integer :: err_tmp

      !WHL - debug
      integer, dimension(:), pointer :: indexToCellID
      integer, pointer :: nCells
      integer, pointer :: config_stats_cell_ID
      integer :: iEdge, iEdgeOnCell
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell

      !WHL - debug - for test-case diagnostics
      logical, parameter :: dome_test = .true.
      logical, parameter :: circular_shelf_test = .false.
      integer :: nCellsPerRow
      integer :: nRows
      integer :: i, k, iRow, iCell, iLayer

      if (circular_shelf_test) then
         nCellsPerRow = 40
         nRows = 46
      elseif (dome_test) then
         nCellsPerRow = 30
         nRows = 34
      endif

      err = 0
      err_tmp = 0

      ! Decide whether to advect tracers
      ! If not, then advect thickness only

      if (present(advectTracersIn)) then
         advectTracers = advectTracersIn
      else
         advectTracers = .false.
      endif

      ! get dimensions
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      ! get arrays from the geometry pool
      call mpas_pool_get_array(geometryPool, 'sfcMassBal', sfcMassBal)
      call mpas_pool_get_array(geometryPool, 'basalMassBal', basalMassBal)
      call mpas_pool_get_array(geometryPool, 'groundedBasalMassBal', groundedBasalMassBal)
      call mpas_pool_get_array(geometryPool, 'floatingBasalMassBal', floatingBasalMassBal)
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(geometryPool, 'layerThicknessEdge', layerThicknessEdge)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)

      ! get arrays from the velocity pool
      call mpas_pool_get_array(velocityPool, 'normalVelocity', normalVelocity)

      ! get arrays from the thermal pool
      call mpas_pool_get_array(thermalPool, 'temperature', temperature)
      call mpas_pool_get_array(thermalPool, 'waterfrac', waterfrac)
      call mpas_pool_get_array(thermalPool, 'enthalpy', enthalpy)

      ! get config variables
      call mpas_pool_get_config(liConfigs, 'config_thickness_advection', config_thickness_advection)
      call mpas_pool_get_config(liConfigs, 'config_print_thickness_advection_info', config_print_thickness_advection_info)
      call mpas_pool_get_config(liConfigs, 'config_ice_density', config_ice_density)

      !WHL - debug
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
      call mpas_pool_get_config(liConfigs, 'config_stats_cell_ID', config_stats_cell_ID)

      ! get fields from the scratch pool
      ! Note: The old/advectedTracersNew fields are scratch fields containing only the tracers that need to be advected.
      !       For example, they will contain temperature or enthalpy (depending on config_thermal_solver) but not both.

      call mpas_pool_get_field(scratchPool, 'workTracerLevelCell', advectedTracersOldField)
      call mpas_allocate_scratch_field(advectedTracersOldField, .true.)
      advectedTracersOld => advectedTracersOldField % array

      call mpas_pool_get_field(scratchPool, 'workTracerLevelCell2', advectedTracersNewField)
      call mpas_allocate_scratch_field(advectedTracersNewField, .true.)
      advectedTracersNew => advectedTracersNewField % array

      !-----------------------------------------------------------------
      ! Horizontal transport of thickness and tracers
      !-----------------------------------------------------------------

      select case (trim(config_thickness_advection))

      case ('fo')

         if (config_print_thickness_advection_info) then
            write(stderrUnit,*) 'Using first-order upwind for thickness and tracer advection'
         endif
         
         !WHL - debug
         if (dome_test .and. config_print_thickness_advection_info) then
            write(stderrUnit,*) ' '
            write(stderrUnit,*) 'Dome test diagnostics:'
            write(stderrUnit,*) ' '
            iLayer = 1
!            write(stderrUnit,*) 'Initial thickness, layer =', iLayer
            do iRow = nRows, 1, -1
               if (mod(iRow,2) == 0) then  ! indent for even-numbered rows
!                  write(stderrUnit,'(a3)',advance='no') '    '
               endif
               do i = nCellsPerRow/2 - 2, nCellsPerRow
                  iCell = (iRow-1)*nCellsPerRow + i
!                  write(stderrUnit,'(f8.3)',advance='no') oldLayerThickness(iLayer,iCell)
               enddo
!               write(stderrUnit,*) ' '
            enddo
         endif

         if (advectTracers) then

            ! Copy the old tracer values into the advectedTracersOld array.
            ! This requires setting a tracer pointer to either temperature or enthalpy,
            !  depending on the value of config_thermal_solver.
            ! Currently (as of Nov. 2015), this is the only tracer transported, but others
            !  may be added later.

            call tracer_setup(&
                 meshPool,          &
                 geometryPool,      &
                 thermalPool,       &
                 advectedTracersOld)

            if (dome_test .and. config_print_thickness_advection_info) then
!               write(stderrUnit,*) ' '
!               write(stderrUnit,*) 'Initial tracer 1, layer =', iLayer
               iLayer = 1
               do iRow = nRows, 1, -1
                  if (mod(iRow,2) == 0) then  ! indent for even-numbered rows
!                     write(stderrUnit,'(a4)',advance='no') '    '
                  endif
                  do i = nCellsPerRow/2 - 2, nCellsPerRow
                     iCell = (iRow-1)*nCellsPerRow + i
!                     write(stderrUnit,'(f10.4)',advance='no') advectedTracersOld(1,iLayer,iCell)
                  enddo
!                  write(stderrUnit,*) ' '
               enddo
               
            endif

         else   ! no tracer advection; pass a tracer array full of zeroes

            advectedTracersOld(:,:,:) = 0.0_RKIND

         endif

         ! Transport thickness and tracers using a first-order upwind scheme
         ! Note:  For the enthalpy scheme, temperature and waterfrac are the primary prognostic
         !        variables to be updated, but enthalpy is the advected tracer (for reasons of
         !        energy conservation).

         if (config_print_thickness_advection_info) then
            do iCell = 1, nCells
               if (indexToCellID(iCell) == config_stats_cell_ID) then
                  write(stderrUnit,*) 'Before thickness advection, iCell (local, global), thickness:', &
                       iCell, indexToCellID(iCell), sum(layerThicknessOld(:,iCell))
                  do iEdgeOnCell = 1, nEdgesOnCell(iCell)
                     iEdge = edgesOnCell(iEdgeOnCell,iCell)
                     write(stderrUnit,*) 'iEdge, normalVelocity:', iEdge, normalVelocity(1,iEdge)
                  enddo
               endif
            enddo
         endif

         call advect_thickness_tracers_upwind(&
              dt,                      &
              meshPool,                &
              normalVelocity,          &
              layerThicknessEdge,      & 
              layerThicknessOld,       &
              advectedTracersOld,      &
              layerThicknessNew,       &
              advectedTracersNew,      &
              allowableDt,             &
              err)

         if (config_print_thickness_advection_info) then
            do iCell = 1, nCells
               if (indexToCellID(iCell) == config_stats_cell_ID) then
                  write(stderrUnit,*) ' '
                  write(stderrUnit,*) 'After thickness advection, iCell, new layer thickness:'
                  do k = 1, nVertLevels
                     write(stderrUnit,*) iCell, layerThicknessNew(k,iCell)
                  enddo
               endif
            enddo
         endif

         !WHL - debug
         if (dome_test .and. config_print_thickness_advection_info) then
!            write(stderrUnit,*) ' '
            iLayer = 1
!            write(stderrUnit,*) 'New thickness, layer =', iLayer
            do iRow = nRows, 1, -1
               if (mod(iRow,2) == 0) then  ! indent for even-numbered rows
!                  write(stderrUnit,'(a3)',advance='no') '    '
               endif
               do i = nCellsPerRow/2 - 2, nCellsPerRow
                  iCell = (iRow-1)*nCellsPerRow + i
!                  write(stderrUnit,'(f8.3)',advance='no') newLayerThickness(iLayer,iCell)
               enddo
!               write(stderrUnit,*) ' '
            enddo
!            write(stderrUnit,*) ' '
!            write(stderrUnit,*) 'New tracer 1, layer =', iLayer
            iLayer = 1
            do iRow = nRows, 1, -1
               if (mod(iRow,2) == 0) then  ! indent for even-numbered rows
!                  write(stderrUnit,'(a4)',advance='no') '    '
               endif
               do i = nCellsPerRow/2 - 2, nCellsPerRow
                  iCell = (iRow-1)*nCellsPerRow + i
!                  write(stderrUnit,'(f10.4)',advance='no') advectedTracersNew(1,iLayer,iCell)
               enddo
!               write(stderrUnit,*) ' '
            enddo
         endif

         !-----------------------------------------------------------------
         ! Add the surface and basal mass balance to the layer thickness
         !-----------------------------------------------------------------

         !TODO - Will need to treat temperature/enthalpy in a more complex way to conserve energy
         !        when mass is added at the upper and lower surfaces.
         !       For now, we are effectively assuming that ice is added at the surface with the top-layer temperature
         !        and at the bed with the bottom-layer temperature.

         ! Combine various basal mass balance fields based on mask
         ! Grounded and floating basal mass balance should come from the temperature solver.

         ! TODO: more complicated treatment at GL?

         where ( li_mask_is_grounded_ice(cellMask) )

            basalMassBal = groundedBasalMassBal

         elsewhere ( li_mask_is_floating_ice(cellMask) )

            ! Currently, floating and grounded ice are mutually exclusive.
            ! This could change if the GL is parameterized, in which case this logic may need adjustment.
            basalMassBal = floatingBasalMassBal

         elsewhere ( .not. (li_mask_is_ice(cellMask) ) )

            ! We don't allow a positive basal mass balance where ice is not already present.
            basalMassBal = 0.0_RKIND

         end where

         ! Add surface mass balance to new top layer thickness
         ! TODO: Need to decide how to deal with negative SMB that eliminates top layer or all ice (check for negative thickness?)

         layerThicknessNew(1,:) = layerThicknessNew(1,:) + (sfcMassBal/config_ice_density)*dt

         ! Add basal mass balance to new bottom layer thickness
         ! TODO: Need to decide how to deal with negative BMB that eliminates top layer or all ice (check for negative thickness?)

         layerThicknessNew(nVertLevels,:) = layerThicknessNew(nVertLevels,:) + (basalMassBal/config_ice_density)*dt

         ! Update the thickness and cellMask before the vertical remap calculation

         thickness = sum(LayerThicknessNew, 1)

         call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
         
         ! Remap tracers to the standard vertical sigma coordinate
         ! Note: If tracers are not being advected, then this subroutine simply restores the
         !       layer thickness to sigma coordinate values.

         call li_vertical_remap(thickness, cellMask, meshPool, layerThicknessNew, advectedTracersNew, err_tmp)
         err = ior(err, err_tmp)

         if (advectTracers) then

            ! Copy the advectedTracersNew values into the thermal tracer arrays
            ! (temperature, waterfrac, enthalpy)

            call tracer_finish(&
                 meshPool,              &
                 geometryPool,          &
                 thermalPool,           &
                 advectedTracersNew)

         endif

      case ('none')

         ! Do nothing

      case default

         write(stderrUnit,*) trim(config_thickness_advection), ' is not a valid option for thickness/tracer advection.'
         err_tmp = 1

      end select

      err = ior(err,err_tmp)

      ! clean up
      call mpas_deallocate_scratch_field(advectedTracersOldField, .true.)
      call mpas_deallocate_scratch_field(advectedTracersNewField, .true.)

      ! === error check
      if (err > 0) then
         write (stderrUnit,*) "An error has occurred in li_advection_thickness_tracers"
      endif

   !--------------------------------------------------------------------
    end subroutine li_advection_thickness_tracers

!***********************************************************************
! Private subroutines:
!***********************************************************************

!***********************************************************************
!
!  subroutine tracer_setup
!
!> \brief   Assemble a 3D array for tracer advection
!> \author  William Lipscomb
!> \date    November 2015
!> \details 
!>  This routine assembles a 3D array for tracer advection.
!>  Each tracer in the array is transported conservatively, along
!>  with the layer thickness field.
!
!-----------------------------------------------------------------------

    subroutine tracer_setup(&
         meshPool,       &
         geometryPool,   &
         thermalPool,    &
         advectedTracers)

      use li_thermal, only: li_temperature_to_enthalpy

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
           meshPool             !< Input: mesh information

      type (mpas_pool_type), intent(in) :: &
           geometryPool         !< Input: geometry information

      type (mpas_pool_type), intent(inout) :: &
           thermalPool          !< Input: temperature/enthalpy information

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(out) :: &
           advectedTracers  ! tracers to be advected
                            ! dimension 1 = maxTracers, 2 = nVertLevels, 3 = nCells

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      character (len=StrKIND), pointer :: &
           config_thermal_solver

      integer, pointer :: &
           nCellsSolve                ! number of locally owned cells

      real (kind=RKIND), dimension(:,:), pointer :: &
           temperature,           & ! interior ice temperature
           waterfrac,             & ! interior water fraction
           enthalpy                 ! interior ice enthalpy

      real (kind=RKIND), dimension(:), pointer :: &
           layerCenterSigma           ! sigma coordinate at midpoint of each layer

      real (kind=RKIND), dimension(:), pointer :: &
           thickness                  ! ice thickness

      integer :: iCell, iTracer

      ! get dimensions
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

      ! get arrays from mesh pool
      call mpas_pool_get_array(meshPool, 'layerCenterSigma', layerCenterSigma)

      ! get arrays from geometry pool
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)

      ! get arrays from thermal pool
      call mpas_pool_get_array(thermalPool, 'temperature', temperature)
      call mpas_pool_get_array(thermalPool, 'waterfrac', waterfrac)
      call mpas_pool_get_array(thermalPool, 'enthalpy', enthalpy)

      ! get config variables
      call mpas_pool_get_config(liConfigs, 'config_thermal_solver', config_thermal_solver)

      ! initialize
      advectedTracers(:,:,:) = 0.0_RKIND

      ! Note: For the enthalpy solver, it is necessary to transport enthalpy (rather than
      !       temperature and waterfrac separately) in order to conserve energy.
      ! Note: Temperature and waterfrac must be up to date (including halo cells)
      !        before calling this subroutine.

      iTracer = 1    ! temperature or enthalpy is the first index of the tracers array

      if (trim(config_thermal_solver) == 'enthalpy') then  ! advect enthalpy

         ! Rather than assume that enthalpy is up to date, recompute it from temperature and waterfrac

         do iCell = 1, nCellsSolve

            call li_temperature_to_enthalpy(&
                 layerCenterSigma,      &
                 thickness(iCell),      &
                 temperature(:,iCell),  &
                 waterfrac(:,iCell),    &
                 enthalpy(:,iCell))

         enddo

         advectedTracers(iTracer,:,:) = enthalpy(:,:)

      elseif (trim(config_thermal_solver) == 'temperature') then  ! advect temperature

         advectedTracers(iTracer,:,:) = temperature(:,:)

      endif

      ! Note: Other tracers (e.g., ice age or damage) can be added here as needed.
      !       May need to increase maxTracers in the Registry.

    end subroutine tracer_setup

!***********************************************************************
!
!  subroutine tracer_finish
!
!> \brief   Copy the new tracer values into prognostic arrays
!> \author  William Lipscomb
!> \date    November 2015
!> \details
!>  Given the new tracers after advection, this routine copies values
!>  into the appropriate prognostic arrays (e.g, temperature or enthalpy).
!
!-----------------------------------------------------------------------

    subroutine tracer_finish(&
         meshPool,           &
         geometryPool,       &
         thermalPool,        &
         advectedTracers)

      use li_thermal, only: li_enthalpy_to_temperature

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
           advectedTracers      !< Input: advected tracers
                                !         dimension 1 = maxTracers, 2 = nVertLevels, 3 = nCells

      type (mpas_pool_type), intent(in) :: &
           meshPool             !< Input: mesh information

      type (mpas_pool_type), intent(in) :: &
           geometryPool         !< Input: geometry information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: &
           thermalPool          !< Input/output: temperature/enthalpy information

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      character (len=StrKIND), pointer :: &
           config_thermal_solver

      integer, pointer :: &
           nCellsSolve                ! number of locally owned cells

      real (kind=RKIND), dimension(:), pointer :: &
           layerCenterSigma           ! sigma coordinate at midpoint of each layer

      real (kind=RKIND), dimension(:), pointer :: &
           thickness                  ! ice thickness

      real (kind=RKIND), dimension(:,:), pointer :: &
           temperature,             & ! interior ice temperature
           waterfrac,               & ! water fraction
           enthalpy                   ! interior ice enthalpy

      integer :: iCell, iTracer

      ! get dimensions
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

      ! get arrays from mesh pool
      call mpas_pool_get_array(meshPool, 'layerCenterSigma', layerCenterSigma)

      ! get arrays from geometry pool
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)

      ! get arrays from thermal pool
      call mpas_pool_get_array(thermalPool, 'temperature', temperature)
      call mpas_pool_get_array(thermalPool, 'waterfrac', waterfrac)
      call mpas_pool_get_array(thermalPool, 'enthalpy', enthalpy)

      ! get config variables
      call mpas_pool_get_config(liConfigs, 'config_thermal_solver', config_thermal_solver)

      iTracer = 1

      if (trim(config_thermal_solver) == 'enthalpy') then

         ! update the enthalpy
         enthalpy(:,:) = advectedTracers(iTracer,:,:)

         ! given the enthalpy, compute the temperature and water fraction

         do iCell = 1, nCellsSolve

            call li_enthalpy_to_temperature(&
                 layerCenterSigma,                           &
                 thickness(iCell),                           &
                 enthalpy(:,iCell),     &
                 temperature(:,iCell),  &
                 waterfrac(:,iCell))

         enddo

      elseif (trim(config_thermal_solver) == 'temperature') then

         ! update the temperature
         temperature(:,:) = advectedTracers(iTracer,:,:)

      endif

      ! Note: Other tracers (e.g., ice age or damage) can be added here as needed.
      !       May need to increase maxTracers in the Registry.

    end subroutine tracer_finish

!***********************************************************************

!  routine advect_thickness_tracers_upwind
!
!> \brief   Advect thickness and tracers using a first-order upwind scheme.
!> \author  William Lipscomb
!> \date    November 2015
!> \details 
!>  This routine computes new values of the thickness and tracers in each ice layer
!>  under horizontal advection using a first-order upwind scheme.
!>  Based on subroutine tend_layerThickness_fo_upwind by Matthew Hoffman
!
!-----------------------------------------------------------------------

    subroutine advect_thickness_tracers_upwind(&
         dt,                     &
         meshPool,               &
         normalVelocity,         &
         layerThicknessEdge,     &
         oldLayerThickness,      &
         oldTracers,             &
         newlayerThickness,      &
         newTracers,             &
         minOfMaxAllowableDt,    &
         err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(in) :: &
           dt                    !< Input: time step (s)

      type (mpas_pool_type), intent(in) :: &
           meshPool              !< Input: mesh information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
           normalVelocity,     & !< Input: normal velocity on cell edges
           layerThicknessEdge    !< Input: layer thickness on upstream edge

      real (kind=RKIND), dimension(:,:), intent(in) :: &
           oldLayerThickness     !< Input: old layer thickness

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
           oldTracers            !< Input: old tracer values

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
           newLayerThickness       !< Output: new layer thickness

      real (kind=RKIND), dimension(:,:,:), intent(out) :: &
           newTracers              !< Output: new tracer values

      real (kind=RKIND), intent(out) :: &
           minOfMaxAllowableDt     !< Output: maximum allowable dt based on CFL condition

      integer, intent(out) :: &
           err                     !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer, pointer :: &
           nCellsSolve,     & ! number of locally owned cells
           nVertLevels        ! number of vertical layers

      integer, dimension(:), pointer :: &
           nEdgesOnCell       ! number of edges on each cell

      integer, dimension(:,:), pointer :: &
           cellsOnEdge,     & ! index for 2 cells on each edge
           edgesOnCell,     & ! index for edges on each cell
           edgeSignOnCell     ! sign of edges on each cell

      real (kind=RKIND), dimension(:), pointer :: &
           dvEdge,          & ! distance between vertices on each edge
           dcEdge,          & ! distance between cell centers on each edge
           areaCell           ! area of each cell

      logical, pointer :: &
           config_print_thickness_advection_info

      integer, dimension(:), pointer :: &
           indexToCellID,   & ! global index for local cells
           indexToEdgeID      ! global index for local edges

      integer, pointer :: &
           config_stats_cell_ID      ! global index of diagnostic cell

      real (kind=RKIND) :: &
           invAreaCell,            & ! 1.0/areaCell
           maxAllowableDt,         & ! max allowable dt based on advective CFL condition
           layerNormalVelocity,    & ! normal velocity averaged from interfaces to layer midpoints
           thicknessFluxEdge,      & ! thickness flux on an edge
           thicknessTendency,      & ! net thickness tendency for a cell
           newThickness              ! new layer thickness

      integer :: iEdge, iCell, iCell1, iCell2, iEdgeOnCell, k

      integer :: nTracers

      real (kind=RKIND), dimension(:), allocatable :: &
           tracersEdge,             & ! upstream tracer values on an edge
           thicknessTracerTendency, & ! net thickness*tracer tendency for a cell
           newThicknessTracers        ! new values of thickness*tracer

      real (kind=RKIND), parameter :: bigNumber = 1.0e16_RKIND  ! This is ~300 million years in seconds, but it is small enough not to overflow

      real(kind=RKIND) :: velSign     ! = 1.0_RKIND or -1.0_RKIND depending on sign of velocity

      integer :: err_tmp

      ! Variables for optional conservation check
      ! Note: This conservation check workw on one process only.
      !       For multiprocess runs, the conservation check would require global sums
      !        and would need to be outside a block loop.
      logical, parameter :: checkConservation = .false.
      real (kind=RKIND) :: initVolumeSum, finalVolumeSum, difference
      real (kind=RKIND) :: initEnergySum, finalEnergySum
      real (kind=RKIND), parameter :: eps11 = 1.e-11_RKIND

      err = 0

      nTracers = size(oldTracers,1)

      ! allocate temporary arrays
      allocate(tracersEdge(nTracers))
      allocate(thicknessTracerTendency(nTracers))
      allocate(newThicknessTracers(nTracers))

      ! get dimensions
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      ! get mesh arrays
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
      call mpas_pool_get_array(meshPool, 'indexToEdgeID', indexToEdgeID)

      ! get config variables
      call mpas_pool_get_config(liConfigs, 'config_print_thickness_advection_info', config_print_thickness_advection_info)
      call mpas_pool_get_config(liConfigs, 'config_stats_cell_ID', config_stats_cell_ID)

      ! initialize output variables
      minOfMaxAllowableDt = bigNumber
      newLayerThickness(:,:) = 0.0_RKIND
      newTracers(:,:,:) = 0.0_RKIND

      if (checkConservation) then

         ! compute initial ice volume (= area*thickness) and energy (= area*thickness*tracer1)
         initVolumeSum = 0.0_RKIND
         initEnergySum = 0.0_RKIND
         do iCell = 1, nCellsSolve
            initVolumeSum = initVolumeSum + areaCell(iCell) * sum(oldLayerThickness(:,iCell))
            initEnergySum = initEnergySum + areaCell(iCell) * sum(oldLayerThickness(:,iCell)*oldTracers(1,:,iCell))
         enddo

      endif

      ! Note: This loop structure (nCells loop outside nEdgesOnCell loop) results in double calculation of fluxes
      !       across each edge. Upwind advection is so cheap, however, that the extra cost is minimal.
      ! TODO: Reorganize as a loop over edges, and then cells?

      ! loop over locally owned cells
      do iCell = 1, nCellsSolve

         invAreaCell = 1.0_RKIND / areaCell(iCell)

         if (indexToCellID(iCell) == config_stats_cell_ID .and. config_print_thickness_advection_info) then
            write(stderrUnit,*) 'In advect_thickness_tracer, iCell =', iCell
            write(stderrUnit,*)
            write(stderrUnit,*) 'k, iEdgeOnCell, layerNormalVelocity, layerThicknessEdge, thicknessFluxEdge, thicknessTendency' 
         endif

         ! loop over layers
         do k = 1, nVertLevels

            ! initialize the tendencies for this layer
            thicknessTendency = 0.0_RKIND
            thicknessTracerTendency(:) = 0.0_RKIND

            ! compute fluxes for each edge of the cell
            do iEdgeOnCell = 1, nEdgesOnCell(iCell)

               iEdge = edgesOnCell(iEdgeOnCell, iCell)
               iCell1 = cellsOnEdge(1,iEdge)
               iCell2 = cellsOnEdge(2,iEdge)

               ! Average normal velocities from layer interfaces to layer midpoints for advection
               ! TODO - This may make more sense to calculate as a 3d field somewhere else
               !        so it can also be used for tracer advection, output visualization, etc.
               layerNormalVelocity = 0.5_RKIND * (normalVelocity(k,iEdge) + normalVelocity(k+1,iEdge))

               ! Check for potential CFL violation
               if (abs(layerNormalVelocity) > 0.0_RKIND) then
                  maxAllowableDt = (0.5_RKIND * dcEdge(iEdge)) / abs(layerNormalVelocity)
               else
                  maxAllowableDt = bigNumber
               endif

               if ( maxAllowableDt < dt ) then
                  write(stderrUnit,*) 'CFL violation at level, edge ID (local, global):', k, iEdge, indexToEdgeID(iEdge)
                  err = err + 1
               endif
               minOfMaxAllowableDt = min(minOfMaxAllowableDt, maxAllowableDt)

               velSign = sign(1.0_RKIND, layerNormalVelocity)

               ! Assign the upwind tracer values to each edge
               ! Note: layerThicknessEdge is passed in and does not need to be computed here

               tracersEdge(:) = max(velSign * oldTracers(:,k,iCell1),   &
                                    velSign * (-1.0_RKIND) * oldTracers(:,k,iCell2))
               
               ! Increment the thickness and thickness*tracer tendencies

               thicknessFluxEdge = layerNormalVelocity * dvEdge(iEdge) * layerThicknessEdge(k,iEdge)
               thicknessTendency = thicknessTendency &
                                 + edgeSignOnCell(iEdgeOnCell,iCell) * thicknessFluxEdge * invAreaCell
               thicknessTracerTendency(:) = thicknessTracerTendency(:) &
                                 + edgeSignOnCell(iEdgeOnCell,iCell) * thicknessFluxEdge * invAreaCell * tracersEdge(:)

               !WHL - debug
               if (indexToCellID(iCell) == config_stats_cell_ID .and. k==1 .and. config_print_thickness_advection_info) then
                  write(stderrUnit,*) k, iEdgeOnCell, layerNormalVelocity, layerThicknessEdge(k,iEdge), thicknessFluxEdge, thicknessTendency 
               endif

            enddo    ! iEdgeOnCell

            ! Compute the new layer thickness and tracers

            newLayerThickness(k,iCell) = oldLayerThickness(k,iCell) + thicknessTendency * dt

            newThicknessTracers(:) = oldLayerThickness(k,iCell)*oldTracers(:,k,iCell) + thicknessTracerTendency(:) * dt

            if (newLayerThickness(k,iCell) > 0.0_RKIND) then
               newTracers(:,k,iCell) = newThicknessTracers(:) / newLayerThickness(k,iCell)
            else
               newTracers(:,k,iCell) = 0.0_RKIND
            endif

         enddo     ! k

      enddo   ! iCell


      if (checkConservation) then

         ! compute final ice volume (= area*thickness) and energy (= area*thickness*tracer1)
         finalVolumeSum = 0.0_RKIND
         finalEnergySum = 0.0_RKIND
         do iCell = 1, nCellsSolve
            finalVolumeSum = finalVolumeSum + areaCell(iCell) * sum(newLayerThickness(:,iCell))
            finalEnergySum = finalEnergySum + areaCell(iCell) * sum(newLayerThickness(:,iCell)*newTracers(1,:,iCell))
         enddo

         if (config_print_thickness_advection_info) then
            write(stderrUnit,*) 'init volume, final volume, difference:', initVolumeSum, finalVolumeSum, finalVolumeSum - initVolumeSum
            write(stderrUnit,*) 'init energy, final energy, difference:', initEnergySum, finalEnergySum, finalEnergySum - initEnergySum
         endif

         !TODO - Make these fatal errors instead of warnings?
         difference = abs(finalVolumeSum - initVolumeSum)
         if (difference/initVolumeSum > eps11) then
            write(stderrUnit,*) 'WARNING: upwind advection, mass conservation error'
            write(stderrUnit,*) 'init volume, final volume, difference:', initVolumeSum, finalVolumeSum, finalVolumeSum - initVolumeSum
         endif

         difference = abs(finalEnergySum - initEnergySum)
         if (difference/initEnergySum > eps11) then
            write(stderrUnit,*) 'WARNING: upwind advection, mass*tracer conservation error'
            write(stderrUnit,*) 'init energy, final energy, difference:', initEnergySum, finalEnergySum, finalEnergySum - initEnergySum
         endif

      endif

      ! clean up
      deallocate(tracersEdge)
      deallocate(thicknessTracerTendency)
      deallocate(newThicknessTracers)

      ! === error check
      if (err > 0) then
          write (stderrUnit,*) "An error has occurred in advect_thickness_tracers_upwind."
      endif

   !--------------------------------------------------------------------

    end subroutine advect_thickness_tracers_upwind

!-----------------------------------------------------------------------
!
!  subroutine li_vertical_remap
!
!> \brief   Vertical remapping of thickness and tracers
!> \author  Matt Hoffman
!> \date    October 2013; revised November 2015
!> \details
!>  This routine performs vertical remapping of thickness and tracers from one vertical
!>  coordinate system to another, as is required for our sigma coordinate system.
!>  The remapping is first-order accurate.
!>  This uses code from the CISM glissade_transport.F90 module written by Bill Lipscomb.
!>  I have altered the array structures to work with MPAS and refactored it.
!>  It now does all calculations column-wise, so it can be vectorized using
!>  OpenMP over either blocks or cells.
!>  I have dealt with 0 thickness in non-ice cells using eps and the cellMask,
!>  rather than using if/where-statements.
!
!-----------------------------------------------------------------------
   subroutine li_vertical_remap(thickness, cellMask, meshPool, layerThickness, tracers, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: LI mesh information

      real(kind=RKIND), dimension(:), intent(in) :: &
         thickness         !< Input: ice thickness

      integer, dimension(:), intent(in) :: &
         cellMask          !< Input: mask for cells (needed for determining presence/absence of ice)

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         layerThickness    !< Input/Output: thickness of layers (to be updated)

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tracers           !< Input/Output: tracer values (to be updated)

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

      ! pointers to mesh arrays
      real (kind=RKIND), dimension(:), pointer :: layerThicknessFractions, layerInterfaceSigma

      ! local arrays
      real (kind=RKIND), dimension(:), allocatable :: layerInterfaceSigma_Input
      real (kind=RKIND), dimension(:,:), allocatable :: hTsum

      ! counters, mesh variables, index variables
      integer, pointer :: nCells, nVertLevels
      integer :: nTracers, iCell, k, k1, k2, nt

      ! stuff for making calculations
      real(kind=RKIND) :: thisThk, zhi, zlo, hOverlap

      ! variables for optional conservation check
      ! set to 'true' for now, since it is inexpensive and might catch problems
      logical, parameter :: checkConservation = .true.
      real (kind=RKIND) :: initEnergySum, finalEnergySum, difference
      real (kind=RKIND), parameter :: eps11 = 1.e-11_RKIND
      real(kind=RKIND), parameter :: eps30 = 1.0d-30  ! Make a tiny value so we don't corrupt any cells with ice 

      !WHL - debug
      integer, dimension(:), pointer :: indexToCellID
      integer, pointer :: config_stats_cell_ID

      err = 0

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      nTracers = size(tracers, 1)

      call mpas_pool_get_array(meshPool, 'layerThicknessFractions', layerThicknessFractions)
      call mpas_pool_get_array(meshPool, 'layerInterfaceSigma', layerInterfaceSigma)

      !WHL - debug
      call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
      call mpas_pool_get_config(liConfigs, 'config_stats_cell_ID', config_stats_cell_ID)

      allocate(layerInterfaceSigma_Input(nVertLevels+1))
      allocate(hTsum(nTracers, nVertLevels))

      !TODO  We could avoid the eps30 by inserting 'if (thickness(iCell) > 0.0_RKIND'
      !      Not sure how important it is to avoid the 'if'

      ! loop over cells
      do iCell = 1, nCells

         if (checkConservation) then   ! compute sum of layerThickness*tracer1 in the column
            initEnergySum = sum(layerThickness(:,iCell)*tracers(1,:,iCell))
         endif

         thisThk = thickness(iCell)

         ! *** Calculate vertical sigma coordinates of each layer interface for the input non-sigma state
         !     (we already have the desired new sigma-based state as layerInterfaceSigma)

         layerInterfaceSigma_Input(1) = 0.0_RKIND
         do k = 2, nVertLevels
            layerInterfaceSigma_Input(k) = layerInterfaceSigma_Input(k-1) + layerThickness(k-1, iCell) / (thisThk + eps30)
         end do
         layerInterfaceSigma_Input(nVertLevels+1) = 1.0_RKIND
         
         !WHL - debug
         if (indexToCellID(iCell) == config_stats_cell_ID) then
!            write(stderrUnit,*) ' '
!            write(stderrUnit,*) 'thisThk:', thisThk
!            write(stderrUnit,*) 'Start layer 1:', layerThickness(1,iCell)
         endif

         ! *** Compute new layer thicknesses (layerInterfaceSigma coordinates)
         
         do k = 1, nVertLevels
            layerThickness(k,iCell) = layerThicknessFractions(k) * (thisThk + eps30)
         end do

         !WHL - debug
         if (indexToCellID(iCell) == config_stats_cell_ID) then
!            write(stderrUnit,*) 'New layer 1:', layerThickness(1,iCell)
         endif

         ! Check for conservation of mass.  Put any residual in the top layer.
         layerThickness(1,iCell) = layerThickness(1,iCell) + (thisThk - sum(layerThickness(:,iCell)) )

         !WHL - debug
         if (indexToCellID(iCell) == config_stats_cell_ID) then
!            write(stderrUnit,*) 'Final layer 1:', layerThickness(1,iCell)
         endif

         ! TODO This conservation check may make layerThicknesses inconsistent with the sigma levels (which are used below in tracer remapping).

         !-----------------------------------------------------------------
         ! Compute sum of h*T for each new layer (k2) by integrating
         ! over the regions of overlap with old layers (k1).
         ! Note: It might be worth trying a more efficient
         !       search algorithm if the number of layers is large.
         !       This algorithm scales as nlyr^2.
         !-----------------------------------------------------------------
         
         do k2 = 1, nVertLevels
            hTsum(:,k2) = 0.0_RKIND
            do k1 = 1, nVertLevels
               do nt = 1, nTracers
                  zhi = min (layerInterfaceSigma_Input(k1+1), layerInterfaceSigma(k2+1))
                  zlo = max (layerInterfaceSigma_Input(k1), layerInterfaceSigma(k2))
                  hOverlap = max (zhi-zlo, 0.0_RKIND) * thisThk
                  hTsum(nt,k2) = htsum(nt,k2) + tracers(nt,k1,iCell) * hOverlap
               enddo         ! nt
            enddo            ! k1
         enddo               ! k2
         
         !-----------------------------------------------------------------
         ! Compute tracer values in new layers
         !-----------------------------------------------------------------
         
         do k = 1, nVertLevels
            do nt = 1, nTracers
               tracers(nt,k,iCell) = hTsum(nt,k) / (layerThickness(k, iCell) + eps30)
               layerThickness(k,iCell) = layerThickness(k,iCell) * real(li_mask_is_ice_int(cellMask(iCell)), RKIND)
               ! Note: the last term is used to multiply by 0 where there is no ice in order to clear out garbage introduced by eps30 
            enddo         ! nt
         enddo            ! k
         
         if (checkConservation) then   ! compute new sum of layerThickness*tracer1 in the column
            finalEnergySum = sum(layerThickness(:,iCell)*tracers(1,:,iCell))
            difference = abs(finalEnergySum - initEnergySum)
            if (initEnergySum > eps11) then
               if (difference/initEnergySum > eps11) then
                  write(stderrUnit,*) 'WARNING: vertical_remap, mass*tracer conservation error, iCell =', iCell
                  write(stderrUnit,*) 'init energy, final energy, difference:', initEnergySum, finalEnergySum, finalEnergySum - initEnergySum
               endif
            endif
         endif
         
      enddo ! nCellsSolve

      ! clean up
      deallocate(layerInterfaceSigma_Input)
      deallocate(hTsum)

   end subroutine li_vertical_remap

!***********************************************************************

 end module li_advection

