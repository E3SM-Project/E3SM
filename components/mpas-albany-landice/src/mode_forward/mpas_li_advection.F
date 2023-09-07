! Copyright (c) 2013-2018,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_advection
!
!> \brief MPAS land ice advection solver
!> \author William Lipscomb, Trevor Hillebrand
!> \date   December 2015, updated Sept 2023 to include FCT
!> \details
!>  This module contains the routines for advecting ice thickness
!>  and tracers for land ice.
!-----------------------------------------------------------------------

module li_advection

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_log

   use li_setup
   use li_mask
   use li_constants
   use li_config
   use li_mesh

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
   public :: li_advection_thickness_tracers, &
             li_layer_normal_velocity, &
             li_update_geometry

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

 contains


!***********************************************************************
!
!  subroutine li_advection_thickness_tracers
!
!> \brief   Advection for ice thickness and tracers
!> \author  William Lipscomb, Trevor Hillebrand
!> \date    November 2015, updated Sept 2023 to include FCT
!> \details
!>  This routine (1) computes new values of ice thickness and tracers under
!>  horizontal advection, (2) applies surface and basal mass balance
!>  terms, and (3) remaps tracers back onto the standard vertical grid.
!>  Based on subroutine li_tendency_thickness by Matthew Hoffman
!>  Note: The thickness and tracer fields must have had halo updates
!>        before calling this subroutine.
!-----------------------------------------------------------------------

   subroutine li_advection_thickness_tracers(&
        dt,                     &
        meshPool,               &
        velocityPool,           &
        geometryPool,           &
        thermalPool,            &
        scratchPool,            &
        err,                    &
        advectTracersIn)

     use li_setup, only: li_calculate_layerThickness
     use li_tracer_advection_fct, only: li_tracer_advection_fct_tend
     use li_tracer_advection_fct_shared
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(in) :: &
           dt                     !< Input: time step (s)

      type (mpas_pool_type), intent(in) :: &
           meshPool               !< Input: mesh information

      logical, intent(in), optional :: &
           advectTracersIn        !< Input: if true, then advect tracers as well as thickness
                                  !TODO: Default is false; change to true?

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

      integer, intent(out) :: &
           err                    !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer, pointer :: nVertLevels

      real (kind=RKIND), dimension(:), pointer :: &
           thickness,             & ! ice thickness (updated in this subroutine)
           bedTopography,         & ! bed topography
           sfcMassBal,            & ! surface mass balance (potential forcing)
           sfcMassBalApplied,     & ! surface mass balance (actually applied)
           groundedSfcMassBalApplied,     & ! surface mass balance on grounded locations (actually applied)
           basalMassBal,          & ! basal mass balance
           groundedBasalMassBal,  & ! basal mass balance for grounded ice
           floatingBasalMassBal,  & ! basal mass balance for floating ice
           basalMassBalApplied,          & ! basal mass balance limited by ice thickness
           groundedBasalMassBalApplied,  & ! basal mass balance for grounded ice limited by ice thickness
           floatingBasalMassBalApplied,  & ! basal mass balance for floating ice limited by ice thickness
           dynamicThickening,     & ! dynamic thickening rate
           groundedToFloatingThickness, & ! thickness changing from grounded to floating or vice versa
           fluxAcrossGroundingLine, &     ! magnitude of flux across GL
           fluxAcrossGroundingLineOnCells ! magnitude of flux across GL, calculated on cells

      real (kind=RKIND), dimension(:,:), pointer :: &
           temperature,           & ! interior ice temperature
           waterFrac,             & ! interior water fraction
           drainedInternalMeltRate,      & ! excess internal melt drained to the bed
           enthalpy                 ! interior ice enthalpy

      real (kind=RKIND), dimension(:,:), pointer :: &
           layerNormalVelocity,   & ! normal component of velocity on cell edges at layer midpoints
           layerThickness,        & ! thickness of each layer
           layerThicknessOld,     & ! old layer thickness
           layerThicknessEdge,    & ! layer thickness on upstream edge of cell
           normalVelocity           ! horizontal velocity on interfaces

      real (kind=RKIND), dimension(:,:,:), pointer :: &
           advectedTracers,       & ! values of advected tracers
           advectedTracersOld       ! old values of advected tracers

      real (kind=RKIND), dimension(:,:), pointer :: &
           surfaceTracers,        & ! tracer values for new ice at upper surface
           basalTracers             ! tracer values for new ice at lower surface

      type (field3DReal), pointer :: &
           advectedTracersField,  & ! scratch field containing values of advected tracers
           advectedTracersOldField  ! scratch field containing old values of advected tracers

      type (field2DReal), pointer :: &
           layerThicknessOldField   ! scratch field containing old values of layer thickness

      type (field2DReal), pointer :: &
           surfaceTracersField,   & ! scratch field containing values of surface tracers
           basalTracersField        ! scratch field containing values of basal tracers

      type (field1DInteger), pointer :: &
           cellMaskTemporaryField, & ! scratch field containing old values of cellMask
           thermalCellMaskField

      ! Allocatable arrays need for flux-corrected transport advection
      real (kind=RKIND), dimension(:,:,:), allocatable :: tend
      real (kind=RKIND), dimension(:,:,:), allocatable :: activeTracerHorizontalAdvectionEdgeFlux

      integer, dimension(:), pointer :: &
           cellMask,              & ! integer bitmask for cells
           edgeMask,              & ! integer bitmask for edges
           thermalCellMask

      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND), dimension(:), pointer :: dvEdge

      character (len=StrKIND), pointer :: &
           config_thickness_advection   ! method for advecting thickness and tracers

      logical, pointer :: &
           config_print_thickness_advection_info, &  !TODO - change to config_print_advection_info?
           config_print_thermal_info, &
           config_thermal_calculate_bmb

      real (kind=RKIND), pointer :: &
           config_ice_density,    & ! rhoi
           config_ocean_density,  & ! rhoo
           config_sea_level,      & ! sea level relative to z = 0
           config_thermal_thickness, &
           config_advection_coef_3rd_order

      integer, pointer :: config_num_halos

      logical :: advectTracers     ! if true, then advect tracers as well as thickness

      integer :: iCell1, iCell2, theGroundedCell
      integer, parameter :: horizAdvOrder = 4
      integer, parameter :: vertAdvOrder = 1

      real(kind=RKIND) :: GLfluxSign, thicknessFluxEdge

      integer :: err_tmp
      integer :: nTracers
      integer, pointer :: maxEdges2

      !WHL - debug
      integer, dimension(:), pointer :: indexToCellID, indexToEdgeID
      integer, pointer :: nCells, nEdges
      integer, pointer :: config_stats_cell_ID
      integer :: iEdge, iEdgeOnCell
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell

      !WHL - debug - for test-case diagnostics
      logical, parameter :: dome_test = .false.
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
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
      call mpas_pool_get_array(geometryPool, 'sfcMassBal', sfcMassBal)
      call mpas_pool_get_array(geometryPool, 'sfcMassBalApplied', sfcMassBalApplied)
      call mpas_pool_get_array(geometryPool, 'groundedSfcMassBalApplied', groundedSfcMassBalApplied)
      call mpas_pool_get_array(geometryPool, 'basalMassBal', basalMassBal)
      call mpas_pool_get_array(geometryPool, 'groundedBasalMassBal', groundedBasalMassBal)
      call mpas_pool_get_array(geometryPool, 'floatingBasalMassBal', floatingBasalMassBal)
      call mpas_pool_get_array(geometryPool, 'basalMassBalApplied', basalMassBalApplied)
      call mpas_pool_get_array(geometryPool, 'groundedBasalMassBalApplied', groundedBasalMassBalApplied)
      call mpas_pool_get_array(geometryPool, 'floatingBasalMassBalApplied', floatingBasalMassBalApplied)
      call mpas_pool_get_array(geometryPool, 'layerThickness', layerThickness)
      call mpas_pool_get_array(geometryPool, 'layerThicknessEdge', layerThicknessEdge)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_array(geometryPool, 'edgeMask', edgeMask)
      call mpas_pool_get_array(geometryPool, 'dynamicThickening', dynamicThickening)
      call mpas_pool_get_array(geometryPool, 'groundedToFloatingThickness', groundedToFloatingThickness)

      ! get arrays from the velocity pool
      call mpas_pool_get_array(velocityPool, 'layerNormalVelocity', layerNormalVelocity)
      call mpas_pool_get_array(velocityPool, 'normalVelocity', normalVelocity)
      call mpas_pool_get_array(velocityPool, 'fluxAcrossGroundingLine', fluxAcrossGroundingLine)
      call mpas_pool_get_array(velocityPool, 'fluxAcrossGroundingLineOnCells', fluxAcrossGroundingLineOnCells)

      ! get arrays from the thermal pool
      call mpas_pool_get_array(thermalPool, 'temperature', temperature)
      call mpas_pool_get_array(thermalPool, 'waterFrac', waterFrac)
      call mpas_pool_get_array(thermalPool, 'enthalpy', enthalpy)
      call mpas_pool_get_array(thermalPool, 'drainedInternalMeltRate', drainedInternalMeltRate)

      ! get config variables
      call mpas_pool_get_config(liConfigs, 'config_thickness_advection', config_thickness_advection)
      call mpas_pool_get_config(liConfigs, 'config_print_thickness_advection_info', config_print_thickness_advection_info)
      call mpas_pool_get_config(liConfigs, 'config_print_thermal_info', config_print_thermal_info)
      call mpas_pool_get_config(liConfigs, 'config_thermal_calculate_bmb', config_thermal_calculate_bmb)
      call mpas_pool_get_config(liConfigs, 'config_thermal_thickness', config_thermal_thickness)
      call mpas_pool_get_config(liConfigs, 'config_ice_density', config_ice_density)
      call mpas_pool_get_config(liConfigs, 'config_ocean_density', config_ocean_density)
      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
      call mpas_pool_get_config(liConfigs, 'config_num_halos', config_num_halos)

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'maxEdges2', maxEdges2)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
      call mpas_pool_get_array(meshPool, 'indexToEdgeID', indexToEdgeID)
      call mpas_pool_get_config(liConfigs, 'config_stats_cell_ID', config_stats_cell_ID)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)

      ! get fields from the scratch pool
      ! Note: The advectedTracers field is a scratch field containing only the tracers that need to be advected.
      !       For example, it contains temperature or enthalpy (depending on config_thermal_solver) but not both.

      call mpas_pool_get_field(scratchPool, 'workTracerLevelCell', advectedTracersField)
      call mpas_allocate_scratch_field(advectedTracersField, .true.)
      advectedTracers => advectedTracersField % array

      call mpas_pool_get_field(scratchPool, 'workLevelCell', layerThicknessOldField)
      call mpas_allocate_scratch_field(layerThicknessOldField, .true.)
      layerThicknessOld => layerThicknessOldField % array

      call mpas_pool_get_field(scratchPool, 'workTracerLevelCell2', advectedTracersOldField)
      call mpas_allocate_scratch_field(advectedTracersOldField, .true.)
      advectedTracersOld => advectedTracersOldField % array

      call mpas_pool_get_field(scratchPool, 'workTracerCell', surfaceTracersField)
      call mpas_allocate_scratch_field(surfaceTracersField, .true.)
      surfaceTracers => surfaceTracersField % array

      call mpas_pool_get_field(scratchPool, 'workTracerCell2', basalTracersField)
      call mpas_allocate_scratch_field(basalTracersField, .true.)
      basalTracers => basalTracersField % array

      call mpas_pool_get_field(geometryPool, 'cellMaskTemporary', cellMaskTemporaryField)
      call mpas_allocate_scratch_field(cellMaskTemporaryField, .true.)

      call mpas_pool_get_field(scratchPool, 'iceCellMask', thermalCellMaskField)
      call mpas_allocate_scratch_field(thermalCellMaskField, .true.)
      thermalCellMask => thermalCellMaskField % array

      ! define max number of advection cells as max number of edges*2
      nAdvCellsMax = maxEdges2

      ! given the old thickness, compute the thickness in each layer
      call li_calculate_layerThickness(meshPool, thickness, layerThickness)

      ! update masks
      call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
      err = ior(err, err_tmp)

      ! save old copycellMask for determining cells changing from grounded to floating and vice versa
      cellMaskTemporaryField % array(:) = cellMask(:)

      !-----------------------------------------------------------------
      ! Horizontal transport of thickness and tracers
      !-----------------------------------------------------------------

      if ( (trim(config_thickness_advection) .eq. 'fo') .or. &
          (trim(config_thickness_advection) .eq. 'fct') ) then

         if (config_print_thickness_advection_info) then
            call mpas_log_write('Using first-order upwind for thickness advection')
            call mpas_log_write('advectTracers = $l', logicArgs=(/advectTracers/))
            if (dome_test) then
               call mpas_log_write(' ')
               call mpas_log_write('Dome test diagnostics:')
            endif
         endif


         if ( (advectTracers) .or. (trim(config_thickness_advection) .eq. 'fct') ) then

            ! Copy the old tracer values into the advectedTracers array.
            ! This requires setting a tracer pointer to either temperature or enthalpy,
            !  depending on the value of config_thermal_solver.
            ! Currently (as of Nov. 2015), this is the only tracer transported, but others
            !  may be added later.

            call tracer_setup(&
                 meshPool,          &
                 geometryPool,      &
                 thermalPool,       &
                 advectedTracers,   &
                 surfaceTracers,    &
                 basalTracers,      &
                 nTracers)

         else   ! no tracer advection; pass a tracer array full of zeroes

            advectedTracers(:,:,:) = 0.0_RKIND

         endif

         if ((trim(config_thickness_advection) .eq. 'fct') .or. &
             (trim(config_tracer_advection) .eq. 'fct' ) ) then
            allocate(tend(nTracers,nVertLevels,nCells+1))
            tend(:,:,:) = 0.0_RKIND
            allocate(activeTracerHorizontalAdvectionEdgeFlux(nTracers,nVertLevels+1,nEdges+1))
            activeTracerHorizontalAdvectionEdgeFlux(:,:,:) = 0.0_RKIND
         endif

         ! Transport thickness and tracers
         ! Note:  For the enthalpy scheme, temperature and waterFrac are the primary prognostic
         !        variables to be updated, but enthalpy is the advected tracer (for reasons of
         !        energy conservation).

         if (config_print_thickness_advection_info) then

            do iCell = 1, nCells
               if (indexToCellID(iCell) == config_stats_cell_ID) then
                  call mpas_log_write('Before thickness advection, iCell (local=$i, global=$i), thickness=$r, layer 1 tracer:=$r',&
                       intArgs=(/iCell, indexToCellID(iCell)/), realArgs=(/thickness(iCell), advectedTracers(1,1,iCell)/))
                  do iEdgeOnCell = 1, nEdgesOnCell(iCell)
                     iEdge = edgesOnCell(iEdgeOnCell,iCell)
                     call mpas_log_write('iEdge (local=$i, global=$i), &
                             layerNormalVelocity=$r:', &
                             intArgs=(/iEdge, indexToEdgeID(iEdge)/), &
                             realArgs=(/layerNormalVelocity(1,iEdge)/))
                  enddo
               endif
            enddo
         endif

         ! copy old (input) values of layer thickness and tracers to scratch arrays
         layerThicknessOld(:,:) = layerThickness(:,:)
         advectedTracersOld(:,:,:) = advectedTracers(:,:,:)

         ! compute new values of layer thickness and tracers
         if ((trim(config_thickness_advection) .eq. 'fo') .and. &
             ( (trim(config_tracer_advection) .eq. 'fo') .or. &
               (trim(config_tracer_advection) .eq. 'none') ) ) then
             call advect_thickness_tracers_upwind(&
                  dt,                      &
                  meshPool,                &
                  layerNormalVelocity,     &
                  layerThicknessEdge,      &
                  layerThicknessOld,       &
                  advectedTracersOld,      &
                  cellMask,                &
                  edgeMask,                &
                  layerThickness,          &
                  advectedTracers,         &
                  err)
         elseif ((trim(config_thickness_advection) .eq. 'fo') .and. &
                  trim(config_tracer_advection) .eq. 'fct') then
             call advect_thickness_tracers_upwind(&
                  dt,                      &
                  meshPool,                &
                  layerNormalVelocity,     &
                  layerThicknessEdge,      &
                  layerThicknessOld,       &
                  advectedTracersOld,      &
                  cellMask,                &
                  edgeMask,                &
                  layerThickness,          &
                  advectedTracers,         &
                  err)
             ! Reset tracers after fo advection for fct advection
             advectedTracers(:,:,:) = advectedTracersOld(:,:,:)

             ! Pass FO upwind normalThicknessFlux (layerThicknessEdge * layerNormalVelocity)
             ! to fct tracer routine
             call li_tracer_advection_fct_tend(&
                  tend, advectedTracers, layerThicknessOld, &
                  layerThicknessEdge * layerNormalVelocity, 0 * normalVelocity, dt,    &
                  nTracers, activeTracerHorizontalAdvectionEdgeFlux, computeBudgets=.false.)!{{{
         elseif (trim(config_thickness_advection) .eq. 'fct') then
             ! Call fct routine for thickness first, and use activeTracerHorizontalAdvectionEdgeFlux
             ! returned by that call as normalThicknessFlux for call to tracer fct
             call li_tracer_advection_fct_tend(&
                  tend(nTracers:,:,:), advectedTracers(nTracers:,:,:), layerThicknessOld * 0.0_RKIND + 1.0_RKIND, &
                  layerNormalVelocity, 0.0_RKIND * normalVelocity, dt,    &
                  1, activeTracerHorizontalAdvectionEdgeFlux(nTracers:,:,:), computeBudgets=.false.)
             ! layerThickness is last tracer. However, for some reason
             ! this: layerThickness(:,:) = advectedTracers(nTracers,:,:) does not conserve mass!
             ! This does conserve mass:
             layerThickness(:,:) = layerThickness(:,:) + tend(nTracers,:,:) * dt

             if (trim(config_tracer_advection) .eq. 'fct') then
                ! Call fct for tracers, using activeTracerHorizontalAdvectionEdgeFlux
                ! from fct thickness advection as normalThicknessFlux
                call li_tracer_advection_fct_tend(&
                     tend(1:nTracers-1,:,:), advectedTracers(1:nTracers-1,:,:), layerThicknessOld, &
                     activeTracerHorizontalAdvectionEdgeFlux(nTracers,:,:), 0.0_RKIND * normalVelocity, dt,    &
                     nTracers-1, computeBudgets=.false.)
             elseif (trim(config_tracer_advection) .eq. 'none') then
                ! do nothing
             else
                err_tmp = 1
                call mpas_log_write(trim(config_tracer_advection) // &
                                    ' tracer advection is not currently supported with fct thickness advection.', MPAS_LOG_ERR)
             endif
         else
            err_tmp = 1
            call mpas_log_write("config_thickness_advection = " // trim(config_thickness_advection) // &
                                ", config_tracer_advection = " // trim(config_tracer_advection) // &
                                " is not a supported combination.", MPAS_LOG_ERR)
         endif

         if (config_print_thickness_advection_info) then
            do iCell = 1, nCells
               if (indexToCellID(iCell) == config_stats_cell_ID) then
                  call mpas_log_write(' ')
                  call mpas_log_write('After thickness advection, iCell (local, global), new layer thickness, tracer 1:')
                  do k = 1, nVertLevels
                     call mpas_log_write("$i $i $r $r", intArgs=(/iCell, indexToCellID(iCell)/), &
                         realArgs=(/layerThickness(k,iCell), advectedTracers(1,k,iCell)/))
                  enddo
               endif
            enddo
         endif


         ! Calculate dynamicThickening (layerThickness is updated by advection at this point, while thickness is still old)
         dynamicThickening = (sum(layerThickness, 1) - thickness) / dt * scyr  ! units of m/yr


         ! Update the thickness and cellMask before applying the mass balance.
         ! The update is needed because the SMB and BMB depend on whether ice is present.
         thickness = sum(layerThickness, 1)
         call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)



         !-----------------------------------------------------------------
         ! Add the surface and basal mass balance to the layer thickness
         !-----------------------------------------------------------------

         ! Combine various basal mass balance fields based on mask.
         ! Grounded and floating basal mass balance should come from the thermal solver.

         ! TODO: more complicated treatment at GL?

         do iCell = 1, nCells
            if (li_mask_is_grounded_ice(cellMask(iCell))) then
               basalMassBal(iCell) = groundedBasalMassBal(iCell)
            elseif (li_mask_is_floating_ice(cellMask(iCell))) then
               ! Currently, floating and grounded ice are mutually exclusive.
               ! This could change if the GL is parameterized, in which case this logic may need adjustment.
               basalMassBal(iCell) = floatingBasalMassBal(iCell)
            elseif ( .not. (li_mask_is_ice(cellMask(iCell)))) then
               ! We don't allow a positive basal mass balance where ice is not already present.
               basalMassBal(iCell) = 0.0_RKIND
            endif
         enddo

         ! It is possible that excess internal melting was computed and assigned
         ! to the drainedInternalMeltRate array in mpas_li_thermal.F.  If so, then add it to basalMassBal.
         ! floatingBasalMassBal should never be altered because it is an input variable.
         ! Note: Subroutine basal_melt_floating_ice should be called earlier in the time step, before adding this term.
                  ! calculate a mask to identify ice that is thick enough to be thermally active
         if (config_thermal_calculate_bmb) then
            do iCell = 1, nCells
               if (thickness(iCell) > config_thermal_thickness) then
                  thermalCellMask(iCell) = 1
                  do k = 1, nVertLevels
                     basalMassBal(iCell) = basalMassBal(iCell) - drainedInternalMeltRate(k, iCell)
                  enddo
               else
                  thermalCellMask(iCell) = 0
               endif
            enddo
         endif

         call apply_mass_balance(&
              dt,                  &
              config_ice_density,  &
              cellMask,            &
              bedTopography,       &
              sfcMassBal,          &
              sfcMassBalApplied,   &
              groundedSfcMassBalApplied,   &
              basalMassBal,        &
              basalMassBalApplied, &
              groundedBasalMassBalApplied, &
              floatingBasalMassBalApplied, &
              surfaceTracers,      &
              basalTracers,        &
              layerThickness,      &
              advectedTracers)

         ! Update the thickness and cellMask before the vertical remap calculation.
         ! Note: The thickness field (rather than layerThickness) is the primary prognostic.
         !       It should be returned from each subroutine that alters the thickness, or
         !       updated immediately afterward.

         thickness = sum(LayerThickness, 1)

         !WHL - debug - Get mass balance for test cell
         if (config_print_thickness_advection_info) then
            do iCell = 1, nCells
               if (indexToCellID(iCell) == config_stats_cell_ID) then
                  call mpas_log_write(' ')
                  call mpas_log_write('After apply_mass_balance, indexToCellID=$i, thickness=$r', intArgs=(/indexToCellID(iCell)/), &
                          realArgs=(/thickness(iCell)/) )
                  call mpas_log_write('cellMask=$i, is ice=$l, is grounded=$l, is floating=$l', &
                       intArgs=(/cellMask(iCell)/), logicArgs=(/li_mask_is_ice(cellMask(iCell)), &
                       li_mask_is_grounded_ice(cellMask(iCell)), li_mask_is_floating_ice(cellMask(iCell)) /) )
                  call mpas_log_write('basalMassBalApplied=$r, grounded=$r, floating=$r', realArgs=(/ basalMassBalApplied(iCell)*31536000./917.,&
                          groundedBasalMassBalApplied(iCell)*31536000./917., floatingBasalMassBalApplied(iCell)*31536000./917. /) )
               endif
            enddo
         endif

         call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
         err = ior(err, err_tmp)

         ! Calculate volume converted from grounded to floating
         ! This needs to be determined after SMB/BMB are applied because those can change floating/grounded state
         call grounded_to_floating(cellMaskTemporaryField % array, cellMask, thickness, groundedToFloatingThickness, nCells)

         ! Calculate flux across grounding line
         ! Do this after new thickness & mask have been calculated, including SMB/BMB
         fluxAcrossGroundingLine(:) = 0.0_RKIND
         fluxAcrossGroundingLineOnCells(:) = 0.0_RKIND
         do iEdge = 1, nEdges
            if (li_mask_is_grounding_line(edgeMask(iEdge))) then
               iCell1 = cellsOnEdge(1,iEdge)
               iCell2 = cellsOnEdge(2,iEdge)
               if (li_mask_is_grounded_ice(cellMask(iCell1))) then
                  GLfluxSign = 1.0_RKIND ! edge sign convention is positive from iCell1 to iCell2 on an edge
                  theGroundedCell = iCell1
               else
                  GLfluxSign = -1.0_RKIND
                  theGroundedCell = iCell2
               endif
               do k = 1, nVertLevels
                  thicknessFluxEdge = layerNormalVelocity(k,iEdge) * dvEdge(iEdge) * layerThicknessEdge(k,iEdge)
                  fluxAcrossGroundingLine(iEdge) = fluxAcrossGroundingLine(iEdge) + GLfluxSign * thicknessFluxEdge / dvEdge(iEdge)
               enddo
               ! assign to grounded cell in fluxAcrossGroundingLineOnCells
               if (thickness(theGroundedCell) <= 0.0_RKIND) then
                  ! This should never be the case, but checking to avoid possible divide by zero
                  call mpas_log_write("thickness at a grounding line is unexepectedly <=0", MPAS_LOG_ERR)
                  err = ior(err, 1)
                  return
               endif
               fluxAcrossGroundingLineOnCells(theGroundedCell) = fluxAcrossGroundingLineOnCells(theGroundedCell) + &
                       fluxAcrossGroundingLine(iEdge) / thickness(theGroundedCell) * config_ice_density  ! adjust to correct units
            endif
         enddo ! edges


         ! Remap tracers to the standard vertical sigma coordinate
         ! Note: If tracers are not being advected, then this subroutine simply restores the
         !       layer thickness to sigma coordinate values.

         call vertical_remap(thickness, cellMask, meshPool, layerThickness, advectedTracers, err_tmp)
         err = ior(err, err_tmp)

         if (config_print_thickness_advection_info) then
            do iCell = 1, nCells
               if (indexToCellID(iCell) == config_stats_cell_ID) then
                  call mpas_log_write(' ')
                  call mpas_log_write('After vertical remap, iCell, new layer thickness, tracer 1:')
                  do k = 1, nVertLevels
                     call mpas_log_write("$i $r $r", intArgs=(/iCell/), realArgs=(/layerThickness(k,iCell), &
                             advectedTracers(1,k,iCell)/))
                  enddo
               endif
            enddo
         endif

         if (advectTracers) then

            ! Copy the advectedTracersNew values into the thermal tracer arrays
            ! (temperature, waterFrac, enthalpy)

            call tracer_finish(&
                 meshPool,              &
                 geometryPool,          &
                 thermalPool,           &
                 advectedTracers)
         endif

      elseif (trim(config_thickness_advection) .eq. 'none') then

         ! Do nothing

      else

         call mpas_log_write(trim(config_thickness_advection) // ' is not a valid option for thickness/tracer advection.', &
                 MPAS_LOG_ERR)
         err_tmp = 1

      end if ! config_thickness_advection

      err = ior(err,err_tmp)

      ! Deallocate arrays for fct
      if ( (trim(config_thickness_advection) .eq. 'fct') .or. &
           (trim(config_tracer_advection) .eq. 'fct') ) then
         deallocate( nAdvCellsForEdge, &
                     advCellsForEdge, &
                     advCoefs, &
                     advCoefs3rd, &
                     advMaskHighOrder, &
                     advMask2ndOrder, &
                     tend, &
                     activeTracerHorizontalAdvectionEdgeFlux)
      endif

      ! clean up
      call mpas_deallocate_scratch_field(advectedTracersField, .true.)
      call mpas_deallocate_scratch_field(advectedTracersOldField, .true.)
      call mpas_deallocate_scratch_field(layerThicknessOldField, .true.)
      call mpas_deallocate_scratch_field(basalTracersField, .true.)
      call mpas_deallocate_scratch_field(surfaceTracersField, .true.)
      call mpas_deallocate_scratch_field(cellMaskTemporaryField, .true.)
      call mpas_deallocate_scratch_field(thermalCellMaskField, .true.)

      ! === error check
      if (err > 0) then
         call mpas_log_write("An error has occurred in li_advection_thickness_tracers", MPAS_LOG_ERR)
      endif

   !--------------------------------------------------------------------
    end subroutine li_advection_thickness_tracers


!***********************************************************************
!
!  subroutine li_update_geometry
!
!> \brief   Calculates the lower and upper surface
!> \author  Matt Hoffman
!> \date    23 May 2017
!> \details
!>   This routine calculates lowerSurface and upperSurface from thickness
!>   cellMask, and bedTopography.  It is a standalone routine so it can be called
!>   consistently from multiple places as needed.
!-----------------------------------------------------------------------

   subroutine li_update_geometry(geometryPool)

     !-----------------------------------------------------------------
     ! input variables
     !-----------------------------------------------------------------

     !-----------------------------------------------------------------
     ! input/output variables
     !-----------------------------------------------------------------
     type (mpas_pool_type), intent(inout) :: geometryPool !< Input/Output: geometry object

     !-----------------------------------------------------------------
     ! output variables
     !-----------------------------------------------------------------


     !-----------------------------------------------------------------
     ! local variables
     !-----------------------------------------------------------------
     real (kind=RKIND), pointer :: config_sea_level, config_ice_density, config_ocean_density
     real (kind=RKIND), dimension(:), pointer :: thickness, upperSurface, &
        lowerSurface, bedTopography
     integer, dimension(:), pointer :: cellMask
     integer, pointer :: nCells
     integer :: iCell

     call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
     call mpas_pool_get_config(liConfigs, 'config_ice_density', config_ice_density)
     call mpas_pool_get_config(liConfigs, 'config_ocean_density', config_ocean_density)

     call mpas_pool_get_dimension(geometryPool, 'nCells', nCells)

     call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
     call mpas_pool_get_array(geometryPool, 'thickness', thickness)
     call mpas_pool_get_array(geometryPool, 'upperSurface', upperSurface)
     call mpas_pool_get_array(geometryPool, 'lowerSurface', lowerSurface)
     call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)


     ! Lower surface is based on floatation for floating ice.
     ! For grounded ice (and non-ice areas) it is the bed.
     where ( li_mask_is_floating_ice(cellMask) )
        lowerSurface = config_sea_level - thickness * (config_ice_density / config_ocean_density)
     elsewhere
        lowerSurface = bedTopography
     end where
     ! Make sure lowerSurface calculation is reasonable.  This check could be deleted once this has been throroughly tested.
     !do iCell = 1, nCells
     !   if (lowerSurface(iCell) < bedTopography(iCell)) then
     !      call mpas_log_write('lowerSurface less than bedTopography at cell:i $i', intArgs=(/iCell/))
     !      err = 1
     !   endif
     !end do

     ! Upper surface is the lower surface plus the thickness
     upperSurface(:) = lowerSurface(:) + thickness(:)

   !--------------------------------------------------------------------
   end subroutine li_update_geometry



!***********************************************************************
! Private subroutines:
!***********************************************************************

!***********************************************************************
!
!  subroutine apply_mass_balance
!
!> \brief   Apply surface and basal mass balance
!> \author  William Lipscomb
!> \date    February 2016
!> \details
!>  This routine applies the surface and basal mass balance, with a
!>  conservative treatment of tracers such as temperature.
!
!-----------------------------------------------------------------------

    subroutine apply_mass_balance(&
              dt,                  &
              rhoi,                &
              cellMask,            &
              bedTopography,       &
              sfcMassBal,          &
              sfcMassBalApplied,   &
              groundedSfcMassBalApplied,   &
              basalMassBal,        &
              basalMassBalApplied, &
              groundedBasalMassBalApplied, &
              floatingBasalMassBalApplied, &
              surfaceTracers,      &
              basalTracers,        &
              layerThickness,      &
              advectedTracers)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(in) ::  &
           dt,                 & !< Input: time step (s)
           rhoi                  !< Input: ice density (kg/m^3)

      integer, dimension(:), intent(in) :: cellMask !< Input: mask on cells

      real (kind=RKIND), dimension(:), intent(in) ::  &
           bedTopography,      & !< Input: bed elevation (m)
           sfcMassBal,         & !< Input: surface mass balance (kg/m^2/s)
           basalMassBal          !< Input: basal mass balance (kg/m^2/s)

      real(kind=RKIND), dimension(:,:), intent(in) :: &
           surfaceTracers,     & !< Input: tracer values of new ice at upper surface
           basalTracers          !< Input: tracer values of new ice at lower surface

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      real(kind=RKIND), dimension(:,:), intent(inout) :: &
           layerThickness        !< Input/Output: ice thickness in each layer (m)

      real(kind=RKIND), dimension(:,:,:), intent(inout) :: &
           advectedTracers       !< Input/Output: tracer values in each layer

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real(kind=RKIND), dimension(:), intent(out) :: &
           sfcMassBalApplied !< Output: surface mass balance actually applied on this time step (kg/m^2/s)

      real(kind=RKIND), dimension(:), intent(out) :: &
           groundedSfcMassBalApplied !< Output: surface mass balance actually applied to grounded ice on this time step (kg/m^2/s)

      real(kind=RKIND), dimension(:), intent(out) :: &
           basalMassBalApplied,         & !< Output: basal mass balance actually applied on this time step (kg/m^2/s)
           groundedBasalMassBalApplied, & !< Output: basal mass balance actually applied to grounded ice on this time step (kg/m^2/s)
           floatingBasalMassBalApplied    !< Output: basal mass balance actually applied to floating ice on this time step (kg/m^2/s)
      ! local variables

      real (kind=RKIND) ::  &
           sfcAccum, basalAccum,   & ! surface and basal accumulation (m), positive for ice gain
           sfcAblat, basalAblat      ! surface and basal ablation (m), positive for ice loss

      real (kind=RKIND), dimension(:), allocatable ::  &
           thckTracerProducts        ! thickness-tracer products

      integer :: nCells     ! number of cells
      integer :: nLayers    ! number of layers
      integer :: nTracers   ! number of tracers
      real (kind=RKIND), pointer :: config_sea_level ! sea level relative to z = 0
      logical, pointer :: config_zero_sfcMassBalApplied_over_bare_land

      integer :: iCell, iLayer, iTracer

      nLayers = size(layerThickness,1)
      nCells = size(layerThickness,2)
      nTracers = size(advectedTracers,1)

      allocate(thckTracerProducts(nTracers))

      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
      call mpas_pool_get_config(liConfigs, 'config_zero_sfcMassBalApplied_over_bare_land', &
                                config_zero_sfcMassBalApplied_over_bare_land)

      ! apply surface mass balance
      ! If positive, then add the SMB to the top layer, conserving mass*tracer products
      ! If negative, then melt from the top down until the melting term is used up or the ice is gone

      ! Initialize applied SMB and BMB fields
      sfcMassBalApplied(:) = sfcMassBal(:)
      basalMassBalApplied(:) = basalMassBal(:)

      do iCell = 1, nCells

         ! initialize accumulation/ablation terms

         sfcAccum = 0.0_RKIND
         sfcAblat = 0.0_RKIND
         basalAccum = 0.0_RKIND
         basalAblat = 0.0_RKIND

         ! apply surface mass balance

         ! Zero out any positive surface mass balance for ice-free ocean cells
         if (sfcMassBalApplied(iCell) > 0.0_RKIND .and. &
            bedTopography(iCell) < config_sea_level .and. .not.(li_mask_is_ice(cellMask(iCell)) ) ) then
            sfcMassBalApplied(iCell) = 0.0_RKIND
         end if

         ! Zero out positive surface mass balance for bare land cells
         if (config_zero_sfcMassBalApplied_over_bare_land) then
            if (sfcMassBalApplied(iCell) > 0.0_RKIND .and. &
               bedTopography(iCell) >= config_sea_level .and. .not.(li_mask_is_ice(cellMask(iCell)) ) ) then
               sfcMassBalApplied(iCell) = 0.0_RKIND
            endif
         endif

         ! surface accumulation
         if (sfcMassBalApplied(iCell) > 0.0_RKIND) then

            ! modify tracers conservatively in top layer

            sfcAccum = sfcMassBalApplied(iCell) * dt / rhoi

            ! compute mass-tracer products in top layer
            thckTracerProducts(:) = layerThickness(1,iCell)*advectedTracers(:,1,iCell)  &
                                  + sfcAccum * surfaceTracers(:,iCell)

            ! new thickness in top layer
            layerThickness(1,iCell) = layerThickness(1,iCell) + sfcAccum

            ! new tracers in top layer
            advectedTracers(:,1,iCell) = thckTracerProducts(:) / layerThickness(1,iCell)

         ! surface ablation
         elseif (sfcMassBalApplied(iCell) < 0.0_RKIND) then

            sfcAblat = -1.0_RKIND * sfcMassBalApplied(iCell) * dt /rhoi    ! positive for melting

            do iLayer = 1, nLayers
               if (sfcAblat > layerThickness(iLayer,iCell)) then   ! melt the entire layer
                  sfcAblat = sfcAblat - layerThickness(iLayer,iCell)
                  layerThickness(iLayer,iCell) = 0.0_RKIND
                  advectedTracers(:,iLayer,iCell) = 0.0_RKIND
               else   ! melt part of the layer
                  layerThickness(iLayer,iCell) = layerThickness(iLayer,iCell) - sfcAblat
                  sfcAblat = 0.0_RKIND
                  exit
               endif
            enddo

            ! Adjust applied SMB to only be the amount actually applied above (this will include ice-free cells)
            if (sfcAblat > 0.0_RKIND) then
               sfcMassBalApplied(iCell) = sfcMassBalApplied(iCell) + sfcAblat * rhoi / dt
               !TODO - If remaining sfcAblat > 0, then keep track of it to conserve energy (?)
            endif


         endif   ! sfcMassBal > 0

         ! apply basal mass balance

         if (basalMassBalApplied(iCell) > 0.0_RKIND) then

            ! basal freeze-on
            ! modify tracers conservatively in top layer

            basalAccum = basalMassBalApplied(iCell) * dt / rhoi

            ! compute mass-tracer products in bottom layer
            thckTracerProducts(:) = layerThickness(nLayers,iCell)*advectedTracers(:,nLayers,iCell)  &
                                  + basalAccum * basalTracers(:,iCell)

            ! new thickness in top layer
            layerThickness(nLayers,iCell) = layerThickness(nLayers,iCell) + basalAccum

            ! new tracers in top layer
            advectedTracers(:,nLayers,iCell) = thckTracerProducts(:) / layerThickness(nLayers,iCell)

         elseif (basalMassBalApplied(iCell) < 0.0_RKIND) then

            ! surface ablation from the bottom up

            basalAblat = -basalMassBalApplied(iCell) * dt /rhoi    ! positive for melting

            do iLayer = nLayers, 1, -1
               if (basalAblat > layerThickness(iLayer,iCell)) then   ! melt the entire layer
                  basalAblat = basalAblat - layerThickness(iLayer,iCell)
                  layerThickness(iLayer,iCell) = 0.0_RKIND
                  advectedTracers(:,iLayer,iCell) = 0.0_RKIND
               else   ! melt part of the layer
                  layerThickness(iLayer,iCell) = layerThickness(iLayer,iCell) - basalAblat
                  basalAblat = 0.0_RKIND
                  exit
               endif
            enddo

            if (basalAblat > 0.0_RKIND) then
               basalMassBalApplied(iCell) = basalMassBalApplied(iCell) + basalAblat * rhoi / dt
               !TODO - If remaining basalAblat > 0, then keep track of it to conserve energy
            endif

         endif   ! basalMassBalApplied > 0

      enddo   ! iCell

      ! Separate grounded and floating components, as necessary.
      where (li_mask_is_grounded_ice(cellMask) .or. bedTopography > config_sea_level)
            groundedSfcMassBalApplied = sfcMassBalApplied
      elsewhere
            groundedSfcMassBalApplied = 0.0_RKIND
      end where

      do iCell = 1, nCells
         if (li_mask_is_grounded_ice(cellMask(iCell))) then
            groundedBasalMassBalApplied(iCell) = basalMassBalApplied(iCell)
            floatingBasalMassBalApplied(iCell) = 0.0_RKIND
         elseif (li_mask_is_floating_ice(cellMask(iCell))) then
            floatingBasalMassBalApplied(iCell) = basalMassBalApplied(iCell)
            groundedBasalMassBalApplied(iCell) = 0.0_RKIND
         else
            groundedBasalMassBalApplied(iCell) = 0.0_RKIND
            floatingBasalMassBalApplied(iCell) = 0.0_RKIND
         endif
      enddo

      deallocate(thckTracerProducts)

    end subroutine apply_mass_balance

!***********************************************************************
!
!  subroutine tracer_setup
!
!> \brief   Assemble a 3D array for tracer advection
!> \author  William Lipscomb, Trevor Hillebrand
!> \date    November 2015, updated Sept 2023 for FCT
!> \details
!>  This routine assembles a 3D array for tracer advection.
!>  Each tracer in the array is transported conservatively, along
!>  with the layer thickness field.
!
!-----------------------------------------------------------------------

    subroutine tracer_setup(&
         meshPool,         &
         geometryPool,     &
         thermalPool,      &
         advectedTracers,  &
         surfaceTracers,   &
         basalTracers,     &
         nTracers)

      use li_thermal, only: li_temperature_to_enthalpy_kelvin
      use li_tracer_advection_fct_shared
      use li_tracer_advection_fct
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
           advectedTracers   ! tracers to be advected
                             ! dimension 1 = maxTracers, 2 = nVertLevels, 3 = nCells

      real (kind=RKIND), dimension(:,:), intent(out) :: &
           surfaceTracers, & ! tracer values for new ice at upper surface
           basalTracers      ! tracer values for new ice at lower surface
                             ! dimension 1 = maxTracers, 2 = nCells

      integer, intent(out)  :: nTracers
      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      character (len=StrKIND), pointer :: &
           config_thermal_solver, &
           config_thickness_advection

      logical, pointer :: &
           config_calculate_damage

      integer, pointer :: &
           nCells                     ! number of cells

      integer, pointer :: &
           nVertLevels      ! number of vertical layers

      real (kind=RKIND), dimension(:,:), pointer :: &
           temperature,             & ! interior ice temperature
           waterFrac,               & ! interior water fraction
           enthalpy,                & ! interior ice enthalpy
           layerThickness

      real (kind=RKIND), dimension(:), pointer :: &
           surfaceAirTemperature,   & ! surface air temperature
           basalTemperature           ! basal ice temperature

      real (kind=RKIND), dimension(:), pointer :: &
           layerCenterSigma           ! sigma coordinate at midpoint of each layer

      real (kind=RKIND), dimension(:), pointer :: &
           thickness                  ! ice thickness

      real (kind=RKIND), dimension(:), pointer :: &
           damage, & ! damage
           passiveTracer2d

      real (kind=RKIND), pointer ::  &
           config_ice_density         ! ice density

      integer :: iCell, iTracer, k, err, err1, err2

      err  = 0
      err1 = 0
      err2 = 0

      ! get dimensions
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      ! get arrays from mesh pool
      call mpas_pool_get_array(meshPool, 'layerCenterSigma', layerCenterSigma)

      ! get arrays from geometry pool
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(geometryPool, 'layerThickness', layerThickness)
      call mpas_pool_get_array(geometryPool, 'damage', damage)
      call mpas_pool_get_array(geometryPool, 'passiveTracer2d', passiveTracer2d)

      ! get arrays from thermal pool
      call mpas_pool_get_array(thermalPool, 'temperature', temperature)
      call mpas_pool_get_array(thermalPool, 'waterFrac', waterFrac)
      call mpas_pool_get_array(thermalPool, 'enthalpy', enthalpy)
      call mpas_pool_get_array(thermalPool, 'surfaceAirTemperature', surfaceAirTemperature)
      call mpas_pool_get_array(thermalPool, 'basalTemperature', basalTemperature)

      ! get config variables
      call mpas_pool_get_config(liConfigs, 'config_thermal_solver', config_thermal_solver)
      call mpas_pool_get_config(liConfigs, 'config_calculate_damage', config_calculate_damage)
      call mpas_pool_get_config(liConfigs, 'config_ice_density', config_ice_density)
      call mpas_pool_get_config(liConfigs, 'config_thickness_advection', config_thickness_advection)

      ! initialize
      advectedTracers(:,:,:) = 0.0_RKIND

      ! Notes:
      ! (1) For the enthalpy solver, it is necessary to transport enthalpy (rather than
      !     temperature and waterFrac separately) in order to conserve energy.
      ! (2) Temperature and waterFrac must be up to date (including halo cells)
      !     before calling this subroutine.
      ! (3) Surface and basal temperature (or enthalpy) are not transported, but their
      !     values are applied to new accumulation at either surface.
      ! (4) In most cases the surface temperature is equal to min(surfaceAirTemperature, 273.15),
      !     so we could set surfaceTracers to surfaceTemperature.  But if a positive SMB
      !     is applied to a previously ice-free cell, then surfaceTemperature has not yet been
      !     initialized to the correct value, so we explicitly use min(surfaceAirTemperature, 273.15).

      iTracer = 0    ! initialize the tracer index

      if (trim(config_thermal_solver) == 'enthalpy') then  ! advect enthalpy

         ! Rather than assume that enthalpy is up to date, recompute it from temperature and waterFrac

         do iCell = 1, nCells

            call li_temperature_to_enthalpy_kelvin(&
                 layerCenterSigma,      &
                 thickness(iCell),      &
                 temperature(:,iCell),  &
                 waterFrac(:,iCell),    &
                 enthalpy(:,iCell))

         enddo

         ! increment the tracer index
         iTracer = iTracer + 1

         ! copy enthalpy to the tracer array
         advectedTracers(iTracer,:,:) = enthalpy(:,:)

         ! set enthalpy of new ice at upper and lower surfaces
         ! convert temperature to enthalpy assuming waterFrac = 0 for new ice
         surfaceTracers(iTracer,:) = min(surfaceAirTemperature(:), kelvin_to_celsius) * config_ice_density*cp_ice
         basalTracers(iTracer,:) = basalTemperature(:) * config_ice_density*cp_ice

      elseif (trim(config_thermal_solver) == 'temperature') then  ! advect temperature

         ! increment the tracer index
         iTracer = iTracer + 1

         ! copy temperature to the tracer array
         advectedTracers(iTracer,:,:) = temperature(:,:)

         ! set enthalpy of new ice at upper and lower surfaces
         surfaceTracers(iTracer,:) = min(surfaceAirTemperature(:), kelvin_to_celsius)
         basalTracers(iTracer,:) = basalTemperature(:)

      endif

      if (config_calculate_damage) then

         ! increment the tracer index
         iTracer = iTracer + 1

         ! copy 2d damage to all vertical levels in the tracer array
         do k = 1, nVertLevels
            advectedTracers(iTracer,k,:) = damage(:)
         enddo

         ! set damage of new ice at upper and lower surfaces to 0
         ! TODO: Assess if this is the most appropriate option.  Other options would be setting damage equal to the existing value
         ! or the Nye zero stress value.  BISICLES implementation of the Bassis & Ma damage model sets damage associated with
         ! surface and basal accumulation as we are doing now.
         surfaceTracers(iTracer,:) = 0.0_RKIND
         basalTracers(iTracer,:) = 0.0_RKIND

      endif

      ! always advect passiveTracer2d
      iTracer = iTracer + 1
      ! copy 2d passiveTracer2d to all vertical levels in the tracer array
      do k = 1, nVertLevels
         advectedTracers(iTracer,k,:) = passiveTracer2d(:)
      enddo
      surfaceTracers(iTracer,:) = 0.0_RKIND
      basalTracers(iTracer,:) = 0.0_RKIND
 
      if (trim(config_thickness_advection) == 'fct') then
         ! increment the tracer index
         iTracer = iTracer + 1

         ! copy 2d thickness to all vertical levels in the tracer array
         advectedTracers(iTracer,:,:) = layerThickness(:,:)
         surfaceTracers(iTracer,:) = 0.0_RKIND
         basalTracers(iTracer,:) = 0.0_RKIND

      endif

      nTracers = iTracer
      ! Note: Other tracers (e.g., ice age) can be added here as needed.
      !       May need to increase maxTracers in the Registry.
      if ( (trim(config_tracer_advection) == 'fct') .or. &
           (trim(config_thickness_advection) == 'fct') ) then
         call li_tracer_advection_fct_shared_init(geometryPool, err1)
         call li_tracer_advection_fct_init(err2)

         if (err1 /= 0 .or. err2 /= 0) then
            err = 1
            call mpas_log_write(                                 &
               'Error encountered during fct tracer advection init', &
               MPAS_LOG_ERR, masterOnly=.true.)
         endif
      endif

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

      use li_thermal, only: li_enthalpy_to_temperature_kelvin

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

      logical, pointer :: &
           config_calculate_damage

      integer, pointer :: &
           nCellsSolve                ! number of locally owned cells

      real (kind=RKIND), dimension(:), pointer :: &
           layerCenterSigma           ! sigma coordinate at midpoint of each layer

      real (kind=RKIND), dimension(:), pointer :: &
           thickness                  ! ice thickness

      real (kind=RKIND), dimension(:), pointer :: &
           damage, & ! damage
           passiveTracer2d

      real (kind=RKIND), dimension(:,:), pointer :: &
           temperature,             & ! interior ice temperature
           waterFrac,               & ! water fraction
           enthalpy                   ! interior ice enthalpy

      real (kind=RKIND), dimension(:), pointer :: layerThicknessFractions

      integer :: iCell, iTracer

      ! get dimensions
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

      ! get arrays from mesh pool
      call mpas_pool_get_array(meshPool, 'layerCenterSigma', layerCenterSigma)
      call mpas_pool_get_array(meshPool, 'layerThicknessFractions', layerThicknessFractions)

      ! get arrays from geometry pool
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(geometryPool, 'damage', damage)
      call mpas_pool_get_array(geometryPool, 'passiveTracer2d', passiveTracer2d)

      ! get arrays from thermal pool
      call mpas_pool_get_array(thermalPool, 'temperature', temperature)
      call mpas_pool_get_array(thermalPool, 'waterFrac', waterFrac)
      call mpas_pool_get_array(thermalPool, 'enthalpy', enthalpy)

      ! get config variables
      call mpas_pool_get_config(liConfigs, 'config_thermal_solver', config_thermal_solver)
      call mpas_pool_get_config(liConfigs, 'config_calculate_damage', config_calculate_damage)

      iTracer = 0

      if (trim(config_thermal_solver) == 'enthalpy') then

         iTracer = iTracer + 1

         ! update the enthalpy
         enthalpy(:,:) = advectedTracers(iTracer,:,:)

         ! given the enthalpy, compute the temperature and water fraction

         do iCell = 1, nCellsSolve

            call li_enthalpy_to_temperature_kelvin(&
                 layerCenterSigma,                           &
                 thickness(iCell),                           &
                 enthalpy(:,iCell),     &
                 temperature(:,iCell),  &
                 waterFrac(:,iCell))

         enddo

      elseif (trim(config_thermal_solver) == 'temperature') then

         iTracer = iTracer + 1

         ! update the temperature
         temperature(:,:) = advectedTracers(iTracer,:,:)

      endif

      if (config_calculate_damage) then

         iTracer = iTracer + 1

         do iCell = 1, nCellsSolve
            damage(iCell) = sum(advectedTracers(iTracer,:,iCell) * layerThicknessFractions)  ! thickness-weighted average
         enddo
      endif

      ! Always advect passiveTracer2d
      iTracer = iTracer + 1
      do iCell = 1, nCellsSolve
         passiveTracer2d(iCell) = sum(advectedTracers(iTracer,:,iCell) * layerThicknessFractions)  ! thickness-weighted average
      enddo

      ! Note: Other tracers (e.g., ice age) can be added here as needed.
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
         layerNormalVelocity,    &
         layerThicknessEdge,     &
         layerThicknessOld,      &
         tracersOld,             &
         cellMask,               &
         edgeMask,               &
         layerThicknessNew,      &
         tracersNew,             &
         err,                    &
         advectTracersIn)

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
           layerNormalVelocity,& !< Input: normal velocity averaged from interfaces to layer midpoints
           layerThicknessEdge    !< Input: layer thickness on upstream edge

      logical, intent(in), optional :: &
           advectTracersIn       !< Input: if true, then advect tracers as well as thickness
                                 !  Default (if not passed in) is to advect tracers

      real (kind=RKIND), dimension(:,:), intent(in) :: &
           layerThicknessOld     !< Input: layer thickness

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
           tracersOld            !< Input: tracer values

      integer, dimension(:), intent(in) :: &
           cellMask              !< Input: mask on cells

      integer, dimension(:), intent(in) :: &
           edgeMask              !< Input: mask on edges

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
           layerThicknessNew     !< Output: layer thickness

      real (kind=RKIND), dimension(:,:,:), intent(out) :: &
           tracersNew            !< Output: tracer values

      integer, intent(out) :: &
           err                   !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer, pointer :: &
           nCellsSolve,     & ! number of locally owned cells
           nVertLevels,     & ! number of vertical layers
           nEdges             ! number of edges

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
           thicknessFluxEdge,      & ! thickness flux on an edge
           thicknessTendency,      & ! net thickness tendency for a cell
           newThickness              ! new layer thickness

      integer :: iEdge, iCell, iCell1, iCell2, iEdgeOnCell, k

      integer :: nTracers            ! number of tracers

      logical :: advectTracers       ! if true, then advect tracers as well as thickness

      real (kind=RKIND), dimension(:), allocatable :: &
           tracersEdge,             &! upstream tracer values on an edge
           thicknessTracerTendency, &! net thickness*tracer tendency for a cell
           newThicknessTracers       ! new values of thickness*tracer

      real (kind=RKIND), parameter :: bigNumber = 1.0e11_RKIND
      ! This is ~5000 years in seconds, but it is small enough not to overflow

      real(kind=RKIND) :: velSign     ! = 1.0_RKIND or -1.0_RKIND depending on sign of velocity
      real(kind=RKIND) :: GLfluxSign

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

      if (present(advectTracersIn)) then
         advectTracers = advectTracersIn
      else
         advectTracers = .true.
      endif

      nTracers = size(tracersOld,1)

      ! allocate temporary arrays
      allocate(tracersEdge(nTracers))
      allocate(thicknessTracerTendency(nTracers))
      allocate(newThicknessTracers(nTracers))

      ! get dimensions
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
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

      if (checkConservation) then

         ! compute initial ice volume (= area*thickness) and energy (= area*thickness*tracer1)
         initVolumeSum = 0.0_RKIND
         initEnergySum = 0.0_RKIND
         do iCell = 1, nCellsSolve
            initVolumeSum = initVolumeSum + areaCell(iCell) * sum(layerThicknessOld(:,iCell))
            initEnergySum = initEnergySum + areaCell(iCell) * sum(layerThicknessOld(:,iCell)*tracersOld(1,:,iCell))
         enddo

      endif


      ! Note: This loop structure (nCells loop outside nEdgesOnCell loop) results in double calculation of fluxes
      !       across each edge. But upwind advection is cheap, so the extra cost is minimal.

      ! loop over locally owned cells
      do iCell = 1, nCellsSolve

         invAreaCell = 1.0_RKIND / areaCell(iCell)

         if (indexToCellID(iCell) == config_stats_cell_ID .and. config_print_thickness_advection_info) then
            call mpas_log_write('In advect_thickness_tracer, iCell (local=$i, global=$i) =', &
               intArgs=(/iCell, indexToCellID(iCell)/))
            call mpas_log_write('k, iEdgeOnCell, layerNormalVelocity, layerThicknessEdge, ' // &
               'thicknessFluxEdge, thicknessTendency, thicknessTracerTendency:')
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

               ! Increment the thickness and thickness*tracer tendencies

               thicknessFluxEdge = layerNormalVelocity(k,iEdge) * dvEdge(iEdge) * layerThicknessEdge(k,iEdge)
               thicknessTendency = thicknessTendency &
                                 + edgeSignOnCell(iEdgeOnCell,iCell) * thicknessFluxEdge * invAreaCell

               if (advectTracers) then

                  velSign = sign(1.0_RKIND, layerNormalVelocity(k,iEdge))

                  ! Assign the upwind tracer values to each edge
                  ! Note: layerThicknessEdge is passed in and does not need to be computed here

                  tracersEdge(:) = max(velSign * tracersOld(:,k,iCell1),   &
                                       velSign * (-1.0_RKIND) * tracersOld(:,k,iCell2))

                  thicknessTracerTendency(:) = thicknessTracerTendency(:) &
                                             + edgeSignOnCell(iEdgeOnCell,iCell) * thicknessFluxEdge * invAreaCell * tracersEdge(:)

               endif   ! advectTracers

               !WHL - debug
               if (indexToCellID(iCell) == config_stats_cell_ID .and. k==1 .and. config_print_thickness_advection_info) then
                  call mpas_log_write("$i $i $r $r $r $r $r", &
                        intArgs=(/k, iEdgeOnCell/), realArgs=(/layerNormalVelocity(k,iEdge), &
                        layerThicknessEdge(k,iEdge), thicknessFluxEdge, thicknessTendency, thicknessTracerTendency(1)/) )
               endif

            enddo      ! iEdgeOnCell

            ! Compute the new layer thickness and tracers

            if (advectTracers) then

               ! new thickness*tracer products
               newThicknessTracers(:) = layerThicknessOld(k,iCell)*tracersOld(:,k,iCell) + thicknessTracerTendency(:) * dt

               ! new layer thickness
               layerThicknessNew(k,iCell) = layerThicknessOld(k,iCell) + thicknessTendency * dt

               ! new tracers
               if (layerThicknessNew(k,iCell) > 0.0_RKIND) then
                  tracersNew(:,k,iCell) = newThicknessTracers(:) / layerThicknessNew(k,iCell)
               else
                  tracersNew(:,k,iCell) = 0.0_RKIND
               endif

            else   ! advecting thickness only

               layerThicknessNew(k,iCell) = layerThicknessOld(k,iCell) + thicknessTendency * dt

            endif  ! advectTracers

         enddo     ! k

      enddo   ! iCell


      if (checkConservation) then

         ! compute final ice volume (= area*thickness) and energy (= area*thickness*tracer1)
         finalVolumeSum = 0.0_RKIND
         finalEnergySum = 0.0_RKIND
         do iCell = 1, nCellsSolve
            finalVolumeSum = finalVolumeSum + areaCell(iCell) * sum(layerThicknessNew(:,iCell))
            finalEnergySum = finalEnergySum + areaCell(iCell) * sum(layerThicknessNew(:,iCell)*tracersNew(1,:,iCell))
         enddo

         if (config_print_thickness_advection_info) then
            call mpas_log_write('init volume, final volume, difference: $r $r $r' , &
               realArgs=(/initVolumeSum, finalVolumeSum, finalVolumeSum - initVolumeSum/))
            call mpas_log_write('init energy, final energy, difference: $r $r $r', &
               realArgs=(/initEnergySum, finalEnergySum, finalEnergySum - initEnergySum/))
         endif

         !TODO - Make these fatal errors instead of warnings?
         difference = abs(finalVolumeSum - initVolumeSum)
         if (difference/initVolumeSum > eps11) then
            call mpas_log_write('upwind advection, mass conservation error', MPAS_LOG_WARN)
            call mpas_log_write('init volume, final volume, difference: $r $r $r', &
               realArgs=(/initVolumeSum, finalVolumeSum, finalVolumeSum - initVolumeSum/))
         endif

         difference = abs(finalEnergySum - initEnergySum)
         if (difference/initEnergySum > eps11) then
            call mpas_log_write('upwind advection, mass*tracer conservation error', MPAS_LOG_WARN)
            call mpas_log_write('init energy, final energy, difference: $r $r $r', &
               realArgs=(/initEnergySum, finalEnergySum, finalEnergySum - initEnergySum/))
         endif

      endif

      ! clean up
      deallocate(tracersEdge)
      deallocate(thicknessTracerTendency)
      deallocate(newThicknessTracers)

      ! === error check
      if (err > 0) then
          call mpas_log_write("An error has occurred in advect_thickness_tracers_upwind.", MPAS_LOG_ERR)
      endif

   !--------------------------------------------------------------------

    end subroutine advect_thickness_tracers_upwind

!***********************************************************************

!  routine li_layer_normal_velocity
!
!> \brief   Compute layer normal velocity and diagnose the advective CFL limit
!> \author  William Lipscomb
!> \date    January 2016
!> \details
!>  This routine computed layer normal velocities on each edge in preparation
!>  for advection. Given these velocities, it diagnoses the advective CFL limit
!>  for the edges on this block. The logic here was originally part of
!>  subroutine advect_thickness_tracers_upwind, but was moved to a separate
!>  subroutine so that an adaptive timestep could be diagnosed at the start
!>  of the time step.
!
!-----------------------------------------------------------------------

    subroutine li_layer_normal_velocity(&
         meshPool,               &
         normalVelocity,         &
         edgeMask,               &
         layerNormalVelocity,    &
         minOfMaxAllowableDt,    &
         err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
           meshPool              !< Input: mesh information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
           normalVelocity        !< Input: normal velocity on cell edges

      integer, dimension(:), intent(in) :: &
           edgeMask  !< Input: normal velocity on cell edges

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
           layerNormalVelocity     !< Output: normal velocity on cell edges, averaged to layer midpoints

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
           nEdgesSolve,     & ! dnumber of locally owned edges
           nVertLevels        ! number of vertical layers

      real (kind=RKIND), dimension(:), pointer :: &
           dcEdge             ! distance between cell centers on each edge

      logical, pointer :: &
           config_print_thickness_advection_info

      integer, dimension(:), pointer :: &
           indexToEdgeID      ! global index for edges

      real (kind=RKIND) :: &
           maxAllowableDt     ! max allowable dt based on advective CFL condition

      integer :: iEdge, k

      real (kind=RKIND), parameter :: bigNumber = 1.0e11_RKIND
      ! This is ~5000 years in seconds, but is small enough not to overflow

      err = 0

      ! get dimensions
      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      ! get mesh arrays
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'indexToEdgeID', indexToEdgeID)

      ! get config variables
      call mpas_pool_get_config(liConfigs, 'config_print_thickness_advection_info', config_print_thickness_advection_info)

      ! initialize output variables
      minOfMaxAllowableDt = bigNumber

      ! loop over local edges
      do iEdge = 1, nEdgesSolve

         if (li_mask_is_dynamic_ice(edgeMask(iEdge))) then

         ! loop over layers
         do k = 1, nVertLevels

            ! average normal velocities from layer interfaces to layer midpoints for advection
            layerNormalVelocity(k,iEdge) = 0.5_RKIND * (normalVelocity(k,iEdge) + normalVelocity(k+1,iEdge))

            ! Check for potential CFL violation
            if (abs(layerNormalVelocity(k,iEdge)) > 0.0_RKIND) then
               maxAllowableDt = (0.5_RKIND * dcEdge(iEdge)) / abs(layerNormalVelocity(k,iEdge))
            else
               maxAllowableDt = bigNumber
            endif

            minOfMaxAllowableDt = min(minOfMaxAllowableDt, maxAllowableDt)

         enddo     ! k

         else
            layerNormalVelocity(:,iEdge) = 0.0_RKIND
         endif

      enddo   ! iEdge

      ! === error check
      if (err > 0) then
          call mpas_log_write("An error has occurred in li_layer_normal_velocity.", MPAS_LOG_ERR)
      endif

   !--------------------------------------------------------------------

    end subroutine li_layer_normal_velocity

!-----------------------------------------------------------------------
!
!  subroutine vertical_remap
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
!
!-----------------------------------------------------------------------
   subroutine vertical_remap(thickness, cellMask, meshPool, layerThickness, tracers, err)

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
      integer, dimension(:), pointer :: indexToCellID

      ! local arrays
      real (kind=RKIND), dimension(:), allocatable :: layerInterfaceSigma_Input
      real (kind=RKIND), dimension(:,:), allocatable :: hTsum

      ! counters, mesh variables, index variables
      integer, pointer :: nCells, nVertLevels
      integer :: nTracers, iCell, k, k1, k2, nt

      ! stuff for making calculations
      real(kind=RKIND) :: zhi, zlo, hOverlap

      ! variables for optional conservation check
      ! set to 'true' for now, since it is inexpensive and might catch problems
      logical, parameter :: checkConservation = .true.
      real (kind=RKIND) :: initEnergySum, finalEnergySum, difference
      real (kind=RKIND), parameter :: eps11 = 1.e-11_RKIND

      err = 0

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      nTracers = size(tracers, 1)

      call mpas_pool_get_array(meshPool, 'layerThicknessFractions', layerThicknessFractions)
      call mpas_pool_get_array(meshPool, 'layerInterfaceSigma', layerInterfaceSigma)
      call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)

      allocate(layerInterfaceSigma_Input(nVertLevels+1))
      allocate(hTsum(nTracers, nVertLevels))

      ! loop over cells
      do iCell = 1, nCells

         if (checkConservation) then   ! compute sum of layerThickness*tracer1 in the column
            initEnergySum = sum(layerThickness(:,iCell)*tracers(1,:,iCell))
         endif

         if (thickness(iCell) > 0.0_RKIND) then

            ! *** Calculate vertical sigma coordinates of each layer interface for the input non-sigma state
            !     (we already have the desired new sigma-based state as layerInterfaceSigma)

            layerInterfaceSigma_Input(1) = 0.0_RKIND
            do k = 2, nVertLevels
               layerInterfaceSigma_Input(k) = layerInterfaceSigma_Input(k-1) + layerThickness(k-1, iCell) / thickness(iCell)
            end do
            layerInterfaceSigma_Input(nVertLevels+1) = 1.0_RKIND

            ! *** Compute new layer thicknesses (layerInterfaceSigma coordinates)

            do k = 1, nVertLevels
               layerThickness(k,iCell) = layerThicknessFractions(k) * thickness(iCell)
            end do

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
                     hOverlap = max (zhi-zlo, 0.0_RKIND) * thickness(iCell)
                     hTsum(nt,k2) = htsum(nt,k2) + tracers(nt,k1,iCell) * hOverlap
                  enddo         ! nt
               enddo            ! k1
            enddo               ! k2

            !-----------------------------------------------------------------
            ! Compute tracer values in new layers
            !-----------------------------------------------------------------

            do k = 1, nVertLevels
               do nt = 1, nTracers
                  tracers(nt,k,iCell) = hTsum(nt,k) / (layerThickness(k, iCell))
               enddo         ! nt
            enddo            ! k

         endif    ! thickness > 0.0

         if (checkConservation) then   ! compute new sum of layerThickness*tracer1 in the column
            finalEnergySum = sum(layerThickness(:,iCell)*tracers(1,:,iCell))
            difference = abs(finalEnergySum - initEnergySum)
            if (initEnergySum > eps11) then
               if (difference/initEnergySum > eps11) then
                  call mpas_log_write('vertical_remap, mass*tracer conservation error, indexToCellID = $i', &
                     MPAS_LOG_WARN, intArgs=(/indexToCellID(iCell)/))
                  call mpas_log_write('init energy, final energy, difference: $r $r $r', &
                     realArgs=(/initEnergySum, finalEnergySum, finalEnergySum - initEnergySum/))
               endif
            endif
         endif

      enddo ! nCellsSolve

      ! clean up
      deallocate(layerInterfaceSigma_Input)
      deallocate(hTsum)

    end subroutine vertical_remap

    subroutine grounded_to_floating(cellMaskOrig, cellMaskNew, thicknessNew, groundedToFloatingThickness, nCells)
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      integer, dimension(:), intent(in) :: &
         cellMaskOrig          !< Input: mask for cells before advection

      integer, dimension(:), intent(in) :: &
         cellMaskNew          !< Input: mask for cells after advection

      real(kind=RKIND), dimension(:), intent(in) :: &
         thicknessNew         !< Input: ice thickness after advection

      integer, pointer, intent(in) :: &
         nCells
      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      real(kind=RKIND), dimension(:), intent(out) :: &
         groundedToFloatingThickness         !< Input: ice thickness

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      integer :: iCell

      groundedToFloatingThickness(:) = 0.0_RKIND

      do iCell = 1, nCells
        ! find locations that had grounded ice before but now have floating ice
        if (li_mask_is_grounded_ice(cellMaskOrig(iCell)) .and. &
            li_mask_is_floating_ice(cellMaskNew(iCell)) ) then
           groundedToFloatingThickness(iCell) = thicknessNew(iCell)
        ! find locations that had floating ice before but now have grounded ice
        else if (li_mask_is_floating_ice(cellMaskOrig(iCell)) .and. &
                 li_mask_is_grounded_ice(cellMaskNew(iCell)) ) then
           groundedToFloatingThickness(iCell) = -1.0_RKIND * thicknessNew(iCell)
        endif
      enddo

    end subroutine grounded_to_floating

!***********************************************************************

 end module li_advection

