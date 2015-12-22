! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_tendency
!
!> \brief MPAS land ice tendency driver
!> \author Matt Hoffman
!> \date   17 April 2011
!> \details
!>  This module contains the routines for computing
!>  various tendencies for land ice. As well as routines
!>  for computing diagnostic variables.
!>  Note: This module will likely be deprecated in the near future
!>  (WHL, Dec. 2015)
!-----------------------------------------------------------------------

module li_tendency

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
   public :: li_tendency_thickness

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   !WHL - debug
   integer, parameter :: iTest = 485

 contains


!***********************************************************************
!
!  subroutine li_tendency_thickness
!
!> \brief   Computes tendency term from horizontal advection of thickness
!> \author  Matt Hoffman
!> \date    16 April 2012
!> \details 
!>  This routine computes the horizontal advection tendency for
!>  thickness based on current state and user choices of forcings. Based on
!>  ocn_thick_hadv_tend in the ocean core.
!
!-----------------------------------------------------------------------

    subroutine li_tendency_thickness(meshPool, velocityPool, geometryPool, layerThicknessTendency, dt, dminfo, allowableDt, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      type (mpas_pool_type), intent(in) :: &
         velocityPool      !< Input: velocity information

      real (kind=RKIND), intent(in) :: &
         dt       !< Input: dt

      type (dm_info), pointer, intent(in) :: &
         dminfo      !< Input: domain info

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: &
         geometryPool         !< Input: geometry information to be updated

      real (kind=RKIND), dimension(:,:), pointer, intent(inout) :: &
         layerThicknessTendency         !< Input/Output: layer thickness tendency

!!!      real (kind=RKIND), dimension(:), pointer, intent(inout) :: &
!!!         iceArea_tend         !< Input/Output: ice area tendency


      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag
      real (kind=RKIND), intent(out) :: allowableDt

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      integer, pointer :: nVertLevels
      real (kind=RKIND), dimension(:), pointer :: sfcMassBal
      real (kind=RKIND), dimension(:), pointer :: basalMassBal
      real (kind=RKIND), dimension(:), pointer :: groundedBasalMassBal
      real (kind=RKIND), dimension(:), pointer :: floatingBasalMassBal
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity, layerThicknessEdge
      integer, dimension(:), pointer :: cellMask, edgeMask
      character (len=StrKIND), pointer :: config_thickness_advection
      real (kind=RKIND), pointer :: config_ice_density
      integer :: err_tmp

      !WHL - debug
      integer, dimension(:), pointer :: indexToCellID
      integer, pointer :: nCells
      integer, pointer :: config_stats_cell_ID
      logical, pointer :: config_print_thickness_advection_info
      integer :: iEdge, iEdgeOnCell
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell

      !WHL - debug - for test-case diagnostics
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness
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

      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(geometryPool, 'sfcMassBal', sfcMassBal)
      call mpas_pool_get_array(geometryPool, 'basalMassBal', basalMassBal)
      call mpas_pool_get_array(geometryPool, 'groundedBasalMassBal', groundedBasalMassBal)
      call mpas_pool_get_array(geometryPool, 'floatingBasalMassBal', floatingBasalMassBal)
      ! Assuming tendency will always be calculated using time level 1!
      call mpas_pool_get_array(velocityPool, 'normalVelocity', normalVelocity)
      call mpas_pool_get_array(geometryPool, 'layerThicknessEdge', layerThicknessEdge)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_array(geometryPool, 'edgeMask', edgeMask)
      call mpas_pool_get_config(liConfigs, 'config_thickness_advection', config_thickness_advection)
      call mpas_pool_get_config(liConfigs, 'config_ice_density', config_ice_density)

      !WHL - debug
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
      call mpas_pool_get_config(liConfigs, 'config_stats_cell_ID', config_stats_cell_ID)
      call mpas_pool_get_config(liConfigs, 'config_print_thickness_advection_info', config_print_thickness_advection_info)

     ! 0 tendency
      layerThicknessTendency = 0.0_RKIND
!!!     iceArea_tend = 0.0_RKIND

      select case (trim(config_thickness_advection))
      case ('fo')  !===================================================

         if (config_print_thickness_advection_info) then
            write(stderrUnit,*) 'Using first-order upwind for thickness advection'
            !WHL - debug
            write(stderrUnit,*) 'itest, cellMask:', itest, cellMask(itest)
         endif

!!!        ! Alternate call to calculate thickness tendency for entire ice column (deprecated)
!!!        !call li_tend_h_fo_upwind(mesh, state % normalVelocity % array, &
!!!        !   state % layerThickness % array, state % thickness % array, thickness_tend, dt, err)

         !WHL - debug
         if (dome_test .and. config_print_thickness_advection_info) then
            write(stderrUnit,*) ' '
            write(stderrUnit,*) 'Dome test diagnostics:'
            write(stderrUnit,*) ' '
            call mpas_pool_get_array(geometryPool, 'layerThickness', layerThickness, timeLevel=2)  ! timeLevel 2 is value from end of previous step
            iLayer = 1
!            write(stderrUnit,*) 'Initial thickness, layer =', iLayer
            do iRow = nRows, 1, -1
               if (mod(iRow,2) == 0) then  ! indent for even-numbered rows
!                  write(stderrUnit,'(a3)',advance='no') '    '
               endif
               do i = nCellsPerRow/2 - 2, nCellsPerRow
                  iCell = (iRow-1)*nCellsPerRow + i
!                  write(stderrUnit,'(f8.3)',advance='no') layerThickness(iLayer,iCell)
               enddo
!               write(stderrUnit,*) ' '
            enddo
         endif

         !WHL - debug
         if (config_print_thickness_advection_info) then
            do iCell = 1, nCells
!               write(stderrUnit,*) 'iCell:', iCell
               if (indexToCellID(iCell) == config_stats_cell_ID) then
                  write(stderrUnit,*) 'Before thickness advection, iCell (local, global), thickness:', &
                       iCell, indexToCellID(iCell), sum(layerThickness(:,iCell))
                  do iEdgeOnCell = 1, nEdgesOnCell(iCell)
                     iEdge = edgesOnCell(iEdgeOnCell,iCell)
                     write(stderrUnit,*) 'iEdge, normalVelocity:', iEdge, normalVelocity(1,iEdge)
                  enddo
               endif
            enddo
         endif

        call tend_layerthickness_fo_upwind(meshPool, normalVelocity, &
             layerThicknessEdge, edgeMask, layerThicknessTendency, dt, allowableDt, err_tmp)

         !WHL - debug
         if (config_print_thickness_advection_info) then
            do iCell = 1, nCells
               if (indexToCellID(iCell) == config_stats_cell_ID) then
                  write(stderrUnit,*) ' '
                  write(stderrUnit,*) 'After thickness advection, iCell, new layer thickness:'
                  do k = 1, nVertLevels
                     write(stderrUnit,*) iCell, layerThickness(k,iCell) + layerThicknessTendency(k,iCell)*dt
                  enddo
               endif
            enddo
         endif

         !WHL - debug
         if (dome_test) then
!            write(stderrUnit,*) ' '
            iLayer = 1
!            write(stderrUnit,*) 'New thickness, layer =', iLayer
            do iRow = nRows, 1, -1
               if (mod(iRow,2) == 0) then  ! indent for even-numbered rows
!                  write(stderrUnit,'(a3)',advance='no') '    '
               endif
               do i = nCellsPerRow/2 - 2, nCellsPerRow
                  iCell = (iRow-1)*nCellsPerRow + i
!                  write(stderrUnit,'(f8.3)',advance='no') layerThickness(iLayer,iCell) + layerThicknessTendency(iLayer,iCell)*dt
               enddo
!               write(stderrUnit,*) ' '
            enddo
!            write(stderrUnit,*) ' '
!            write(stderrUnit,*) 'layerThicknessTendency*dt, layer =', iLayer
            do iRow = nRows, 1, -1
               if (mod(iRow,2) == 0) then  ! indent for even-numbered rows
!                  write(stderrUnit,'(a3)',advance='no') '    '
               endif
               do i = nCellsPerRow/2 - 2, nCellsPerRow
                  iCell = (iRow-1)*nCellsPerRow + i
!                  write(stderrUnit,'(f8.3)',advance='no') layerThicknessTendency(iLayer,iCell)*dt
               enddo
!               write(stderrUnit,*) ' '
            enddo
         endif


        !!! Experimental implementation of using FCT for thickness tendency calculation
        ! This requires additional work, including a separate tendency calculation for 
        ! advancing cells using FO-Upwind.
     !case ('FCT')  !===================================================
        ! MJH TEMP TRYING IT WITH THICKNESS
        !      stateOld % sup_thickness % array(1,:,:) = stateOld % layerThickness % array  !!!!!  MJH TEMP
        !      block % tend % sup_thickness % array = 0.0
!!! NOTE FOR FUTURE: These two calls to tracer_advection have changed since this code was written - they are now in framework and the calling signature has changed!
        !         call mpas_tracer_advection_tend(stateOld % sup_thickness % array, uh , &
        !                  wTop,   stateOld % sup_thickness % array(1,:,:), stateOld % sup_thickness % array(1,:,:), &
        !                  dt/SecondsInYear  , mesh,   0.0*layerThickness_tend,    block % tend % sup_thickness % array)
 
        ! Doug suggested actually doing it this way:
        !         call mpas_tracer_advection_tend(stateOld % sup_thickness % array, normalVelocity , &
        !                  wTop,   stateOld % sup_thickness % array(1,:,:) * 0.0 + 1.0 , stateOld % sup_thickness % array(1,:,:) * 0.0 + 1.0, &
        !                  dt/SecondsInYear  , mesh,   0.0*layerThickness_tend,    block % tend % sup_thickness % array)

        ! where (stateOld % sup_thickness % array .gt. 1.0e-9)
        !  block % tend % sup_thickness % array = block % tend % sup_thickness % array / stateOld % sup_thickness % array
        ! else where 
        !  block % tend % sup_thickness % array = 0.0
        ! end where
         ! assign to the thickness tend i actually use
        ! layerThickness_tend = block % tend % sup_thickness % array(1,:,:)

      case ('none')  !===================================================
        ! Do nothing
      case default  !===================================================
        write(stderrUnit,*) trim(config_thickness_advection), ' is not a valid thickness advection option.'
        err_tmp = 1
      end select  !===================================================
      err = ior(err,err_tmp)

      ! Add the mass balance to the tendencies
      select case (trim(config_thickness_advection))
      case ('None')  !===================================================
        ! Do nothing - don't add the MB
      case default

        ! Combine various basal mass balnce fields based on mask
        ! 1. make adjustments for where the ice is grounded and floating.
        !    TODO: more complicated treatment at GL?
        where ( li_mask_is_grounded_ice(cellMask) )
           ! Apply marineBasalMassBal to floating ice only.
           basalMassBal = groundedBasalMassBal
        elsewhere ( li_mask_is_floating_ice(cellMask) )
           ! Currently, floating and grounded ice are mutually exclusive.
           ! This could change if the GL is parameterized, in which case this logic may need adjustment.
           ! Grounded BMB should come from the temperature solver.
           basalMassBal = floatingBasalMassBal
        elsewhere ( .not. (li_mask_is_ice(cellMask) ) )
           ! We don't allow a positive basal mass balance where ice is not already present.
           basalMassBal = 0.0_RKIND
        end where

        ! Add surface mass balance to tendency
        ! TODO: Need to decide how to deal with negative SMB that eliminates top layer or all ice (check for negative thickness?)
        layerThicknessTendency(1,:) = layerThicknessTendency(1,:) + sfcMassBal / config_ice_density  ! (tendency in meters per year)
        ! TODO THIS MIGHT RESULT IN  NEGATIVE LAYER THICKNESS!

        ! Add basal mass balance to tendency
        ! TODO: Need to decide how to deal with negative BMB that eliminates top layer or all ice (check for negative thickness?)
        layerThicknessTendency(nVertLevels,:) = layerThicknessTendency(nVertLevels,:) &
                                    + basalMassBal / config_ice_density    ! (tendency in meters per year)
        ! TODO THIS MIGHT RESULT IN  NEGATIVE LAYER THICKNESS!

      end select

      ! === error check
      if (err > 0) then
         write (stderrUnit,*) "An error has occurred in li_tendency_thickness."
      endif

      !--------------------------------------------------------------------
    end subroutine li_tendency_thickness

!***********************************************************************
! Private subroutines:
!***********************************************************************


!***********************************************************************
!
!  subroutine tend_layerThickness_fo_upwind
!
!> \brief   Computes tendency term from horizontal advection of thickness layers
!> \author  Matt Hoffman
!> \date    16 April 2012
!> \details 
!>  This routine computes the horizontal advection tendency for each
!>  thickness layer using first-order upwinding. Based on
!>  ocn_thick_hadv_tend in the ocean core.  This is an alternative to lice_tend_h
!>  that calculates the tendency for each layer, which would then need to be
!>  added up to calculate the change in thickness.  The two methods yield identical
!>  results.  
!
!-----------------------------------------------------------------------
   subroutine tend_layerThickness_fo_upwind(meshPool, normalVelocity, layerThicknessEdge, edgeMask, tend, dt, MinOfMaxAllowableDt, err)

      use mpas_timekeeping

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity    !< Input: velocity

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThicknessEdge    !< Input: thickness of each layer on edges

      integer, dimension(:), intent(in) :: &
         edgeMask    !< Input: mask on edges

      real (kind=RKIND) :: dt

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend         !< Input/Output: velocity tendency

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag
      real (kind=RKIND) :: MinOfMaxAllowableDt

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer, pointer :: nCells, nVertLevels
      integer, dimension(:,:), pointer :: edgesOnCell, edgeSignOnCell
      integer, dimension(:), pointer :: nEdgesOnCell
      real (kind=RKIND), dimension(:), pointer :: dvEdge, areaCell, dcEdge
      real (kind=RKIND) :: invAreaCell, flux, maxAllowableDt, layerNormalVelocity
      integer :: iEdge, iCell, i, k
      real (kind=RKIND), parameter :: bigNumber = 1.0e16_RKIND  ! This is ~300 million years in seconds, but it is small enough not too overflow
      integer :: err_tmp

      ! Only needed for optional check for mass conservation
      !real (kind=RKIND) :: tendVolSum

      !WHL - debug
      logical, pointer :: config_print_thickness_advection_info

      err = 0

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

      !WHL - debug
      call mpas_pool_get_config(liConfigs, 'config_print_thickness_advection_info', config_print_thickness_advection_info)

      MinOfMaxAllowableDt = bigNumber

      do iCell = 1, nCells
        invAreaCell = 1.0 / areaCell(iCell)

        !WHL - debug
         if (iCell == itest .and. config_print_thickness_advection_info) then
            write(stderrUnit,*) 'In tend_layerThickness, iCell =', iCell
            write(stderrUnit,*)
            write(stderrUnit,*) 'k, iEdgeOnCell, layerNormalVelocity, layerThicknessEdge, thicknessFluxEdge, thicknessTendency' 
         endif

        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)

          do k = 1, nVertLevels

            ! Average native velocities from layer interfaces to layer midpoints for advection
            ! TODO This may make more sense to calculate as a 3d field somewhere else so it can also be used for tracer advection, output visualization, etc.
            layerNormalVelocity = 0.5_RKIND * (normalVelocity(k, iEdge) + normalVelocity(k+1, iEdge))

            if (abs(layerNormalVelocity) > 0.0_RKIND) then
               maxAllowableDt = (0.5_RKIND * dcEdge(iEdge)) / abs(layerNormalVelocity)
            else
               maxAllowableDt = bigNumber
            endif
            if ( maxAllowableDt < dt ) then
                 !write(stderrUnit,*) 'CFL violation at level, edge', k, iEdge                      
                 err = err + 1
            endif
            MinOfMaxAllowableDt = min(MinOfMaxAllowableDt, maxAllowableDt)

            flux = layerNormalVelocity * dvEdge(iEdge) * layerThicknessEdge(k, iEdge)
            tend(k, iCell) = tend(k, iCell) + edgeSignOnCell(i, iCell) * flux * invAreaCell

            !WHL - debug
            if (iCell == itest .and. k==1 .and. config_print_thickness_advection_info) then
               write(stderrUnit,*) k, i, layerNormalVelocity, layerThicknessEdge(k,iEdge), flux, tend(k,iCell)
            endif

          end do
        end do
      end do


      ! Optional check for mass conservation
      !tendVolSum = 0.0_RKIND
      !do iCell=1, mesh % nCells
      !     tendVolSum  = tendVolSum + sum(tend(:,iCell)) * areaCell(iCell)
      !end do
      !print *,'SUM OF VOLUME TENDENCY  ======', tendVolSum

      ! === error check
      if (err > 0) then
          write (stderrUnit,*) "An error has occurred in tend_layerThickness_fo_upwind."
      endif

   !--------------------------------------------------------------------

   end subroutine tend_layerThickness_fo_upwind!}}}

!***********************************************************************

end module li_tendency

