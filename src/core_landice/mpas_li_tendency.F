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
!> \version SVN:$Id:$
!> \details
!>  This module contains the routines for computing
!>  various tendencies for land ice. As well as routines
!>  for computing diagnostic variables.
!
!-----------------------------------------------------------------------

module li_tendency

   use mpas_grid_types
   use mpas_configure
   use mpas_constants
   use mpas_dmpar
!!!   use li_mask
!!!   use mpas_tracer_advection

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
   public :: li_tendency_thickness, &
             li_tendency_tracers,   &
             li_apply_calving

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------



   contains


!***********************************************************************
!
!  subroutine li_tendency_thickness
!
!> \brief   Computes tendency term from horizontal advection of thickness
!> \author  Matt Hoffman
!> \date    16 April 2012
!> \version SVN:$Id$
!> \details 
!>  This routine computes the horizontal advection tendency for
!>  thickness based on current state and user choices of forcings. Based on
!>  ocn_thick_hadv_tend in the ocean core.
!
!-----------------------------------------------------------------------

   subroutine li_tendency_thickness(mesh, state, layerThickness_tend, dt, dminfo, allowableDt, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         mesh          !< Input: mesh information

      real (kind=RKIND), intent(in) :: &
         dt       !< Input: dt

      type (dm_info), pointer, intent(in) :: &
         dminfo      !< Input: domain info

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (state_type), intent(inout) :: &
         state         !< Input: state to use to calculate tendency (old time level)
         ! Note: state needs to be inout (rather than just in) so that adjust_marine_boundary_fluxes can modify it.

      real (kind=RKIND), dimension(:,:), pointer, intent(inout) :: &
         layerThickness_tend         !< Input/Output: layer thickness tendency

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
      integer nVertLevels
      real (kind=RKIND), dimension(:), pointer :: sfcMassBal
!!!      real (kind=RKIND), dimension(:), pointer :: iceArea, areaCell, marineBasalMassBal
      integer, dimension(:), pointer :: cellMask
      integer :: err_tmp


      err = 0
      err_tmp = 0

     nVertLevels = mesh % nVertLevels
     sfcMassBal => mesh % sfcMassBal % array
!!!     marineBasalMassBal => mesh % marineBasalMassBal % array
!!!     iceArea => state % iceArea % array
!!!     areaCell => mesh % areaCell % array
     cellMask => state % cellMask % array


     ! 0 tendency
     layerThickness_tend = 0.0_RKIND
!!!     iceArea_tend = 0.0_RKIND

     select case (config_thickness_advection)
     case ('fo')  !===================================================
        !print *,'Using FO Upwind for thickness advection'

!!!        ! Alternate call to calculate thickness tendency for entire ice column (deprecated)
!!!        !call li_tend_h_fo_upwind(mesh, state % normalVelocity % array, &
!!!        !   state % layerThickness % array, state % thickness % array, thickness_tend, dt, err)

        call tend_layerthickness_fo_upwind(mesh, state % normalVelocity % array, &
             state % layerThicknessEdge % array, state % edgeMask % array, layerThickness_tend, dt, allowableDt, err_tmp)       


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
        write(0,*) trim(config_thickness_advection), ' is not a valid thickness advection option.'
        err_tmp = 1
     end select  !===================================================
     err = ior(err,err_tmp)


     ! Add the MB to the tendencies
     select case (config_thickness_advection)
     case ('None')  !===================================================
        ! Do nothing - don't add the MB
     case default
! Commenting BMB out for now.
!!!        ! Make some potential adjustments to BMB before applying them.
!!!        !    It's ok to overwrite the values with 0's here, because each time step
!!!        !    we get a fresh copy of the array from the annual_forcing subroutine.
!!!        ! 1. make adjustments for where the ice is grounded and floating.
!!!        !    TODO: more complicated treatment at GL?
!!!        where ( li_mask_is_grounded_ice(cellMask) )
!!!           ! Apply marineBasalMassBal to floating ice only.  
!!!           marineBasalMassBal = 0.0_RKIND
!!!        elsewhere ( li_mask_is_floating_ice(cellMask) )
!!!           ! Currently, floating and grounded ice are mutually exclusive.
!!!           ! This could change if the GL is parameterized, in which case this logic may need adjustment.
!!!           ! Grounded BMB should come from the temperature solver.
!!!           ! < PLACEHOLDER >
!!!        elsewhere ( .not. (li_mask_is_ice(cellMask) )
!!!              ! We don't allow a positive BMB where ice is not already present.
!!!              mesh % marineBasalMassBal % array = 0.0_RKIND
!!!        end where

        ! Add surface mass balance to tendency
        ! TODO: Need to decide how to deal with negative SMB that eliminates top layer or all ice (check for negative thickness?)
        layerThickness_tend(1,:) = layerThickness_tend(1,:) + sfcMassBal / config_ice_density  ! (tendency in meters per year)
        ! TODO THIS MIGHT RESULT IN  NEGATIVE LAYER THICKNESS!

!!!        ! Add basal mass balance to tendency
!!!        ! TODO: Need to decide how to deal with negative BMB that eliminates top layer or all ice (check for negative thickness?)
!!!        layerThickness_tend(nVertLevels,:) = layerThickness_tend(nVertLevels,:) &
!!!                                    + mesh % marineBasalMassBal % array    ! (tendency in meters per year)
!!!        ! TODO Add in grounded ice basal mass balance once temperature diffusion is calculated
!!!        ! TODO THIS MIGHT RESULT IN  NEGATIVE LAYER THICKNESS!
     end select

     ! === error check
     if (err > 0) then
         write (0,*) "An error has occurred in li_tendency_thickness."
     endif

   !--------------------------------------------------------------------
   end subroutine li_tendency_thickness



!***********************************************************************
!
!  subroutine li_tendency_tracers
!
!> \brief   Computes tendency term from horizontal advection of thickness
!> \author  Matt Hoffman
!> \date    16 April 2012
!> \version SVN:$Id$
!> \details 
!>  This routine computes the horizontal advection tendency for
!>  thickness based on current state and user choices of forcings. Based on
!>  ocn_thick_hadv_tend in the ocean core.
!
!-----------------------------------------------------------------------

   subroutine li_tendency_tracers(mesh, state, layerThickness_tend, tracer_tendency, dt, dminfo, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         mesh          !< Input: mesh information

      type (state_type), intent(in) :: &
         state         !< Input: state to use to calculate tendency

      real (kind=RKIND), dimension(:,:), pointer, intent(in) :: &
         layerThickness_tend         !< Input/Output: layer thickness tendency     

      real (kind=RKIND), intent(in) :: &
         dt       !< Input: dt

      type (dm_info), pointer, intent(in) :: &
         dminfo      !< Input: domain info

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), pointer, intent(inout) :: &
         tracer_tendency         !< Input/Output: tracers tendency

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
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity, layerThickness, layerThicknessEdge
      real (kind=RKIND), dimension(:,:), allocatable :: wTop, uh
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers
      integer :: iEdge, k

! Commenting out this entire subroutine for now.
!!!     ! 0 tendency
!!!     tracer_tendency = 0.0_RKIND

!!!     select case (config_tracer_advection)
!!!     case ('FO') !===================================================
!!!         ! Assign pointers, etc.
!!!         layerThicknessEdge => state % layerThicknessEdge % array
!!!         normalVelocity => state % normalVelocity % array
!!!         tracers => state % tracers % array
!!!         allocate(uh(mesh % nVertLevels, mesh % nEdges + 1))

!!!         ! FO code needs u*h on edges - this could be calculated elsewhere and saved, potentially (e.g. it's already calculated locally in the FO upwind thickness routine)
!!!         ! layerThicknessEdge and normalVelocity should match those used in thickness advection!
!!!         do iEdge = 1, mesh % nEdges
!!!           do k = 1, mesh % nVertLevels
!!!               uh(k, iEdge) = normalVelocity(k, iEdge) * layerThicknessEdge(k, iEdge)
!!!           end do
!!!         end do

!!!         ! Call the FO!  
!!!         call tracer_advection_tend_fo(tracers, uh, mesh, tracer_tendency, err)

!!!         deallocate(uh)

!!!     case ('FCT')  !===================================================
!!!         ! Assign pointers, etc.
!!!         layerThickness => state % layerThickness % array
!!!         layerThicknessEdge => state % layerThicknessEdge % array
!!!         normalVelocity => state % normalVelocity % array
!!!         tracers => state % tracers % array
!!!         allocate(wTop(mesh % nVertLevels + 1, mesh % nCells + 1))
!!!         allocate(uh(mesh % nVertLevels, mesh % nEdges + 1))

!!!         ! Setup 0 vertical velocity field - vertical velocity not needed with sigma coordinates
!!!         wTop = 0.0_RKIND

!!!         ! prepare for FCT call

!!!         ! This is a kluge because the FCT code will divide by zero if thickness is zero \todo: deal with more elegantly!
!!!         where (layerThickness == 0.0_RKIND)
!!!             layerThickness = 1.0e-12_RKIND
!!!         end where

!!!         ! FCT code needs u*h on edges - this could be calculated elsewhere and saved, potentially (e.g. it's already calculated locally in the FO upwind thickness routine)
!!!         ! layerThicknessEdge and normalVelocity should match those used in thickness advection!
!!!         do iEdge = 1, mesh % nEdges
!!!           do k = 1, mesh % nVertLevels
!!!               uh(k, iEdge) = normalVelocity(k, iEdge) * layerThicknessEdge(k, iEdge)
!!!           end do
!!!         end do

!!!         !print *,'uh', maxval(uh), minval(uh)
!!!         !print *,'wTop', maxval(wTop), minval(wTop)
!!!         !print *,'layerThicknessOld', maxval(layerThicknessOld), minval(layerThicknessOld)
!!!         !print *,'layerThickness_tend', maxval(layerThickness_tend), minval(layerThickness_tend)
!!!         !print *,'  temp min/max:', minval(stateOld % tracers % array(stateOld%index_temperature,:,1:mesh%nCells)), maxval(stateOld % tracers % array(stateOld%index_temperature,:,1:mesh%nCells))

!!!         ! Call the FCT!  (this will likely move from the ocean core to operators eventually)
!!!         call mpas_ocn_tracer_advection_tend(tracers, uh , &
!!!                  wTop, layerThickness, layerThickness, dt/SecondsInYear,  &
!!!                  mesh, layerThickness_tend, tracer_tendency)
!!!         !print *,'tracer_tendency', maxval(tracer_tendency), minval(tracer_tendency)

!!!         ! Set thickness back to 0 where needed.  This is a kluge because the FCT code will divide by zero if thickness is zero \todo: deal with more elegantly!
!!!         where (layerThickness == 1.0e-12_RKIND)
!!!             layerThickness = 0.0_RKIND
!!!         end where


!!!         deallocate(wTop)
!!!         deallocate(uh)

!!!     case ('None')  !===================================================
!!!        ! Do nothing
!!!     case default  !===================================================
!!!        write(0,*) trim(config_tracer_advection), ' is not a valid tracer advection option.'
!!!        call mpas_dmpar_abort(dminfo)
!!!     end select  !===================================================

   !--------------------------------------------------------------------
   end subroutine li_tendency_tracers



!***********************************************************************
!
!  subroutine li_apply_calving
!
!> \brief   Applies a calving 'law' to any marine-terminating ice
!> \author  Matt Hoffman
!> \date    19 February 2013
!> \version SVN:$Id$
!> \details 
!> This routine is a driver for applying various calving laws to any
!> marine-terminating ice.
!
!-----------------------------------------------------------------------
   subroutine li_apply_calving(mesh, state, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         mesh          !< Input: mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (state_type), intent(inout) :: &
         state         !< Input/Output: state for which to update diagnostic variables

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
      integer, dimension(:), pointer :: cellMask, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell
      real (kind=RKIND), dimension(:), pointer :: thickness, iceArea, areaCell
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers
      integer :: iCell, nCells, numCalvedCells, nThickEdges, iEdge
      real (kind=RKIND) :: physicalThickness

! Commenting out this entire subroutine for now.
!!!      areaCell => mesh % areaCell % array
!!!      nEdgesOnCell => mesh % nEdgesOnCell % array
!!!      cellsOnCell => mesh % cellsOnCell % array

!!!      cellMask => state % cellMask % array
!!!      thickness => state % thickness % array
!!!      iceArea => state % iceArea % array
!!!      layerThickness => state % layerThickness % array
!!!      tracers => state % tracers % array

!!!      nCells = mesh % nCells

!!!      select case (config_calving_law)

!!!      !====================================
!!!      !case ('ocean-kill')

!!!      !====================================
!!!      case ('critical-thickness')
!!!         ! This case loops around the floating ice margin, removing any ice below
!!!         ! a critical thickness.  Multiple loops are needed until the entire region
!!!         ! of too-thin ice has been found.
!!!         ! This could be moved to its own subroutine...
!!!         numCalvedCells = 0

!!!         ! Loop across the entire domain multiple times until the inward progression of too-thin ice has been totally mapped
!!!         do 

!!!         ! Update the mask before we look for the margin, since the margin may have advanced.
!!!         call li_calculate_mask(mesh, state, err)

!!!         do iCell = 1, nCells
!!!            if ( (MASK_IS_FLOATING(cellMask(iCell))) .and. (MASK_IS_MARGIN(cellMask(iCell))) ) then
!!!               ! Convert to physical thickness - only needed if CFBC is on but this is always defined.
!!!               physicalThickness = thickness(iCell) * areaCell(iCell) / iceArea(iCell)  ! iceArea should always be > 0 since we are only checking on ice cells.
!!!               if (iceArea(iCell) == 0.0) then
!!!                  write(0,*) 'iceArea is 0 on a cell with ice!'
!!!                  err = 1
!!!               endif
!!!               if (physicalThickness <= config_calving_critical_thickness) then
!!!                  numCalvedCells = numCalvedCells + 1
!!!                  ! TODO tally a sum of mass lost to calving
!!!                  ! Make this cell ice-free and set tracers there to 0
!!!                  thickness(iCell) = 0.0_RKIND
!!!                  layerThickness(:,iCell) = 0.0_RKIND
!!!                  tracers(:,:,iCell) = 0.0_RKIND
!!!                  iceArea(iCell) = 0.0_RKIND
!!!               end if
!!!            end if
!!!         end do

!!!         if (numCalvedCells == 0) then 
!!!            ! if we haven't calved any cells this time though the loop, then we have mapped the extent of the too-thin ice
!!!            exit
!!!         else
!!!            numCalvedCells = 0
!!!         end  if

!!!         end do  ! start the search again

!!!         ! clean up - remove any hanging pieces that will not be dynamically active anymore. The mask should already be updated with current ice extent so we don't need to do it again.
!!!         do iCell = 1, nCells
!!!            if ( (MASK_IS_FLOATING(cellMask(iCell))) .and. (MASK_IS_MARGIN(cellMask(iCell))) ) then
!!!               nThickEdges = 0
!!!               do iEdge = 1, nEdgesOnCell(iCell)
!!!                  if ( MASK_IS_THICK_ICE(cellMask(cellsOnCell(iEdge, iCell))) ) then
!!!                     nThickEdges = nThickEdges + 1
!!!                  endif
!!!               end do
!!!               if (nThickEdges == 0) then
!!!                  ! Make this cell ice-free and set tracers there to 0
!!!                  thickness(iCell) = 0.0_RKIND
!!!                  layerThickness(:,iCell) = 0.0_RKIND
!!!                  tracers(:,:,iCell) = 0.0_RKIND
!!!                  iceArea(iCell) = 0.0_RKIND
!!!               end if
!!!            end if
!!!         end do

!!!         ! Update the mask one last time (probably not needed because li_diagnostic solve is coming up)
!!!         call li_calculate_mask(mesh, state, err)

!!!      !====================================
!!!      !case ('eigencalving')

!!!      end select

   !--------------------------------------------------------------------
   end subroutine li_apply_calving



!***********************************************************************
!***********************************************************************
! Private subroutines:
!***********************************************************************
!***********************************************************************



!***********************************************************************
!
!  subroutine tend_layerThickness_fo_upwind
!
!> \brief   Computes tendency term from horizontal advection of thickness layers
!> \author  Matt Hoffman
!> \date    16 April 2012
!> \version SVN:$Id$
!> \details 
!>  This routine computes the horizontal advection tendency for each
!>  thickness layer using first-order upwinding. Based on
!>  ocn_thick_hadv_tend in the ocean core.  This is an alternative to lice_tend_h
!>  that calculates the tendency for each layer, which would then need to be
!>  added up to calculate the change in thickness.  The two methods yield identical
!>  results.  
!
!-----------------------------------------------------------------------
   subroutine tend_layerThickness_fo_upwind(mesh, normalVelocity, layerThicknessEdge, edgeMask, tend, dt, MinOfMaxAllowableDt, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         mesh          !< Input: mesh information

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

      integer :: nCells, nVertLevels, iEdge, iCell, i, k
      integer, dimension(:,:), pointer :: edgesOnCell, edgeSignOnCell
      integer, dimension(:), pointer :: nEdgesOnCell
      real (kind=RKIND) :: invAreaCell, flux, maxAllowableDt
      real (kind=RKIND), dimension(:), pointer :: dvEdge, areaCell, dcEdge

      ! Only needed for optional check for mass conservation
      !real (kind=RKIND) :: tendVolSum

      err = 0

      nCells = mesh % nCells
      nVertLevels = mesh % nVertLevels
      nEdgesOnCell => mesh % nEdgesOnCell % array
      edgesOnCell => mesh % edgesOnCell % array
      edgeSignOnCell => mesh % edgeSignOnCell % array
      dvEdge => mesh % dvEdge % array
      dcEdge => mesh % dcEdge % array
      areaCell => mesh % areaCell % array

      MinOfMaxAllowableDt = 1.0e36_RKIND

      do iCell = 1, nCells
        invAreaCell = 1.0 / areaCell(iCell)
        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)
          do k = 1, nVertLevels

            if (abs(normalVelocity(k, iEdge)) > 0.0_RKIND) then
               maxAllowableDt = (0.5_RKIND * dcEdge(iEdge)) / abs(normalVelocity(k, iEdge))  ! in years 
            else
               maxAllowableDt = 1.0e36_RKIND
            endif
            if ( maxAllowableDt < dt ) then
                 !write(0,*) 'CFL violation at level, edge', k, iEdge                      
                 err = err + 1
            endif
            MinOfMaxAllowableDt = min(MinOfMaxAllowableDt, maxAllowableDt)

            flux = normalVelocity(k, iEdge) * dvEdge(iEdge) * layerThicknessEdge(k, iEdge)
            tend(k, iCell) = tend(k, iCell) + edgeSignOnCell(i, iCell) * flux * invAreaCell
          end do
        end do
      end do

      if (err > 0) then
         write(0,*) 'CFL violation on this processor on ', err, ' level-edges!  Maximum allowable time step (seconds) for this processor is ', MinOfMaxAllowableDt
         err = 1
      endif

      if (config_print_thickness_advection_info .eqv. .true.) then
         write(6,*) '  Maximum allowable time step (seconds) on this processor is ~', MinOfMaxAllowableDt
      endif

      ! Optional check for mass conservation
      !tendVolSum = 0.0_RKIND
      !do iCell=1, mesh % nCells
      !     tendVolSum  = tendVolSum + sum(tend(:,iCell)) * areaCell(iCell)
      !end do
      !print *,'SUM OF VOLUME TENDENCY  ======', tendVolSum

      ! === error check
      if (err > 0) then
          write (0,*) "An error has occurred in tend_layerThickness_fo_upwind."
      endif

   !--------------------------------------------------------------------

   end subroutine tend_layerThickness_fo_upwind!}}}



!***********************************************************************
!
!  subroutine tracer_advection_tend_fo
!
!> \brief   tracer_advection_tend_fo
!> \author  Matt Hoffman
!> \date    06 February 2013
!> \version SVN:$Id$
!> \details 
!>  

!
!-----------------------------------------------------------------------
   subroutine tracer_advection_tend_fo(tracers, uh, mesh, tracer_tendency, err)
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         mesh          !< Input: mesh information

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         tracers    !< Input: 

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         uh    !< Input: 

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tracer_tendency   !< Input: layer tracer tendencies

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
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell, cellsOnEdge
      real (kind=RKIND), dimension(:), pointer :: areaCell, dvEdge
      ! counters, mesh variables, index variables
      integer :: nTracers, nCells, nVertLevels
      integer :: iTracer, iCell, iEdge, thisEdge, k
      ! stuff for making calculations
      real (kind=RKIND) :: q, flux

! Commenting out this entire subroutine for now.
!!!      nCells = mesh % nCells
!!!      nVertLevels = mesh % nVertLevels
!!!      nTracers = size(tracers, 1)

!!!      nEdgesOnCell => mesh % nEdgesOnCell % array
!!!      edgesOnCell => mesh % edgesOnCell % array
!!!      cellsOnEdge => mesh % cellsOnEdge % array
!!!      areaCell => mesh % areaCell % array
!!!      dvEdge => mesh % dvEdge % array

!!!      do iTracer = 1, nTracers
!!!      do iCell = 1, nCells
!!!         do iEdge = 1, nEdgesOnCell(iCell)
!!!            ! What is this edge's index?
!!!            thisEdge = edgesOnCell(iEdge, iCell)
!!!            do k = 1, nVertLevels
!!!               ! Find q-edge using in an upwind sense
!!!               if (cellsOnEdge(1, thisEdge) == iCell ) then
!!!                  if ( uh(k,thisEdge) < 0.0_RKIND) then
!!!                     ! iCell is the downwind cell so use the other cell for q
!!!                     q = tracers(iTracer, k, cellsOnEdge(2, thisEdge) )
!!!                     flux = -1.0_RKIND * dvEdge(thisEdge) * uh(k, thisEdge) * q
!!!                  else
!!!                     ! iCell is up the upwind cell or the flux is 0
!!!                     q = tracers(iTracer, k, iCell)
!!!                     flux = -1.0_RKIND * dvEdge(thisEdge) * uh(k, thisEdge) * q
!!!                  end if
!!!               else   !  if (cellsOnEdge(2, thisEdge) == iCell ) then
!!!                  if ( uh(k,thisEdge) > 0.0_RKIND) then
!!!                     ! iCell is the downwind cell so use the other cell for q
!!!                     q = tracers(iTracer, k, cellsOnEdge(1, thisEdge) )
!!!                     flux = dvEdge(thisEdge) * uh(k, thisEdge) * q
!!!                  else
!!!                     ! iCell is up the upwind cell or the flux is 0
!!!                     q = tracers(iTracer, k, iCell)
!!!                     flux = dvEdge(thisEdge) * uh(k, thisEdge) * q
!!!                  end if
!!!               end if

!!!               tracer_tendency(iTracer, k, iCell) = tracer_tendency(iTracer, k, iCell) + flux / areaCell(iCell) 
!!!            end do ! vert levs
!!!         end do ! edges
!!!      end do ! cells
!!!      end do ! tracers

   end subroutine tracer_advection_tend_fo



! This subroutine calculates the tendency for the entire thickness column rather than layer-by-layer.
! Tracer advection requires knowing thickness tendencies layer-by-layer, so this has been
! deprecated.
!!***********************************************************************
!!
!!  subroutine li_tend_h_fo_upwind
!!
!!> \brief   Computes tendency term from horizontal advection of thickness
!!> \author  Matt Hoffman
!!> \date    16 April 2012
!!> \version SVN:$Id$
!!> \details 
!!>  This routine computes the horizontal advection tendency for
!!>  thickness based on current state and user choices of forcings. Based on
!!>  ocn_thick_hadv_tend in the ocean core.
!!
!!-----------------------------------------------------------------------

!   subroutine li_tend_h_fo_upwind(mesh, normalVelocity, layerThickness, thickness, thickness_tend, dt, err)!{{{

!      !-----------------------------------------------------------------
!      !
!      ! input variables
!      !
!      !-----------------------------------------------------------------

!      real (kind=RKIND), dimension(:,:), intent(in) :: &
!         normalVelocity    !< Input: velocity

!      real (kind=RKIND), dimension(:,:), intent(in) :: &
!         layerThickness    !< Input: thickness of each layer

!      real (kind=RKIND), dimension(:), intent(in) :: &
!         thickness    !< Input: thickness

!      type (mesh_type), intent(in) :: &
!         mesh          !< Input: mesh information

!      real (kind=RKIND) :: dt

!      !-----------------------------------------------------------------
!      !
!      ! input/output variables
!      !
!      !-----------------------------------------------------------------

!      real (kind=RKIND), dimension(:), intent(inout) :: &
!         thickness_tend        !< Input/Output: thickness tendency

!      !-----------------------------------------------------------------
!      !
!      ! output variables
!      !
!      !-----------------------------------------------------------------

!      integer, intent(out) :: err !< Output: error flag

!      !-----------------------------------------------------------------
!      !
!      ! local variables
!      !
!      !-----------------------------------------------------------------

!      integer :: iEdge, nEdges, cell1, cell2, nVertLevels, k, CellUpwind, CellDownwind

!      integer, dimension(:,:), pointer :: cellsOnEdge

!      real (kind=RKIND) :: flux, VelSign, h_edge, ubar, maxAllowableDt
!      real (kind=RKIND), dimension(:), pointer :: dvEdge, areaCell, dcEdge

!      err = 0
!      maxAllowableDt = 1.0e36

!      nEdges = mesh % nEdges
!      nVertLevels = mesh % nVertLevels
!      cellsOnEdge => mesh % cellsOnEdge % array
!      dvEdge => mesh % dvEdge % array
!      dcEdge => mesh % dcEdge % array
!      areaCell => mesh % areaCell % array

!      !print *,'Max velocity magn.:', maxval(abs(normalVelocity))

!      ! Zero the tendency before accumulating
!      thickness_tend = 0.0

!      do iEdge=1,nEdges
!         cell1 = cellsOnEdge(1,iEdge)
!         cell2 = cellsOnEdge(2,iEdge)

!         VelSign = sign(1.0, normalVelocity(1, iEdge))
!         if (VelSign .gt. 0.0) then
!            CellUpwind = cell1
!            CellDownwind = cell2
!         else
!            CellUpwind = cell2
!            CellDownwind = cell1
!         endif
!         if (thickness(cellUpwind) .gt. 0.0) then  ! Don't calculate for non-ice cells - would result in divide by 0
!                 h_edge = thickness(CellUpwind)
!                 flux = 0.0
!                 do k=1, nVertLevels
!                    ! Calculate thickness averaged velocity  -  this make be calculated externally and passed in
!                    flux = flux + layerThickness(k, cellUpwind) * abs(normalVelocity(k, iEdge))
!                 enddo
!                 ubar = flux / thickness(cellUpwind)
!                 if ( (abs(ubar) * dt/SecondsInYear) .gt. (0.5 * dcEdge(iEdge))) then
!                      !maxAllowableDt = min(maxAllowableDt, (0.5 * dcEdge)/ubar )
!                      write(0,*) 'CFL violation at edge', iEdge
!                      err = err + 1
!                 endif
!                 thickness_tend(cellUpwind) = thickness_tend(cellUpwind) - flux * dvEdge(iEdge) / areaCell(cellUpwind)
!                 thickness_tend(cellDownwind) = thickness_tend(cellDownwind) + flux * dvEdge(iEdge) / areaCell(cellDownwind)
!         endif
!      end do

!      if (err .gt. 0) then
!         write(6,*) 'CFL violation at ', err, ' edges!  Maximum time step should be ', maxAllowableDt
!         err = 1
!      endif
!   !--------------------------------------------------------------------

!   end subroutine li_tend_h_fo_upwind !}}}



end module li_tendency

