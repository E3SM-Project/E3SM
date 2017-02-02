! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_velocity_simple
!
!> \MPAS land-ice simple velocity driver
!> \author William Lipscomb
!> \date   October 2015
!> \details
!>  This module contains the routines for calculating simple velocity fields
!>  (e.g., uniform in x direction, radially symmetric).
!>
!
!-----------------------------------------------------------------------

module li_velocity_simple

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use li_mask
   use li_setup
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
   public :: li_velocity_simple_init, &
             li_velocity_simple_finalize, &
             li_velocity_simple_block_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------



!***********************************************************************

contains

!***********************************************************************
!
!  routine li_velocity_simple_init
!
!> \brief   Initializes simple velocity
!> \author  William Lipscomb
!> \date    October 2015
!> \details
!>  This routine initializes the simple velocity cases.
!
!-----------------------------------------------------------------------

   subroutine li_velocity_simple_init(domain, err)

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

   !--------------------------------------------------------------------

    end subroutine li_velocity_simple_init



!***********************************************************************
!
!  routine li_velocity_simple_block_init
!
!> \brief   Initializes blocks for simple velocity
!> \author  William Lipscomb
!> \date    October 2015
!> \details
!>  This routine initializes each block with a simple velocity field
!>  (uniform velocity in a straight line, or radially symmetric).
!>  NOTE: This subroutine assumes flow in a plane with all z components = 0.
!
!-----------------------------------------------------------------------

   subroutine li_velocity_simple_block_init(block, err)

     use mpas_vector_reconstruction

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
      type (block_type), intent(inout) :: &
           block          !< Input/Output: block object

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

      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: velocityPool
      type (mpas_pool_type), pointer :: scratchPool

      integer, pointer :: nCells, nEdges
      integer, pointer :: nCellsSolve, nEdgesSolve
      integer, pointer :: nVertInterfaces
      integer, pointer :: config_stats_cell_ID

      character(len=StrKind), pointer :: config_simple_velocity_type

      !NOTE: Assume a planar mesh, so z coordinates are not needed
      real (kind=RKIND), dimension(:), pointer :: &
           xCell, yCell,   &       ! cell center coordinates
           xEdge, yEdge,   &       ! edge midpoint coordinates
           dcEdge                  ! distance between the 2 cell centers on each side of an edge

      ! prescribed velocity at cell centers
      type (field1dReal), pointer :: uVelocityXField
      type (field1dReal), pointer :: uVelocityYField
      real (kind=RKIND), dimension(:), pointer :: uVelocityX
      real (kind=RKIND), dimension(:), pointer :: uVelocityY

      real (kind=RKIND), dimension(:,:), pointer :: &
           normalVelocityInitial,                         & ! normal component of velocity on edges
           uReconstructX, uReconstructY, uReconstructZ,   & ! x/y/z velocity components at cell center
           uReconstructZonal, uReconstructMeridional        ! zonal and meridional velocity components at cell center

      integer, dimension(:,:), pointer :: cellsOnEdge     ! indices for the 2 cells on each edge

      real (kind=RKIND), dimension(2) :: unitNormalVector   ! x/y components of normal vector on an edge

      real (kind=RKIND) :: magnitude, radius, speed, xDiff, yDiff

      real (kind=RKIND) :: uEdgeX, uEdgeY   ! x/y components of velocity at edge midpoints

      integer :: err_tmp

      integer :: iLevel, iEdge, iCell, iCell1, iCell2

      real (kind=RKIND), parameter :: flowSpeed = 1000._RKIND/scyr  ! flow speed (m/s)
                                                                    ! applies to uniform flow
      real (kind=RKIND), parameter :: flowTheta = 0.0_RKIND         ! direction of flow (0 < theta < 2*pi)
                                                                    ! applied to uniform flow (not radial)

      !Note: For radial flow, the user may want to reset these parameters
      real (kind=RKIND), parameter :: flowGradient = 1.2e-3_RKIND/scyr   ! du/dr for radial flow

      real (kind=RKIND), parameter :: xCenter = 0.0_RKIND   ! x coordinate of center of radial flow
      real (kind=RKIND), parameter :: yCenter = 0.0_RKIND   ! y coordinate of center of radial flow

      !WHL - debug diagnostics only
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell  ! index for each edge on a cell
      real(kind=RKIND), dimension(:,:,:), pointer ::  &
           coeffsReconstruct        ! coefficients for reconstructing edge-based fields at cell centers
      integer :: iEdgeOnCell

      logical, parameter :: velocity_simple_bug_check = .false.

      !--------------------------------
      !WHL - optional thickness, SMB and topography tweaking for the radial velocity field and circular-shelf test case
      !TODO - Remove these options after testing

      character(len=StrKIND), pointer :: config_calving

      logical, parameter :: radialMelting = .false.
!!      logical, parameter :: radialMelting = .true.
      type (mpas_pool_type), pointer :: geometryPool
      real (kind=RKIND), dimension(:), pointer :: thickness
      real (kind=RKIND), dimension(:), pointer :: sfcMassBal
      real (kind=RKIND), dimension(:), pointer :: bedTopography

      real (kind=RKIND), parameter :: &
           maxRadius = 21000.0_RKIND   ! ice radius (m) for circular shelf problem

      real (kind=RKIND), parameter :: &
           maxMelt = 100.0_RKIND * 910.0_RKIND / scyr   ! max melt rate, kg/m2/s (converted from 100 m/yr)
                                                        ! for radial melting option
      real (kind=RKIND), parameter :: &
           spikeTopography = -880.0_RKIND    ! elevation of spike that grounds the ice
                                             ! for config_calving = 'floating'

      integer, parameter :: ncellsPerRow = 40
      integer, parameter :: nRows = 46
      integer :: i, iRow
      !--------------------------------

      ! No block init needed.
      err = 0
      err_tmp = 0

      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

      ! Set needed variables and pointers

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_dimension(meshPool, 'nVertInterfaces', nVertInterfaces)

      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
      call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)

      call mpas_pool_get_array(velocityPool, 'normalVelocityInitial', normalVelocityInitial)
      call mpas_pool_get_array(velocityPool, 'uReconstructX', uReconstructX)
      call mpas_pool_get_array(velocityPool, 'uReconstructY', uReconstructY)
      call mpas_pool_get_array(velocityPool, 'uReconstructZ', uReconstructZ)
      call mpas_pool_get_array(velocityPool, 'uReconstructZonal', uReconstructZonal)
      call mpas_pool_get_array(velocityPool, 'uReconstructMeridional', uReconstructMeridional)

      call mpas_pool_get_field(scratchPool, 'workCell', uVelocityXField)
      call mpas_allocate_scratch_field(uVelocityXField, .true.)
      uVelocityX => uVelocityXField % array

      call mpas_pool_get_field(scratchPool, 'workCell2', uVelocityYField)
      call mpas_allocate_scratch_field(uVelocityYField, .true.)
      uVelocityY => uVelocityYField % array

      call mpas_pool_get_config(liConfigs, 'config_stats_cell_ID', config_stats_cell_ID)
      call mpas_pool_get_config(liConfigs, 'config_simple_velocity_type', config_simple_velocity_type)

      !WHL - debug diagnostics only
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'coeffs_reconstruct', coeffsReconstruct)

      uVelocityX(:) = 0.0_RKIND
      uVelocityY(:) = 0.0_RKIND

      ! prescribe the x and y velocity components at cell centers (with no vertical variation)

      if (trim(config_simple_velocity_type) == 'uniform') then

         uVelocityX(:) = flowSpeed * cos(flowTheta)
         uVelocityY(:) = flowSpeed * sin(flowTheta)

      elseif (trim(config_simple_velocity_type) == 'radial') then

         do iCell = 1, nCells
            xDiff = xCell(iCell) - xCenter
            yDiff = yCell(iCell) - yCenter
            radius = sqrt(xDiff**2 + yDiff**2)
            if (radius > 0.0_RKIND) then
               speed = flowGradient * radius
               uVelocityX(iCell) = speed * xDiff/radius
               uVelocityY(iCell) = speed * yDiff/radius
            else
               uVelocityX(iCell) = 0.0_RKIND
               uVelocityY(iCell) = 0.0_RKIND
            endif
         enddo

      endif

      ! given the velocity components at cell centers, compute the normal velocity component on edges

      normalVelocityInitial(:,:) = 0.0_RKIND

      do iEdge = 1, nEdgesSolve

         iLevel = 1
         iCell1 = cellsOnEdge(1,iEdge)
         iCell2 = cellsOnEdge(2,iEdge)

         ! average the velocity from the neighboring cells to the edge
         uEdgeX = 0.5_RKIND * (uVelocityX(iCell1) + uVelocityX(iCell2))
         uEdgeY = 0.5_RKIND * (uVelocityY(iCell1) + uVelocityY(iCell2))

         ! Compute the components of the normal vector on the edge

         unitNormalVector(1) = xEdge(iEdge) - xCell(iCell1)
         unitNormalVector(2) = yEdge(iEdge) - yCell(iCell1)
         magnitude = sqrt(unitNormalVector(1)**2 + unitNormalVector(2)**2)

         ! Note: The magnitude should be dcEdge/2.
         !       But this may not be the case for edges at the border of a periodic domain;
         !        for these cells, the magnitude may be comparable to the domain size.
         !       For such edges, create the normal vector from iCell2 instead.
         if (magnitude > dcEdge(iEdge)) then
!            write(stderrUnit,*) 'Use iCell2 instead: iEdge, iCell1, iCell2, magnitude, dcEdge =', &
!                iEdge, iCell1, iCell2, magnitude, dcEdge(iEdge)
            unitNormalVector(1) = -(xEdge(iEdge) - xCell(iCell2))
            unitNormalVector(2) = -(yEdge(iEdge) - yCell(iCell2))
            magnitude = sqrt(unitNormalVector(1)**2 + unitNormalVector(2)**2)
         endif

         unitNormalVector(:) = unitNormalVector(:)/magnitude

         ! Compute the dot product of the velocity with the normal vector
         ! Set to the same value everywhere in the column
         normalVelocityInitial(:,iEdge) = uEdgeX*unitNormalVector(1) + uEdgeY*unitNormalVector(2)

      enddo   ! iEdge

      !WHL - debug
      iCell = config_stats_cell_ID
      iLevel = 1
      write(stderrUnit,*) ' '
      write(stderrUnit,*) 'Prescribed velocity, iCell, uvel, vvel (m/yr):', iCell, uVelocityX(iCell)*scyr, uVelocityY(iCell)*scyr
      write(stderrUnit,*) 'xCell, yCell:', xCell(iCell), yCell(iCell)
      write(stderrUnit,*) 'iEdgeOnCell, cellsOnEdge, normalVelocity:'
      do iEdgeOnCell = 1, nEdgesOnCell(iCell)
         iEdge = edgesOnCell(iEdgeOnCell,iCell)
         write(stderrUnit,*) iEdgeOnCell, cellsOnEdge(:,iEdge), normalVelocityInitial(iLevel,iEdge)*scyr
      enddo

      if (velocity_simple_bug_check) then

         ! Make sure we can recover the cell-center velocity to a good approximation
         call mpas_reconstruct(meshPool, normalVelocityInitial,               &
                               uReconstructX, uReconstructY, uReconstructZ, &
                               uReconstructZonal, uReconstructMeridional )

         ! Loop over cells, comparing the reconstructed velocity to the prescribed velocity
         ! Note: Currently, the reconstruction coefficients are not correct for cells at the edge of a periodic domain,
         !        so errors will be generated even though the normal velocities are correct.
         !       For this reason I have commented out the warning messages.

         do iCell = 1, nCellsSolve

            speed = sqrt(uVelocityX(iCell)**2 + uVelocityY(iCell)**2)

            if (iCell == config_stats_cell_ID) then
               iLevel = 1
               write(stderrUnit,*) ' '
               write(stderrUnit,*) 'Velocity reconstruction, iCell =', iCell
               write(stderrUnit,*) 'Initial velocity:', uVelocityX(iCell), uVelocityY(iCell)
               write(stderrUnit,*) 'Reconstructed velocity:', uReconstructX(iLevel,iCell), uReconstructY(iLevel,iCell)
               write(stderrUnit,*) 'Reconstruction coefficients:'
               write(stderrUnit,*) ' '
               do iEdgeOnCell = 1, nEdgesOnCell(iCell)
                  write(stderrUnit,*) iEdgeOnCell, coeffsReconstruct(:,iEdgeOnCell,iCell)
               enddo
               write(stderrUnit,*) ' '
               write(stderrUnit,*) 'iEdgeOnCell, cellsOnEdge, normalVelocity:'
               do iEdgeOnCell = 1, nEdgesOnCell(iCell)
                  iEdge = edgesOnCell(iEdgeOnCell,iCell)
                  write(stderrUnit,*) iEdgeOnCell, cellsOnEdge(:,iEdge), normalVelocityInitial(iLevel,iEdge)
               enddo
            endif

            iLevel = 1   ! Check at one level only, since the velocity is vertically uniform

            if (abs(uReconstructX(iLevel,iCell) - uVelocityX(iCell)) > 1.e-8_RKIND*speed .or.  &
                abs(uReconstructY(iLevel,iCell) - uVelocityY(iCell)) > 1.e-8_RKIND*speed) then

               xDiff = abs(uReconstructX(iLevel,iCell) - uVelocityX(iCell)) / speed
               yDiff = abs(uReconstructY(iLevel,iCell) - uVelocityY(iCell)) / speed

               write(stderrUnit,*) ' '
               write(stderrUnit,*) 'WARNING: Reconstructed velocity not equal to uniform velocity, iCell, xDiff, yDiff=', &
                  iCell, xDiff, yDiff
               write(stderrUnit,*) 'Prescribed velocity:   ', uVelocityX(iCell), uVelocityY(iCell)
               write(stderrUnit,*) 'Reconstructed velocity:', uReconstructX(iLevel,iCell), uReconstructY(iLevel,iCell)
               write(stderrUnit,*) 'iEdgeOnCell, cellsOnEdge, normal velocity:'
               do iEdgeOnCell = 1, nEdgesOnCell(iCell)
                  iEdge = edgesOnCell(iEdgeOnCell,iCell)
                  write(stderrUnit,*) iEdgeOnCell, cellsOnEdge(:,iEdge), normalVelocityInitial(iLevel,iEdge)
               enddo
               err = 1
            endif

         enddo   ! iCell

      endif   ! bug check


      !--------------------------------
      !TODO - Remove these options after testing the calving scheme

      call mpas_pool_get_config(liConfigs, 'config_calving', config_calving)

      if (radialMelting) then   ! force the calving front to retreat

         write(stderrUnit,*) 'Setting up radially symmetric melting'
         write(stderrUnit,*) 'Melt rate at periphery (m/yr) =', maxMelt * scyr / 910.0_RKIND

         ! Zero out the normal velocities, since we are testing ice retreat
         normalVelocityInitial(:,:) = 0.0_RKIND

         ! Set the thickness and melt rate
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(geometryPool, 'sfcMassBal', sfcMassBal)

         do iCell = 1, nCells
            xDiff = xCell(iCell) - xCenter
            yDiff = yCell(iCell) - yCenter
            radius = sqrt(xDiff**2 + yDiff**2)
            ! set thickness to taper away from the center
            thickness(iCell) = thickness(iCell) * (1.0_RKIND - radius/maxRadius)
            ! set melting to increase away from the center
            sfcMassBal(iCell) = maxMelt * (-radius/maxRadius)
         enddo

         write(stderrUnit,*) ' '
         write(stderrUnit,*) 'thickness (m):'
         do iRow = nRows, 1, -1
            if (mod(iRow,2) == 0) then  ! indent for even-numbered rows
               write(stderrUnit,'(a3)',advance='no') '    '
            endif
            do i = nCellsPerRow/2 - 2, nCellsPerRow
               iCell = (iRow-1)*nCellsPerRow + i
               write(stderrUnit,'(f8.2)',advance='no') thickness(iCell)
            enddo
            write(stderrUnit,*) ' '
         enddo

         write(stderrUnit,*) ' '
         write(stderrUnit,*) 'melt rate (m/yr):'
         do iRow = nRows, 1, -1
            if (mod(iRow,2) == 0) then  ! indent for even-numbered rows
               write(stderrUnit,'(a3)',advance='no') '    '
            endif
            do i = nCellsPerRow/2 - 2, nCellsPerRow
               iCell = (iRow-1)*nCellsPerRow + i
               write(stderrUnit,'(f8.2)',advance='no') -sfcMassBal(iCell)*scyr/910.0_RKIND
            enddo
            write(stderrUnit,*) ' '
         enddo

      endif   ! radialMelting

      if (trim(config_calving) == 'topographic_threshold') then

         write(stderrUnit,*) 'Setting topography to drop off at periphery'

         ! Set the bed topography to drop off near the periphery of the ice
         ! so as to check that the calving topographic threshold option is working.
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)

         do iCell = 1, nCells
            xDiff = xCell(iCell) - xCenter
            yDiff = yCell(iCell) - yCenter
            radius = sqrt(xDiff**2 + yDiff**2)
            if (radius > 0.9_RKIND*maxRadius) then   ! close to the edge
               bedTopography(iCell) = bedTopography(iCell) * radius/(0.9_RKIND*maxRadius)
            endif
         enddo

         write(stderrUnit,*) ' '
         write(stderrUnit,*) 'bedTopography (m):'
         do iRow = nRows, 1, -1
            if (mod(iRow,2) == 0) then  ! indent for even-numbered rows
               write(stderrUnit,'(a3)',advance='no') '    '
            endif
            do i = nCellsPerRow/2 - 2, nCellsPerRow
               iCell = (iRow-1)*nCellsPerRow + i
               write(stderrUnit,'(f8.2)',advance='no') -bedTopography(iCell)
            enddo
            write(stderrUnit,*) ' '
         enddo

      elseif (trim(config_calving) == 'floating') then

         write(stderrUnit,*) 'Setting topography to be mostly grounded'
         write(stderrUnit,*) 'Spike depth =', spikeTopography

         ! Put in a large spike that grounds most of the ice, but leaves the peripheral ice floating,
         ! so as to check that the calving no-float option is working.
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)

         do iCell = 1, nCells
            xDiff = xCell(iCell) - xCenter
            yDiff = yCell(iCell) - yCenter
            radius = sqrt(xDiff**2 + yDiff**2)
            if (radius < 0.9_RKIND*maxRadius) then   ! inner part of ice shelf
               bedTopography(iCell) = spikeTopography
            endif
         enddo

         write(stderrUnit,*) ' '
         write(stderrUnit,*) 'bed depth (m):'
         do iRow = nRows, 1, -1
            if (mod(iRow,2) == 0) then  ! indent for even-numbered rows
               write(stderrUnit,'(a3)',advance='no') '    '
            endif
            do i = nCellsPerRow/2 - 2, nCellsPerRow
               iCell = (iRow-1)*nCellsPerRow + i
               write(stderrUnit,'(f8.2)',advance='no') -bedTopography(iCell)
            enddo
            write(stderrUnit,*) ' '
         enddo

      endif   ! config_calving
      !--------------------------------

      ! clean up
      call mpas_allocate_scratch_field(uVelocityXField, .true.)
      call mpas_allocate_scratch_field(uVelocityYField, .true.)

      ! === error check
      if (err > 0) then
          write (stderrUnit,*) "An error has occurred in li_velocity_uniform_init."
      endif

   !--------------------------------------------------------------------
    end subroutine li_velocity_simple_block_init


!***********************************************************************
!
!  routine li_velocity_simple_finalize
!
!> \brief   finalizes simple velocity
!> \author  William Lipscomb
!> \date    October 2015
!> \details
!>  This routine finalizes the simple velocity cases.
!
!-----------------------------------------------------------------------

   subroutine li_velocity_simple_finalize(err)

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


   !--------------------------------------------------------------------

    end subroutine li_velocity_simple_finalize



   ! private subroutines




!***********************************************************************

  end module li_velocity_simple

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
