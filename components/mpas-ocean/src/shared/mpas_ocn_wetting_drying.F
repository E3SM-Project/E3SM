! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_wetting_drying
!
!> \brief MPAS ocean wetting and drying
!> \author Phillip J. Wolfram
!> \date   03/09/2018
!> \details
!>  This module contains the wetting and drying algorithms
!>  to prevent cell thicknesses from becoming too small.
!>
!
!-----------------------------------------------------------------------

module ocn_wetting_drying

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use mpas_io_units

   use ocn_constants
   use ocn_config
   use ocn_diagnostics
   use ocn_diagnostics_variables
   use ocn_gm

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

   public :: ocn_wetting_drying_verify, ocn_prevent_drying_rk4

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------
   real(kind=RKIND), parameter :: eps = 1.0e-12_RKIND


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_wetting_drying_verify
!
!> \brief   Verifies that cells are not too dry.
!> \author  Phillip J. Wolfram
!> \date    03/09/2018
!> \details
!>  This routine checks that the minimum thickness in a cell is not
!>  too small.
!
!-----------------------------------------------------------------------

   subroutine ocn_wetting_drying_verify( block , minHeight, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (block_type), intent(in) :: block !< block for computation
      real (kind=RKIND), intent(in) :: minHeight

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

      type (mpas_pool_type), pointer :: statePool, meshPool, tendPool
      integer, dimension(:), pointer :: minLevelCell, maxLevelCell
      real (kind=RKIND), dimension(:), pointer :: sshSubcycleNew
      real (kind=RKIND), dimension(:), pointer :: bottomDepth
      integer, pointer :: nCellsSolve
      integer :: iCell, k
      integer :: debugUnit
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessCur
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessNew
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessTend
      real (kind=RKIND), dimension(:), pointer :: lonCell, latCell
      real (kind=RKIND) :: minThickness, layerThick
      character (len=StrKIND) :: debugFilename

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing coefficients
      !
      !-----------------------------------------------------------------

      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'tend', tendPool)

      call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, timeLevel=1)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, timeLevel=2)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)
      call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleNew, 2)
      call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      call mpas_pool_get_array(meshPool, 'latCell', latCell)

      err = 0

      call mpas_log_write( 'Verifying that cells are not dry...')

      ! check to make sure that there is no layer that is too dry
      minThickness = +1.0E34
      do iCell = 1, nCellsSolve
        do k = minLevelCell(iCell), maxLevelCell(iCell)
          ! use ssh as a proxy too for baroclinic mode
          if (trim(config_time_integrator) == 'split_explicit' .or. trim(config_time_integrator) == 'split_implicit') then
            layerThick = min(layerThicknessNew(k, iCell), (sshSubcycleNew(iCell)+bottomDepth(iCell))/maxLevelCell(iCell))
          else
            layerThick = layerThicknessNew(k, iCell)
          end if
          minThickness = min(minThickness, layerThick)
          if ( layerThick <  minHeight ) then
            ! report error if layerThickness constraint is violated
            if (err == 0) then
              debugFilename = ocn_build_log_filename('mpas_ocean_layerThickness_block_stats_', block % blockID)
              call mpas_new_unit(debugUnit)
              open(unit=debugUnit, file=debugFilename, form='formatted', status='unknown')
            end if
            err = iCell
            write(debugUnit, '(A, I5, A, I5, A, ES14.7, A, ES14.7, A, ES14.7, A)') &
              'ERROR: layerThickness too small at iCell=', iCell, ' k = ', k, &
              ', lon= ', lonCell(iCell), ', lat= ', latCell(iCell), &
              ' with thickness of ', layerThick , '.'
          end if
        end do
      end do

      call mpas_log_write('Minimum thickness is $r.', realArgs=(/minThickness/))

      if ( err > 0) then
        ! end the simulation if layerThickness constraint is violated
        call mpas_log_write( 'Done, some cells have dried.')
        flush(debugUnit)
        close(debugUnit)
        call mpas_release_unit(debugUnit)
        call mpas_log_write( 'ERROR: Layer thickness smaller than $r, see ' // debugFilename , &
                             realArgs=(/ minHeight /), messageType=MPAS_LOG_CRIT)
      else
        call mpas_log_write( 'Done verifying that cells are wet.')
      end if

   end subroutine ocn_wetting_drying_verify !}}}


!***********************************************************************
!
!  routine ocn_prevent_drying_rk4
!
!> \brief   Prevents velocity tendency from causing cell drying
!> \author  Phillip J. Wolfram
!> \date    03/20/2018
!> \details
!>  This routine modifies velocity tendency to prevent cell drying.
!
!-----------------------------------------------------------------------

   subroutine ocn_prevent_drying_rk4(block, dt, rkSubstepWeight, config_zero_drying_velocity, err) !{{{

     !-----------------------------------------------------------------
     !
     ! input variables
     !
     !-----------------------------------------------------------------

     type (block_type), intent(in) :: block
     real (kind=RKIND), intent(in) :: dt
     real (kind=RKIND), intent(in) :: rkSubstepWeight
     logical, pointer :: config_zero_drying_velocity

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


      type (mpas_pool_type), pointer :: tendPool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: provisStatePool
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessCur
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessProvis
      real (kind=RKIND), dimension(:, :), pointer :: normalVelocity

      integer, dimension(:), pointer :: minLevelEdgeTop
      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:), pointer :: maxLevelEdgeBot
      integer, pointer :: nEdges
      integer :: iEdge, k

      err = 0

     call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
     call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
     call mpas_pool_get_subpool(block % structs, 'state', statePool)
     call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)

     call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, 1)
     ! use thickness at n because constraint is h_n + dt*T_h > h_min
     call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
     call mpas_pool_get_array(provisStatePool, 'layerThickness', layerThicknessProvis, 1)

     call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
     call mpas_pool_get_array(meshPool, 'minLevelEdgeTop', minLevelEdgeTop)
     call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
     call mpas_pool_get_array(meshPool, 'maxLevelEdgeBot', maxLevelEdgeBot)

     !$omp parallel
     !$omp do schedule(runtime)
     do iEdge = 1, nEdges
       wettingVelocity(:, iEdge) = 0.0_RKIND
     end do
     !$omp end do
     !$omp end parallel

     ! ensure cells stay wet by selectively damping cells with a damping tendency to make sure tendency doesn't dry cells

     call ocn_wetting_drying_wettingVelocity(meshPool, layerThickEdgeFlux, layerThicknessCur, layerThicknessProvis, &
                                             normalTransportVelocity, rkSubstepWeight, wettingVelocity, err)

     ! prevent drying from happening with selective wettingVelocity
     !$omp parallel
     !$omp do schedule(runtime) private(k)
     do iEdge = 1, nEdges
       do k = minLevelEdgeTop(iEdge), maxLevelEdgeBot(iEdge)
         if (abs(normalTransportVelocity(k,iEdge) + wettingVelocity(k,iEdge)) < eps)  then
           ! prevent spurious flux for close to zero values
           normalTransportVelocity(k, iEdge) = 0.0_RKIND
           normalVelocity(k, iEdge) = 0.0_RKIND
         else if (abs(normalTransportVelocity(k,iEdge) + wettingVelocity(k,iEdge)) <= abs(normalTransportVelocity(k,iEdge))) then
           normalTransportVelocity(k, iEdge) = normalTransportVelocity(k, iEdge) + wettingVelocity(k, iEdge)
           normalVelocity(k, iEdge) = normalVelocity(k, iEdge) + wettingVelocity(k, iEdge)
         end if

         if (abs(wettingVelocity(k, iEdge)) > 0.0_RKIND .and. config_zero_drying_velocity) then
           normalTransportVelocity(k, iEdge) = 0.0_RKIND
           normalVelocity(k, iEdge) = 0.0_RKIND
         end if

       end do
     end do
     !$omp end do
     !$omp end parallel

   end subroutine ocn_prevent_drying_rk4 !}}}


!***********************************************************************
!
!  routine ocn_wetting_drying_wettingVelocity
!
!> \brief   Computes velocity to prevent cell drying
!> \author  Phillip J. Wolfram
!> \date    03/19/2018
!> \details
!>  This routine adds wetting velocity opposed to drying motion
!>  to prevent cells from drying.
!
!-----------------------------------------------------------------------
   subroutine ocn_wetting_drying_wettingVelocity(meshPool, layerThickEdgeFlux, layerThicknessCur, layerThicknessProvis, &

       normalVelocity, dt, wettingVelocity, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool           !< Input: horizonal mesh information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThicknessCur    !< Input: layer thickness at old time

       real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThicknessProvis    !< Input: provisional layer thickness

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickEdgeFlux  !< Input: flux-related layerThickness at an edge

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity     !< Input: transport

      real (kind=RKIND), intent(in) :: &
         dt     !< Input: time step

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         wettingVelocity          !< Input/Output: velocity wettingVelocityency

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

      integer :: iEdge, iCell, k, i
      integer, pointer :: nVertLevels, nCells
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:), pointer :: minLevelCell, maxLevelCell
      integer, dimension(:), pointer :: minLevelEdgeTop, maxLevelEdgeBot
      integer, dimension(:,:), pointer :: edgesOnCell
      integer, dimension(:,:), pointer :: edgeSignOnCell

      real (kind=RKIND) :: divOutFlux
      real (kind=RKIND) :: invAreaCell
      real (kind=RKIND) :: layerThickness
      real (kind=RKIND), dimension(:), pointer :: dvEdge
      real (kind=RKIND), dimension(:), pointer :: areaCell

      err = 0

      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
      call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'minLevelEdgeTop', minLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeBot', maxLevelEdgeBot)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      ! need predicted transport velocity to limit drying flux
      !$omp parallel
      !$omp do schedule(runtime) private(invAreaCell, i, iEdge, k, divOutFlux, layerThickness)
      do iCell = 1, nCells
        invAreaCell = 1.0_RKIND / areaCell(iCell)
        ! can switch with maxLevelEdgeBot(iEdge)
        do k = minLevelCell(iCell), maxLevelCell(iCell)
          divOutFlux = 0.0_RKIND
          layerThickness = min(layerThicknessProvis(k, iCell), layerThicknessCur(k, iCell))
          do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            if (k <= maxLevelEdgeBot(iEdge) .or. k >= minLevelEdgeTop(iEdge)) then
              ! only consider divergence flux leaving the cell
              if ( normalVelocity(k, iEdge) * edgeSignOnCell(i, iCell) < 0.0_RKIND ) then
                divOutFlux = divOutFlux + normalVelocity(k, iEdge) * edgeSignOnCell(i, iCell) &
                  * layerThickEdgeFlux(k, iEdge) * dvEdge(iEdge)  * invAreaCell
              end if
            end if
          end do

          ! if layer thickness is too small, limit divergence flux outwards with opposite velocity
          if ((layerThickness + dt*divOutFlux ) <= (config_drying_min_cell_height + config_drying_safety_height))  then
            ! limit divOutFlux out of cell to keep it wet
            divOutFlux = abs(divOutFlux)
            divOutFlux = (layerThickness - (config_drying_min_cell_height + eps)) / (dt*divOutFlux + eps)

            do i = 1, nEdgesOnCell(iCell)
              iEdge = edgesOnCell(i, iCell)
              if (k <= maxLevelEdgeBot(iEdge) .or. k >= minLevelEdgeTop(iEdge)) then
                if ( normalVelocity(k, iEdge) * edgeSignOnCell(i, iCell) <= 0.0_RKIND ) then
                  ! each outgoing velocity is penalized (but not the incoming, wetting velocities)
                  ! square the fractional term to make values near zero go to zero much quicker (to prevent threshold from being hit)
                  wettingVelocity(k, iEdge) =  - (min(max(0.0_RKIND, 1.0_RKIND - (divOutFlux*divOutFlux)), 1.0_RKIND)) * normalVelocity(k, iEdge)
                  ! just go with simple boolean approach for zero wetting velocity for debugging purposes
                  if (config_zero_drying_velocity) then
                    wettingVelocity(k, iEdge) = 1.0_RKIND
                  end if
                end if
              end if
            end do
          end if

        end do
      end do
      !$omp end do
      !$omp end parallel

   end subroutine ocn_wetting_drying_wettingVelocity !}}}


end module ocn_wetting_drying

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! vim: foldmethod=marker
