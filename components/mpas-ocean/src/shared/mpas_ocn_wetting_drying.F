! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
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
   use ocn_mesh

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
      real (kind=RKIND), dimension(:), pointer :: sshSubcycleNew
      integer, pointer :: nCellsSolve
      integer :: iCell, k
      integer :: debugUnit
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessCur
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessNew
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessTend
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
      call mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)
      call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleNew, 2)

      err = 0

      call mpas_log_write( 'Verifying that cells are not dry...')

      ! check to make sure that there is no layer that is too dry
      minThickness = +1.0E34
      do iCell = 1, nCellsSolve
        do k = minLevelCell(iCell), maxLevelCell(iCell)
          ! use ssh as a proxy too for baroclinic mode
          ! Note: wetting-drying currently not supported for either of these time integration methods
          if (trim(config_time_integrator) == 'split_explicit' .or. &
              trim(config_time_integrator) == 'split_implicit' .or. &
              trim(config_time_integrator) == 'split_explicit_ab2') then
            layerThick = min(layerThicknessNew(k, iCell), &
                             (sshSubcycleNew(iCell) + bottomDepth(iCell))/maxLevelCell(iCell))
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
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: provisStatePool
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessCur
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessProvis
      real (kind=RKIND), dimension(:, :), pointer :: normalVelocity

      integer :: iEdge, k

      err = 0

     call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
     call mpas_pool_get_subpool(block % structs, 'state', statePool)
     call mpas_pool_get_subpool(block % structs, 'provis_state', provisStatePool)

     call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, 1)
     ! use thickness at n because constraint is h_n + dt*T_h > h_min
     call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
     call mpas_pool_get_array(provisStatePool, 'layerThickness', layerThicknessProvis, 1)


     !$omp parallel
     !$omp do schedule(runtime)
     do iEdge = 1, nEdgesAll
       wettingVelocityFactor(:, iEdge) = 0.0_RKIND
     end do
     !$omp end do
     !$omp end parallel

     ! ensure cells stay wet by selectively damping cells with a damping tendency to make 
     ! sure tendency doesn't dry cells

     call ocn_wetting_drying_wettingVelocity(layerThickEdgeFlux, layerThicknessCur, layerThicknessProvis, &
                                             normalTransportVelocity, rkSubstepWeight, wettingVelocityFactor, err)

     ! prevent drying from happening with selective wettingVelocityFactor
     if (config_zero_drying_velocity) then
       !$omp parallel
       !$omp do schedule(runtime) private(k)
       do iEdge = 1, nEdgesAll
         do k = minLevelEdgeTop(iEdge), maxLevelEdgeBot(iEdge)

           if (abs(wettingVelocityFactor(k, iEdge)) > 0.0_RKIND) then
             normalTransportVelocity(k, iEdge) = (1.0_RKIND - &
               wettingVelocityFactor(k, iEdge)) * normalTransportVelocity(k, iEdge)
             normalVelocity(k, iEdge) = (1.0_RKIND - &
               wettingVelocityFactor(k, iEdge)) * normalVelocity(k, iEdge)
           end if

         end do
       end do
       !$omp end do
       !$omp end parallel
     end if

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
   subroutine ocn_wetting_drying_wettingVelocity(layerThickEdgeFlux, layerThicknessCur, layerThicknessProvis, &

       normalVelocity, dt, wettingVelocityFactor, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

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
         wettingVelocityFactor          !< Input/Output: velocity wettingVelocityFactor

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

      integer :: cell1, cell2, iEdge, iCell, k, i

      real (kind=RKIND) :: divFlux, divOutFlux
      real (kind=RKIND) :: layerThickness
      real (kind=RKIND) :: hCrit, hRampMin, hEdgeTotal

      character (len=100) :: log_string

      err = 0

      hCrit = config_drying_min_cell_height

      if (.not. config_zero_drying_velocity) return

      ! need predicted transport velocity to limit drying flux
      !$omp parallel
      !$omp do schedule(runtime) private(i, iEdge, k, divOutFlux, layerThickness)
      do iCell = 1, nCellsAll
        do k = minLevelCell(iCell), maxLevelCell(iCell)
          divOutFlux = 0.0_RKIND
          layerThickness = min(layerThicknessProvis(k, iCell), layerThicknessCur(k, iCell))
          do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            if (k <= maxLevelEdgeTop(iEdge) .and. k >= minLevelEdgeBot(iEdge)) then
              ! only consider divergence flux leaving the cell
              if ( normalVelocity(k, iEdge) * edgeSignOnCell(i, iCell) < 0.0_RKIND ) then
                divOutFlux = divOutFlux + &
                             normalVelocity(k, iEdge) * edgeSignOnCell(i, iCell) * &
                             layerThickEdgeFlux(k, iEdge) * dvEdge(iEdge) * &
                             invAreaCell(iCell)
              end if
            end if
          end do
          layerThickness = layerThickness + dt * divOutFlux

          ! if layer thickness is too small, limit divergence flux outwards with
          ! opposite velocity
          if (layerThickness <= &
              hCrit + config_drying_safety_height) then
            do i = 1, nEdgesOnCell(iCell)
              iEdge = edgesOnCell(i, iCell)
              if (k <= maxLevelEdgeBot(iEdge) .and. k >= minLevelEdgeTop(iEdge)) then
                if ( normalVelocity(k, iEdge) * edgeSignOnCell(i, iCell) <= 0.0_RKIND ) then
                  wettingVelocityFactor(k, iEdge) = 1.0_RKIND
                end if
              end if
            end do
          elseif (config_zero_drying_velocity_ramp .and. &
                 (layerThickness > &
                 hCrit + config_drying_safety_height) .and. &
                 (layerThickness <= config_zero_drying_velocity_ramp_hmax)) then

            hRampMin = config_zero_drying_velocity_ramp_hmin
            ! Following O'Dea et al. (2021), if total upwinded wct is less than
            ! 2*critical thickness, apply damping at each edge
            do i = 1, nEdgesOnCell(iCell)
              iEdge = edgesOnCell(i, iCell)
              if (k <= maxLevelEdgeBot(iEdge) .and. k >= minLevelEdgeTop(iEdge)) then
                if ( normalVelocity(k, iEdge) * edgeSignOnCell(i, iCell) <= 0.0_RKIND ) then
                  wettingVelocityFactor(k, iEdge) = 1.0_RKIND - &
                    tanh(50.0_RKIND * (layerThickness - hRampMin)/hRampMin)
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
