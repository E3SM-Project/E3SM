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
   use ocn_subgrid
   use ocn_vel_pressure_grad

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
   public :: ocn_wetting_velocity_factor_on_cell_edges

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

   subroutine ocn_prevent_drying_rk4(domain, block, dt, rkSubstepWeight, config_zero_drying_velocity, err) !{{{

     !-----------------------------------------------------------------
     !
     ! input variables
     !
     !-----------------------------------------------------------------

     type (domain_type), intent(inout) :: domain
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

     call ocn_wetting_drying_wettingVelocity(domain, layerThickEdgeFlux, layerThicknessCur, layerThicknessProvis, &
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

   subroutine ocn_wetting_drying_wettingVelocity(domain, layerThickEdgeFlux, layerThicknessCur, layerThicknessProvis, &
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

      type (domain_type), optional, intent(inout) :: &
         domain !< Input/Output: domain information, needed for halo exchange


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

      real (kind=RKIND) :: divOutFlux
      real (kind=RKIND) :: columnThickness
      real (kind=RKIND) :: hCrit, hRampMin, hRampMax, hEdgeTotal

      character (len=100) :: log_string

      integer:: cellDummy(2), cellCur, CellNei, sshCur, sshNei
      real (kind=RKIND), dimension(:), pointer :: sshCell
      real (kind=RKIND), dimension(:, :), pointer :: layerThicknessCell


      err = 0

      if (.not. config_zero_drying_velocity) return

      if (config_use_ssh_gradient_wetting_drying .and. &
          pGradType /= pGradTypeSSHgrad) then

         call mpas_log_write("config_use_ssh_gradient_wetting_drying requires " // &
                             "config_pressure_gradient_type = 'ssh_gradient'", MPAS_LOG_CRIT)
      endif

      hRampMin = config_zero_drying_velocity_ramp_hmin
      hRampMax = config_zero_drying_velocity_ramp_hmax
      hCrit = config_drying_min_cell_height + config_drying_safety_height

      layerThicknessCell => layerThicknessCellWetDry
      ! need predicted transport velocity to limit drying flux
      !$omp parallel
      !$omp do schedule(runtime) private(i, iEdge, k, divOutFlux)
      do iCell = 1, nCellsAll
        do k = minLevelCell(iCell), maxLevelCell(iCell)
          divOutFlux = 0.0_RKIND
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

          layerThicknessCell(k, iCell) = min(layerThicknessProvis(k, iCell), layerThicknessCur(k, iCell)) + dt * divOutFlux
          if ( .not. config_use_ssh_gradient_wetting_drying ) then
             call ocn_wetting_velocity_factor_on_cell_edges( &
                      wettingVelocityFactor, layerThicknessCell(k, iCell), normalVelocity, iCell, k)
          end if
        end do
      end do
      !$omp end do
      !$omp end parallel

      call mpas_dmpar_field_halo_exch(domain, 'layerThicknessCellWetDry')


      ! This wetting and drying implementation is based on the one above.
      ! However, it has special considerations for edges with normal velocities
      ! that have been previously zeroed to limit drying (or are zero initially),
      ! but have ssh gradients that would cause flow into the current cell,
      ! i.e., the neighboring cell's ssh is higher.
      ! There are two cases that require special treatment:
      !  1) In the case that the neighbor cell is predicted to be wet, flow should
      !     be allowed into the current cell.
      !  2) In the case that the neighor cell is predicted to be dry, flow should
      !     be prevented completely (and not ramped).
      ! Currently, this only works for single layer configurations with the 'ssh_gradient'
      ! pressure gradient option
      if (config_use_ssh_gradient_wetting_drying ) then
        ! Compute predicted ssh
        sshCell => sshCellWetDry
        sshCell = 0.0_RKIND

        wettingVelocityFactor = 0.0_RKIND

        if ( config_use_subgrid_wetting_drying ) then
          do iCell = 1, nCellsAll
            k = 1
            call ocn_subgrid_ssh_lookup(layerThicknessCell(k, iCell), &
                                        subgridWetVolumeCellTable(:,iCell), &
                                        subgridSshCellTableRange(:,iCell), &
                                        bottomDepth(iCell), &
                                        subgridCellBathymetryMin(iCell), &
                                        sshCell(iCell))
          enddo
        else
          !$omp parallel
          !$omp do schedule(runtime) private(k, columnThickness)
          do iCell = 1, nCellsAll
            columnThickness = 0.0_RKIND
            do k = minLevelCell(iCell), maxLevelCell(iCell)
              columnThickness = columnThickness + layerThicknessCell(k, iCell)
            enddo
            sshCell(iCell) = columnThickness - bottomDepth(iCell)
          enddo
          !$omp end do
          !$omp end parallel
        end if

        !$omp parallel
        !$omp do schedule(runtime) private(k, i, iEdge, cellDummy, cellNei, cellCur, sshCur, sshNei)
        do iCell = 1, nCellsAll
          do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = minLevelCell(iCell), maxLevelCell(iCell)
              if (k <= maxLevelEdgeTop(iEdge) .and. k >= minLevelEdgeBot(iEdge) .and. &
                 layerThicknessCell(k, iCell) <= hCrit ) then
                ! If edge velocity has been zeroed, check neighbor and currnt cell ssh.
                ! Allow flow into current cell if neighbor cell is wet and the ssh gradient
                ! would cause flow from neighbor cell into current cell.
                if ( normalVelocity(k,iEdge)*edgeSignOnCell(i, iCell) < 0.0_RKIND ) then

                  wettingVelocityFactor(k,iEdge) = 1.0_RKIND

                elseif ( normalVelocity(k,iEdge)*edgeSignOnCell(i, iCell) == 0.0_RKIND ) then

                  cellDummy(1:2) = cellsOnEdge(1:2,iEdge)
                  cellNei = merge( cellDummy(2), cellDummy(1), iCell == cellDummy(1) )
                  cellCur = iCell
                  if ( cellNei > nCellsAll ) cellNei = cellCur
                  sshCur = sshCell(cellCur)
                  sshNei = sshCell(cellNei)

                  wettingVelocityFactor(k,iEdge) = 1.0_RKIND

                  ! if the neigbor cell is anticipated wet
                  if ( layerThicknessCell(k, cellNei) >  hCrit ) then
                    if ( sshCur < sshNei ) then
                       wettingVelocityFactor(k,iEdge) = 0.0_RKIND
                    end if
                  end if

                end if ! velocity check
              endif ! k check and thickness check
            enddo ! k
          enddo ! i
        end do ! iCell
        !$omp end do
        !$omp end parallel

        !$omp parallel
        !$omp do schedule(runtime) private(k, i, iEdge, cellDummy, cellNei, cellCur, sshCur, sshNei)
        do iCell = 1, nCellsAll
          do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = minLevelCell(iCell), maxLevelCell(iCell)
              if (k <= maxLevelEdgeTop(iEdge) .and. k >= minLevelEdgeBot(iEdge) .and. &
                  config_zero_drying_velocity_ramp .and. &
                  (layerThicknessCell(k, iCell) > hCrit) .and. &
                  (layerThicknessCell(k, iCell) <= hRampMax) ) then

                if ( normalVelocity(k, iEdge) * edgeSignOnCell(i, iCell) <= 0.0_RKIND ) then
                  wettingVelocityFactor(k, iEdge) = 1.0_RKIND - &
                      tanh(50.0_RKIND * (layerThicknessCell(k, iCell) - hRampMin)/hRampMin)
                endif

                ! If edge velocity has been zeroed, check neighbor and current cell ssh.
                ! Prevent flow out of cell if neighbor cell is dry and the ssh gradient
                ! would cause flow from neighbor cell into current cell.
                if ( normalVelocity(k, iEdge) * edgeSignOnCell(i,iCell) == 0.0_RKIND ) then

                  cellDummy(1:2) = cellsOnEdge(1:2,iEdge)
                  cellNei = merge( cellDummy(2), cellDummy(1), iCell == cellDummy(1) )
                  cellCur = iCell
                  if ( cellNei > nCellsAll ) cellNei = cellCur
                  sshCur = sshCell(cellCur)
                  sshNei = sshCell(cellNei)

                  ! if the neigbor cell is anticipated dry
                  if ( layerThicknessCell(k, cellNei) <=  hCrit ) then
                     if ( sshCur < sshNei ) then
                         wettingVelocityFactor(k,iEdge) = 1.0_RKIND
                     end if
                  end if
                end if ! velocity check
              end if ! thickness check
            enddo ! k
          enddo ! iEdge
        end do ! iCell
        !$omp end do
        !$omp end parallel

      end if

   end subroutine ocn_wetting_drying_wettingVelocity !}}}

!***********************************************************************
!
!  routine ocn_wetting_velocity_factor_on_cell_edges
!
!> \brief   Computes wettingVelocityFactor at edges of iCell
!> \author  Giacomo Capodaglio
!> \date    09/06/2023
!> \details
!>  This routine computes wettingVelocityFactor at the edges of iCell
!
!-----------------------------------------------------------------------
   subroutine ocn_wetting_velocity_factor_on_cell_edges(wettingVelocityFactor, layerThicknessTmp, normalVelocity, iCell, k)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(in) :: &
         layerThicknessTmp    !< Input: layer thickness to use for limiting

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity       !< Input: transport

      integer, intent(in) :: &
         iCell                !< Input: wettingVelocityFactor is computed at the edges of iCell

      integer, intent(in) :: &
         k                    !< Input: wettingVelocityFactor is computed at layer k

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         wettingVelocityFactor          !< Input/Output: velocity wettingVelocityFactor

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: i, iEdge

      real (kind=RKIND) :: hCrit, hRampMin

      hCrit = config_drying_min_cell_height

      ! if layer thickness is too small, limit divergence flux outwards with
      ! opposite velocity
      if (layerThicknessTmp <= &
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
              (layerThicknessTmp > &
              hCrit + config_drying_safety_height) .and. &
              (layerThicknessTmp <= config_zero_drying_velocity_ramp_hmax)) then

         hRampMin = config_zero_drying_velocity_ramp_hmin
         ! Following O'Dea et al. (2021), if total upwinded wct is less than
         ! 2*critical thickness, apply damping at each edge
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            if (k <= maxLevelEdgeBot(iEdge) .and. k >= minLevelEdgeTop(iEdge)) then
               if ( normalVelocity(k, iEdge) * edgeSignOnCell(i, iCell) <= 0.0_RKIND ) then
                  wettingVelocityFactor(k, iEdge) = 1.0_RKIND - &
                  tanh(50.0_RKIND * (layerThicknessTmp - hRampMin)/hRampMin)
               end if
            end if
         end do

      end if

   end subroutine ocn_wetting_velocity_factor_on_cell_edges!}}}

end module ocn_wetting_drying

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! vim: foldmethod=marker
