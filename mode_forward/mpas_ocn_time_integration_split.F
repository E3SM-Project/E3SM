! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_split
!
!> \brief MPAS ocean split explicit time integration scheme
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains the routine for the split explicit
!>  time integration scheme
!
!-----------------------------------------------------------------------


module ocn_time_integration_split

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar
   use mpas_vector_reconstruction
   use mpas_spline_interpolation
   use mpas_timer
   use mpas_threading
   use mpas_timekeeping

   use ocn_tendency
   use ocn_diagnostics
   use ocn_gm

   use ocn_equation_of_state
   use ocn_vmix
   use ocn_time_average_coupled

   use ocn_effective_density_in_land_ice

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

   public :: ocn_time_integrator_split, ocn_time_integration_split_init

   character (len=*), parameter :: subcycleGroupName = 'subcycleFields'
   character (len=*), parameter :: finalBtrGroupName = 'finalBtrFields'
   integer :: nBtrSubcycles

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_split
!
!> \brief MPAS ocean split explicit time integration scheme
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This routine integrates a master time step (dt) using a
!>  split explicit time integrator.
!
!-----------------------------------------------------------------------

    subroutine ocn_time_integrator_split(domain, dt)!{{{
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Advance model state forward in time by the specified time step using
    !   Split_Explicit timestepping scheme
    !
    ! Input: domain - current model state in time level 1 (e.g., time_levs(1)state%h(:,:))
    !                 plus mesh meta-data
    ! Output: domain - upon exit, time level 2 (e.g., time_levs(2)%state%h(:,:)) contains
    !                  model state advanced forward in time by dt seconds
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt

      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: verticalMeshPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: tendPool
      type (mpas_pool_type), pointer :: tracersTendPool
      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: swForcingPool

      type (dm_info) :: dminfo
      integer :: iCell, i,k,j, iEdge, cell1, cell2, split_explicit_step, split, &
                 eoe, oldBtrSubcycleTime, newBtrSubcycleTime, uPerpTime, BtrCorIter, &
                 stage1_tend_time
      integer, dimension(:), allocatable :: n_bcl_iter
      type (block_type), pointer :: block
      real (kind=RKIND) :: normalThicknessFluxSum, thicknessSum, flux, sshEdge, hEdge1, &
                 CoriolisTerm, normalVelocityCorrection, temp, temp_h, coef, barotropicThicknessFlux_coeff, sshCell1, sshCell2
      integer :: useVelocityCorrection, err
      real (kind=RKIND), dimension(:,:), pointer :: &
                 vertViscTopOfEdge, vertDiffTopOfCell
      real (kind=RKIND), dimension(:,:,:), pointer :: tracersGroup
      real (kind=RKIND), dimension(:), allocatable:: uTemp
      real (kind=RKIND), dimension(:), pointer :: btrvel_temp
      type (field1DReal), pointer :: btrvel_tempField
      real (kind=RKIND), dimension(:,:), allocatable:: tracersTemp

      integer :: tsIter
      integer :: edgeHaloComputeCounter, cellHaloComputeCounter

      ! Config options
      character (len=StrKIND), pointer :: config_time_integrator
      integer, pointer :: config_n_bcl_iter_mid, config_n_bcl_iter_beg, config_n_bcl_iter_end
      integer, pointer :: config_n_ts_iter, config_btr_subcycle_loop_factor
      integer, pointer :: config_n_btr_cor_iter, config_num_halos
      logical, pointer :: config_use_standardGM

      logical, pointer :: config_use_freq_filtered_thickness, config_btr_solve_SSH2, config_filter_btr_mode
      logical, pointer :: config_vel_correction, config_prescribe_velocity, config_prescribe_thickness
      logical, pointer :: config_use_cvmix_kpp
      logical, pointer :: config_use_tracerGroup
      character (len=StrKIND), pointer :: config_land_ice_flux_mode

      real (kind=RKIND), pointer :: config_mom_del4, config_btr_gam1_velWt1, config_btr_gam2_SSHWt1
      real (kind=RKIND), pointer :: config_btr_gam3_velWt2

      ! Dimensions
      integer :: nCells, nEdges
      integer, pointer :: nCellsPtr, nEdgesPtr, nVertLevels, num_tracersGroup, startIndex, endIndex
      integer, pointer :: indexTemperature, indexSalinity
      integer, pointer :: indexSurfaceVelocityZonal, indexSurfaceVelocityMeridional
      integer, pointer :: indexSSHGradientZonal, indexSSHGradientMeridional
      integer, dimension(:), pointer :: nCellsArray, nEdgesArray

      ! Mesh array pointers
      integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeTop, nEdgesOnEdge, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, edgeMask, edgesOnEdge
      integer, dimension(:,:), pointer :: edgesOnCell, edgeSignOnCell

      real (kind=RKIND), dimension(:), pointer :: dcEdge, fEdge, bottomDepth, refBottomDepthTopOfCell
      real (kind=RKIND), dimension(:), pointer :: dvEdge, areaCell
      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge

      ! State Array Pointers
      real (kind=RKIND), dimension(:), pointer :: sshSubcycleCur, sshSubcycleNew
      real (kind=RKIND), dimension(:), pointer :: normalBarotropicVelocitySubcycleCur, normalBarotropicVelocitySubcycleNew
      real (kind=RKIND), dimension(:), pointer :: sshCur, sshNew
      real (kind=RKIND), dimension(:), pointer :: normalBarotropicVelocityCur, normalBarotropicVelocityNew
      real (kind=RKIND), dimension(:,:), pointer :: normalBaroclinicVelocityCur, normalBaroclinicVelocityNew
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocityCur, normalVelocityNew
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessCur, layerThicknessNew
      real (kind=RKIND), dimension(:,:), pointer :: highFreqThicknessCur, highFreqThicknessNew
      real (kind=RKIND), dimension(:,:), pointer :: lowFreqDivergenceCur, lowFreqDivergenceNew
      real (kind=RKIND), dimension(:,:,:), pointer :: tracersGroupCur, tracersGroupNew

      ! Tend Array Pointers
      real (kind=RKIND), dimension(:), pointer :: sshTend
      real (kind=RKIND), dimension(:,:), pointer :: highFreqThicknessTend
      real (kind=RKIND), dimension(:,:), pointer :: lowFreqDivergenceTend
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocityTend, layerThicknessTend
      real (kind=RKIND), dimension(:,:,:), pointer :: tracersGroupTend, activeTracersTend

      ! Diagnostics Array Pointers
      real (kind=RKIND), dimension(:), pointer :: barotropicForcing, barotropicThicknessFlux
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessEdge, normalTransportVelocity, normalGMBolusVelocity
      real (kind=RKIND), dimension(:,:), pointer :: vertAleTransportTop
      real (kind=RKIND), dimension(:,:), pointer :: velocityX, velocityY, velocityZ
      real (kind=RKIND), dimension(:,:), pointer :: velocityZonal, velocityMeridional
      real (kind=RKIND), dimension(:), pointer :: gradSSH
      real (kind=RKIND), dimension(:), pointer :: gradSSHX, gradSSHY, gradSSHZ
      real (kind=RKIND), dimension(:), pointer :: gradSSHZonal, gradSSHMeridional
      real (kind=RKIND), dimension(:,:), pointer :: surfaceVelocity, SSHGradient

      ! Diagnostics Field Pointers
      type (field2DReal), pointer :: normalizedRelativeVorticityEdgeField, divergenceField, relativeVorticityField
      type (field1DReal), pointer :: barotropicThicknessFluxField, boundaryLayerDepthField, effectiveDensityField

      ! State/Tend Field Pointers
      type (field1DReal), pointer :: normalBarotropicVelocitySubcycleField, sshSubcycleField
      type (field2DReal), pointer :: highFreqThicknessField, lowFreqDivergenceField
      type (field2DReal), pointer :: normalBaroclinicVelocityField, layerThicknessField
      type (field2DReal), pointer :: normalVelocityField
      type (field3DReal), pointer :: tracersGroupField

      ! tracer iterators
      type (mpas_pool_iterator_type) :: groupItr
      character (len=StrKIND) :: modifiedGroupName
      character (len=StrKIND) :: configName
      integer :: threadNum

      call mpas_timer_start("se timestep")

      call mpas_pool_get_config(domain % configs, 'config_n_bcl_iter_beg', config_n_bcl_iter_beg)
      call mpas_pool_get_config(domain % configs, 'config_n_bcl_iter_mid', config_n_bcl_iter_mid)
      call mpas_pool_get_config(domain % configs, 'config_n_bcl_iter_end', config_n_bcl_iter_end)
      call mpas_pool_get_config(domain % configs, 'config_n_ts_iter', config_n_ts_iter)
      call mpas_pool_get_config(domain % configs, 'config_btr_subcycle_loop_factor', config_btr_subcycle_loop_factor)
      call mpas_pool_get_config(domain % configs, 'config_btr_gam1_velWt1', config_btr_gam1_velWt1)
      call mpas_pool_get_config(domain % configs, 'config_btr_gam3_velWt2', config_btr_gam3_velWt2)
      call mpas_pool_get_config(domain % configs, 'config_btr_solve_SSH2', config_btr_solve_SSH2)
      call mpas_pool_get_config(domain % configs, 'config_n_btr_cor_iter', config_n_btr_cor_iter)
      call mpas_pool_get_config(domain % configs, 'config_btr_gam2_SSHWt1', config_btr_gam2_SSHWt1)
      call mpas_pool_get_config(domain % configs, 'config_filter_btr_mode', config_filter_btr_mode)

      call mpas_pool_get_config(domain % configs, 'config_mom_del4', config_mom_del4)
      call mpas_pool_get_config(domain % configs, 'config_use_freq_filtered_thickness', config_use_freq_filtered_thickness)
      call mpas_pool_get_config(domain % configs, 'config_time_integrator', config_time_integrator)
      call mpas_pool_get_config(domain % configs, 'config_vel_correction', config_vel_correction)

      call mpas_pool_get_config(domain % configs, 'config_prescribe_velocity', config_prescribe_velocity)
      call mpas_pool_get_config(domain % configs, 'config_prescribe_thickness', config_prescribe_thickness)

      call mpas_pool_get_config(domain % configs, 'config_prescribe_velocity', config_prescribe_velocity)
      call mpas_pool_get_config(domain % configs, 'config_prescribe_thickness', config_prescribe_thickness)

      call mpas_pool_get_config(domain % configs, 'config_use_standardGM', config_use_standardGM)
      call mpas_pool_get_config(domain % configs, 'config_use_cvmix_kpp', config_use_cvmix_kpp)
      call mpas_pool_get_config(domain % configs, 'config_land_ice_flux_mode', config_land_ice_flux_mode)

      call mpas_pool_get_config(domain % configs, 'config_num_halos', config_num_halos)

      allocate(n_bcl_iter(config_n_ts_iter))

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      !  Prep variables before first iteration
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call mpas_timer_start("se prep")
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCellsPtr)
         call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
         call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)

         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

         call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', normalBaroclinicVelocityCur, 1)
         call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityCur, 1)
         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)

         call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', normalBaroclinicVelocityNew, 2)
         call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityNew, 2)
         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)

         call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
         call mpas_pool_get_array(statePool, 'ssh', sshNew, 2)

         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)

         call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessCur, 1)
         call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessNew, 2)

         call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceCur, 1)
         call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceNew, 2)

         call mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', vertAleTransportTop)

         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
! mrp
         call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)

         nCells = nCellsPtr
         nEdges = nEdgesPtr

         ! Initialize * variables that are used to compute baroclinic tendencies below.

         !$omp do schedule(runtime) private(k)
         do iEdge = 1, nEdges
            do k = 1, nVertLevels !maxLevelEdgeTop % array(iEdge)

               ! The baroclinic velocity needs be recomputed at the beginning of a
               ! timestep because the implicit vertical mixing is conducted on the
               ! total u.  We keep normalBarotropicVelocity from the previous timestep.
               ! Note that normalBaroclinicVelocity may now include a barotropic component, because the
               ! weights layerThickness have changed.  That is OK, because the barotropicForcing variable
               ! subtracts out the barotropic component from the baroclinic.
               normalBaroclinicVelocityCur(k,iEdge) = normalVelocityCur(k,iEdge) - normalBarotropicVelocityCur(iEdge)

               normalVelocityNew(k,iEdge) = normalVelocityCur(k,iEdge)

               normalBaroclinicVelocityNew(k,iEdge) = normalBaroclinicVelocityCur(k,iEdge)
            end do
         end do
         !$omp end do

         !$omp do schedule(runtime) private(k)
         do iCell = 1, nCells
            sshNew(iCell) = sshCur(iCell)
            do k = 1, maxLevelCell(iCell)
               layerThicknessNew(k,iCell) = layerThicknessCur(k,iCell)
               ! set vertAleTransportTop to zero for stage 1 velocity tendency, first time through.
               vertAleTransportTop(k,iCell) = 0.0_RKIND
            end do
         end do
         !$omp end do


         threadnum = mpas_threading_get_thread_num()

         call mpas_pool_begin_iteration(tracersPool)
         do while ( mpas_pool_get_next_member(tracersPool, groupItr))
            if ( groupItr % memberType == MPAS_POOL_FIELD ) then
               call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupCur, 1)
               call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupNew, 2)

               if ( associated(tracersGroupCur) .and. associated(tracersGroupNew) ) then
                  !$omp do schedule(runtime) private(k)
                  do iCell = 1, nCells
                     do k = 1, maxLevelCell(iCell)
                        tracersGroupNew(:,k,iCell) = tracersGroupCur(:,k,iCell)
                     end do
                  end do
                  !$omp end do
               end if
            end if
         end do


         if (associated(highFreqThicknessNew)) then
            !$omp do schedule(runtime)
            do iCell = 1, nCells
               highFreqThicknessNew(:, iCell) = highFreqThicknessCur(:, iCell)
            end do
            !$omp end do
         end if

         if (associated(lowFreqDivergenceNew)) then
            !$omp do schedule(runtime)
            do iCell = 1, nCells
               lowFreqDivergenceNew(:, iCell) = lowFreqDivergenceCur(:, iCell)
            end do
            !$omp end do
         endif

         block => block % next
      end do

      call mpas_timer_stop("se prep")
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! BEGIN large iteration loop
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      n_bcl_iter = config_n_bcl_iter_mid
      n_bcl_iter(1) = config_n_bcl_iter_beg
      n_bcl_iter(config_n_ts_iter) = config_n_bcl_iter_end

      do split_explicit_step = 1, config_n_ts_iter
         call mpas_timer_start('se loop')

         stage1_tend_time = min(split_explicit_step,2)

         call mpas_pool_get_subpool(domain % blocklist % structs, 'diagnostics', diagnosticsPool)

         call mpas_threading_barrier()
         ! ---  update halos for diagnostic ocean boundayr layer depth
         if (config_use_cvmix_kpp) then
            call mpas_timer_start("se halo diag obd")
            call mpas_dmpar_field_halo_exch(domain, 'boundaryLayerDepth')
            call mpas_timer_stop("se halo diag obd")
         end if

         ! ---  update halos for diagnostic variables
         call mpas_timer_start("se halo diag")

         call mpas_dmpar_field_halo_exch(domain, 'normalizedRelativeVorticityEdge')
         if (config_mom_del4 > 0.0) then
           call mpas_dmpar_field_halo_exch(domain, 'divergence')
           call mpas_dmpar_field_halo_exch(domain, 'relativeVorticity')
         end if
         call mpas_timer_stop("se halo diag")

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         !  Stage 1: Baroclinic velocity (3D) prediction, explicit with long timestep
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (config_use_freq_filtered_thickness) then
            call mpas_timer_start("se freq-filtered-thick computations")

            block => domain % blocklist
            do while (associated(block))
               call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
               call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
               call mpas_pool_get_subpool(block % structs, 'state', statepool)
               call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
               call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
               call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

               call ocn_tend_freq_filtered_thickness(tendPool, statePool, diagnosticsPool, meshPool, stage1_tend_time)
               block => block % next
            end do
            call mpas_timer_stop("se freq-filtered-thick computations")

            call mpas_threading_barrier()

            call mpas_timer_start("se freq-filtered-thick halo update")

            call mpas_dmpar_field_halo_exch(domain, 'tendHighFreqThickness')
            call mpas_dmpar_field_halo_exch(domain, 'tendLowFreqDivergence')

            call mpas_timer_stop("se freq-filtered-thick halo update")
            call mpas_threading_barrier()

            block => domain % blocklist
            do while (associated(block))
               call mpas_pool_get_dimension(block % dimensions, 'nCells', nCellsPtr)
               call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)

               call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
               call mpas_pool_get_subpool(block % structs, 'state', statePool)
               call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
               call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
               call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)

               call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

               call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessCur, 1)
               call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessNew, 2)

               call mpas_pool_get_array(tendPool, 'highFreqThickness', highFreqThicknessTend)

               nCells = nCellsPtr

               !$omp do schedule(runtime) private(k)
               do iCell = 1, nCells
                  do k = 1, maxLevelCell(iCell)
                     ! this is h^{hf}_{n+1}
                     highFreqThicknessNew(k,iCell) = highFreqThicknessCur(k,iCell) + dt * highFreqThicknessTend(k,iCell)
                  end do
               end do
               !$omp end do

               block => block % next
            end do

         endif

         ! compute velocity tendencies, T(u*,w*,p*)
         call mpas_timer_start("se bcl vel")

         call mpas_timer_start('se bcl vel tend')
         block => domain % blocklist
         do while (associated(block))
           call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
           call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
           call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
           call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
           call mpas_pool_get_subpool(block % structs, 'state', statePool)
           call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
           call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
           call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
           call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)

           call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
           call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, stage1_tend_time)
           call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)

           call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessNew, 2)

           call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)

           call ocn_tend_vel(tendPool, statePool, forcingPool, diagnosticsPool, meshPool, scratchPool, stage1_tend_time)

           block => block % next
         end do
         call mpas_timer_stop('se bcl vel tend')

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! BEGIN baroclinic iterations on linear Coriolis term
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         do j=1,n_bcl_iter(split_explicit_step)

            ! Use this G coefficient to avoid an if statement within the iEdge loop.
            if (trim(config_time_integrator) == 'unsplit_explicit') then
               split = 0
            elseif (trim(config_time_integrator) == 'split_explicit') then
               split = 1
            endif

            call mpas_timer_start('bcl iters on linear Coriolis')
            block => domain % blocklist
            do while (associated(block))
               call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
               call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
               call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)

               call mpas_pool_get_subpool(block % structs, 'state', statePool)
               call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
               call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
               call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
               call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
               call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

               call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
               call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
               call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)

               call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
               call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', normalBaroclinicVelocityCur, 1)
               call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', normalBaroclinicVelocityNew, 2)
               call mpas_pool_get_array(statePool, 'ssh', sshNew, 2)

               call mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)

               call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
               call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)

               ! Only need to loop over the 1 halo, since there is a halo exchange immediately after this computation.
               nEdges = nEdgesArray( 1 )

               ! Put f*normalBaroclinicVelocity^{perp} in normalVelocityNew as a work variable
               call ocn_fuperp(statePool, meshPool, 2)

               allocate(uTemp(nVertLevels))

               !$omp do schedule(runtime) private(cell1, cell2, k, normalThicknessFluxSum, thicknessSum)
               do iEdge = 1, nEdges
                  cell1 = cellsOnEdge(1,iEdge)
                  cell2 = cellsOnEdge(2,iEdge)

                  uTemp = 0.0_RKIND  ! could put this after with uTemp(maxleveledgetop+1:nvertlevels)=0
                  do k = 1, maxLevelEdgeTop(iEdge)

                     ! normalBaroclinicVelocityNew = normalBaroclinicVelocityOld + dt*(-f*normalBaroclinicVelocityPerp
                     !                             + T(u*,w*,p*) + g*grad(SSH*) )
                     ! Here uNew is a work variable containing -fEdge(iEdge)*normalBaroclinicVelocityPerp(k,iEdge)
                      uTemp(k) = normalBaroclinicVelocityCur(k,iEdge) &
                         + dt * (normalVelocityTend(k,iEdge) &
                         + normalVelocityNew(k,iEdge) &  ! this is f*normalBaroclinicVelocity^{perp}
                         + split * gravity * (  sshNew(cell2) - sshNew(cell1) ) &
                          / dcEdge(iEdge) )
                  enddo

                  ! thicknessSum is initialized outside the loop because on land boundaries
                  ! maxLevelEdgeTop=0, but I want to initialize thicknessSum with a
                  ! nonzero value to avoid a NaN.
                  normalThicknessFluxSum = layerThicknessEdge(1,iEdge) * uTemp(1)
                  thicknessSum  = layerThicknessEdge(1,iEdge)

                  do k = 2, maxLevelEdgeTop(iEdge)
                     normalThicknessFluxSum = normalThicknessFluxSum + layerThicknessEdge(k,iEdge) * uTemp(k)
                     thicknessSum  =  thicknessSum + layerThicknessEdge(k,iEdge)
                  enddo
                  barotropicForcing(iEdge) = split * normalThicknessFluxSum / thicknessSum / dt


                  do k = 1, maxLevelEdgeTop(iEdge)
                     ! These two steps are together here:
                     !{\bf u}'_{k,n+1} = {\bf u}'_{k,n} - \Delta t {\overline {\bf G}}
                     !{\bf u}'_{k,n+1/2} = \frac{1}{2}\left({\bf u}^{'}_{k,n} +{\bf u}'_{k,n+1}\right)
                     ! so that normalBaroclinicVelocityNew is at time n+1/2
                     normalBaroclinicVelocityNew(k,iEdge) = 0.5_RKIND*( &
                       normalBaroclinicVelocityCur(k,iEdge) + uTemp(k) - dt * barotropicForcing(iEdge))

                  enddo

               enddo ! iEdge
               !$omp end do

               deallocate(uTemp)

               block => block % next
            end do

            call mpas_threading_barrier()

            call mpas_timer_start("se halo normalBaroclinicVelocity")
            call mpas_dmpar_field_halo_exch(domain, 'normalBaroclinicVelocity', timeLevel=2)
            call mpas_timer_stop("se halo normalBaroclinicVelocity")

            call mpas_timer_stop('bcl iters on linear Coriolis')

         end do  ! do j=1,config_n_bcl_iter

         call mpas_timer_start('se halo barotropicForcing')
         call mpas_dmpar_field_halo_exch(domain, 'barotropicForcing')
         call mpas_timer_stop('se halo barotropicForcing')

         call mpas_timer_stop("se bcl vel")
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! END baroclinic iterations on linear Coriolis term
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         !  Stage 2: Barotropic velocity (2D) prediction, explicitly subcycled
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         call mpas_timer_start("se btr vel")

         oldBtrSubcycleTime = 1
         newBtrSubcycleTime = 2

         if (trim(config_time_integrator) == 'unsplit_explicit') then

            call mpas_timer_start('btr vel ue')
            block => domain % blocklist
            do while (associated(block))
               call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
               call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
               call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)

               call mpas_pool_get_subpool(block % structs, 'state', statePool)
               call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
               call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
               call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

               call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityNew, 2)
               call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
               call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', normalBaroclinicVelocityNew, 2)

               call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)
               call mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', normalGMBolusVelocity)

               call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)

               nEdges = nEdgesPtr

               ! For Split_Explicit unsplit, simply set normalBarotropicVelocityNew=0, normalBarotropicVelocitySubcycle=0, and
               ! uNew=normalBaroclinicVelocityNew

               !$omp do schedule(runtime) private(k)
               do iEdge = 1, nEdges
                  normalBarotropicVelocityNew(iEdge) = 0.0_RKIND
                  do k = 1, nVertLevels
                     normalVelocityNew(k, iEdge)  = normalBaroclinicVelocityNew(k, iEdge)

                     ! normalTransportVelocity = normalBaroclinicVelocity + normalGMBolusVelocity
                     ! This is u used in advective terms for layerThickness and tracers
                     ! in tendency calls in stage 3.
                     normalTransportVelocity(k,iEdge) = edgeMask(k,iEdge) &
                           *( normalBaroclinicVelocityNew(k,iEdge) + normalGMBolusVelocity(k,iEdge) )

                  enddo
               end do  ! iEdge
               !$omp end do

               block => block % next
            end do  ! block
            call mpas_timer_stop('btr vel ue')

         elseif (trim(config_time_integrator) == 'split_explicit') then

            ! Initialize variables for barotropic subcycling
            call mpas_timer_start('btr vel se init')
            block => domain % blocklist
            do while (associated(block))
               call mpas_pool_get_dimension(block % dimensions, 'nCells', nCellsPtr)
               call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
               call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
               call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

               call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
               call mpas_pool_get_subpool(block % structs, 'state', statePool)
               call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

               call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)
               call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux', barotropicThicknessFlux)

               call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
               call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, oldBtrSubcycleTime)
               call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleCur, &
                                      oldBtrSubcycleTime)
               call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityCur, 1)
               call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityNew, 2)

               nCells = nCellsPtr
               nEdges = nEdgesPtr

               if (config_filter_btr_mode) then
                  !$omp do schedule(runtime)
                  do iEdge = 1, nEdges
                     barotropicForcing(iEdge) = 0.0_RKIND
                  end do
                  !$omp end do
               endif

               !$omp do schedule(runtime)
               do iCell = 1, nCells
                  ! sshSubcycleOld = sshOld
                  sshSubcycleCur(iCell) = sshCur(iCell)
               end do
               !$omp end do

               !$omp do schedule(runtime)
               do iEdge = 1, nEdges

                  ! normalBarotropicVelocitySubcycleOld = normalBarotropicVelocityOld
                  normalBarotropicVelocitySubcycleCur(iEdge) = normalBarotropicVelocityCur(iEdge)

                  ! normalBarotropicVelocityNew = BtrOld  This is the first for the summation
                  normalBarotropicVelocityNew(iEdge) = normalBarotropicVelocityCur(iEdge)

                  ! barotropicThicknessFlux = 0
                  barotropicThicknessFlux(iEdge) = 0.0_RKIND
               end do
               !$omp end do

               block => block % next
            end do  ! block
            call mpas_timer_stop('btr vel se init')

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! BEGIN Barotropic subcycle loop
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! Allocate subcycled scratch fields before starting subcycle loop
            call mpas_pool_get_subpool(domain % blocklist % structs, 'scratch', scratchPool)
            call mpas_pool_get_field(scratchPool, 'btrvel_temp', btrvel_tempField)
            call mpas_allocate_scratch_field(btrvel_tempField, .false.)

            call mpas_threading_barrier()

            call mpas_timer_start('btr se subcycle loop')
            do j = 1, nBtrSubcycles * config_btr_subcycle_loop_factor
               cellHaloComputeCounter = config_num_halos
               edgeHaloComputeCounter = config_num_halos + 1

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! Barotropic subcycle: VELOCITY PREDICTOR STEP
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               if (config_btr_gam1_velWt1 > 1.0e-12_RKIND) then  ! only do this part if it is needed in next SSH solve
                  uPerpTime = oldBtrSubcycleTime

                  block => domain % blocklist
                  do while (associated(block))
                     call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
                     call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

                     call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
                     call mpas_pool_get_subpool(block % structs, 'state', statePool)
                     call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
                     call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

                     call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
                     call mpas_pool_get_array(meshPool, 'nEdgesOnEdge', nEdgesOnEdge)
                     call mpas_pool_get_array(meshPool, 'edgesOnEdge', edgesOnEdge)
                     call mpas_pool_get_array(meshPool, 'weightsOnEdge', weightsOnEdge)
                     call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
                     call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
                     call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)

                     call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleCur, &
                                              uPerpTime)
                     call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleNew, &
                                              newBtrSubcycleTime)
                     call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, oldBtrSubcycleTime)

                     call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)

                     nEdges = nEdgesPtr
                     nEdges = nEdgesArray( edgeHaloComputeCounter )

                     !$omp do schedule(runtime) private(cell1, cell2, CoriolisTerm, i, eoe)
                     do iEdge = 1, nEdges

                        cell1 = cellsOnEdge(1,iEdge)
                        cell2 = cellsOnEdge(2,iEdge)

                        ! Compute the barotropic Coriolis term, -f*uPerp
                        CoriolisTerm = 0.0_RKIND
                        do i = 1, nEdgesOnEdge(iEdge)
                           eoe = edgesOnEdge(i,iEdge)
                           CoriolisTerm = CoriolisTerm + weightsOnEdge(i,iEdge) &
                                        * normalBarotropicVelocitySubcycleCur(eoe) * fEdge(eoe)
                        end do

                        ! normalBarotropicVelocityNew = normalBarotropicVelocityOld + dt/J*(-f*normalBarotropicVelocityoldPerp
                        !                             - g*grad(SSH) + G)
                        normalBarotropicVelocitySubcycleNew(iEdge) &
                          = (normalBarotropicVelocitySubcycleCur(iEdge) &
                          + dt / nBtrSubcycles * (CoriolisTerm - gravity &
                          * (sshSubcycleCur(cell2) - sshSubcycleCur(cell1) ) &
                          / dcEdge(iEdge) + barotropicForcing(iEdge))) * edgeMask(1, iEdge)
                     end do
                     !$omp end do

                     block => block % next
                  end do  ! block
              endif ! config_btr_gam1_velWt1>1.0e-12

              ! 1 Halo from edges is corrupted here, so reduce the edge halo layers by 1
              edgeHaloComputeCounter = edgeHaloComputeCounter - 1

              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ! Barotropic subcycle: SSH PREDICTOR STEP
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              block => domain % blocklist
              do while (associated(block))
                call mpas_pool_get_dimension(block % dimensions, 'nCells', nCellsPtr)
                call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
                call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
                call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

                call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
                call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
                call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
                call mpas_pool_get_subpool(block % structs, 'state', statePool)
                call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
                call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

                call mpas_pool_get_array(tendPool, 'ssh', sshTend)

                call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
                call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
                call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
                call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
                call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
                call mpas_pool_get_array(meshPool, 'refBottomDepthTopOfCell', refBottomDepthTopOfCell)
                call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
                call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
                call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

                call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, oldBtrSubcycleTime)
                call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleNew, newBtrSubcycleTime)
                call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleCur, &
                                            oldBtrSubcycleTime)
                call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleNew, &
                                            newBtrSubcycleTime)

                call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux', barotropicThicknessFlux)

                nCells = nCellsPtr
                nEdges = nEdgesPtr

                nCells = nCellsArray( cellHaloComputeCounter )
                nEdges = nEdgesArray( edgeHaloComputeCounter )

                if (config_btr_solve_SSH2) then
                   ! If config_btr_solve_SSH2=.true., then do NOT accumulate barotropicThicknessFlux in this SSH predictor
                   ! section, because it will be accumulated in the SSH corrector section.
                   barotropicThicknessFlux_coeff = 0.0_RKIND
                else
                   ! otherwise, DO accumulate barotropicThicknessFlux in this SSH predictor section
                   barotropicThicknessFlux_coeff = 1.0_RKIND
                endif

                ! config_btr_gam1_velWt1 sets the forward weighting of velocity in the SSH computation
                ! config_btr_gam1_velWt1=  1     flux = normalBarotropicVelocityNew*H
                ! config_btr_gam1_velWt1=0.5     flux = 1/2*(normalBarotropicVelocityNew+normalBarotropicVelocityOld)*H
                ! config_btr_gam1_velWt1=  0     flux = normalBarotropicVelocityOld*H

                !$omp do schedule(runtime) private(i, iEdge, cell1, cell2, sshEdge, thicknessSum, flux)
                do iCell = 1, nCells
                  sshTend(iCell) = 0.0_RKIND
                  do i = 1, nEdgesOnCell(iCell)
                    iEdge = edgesOnCell(i, iCell)

                    cell1 = cellsOnEdge(1, iEdge)
                    cell2 = cellsOnEdge(2, iEdge)

                    sshEdge = 0.5_RKIND * (sshSubcycleCur(cell1) + sshSubcycleCur(cell2) )

                   ! method 0: orig, works only without pbc:
                   !thicknessSum = sshEdge + refBottomDepthTopOfCell(maxLevelEdgeTop(iEdge)+1)

                   ! method 1, matches method 0 without pbcs, works with pbcs.
                   thicknessSum = sshEdge + min(bottomDepth(cell1), bottomDepth(cell2))

                   ! method 2: may be better than method 1.
                   ! Take average  of full thickness at two neighboring cells.
                   !thicknessSum = sshEdge + 0.5 *( bottomDepth(cell1) + bottomDepth(cell2) )


                    flux = ((1.0-config_btr_gam1_velWt1) * normalBarotropicVelocitySubcycleCur(iEdge) &
                           + config_btr_gam1_velWt1 * normalBarotropicVelocitySubcycleNew(iEdge)) &
                           * thicknessSum

                    sshTend(iCell) = sshTend(iCell) + edgeSignOncell(i, iCell) * flux &
                           * dvEdge(iEdge)

                  end do

                  ! SSHnew = SSHold + dt/J*(-div(Flux))
                  sshSubcycleNew(iCell) = sshSubcycleCur(iCell) + dt / nBtrSubcycles * sshTend(iCell) / areaCell(iCell)
                end do
                !$omp end do

                !$omp do schedule(runtime) private(cell1, cell2, sshEdge, thicknessSum, flux)
                do iEdge = 1, nEdges
                   cell1 = cellsOnEdge(1,iEdge)
                   cell2 = cellsOnEdge(2,iEdge)

                   sshEdge = 0.5_RKIND * (sshSubcycleCur(cell1) &
                             + sshSubcycleCur(cell2) )

                   ! method 0: orig, works only without pbc:
                   !thicknessSum = sshEdge + refBottomDepthTopOfCell(maxLevelEdgeTop(iEdge)+1)

                   ! method 1, matches method 0 without pbcs, works with pbcs.
                   thicknessSum = sshEdge + min(bottomDepth(cell1), bottomDepth(cell2))

                   ! method 2: may be better than method 1.
                   ! take average  of full thickness at two neighboring cells
                   !thicknessSum = sshEdge + 0.5 *(  bottomDepth(cell1) &
                   !                       + bottomDepth(cell2) )

                   flux = ((1.0-config_btr_gam1_velWt1) * normalBarotropicVelocitySubcycleCur(iEdge) &
                          + config_btr_gam1_velWt1 * normalBarotropicVelocitySubcycleNew(iEdge)) &
                          * thicknessSum

                   barotropicThicknessFlux(iEdge) = barotropicThicknessFlux(iEdge) + barotropicThicknessFlux_coeff * flux
                end do
                !$omp end do

                block => block % next
              end do  ! block

              ! 1 cell halo layer is now corrupted, so remove one from computing on.
              cellHaloComputeCounter = cellHaloComputeCounter - 1

              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ! Barotropic subcycle: VELOCITY CORRECTOR STEP
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

              ! 1 edge halo is already corrupted from the predictor step.
              do BtrCorIter = 1, config_n_btr_cor_iter
                uPerpTime = newBtrSubcycleTime

                block => domain % blocklist
                do while (associated(block))
                   call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
                   call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

                   call mpas_pool_get_subpool(block % structs, 'state', statePool)
                   call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
                   call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
                   call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
                   call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

                   call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleCur, &
                                            oldBtrSubcycleTime)
                   call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleNew, &
                                            newBtrSubcycleTime)
                   call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, oldBtrSubcycleTime)
                   call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleNew, newBtrSubcycleTime)

                   call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
                   call mpas_pool_get_array(meshPool, 'nEdgesOnEdge', nEdgesOnEdge)
                   call mpas_pool_get_array(meshPool, 'edgesOnEdge', edgesOnEdge)
                   call mpas_pool_get_array(meshPool, 'weightsOnEdge', weightsOnEdge)
                   call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
                   call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
                   call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)

                   call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)

                   call mpas_pool_get_field(scratchPool, 'btrvel_temp', btrvel_tempField)
                   btrvel_temp => btrvel_tempField % array

                   ! Need to initialize btr_vel_temp over the one more halo than we're computing over
                   nEdges = nEdgesPtr

                   nEdges = nEdgesArray( min(edgeHaloComputeCounter + 1, config_num_halos + 1) )

                   !$omp do schedule(runtime)
                   do iEdge = 1, nEdges+1
                      btrvel_temp(iEdge) = normalBarotropicVelocitySubcycleNew(iEdge)
                   end do
                   !$omp end do

                   nEdges = nEdgesArray( edgeHaloComputeCounter )

                   !$omp do schedule(runtime) private(cell1, cell2, eoe, CoriolisTerm, i, sshCell1, sshCell2)
                   do iEdge = 1, nEdges
                     cell1 = cellsOnEdge(1,iEdge)
                     cell2 = cellsOnEdge(2,iEdge)

                     ! Compute the barotropic Coriolis term, -f*uPerp
                     CoriolisTerm = 0.0_RKIND
                     do i = 1, nEdgesOnEdge(iEdge)
                        eoe = edgesOnEdge(i,iEdge)
                        CoriolisTerm = CoriolisTerm + weightsOnEdge(i,iEdge) &
                             !* normalBarotropicVelocitySubcycleNew(eoe) &
                             * btrvel_temp(eoe) * fEdge(eoe)
                     end do

                     ! In this final solve for velocity, SSH is a linear
                     ! combination of SSHold and SSHnew.
                     sshCell1 = (1-config_btr_gam2_SSHWt1) * sshSubcycleCur(cell1) + config_btr_gam2_SSHWt1 * sshSubcycleNew(cell1)
                     sshCell2 = (1-config_btr_gam2_SSHWt1) * sshSubcycleCur(cell2) + config_btr_gam2_SSHWt1 * sshSubcycleNew(cell2)

                     ! normalBarotropicVelocityNew = normalBarotropicVelocityOld + dt/J*(-f*normalBarotropicVelocityoldPerp
                     !                             - g*grad(SSH) + G)
                     normalBarotropicVelocitySubcycleNew(iEdge) = (normalBarotropicVelocitySubcycleCur(iEdge) &
                         + dt / nBtrSubcycles *(CoriolisTerm - gravity *(sshCell2 - sshCell1) / dcEdge(iEdge) &
                         + barotropicForcing(iEdge))) * edgeMask(1,iEdge)
                   end do
                   !$omp end do

                   block => block % next
                end do  ! block

                if ( edgeHaloComputeCounter == 1 .and. BtrCorIter /= config_n_btr_cor_iter ) then
                   edgeHaloComputeCounter = config_num_halos + 1
                   call mpas_timer_start("se halo btr vel corr")
                   call mpas_dmpar_field_halo_exch(domain, 'normalBarotropicVelocitySubcycle', timeLevel=newBtrSubcycleTime)
                   call mpas_timer_stop("se halo btr vel corr")
                else
                   edgeHaloComputeCounter = edgeHaloComputeCounter - 1
                end if

              end do !do BtrCorIter=1,config_n_btr_cor_iter

              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ! Barotropic subcycle: SSH CORRECTOR STEP
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              if (config_btr_solve_SSH2) then

                block => domain % blocklist
                do while (associated(block))
                   call mpas_pool_get_dimension(block % dimensions, 'nCells', nCellsPtr)
                   call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
                   call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
                   call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

                   call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
                   call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
                   call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
                   call mpas_pool_get_subpool(block % structs, 'state', statePool)
                   call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
                   call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

                   call mpas_pool_get_array(tendPool, 'ssh', sshTend)

                   call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
                   call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
                   call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
                   call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
                   call mpas_pool_get_array(meshPool, 'refBottomDepthTopOfCell', refBottomDepthTopOfCell)
                   call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
                   call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
                   call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)

                   call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, oldBtrSubcycleTime)
                   call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleNew, newBtrSubcycleTime)
                   call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleCur, &
                                            oldBtrSubcycleTime)
                   call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleNew, &
                                            newBtrSubcycleTime)

                   call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux', barotropicThicknessFlux)

                   nCells = nCellsPtr
                   nEdges = nEdgesPtr

                   nCells = nCellsArray( cellHaloComputeCounter )
                   nEdges = nEdgesArray( edgeHaloComputeCounter )

                   ! config_btr_gam3_velWt2 sets the forward weighting of velocity in the SSH computation
                   ! config_btr_gam3_velWt2=  1     flux = normalBarotropicVelocityNew*H
                   ! config_btr_gam3_velWt2=0.5     flux = 1/2*(normalBarotropicVelocityNew+normalBarotropicVelocityOld)*H
                   ! config_btr_gam3_velWt2=  0     flux = normalBarotropicVelocityOld*H

                   !$omp do schedule(runtime) private(i, iEdge, cell1, cell2, sshCell1, sshCell2, sshEdge, thicknessSum, flux)
                   do iCell = 1, nCells
                     sshTend(iCell) = 0.0_RKIND
                     do i = 1, nEdgesOnCell(iCell)
                       iEdge = edgesOnCell(i, iCell)

                       cell1 = cellsOnEdge(1,iEdge)
                       cell2 = cellsOnEdge(2,iEdge)

                       ! SSH is a linear combination of SSHold and SSHnew.
                       sshCell1 = (1-config_btr_gam2_SSHWt1)* sshSubcycleCur(cell1) &
                                 +   config_btr_gam2_SSHWt1 * sshSubcycleNew(cell1)
                       sshCell2 = (1-config_btr_gam2_SSHWt1)* sshSubcycleCur(cell2) &
                                 +   config_btr_gam2_SSHWt1 * sshSubcycleNew(cell2)

                       sshEdge = 0.5_RKIND * (sshCell1 + sshCell2)

                      ! method 0: orig, works only without pbc:
                      !thicknessSum = sshEdge + refBottomDepthTopOfCell(maxLevelEdgeTop(iEdge)+1)

                      ! method 1, matches method 0 without pbcs, works with pbcs.
                      thicknessSum = sshEdge + min(bottomDepth(cell1), bottomDepth(cell2))

                      ! method 2: may be better than method 1.
                      ! take average  of full thickness at two neighboring cells
                      !thicknessSum = sshEdge + 0.5 *( bottomDepth(cell1) + bottomDepth (cell2) )

                       flux = ((1.0-config_btr_gam3_velWt2) * normalBarotropicVelocitySubcycleCur(iEdge) &
                              + config_btr_gam3_velWt2 * normalBarotropicVelocitySubcycleNew(iEdge)) &
                              * thicknessSum

                       sshTend(iCell) = sshTend(iCell) + edgeSignOnCell(i, iCell) * flux &
                              * dvEdge(iEdge)

                     end do

                     ! SSHnew = SSHold + dt/J*(-div(Flux))
                     sshSubcycleNew(iCell) = sshSubcycleCur(iCell) &
                          + dt / nBtrSubcycles * sshTend(iCell) / areaCell(iCell)
                   end do
                   !$omp end do

                   !$omp do schedule(runtime) private(cell1, cell2, sshCell1, sshCell2, sshEdge, thicknessSum, flux)
                   do iEdge = 1, nEdges
                      cell1 = cellsOnEdge(1,iEdge)
                      cell2 = cellsOnEdge(2,iEdge)

                      ! SSH is a linear combination of SSHold and SSHnew.
                      sshCell1 = (1-config_btr_gam2_SSHWt1)* sshSubcycleCur(cell1) + config_btr_gam2_SSHWt1 * sshSubcycleNew(cell1)
                      sshCell2 = (1-config_btr_gam2_SSHWt1)* sshSubcycleCur(cell2) + config_btr_gam2_SSHWt1 * sshSubcycleNew(cell2)
                      sshEdge = 0.5_RKIND * (sshCell1 + sshCell2)

                      ! method 0: orig, works only without pbc:
                      !thicknessSum = sshEdge + refBottomDepthTopOfCell(maxLevelEdgeTop(iEdge)+1)

                      ! method 1, matches method 0 without pbcs, works with pbcs.
                      thicknessSum = sshEdge + min(bottomDepth(cell1), bottomDepth(cell2))

                      ! method 2, better, I think.
                      ! take average  of full thickness at two neighboring cells
                      !thicknessSum = sshEdge + 0.5 *( bottomDepth(cell1) + bottomDepth(cell2) )

                      flux = ((1.0-config_btr_gam3_velWt2) * normalBarotropicVelocitySubcycleCur(iEdge) &
                             + config_btr_gam3_velWt2 * normalBarotropicVelocitySubcycleNew(iEdge)) &
                             * thicknessSum

                      barotropicThicknessFlux(iEdge) = barotropicThicknessFlux(iEdge) + flux
                   end do
                   !$omp end do

                   block => block % next
                end do  ! block
               endif ! config_btr_solve_SSH2

               !   boundary update on SSHnew
               call mpas_timer_start("se halo subcycle")
               call mpas_dmpar_exch_group_create(domain, subcycleGroupName)
               call mpas_dmpar_exch_group_add_field(domain, subcycleGroupName, 'sshSubcycle', timeLevel=newBtrSubcycleTime)
               call mpas_dmpar_exch_group_add_field(domain, subcycleGroupName, 'normalBarotropicVelocitySubcycle', timeLevel=newBtrSubcycleTime)

               call mpas_threading_barrier()
               call mpas_dmpar_exch_group_full_halo_exch(domain, subcycleGroupName)

               call mpas_dmpar_exch_group_destroy(domain, subcycleGroupName)
               call mpas_timer_stop("se halo subcycle")

               ! Reset the halo counters
               edgeHaloComputeCounter = config_num_halos + 1
               cellHaloComputeCounter = config_num_halos

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! Barotropic subcycle: Accumulate running sums, advance timestep pointers
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               block => domain % blocklist
               do while (associated(block))
                  call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
                  call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

                  call mpas_pool_get_subpool(block % structs, 'state', statePool)
                  call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

                  call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityNew, 2)
                  call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleNew, &
                                           newBtrSubcycleTime)

                  ! normalBarotropicVelocityNew = normalBarotropicVelocityNew + normalBarotropicVelocitySubcycleNEW
                  ! This accumulates the sum.
                  ! If the Barotropic Coriolis iteration is limited to one, this could
                  ! be merged with the above code.

                  nEdges = nEdgesPtr

                  !$omp do schedule(runtime)
                  do iEdge = 1, nEdges
                       normalBarotropicVelocityNew(iEdge) = normalBarotropicVelocityNew(iEdge) &
                                                          + normalBarotropicVelocitySubcycleNew(iEdge)
                  end do  ! iEdge
                  !$omp end do

                  block => block % next
               end do  ! block

               ! advance time pointers
               oldBtrSubcycleTime = mod(oldBtrSubcycleTime,2)+1
               newBtrSubcycleTime = mod(newBtrSubcycleTime,2)+1

            end do ! j=1,nBtrSubcycles
            call mpas_timer_stop('btr se subcycle loop')

            call mpas_threading_barrier()
            call mpas_pool_get_subpool(domain % blocklist % structs, 'scratch', scratchPool)
            call mpas_pool_get_field(scratchPool, 'btrvel_temp', btrvel_tempField)
            call mpas_deallocate_scratch_field(btrvel_tempField, .false.)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! END Barotropic subcycle loop
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! Normalize Barotropic subcycle sums: ssh, normalBarotropicVelocity, and F
            call mpas_timer_start('btr se norm')
            block => domain % blocklist
            do while (associated(block))
               call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
               call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

               call mpas_pool_get_subpool(block % structs, 'state', statePool)
               call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
               call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

               call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityNew, 2)

               call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux', barotropicThicknessFlux)

               nEdges = nEdgesPtr

               nEdges = nEdgesArray(1)

               !$omp do schedule(runtime)
               do iEdge = 1, nEdges
                  barotropicThicknessFlux(iEdge) = barotropicThicknessFlux(iEdge) &
                      / (nBtrSubcycles * config_btr_subcycle_loop_factor)

                  normalBarotropicVelocityNew(iEdge) = normalBarotropicVelocityNew(iEdge) &
                     / (nBtrSubcycles * config_btr_subcycle_loop_factor + 1)
               end do
               !$omp end do

               block => block % next
            end do  ! block
            call mpas_timer_stop('btr se norm')

            call mpas_threading_barrier()

            ! boundary update on F
            call mpas_timer_start("se halo F and btr vel")
            call mpas_dmpar_exch_group_create(domain, finalBtrGroupName)

            call mpas_dmpar_exch_group_add_field(domain, finalBtrGroupName, 'barotropicThicknessFlux')
            call mpas_dmpar_exch_group_add_field(domain, finalBtrGroupName, 'normalBarotropicVelocity', timeLevel=2)

            call mpas_threading_barrier()
            call mpas_dmpar_exch_group_full_halo_exch(domain, finalBtrGroupName)

            call mpas_dmpar_exch_group_destroy(domain, finalBtrGroupName)
            call mpas_timer_stop("se halo F and btr vel")

            call mpas_threading_barrier()

            ! Check that you can compute SSH using the total sum or the individual increments
            ! over the barotropic subcycles.
            ! efficiency: This next block of code is really a check for debugging, and can
            ! be removed later.
            call mpas_timer_start('btr se ssh verif')
            block => domain % blocklist
            do while (associated(block))
               call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
               call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
               call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)

               call mpas_pool_get_subpool(block % structs, 'state', statePool)
               call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
               call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
               call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

               call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityNew, 2)
               call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', normalBaroclinicVelocityNew, 2)

               call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)
               call mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', normalGMBolusVelocity)
               call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
               call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux', barotropicThicknessFlux)

               call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
               call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)

               nEdges = nEdgesPtr

               nEdges = nEdgesArray( config_num_halos )

               allocate(uTemp(nVertLevels))

               ! Correction velocity    normalVelocityCorrection = (Flux - Sum(h u*))/H
               ! or, for the full latex version:
               !{\bf u}^{corr} = \left( {\overline {\bf F}}
               !  - \sum_{k=1}^{N^{edge}} h_{k,*}^{edge}  {\bf u}_k^{avg} \right)
               ! \left/ \sum_{k=1}^{N^{edge}} h_{k,*}^{edge}   \right.

               if (config_vel_correction) then
                  useVelocityCorrection = 1
               else
                  useVelocityCorrection = 0
               endif

               !$omp do schedule(runtime) private(k, normalThicknessFluxSum, thicknessSum, normalVelocityCorrection)
               do iEdge = 1, nEdges

                  ! velocity for normalVelocityCorrectionection is normalBarotropicVelocity + normalBaroclinicVelocity + uBolus
                  uTemp(:) = normalBarotropicVelocityNew(iEdge) + normalBaroclinicVelocityNew(:,iEdge) &
                           + normalGMBolusVelocity(:,iEdge)

                  ! thicknessSum is initialized outside the loop because on land boundaries
                  ! maxLevelEdgeTop=0, but I want to initialize thicknessSum with a
                  ! nonzero value to avoid a NaN.
                  normalThicknessFluxSum = layerThicknessEdge(1,iEdge) * uTemp(1)
                  thicknessSum  = layerThicknessEdge(1,iEdge)

                  do k = 2, maxLevelEdgeTop(iEdge)
                     normalThicknessFluxSum = normalThicknessFluxSum + layerThicknessEdge(k,iEdge) * uTemp(k)
                     thicknessSum  =  thicknessSum + layerThicknessEdge(k,iEdge)
                  enddo

                  normalVelocityCorrection = useVelocityCorrection * (( barotropicThicknessFlux(iEdge) - normalThicknessFluxSum) &
                                           / thicknessSum)

                  do k = 1, nVertLevels

                     ! normalTransportVelocity = normalBarotropicVelocity + normalBaroclinicVelocity + normalGMBolusVelocity
                     !                         + normalVelocityCorrection
                     ! This is u used in advective terms for layerThickness and tracers
                     ! in tendency calls in stage 3.
!mrp note: in QC version, there is an if (config_use_standardGM) on adding normalGMBolusVelocity
! I think it is not needed because normalGMBolusVelocity=0 when GM not on.
                     normalTransportVelocity(k,iEdge) &
                           = edgeMask(k,iEdge) &
                           *( normalBarotropicVelocityNew(iEdge) + normalBaroclinicVelocityNew(k,iEdge) &
                           + normalGMBolusVelocity(k,iEdge) + normalVelocityCorrection )
                  enddo

               end do ! iEdge
               !$omp end do

               deallocate(uTemp)

               block => block % next
            end do  ! block
            call mpas_timer_stop('btr se ssh verif')

         endif ! split_explicit

         call mpas_timer_stop("se btr vel")

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         !  Stage 3: Tracer, density, pressure, vertical velocity prediction
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! Thickness tendency computations and thickness halo updates are completed before tracer
         ! tendency computations to allow monotonic advection.
         call mpas_timer_start('se thick tend')
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
            call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
            call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)

            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
            call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
            call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessNew, 2)

            call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
            call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)
            call mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', vertAleTransportTop)

            ! compute vertAleTransportTop.  Use normalTransportVelocity for advection of layerThickness and tracers.
            ! Use time level 1 values of layerThickness and layerThicknessEdge because
            ! layerThickness has not yet been computed for time level 2.
            call mpas_timer_start('thick vert trans vel top')
            if (associated(highFreqThicknessNew)) then
               call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
                 layerThicknessCur, layerThicknessEdge, normalTransportVelocity, &
                 sshCur, dt, vertAleTransportTop, err, highFreqThicknessNew)
            else
               call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, scratchPool, &
                 layerThicknessCur, layerThicknessEdge, normalTransportVelocity, &
                 sshCur, dt, vertAleTransportTop, err)
            endif
            call mpas_timer_stop('thick vert trans vel top')

            call ocn_tend_thick(tendPool, forcingPool, diagnosticsPool, meshPool)

            block => block % next
         end do
         call mpas_timer_stop('se thick tend')

         call mpas_threading_barrier()

         ! update halo for thickness tendencies
         call mpas_timer_start("se halo thickness")

         call mpas_dmpar_field_halo_exch(domain, 'tendLayerThickness')

         call mpas_timer_stop("se halo thickness")

         call mpas_threading_barrier()

         call mpas_timer_start('se tracer tend', .false.)
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
            call mpas_pool_get_subpool(block % structs, 'shortwave', swForcingPool)
            call ocn_tend_tracer(tendPool, statePool, forcingPool, diagnosticsPool, meshPool, swForcingPool, scratchPool, dt, 2)

            block => block % next
         end do
         call mpas_timer_stop('se tracer tend')

         call mpas_threading_barrier()

         ! update halo for tracer tendencies
         call mpas_timer_start("se halo tracers")
         call mpas_pool_get_subpool(domain % blocklist % structs, 'tend', tendPool)
         call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)

         call mpas_pool_begin_iteration(tracersTendPool)
         do while ( mpas_pool_get_next_member(tracersTendPool, groupItr) )
            if ( groupItr % memberType == MPAS_POOL_FIELD ) then
               call mpas_dmpar_field_halo_exch(domain, groupItr % memberName)
            end if
         end do
         call mpas_timer_stop("se halo tracers")

         call mpas_threading_barrier()

         call mpas_timer_start('se loop fini')
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_dimension(block % dimensions, 'nCells', nCellsPtr)
            call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
            call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
            call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
            call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)

            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
            call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
            call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

            call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
            call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)
            call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)

            call mpas_pool_get_array(tracersPool, 'activeTracers', tracersGroupCur, 1)
            call mpas_pool_get_array(tracersPool, 'activeTracers', tracersGroupNew, 2)
            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)
            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
            call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessCur, 1)
            call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessNew, 2)
            call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceCur, 1)
            call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceNew, 2)
            call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityCur, 1)
            call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityNew, 2)
            call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', normalBaroclinicVelocityCur, 1)
            call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', normalBaroclinicVelocityNew, 2)

            call mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)
            call mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)
            call mpas_pool_get_array(tendPool, 'highFreqThickness', highFreqThicknessTend)
            call mpas_pool_get_array(tendPool, 'lowFreqDivergence', lowFreqDivergenceTend)

            call mpas_pool_get_array(tracersTendPool, 'activeTracersTend', activeTracersTend)

            nCells = nCellsPtr
            nEdges = nEdgesPtr

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            !  If iterating, reset variables for next iteration
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (split_explicit_step < config_n_ts_iter) then

               ! Get indices for dynamic tracers (Includes T&S).
               call mpas_pool_get_dimension(tracersPool, 'activeGRP_start', startIndex)
               call mpas_pool_get_dimension(tracersPool, 'activeGRP_end', endIndex)

               ! Only need T & S for earlier iterations,
               ! then all the tracers needed the last time through.

               !$omp do schedule(runtime) private(k, temp_h, temp, i)
               do iCell = 1, nCells
                  ! sshNew is a pointer, defined above.
                  do k = 1, maxLevelCell(iCell)

                     ! this is h_{n+1}
                     temp_h = layerThicknessCur(k,iCell) + dt * layerThicknessTend(k,iCell)

                     ! this is h_{n+1/2}
                     layerThicknessNew(k,iCell) = 0.5*( layerThicknessCur(k,iCell) + temp_h)

                     do i = startIndex, endIndex
                        ! This is Phi at n+1
                        temp = ( tracersGroupCur(i,k,iCell) * layerThicknessCur(k,iCell) + dt * activeTracersTend(i,k,iCell)) &
                             / temp_h

                        ! This is Phi at n+1/2
                        tracersGroupNew(i,k,iCell) = 0.5_RKIND * ( tracersGroupCur(i,k,iCell) + temp )
                     end do
                  end do
               end do ! iCell
               !$omp end do

               if (config_use_freq_filtered_thickness) then
                  !$omp do schedule(runtime) private(k, temp)
                  do iCell = 1, nCells
                     do k = 1, maxLevelCell(iCell)

                        ! h^{hf}_{n+1} was computed in Stage 1

                        ! this is h^{hf}_{n+1/2}
                        highFreqThicknessnew(k,iCell) = 0.5_RKIND * (highFreqThicknessCur(k,iCell) + highFreqThicknessNew(k,iCell))

                        ! this is D^{lf}_{n+1}
                        temp = lowFreqDivergenceCur(k,iCell) &
                         + dt * lowFreqDivergenceTend(k,iCell)

                        ! this is D^{lf}_{n+1/2}
                        lowFreqDivergenceNew(k,iCell) = 0.5_RKIND * (lowFreqDivergenceCur(k,iCell) + temp)
                     end do
                  end do
                  !$omp end do
               end if

               !$omp do schedule(runtime) private(k)
               do iEdge = 1, nEdges

                  do k = 1, nVertLevels

                     ! u = normalBarotropicVelocity + normalBaroclinicVelocity
                     ! here normalBaroclinicVelocity is at time n+1/2
                     ! This is u used in next iteration or step
                     normalVelocityNew(k,iEdge) = edgeMask(k,iEdge) * ( normalBarotropicVelocityNew(iEdge) &
                                                + normalBaroclinicVelocityNew(k,iEdge) )

                  enddo

               end do ! iEdge
               !$omp end do

               ! Efficiency note: We really only need this to compute layerThicknessEdge, density, pressure, and SSH
               ! in this diagnostics solve.
               call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, tracersPool, 2)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            !  If large iteration complete, compute all variables at time n+1
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            elseif (split_explicit_step == config_n_ts_iter) then

               !$omp do schedule(runtime) private(k)
               do iCell = 1, nCells
                  do k = 1, maxLevelCell(iCell)
                     ! this is h_{n+1}
                     layerThicknessNew(k,iCell) = layerThicknessCur(k,iCell) + dt * layerThicknessTend(k,iCell)
                  end do
               end do
               !$omp end do

               call mpas_pool_begin_iteration(tracersPool)
               do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
                  if ( groupItr % memberType == MPAS_POOL_FIELD ) then
                     configName = 'config_use_' // trim(groupItr % memberName)
                     call mpas_pool_get_config(domain % configs, configName, config_use_tracerGroup)

                     if ( config_use_tracerGroup ) then
                        call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupCur, 1)
                        call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroupNew, 2)

                        modifiedGroupName = trim(groupItr % memberName) // 'Tend'
                        call mpas_pool_get_array(tracersTendPool, modifiedGroupName, tracersGroupTend)

                        !$omp do schedule(runtime) private(k)
                        do iCell = 1, nCells
                           do k = 1, maxLevelCell(iCell)
                              tracersGroupNew(:,k,iCell) = (tracersGroupCur(:,k,iCell) * layerThicknessCur(k,iCell) + dt &
                                                         * tracersGroupTend(:,k,iCell) ) / layerThicknessNew(k,iCell)
!mrp limit Salinity
                              tracersGroupNew(indexSalinity,k,iCell) = max(0.001_RKIND, tracersGroupNew(indexSalinity,k,iCell))
                           end do
                        end do
                        !$omp end do
                     end if
                  end if
               end do

               if (config_use_freq_filtered_thickness) then
                  !$omp do schedule(runtime) private(k)
                  do iCell = 1, nCells
                     do k = 1, maxLevelCell(iCell)

                        ! h^{hf}_{n+1} was computed in Stage 1

                        ! this is D^{lf}_{n+1}
                        lowFreqDivergenceNew(k,iCell) = lowFreqDivergenceCur(k,iCell) + dt * lowFreqDivergenceTend(k,iCell)
                     end do
                  end do
                  !$omp end do
               end if

               ! Recompute final u to go on to next step.
               ! u_{n+1} = normalBarotropicVelocity_{n+1} + normalBaroclinicVelocity_{n+1}
               ! Right now normalBaroclinicVelocityNew is at time n+1/2, so back compute to get normalBaroclinicVelocity
               !   at time n+1 using normalBaroclinicVelocity_{n+1/2} = 1/2*(normalBaroclinicVelocity_n + u_Bcl_{n+1})
               ! so the following lines are
               ! u_{n+1} = normalBarotropicVelocity_{n+1} + 2*normalBaroclinicVelocity_{n+1/2} - normalBaroclinicVelocity_n
               ! note that normalBaroclinicVelocity is recomputed at the beginning of the next timestep due to Imp Vert mixing,
               ! so normalBaroclinicVelocity does not have to be recomputed here.

               !$omp do schedule(runtime) private(k)
               do iEdge = 1, nEdges
                  do k = 1, maxLevelEdgeTop(iEdge)
                     normalVelocityNew(k,iEdge) = normalBarotropicVelocityNew(iEdge) + 2 * normalBaroclinicVelocityNew(k,iEdge) &
                                                - normalBaroclinicVelocityCur(k,iEdge)
                  end do
               end do ! iEdges
               !$omp end do

            endif ! split_explicit_step

            block => block % next
         end do

         call mpas_timer_stop('se loop fini')
         call mpas_timer_stop('se loop')

      end do  ! split_explicit_step = 1, config_n_ts_iter
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! END large iteration loop
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call mpas_timer_start("se implicit vert mix")

      block => domain % blocklist
      do while(associated(block))
        call mpas_pool_get_subpool(block % structs, 'state', statePool)
        call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
        call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
        call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

        ! Call ocean diagnostic solve in preparation for vertical mixing.  Note
        ! it is called again after vertical mixing, because u and tracers change.
        ! For Richardson vertical mixing, only density, layerThicknessEdge, and kineticEnergyCell need to
        ! be computed.  For kpp, more variables may be needed.  Either way, this
        ! could be made more efficient by only computing what is needed for the
        ! implicit vmix routine that follows.
        call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, tracersPool, 2)

        ! Compute normalGMBolusVelocity; it will be added to the baroclinic modes in Stage 2 above.
        if (config_use_standardGM) then
           call ocn_gm_compute_Bolus_velocity(diagnosticsPool, meshPool, scratchPool)
        end if
        call ocn_vmix_implicit(dt, meshPool, diagnosticsPool, statePool, forcingPool, scratchPool, err, 2)

        block => block % next
      end do

      call mpas_threading_barrier()

      ! Update halo on u and tracers, which were just updated for implicit vertical mixing.  If not done,
      ! this leads to lack of volume conservation.  It is required because halo updates in stage 3 are only
      ! conducted on tendencies, not on the velocity and tracer fields.  So this update is required to
      ! communicate the change due to implicit vertical mixing across the boundary.
      call mpas_timer_start('se vmix halos')
      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)


      call mpas_timer_start('se vmix halos normalVelFld')
      call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', timeLevel=2)
      call mpas_timer_stop('se vmix halos normalVelFld')

      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
            call mpas_dmpar_field_halo_exch(domain, groupItr % memberName, timeLevel=2)
         end if
      end do
      call mpas_timer_stop('se vmix halos')

      call mpas_timer_stop("se implicit vert mix")

      call mpas_threading_barrier()

      call mpas_timer_start('se fini')
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCellsPtr)
         call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
         call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)

         call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)
         call mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', normalGMBolusVelocity)
         call mpas_pool_get_array(diagnosticsPool, 'velocityX', velocityX)
         call mpas_pool_get_array(diagnosticsPool, 'velocityY', velocityY)
         call mpas_pool_get_array(diagnosticsPool, 'velocityZ', velocityZ)
         call mpas_pool_get_array(diagnosticsPool, 'velocityZonal', velocityZonal)
         call mpas_pool_get_array(diagnosticsPool, 'velocityMeridional', velocityMeridional)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSH', gradSSH)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHX', gradSSHX)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHY', gradSSHY)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHZ', gradSSHZ)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHZonal', gradSSHZonal)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHMeridional', gradSSHMeridional)

         call mpas_pool_get_array(diagnosticsPool, 'surfaceVelocity', surfaceVelocity)
         call mpas_pool_get_array(diagnosticsPool, 'SSHGradient', SSHGradient)

         call mpas_pool_get_dimension(diagnosticsPool, 'index_surfaceVelocityZonal', indexSurfaceVelocityZonal)
         call mpas_pool_get_dimension(diagnosticsPool, 'index_surfaceVelocityMeridional', indexSurfaceVelocityMeridional)
         call mpas_pool_get_dimension(diagnosticsPool, 'index_SSHGradientZonal', indexSSHGradientZonal)
         call mpas_pool_get_dimension(diagnosticsPool, 'index_SSHGradientMeridional', indexSSHGradientMeridional)

         nCells = nCellsPtr
         nEdges = nEdgesPtr

         if (config_prescribe_velocity) then
            !$omp do schedule(runtime)
            do iEdge = 1, nEdges
               normalVelocityNew(:, iEdge) = normalVelocityCur(:, iEdge)
            end do
            !$omp end do
         end if

         if (config_prescribe_thickness) then
            !$omp do schedule(runtime)
            do iCell = 1, nCells
               layerThicknessNew(:, iCell) = layerThicknessCur(:, iCell)
            end do
            !$omp end do
         end if

         call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, tracersPool, 2)

         ! Update the effective desnity in land ice if we're coupling to land ice
         call ocn_effective_density_in_land_ice_update(meshPool, forcingPool, statePool, scratchPool, err)

         ! Compute normalGMBolusVelocity; it will be added to normalVelocity in Stage 2 of the next cycle.
         if (config_use_standardGM) then
            call ocn_gm_compute_Bolus_velocity(diagnosticsPool, meshPool, scratchPool)
         end if

         call mpas_timer_start('se final mpas reconstruct', .false.)
         call mpas_threading_barrier()
         !$omp master
         call mpas_reconstruct(meshPool, normalVelocityNew,  &
                          velocityX, velocityY, velocityZ,   &
                          velocityZonal, velocityMeridional, &
                          includeHalos = .true.)

         call mpas_reconstruct(meshPool, gradSSH,          &
                          gradSSHX, gradSSHY, gradSSHZ,    &
                          gradSSHZonal, gradSSHMeridional, &
                          includeHalos = .true.)
         !$omp end master
         call mpas_threading_barrier()
         call mpas_timer_stop('se final mpas reconstruct')

         !$omp do schedule(runtime)
         do iCell = 1, nCells
            surfaceVelocity(indexSurfaceVelocityZonal, iCell) = velocityZonal(1, iCell)
            surfaceVelocity(indexSurfaceVelocityMeridional, iCell) = velocityMeridional(1, iCell)

            SSHGradient(indexSSHGradientZonal, iCell) = gradSSHZonal(iCell)
            SSHGradient(indexSSHGradientMeridional, iCell) = gradSSHMeridional(iCell)
         end do
         !$omp end do

         call ocn_time_average_coupled_accumulate(diagnosticsPool, statePool, forcingPool, 2)
         call mpas_threading_barrier()

         if (config_use_standardGM) then
            call ocn_reconstruct_gm_vectors(diagnosticsPool, meshPool)
         end if

         block => block % next
      end do

      if (trim(config_land_ice_flux_mode) == 'coupled') then
         call mpas_timer_start("se effective density halo")
         call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
         call mpas_pool_get_field(statePool, 'effectiveDensityInLandIce', effectiveDensityField, 2)
         call mpas_dmpar_exch_halo_field(effectiveDensityField)
         call mpas_timer_stop("se effective density halo")
      end if

      call mpas_timer_stop('se fini')
      call mpas_timer_stop("se timestep")

      deallocate(n_bcl_iter)

   end subroutine ocn_time_integrator_split!}}}

!***********************************************************************
!
!  routine ocn_time_integration_split_init
!
!> \brief   Initialize split-explicit time stepping within MPAS-Ocean core
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  This routine initializes variables required for the split-explicit time
!>  stepper.
!
!-----------------------------------------------------------------------
   subroutine ocn_time_integration_split_init(domain)!{{{
   ! Initialize splitting variables

      type (domain_type), intent(inout) :: domain

      integer :: i, iCell, iEdge, iVertex, k
      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: statePool, meshPool, tracersPool

      integer :: iTracer, cell, cell1, cell2
      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: cellsOnEdge
      real (kind=RKIND) :: normalThicknessFluxSum, layerThicknessSum, layerThicknessEdge1
      real (kind=RKIND), dimension(:), pointer :: refBottomDepth, normalBarotropicVelocity

      real (kind=RKIND), dimension(:,:), pointer :: layerThickness
      real (kind=RKIND), dimension(:,:), pointer :: normalBaroclinicVelocity, normalVelocity
      integer, pointer :: nVertLevels, nCells, nEdges
      character (len=StrKIND), pointer :: config_time_integrator, config_btr_dt, config_dt
      logical, pointer :: config_filter_btr_mode, config_do_restart

      type (mpas_time_type) :: nowTime
      type (mpas_timeInterval_type) :: fullTimeStep, barotropicTimeStep, remainder, zeroInterval

      integer :: iErr

      call mpas_pool_get_config(domain % configs, 'config_do_restart', config_do_restart)

      ! Determine the number of barotropic subcycles based on the ratio of time steps
      call mpas_pool_get_config(domain % configs, 'config_time_integrator', config_time_integrator)
      call mpas_pool_get_config(domain % configs, 'config_btr_dt', config_btr_dt)
      call mpas_pool_get_config(domain % configs, 'config_dt', config_dt)

      nowTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
      call mpas_set_timeInterval( zeroInterval, S=0 )

      call mpas_set_timeInterval( fullTimeStep , timeString=config_dt )
      call mpas_set_timeInterval( barotropicTimeStep, timeString=config_btr_dt )

      call mpas_interval_division( nowTime, fullTimeStep, barotropicTimeStep, nBtrSubcycles, remainder )

      if ( remainder > zeroInterval ) then
         nBtrSubcycles = nBtrSubcycles + 1
      end if

      if (trim(config_time_integrator) == 'split_explicit') then
         write(stderrUnit, *) ''
         write(stderrUnit, *) '*********************************************************************************************************'
         write(stderrUnit, *) 'INFO: The split explicit time integration is configured to use: ', nBtrSubcycles, ' barotropic subcycles'
         write(stderrUnit, *) '*********************************************************************************************************'
         write(stderrUnit, *) ''
      end if

      if ( .not. config_do_restart ) then
         ! Initialize z-level mesh variables from h, read in from input file.
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_config(block % configs, 'config_time_integrator', config_time_integrator)
            call mpas_pool_get_config(block % configs, 'config_filter_btr_mode', config_filter_btr_mode)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

            call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
            call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
            call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)

            call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, 1)
            call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocity, 1)
            call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', normalBaroclinicVelocity, 1)

            call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
            call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
            call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)

            ! Compute barotropic velocity at first timestep
            ! This is only done upon start-up.
            if (trim(config_time_integrator) == 'unsplit_explicit') then
               call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocity)

               do iEdge = 1, nEdges
                  normalBarotropicVelocity(iEdge) = 0.0_RKIND
                  normalBaroclinicVelocity(:, iEdge) = normalVelocity(:, iEdge)
               end do

            elseif (trim(config_time_integrator) == 'split_explicit') then

               write(stderrUnit, *) ''
               write(stderrUnit, *) '*********************************************************************************************************'
               write(stderrUnit, *) 'INFO: The split explicit time integration is configured to use: ', nBtrSubcycles, ' barotropic subcycles'
               write(stderrUnit, *) '*********************************************************************************************************'
               write(stderrUnit, *) ''

               if (config_filter_btr_mode) then
                  do iCell = 1, nCells
                     layerThickness(1,iCell) = refBottomDepth(1)
                  enddo
               endif

               do iEdge = 1, nEdges
                  cell1 = cellsOnEdge(1,iEdge)
                  cell2 = cellsOnEdge(2,iEdge)

                  ! normalBarotropicVelocity = sum(h*u)/sum(h) on each edge
                  ! ocn_diagnostic_solve has not yet been called, so compute hEdge
                  ! just for this edge.

                  ! thicknessSum is initialized outside the loop because on land boundaries
                  ! maxLevelEdgeTop=0, but I want to initialize thicknessSum with a
                  ! nonzero value to avoid a NaN.
                  layerThicknessEdge1 = 0.5_RKIND*( layerThickness(1,cell1) + layerThickness(1,cell2) )
                  normalThicknessFluxSum = layerThicknessEdge1 * normalVelocity(1,iEdge)
                  layerThicknessSum = layerThicknessEdge1

                  do k=2, maxLevelEdgeTop(iEdge)
                     ! ocn_diagnostic_solve has not yet been called, so compute hEdge
                     ! just for this edge.
                     layerThicknessEdge1 = 0.5_RKIND*( layerThickness(k,cell1) + layerThickness(k,cell2) )

                     normalThicknessFluxSum = normalThicknessFluxSum &
                        + layerThicknessEdge1 * normalVelocity(k,iEdge)
                     layerThicknessSum = layerThicknessSum + layerThicknessEdge1

                  enddo
                  normalBarotropicVelocity(iEdge) = normalThicknessFluxSum / layerThicknessSum

                  ! normalBaroclinicVelocity(k,iEdge) = normalVelocity(k,iEdge) - normalBarotropicVelocity(iEdge)
                  do k = 1, maxLevelEdgeTop(iEdge)
                     normalBaroclinicVelocity(k,iEdge) = normalVelocity(k,iEdge) - normalBarotropicVelocity(iEdge)
                  enddo

                  ! normalBaroclinicVelocity=0, normalVelocity=0 on land cells
                  do k = maxLevelEdgeTop(iEdge)+1, nVertLevels
                     normalBaroclinicVelocity(k,iEdge) = 0.0_RKIND
                     normalVelocity(k,iEdge) = 0.0_RKIND
                  enddo
               enddo

               if (config_filter_btr_mode) then
                  ! filter normalBarotropicVelocity out of initial condition

                   normalVelocity(:,:) = normalBaroclinicVelocity(:,:)
                   normalBarotropicVelocity(:) = 0.0_RKIND

               endif

            endif

         block => block % next
         end do
      end if

   end subroutine ocn_time_integration_split_init!}}}

end module ocn_time_integration_split

! vim: foldmethod=marker
