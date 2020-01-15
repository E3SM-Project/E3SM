! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_si
!
!> \brief MPAS ocean semi-implicit time integration scheme
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011 (split explicit base code)
!
!> \author Hyun-Gyu Kang (Oak Ridge National Laboratory)
!> \date   September 2019 (semi-implicit code)
!> \details
!>  This module contains the routine for the semi-implicit
!>  time integration scheme based on the split-explicit code.
!>  Only stage 2 (barotropic mode) is changed from the explicit
!>  subcycling scheme to the semi-implicit scheme.
!-----------------------------------------------------------------------


module ocn_time_integration_si

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar
   use mpas_vector_reconstruction
   use mpas_spline_interpolation
   use mpas_timer
   use mpas_threading
   use mpas_timekeeping
   use mpas_log

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

   public :: ocn_time_integrator_si, ocn_time_integration_si_init
   public :: ocn_time_integrator_si_preconditioner

   character (len=*), parameter :: iterGroupName = 'iterFields'
   character (len=*), parameter :: finalBtrGroupName = 'finalBtrFields'
   integer :: nBtrSubcycles

   ! Global variables for the semi-implicit time stepper ---------------
   real (kind=RKIND), allocatable,dimension(:,:) :: prec_ivmat
   real (kind=RKIND) :: total_num_cells,mean_num_cells,area_mean
   integer :: nPrecVec

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_si
!
!> \brief MPAS ocean split explicit time integration scheme
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011 (split explicit base code)
!> \author Hyun-Gyu Kang (Oak Ridge National Laboratory)
!> \date   JAN 2019 (semi-implicit code)
!> \details
!>  This routine integrates a master time step (dt) using a
!>  semi-implicit time integrator.
!
!-----------------------------------------------------------------------

    subroutine ocn_time_integrator_si(domain, dt)!{{{
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
      integer :: iCell, i,k,j, iEdge, cell1, cell2, split_implicit_step, split, &
                 eoe, oldBtrSubcycleTime, newBtrSubcycleTime, uPerpTime, BtrCorIter, &
                 stage1_tend_time,iter,iter_u,iter_max
      integer, dimension(:), allocatable :: n_bcl_iter
      type (block_type), pointer :: block

      real (kind=RKIND) :: normalThicknessFluxSum, thicknessSum,thicknessSumCur,thicknessSumLag,thicknessSumMid,  &
                           flux,flux1,flux2, sshEdge,sshEdgeCur,sshEdgeLag,sshEdgeMid, hEdge1, &
                           CoriolisTerm, normalVelocityCorrection, temp, temp_h, coef, barotropicThicknessFlux_coeff, &
                           sshDiffCur,sshDiffNew,sshDiffLag,sshDiffMid
      real (kind=RKIND) :: fluxb1,fluxb2,fluxAx,sshTendb1,sshTendb2,sshTendAx

      integer :: useVelocityCorrection, err
      real (kind=RKIND), dimension(:,:), pointer :: &
                 vertViscTopOfEdge, vertDiffTopOfCell
      real (kind=RKIND), dimension(:,:,:), pointer :: tracersGroup
      real (kind=RKIND), dimension(:), allocatable:: uTemp
      real (kind=RKIND), dimension(:), pointer :: btrvel_temp
      type (field1DReal), pointer :: btrvel_tempField
      logical :: activeTracersOnly ! if true only compute tendencies for active tracers  
      integer :: tsIter
      integer :: edgeHaloComputeCounter, cellHaloComputeCounter
      integer :: neededHalos

      ! Config options
      character (len=StrKIND), pointer :: config_time_integrator
      integer, pointer :: config_n_bcl_iter_mid, config_n_bcl_iter_beg, config_n_bcl_iter_end
      integer, pointer :: config_n_ts_iter, config_btr_subcycle_loop_factor
      integer, pointer :: config_n_btr_cor_iter, config_num_halos
      logical, pointer :: config_use_standardGM

      logical, pointer :: config_use_freq_filtered_thickness, config_btr_solve_SSH2, config_filter_btr_mode
      logical, pointer :: config_vel_correction, config_prescribe_velocity, config_prescribe_thickness
      logical, pointer :: config_disable_thick_all_tend
      logical, pointer :: config_disable_vel_all_tend
      logical, pointer :: config_disable_tr_all_tend
      logical, pointer :: config_use_cvmix_kpp
      logical, pointer :: config_use_tracerGroup
      logical, pointer :: config_compute_active_tracer_budgets
      logical, pointer :: config_use_tidal_potential_forcing

      character (len=StrKIND), pointer :: config_land_ice_flux_mode

      real (kind=RKIND), pointer :: config_mom_del4, config_btr_gam1_velWt1, config_btr_gam2_SSHWt1
      real (kind=RKIND), pointer :: config_btr_gam3_velWt2
      real (kind=RKIND), pointer :: config_self_attraction_and_loading_beta

      ! Configs for the semi-implicit solver 
      character (len=StrKIND), pointer :: config_btr_si_preconditioner
      real (kind=RKIND), pointer :: config_btr_si_tolerance
      integer, pointer :: config_n_btr_si_outer_iter

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
      real (kind=RKIND), dimension(:), pointer :: sshSubcycleCurWithTides, sshSubcycleNewWithTides
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
      real (kind=RKIND), dimension(:), pointer :: tidalPotentialEta


      ! Semi-implicit Array Pointers
      real (kind=RKIND), dimension(:), pointer :: CGvec_r0,CGvec_r00,CGvec_r1,CGvec_rh0,CGvec_rh1,CGvec_ph0,CGvec_ph1
      real (kind=RKIND), dimension(:), pointer :: CGvec_v0,CGvec_v1,CGvec_s0,CGvec_s1,CGvec_sh0,CGvec_sh1
      real (kind=RKIND), dimension(:), pointer :: CGvec_w0,CGvec_w1,CGvec_wh0,CGvec_wh1
      real (kind=RKIND), dimension(:), pointer :: CGvec_t0,CGvec_t1,CGvec_q0,CGvec_q1,CGvec_qh0
      real (kind=RKIND), dimension(:), pointer :: CGvec_z0,CGvec_z1,CGvec_zh0,CGvec_zh1
      real (kind=RKIND), dimension(:), pointer :: CGvec_y0
      real (kind=RKIND), dimension(:), pointer :: barotropicCoriolisTerm
      ! Semi-implicit variables
      real (kind=RKIND), dimension(5) :: CGcst_allreduce,CGcst_allreduce_global
      real (kind=RKIND) :: CGcst_r00r0,CGcst_r00w0,CGcst_r00r1,CGcst_r00w1
      real (kind=RKIND) :: CGcst_r00r0_global,CGcst_r00w0_global,CGcst_r00r1_global,CGcst_r00w1_global

      real (kind=RKIND) :: CGcst_r0r0,CGcst_r1r1
      real (kind=RKIND) :: CGcst_r0r0_global,CGcst_r1r1_global

      real (kind=RKIND) :: CGcst_r00s0,CGcst_r00z0,CGcst_q0y0,CGcst_y0y0
      real (kind=RKIND) :: CGcst_r00s0_global,CGcst_r00z0_global,CGcst_q0y0_global,CGcst_y0y0_global

      real (kind=RKIND) :: CGcst_alpha0,CGcst_alpha1,CGcst_beta0,CGcst_beta1,CGcst_omega0,CGcst_omega1
      real (kind=RKIND) :: sshCurArea,sshLagArea,sshTendA,sshTendB,sshTendC,temp1,temp2,crit_main,resid,wgt

      ! Diagnostics Field Pointers
      type (field2DReal), pointer :: normalizedRelativeVorticityEdgeField, divergenceField, relativeVorticityField
      type (field1DReal), pointer :: barotropicThicknessFluxField, boundaryLayerDepthField, effectiveDensityField
      ! tracer tendencies brought in here to normalize by new layer thickness
      real (kind=RKIND), dimension(:,:,:), pointer :: &
        activeTracerHorizontalAdvectionTendency,      &
        activeTracerVerticalAdvectionTendency,        &
        activeTracerSurfaceFluxTendency,              &
        activeTracerNonLocalTendency,                 &
        activeTracerHorizontalAdvectionEdgeFlux

      real (kind=RKIND), dimension(:,:), pointer :: &
        temperatureShortWaveTendency
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

      integer :: temp_mask

      real(kind=RKIND) :: alpha1,alpha2
      logical :: first_step = .true.


      dminfo = domain % dminfo

      call mpas_timer_start("si timestep")

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
      call mpas_pool_get_config(domain % configs, 'config_disable_vel_all_tend', config_disable_vel_all_tend)
      call mpas_pool_get_config(domain % configs, 'config_disable_thick_all_tend', config_disable_thick_all_tend)
      call mpas_pool_get_config(domain % configs, 'config_disable_tr_all_tend', config_disable_tr_all_tend)
      call mpas_pool_get_config(domain % configs, 'config_use_tidal_potential_forcing', config_use_tidal_potential_forcing)
      call mpas_pool_get_config(domain % configs, 'config_self_attraction_and_loading_beta',config_self_attraction_and_loading_beta)

      call mpas_pool_get_config(domain % configs, 'config_prescribe_velocity', config_prescribe_velocity)
      call mpas_pool_get_config(domain % configs, 'config_prescribe_thickness', config_prescribe_thickness)

      call mpas_pool_get_config(domain % configs, 'config_prescribe_velocity', config_prescribe_velocity)
      call mpas_pool_get_config(domain % configs, 'config_prescribe_thickness', config_prescribe_thickness)

      call mpas_pool_get_config(domain % configs, 'config_use_standardGM', config_use_standardGM)
      call mpas_pool_get_config(domain % configs, 'config_use_cvmix_kpp', config_use_cvmix_kpp)
      call mpas_pool_get_config(domain % configs, 'config_land_ice_flux_mode', config_land_ice_flux_mode)

      call mpas_pool_get_config(domain % configs, 'config_num_halos', config_num_halos)

      call mpas_pool_get_config(domain % configs, 'config_compute_active_tracer_budgets', config_compute_active_tracer_budgets)


      ! configurations for the semi-implicit barotropic mode solver
      call mpas_pool_get_config(domain % configs, 'config_btr_si_preconditioner', config_btr_si_preconditioner)
      call mpas_pool_get_config(domain % configs, 'config_btr_si_tolerance'     , config_btr_si_tolerance     )
      call mpas_pool_get_config(domain % configs, 'config_n_btr_si_outer_iter'  , config_n_btr_si_outer_iter  )

      ! Tolerance for main iteration -- !
      crit_main  = config_btr_si_tolerance * area_mean

      !--- Backward-Euler at the 1st time step only to remove spurious oscillation caused by CN method
      if ( first_step ) then
         first_step = .false.
         alpha1 = 1.0
         alpha2 = 0.0
      else
         alpha1 = 0.5
         alpha2 = 0.5
      endif

      allocate(n_bcl_iter(config_n_ts_iter))

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      !  Prep variables before first iteration
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call mpas_timer_start("si prep")
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

      call mpas_timer_stop("si prep")

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! BEGIN large iteration loop
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      n_bcl_iter = config_n_bcl_iter_mid
      n_bcl_iter(1) = config_n_bcl_iter_beg
      n_bcl_iter(config_n_ts_iter) = config_n_bcl_iter_end

      do split_implicit_step = 1, config_n_ts_iter

         if (config_disable_thick_all_tend .and. config_disable_vel_all_tend .and. config_disable_tr_all_tend) then
           exit ! don't compute in loop meant to update velocity, thickness, and tracers
         end if

         call mpas_timer_start('si loop')

         stage1_tend_time = min(split_implicit_step,2)

         call mpas_pool_get_subpool(domain % blocklist % structs, 'diagnostics', diagnosticsPool)

         ! ---  update halos for diagnostic ocean boundary layer depth
         if (config_use_cvmix_kpp) then
            call mpas_timer_start("si halo diag obd")
            call mpas_dmpar_field_halo_exch(domain, 'boundaryLayerDepth')
            call mpas_timer_stop("si halo diag obd")
         end if

         ! ---  update halos for diagnostic variables
         call mpas_timer_start("si halo diag")

         call mpas_dmpar_field_halo_exch(domain, 'normalizedRelativeVorticityEdge')
         if (config_mom_del4 > 0.0_RKIND) then
           call mpas_dmpar_field_halo_exch(domain, 'divergence')
           call mpas_dmpar_field_halo_exch(domain, 'relativeVorticity')
         end if
         call mpas_timer_stop("si halo diag")

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         !  Stage 1: Baroclinic velocity (3D) prediction, explicit with long timestep
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (config_use_freq_filtered_thickness) then
            call mpas_timer_start("si freq-filtered-thick computations")


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
            call mpas_timer_stop("si freq-filtered-thick computations")


            call mpas_timer_start("si freq-filtered-thick halo update")

            call mpas_dmpar_field_halo_exch(domain, 'tendHighFreqThickness')
            call mpas_dmpar_field_halo_exch(domain, 'tendLowFreqDivergence')

            call mpas_timer_stop("si freq-filtered-thick halo update")

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
         call mpas_timer_start("si bcl vel")

         call mpas_timer_start('si bcl vel tend')
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

           call ocn_tend_vel(tendPool, statePool, forcingPool, diagnosticsPool, meshPool, scratchPool, stage1_tend_time, dt)

           block => block % next
         end do
         call mpas_timer_stop('si bcl vel tend')

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! BEGIN baroclinic iterations on linear Coriolis term
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         do j=1,n_bcl_iter(split_implicit_step)

            ! Use this G coefficient to avoid an if statement within the iEdge loop.
            split = 1

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

               !$omp do schedule(runtime)
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

            call mpas_timer_start("si halo normalBaroclinicVelocity")
            call mpas_dmpar_field_halo_exch(domain, 'normalBaroclinicVelocity', timeLevel=2)
            call mpas_timer_stop("si halo normalBaroclinicVelocity")

            call mpas_timer_stop('bcl iters on linear Coriolis')

         end do  ! do j=1,config_n_bcl_iter

         call mpas_timer_start('si halo barotropicForcing')
         call mpas_dmpar_field_halo_exch(domain, 'barotropicForcing')
         call mpas_timer_stop('si halo barotropicForcing')

         call mpas_timer_stop("si bcl vel")
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! END baroclinic iterations on linear Coriolis term
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         !  Stage 2: Barotropic velocity (2D) prediction, semi-implicit
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         ! Stage 2.1 : Preparation of varibales before outer two iterations
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         call mpas_timer_start("si btr vel")

         cellHaloComputeCounter = config_num_halos
         edgeHaloComputeCounter = config_num_halos + 1

         if ( split_implicit_step == 1 ) then
   
            ! Initialize variables for barotropic subcycling
            call mpas_timer_start('btr vel si init')

            block => domain % blocklist
            do while (associated(block))
               call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
               call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
   
               call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
               call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
               call mpas_pool_get_subpool(block % structs, 'state', statePool)
               call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
               call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
               call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
   
               call mpas_pool_get_array(meshPool, 'nEdgesOnEdge', nEdgesOnEdge)
               call mpas_pool_get_array(meshPool, 'edgesOnEdge', edgesOnEdge)
               call mpas_pool_get_array(meshPool, 'weightsOnEdge', weightsOnEdge)
               call mpas_pool_get_array(meshPool, 'fEdge' , fEdge)
   
               call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)
               call mpas_pool_get_array(diagnosticsPool, 'barotropicCoriolisTerm',barotropicCoriolisTerm)
   
               call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
               call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, 1)
               call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleCur, 1)
               call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityCur, 1)
   
               nCells = nCellsArray(1)
               nEdges = nEdgesArray(1)
   
               if (config_filter_btr_mode) then
                  !$omp do schedule(runtime)
                  do iEdge = 1, nEdges
                     barotropicForcing(iEdge) = 0.0_RKIND
                  end do
                  !$omp end do
               endif
   
               do iCell = 1, nCells
                  ! sshSubcycleOld = sshOld
                  sshSubcycleCur(iCell) = sshCur(iCell)
               end do

               do iEdge = 1, nEdges
                  normalBarotropicVelocitySubcycleCur(iEdge) = normalBarotropicVelocityCur(iEdge)

                  ! Compute the barotropic Coriolis term, -f*uPerp
                  barotropicCoriolisTerm(iEdge) = 0.0_RKIND

                  do i = 1, nEdgesOnEdge(iEdge)
                     eoe = edgesOnEdge(i,iEdge)
                     barotropicCoriolisTerm(iEdge) =    barotropicCoriolisTerm(iEdge) + weightsOnEdge(i,iEdge) &
                                                      * normalBarotropicVelocityCur(eoe) * fEdge(eoe)
                  end do ! i
               end do ! iEdge
   
               block => block % next
            end do  ! block
   
            call mpas_timer_start("si halo init")
            call mpas_dmpar_exch_group_create(domain, iterGroupName)
            call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'barotropicCoriolisTerm', 1 )
            call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'normalBarotropicVelocitySubcycle', 1 )
            call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'sshSubcycle', 1 )
            call mpas_threading_barrier()
            call mpas_dmpar_exch_group_full_halo_exch(domain, iterGroupName)
            call mpas_dmpar_exch_group_destroy(domain, iterGroupName)
            call mpas_timer_stop("si halo init")
   
            call mpas_timer_stop('btr vel si init')

         else ! split_implicit_step

            call mpas_timer_start('btr vel si init')

            block => domain % blocklist
            do while (associated(block))
               call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
               call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
   
               call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
               call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
               call mpas_pool_get_subpool(block % structs, 'state', statePool)
               call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
               call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
               call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
   
               call mpas_pool_get_array(meshPool, 'nEdgesOnEdge', nEdgesOnEdge)
               call mpas_pool_get_array(meshPool, 'edgesOnEdge', edgesOnEdge)
               call mpas_pool_get_array(meshPool, 'weightsOnEdge', weightsOnEdge)
               call mpas_pool_get_array(meshPool, 'fEdge' , fEdge)
   
               call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)
               call mpas_pool_get_array(diagnosticsPool, 'barotropicCoriolisTerm',barotropicCoriolisTerm)
   
               call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
               call mpas_pool_get_array(statePool, 'ssh', sshNew, 2)
               call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, 1)
               call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityCur, 1)
               call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleCur, 1)
   
               nCells = nCellsArray(1)
               nEdges = nEdgesArray(2)
   
               if (config_filter_btr_mode) then
                  !$omp do schedule(runtime)
                  do iEdge = 1, nEdges
                     barotropicForcing(iEdge) = 0.0_RKIND
                  end do
                  !$omp end do
               endif
   
               do iCell = 1, nCells
                          sshCur(iCell) = 0.5*sshNew(iCell) + 0.5*sshCur(iCell)
                  sshSubcycleCur(iCell) = sshCur(iCell)
               end do

               do iEdge = 1, nEdges
                  normalBarotropicVelocityCur(iEdge) = normalBarotropicVelocitySubcycleCur(iEdge)

                  ! Compute the barotropic Coriolis term, -f*uPerp
                  barotropicCoriolisTerm(iEdge) = 0.0_RKIND

                  do i = 1, nEdgesOnEdge(iEdge)
                     eoe = edgesOnEdge(i,iEdge)
                     barotropicCoriolisTerm(iEdge) =    barotropicCoriolisTerm(iEdge) + weightsOnEdge(i,iEdge) &
                                                       * normalBarotropicVelocitySubcycleCur(eoe) * fEdge(eoe)
                  end do ! i
               end do ! iEdge
   
               block => block % next
            end do  ! block


            ! Subtract tidal potential from ssh, if needed
            !   Subtract the tidal potential from the current subcycle ssh and store and a work array.
            !   Then point sshSubcycleCur to the work array so the tidal potential terms are included
            !   in the grad operator inside the edge loop.
            if (config_use_tidal_potential_forcing) then

               block => domain % blocklist
               do while (associated(block))
                  call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdgesPtr)
                  call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

                  call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
                  call mpas_pool_get_subpool(block % structs, 'state', statePool)
                  call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
                  call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
                  call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)

                  call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, 1)
                  call mpas_pool_get_array(forcingPool, 'sshSubcycleCurWithTides', sshSubcycleCurWithTides)
                  call mpas_pool_get_array(forcingPool, 'tidalPotentialEta', tidalPotentialEta)
                  call mpas_pool_get_dimension(block % dimensions, 'nCells', nCellsPtr)

                  nCells = nCellsPtr
                  do iCell = 1, nCells
                    sshSubcycleCurWithTides(iCell) = sshSubcycleCur(iCell) - tidalPotentialEta(iCell) &
                                                   - config_self_attraction_and_loading_beta * sshSubcycleCur(iCell)
                  end do

                  call mpas_pool_get_array(forcingPool, 'sshSubcycleCurWithTides', sshSubcycleCur)

                  block => block % next
               end do  ! block

            endif !config_use_tidal_potential_forcing

            call mpas_timer_start("si halo init")
            call mpas_dmpar_exch_group_create(domain, iterGroupName)
            call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'barotropicCoriolisTerm', 1 )
            call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'normalBarotropicVelocity', 1 )
            call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'sshSubcycle', 1 )
            call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'ssh', 1 )
            call mpas_threading_barrier()
            call mpas_dmpar_exch_group_full_halo_exch(domain, iterGroupName)
            call mpas_dmpar_exch_group_destroy(domain, iterGroupName)
            call mpas_timer_stop("si halo init")

            call mpas_timer_stop('btr vel si init')

         endif ! split_implicit_step



         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         ! Stage 2.2 : Compute initial residual
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         call mpas_timer_start("si btr residual")

         ! SpMV -----------------------------------------------------------------------------------!

         block => domain % blocklist
         do while (associated(block))

            call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
            call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

            call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
            call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
            call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
            call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
            call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
            call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
            call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
            call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

            call mpas_pool_get_array(statePool, 'ssh', sshCur,1)
            call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityCur,1)

            call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)
            call mpas_pool_get_array(diagnosticsPool, 'barotropicCoriolisTerm',barotropicCoriolisTerm)
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_r0', CGvec_r0)
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_r00', CGvec_r00)

            nCells = nCellsArray( 1 )
            nEdges = nEdgesArray( 2 )

            do iCell = 1, nCells
               sshTendb1 = 0.0_RKIND
               sshTendb2 = 0.0_RKIND
               sshTendAx = 0.0_RKIND

               do i = 1, nEdgesOnCell(iCell)
                  iEdge = edgesOnCell(i, iCell)

                  cell1 = cellsOnEdge(1, iEdge)
                  cell2 = cellsOnEdge(2, iEdge)

                  ! Interpolation sshEdge
                  sshEdgeCur = 0.5_RKIND * (sshCur(cell1) + sshCur(cell2))

                  ! method 1, matches method 0 without pbcs, works with pbcs.
                  thicknessSumCur = sshEdgeCur + min(bottomDepth(cell1), bottomDepth(cell2))

                  ! nabla (ssh^0)
                  sshDiffCur = (  sshCur(cell2) -   sshCur(cell1)) / dcEdge(iEdge)

                  fluxb1 = thicknessSumCur * normalBarotropicVelocityCur(iEdge)
                  fluxb2 = thicknessSumCur * (alpha2*gravity*sshDiffCur + (-barotropicCoriolisTerm(iEdge)-barotropicForcing(iEdge)))
                  fluxAx = thicknessSumCur * sshDiffCur
  
                  sshTendb1 = sshTendb1 + edgeSignOnCell(i, iCell) * fluxb1 * dvEdge(iEdge)
                  sshTendb2 = sshTendb2 + edgeSignOnCell(i, iCell) * fluxb2 * dvEdge(iEdge)
                  sshTendAx = sshTendAx + edgeSignOnCell(i, iCell) * fluxAx * dvEdge(iEdge)
               end do ! i

               sshTendb1 = (1.0_RKIND/(gravity*dt*alpha1**2.0)) * sshTendb1
               sshTendb2 = (1.0_RKIND/(gravity   *alpha1     )) * sshTendb2

               sshCurArea = (1.0_RKIND/(gravity*dt**2.0*alpha1**2.0)) *   sshCur(iCell) * areaCell(iCell)

               CGvec_r0(iCell) = (-sshCurArea - sshTendb1 + sshTendb2)   &
                                -(-sshCurArea - sshTendAx) 
               CGvec_r00(iCell) = CGvec_r0(iCell)
            end do ! iCell

            block => block % next
         end do  ! block


         ! Preconditioning ------------------------------------------------------------------------!

         if ( trim(config_btr_si_preconditioner) == 'ras' ) then
            call mpas_timer_start("si halo r0")
            call mpas_dmpar_exch_group_create(domain, iterGroupName)
            call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'CGvec_r0', 1)
            call mpas_threading_barrier()
            call mpas_dmpar_exch_group_full_halo_exch(domain, iterGroupName)
            call mpas_dmpar_exch_group_destroy(domain, iterGroupName)
            call mpas_timer_stop("si halo r0")
         end if

         block => domain % blocklist
         do while (associated(block))

            call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
            call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

            call mpas_pool_get_array(diagnosticsPool, 'CGvec_r0', CGvec_r0)
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_rh0', CGvec_rh0)

            nCells = nCellsArray( 1 )
            nEdges = nEdgesArray( 2 )

            if ( trim(config_btr_si_preconditioner) == 'ras' ) then
               ! RAS preconditioning: Use BLAS for symmetric matrix-vector multiplication
               call DSYMV('U', nPrecVec, 1.0_RKIND, prec_ivmat, nPrecVec, CGvec_r0(1:nPrecVec), 1, 0.0_RKIND, CGvec_rh0(1:nPrecVec), 1)

            elseif ( trim(config_btr_si_preconditioner) == 'block_jacobi' ) then
               ! Block-Jacobi preconditioning: Use BLAS for symmetric matrix-vector multiplication
               call DSYMV('U', nPrecVec, 1.0_RKIND, prec_ivmat, nPrecVec, CGvec_r0(1:nPrecVec), 1, 0.0_RKIND, CGvec_rh0(1:nPrecVec), 1)

            elseif ( trim(config_btr_si_preconditioner) == 'jacobi' ) then
               ! Jacobi preconditioning 
               CGvec_rh0(1:nPrecVec) = CGvec_r0(1:nPrecVec) * prec_ivmat(1:nPrecVec,1)

            elseif ( trim(config_btr_si_preconditioner) == 'none' ) then
               ! No preconditioning 
               CGvec_rh0(1:nPrecVec) = CGvec_r0(1:nPrecVec)
            end if

            block => block % next
         end do  ! block

         call mpas_timer_start("si halo r0")
         call mpas_dmpar_exch_group_create(domain, iterGroupName)
         call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'CGvec_rh0', 1)
         call mpas_threading_barrier()
         call mpas_dmpar_exch_group_full_halo_exch(domain, iterGroupName)
         call mpas_dmpar_exch_group_destroy(domain, iterGroupName)
         call mpas_timer_stop("si halo r0")


         ! SpMV -----------------------------------------------------------------------------------!

         block => domain % blocklist
         do while (associated(block))

            call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
            call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
 
            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
 
            call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
            call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
            call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
            call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
            call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
            call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
            call mpas_pool_get_array(meshPool, 'refBottomDepthTopOfCell', refBottomDepthTopOfCell)
            call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
            call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
            call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
 
            call mpas_pool_get_array(statePool, 'ssh', sshCur,1)
 
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_r0' , CGvec_r0 )
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_r00', CGvec_r00)
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_rh0', CGvec_rh0)
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_w0' , CGvec_w0)
 
            nCells = nCellsArray( 1 )
            nEdges = nEdgesArray( 2 )
 
            CGcst_r00r0 = 0.0_RKIND
            CGcst_r00w0 = 0.0_RKIND
             
            do iCell = 1, nCells
 
               sshTendAx = 0.0_RKIND
  
               do i = 1, nEdgesOnCell(iCell)
                  iEdge = edgesOnCell(i, iCell)
   
                  cell1 = cellsOnEdge(1, iEdge)
                  cell2 = cellsOnEdge(2, iEdge)
   
                  ! Interpolation sshEdge
                  sshEdgeCur = 0.5_RKIND * (sshCur(cell1) + sshCur(cell2))
   
                  ! method 1, matches method 0 without pbcs, works with pbcs.
                  thicknessSumCur = sshEdgeCur + min(bottomDepth(cell1), bottomDepth(cell2))
   
                  ! nabla (ssh^0)
                  sshDiffCur = (CGvec_rh0(cell2)- CGvec_rh0(cell1)) / dcEdge(iEdge)
                      fluxAx = thicknessSumCur * sshDiffCur
                   sshTendAx = sshTendAx + edgeSignOnCell(i, iCell) * fluxAx * dvEdge(iEdge)
   
               end do ! i
  
               sshCurArea = (1.0_RKIND/(gravity*dt**2.0*alpha1**2.0)) * CGvec_rh0(iCell) * areaCell(iCell)
                  
               CGvec_w0(iCell) = -sshCurArea - sshTendAx 
  
               CGcst_r00r0 = CGcst_r00r0 + CGvec_r00(iCell) * CGvec_r0(iCell)
               CGcst_r00w0 = CGcst_r00w0 + CGvec_r00(iCell) * CGvec_w0(iCell)
 
            end do ! iCell
 
            block => block % next
         end do  ! block


         ! Preconditioning ------------------------------------------------------------------------!

         if ( trim(config_btr_si_preconditioner) == 'ras' ) then
            call mpas_timer_start("si halo r0")
            call mpas_dmpar_exch_group_create(domain, iterGroupName)
            call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'CGvec_w0', 1)
            call mpas_threading_barrier()
            call mpas_dmpar_exch_group_full_halo_exch(domain, iterGroupName)
            call mpas_dmpar_exch_group_destroy(domain, iterGroupName)
            call mpas_timer_stop("si halo r0")
         end if 

         block => domain % blocklist
         do while (associated(block))

            call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
            call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
 
            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
 
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_w0' , CGvec_w0)
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_wh0', CGvec_wh0)
 
            nCells = nCellsArray( 1 )
            nEdges = nEdgesArray( 2 )
 
            if ( trim(config_btr_si_preconditioner) == 'ras' ) then
               ! RAS preconditioning: Use BLAS for symmetric matrix-vector multiplication
               call DSYMV('U', nPrecVec, 1.0_RKIND, prec_ivmat, nPrecVec, CGvec_w0(1:nPrecVec), 1, 0.0_RKIND, CGvec_wh0(1:nPrecVec), 1)
 
            elseif ( trim(config_btr_si_preconditioner) == 'block_jacobi' ) then
               ! Block-Jacobi preconditioning: Use BLAS for symmetric matrix-vector multiplication
               call DSYMV('U', nPrecVec, 1.0_RKIND, prec_ivmat, nPrecVec, CGvec_w0(1:nPrecVec), 1, 0.0_RKIND, CGvec_wh0(1:nPrecVec), 1)
 
            elseif ( trim(config_btr_si_preconditioner) == 'jacobi' ) then
               ! Jacobi preconditioning 
               CGvec_wh0(1:nPrecVec) = CGvec_w0(1:nPrecVec) * prec_ivmat(1:nPrecVec,1)
 
            elseif ( trim(config_btr_si_preconditioner) == 'none' ) then
               ! No preconditioning 
               CGvec_wh0(1:nPrecVec) = CGvec_w0(1:nPrecVec)
            end if
 
            block => block % next
         end do  ! block

         call mpas_timer_start("si halo r0")
         call mpas_dmpar_exch_group_create(domain, iterGroupName)
         call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'CGvec_wh0', 1)
         call mpas_threading_barrier()
         call mpas_dmpar_exch_group_full_halo_exch(domain, iterGroupName)
         call mpas_dmpar_exch_group_destroy(domain, iterGroupName)
         call mpas_timer_stop("si halo r0")


         ! SpMV -----------------------------------------------------------------------------------!

         block => domain % blocklist
         do while (associated(block))

            call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
            call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
 
            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
 
            call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
            call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
            call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
            call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
            call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
            call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
            call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
            call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
 
            call mpas_pool_get_array(statePool, 'ssh', sshCur,1)
 
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_wh0', CGvec_wh0)
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_t0' , CGvec_t0)
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_ph0', CGvec_ph0)
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_v0' , CGvec_v0)
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_s0' , CGvec_s0)
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_sh0', CGvec_sh0)
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_z0' , CGvec_z0)
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_zh0', CGvec_zh0)
 
            nCells = nCellsArray( 1 )
            nEdges = nEdgesArray( 2 )
 
            do iCell = 1, nCells
 
               sshTendAx = 0.0_RKIND
  
               do i = 1, nEdgesOnCell(iCell)

                  iEdge = edgesOnCell(i, iCell)
   
                  cell1 = cellsOnEdge(1, iEdge)
                  cell2 = cellsOnEdge(2, iEdge)
   
                  ! Interpolation sshEdge
                  sshEdgeCur = 0.5_RKIND * (sshCur(cell1) + sshCur(cell2))
   
                  ! method 1, matches method 0 without pbcs, works with pbcs.
                  thicknessSumCur = sshEdgeCur + min(bottomDepth(cell1), bottomDepth(cell2))
   
                  ! nabla (ssh^0)
                  sshDiffCur = (CGvec_wh0(cell2)- CGvec_wh0(cell1)) / dcEdge(iEdge)
   
                  !--------------------------------------------------------------!
                     fluxAx = thicknessSumCur * sshDiffCur
     
                  sshTendAx = sshTendAx + edgeSignOnCell(i, iCell) * fluxAx * dvEdge(iEdge)
                  !--------------------------------------------------------------!
   
               end do ! i
 
               sshCurArea = (1.0_RKIND/(gravity*dt**2.0*alpha1**2.0)) * CGvec_wh0(iCell) * areaCell(iCell)
                
               CGvec_t0(iCell) = -sshCurArea - sshTendAx 
            
               CGvec_ph0(iCell) = 0.0_RKIND
               CGvec_v0(iCell)  = 0.0_RKIND
               CGvec_sh0(iCell) = 0.0_RKIND
               CGvec_z0(iCell)  = 0.0_RKIND
               CGvec_zh0(iCell) = 0.0_RKIND
               CGvec_v0(iCell)  = 0.0_RKIND
               CGvec_s0(iCell)  = 0.0_RKIND
                 
            end do ! iCell
 
            block => block % next
         end do  ! block


         CGcst_allreduce(:) = 0.0_RKIND
         CGcst_allreduce(1) = CGcst_r00r0
         CGcst_allreduce(2) = CGcst_r00w0

         ! Global sum across CPUs
         call mpas_timer_start("si reduction r0")
         call mpas_dmpar_sum_real_array(dminfo, 5, CGcst_allreduce, CGcst_allreduce_global)
         call mpas_timer_stop ("si reduction r0")

         CGcst_r00r0_global = CGcst_allreduce_global(1)
         CGcst_r00w0_global = CGcst_allreduce_global(2)

         CGcst_alpha0 = CGcst_r00r0_global / CGcst_r00w0_global
         CGcst_beta0  = 0.0_RKIND
         CGcst_omega0 = 0.0_RKIND

         call mpas_timer_stop("si btr residual")


         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         ! Stage 2.3 : Outer iterations - lagged values are sufficiently up to date
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         call mpas_timer_start("si btr iteration")

         !**************************************************************!
         do iter = 1,config_n_btr_si_outer_iter
         !**************************************************************!

            block => domain % blocklist
            do while (associated(block))
               call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
               call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
    
               call mpas_pool_get_subpool(block % structs, 'tend'       , tendPool       )
               call mpas_pool_get_subpool(block % structs, 'mesh'       , meshPool       )
               call mpas_pool_get_subpool(block % structs, 'state'      , statePool      )
               call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
    
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_r0' , CGvec_r0 )
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_r00', CGvec_r00)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_rh0', CGvec_rh0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_w0' , CGvec_w0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_wh0', CGvec_wh0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_t0' , CGvec_t0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_ph0', CGvec_ph0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_ph1', CGvec_ph1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_v0' , CGvec_v0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_s0' , CGvec_s0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_s1' , CGvec_s1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_sh0', CGvec_sh0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_sh1', CGvec_sh1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_q0' , CGvec_q0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_qh0', CGvec_qh0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_z0' , CGvec_z0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_z1' , CGvec_z1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_zh0', CGvec_zh0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_y0' , CGvec_y0)
    
               nCells = nCellsArray(1)
               nEdges = nEdgesArray(2)
   
               do iCell = 1, nCells
                  CGvec_ph1(iCell) = CGvec_rh0(iCell) + CGcst_beta0  * (CGvec_ph0(iCell)-CGcst_omega0*CGvec_sh0(iCell))
                  CGvec_s1(iCell)  = CGvec_w0(iCell)  + CGcst_beta0  * (CGvec_s0(iCell)-CGcst_omega0*CGvec_z0(iCell))
    	          CGvec_sh1(iCell) = CGvec_wh0(iCell) + CGcst_beta0  * (CGvec_sh0(iCell)-CGcst_omega0*CGvec_zh0(iCell))
    	          CGvec_z1(iCell)  = CGvec_t0(iCell)  + CGcst_beta0  * (CGvec_z0(iCell)-CGcst_omega0*CGvec_v0(iCell))
    	          CGvec_q0(iCell)  = CGvec_r0(iCell)  - CGcst_alpha0 * CGvec_s1(iCell)
                  CGvec_qh0(iCell) = CGvec_rh0(iCell) - CGcst_alpha0 * CGvec_sh1(iCell)
    	          CGvec_y0(iCell)  = CGvec_w0(iCell)  - CGcst_alpha0 * CGvec_z1(iCell)
               end do ! iCell
   
   
               ! Begin reduction --------------------------------------------------------------------!
   
               CGcst_q0y0  = 0.0_RKIND
               CGcst_y0y0  = 0.0_RKIND
               CGcst_r00r0 = 0.0_RKIND
   
               do iCell = 1,nCells
                  CGcst_q0y0  = CGcst_q0y0  + CGvec_q0(iCell)  * CGvec_y0(iCell)               
                  CGcst_y0y0  = CGcst_y0y0  + CGvec_y0(iCell)  * CGvec_y0(iCell)               
                  CGcst_r00r0 = CGcst_r00r0 + CGvec_r00(iCell) * CGvec_r0(iCell)
               end do
   
               block => block % next
            end do  ! block
   

            CGcst_allreduce(:) = 0.0_RKIND
            CGcst_allreduce(1) = CGcst_q0y0
            CGcst_allreduce(2) = CGcst_y0y0
            CGcst_allreduce(3) = CGcst_r00r0
   
            ! Global sum across CPUs
            call mpas_timer_start("si reduction iter")
            call mpas_dmpar_sum_real_array(dminfo, 5, CGcst_allreduce, CGcst_allreduce_global)
            call mpas_timer_stop("si reduction iter")
   
            CGcst_q0y0_global  = CGcst_allreduce_global(1)
            CGcst_y0y0_global  = CGcst_allreduce_global(2)
            CGcst_r00r0_global = CGcst_allreduce_global(3)
   
   
            ! Preconditioning ------------------------------------------------------------------------!
   
            if ( trim(config_btr_si_preconditioner) == 'ras' ) then
               call mpas_timer_start("si halo iter")
               call mpas_dmpar_exch_group_create(domain, iterGroupName)
               call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'CGvec_z1', 1)
               call mpas_threading_barrier()
               call mpas_dmpar_exch_group_full_halo_exch(domain, iterGroupName)
               call mpas_dmpar_exch_group_destroy(domain, iterGroupName)
               call mpas_timer_stop("si halo iter")
            end if
   
   
            block => domain % blocklist
            do while (associated(block))
   
               call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
               call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
    
               call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
               call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
               call mpas_pool_get_subpool(block % structs, 'state', statePool)
               call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
    
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_z1' , CGvec_z1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_zh1', CGvec_zh1)
    
               nCells = nCellsArray(1)
               nEdges = nEdgesArray(2)

               if ( trim(config_btr_si_preconditioner) == 'ras' ) then
                  ! RAS preconditioning: Use BLAS for symmetric matrix-vector multiplication
                  call DSYMV('U', nPrecVec, 1.0_RKIND, prec_ivmat, nPrecVec, CGvec_z1(1:nPrecVec), 1, 0.0_RKIND, CGvec_zh1(1:nPrecVec), 1)
   
               elseif ( trim(config_btr_si_preconditioner) == 'block_jacobi' ) then
                  ! Block-Jacobi preconditioning: Use BLAS for symmetric matrix-vector multiplication
                  call DSYMV('U', nPrecVec, 1.0_RKIND, prec_ivmat, nPrecVec, CGvec_z1(1:nPrecVec), 1, 0.0_RKIND, CGvec_zh1(1:nPrecVec), 1)
   
               elseif ( trim(config_btr_si_preconditioner) == 'jacobi' ) then
                  ! Jacobi preconditioning 
                  CGvec_zh1(1:nPrecVec) = CGvec_z1(1:nPrecVec) * prec_ivmat(1:nPrecVec,1)
   
               elseif ( trim(config_btr_si_preconditioner) == 'none' ) then
                  ! No preconditioning 
                  CGvec_zh1(1:nPrecVec) = cgvec_z1(1:nprecvec)
               end if
   
               block => block % next
            end do  ! block
   
            call mpas_timer_start("si halo iter")
            call mpas_dmpar_exch_group_create(domain, iterGroupName)
            call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'CGvec_zh1', 1)
            call mpas_threading_barrier()
            call mpas_dmpar_exch_group_full_halo_exch(domain, iterGroupName)
            call mpas_dmpar_exch_group_destroy(domain, iterGroupName)
            call mpas_timer_stop("si halo iter")
   
   
            ! SpMV -----------------------------------------------------------------------------------!

            block => domain % blocklist
            do while (associated(block))
   
               call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
               call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
    
               call mpas_pool_get_subpool(block % structs, 'tend'       , tendPool       )
               call mpas_pool_get_subpool(block % structs, 'mesh'       , meshPool       )
               call mpas_pool_get_subpool(block % structs, 'state'      , statePool      )
               call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
    
               call mpas_pool_get_array(meshPool, 'nEdgesOnCell',            nEdgesOnCell           )
               call mpas_pool_get_array(meshPool, 'edgesOnCell',             edgesOnCell            )
               call mpas_pool_get_array(meshPool, 'cellsOnEdge',             cellsOnEdge            )
               call mpas_pool_get_array(meshPool, 'dcEdge',                  dcEdge                 )
               call mpas_pool_get_array(meshPool, 'bottomDepth',             bottomDepth            )
               call mpas_pool_get_array(meshPool, 'edgeSignOnCell',          edgeSignOnCell         )
               call mpas_pool_get_array(meshPool, 'dvEdge',                  dvEdge                 )
               call mpas_pool_get_array(meshPool, 'areaCell',                areaCell               )
    
               call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
               call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, 1)
               call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleNew, 2)
    
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_r1' , CGvec_r1 )
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_rh1', CGvec_rh1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_w1' , CGvec_w1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_wh0', CGvec_wh0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_t0' , CGvec_t0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_ph1', CGvec_ph1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_v1' , CGvec_v1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_q0' , CGvec_q0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_qh0', CGvec_qh0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_zh1', CGvec_zh1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_y0' , CGvec_y0)
    
               nCells = nCellsArray(1)
               nEdges = nEdgesArray(2)
   
               do iCell = 1, nCells
                  sshTendAx = 0.0_RKIND
   
                  do i = 1, nEdgesOnCell(iCell)
                     iEdge = edgesOnCell(i, iCell)
     
                     cell1 = cellsOnEdge(1, iEdge)
                     cell2 = cellsOnEdge(2, iEdge)
     
                     ! Interpolation sshEdge
                     sshEdgeLag = 0.5_RKIND * (sshCur(cell1) + sshCur(cell2))
     
                     ! method 1, matches method 0 without pbcs, works with pbcs.
                     thicknessSumLag = sshEdgeLag + min(bottomDepth(cell1), bottomDepth(cell2))
     
                     ! nabla (ssh^0)
                     sshDiffNew = (CGvec_zh1(cell2)-CGvec_zh1(cell1)) / dcEdge(iEdge)
     
                     fluxAx = thicknessSumLag * sshDiffNew
     
                     sshTendAx = sshTendAx + edgeSignOnCell(i, iCell) * fluxAx * dvEdge(iEdge)
                  end do ! i
                   
                  sshLagArea = (1.0_RKIND/(gravity*dt**2.0*alpha1**2.0)) * CGvec_zh1(iCell) * areaCell(iCell)
   
                  CGvec_v1(iCell) = -sshLagArea - sshTendAx
       
               end do ! iCell
      
               ! End reduction -----------------------------------------------------------------------!
   
               ! Omega1
               CGcst_omega0 = CGcst_q0y0_global / CGcst_y0y0_global
   
               do iCell = 1,nCells
                  sshSubcycleNew(iCell) = sshSubcycleCur(iCell) + CGcst_alpha0 * CGvec_ph1(iCell) &
                                                                + CGcst_omega0 * CGvec_qh0(iCell)
                  CGvec_r1(iCell)  = CGvec_q0(iCell)  - CGcst_omega0 * CGvec_y0(iCell)
                  CGvec_rh1(iCell) = CGvec_qh0(iCell) - CGcst_omega0 * (CGvec_wh0(iCell)-CGcst_alpha0*CGvec_zh1(iCell))
                  CGvec_w1(iCell)  = CGvec_y0(iCell)  - CGcst_omega0 * (CGvec_t0(iCell)-CGcst_alpha0*CGvec_v1(iCell))
               end do
   
               block => block % next
            end do  ! block
   
   
            ! Begin reduction ------------------------------------------------------------------------!
   
            CGcst_r00r1 = 0.0_RKIND
            CGcst_r00w1 = 0.0_RKIND
            CGcst_r00s0 = 0.0_RKIND
            CGcst_r00z0 = 0.0_RKIND
            CGcst_r1r1  = 0.0_RKIND
   
            do iCell = 1,nCells
               CGcst_r00r1 = CGcst_r00r1 + CGvec_r00(iCell) * CGvec_r1(iCell)               
               CGcst_r00w1 = CGcst_r00w1 + CGvec_r00(iCell) * CGvec_w1(iCell)               
               CGcst_r00s0 = CGcst_r00s0 + CGvec_r00(iCell) * CGvec_s1(iCell) ! s1
               CGcst_r00z0 = CGcst_r00z0 + CGvec_r00(iCell) * CGvec_z1(iCell) ! z1
               CGcst_r1r1  = CGcst_r1r1  + CGvec_r1(iCell)  * CGvec_r1(iCell) 
            end do
   
            CGcst_allreduce(1) = CGcst_r00r1
            CGcst_allreduce(2) = CGcst_r00w1
            CGcst_allreduce(3) = CGcst_r00s0
            CGcst_allreduce(4) = CGcst_r00z0
            CGcst_allreduce(5) = CGcst_r1r1
   
            ! Global sum across CPUs
            call mpas_timer_start("si reduction iter")
            call mpas_dmpar_sum_real_array(dminfo, 5, CGcst_allreduce, CGcst_allreduce_global)
            call mpas_timer_stop("si reduction iter")
   
            CGcst_r00r1_global = CGcst_allreduce_global(1)
            CGcst_r00w1_global = CGcst_allreduce_global(2)
            CGcst_r00s0_global = CGcst_allreduce_global(3)
            CGcst_r00z0_global = CGcst_allreduce_global(4)
                         resid = CGcst_allreduce_global(5)
   
   
            ! Preconditioning ------------------------------------------------------------------------!
   
            if ( trim(config_btr_si_preconditioner) == 'ras' ) then
               call mpas_timer_start("si halo iter")
               call mpas_dmpar_exch_group_create(domain, iterGroupName)
               call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'CGvec_w1', 1)
               call mpas_threading_barrier()
               call mpas_dmpar_exch_group_full_halo_exch(domain, iterGroupName)
               call mpas_dmpar_exch_group_destroy(domain, iterGroupName)
               call mpas_timer_stop("si halo iter")
            end if
   
   
            block => domain % blocklist
            do while (associated(block))
   
               call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
               call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
     
               call mpas_pool_get_subpool(block % structs, 'tend'       , tendPool       )
               call mpas_pool_get_subpool(block % structs, 'mesh'       , meshPool       )
               call mpas_pool_get_subpool(block % structs, 'state'      , statePool      )
               call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
   
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_w1' , CGvec_w1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_wh1', CGvec_wh1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_v0' , CGvec_v0)
   
               nCells = nCellsArray(1)
               nEdges = nEdgesArray(2)
   
               if ( trim(config_btr_si_preconditioner) == 'ras' ) then
                  ! RAS preconditioning: Use BLAS for symmetric matrix-vector multiplication
                  call DSYMV('U', nPrecVec, 1.0_RKIND, prec_ivmat, nPrecVec, CGvec_w1(1:nPrecVec), 1, 0.0_RKIND, CGvec_wh1(1:nPrecVec), 1)
   
               elseif ( trim(config_btr_si_preconditioner) == 'block_jacobi' ) then
                  ! Block-Jacobi preconditioning: Use BLAS for symmetric matrix-vector multiplication
                  call DSYMV('U', nPrecVec, 1.0_RKIND, prec_ivmat, nPrecVec, CGvec_w1(1:nPrecVec), 1, 0.0_RKIND, CGvec_wh1(1:nPrecVec), 1)
               elseif ( trim(config_btr_si_preconditioner) == 'jacobi' ) then
                  ! Jacobi preconditioning 
                  CGvec_wh1(1:nPrecVec) = CGvec_w1(1:nPrecVec) * prec_ivmat(1:nPrecVec,1)
   
               elseif ( trim(config_btr_si_preconditioner) == 'none' ) then
                  ! No preconditioning 
                  CGvec_wh1(1:nPrecVec) = CGvec_w1(1:nPrecVec)
               end if
   
               block => block % next
            end do  ! block
   
            call mpas_timer_start("si halo iter")
            call mpas_dmpar_exch_group_create(domain, iterGroupName)
            call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'CGvec_wh1', 1)
            call mpas_threading_barrier()
            call mpas_dmpar_exch_group_full_halo_exch(domain, iterGroupName)
            call mpas_dmpar_exch_group_destroy(domain, iterGroupName)
            call mpas_timer_stop("si halo iter")
   
   
            ! SpMV -----------------------------------------------------------------------------------!
   
            block => domain % blocklist
            do while (associated(block))
   
               call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
               call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
   
               call mpas_pool_get_subpool(block % structs, 'tend'       , tendPool       )
               call mpas_pool_get_subpool(block % structs, 'mesh'       , meshPool       )
               call mpas_pool_get_subpool(block % structs, 'state'      , statePool      )
               call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
   
               call mpas_pool_get_array(meshPool, 'nEdgesOnCell',            nEdgesOnCell           )
               call mpas_pool_get_array(meshPool, 'edgesOnCell',             edgesOnCell            )
               call mpas_pool_get_array(meshPool, 'cellsOnEdge',             cellsOnEdge            )
               call mpas_pool_get_array(meshPool, 'dcEdge',                  dcEdge                 )
               call mpas_pool_get_array(meshPool, 'bottomDepth',             bottomDepth            )
               call mpas_pool_get_array(meshPool, 'edgeSignOnCell',          edgeSignOnCell         )
               call mpas_pool_get_array(meshPool, 'dvEdge',                  dvEdge                 )
               call mpas_pool_get_array(meshPool, 'areaCell',                areaCell               )
   
               call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
   
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_r0' , CGvec_r0 )
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_r1' , CGvec_r1 )
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_rh0', CGvec_rh0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_rh1', CGvec_rh1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_w0' , CGvec_w0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_w1' , CGvec_w1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_wh0', CGvec_wh0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_wh1', CGvec_wh1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_t0' , CGvec_t0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_t1' , CGvec_t1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_ph0', CGvec_ph0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_ph1', CGvec_ph1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_v0' , CGvec_v0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_v1' , CGvec_v1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_s0' , CGvec_s0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_s1' , CGvec_s1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_sh0', CGvec_sh0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_sh1', CGvec_sh1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_z0' , CGvec_z0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_z1' , CGvec_z1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_zh0', CGvec_zh0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_zh1', CGvec_zh1)
   
               nCells = nCellsArray(1)
               nEdges = nEdgesArray(2)
   
               do iCell = 1, nCells
                  sshTendAx = 0.0_RKIND
   
                  do i = 1, nEdgesOnCell(iCell)
                     iEdge = edgesOnCell(i, iCell)
     
                     cell1 = cellsOnEdge(1, iEdge)
                     cell2 = cellsOnEdge(2, iEdge)
     
                     ! Interpolation sshEdge
                     sshEdgeLag = 0.5_RKIND * (sshCur(cell1) + sshCur(cell2))
     
                     ! method 1, matches method 0 without pbcs, works with pbcs.
                     thicknessSumLag = sshEdgeLag + min(bottomDepth(cell1), bottomDepth(cell2))
     
                     ! nabla (ssh^0)
                     sshDiffNew = (CGvec_wh1(cell2)-CGvec_wh1(cell1)) / dcEdge(iEdge)
     
                     fluxAx = thicknessSumLag * sshDiffNew
    
                     sshTendAx = sshTendAx + edgeSignOnCell(i, iCell) * fluxAx * dvEdge(iEdge)
                  end do ! i
                   
                  sshLagArea = (1.0_RKIND/(gravity*dt**2.0*alpha1**2.0)) * CGvec_wh1(iCell) * areaCell(iCell)
   
                  CGvec_t1(iCell) = -sshLagArea - sshTendAx
               end do ! iCell
          
               ! End reduction -----------------------------------------------------------------------!
   
               CGcst_beta0 = (CGcst_alpha0/CGcst_omega0) * CGcst_r00r1_global / CGcst_r00r0_global
               CGcst_alpha1 =  CGcst_r00r1_global / (  CGcst_r00w1_global + CGcst_beta0 * CGcst_r00s0_global  &
                                                     - CGcst_beta0 * CGcst_omega0 * CGcst_r00z0_global       )
   
               do iCell = 1,nCells
                  CGvec_r0(iCell) = CGvec_r1(iCell)
                  CGvec_s0(iCell) = CGvec_s1(iCell)
                  CGvec_z0(iCell) = CGvec_z1(iCell)
                  CGvec_w0(iCell) = CGvec_w1(iCell)
                  CGvec_t0(iCell) = CGvec_t1(iCell)
                  CGvec_v0(iCell) = CGvec_v1(iCell)
    
                  CGvec_rh0(iCell) = CGvec_rh1(iCell)
                  CGvec_sh0(iCell) = CGvec_sh1(iCell)
                  CGvec_ph0(iCell) = CGvec_ph1(iCell)
                  CGvec_wh0(iCell) = CGvec_wh1(iCell)
                  CGvec_zh0(iCell) = CGvec_zh1(iCell)
    
                  sshSubcycleCur(iCell) = sshSubcycleNew(iCell)
                          sshNew(iCell) = sshSubcycleNew(iCell)
               end do ! iCell
       	
               CGcst_r00r0_global = CGcst_r00r1_global
               CGcst_alpha0       = CGcst_alpha1
   
               block => block % next
            end do  ! block

         !**************************************************************!
         end do ! do iter = 2
         !**************************************************************!

         !   boundary update on SSHnew
         call mpas_timer_start("si halo iter")
         call mpas_dmpar_exch_group_create(domain, iterGroupName)
         call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'sshSubcycle', timeLevel=1)
         call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'ssh', timeLevel=2)
         call mpas_threading_barrier()
         call mpas_dmpar_exch_group_full_halo_exch(domain, iterGroupName)
         call mpas_dmpar_exch_group_destroy(domain, iterGroupName)
         call mpas_timer_stop("si halo iter")

         call mpas_timer_stop("si btr iteration")


         call mpas_timer_start("si btr vel update")
         ! Update normalBarotropicVelocitySubcycle ------------------------------------------------!
         
         block => domain % blocklist
         do while (associated(block))

            call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
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

            call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleCur, 1)
            call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleNew, 2)
            call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityCur,1)
            call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
            call mpas_pool_get_array(statePool, 'ssh', sshNew, 2)

            call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)
            call mpas_pool_get_array(diagnosticsPool, 'barotropicCoriolisTerm',barotropicCoriolisTerm)

            nCells = nCellsArray( cellHaloComputeCounter )
            nEdges = nEdgesArray( edgeHaloComputeCounter )

            !$omp do schedule(runtime) private(cell1, cell2, CoriolisTerm, i, eoe)
            do iEdge = 1, nEdges

               temp_mask = edgeMask(1, iEdge)

               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)

               ! Compute the barotropic Coriolis term, -f*uPerp
               CoriolisTerm = 0.0_RKIND

               do i = 1, nEdgesOnEdge(iEdge)
                  eoe = edgesOnEdge(i,iEdge)
                  wgt = weightsOnEdge(i,iEdge)
                  CoriolisTerm = CoriolisTerm + wgt  &
                               * normalBarotropicVelocityCur(eoe) * fEdge(eoe)
               end do
                    barotropicCoriolisTerm(iEdge) = CoriolisTerm

               normalBarotropicVelocitySubcycleNew(iEdge) &
                  = temp_mask &
                  * (normalBarotropicVelocityCur(iEdge) &
                  - dt * (-barotropicCoriolisTerm(iEdge) + gravity &
                  * ( (alpha1*sshNew(cell2)+alpha2*sshCur(cell2)) - (alpha1*sshNew(cell1)+alpha2*sshCur(cell1))) &
                  / dcEdge(iEdge) - barotropicForcing(iEdge)))
            end do ! iEdge

            do iEdge = 1, nEdges
               normalBarotropicVelocitySubcycleCur(iEdge) =  alpha1* normalBarotropicVelocitySubcycleNew(iEdge)  &
                                                           + alpha2* normalBarotropicVelocityCur(iEdge)        
            end do ! iEdge

            block => block % next
         end do  ! block

         call mpas_timer_start("si halo btr vel")
         call mpas_dmpar_exch_group_create(domain, iterGroupName)
         call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'normalBarotropicVelocitySubcycle',1)
         call mpas_threading_barrier()
         call mpas_dmpar_exch_group_full_halo_exch(domain, iterGroupName)
         call mpas_dmpar_exch_group_destroy(domain, iterGroupName)
         call mpas_timer_stop("si halo btr vel")
                
         
         ! Update barotropicCoriolisTerm ----------------------------------------------------------!
 
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
            call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
            call mpas_pool_get_subpool(block % structs, 'tend'       , tendPool       )

            call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
            call mpas_pool_get_array(meshPool, 'nEdgesOnEdge', nEdgesOnEdge)
            call mpas_pool_get_array(meshPool, 'edgesOnEdge', edgesOnEdge)
            call mpas_pool_get_array(meshPool, 'weightsOnEdge', weightsOnEdge)
            call mpas_pool_get_array(meshPool, 'fEdge', fEdge)

            call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleCur, 1)
            call mpas_pool_get_array(diagnosticsPool, 'barotropicCoriolisTerm', barotropicCoriolisTerm)

            nCells = nCellsArray(1)
            nEdges = nCellsArray(2)

            !$omp do schedule(runtime) private(cell1, cell2, CoriolisTerm, i, eoe)
            do iEdge = 1, nEdges

               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)

               ! Compute the barotropic Coriolis term, -f*uPerp
               temp1 = 0.0_RKIND
               do i = 1, nEdgesOnEdge(iEdge)
                  eoe = edgesOnEdge(i,iEdge)
                  temp1 = temp1 + weightsOnEdge(i,iEdge)  &
                                * normalBarotropicVelocitySubcycleCur(eoe) * fEdge(eoe)
               end do

               barotropicCoriolisTerm(iEdge) = temp1
            end do
            !$omp end do

            block => block % next
         end do  ! block

         call mpas_timer_start("si halo btr vel")
         call mpas_dmpar_exch_group_create(domain, iterGroupName)
         call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'barotropicCoriolisTerm', 1 )
         call mpas_threading_barrier()
         call mpas_dmpar_exch_group_full_halo_exch(domain, iterGroupName)
         call mpas_dmpar_exch_group_destroy(domain, iterGroupName)
         call mpas_timer_stop("si halo btr vel")

         call mpas_timer_stop("si btr vel update")


         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         ! Stage 2.4 : Recompute initial residual with lagged values
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         call mpas_timer_start("si btr residual")

         ! SpMV -----------------------------------------------------------------------------------!

         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
            call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
 
            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
 
            call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
            call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
            call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
            call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
            call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
            call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
            call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
            call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
 
            call mpas_pool_get_array(statePool, 'ssh', sshCur,1)
            call mpas_pool_get_array(statePool, 'ssh', sshNew,2)
            call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityCur,1)
 
            call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)
            call mpas_pool_get_array(diagnosticsPool, 'barotropicCoriolisTerm',barotropicCoriolisTerm)
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_r0', CGvec_r0)
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_r00', CGvec_r00)
 
            nCells = nCellsArray( 1 )
            nEdges = nEdgesArray( 2 )
 
            do iCell = 1, nCells
               sshTendb1 = 0.0_RKIND
               sshTendb2 = 0.0_RKIND
               sshTendAx = 0.0_RKIND

               do i = 1, nEdgesOnCell(iCell)
                  iEdge = edgesOnCell(i, iCell)
   
                  cell1 = cellsOnEdge(1, iEdge)
                  cell2 = cellsOnEdge(2, iEdge)
   
                  ! Interpolation sshEdge
                  sshEdgeCur = 0.5_RKIND * (sshCur(cell1) + sshCur(cell2))
                  sshEdgeLag = 0.5_RKIND * (sshNew(cell1) + sshNew(cell2))
                  sshEdgeMid = 0.5_RKIND * (sshEdgeLag + sshEdgeCur)
   
                  ! method 1, matches method 0 without pbcs, works with pbcs.
                  thicknessSumCur = sshEdgeCur + min(bottomDepth(cell1), bottomDepth(cell2))
                  thicknessSumLag = sshEdgeLag + min(bottomDepth(cell1), bottomDepth(cell2))
                  thicknessSumMid = sshEdgeMid + min(bottomDepth(cell1), bottomDepth(cell2))
   
                  ! nabla (ssh^0)
                  sshDiffCur = (sshCur(cell2)-sshCur(cell1)) / dcEdge(iEdge)
                  sshDiffLag = (sshNew(cell2)-sshNew(cell1)) / dcEdge(iEdge)
   
                  fluxb1 = thicknessSumMid * normalBarotropicVelocityCur(iEdge)
                  fluxb2 = thicknessSumLag * (alpha2*gravity*sshDiffCur + (-barotropicCoriolisTerm(iEdge)-barotropicForcing(iEdge)) )
                  fluxAx = thicknessSumLag * sshDiffLag
     
                  sshTendb1 = sshTendb1 + edgeSignOnCell(i, iCell) * fluxb1 * dvEdge(iEdge)
                  sshTendb2 = sshTendb2 + edgeSignOnCell(i, iCell) * fluxb2 * dvEdge(iEdge)
                  sshTendAx = sshTendAx + edgeSignOnCell(i, iCell) * fluxAx * dvEdge(iEdge)
               end do ! i
 
               sshTendb1 = (1.0_RKIND/(gravity*dt*alpha1**2.0)) * sshTendb1
               sshTendb2 = (1.0_RKIND/(gravity   *alpha1     )) * sshTendb2

               sshCurArea = (1.0_RKIND/(gravity*dt**2.0*alpha1**2.0)) * sshCur(iCell) * areaCell(iCell)
               sshLagArea = (1.0_RKIND/(gravity*dt**2.0*alpha1**2.0)) * sshNew(iCell) * areaCell(iCell)

               CGvec_r0(iCell) = (-sshCurArea - sshTendb1 + sshTendb2)   &
                                -(-sshLagArea - sshTendAx) 
               CGvec_r00(iCell) = CGvec_r0(iCell)
            end do ! iCell
 
            block => block % next
         end do  ! block


         ! Preconditioning ------------------------------------------------------------------------!

         if ( trim(config_btr_si_preconditioner) == 'ras' ) then
            call mpas_timer_start("si halo r0")
            call mpas_dmpar_exch_group_create(domain, iterGroupName)
            call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'CGvec_r0', 1)
            call mpas_threading_barrier()
            call mpas_dmpar_exch_group_full_halo_exch(domain, iterGroupName)
            call mpas_dmpar_exch_group_destroy(domain, iterGroupName)
            call mpas_timer_stop("si halo r0")
         end if

         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
            call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

            call mpas_pool_get_array(diagnosticsPool, 'CGvec_r0', CGvec_r0)
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_r00', CGvec_r00)
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_rh0', CGvec_rh0)

            nCells = nCellsArray( 1 )
            nEdges = nEdgesArray( 2 )

            if ( trim(config_btr_si_preconditioner) == 'ras' ) then
               ! RAS preconditioning: Use BLAS for symmetric matrix-vector multiplication
               call DSYMV('U', nPrecVec, 1.0_RKIND, prec_ivmat, nPrecVec, CGvec_r0(1:nPrecVec), 1, 0.0_RKIND, CGvec_rh0(1:nPrecVec), 1)

            elseif ( trim(config_btr_si_preconditioner) == 'block_jacobi' ) then
               ! Block-Jacobi preconditioning: Use BLAS for symmetric matrix-vector multiplication
               call DSYMV('U', nPrecVec, 1.0_RKIND, prec_ivmat, nPrecVec, CGvec_r0(1:nPrecVec), 1, 0.0_RKIND, CGvec_rh0(1:nPrecVec), 1)

            elseif ( trim(config_btr_si_preconditioner) == 'jacobi' ) then
               ! Jacobi preconditioning
               CGvec_rh0(1:nPrecVec) = CGvec_r0(1:nPrecVec) * prec_ivmat(1:nPrecVec,1)

            elseif ( trim(config_btr_si_preconditioner) == 'none' ) then
               ! No preconditioning
               CGvec_rh0(1:nPrecVec) = CGvec_r0(1:nPrecVec)

            end if

            block => block % next
         end do  ! block

         call mpas_timer_start("si halo r0")
         call mpas_dmpar_exch_group_create(domain, iterGroupName)
         call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'CGvec_rh0', 1)
         call mpas_threading_barrier()
         call mpas_dmpar_exch_group_full_halo_exch(domain, iterGroupName)
         call mpas_dmpar_exch_group_destroy(domain, iterGroupName)
         call mpas_timer_stop("si halo r0")


         ! SpMV -----------------------------------------------------------------------------------!

         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
            call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
 
            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
 
            call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
            call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
            call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
            call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
            call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
            call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
            call mpas_pool_get_array(meshPool, 'refBottomDepthTopOfCell', refBottomDepthTopOfCell)
            call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
            call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
            call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
 
            call mpas_pool_get_array(statePool, 'ssh', sshNew,2)
 
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_r0' , CGvec_r0 )
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_r00', CGvec_r00)
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_rh0', CGvec_rh0)
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_w0' , CGvec_w0)
 
            nCells = nCellsArray( 1 )
            nEdges = nEdgesArray( 2 )
 
            CGcst_r00r0 = 0.0_RKIND
            CGcst_r00w0 = 0.0_RKIND
 
            do iCell = 1, nCells
 
               sshTendAx = 0.0_RKIND
 
               do i = 1, nEdgesOnCell(iCell)
                  iEdge = edgesOnCell(i, iCell)
 
                  cell1 = cellsOnEdge(1, iEdge)
                  cell2 = cellsOnEdge(2, iEdge)
 
                  ! Interpolation sshEdge
                  sshEdgeLag = 0.5_RKIND * (sshNew(cell1) + sshNew(cell2))
 
                  ! method 1, matches method 0 without pbcs, works with pbcs.
                  thicknessSumLag = sshEdgeLag + min(bottomDepth(cell1), bottomDepth(cell2))
 
                  ! nabla (ssh^0)
                  sshDiffLag = (CGvec_rh0(cell2)- CGvec_rh0(cell1)) / dcEdge(iEdge)
 
                  fluxAx = thicknessSumLag * sshDiffLag
 
                  sshTendAx = sshTendAx + edgeSignOnCell(i, iCell) * fluxAx * dvEdge(iEdge)
 
               end do ! i
 
               sshCurArea = (1.0_RKIND/(gravity*dt**2.0*alpha1**2.0)) * CGvec_rh0(iCell) * areaCell(iCell)
 
               CGvec_w0(iCell) = -sshCurArea - sshTendAx
 
               CGcst_r00r0 = CGcst_r00r0 + CGvec_r00(iCell) * CGvec_r0(iCell)
               CGcst_r00w0 = CGcst_r00w0 + CGvec_r00(iCell) * CGvec_w0(iCell)
 
            end do ! iCell
 
            block => block % next
         end do  ! block


         ! Preconditioning ------------------------------------------------------------------------!

         if ( trim(config_btr_si_preconditioner) == 'ras' ) then
            call mpas_timer_start("si halo r0")
            call mpas_dmpar_exch_group_create(domain, iterGroupName)
            call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'CGvec_w0', 1)
            call mpas_threading_barrier()
            call mpas_dmpar_exch_group_full_halo_exch(domain, iterGroupName)
            call mpas_dmpar_exch_group_destroy(domain, iterGroupName)
            call mpas_timer_stop("si halo r0")
         end if

         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
            call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
 
            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
 
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_w0' , CGvec_w0)
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_wh0', CGvec_wh0)
 
            nCells = nCellsArray( 1 )
            nEdges = nEdgesArray( 2 )
 
            if ( trim(config_btr_si_preconditioner) == 'ras' ) then
               ! RAS preconditioning: Use BLAS for symmetric matrix-vector multiplication
               call DSYMV('U', nPrecVec, 1.0_RKIND, prec_ivmat, nPrecVec, CGvec_w0(1:nPrecVec), 1, 0.0_RKIND, CGvec_wh0(1:nPrecVec), 1)
 
            elseif ( trim(config_btr_si_preconditioner) == 'block_jacobi' ) then
               ! Block-Jacobi preconditioning: Use BLAS for symmetric matrix-vector multiplication
               call DSYMV('U', nPrecVec, 1.0_RKIND, prec_ivmat, nPrecVec, CGvec_w0(1:nPrecVec), 1, 0.0_RKIND, CGvec_wh0(1:nPrecVec), 1)
 
            elseif ( trim(config_btr_si_preconditioner) == 'jacobi' ) then
               ! Jacobi preconditioning
               CGvec_wh0(1:nPrecVec) = CGvec_w0(1:nPrecVec) * prec_ivmat(1:nPrecVec,1)
 
            elseif ( trim(config_btr_si_preconditioner) == 'none' ) then
               ! No preconditioning
               CGvec_wh0(1:nPrecVec) = CGvec_w0(1:nPrecVec)
            end if

           block => block % next
         end do  ! block

         call mpas_timer_start("si halo r0")
         call mpas_dmpar_exch_group_create(domain, iterGroupName)
         call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'CGvec_wh0', 1)
         call mpas_threading_barrier()
         call mpas_dmpar_exch_group_full_halo_exch(domain, iterGroupName)
         call mpas_dmpar_exch_group_destroy(domain, iterGroupName)
         call mpas_timer_stop("si halo r0")


         ! SpMV -----------------------------------------------------------------------------------!

         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
            call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
 
            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
 
            call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
            call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
            call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
            call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
            call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
            call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
            call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
            call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
 
            call mpas_pool_get_array(statePool, 'ssh', sshNew,2)
 
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_wh0', CGvec_wh0)
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_t0' , CGvec_t0)
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_ph0', CGvec_ph0)
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_v0' , CGvec_v0)
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_s0' , CGvec_s0)
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_sh0', CGvec_sh0)
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_z0' , CGvec_z0)
            call mpas_pool_get_array(diagnosticsPool, 'CGvec_zh0', CGvec_zh0)
 
            nCells = nCellsArray( 1 )
            nEdges = nEdgesArray( 2 )

            do iCell = 1, nCells
 
               sshTendAx = 0.0_RKIND
  
               do i = 1, nEdgesOnCell(iCell)
                  iEdge = edgesOnCell(i, iCell)
   
                  cell1 = cellsOnEdge(1, iEdge)
                  cell2 = cellsOnEdge(2, iEdge)
   
                  ! Interpolation sshEdge
                  sshEdgeLag = 0.5_RKIND * (sshNew(cell1) + sshNew(cell2))
   
                  ! method 1, matches method 0 without pbcs, works with pbcs.
                  thicknessSumLag = sshEdgeLag + min(bottomDepth(cell1), bottomDepth(cell2))
   
                  ! nabla (ssh^0)
                  sshDiffCur = (CGvec_wh0(cell2)- CGvec_wh0(cell1)) / dcEdge(iEdge)
   
                  fluxAx = thicknessSumLag * sshDiffCur
   
                  sshTendAx = sshTendAx + edgeSignOnCell(i, iCell) * fluxAx * dvEdge(iEdge)
               end do ! i
 
               sshCurArea = (1.0_RKIND/(gravity*dt**2.0*alpha1**2.0)) * CGvec_wh0(iCell) * areaCell(iCell)

               CGvec_t0(iCell) = -sshCurArea - sshTendAx

               CGvec_ph0(iCell) = 0.0_RKIND
               CGvec_sh0(iCell) = 0.0_RKIND
               CGvec_z0(iCell)  = 0.0_RKIND
               CGvec_zh0(iCell) = 0.0_RKIND
               CGvec_v0(iCell)  = 0.0_RKIND
               CGvec_s0(iCell)  = 0.0_RKIND
 
            end do ! iCell
 
            block => block % next
         end do  ! block

         CGcst_allreduce(:) = 0.0_RKIND
         CGcst_allreduce(1) = CGcst_r00r0
         CGcst_allreduce(2) = CGcst_r00w0

         ! Global sum across CPU 
         call mpas_timer_start("si reduction r0")
         call mpas_dmpar_sum_real_array(dminfo, 5, CGcst_allreduce, CGcst_allreduce_global)
         call mpas_timer_stop("si reduction r0")

         CGcst_r00r0_global = CGcst_allreduce_global(1)
         CGcst_r00w0_global = CGcst_allreduce_global(2)

         CGcst_alpha0 = CGcst_r00r0_global / CGcst_r00w0_global
         CGcst_beta0  = 0.0_RKIND
         CGcst_omega0 = 0.0_RKIND

         call mpas_timer_stop("si btr residual")


         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         ! Stage 2.5 : Main iterations
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         call mpas_timer_start("si btr iteration")

         iter = 0
         resid = (crit_main+100.0)**2.0

         !**************************************************************!
         do while ( dsqrt(resid) > crit_main ) 
         !**************************************************************!

            iter = iter + 1
   
            block => domain % blocklist
            do while (associated(block))
               call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
               call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
     
               call mpas_pool_get_subpool(block % structs, 'tend'       , tendPool       )
               call mpas_pool_get_subpool(block % structs, 'mesh'       , meshPool       )
               call mpas_pool_get_subpool(block % structs, 'state'      , statePool      )
               call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

               call mpas_pool_get_subpool(statePool, 'tracers'    , tracersPool    )
               call mpas_pool_get_subpool(tendPool , 'tracersTend', tracersTendPool)
     
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_r0' , CGvec_r0 )
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_r00', CGvec_r00)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_rh0', CGvec_rh0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_w0' , CGvec_w0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_wh0', CGvec_wh0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_t0' , CGvec_t0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_ph0', CGvec_ph0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_ph1', CGvec_ph1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_v0' , CGvec_v0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_s0' , CGvec_s0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_s1' , CGvec_s1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_sh0', CGvec_sh0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_sh1', CGvec_sh1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_q0' , CGvec_q0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_qh0', CGvec_qh0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_z0' , CGvec_z0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_z1' , CGvec_z1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_zh0', CGvec_zh0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_y0' , CGvec_y0)
     
               nCells = nCellsArray(1)
               nEdges = nEdgesArray(2)
     
               do iCell = 1, nCells
                  CGvec_ph1(iCell) = CGvec_rh0(iCell) + CGcst_beta0  * (CGvec_ph0(iCell)-CGcst_omega0*CGvec_sh0(iCell))
                  CGvec_s1(iCell)  = CGvec_w0(iCell)  + CGcst_beta0  * (CGvec_s0(iCell)-CGcst_omega0*CGvec_z0(iCell))
                  CGvec_sh1(iCell) = CGvec_wh0(iCell) + CGcst_beta0  * (CGvec_sh0(iCell)-CGcst_omega0*CGvec_zh0(iCell))
                  CGvec_z1(iCell)  = CGvec_t0(iCell)  + CGcst_beta0  * (CGvec_z0(iCell)-CGcst_omega0*CGvec_v0(iCell))
                  CGvec_q0(iCell)  = CGvec_r0(iCell)  - CGcst_alpha0 * CGvec_s1(iCell)
                  CGvec_qh0(iCell) = CGvec_rh0(iCell) - CGcst_alpha0 * CGvec_sh1(iCell)
                  CGvec_y0(iCell)  = CGvec_w0(iCell)  - CGcst_alpha0 * CGvec_z1(iCell)
               end do ! iCell
   
   
               ! Begin reduction ----------------------------------------------------------------------!
   
               CGcst_q0y0  = 0.0_RKIND
               CGcst_y0y0  = 0.0_RKIND
               CGcst_r00r0 = 0.0_RKIND
   
               do iCell = 1,nCells
                  CGcst_q0y0  = CGcst_q0y0  + CGvec_q0(iCell)  * CGvec_y0(iCell)               
                  CGcst_y0y0  = CGcst_y0y0  + CGvec_y0(iCell)  * CGvec_y0(iCell)               
                  CGcst_r00r0 = CGcst_r00r0 + CGvec_r00(iCell) * CGvec_r0(iCell)
               end do
   
               block => block % next
            end do  ! block
   
            CGcst_allreduce(:) = 0.0_RKIND
            CGcst_allreduce(1) = CGcst_q0y0
            CGcst_allreduce(2) = CGcst_y0y0
            CGcst_allreduce(3) = CGcst_r00r0
   
            ! Global sum across CPUs
            call mpas_timer_start("si reduction iter")
            call mpas_dmpar_sum_real_array(dminfo, 5, CGcst_allreduce, CGcst_allreduce_global)
            call mpas_timer_stop("si reduction iter")
   
            CGcst_q0y0_global  = CGcst_allreduce_global(1)
            CGcst_y0y0_global  = CGcst_allreduce_global(2)
            CGcst_r00r0_global = CGcst_allreduce_global(3)
   
   
            ! Preconditioning ------------------------------------------------------------------------!
   
            if ( trim(config_btr_si_preconditioner) == 'ras' ) then
               call mpas_timer_start("si halo iter")
               call mpas_dmpar_exch_group_create(domain, iterGroupName)
               call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'CGvec_z1', 1)
               call mpas_threading_barrier()
               call mpas_dmpar_exch_group_full_halo_exch(domain, iterGroupName)
               call mpas_dmpar_exch_group_destroy(domain, iterGroupName)
               call mpas_timer_stop("si halo iter")
            end if
   
            block => domain % blocklist
            do while (associated(block))
   
               call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
               call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
    
               call mpas_pool_get_subpool(block % structs, 'tend'       , tendPool       )
               call mpas_pool_get_subpool(block % structs, 'mesh'       , meshPool       )
               call mpas_pool_get_subpool(block % structs, 'state'      , statePool      )
               call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
    
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_z1' , CGvec_z1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_zh1', CGvec_zh1)
    
               nCells = nCellsArray(1)
               nEdges = nEdgesArray(2)
   
               if ( trim(config_btr_si_preconditioner) == 'ras' ) then
                  ! RAS preconditioning: Use BLAS for symmetric matrix-vector multiplication
                  call DSYMV('U', nPrecVec, 1.0_RKIND, prec_ivmat, nPrecVec, CGvec_z1(1:nPrecVec), 1, 0.0_RKIND, CGvec_zh1(1:nPrecVec), 1)
   
               elseif ( trim(config_btr_si_preconditioner) == 'block_jacobi' ) then
                  ! Block-Jacobi preconditioning: Use BLAS for symmetric matrix-vector multiplication
                  call DSYMV('U', nPrecVec, 1.0_RKIND, prec_ivmat, nPrecVec, CGvec_z1(1:nPrecVec), 1, 0.0_RKIND, CGvec_zh1(1:nPrecVec), 1)
   
               elseif ( trim(config_btr_si_preconditioner) == 'jacobi' ) then
                  ! Jacobi preconditioning 
                  CGvec_zh1(1:nPrecVec) = CGvec_z1(1:nPrecVec) * prec_ivmat(1:nPrecVec,1)
   
               elseif ( trim(config_btr_si_preconditioner) == 'none' ) then
                  ! No preconditioning 
                  CGvec_zh1(1:nPrecVec) = CGvec_z1(1:nPrecVec)
               end if
   
               block => block % next
            end do  ! block
   
            call mpas_timer_start("si halo iter")
            call mpas_dmpar_exch_group_create(domain, iterGroupName)
            call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'CGvec_zh1', 1)
            call mpas_threading_barrier()
            call mpas_dmpar_exch_group_full_halo_exch(domain, iterGroupName)
            call mpas_dmpar_exch_group_destroy(domain, iterGroupName)
            call mpas_timer_stop("si halo iter")
   
   
            ! SpMV -----------------------------------------------------------------------------------!
   
            block => domain % blocklist
            do while (associated(block))
   
               call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
               call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
    
               call mpas_pool_get_subpool(block % structs, 'tend'       , tendPool       )
               call mpas_pool_get_subpool(block % structs, 'mesh'       , meshPool       )
               call mpas_pool_get_subpool(block % structs, 'state'      , statePool      )
               call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
    
               call mpas_pool_get_array(meshPool, 'nEdgesOnCell',            nEdgesOnCell           )
               call mpas_pool_get_array(meshPool, 'edgesOnCell',             edgesOnCell            )
               call mpas_pool_get_array(meshPool, 'cellsOnEdge',             cellsOnEdge            )
               call mpas_pool_get_array(meshPool, 'dcEdge',                  dcEdge                 )
               call mpas_pool_get_array(meshPool, 'bottomDepth',             bottomDepth            )
               call mpas_pool_get_array(meshPool, 'edgeSignOnCell',          edgeSignOnCell         )
               call mpas_pool_get_array(meshPool, 'dvEdge',                  dvEdge                 )
               call mpas_pool_get_array(meshPool, 'areaCell',                areaCell               )
    
               call mpas_pool_get_array(statePool, 'ssh', sshNew, 2)
               call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, 1)
               call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleNew, 2)
    
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_r1' , CGvec_r1 )
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_rh1', CGvec_rh1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_w1' , CGvec_w1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_wh0', CGvec_wh0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_t0' , CGvec_t0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_ph1', CGvec_ph1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_v1' , CGvec_v1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_q0' , CGvec_q0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_qh0', CGvec_qh0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_zh1', CGvec_zh1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_y0' , CGvec_y0)
    
               nCells = nCellsArray(1)
               nEdges = nEdgesArray(2)
    
               do iCell = 1, nCells
                  sshTendAx = 0.0_RKIND
   
                  do i = 1, nEdgesOnCell(iCell)
                     iEdge = edgesOnCell(i, iCell)
      
                     cell1 = cellsOnEdge(1, iEdge)
                     cell2 = cellsOnEdge(2, iEdge)
      
                     ! Interpolation sshEdge
                     sshEdgeLag = 0.5_RKIND * (sshNew(cell1) + sshNew(cell2))
      
                     ! method 1, matches method 0 without pbcs, works with pbcs.
                     thicknessSumLag = sshEdgeLag + min(bottomDepth(cell1), bottomDepth(cell2))
      
                     ! nabla (ssh^0)
                     sshDiffNew = (CGvec_zh1(cell2)-CGvec_zh1(cell1)) / dcEdge(iEdge)
      
                     fluxAx = thicknessSumLag * sshDiffNew
     
                     sshTendAx = sshTendAx + edgeSignOnCell(i, iCell) * fluxAx * dvEdge(iEdge)
                  end do ! i
                   
                  sshLagArea = (1.0_RKIND/(gravity*dt**2.0*alpha1**2.0)) * CGvec_zh1(iCell) * areaCell(iCell)
    
                  CGvec_v1(iCell) = -sshLagArea - sshTendAx
            
               end do ! iCell
           
               ! End reduction ------------------------------------------------------------------------!
   
               CGcst_omega0 = CGcst_q0y0_global / CGcst_y0y0_global
   
               do iCell = 1,nCells
                  sshSubcycleNew(iCell) = sshSubcycleCur(iCell) + CGcst_alpha0 * CGvec_ph1(iCell) &
                                                                + CGcst_omega0 * CGvec_qh0(iCell)
                  CGvec_r1(iCell)  = CGvec_q0(iCell)  - CGcst_omega0 * CGvec_y0(iCell)
                  CGvec_rh1(iCell) = CGvec_qh0(iCell) - CGcst_omega0 * (CGvec_wh0(iCell)-CGcst_alpha0*CGvec_zh1(iCell))
                  CGvec_w1(iCell)  = CGvec_y0(iCell)  - CGcst_omega0 * (CGvec_t0(iCell)-CGcst_alpha0*CGvec_v1(iCell))
               end do
   
               block => block % next
            end do  ! block
   
   
            ! Begin reduction ------------------------------------------------------------------------!
   
            CGcst_r00r1 = 0.0_RKIND
            CGcst_r00w1 = 0.0_RKIND
            CGcst_r00s0 = 0.0_RKIND
            CGcst_r00z0 = 0.0_RKIND
            CGcst_r1r1  = 0.0_RKIND
   
            do iCell = 1,nCells
               CGcst_r00r1 = CGcst_r00r1 + CGvec_r00(iCell) * CGvec_r1(iCell)               
               CGcst_r00w1 = CGcst_r00w1 + CGvec_r00(iCell) * CGvec_w1(iCell)               
               CGcst_r00s0 = CGcst_r00s0 + CGvec_r00(iCell) * CGvec_s1(iCell) ! s1
               CGcst_r00z0 = CGcst_r00z0 + CGvec_r00(iCell) * CGvec_z1(iCell) ! z1
               CGcst_r1r1  = CGcst_r1r1  + CGvec_r1(iCell)  * CGvec_r1(iCell) 
            end do
   
            CGcst_allreduce(1) = CGcst_r00r1
            CGcst_allreduce(2) = CGcst_r00w1
            CGcst_allreduce(3) = CGcst_r00s0
            CGcst_allreduce(4) = CGcst_r00z0
            CGcst_allreduce(5) = CGcst_r1r1
   
            ! Global sum across CPUs
            call mpas_timer_start("si reduction iter")
            call mpas_dmpar_sum_real_array(dminfo, 5, CGcst_allreduce, CGcst_allreduce_global)
            call mpas_timer_stop("si reduction iter")
   
            CGcst_r00r1_global = CGcst_allreduce_global(1)
            CGcst_r00w1_global = CGcst_allreduce_global(2)
            CGcst_r00s0_global = CGcst_allreduce_global(3)
            CGcst_r00z0_global = CGcst_allreduce_global(4)
                         resid = CGcst_allreduce_global(5)
   
            ! Preconditioning ------------------------------------------------------------------------!
   
            if ( trim(config_btr_si_preconditioner) == 'ras' ) then
               call mpas_timer_start("si halo iter")
               call mpas_dmpar_exch_group_create(domain, iterGroupName)
               call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'CGvec_w1', 1)
               call mpas_threading_barrier()
               call mpas_dmpar_exch_group_full_halo_exch(domain, iterGroupName)
               call mpas_dmpar_exch_group_destroy(domain, iterGroupName)
               call mpas_timer_stop("si halo iter")
            end if
   
   
            block => domain % blocklist
            do while (associated(block))
   
               call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
               call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
    
               call mpas_pool_get_subpool(block % structs, 'tend'       , tendPool       )
               call mpas_pool_get_subpool(block % structs, 'mesh'       , meshPool       )
               call mpas_pool_get_subpool(block % structs, 'state'      , statePool      )
               call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
    
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_w1' , CGvec_w1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_wh1', CGvec_wh1)
    
               nCells = nCellsArray(1)
               nEdges = nEdgesArray(2)
   
               if ( trim(config_btr_si_preconditioner) == 'ras' ) then
                  ! RAS preconditioning: Use BLAS for symmetric matrix-vector multiplication
                  call DSYMV('U', nPrecVec, 1.0_RKIND, prec_ivmat, nPrecVec, CGvec_w1(1:nPrecVec), 1, 0.0_RKIND, CGvec_wh1(1:nPrecVec), 1)
   
               elseif ( trim(config_btr_si_preconditioner) == 'block_jacobi' ) then
                  ! Block-Jacobi preconditioning: Use BLAS for symmetric matrix-vector multiplication
                  call DSYMV('U', nPrecVec, 1.0_RKIND, prec_ivmat, nPrecVec, CGvec_w1(1:nPrecVec), 1, 0.0_RKIND, CGvec_wh1(1:nPrecVec), 1)
               elseif ( trim(config_btr_si_preconditioner) == 'jacobi' ) then
                  ! Jacobi preconditioning 
                  CGvec_wh1(1:nPrecVec) = CGvec_w1(1:nPrecVec) * prec_ivmat(1:nPrecVec,1)
   
               elseif ( trim(config_btr_si_preconditioner) == 'none' ) then
                  ! No preconditioning 
                  CGvec_wh1(1:nPrecVec) = CGvec_w1(1:nPrecVec)
               end if
   
               block => block % next
            end do  ! block
   
            call mpas_timer_start("si halo iter")
            call mpas_dmpar_exch_group_create(domain, iterGroupName)
            call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'CGvec_wh1', 1)
            call mpas_threading_barrier()
            call mpas_dmpar_exch_group_full_halo_exch(domain, iterGroupName)
            call mpas_dmpar_exch_group_destroy(domain, iterGroupName)
            call mpas_timer_stop("si halo iter")
   
   
            ! SpMV -----------------------------------------------------------------------------------!
   
            block => domain % blocklist
            do while (associated(block))
   
               call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
               call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
   
               call mpas_pool_get_subpool(block % structs, 'tend'       , tendPool       )
               call mpas_pool_get_subpool(block % structs, 'mesh'       , meshPool       )
               call mpas_pool_get_subpool(block % structs, 'state'      , statePool      )
               call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
   
               call mpas_pool_get_array(meshPool, 'nEdgesOnCell',            nEdgesOnCell           )
               call mpas_pool_get_array(meshPool, 'edgesOnCell',             edgesOnCell            )
               call mpas_pool_get_array(meshPool, 'cellsOnEdge',             cellsOnEdge            )
               call mpas_pool_get_array(meshPool, 'dcEdge',                  dcEdge                 )
               call mpas_pool_get_array(meshPool, 'bottomDepth',             bottomDepth            )
               call mpas_pool_get_array(meshPool, 'edgeSignOnCell',          edgeSignOnCell         )
               call mpas_pool_get_array(meshPool, 'dvEdge',                  dvEdge                 )
               call mpas_pool_get_array(meshPool, 'areaCell',                areaCell               )
   
               call mpas_pool_get_array(statePool, 'ssh', sshNew, 2)
               call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, 1)
               call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleNew, 2)
   
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_r0' , CGvec_r0 )
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_r1' , CGvec_r1 )
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_rh0', CGvec_rh0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_rh1', CGvec_rh1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_w0' , CGvec_w0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_w1' , CGvec_w1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_wh0', CGvec_wh0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_wh1', CGvec_wh1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_t0' , CGvec_t0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_t1' , CGvec_t1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_ph0', CGvec_ph0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_ph1', CGvec_ph1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_v0' , CGvec_v0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_v1' , CGvec_v1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_s0' , CGvec_s0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_s1' , CGvec_s1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_sh0', CGvec_sh0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_sh1', CGvec_sh1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_z0' , CGvec_z0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_z1' , CGvec_z1)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_zh0', CGvec_zh0)
               call mpas_pool_get_array(diagnosticsPool, 'CGvec_zh1', CGvec_zh1)
   
               nCells = nCellsArray(1)
               nEdges = nEdgesArray(2)
   
               do iCell = 1, nCells
                  sshTendAx = 0.0_RKIND
   
                  do i = 1, nEdgesOnCell(iCell)
                     iEdge = edgesOnCell(i, iCell)
   
                     cell1 = cellsOnEdge(1, iEdge)
                     cell2 = cellsOnEdge(2, iEdge)
   
                     ! Interpolation sshEdge
                     sshEdgeLag = 0.5_RKIND * (sshNew(cell1) + sshNew(cell2))
   
                     ! method 1, matches method 0 without pbcs, works with pbcs.
                     thicknessSumLag = sshEdgeLag + min(bottomDepth(cell1), bottomDepth(cell2))
   
                     ! nabla (ssh^0)
                     sshDiffNew = (CGvec_wh1(cell2)-CGvec_wh1(cell1)) / dcEdge(iEdge)
   
                     fluxAx = thicknessSumLag * sshDiffNew
   
                     sshTendAx = sshTendAx + edgeSignOnCell(i, iCell) * fluxAx * dvEdge(iEdge)
                  end do ! i
                     
                  sshLagArea = (1.0_RKIND/(gravity*dt**2.0*alpha1**2.0)) * CGvec_wh1(iCell) * areaCell(iCell)
   
                  CGvec_t1(iCell) = -sshLagArea - sshTendAx
               end do ! iCell
          
               ! End reduction -----------------------------------------------------------------------!
   
               CGcst_beta0 = (CGcst_alpha0/CGcst_omega0) * CGcst_r00r1_global / CGcst_r00r0_global
               CGcst_alpha1 =  CGcst_r00r1_global / (  CGcst_r00w1_global + CGcst_beta0 * CGcst_r00s0_global  &
                                                     - CGcst_beta0 * CGcst_omega0 * CGcst_r00z0_global       )
   
               do iCell = 1,nCells
   
                  CGvec_r0(iCell) = CGvec_r1(iCell)
                  CGvec_s0(iCell) = CGvec_s1(iCell)
                  CGvec_z0(iCell) = CGvec_z1(iCell)
                  CGvec_w0(iCell) = CGvec_w1(iCell)
                  CGvec_t0(iCell) = CGvec_t1(iCell)
                  CGvec_v0(iCell) = CGvec_v1(iCell)
   
                  CGvec_rh0(iCell) = CGvec_rh1(iCell)
                  CGvec_sh0(iCell) = CGvec_sh1(iCell)
                  CGvec_ph0(iCell) = CGvec_ph1(iCell)
                  CGvec_wh0(iCell) = CGvec_wh1(iCell)
                  CGvec_zh0(iCell) = CGvec_zh1(iCell)
   
                  sshSubcycleCur(iCell) = sshSubcycleNew(iCell)
               end do ! iCell
   	    
               CGcst_r00r0_global = CGcst_r00r1_global
               CGcst_alpha0       = CGcst_alpha1
   
               block => block % next
            end do  ! block

            if ( iter > int(mean_num_cells*5) ) then
               call mpas_log_write('******************************************************')
               call mpas_log_write('Iteration number exceeds Max. #iteration: PROGRAM STOP',MPAS_LOG_CRIT)
               call mpas_log_write('Current #Iteration $i = ', intArgs=(/ iter /) )
               call mpas_log_write('Max.    #Iteration $i = ', intArgs=(/ int(mean_num_cells*5) /) )
               call mpas_log_write('******************************************************')
            endif

         !**************************************************************!
         end do ! do iter
         !**************************************************************!

         ! boundary update on SSHnew
         call mpas_timer_start("si halo iter")
         call mpas_dmpar_exch_group_create(domain, iterGroupName)
         call mpas_dmpar_exch_group_add_field(domain, iterGroupName, 'sshSubcycle', timeLevel=1)
         call mpas_threading_barrier()
         call mpas_dmpar_exch_group_full_halo_exch(domain, iterGroupName)
         call mpas_dmpar_exch_group_destroy(domain, iterGroupName)
         call mpas_timer_stop("si halo iter")

         call mpas_timer_stop("si btr iteration")


         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         ! Stage 2.6 : Barotropic velocity update
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         call mpas_timer_start("si btr vel update")

         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
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
            call mpas_pool_get_array(meshPool, 'bottomDepth',             bottomDepth            )

            call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
            call mpas_pool_get_array(statePool, 'ssh', sshNew, 2)
            call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, 1)
            call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleCur, 1)
            call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleNew, 2)
            call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityCur,1)
            call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityNew,2)

            call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)
            call mpas_pool_get_array(diagnosticsPool, 'barotropicCoriolisTerm',barotropicCoriolisTerm)
            call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux', barotropicThicknessFlux)

            nCells = nCellsArray(1)
            nEdges = nEdgesArray(1)

            do iEdge = 1, nEdges

               temp_mask = edgeMask(1, iEdge)

               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)

               ! Compute the barotropic Coriolis term, -f*uPerp
               CoriolisTerm = 0.0_RKIND

               do i = 1, nEdgesOnEdge(iEdge)
                  eoe = edgesOnEdge(i,iEdge)
                  wgt = weightsOnEdge(i,iEdge)
                  CoriolisTerm = Coriolisterm + wgt  &
                               * normalBarotropicVelocitySubcycleCur(eoe) * fEdge(eoe)
               end do
               barotropicCoriolisTerm(iEdge) = CoriolisTerm

               normalBarotropicVelocitySubcycleNew(iEdge) &
                  = temp_mask &
                  * (normalBarotropicVelocityCur(iEdge) &
                  - (dt) * (-barotropicCoriolisTerm(iEdge) + gravity &
                  * ( (alpha1*sshSubcycleCur(cell2)+alpha2*sshCur(cell2)) - (alpha1*sshSubcycleCur(cell1)+alpha2*sshCur(cell1)) ) &
                  / dcEdge(iEdge) - barotropicForcing(iEdge)))

            end do ! iEdge

            do iEdge = 1, nEdges
               normalBarotropicVelocitySubcycleCur(iEdge) =   0.5*normalBarotropicVelocitySubcycleNew(iEdge)  &
                                                            + 0.5*normalBarotropicVelocityCur(iEdge)
            end do ! iEdge


            do iEdge = 1, nEdges
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)

               sshEdge =  alpha2*( 0.5_RKIND * (sshCur(cell1) + sshCur(cell2))) &
                        + alpha1*( 0.5_RKIND * (sshSubcycleCur(cell1) + sshSubcycleCur(cell2)))

               thicknessSum = sshEdge + min(bottomDepth(cell1), bottomDepth(cell2))

               flux =normalBarotropicVelocitySubcycleCur(iEdge) * thicknessSum

               barotropicThicknessFlux(iEdge) = flux

               normalBarotropicVelocityNew(iEdge) = normalBarotropicVelocitySubcycleCur(iEdge)
            end do ! iEdge


            do iCell = 1,nCells
               sshNew(iCell) = sshSubcycleCur(iCell) 
            end do

            ! #Iteration  check
!           if ( split_implicit_step == 2 ) then
!           call mpas_log_write( 'Iteration $i', intArgs=(/ iter /) )
!           endif

            block => block % next
         end do  ! block

         ! boundary update on F
         call mpas_timer_start("si halo btr vel")
         call mpas_dmpar_exch_group_create(domain, finalBtrGroupName)
         call mpas_dmpar_exch_group_add_field(domain, finalBtrGroupName, 'barotropicThicknessFlux')
         call mpas_dmpar_exch_group_add_field(domain, finalBtrGroupName, 'normalBarotropicVelocity',2)
         call mpas_dmpar_exch_group_add_field(domain, finalBtrGroupName, 'ssh',2)
         call mpas_threading_barrier()
         call mpas_dmpar_exch_group_full_halo_exch(domain, finalBtrGroupName)
         call mpas_dmpar_exch_group_destroy(domain, finalBtrGroupName)
         call mpas_timer_stop("si halo btr vel")
         call mpas_threading_barrier()

         call mpas_timer_stop("si btr vel update")

         call mpas_threading_barrier()

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         ! u correct and transport velocity
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! Check that you can compute SSH using the total sum or the individual increments
         ! over the barotropic subcycles.
         ! efficiency: This next block of code is really a check for debugging, and can
         ! be removed later.
         call mpas_timer_start('btr si ssh verif')

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

            !$omp do schedule(runtime)
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
         call mpas_timer_stop('btr si ssh verif')

         call mpas_timer_stop("si btr vel")


         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         !  Stage 3: Tracer, density, pressure, vertical velocity prediction
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! only compute tendencies for active tracers on last large iteration
         if (split_implicit_step < config_n_ts_iter) then
            activeTracersOnly = .true.
         else
            activeTracersOnly = .false.
         endif

         ! Thickness tendency computations and thickness halo updates are completed before tracer
         ! tendency computations to allow monotonic advection.
         call mpas_timer_start('si thick tend')
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
         call mpas_timer_stop('si thick tend')

         ! update halo for thickness tendencies
         call mpas_timer_start("si halo thickness")

         call mpas_dmpar_field_halo_exch(domain, 'tendLayerThickness')

         call mpas_timer_stop("si halo thickness")

         call mpas_timer_start('si tracer tend', .false.)
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
            call ocn_tend_tracer(tendPool, statePool, forcingPool, diagnosticsPool, meshPool, swForcingPool, scratchPool, &
                    dt, activeTracersOnly, 2)

            block => block % next
         end do
         call mpas_timer_stop('si tracer tend')

         ! update halo for tracer tendencies
         call mpas_timer_start("si halo tracers")
         call mpas_pool_get_subpool(domain % blocklist % structs, 'tend', tendPool)
         call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)

         call mpas_pool_begin_iteration(tracersTendPool)
         do while ( mpas_pool_get_next_member(tracersTendPool, groupItr) )
            if ( groupItr % memberType == MPAS_POOL_FIELD ) then
               ! Only compute tendencies for active tracers if activeTracersOnly flag is true.
               if ( .not.activeTracersOnly .or. trim(groupItr % memberName)=='activeTracersTend') then
                  call mpas_dmpar_field_halo_exch(domain, groupItr % memberName)
               end if
            end if
         end do
         call mpas_timer_stop("si halo tracers")

         call mpas_timer_start('si loop fini')
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
            if (split_implicit_step < config_n_ts_iter) then

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
            elseif (split_implicit_step == config_n_ts_iter) then

               !$omp do schedule(runtime) private(k)
               do iCell = 1, nCells
                  do k = 1, maxLevelCell(iCell)
                     ! this is h_{n+1}
                     layerThicknessNew(k,iCell) = layerThicknessCur(k,iCell) + dt * layerThicknessTend(k,iCell)
                  end do
               end do
               !$omp end do

               if (config_compute_active_tracer_budgets) then
                  call mpas_pool_get_array(diagnosticsPool,'activeTracerHorizontalAdvectionTendency', &
                          activeTracerHorizontalAdvectionTendency)
                  call mpas_pool_get_array(diagnosticspool,'activeTracerVerticalAdvectionTendency', &
                          activeTracerVerticalAdvectionTendency)
                  call mpas_pool_get_array(diagnosticsPool,'activeTracerSurfaceFluxTendency',activeTracerSurfaceFluxTendency)
                  call mpas_pool_get_array(diagnosticsPool,'temperatureShortWaveTendency',temperatureShortWaveTendency)
                  call mpas_pool_get_array(diagnosticsPool,'activeTracerNonLocalTendency',activeTracerNonLocalTendency)
                  call mpas_pool_get_array(diagnosticsPool,'activeTracerHorizontalAdvectionEdgeFlux', &
                          activeTracerHorizontalAdvectionEdgeFlux)
                  call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)

                  !$omp do schedule(runtime) private(k)
                  do iEdge = 1, nEdges
                     do k= 1, maxLevelEdgeTop(iEdge)
                        activeTracerHorizontalAdvectionEdgeFlux(:,k,iEdge) = &
                          activeTracerHorizontalAdvectionEdgeFlux(:,k,iEdge) / &
                          layerThicknessEdge(k,iEdge)
                     enddo
                  enddo
                  !$omp end do

                  !$omp do schedule(runtime) private(k)
                  do iCell = 1, nCells
                     do k= 1, maxLevelCell(iCell)
                        activeTracerHorizontalAdvectionTendency(:,k,iCell) = &
                           activeTracerHorizontalAdvectionTendency(:,k,iCell) / &
                           layerThicknessNew(k,iCell)

                        activeTracerVerticalAdvectionTendency(:,k,iCell) = &
                           activeTracerVerticalAdvectionTendency(:,k,iCell) / &
                           layerThicknessNew(k,iCell)

                        activeTracerSurfaceFluxTendency(:,k,iCell) = &
                           activeTracerSurfaceFluxTendency(:,k,iCell) / &
                           layerThicknessNew(k,iCell)

                        temperatureShortWaveTendency(k,iCell) = &
                           temperatureShortWaveTendency(k,iCell) / &
                           layerThicknessNew(k,iCell)

                        activeTracerNonLocalTendency(:,k,iCell) = &
                           activeTracerNonLocalTendency(:,k,iCell) / &
                           layerThicknessNew(k,iCell)
                     end do
                  end do
                  !$omp end do
               endif

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
                           end do
                        end do
                        !$omp end do

                        ! limit salinity in separate loop
                        if ( trim(groupItr % memberName) == 'activeTracers' ) then
                           !$omp do schedule(runtime) private(k)
                           do iCell = 1, nCells
                              do k = 1, maxLevelCell(iCell)
                                 tracersGroupNew(indexSalinity,k,iCell) = max(0.001_RKIND, tracersGroupNew(indexSalinity,k,iCell))
                              end do
                           end do
                           !$omp end do
                        end if

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

            endif ! split_implicit_step

            block => block % next
         end do

         call mpas_timer_stop('si loop fini')
         call mpas_timer_stop('si loop')

      end do  ! split_implicit_step = 1, config_n_ts_iter
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! END large iteration loop
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call mpas_timer_start("si implicit vert mix")

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

        block => block % next
      end do
      call mpas_threading_barrier()

      call mpas_dmpar_field_halo_exch(domain, 'surfaceFrictionVelocity')

      block => domain % blocklist
      do while(associated(block))
        call mpas_pool_get_subpool(block % structs, 'state', statePool)
        call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
        call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
        call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

        ! Compute normalGMBolusVelocity; it will be added to the baroclinic modes in Stage 2 above.
        if (config_use_standardGM) then
           call ocn_gm_compute_Bolus_velocity(diagnosticsPool, meshPool, scratchPool)
        end if
        call ocn_vmix_implicit(dt, meshPool, diagnosticsPool, statePool, forcingPool, scratchPool, err, 2)

        block => block % next
      end do

      ! Update halo on u and tracers, which were just updated for implicit vertical mixing.  If not done,
      ! this leads to lack of volume conservation.  It is required because halo updates in stage 3 are only
      ! conducted on tendencies, not on the velocity and tracer fields.  So this update is required to
      ! communicate the change due to implicit vertical mixing across the boundary.
      call mpas_timer_start('si vmix halos')
      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)


      call mpas_timer_start('si vmix halos normalVelFld')
      call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', timeLevel=2)
      call mpas_timer_stop('si vmix halos normalVelFld')

      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
            call mpas_dmpar_field_halo_exch(domain, groupItr % memberName, timeLevel=2)
         end if
      end do
      call mpas_timer_stop('si vmix halos')

      call mpas_timer_stop("si implicit vert mix")

      call mpas_timer_start('si fini')
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

         call mpas_timer_start('si final mpas reconstruct', .false.)

         call mpas_reconstruct(meshPool, normalVelocityNew,  &
                          velocityX, velocityY, velocityZ,   &
                          velocityZonal, velocityMeridional, &
                          includeHalos = .true.)

         call mpas_reconstruct(meshPool, gradSSH,          &
                          gradSSHX, gradSSHY, gradSSHZ,    &
                          gradSSHZonal, gradSSHMeridional, &
                          includeHalos = .true.)

         call mpas_timer_stop('si final mpas reconstruct')

         !$omp do schedule(runtime)
         do iCell = 1, nCells
            surfaceVelocity(indexSurfaceVelocityZonal, iCell) = velocityZonal(1, iCell)
            surfaceVelocity(indexSurfaceVelocityMeridional, iCell) = velocityMeridional(1, iCell)

            SSHGradient(indexSSHGradientZonal, iCell) = gradSSHZonal(iCell)
            SSHGradient(indexSSHGradientMeridional, iCell) = gradSSHMeridional(iCell)
         end do
         !$omp end do

         call ocn_time_average_coupled_accumulate(diagnosticsPool, statePool, forcingPool, 2)

         if (config_use_standardGM) then
            call ocn_reconstruct_gm_vectors(diagnosticsPool, meshPool)
         end if

         block => block % next
      end do

      if (trim(config_land_ice_flux_mode) == 'coupled') then
         call mpas_timer_start("si effective density halo")
         call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
         call mpas_pool_get_field(statePool, 'effectiveDensityInLandIce', effectiveDensityField, 2)
         call mpas_dmpar_exch_halo_field(effectiveDensityField)
         call mpas_timer_stop("se effective density halo")
      end if

      call mpas_timer_stop('si fini')
      call mpas_timer_stop("si timestep")

      deallocate(n_bcl_iter)

   end subroutine ocn_time_integrator_si!}}}

!***********************************************************************
!
!  routine ocn_time_integration_si_init
!
!> \brief   Initialize split-explicit time stepping within MPAS-Ocean core
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  This routine initializes variables required for the split-explicit time
!>  stepper.
!
!-----------------------------------------------------------------------
   subroutine ocn_time_integration_si_init(domain)!{{{
   ! Initialize splitting variables

      type (domain_type), intent(inout) :: domain

      integer :: i, iCell, iEdge, iVertex, k
      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: statePool, meshPool, tracersPool
      type (dm_info) :: dminfo

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
      integer (kind=I8KIND) :: nBtrSubcyclesI8

      integer, dimension(:), pointer :: nCellsArray, nEdgesArray
      real (kind=RKIND) :: local_num_cells,sum1,total_area_sum,local_area_sum
      real (kind=RKIND), dimension(:), pointer :: areaCell
      integer :: ncpus

      dminfo = domain % dminfo
      ncpus = dminfo % nprocs

      call mpas_pool_get_config(domain % configs, 'config_do_restart', config_do_restart)

      ! Determine the number of barotropic subcycles based on the ratio of time steps
      call mpas_pool_get_config(domain % configs, 'config_time_integrator', config_time_integrator)
      call mpas_pool_get_config(domain % configs, 'config_btr_dt', config_btr_dt)
      call mpas_pool_get_config(domain % configs, 'config_dt', config_dt)

      nowTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
      call mpas_set_timeInterval( zeroInterval, S=0 )

      call mpas_set_timeInterval( fullTimeStep , timeString=config_dt )
      call mpas_set_timeInterval( barotropicTimeStep, timeString=config_btr_dt )

      ! transfer to I8 for division step
      nBtrSubcyclesI8 = nBtrSubcycles
      call mpas_interval_division( nowTime, fullTimeStep, barotropicTimeStep, nBtrSubcyclesI8, remainder )
      nBtrSubcycles = nBtrSubcyclesI8

      if ( remainder > zeroInterval ) then
         nBtrSubcycles = nBtrSubcycles + 1
      end if

      call mpas_log_write( '*******************************************************************************')
      call mpas_log_write( 'The semi-implicit time integration is configured')
      call mpas_log_write( '*******************************************************************************')

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
            call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
            call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)

            call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
            call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, 1)
            call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocity, 1)
            call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', normalBaroclinicVelocity, 1)

            call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
            call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
            call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
            call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

           

            ! Compute barotropic velocity at first timestep
            ! This is only done upon start-up.
            if (trim(config_time_integrator) == 'semi_implicit') then

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

              !--- RMS of mean area
              local_num_cells = nCellsArray(1)
              call mpas_dmpar_sum_real(dminfo,local_num_cells,total_num_cells)

              local_area_sum = 0.0_RKIND
              do iCell = 1,nCellsArray(1)
                local_area_sum = local_area_sum + areaCell(iCell)**2.0
              end do
              call mpas_dmpar_sum_real(dminfo,local_area_sum,total_area_sum)

              area_mean = dsqrt(total_area_sum / total_num_cells)
              mean_num_cells = total_num_cells/ncpus

         block => block % next
         end do
      end if

   end subroutine ocn_time_integration_si_init!}}}

!***********************************************************************
!
!  routine ocn_time_integration_si_preconditioner
!
!> \brief   Construct a Block-Jacobi preconditioning matrix
!> \author  Hyun-Gyu Kang (Oak Ridge National Laboratory)
!> \date    September 2019
!> \details
!>  This routine constructs a Block-Jacobi preconditioner for the
!>  split-implicit time stepper.
!
!-----------------------------------------------------------------------
   subroutine ocn_time_integrator_si_preconditioner(domain, dt)!{{{

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND) , intent(in)    :: dt

      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: tendPool
      type (mpas_pool_type), pointer :: tracersTendPool

      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      real (kind=RKIND) :: thicknessSum,fluxAx
      integer :: iCell, i,k,j, iEdge, cell1, cell2
      character (len=StrKIND), pointer :: config_btr_si_preconditioner

      ! Dimensions
      integer :: nCells, nEdges
      integer, pointer :: nCellsPtr, nEdgesPtr
      integer, dimension(:), pointer :: nCellsArray, nEdgesArray

      ! Mesh array pointers
      integer, dimension(:)  , pointer :: maxLevelCell, maxLevelEdgeTop, nEdgesOnEdge, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, edgeMask, edgesOnEdge
      integer, dimension(:,:), pointer :: edgesOnCell, edgeSignOnCell

      real (kind=RKIND), dimension(:)  , pointer :: dcEdge,bottomDepth
      real (kind=RKIND), dimension(:)  , pointer :: dvEdge,areaCell
      real (kind=RKIND), dimension(:)  , pointer :: CGvec_u1
      real (kind=RKIND) :: temp1,temp2
      !------------------------------------------!
      integer,dimension(:)  ,pointer     :: globalCellId
      integer,dimension(:)  ,allocatable :: numCellRank,numCellRankA2,myGlobalCellId
      integer,dimension(:,:),allocatable :: rankCellId,myLocalStartEndId,iTmpArray
      integer :: rank,ncpus,mpi_ierr
      integer :: local_num_cells,total_num_cells,itmp1,itmp2
      integer :: maxNumCell,maxNumCellA2,isum1,isum2
      integer :: local_start,local_end,nCellsA2,nCellsA3
      !include 'mpif.h'
      !------------------------------------------!

      dminfo = domain % dminfo
        rank = dminfo % my_proc_id
       ncpus = dminfo % nprocs

      call mpas_pool_get_config(domain % configs, 'config_btr_si_preconditioner', config_btr_si_preconditioner)


      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_dimension(block % dimensions, 'nCells'     , nCellsPtr  )
         call mpas_pool_get_dimension(block % dimensions, 'nEdges'     , nEdgesPtr  )
         call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
         call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
 
         call mpas_pool_get_subpool(block % structs, 'tend'       , tendPool       )
         call mpas_pool_get_subpool(block % structs, 'mesh'       , meshPool       )
         call mpas_pool_get_subpool(block % structs, 'state'      , statePool      )
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
 
         call mpas_pool_get_subpool(statePool, 'tracers'    , tracersPool    )
         call mpas_pool_get_subpool(tendPool , 'tracersTend', tracersTendPool)
 
         call mpas_pool_get_array(meshPool, 'nEdgesOnCell',    nEdgesOnCell    )
         call mpas_pool_get_array(meshPool, 'edgesOnCell',     edgesOnCell     )
         call mpas_pool_get_array(meshPool, 'cellsOnEdge',     cellsOnEdge     )
         call mpas_pool_get_array(meshPool, 'dcEdge',          dcEdge          )
         call mpas_pool_get_array(meshPool, 'bottomDepth',     bottomDepth     )
         call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop )
         call mpas_pool_get_array(meshPool, 'edgeSignOnCell',  edgeSignOnCell  )
         call mpas_pool_get_array(meshPool, 'dvEdge',          dvEdge          )
         call mpas_pool_get_array(meshPool, 'areaCell',        areaCell        )
         call mpas_pool_get_array(meshPool, 'nEdgesOnEdge',    nEdgesOnEdge    )
         call mpas_pool_get_array(meshPool, 'edgesOnEdge',     edgesOnEdge     )
         call mpas_pool_get_array(meshPool, 'indexToCellID',   globalCellId    )
 
         nCells   = nCellsArray(1)
         nEdges   = nEdgesArray(1)
         nCellsA2 = nCellsArray(2)
         nCellsA3 = nCellsArray(3)
 
         call mpas_log_write('   config_btr_si_preconditioner: ' // trim(config_btr_si_preconditioner))
 
         local_num_cells = nCellsArray(1)
         call mpas_dmpar_sum_int(dminfo,local_num_cells,total_num_cells)


         ! Restricted Additive Schwarz preconditioner ---------------------------------------------!
         if ( trim(config_btr_si_preconditioner) == 'ras' ) then
 
!           nPrecVec = nCellsA2
            nPrecVec = nCellsA3
  
            allocate(prec_ivmat(1:nPrecVec,1:nPrecVec))
                     prec_ivmat(:,:) = 0.0_RKIND
  
            do iCell = 1, nPrecVec
  
               do i = 1, nEdgesOnCell(iCell)
                  iEdge = edgesOnCell(i, iCell)
                  cell1 = cellsOnEdge(1, iEdge)
                  cell2 = cellsOnEdge(2, iEdge)
    
                  ! method 1, matches method 0 without pbcs, works with pbcs.
                  thicknessSum = min(bottomDepth(cell1), bottomDepth(cell2))
                        fluxAx = edgeSignOnCell(i,iCell)*dvEdge(iEdge)*thicknessSum / dcEdge(iEdge)
    
                  !-------------------------------------------------------------------------------!
                  if ( globalCellId(cell1) > 0 .and. globalCellId(cell1) < total_num_cells+1 ) then
                     if ( cell1 <= nPrecVec) then
                        prec_ivmat(iCell,cell1) = prec_ivmat(iCell,cell1) + fluxAx
                     endif
                  endif
    
                  if ( globalCellId(cell2) > 0 .and. globalCellId(cell2) < total_num_cells+1 ) then
                     if ( cell2 <= nPrecVec) then
                        prec_ivmat(iCell,cell2) = prec_ivmat(iCell,cell2) - fluxAx
                     endif
                  endif
                  !-------------------------------------------------------------------------------!
   
               end do ! i
  
               prec_ivmat(iCell,iCell) = prec_ivmat(iCell,iCell) - (4.0_RKIND/(gravity*dt**2.0)) * areaCell(iCell)
  
            end do ! iCell
  
           ! Preconditioning matrix inversion
            call inverse(nPrecVec,prec_ivmat)
  

         ! Block-Jacobi preconditioner ------------------------------------------------------------!
         elseif ( trim(config_btr_si_preconditioner) == 'block_jacobi' ) then
 
            nPrecVec = nCells ! length of preconditioning vector
  
            allocate(prec_ivmat(1:nPrecVec,1:nPrecVec))
                     prec_ivmat(:,:) = 0.0_RKIND
  
            do iCell = 1, nPrecVec
  
               do i = 1, nEdgesOnCell(iCell)
                  iEdge = edgesOnCell(i, iCell)
                  cell1 = cellsOnEdge(1, iEdge)
                  cell2 = cellsOnEdge(2, iEdge)
    
                  ! method 1, matches method 0 without pbcs, works with pbcs.
                  thicknessSum = min(bottomDepth(cell1), bottomDepth(cell2))
                        fluxAx = edgeSignOnCell(i,iCell)*dvEdge(iEdge)*thicknessSum / dcEdge(iEdge)
    
                  !-------------------------------------------------------------------------------!
                  if ( globalCellId(cell1) > 0 .and. globalCellId(cell1) < total_num_cells+1 ) then
                     if ( cell1 <= nPrecVec) then
                        prec_ivmat(iCell,cell1) = prec_ivmat(iCell,cell1) + fluxAx
                     endif
                  endif
    
                  if ( globalCellId(cell2) > 0 .and. globalCellId(cell2) < total_num_cells+1 ) then
                     if ( cell2 <= nPrecVec) then
                        prec_ivmat(iCell,cell2) = prec_ivmat(iCell,cell2) - fluxAx
                     endif
                  endif
                  !-------------------------------------------------------------------------------!
               end do ! i
   
               prec_ivmat(iCell,iCell) = prec_ivmat(iCell,iCell) - (4.0_RKIND/(gravity*dt**2.0))*areaCell(iCell)
  
            end do ! iCell
  
            ! Preconditioning matrix inversion
            call inverse(nPrecVec,prec_ivmat)
 

         ! Jacobi preconditioner ------------------------------------------------------------------!
         else if ( trim(config_btr_si_preconditioner) == 'jacobi' ) then
 
            nPrecVec = nCells ! length of preconditioning vector
    
            allocate(prec_ivmat(1:nPrecVec,1))
                     prec_ivmat(:,:) = 0.0_RKIND
    
            do iCell = 1, nPrecVec
  
               do i = 1, nEdgesOnCell(iCell)
                  iEdge = edgesOnCell(i, iCell)
                  cell1 = cellsOnEdge(1, iEdge)
                  cell2 = cellsOnEdge(2, iEdge)
  
                  ! method 1, matches method 0 without pbcs, works with pbcs.
                  thicknessSum =  min(bottomDepth(cell1), bottomDepth(cell2))
  
                  !-------------------------------------------------------------------------!
                  fluxAx = edgeSignOnCell(i,iCell)*dvEdge(iEdge)*thicknessSum / dcEdge(iEdge)
  
                  if (cell1 == iCell) then
                     prec_ivmat(iCell,1) = prec_ivmat(iCell,1) + fluxAx  ! reversed sign
                  elseif ( cell2 == iCell) then
                     prec_ivmat(iCell,1) = prec_ivmat(iCell,1) - fluxAx  ! reversed sign
                  endif
                  !-------------------------------------------------------------------------!
               end do ! i
  
               temp1 = prec_ivmat(iCell,1) - (4.0_RKIND/(gravity*dt**2.0))*areaCell(iCell)
  
               prec_ivmat(iCell,1) = 1.0_RKIND / temp1
  
            end do ! iCell
   
   
         ! No preconditioner ----------------------------------------------------------------------!
         else if ( trim(config_btr_si_preconditioner) == 'none' ) then
 
            nPrecVec = nCells ! length of preconditioning vector
 
            allocate(prec_ivmat(1,1))
            prec_ivmat(:,:) = 1.0_RKIND ! This array is not used only for 'none'.
 
         else
 
            call mpas_log_write('Incorrect choice for config_btr_si_preconditioner: ' // trim(config_btr_si_preconditioner) // &
                                '   choices are: ras, block_jacobi, jacobi, none',MPAS_LOG_CRIT)
 
         endif ! config_btr_si_preconditioner
 
         block => block % next
      end do  ! block

   end subroutine ocn_time_integrator_si_preconditioner !}}}

end module ocn_time_integration_si

!***********************************************************************
!
!  routine inverse
!
!> \brief   Inverte a square matrix using the LAPACK
!> \author  Hyun-Gyu Kang (Oak Ridge National Laboratory)
!> \date    September 2019
!> \details
!>  This routine inverts a square matrix.
!>  Matrix (A) is destroyed and changed to its inversion.
!
!----------------------------------------------------------------------
subroutine inverse (n,Ainv) !{{{
    implicit none
    real(kind=8),dimension(n,n) :: Ainv
    real(kind=8),dimension(n)   :: work
    integer     ,dimension(n)   :: ipiv
    integer                     :: n,info,i,j

    ! (-) sign for positivie definite matrix
    Ainv = -Ainv

    ! DPOTRF computes the Cholesky factorizaton of a symmetric positive
    ! defnitie matrix.
    call DPOTRF('U',n,Ainv,n,info)
!   call DGETRF(n,n,Ainv,n,ipiv,info) ! For unsymmetric matrix

    if (info.lt.0) stop 'ERROR in inverse routine: Matrix is numerically singular!'
    if (info.gt.0) stop 'ERROR in inverse routine: Matrix is not positive definite!'

    ! DPOTRI computes the inverse of a symmetric positive definite matrix,
    ! using the Cholesky factorization computed by DPOTRF
    call DPOTRI('U',n,Ainv,n,info)
!   call DGETRI(n,Ainv,n,ipiv,work,n,info) ! For unsymmetric matrix

    Ainv = -Ainv ! (-) sign to restore a matrix sign
    if (info.ne.0) stop 'ERROR in inverse routine: Matrix inversion failed!'
end subroutine inverse !}}}

! vim: foldmethod=marker
