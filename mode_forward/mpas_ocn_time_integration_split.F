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

   use mpas_grid_types
   use mpas_constants
   use mpas_dmpar
   use mpas_vector_reconstruction
   use mpas_spline_interpolation
   use mpas_timer

   use ocn_tendency
   use ocn_diagnostics

   use ocn_equation_of_state
   use ocn_vmix
   use ocn_time_average
   use ocn_time_average_coupled

   use ocn_sea_ice

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

   public :: ocn_time_integrator_split

   type (timer_node), pointer :: timer_main, timer_prep, timer_bcl_vel, timer_btr_vel, timer_diagnostic_update, timer_implicit_vmix, &
                                 timer_halo_diagnostic, timer_halo_normalBarotropicVelocity, timer_halo_ssh, timer_halo_f, timer_halo_thickness, & 
                                 timer_halo_tracers, timer_halo_normalBaroclinicVelocity

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_split
!
!> \brief MPAS ocean split explicit time integration scheme
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This routine integrates a single time step (dt) using a
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
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: verticalMeshPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: tendPool
      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: averagePool
      type (mpas_pool_type), pointer :: scratchPool

      type (dm_info) :: dminfo
      integer :: iCell, i,k,j, iEdge, cell1, cell2, split_explicit_step, split, &
                 eoe, oldBtrSubcycleTime, newBtrSubcycleTime, uPerpTime, BtrCorIter, &
                 stage1_tend_time
      integer, dimension(:), allocatable :: n_bcl_iter
      type (block_type), pointer :: block
      real (kind=RKIND) :: normalThicknessFluxSum, thicknessSum, flux, sshEdge, hEdge1, &
                 CoriolisTerm, uCorr, temp, temp_h, coef, barotropicThicknessFlux_coeff, sshCell1, sshCell2
      integer :: ucorr_coef, err
      real (kind=RKIND), dimension(:,:), pointer :: &
                 vertViscTopOfEdge, vertDiffTopOfCell
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers
      real (kind=RKIND), dimension(:), allocatable:: uTemp
      real (kind=RKIND), dimension(:,:), allocatable:: tracersTemp

      integer :: tsIter

      ! Config options
      character (len=StrKIND), pointer :: config_time_integrator
      integer, pointer :: config_n_bcl_iter_mid, config_n_bcl_iter_beg, config_n_bcl_iter_end
      integer, pointer :: config_n_ts_iter, config_btr_subcycle_loop_factor, config_n_btr_subcycles
      integer, pointer :: config_n_btr_cor_iter

      logical, pointer :: config_use_freq_filtered_thickness, config_btr_solve_SSH2, config_filter_btr_mode
      logical, pointer :: config_vel_correction, config_prescribe_velocity, config_prescribe_thickness

      real (kind=RKIND), pointer :: config_mom_del4, config_btr_gam1_velWt1, config_btr_gam2_SSHWt1
      real (kind=RKIND), pointer :: config_btr_gam3_velWt2

      ! Dimensions
      integer, pointer :: nCells, nEdges, nVertLevels, num_tracers, startIndex, endIndex
      integer, pointer :: indexTemperature, indexSalinity
      integer, pointer :: indexZonalSurfaceVelocity, indexMeridionalSurfaceVelocity
      integer, pointer :: indexZonalSSHGradient, indexMeridionalSSHGradient

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
      real (kind=RKIND), dimension(:,:,:), pointer :: tracersCur, tracersNew

      ! Tend Array Pointers
      real (kind=RKIND), dimension(:), pointer :: sshTend
      real (kind=RKIND), dimension(:,:), pointer :: highFreqThicknessTend
      real (kind=RKIND), dimension(:,:), pointer :: lowFreqDivergenceTend
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocityTend, layerThicknessTend
      real (kind=RKIND), dimension(:,:,:), pointer :: tracersTend

      ! Diagnostics Array Pointers
      real (kind=RKIND), dimension(:), pointer :: barotropicForcing, barotropicThicknessFlux
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessEdge, uTransport, uBolusGM
      real (kind=RKIND), dimension(:,:), pointer :: vertTransportVelocityTop
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocityX, normalVelocityY, normalVelocityZ
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocityZonal, normalVelocityMeridional
      real (kind=RKIND), dimension(:,:), pointer :: gradSSH
      real (kind=RKIND), dimension(:,:), pointer :: gradSSHX, gradSSHY, gradSSHZ
      real (kind=RKIND), dimension(:,:), pointer :: gradSSHZonal, gradSSHMeridional
      real (kind=RKIND), dimension(:,:), pointer :: surfaceVelocity, SSHGradient

      ! Forcing Array Pointer
      real (kind=RKIND), dimension(:), pointer :: seaIceEnergy

      ! Diagnostics Field Pointers
      type (field2DReal), pointer :: normalizedRelativeVorticityEdgeField, divergenceField, relativeVorticityField
      type (field1DReal), pointer :: barotropicThicknessFluxField

      ! State/Tend Field Pointers
      type (field1DReal), pointer :: normalBarotropicVelocitySubcycleField, sshSubcycleField
      type (field2DReal), pointer :: highFreqThicknessField, lowFreqDivergenceField
      type (field2DReal), pointer :: normalBaroclinicVelocityField, layerThicknessField
      type (field2DReal), pointer :: normalVelocityField
      type (field3DReal), pointer :: tracersField

      call mpas_timer_start("se timestep", .false., timer_main)

      call mpas_pool_get_config(domain % configs, 'config_n_bcl_iter_beg', config_n_bcl_iter_beg)
      call mpas_pool_get_config(domain % configs, 'config_n_bcl_iter_mid', config_n_bcl_iter_mid)
      call mpas_pool_get_config(domain % configs, 'config_n_bcl_iter_end', config_n_bcl_iter_end)
      call mpas_pool_get_config(domain % configs, 'config_n_ts_iter', config_n_ts_iter)
      call mpas_pool_get_config(domain % configs, 'config_n_btr_subcycles', config_n_btr_subcycles)
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

      allocate(n_bcl_iter(config_n_ts_iter))

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      !  Prep variables before first iteration
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call mpas_timer_start("se prep", .false., timer_prep)
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
         call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)

         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

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

         call mpas_pool_get_array(statePool, 'tracers', tracersCur, 1)
         call mpas_pool_get_array(statePool, 'tracers', tracersNew, 2)

         call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessCur, 1)
         call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessNew, 2)

         call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceCur, 1)
         call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergenceNew, 2)

         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

         ! Initialize * variables that are used to compute baroclinic tendencies below.
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

               ! DWJ-POOL What's this for?
!                block % diagnostics % layerThicknessEdge % array(k,iEdge) &
!              = block % diagnostics % layerThicknessEdge % array(k,iEdge)
            end do 
         end do 

         sshNew(:) = sshCur(:)

         do iCell = 1, nCells  
            do k = 1, maxLevelCell(iCell)
               layerThicknessNew(k,iCell) = layerThicknessCur(k,iCell)

               tracersNew(:,k,iCell) = tracersCur(:,k,iCell) 
            end do
         end do

         if (associated(highFreqThicknessNew)) then
            highFreqThicknessNew(:,:) = highFreqThicknessCur(:,:)
         end if

         if (associated(lowFreqDivergenceNew)) then
            lowFreqDivergenceNew(:,:) = lowFreqDivergenceCur(:,:)
         endif

         block => block % next
      end do

      call mpas_timer_stop("se prep", timer_prep)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! BEGIN large iteration loop 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      n_bcl_iter = config_n_bcl_iter_mid
      n_bcl_iter(1) = config_n_bcl_iter_beg
      n_bcl_iter(config_n_ts_iter) = config_n_bcl_iter_end

      do split_explicit_step = 1, config_n_ts_iter
         stage1_tend_time = min(split_explicit_step,2)

         ! ---  update halos for diagnostic variables
         call mpas_timer_start("se halo diag", .false., timer_halo_diagnostic)

         call mpas_pool_get_subpool(domain % blocklist % structs, 'diagnostics', diagnosticsPool)

         call mpas_pool_get_field(diagnosticsPool, 'normalizedRelativeVorticityEdge', normalizedRelativeVorticityEdgeField)
         call mpas_pool_get_field(diagnosticsPool, 'divergence', divergenceField)
         call mpas_pool_get_field(diagnosticsPool, 'relativeVorticity', relativeVorticityField)
         call mpas_dmpar_exch_halo_field(normalizedRelativeVorticityEdgeField)
         if (config_mom_del4 > 0.0) then
           call mpas_dmpar_exch_halo_field(divergenceField)
           call mpas_dmpar_exch_halo_field(relativeVorticityField)
         end if
         call mpas_timer_stop("se halo diag", timer_halo_diagnostic)

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
               call mpas_pool_get_subpool(block % structs, 'state', statepool)
               call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
               call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

               call ocn_tend_freq_filtered_thickness(tendPool, statePool, diagnosticsPool, meshPool, stage1_tend_time)
               block => block % next
            end do
            call mpas_timer_stop("se freq-filtered-thick computations")

            call mpas_timer_start("se freq-filtered-thick halo update")
            call mpas_pool_get_subpool(domain % blocklist % structs, 'tend', tendPool)

            call mpas_pool_get_field(tendPool, 'highFreqThickness', highFreqThicknessField)
            call mpas_pool_get_field(tendPool, 'lowFreqDivergence', lowFreqDivergenceField)

            call mpas_dmpar_exch_halo_field(highFreqThicknessField)
            call mpas_dmpar_exch_halo_field(lowFreqDivergenceField)
            call mpas_timer_stop("se freq-filtered-thick halo update")

            block => domain % blocklist
            do while (associated(block))
               call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)

               call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
               call mpas_pool_get_subpool(block % structs, 'state', statePool)
               call mpas_pool_get_subpool(block % structs, 'tend', tendPool)

               call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

               call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessCur, 1)
               call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessNew, 2)

               call mpas_pool_get_array(tendPool, 'highFreqThickness', highFreqThicknessTend)

               do iCell = 1, nCells
                  do k = 1, maxLevelCell(iCell)
                     ! this is h^{hf}_{n+1}
                     highFreqThicknessNew(k,iCell) = highFreqThicknessCur(k,iCell) + dt * highFreqThicknessTend(k,iCell) 
                  end do
               end do
               block => block % next
            end do

         endif


         ! compute velocity tendencies, T(u*,w*,p*)
         call mpas_timer_start("se bcl vel", .false., timer_bcl_vel)

         block => domain % blocklist
         do while (associated(block))
           call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
           call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
           call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
           call mpas_pool_get_subpool(block % structs, 'state', statePool)
           call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
           call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
           call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)

           call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
           call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, stage1_tend_time)
           call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)

           call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessNew, 2)

           call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
           call mpas_pool_get_array(diagnosticsPool, 'vertTransportVelocityTop', vertTransportVelocityTop)

           ! compute vertTransportVelocityTop.  Use u (rather than uTransport) for momentum advection.
           ! Use the most recent time level available.
           call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, &
              layerThicknessCur, layerThicknessEdge, normalVelocityCur, &
              sshCur, highFreqThicknessNew, dt, vertTransportVelocityTop, err)

            call ocn_tend_vel(tendPool, statePool, forcingPool, diagnosticsPool, meshPool, scratchPool, stage1_tend_time)

            block => block % next
         end do

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

            block => domain % blocklist
            do while (associated(block))
               call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
               call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)

               call mpas_pool_get_subpool(block % structs, 'state', statePool)
               call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
               call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
               call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

               call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
               call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
               call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)

               call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', normalBaroclinicVelocityCur, 1)
               call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', normalBaroclinicVelocityNew, 2)
               call mpas_pool_get_array(statePool, 'ssh', sshNew, 2)

               call mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)

               call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
               call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)

               allocate(uTemp(nVertLevels))

               ! Put f*normalBaroclinicVelocity^{perp} in uNew as a work variable
               call ocn_fuperp(statePool, meshPool, 2)

               do iEdge = 1, nEdges
                  cell1 = cellsOnEdge(1,iEdge)
                  cell2 = cellsOnEdge(2,iEdge)

                  uTemp = 0.0  ! could put this after with uTemp(maxleveledgetop+1:nvertlevels)=0
                  do k = 1, maxLevelEdgeTop(iEdge)

                     ! normalBaroclinicVelocityNew = normalBaroclinicVelocityOld + dt*(-f*normalBaroclinicVelocityPerp + T(u*,w*,p*) + g*grad(SSH*) )
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
                     normalBaroclinicVelocityNew(k,iEdge) = 0.5*( &
                       normalBaroclinicVelocityCur(k,iEdge) + uTemp(k) - dt * barotropicForcing(iEdge))

                  enddo
 
               enddo ! iEdge

               deallocate(uTemp)

               block => block % next
            end do

            call mpas_timer_start("se halo normalBaroclinicVelocity", .false., timer_halo_normalBaroclinicVelocity)
            call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
            call mpas_pool_get_field(statePool, 'normalBaroclinicVelocity', normalBaroclinicVelocityField, 2)

            call mpas_dmpar_exch_halo_field(normalBaroclinicVelocityField)
            call mpas_timer_stop("se halo normalBaroclinicVelocity", timer_halo_normalBaroclinicVelocity)

         end do  ! do j=1,config_n_bcl_iter

         call mpas_timer_stop("se bcl vel", timer_bcl_vel)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! END baroclinic iterations on linear Coriolis term
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         !  Stage 2: Barotropic velocity (2D) prediction, explicitly subcycled
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         call mpas_timer_start("se btr vel", .false., timer_btr_vel)

         oldBtrSubcycleTime = 1
         newBtrSubcycleTime = 2

         if (trim(config_time_integrator) == 'unsplit_explicit') then

            block => domain % blocklist
            do while (associated(block))
               call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
               call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)

               call mpas_pool_get_subpool(block % structs, 'state', statePool)
               call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
               call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

               call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityNew, 2)
               call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
               call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', normalBaroclinicVelocityNew, 2)

               call mpas_pool_get_array(diagnosticsPool, 'uTransport', uTransport)
               call mpas_pool_get_array(diagnosticsPool, 'uBolusGM', uBolusGM)

               call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)

               ! For Split_Explicit unsplit, simply set normalBarotropicVelocityNew=0, normalBarotropicVelocitySubcycle=0, and uNew=normalBaroclinicVelocityNew
               normalBarotropicVelocityNew(:) = 0.0

               normalVelocityNew(:,:)  = normalBaroclinicVelocityNew(:,:) 

               do iEdge = 1, nEdges
                  do k = 1, nVertLevels

                     ! uTranport = normalBaroclinicVelocity + uBolus 
                     ! This is u used in advective terms for layerThickness and tracers 
                     ! in tendency calls in stage 3.
                     uTransport(k,iEdge) = edgeMask(k,iEdge) &
                           *( normalBaroclinicVelocityNew(k,iEdge) + uBolusGM(k,iEdge) )

                  enddo
               end do  ! iEdge
   
               block => block % next
            end do  ! block

         elseif (trim(config_time_integrator) == 'split_explicit') then

            ! Initialize variables for barotropic subcycling
            block => domain % blocklist
            do while (associated(block))
               call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
               call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)

               call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
               call mpas_pool_get_subpool(block % structs, 'state', statePool)

               call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)
               call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux', barotropicThicknessFlux)

               call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
               call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, oldBtrSubcycleTime)
               call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleCur, oldBtrSubcycleTime)
               call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityCur, 1)
               call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityNew, 2)

               if (config_filter_btr_mode) then
                  barotropicForcing(:) = 0.0
               endif

               do iCell = 1, nCells
                  ! sshSubcycleOld = sshOld  
                  sshSubcycleCur(iCell) = sshCur(iCell)  
               end do

               do iEdge = 1, nEdges

                  ! normalBarotropicVelocitySubcycleOld = normalBarotropicVelocityOld 
                  normalBarotropicVelocitySubcycleCur(iEdge) = normalBarotropicVelocityCur(iEdge) 

                  ! normalBarotropicVelocityNew = BtrOld  This is the first for the summation
                  normalBarotropicVelocityNew(iEdge) = normalBarotropicVelocityCur(iEdge) 

                  ! barotropicThicknessFlux = 0  
                  barotropicThicknessFlux(iEdge) = 0.0
               end do

               block => block % next
            end do  ! block

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! BEGIN Barotropic subcycle loop
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do j = 1, config_n_btr_subcycles * config_btr_subcycle_loop_factor

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! Barotropic subcycle: VELOCITY PREDICTOR STEP
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               if (config_btr_gam1_velWt1 > 1.0e-12) then  ! only do this part if it is needed in next SSH solve
                  uPerpTime = oldBtrSubcycleTime

                  block => domain % blocklist
                  do while (associated(block))
                     call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)

                     call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
                     call mpas_pool_get_subpool(block % structs, 'state', statePool)
                     call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

                     call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
                     call mpas_pool_get_array(meshPool, 'nEdgesOnEdge', nEdgesOnEdge)
                     call mpas_pool_get_array(meshPool, 'edgesOnEdge', edgesOnEdge)
                     call mpas_pool_get_array(meshPool, 'weightsOnEdge', weightsOnEdge)
                     call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
                     call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
                     call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)

                     call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleCur, uPerpTime)
                     call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleNew, newBtrSubcycleTime)
                     call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, oldBtrSubcycleTime)

                     call mpas_pool_get_array(diagnosticsPool, 'barotropicForcing', barotropicForcing)

                     do iEdge = 1, nEdges

                        cell1 = cellsOnEdge(1,iEdge)
                        cell2 = cellsOnEdge(2,iEdge)

                        ! Compute the barotropic Coriolis term, -f*uPerp
                        CoriolisTerm = 0.0
                        do i = 1, nEdgesOnEdge(iEdge)
                           eoe = edgesOnEdge(i,iEdge)
                           CoriolisTerm = CoriolisTerm + weightsOnEdge(i,iEdge) &
                                        * normalBarotropicVelocitySubcycleCur(eoe) * fEdge(eoe)
                        end do
      
                        ! normalBarotropicVelocityNew = normalBarotropicVelocityOld + dt/J*(-f*normalBarotropicVelocityoldPerp - g*grad(SSH) + G)
                        normalBarotropicVelocitySubcycleNew(iEdge) &
                          = (normalBarotropicVelocitySubcycleCur(iEdge) &
                          + dt / config_n_btr_subcycles * (CoriolisTerm - gravity &
                          * (sshSubcycleCur(cell2) - sshSubcycleCur(cell1) ) &
                          / dcEdge(iEdge) + barotropicForcing(iEdge))) * edgeMask(1, iEdge)
                     end do

                     block => block % next
                  end do  ! block

                !   boundary update on normalBarotropicVelocityNew
                call mpas_timer_start("se halo normalBarotropicVelocity", .false., timer_halo_normalBarotropicVelocity)
                call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)

                call mpas_pool_get_field(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleField, newBtrSubcycleTime)
                call mpas_dmpar_exch_halo_field(normalBarotropicVelocitySubcycleField)
                call mpas_timer_stop("se halo normalBarotropicVelocity", timer_halo_normalBarotropicVelocity)
              endif ! config_btr_gam1_velWt1>1.0e-12

              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ! Barotropic subcycle: SSH PREDICTOR STEP 
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              block => domain % blocklist
              do while (associated(block))
                call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
                call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)

                call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
                call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
                call mpas_pool_get_subpool(block % structs, 'state', statePool)
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
                call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleCur, oldBtrSubcycleTime)
                call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleNew, newBtrSubcycleTime)

                call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux', barotropicThicknessFlux)
      
                sshTend(:) = 0.0
      
                if (config_btr_solve_SSH2) then
                   ! If config_btr_solve_SSH2=.true., then do NOT accumulate barotropicThicknessFlux in this SSH predictor 
                   ! section, because it will be accumulated in the SSH corrector section.
                   barotropicThicknessFlux_coeff = 0.0
                else
                   ! otherwise, DO accumulate barotropicThicknessFlux in this SSH predictor section
                   barotropicThicknessFlux_coeff = 1.0
                endif
      
                ! config_btr_gam1_velWt1 sets the forward weighting of velocity in the SSH computation
                ! config_btr_gam1_velWt1=  1     flux = normalBarotropicVelocityNew*H
                ! config_btr_gam1_velWt1=0.5     flux = 1/2*(normalBarotropicVelocityNew+normalBarotropicVelocityOld)*H
                ! config_btr_gam1_velWt1=  0     flux = normalBarotropicVelocityOld*H

                do iCell = 1, nCells
                  do i = 1, nEdgesOnCell(iCell)
                    iEdge = edgesOnCell(i, iCell)

                    cell1 = cellsOnEdge(1, iEdge)
                    cell2 = cellsOnEdge(2, iEdge)

                    sshEdge = 0.5 * (sshSubcycleCur(cell1) + sshSubcycleCur(cell2) )

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
                end do

                do iEdge = 1, nEdges
                   cell1 = cellsOnEdge(1,iEdge)
                   cell2 = cellsOnEdge(2,iEdge)

                   sshEdge = 0.5 * (sshSubcycleCur(cell1) &
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
      
                ! SSHnew = SSHold + dt/J*(-div(Flux))
                do iCell = 1, nCells 
                   sshSubcycleNew(iCell) = sshSubcycleCur(iCell) + dt / config_n_btr_subcycles * sshTend(iCell) / areaCell(iCell)
                end do
      
                block => block % next
              end do  ! block
      
              !   boundary update on SSHnew
              call mpas_timer_start("se halo ssh", .false., timer_halo_ssh)
              call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)

              call mpas_pool_get_field(statePool, 'sshSubcycle', sshSubcycleField, newBtrSubcycleTime)
              call mpas_dmpar_exch_halo_field(sshSubcycleField)
              call mpas_timer_stop("se halo ssh", timer_halo_ssh)
      
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ! Barotropic subcycle: VELOCITY CORRECTOR STEP
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              do BtrCorIter = 1, config_n_btr_cor_iter
                uPerpTime = newBtrSubcycleTime
      
                block => domain % blocklist
                do while (associated(block))
                   call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)

                   call mpas_pool_get_subpool(block % structs, 'state', statePool)
                   call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
                   call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

                   call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleCur, oldBtrSubcycleTime)
                   call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleNew, newBtrSubcycleTime)
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

                   allocate(utemp(nEdges+1))

                   uTemp(:) = normalBarotropicVelocitySubcycleNew(:)
                   do iEdge = 1, nEdges 
                     cell1 = cellsOnEdge(1,iEdge)
                     cell2 = cellsOnEdge(2,iEdge)
      
                     ! Compute the barotropic Coriolis term, -f*uPerp
                     CoriolisTerm = 0.0
                     do i = 1, nEdgesOnEdge(iEdge)
                        eoe = edgesOnEdge(i,iEdge)
                        CoriolisTerm = CoriolisTerm + weightsOnEdge(i,iEdge) &
                             !* normalBarotropicVelocitySubcycleNew(eoe) &
                             * uTemp(eoe) * fEdge(eoe) 
                     end do
      
                     ! In this final solve for velocity, SSH is a linear
                     ! combination of SSHold and SSHnew.
                     sshCell1 = (1-config_btr_gam2_SSHWt1) * sshSubcycleCur(cell1) + config_btr_gam2_SSHWt1 * sshSubcycleNew(cell1)
                     sshCell2 = (1-config_btr_gam2_SSHWt1) * sshSubcycleCur(cell2) + config_btr_gam2_SSHWt1 * sshSubcycleNew(cell2)
    
                     ! normalBarotropicVelocityNew = normalBarotropicVelocityOld + dt/J*(-f*normalBarotropicVelocityoldPerp - g*grad(SSH) + G)
                     normalBarotropicVelocitySubcycleNew(iEdge) = (normalBarotropicVelocitySubcycleCur(iEdge) & 
                         + dt / config_n_btr_subcycles *(CoriolisTerm - gravity *(sshCell2 - sshCell1) / dcEdge(iEdge) &
                         + barotropicForcing(iEdge))) * edgeMask(1,iEdge)
                   end do
                   deallocate(uTemp)
      
                   block => block % next
                end do  ! block
      
                !   boundary update on normalBarotropicVelocityNew
                call mpas_timer_start("se halo normalBarotropicVelocity", .false., timer_halo_normalBarotropicVelocity)
                call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)

                call mpas_pool_get_field(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleField, newBtrSubcycleTime)

                call mpas_dmpar_exch_halo_field(normalBarotropicVelocitySubcycleField)
                call mpas_timer_stop("se halo normalBarotropicVelocity", timer_halo_normalBarotropicVelocity)
              end do !do BtrCorIter=1,config_n_btr_cor_iter
      
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ! Barotropic subcycle: SSH CORRECTOR STEP
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              if (config_btr_solve_SSH2) then
      
                block => domain % blocklist
                do while (associated(block))
                   call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
                   call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)

                   call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
                   call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
                   call mpas_pool_get_subpool(block % structs, 'state', statePool)
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
                   call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleCur, oldBtrSubcycleTime)
                   call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleNew, newBtrSubcycleTime)

                   call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux', barotropicThicknessFlux)
                
                   sshTend(:) = 0.0
      
                   ! config_btr_gam3_velWt2 sets the forward weighting of velocity in the SSH computation
                   ! config_btr_gam3_velWt2=  1     flux = normalBarotropicVelocityNew*H
                   ! config_btr_gam3_velWt2=0.5     flux = 1/2*(normalBarotropicVelocityNew+normalBarotropicVelocityOld)*H
                   ! config_btr_gam3_velWt2=  0     flux = normalBarotropicVelocityOld*H

                   do iCell = 1, nCells
                     do i = 1, nEdgesOnCell(iCell)
                       iEdge = edgesOnCell(i, iCell)

                       cell1 = cellsOnEdge(1,iEdge)
                       cell2 = cellsOnEdge(2,iEdge)

                       ! SSH is a linear combination of SSHold and SSHnew.
                       sshCell1 = (1-config_btr_gam2_SSHWt1)* sshSubcycleCur(cell1) &
                                 +   config_btr_gam2_SSHWt1 * sshSubcycleNew(cell1)
                       sshCell2 = (1-config_btr_gam2_SSHWt1)* sshSubcycleCur(cell2) &
                                 +   config_btr_gam2_SSHWt1 * sshSubcycleNew(cell2)
 
                       sshEdge = 0.5 * (sshCell1 + sshCell2)

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
                   end do

                   do iEdge = 1, nEdges
                      cell1 = cellsOnEdge(1,iEdge)
                      cell2 = cellsOnEdge(2,iEdge)
      
                      ! SSH is a linear combination of SSHold and SSHnew.
                      sshCell1 = (1-config_btr_gam2_SSHWt1)* sshSubcycleCur(cell1) + config_btr_gam2_SSHWt1 * sshSubcycleNew(cell1)
                      sshCell2 = (1-config_btr_gam2_SSHWt1)* sshSubcycleCur(cell2) + config_btr_gam2_SSHWt1 * sshSubcycleNew(cell2)
                      sshEdge = 0.5 * (sshCell1 + sshCell2)

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
      
                   ! SSHnew = SSHold + dt/J*(-div(Flux))
                   do iCell = 1, nCells 
                      sshSubcycleNew(iCell) = sshSubcycleCur(iCell) & 
                           + dt / config_n_btr_subcycles * sshTend(iCell) / areaCell(iCell)
                   end do
      
                   block => block % next
                end do  ! block
      
                !   boundary update on SSHnew
                call mpas_timer_start("se halo ssh", .false., timer_halo_ssh)
                call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)

                call mpas_pool_get_field(statePool, 'sshSubcycle', sshSubcycleField)

                call mpas_dmpar_exch_halo_field(sshSubcycleField)
                call mpas_timer_stop("se halo ssh", timer_halo_ssh)
               endif ! config_btr_solve_SSH2
      
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! Barotropic subcycle: Accumulate running sums, advance timestep pointers
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
               block => domain % blocklist
               do while (associated(block))
                  call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)

                  call mpas_pool_get_subpool(block % structs, 'state', statePool)

                  call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityNew, 2)
                  call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', normalBarotropicVelocitySubcycleNew, newBtrSubcycleTime)
      
                  ! normalBarotropicVelocityNew = normalBarotropicVelocityNew + normalBarotropicVelocitySubcycleNEW
                  ! This accumulates the sum.
                  ! If the Barotropic Coriolis iteration is limited to one, this could 
                  ! be merged with the above code.
                  do iEdge = 1, nEdges 
                       normalBarotropicVelocityNew(iEdge) = normalBarotropicVelocityNew(iEdge) + normalBarotropicVelocitySubcycleNew(iEdge)  
                  end do  ! iEdge
                  block => block % next
               end do  ! block
      
               ! advance time pointers
               oldBtrSubcycleTime = mod(oldBtrSubcycleTime,2)+1
               newBtrSubcycleTime = mod(newBtrSubcycleTime,2)+1
      
            end do ! j=1,config_n_btr_subcycles
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! END Barotropic subcycle loop
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! Normalize Barotropic subcycle sums: ssh, normalBarotropicVelocity, and F
            block => domain % blocklist
            do while (associated(block))
               call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)

               call mpas_pool_get_subpool(block % structs, 'state', statePool)
               call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

               call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityNew, 2)

               call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux', barotropicThicknessFlux)
      
               do iEdge = 1, nEdges
                  barotropicThicknessFlux(iEdge) = barotropicThicknessFlux(iEdge) &
                      / (config_n_btr_subcycles * config_btr_subcycle_loop_factor)
      
                  normalBarotropicVelocityNew(iEdge) = normalBarotropicVelocityNew(iEdge) & 
                     / (config_n_btr_subcycles * config_btr_subcycle_loop_factor + 1)
               end do
      
               block => block % next
            end do  ! block
      
      
            ! boundary update on F
            call mpas_timer_start("se halo F", .false., timer_halo_f)
            call mpas_pool_get_subpool(domain % blocklist % structs, 'diagnostics', diagnosticsPool)

            call mpas_pool_get_field(diagnosticsPool, 'barotropicThicknessFlux', barotropicThicknessFluxField)

            call mpas_dmpar_exch_halo_field(barotropicThicknessFluxField)
            call mpas_timer_stop("se halo F", timer_halo_f)


            ! Check that you can compute SSH using the total sum or the individual increments
            ! over the barotropic subcycles.
            ! efficiency: This next block of code is really a check for debugging, and can 
            ! be removed later.
            block => domain % blocklist
            do while (associated(block))
               call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
               call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)

               call mpas_pool_get_subpool(block % structs, 'state', statePool)
               call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
               call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

               call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocityNew, 2)
               call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', normalBaroclinicVelocityNew, 2)

               call mpas_pool_get_array(diagnosticsPool, 'uTransport', uTransport)
               call mpas_pool_get_array(diagnosticsPool, 'uBolusGM', uBolusGM)
               call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
               call mpas_pool_get_array(diagnosticsPool, 'barotropicThicknessFlux', barotropicThicknessFlux)

               call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
               call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)

               allocate(uTemp(nVertLevels))

               ! Correction velocity    uCorr = (Flux - Sum(h u*))/H
               ! or, for the full latex version:
               !{\bf u}^{corr} = \left( {\overline {\bf F}} 
               !  - \sum_{k=1}^{N^{edge}} h_{k,*}^{edge}  {\bf u}_k^{avg} \right)
               ! \left/ \sum_{k=1}^{N^{edge}} h_{k,*}^{edge}   \right. 

               if (config_vel_correction) then
                  ucorr_coef = 1
               else
                  ucorr_coef = 0
               endif

               do iEdge = 1, nEdges

                  ! velocity for uCorrection is normalBarotropicVelocity + normalBaroclinicVelocity + uBolus
                  uTemp(:) = normalBarotropicVelocityNew(iEdge) + normalBaroclinicVelocityNew(:,iEdge) + uBolusGM(:,iEdge)

                  ! thicknessSum is initialized outside the loop because on land boundaries 
                  ! maxLevelEdgeTop=0, but I want to initialize thicknessSum with a 
                  ! nonzero value to avoid a NaN.
                  normalThicknessFluxSum = layerThicknessEdge(1,iEdge) * uTemp(1)
                  thicknessSum  = layerThicknessEdge(1,iEdge)

                  do k = 2, maxLevelEdgeTop(iEdge)
                     normalThicknessFluxSum = normalThicknessFluxSum + layerThicknessEdge(k,iEdge) * uTemp(k)
                     thicknessSum  =  thicknessSum + layerThicknessEdge(k,iEdge)
                  enddo

                  uCorr =   ucorr_coef*(( barotropicThicknessFlux(iEdge) - normalThicknessFluxSum)/thicknessSum)

                  do k = 1, nVertLevels

                     ! uTranport = normalBarotropicVelocity + normalBaroclinicVelocity + uBolus + uCorrection
                     ! This is u used in advective terms for layerThickness and tracers 
                     ! in tendency calls in stage 3.
                     uTransport(k,iEdge) &
                           = edgeMask(k,iEdge) &
                           *( normalBarotropicVelocityNew(iEdge) + normalBaroclinicVelocityNew(k,iEdge) &
                           + uBolusGM(k,iEdge) + uCorr )
                  enddo

               end do ! iEdge

               deallocate(uTemp)

               block => block % next
            end do  ! block

         endif ! split_explicit  

         call mpas_timer_stop("se btr vel", timer_btr_vel)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         !  Stage 3: Tracer, density, pressure, vertical velocity prediction
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! Thickness tendency computations and thickness halo updates are completed before tracer 
         ! tendency computations to allow monotonic advection.
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)

            call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
            call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
            call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThicknessNew, 2)

            call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
            call mpas_pool_get_array(diagnosticsPool, 'uTransport', uTransport)
            call mpas_pool_get_array(diagnosticsPool, 'vertTransportVelocityTop', vertTransportVelocityTop)

            ! compute vertTransportVelocityTop.  Use uTransport for advection of layerThickness and tracers.
            ! Use time level 1 values of layerThickness and layerThicknessEdge because 
            ! layerThickness has not yet been computed for time level 2.
            call ocn_vert_transport_velocity_top(meshPool, verticalMeshPool, &
               layerThicknessCur, layerThicknessEdge, uTransport, &
               sshCur, highFreqThicknessNew, dt, vertTransportVelocityTop, err)

            call ocn_tend_thick(tendPool, forcingPool, diagnosticsPool, meshPool)

            block => block % next
         end do

         ! update halo for thickness tendencies
         call mpas_timer_start("se halo thickness", .false., timer_halo_thickness)
         call mpas_pool_get_subpool(domain % blocklist % structs, 'tend', tendPool)

         call mpas_pool_get_field(tendPool, 'layerThickness', layerThicknessField)

         call mpas_dmpar_exch_halo_field(layerThicknessField)
         call mpas_timer_stop("se halo thickness", timer_halo_thickness)

         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call ocn_tend_tracer(tendPool, statePool, forcingPool, diagnosticsPool, meshPool, dt, 2)

            block => block % next
         end do

         ! update halo for tracer tendencies
         call mpas_timer_start("se halo tracers", .false., timer_halo_tracers)
         call mpas_pool_get_subpool(domain % blocklist % structs, 'tend', tendPool)

         call mpas_pool_get_field(tendPool, 'tracers', tracersField)

         call mpas_dmpar_exch_halo_field(tracersField)
         call mpas_timer_stop("se halo tracers", timer_halo_tracers)

         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
            call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
            call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)

            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
            call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
            call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
            call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

            call mpas_pool_get_dimension(statePool, 'num_tracers', num_tracers)

            call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
            call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)

            call mpas_pool_get_array(statePool, 'tracers', tracersCur, 1)
            call mpas_pool_get_array(statePool, 'tracers', tracersNew, 2)
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

            call mpas_pool_get_array(tendPool, 'tracers', tracersTend)
            call mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)
            call mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)
            call mpas_pool_get_array(tendPool, 'highFreqThickness', highFreqThicknessTend)
            call mpas_pool_get_array(tendPool, 'lowFreqDivergence', lowFreqDivergenceTend)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            !  If iterating, reset variables for next iteration
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (split_explicit_step < config_n_ts_iter) then

               ! Get indices for dynamic tracers (Includes T&S).
               call mpas_pool_get_dimension(statePool, 'dynamics_start', startIndex)
               call mpas_pool_get_dimension(statePool, 'dynamics_end', endIndex)

               ! Only need T & S for earlier iterations,
               ! then all the tracers needed the last time through.
               do iCell = 1, nCells
                  ! sshNew is a pointer, defined above.
                  do k = 1, maxLevelCell(iCell)

                     ! this is h_{n+1}
                     temp_h = layerThicknessCur(k,iCell) + dt * layerThicknessTend(k,iCell) 

                     ! this is h_{n+1/2}
                     layerThicknessNew(k,iCell) = 0.5*( layerThicknessCur(k,iCell) + temp_h)

                     do i = startIndex, endIndex
                        ! This is Phi at n+1
                        temp = ( tracersCur(i,k,iCell) * layerThicknessCur(k,iCell) + dt * tracersTend(i,k,iCell)) / temp_h
  
                        ! This is Phi at n+1/2
                        tracersNew(i,k,iCell) = 0.5 * ( tracersCur(i,k,iCell) + temp )
                     end do
                  end do
               end do ! iCell

               if (config_use_freq_filtered_thickness) then
                  do iCell = 1, nCells
                     do k = 1, maxLevelCell(iCell)

                        ! h^{hf}_{n+1} was computed in Stage 1

                        ! this is h^{hf}_{n+1/2}
                        highFreqThicknessnew(k,iCell) = 0.5 * (highFreqThicknessCur(k,iCell) + highFreqThicknessNew(k,iCell))

                        ! this is D^{lf}_{n+1}
                        temp = lowFreqDivergenceCur(k,iCell) &
                         + dt * lowFreqDivergenceTend(k,iCell) 

                        ! this is D^{lf}_{n+1/2}
                        lowFreqDivergenceNew(k,iCell) = 0.5 * (lowFreqDivergenceCur(k,iCell) + temp)
                     end do
                  end do
               end if

               do iEdge = 1, nEdges

                  do k = 1, nVertLevels

                     ! u = normalBarotropicVelocity + normalBaroclinicVelocity 
                     ! here normalBaroclinicVelocity is at time n+1/2
                     ! This is u used in next iteration or step
                     normalVelocityNew(k,iEdge) = edgeMask(k,iEdge) * ( normalBarotropicVelocityNew(iEdge) + normalBaroclinicVelocityNew(k,iEdge) )

                  enddo

               end do ! iEdge

               ! Efficiency note: We really only need this to compute layerThicknessEdge, density, pressure, and SSH 
               ! in this diagnostics solve.
               call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, 2)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            !  If large iteration complete, compute all variables at time n+1
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            elseif (split_explicit_step == config_n_ts_iter) then

               do iCell = 1, nCells
                  do k = 1, maxLevelCell(iCell)

                     ! this is h_{n+1}
                     layerThicknessNew(k,iCell) = layerThicknessCur(k,iCell) + dt * layerThicknessTend(k,iCell) 

                     ! This is Phi at n+1
                     do i = 1, num_tracers
                        tracersNew(i,k,iCell) = (tracersCur(i,k,iCell) * layerThicknessCur(k,iCell) + dt * tracersTend(i,k,iCell) ) &
                         / layerThicknessNew(k,iCell)

                     enddo
                  end do
               end do

               if (config_use_freq_filtered_thickness) then
                  do iCell = 1, nCells
                     do k = 1, maxLevelCell(iCell)

                        ! h^{hf}_{n+1} was computed in Stage 1

                        ! this is D^{lf}_{n+1}
                        lowFreqDivergenceNew(k,iCell) = lowFreqDivergenceCur(k,iCell) + dt * lowFreqDivergenceTend(k,iCell) 
                     end do
                  end do
               end if

               ! Recompute final u to go on to next step.
               ! u_{n+1} = normalBarotropicVelocity_{n+1} + normalBaroclinicVelocity_{n+1} 
               ! Right now normalBaroclinicVelocityNew is at time n+1/2, so back compute to get normalBaroclinicVelocity at time n+1
               !   using normalBaroclinicVelocity_{n+1/2} = 1/2*(normalBaroclinicVelocity_n + u_Bcl_{n+1})
               ! so the following lines are
               ! u_{n+1} = normalBarotropicVelocity_{n+1} + 2*normalBaroclinicVelocity_{n+1/2} - normalBaroclinicVelocity_n
               ! note that normalBaroclinicVelocity is recomputed at the beginning of the next timestep due to Imp Vert mixing,
               ! so normalBaroclinicVelocity does not have to be recomputed here.
      
               do iEdge = 1, nEdges
                  do k = 1, maxLevelEdgeTop(iEdge)
                     normalVelocityNew(k,iEdge) = normalBarotropicVelocityNew(iEdge) + 2 * normalBaroclinicVelocityNew(k,iEdge) - normalBaroclinicVelocityCur(k,iEdge)
                  end do
               end do ! iEdges

            endif ! split_explicit_step

            block => block % next
         end do



      end do  ! split_explicit_step = 1, config_n_ts_iter
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! END large iteration loop 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Perform Sea Ice Formation Adjustment
      block => domain % blocklist
      do while(associated(block))
        call mpas_pool_get_subpool(block % structs, 'state', statePool)
        call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
        call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

        call mpas_pool_get_dimension(statePool, 'index_temperature', indexTemperature)
        call mpas_pool_get_dimension(statePool, 'index_salinity', indexSalinity)

        call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)
        call mpas_pool_get_array(statePool, 'tracers', tracersNew, 2)

        call mpas_pool_get_array(forcingPool, 'seaIceEnergy', seaIceEnergy)

        call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, 2)
        call ocn_sea_ice_formation(meshPool, indexTemperature, indexSalinity, layerThicknessNew, &
                                   tracersNew, seaIceEnergy, err)
        block => block % next
      end do

      call mpas_timer_start("se implicit vert mix")
      block => domain % blocklist
      do while(associated(block))
        call mpas_pool_get_subpool(block % structs, 'state', statePool)
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
        call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, 2)

        call ocn_vmix_implicit(dt, meshPool, diagnosticsPool, statePool, err, 2)

        block => block % next
      end do

      ! Update halo on u and tracers, which were just updated for implicit vertical mixing.  If not done, 
      ! this leads to lack of volume conservation.  It is required because halo updates in stage 3 are only
      ! conducted on tendencies, not on the velocity and tracer fields.  So this update is required to 
      ! communicate the change due to implicit vertical mixing across the boundary.
      call mpas_timer_start("se implicit vert mix halos")
      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)

      call mpas_pool_get_field(statePool, 'normalVelocity', normalVelocityField, 2)
      call mpas_pool_get_field(statePool, 'tracers', tracersField, 2)

      call mpas_dmpar_exch_halo_field(normalVelocityField)
      call mpas_dmpar_exch_halo_field(tracersField)
      call mpas_timer_stop("se implicit vert mix halos")

      call mpas_timer_stop("se implicit vert mix")

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'average', averagePool)

         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessCur, 1)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, 2)

         call mpas_pool_get_array(diagnosticsPool, 'uTransport', uTransport)
         call mpas_pool_get_array(diagnosticsPool, 'uBolusGM', uBolusGM)
         call mpas_pool_get_array(diagnosticsPool, 'normalVelocityX', normalVelocityX)
         call mpas_pool_get_array(diagnosticsPool, 'normalVelocityY', normalVelocityY)
         call mpas_pool_get_array(diagnosticsPool, 'normalVelocityZ', normalVelocityZ)
         call mpas_pool_get_array(diagnosticsPool, 'normalVelocityZonal', normalVelocityZonal)
         call mpas_pool_get_array(diagnosticsPool, 'normalVelocityMeridional', normalVelocityMeridional)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSH', gradSSH)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHX', gradSSHX)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHY', gradSSHY)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHZ', gradSSHZ)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHZonal', gradSSHZonal)
         call mpas_pool_get_array(diagnosticsPool, 'gradSSHMeridional', gradSSHMeridional)

         call mpas_pool_get_array(diagnosticsPool, 'surfaceVelocity', surfaceVelocity)
         call mpas_pool_get_array(diagnosticsPool, 'SSHGradient', SSHGradient)

         call mpas_pool_get_dimension(diagnosticsPool, 'index_zonalSurfaceVelocity', indexZonalSurfaceVelocity)
         call mpas_pool_get_dimension(diagnosticsPool, 'index_meridionalSurfaceVelocity', indexMeridionalSurfaceVelocity)
         call mpas_pool_get_dimension(diagnosticsPool, 'index_zonalSSHGradient', indexZonalSSHGradient)
         call mpas_pool_get_dimension(diagnosticsPool, 'index_meridionalSSHGradient', indexMeridionalSSHGradient)

         if (config_prescribe_velocity) then
            normalVelocityNew(:,:) = normalVelocityCur(:,:)
         end if

         if (config_prescribe_thickness) then
            layerThicknessNew(:,:) = layerThicknessCur(:,:)
         end if

         call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, 2)

         ! Compute velocity transport, used in advection terms of layerThickness and tracer tendency
         uTransport(:,:) = normalVelocityNew(:,:) + uBolusGM(:,:)

         call mpas_reconstruct(meshPool, normalVelocityNew,                  &
                          normalVelocityX, normalVelocityY, normalVelocityZ, &
                          normalVelocityZonal, normalVelocityMeridional      &
                         )

         call mpas_reconstruct(meshPool, gradSSH,         &
                          gradSSHX, gradSSHY, gradSSHZ,   &
                          gradSSHZonal, gradSSHMeridional &
                         )

         surfaceVelocity(indexZonalSurfaceVelocity, :) = normalVelocityZonal(1, :)
         surfaceVelocity(indexMeridionalSurfaceVelocity, :) = normalVelocityMeridional(1, :)

         SSHGradient(indexZonalSSHGradient, :) = gradSSHZonal(1, :)
         SSHGradient(indexMeridionalSSHGradient, :) = gradSSHMeridional(1, :)

         call ocn_time_average_accumulate(averagePool, statePool, diagnosticsPool, 2)
         call ocn_time_average_coupled_accumulate(diagnosticsPool, forcingPool)

         block => block % next
      end do

      call mpas_timer_stop("se timestep", timer_main)

      deallocate(n_bcl_iter)

   end subroutine ocn_time_integrator_split!}}}

end module ocn_time_integration_split

! vim: foldmethod=marker
