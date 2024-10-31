! Copyright (c) 2016-2018,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_SGH_driver
!
!> \brief MPAS land ice SGH primary routines
!> \author Matt Hoffman
!> \date   27 June 2016
!> \details
!>  This module contains the main driver routines for
!>  for subglacial hydro.
!
!-----------------------------------------------------------------------

module li_subglacial_hydro

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_timer
   use mpas_log

   use li_setup
   use li_constants
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

   public :: li_SGH_init, &
             li_SGH_solve, &
             li_SGH_finalize

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------
   ! Minimum gradMagPhiBaseEdge and gradMagPhiEdge allowed before all dependent variables are zeroed out
   real(kind=RKIND), parameter :: SMALL_GRADPHI = 1.0e-6_RKIND

   !Minimum outflowing hydropotential slope applied at grounding line
   real(kind=RKIND), parameter :: MIN_PHISLOPE_GL = 1e-10_RKIND

!***********************************************************************
   contains


!***********************************************************************
!
!  routine li_SGH_init
!
!> \brief   Initialize SGH
!> \author  Matt Hoffman
!> \date    27 June 2016
!> \details
!>  This routine initializes the subglacial hydro model
!-----------------------------------------------------------------------
   subroutine li_SGH_init(domain, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      ! Pools pointers
      logical, pointer :: config_SGH
      logical, pointer :: config_ocean_connection_N
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: hydroPool
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: velocityPool
      real (kind=RKIND), pointer :: deltatSGH
      real (kind=RKIND), dimension(:), pointer :: waterThickness
      real (kind=RKIND), dimension(:), pointer :: tillWaterThickness
      real (kind=RKIND), dimension(:), pointer :: waterPressure
      real (kind=RKIND), dimension(:), pointer :: thickness
      real (kind=RKIND), dimension(:), pointer :: bedTopography
      real (kind=RKIND), dimension(:), pointer :: iceThicknessHydro
      real (kind=RKIND), dimension(:), pointer :: hydropotential
      integer, dimension(:), pointer :: cellMask
      real (kind=RKIND), pointer :: tillMax
      real (kind=RKIND), pointer :: rhoi, rhoo
      logical, pointer :: config_do_restart
      real (kind=RKIND), pointer :: config_sea_level
      integer, pointer :: config_num_halos
      integer :: err_tmp


      err = 0
      err_tmp = 0

      call mpas_pool_get_config(liConfigs, 'config_SGH', config_SGH)
      if (.not. config_SGH) then
         call mpas_pool_get_config(liConfigs, 'config_ocean_connection_N', config_ocean_connection_N)
         if (config_ocean_connection_N) then
            call ocean_connection_N(domain)
         endif
         ! If SGH is not active, skip everything
         return
      endif


      call mpas_timer_start("hydro init")

      call mpas_log_write('Beginning subglacial hydro init.')

      ! Check number of halos
      call mpas_pool_get_config(liConfigs, 'config_num_halos', config_num_halos)
      if (config_num_halos < 1) then
         call mpas_log_write("Subglacial hydrology requires config_num_halos >= 1", MPAS_LOG_ERR)
         err = ior(err, 1)
      endif

      call mpas_pool_get_config(liConfigs, 'config_do_restart', config_do_restart)
      call mpas_pool_get_config(liConfigs, 'config_SGH_till_max', tillMax)
      call mpas_pool_get_config(liConfigs, 'config_ice_density', rhoi)
      call mpas_pool_get_config(liConfigs, 'config_ocean_density', rhoo)
      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)

         call mpas_pool_get_array(hydroPool, 'deltatSGH', deltatSGH)

         if (.not. config_do_restart) then
            ! On cold start, set initial timestep to a small value.
            ! On a restart we will use the value from the last subcycle of the previous
            ! master timestep to make restarts BFB.
            ! Note the value here is only used for the first update of thetill model;
            ! the sheet and channel models calculate the adaptive timestep, but that is not
            ! available for the till model as the code is currently organized.
            ! TODO: Move till update until after adaptive timestep has been set
            deltatSGH = 1.0e-4_RKIND ! in seconds
         endif

         ! Mask needs to be initialized for pressure calcs to be correct
         call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
         err = ior(err, err_tmp)

         call calc_hydro_mask(domain)

         ! remove invalid values - not necessary on restart, but shouldn't hurt
         call mpas_pool_get_array(hydroPool, 'waterThickness', waterThickness)
         waterThickness = max(0.0_RKIND, waterThickness)

         call mpas_pool_get_array(hydroPool, 'tillWaterThickness', tillWaterThickness)
         tillWaterThickness = max(0.0_RKIND, tillWaterThickness)
         tillWaterThickness = min(tillMax, tillWaterThickness)

         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(hydroPool, 'iceThicknessHydro', iceThicknessHydro) 
         call calc_iceThicknessHydro(block, err_tmp) !adjust ice thickness along boundaries
         err = ior(err,err_tmp)
         block => block % next
      end do

      !update halo for iceThicknessHydro
      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'iceThicknessHydro')
      call mpas_timer_stop("halo updates")
      
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)

         call mpas_pool_get_array(hydroPool, 'waterPressure', waterPressure)
         call mpas_pool_get_array(hydroPool, 'hydropotential', hydropotential)
         call mpas_pool_get_array(hydroPool, 'iceThicknessHydro', iceThicknessHydro) 
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
         
         waterPressure = max(0.0_RKIND, waterPressure)
         where (li_mask_is_grounded_ice(cellMask))
            waterPressure = min(waterPressure, rhoi * gravity * iceThicknessHydro)
         end where
         
         ! set pressure and hydropotential correctly on ice-free land and in ocean
         where ((.not. (li_mask_is_grounded_ice(cellMask))) .and. (bedTopography > config_sea_level))
            waterPressure = 0.0_RKIND
            hydropotential = rho_water * gravity * bedTopography
         elsewhere ((.not. (li_mask_is_grounded_ice(cellMask))) .and. (bedTopography <= config_sea_level))
            waterPressure = gravity * rhoo * (config_sea_level - bedTopography)
            hydropotential = 0.0_RKIND
         end where
        
         ! Initialize diagnostic pressure variables
         call calc_pressure_diag_vars(block, err_tmp)
         err = ior(err, err_tmp)

         !smooth water pressure for calculation of waterPressureSlopeNormal
         call calc_waterPressureSmooth(block, err_tmp) 
         err = ior(err,err_tmp)

         block => block % next
      end do

      !updates halos for waterPressureSmooth
      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'waterPressureSmooth')
      call mpas_timer_stop("halo updates")
     
      ! === error check
      if (err > 0) then
          call mpas_log_write("An error has occurred in li_SGH_init.", MPAS_LOG_ERR)
      endif

      call mpas_timer_stop("hydro init")

   !--------------------------------------------------------------------
   end subroutine li_SGH_init



!***********************************************************************
!
!  routine li_SGH_solve
!
!> \brief   Solve and update SGH for current time step
!> \author  Matt Hoffman
!> \date    27 June 2016
!> \details
!>  This routine solves and updates the subglacial hydro model
!>  for the current ISM time step.
!-----------------------------------------------------------------------
   subroutine li_SGH_solve(domain, err)
      use mpas_vector_reconstruction
      use li_diagnostic_vars

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      ! Pools pointers
      logical, pointer :: config_SGH
      logical, pointer :: config_ocean_connection_N
      logical, pointer :: config_SGH_chnl_active
      character (len=StrKIND), pointer :: config_SGH_basal_melt
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: hydroPool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: thermalPool
      type (mpas_pool_type), pointer :: velocityPool
      real (kind=RKIND), dimension(:), pointer :: thickness
      real (kind=RKIND), dimension(:,:), pointer :: temperature, flowParamA
      real (kind=RKIND), dimension(:), pointer :: Wtill, WtillOld
      real (kind=RKIND), dimension(:), pointer :: basalMeltInput
      real (kind=RKIND), dimension(:), pointer :: groundedBasalMassBal
      real (kind=RKIND), dimension(:), pointer :: basalHeatFlux
      real (kind=RKIND), dimension(:), pointer :: basalFrictionFlux
      real (kind=RKIND), dimension(:), pointer :: externalWaterInput
      real (kind=RKIND), dimension(:), pointer :: divergence
      real (kind=RKIND), dimension(:), pointer :: waterFlux
      real (kind=RKIND), dimension(:), pointer :: waterThickness
      real (kind=RKIND), dimension(:), pointer :: waterThicknessOld
      real (kind=RKIND), dimension(:), pointer :: waterThicknessTendency
      real (kind=RKIND), dimension(:), pointer :: divergenceChannel
      real (kind=RKIND), dimension(:), pointer :: channelAreaChangeCell
      real (kind=RKIND), dimension(:), pointer :: channelMeltInputCell
      real (kind=RKIND), dimension(:), pointer :: dvEdge
      real (kind=RKIND), dimension(:), pointer :: areaCell
      real (kind=RKIND), dimension(:), pointer :: waterVelocity
      real (kind=RKIND), dimension(:), pointer :: waterVelocityCellX
      real (kind=RKIND), dimension(:), pointer :: waterVelocityCellY
      real (kind=RKIND), dimension(:), pointer :: effectivePressureSGH
      real (kind=RKIND), dimension(:), pointer :: effectivePressure
      real (kind=RKIND), dimension(:), allocatable :: cellJunk
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer, dimension(:,:), pointer :: edgeSignOnCell
      integer, dimension(:), pointer :: cellMask
      real (kind=RKIND), pointer :: deltatSGH
      real (kind=RKIND), pointer :: masterDeltat
      real (kind=RKIND), pointer :: Cd
      real (kind=RKIND), pointer :: tillMax
      integer, pointer :: nCellsSolve
      integer, pointer :: nCells
      integer :: iCell, iEdge, iEdgeOnCell
      real (kind=RKIND) :: timeLeft ! subcycling time remaining until master dt is reached
      real (kind=RKIND) :: deltatSGHaccumulated
      integer :: numSubCycles ! number of subcycles
      integer :: err_tmp

      err = 0
      err_tmp = 0

      call mpas_pool_get_config(liConfigs, 'config_SGH', config_SGH)
      if (.not. config_SGH) then
         call mpas_pool_get_config(liConfigs, 'config_ocean_connection_N', config_ocean_connection_N)
         if (config_ocean_connection_N) then
            call ocean_connection_N(domain)
         endif
         ! If SGH is not active, skip everything
         return
      endif

      call calc_hydro_mask(domain)

      call mpas_log_write('Beginning subglacial hydro solve.')
      call mpas_pool_get_config(liConfigs, 'config_SGH_chnl_active', config_SGH_chnl_active)
      call mpas_pool_get_config(liConfigs, 'config_SGH_till_drainage', Cd)
      call mpas_pool_get_config(liConfigs, 'config_SGH_till_max', tillMax)
      call mpas_pool_get_config(liConfigs, 'config_SGH_basal_melt', config_SGH_basal_melt)

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'thermal', thermalPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_array(thermalPool, 'temperature', temperature)
         call mpas_pool_get_array(velocityPool, 'flowParamA', flowParamA)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)

         call calc_iceThicknessHydro(block, err_tmp)
         err = ior(err, err_tmp)

         call li_calculate_flowParamA(meshPool, temperature, thickness, flowParamA, err_tmp)
         err = ior(err, err_tmp)

         block => block % next
      end do

      !update halo for iceThicknessHydro
      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'iceThicknessHydro')
      call mpas_timer_stop("halo updates")
      
      ! initialize while loop
      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)  ! can get from any block
      call mpas_pool_get_array(meshPool, 'deltat', masterDeltat)
      timeLeft = masterDeltat  ! in seconds
      numSubCycles = 0

      block => domain % blocklist
      do while (associated(block))
        
         call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
         call mpas_pool_get_array(hydroPool, 'effectivePressure', effectivePressure)
         effectivePressure = 0.0_RKIND
         deltatSGHaccumulated = 0.0_RKIND

         block => block % next
      end do


      ! =============
      ! =============
      ! =============
      ! subcycle hydro model until master dt is reached
      ! =============
      ! =============
      ! =============
      timecycling: do while (timeLeft > 0.0_RKIND)
      numSubCycles = numSubCycles + 1


      ! =============
      ! Calculate time-varying forcing, as needed
      ! =============
      block => domain % blocklist
      do while (associated(block))

         ! Decide where the basal melt input comes from
         ! Either prescribed, from the temperature model, or calculated here
         if (trim(config_SGH_basal_melt) == 'file') then
             ! do nothing
         elseif (trim(config_SGH_basal_melt) == 'thermal') then
             call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
             call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
             call mpas_pool_get_array(hydroPool, 'basalMeltInput', basalMeltInput)
             call mpas_pool_get_array(geometryPool,'groundedBasalMassBal',groundedBasalMassBal)
             basalMeltInput = -1.0_RKIND * groundedBasalMassBal ! TODO: Ensure positive flux?
         elseif (trim(config_SGH_basal_melt) == 'basal_heat') then
             call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
             call mpas_pool_get_subpool(block % structs, 'thermal', thermalPool)
             call mpas_pool_get_array(hydroPool, 'basalMeltInput', basalMeltInput)
             call mpas_pool_get_array(thermalPool, 'basalHeatFlux', basalHeatFlux)
             call mpas_pool_get_array(thermalPool, 'basalFrictionFlux', basalFrictionFlux)
             basalMeltInput = (basalHeatFlux + basalFrictionFlux) / latent_heat_ice
         else
            call mpas_log_write("Unknown value provided for config_SGH_basal_melt: " // trim(config_SGH_basal_melt), MPAS_LOG_ERR)
            err = ior(err, 1)
         endif

         ! SHMIP forcing will override basalMeltInput.
         call shmip_timevarying_forcing(block, err_tmp)
         err = ior(err, err_tmp)

         block => block % next
      end do

      ! Perform halo update, if needed
      if ((trim(config_SGH_basal_melt) == 'thermal') .or. (trim(config_SGH_basal_melt) == 'basal_heat')) then
         call mpas_timer_start("halo updates")
         ! (groundedBasalMassBal will not have had a halo update,
         !  because a halo update for it is unneeded except for the SGH model.)
         call mpas_dmpar_field_halo_exch(domain, 'basalMeltInput')
         call mpas_timer_stop("halo updates")
      endif


      ! =============
      ! Update till water layer thickness
      ! =============
      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)

         call mpas_pool_get_array(hydroPool, 'tillWaterThickness', Wtill)
         call mpas_pool_get_array(hydroPool, 'tillWaterThicknessOld', WtillOld)
         call mpas_pool_get_array(hydroPool, 'deltatSGH', deltatSGH)
         call mpas_pool_get_array(hydroPool, 'basalMeltInput', basalMeltInput)
         call mpas_pool_get_array(hydroPool, 'externalWaterInput', externalWaterInput)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)

         WtillOld = Wtill
         Wtill = Wtill + deltatSGH * ( (basalMeltInput + externalWaterInput) / rho_water - Cd)
         Wtill = Wtill * li_mask_is_grounded_ice_int(cellMask)  ! zero Wtill in non-grounded locations
         Wtill = min(Wtill, tillmax)
         Wtill = max(0.0_RKIND, Wtill)

         block => block % next
      end do


      ! =============
      ! Calculate edge quantities and advective fluxes
      ! =============
      block => domain % blocklist
      do while (associated(block))

         call calc_edge_quantities(block, err_tmp)
         err = ior(err, err_tmp)

         block => block % next
      end do
      ! Update halos on edge quantities
      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'waterFlux')
      ! intermediate fields will be out of date, but will be correct in output files
      ! waterVelocity needs to be updated in order to get waterVelocityCellX/Y correct
      call mpas_dmpar_field_halo_exch(domain, 'waterVelocity')
      call mpas_timer_stop("halo updates")

      ! Calculate reconstructed velocities on cell centers for viz
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_array(hydroPool, 'waterVelocity', waterVelocity)
         call mpas_pool_get_array(hydroPool, 'waterVelocityCellX', waterVelocityCellX)
         call mpas_pool_get_array(hydroPool, 'waterVelocityCellY', waterVelocityCellY)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         allocate(cellJunk(nCells+1))
         call mpas_reconstruct(meshPool, waterVelocity, waterVelocityCellX, waterVelocityCellY, &
              cellJunk, cellJunk, cellJunk)
         deallocate(cellJunk)
         waterVelocityCellX = waterVelocityCellX * li_mask_is_grounded_ice_int(cellMask)  ! zero in non-grounded locations
         waterVelocityCellY = waterVelocityCellY * li_mask_is_grounded_ice_int(cellMask)  ! zero in non-grounded locations

         block => block % next
      end do


      ! =============
      ! Update Channel fields
      ! =============
      if (config_SGH_chnl_active) then
         block => domain % blocklist
         do while (associated(block))

            call update_channel(block, err_tmp)
            err = ior(err, err_tmp)

            block => block % next
         end do
         ! Update halos on channel
         call mpas_timer_start("halo updates")
         call mpas_dmpar_field_halo_exch(domain, 'channelChangeRate')
         call mpas_dmpar_field_halo_exch(domain, 'channelDischarge')
         call mpas_dmpar_field_halo_exch(domain, 'channelMelt')
         call mpas_dmpar_field_halo_exch(domain, 'channelPressureFreeze')
         call mpas_timer_stop("halo updates")
      endif


      ! =============
      ! Calculate adaptive time step
      ! =============
      call check_timestep(domain, timeLeft, numSubCycles, err_tmp)
      err = ior(err, err_tmp)
      if (err > 0) then
         exit timecycling
      endif


      ! =============
      ! Compute flux divergence
      ! =============
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_array(hydroPool, 'divergence', divergence)
         call mpas_pool_get_array(hydroPool, 'waterFlux', waterFlux)
         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
         call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
         call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
         call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)

         divergence(:) = 0.0_RKIND  ! zero div before starting
         ! loop over locally owned cells
         do iCell = 1, nCellsSolve
            if (li_mask_is_grounded_ice(cellMask(iCell))) then  ! can skip for non-grounded ice
               ! compute fluxes for each edge of the cell
               do iEdgeOnCell = 1, nEdgesOnCell(iCell)
                  iEdge = edgesOnCell(iEdgeOnCell, iCell)
                  ! add on advective & diffusive fluxes
                  divergence(iCell) = divergence(iCell) - waterFlux(iEdge) * dvEdge(iEdge) * edgeSignOnCell(iEdgeOnCell, iCell)
               end do ! edges
            end if
         end do ! cells
         divergence(1:nCellsSolve) = divergence(1:nCellsSolve) / areaCell(1:nCellsSolve)

         block => block % next
      end do
      ! Update halos on divergence
      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'divergence')
      call mpas_timer_stop("halo updates")


      ! =============
      ! Update channel area now that we have time step
      ! =============
      if (config_SGH_chnl_active) then
         block => domain % blocklist
         do while (associated(block))

            call evolve_channel(block, err_tmp)
            err = ior(err, err_tmp)

            block => block % next
         end do
         call mpas_timer_start("halo updates")
         call mpas_dmpar_field_halo_exch(domain, 'divergenceChannel')
         call mpas_dmpar_field_halo_exch(domain, 'channelAreaChangeCell')
         call mpas_dmpar_field_halo_exch(domain, 'channelMeltInputCell')
         call mpas_dmpar_field_halo_exch(domain, 'channelArea')
         call mpas_timer_stop("halo updates")
      endif

      ! =============
      ! Calculate total grounding line discharges
      ! =============
      block => domain % blocklist
      do while (associated(block))

         call calc_gl_totals(block, err_tmp)
         err = ior(err, err_tmp)

         block => block % next
      end do

      ! =============
      ! Update water layer thickness
      ! =============
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
         call mpas_pool_get_array(hydroPool, 'waterThickness', waterThickness)
         call mpas_pool_get_array(hydroPool, 'waterThicknessOld', waterThicknessOld)
         call mpas_pool_get_array(hydroPool, 'waterThicknessTendency', waterThicknessTendency)
         call mpas_pool_get_array(hydroPool, 'tillWaterThickness', Wtill)
         call mpas_pool_get_array(hydroPool, 'tillWaterThicknessOld', WtillOld)
         call mpas_pool_get_array(hydroPool, 'deltatSGH', deltatSGH)
         call mpas_pool_get_array(hydroPool, 'basalMeltInput', basalMeltInput)
         call mpas_pool_get_array(hydroPool, 'externalWaterInput', externalWaterInput)
         call mpas_pool_get_array(hydroPool, 'divergence', divergence)
         call mpas_pool_get_array(hydroPool, 'divergenceChannel', divergenceChannel)
         call mpas_pool_get_array(hydroPool, 'channelAreaChangeCell', channelAreaChangeCell)
         call mpas_pool_get_array(hydroPool, 'channelMeltInputCell', channelMeltInputCell)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)

         waterThicknessOld = waterThickness
         waterThickness = waterThicknessOld + deltatSGH * ( &
             (basalMeltInput + externalWaterInput) / rho_water &
             + channelMeltInputCell &
             - divergence  &
             - divergenceChannel - channelAreaChangeCell  &
             - (Wtill - WtillOld) / deltatSGH )
         waterThickness = waterThickness * li_mask_is_grounded_ice_int(cellMask)  ! zero in non-grounded locations
         waterThickness = max(0.0_RKIND, waterThickness)
         divergence = divergence * li_mask_is_grounded_ice_int(cellMask)  ! zero in non-grounded locations for more convenient viz
         waterThicknessTendency = (waterThickness - waterThicknessOld) / deltatSGH

         block => block % next
      end do

      ! =============
      ! Calculate pressure field
      ! =============
      ! Note: In Bueler and van Pelt, pressure is updated before waterThickness
      ! (step vii on page 1625).  It does not matter which order they are calculated
      ! because the update calculations are formulated on intermediate variables
      ! such that neither depend directly on the other in the actual update.
      ! However, calc_pressure also calculates the full hydropotential (including
      ! the water thickness elevation head) which *is* a direct function of water thickness,
      ! making the ordering of that calculation important.
      ! Because the pressure variables calculated here are meant to be used on the
      ! following timestep (forward Euler), the full hydropotential should be using
      ! the updated waterThickness, not the old one.  To achieve that, calc_pressure
      ! is called here, after waterThickness has been updated.
      ! Note that the full hydropotential is only used by the channel model, so this
      ! ordering choice is only needed to support channels.  Without it, runs with
      ! channels would use an out of date waterThickness in the hydropotential and
      ! do not restart correctly due to the order of operations mismatch.
      block => domain % blocklist
      do while (associated(block))

         call calc_pressure(block, err_tmp)
         err = ior(err, err_tmp)
        
         call calc_waterPressureSmooth(block, err_tmp) !compute smoothed version of waterPressure
         err = ior(err,err_tmp)

         block => block % next
      end do

      !updates halos for waterPressure
      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'waterPressureSmooth')
      call mpas_timer_stop("halo updates")

      !Average effectivePressureSGH over coupling interval for use in dynamics model 
      block => domain % blocklist
      do while (associated(block))
        
         call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
         call mpas_pool_get_array(hydroPool, 'effectivePressureSGH', effectivePressureSGH)
         call mpas_pool_get_array(hydroPool, 'effectivePressure', effectivePressure)
         call mpas_pool_get_array(hydroPool, 'deltatSGH', deltatSGH)

         effectivePressure = (effectivePressure * deltatSGHaccumulated + effectivePressureSGH * deltatSGH) / (deltatSGHaccumulated + deltatSGH)
   
         block => block % next
      end do

      deltatSGHaccumulated = deltatSGHaccumulated + deltatSGH
      
      ! =============
      ! =============
      ! =============
      end do timecycling ! while timeLeft>0
      ! =============
      ! =============
      ! =============

      call mpas_log_write("Hydro model subcycled $i times.", intArgs=(/numSubCycles/))


      ! === error check
      if (err > 0) then
          call mpas_log_write("An error has occurred in li_SGH_solve.", MPAS_LOG_ERR)
      endif

   !--------------------------------------------------------------------
   end subroutine li_SGH_solve



!***********************************************************************
!
!  routine li_SGH_finalize
!
!> \brief   Finalize SGH
!> \author  Matt Hoffman
!> \date    27 June 2016
!> \details
!>  This routine finalizes the subglacial hydro model
!-----------------------------------------------------------------------
   subroutine li_SGH_finalize(domain, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      ! Pools pointers
      logical, pointer :: config_SGH
      type (block_type), pointer :: block
      integer :: err_tmp

      err = 0
      err_tmp = 0

      call mpas_pool_get_config(liConfigs, 'config_SGH', config_SGH)
      if (.not. config_SGH) then
         return
      endif

      block => domain % blocklist
      do while (associated(block))

         block => block % next
      end do

      ! === error check
      if (err > 0) then
          call mpas_log_write("An error has occurred in li_SGH_finalize.", MPAS_LOG_ERR)
      endif

   !--------------------------------------------------------------------
   end subroutine li_SGH_finalize



   !--------------------------------------------------------------------
   !--------------------------------------------------------------------
   ! Local routines
   !--------------------------------------------------------------------
   !--------------------------------------------------------------------



!***********************************************************************
!
!  routine calc_edge_quantities
!
!> \brief   Calculate SGH fields on edges
!> \author  Matt Hoffman
!> \date    27 June 2016
!> \details
!>  This routine calculates needed SGH fields on edges
!-----------------------------------------------------------------------
   subroutine calc_edge_quantities(block, err)

      use mpas_geometry_utils, only: mpas_cells_to_points_using_baryweights
      use mpas_vector_operations, only: mpas_tangential_vector_1d
      use li_setup, only: li_cells_to_vertices_1dfield_using_kiteAreas

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (block_type), intent(inout) :: block    !< Input/Output: block object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      ! Pools pointers
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: hydroPool
      real (kind=RKIND), dimension(:), pointer :: bedTopography
      real (kind=RKIND), dimension(:), pointer :: thickness
      real (kind=RKIND), dimension(:), pointer :: hydropotentialBase
      real (kind=RKIND), dimension(:), pointer :: hydropotential
      real (kind=RKIND), dimension(:), pointer :: hydropotentialBaseVertex
      real (kind=RKIND), dimension(:), pointer :: hydropotentialVertex
      real (kind=RKIND), dimension(:), pointer :: waterPressure
      real (kind=RKIND), dimension(:), pointer :: waterThicknessEdge
      real (kind=RKIND), dimension(:), pointer :: waterThicknessEdgeUpwind
      real (kind=RKIND), dimension(:), pointer :: waterThickness
      real (kind=RKIND), dimension(:), pointer :: hydropotentialBaseSlopeNormal
      real (kind=RKIND), dimension(:), pointer :: hydropotentialSlopeNormal
      real (kind=RKIND), dimension(:), pointer :: waterPressureSlopeNormal
      real (kind=RKIND), dimension(:), pointer :: hydropotentialBaseSlopeTangent
      real (kind=RKIND), dimension(:), pointer :: hydropotentialSlopeTangent
      real (kind=RKIND), dimension(:), pointer :: gradMagPhiBaseEdge
      real (kind=RKIND), dimension(:), pointer :: gradMagPhiEdge
      real (kind=RKIND), dimension(:), pointer :: effectiveConducEdge
      real (kind=RKIND), dimension(:), pointer :: diffusivity
      real (kind=RKIND), dimension(:), pointer :: dcEdge
      real (kind=RKIND), dimension(:), pointer :: dvEdge
      real (kind=RKIND), dimension(:), pointer :: waterVelocity
      real (kind=RKIND), dimension(:), pointer :: waterFlux
      real (kind=RKIND), dimension(:), pointer :: waterFluxAdvec
      real (kind=RKIND), dimension(:), pointer :: waterFluxDiffu
      real (kind=RKIND), dimension(:), pointer :: waterPressureSmooth
      integer, dimension(:), pointer :: hydroMarineMarginMask
      integer, dimension(:), pointer :: hydroTerrestrialMarginMask
      integer, dimension(:), pointer :: waterFluxMask
      integer, dimension(:,:), pointer :: edgeSignOnCell
      integer, dimension(:), pointer :: cellMask
      integer, dimension(:), pointer :: edgeMask
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer, dimension(:,:), pointer :: edgesOnCell
      integer, dimension(:,:), pointer :: verticesOnEdge
      integer, dimension(:,:), pointer :: baryCellsOnVertex
      real (kind=RKIND), dimension(:,:), pointer :: baryWeightsOnVertex
      real (kind=RKIND), pointer :: alpha, beta
      real (kind=RKIND), pointer :: conduc_coeff
      real (kind=RKIND), pointer :: conduc_coeff_drowned
      real (kind=RKIND), pointer :: bedRoughMax
      real (kind=RKIND) :: conduc_coeff_wtd
      character (len=StrKIND), pointer :: config_SGH_tangent_slope_calculation
      real (kind=RKIND), pointer :: config_sea_level
      real (kind=RKIND), pointer :: rhoo
      integer, pointer :: nEdges
      integer, pointer :: nCells
      integer, pointer :: nVertices
      integer :: iEdge, cell1, cell2
      integer :: i, j, iVertex, iCell
      real (kind=RKIND) :: velSign
      integer :: numGroundedCells
      real(kind=RKIND), parameter :: SMALL_CONDUC = 1.0e-30_RKIND
      integer :: err_tmp
      
      err = 0
      err_tmp = 0

      ! Get pools things
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
      call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)

      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)

      call mpas_pool_get_config(liConfigs, 'config_SGH_alpha', alpha)
      call mpas_pool_get_config(liConfigs, 'config_SGH_beta', beta)
      call mpas_pool_get_config(liConfigs, 'config_SGH_conduc_coeff', conduc_coeff)
      call mpas_pool_get_config(liConfigs, 'config_SGH_conduc_coeff_drowned', conduc_coeff_drowned)
      call mpas_pool_get_config(liConfigs, 'config_SGH_bed_roughness_max', bedRoughMax)
      call mpas_pool_get_config(liConfigs, 'config_SGH_tangent_slope_calculation', config_SGH_tangent_slope_calculation)
      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
      call mpas_pool_get_config(liConfigs, 'config_ocean_density', rhoo)
      
      call mpas_pool_get_array(hydroPool, 'waterThickness', waterThickness)
      call mpas_pool_get_array(hydroPool, 'waterPressure', waterPressure)
      call mpas_pool_get_array(hydroPool, 'hydropotentialBase', hydropotentialBase)
      call mpas_pool_get_array(hydroPool, 'hydropotential', hydropotential)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(hydroPool, 'waterThicknessEdge', waterThicknessEdge)
      call mpas_pool_get_array(hydroPool, 'waterThicknessEdgeUpwind', waterThicknessEdgeUpwind)
      call mpas_pool_get_array(hydroPool, 'hydropotentialBaseSlopeNormal', hydropotentialBaseSlopeNormal)
      call mpas_pool_get_array(hydroPool, 'hydropotentialSlopeNormal', hydropotentialSlopeNormal)
      call mpas_pool_get_array(hydroPool, 'waterPressureSlopeNormal', waterPressureSlopeNormal)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
      call mpas_pool_get_array(hydroPool, 'hydropotentialBaseSlopeTangent', hydropotentialBaseSlopeTangent)
      call mpas_pool_get_array(hydroPool, 'hydropotentialSlopeTangent', hydropotentialSlopeTangent)
      call mpas_pool_get_array(hydroPool, 'gradMagPhiBaseEdge', gradMagPhiBaseEdge)
      call mpas_pool_get_array(hydroPool, 'gradMagPhiEdge', gradMagPhiEdge)
      call mpas_pool_get_array(hydroPool, 'effectiveConducEdge', effectiveConducEdge)
      call mpas_pool_get_array(hydroPool, 'diffusivity', diffusivity)
      call mpas_pool_get_array(hydroPool, 'waterVelocity', waterVelocity)
      call mpas_pool_get_array(hydroPool, 'waterFlux', waterFlux)
      call mpas_pool_get_array(hydroPool, 'waterFluxAdvec', waterFluxAdvec)
      call mpas_pool_get_array(hydroPool, 'waterFluxDiffu', waterFluxDiffu)
      call mpas_pool_get_array(hydroPool, 'waterFluxMask', waterFluxMask)
      call mpas_pool_get_array(hydroPool, 'hydroMarineMarginMask', hydroMarineMarginMask)
      call mpas_pool_get_array(hydroPool, 'hydroTerrestrialMarginMask', hydroTerrestrialMarginMask)
      call mpas_pool_get_array(geometryPool, 'edgeMask', edgeMask)
      call mpas_pool_get_array(hydroPool, 'waterPressureSmooth', waterPressureSmooth) 
      
      do iEdge = 1, nEdges
         cell1 = cellsOnEdge(1, iEdge)
         cell2 = cellsOnEdge(2, iEdge)

         !waterThicknessEdge(iEdge) = 0.5_RKIND * ( waterThickness(cell1) + waterThickness(cell2) )
         ! This version ignores the thickness where there is no grounded ice (one-sided average at margin)
         numGroundedCells = li_mask_is_grounded_ice_int(cellMask(cell1)) + li_mask_is_grounded_ice_int(cellMask(cell2))
         if (numGroundedCells > 0) then
            ! Assuming here that waterThickness has been zeroed in non-grounded locations
            waterThicknessEdge(iEdge) = ( waterThickness(cell1) + waterThickness(cell2) ) / real(numGroundedCells)
         else
            waterThicknessEdge(iEdge) = 0.0_RKIND
         endif

         hydropotentialBaseSlopeNormal(iEdge) = (hydropotentialBase(cell2) - hydropotentialBase(cell1)) / dcEdge(iEdge)
         hydropotentialSlopeNormal(iEdge) = (hydropotential(cell2) - hydropotential(cell1)) / dcEdge(iEdge)

         waterPressureSlopeNormal(iEdge) = (waterPressureSmooth(cell2) - waterPressureSmooth(cell1)) / dcEdge(iEdge)
      end do

      ! At terrestrial margin, ignore the downslope bed topography gradient.  Including it can lead to unrealistically large
      ! hydropotential gradients and unstable channel growth.
      ! The hydropotential at the terrestrial margin should be determined by the geometry
      ! at the edge of the cell in a 1-sided sense. For hydropotentialBase = rho*g*Zb + Pw, this means hydropotentialBaseSlopeNormal
      ! at the terrestrial margin is equal to Pw/dcEdge. 
      ! This one-sided implementation also creates outflowing conditions at terrestrial boundary
      do iEdge = 1, nEdges
         if (hydroTerrestrialMarginMask(iEdge) == 1) then
            cell1 = cellsOnEdge(1, iEdge)
            cell2 = cellsOnEdge(2, iEdge)
            if (li_mask_is_grounded_ice(cellMask(cell1))) then ! cell2 is the icefree cell - replace phi there with cell1 Phig
               
               hydropotentialBaseSlopeNormal(iEdge) = - waterPressure(cell1) / dcEdge(iEdge)
               hydropotentialSlopeNormal(iEdge) = - (rho_water * gravity * waterThickness(cell1) + waterPressure(cell1)) / dcEdge(iEdge)

            else ! cell1 is the icefree cell - replace phi there with cell2 Phig
                  
               hydropotentialBaseSlopeNormal(iEdge) = waterPressure(cell2) / dcEdge(iEdge) 
               hydropotentialSlopeNormal(iEdge) = (rho_water * gravity * waterThickness(cell2) + waterPressure(cell2)) / dcEdge(iEdge)

            endif ! which cell is icefree
         endif ! if edge of grounded ice
      end do

      ! Disallow inflow from the marine margin by imposing a minimum outflowing hydropotential gradient at the grounding line.
      do iEdge = 1, nEdges
         if ( hydroMarineMarginMask(iEdge) == 1) then
            cell1 = cellsOnEdge(1, iEdge)
            cell2 = cellsOnEdge(2, iEdge)
            if (li_mask_is_grounded_ice(cellMask(cell1))) then ! cell2 is the cell outside the hydro domain
               
               if (hydropotentialBaseSlopeNormal(iEdge) > -MIN_PHISLOPE_GL) then
                  hydropotentialBaseSlopeNormal(iEdge) = -MIN_PHISLOPE_GL
               endif
               if (hydropotentialSlopeNormal(iEdge) > -MIN_PHISLOPE_GL) then
                  hydropotentialSlopeNormal(iEdge) = -MIN_PHISLOPE_GL
               endif

            else ! cell1 is the cell outside the hydro domain

               if (hydropotentialBaseSlopeNormal(iEdge) < MIN_PHISLOPE_GL) then
                  hydropotentialBaseSlopeNormal(iEdge) = MIN_PHISLOPE_GL
               endif
               if (hydropotentialSlopeNormal(iEdge) < MIN_PHISLOPE_GL) then
                  hydropotentialSlopeNormal(iEdge) = MIN_PHISLOPE_GL
               endif

            endif ! which cell is icefree
         endif ! if edge of grounded ice
      end do

      ! zero gradients along zero flux boundaries
      do iEdge = 1, nEdges
         if (waterFluxMask(iEdge) == 2) then
            hydropotentialBaseSlopeNormal(iEdge) = 0.0_RKIND
            hydropotentialSlopeNormal(iEdge) = 0.0_RKIND
            waterPressureSlopeNormal(iEdge) = 0.0_RKIND
         endif
      end do

      ! Calculate tangent slope of hydropotential and hydropotentialBase - three possible methods to consider

      ! Calculate hydropotentialBaseVertex if needed
      call mpas_pool_get_array(hydroPool, 'hydropotentialBaseVertex', hydropotentialBaseVertex)
         ! < this array could be protected by logic if desired
      ! calculate hydropotentialVertex if needed
      call mpas_pool_get_array(hydroPool, 'hydropotentialVertex', hydropotentialVertex)
         ! < this array could be protected by logic if desired
      select case (trim(config_SGH_tangent_slope_calculation))
      case ('from_vertex_barycentric')
         call mpas_pool_get_array(meshPool, 'baryCellsOnVertex', baryCellsOnVertex)
         call mpas_pool_get_array(meshPool, 'baryWeightsOnVertex', baryWeightsOnVertex)
         call mpas_cells_to_points_using_baryweights(meshPool, baryCellsOnVertex(:, 1:nVertices), &
              baryWeightsOnVertex(:, 1:nVertices), hydropotentialBase, hydropotentialBaseVertex(1:nVertices), err_tmp)
              err = ior(err, err_tmp)
         call mpas_cells_to_points_using_baryweights(meshPool, baryCellsOnVertex(:, 1:nVertices), &
              baryWeightsOnVertex(:, 1:nVertices), hydropotential, hydropotentialVertex(1:nVertices), err_tmp)
      case ('from_vertex_barycentric_kiteareas')
         call li_cells_to_vertices_1dfield_using_kiteAreas(meshPool, hydropotentialBase, hydropotentialBaseVertex)
         call li_cells_to_vertices_1dfield_using_kiteAreas(meshPool, hydropotential, hydropotentialVertex)
      end select

      ! Now perform tangent slope calculation based on method chosen
      ! For the two vertex based methods, edges with vertices on the boundary of the mesh will have invalid values
      ! For the normal slope method, any edges on cells that are at the edge of the mesh will have
      ! contaminated values.
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      select case (trim(config_SGH_tangent_slope_calculation))
      case ('from_vertex_barycentric', 'from_vertex_barycentric_kiteareas')
         do iEdge = 1, nEdges
            ! Only calculate slope for edges that have ice on at least one side.
            if ( li_mask_is_ice(edgeMask(iEdge)) ) then
               hydropotentialBaseSlopeTangent(iEdge) = ( hydropotentialBaseVertex(verticesOnEdge(1,iEdge)) -  &
                     hydropotentialBaseVertex(verticesOnEdge(2,iEdge)) ) / dvEdge(iEdge)
               hydropotentialSlopeTangent(iEdge) = ( hydropotentialVertex(verticesOnEdge(1,iEdge)) -  &
                     hydropotentialVertex(verticesOnEdge(2,iEdge)) ) / dvEdge(iEdge)
            else
               hydropotentialBaseSlopeTangent(iEdge) = 0.0_RKIND
               hydropotentialSlopeTangent(iEdge) = 0.0_RKIND
            endif
            ! zero tangent slope at waterFluxMask==2
            if (waterFluxMask(iEdge) == 2) then
                hydropotentialBaseSlopeTangent(iEdge) = 0.0_RKIND
                hydropotentialSlopeTangent(iEdge) = 0.0_RKIND
             endif
         end do  ! edges
      case ('from_normal_slope')
         ! Do first with hydropotentialBase
         call mpas_tangential_vector_1d(hydropotentialBaseSlopeNormal, meshPool, &
         includeHalo=.false., tangentialVector=hydropotentialBaseSlopeTangent)
         ! Repeat for full hydropotential
         call mpas_tangential_vector_1d(hydropotentialSlopeNormal, meshPool, &
         includeHalo=.false., tangentialVector=hydropotentialSlopeTangent)
         ! ensure that edges that don't have ice on at least one side have tangent slope set to zero
         do iEdge = 1, nEdges
            if ( .not. li_mask_is_ice(edgeMask(iEdge)) ) then
               hydropotentialBaseSlopeTangent(iEdge) = 0.0_RKIND
               hydropotentialSlopeTangent(iEdge) = 0.0_RKIND
            endif
            if (waterFluxMask(iEdge) == 2) then
               hydropotentialBaseSlopeTangent(iEdge) = 0.0_RKIND
               hydropotentialSlopeTangent(iEdge) = 0.0_RKIND
            endif

         end do  ! edges
      case default
         call mpas_log_write('Invalid value for config_SGH_tangent_slope_calculation.', MPAS_LOG_ERR)
         err = 1
      end select

      ! calculate magnitude of gradient of Phi
      gradMagPhiEdge = sqrt(hydropotentialSlopeNormal**2 + hydropotentialSlopeTangent**2)

      ! calculate magnitude of gradient of hydropotentialBase
      gradMagPhiBaseEdge = sqrt(hydropotentialBaseSlopeNormal**2 +&
      hydropotentialBaseSlopeTangent**2)

      ! calculate effective conductivity on edges
      if (conduc_coeff_drowned > 0.0_RKIND) then
         ! Use a thickness weighted conductivity coeff. when water thickness exceeds bump height
         do iEdge = 1, nEdges
            if (gradMagPhiBaseEdge(iEdge) < SMALL_GRADPHI) then
               effectiveConducEdge(iEdge) = 0.0_RKIND
            else
               conduc_coeff_wtd = (conduc_coeff * min(waterThicknessEdge(iEdge), bedRoughMax) + &
                  conduc_coeff_drowned * max(waterThicknessEdge(iEdge) - bedRoughMax, 0.0_RKIND)) / &
                  (waterThicknessEdge(iEdge) + 1.0e-16_RKIND)  ! Regularization only applies where value doesn't matter

               effectiveConducEdge(iEdge) = conduc_coeff_wtd * waterThicknessEdge(iEdge)**(alpha-1.0_RKIND) * &
                  gradMagPhiBaseEdge(iEdge)**(beta - 2.0_RKIND)
            endif
         end do
      else
         do iEdge = 1, nEdges
            ! Just use a single conductivity coeff.
            if (gradMagPhiBaseEdge(iEdge) < SMALL_GRADPHI) then
               effectiveConducEdge(iEdge) = 0.0_RKIND
            else
               effectiveConducEdge(iEdge) = conduc_coeff * waterThicknessEdge(iEdge)**(alpha-1.0_RKIND) *&
                  gradMagPhiBaseEdge(iEdge)**(beta - 2.0_RKIND)
            endif
         enddo
      endif

      where (effectiveConducEdge < SMALL_CONDUC)
              effectiveConducEdge = 0.0_RKIND
      end where

      ! calculate diffusivity on edges
      diffusivity(:) = rho_water * gravity * effectiveConducEdge(:) * waterThicknessEdge(:)

      do iEdge = 1, nEdges
         cell1 = cellsOnEdge(1, iEdge)
         cell2 = cellsOnEdge(2, iEdge)
         waterVelocity(iEdge) = -1.0_RKIND * effectiveConducEdge(iEdge) * hydropotentialBaseSlopeNormal(iEdge)
         velSign = sign(1.0_RKIND, waterVelocity(iEdge))
         waterThicknessEdgeUpwind(iEdge) = max(velSign * waterThickness(cell1),   &
                     velSign * (-1.0_RKIND) * waterThickness(cell2))

         ! advective flux
         waterFluxAdvec(iEdge) = waterVelocity(iEdge) * waterThicknessEdgeUpwind(iEdge)

         ! diffusive flux
         ! Note: Water thickness at the GL is undefined.  I don't think we want to assume it's 0, and
         ! we also don't want to assume it's the ocean water column height, because that would imply
         ! a diffusive flux inward, which is undesirable. So disabling diffusion at the GL.
         if (hydroMarineMarginMask(iEdge) == 1) then
            diffusivity(iEdge) = 0.0_RKIND
            waterFluxDiffu(iEdge) = 0.0_RKIND
         else
            waterFluxDiffu(iEdge) = -1.0_RKIND * diffusivity(iEdge) * (waterThickness(cell2) - waterThickness(cell1)) &
               / dcEdge(iEdge)
         endif
      end do
      where (waterFluxMask == 2)
         diffusivity = 0.0_RKIND
         waterFluxAdvec = 0.0_RKIND
         waterFluxDiffu = 0.0_RKIND
         waterVelocity = 0.0_RKIND
      end where
      where (waterThicknessEdgeUpwind == 0.0_RKIND)  ! if no water to give, should have no flux!
         waterFluxAdvec = 0.0_RKIND ! Should already be 0
         waterFluxDiffu = 0.0_RKIND ! Might not be 0
         waterVelocity = 0.0_RKIND
      end where

      waterFlux(:) = waterFluxAdvec(:) + waterFluxDiffu(:)

   !--------------------------------------------------------------------
   end subroutine calc_edge_quantities




!***********************************************************************
!
!  routine check_timestep
!
!> \brief   Calculate SGH timesteps and check current timestep
!> \author  Matt Hoffman
!> \date    7 July 2016
!> \details
!>  This routine calculates the three timesteps associated with the
!>  SGH solver and compares them to the current model timestep.
!-----------------------------------------------------------------------
   subroutine check_timestep(domain, timeLeft, numSubCycles, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      integer, intent(in) :: numSubCycles  !< Input: number of subcycles taken so far

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain    !< Input/Output: domain object
      real (kind=RKIND), intent(inout) :: timeLeft  !< Input/Output: time remaining for subcycling (seconds)

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      ! Pools pointers
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: hydroPool
      real (kind=RKIND), dimension(:), pointer :: waterVelocity
      real (kind=RKIND), dimension(:), pointer :: channelVelocity
      real (kind=RKIND), dimension(:), pointer :: diffusivity
      real (kind=RKIND), dimension(:), pointer :: channelDiffusivity
      real (kind=RKIND), dimension(:), pointer :: dcEdge
      integer, dimension(:), pointer :: indexToEdgeID
      real (kind=RKIND), pointer :: deltatSGH
      real (kind=RKIND), pointer :: deltatSGHadvec
      real (kind=RKIND), pointer :: deltatSGHdiffu
      real (kind=RKIND), pointer :: deltatSGHpressure
      real (kind=RKIND), pointer :: deltatSGHadvecChannel
      real (kind=RKIND), pointer :: deltatSGHdiffuChannel
      real (kind=RKIND), pointer :: porosity
      type (block_type), pointer :: block
      real (kind=RKIND), pointer :: deltat
      integer, pointer :: nEdgesSolve
      logical, pointer :: config_SGH_chnl_active, config_SGH_chnl_include_DCFL
      ! in the following variables, "Proc" indicates the value on the current processor,
      ! and "Block" indicates value on current block
      real (kind=RKIND) :: dtSGHadvecBlock, dtSGHadvecProc
      real (kind=RKIND) :: dtSGHdiffuBlock, dtSGHdiffuProc
      real (kind=RKIND) :: dtSGHpressureBlock, dtSGHpressureProc
      real (kind=RKIND) :: dtSGHadvecChanBlock, dtSGHadvecChanProc
      real (kind=RKIND) :: dtSGHdiffuChanBlock, dtSGHdiffuChanProc
      integer :: err_tmp
      real(kind=RKIND), parameter :: bigNumber = 1.0e36_RKIND
      real(kind=RKIND), pointer :: CFLfraction
      real(kind=RKIND), pointer :: maxDt
      real(kind=RKIND) :: proposedDt
      real(kind=RKIND) :: masterDt
      real(kind=RKIND), dimension(5) :: localMinValues, reducedMinValues
      integer, dimension(5) :: localLocValues, reducedLocValues
      integer :: iEdge
      real(kind=RKIND), save :: minDtSoFar
      integer, save :: numDtLess1
      logical :: printDtInfo


      err = 0
      err_tmp = 0

      printDtInfo = .false.
      if (numSubCycles == 1) then
         minDtSoFar = 1.0e36_RKIND
         numDtLess1 = 0
      endif

      call mpas_pool_get_config(liConfigs, 'config_SGH_englacial_porosity', porosity)
      call mpas_pool_get_config(liConfigs, 'config_SGH_chnl_active', config_SGH_chnl_active)

      dtSGHadvecProc = bigNumber
      dtSGHdiffuProc = bigNumber
      dtSGHpressureProc = bigNumber
      dtSGHadvecChanProc = bigNumber
      dtSGHdiffuChanProc = bigNumber

      ! ---
      ! Find local (block) limiting dt for various CFL conditions
      ! ---
      block => domain % blocklist
      do while (associated(block))

         ! Get pools things
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)

         call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)

         call mpas_pool_get_array(hydroPool, 'deltatSGH', deltatSGH)
         call mpas_pool_get_array(hydroPool, 'deltatSGHadvec', deltatSGHadvec)
         call mpas_pool_get_array(hydroPool, 'deltatSGHdiffu', deltatSGHdiffu)
         call mpas_pool_get_array(hydroPool, 'deltatSGHpressure', deltatSGHpressure)
         call mpas_pool_get_array(hydroPool, 'deltatSGHadvecChannel', deltatSGHadvecChannel)
         call mpas_pool_get_array(hydroPool, 'deltatSGHdiffuChannel', deltatSGHdiffuChannel)
         call mpas_pool_get_array(hydroPool, 'waterVelocity', waterVelocity)
         call mpas_pool_get_array(hydroPool, 'channelVelocity', channelVelocity)
         call mpas_pool_get_array(hydroPool, 'diffusivity', diffusivity)
         call mpas_pool_get_array(hydroPool, 'channelDiffusivity', channelDiffusivity)
         call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
         call mpas_pool_get_array(meshPool, 'indexToEdgeID', indexToEdgeID)

         dtSGHadvecBlock = bigNumber
         dtSGHdiffuBlock = bigNumber
         dtSGHpressureBlock = bigNumber
         dtSGHadvecChanBlock = bigNumber
         dtSGHdiffuChanBlock = bigNumber

         do iEdge = 1, nEdgesSolve
            ! Calculate advective CFL-limited time step
            if (abs(waterVelocity(iEdge)) > 0.0) then
               dtSGHadvecBlock = min(dtSGHadvecBlock, 0.5_RKIND * dcEdge(iEdge) / abs(waterVelocity(iEdge)))
            endif

            if (diffusivity(iEdge) > 0.0) then
               ! Calculate diffusive CFL-limited time step
               dtSGHdiffuBlock = min(dtSGHdiffuBlock, 0.25_RKIND * dcEdge(iEdge)**2 / diffusivity(iEdge))
               ! Calculate pressure limited time step
               dtSGHpressureBlock = min(dtSGHpressureBlock, 1.0_RKIND * porosity * dcEdge(iEdge)**2 &
                      / (2.0_RKIND * diffusivity(iEdge)))
            endif

            if (config_SGH_chnl_active) then
               if (abs(channelVelocity(iEdge)) > 0.0) then
                  ! Calculate channel advection limited time step
                  dtSGHadvecChanBlock = min(dtSGHadvecChanBlock, 0.5_RKIND * dcEdge(iEdge) / (abs(channelVelocity(iEdge))))
               endif
               ! Calculate channel diffusion limited time step
               if (channelDiffusivity(iEdge) > 0.0) then
                  dtSGHdiffuChanBlock = min(dtSGHdiffuChanBlock, 0.25_RKIND * dcEdge(iEdge)**2 / channelDiffusivity(iEdge))
               endif
            endif
         enddo

         dtSGHadvecProc = min(dtSGHadvecProc, dtSGHadvecBlock)
         dtSGHdiffuProc = min(dtSGHdiffuProc, dtSGHdiffuBlock)
         dtSGHpressureProc = min(dtSGHpressureProc, dtSGHpressureBlock)
         dtSGHadvecChanProc = min(dtSGHadvecChanProc, dtSGHadvecChanBlock)
         dtSGHdiffuChanProc = min(dtSGHdiffuChanProc, dtSGHdiffuChanBlock)

         ! Master deltat is needed below, so grab it in this block loop
         call mpas_pool_get_array(meshPool, 'deltat', deltat)

         block => block % next
      end do


      ! ---
      ! reduce across procs
      ! ---
      localMinValues(1) = dtSGHadvecProc
      localMinValues(2) = dtSGHdiffuProc
      localMinValues(3) = dtSGHpressureProc
      if (config_SGH_chnl_active) then
         localMinValues(4) = dtSGHadvecChanProc
         localMinValues(5) = dtSGHdiffuChanProc
         call mpas_timer_start("global reduce")
         call mpas_dmpar_min_real_array(domain % dminfo, 5, localMinValues, reducedMinValues)
         call mpas_timer_stop("global reduce")
         deltatSGHadvecChannel = reducedMinValues(4)
         deltatSGHdiffuChannel = reducedMinValues(5)
      else
         call mpas_timer_start("global reduce")
         call mpas_dmpar_min_real_array(domain % dminfo, 3, localMinValues(1:3), reducedMinValues(1:3))
         call mpas_timer_stop("global reduce")
      endif
      deltatSGHadvec = reducedMinValues(1)
      deltatSGHdiffu = reducedMinValues(2)
      deltatSGHpressure = reducedMinValues(3)


      ! ---
      ! Find and set the new hydro subcycling dt
      ! ---
      call mpas_pool_get_config(liConfigs, 'config_SGH_adaptive_timestep_fraction', CFLfraction)
      call mpas_pool_get_config(liConfigs, 'config_SGH_max_adaptive_timestep', maxDt)
      call mpas_pool_get_config(liConfigs, 'config_SGH_chnl_include_DCFL', config_SGH_chnl_include_DCFL)

      ! Find smallest of 3 or 5 limiting time steps
      proposedDt = min(deltatSGHadvec, deltatSGHdiffu, deltatSGHpressure)
      if (config_SGH_chnl_active) then
         proposedDt = min(proposedDt, deltatSGHadvecChannel)
         if (config_SGH_chnl_include_DCFL) then
            proposedDt = min(proposedDt, deltatSGHdiffuChannel)
         endif
      endif
      proposedDt = proposedDt * CFLfraction

      if (proposedDt < 1.0_RKIND) then
         numDtLess1 = numDtLess1 + 1
      endif
      minDtSoFar = min(minDtSoFar, proposedDt)

      ! Write out dt stats every now and then
      if (mod(numSubCycles, 2000)==0) then
         call mpas_log_write("SGH subcycle #$i. Time left=$r s.", intArgs=(/numSubCycles/), realArgs=(/timeLeft/))
         printDtInfo=.true.
      endif

      ! Write out dt info on the final subcycle
      if (proposedDt >= timeLeft) then
         call mpas_log_write("SGH completed subcycling with $i timesteps.", intArgs=(/numSubCycles/))
         printDtInfo=.true.
      endif

      if (proposedDt < 1.0E-2_RKIND) then
         call mpas_log_write("CFL conditions are limiting subglacial hydrology timestep to <1e-2s ($r) on subcycle $i.", &
                 MPAS_LOG_WARN, realArgs=(/proposedDt/), intArgs=(/numSubCycles/))
         printDtInfo=.true.
      endif

      if (proposedDt < 1.0E-4_RKIND) then
         call mpas_log_write("CFL conditions are limiting subglacial hydrology timestep to <1e-4s.", MPAS_LOG_ERR)
         printDtInfo=.true.
         err = ior(err, 1)
      endif

      if (printDtInfo) then
         ! Determine location of limiting Dt
         ! NOTE: ignoring multiple blocks here.
         if (localMinValues(1) == reducedMinValues(1)) then ! limiting edge is on this proc
            localLocValues(1) = indexToEdgeID(minloc(dcEdge(1:nEdgesSolve) / (abs(waterVelocity(1:nEdgesSolve)) &
               + 1.0e-12_RKIND), dim=1))
         else
            localLocValues(1) = -1
         endif

         if (localMinValues(2) == reducedMinValues(2)) then ! limiting edge is on this proc
            localLocValues(2) = indexToEdgeID(minloc(dcEdge(1:nEdgesSolve)**2 / (diffusivity(1:nEdgesSolve) + 1.0e-12_RKIND), &
                    dim=1))
         else
            localLocValues(2) = -1
         endif

         if (localMinValues(3) == reducedMinValues(3)) then ! limiting edge is on this proc
            localLocValues(3) = indexToEdgeID(minloc(porosity * dcEdge(1:nEdgesSolve)**2 &
                      / (2.0_RKIND * diffusivity(1:nEdgesSolve) + 1.0e-12_RKIND), dim=1))
         else
            localLocValues(3) = -1
         endif

         if (config_SGH_chnl_active) then
            if (localMinValues(4) == reducedMinValues(4)) then ! limiting edge is on this proc
               localLocValues(4) = indexToEdgeID(minloc(dcEdge(1:nEdgesSolve) / (abs(channelVelocity(1:nEdgesSolve)) &
                  + 1.0e-12_RKIND), dim=1))
            else
               localLocValues(4) = -1
            endif

            if (localMinValues(5) == reducedMinValues(5)) then ! limiting edge is on this proc
               localLocValues(5) = indexToEdgeID(minloc(dcEdge(1:nEdgesSolve)**2 / (channelDiffusivity(1:nEdgesSolve) + &
                  1.0e-12_RKIND), dim=1))
            else
               localLocValues(5) = -1
            endif

            ! global reduce of locations
            call mpas_timer_start("global reduce")
            call mpas_dmpar_max_int_array(domain % dminfo, 5, localLocValues, reducedLocValues)
            call mpas_timer_stop("global reduce")
         else
            call mpas_timer_start("global reduce")
            call mpas_dmpar_max_int_array(domain % dminfo, 3, localLocValues, reducedLocValues)
            call mpas_timer_stop("global reduce")
         endif


         if (config_SGH_chnl_active) then
            call mpas_log_write("deltatSGH: used=$r, prev=$r, advec=$r, diffu=$r, pressure=$r, advecChannel=$r, diffuChannel=$r",&
              realArgs=(/proposedDt, deltatSGH, deltatSGHadvec, deltatSGHdiffu, deltatSGHpressure, deltatSGHadvecChannel, &
                         deltatSGHdiffuChannel/))
            call mpas_log_write("Limiting edge IDs: advec=$i, diffu=$i, pressure=$i, advecChannel=$i, diffuChannel=$i", &
               intArgs=reducedLocValues(1:5))
         else
            call mpas_log_write("deltatSGH: used=$r, previous=$r, advec=$r, diffu=$r, pressure=$r", &
              realArgs=(/proposedDt, deltatSGH, deltatSGHadvec, deltatSGHdiffu, deltatSGHpressure/))
            call mpas_log_write("Limiting edge IDs: advec=$i, diffu=$i, pressure=$i", &
               intArgs=reducedLocValues(1:3))
         endif
         call mpas_log_write("  min hydro dt used=$r; There were $i subcyles with dt<1s.", &
                 realArgs=(/minDtSoFar/), intArgs=(/numDtLess1/))
      endif

      ! Don't exceed the maximum allowed hydro time step
      proposedDt = min(proposedDt, maxDt)
      ! Don't let SGH time step exceed time left in the master model dt
      proposedDt = min(proposedDt, timeLeft)

      !call mpas_log_write('  Setting SGH time step to: $r seconds', realArgs=(/proposedDt/))

      timeLeft = timeLeft - proposedDt
      !call mpas_log_write("dt=$r, new TIMELEFT=$r", realArgs=(/proposedDt, timeLeft/))


      ! ---
      ! Assign new time step to all blocks
      ! ---
      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
         call mpas_pool_get_array(hydroPool, 'deltatSGH', deltatSGH)
         deltatSGH = proposedDt

         block => block % next
      end do


      ! compare timesteps
      if (deltatSGH > deltatSGHadvec) then
         call mpas_log_write("deltatSGH > deltatSGHadvec  $r > $r", MPAS_LOG_WARN, realArgs=(/deltatSGH, deltatSGHadvec/))
      endif
      if (deltatSGH > deltatSGHdiffu) then
         call mpas_log_write("deltatSGH > deltatSGHdiffu  $r > $r", MPAS_LOG_WARN, realArgs=(/deltatSGH, deltatSGHdiffu/))
      endif
      if (deltatSGH > deltatSGHpressure) then
         call mpas_log_write("deltatSGH > deltatSGHpressure  $r > $r", MPAS_LOG_WARN, realArgs=(/deltatSGH, deltatSGHpressure/))
      endif
      if (config_SGH_chnl_active .and. (deltatSGH > deltatSGHadvecChannel)) then
         call mpas_log_write("deltatSGH > deltatSGHadvecChannel  $r > $r", MPAS_LOG_WARN, &
            realArgs=(/deltatSGH, deltatSGHadvecChannel/))
      endif
      if (config_SGH_chnl_active .and. config_SGH_chnl_include_DCFL .and. (deltatSGH > deltatSGHdiffuChannel)) then
         call mpas_log_write("deltatSGH > deltatSGHdiffuChannel  $r > $r", MPAS_LOG_WARN, &
            realArgs=(/deltatSGH, deltatSGHdiffuChannel/))
      endif


   !--------------------------------------------------------------------
   end subroutine check_timestep




!***********************************************************************
!
!  routine calc_pressure
!
!> \brief   Calculate SGH water pressure
!> \author  Matt Hoffman
!> \date    5 July 2016
!> \details
!>  This routine calculates SGH water pressure
!-----------------------------------------------------------------------
   subroutine calc_pressure(block, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (block_type), intent(inout) :: block    !< Input/Output: block object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      ! Pools pointers
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: hydroPool
      type (mpas_pool_type), pointer :: velocityPool
      real (kind=RKIND), dimension(:), pointer :: waterPressure
      real (kind=RKIND), dimension(:), pointer :: waterPressureOld
      real (kind=RKIND), dimension(:), pointer :: waterPressureTendency
      real (kind=RKIND), dimension(:), pointer :: waterThickness
      real (kind=RKIND), dimension(:), pointer :: effectivePressureSGH
      real (kind=RKIND), dimension(:), pointer :: zeroOrderSum
      real (kind=RKIND), dimension(:), pointer :: closingRate
      real (kind=RKIND), dimension(:), pointer :: openingRate
      real (kind=RKIND), dimension(:), pointer :: basalMeltInput
      real (kind=RKIND), dimension(:), pointer :: externalWaterInput
      real (kind=RKIND), dimension(:), pointer :: Wtill, WtillOld
      real (kind=RKIND), dimension(:), pointer :: divergence
      real (kind=RKIND), dimension(:), pointer :: basalSpeed
      real (kind=RKIND), dimension(:,:), pointer :: flowParamA
      real (kind=RKIND), dimension(:), pointer :: thickness
      real (kind=RKIND), dimension(:), pointer :: divergenceChannel
      real (kind=RKIND), dimension(:), pointer :: channelAreaChangeCell
      real (kind=RKIND), dimension(:), pointer :: bedTopography
      real (kind=RKIND), dimension(:), pointer :: iceThicknessHydro      
      integer, dimension(:), pointer :: hydroMarineMarginMask
      integer, dimension(:), pointer :: cellMask
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell
      real (kind=RKIND), pointer :: deltatSGH
      real (kind=RKIND), pointer :: bedRough, bedRoughMax
      real (kind=RKIND), pointer :: rhoi
      real (kind=RKIND), pointer :: creepCoeff
      real (kind=RKIND), pointer :: porosity
      integer, pointer :: nVertLevels
      integer, pointer :: nCells
      character (len=StrKIND), pointer :: config_SGH_pressure_calc
      real (kind=RKIND), pointer :: config_sea_level
      real (kind=RKIND), pointer :: rhoo
      integer :: err_tmp
      integer :: iCell
      integer :: iEdge
      logical :: onMarineMargin

      err = 0
      err_tmp = 0


      ! Get pools things
      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
      call mpas_pool_get_config(liConfigs, 'config_ocean_density', rhoo)

      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
      call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
      call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)

      call mpas_pool_get_config(liConfigs, 'config_ice_density', rhoi)
      call mpas_pool_get_config(liConfigs, 'config_SGH_bed_roughness', bedRough)
      call mpas_pool_get_config(liConfigs, 'config_SGH_bed_roughness_max', bedRoughMax)
      call mpas_pool_get_config(liConfigs, 'config_SGH_creep_coefficient', creepCoeff)
      call mpas_pool_get_config(liConfigs, 'config_SGH_englacial_porosity', porosity)
      call mpas_pool_get_config(liConfigs, 'config_SGH_pressure_calc', config_SGH_pressure_calc)

      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)

      call mpas_pool_get_array(hydroPool, 'effectivePressureSGH', effectivePressureSGH)
      call mpas_pool_get_array(hydroPool, 'waterPressure', waterPressure)
      call mpas_pool_get_array(hydroPool, 'waterPressureOld', waterPressureOld)
      call mpas_pool_get_array(hydroPool, 'waterPressureTendency', waterPressureTendency)
      call mpas_pool_get_array(hydroPool, 'waterThickness', waterThickness)
      call mpas_pool_get_array(hydroPool, 'zeroOrderSum', zeroOrderSum)
      call mpas_pool_get_array(hydroPool, 'closingRate', closingRate)
      call mpas_pool_get_array(hydroPool, 'openingRate', openingRate)
      call mpas_pool_get_array(hydroPool, 'deltatSGH', deltatSGH)
      call mpas_pool_get_array(hydroPool, 'basalMeltInput', basalMeltInput)
      call mpas_pool_get_array(hydroPool, 'externalWaterInput', externalWaterInput)
      call mpas_pool_get_array(hydroPool, 'tillWaterThickness', Wtill)
      call mpas_pool_get_array(hydroPool, 'tillWaterThicknessOld', WtillOld)
      call mpas_pool_get_array(hydroPool, 'divergence', divergence)
      call mpas_pool_get_array(hydroPool, 'divergenceChannel', divergenceChannel)
      call mpas_pool_get_array(hydroPool, 'channelAreaChangeCell', channelAreaChangeCell)
      call mpas_pool_get_array(hydroPool, 'hydroMarineMarginMask', hydroMarineMarginMask)
      call mpas_pool_get_array(velocityPool, 'basalSpeed', basalSpeed)
      call mpas_pool_get_array(velocityPool, 'flowParamA', flowParamA)
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
      call mpas_pool_get_array(hydroPool, 'iceThicknessHydro', iceThicknessHydro)

      openingRate(:) = bedRough * basalSpeed(:) * (bedRoughMax - waterThickness(:))
      !openingRate(:) = bedRough * basalSpeed(:) * (bedRoughMax - waterThickness(:)) + &
      !   basalMeltInput(:) / rhoi  ! Hewitt 2011 opening
      openingRate = max(0.0_RKIND, openingRate)

      closingRate(:) = creepCoeff * flowParamA(nVertLevels, :) * &
              effectivePressureSGH(:)**3 * waterThickness(:)
!      closingRate(:) = waterThickness(:) * effectivePressureSGH(:) / 1.0e13_RKIND
!          ! Hewitt 2011 creep closure form.  Denominator is ice viscosity

      zeroOrderSum = closingRate - openingRate + (basalMeltInput + externalWaterInput) / rho_water - &
         (Wtill - WtillOld) / deltatSGH

      waterPressureOld = waterPressure

      select case (trim(config_SGH_pressure_calc))
      case ('cavity')

         where (li_mask_is_grounded_ice(cellMask))        
            waterPressure = (zeroOrderSum - divergence - divergenceChannel - channelAreaChangeCell) * &
                   rho_water * gravity * deltatSGH / porosity + waterPressureOld
         elsewhere ((.not. (li_mask_is_grounded_ice(cellMask))) .and. (bedTopography > config_sea_level))
            waterPressure = 0.0_RKIND
         elsewhere ! should evaluate to ((.not. (li_mask_is_grounded_ice(cellMask))) .and. (bedTopography <= config_sea_level))
            waterPressure = gravity * rhoo * (config_sea_level - bedTopography)
         end where

      case ('overburden')

         where (li_mask_is_grounded_ice(cellMask))         
            waterPressure = rhoi * gravity * iceThicknessHydro
         elsewhere ((.not. (li_mask_is_grounded_ice(cellMask))) .and. (bedTopography > config_sea_level))
            waterPressure = 0.0_RKIND
         elsewhere ! should evaluate to ((.not. (li_mask_is_grounded_ice(cellMask))) .and. (bedTopography <= config_sea_level))
            waterPressure = gravity * rhoo * (config_sea_level - bedTopography)
         end where

      case default
         call mpas_log_write("Invalid option specified for config_SGH_pressure_calc:" // config_SGH_pressure_calc, MPAS_LOG_ERR)
         err = ior(err, 1)
      end select

      waterPressure = max(0.0_RKIND, waterPressure)
      where (li_mask_is_grounded_ice(cellMask))
         waterPressure = min(waterPressure, rhoi * gravity * iceThicknessHydro)
      end where
      
      do iCell = 1, nCells
        onMarineMargin = .false.
        do iEdge = 1, nEdgesOnCell(iCell)
           if (hydroMarineMarginMask(edgesOnCell(iEdge, iCell)) == 1) then
              onMarineMargin = .true.
              exit
           endif
        enddo
        if (onMarineMargin) then
           ! At marine margin, don't let water pressure fall below ocean pressure
           ! TODO: Not sure if this should include the water layer thickness term.  Leaving it off.
           if (waterPressure(iCell) < rho_water * gravity * (config_sea_level - bedTopography(iCell))) then
               waterPressure(iCell) = rho_water * gravity * (config_sea_level - bedTopography(iCell))
           endif
        endif
      enddo

      waterPressureTendency = (waterPressure - waterPressureOld) / deltatSGH

      call calc_pressure_diag_vars(block, err_tmp)
      err = ior(err, err_tmp)

   !--------------------------------------------------------------------
   end subroutine calc_pressure


!***********************************************************************
!
!  routine calc_pressure_diag_vars
!
!> \brief   Calculate SGH diagnostic variables related to  pressure
!> \author  Matt Hoffman
!> \date    5 July 2016
!> \details
!>  This routine calculates variables related to water pressure
!-----------------------------------------------------------------------
   subroutine calc_pressure_diag_vars(block, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (block_type), intent(inout) :: block    !< Input/Output: block object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      ! Pools pointers
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: hydroPool
      real (kind=RKIND), pointer :: rhoi, rhoo
      real (kind=RKIND), dimension(:), pointer :: waterPressure
      real (kind=RKIND), dimension(:), pointer :: bedTopography
      real (kind=RKIND), dimension(:), pointer :: hydropotentialBase
      real (kind=RKIND), dimension(:), pointer :: waterThickness
      real (kind=RKIND), dimension(:), pointer :: hydropotential
      real (kind=RKIND), dimension(:), pointer :: effectivePressureSGH
      real (kind=RKIND), dimension(:), pointer :: iceThicknessHydro
      integer, dimension(:), pointer :: cellMask
      real (kind=RKIND), pointer :: config_sea_level

      err = 0

      ! Get pools things
      call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
      call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)

      call mpas_pool_get_config(liConfigs, 'config_ice_density', rhoi)
      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
      call mpas_pool_get_config(liConfigs, 'config_ocean_density', rhoo)

      call mpas_pool_get_array(hydroPool, 'effectivePressureSGH', effectivePressureSGH)
      call mpas_pool_get_array(hydroPool, 'waterPressure', waterPressure)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
      call mpas_pool_get_array(hydroPool, 'hydropotentialBase', hydropotentialBase)
      call mpas_pool_get_array(hydroPool, 'waterThickness', waterThickness)
      call mpas_pool_get_array(hydroPool, 'hydropotential', hydropotential)
      call mpas_pool_get_array(hydroPool, 'iceThicknessHydro', iceThicknessHydro)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)

      effectivePressureSGH = rhoi * gravity * iceThicknessHydro - waterPressure
         ! < this should evalute to 0 for floating ice if Pw set correctly there.
      where (.not. (li_mask_is_grounded_ice(cellMask)))
         effectivePressureSGH = 0.0_RKIND  ! zero effective pressure where no ice to avoid confusion
      end where
    
      hydropotentialBase = rho_water * gravity * bedTopography + waterPressure
      where ((.not. li_mask_is_grounded_ice(cellMask)) .and. (bedTopography <= config_sea_level)) 
         hydropotentialBase = 0.0_RKIND !zero hydropotential in ocean
      elsewhere ((.not. li_mask_is_grounded_ice(cellMask)) .and. (bedTopography > config_sea_level))
         hydropotentialBase = rho_water * gravity * bedTopography ! for completeness, but won't matter with one-side slope calculations on terrestrial boundaries
      end where      

      ! hydropotential with water thickness
      hydropotential = hydropotentialBase + rho_water * gravity * waterThickness

   !--------------------------------------------------------------------
   end subroutine calc_pressure_diag_vars



!***********************************************************************
!
!  routine update_channel
!
!> \brief   Calculate SGH channel area
!> \author  Matt Hoffman
!> \date    28 July 2016
!> \details
!>  This routine updates the channel area in the subglacial hydrology model.
!>  It uses the conduit space evolution equation.
!-----------------------------------------------------------------------
   subroutine update_channel(block, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (block_type), intent(inout) :: block    !< Input/Output: block object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      ! Pools pointers
      type (mpas_pool_type), pointer :: hydroPool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: velocityPool
      type (mpas_pool_type), pointer :: geometryPool
      real (kind=RKIND), pointer :: Kc
      real (kind=RKIND), pointer :: alpha_c
      real (kind=RKIND), pointer :: beta_c
      real (kind=RKIND), pointer :: creep_coeff
      real (kind=RKIND), pointer :: rhoi
      real (kind=RKIND), pointer :: config_SGH_incipient_channel_width
      logical, pointer :: config_SGH_include_pressure_melt
      real (kind=RKIND), dimension(:), pointer :: channelArea
      real (kind=RKIND), dimension(:), pointer :: channelMelt
      real (kind=RKIND), dimension(:), pointer :: channelPressureFreeze
      real (kind=RKIND), dimension(:), pointer :: channelDischarge
      real (kind=RKIND), dimension(:), pointer :: channelVelocity
      real (kind=RKIND), dimension(:), pointer :: gradMagPhiEdge
      real (kind=RKIND), dimension(:), pointer :: waterFlux
      real (kind=RKIND), dimension(:), pointer :: hydropotentialSlopeNormal
      real (kind=RKIND), dimension(:), pointer :: waterPressureSlopeNormal
      real (kind=RKIND), dimension(:), pointer :: channelOpeningRate
      real (kind=RKIND), dimension(:), pointer :: channelClosingRate
      real (kind=RKIND), dimension(:), pointer :: channelChangeRate
      real (kind=RKIND), dimension(:), pointer :: flowParamAChannel
      real (kind=RKIND), dimension(:), pointer :: channelEffectivePressure
      real (kind=RKIND), dimension(:), pointer :: effectivePressureSGH
      real (kind=RKIND), dimension(:), pointer :: channelDiffusivity
      real (kind=RKIND), dimension(:), pointer :: waterThicknessEdgeUpwind
      integer, dimension(:), pointer :: waterFluxMask
      integer, dimension(:), pointer :: hydroMarineMarginMask
      integer, dimension(:), pointer :: edgeMask
      real (kind=RKIND), dimension(:,:), pointer :: flowParamA
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer, pointer :: nVertLevels
      integer, pointer :: nEdgesSolve
      integer :: iEdge, cell1, cell2

      err = 0

      ! Get pools things
      call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
      call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
      call mpas_pool_get_config(liConfigs, 'config_ice_density', rhoi)
      call mpas_pool_get_config(liConfigs, 'config_SGH_chnl_conduc_coeff', Kc)
      call mpas_pool_get_config(liConfigs, 'config_SGH_chnl_alpha', alpha_c)
      call mpas_pool_get_config(liConfigs, 'config_SGH_chnl_beta', beta_c)
      call mpas_pool_get_config(liConfigs, 'config_SGH_chnl_creep_coefficient', creep_coeff)
      call mpas_pool_get_config(liConfigs, 'config_SGH_incipient_channel_width', config_SGH_incipient_channel_width)
      call mpas_pool_get_config(liConfigs, 'config_SGH_include_pressure_melt', config_SGH_include_pressure_melt)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_array(hydroPool, 'channelArea', channelArea)
      call mpas_pool_get_array(hydroPool, 'channelMelt', channelMelt)
      call mpas_pool_get_array(hydroPool, 'channelPressureFreeze', channelPressureFreeze)
      call mpas_pool_get_array(hydroPool, 'channelDischarge', channelDischarge)
      call mpas_pool_get_array(hydroPool, 'channelVelocity', channelVelocity)
      call mpas_pool_get_array(hydroPool, 'gradMagPhiEdge', gradMagPhiEdge)
      call mpas_pool_get_array(hydroPool, 'hydropotentialSlopeNormal', hydropotentialSlopeNormal)
      call mpas_pool_get_array(hydroPool, 'waterPressureSlopeNormal', waterPressureSlopeNormal)
      call mpas_pool_get_array(hydroPool, 'waterFlux', waterFlux)
      call mpas_pool_get_array(hydroPool, 'channelOpeningRate', channelOpeningRate)
      call mpas_pool_get_array(hydroPool, 'channelClosingRate', channelClosingRate)
      call mpas_pool_get_array(hydroPool, 'channelChangeRate', channelChangeRate)
      call mpas_pool_get_array(hydroPool, 'flowParamAChannel', flowParamAChannel)
      call mpas_pool_get_array(hydroPool, 'channelEffectivePressure', channelEffectivePressure)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(velocityPool, 'flowParamA', flowParamA)
      call mpas_pool_get_array(hydroPool, 'effectivePressureSGH', effectivePressureSGH)
      call mpas_pool_get_array(hydroPool, 'waterFluxMask', waterFluxMask)
      call mpas_pool_get_array(hydroPool, 'hydroMarineMarginMask', hydroMarineMarginMask)
      call mpas_pool_get_array(hydroPool, 'channelDiffusivity', channelDiffusivity)
      call mpas_pool_get_array(hydroPool, 'waterThicknessEdgeUpwind', waterThicknessEdgeUpwind)
      call mpas_pool_get_array(geometryPool, 'edgeMask', edgeMask)
      ! Calculate terms needed for opening (melt) rate

      where(gradMagPhiEdge < SMALL_GRADPHI)
         channelDischarge(:) = 0.0_RKIND
      elsewhere
         channelDischarge = -1.0_RKIND * Kc * channelArea**alpha_c * gradMagPhiEdge**(beta_c - 2.0_RKIND) * &
            hydropotentialSlopeNormal
      end where

      where (waterFluxMask == 2)
         channelDischarge = 0.0_RKIND
         channelArea = 0.0_RKIND
      end where

      ! Note: an edge with only one grounded cell neighbor is called floating, so this logic retains channel vars
      ! on those edges to allow channel discharge across GL
      where (.not. ( (li_mask_is_grounded_ice(edgeMask)) .or. (hydroMarineMarginMask==1) ) )
         channelArea = 0.0_RKIND
         channelDischarge = 0.0_RKIND
      end where

      ! Similar to waterFlux, shut off channelDischarge if no water upstream. Prevents self-sustaining channels in the absence of distributed water
      where (waterThicknessEdgeUpwind == 0.0_RKIND)
         channelArea = 0.0_RKIND
         channelDischarge = 0.0_RKIND
      end where

      channelVelocity = channelDischarge / (channelArea + 1.0e-12_RKIND)

      ! diffusivity used only to limit channel dt right now
      where(gradMagPhiEdge < SMALL_GRADPHI)
         channelDiffusivity = 0.0_RKIND
      elsewhere
         channelDiffusivity = abs(rho_water * gravity * channelArea *  &
            Kc * channelArea**(alpha_c - 1.0_RKIND) * gradMagPhiEdge**(beta_c - 2.0_RKIND))
      end where

      channelMelt = (abs(channelDischarge * hydropotentialSlopeNormal) &  ! channel dissipation
                  +  abs(waterFlux * hydropotentialSlopeNormal * config_SGH_incipient_channel_width) & !some sheet dissipation
                  ) / latent_heat_ice
      channelPressureFreeze = -1.0_RKIND * iceMeltingPointPressureDependence * cp_freshwater * rho_water * &
         (channelDischarge + waterFlux * config_SGH_incipient_channel_width) &
         * waterPressureSlopeNormal / latent_heat_ice

      if (config_SGH_include_pressure_melt) then
         channelOpeningRate = (channelMelt - channelPressureFreeze) / rhoi
      else
         channelOpeningRate = channelMelt / rhoi
      endif

      ! Calculate terms needed for closing (creep) rate
      ! Need cell center quantities on edges
      do iEdge = 1, nEdgesSolve
         cell1 = cellsOnEdge(1, iEdge)
         cell2 = cellsOnEdge(2, iEdge)

         ! Not sure if these ought to be upwind average, but using centered
         flowParamAChannel(iEdge) = 0.5_RKIND * (flowParamA(nVertLevels, cell1) + flowParamA(nVertLevels, cell2) )
         channelEffectivePressure(iEdge) = 0.5_RKIND * (effectivePressureSGH(cell1) + effectivePressureSGH(cell2))
      end do
      channelClosingRate(:) = creep_coeff * channelArea(:) * flowParamAChannel(:) * channelEffectivePressure(:)**3

      where (waterFluxMask == 2)
         channelOpeningRate = 0.0_RKIND
         channelClosingRate = 0.0_RKIND
      end where

      channelChangeRate = channelOpeningRate - channelClosingRate

   !--------------------------------------------------------------------
   end subroutine update_channel


!***********************************************************************
!
!  routine evolve_channel
!
!> \brief   Evolve SGH channel area, calculate changes on cells
!> \author  Matt Hoffman
!> \date    28 July 2016
!> \details
!>  This routine updates the channel area in the subglacial hydrology model.
!>  It uses the conduit space evolution equation.
!-----------------------------------------------------------------------
   subroutine evolve_channel(block, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (block_type), intent(inout) :: block    !< Input/Output: block object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      ! Pools pointers
      type (mpas_pool_type), pointer :: hydroPool
      type (mpas_pool_type), pointer :: meshPool
      real (kind=RKIND), dimension(:), pointer :: channelArea
      real (kind=RKIND), dimension(:), pointer :: channelDischarge
      real (kind=RKIND), dimension(:), pointer :: channelMelt, channelPressureFreeze
      real (kind=RKIND), dimension(:), pointer :: channelChangeRate
      real (kind=RKIND), dimension(:), pointer :: divergenceChannel
      real (kind=RKIND), dimension(:), pointer :: channelAreaChangeCell
      real (kind=RKIND), dimension(:), pointer :: channelMeltInputCell
      real (kind=RKIND), pointer :: deltatSGH
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer, dimension(:,:), pointer :: edgeSignOnCell
      real (kind=RKIND), dimension(:), pointer :: areaCell
      real (kind=RKIND), dimension(:), pointer :: dcEdge
      integer, pointer :: nCellsSolve
      integer :: iCell, iEdgeOnCell, iEdge

      call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_array(hydroPool, 'deltatSGH', deltatSGH)
      call mpas_pool_get_array(hydroPool, 'channelArea', channelArea)
      call mpas_pool_get_array(hydroPool, 'channelDischarge', channelDischarge)
      call mpas_pool_get_array(hydroPool, 'channelMelt', channelMelt)
      call mpas_pool_get_array(hydroPool, 'channelPressureFreeze', channelPressureFreeze)
      call mpas_pool_get_array(hydroPool, 'channelChangeRate', channelChangeRate)
      call mpas_pool_get_array(hydroPool, 'divergenceChannel', divergenceChannel)
      call mpas_pool_get_array(hydroPool, 'channelAreaChangeCell', channelAreaChangeCell)
      call mpas_pool_get_array(hydroPool, 'channelMeltInputCell', channelMeltInputCell)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)

      ! Calculate flux divergence on cells and channel area change on cells
      divergenceChannel(:) = 0.0_RKIND  ! zero div before accumulating
      channelAreaChangeCell(:) = 0.0_RKIND  ! zero before accumulating
      channelMeltInputCell(:) = 0.0_RKIND  ! zero before accumulating
      ! loop over locally owned cells
      do iCell = 1, nCellsSolve
         ! TODO: could limit to grounded cells only
         ! compute fluxes for each edge of the cell
         do iEdgeOnCell = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(iEdgeOnCell, iCell)
            ! add on advective & diffusive fluxes
            divergenceChannel(iCell) = divergenceChannel(iCell) - channelDischarge(iEdge) * edgeSignOnCell(iEdgeOnCell, iCell)
            channelAreaChangeCell(iCell) = channelChangeRate(iEdge) * dcEdge(iEdge) * 0.5_RKIND + channelAreaChangeCell(iCell)
               ! < only half of channel is in this cell
            channelMeltInputCell(iCell) = 0.5_RKIND * (channelMelt(iEdge) - channelPressureFreeze(iEdge)) * dcEdge(iEdge) / rho_water + channelMeltInputCell(iCell)
         end do ! edges

      end do ! cells
      divergenceChannel(1:nCellsSolve) = divergenceChannel(1:nCellsSolve) / areaCell(1:nCellsSolve)
      channelAreaChangeCell(1:nCellsSolve) = channelAreaChangeCell(1:nCellsSolve) / areaCell(1:nCellsSolve)
      channelMeltInputCell(1:nCellsSolve) = channelMeltInputCell(1:nCellsSolve) / areaCell(1:nCellsSolve)

      ! Now update channel area
      channelArea = channelChangeRate * deltatSGH + channelArea
      ! If sheet dissipation contributes to channel, there should be no need for a minimum channel size
      channelArea = max(1.0e-8_RKIND, channelArea)  ! make some tiny value when it goes negative

   !--------------------------------------------------------------------
   end subroutine evolve_channel


!***********************************************************************
!
!  routine shmip_timevarying_forcing
!
!> \brief   Calculate time-varying forcings for SHMIP experiments
!> \author  Matt Hoffman
!> \date    18 January 2017
!> \details
!>  This routine calculates time-varying forcings for the Subglacial
!>  Hydrology Model Intercomparison Project's Experiments C and D.
!>  I chose to include these in model source code rather than through
!>  a forcing module because 1) this allows the SHMIP analytic forcing
!>  functions to be reproduced exactly and 2) with the short hydro model
!>  time steps, using a forcing i/o module would require a lot of file
!>  reading during each run.
!-----------------------------------------------------------------------
   subroutine shmip_timevarying_forcing(block, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (block_type), intent(inout) :: block    !< Input/Output: block object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      ! Pools pointers
      type (mpas_pool_type), pointer :: hydroPool
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: meshPool
      character (len=StrKIND), pointer :: config_SGH_shmip_forcing
      real(kind=RKIND), pointer :: daysSinceStart
      real (kind=RKIND), dimension(:), pointer :: externalWaterInput
      real (kind=RKIND), dimension(:), pointer :: basalMeltInput
      real (kind=RKIND), dimension(:), pointer :: upperSurface
      real (kind=RKIND), dimension(:), pointer :: areaCell
      ! Test case specific variables
      real(kind=RKIND) :: ra ! relative amplitude for versions of C experiment
      real(kind=RKIND) :: temperature0  ! temperature at 0m elevation for D experiment
      real(kind=RKIND) :: DT ! temperature offset for versions of D experiment

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_SGH_shmip_forcing', config_SGH_shmip_forcing)

      if (config_SGH_shmip_forcing=='none') then
         ! Default is to do nothing and skip rest of routine
         return
      endif

      call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
      call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_array(meshPool, 'daysSinceStart', daysSinceStart)
      call mpas_pool_get_array(hydroPool, 'externalWaterInput', externalWaterInput)

      select case (config_SGH_shmip_forcing(1:1))  ! Check on first character only
      case ('C')
         ! Note: Case C assumes that moulin locations are input via externalWaterInput
         ! and have a very small value that will not affect the sum of
         ! externalWaterInput+basalMeltInput.
         ! The actual forcing is stored in basalMeltInput.  Note this is the moulin
         ! and the basal contributions combined, unlike in the other tests
         ! where I keep those separate.

         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_array(hydroPool, 'basalMeltInput', basalMeltInput)

         select case (config_SGH_shmip_forcing(2:2)) ! get experiment number
         case ('1')
            ra = 0.25_RKIND
         case ('2')
            ra = 0.5_RKIND
         case ('3')
            ra = 1.0_RKIND
         case ('4')
            ra = 2.0_RKIND
         case default
            call mpas_log_write("Unknown value for config_SGH_shmip_forcing:" // config_SGH_shmip_forcing, MPAS_LOG_ERR)
            err = ior(err, 1)
         end select

         ! basalMeltInput calculated in SHMIP m/s units here
         where (externalWaterInput > 0.0_RKIND)
            ! 0.9/areaCell is the value being input into each moulin in B5
            basalMeltInput = 0.9_RKIND / areaCell * (1.0_RKIND - ra * sin(2.0_RKIND * pii * daysSinceStart))
         elsewhere
            ! where the input file does not identify a moulin, set source term to 0
            basalMeltInput = 0.0_RKIND
         end where
         ! instructions say to zero negative values
         where (basalMeltInput < 0.0_RKIND)
            basalMeltInput = 0.0_RKIND
         end where
         ! Add in basal melting contribution specified by instructions AND convert from m/s to kg/m2/s
         basalMeltInput = (basalMeltInput + 7.93e-11_RKIND) * 1000.0_RKIND

      case ('D')
         select case (config_SGH_shmip_forcing(2:2)) ! get experiment number
         case ('1')
            DT = -4.0_RKIND
         case ('2')
            DT = -2.0_RKIND
         case ('3')
            DT = 0.0_RKIND
         case ('4')
            DT = 2.0_RKIND
         case ('5')
            DT = 4.0_RKIND
         case default
            call mpas_log_write("Unknown value for config_SGH_shmip_forcing:" // config_SGH_shmip_forcing, MPAS_LOG_ERR)
            err = ior(err, 1)
         end select

         call mpas_pool_get_array(geometryPool, 'upperSurface', upperSurface)
         temperature0 = -16.0_RKIND * cos(2.0_RKIND * pii / 365.0 * daysSinceStart) - 5.0_RKIND + DT
         externalWaterInput(:) = (upperSurface(:) * (-0.0075_RKIND) + temperature0) * (0.01_RKIND / 86400.0_RKIND) * 1000.0_RKIND
           ! < the 1000 factor converts from m/s to kg/m2/s
         where (externalWaterInput < 0.0_RKIND)
            externalWaterInput = 0.0_RKIND
         end where

      case default
         call mpas_log_write("Unknown value for config_SGH_shmip_forcing:" // config_SGH_shmip_forcing, MPAS_LOG_ERR)
         err = ior(err, 1)
      end select


   !--------------------------------------------------------------------
   end subroutine shmip_timevarying_forcing



!***********************************************************************
!
!  routine ocean_connection_N
!
!> \brief   Calculate effective pressure assuming perfect ocean connection
!> \author  Matt Hoffman
!> \date    02 June 2020
!> \details
!>  This routine calculates effective pressure assuming a perfect hydrologic
!>  connection with the ocean.  It can be used as a simple alternative to the
!>  full subglacial hydrology.
!-----------------------------------------------------------------------
   subroutine ocean_connection_N(domain)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: hydroPool
      type (mpas_pool_type), pointer :: geometryPool
      real (kind=RKIND), dimension(:), pointer :: bedTopography
      real (kind=RKIND), dimension(:), pointer :: effectivePressure
      real (kind=RKIND), dimension(:), pointer :: thickness
      real (kind=RKIND), pointer :: rhoi, rhoo

      ! Calculate N assuming perfect ocean connection
      call mpas_log_write('Calculating N assuming perfect ocean connection.')

      call mpas_pool_get_config(liConfigs, 'config_ice_density', rhoi)
      call mpas_pool_get_config(liConfigs, 'config_ocean_density', rhoo)
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
         call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
         call mpas_pool_get_array(hydroPool, 'effectivePressure', effectivePressure)
         call mpas_pool_get_array(hydroPool, 'thickness', thickness)
         
         effectivePressure = rhoi * gravity * thickness - rhoi * gravity * max(0.0_RKIND,  -1.0_RKIND * rhoo/rhoi * bedTopography)
         effectivePressure = max(effectivePressure, 0.0_RKIND) ! This is just to zero out N in the open ocean to avoid confusion

         block => block % next
      end do

   !--------------------------------------------------------------------
   end subroutine ocean_connection_N


!***********************************************************************
!
!  routine calc_hydro_mask
!
!> \brief   Calculate the boundaries of the active hydrology domain
!> \author  Matt Hoffman
!> \date    24 October 2022
!> \details
!>  This routine calculates a mask of the boundaries of the active hydrology domain.
!>  If there is no waterFluxMask set around domain boundaries, then calc_hydro_mask creates one.
!-----------------------------------------------------------------------
   subroutine calc_hydro_mask(domain)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: hydroPool
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: meshPool
      real (kind=RKIND), dimension(:), pointer :: bedTopography
      integer, dimension(:), pointer :: hydroMarineMarginMask
      integer, dimension(:), pointer :: hydroTerrestrialMarginMask
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer, dimension(:), pointer :: cellMask
      integer, dimension(:), pointer :: waterFluxMask
      integer, pointer :: nEdgesSolve
      integer, pointer :: nCells
      integer :: cell1, cell2, iEdge
      real (kind=RKIND), pointer :: config_sea_level
      integer :: wfmWarning

      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_array(hydroPool, 'hydroMarineMarginMask', hydroMarineMarginMask)
         call mpas_pool_get_array(hydroPool, 'hydroTerrestrialMarginMask', hydroTerrestrialMarginMask)
         call mpas_pool_get_array(hydroPool, 'waterFluxMask', waterFluxMask)
         call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
         call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

         hydroMarineMarginMask(:) = 0
         hydroTerrestrialMarginMask(:) = 0
         wfmWarning = 0

         do iEdge = 1, nEdgesSolve
            cell1 = cellsOnEdge(1, iEdge)
            cell2 = cellsOnEdge(2, iEdge)

            !ensure no-flow conditions on edges of domain if not set up
            !in input file
            if ( ((cell1 == nCells + 1) .or. (cell2 == nCells +1 )) &
               .and. (waterFluxMask(iEdge) .ne. 2)) then
                 waterFluxMask(iEdge) = 2
                 wfmWarning = 1
            endif

            ! hydroMarineMarginMask: We are looking for edges with 1 edge grounded ice and the other edge floating ice or open ocean
            ! Exclude no-flow boundaries
            if ( (li_mask_is_grounded_ice(cellMask(cell1))) .and. &
                 (li_mask_is_floating_ice(cellMask(cell2)) .or. &
                 ((bedTopography(cell2) < config_sea_level) .and. (.not. li_mask_is_ice(cellMask(cell2)))) ) &
                  .and. (waterFluxMask(iEdge) .ne. 2) ) then
               hydroMarineMarginMask(iEdge) = 1
            elseif ( (li_mask_is_grounded_ice(cellMask(cell2))) .and. &
                     (li_mask_is_floating_ice(cellMask(cell1)) .or. &
                     ((bedTopography(cell1) < config_sea_level) .and. (.not. li_mask_is_ice(cellMask(cell1)))) ) &
                     .and. (waterFluxMask(iEdge) .ne. 2) ) then
               hydroMarineMarginMask(iEdge) = 1
            endif

            ! hydroTerrestrialMarginMask: Look for edges with 1 cell on grounding ice and the other cell on land without ice
            ! Exclude no-flow boundaries
            if ((li_mask_is_grounded_ice(cellMask(cell1))) .and. ( .not. li_mask_is_ice(cellMask(cell2))) &
               .and. (bedTopography(cell2) >= config_sea_level) &
               .and. (waterFluxMask(iEdge) .ne. 2) ) then
               hydroTerrestrialMarginMask(iEdge) = 1
            elseif ((li_mask_is_grounded_ice(cellMask(cell2))) .and. ( .not. li_mask_is_ice(cellMask(cell1))) &
               .and. (bedTopography(cell1) >= config_sea_level) &
               .and. (waterFluxMask(iEdge) .ne. 2) ) then
               hydroTerrestrialMarginMask(iEdge) = 1
            endif
         enddo

         block => block % next
      end do

      if (wfmWarning == 1) then
         call mpas_log_write('Changing waterFluxMask to enforce no-flow conditions at domain boundaries', MPAS_LOG_WARN)
      endif
      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'hydroMarineMarginMask')
      call mpas_dmpar_field_halo_exch(domain, 'hydroTerrestrialMarginMask')
      call mpas_dmpar_field_halo_exch(domain, 'waterFluxMask')
      call mpas_timer_stop("halo updates")

   !--------------------------------------------------------------------
   end subroutine calc_hydro_mask
   

!***********************************************************************
!
!  routine calc_iceThicknessHydro
!
!> \brief   Calculate version of ice thickness used by hydrology model
!> \author  Alex Hager
!> \date    7 June 2023
!> \details
!> This routine calculates a modified ice thickness that is altered to 
!> avoid local minima in hydropotential
!-----------------------------------------------------------------------
   subroutine calc_iceThicknessHydro(block, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (block_type), intent(inout) :: block    !< Input/Output: block object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      ! Pools pointers
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: hydroPool
      type (mpas_pool_type), pointer :: meshPool
      real (kind=RKIND), dimension(:), pointer :: thickness
      real (kind=RKIND), dimension(:), pointer :: bedTopography
      real (kind=RKIND), dimension(:), pointer :: iceThicknessHydro
      real (kind=RKIND), dimension(:), pointer :: upperSurface
      integer, dimension(:,:), pointer :: cellsOnCell
      integer, dimension(:), pointer :: cellMask
      real (kind=RKIND), dimension(:), pointer :: areaCell
      integer, pointer :: nCells
      integer, dimension(:), pointer :: nEdgesOnCell
      integer :: iCell
      integer :: jCell
      integer :: iNeighbor
      real (kind=RKIND) :: maxNeighborHeight
      real (kind=RKIND) :: minNeighborHeight
      real (kind=RKIND) :: meanNeighborHeight
      real (kind=RKIND) :: totalNeighborHeight
      real (kind=RKIND) :: totalArea
      real (kind=RKIND), parameter :: bigValue = 1.0e6_RKIND
      real (kind=RKIND), parameter :: bigNegativeValue = -1.0e6_RKIND
      logical, pointer :: config_SGH_use_iceThicknessHydro
      err = 0

      ! Get pools things
      call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
      call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_array(hydroPool, 'iceThicknessHydro', iceThicknessHydro)
      call mpas_pool_get_array(geometryPool, 'upperSurface', upperSurface)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_config(liConfigs, 'config_SGH_use_iceThicknessHydro', config_SGH_use_iceThicknessHydro)

      iceThicknessHydro = thickness 

      !iceThicknessHydro is equivalent to thickness if
      !config_SGH_use_iceThicknessHydro is false. Default is to alter
      !ice thickness

      if (config_SGH_use_iceThicknessHydro) then
         do iCell = 1, nCells

            if (li_mask_is_grounded_ice(cellMask(iCell))) then !identify grounded ice
               
                  maxNeighborHeight = bigNegativeValue !initialize
                  minNeighborHeight = bigValue
                  
                  totalNeighborHeight = 0.0_RKIND
                  totalArea = 0.0_RKIND

                  do jCell = 1, nEdgesOnCell(iCell)
                  
                     iNeighbor = cellsOnCell(jCell,iCell)

                     !Only include neighbor cell in averaging if it contains grounded ice
                     if ((li_mask_is_grounded_ice(cellMask(iNeighbor)))) then

                         minNeighborHeight = min(minNeighborHeight, upperSurface(iNeighbor))

                         maxNeighborHeight = max(maxNeighborHeight, upperSurface(iNeighbor))

                         totalNeighborHeight = totalNeighborHeight + upperSurface(iNeighbor)*areaCell(iNeighbor)

                         totalArea = totalArea + areaCell(iNeighbor)
                         
                     endif
                  end do
                  
                  if ((totalNeighborHeight == 0.0_RKIND) .or. (totalArea == 0.0_RKIND)) then
                     meanNeighborHeight = upperSurface(iCell) !no smoothing in single-cell islands
                  else
                     meanNeighborHeight = totalNeighborHeight / totalArea !area-weighted average height
                  endif

                  !only adjust surface elevation if cell is local minimum or
                  !maximum
                  if ((upperSurface(iCell) < minNeighborHeight) &
                          .and. (minNeighborHeight < bigValue)) then

                     iceThicknessHydro(iCell) = meanNeighborHeight - bedTopography(iCell)

                  elseif ((upperSurface(iCell) > maxNeighborHeight) &
                          .and. (maxNeighborHeight > bigNegativeValue)) then

                     iceThicknessHydro(iCell) = meanNeighborHeight - bedTopography(iCell)

                  endif

                  !avoid negatives, default to original thickness
                  if (iceThicknessHydro(iCell) < 0.0_RKIND) then
                     iceThicknessHydro(iCell) = thickness(iCell)
                  endif
            endif
         enddo
      endif

   end subroutine calc_iceThicknessHydro

!***********************************************************************
!
!  routine calc_gl_total
!
!> \brief   Calculate total grounding line discharge on
!           adjacent ocean cell
!> \author  Alex Hager
!> \date    27 March 2024
!> \details Find the total amount of freshwater entering the first ocean cell from the grounding line.
!-----------------------------------------------------------------------
   subroutine calc_gl_totals(block, err)
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (block_type), intent(inout) :: block    !< Input/Output: block object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (mpas_pool_type), pointer :: hydroPool
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: meshPool

      real (kind=RKIND), dimension(:), pointer :: distGroundingLineDischargeCell
      real (kind=RKIND), dimension(:), pointer :: chnlGroundingLineDischargeCell
      real (kind=RKIND), dimension(:), pointer :: totalGroundingLineDischargeCell
      real (kind=RKIND), dimension(:), pointer :: channelDischarge
      real (kind=RKIND), dimension(:), pointer :: waterFlux
      real (kind=RKIND), dimension(:), pointer :: dvEdge
      integer iEdge, cell1, cell2
      integer, pointer :: nEdgesSolve
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer, dimension(:), pointer :: hydroMarineMarginMask
      integer, dimension(:), pointer :: cellMask
      real (kind=RKIND) :: distGroundingLineDischargeEdge
      real (kind=RKIND) :: chnlGroundingLineDischargeEdge
      real (kind=RKIND) :: totalGroundingLineDischargeEdge

      call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
      call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

      call mpas_pool_get_array(hydroPool, 'distGroundingLineDischargeCell', distGroundingLineDischargeCell)
      call mpas_pool_get_array(hydroPool, 'chnlGroundingLineDischargeCell', chnlGroundingLineDischargeCell)
      call mpas_pool_get_array(hydroPool, 'totalGroundingLineDischargeCell', totalGroundingLineDischargeCell)
      call mpas_pool_get_array(hydroPool, 'channelDischarge', channelDischarge)
      call mpas_pool_get_array(hydroPool, 'waterFlux', waterFlux)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(hydroPool, 'hydroMarineMarginMask', hydroMarineMarginMask)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)

      distGroundingLineDischargeCell(:) = 0.0_RKIND
      chnlGroundingLineDischargeCell(:) = 0.0_RKIND
      totalGroundingLineDischargeCell(:) = 0.0_RKIND

      do iEdge = 1, nEdgesSolve

         if (hydroMarineMarginMask(iEdge) == 1) then

            cell1 = cellsOnEdge(1, iEdge)
            cell2 = cellsOnEdge(2, iEdge)

            !calculate totals at each grounding line edge
            distGroundingLineDischargeEdge = abs(waterFlux(iEdge) * dvEdge(iEdge))
            chnlGroundingLineDischargeEdge = abs(channelDischarge(iEdge))
            totalGroundingLineDischargeEdge = distGroundingLineDischargeEdge + chnlGroundingLineDischargeEdge

            !Assign sum of edge totals to adjacent ungrounded cell
            if (li_mask_is_grounded_ice(cellMask(cell1))) then
               distGroundingLineDischargeCell(cell2) = distGroundingLineDischargeCell(cell2) + distGroundingLineDischargeEdge
               chnlGroundingLineDischargeCell(cell2) = chnlGroundingLineDischargeCell(cell2) + chnlGroundingLineDischargeEdge
               totalGroundingLineDischargeCell(cell2) = totalGroundingLineDischargeCell(cell2) + totalGroundingLineDischargeEdge
            elseif (li_mask_is_grounded_ice(cellMask(cell2))) then
               distGroundingLineDischargeCell(cell1) = distGroundingLineDischargeCell(cell1) + distGroundingLineDischargeEdge
               chnlGroundingLineDischargeCell(cell1) = chnlGroundingLineDischargeCell(cell1) + chnlGroundingLineDischargeEdge
               totalGroundingLineDischargeCell(cell1) = totalGroundingLineDischargeCell(cell1) + totalGroundingLineDischargeEdge
            endif
         endif
      enddo
   end subroutine calc_gl_totals

!***********************************************************************
!
!  routine calc_waterPressureSmooth
!
!> \brief   Calculate smoothed version of waterPressure used for calculation
!>          of waterPressureSlopeNormal
!> \author  Alex Hager
!> \date    29 February 2024
!> \details Creates a smoothed version of waterPressure, called waterPressureSmooth,
!>          that is used for calculation of waterPressureSlopeNormal.
!>          This is necessary to increase stability with channelPressureFreeze
!>          with spatially variable bedTopography and upperSurface. waterPressureSmooth 
!>          an area-weighted average of the current cells and its neighbors. Possible to
!>          perform multiple iterations of smoothing by adjusting
!>          config_SGH_iter_smooth_waterPressureSlopeNormal
!-----------------------------------------------------------------------
   subroutine calc_waterPressureSmooth(block,err)
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (block_type), intent(inout) :: block    !< Input/Output: block object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: hydroPool
      real (kind=RKIND), dimension(:), pointer :: bedTopography
      real (kind=RKIND), dimension(:), pointer :: waterPressure
      real (kind=RKIND), dimension(:), pointer :: waterPressureSmooth
      real (kind=RKIND), dimension(:), pointer :: areaCell
      real (kind=RKIND), dimension(:), allocatable :: pressureField
      real(kind=RKIND) :: totalPressure
      real(kind=RKIND) :: totalArea
      integer, dimension(:), pointer :: cellMask
      integer, dimension(:,:), pointer :: cellsOnCell
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, pointer :: nTimesSmooth
      integer, pointer :: nCells
      integer :: iCell, jEdge
      integer :: iNeighbor 
      integer :: iter

      err = 0

      ! Get pools things
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
      call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_config(liConfigs, 'config_SGH_iter_smooth_waterPressureSlopeNormal', nTimesSmooth)
      call mpas_pool_get_array(hydroPool, 'waterPressure', waterPressure)
      call mpas_pool_get_array(hydroPool, 'waterPressureSmooth', waterPressureSmooth)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      allocate(pressureField(nCells+1))

      waterPressureSmooth = waterPressure
      pressureField = waterPressure

      do iter = 1, nTimesSmooth !May need to iterate to smooth properly

         do iCell = 1, nCells
            if (li_mask_is_grounded_ice(cellMask(iCell))) then !smooth over all grounded cells

               totalPressure = pressureField(iCell)*areaCell(iCell) 
               totalArea = areaCell(iCell)
          
               do jEdge = 1, nEdgesOnCell(iCell)
                  iNeighbor = cellsOnCell(jEdge,iCell)

                  if ((li_mask_is_grounded_ice(cellMask(iNeighbor)))) then !only include GROUNDED neighboring cells in smoothing 

                     totalPressure = totalPressure + pressureField(iNeighbor)*areaCell(iNeighbor)

                     totalArea = totalArea + areaCell(iNeighbor)
                  endif
               end do
          
               waterPressureSmooth(iCell) = totalPressure / totalArea !area-weighted average pressure
            endif
         end do
         pressureField = waterPressureSmooth
      end do
      
      deallocate(pressureField)
   end subroutine calc_waterPressureSmooth

end module li_subglacial_hydro
