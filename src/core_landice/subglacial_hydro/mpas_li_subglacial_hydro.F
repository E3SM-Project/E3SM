! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
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
   use mpas_constants, only: gravity
   use mpas_timer
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
      logical, pointer :: config_SGH_adaptive_timestep
      logical, pointer :: config_adaptive_timestep
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: hydroPool
      type (mpas_pool_type), pointer :: geometryPool
      real (kind=RKIND), pointer :: deltat
      real (kind=RKIND), pointer :: deltatSGH
      real (kind=RKIND), dimension(:), pointer :: waterThickness
      real (kind=RKIND), dimension(:), pointer :: tillWaterThickness
      real (kind=RKIND), dimension(:), pointer :: waterPressure
      real (kind=RKIND), dimension(:), pointer :: thickness
      real (kind=RKIND), pointer :: tillMax
      real (kind=RKIND), pointer :: rhoi
      integer, pointer :: config_num_halos
      integer :: err_tmp


      err = 0
      err_tmp = 0

      call mpas_pool_get_config(liConfigs, 'config_SGH', config_SGH)
      if (.not. config_SGH) then
         ! If SGH is not active, skip everything
         return
      endif

      call mpas_timer_start("hydro init")

      write(stdoutUnit,*) 'Beginning subglacial hydro init.'

      ! Check number of halos
      call mpas_pool_get_config(liConfigs, 'config_num_halos', config_num_halos)
      if (config_num_halos < 1) then
         write (stdErrUnit,*) "ERROR: Subglacial hydrology requires config_num_halos >= 1"
         err = ior(err, 1)
      endif

      ! Check for valid time steppers
      call mpas_pool_get_config(liConfigs, 'config_SGH_adaptive_timestep', config_SGH_adaptive_timestep)
      call mpas_pool_get_config(liConfigs, 'config_adaptive_timestep', config_adaptive_timestep)
      if (config_adaptive_timestep .and. config_SGH_adaptive_timestep) then
         write (stdErrUnit,*) "ERROR: Subglacial hydrology cannot use its own adaptive time stepper if primary time stepper is adaptive.  This may be relaxed in the future."
         err = ior(err, 1)
      endif

      call mpas_pool_get_config(liConfigs, 'config_SGH_till_max', tillMax)
      call mpas_pool_get_config(liConfigs, 'config_ice_density', rhoi)

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)

         call mpas_pool_get_array(meshPool, 'deltat', deltat)
         call mpas_pool_get_array(hydroPool, 'deltatSGH', deltatSGH)

         if (config_SGH_adaptive_timestep) then
            deltatSGH = 1.0_RKIND  ! Until init is done properly, make this tiny
         else
            ! For now, force hydro dt to be master dt
            deltatSGH = deltat
         endif

         ! remove invalid values - not necessary on restart, but shouldn't hurt
         call mpas_pool_get_array(hydroPool, 'waterThickness', waterThickness)
         waterThickness = max(0.0_RKIND, waterThickness)

         call mpas_pool_get_array(hydroPool, 'tillWaterThickness', tillWaterThickness)
         tillWaterThickness = max(0.0_RKIND, tillWaterThickness)
         tillWaterThickness = min(tillMax, tillWaterThickness)

         call mpas_pool_get_array(hydroPool, 'waterPressure', waterPressure)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         waterPressure = max(0.0_RKIND, waterPressure)
         waterPressure = min(waterPressure, rhoi * gravity * thickness)

         ! Initialize diagnostic pressure variables
         call calc_pressure_diag_vars(block, err_tmp)
         err = ior(err, err_tmp)

         block => block % next
      end do

      ! === error check
      if (err > 0) then
          write (stderrUnit,*) "An error has occurred in li_SGH_init."
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
      logical, pointer :: config_SGH_chnl_active
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: hydroPool
      type (mpas_pool_type), pointer :: meshPool
      real (kind=RKIND), dimension(:), pointer :: Wtill, WtillOld
      real (kind=RKIND), dimension(:), pointer :: basalMeltInput
      real (kind=RKIND), dimension(:), pointer :: externalWaterInput
      real (kind=RKIND), dimension(:), pointer :: divergence
      real (kind=RKIND), dimension(:), pointer :: waterFlux
      real (kind=RKIND), dimension(:), pointer :: waterThickness
      real (kind=RKIND), dimension(:), pointer :: waterThicknessOld
      real (kind=RKIND), dimension(:), pointer :: waterThicknessTendency
      real (kind=RKIND), dimension(:), pointer :: divergenceChannel
      real (kind=RKIND), dimension(:), pointer :: channelAreaChangeCell
      real (kind=RKIND), dimension(:), pointer :: dvEdge
      real (kind=RKIND), dimension(:), pointer :: areaCell
      real (kind=RKIND), dimension(:), pointer :: waterVelocity
      real (kind=RKIND), dimension(:), pointer :: waterVelocityCellX
      real (kind=RKIND), dimension(:), pointer :: waterVelocityCellY
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
      integer :: numSubCycles ! number of subcycles
      integer :: err_tmp



      err = 0
      err_tmp = 0

      call mpas_pool_get_config(liConfigs, 'config_SGH', config_SGH)
      if (.not. config_SGH) then
         ! If SGH is not active, skip everything
         return
      endif

      write(stdoutUnit,*) 'Beginning subglacial hydro solve.'
      call mpas_pool_get_config(liConfigs, 'config_SGH_chnl_active', config_SGH_chnl_active)
      call mpas_pool_get_config(liConfigs, 'config_SGH_till_drainage', Cd)
      call mpas_pool_get_config(liConfigs, 'config_SGH_till_max', tillMax)



      ! initialize while loop
      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)  ! can get from any block
      call mpas_pool_get_array(meshPool, 'deltat', masterDeltat)
      timeLeft = masterDeltat  ! in seconds
      numSubCycles = 0
      ! =============
      ! =============
      ! =============
      ! subcycle hydro model until master dt is reached
      ! =============
      ! =============
      ! =============
      do while (timeLeft > 0.0_RKIND)
      numSubCycles = numSubCycles + 1


      ! =============
      ! Calculate SHMIP time-varying forcing, if requested
      ! =============
      block => domain % blocklist
      do while (associated(block))
         call shmip_timevarying_forcing(block, err_tmp)
         err = ior(err, err_tmp)

         block => block % next
      end do


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
         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_array(hydroPool, 'waterVelocity', waterVelocity)
         call mpas_pool_get_array(hydroPool, 'waterVelocityCellX', waterVelocityCellX)
         call mpas_pool_get_array(hydroPool, 'waterVelocityCellY', waterVelocityCellY)
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
         call mpas_timer_stop("halo updates")
      endif


      ! =============
      ! Calculate time step if being adaptive, otherwise check for CFL
      ! =============
      call check_timestep(domain, timeLeft, err_tmp)
      err = ior(err, err_tmp)


      ! =============
      ! Compute flux divergence
      ! =============
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
         call mpas_pool_get_array(hydroPool, 'divergence', divergence)
         call mpas_pool_get_array(hydroPool, 'waterFlux', waterFlux)
         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
         call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
         call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
         call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)

         divergence(:) = 0.0_RKIND  ! zero div before starting
         ! loop over locally owned cells
         do iCell = 1, nCellsSolve
            ! compute fluxes for each edge of the cell
            do iEdgeOnCell = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(iEdgeOnCell, iCell)
               ! add on advective & diffusive fluxes
               divergence(iCell) = divergence(iCell) - waterFlux(iEdge) * dvEdge(iEdge) * edgeSignOnCell(iEdgeOnCell, iCell)
            end do ! edges
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
         call mpas_dmpar_field_halo_exch(domain, 'channelArea')
         call mpas_timer_stop("halo updates")
      endif


      ! =============
      ! Calculate pressure field
      ! =============
      block => domain % blocklist
      do while (associated(block))

         call calc_pressure(block, err_tmp)
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
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)

         waterThicknessOld = waterThickness
         waterThickness = waterThicknessOld + deltatSGH * ( (basalMeltInput + externalWaterInput) / rho_water - divergence  &
             - divergenceChannel - channelAreaChangeCell  &
             - (Wtill - WtillOld) / deltatSGH)
         waterThickness = waterThickness * li_mask_is_grounded_ice_int(cellMask)  ! zero in non-grounded locations
         waterThickness = max(0.0_RKIND, waterThickness)
         divergence = divergence * li_mask_is_grounded_ice_int(cellMask)  ! zero in non-grounded locations for more convenient viz
         waterThicknessTendency = (waterThickness - waterThicknessOld) / deltatSGH

         block => block % next
      end do


      ! =============
      ! Update model time
      ! =============


      ! =============
      ! =============
      ! =============
      end do  ! while timeLeft>0
      ! =============
      ! =============
      ! =============

      write(stdErrUnit,*) "Hydro model subcycled n times:", numSubCycles


      ! === error check
      if (err > 0) then
          write (stderrUnit,*) "An error has occurred in li_SGH_solve."
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
          write (stderrUnit,*) "An error has occurred in li_SGH_finalize."
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
      real (kind=RKIND), dimension(:), pointer :: hydropotentialBase
      real (kind=RKIND), dimension(:), pointer :: hydropotentialBaseVertex
      real (kind=RKIND), dimension(:), pointer :: waterPressure
      real (kind=RKIND), dimension(:), pointer :: waterThicknessEdge
      real (kind=RKIND), dimension(:), pointer :: waterThicknessEdgeUpwind
      real (kind=RKIND), dimension(:), pointer :: waterThickness
      real (kind=RKIND), dimension(:), pointer :: hydropotentialBaseSlopeNormal
      real (kind=RKIND), dimension(:), pointer :: waterPressureSlopeNormal
      real (kind=RKIND), dimension(:), pointer :: hydropotentialBaseSlopeTangent
      real (kind=RKIND), dimension(:), pointer :: gradMagPhiEdge
      real (kind=RKIND), dimension(:), pointer :: effectiveConducEdge
      real (kind=RKIND), dimension(:), pointer :: diffusivity
      real (kind=RKIND), dimension(:), pointer :: dcEdge
      real (kind=RKIND), dimension(:), pointer :: dvEdge
      real (kind=RKIND), dimension(:), pointer :: waterVelocity
      real (kind=RKIND), dimension(:), pointer :: waterFlux
      real (kind=RKIND), dimension(:), pointer :: waterFluxAdvec
      real (kind=RKIND), dimension(:), pointer :: waterFluxDiffu
      integer, dimension(:), pointer :: waterFluxMask
      integer, dimension(:,:), pointer :: edgeSignOnCell
      integer, dimension(:), pointer :: cellMask
      integer, dimension(:), pointer :: edgeMask
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer, dimension(:,:), pointer :: verticesOnEdge
      integer, dimension(:,:), pointer :: baryCellsOnVertex
      real (kind=RKIND), dimension(:,:), pointer :: baryWeightsOnVertex
      real (kind=RKIND), pointer :: alpha, beta
      real (kind=RKIND), pointer :: conduc_coeff
      character (len=StrKIND), pointer :: config_SGH_tangent_slope_calculation
      integer, pointer :: nEdgesSolve
      integer, pointer :: nCells
      integer, pointer :: nVertices
      integer :: iEdge, cell1, cell2
      real (kind=RKIND) :: velSign
      integer :: numGroundedCells
      integer :: err_tmp

      err = 0
      err_tmp = 0


      ! Get pools things
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
      call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)

      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)

      call mpas_pool_get_config(liConfigs, 'config_SGH_alpha', alpha)
      call mpas_pool_get_config(liConfigs, 'config_SGH_beta', beta)
      call mpas_pool_get_config(liConfigs, 'config_SGH_conduc_coeff', conduc_coeff)
      call mpas_pool_get_config(liConfigs, 'config_SGH_tangent_slope_calculation', config_SGH_tangent_slope_calculation)

      call mpas_pool_get_array(hydroPool, 'waterThickness', waterThickness)
      call mpas_pool_get_array(hydroPool, 'waterPressure', waterPressure)
      call mpas_pool_get_array(hydroPool, 'hydropotentialBase', hydropotentialBase)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
      call mpas_pool_get_array(hydroPool, 'waterThicknessEdge', waterThicknessEdge)
      call mpas_pool_get_array(hydroPool, 'waterThicknessEdgeUpwind', waterThicknessEdgeUpwind)
      call mpas_pool_get_array(hydroPool, 'hydropotentialBaseSlopeNormal', hydropotentialBaseSlopeNormal)
      call mpas_pool_get_array(hydroPool, 'waterPressureSlopeNormal', waterPressureSlopeNormal)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
      call mpas_pool_get_array(hydroPool, 'hydropotentialBaseSlopeTangent', hydropotentialBaseSlopeTangent)
      call mpas_pool_get_array(hydroPool, 'gradMagPhiEdge', gradMagPhiEdge)
      call mpas_pool_get_array(hydroPool, 'effectiveConducEdge', effectiveConducEdge)
      call mpas_pool_get_array(hydroPool, 'diffusivity', diffusivity)
      call mpas_pool_get_array(hydroPool, 'waterVelocity', waterVelocity)
      call mpas_pool_get_array(hydroPool, 'waterFlux', waterFlux)
      call mpas_pool_get_array(hydroPool, 'waterFluxAdvec', waterFluxAdvec)
      call mpas_pool_get_array(hydroPool, 'waterFluxDiffu', waterFluxDiffu)
      call mpas_pool_get_array(hydroPool, 'waterFluxMask', waterFluxMask)


      do iEdge = 1, nEdgesSolve
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
         waterPressureSlopeNormal(iEdge) = (waterPressure(cell2) - waterPressure(cell1)) / dcEdge(iEdge)
      end do


      ! Calculate tangent slope of hydropotentialBase - three possible methods to consider

      ! Calculate hydropotentialBaseVertex if needed
      call mpas_pool_get_array(hydroPool, 'hydropotentialBaseVertex', hydropotentialBaseVertex)  ! this array could be protected by logic if desired
      select case (trim(config_SGH_tangent_slope_calculation))
      case ('from_vertex_barycentric')
         call mpas_pool_get_array(meshPool, 'baryCellsOnVertex', baryCellsOnVertex)
         call mpas_pool_get_array(meshPool, 'baryWeightsOnVertex', baryWeightsOnVertex)
         call mpas_cells_to_points_using_baryweights(meshPool, baryCellsOnVertex(:, 1:nVertices), &
              baryWeightsOnVertex(:, 1:nVertices), hydropotentialBase, hydropotentialBaseVertex(1:nVertices), err_tmp)
              err = ior(err, err_tmp)
      case ('from_vertex_barycentric_kiteareas')
         call li_cells_to_vertices_1dfield_using_kiteAreas(meshPool, hydropotentialBase, hydropotentialBaseVertex)
      end select

      ! Now perform tangent slope calculation based on method chosen
      select case (trim(config_SGH_tangent_slope_calculation))
      case ('from_vertex_barycentric', 'from_vertex_barycentric_kiteareas')
         call mpas_pool_get_array(geometryPool, 'edgeMask', edgeMask)
         call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
         call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
         do iEdge = 1, nEdgesSolve
            ! Only calculate slope for edges that have ice on at least one side.
            if ( li_mask_is_ice(edgeMask(iEdge)) ) then
               hydropotentialBaseSlopeTangent(iEdge) = ( hydropotentialBaseVertex(verticesOnEdge(1,iEdge)) -  &
                     hydropotentialBaseVertex(verticesOnEdge(2,iEdge)) ) / dvEdge(iEdge)
            else
               hydropotentialBaseSlopeTangent(iEdge) = 0.0_RKIND
            endif
         end do  ! edges
      case ('from_normal_slope')
         call mpas_tangential_vector_1d(hydropotentialBaseSlopeNormal, meshPool, &
                  includeHalo=.false., tangentialVector=hydropotentialBaseSlopeTangent)
      case default
         write (stdErrUnit,*) 'Error: Invalid value for config_SGH_tangent_slope_calculation.'
         err = 1
      end select

      ! calculate magnitude of gradient of Phi
      gradMagPhiEdge = sqrt(hydropotentialBaseSlopeNormal**2 + hydropotentialBaseSlopeTangent**2)

      ! calculate effective conductivity on edges
      ! OLD: USE REGULARIZATION:  effectiveConducEdge(:) = conduc_coeff * waterThicknessEdge(:)**(alpha-1.0_RKIND) * (gradMagPhiEdge(:)+1.0e-10_RKIND)**(beta - 2.0_RKIND)   ! 1e-10 used for regularization
      ! Do not calculate the conductivity where it is tiny to avoid blowups
      where(gradMagPhiEdge < 0.01_RKIND)
         effectiveConducEdge(:) = 0.0_RKIND
      elsewhere
         effectiveConducEdge(:) = conduc_coeff * waterThicknessEdge(:)**(alpha-1.0_RKIND) * gradMagPhiEdge(:)**(beta - 2.0_RKIND)
      end where

      ! calculate diffusivity on edges
      diffusivity(:) = rho_water * gravity * effectiveConducEdge(:) * waterThicknessEdge(:)

      do iEdge = 1, nEdgesSolve
         cell1 = cellsOnEdge(1, iEdge)
         cell2 = cellsOnEdge(2, iEdge)
         waterVelocity(iEdge) = -1.0_RKIND * effectiveConducEdge(iEdge) * hydropotentialBaseSlopeNormal(iEdge)
         velSign = sign(1.0_RKIND, waterVelocity(iEdge))
         waterThicknessEdgeUpwind(iEdge) = max(velSign * waterThickness(cell1),   &
                     velSign * (-1.0_RKIND) * waterThickness(cell2))

         ! advective flux
         waterFluxAdvec(iEdge) = waterVelocity(iEdge) * waterThicknessEdgeUpwind(iEdge)

         ! diffusive flux
         numGroundedCells = li_mask_is_grounded_ice_int(cellMask(cell1)) + li_mask_is_grounded_ice_int(cellMask(cell2))
         if (numGroundedCells < 2) then
            waterFluxDiffu(iEdge) = 0.0_RKIND
         else
            waterFluxDiffu(iEdge) = -1.0_RKIND * diffusivity(iEdge) * (waterThickness(cell2) - waterThickness(cell1)) / dcEdge(iEdge)
         endif
      end do
      where (waterFluxMask == 2)
         waterFluxAdvec = 0.0_RKIND
         waterFluxDiffu = 0.0_RKIND
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
   subroutine check_timestep(domain, timeLeft, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

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
!        type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: hydroPool
      real (kind=RKIND), dimension(:), pointer :: waterVelocity
      real (kind=RKIND), dimension(:), pointer :: channelVelocity
      real (kind=RKIND), dimension(:), pointer :: diffusivity
      real (kind=RKIND), dimension(:), pointer :: channelDiffusivity
      real (kind=RKIND), dimension(:), pointer :: dcEdge
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
      logical, pointer :: config_SGH_adaptive_timestep
      logical, pointer :: config_SGH_chnl_active
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


      err = 0
      err_tmp = 0

      call mpas_pool_get_config(liConfigs, 'config_SGH_englacial_porosity', porosity)
      call mpas_pool_get_config(liConfigs, 'config_SGH_adaptive_timestep', config_SGH_adaptive_timestep)
      call mpas_pool_get_config(liConfigs, 'config_SGH_chnl_active', config_SGH_chnl_active)

      dtSGHadvecProc = bigNumber
      dtSGHdiffuProc = bigNumber
      dtSGHpressureProc = bigNumber
      dtSGHadvecChanProc = bigNumber
      dtSGHdiffuChanProc = bigNumber

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

         ! Calculate advective CFL-limited time step
         dtSGHadvecBlock = 0.5_RKIND * minval(dcEdge(1:nEdgesSolve) / (abs(waterVelocity(1:nEdgesSolve)) + 1.0e-12_RKIND))  ! regularize
         dtSGHadvecProc = min(dtSGHadvecProc, dtSGHadvecBlock)

         ! Calculate diffusive CFL-limited time step
         dtSGHdiffuBlock = 0.25_RKIND * minval(dcEdge(1:nEdgesSolve)**2 / (diffusivity(1:nEdgesSolve) + 1.0e-12_RKIND))
         dtSGHdiffuProc = min(dtSGHdiffuProc, dtSGHdiffuBlock)

         ! Calculate pressure limited time step
         dtSGHpressureBlock = 1.0_RKIND * minval(porosity * dcEdge(1:nEdgesSolve)**2 &
                   / (2.0_RKIND * diffusivity(1:nEdgesSolve) + 1.0e-12_RKIND))
         dtSGHpressureProc = min(dtSGHpressureProc, dtSGHpressureBlock)

         if (config_SGH_chnl_active) then
            ! Calculate channel advection limited time step
            dtSGHadvecChanBlock = 0.5_RKIND * minval(dcEdge(1:nEdgesSolve) / (abs(channelVelocity(1:nEdgesSolve)) + 1.0e-12_RKIND))  ! regularize
            dtSGHadvecChanProc = min(dtSGHadvecChanProc, dtSGHadvecChanBlock)
            ! Calculate channel diffusion limited time step
            dtSGHdiffuChanBlock = 0.25_RKIND * minval(dcEdge(1:nEdgesSolve)**2 / (channelDiffusivity(1:nEdgesSolve) + 1.0e-12_RKIND))
            dtSGHdiffuChanProc = min(dtSGHdiffuChanProc, dtSGHdiffuChanBlock)
         endif

         ! Master deltat is needed below, so grab it in this block loop
         call mpas_pool_get_array(meshPool, 'deltat', deltat)

         block => block % next
      end do

      ! reduce across procs
      ! TODO Pack these into array so only need 1 comm?
      call mpas_timer_start("global reduce")
      call mpas_dmpar_min_real(domain % dminfo, dtSGHadvecProc, deltatSGHadvec)
      call mpas_dmpar_min_real(domain % dminfo, dtSGHdiffuProc, deltatSGHdiffu)
      call mpas_dmpar_min_real(domain % dminfo, dtSGHpressureProc, deltatSGHpressure)
      if (config_SGH_chnl_active) then
         call mpas_dmpar_min_real(domain % dminfo, dtSGHadvecChanProc, deltatSGHadvecChannel)
         call mpas_dmpar_min_real(domain % dminfo, dtSGHdiffuChanProc, deltatSGHdiffuChannel)
      endif
      call mpas_timer_stop("global reduce")


      ! Set time step
      if (config_SGH_adaptive_timestep) then
         call mpas_pool_get_config(liConfigs, 'config_SGH_adaptive_timestep_fraction', CFLfraction)
         call mpas_pool_get_config(liConfigs, 'config_SGH_max_adaptive_timestep', maxDt)

         ! Find smallest of 3 or 5 limiting time steps
         proposedDt = min(deltatSGHadvec, deltatSGHdiffu, deltatSGHpressure)
         if (config_SGH_chnl_active) then
            proposedDt = min(proposedDt, deltatSGHadvecChannel, deltatSGHdiffuChannel)
         endif
         proposedDt = proposedDt * CFLfraction


         ! Write out dt info on the final subcycle
         if (proposedDt >= timeLeft) then
            if (config_SGH_chnl_active) then
               write(stdErrUnit,*) "deltatSGH, deltatSGHadvec, deltatSGHdiffu, deltatSGHpressure, deltatSGHadvecChannel, deltatSGHdiffuChannel:"
               write(stdErrUnit,*) deltatSGH, deltatSGHadvec, deltatSGHdiffu, deltatSGHpressure, deltatSGHadvecChannel, deltatSGHdiffuChannel
            else
               write(stdErrUnit,*) "deltatSGH, deltatSGHadvec, deltatSGHdiffu, deltatSGHpressure:"
               write(stdErrUnit,*) deltatSGH, deltatSGHadvec, deltatSGHdiffu, deltatSGHpressure
            endif
         endif


         ! Don't exceed the maximum allowed hydro time step
         proposedDt = min(proposedDt, maxDt)
         ! Don't let SGH time step exceed time left in the master model dt
         proposedDt = min(proposedDt, timeLeft)

         !write(stdOutUnit,*) '  Setting SGH time step (seconds) to:', proposedDt

         timeLeft = timeLeft - proposedDt
         !write (stdErrUnit,*) "dt, new TIMELEFT", proposedDt, timeLeft

      else ! not SGH adaptive
         ! master dt could change if using master adaptive time stepper, so update it in case
         proposedDt = deltat
      endif


      ! Assign new time step to all blocks
      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
         call mpas_pool_get_array(hydroPool, 'deltatSGH', deltatSGH)
         deltatSGH = proposedDt

         block => block % next
      end do


      ! compare timesteps
      if (deltatSGH > deltatSGHadvec) then
         write(stdErrUnit,*) "WARNING: deltatSGH > deltatSGHadvec  ", deltatSGH, deltatSGHadvec
      endif
      if (deltatSGH > deltatSGHdiffu) then
         write(stdErrUnit,*) "WARNING: deltatSGH > deltatSGHdiffu  ", deltatSGH, deltatSGHdiffu
      endif
      if (deltatSGH > deltatSGHpressure) then
         write(stdErrUnit,*) "WARNING: deltatSGH > deltatSGHpressure  ", deltatSGH, deltatSGHpressure
      endif
      if (config_SGH_chnl_active .and. (deltatSGH > deltatSGHadvecChannel)) then
         write(stdErrUnit,*) "WARNING: deltatSGH > deltatSGHadvecChannel  ", deltatSGH, deltatSGHadvecChannel
      endif
      if (config_SGH_chnl_active .and. (deltatSGH > deltatSGHdiffuChannel)) then
         write(stdErrUnit,*) "WARNING: deltatSGH > deltatSGHdiffuChannel  ", deltatSGH, deltatSGHdiffuChannel
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
      real (kind=RKIND), dimension(:), pointer :: effectivePressure
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
      integer, dimension(:), pointer :: cellMask
      real (kind=RKIND), pointer :: deltatSGH
      real (kind=RKIND), pointer :: bedRough, bedRoughMax
      real (kind=RKIND), pointer :: rhoi
      real (kind=RKIND), pointer :: creepCoeff
      real (kind=RKIND), pointer :: porosity
      integer, pointer :: nVertLevels
      character (len=StrKIND), pointer :: config_SGH_pressure_calc
      integer :: err_tmp

      err = 0
      err_tmp = 0


      ! Get pools things
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

      call mpas_pool_get_array(hydroPool, 'effectivePressure', effectivePressure)
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
      call mpas_pool_get_array(velocityPool, 'basalSpeed', basalSpeed)
      call mpas_pool_get_array(velocityPool, 'flowParamA', flowParamA)
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)

      openingRate(:) = bedRough * basalSpeed(:) * (bedRoughMax - waterThickness(:))
      !openingRate(:) = bedRough * basalSpeed(:) * (bedRoughMax - waterThickness(:)) + basalMeltInput(:) / rhoi  ! Hewitt 2011 opening
      openingRate = max(0.0_RKIND, openingRate)

      closingRate(:) = creepCoeff * flowParamA(nVertLevels, :) * effectivePressure(:)**3 * waterThickness(:)
!      closingRate(:) = waterThickness(:) * effectivePressure(:) / 1.0e13_RKIND  ! Hewitt 2011 creep closure form.  Denominator is ice viscosity

      zeroOrderSum = closingRate - openingRate + (basalMeltInput + externalWaterInput) / rho_water - (Wtill - WtillOld) / deltatSGH

      waterPressureOld = waterPressure

      select case (trim(config_SGH_pressure_calc))
      case ('cavity')

         where (li_mask_is_floating_ice(cellMask))
            waterPressure = rhoi * gravity * thickness
         elsewhere (.not. li_mask_is_ice(cellMask))
            waterPressure = 0.0_RKIND
         elsewhere
            waterPressure = (zeroOrderSum - divergence - divergenceChannel - channelAreaChangeCell) * rho_water * gravity * deltatSGH / porosity + waterPressureOld
         end where

      case ('overburden')
         where (li_mask_is_floating_ice(cellMask))
            waterPressure = rhoi * gravity * thickness
         elsewhere (.not. li_mask_is_ice(cellMask))
            waterPressure = 0.0_RKIND
         elsewhere
            waterPressure = rhoi * gravity * thickness
         end where

      case default
         write(stdErrUnit,*) "ERROR: Invalid option specified for config_SGH_pressure_calc:", config_SGH_pressure_calc
         err = ior(err, 1)
      end select

      waterPressure = max(0.0_RKIND, waterPressure)
      waterPressure = min(waterPressure, rhoi * gravity * thickness)
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
      real (kind=RKIND), pointer :: rhoi
      real (kind=RKIND), dimension(:), pointer :: thickness
      real (kind=RKIND), dimension(:), pointer :: waterPressure
      real (kind=RKIND), dimension(:), pointer :: bedTopography
      real (kind=RKIND), dimension(:), pointer :: hydropotentialBase
      real (kind=RKIND), dimension(:), pointer :: waterThickness
      real (kind=RKIND), dimension(:), pointer :: hydropotential
      real (kind=RKIND), dimension(:), pointer :: effectivePressure

      err = 0

      ! Get pools things
      call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
      call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)

      call mpas_pool_get_config(liConfigs, 'config_ice_density', rhoi)

      call mpas_pool_get_array(hydroPool, 'effectivePressure', effectivePressure)
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(hydroPool, 'waterPressure', waterPressure)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
      call mpas_pool_get_array(hydroPool, 'hydropotentialBase', hydropotentialBase)
      call mpas_pool_get_array(hydroPool, 'waterThickness', waterThickness)
      call mpas_pool_get_array(hydroPool, 'hydropotential', hydropotential)

      effectivePressure = rhoi * gravity * thickness - waterPressure
      ! This is the hydropotential without the water thickness piece
      hydropotentialBase = rho_water * gravity * bedTopography + waterPressure
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
!!      type (mpas_pool_type), pointer :: geometryPool
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
      real (kind=RKIND), dimension(:), pointer :: hydropotentialBaseSlopeNormal
      real (kind=RKIND), dimension(:), pointer :: waterPressureSlopeNormal
      real (kind=RKIND), dimension(:), pointer :: channelOpeningRate
      real (kind=RKIND), dimension(:), pointer :: channelClosingRate
      real (kind=RKIND), dimension(:), pointer :: channelChangeRate
      real (kind=RKIND), dimension(:), pointer :: flowParamAChannel
      real (kind=RKIND), dimension(:), pointer :: channelEffectivePressure
      real (kind=RKIND), dimension(:), pointer :: effectivePressure
      real (kind=RKIND), dimension(:), pointer :: channelDiffusivity
      integer, dimension(:), pointer :: waterFluxMask
      integer, dimension(:), pointer :: edgeMask
      real (kind=RKIND), dimension(:,:), pointer :: flowParamA
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer, pointer :: nVertLevels

      integer, pointer :: nEdgesSolve
      integer :: iEdge, cell1, cell2


      err = 0

      ! Get pools things
      call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
!      call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
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

      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(hydroPool, 'channelArea', channelArea)
      call mpas_pool_get_array(hydroPool, 'channelMelt', channelMelt)
      call mpas_pool_get_array(hydroPool, 'channelPressureFreeze', channelPressureFreeze)
      call mpas_pool_get_array(hydroPool, 'channelDischarge', channelDischarge)
      call mpas_pool_get_array(hydroPool, 'channelVelocity', channelVelocity)
      call mpas_pool_get_array(hydroPool, 'gradMagPhiEdge', gradMagPhiEdge)
      call mpas_pool_get_array(hydroPool, 'hydropotentialBaseSlopeNormal', hydropotentialBaseSlopeNormal)
      call mpas_pool_get_array(hydroPool, 'waterPressureSlopeNormal', waterPressureSlopeNormal)
      call mpas_pool_get_array(hydroPool, 'waterFlux', waterFlux)
      call mpas_pool_get_array(hydroPool, 'channelOpeningRate', channelOpeningRate)
      call mpas_pool_get_array(hydroPool, 'channelClosingRate', channelClosingRate)
      call mpas_pool_get_array(hydroPool, 'channelChangeRate', channelChangeRate)
      call mpas_pool_get_array(hydroPool, 'flowParamAChannel', flowParamAChannel)
      call mpas_pool_get_array(hydroPool, 'channelEffectivePressure', channelEffectivePressure)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(velocityPool, 'flowParamA', flowParamA)
      call mpas_pool_get_array(hydroPool, 'effectivePressure', effectivePressure)
      call mpas_pool_get_array(hydroPool, 'waterFluxMask', waterFluxMask)
      call mpas_pool_get_array(hydroPool, 'channelDiffusivity', channelDiffusivity)
      call mpas_pool_get_array(geometryPool, 'edgeMask', edgeMask)

      ! Calculate terms needed for opening (melt) rate

      where(gradMagPhiEdge < 0.01_RKIND)
         channelDischarge(:) = 0.0_RKIND
      elsewhere
         channelDischarge = -1.0_RKIND * Kc * channelArea**alpha_c * gradMagPhiEdge**(beta_c - 2.0_RKIND) * hydropotentialBaseSlopeNormal
      end where

      where (waterFluxMask == 2)
         channelDischarge = 0.0_RKIND
         channelArea = 0.0_RKIND
      end where

      where (.not. li_mask_is_ice(edgeMask))
         channelArea = 0.0_RKIND
         channelDischarge = 0.0_RKIND
      end where

      channelVelocity = channelDischarge / (channelArea + 1.0e-12_RKIND)

      ! diffusivity used only to limit channel dt right now
      where(gradMagPhiEdge < 0.01_RKIND)
         channelDiffusivity = 0.0_RKIND
      elsewhere
         channelDiffusivity = abs(rho_water * gravity * channelArea *  &
            Kc * channelArea**(alpha_c - 1.0_RKIND) * gradMagPhiEdge**(beta_c - 2.0_RKIND))
      end where

      channelMelt = (abs(channelDischarge * hydropotentialBaseSlopeNormal) &  ! channel dissipation
                  +  abs(waterFlux * hydropotentialBaseSlopeNormal * config_SGH_incipient_channel_width)        &  ! some sheet dissipation
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
         flowParamAChannel(iEdge) = 0.5_RKIND * ( flowParamA(nVertLevels, cell1) + flowParamA(nVertLevels, cell2) )
         channelEffectivePressure(iEdge) = 0.5_RKIND * (effectivePressure(cell1) + effectivePressure(cell2))
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
      real (kind=RKIND), dimension(:), pointer :: channelChangeRate
      real (kind=RKIND), dimension(:), pointer :: divergenceChannel
      real (kind=RKIND), dimension(:), pointer :: channelAreaChangeCell
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
      call mpas_pool_get_array(hydroPool, 'channelChangeRate', channelChangeRate)
      call mpas_pool_get_array(hydroPool, 'divergenceChannel', divergenceChannel)
      call mpas_pool_get_array(hydroPool, 'channelAreaChangeCell', channelAreaChangeCell)
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
      ! loop over locally owned cells
      do iCell = 1, nCellsSolve
         ! compute fluxes for each edge of the cell
         do iEdgeOnCell = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(iEdgeOnCell, iCell)
            ! add on advective & diffusive fluxes
            divergenceChannel(iCell) = divergenceChannel(iCell) - channelDischarge(iEdge) * edgeSignOnCell(iEdgeOnCell, iCell)
            channelAreaChangeCell(iCell) = channelChangeRate(iEdge) * dcEdge(iEdge) * 0.5_RKIND  ! only half of channel is in this cell
         end do ! edges
      end do ! cells
      divergenceChannel(1:nCellsSolve) = divergenceChannel(1:nCellsSolve) / areaCell(1:nCellsSolve)
      channelAreaChangeCell(1:nCellsSolve) = channelAreaChangeCell(1:nCellsSolve) / areaCell(1:nCellsSolve)

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
            write(stdErrUnit,*) "ERROR: Unknown value for config_SGH_shmip_forcing!"
            err = ior(err, 1)
         end select
         where (externalWaterInput > 0.0_RKIND)
            basalMeltInput = 0.9 / areaCell * (1.0_RKIND - ra * sin(2.0_RKIND * pii * daysSinceStart))
         elsewhere
            basalMeltInput = 0.0_RKIND
         end where
         where (basalMeltInput < 0.0_RKIND)
            basalMeltInput = 0.0_RKIND
         end where
         basalMeltInput = (basalMeltInput + 7.93e-11_RKIND) / 1000.0_RKIND

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
            write(stdErrUnit,*) "ERROR: Unknown value for config_SGH_shmip_forcing!"
            err = ior(err, 1)
         end select

         call mpas_pool_get_array(geometryPool, 'upperSurface', upperSurface)
         temperature0 = -16.0_RKIND * cos(2.0_RKIND * pii / 365.0 * daysSinceStart) - 5.0_RKIND + DT
         externalWaterInput(:) = (upperSurface(:) * (-0.0075_RKIND) + temperature0) * (0.01_RKIND / 86400.0_RKIND) * 1000.0_RKIND  ! the 1000 factor converts from m/s to kg/m2/s
         where (externalWaterInput < 0.0_RKIND)
            externalWaterInput = 0.0_RKIND
         end where

      case default
         write(stdErrUnit,*) "ERROR: Unknown value for config_SGH_shmip_forcing!"
         err = ior(err, 1)
      end select


   !--------------------------------------------------------------------
   end subroutine shmip_timevarying_forcing


end module li_subglacial_hydro
