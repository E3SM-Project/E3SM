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
   use mpas_constants
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
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: hydroPool
      real (kind=RKIND), pointer :: deltat
      real (kind=RKIND), pointer :: deltatSGH
      integer :: err_tmp

      err = 0
      err_tmp = 0

      call mpas_pool_get_config(liConfigs, 'config_SGH', config_SGH)
      if (.not. config_SGH) then
         ! If SGH is not active, skip everything
         return
      endif

      write(stdoutUnit,*) 'Beginning subglacial hydro init.'
!      call mpas_pool_get_config(liConfigs, 'config_time_integration', config_time_integration)
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
         call mpas_pool_get_array(meshPool, 'deltat', deltat)
         call mpas_pool_get_array(hydroPool, 'deltatSGH', deltatSGH)

         ! For now, force hydro dt to be master dt
         deltatSGH = deltat

         block => block % next
      end do

      ! === error check
      if (err > 0) then
          write (stderrUnit,*) "An error has occurred in li_SGH_init."
      endif

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
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: hydroPool
      type (mpas_pool_type), pointer :: meshPool
      real (kind=RKIND), dimension(:), pointer :: Wtill, WtillOld
      real (kind=RKIND), dimension(:), pointer :: meltInput
      real (kind=RKIND), dimension(:), pointer :: divergence
      real (kind=RKIND), dimension(:), pointer :: diffusivity
      real (kind=RKIND), dimension(:), pointer :: waterFlux
      real (kind=RKIND), dimension(:), pointer :: waterThickness
      real (kind=RKIND), dimension(:), pointer :: waterThicknessOld
      real (kind=RKIND), dimension(:), pointer :: dcEdge
      real (kind=RKIND), dimension(:), pointer :: dvEdge
      real (kind=RKIND), dimension(:), pointer :: areaCell
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer, dimension(:,:), pointer :: edgeSignOnCell
      integer, dimension(:), pointer :: cellMask
      real (kind=RKIND), pointer :: deltatSGH
      real (kind=RKIND), pointer :: deltat
      real (kind=RKIND), pointer :: Cd
      real (kind=RKIND), pointer :: tillMax
      integer, pointer :: nCellsSolve
      integer :: iCell, iEdge, iEdgeOnCell
      integer :: cell1, cell2
      integer :: err_tmp



      err = 0
      err_tmp = 0

      call mpas_pool_get_config(liConfigs, 'config_SGH', config_SGH)
      if (.not. config_SGH) then
         ! If SGH is not active, skip everything
         return
      endif

      write(stdoutUnit,*) 'Beginning subglacial hydro solve.'
      call mpas_pool_get_config(liConfigs, 'config_SGH_till_drainage', Cd)
      call mpas_pool_get_config(liConfigs, 'config_SGH_till_max', tillMax)

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
         call mpas_pool_get_array(hydroPool, 'meltInput', meltInput)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)

         WtillOld = Wtill
         Wtill = Wtill + deltatSGH * (meltInput / rho_water - Cd)
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
      call mpas_dmpar_field_halo_exch(domain, 'diffusivity')
      call mpas_dmpar_field_halo_exch(domain, 'waterFlux')
      ! intermediate fields will be out of date, but will be correct in output files
      call mpas_timer_stop("halo updates")


      ! =============
      ! Calculate time step if being adaptive, otherwise check for CFL
      ! =============
      call check_timestep(domain, err_tmp)
      err = ior(err, err_tmp)


      ! =============
      ! Compute flux divergence
      ! =============
      !call calc_flux_div(divergence, waterFlux, diffusivity, dcEdge)

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
         call mpas_pool_get_array(hydroPool, 'divergence', divergence)
         call mpas_pool_get_array(hydroPool, 'diffusivity', diffusivity)
         call mpas_pool_get_array(hydroPool, 'waterFlux', waterFlux)
         call mpas_pool_get_array(hydroPool, 'waterThickness', waterThickness)
         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
         call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
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
               cell1 = cellsOnEdge(1, iEdge)
               cell2 = cellsOnEdge(2, iEdge)
               ! add on advective flux
               divergence(iCell) = divergence(iCell) - waterFlux(iEdge) * edgeSignOnCell(iEdgeOnCell, iCell)
               ! add on diffusive flux
               divergence(iCell) = divergence(iCell) + diffusivity(iEdge) * &
                 (waterThickness(cell1) * edgeSignOnCell(iEdgeOnCell, iCell) - &
                  waterThickness(cell2) * edgeSignOnCell(iEdgeOnCell, iCell)) / dcEdge(iEdge) * dvEdge(iEdge)
            end do ! edges
         end do ! cells
         divergence(1:nCellsSolve) = divergence(1:nCellsSolve) / areaCell(1:nCellsSolve)
         block => block % next
      end do


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
         call mpas_pool_get_array(hydroPool, 'tillWaterThickness', Wtill)
         call mpas_pool_get_array(hydroPool, 'tillWaterThicknessOld', WtillOld)
         call mpas_pool_get_array(hydroPool, 'deltatSGH', deltatSGH)
         call mpas_pool_get_array(hydroPool, 'meltInput', meltInput)
         call mpas_pool_get_array(hydroPool, 'divergence', divergence)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)

         waterThicknessOld = waterThickness
         waterThickness = waterThicknessOld + deltatSGH * (meltInput / rho_water - divergence  &
             - (Wtill - WtillOld) / deltatSGH)
         waterThickness = waterThickness * li_mask_is_grounded_ice_int(cellMask)  ! zero in non-grounded locations
         waterThickness = max(0.0_RKIND, waterThickness)
         divergence = divergence * li_mask_is_grounded_ice_int(cellMask)  ! zero in non-grounded locations for more convenient viz

         block => block % next
      end do


      ! =============
      ! Update model time
      ! =============


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
!         ! Assign the time stamp for this time step
!         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
!         call mpas_pool_get_array(meshPool, 'xtime', xtime)
!         xtime = timeStamp
!
!         call mpas_pool_get_array(meshPool, 'deltat', deltat_output)
!         deltat_output = dtSeconds
!
!         ! compute time since start of simulation, in days
!         call mpas_pool_get_array(meshPool, 'simulationStartTime', simulationStartTime)
!         call mpas_pool_get_array(meshPool, 'daysSinceStart',daysSinceStart)
!         call mpas_set_time(simulationStartTime_timeType, dateTimeString=simulationStartTime)
!         call mpas_get_timeInterval(currTime - simulationStartTime_timeType, dt=daysSinceStart)
!         daysSinceStart = daysSinceStart / seconds_per_day

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

      use mpas_vector_operations, only: mpas_tangential_vector_1d
      use mpas_vector_reconstruction

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
      real (kind=RKIND), dimension(:), pointer :: waterPressure
      real (kind=RKIND), dimension(:), pointer :: hydropotentialBase
      real (kind=RKIND), dimension(:), pointer :: waterThicknessEdge
      real (kind=RKIND), dimension(:), pointer :: waterThicknessEdgeUpwind
      real (kind=RKIND), dimension(:), pointer :: waterThickness
      real (kind=RKIND), dimension(:), pointer :: hydropotentialBaseSlopeNormal
      real (kind=RKIND), dimension(:), pointer :: hydropotentialBaseSlopeTangent
      real (kind=RKIND), dimension(:), pointer :: gradMagPhiEdge
      real (kind=RKIND), dimension(:), pointer :: effectiveConducEdge
      real (kind=RKIND), dimension(:), pointer :: diffusivity
      real (kind=RKIND), dimension(:), pointer :: dcEdge
      real (kind=RKIND), dimension(:), pointer :: dvEdge
      real (kind=RKIND), dimension(:), pointer :: waterVelocity
      real (kind=RKIND), dimension(:), pointer :: waterFlux
      real (kind=RKIND), dimension(:), pointer :: waterVelocityCellX
      real (kind=RKIND), dimension(:), pointer :: waterVelocityCellY
      integer, dimension(:), pointer :: cellMask
      real (kind=RKIND), dimension(:), allocatable :: cellJunk
      integer, dimension(:,:), pointer :: cellsOnEdge
      real (kind=RKIND), pointer :: alpha, beta
      real (kind=RKIND), pointer :: conduc_coeff
      character (len=StrKIND), pointer :: config_SGH_tangent_slope_calculation
      integer, pointer :: nEdgesSolve
      integer, pointer :: nCells
      integer :: iEdge, cell1, cell2
      real (kind=RKIND) :: velSign
      integer :: err_tmp

      err = 0
      err_tmp = 0


      ! Get pools things
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
      call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)

      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      call mpas_pool_get_config(liConfigs, 'config_SGH_alpha', alpha)
      call mpas_pool_get_config(liConfigs, 'config_SGH_beta', beta)
      call mpas_pool_get_config(liConfigs, 'config_SGH_conduc_coeff', conduc_coeff)
      call mpas_pool_get_config(liConfigs, 'config_SGH_tangent_slope_calculation', config_SGH_tangent_slope_calculation)

      call mpas_pool_get_array(hydroPool, 'waterThickness', waterThickness)
      call mpas_pool_get_array(hydroPool, 'hydropotentialBase', hydropotentialBase)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
      call mpas_pool_get_array(hydroPool, 'waterPressure', waterPressure)
      call mpas_pool_get_array(hydroPool, 'waterThicknessEdge', waterThicknessEdge)
      call mpas_pool_get_array(hydroPool, 'waterThicknessEdgeUpwind', waterThicknessEdgeUpwind)
      call mpas_pool_get_array(hydroPool, 'hydropotentialBaseSlopeNormal', hydropotentialBaseSlopeNormal)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(hydroPool, 'hydropotentialBaseSlopeTangent', hydropotentialBaseSlopeTangent)
      call mpas_pool_get_array(hydroPool, 'gradMagPhiEdge', gradMagPhiEdge)
      call mpas_pool_get_array(hydroPool, 'effectiveConducEdge', effectiveConducEdge)
      call mpas_pool_get_array(hydroPool, 'diffusivity', diffusivity)
      call mpas_pool_get_array(hydroPool, 'waterVelocity', waterVelocity)
      call mpas_pool_get_array(hydroPool, 'waterFlux', waterFlux)
      call mpas_pool_get_array(hydroPool, 'waterVelocityCellX', waterVelocityCellX)
      call mpas_pool_get_array(hydroPool, 'waterVelocityCellY', waterVelocityCellY)


      ! This is the hydropotential without the water thickness piece
      hydropotentialBase = rho_water * gravity * bedTopography + waterPressure

      do iEdge = 1, nEdgesSolve
         cell1 = cellsOnEdge(1, iEdge)
         cell2 = cellsOnEdge(2, iEdge)

         waterThicknessEdge(iEdge) = 0.5_RKIND * ( waterThickness(cell1) + waterThickness(cell2) )

         hydropotentialBaseSlopeNormal(iEdge) = (hydropotentialBase(cell1) - hydropotentialBase(cell2)) / dcEdge(iEdge)
      end do


      ! Calculate tangent slope of hydropotentialBase
      select case (trim(config_SGH_tangent_slope_calculation))
      case ('from_vertex_barycentric', 'from_vertex_barycentric_kiteareas')
         do iEdge = 1, nEdgesSolve
!            ! Only calculate slope for edges that have ice on at least one side.
!            if ( li_mask_is_dynamic_ice(edgeMask(iEdge)) ) then
!               hydropotentialBaseSlopeTangent(iEdge) = ( upperSurfaceVertex(verticesOnEdge(1,iEdge)) -  &
!                     upperSurfaceVertex(verticesOnEdge(2,iEdge)) ) / dvEdge(iEdge)
!            else
!               hydropotentialBaseSlopeTangent(iEdge) = 0.0_RKIND
!            endif
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
      effectiveConducEdge(:) = conduc_coeff * waterThicknessEdge(:)**(alpha-1.0_RKIND) * (gradMagPhiEdge(:)+1.0e-10_RKIND)**(beta - 2.0_RKIND)

      ! calculate diffusivity on edges
      diffusivity(:) = rho_water * gravity * effectiveConducEdge(:) * waterThicknessEdge(:)

      do iEdge = 1, nEdgesSolve
         cell1 = cellsOnEdge(1, iEdge)
         cell2 = cellsOnEdge(2, iEdge)
         waterVelocity(iEdge) = 1.0_RKIND / dcEdge(iEdge) * effectiveConducEdge(iEdge) * &
             (hydropotentialBase(cell1) - hydropotentialBase(cell2))
         velSign = sign(1.0_RKIND, waterVelocity(iEdge))
         waterThicknessEdgeUpwind(iEdge) = max(velSign * waterThickness(cell1),   &
                     velSign * (-1.0_RKIND) * waterThickness(cell2))
         waterFlux(iEdge) = waterVelocity(iEdge) * dvEdge(iEdge) * waterThicknessEdgeUpwind(iEdge)
      end do

      ! Calculate reconstructed velocities on cell centers for viz
      allocate(cellJunk(nCells+1))
      call mpas_reconstruct(meshPool, waterVelocity, waterVelocityCellX, waterVelocityCellY, &
              cellJunk, cellJunk, cellJunk)
      deallocate(cellJunk)
      waterVelocityCellX = waterVelocityCellX * li_mask_is_grounded_ice_int(cellMask)  ! zero in non-grounded locations
      waterVelocityCellY = waterVelocityCellY * li_mask_is_grounded_ice_int(cellMask)  ! zero in non-grounded locations

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
   subroutine check_timestep(domain, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain    !< Input/Output: domain object

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
      real (kind=RKIND), dimension(:), pointer :: diffusivity
      real (kind=RKIND), dimension(:), pointer :: dcEdge
      real (kind=RKIND), pointer :: deltatSGH
      real (kind=RKIND), pointer :: deltatSGHadvec
      real (kind=RKIND), pointer :: deltatSGHdiffu
      real (kind=RKIND), pointer :: deltatSGHpressure
      real (kind=RKIND), pointer :: porosity
      type (block_type), pointer :: block
      integer, pointer :: nEdgesSolve
      real (kind=RKIND) :: dtSGHadvecBlock, dtSGHadvecProc
      real (kind=RKIND) :: dtSGHdiffuBlock, dtSGHdiffuProc
      real (kind=RKIND) :: dtSGHpressureBlock, dtSGHpressureProc
      integer :: err_tmp
      real(kind=RKIND) :: bigNumber
      parameter (bigNumber = 1.0e36_RKIND)

      err = 0
      err_tmp = 0

      call mpas_pool_get_config(liConfigs, 'config_SGH_englacial_porosity', porosity)

      dtSGHadvecProc = bigNumber
      dtSGHdiffuProc = bigNumber
      dtSGHpressureProc = bigNumber

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
         call mpas_pool_get_array(hydroPool, 'waterVelocity', waterVelocity)
         call mpas_pool_get_array(hydroPool, 'diffusivity', diffusivity)
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

         block => block % next
      end do

      call mpas_dmpar_min_real(domain % dminfo, dtSGHadvecProc, deltatSGHadvec)
      call mpas_dmpar_min_real(domain % dminfo, dtSGHdiffuProc, deltatSGHdiffu)
      call mpas_dmpar_min_real(domain % dminfo, dtSGHpressureProc, deltatSGHpressure)

      ! optional: always write out for debugging
      write(stdErrUnit,*) "deltatSGH, deltatSGHadvec, deltatSGHdiffu, deltatSGHpressure:", deltatSGH, deltatSGHadvec, deltatSGHdiffu, deltatSGHpressure

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
!      real (kind=RKIND), dimension(:), pointer :: bedTopography
      real (kind=RKIND), dimension(:), pointer :: waterPressure
      real (kind=RKIND), dimension(:), pointer :: waterPressureOld
      real (kind=RKIND), dimension(:), pointer :: waterThickness
      real (kind=RKIND), dimension(:), pointer :: effectivePressure
      real (kind=RKIND), dimension(:), pointer :: zeroOrderSum
      real (kind=RKIND), dimension(:), pointer :: closingRate
      real (kind=RKIND), dimension(:), pointer :: openingRate
      real (kind=RKIND), dimension(:), pointer :: meltInput
      real (kind=RKIND), dimension(:), pointer :: Wtill, WtillOld
      real (kind=RKIND), dimension(:), pointer :: divergence
      real (kind=RKIND), dimension(:), pointer :: basalSpeed
      real (kind=RKIND), dimension(:,:), pointer :: flowParamA
      real (kind=RKIND), dimension(:), pointer :: bedTopography
      real (kind=RKIND), dimension(:), pointer :: thickness
      integer, dimension(:), pointer :: cellMask
      real (kind=RKIND), pointer :: deltatSGH
      real (kind=RKIND), pointer :: bedRough, bedRoughMax
      real (kind=RKIND), pointer :: rhoi
      real (kind=RKIND), pointer :: creepCoeff
      real (kind=RKIND), pointer :: porosity
      integer, pointer :: nVertLevels
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

      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(hydroPool, 'effectivePressure', effectivePressure)
      call mpas_pool_get_array(hydroPool, 'waterPressure', waterPressure)
      call mpas_pool_get_array(hydroPool, 'waterPressureOld', waterPressureOld)
      call mpas_pool_get_array(hydroPool, 'waterThickness', waterThickness)
      call mpas_pool_get_array(hydroPool, 'zeroOrderSum', zeroOrderSum)
      call mpas_pool_get_array(hydroPool, 'closingRate', closingRate)
      call mpas_pool_get_array(hydroPool, 'openingRate', openingRate)
      call mpas_pool_get_array(hydroPool, 'deltatSGH', deltatSGH)
      call mpas_pool_get_array(hydroPool, 'meltInput', meltInput)
      call mpas_pool_get_array(hydroPool, 'tillWaterThickness', Wtill)
      call mpas_pool_get_array(hydroPool, 'tillWaterThicknessOld', WtillOld)
      call mpas_pool_get_array(hydroPool, 'divergence', divergence)
      call mpas_pool_get_array(velocityPool, 'basalSpeed', basalSpeed)
      call mpas_pool_get_array(velocityPool, 'flowParamA', flowParamA)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)

      openingRate(:) = bedRough * basalSpeed(:) * (bedRoughMax - waterThickness(:))
      openingRate = max(0.0_RKIND, openingRate)

      closingRate(:) = creepCoeff * flowParamA(nVertLevels, :) * effectivePressure(:)**3 * waterThickness(:)

      zeroOrderSum = closingRate - openingRate + meltInput / rho_water - (Wtill - WtillOld) / deltatSGH

      waterPressureOld = waterPressure

      where (li_mask_is_floating_ice(cellMask))
         waterPressure = rhoi * gravity * thickness
      elsewhere (.not. li_mask_is_ice(cellMask))
         waterPressure = 0.0_RKIND
      elsewhere
         waterPressure = (zeroOrderSum - divergence) * rho_water * gravity * deltatSGH / porosity + waterPressureOld
      end where
      waterPressure = max(0.0_RKIND, waterPressure)
      waterPressure = min(waterPressure, rhoi * gravity * thickness)

      effectivePressure = rhoi * gravity * thickness - waterPressure

   !--------------------------------------------------------------------
   end subroutine calc_pressure


end module li_subglacial_hydro
