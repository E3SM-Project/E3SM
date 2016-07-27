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
      real (kind=RKIND), dimension(:), pointer :: basalMeltInput
      real (kind=RKIND), dimension(:), pointer :: externalWaterInput
      real (kind=RKIND), dimension(:), pointer :: divergence
      real (kind=RKIND), dimension(:), pointer :: waterFlux
      real (kind=RKIND), dimension(:), pointer :: waterThickness
      real (kind=RKIND), dimension(:), pointer :: waterThicknessOld
      real (kind=RKIND), dimension(:), pointer :: waterThicknessTendency
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
      call mpas_timer_stop("halo updates")


      ! =============
      ! Calculate time step if being adaptive, otherwise check for CFL
      ! =============
      call check_timestep(domain, err_tmp)
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
               cell1 = cellsOnEdge(1, iEdge)
               cell2 = cellsOnEdge(2, iEdge)
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
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)

         waterThicknessOld = waterThickness
         waterThickness = waterThicknessOld + deltatSGH * ( (basalMeltInput + externalWaterInput) / rho_water - divergence  &
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
      real (kind=RKIND), dimension(:), pointer :: waterVelocity
      real (kind=RKIND), dimension(:), pointer :: waterFlux
      real (kind=RKIND), dimension(:), pointer :: waterFluxAdvec
      real (kind=RKIND), dimension(:), pointer :: waterFluxDiffu
      real (kind=RKIND), dimension(:), pointer :: waterVelocityCellX
      real (kind=RKIND), dimension(:), pointer :: waterVelocityCellY
      integer, dimension(:,:), pointer :: edgeSignOnCell
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

      call mpas_pool_get_config(liConfigs, 'config_SGH_alpha', alpha)
      call mpas_pool_get_config(liConfigs, 'config_SGH_beta', beta)
      call mpas_pool_get_config(liConfigs, 'config_SGH_conduc_coeff', conduc_coeff)
      call mpas_pool_get_config(liConfigs, 'config_SGH_tangent_slope_calculation', config_SGH_tangent_slope_calculation)

      call mpas_pool_get_array(hydroPool, 'waterThickness', waterThickness)
      call mpas_pool_get_array(hydroPool, 'hydropotentialBase', hydropotentialBase)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
      call mpas_pool_get_array(hydroPool, 'waterThicknessEdge', waterThicknessEdge)
      call mpas_pool_get_array(hydroPool, 'waterThicknessEdgeUpwind', waterThicknessEdgeUpwind)
      call mpas_pool_get_array(hydroPool, 'hydropotentialBaseSlopeNormal', hydropotentialBaseSlopeNormal)
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
      call mpas_pool_get_array(hydroPool, 'waterVelocityCellX', waterVelocityCellX)
      call mpas_pool_get_array(hydroPool, 'waterVelocityCellY', waterVelocityCellY)



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
      end do


      ! Calculate tangent slope of hydropotentialBase
      select case (trim(config_SGH_tangent_slope_calculation))
      case ('from_vertex_barycentric', 'from_vertex_barycentric_kiteareas')
         ! TODO: add these options if desired
         write (stdErrUnit,*) "Error: 'from_vertex_barycentric' and 'from_vertex_barycentric_kiteareas' are not yet supported for the config_SGH_tangent_slope_calculation option."
         err = 1
!         do iEdge = 1, nEdgesSolve
!            ! Only calculate slope for edges that have ice on at least one side.
!            if ( li_mask_is_dynamic_ice(edgeMask(iEdge)) ) then
!               hydropotentialBaseSlopeTangent(iEdge) = ( upperSurfaceVertex(verticesOnEdge(1,iEdge)) -  &
!                     upperSurfaceVertex(verticesOnEdge(2,iEdge)) ) / dvEdge(iEdge)
!            else
!               hydropotentialBaseSlopeTangent(iEdge) = 0.0_RKIND
!            endif
!         end do  ! edges
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
      waterFlux(:) = waterFluxAdvec(:) + waterFluxDiffu(:)


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
      real (kind=RKIND), pointer :: deltat
      integer, pointer :: nEdgesSolve
      logical, pointer :: config_SGH_adaptive_timestep
      ! in the following variables, "Proc" indicates the value on the current processor,
      ! and "Block" indicates value on current block
      real (kind=RKIND) :: dtSGHadvecBlock, dtSGHadvecProc
      real (kind=RKIND) :: dtSGHdiffuBlock, dtSGHdiffuProc
      real (kind=RKIND) :: dtSGHpressureBlock, dtSGHpressureProc
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

         ! Master deltat is needed below, so grab it in this block loop
         call mpas_pool_get_array(meshPool, 'deltat', deltat)

         block => block % next
      end do

      ! reduce across procs
      call mpas_dmpar_min_real(domain % dminfo, dtSGHadvecProc, deltatSGHadvec)
      call mpas_dmpar_min_real(domain % dminfo, dtSGHdiffuProc, deltatSGHdiffu)
      call mpas_dmpar_min_real(domain % dminfo, dtSGHpressureProc, deltatSGHpressure)


      ! Set time step
      if (config_SGH_adaptive_timestep) then
         call mpas_pool_get_config(liConfigs, 'config_SGH_adaptive_timestep_fraction', CFLfraction)
         call mpas_pool_get_config(liConfigs, 'config_SGH_max_adaptive_timestep', maxDt)

         ! Find smallest of 3 limiting time steps
         proposedDt = min(deltatSGHadvec, deltatSGHdiffu, deltatSGHpressure)
         proposedDt = proposedDt * CFLfraction

         ! Don't let SGH time step exceed master time step - technically need to check the interval til completion, because we might already be partway through the deltat
         proposedDt = min(proposedDt, maxDt)

         write(stdOutUnit,*) '  Setting SGH time step (seconds) to:', proposedDt

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

         ! For now, hijack master dt. TODO subcycle hydro under master dt
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_array(meshPool, 'deltat', deltat)
         deltat = proposedDt


         block => block % next
      end do


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
            waterPressure = (zeroOrderSum - divergence) * rho_water * gravity * deltatSGH / porosity + waterPressureOld
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


end module li_subglacial_hydro
