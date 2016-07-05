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

module li_SGH_driver

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar
   use mpas_timekeeping
   use li_setup
   use li_constants

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
             li_SGH_final

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
      logical, pointer :: config_subglacial_hydro
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: hydroPool
      real (kind=RKIND), pointer :: deltat
      real (kind=RKIND), pointer :: deltatSGH
      integer :: err_tmp

      err = 0
      err_tmp = 0

      call mpas_pool_get_config(liConfigs, 'config_subglacial_hydro', config_subglacial_hydro)
      if (.not. config_subglacial_hydro) then
         ! If SGH is not active, skip everything
         return
      endif

!      call mpas_pool_get_config(liConfigs, 'config_time_integration', config_time_integration)
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
         call mpas_pool_get_array(meshPool, 'deltat', deltat)
         call mpas_pool_get_array(meshPool, 'deltatSGH', deltatSGH)

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

      use li_mask

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
      logical, pointer :: config_subglacial_hydro
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: hydroPool
      type (mpas_pool_type), pointer :: meshPool
      real (kind=RKIND), dimension(:), pointer :: Wtill
      real (kind=RKIND), dimension(:), pointer :: meltInput
      real (kind=RKIND), dimension(:), pointer :: divergence
      real (kind=RKIND), dimension(:), pointer :: diffusivity
      real (kind=RKIND), dimension(:), pointer :: waterFlux
      real (kind=RKIND), dimension(:), pointer :: waterThickness
      real (kind=RKIND), dimension(:), pointer :: dcEdge
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer, dimension(:,:), pointer :: edgeSignOnCell
      integer, dimension(:), pointer :: cellMask
      real (kind=RKIND), pointer :: deltatSGH
      real (kind=RKIND), pointer :: deltat
      real (kind=RKIND), pointer :: rhow
      real (kind=RKIND), pointer :: Cd
      integer, pointer :: nCellsSolve
      integer :: iCell, cell1, cell2, iEdge, iEdgeOnCell
      integer :: err_tmp



      err = 0
      err_tmp = 0

      call mpas_pool_get_config(liConfigs, 'config_subglacial_hydro', config_subglacial_hydro)
      if (.not. config_subglacial_hydro) then
         ! If SGH is not active, skip everything
         return
      endif

      call mpas_pool_get_config(liConfigs, 'config_water_density', rhow)
      call mpas_pool_get_config(liConfigs, 'config_SGH_till_drainage', Cd)

      ! =============
      ! Update till water layer thickness
      ! =============
      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)

         call mpas_pool_get_array(hydroPool, 'tillWaterThickness', Wtill)
         call mpas_pool_get_array(hydroPool, 'deltatSGH', deltatSGH)
         call mpas_pool_get_array(hydroPool, 'meltInput', meltInput)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)

         Wtill = Wtill + deltatSGH * (meltInput / rhow - Cd)
         Wtill = Wtill * li_mask_is_grounded_ice_int(cellMask)  ! zero Wtill in non-grounded locations

         block => block % next
      end do


      ! =============
      ! Calculate edge quantities
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
      !call check_time_step()

      ! =============
      ! Compute advective fluxes
      ! =============
      block => domain % blocklist
      do while (associated(block))

         !call calc_advective(waterFlux, waterVelo, effectiveConducEdge, dcEdge, dvEdge, hydropotentialBase)
         err = ior(err, err_tmp)

         block => block % next
      end do

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
         call mpas_pool_get_dimension(meshPool, 'nCellssSolve', nCellsSolve)
         call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
         call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
         call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
         call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
         call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)

         divergence(:) = 0.0_RKIND  ! zero div before starting
         ! loop over locally owned cells
         do iCell = 1, nCellsSolve
            ! compute fluxes for each edge of the cell
            do iEdgeOnCell = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(iEdgeOnCell, iCell)
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               ! add on advective flux
               divergence(iCell) = divergence(iCell) + waterFlux(iEdge)
               ! add on diffusive flux
               divergence(iCell) = divergence(iCell) + diffusivity(iEdge) * &
                 (waterThickness(cell1) * edgeSignOnCell(iEdge, iCell) - &
                  waterThickness(cell2) * edgeSignOnCell(iEdge, iCell)) / dcEdge(iEdge)
            end do ! edges
         end do ! cells

         block => block % next
      end do


      ! =============
      ! Calculate pressure field
      ! =============
      !call calc_pressure()

      ! =============
      ! Update water layer thickness
      ! =============

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
!  routine li_SGH_final
!
!> \brief   Finalize SGH
!> \author  Matt Hoffman
!> \date    27 June 2016
!> \details
!>  This routine finalizes the subglacial hydro model
!-----------------------------------------------------------------------
   subroutine li_SGH_final(domain, err)

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
      logical, pointer :: config_subglacial_hydro
      type (block_type), pointer :: block
      integer :: err_tmp

      err = 0
      err_tmp = 0

      call mpas_pool_get_config(liConfigs, 'config_subglacial_hydro', config_subglacial_hydro)
      if (.not. config_subglacial_hydro) then
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
          write (stderrUnit,*) "An error has occurred in li_SGH_final."
      endif

   !--------------------------------------------------------------------
   end subroutine li_SGH_final



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

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (block_type), pointer :: block    !< Input/Output: block object

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
      real (kind=RKIND), dimension(:,:), pointer :: cellsOnEdge
      real (kind=RKIND), pointer :: rhow, gravity
      real (kind=RKIND), pointer :: alpha, beta
      real (kind=RKIND), pointer :: conduc_coeff
      character (len=StrKIND), pointer :: config_SGH_tangent_slope_calculation
      integer, pointer :: nEdgesSolve
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

      call mpas_pool_get_config(liConfigs, 'config_water_density', rhow)
      call mpas_pool_get_config(liConfigs, 'config_SGH_alpha', alpha)
      call mpas_pool_get_config(liConfigs, 'config_SGH_beta', beta)
      call mpas_pool_get_config(liConfigs, 'config_SGH_conduc_coeff', conduc_coeff)
      call mpas_pool_get_config(liConfigs, 'config_SGH_tangent_slope_calculation', config_SGH_tangent_slope_calculation)

      call mpas_pool_get_array(hydroPool, 'hydropotentialBase', hydropotentialBase)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
      call mpas_pool_get_array(geometryPool, 'waterPressure', waterPressure)
      call mpas_pool_get_array(hydroPool, 'waterThicknessEdge', waterThicknessEdge)
      call mpas_pool_get_array(hydroPool, 'waterThicknessEdgeUpwind', waterThicknessEdgeUpwind)
      call mpas_pool_get_array(hydroPool, 'hydropotentialBaseSlopeNormal', hydropotentialBaseSlopeNormal)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(hydroPool, 'hydropotentialBaseSlopeTangent', hydropotentialBaseSlopeTangent)
      call mpas_pool_get_array(hydroPool, 'gradMagPhiEdge', gradMagPhiEdge)
      call mpas_pool_get_array(hydroPool, 'effectiveConducEdge', effectiveConducEdge)
      call mpas_pool_get_array(hydroPool, 'diffusivity', diffusivity)
      call mpas_pool_get_array(hydroPool, 'waterVelocity', waterVelocity)
      call mpas_pool_get_array(hydroPool, 'waterFlux', waterFlux)


      ! This is the hydropotential without the water thickness piece
      hydropotentialBase = rhow * gravity * bedTopography + waterPressure

      do iEdge = 1, nEdgesSolve
         cell1 = cellsOnEdge(iEdge, 1)
         cell2 = cellsOnEdge(iEdge, 2)

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
      effectiveConducEdge(:) = conduc_coeff * waterThicknessEdge(:)**(alpha-1.0_RKIND) * gradMagPhiEdge(:)**(beta - 2.0_RKIND)

      ! calculate diffusivity on edges
      diffusivity(:) = rhow * gravity * effectiveConducEdge(:) * waterThicknessEdge(:)

      do iEdge = 1, nEdgesSolve
         cell1 = cellsOnEdge(iEdge, 1)
         cell2 = cellsOnEdge(iEdge, 2)
         waterVelocity(iEdge) = -1.0_RKIND / dcEdge(iEdge) * effectiveConducEdge(iEdge) * &
             (hydropotentialBase(cell1) - hydropotentialBase(cell2))
         velSign = sign(1.0_RKIND, waterVelocity(iEdge))
         waterThicknessEdgeUpwind(iEdge) = max(velSign * waterThickness(cell1),   &
                     velSign * (-1.0_RKIND) * waterThickness(cell2))
         waterFlux(iEdge) = waterVelocity(iEdge) * dvEdge(iEdge) * waterThicknessEdgeUpwind(iEdge)
      end do

   !--------------------------------------------------------------------
   end subroutine calc_edge_quantities




!***********************************************************************
!
!  routine calc_advective_flux
!
!> \brief   Calculate SGH advective fluxes
!> \author  Matt Hoffman
!> \date    5 July 2016
!> \details
!>  This routine calculates needed SGH advective fluxes
!-----------------------------------------------------------------------
!     subroutine calc_advective_flux(block, err)
!  
!        !-----------------------------------------------------------------
!        ! input variables
!        !-----------------------------------------------------------------
!  
!        !-----------------------------------------------------------------
!        ! input/output variables
!        !-----------------------------------------------------------------
!        type (block_type), pointer :: block    !< Input/Output: block object
!  
!        !-----------------------------------------------------------------
!        ! output variables
!        !-----------------------------------------------------------------
!        integer, intent(out) :: err !< Output: error flag
!  
!        !-----------------------------------------------------------------
!        ! local variables
!        !-----------------------------------------------------------------
!        ! Pools pointers
!        type (mpas_pool_type), pointer :: meshPool
!        type (mpas_pool_type), pointer :: geometryPool
!        type (mpas_pool_type), pointer :: hydroPool
!        real (kind=RKIND), dimension(:), pointer :: bedTopography
!        real (kind=RKIND), dimension(:), pointer :: waterPressure
!        real (kind=RKIND), dimension(:), pointer :: hydropotentialBase
!        real (kind=RKIND), dimension(:), pointer :: waterThicknessEdge
!        real (kind=RKIND), dimension(:), pointer :: waterThickness
!        real (kind=RKIND), dimension(:), pointer :: hydropotentialBaseSlopeNormal
!        real (kind=RKIND), dimension(:), pointer :: hydropotentialBaseSlopeTangent
!        real (kind=RKIND), dimension(:), pointer :: gradMagPhiEdge
!        real (kind=RKIND), dimension(:), pointer :: effectiveConducEdge
!        real (kind=RKIND), dimension(:), pointer :: diffusivity
!        real (kind=RKIND), dimension(:), pointer :: dcEdge
!        real (kind=RKIND), dimension(:), pointer :: waterVelocity
!        real (kind=RKIND), dimension(:,:), pointer :: cellsOnEdge
!        real (kind=RKIND), pointer :: rhow, gravity
!        real (kind=RKIND), pointer :: alpha, beta
!        real (kind=RKIND), pointer :: conduc_coeff
!        character (len=StrKIND), pointer :: config_SGH_tangent_slope_calculation
!        integer, pointer :: nEdgesSolve
!        integer :: iEdge, cell1, cell2
!        integer :: err_tmp
!  
!        err = 0
!        err_tmp = 0
!  
!  
!        ! Get pools things
!        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
!        call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
!        call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
!  
!        call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
!  
!        call mpas_pool_get_config(liConfigs, 'config_water_density', rhow)
!        call mpas_pool_get_config(liConfigs, 'config_SGH_alpha', alpha)
!        call mpas_pool_get_config(liConfigs, 'config_SGH_beta', beta)
!        call mpas_pool_get_config(liConfigs, 'config_SGH_conduc_coeff', conduc_coeff)
!        call mpas_pool_get_config(liConfigs, 'config_SGH_tangent_slope_calculation', config_SGH_tangent_slope_calculation)
!  
!        call mpas_pool_get_array(hydroPool, 'hydropotentialBase', hydropotentialBase)
!        call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
!        call mpas_pool_get_array(geometryPool, 'waterPressure', waterPressure)
!        call mpas_pool_get_array(hydroPool, 'waterThicknessEdge', waterThicknessEdge)
!        call mpas_pool_get_array(hydroPool, 'hydropotentialBaseSlopeNormal', hydropotentialBaseSlopeNormal)
!        call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
!        call mpas_pool_get_array(hydroPool, 'hydropotentialBaseSlopeTangent', hydropotentialBaseSlopeTangent)
!        call mpas_pool_get_array(hydroPool, 'gradMagPhiEdge', gradMagPhiEdge)
!        call mpas_pool_get_array(hydroPool, 'effectiveConducEdge', effectiveConducEdge)
!        call mpas_pool_get_array(hydroPool, 'diffusivity', diffusivity)
!        call mpas_pool_get_array(hydroPool, 'waterVelocity', waterVelocity)
!  
!  
!        !-------------
!        select case (trim(config_SGH_advection))
!        !-------------
!        case ("fo")
!  
!           ! Note: This loop structure (nCells loop outside nEdgesOnCell loop) results in double calculation of fluxes
!           !       across each edge. But upwind advection is cheap, so the extra cost is minimal.
!  
!           ! loop over locally owned cells
!           do iCell = 1, nCellsSolve
!              ! compute fluxes for each edge of the cell
!              do iEdgeOnCell = 1, nEdgesOnCell(iCell)
!  
!                 iEdge = edgesOnCell(iEdgeOnCell, iCell)
!                 cell1 = cellsOnEdge(1,iEdge)
!                 cell2 = cellsOnEdge(2,iEdge)
!  
!                 waterFlux = waterVelocity(iEdge) * dvEdge(iEdge) * waterThicknessEdgeUpwind(iEdge)
!              end do ! edges
!           end do ! cells
!  
!        !-------------
!        case ("fct")
!           write (stdErrUnit,*) 'Error: Invalid value for config_SGH_advection.'
!           err = 1
!  
!        !-------------
!        case default
!           write (stdErrUnit,*) 'Error: Invalid value for config_SGH_advection.'
!           err = 1
!        !-------------
!        end select
!        !-------------
!  
!     !--------------------------------------------------------------------
!     end subroutine calc_advective_flux

end module li_SGH_driver
