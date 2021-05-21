! Copyright (c) 2013-2018,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_velocity
!
!> \MPAS land-ice velocity driver
!> \author William Lipscomb
!> \date   10 January 2012
!> \details
!>  This module contains the routines for calling dycores
!>  for calculating ice velocity.
!>
!
!-----------------------------------------------------------------------


module li_velocity

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use mpas_log

   use li_velocity_external
   use li_velocity_simple
   use li_sia
   use li_setup

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

   public :: li_velocity_init, &
             li_velocity_finalize, &
             li_velocity_block_init, &
             li_velocity_solve

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------


!***********************************************************************

contains

!***********************************************************************
!
!  routine li_velocity_init
!
!> \brief   Initializes velocity solver
!> \author  Xylar Asay-Davis
!> \date    18 January 2012
!> \details
!>  This routine initializes the ice velocity solver.
!
!-----------------------------------------------------------------------

   subroutine li_velocity_init(domain, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object

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
      character (len=StrKIND), pointer :: config_velocity_solver
      integer :: err_tmp

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_velocity_solver', config_velocity_solver)

      call mpas_log_write('Using ' // trim(config_velocity_solver) // ' velocity solver.')
      select case (config_velocity_solver)
      case ('none')
         ! Do nothing
      case ('sia')
          call li_sia_init(domain, err)
      case ('simple')
          call li_velocity_simple_init(domain, err)
      case ('L1L2', 'FO', 'Stokes')
          call li_velocity_external_init(domain, err)
      case default
          call mpas_log_write(trim(config_velocity_solver) // ' is not a valid land ice velocity solver option.', MPAS_LOG_ERR)
          err = 1
      end select

      ! === error check
      if (err > 0) then
          call mpas_log_write("An error has occurred in li_velocity_init.", MPAS_LOG_ERR)
      endif

   !--------------------------------------------------------------------
   end subroutine li_velocity_init



!***********************************************************************
!
!  routine li_velocity_block_init
!
!> \brief   Initializes velocity solver
!> \author  William Lipscomb
!> \date    10 January 2012
!> \details
!>  This routine initializes each block of the ice velocity solver.
!
!-----------------------------------------------------------------------
   subroutine li_velocity_block_init(block, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (block_type), intent(inout) :: &
         block          !< Input/Output: block object

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

      character (len=StrKIND), pointer :: config_velocity_solver
      type (mpas_pool_type), pointer :: velocityPool

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_velocity_solver', config_velocity_solver)

      select case (config_velocity_solver)
      case ('none')
         ! Do nothing
      case ('sia')
          call li_sia_block_init(block, err)
      case ('simple')
          call li_velocity_simple_block_init(block, err)
      case ('L1L2', 'FO', 'Stokes')
          call li_velocity_external_block_init(block, err)
      case default
          call mpas_log_write(trim(config_velocity_solver) // ' is not a valid land ice velocity solver option.', MPAS_LOG_ERR)
          err = 1
      end select

      ! === error check
      if (err > 0) then
          call mpas_log_write("An error has occurred in li_velocity_block_init.", MPAS_LOG_ERR)
      endif

   !--------------------------------------------------------------------

   end subroutine li_velocity_block_init



!***********************************************************************
!
!  routine li_velocity_solve
!
!> \brief   Solver for calculating ice velocity
!> \author  William Lipscomb
!> \date    10 January 2012
!> \details
!>  This routine calls velocity solvers.
!
!-----------------------------------------------------------------------
   subroutine li_velocity_solve(domain, solveVelo, err)

      use mpas_vector_reconstruction
      use li_mask
      use li_advection
      use li_thermal

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      logical, intent(in) :: solveVelo !< Input: logical determining if
      !< velocity should be solved, or if an existing solution should be
      !< used to calculate the diagnostic fields related to velocty
      !< (e.g. on a restart)

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object
      ! Note: domain is passed in because halo updates are needed in this routine
      ! and halo updates have to happen outside block loops, which requires domain.

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: thermalPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: velocityPool
      type (mpas_pool_type), pointer :: hydroPool
      ! pointers to get from pools
      character (len=StrKIND), pointer :: config_velocity_solver
      logical, pointer :: config_do_velocity_reconstruction_for_external_dycore
      logical, pointer :: config_print_velocity_cleanup_details
      real (kind=RKIND), pointer :: config_dynamic_thickness
      logical, pointer :: config_adaptive_timestep_include_DCFL
      integer, pointer :: nEdgesSolve
      integer, pointer :: nEdges
      integer, pointer :: nVertInterfaces
      integer, pointer :: nCells
      integer, dimension(:), pointer :: edgeMask, cellMask, vertexMask, vertexMaskOld
      integer, dimension(:,:), pointer :: dirichletVelocityMaskOld, dirichletVelocityMaskNew
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity, normalVelocityInitial
      real (kind=RKIND), dimension(:,:), pointer :: uReconstructX, uReconstructY, uReconstructZ, &
         uReconstructZonal, uReconstructMeridional
      real (kind=RKIND), dimension(:,:), pointer :: edgeNormalVectors
      real (kind=RKIND), dimension(:), pointer :: thickness
      real (kind=RKIND), dimension(:), pointer :: surfaceSpeed, basalSpeed
      real (kind=RKIND), dimension(:), pointer :: xvelmean, yvelmean
      real (kind=RKIND), dimension(:), pointer :: normalSlopeEdge
      real (kind=RKIND), dimension(:), pointer :: dcEdge
      real (kind=RKIND), dimension(:), pointer :: layerInterfaceFractions
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer, dimension(:,:), pointer :: cellsOnVertex
      integer, dimension(:,:), pointer :: verticesOnEdge
      real (kind=RKIND), dimension(:), pointer :: upperSurface
      integer, dimension(:), pointer :: indexToEdgeID
      integer, pointer :: anyDynamicVertexMaskChanged
      integer, pointer :: dirichletMaskChanged
      ! truly local variables
      integer :: cell1, cell2
      integer :: cell3, cell4, thisCell
      integer :: vertex1, vertex2
      integer :: iEdge
      integer :: iCell
      integer :: uphillMarginEdgesFixed
      integer :: err_tmp
      integer :: k
      real (kind=RKIND) :: maxThicknessOnProc, maxThicknessAllProcs
      real (kind=RKIND) :: xVelEdge, yVelEdge
      integer :: blockDynamicVertexMaskChanged, procDynamicVertexMaskChanged
      integer :: blockDirichletMaskChanged, procDirichletMaskChanged

      call mpas_timer_start("velocity solve")

      err_tmp = 0
      err = 0

      call mpas_pool_get_config(liConfigs, 'config_velocity_solver', config_velocity_solver)
      call mpas_pool_get_config(liConfigs, 'config_do_velocity_reconstruction_for_external_dycore', &
         config_do_velocity_reconstruction_for_external_dycore)
      call mpas_pool_get_config(liConfigs, 'config_print_velocity_cleanup_details', config_print_velocity_cleanup_details)
      call mpas_pool_get_config(liConfigs, 'config_dynamic_thickness', config_dynamic_thickness)
      call mpas_pool_get_config(liConfigs, 'config_adaptive_timestep_include_DCFL', config_adaptive_timestep_include_DCFL)

      uphillMarginEdgesFixed = 0


      ! In some situations (first time level of a restart), we don't want to actually solve
      ! velocity again, but we do want to calculate the diagnostic fields related to velocity.
      if (solveVelo) then


      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'thickness')
      call mpas_timer_stop("halo updates")

      ! Update mask just to be safe (might be redundant)
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)

         call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
         call li_update_geometry(geometryPool)

         block => block % next
      end do

      ! Update halos on masks - the outermost cells/edges/vertices may be wrong for mask components that need neighbor information
      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'cellMask')
      call mpas_dmpar_field_halo_exch(domain, 'edgeMask')
      call mpas_dmpar_field_halo_exch(domain, 'vertexMask')
      call mpas_timer_stop("halo updates")

      !   This first block for calculating normalSlopeEdge is also needed for the DCFL calculation
      !   that could occur with any velocity solver
      if ( (trim(config_velocity_solver) == 'sia') .or. (config_adaptive_timestep_include_DCFL) ) then
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
            call mpas_pool_get_array(geometryPool, 'normalSlopeEdge', normalSlopeEdge)
            call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
            call mpas_pool_get_array(geometryPool, 'edgeMask', edgeMask)
            call mpas_pool_get_array(geometryPool, 'upperSurface', upperSurface)
            call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
            call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)

            ! Calculate normal slope
            do iEdge = 1, nEdges
               ! Only calculate slope for edges that have ice on at least one side.
               if ( li_mask_is_dynamic_ice(edgeMask(iEdge)) ) then
                  cell1 = cellsOnEdge(1,iEdge)
                  cell2 = cellsOnEdge(2,iEdge)
                  ! Calculate slope at edge
                  normalSlopeEdge(iEdge) = (upperSurface(cell1) - upperSurface(cell2) ) / dcEdge(iEdge)
               else
                  normalSlopeEdge(iEdge) = 0.0_RKIND
               endif
            end do  ! edges

            block => block % next
         end do
      endif

      if ( (trim(config_velocity_solver) == 'L1L2') .or. &
           (trim(config_velocity_solver) == 'FO') .or. &
           (trim(config_velocity_solver) == 'Stokes') ) then

         ! External solvers may not be able to cope with no ice in the domain, so determine if that is the case
         ! Don't bother checking this with SIA because it requires an extra global reduce
         maxThicknessOnProc = 0.0_RKIND  ! initialize to
         procDynamicVertexMaskChanged = 0
         procDirichletMaskChanged = 0

         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)

            call mpas_pool_get_array(geometryPool, 'thickness', thickness)
            maxThicknessOnProc = max(maxThicknessOnProc, maxval(thickness))

            ! The interface expects an array where 1's are floating edges and 0's are non-floating edges.
            call mpas_pool_get_array(geometryPool, 'edgeMask', edgeMask)
            call mpas_pool_get_array(geometryPool, 'vertexMask', vertexMask, timeLevel=1)

            ! Determine if the vertex mask changed during this time step for this block (needed for external dycores)
            ! TODO:  There may be some aspects of the mask that are ok change for external dycores,
            !        but for now just check the whole thing.
            call mpas_pool_get_array(geometryPool, 'vertexMask', vertexMaskOld, timeLevel=2)
            if ( sum(li_mask_is_dynamic_ice_int(vertexMask) - li_mask_is_dynamic_ice_int(vertexMaskOld)) /= 0 ) then
                blockDynamicVertexMaskChanged = 1
            else
                blockDynamicVertexMaskChanged = 0
            endif
            !print *, 'blockVertexMaskChanged ', blockVertexMaskChanged
            ! Determine if any blocks on this processor had a change to the vertex mask
            procDynamicVertexMaskChanged = max(procDynamicVertexMaskChanged, blockDynamicVertexMaskChanged)
            !print *,'procVertexMaskChanged', procVertexMaskChanged

            ! Also check to see if the Dirichlet b.c. mask has changed
            call mpas_pool_get_array(velocityPool, 'dirichletVelocityMask', dirichletVelocityMaskOld, timeLevel=2)
            call mpas_pool_get_array(velocityPool, 'dirichletVelocityMask', dirichletVelocityMaskNew, timeLevel=1)
            if ( sum(dirichletVelocityMaskNew - dirichletVelocityMaskOld) /= 0 ) then
                blockDirichletMaskChanged = 1
            else
                blockDirichletMaskChanged = 0
            endif
            ! Determine if any blocks on this processor had a change to the vertex mask
            procDirichletMaskChanged = max(procDirichletMaskChanged, blockDirichletMaskChanged)

            ! Set beta for solver to use
            ! (this could potentially be applied to an SIA solver, so it is calculated in this module)
            call calculate_beta(block, err_tmp)
            err = ior(err, err_tmp)


            block => block % next
         end do

         ! check for thickness on any proc
         call mpas_dmpar_max_real(domain % dminfo, maxThicknessOnProc, maxThicknessAllProcs)

         ! Determine if the vertex mask has changed on any processor and store the value for later use
         ! (need to exit the block loop to do so)
         ! TODO Update all blocks with result (if ever support multiple blocks)
         call mpas_pool_get_array(velocityPool, 'anyDynamicVertexMaskChanged', anyDynamicVertexMaskChanged)
         call mpas_dmpar_max_int(domain % dminfo, procDynamicVertexMaskChanged, anyDynamicVertexMaskChanged)
         !print *,'anyDynamicVertexMaskChanged', anyDynamicVertexMaskChanged
         ! Do the same for the Dirichlet b.c. mask
         call mpas_pool_get_array(velocityPool, 'dirichletMaskChanged', dirichletMaskChanged)
         call mpas_dmpar_max_int(domain % dminfo, procDirichletMaskChanged, dirichletMaskChanged)
         !print *,'dirichletMaskChanged', dirichletMaskChanged

      endif



      ! External solvers do not support multiple blocks but the MPAS SIA solver does.
      block => domain % blocklist
      do while (associated(block))
         ! Get variables from pools
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'thermal', thermalPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)

         call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
         call mpas_pool_get_array(velocityPool, 'normalVelocity', normalVelocity)
         call mpas_pool_get_array(geometryPool, 'edgeMask', edgeMask)
         call mpas_pool_get_array(meshPool, 'indexToEdgeID', indexToEdgeID)
         call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
         call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
         call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_array(geometryPool, 'upperSurface', upperSurface)


         ! Solve velocity
         select case (config_velocity_solver)
         case ('none')

            ! Do nothing

         case ('sia')

             call li_sia_solve(meshPool, geometryPool, thermalPool, velocityPool, err_tmp)

         case ('L1L2', 'FO', 'Stokes')

             if (maxThicknessAllProcs < config_dynamic_thickness) then
                ! External dycores may not be able to handle case when there is no ice
                call mpas_pool_get_array(velocityPool, 'uReconstructX', uReconstructX)
                call mpas_pool_get_array(velocityPool, 'uReconstructY', uReconstructY)
                call mpas_pool_get_array(velocityPool, 'uReconstructZ', uReconstructZ)
                normalVelocity = 0.0_RKIND
                uReconstructX = 0.0_RKIND
                uReconstructY = 0.0_RKIND
                uReconstructZ = 0.0_RKIND
                call mpas_log_write("Notice: Skipping velocity solve because there is no dynamic ice in domain.", MPAS_LOG_WARN)
             else
                call li_velocity_external_solve(meshPool, geometryPool, thermalPool, hydroPool, scratchPool, velocityPool, err_tmp)
             endif

         case('simple')

            ! Set the normal velocities to the values computed at initialization
            ! Note: The reason these velocities are reset to the initial values during each time step is that
            !       they may have been altered during the previous time step (e.g., set to zero for non-dynamic edges).
            call mpas_pool_get_array(velocityPool, 'normalVelocityInitial', normalVelocityInitial)
            normalVelocity = normalVelocityInitial

            ! Fix up these velocities by setting them to zero on non-dynamic edges
            do iEdge = 1, nEdgesSolve
               if (.not.(li_mask_is_dynamic_ice(edgeMask(iEdge)))) then
                  normalVelocity(:,iEdge) = 0.0_RKIND
               endif
            enddo

        case default

             call mpas_log_write(trim(config_velocity_solver) // ' is not a valid land ice velocity solver option.', MPAS_LOG_ERR)
             err = 1
             call mpas_timer_stop("velocity solve")
             return

         end select
         err = ior(err, err_tmp)

         ! Some "quality control" of normalVelocity
         do iEdge = 1, nEdgesSolve

            ! Don't allow normalVelocity on edges where an unglaciated cell with
            ! higher elevation neighbors a glaciated cell.  Some velocity solvers
            ! could generate a nonzero velocity on these edges.  In the case of a
            ! velocity directed into the ice sheet, this probably does no harm
            ! for advection because there is no ice to advect in, but it could result
            ! in overly restrictive advective CFL conditions.  In the case of velocity
            ! directed out of the ice sheet, this would result in uphill flow which is
            ! highly unlikely to be physically correct.  (It could be possible in a HO
            ! stress balance where stress transfer 'overrides' the driving stress, but
            ! this seems unlikely to be significant.)  Therefore, always zero velocity
            ! in these situations.
            if ( li_mask_is_dynamic_margin(edgeMask(iEdge)) ) then
               cell1 = cellsOnEdge(1, iEdge)
               cell2 = cellsOnEdge(2, iEdge)
               if ( ( li_mask_is_dynamic_ice(cellMask(cell1)) .and.     &
                      upperSurface(cell2) > upperSurface(cell1) ) .or.  &
                    ( li_mask_is_dynamic_ice(cellMask(cell2)) .and.     &
                      upperSurface(cell1) > upperSurface(cell2) ) ) then
                  if (config_print_velocity_cleanup_details) then
                     call mpas_log_write("Notice: Nonzero velocity has been calculated on an 'uphill' margin edge.  " // &
                        "normalVelocity here has been set to 0.  Location is edge index: $i", intArgs=(/indexToEdgeID(iEdge)/))
                  endif
                  normalVelocity(:, iEdge) = 0.0_RKIND
                  uphillMarginEdgesFixed = uphillMarginEdgesFixed + 1
               endif
            endif

         enddo

         block => block % next
      end do

     if (uphillMarginEdgesFixed > 0) then
          call mpas_log_write("Notice: Nonzero velocity has been calculated on 'uphill' margin edge(s).  normalVelocity has " &
                            // "been set to 0 at these location(s).  Number of edges affected on this processor: $i", &
                            intArgs=(/uphillMarginEdgesFixed/))
      endif

      ! ---
      ! --- update halos on velocity
      ! ---
      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'normalVelocity')
      call mpas_timer_stop("halo updates")

      endif  ! if solveVelo


      ! -- Now that velocity is solved, update diagnostic fields related to velocity --

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'thermal', thermalPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

         ! ---
         ! --- Calculate reconstructed velocities
         ! ---
         ! do this after velocity halo update in case velocities on the 1-halo edge are wrong (depends on velocity solver)
         ! Still do this even if we didn't calculate velocity because on a restart these will be defined at the initial time.
         call mpas_pool_get_dimension(meshPool, 'nVertInterfaces', nVertInterfaces)
         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_array(meshPool, 'layerInterfaceFractions', layerInterfaceFractions)
         call mpas_pool_get_array(velocityPool, 'normalVelocity', normalVelocity)
         call mpas_pool_get_array(velocityPool, 'uReconstructX', uReconstructX)
         call mpas_pool_get_array(velocityPool, 'uReconstructY', uReconstructY)
         call mpas_pool_get_array(velocityPool, 'uReconstructZ', uReconstructZ)
         call mpas_pool_get_array(velocityPool, 'uReconstructZonal', uReconstructZonal)
         call mpas_pool_get_array(velocityPool, 'uReconstructMeridional', uReconstructMeridional)
         call mpas_pool_get_array(velocityPool, 'surfaceSpeed', surfaceSpeed)
         call mpas_pool_get_array(velocityPool, 'basalSpeed', basalSpeed)
         call mpas_pool_get_array(velocityPool, 'xvelmean', xvelmean)
         call mpas_pool_get_array(velocityPool, 'yvelmean', yvelmean)

         ! Velocities need to be reconstructed at cell centers for the native SIA dycore and for prescribed simple velocities.
         ! External dycores return their native velocities at cell center locations,
         ! but these can optionally be overwritten by reconstructed velocities for testing.
         if ( (trim(config_velocity_solver) == 'sia') .or.       &
              (trim(config_velocity_solver) == 'simple') .or.    &
               config_do_velocity_reconstruction_for_external_dycore ) then
            call mpas_reconstruct(meshPool, normalVelocity,               &
                             uReconstructX, uReconstructY, uReconstructZ, &
                             uReconstructZonal, uReconstructMeridional )
         elseif (trim(config_velocity_solver) == 'none') then
            ! The none velocity solver may have data velocity in uReconstructX/Y.
            ! Fill out normalVelocity in case we want to advect with this field
            call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
            call mpas_pool_get_array(meshPool, 'edgeNormalVectors', edgeNormalVectors)
            call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
            do iEdge=1,nEdges
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               do k=1,nVertInterfaces
                 ! average neighboring cell-centered vectors to the edge
                 xVelEdge = 0.5_RKIND*(uReconstructX(k,cell1) + uReconstructX(k,cell2))
                 yVelEdge = 0.5_RKIND*(uReconstructY(k,cell1) + uReconstructY(k,cell2))
                 ! normal component at edge: take dot products with unit vectors at edge
                 normalVelocity(k,iEdge) = xVelEdge * edgeNormalVectors(1,iEdge) + yVelEdge * edgeNormalVectors(2,iEdge)
               enddo
            enddo
         else
            ! For 2-d meshes, these are set by mpas_reconstruct, so set them for HO dycores
            uReconstructZonal = uReconstructX
            uReconstructMeridional = uReconstructY
         end if

         ! ---
         ! --- Remove floating ice where velocity is unrealistically fast
         ! --- (Removed ice is added to calving flux)
         ! --- (TODO: make this a namelist option?)
         ! ---
         call remove_fast_ice(meshPool, uReconstructX, uReconstructY, normalVelocity, geometryPool)

         ! ---
         ! --- Calculate diagnostic speed arrays
         ! ---
         surfaceSpeed = sqrt(uReconstructX(1,:)**2 + uReconstructY(1,:)**2)
         basalSpeed   = sqrt(uReconstructX(nVertInterfaces,:)**2 + uReconstructY(nVertInterfaces,:)**2)

         ! calculate vertical mean velocity
         do iCell = 1, nCells
            xvelmean(iCell) = sum(uReconstructX(:,iCell) * layerInterfaceFractions(:))
            yvelmean(iCell) = sum(uReconstructY(:,iCell) * layerInterfaceFractions(:))
         enddo
         ! ---
         ! --- Calculate strain rates on cell centers
         ! ---
         call calculate_strain_rates_and_stresses(meshPool, geometryPool, thermalPool, scratchPool, velocityPool, err)
         err = ior(err, err_tmp)

         block => block % next
      end do

      ! calculate velocity-related thermal source term fields
      ! (do it here to make restarts easier - this way only the fields need to be restart
      !  variables, rather than the fields that are input to the calculations.)
      if (solveVelo .and. (trim(config_velocity_solver) /= 'none')) then

         ! (on a restart value, do not try to calculate these - instead use the fields in the restart file)
         call li_basal_friction(domain, err_tmp)
         err = ior(err, err_tmp)

         if (trim(config_velocity_solver) == 'sia') then
            ! For FO solver, heatDissipation is passed back from Albany
            call li_heat_dissipation_sia(domain, err_tmp)
            err = ior(err, err_tmp)
         endif
      endif

      ! === error check
      if (err > 0) then
          call mpas_log_write("An error has occurred in li_velocity_solve.", MPAS_LOG_ERR)
      endif

      call mpas_timer_stop("velocity solve")

   !--------------------------------------------------------------------
   end subroutine li_velocity_solve



!***********************************************************************
!
!  routine li_velocity_finalize
!
!> \brief   Finalizes velocity solver
!> \author  Xylar Asay-Davis
!> \date    18 January 2012
!> \details
!>  This routine finalizes the ice velocity solver.
!
!-----------------------------------------------------------------------

   subroutine li_velocity_finalize(domain, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain   !< Input/Output: domain object

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
      character (len=StrKIND), pointer :: config_velocity_solver

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_velocity_solver', config_velocity_solver)

      select case (config_velocity_solver)
      case ('none')
         ! Do nothing
      case ('sia')
          call li_sia_finalize(domain, err)
      case ('simple')
          call li_velocity_simple_finalize(err)
      case ('L1L2', 'FO', 'Stokes')
          call li_velocity_external_finalize(err)
      case default
          call mpas_log_write(trim(config_velocity_solver) // ' is not a valid land ice velocity solver option.', MPAS_LOG_ERR)
          err = 1
          return
      end select

      ! === error check
      if (err > 0) then
          call mpas_log_write("An error has occurred in li_velocity_finalize.", MPAS_LOG_ERR)
      endif

   !--------------------------------------------------------------------
   end subroutine li_velocity_finalize



!***********************************************************************
!  private subroutines
!***********************************************************************


!***********************************************************************
!
!  routine calculate_beta
!
!> \brief  Calulates the beta field to be used by velocity solver
!> \author Matt Hoffman
!> \date   3 October 2013
!> \details
!>  This routine calculates the beta field to be used by the velocity solver.
!>  The resulting field to be used is stored as 'betaSolve'.
!>  Additional methods for calculating beta can be added as they are developed.
!
!-----------------------------------------------------------------------

   subroutine calculate_beta(block, err)

      use li_mask

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (block_type), pointer, intent(inout) :: block

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (mpas_pool_type), pointer :: hydroPool
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: velocityPool
      type (mpas_pool_type), pointer :: thermalPool
      type (mpas_pool_type), pointer :: meshPool
      real (kind=RKIND), dimension(:), pointer :: betaSolve, beta
      real (kind=RKIND), dimension(:), pointer :: basalTemperature, basalPmpTemperature
      integer, dimension(:), pointer :: cellMask
      logical, pointer :: config_use_glp
      logical, pointer :: config_beta_thawed_only
      logical, pointer :: hydroActive
      real (kind=RKIND) :: betaAccum
      integer :: nBetaValues
      integer :: iCell, iCell2
      integer :: neighbor
      integer, pointer :: nCells
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_use_glp', config_use_glp)
      call mpas_pool_get_config(liConfigs, 'config_beta_thawed_only', config_beta_thawed_only)

      call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
      call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
      call mpas_pool_get_subpool(block % structs, 'thermal', thermalPool)
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

      call mpas_pool_get_array(velocityPool, 'betaSolve', betaSolve)
      call mpas_pool_get_array(velocityPool, 'beta', beta)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)

      betaSolve = beta

      if (config_beta_thawed_only) then
         call mpas_pool_get_array(thermalPool, 'basalTemperature', basalTemperature)
         call mpas_pool_get_array(thermalPool, 'basalPmpTemperature', basalPmpTemperature)
         where (basalPmpTemperature - basalTemperature > 0.01_RKIND)
            ! Use a small difference to account for roundoff.  0.01 is value used in thermal module.
            betaSolve = 1.0e6_RKIND  ! no-slip value for frozen bed locations
         end where
      endif

   !--------------------------------------------------------------------
   end subroutine calculate_beta


!***********************************************************************
!
!  routine calculate_strain_rates_and_stresses
!
!> \brief  Calulates strain rates and stresses on cell centers
!> \author Matt Hoffman
!> \date   Feb. 2018
!> \details
!>  This routine calculates strain rates and stresses on cell centers.
!
!-----------------------------------------------------------------------

   subroutine calculate_strain_rates_and_stresses(meshPool, geometryPool, thermalPool, scratchPool, velocityPool, err)

      use li_setup
      use li_diagnostic_vars

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: &
         meshPool                                     !< Input: mesh object

      type (mpas_pool_type), intent(in) :: &
         geometryPool                                 !< Input: geometry object
      type (mpas_pool_type), intent(in) :: &
         thermalPool                                  !< Input: thermal object
      type (mpas_pool_type), intent(in) :: &
         scratchPool                                  !< Input: scratch object

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: &
         velocityPool                                     !< Input: mesh object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      real(kind=RKIND), dimension(:), pointer :: exx, eyy, exy, eTheta, eMax, eMin, dudy, dvdx
      real(kind=RKIND), dimension(:), pointer :: xvelmean, yvelmean
      real(kind=RKIND), dimension(:), pointer :: tauxx, tauyy, tauxy, tauMax, tauMin
      integer :: err_tmp

      real(kind=RKIND), pointer :: config_flowLawExponent

      real (kind=RKIND), dimension(:), pointer :: thickness
      real (kind=RKIND), dimension(:,:), pointer :: temperature
      real (kind=RKIND), dimension(:,:), pointer :: flowParamA
      real (kind=RKIND), dimension(:), pointer :: effectiveViscosity
      real (kind=RKIND), dimension(:), pointer :: stiffnessFactor

      type (field1dReal), pointer :: meanFlowParamAVar
      real (kind=RKIND), dimension(:), pointer :: meanFlowParamA

      integer, pointer :: nVertLevels
      integer, pointer :: nCells

      integer :: iCell

      real(kind=RKIND) :: eEff

      err = 0

      call mpas_pool_get_array(velocityPool, 'xvelmean', xvelmean)
      call mpas_pool_get_array(velocityPool, 'yvelmean', yvelmean)

      call mpas_pool_get_config(liConfigs, 'config_flowLawExponent', config_flowLawExponent)

      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      call mpas_pool_get_array(velocityPool, 'exx', exx)
      call mpas_pool_get_array(velocityPool, 'eyy', eyy)
      call mpas_pool_get_array(velocityPool, 'dudy', dudy)
      call mpas_pool_get_array(velocityPool, 'dvdx', dvdx)
      call mpas_pool_get_array(velocityPool, 'exy', exy)
      call mpas_pool_get_array(velocityPool, 'eTheta', eTheta)
      call mpas_pool_get_array(velocityPool, 'eMax', eMax)
      call mpas_pool_get_array(velocityPool, 'eMin', eMin)

      call mpas_pool_get_array(velocityPool, 'tauxx', tauxx)
      call mpas_pool_get_array(velocityPool, 'tauyy', tauyy)
      call mpas_pool_get_array(velocityPool, 'tauxy', tauxy)
      call mpas_pool_get_array(velocityPool, 'tauMax', tauMax)
      call mpas_pool_get_array(velocityPool, 'tauMin', tauMin)

      call mpas_pool_get_array(thermalPool, 'temperature', temperature)
      call mpas_pool_get_array(velocityPool, 'flowParamA', flowParamA)
      call mpas_pool_get_array(velocityPool, 'stiffnessFactor', stiffnessFactor)

      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(geometryPool, 'effectiveViscosity', effectiveViscosity)

      call mpas_pool_get_field(scratchPool, 'meanFlowParamA', meanFlowParamAVar)
      call mpas_allocate_scratch_field(meanFlowParamAVar, .true.)
      meanFlowParamA => meanFlowParamAVar % array

      ! Calculate strain rates
      call li_compute_gradient_2d(meshPool, xvelmean, exx, dudy, err_tmp)
      err = ior(err, err_tmp)

      call li_compute_gradient_2d(meshPool, yvelmean, dvdx, eyy, err_tmp)
      err = ior(err, err_tmp)

      exy = 0.5_RKIND * (dudy + dvdx)

      ! Calculate principal strain rate angle
      eTheta = 0.5_RKIND * atan( (2.0_RKIND * exy) / (exx - eyy + 1.0e-42_RKIND) )

      ! Calculate principal strain rates
      eMax = 0.5_RKIND * (exx + eyy) + sqrt( (0.5_RKIND * (exx - eyy))**2 + exy**2)
      eMin = 0.5_RKIND * (exx + eyy) - sqrt( (0.5_RKIND * (exx - eyy))**2 + exy**2)

      call li_calculate_flowParamA(meshPool, temperature, thickness, flowParamA, err_tmp)
      err = ior(err, err_tmp)

      ! calculate the depth-averaged flow parameter A
      meanFlowParamA(:) = sum(flowParamA(:,:), dim=1)/REAL(nVertLevels, kind=RKIND)

      ! calculate effective viscosity from strain rate and flow param.
      do iCell = 1, nCells
         eEff = sqrt(exx(iCell)**2 + eyy(iCell)**2 + exx(iCell)*eyy(iCell) + exy(iCell)**2) ! effective strain rate
         if ( (eEff == 0.0_RKIND) .or. (meanFlowParamA(iCell) == 0.0_RKIND) ) then
            effectiveViscosity(iCell) = 0.0_RKIND
         else
            effectiveViscosity(iCell) = &
            0.5_RKIND*stiffnessFactor(iCell)*meanFlowParamA(iCell)**(-1.0_RKIND/config_flowLawExponent)* &
            eEff**((1.0_RKIND-config_flowLawExponent)/config_flowLawExponent)
         endif
      enddo

      ! calculate deviatoric stresses
      tauxx(:) = 2.0_RKIND*effectiveViscosity(:)*exx(:)
      tauyy(:) = 2.0_RKIND*effectiveViscosity(:)*eyy(:)
      tauxy(:) = 2.0_RKIND*effectiveViscosity(:)*exy(:)

      ! calculate principal stresses
      tauMin(:) = (tauxx(:)+tauyy(:))/2.0_RKIND - sqrt(((tauxx(:)-tauyy(:))/2.0_RKIND)**2 + tauxy(:)**2)
      tauMax(:) = (tauxx(:)+tauyy(:))/2.0_RKIND + sqrt(((tauxx(:)-tauyy(:))/2.0_RKIND)**2 + tauxy(:)**2)

      call mpas_deallocate_scratch_field(meanFlowParamAVar, .true.)

   !--------------------------------------------------------------------
   end subroutine calculate_strain_rates_and_stresses

!***********************************************************************
!
!  routine remove_fast_ice
!
!> \brief  Removes ice with a fast surface speed, presumably an iceberg
!> \author Matt Hoffman
!> \date   March 2018
!> \details
!>  Inspired by BISICLES/Dan Martin.  Could/should make an option to disable.
!
!-----------------------------------------------------------------------

   subroutine remove_fast_ice(meshPool, uReconstructX, uReconstructY, normalVelocity, geometryPool)

      use li_mask

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: meshPool

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      real(kind=RKIND), dimension(:,:), intent(inout) :: uReconstructX
      real(kind=RKIND), dimension(:,:), intent(inout) :: uReconstructY
      real(kind=RKIND), dimension(:,:), intent(inout) :: normalVelocity
      type (mpas_pool_type), intent(inout) :: geometryPool

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      real(kind=RKIND), pointer :: config_unrealistic_velocity
      integer, pointer :: nCells
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell
      real(kind=RKIND), dimension(:), pointer :: thickness, calvingThickness
      real(kind=RKIND), dimension(:), pointer :: upperSurface, lowerSurface
      integer, dimension(:), pointer :: cellMask
      integer :: iCell, iEdgeOnCell, iEdge
      integer :: cellCount

      call mpas_pool_get_config(liConfigs, 'config_unrealistic_velocity', config_unrealistic_velocity)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)
      call mpas_pool_get_array(geometryPool, 'upperSurface', upperSurface)
      call mpas_pool_get_array(geometryPool, 'lowerSurface', lowerSurface)


      cellCount = 0
      do iCell = 1, nCells
         if ((li_mask_is_ice(cellMask(iCell))) .and. &
            (sqrt(uReconstructX(1,iCell)**2 + uReconstructY(1,iCell)**2) > config_unrealistic_velocity)) then
         ! "fast" ice that is removed is added to the calving flux
            calvingThickness(iCell) = calvingThickness(iCell) + thickness(iCell)
            thickness(iCell) = 0.0_RKIND
            uReconstructX(:,iCell) = 0.0_RKIND
            uReconstructY(:,iCell) = 0.0_RKIND
            do iEdgeOnCell = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(iEdgeOnCell,iCell)
               normalVelocity(:,iEdge) = 0.0_RKIND
            enddo
            upperSurface(iCell) = 0.0_RKIND
            lowerSurface(iCell) = 0.0_RKIND
            ! TODO Should also zero temperature/tracers
            cellCount = cellCount + 1
         endif
      end do

      if (cellCount > 0) then
         call mpas_log_write("Removed $i cells due to speeds greater than $r m/s.", &
            intArgs=(/cellCount/), realArgs=(/config_unrealistic_velocity/))
      endif

   !--------------------------------------------------------------------
   end subroutine remove_fast_ice

!***********************************************************************

end module li_velocity

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
