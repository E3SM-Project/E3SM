! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
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

      write(stdoutUnit,*) 'Using ', trim(config_velocity_solver), ' dynamical core.'
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
          write(stderrUnit,*) trim(config_velocity_solver), ' is not a valid land ice velocity solver option.'
          err = 1
      end select

      ! === error check
      if (err > 0) then
          write (stderrUnit,*) "An error has occurred in li_velocity_init."
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
          write(stderrUnit,*) trim(config_velocity_solver), ' is not a valid land ice velocity solver option.'
          err = 1
      end select

      ! === error check
      if (err > 0) then
          write (stderrUnit,*) "An error has occurred in li_velocity_block_init."
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
      ! pointers to get from pools
      character (len=StrKIND), pointer :: config_velocity_solver
      logical, pointer :: config_do_velocity_reconstruction_for_external_dycore
      logical, pointer :: config_print_velocity_cleanup_details
      real (kind=RKIND), pointer :: config_dynamic_thickness
      integer, pointer :: nEdgesSolve
      integer, pointer :: nEdges
      integer, pointer :: nVertInterfaces
      integer, dimension(:), pointer :: edgeMask, cellMask
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity, normalVelocityInitial
      real (kind=RKIND), dimension(:,:), pointer :: uReconstructX, uReconstructY, uReconstructZ, uReconstructZonal, uReconstructMeridional
      real (kind=RKIND), dimension(:,:), pointer :: edgeNormalVectors
      real (kind=RKIND), dimension(:), pointer :: thickness
      real (kind=RKIND), dimension(:), pointer :: surfaceSpeed, basalSpeed
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer, dimension(:,:), pointer :: cellsOnVertex
      integer, dimension(:,:), pointer :: verticesOnEdge
      real (kind=RKIND), dimension(:), pointer :: upperSurface
      integer, dimension(:), pointer :: indexToEdgeID
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

      call mpas_timer_start("velocity solve")

      err_tmp = 0
      err = 0

      call mpas_pool_get_config(liConfigs, 'config_velocity_solver', config_velocity_solver)
      call mpas_pool_get_config(liConfigs, 'config_do_velocity_reconstruction_for_external_dycore', config_do_velocity_reconstruction_for_external_dycore)
      call mpas_pool_get_config(liConfigs, 'config_print_velocity_cleanup_details', config_print_velocity_cleanup_details)
      call mpas_pool_get_config(liConfigs, 'config_dynamic_thickness', config_dynamic_thickness)

      uphillMarginEdgesFixed = 0


      ! In some situations (first time level of a restart), we don't want to actually solve
      ! velocity again, but we do want to calculate the diagnostic fields related to velocity.
      if (solveVelo) then

      ! External solvers may not be able to cope with no ice in the domain, so determine if that is the case
      ! Don't bother checking this with SIA because it requires an extra global reduce
      if ( (trim(config_velocity_solver) == 'L1L2') .or. (trim(config_velocity_solver) == 'FO') .or. &
           (trim(config_velocity_solver) == 'Stokes') ) then
         maxThicknessOnProc = 0.0  ! initialize to 0
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
            call mpas_pool_get_array(geometryPool, 'thickness', thickness)
            maxThicknessOnProc = max(maxThicknessOnProc, maxval(thickness))

            block => block % next
         end do

         call mpas_dmpar_max_real(domain % dminfo, maxThicknessOnProc, maxThicknessAllProcs)
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

             call li_sia_solve(meshPool, geometryPool, velocityPool, err_tmp)

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
                write (stderrUnit,*) "Notice: Skipping velocity solve because there is no dynamic ice in domain."
             else
                call li_velocity_external_solve(meshPool, geometryPool, thermalPool, scratchPool, velocityPool, err_tmp)
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

             write(stderrUnit,*) 'Error: ', trim(config_velocity_solver), ' is not a valid land ice velocity solver option.'
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
                     write (stderrUnit,*) "Notice: Nonzero velocity has been calculated on an 'uphill' margin edge.  normalVelocity here has been set to 0.  Location is edge index:", indexToEdgeID(iEdge)
                  endif
                  normalVelocity(:, iEdge) = 0.0_RKIND
                  uphillMarginEdgesFixed = uphillMarginEdgesFixed + 1
               endif
            endif

         enddo

         block => block % next
      end do

     if (uphillMarginEdgesFixed > 0) then
          write (stderrUnit,*) "Notice: Nonzero velocity has been calculated on 'uphill' margin edge(s).  normalVelocity has " &
                            // "been set to 0 at these location(s).  Number of edges affected on this processor:", &
                            uphillMarginEdgesFixed
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

         ! ---
         ! --- Calculate reconstructed velocities
         ! ---
         ! do this after velocity halo update in case velocities on the 1-halo edge are wrong (depends on velocity solver)
         ! Still do this even if we didn't calculate velocity because on a restart these will be defined at the initial time.
         call mpas_pool_get_dimension(meshPool, 'nVertInterfaces', nVertInterfaces)
         call mpas_pool_get_array(velocityPool, 'normalVelocity', normalVelocity)
         call mpas_pool_get_array(velocityPool, 'uReconstructX', uReconstructX)
         call mpas_pool_get_array(velocityPool, 'uReconstructY', uReconstructY)
         call mpas_pool_get_array(velocityPool, 'uReconstructZ', uReconstructZ)
         call mpas_pool_get_array(velocityPool, 'uReconstructZonal', uReconstructZonal)
         call mpas_pool_get_array(velocityPool, 'uReconstructMeridional', uReconstructMeridional)
         call mpas_pool_get_array(velocityPool, 'surfaceSpeed', surfaceSpeed)
         call mpas_pool_get_array(velocityPool, 'basalSpeed', basalSpeed)

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
            ! The none velocity solver may have data velocity in uReconstructX/Y.  Fill out normalVelocity in case we want to advect with this field
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
         ! --- Calculate diagnostic speed arrays
         ! ---
         surfaceSpeed = sqrt(uReconstructX(1,:)**2 + uReconstructY(1,:)**2)
         basalSpeed   = sqrt(uReconstructX(nVertInterfaces,:)**2 + uReconstructY(nVertInterfaces,:)**2)

         block => block % next
      end do

      ! === error check
      if (err > 0) then
          write (stderrUnit,*) "An error has occurred in li_velocity_solve."
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
          write(stdoutUnit,*) trim(config_velocity_solver), ' is not a valid land ice velocity solver option.'
          err = 1
          return
      end select

      ! === error check
      if (err > 0) then
          write (stderrUnit,*) "An error has occurred in li_velocity_finalize."
      endif

   !--------------------------------------------------------------------
   end subroutine li_velocity_finalize



!***********************************************************************

end module li_velocity

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
