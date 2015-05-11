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
   use mpas_configure
   use li_velocity_external
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
   subroutine li_velocity_solve(domain, err)

      use mpas_vector_reconstruction
      use li_mask

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

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
      type (mpas_pool_type), pointer :: velocityPool
      ! pointers to get from pools
      character (len=StrKIND), pointer :: config_velocity_solver
      logical, pointer :: config_do_velocity_reconstruction_for_external_dycore
      integer, pointer :: nEdgesSolve
      integer, pointer :: nVertInterfaces
      integer, dimension(:), pointer :: edgeMask, cellMask
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity, uReconstructX, uReconstructY, uReconstructZ, uReconstructZonal, uReconstructMeridional
      real (kind=RKIND), dimension(:), pointer :: surfaceSpeed, basalSpeed
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer, dimension(:,:), pointer :: cellsOnVertex
      integer, dimension(:,:), pointer :: verticesOnEdge
      real (kind=RKIND), dimension(:), pointer :: upperSurface
      type (field2DReal), pointer :: normalVelocityField
      integer, dimension(:), pointer :: indexToEdgeID
      ! truly local variables
      integer :: cell1, cell2
      integer :: cell3, cell4, thisCell
      integer :: vertex1, vertex2
      integer :: iEdge
      integer :: iCell
      integer :: err_tmp

      call mpas_timer_start("velocity solve")

      err_tmp = 0
      err = 0

      call mpas_pool_get_config(liConfigs, 'config_velocity_solver', config_velocity_solver)
      call mpas_pool_get_config(liConfigs, 'config_do_velocity_reconstruction_for_external_dycore', config_do_velocity_reconstruction_for_external_dycore)

      ! External solvers do not support multiple blocks but the MPAS SIA solver does.
      block => domain % blocklist
      do while (associated(block))
         ! Get variables from pools
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'thermal', thermalPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call mpas_pool_get_config(liConfigs, 'config_velocity_solver', config_velocity_solver)


         ! Solve velocity
         select case (config_velocity_solver)
         case ('none')
            ! Do nothing
         case ('sia')
             call li_sia_solve(meshPool, geometryPool, velocityPool, err_tmp)
         case ('L1L2', 'FO', 'Stokes')
             call li_velocity_external_solve(meshPool, geometryPool, thermalPool, velocityPool, err_tmp)
         case default
             write(stderrUnit,*) 'Error: ', trim(config_velocity_solver), ' is not a valid land ice velocity solver option.'
             err = 1
             call mpas_timer_stop("velocity solve")
             return
         end select
         err = ior(err, err_tmp)


         ! Some "quality control" of normalVelocity
         call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
         call mpas_pool_get_array(velocityPool, 'normalVelocity', normalVelocity)
         call mpas_pool_get_array(geometryPool, 'edgeMask', edgeMask)
         call mpas_pool_get_array(meshPool, 'indexToEdgeID', indexToEdgeID)
         call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
         call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
         call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_array(geometryPool, 'upperSurface', upperSurface)

         do iEdge = 1, nEdgesSolve

            ! Check if the velocity solver has returned a velocity on any non-dynamic edges
            if ( li_mask_is_ice(edgeMask(iEdge)) .and.                    &
                (.not.(li_mask_is_dynamic_ice(edgeMask(iEdge)))) .and.   &
                (maxval(abs(normalVelocity(:,iEdge))) /= 0.0_RKIND)      &
               ) then
               ! There is an edge case where this is ok.  If there are two peninsulas of dynamic ice
               ! with a single 'row' of nondynamic cells between them, the FEM velo solver will likely
               ! calculate a nonzero velocity on an edge that has 0 thickness.  The two FEM elements
               ! neighboring this edge have nonzero thickness everywhere except along this edge, and
               ! so there is no guarantee of zero-velocity on this edge.  Schematically, this looks like:
               !
               !    \   I   /
               !   A |--e--| A
               !    /   I   \
               !
               ! where the lines are edges, and e is the edge with the issue.  I's are inactive cells, and
               ! A's are active cells.  So check for this specific situation before calling this an error.
               cell1 = cellsOnEdge(1, iEdge)
               cell2 = cellsOnEdge(2, iEdge)
               ! Criterion 1: both cells adjacent to edge are inactive
               if ( ( .not. li_mask_is_dynamic_ice(cellMask(cell1)) )   .and.     &
                    ( .not. li_mask_is_dynamic_ice(cellMask(cell2)) ) ) then
                  ! Criterion 2: both remaining cells adjacent to edge's vertices are active
                  vertex1 = verticesOnEdge(1, iEdge)
                  cell3 = -999
                  do iCell = 1, 3
                     thisCell = cellsOnVertex(iCell, vertex1)
                     if ((thisCell /= cell1) .and. (thisCell /= cell2)) then
                        cell3 = thisCell
                        exit  ! we found the remaining cell on the vertex
                     endif
                  enddo
                  vertex2 = verticesOnEdge(2, iEdge)
                  cell4 = -999
                  do iCell = 1, 3
                     thisCell = cellsOnVertex(iCell, vertex2)
                     if ((thisCell /= cell1) .and. (thisCell /= cell2)) then
                        cell4 = thisCell
                        exit  ! we found the remaining cell on the vertex
                     endif
                  enddo
                  if ( (li_mask_is_dynamic_ice(cellMask(cell3))) .and.     &
                       (li_mask_is_dynamic_ice(cellMask(cell4))) ) then
                     write (stderrUnit,*) "Notice: External velocity solver returned a normalVelocity on a non-dynamic edge, but this is ok because the location is in a non-dynamic 'inlet'.  normalVelocity has been set to 0 at this location.  Location is edge index:", indexToEdgeID(iEdge)
                     normalVelocity(:,iEdge) = 0.0_RKIND
                  else
                     write (stderrUnit,*) 'ERROR: VELO ON NON-DYNAMIC EDGE, edge=', indexToEdgeID(iEdge)
                     err_tmp= 1
                     !!!normalVelocity(:,iEdge) = 0.0_RKIND  ! a hack to ignore this error.
                  endif  ! Criterion 2 check
               endif  ! Criterion 1 check
            endif

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
               cell2 = cellsOnEdge(1, iEdge)
               if ( ( li_mask_is_dynamic_ice(cellMask(cell1)) .and.     &
                      upperSurface(cell2) > upperSurface(cell1) ) .or.  &
                    ( li_mask_is_dynamic_ice(cellMask(cell2)) .and.     &
                      upperSurface(cell1) > upperSurface(cell2) ) ) then
                  normalVelocity(:, iEdge) = 0.0_RKIND
               endif
            endif

         enddo

         if (err_tmp == 1) then
            write(stderrUnit,*) 'Error: Velocity has been calculated on non-dynamic edges.  There is a problem with the velocity solver.' !!!  Velocity on those edges have been set to 0, but this should be a fatal error.'
            err = 1
         end if

         block => block % next
      end do


      ! ---
      ! --- update halos on velocity
      ! ---
      call mpas_timer_start("halo updates")
      call mpas_pool_get_subpool(domain % blocklist % structs, 'velocity', velocityPool)
      call mpas_pool_get_field(velocityPool, 'normalVelocity', normalVelocityField)
      call mpas_dmpar_exch_halo_field(normalVelocityField)
      call mpas_timer_stop("halo updates")


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

         ! Native SIA dycore needs to have reconstructed velocities calculated.
         ! External dycores return their native velocities at cell center locations,
         ! but these can optionally be overwritten by reconstructed velocities for testing.
         if ( (trim(config_velocity_solver) == 'sia') .or.    &
               config_do_velocity_reconstruction_for_external_dycore ) then
            call mpas_reconstruct(meshPool, normalVelocity,               &
                             uReconstructX, uReconstructY, uReconstructZ, &
                             uReconstructZonal, uReconstructMeridional )
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
