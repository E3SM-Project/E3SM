! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_ssh_and_landIcePressure
!
!> \brief MPAS ocean initialize matching SSH and land-ice pressure
!> \author Xylar Asay-Davis
!> \date   06/05/2015
!> \details
!>  This module contains the routines for aiding in initializing the
!>  land-ice pressure based on the sea-surface height (SSH)
!>  so that the barotropic pressure-gradient force (PGF) is initially small
!
!-----------------------------------------------------------------------

module ocn_init_ssh_and_landIcePressure

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants

   use ocn_constants
   use ocn_config
   use ocn_diagnostics_variables

   use ocn_equation_of_state

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

   public :: ocn_init_ssh_and_landIcePressure_balance


   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_ssh_and_landIcePressure_balance
!
!> \brief   Compute the balance land-ice pressure given the SSH
!> \author  Xylar Asay-Davis
!> \date    8/8/2016
!> \details
!>  This routine produces an initial guess at land-ice pressure using the
!>  density of the topmost layer as the effective density of seawater within
!>  the land ice. The resulting land-ice pressure is approximately consistent
!>  with the ice draft in the sense that the horizontal pressure-gradient force (HPGF)
!>  should be small at the ocean surface.
!>  ocn_init_vertical_grid() should be called to produce the appropriate
!>  vertical grid before calling this subroutine and activeTracers should be
!>  initialized based on that vertical grid, so that density can be computed
!>  here.

!-----------------------------------------------------------------------

   subroutine ocn_init_ssh_and_landIcePressure_balance(domain, iErr)!{{{

   !--------------------------------------------------------------------

     type (domain_type), intent(inout) :: domain
     integer, intent(out) :: iErr

     type (block_type), pointer :: block_ptr

     type (mpas_pool_type), pointer :: meshPool, forcingPool, statePool, &
                                       verticalMeshPool

     type (mpas_pool_type), pointer :: tracersPool

     real (kind=RKIND), dimension(:), pointer :: ssh

     real (kind=RKIND), dimension(:), pointer :: landIcePressure, landIceDraft, &
                                                 effectiveDensityInLandIce

     integer, dimension(:), pointer :: minLevelCell

     integer, pointer :: nCells

     integer :: iCell

     iErr = 0

     ! compute density

     block_ptr => domain % blocklist
     do while(associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)

       call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

       call ocn_equation_of_state_density(statePool, meshPool, tracersSurfaceValue, &
           nCells, 0, 'relative', density, iErr, &
           timeLevelIn=1)

       if(iErr .ne. 0) then
         call mpas_log_write( 'ocn_equation_of_state_density failed.', MPAS_LOG_CRIT)
         return
       end if

       block_ptr => block_ptr % next
     end do !block_ptr

     block_ptr => domain % blocklist
     do while(associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
       call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
       call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

       call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

       call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)

       call mpas_pool_get_array(forcingPool, 'landIcePressure', landIcePressure)

       call mpas_pool_get_array(forcingPool, 'landIceDraft', landIceDraft)

       call mpas_pool_get_array(statePool, 'ssh', ssh, 1)
       call mpas_pool_get_array(statePool, 'effectiveDensityInLandIce', effectiveDensityInLandIce, 1)

       call ocn_equation_of_state_density(statePool, meshPool, tracersSurfaceValue, &
           nCells, 0, 'relative', density, iErr, &
           timeLevelIn=1)

       if(iErr .ne. 0) then
         call mpas_log_write( 'ocn_equation_of_state_density failed.', MPAS_LOG_CRIT)
         return
       end if

       do iCell = 1, nCells
         if(sshAdjustmentMask(iCell) == 0) then
           ssh(iCell) = 0.0_RKIND
           landIcePressure(iCell) = 0.0_RKIND

           if (associated(effectiveDensityInLandIce)) &
             ! effective density cannot be determined
             effectiveDensityInLandIce(iCell) = 0.0_RKIND
           cycle
         end if

         landIcePressure(iCell) = &
           max(0.0_RKIND, -density(minLevelCell(iCell), iCell)*gravity*ssh(iCell))
         if (associated(effectiveDensityInLandIce)) &
           effectiveDensityInLandIce(iCell) = density(minLevelCell(iCell), iCell)
       end do

       ! copy the SSH into the landIceDraft so we can use it later to remove it when
       ! computing sea-surface tilt
       landIceDraft(:) = ssh(:)

       block_ptr => block_ptr % next
     end do !block_ptr

   !--------------------------------------------------------------------

   end subroutine ocn_init_ssh_and_landIcePressure_balance

!***********************************************************************

end module ocn_init_ssh_and_landIcePressure

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
