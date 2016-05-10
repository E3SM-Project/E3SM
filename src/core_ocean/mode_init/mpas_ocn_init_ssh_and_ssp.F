! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_ssh_and_ssp
!
!> \brief MPAS ocean initialize matching SSH and SSP
!> \author Xylar Asay-Davis
!> \date   06/05/2015
!> \details
!>  This module contains the routines for aiding in initializing the
!>  sea-surface pressure (SSP) based on the sea-surface height (SSH)
!>  so that the barotropic pressure-gradient force (PGF) is initially small
!
!-----------------------------------------------------------------------

module ocn_init_ssh_and_ssp

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants

   use ocn_constants
   use ocn_init_interpolation
   use ocn_init_vertical_grids

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

   public :: ocn_init_ssh_and_ssp_vertical_grid, &
             ocn_init_ssh_and_ssp_balance


   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains


!***********************************************************************
!
!  routine ocn_init_ssh_and_ssp_vertical_grid
!
!> \brief   Initialize z* vertical grid based on SSH
!> \author  Xylar Asay-Davis
!> \date    10/21/2015
!> \details
!>  This routine sets up the vertical grid (layerThickness,
!>  zMid and restingThickness) needed for computing SSH from
!>  SSP or visa versa. bottomDepth, refBottomDepth, maxLevelCell
!>  and modifySSHMask must have been computed by the test case
!>  before calling this routine.  If
!>  config_iterative_init_variable = 'ssp', the test
!>  case must compute refSSH before calling this routine.
!>  modifySSHMask should be set to 1 wherever the ssh or ssp
!>  should be modified for consistency (e.g. under land ice). This
!>  routine will take care of setting up partial bottom cells
!>  by calling ocn_alter_bottomDepth_for_pbcs (except for the
!>  Haney-number-constrained coordinate, which handle thin bottom
!>  cells via the Haney-number constraint.

!-----------------------------------------------------------------------

   subroutine ocn_init_ssh_and_ssp_vertical_grid(domain, iErr)!{{{

   !--------------------------------------------------------------------

     type (domain_type), intent(inout) :: domain
     integer, intent(out) :: iErr

   !--------------------------------------------------------------------

     iErr = 0
     call ocn_init_vertical_grid(domain, iErr, updateOnly=.false.)

   end subroutine ocn_init_ssh_and_ssp_vertical_grid

!***********************************************************************
!
!  routine ocn_init_ssh_and_ssp_balance
!
!> \brief   Compute the balance SSP given the SSH or visa versa
!> \author  Xylar Asay-Davis
!> \date    10/21/2015
!> \details
!>  This routine either updates SSH based on SSP (if config_iterative_init_variable = 'ssh')
!>  or visa versa (if config_iterative_init_variable = 'ssp').  The routine either produces
!>  and initial guess at SSP or SSH (if config_read_ssh_and_ssp_from_stream = .false.)
!>  or it updates SSP or SSH based on the change in SSH over a forward run
!>  (if config_read_ssh_and_ssp_from_stream = .true.).
!>  The SSP and SSH are approximately consistent with one another
!>  in the sense that the horizontal pressure-gradient force (HPGF)
!>  should be small at the ocean surface.
!>  ocn_init_ssh_and_ssp_vertical_grid should be called to produce
!>  the appropriate vertical grid before calling this subroutine.
!>  activeTracers should be initialized based on this vertical grid.
!>  Upon completion, the vertical grid will have been updated (if necessary)
!>  to be consistent with the SSH, in which case the activeTracers will have been
!>  interpolated to the new grid.

!-----------------------------------------------------------------------

   subroutine ocn_init_ssh_and_ssp_balance(domain, iErr)!{{{

   !--------------------------------------------------------------------

     type (domain_type), intent(inout) :: domain
     integer, intent(out) :: iErr

     logical, pointer :: config_read_ssh_and_ssp_from_stream

     character (len=StrKIND), pointer :: config_iterative_init_variable

   !--------------------------------------------------------------------

     iErr = 0

     call mpas_pool_get_config(ocnConfigs, 'config_read_ssh_and_ssp_from_stream', config_read_ssh_and_ssp_from_stream)
     call mpas_pool_get_config(ocnConfigs, 'config_iterative_init_variable', config_iterative_init_variable)


     if(config_read_ssh_and_ssp_from_stream) then
       if(config_iterative_init_variable == 'ssp') then
         call ocn_update_ssp_from_ssh(domain, iErr)

         if(iErr .ne. 0) then
           write(stderrUnit,*) 'ERROR: ocn_update_ssp_from_ssh failed.'
           return
         end if
       end if
     else
       ! In this case, also recompute activeTracers, zMid and layerThickness taking
       ! the ssh into account
       call initial_guess_ssp_ssh(domain, iErr)

       if(iErr .ne. 0) then
         write(stderrUnit,*) 'ERROR: initial_guess_ssp_ssh failed.'
         return
       end if
     end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_ssh_and_ssp_balance

!***********************************************************************
!
! PRIVATE SUBROUTINES
!
!***********************************************************************

!***********************************************************************
!
!  routine ocn_init_vertical_grid
!
!> \brief   Initialize z* vertical grid based on SSH
!> \author  Xylar Asay-Davis
!> \date    10/21/2015
!> \details
!>  This routine sets up the vertical grid (layerThickness,
!>  zMid and restingThickness) needed for computing SSH from
!>  SSP or visa versa. bottomDepth, refBottomDepth and maxLevelCell
!>  must have been computed by the test case before calling this
!>  routine.  If config_iterative_init_variable = 'ssp', the test
!>  case must compute refSSH before calling this routine.  This
!>  routine will take care of setting up partial bottom cells
!>  by calling ocn_alter_bottomDepth_for_pbcs (except for the
!>  Haney-number-constrained coordinate, which handle thin bottom
!>  cells via the Haney-number constraint.
!-----------------------------------------------------------------------

   subroutine ocn_init_vertical_grid(domain, iErr, updateOnly)!{{{

   !--------------------------------------------------------------------

     type (domain_type), intent(inout) :: domain
     integer, intent(out) :: iErr
     logical, intent(in) :: updateOnly

     type (block_type), pointer :: block_ptr

     type (mpas_pool_type), pointer :: meshPool, statePool, diagnosticsPool, verticalMeshPool

     logical, pointer :: config_read_ssh_and_ssp_from_stream , &
                         config_modify_open_ocean_ssh, &
                         config_use_rx1_constraint

     character (len=StrKIND), pointer :: config_iterative_init_variable

     ! Define dimension pointers
     integer, pointer :: nCells, nVertLevels

     ! Define variable pointers
     integer, dimension(:), pointer :: maxLevelCell, modifySSHMask
     real (kind=RKIND), dimension(:), pointer :: refBottomDepth, bottomDepth, deltaSSH, &
                                                 refSSH, ssh
     real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness, zMid

     integer :: iCell

     logical :: initWithSSH, initZStarWithSSH, initZStarWithoutSSH, initRx1WithSSH

   !--------------------------------------------------------------------

     iErr = 0

     call mpas_pool_get_config(ocnConfigs, 'config_read_ssh_and_ssp_from_stream', config_read_ssh_and_ssp_from_stream)
     call mpas_pool_get_config(ocnConfigs, 'config_iterative_init_variable', config_iterative_init_variable)

     if(config_iterative_init_variable .ne. 'ssh' &
        .and. config_iterative_init_variable .ne. 'ssp') then
       iErr = 1
       write(stderrUnit,*) 'ERROR: invalid value for config_iterative_init_variable', trim(config_iterative_init_variable)
       return
     end if
     call mpas_pool_get_config(ocnConfigs, 'config_use_rx1_constraint', config_use_rx1_constraint)

     initWithSSH = updateOnly .or. config_read_ssh_and_ssp_from_stream
     initZStarWithSSH = initWithSSH .and. .not. config_use_rx1_constraint
     initZStarWithoutSSH = .not. initWithSSH
     initRx1WithSSH = initWithSSH .and. config_use_rx1_constraint

     if(.not. updateOnly) then
       ! we haven't computed deltaSSH yet and may need to modify ssh

       block_ptr => domain % blocklist
       do while(associated(block_ptr))
         call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)

         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

         call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
         call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

         call mpas_pool_get_array(statePool, 'ssh', ssh, 1)

         call mpas_pool_get_array(diagnosticsPool, 'refSSH', refSSH)
         call mpas_pool_get_array(diagnosticsPool, 'deltaSSH', deltaSSH)
         call mpas_pool_get_array(diagnosticsPool, 'modifySSHMask', modifySSHMask)

         if(config_read_ssh_and_ssp_from_stream) then
           where(modifySSHMask(:) == 1)
             deltaSSH(:) = ssh(:) - refSSH(:)
           elsewhere
             deltaSSH(:) = 0.0_RKIND
           end where
         end if
         if(config_iterative_init_variable == 'ssp') then
           ! we're updating SSP so we want to keep ssh fixed at refSSH
           ssh(:) = refSSH(:)
         end if

         do iCell = 1, nCells
           if(modifySSHMask(iCell) == 0) then
             ! we don't want the SSH to be modified here, so set it back to zero
             ssh(iCell) = 0.0_RKIND
           end if
         end do !iCell

         !pbcs are handled internally as part of computing the Haney-number-constrained coordinate
         if(.not. initRx1WithSSH) then
           do iCell = 1, nCells
             call ocn_alter_bottomDepth_for_pbcs(bottomDepth(iCell), refBottomDepth, maxLevelCell(iCell), iErr)
             if(iErr .ne. 0) then
               write(stderrUnit,*) 'ERROR: ocn_alter_bottomDepth_for_pbcs failed.'
               return
             end if
           end do
         end if
         block_ptr => block_ptr % next
       end do !block_ptr
     end if

     block_ptr => domain % blocklist
     do while(associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
       call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)

       call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
       call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

       call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
       call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
       call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

       call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
       call mpas_pool_get_array(statePool, 'ssh', ssh, 1)

       call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)

       call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
       call mpas_pool_get_array(diagnosticsPool, 'refSSH', refSSH)
       call mpas_pool_get_array(diagnosticsPool, 'deltaSSH', deltaSSH)
       call mpas_pool_get_array(diagnosticsPool, 'modifySSHMask', modifySSHMask)

       do iCell = 1, nCells
         if(initZStarWithoutSSH) then
           ! We don't know the ssh or ssp yet, and we need tracers on a reference grid to figure it out.
           ! compute restingThickness and reference layerThickness and zMid based on topography with ssh=0
           ! (omitting ssh argument)
           call ocn_compute_layerThickness_zMid_from_bottomDepth(layerThickness(:,iCell),zMid(:,iCell), &
                refBottomDepth,bottomDepth(iCell), &
                maxLevelCell(iCell),nVertLevels,iErr, &
                restingThickness=restingThickness(:,iCell))
         else if(initZStarWithSSH) then
           ! we already know the ssh and ssp we want to use.
           ! compute the layer thicknesses and zMid based on topography and ssh
           call ocn_compute_layerThickness_zMid_from_bottomDepth(layerThickness(:,iCell),zMid(:,iCell), &
                refBottomDepth,bottomDepth(iCell), &
                maxLevelCell(iCell),nVertLevels,iErr, &
                restingThickness=restingThickness(:,iCell), &
                ssh=ssh(iCell))
         end if

         if(iErr .ne. 0) then
           write(stderrUnit,*) 'ERROR: ocn_compute_layerThickness_zMid_from_bottomDepth failed.'
           return
         end if
       end do !iCell

       block_ptr => block_ptr % next
     end do !block_ptr

     if(initRx1WithSSH) then
       ! We already know the ssh and ssp we want to use.
       ! Compute the layer thicknesses and zMid based on topography and ssh.
       ! Use rx1 constraint to recompute the vertical grid.
       call ocn_init_vertical_grid_with_max_rx1(domain, iErr)

       if(iErr .ne. 0) then
         write(stderrUnit,*) 'ERROR: ocn_init_vertical_grid_with_max_rx1 failed.'
         return
       end if

     end if

   end subroutine ocn_init_vertical_grid

!***********************************************************************
!
!  routine initial_guess_ssp_ssh
!
!> \brief   Compute the balance SSP given the SSH or visa versa
!> \author  Xylar Asay-Davis
!> \date    10/12/2015
!> \details
!>  This routine computes the SSH that is approximately consistent
!>  with a given SSP or visa versa, given reference temperature and
!>  salinity fields. activeTracers should be initialized  to reference
!>  T and S and will contain T and S interpolated at zMid on completion.
!>  This subroutine assumes that zMid and/or layerThickness are needed
!>  from computing the reference T and S and have already been computed
!>  with ocn_compute_layerThickness_zMid_from_bottomDepth (omitting the
!>  ssh argument) before calling this routine.  zMid and layerThickness
!>  are recomputed taking the SSH into account in this routine.
!>  The density and bottom pressure are also computed by this routine.

!-----------------------------------------------------------------------

   subroutine initial_guess_ssp_ssh(domain, iErr)!{{{

   !--------------------------------------------------------------------

     type (domain_type), intent(inout) :: domain
     integer, intent(out) :: iErr

     type (block_type), pointer :: block_ptr

     type (mpas_pool_type), pointer :: meshPool, forcingPool, statePool, diagnosticsPool, &
                                       verticalMeshPool, scratchPool

     type (mpas_pool_type), pointer :: tracersPool

     integer, dimension(:), pointer :: maxLevelCell
     real (kind=RKIND), dimension(:), pointer :: ssh

     real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers
     real (kind=RKIND), dimension(:,:), pointer :: layerThickness, zMid

     real (kind=RKIND), dimension(:,:), pointer :: density
     real (kind=RKIND), dimension(:), pointer :: seaSurfacePressure, &
                                                 bottomPressure


     real(kind=RKIND), dimension(:,:), pointer :: origZMid
     integer, dimension(:), pointer :: origMaxLevelCell, modifySSHMask
     type (field2DReal), pointer :: origZMidField
     type (field1DInteger), pointer :: origMaxLevelCellField
     integer, pointer :: nCells, nVertLevels

     character (len=StrKIND), pointer :: config_iterative_init_variable

     integer :: iCell

     iErr = 0

     call mpas_pool_get_config(ocnConfigs, 'config_iterative_init_variable', config_iterative_init_variable)
     call mpas_pool_get_subpool(domain % blocklist % structs, 'scratch', scratchPool)

     call mpas_pool_get_field(scratchPool, 'scratchZMid', origZMidField)
     call mpas_allocate_scratch_field(origZMidField, .false.)
     call mpas_pool_get_field(scratchPool, 'scratchMaxLevelCell', origMaxLevelCellField)
     call mpas_allocate_scratch_field(origMaxLevelCellField, .false.)

     block_ptr => domain % blocklist
     do while(associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
       call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
       call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

       call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
       call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

       call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

       call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

       call mpas_pool_get_array(forcingPool, 'seaSurfacePressure', seaSurfacePressure)

       call mpas_pool_get_array(statePool, 'ssh', ssh, 1)
       call mpas_pool_get_array(diagnosticsPool, 'bottomPressure', bottomPressure)
       call mpas_pool_get_array(diagnosticsPool, 'density', density)
       call mpas_pool_get_array(diagnosticsPool, 'modifySSHMask', modifySSHMask)

       call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
       call mpas_pool_get_array(scratchPool, 'scratchZMid', origZMid)
       call mpas_pool_get_array(scratchPool, 'scratchMaxLevelCell', origMaxLevelCell)

       ! compute the bottom pressure assuming ssh = 0 (the full weight of the water column)
       call compute_density_and_bottom_pressure(meshPool, statePool, diagnosticsPool, scratchPool, &
                                              layerThickness, density, bottomPressure, iErr)
       if(iErr .ne. 0) then
         write(stderrUnit,*) 'ERROR: compute_density_and_bottom_pressure failed.'
         return
       end if

       if(config_iterative_init_variable == 'ssh') then
         do iCell = 1, nCells
           if(modifySSHMask(iCell) == 1) then
             ! compute ssh where pressure equals seaSurfacePressure
             ssh(iCell) = find_z_given_pressure(seaSurfacePressure(iCell), density(:,iCell), &
                                                layerThickness(:,iCell), nVertLevels, maxLevelCell(iCell))
           else
             ssh(iCell) = 0.0_RKIND
           end if
         end do
       end if

       ! save the old zMid for use in tracer inerpolation
       origZMid(:,:) = zMid(:,:)
       origMaxLevelCell(:) = maxLevelCell(:)

       block_ptr => block_ptr % next
     end do !block_ptr

     ! update the vertical grid based on the new ssh
     call ocn_init_vertical_grid(domain, iErr, updateOnly=.true.)

     if(iErr .ne. 0) then
       write(stderrUnit,*) 'ERROR: ocn_init_vertical_grid failed.'
       return
     end if

     block_ptr => domain % blocklist
     do while(associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
       call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
       call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

       call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

       call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
       call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

       call mpas_pool_get_array(forcingPool, 'seaSurfacePressure', seaSurfacePressure)

       call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
       call mpas_pool_get_array(diagnosticsPool, 'bottomPressure', bottomPressure)
       call mpas_pool_get_array(diagnosticsPool, 'density', density)
       call mpas_pool_get_array(diagnosticsPool, 'modifySSHMask', modifySSHMask)

       call mpas_pool_get_array(scratchPool, 'scratchZMid', origZMid)
       call mpas_pool_get_array(scratchPool, 'scratchMaxLevelCell', origMaxLevelCell)

       ! interpolate active tracers to the new zMid
       call interpolate_activeTracers(meshPool, origZMid, zMid, &
                                      origMaxLevelCell, maxLevelCell, &
                                      activeTracers, iErr)
       if(iErr .ne. 0) then
         write(stderrUnit,*) 'ERROR: interpolate_activeTracers failed.'
         return
       end if

       if(config_iterative_init_variable == 'ssp') then
         ! now compute what the weight of each ocean column with the new layer thickness,
         ! stored temporarily in seaSurfacePressure to not require an extra scratch variable
         call compute_density_and_bottom_pressure(meshPool, statePool, diagnosticsPool, scratchPool, &
                                                  layerThickness, density, seaSurfacePressure, iErr)

         if(iErr .ne. 0) then
            write(stderrUnit,*) 'ERROR: compute_density_and_bottom_pressure failed.'
            return
         end if

         ! the SSP is the weight of the full ocean column minus the weight of the column
         ! with the SSH applied (temporarily stored in seaSurfacePressure)
         where(modifySSHMask(:) == 1)
           seaSurfacePressure(:) = bottomPressure(:) - seaSurfacePressure(:)
         elsewhere
           seaSurfacePressure(:) = 0.0_RKIND
         end where
       end if

       block_ptr => block_ptr % next
     end do !block_ptr

     call mpas_deallocate_scratch_field(origZMidField, .false.)
     call mpas_deallocate_scratch_field(origMaxLevelCellField, .false.)

   !--------------------------------------------------------------------

   end subroutine initial_guess_ssp_ssh!}}}

!***********************************************************************
!
!  routine ocn_update_ssp_from_ssh
!
!> \brief   Update the SSP based on a change in SSH
!> \author  Xylar Asay-Davis
!> \date    10/20/2015
!> \details
!>  This routine adds a perturbation to the sea-surface pressure (SSP)
!>  based on the change in the sea-surface height (SSH) over a short
!>  forward run as part of an iterative process for finding a compatible
!>  SSP given a desired SSH.

!-----------------------------------------------------------------------

   subroutine ocn_update_ssp_from_ssh(domain, iErr)!{{{

   !--------------------------------------------------------------------

     type (domain_type), intent(inout) :: domain
     integer, intent(out) :: iErr

     type (block_type), pointer :: block_ptr

     type (mpas_pool_type), pointer :: meshPool, statePool, diagnosticsPool, verticalMeshPool, scratchPool, forcingPool

     real (kind=RKIND), dimension(:), pointer :: deltaSSH

     real (kind=RKIND), dimension(:,:), pointer :: layerThickness

     real (kind=RKIND), dimension(:,:), pointer :: density
     real (kind=RKIND), dimension(:), pointer :: seaSurfacePressure, &
                                                 bottomPressure
     integer, dimension(:), pointer :: modifySSHMask

     integer, pointer :: nCells
     logical, pointer :: config_modify_open_ocean_ssh
     integer :: iCell

     iErr = 0

     call mpas_pool_get_config(ocnConfigs, 'config_modify_open_ocean_ssh', config_modify_open_ocean_ssh)

     block_ptr => domain % blocklist
     do while(associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
       call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)

       call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

       call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

       call mpas_pool_get_array(forcingPool, 'seaSurfacePressure', seaSurfacePressure)

       call mpas_pool_get_array(diagnosticsPool, 'bottomPressure', bottomPressure)
       call mpas_pool_get_array(diagnosticsPool, 'density', density)
       call mpas_pool_get_array(diagnosticsPool, 'deltaSSH', deltaSSH)
       call mpas_pool_get_array(diagnosticsPool, 'modifySSHMask', modifySSHMask)

       ! compute the density and weight of each water column (stored in bottomPressure)
       call compute_density_and_bottom_pressure(meshPool, statePool, diagnosticsPool, scratchPool, &
                                                layerThickness, density, bottomPressure, iErr)

       if(iErr .ne. 0) then
          write(stderrUnit,*) 'ERROR: compute_density_and_bottom_pressure failed.'
          return
       end if

       do iCell = 1,nCells
         if(modifySSHMask(iCell) == 1) then
           ! Moving the SSH up or down by deltaSSH would change the SSP by density(SSH)*g*deltaSSH.
           ! If deltaSSH is positive (moving up), it means the SSP is too small and if deltaSSH
           ! is negative (moving down), it means SSP is too large, the sign of the second term
           ! makes sense.
           seaSurfacePressure(iCell) = max(seaSurfacePressure(iCell) + density(1,iCell)*gravity*deltaSSH(iCell), &
                                           0.0_RKIND)
         else
           ! the SSP should remain zero in the open ocean, though it may be that we are allowing the SSH to evolve here
           seaSurfacePressure(iCell) = 0.0_RKIND
         end if

         ! add the SSP, since bottomPressure only contains the weight of the water column up to now
         bottomPressure(iCell) = bottomPressure(iCell) + seaSurfacePressure(iCell)
       end do

       block_ptr => block_ptr % next
     end do !block_ptr

   !--------------------------------------------------------------------

   end subroutine ocn_update_ssp_from_ssh!}}}

!***********************************************************************
!
!  routine compute_density_and_bottom_pressure
!
!> \brief   Compute bottom pressure from current state
!> \author  Xylar Asay-Davis
!> \date    06/05/2015
!> \details
!>  This routine computes bottom pressure from layerThickness and density
!>  based on T and S from the current state. The bottom pressure with and
!>  without the ssh can be used to determine a seaSurfacePressure that
!>  roughly consistent with the ssh.
!>  bottomDepth, restingThickness, layerThicknesses, zMid, maxLevelCell,
!>  temperature and salinity should have already been initialized before
!>  calling this routine.

!----------------------------------------------------------------------

   subroutine compute_density_and_bottom_pressure(meshPool, statePool, diagnosticsPool, &
                                                  scratchPool, layerThickness, density, &
                                                  bottomPressure, iErr)!{{{

   !--------------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: meshPool, statePool
      type (mpas_pool_type), intent(inout) :: diagnosticsPool, scratchPool
      real (kind=RKIND), dimension(:,:), intent(in) :: layerThickness


      real (kind=RKIND), dimension(:,:), intent(out) :: density
      real (kind=RKIND), dimension(:), intent(out) :: bottomPressure
      integer, intent(out) :: iErr

      ! Define dimension pointers
      integer, pointer :: nCells, nVertLevels

      integer, dimension(:), pointer :: maxLevelCell

      ! Define variable pointers
      integer :: iCell, k


      iErr = 0

      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

      call ocn_equation_of_state_density(statePool, diagnosticsPool, meshPool, scratchPool, 0, 'relative', density, iErr, &
           timeLevelIn=1)

      if(iErr .ne. 0) then
        write(stderrUnit,*) 'ERROR: ocn_equation_of_state_density failed.'
        return
      end if

      do iCell = 1, nCells
        if(maxLevelCell(iCell) <= 0) cycle
        bottomPressure(iCell) = 0.0_RKIND
        do k=1,maxLevelCell(iCell)
          bottomPressure(iCell) = bottomPressure(iCell) &
             + density(k,iCell)*gravity*layerThickness(k,iCell)
        end do
      end do

   !--------------------------------------------------------------------

   end subroutine compute_density_and_bottom_pressure!}}}

!***********************************************************************
!
!  routine interpolate_activeTracers
!
!> \brief   interpolate the active tracers from reference fields
!> \author  Xylar Asay-Davis
!> \date    10/12/2015
!> \details
!>  Perform linear interpolation of T and S from reference fields without
!>  the sea-surface height (SSH) displacement at refZMid to new locations
!>  zMid that take the SSH into account.

!-----------------------------------------------------------------------

   subroutine interpolate_activeTracers(meshPool, inZMid, outZMid, &
                                        inMaxLevelCell, outMaxLevelCell, &
                                        activeTracers, iErr)!{{{

   !--------------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: meshPool
      real (kind=RKIND), dimension(:,:), intent(in) :: inZMid, outZMid
      integer, dimension(:), intent(in) :: inMaxLevelCell, outMaxLevelCell

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: activeTracers
      integer, intent(out) :: iErr

      ! Define dimension pointers
      integer, pointer :: nCells, nVertLevels

      ! Define variable pointers
      integer :: iCell, inKMax, outKMax

      real (kind=RKIND), dimension(:), allocatable :: inTracerColumn, outTracerColumn

      integer :: nTracers, iTracer

      iErr = 0

      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      nTracers = size(activeTracers, dim=1)
      allocate(inTracerColumn(nVertLevels),outTracerColumn(nVertLevels))

      do iCell = 1, nCells
        inKMax = inMaxLevelCell(iCell)
        outKMax = outMaxLevelCell(iCell)
        if((inKMax <= 0) .or. (outKMax <= 0))  cycle

        do iTracer = 1, nTracers
          inTracerColumn(:) = activeTracers(iTracer,:,iCell)
          call ocn_init_interpolation_linear_vert(inZMid(1:inKMax,iCell), &
                                                  inTracerColumn(1:inKMax), &
                                                  inKMax, &
                                                  outZMid(1:outKMax,iCell), &
                                                  outTracerColumn(1:outKMax), &
                                                  outKMax, &
                                                  extrapolate=.true.)
          activeTracers(iTracer,:,iCell) = outTracerColumn(:)
        end do
      end do

      deallocate(inTracerColumn, outTracerColumn)

   !--------------------------------------------------------------------

   end subroutine interpolate_activeTracers!}}}

!***********************************************************************
!
!  funciton find_z_given_pressure
!
!> \brief   interpolate the active tracers from reference fields
!> \author  Xylar Asay-Davis
!> \date    10/13/2015
!> \details
!>  In a column, find the depth at which the hydrostatic pressure reaches a given
!>  value provided a density profile.

!-----------------------------------------------------------------------

   function find_z_given_pressure(pressure, density, layerThickness, nVertLevels, maxLevelCell) result(z)
      real (kind=RKIND), intent(in) :: pressure
      real (kind=RKIND), intent(in), dimension(nVertLevels) :: density, layerThickness
      integer, intent(in) :: nVertLevels, maxLevelCell
      real (kind=RKIND) :: z

      integer :: k
      real (kind=RKIND) :: pressureTop, pressureBot

      pressureTop = 0.0_RKIND
      z = 0.0_RKIND

      if(maxLevelCell <= 0) return

      do k = 1, maxLevelCell
        pressureBot = pressureTop + density(k)*gravity*layerThickness(k)
        if(pressure < pressureBot) then
           ! note: this will simply extrapolate if presssure is negative for some reason
           z = z - (pressure - pressureTop)/(pressureBot - pressureTop)*layerThickness(k)
           return
        end if
        z = z - layerThickness(k)
        pressureTop = pressureBot
     end do

   end function find_z_given_pressure

!***********************************************************************

end module ocn_init_ssh_and_ssp

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
