! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
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

   public :: ocn_init_ssh_and_landIcePressure_vertical_grid, &
             ocn_init_ssh_and_landIcePressure_balance


   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains


!***********************************************************************
!
!  routine ocn_init_ssh_and_landIcePressure_vertical_grid
!
!> \brief   Initialize z* vertical grid based on SSH
!> \author  Xylar Asay-Davis
!> \date    10/21/2015
!> \details
!>  This routine sets up the vertical grid (layerThickness,
!>  zMid and restingThickness) needed for computing SSH from
!>  land-ice pressure or visa versa. bottomDepth, refBottomDepth, maxLevelCell
!>  and modifySSHMask must have been computed by the test case
!>  before calling this routine.  If
!>  config_iterative_init_variable = 'landIcePressure' or 'landIcePressure_from_top_density', the test
!>  case must compute ssh before calling this routine.
!>  modifySSHMask should be set to 1 wherever the ssh or landIcePressure
!>  should be modified for consistency (e.g. under land ice). This
!>  routine will take care of setting up partial bottom cells
!>  by calling ocn_alter_bottomDepth_for_pbcs (except for the
!>  Haney-number-constrained coordinate, which handle thin bottom
!>  cells via the Haney-number constraint.

!-----------------------------------------------------------------------

   subroutine ocn_init_ssh_and_landIcePressure_vertical_grid(domain, iErr)!{{{

   !--------------------------------------------------------------------

     type (domain_type), intent(inout) :: domain
     integer, intent(out) :: iErr

   !--------------------------------------------------------------------

     iErr = 0
     call ocn_init_vertical_grid(domain, updateWithSSH=.false., iErr=iErr)

   end subroutine ocn_init_ssh_and_landIcePressure_vertical_grid

!***********************************************************************
!
!  routine ocn_init_ssh_and_landIcePressure_balance
!
!> \brief   Compute the balance land-ice pressure given the SSH or visa versa
!> \author  Xylar Asay-Davis
!> \date    8/8/2016
!> \details
!>  This routine either updates SSH based on land-ice pressure (if config_iterative_init_variable = 'ssh')
!>  or visa versa (if config_iterative_init_variable = 'landIcePressure' or 'landIcePressure_from_top_density').
!>  The routine produces an initial guess at land-ice pressure or SSH using either the density of the topmost layer
!>  ('landIcePressure_from_top_density') or of all layers above the SSH ('ssh' or 'landIcePressure')
!>  to determine the effective density of seawater within the land ice.
!>  The resulting land-ice pressure and SSH are approximately consistent with one another
!>  in the sense that the horizontal pressure-gradient force (HPGF)
!>  should be small at the ocean surface.
!>  ocn_init_ssh_and_landIcePressure_vertical_grid should be called to produce
!>  the appropriate vertical grid before calling this subroutine.
!>  activeTracers should be initialized based on this vertical grid.
!>  Upon completion, the vertical grid will have been updated
!>  to be consistent with the SSH and the activeTracers will have been
!>  interpolated to the new grid.

!-----------------------------------------------------------------------

   subroutine ocn_init_ssh_and_landIcePressure_balance(domain, iErr)!{{{

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
     real (kind=RKIND), dimension(:), pointer :: landIcePressure, landIceDraft, &
                                                 effectiveDensityInLandIce


     real(kind=RKIND), dimension(:,:), pointer :: origZMid
     integer, dimension(:), pointer :: origMaxLevelCell, modifySSHMask
     type (field2DReal), pointer :: origZMidField
     type (field1DInteger), pointer :: origMaxLevelCellField
     integer, pointer :: nCells, nVertLevels

     character (len=StrKIND), pointer :: config_iterative_init_variable

     integer :: iCell

     iErr = 0

     call mpas_pool_get_config(ocnConfigs, 'config_iterative_init_variable', config_iterative_init_variable)

     ! compute density (needed regardless of config_iterative_init_variable)

     block_ptr => domain % blocklist
     do while(associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
       call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)

       call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

       call mpas_pool_get_array(diagnosticsPool, 'density', density)

       call ocn_equation_of_state_density(statePool, diagnosticsPool, meshPool, &
           scratchPool, nCells, 0, 'relative', density, iErr, &
           timeLevelIn=1)

       if(iErr .ne. 0) then
         write(stderrUnit,*) 'ERROR: ocn_equation_of_state_density failed.'
         return
       end if

       block_ptr => block_ptr % next
     end do !block_ptr

     ! first, handle the simple case where we're going to compute landIcePressure from the density at the top.
     ! In this case, we already computed the correct vertical grid with the ssh, so all we have to do is
     ! landIcePressure from the denisty we just got.
     if(config_iterative_init_variable ==  'landIcePressure_from_top_density') then
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

         call mpas_pool_get_array(forcingPool, 'landIcePressure', landIcePressure)

         call mpas_pool_get_array(forcingPool, 'landIceDraft', landIceDraft)

         call mpas_pool_get_array(statePool, 'ssh', ssh, 1)
         call mpas_pool_get_array(statePool, 'effectiveDensityInLandIce', effectiveDensityInLandIce, 1)
         call mpas_pool_get_array(diagnosticsPool, 'density', density)
         call mpas_pool_get_array(diagnosticsPool, 'modifySSHMask', modifySSHMask)

         call ocn_equation_of_state_density(statePool, diagnosticsPool, meshPool, &
             scratchPool, nCells, 0, 'relative', density, iErr, &
             timeLevelIn=1)

         if(iErr .ne. 0) then
           write(stderrUnit,*) 'ERROR: ocn_equation_of_state_density failed.'
           return
         end if

         do iCell = 1, nCells
           if(modifySSHMask(iCell) == 0) then
             ssh(iCell) = 0.0_RKIND
             landIcePressure(iCell) = 0.0_RKIND

             if (associated(effectiveDensityInLandIce)) &
               ! effective density cannot be determined
               effectiveDensityInLandIce(iCell) = 0.0_RKIND
             cycle
           end if

           landIcePressure(iCell) = max(0.0_RKIND, -density(1,iCell)*gravity*ssh(iCell))
           if (associated(effectiveDensityInLandIce)) &
             effectiveDensityInLandIce(iCell) = density(1,iCell)
         end do

         ! copy the SSH into the landIceDraft so we can use it later to remove it when
         ! computing sea-surface tilt
         landIceDraft(:) = ssh(:)

         block_ptr => block_ptr % next
       end do !block_ptr

       return
     end if

     ! The other cases are more complicated and require interpolating the activeTracers
     ! once the ssh or landIcePressure is determined

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

       call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
       call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

       call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

       call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

       call mpas_pool_get_array(forcingPool, 'landIcePressure', landIcePressure)

       call mpas_pool_get_array(forcingPool, 'landIceDraft', landIceDraft)

       call mpas_pool_get_array(statePool, 'ssh', ssh, 1)
       call mpas_pool_get_array(statePool, 'effectiveDensityInLandIce', effectiveDensityInLandIce, 1)
       call mpas_pool_get_array(diagnosticsPool, 'density', density)
       call mpas_pool_get_array(diagnosticsPool, 'modifySSHMask', modifySSHMask)

       call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
       call mpas_pool_get_array(scratchPool, 'scratchZMid', origZMid)
       call mpas_pool_get_array(scratchPool, 'scratchMaxLevelCell', origMaxLevelCell)


       call ocn_equation_of_state_density(statePool, diagnosticsPool, meshPool, &
           scratchPool, nCells, 0, 'relative', density, iErr, &
           timeLevelIn=1)

       if(iErr .ne. 0) then
         write(stderrUnit,*) 'ERROR: ocn_equation_of_state_density failed.'
         return
       end if

       do iCell = 1, nCells
         if(modifySSHMask(iCell) == 0) then
           ssh(iCell) = 0.0_RKIND
           landIcePressure(iCell) = 0.0_RKIND

           if (associated(effectiveDensityInLandIce)) &
             ! effective density cannot be determined
             effectiveDensityInLandIce(iCell) = 0.0_RKIND

           cycle
         end if

         if(config_iterative_init_variable == 'ssh') then
           ! compute ssh where pressure equals landIcePressure
           ssh(iCell) = find_z_given_pressure(landIcePressure(iCell), density(:,iCell), &
                                              layerThickness(:,iCell), nVertLevels, maxLevelCell(iCell))
         else
           ! compute landIcePressure based on hydrostatic pressure at SSH
           landIcePressure(iCell) = max(0.0_RKIND, find_pressure_given_z(ssh(iCell), density(:,iCell), &
                                              layerThickness(:,iCell), nVertLevels, maxLevelCell(iCell)))
         end if

         if (associated(effectiveDensityInLandIce)) then
           ! the effective density of ocean water in land ice is determined from the land-ice pressure and SSH
           effectiveDensityInLandIce(iCell) = -landIcePressure(iCell)/(gravity*ssh(iCell))
         end if

       end do

       ! save the old zMid for use in tracer inerpolation
       origZMid(:,:) = zMid(:,:)
       origMaxLevelCell(:) = maxLevelCell(:)

       ! copy the SSH into the landIceDraft so we can use it later to remove it when
       ! computing sea-surface tilt
       landIceDraft(:) = ssh(:)

       block_ptr => block_ptr % next
     end do !block_ptr

     ! update the vertical grid based on the new ssh
     call ocn_init_vertical_grid(domain, updateWithSSH=.true., iErr=iErr)

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

       call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)

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

       block_ptr => block_ptr % next
     end do !block_ptr

     call mpas_deallocate_scratch_field(origZMidField, .false.)
     call mpas_deallocate_scratch_field(origMaxLevelCellField, .false.)

   !--------------------------------------------------------------------

   end subroutine ocn_init_ssh_and_landIcePressure_balance

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
!> \date    8/8/2016
!> \details
!>  This routine sets up the vertical grid (layerThickness,
!>  zMid and restingThickness) needed for computing SSH from
!>  land-ice pressure or visa versa. bottomDepth, refBottomDepth and maxLevelCell
!>  must have been computed by the test case before calling this
!>  routine.  If config_iterative_init_variable = 'landIcePressure' or 'landIcePressure_from_top_density', the test
!>  case must compute ssh before calling this routine.  This
!>  routine will take care of setting up partial bottom cells
!>  by calling ocn_alter_bottomDepth_for_pbcs (except for the
!>  Haney-number-constrained coordinate, which handle thin bottom
!>  cells via the Haney-number constraint.
!-----------------------------------------------------------------------

   subroutine ocn_init_vertical_grid(domain, updateWithSSH, iErr)!{{{

   !--------------------------------------------------------------------

     type (domain_type), intent(inout) :: domain
     logical, intent(in) :: updateWithSSH
     integer, intent(out) :: iErr

     type (block_type), pointer :: block_ptr

     type (mpas_pool_type), pointer :: meshPool, statePool, diagnosticsPool, verticalMeshPool

     logical, pointer :: config_use_rx1_constraint

     character (len=StrKIND), pointer :: config_iterative_init_variable

     ! Define dimension pointers
     integer, pointer :: nCells, nVertLevels

     ! Define variable pointers
     integer, dimension(:), pointer :: maxLevelCell
     real (kind=RKIND), dimension(:), pointer :: refBottomDepth, bottomDepth, ssh
     real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness, zMid

     integer :: iCell

     logical :: initWithSSH, initRx1WithSSH

   !--------------------------------------------------------------------

     iErr = 0

     call mpas_pool_get_config(ocnConfigs, 'config_iterative_init_variable', config_iterative_init_variable)

     if(config_iterative_init_variable .ne. 'ssh' &
        .and. config_iterative_init_variable .ne. 'landIcePressure_from_top_density' &
        .and. config_iterative_init_variable .ne. 'landIcePressure') then
       iErr = 1
       write(stderrUnit,*) 'ERROR: invalid value for config_iterative_init_variable', trim(config_iterative_init_variable)
       return
     end if
     call mpas_pool_get_config(ocnConfigs, 'config_use_rx1_constraint', config_use_rx1_constraint)

     ! one reason for using config_iterative_init_variable == 'landIcePressure_from_top_density' is that we can immediately compute
     ! the vertical grid displaced by the ssh
     initWithSSH = updateWithSSH .or. (config_iterative_init_variable == 'landIcePressure_from_top_density')
     initRx1WithSSH = initWithSSH .and. config_use_rx1_constraint

     if(.not. config_use_rx1_constraint .and. .not. updateWithSSH) then
       ! only alter pbcs if we're initializing for the first time (updateWithSSH == .false.) so we haven't done it already
       ! and we're not going to handle bottomDepth another way (with the Haney-number-constrained coordinate)

       block_ptr => domain % blocklist
       do while(associated(block_ptr))
         call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)

         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

         call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
         call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

         do iCell = 1, nCells
           call ocn_alter_bottomDepth_for_pbcs(bottomDepth(iCell), refBottomDepth, maxLevelCell(iCell), iErr)
           if(iErr .ne. 0) then
             write(stderrUnit,*) 'ERROR: ocn_alter_bottomDepth_for_pbcs failed.'
             return
           end if
         end do
         block_ptr => block_ptr % next
       end do !block_ptr
     end if

     if(initRx1WithSSH) then
       ! We already know the ssh and landIcePressure we want to use.
       ! Compute the layer thicknesses and zMid based on topography and ssh.
       ! Use rx1 constraint to recompute the vertical grid.
       call ocn_init_vertical_grid_with_max_rx1(domain, iErr)

       if(iErr .ne. 0) then
         write(stderrUnit,*) 'ERROR: ocn_init_vertical_grid_with_max_rx1 failed.'
         return
       end if

     else
       ! we're initializing to z-star
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

         do iCell = 1, nCells
           if(initWithSSH) then
             ! we already know the ssh and landIcePressure we want to use.
             ! compute the layer thicknesses and zMid based on topography and ssh
             call ocn_compute_layerThickness_zMid_from_bottomDepth(layerThickness(:,iCell),zMid(:,iCell), &
                  refBottomDepth,bottomDepth(iCell), &
                  maxLevelCell(iCell),nVertLevels,iErr, &
                  restingThickness=restingThickness(:,iCell), &
                  ssh=ssh(iCell))
           else
             ! We don't know the ssh or landIcePressure yet, and we need tracers on a reference grid to figure it out.
             ! compute restingThickness and reference layerThickness and zMid based on topography with ssh=0
             ! (omitting ssh argument)
             call ocn_compute_layerThickness_zMid_from_bottomDepth(layerThickness(:,iCell),zMid(:,iCell), &
                  refBottomDepth,bottomDepth(iCell), &
                  maxLevelCell(iCell),nVertLevels,iErr, &
                  restingThickness=restingThickness(:,iCell))
           end if

           if(iErr .ne. 0) then
             write(stderrUnit,*) 'ERROR: ocn_compute_layerThickness_zMid_from_bottomDepth failed.'
             return
           end if
         end do !iCell

         block_ptr => block_ptr % next
       end do !block_ptr
     end if

   end subroutine ocn_init_vertical_grid

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
          outTracerColumn(:) = -9.969209968386869e+36_RKIND
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
!  funciton find_pressure_given_z
!
!> \brief   Determine the pressure at a given depth
!> \author  Xylar Asay-Davis
!> \date    8/8/2016
!> \details
!>  In a column, find the hydrostatic pressure at a given depth with
!>  the given density profile.

!-----------------------------------------------------------------------

   function find_pressure_given_z(z, density, layerThickness, nVertLevels, maxLevelCell) result(pressure)
      real (kind=RKIND), intent(in) :: z
      real (kind=RKIND), intent(in), dimension(nVertLevels) :: density, layerThickness
      integer, intent(in) :: nVertLevels, maxLevelCell
      real (kind=RKIND) :: pressure

      integer :: k
      real (kind=RKIND) :: pressureTop, pressureBot, zTop, zBot

      pressure = 0.0_RKIND
      zTop = 0.0_RKIND

      if(maxLevelCell <= 0) return

      do k = 1, maxLevelCell
        zBot = zTop - layerThickness(k)
        if(z > zBot) then
           ! note: this will simply extrapolate if z is positive for some reason
           pressure = pressure + density(k)*gravity*(zTop - z)
           return
        end if
        pressure = pressure + density(k)*gravity*layerThickness(k)
        zTop = zBot
     end do

   end function find_pressure_given_z

!***********************************************************************
!
!  funciton find_z_given_pressure
!
!> \brief  Find the depth at which pressure has a given value
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

end module ocn_init_ssh_and_landIcePressure

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
