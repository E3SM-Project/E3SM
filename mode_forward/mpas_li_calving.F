! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_calving
!
!> \brief MPAS land ice calving scheme
!> \author William Lipscomb
!> \date   September 2015
!> \details
!>  This module contains several options for calving ice.
!
!-----------------------------------------------------------------------

module li_calving

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar
   use li_setup
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

   public :: do_calving

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------


!***********************************************************************
   contains
!***********************************************************************


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  !  routine do_calving
!
!> \brief MPAS land ice calving scheme
!> \author William Lipscomb
!> \date   September 2015
!> \details
!>  This routine contains several options for calving ice:
!> (0) Do nothing
!> (1) Calve all floating ice
!> (2) Calve ice based on a topographic threshold
!> (3) Calve ice based on an ice thickness threshold
!-----------------------------------------------------------------------

   subroutine do_calving(domain, deltat, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real (kind=RKIND), intent(in) :: deltat  !< Input: time step

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: &
         domain          !< Input/Output: domain object 

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      type (dm_info), pointer :: dminfo
      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: meshPool 
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: scratchPool 
      type (mpas_pool_type), pointer :: velocityPool

      ! calving-relevant config options
      character (len=StrKIND), pointer :: config_calving
      logical, pointer :: config_print_calving_info
      real(kind=RKIND), pointer :: config_calving_topographic_threshold, &
                                   config_calving_thickness_threshold, &
                                   config_calving_timescale, &
                                   config_sea_level

      ! fields for calving masks
      ! The calvingLaw mask is used for the floating and topographic_threshold options
      ! The inactiveMargin and ocean masks are used for the more complex thickness_threshold option

      type (field1dInteger), pointer :: calvingLawMaskField
      integer, dimension(:), pointer :: calvingLawMask     ! = 1 where calving-law criterion is satisfied, else = 0

      type (field1dInteger), pointer :: inactiveMarginMaskField
      integer, dimension(:), pointer :: inactiveMarginMask  ! = 1 for inactive cells (thin or no ice) that have 1 or more active neighbors

      type (field1dInteger), pointer :: oceanMaskField
      integer, dimension(:), pointer :: oceanMask   ! = 1 for cells that are not land and do not have active ice
                                                    ! may include floating cells with inactive ice

      integer, pointer :: nCells

      integer, dimension(:), pointer :: &
           nCellsOnCell,     & ! number of cells that border each cell
           cellMask            ! bit mask describing whether ice is floating, dynamically active, etc.
      
      integer, dimension(:), pointer :: &
           indexToCellID       ! list of global cell IDs

      integer, dimension(:,:), pointer :: &
           cellsOnCell         ! list of cells that neighbor each cell

      real (kind=RKIND) ::  &
           calvingFraction    ! fraction of ice that calves in each column; depends on calving_timescale

      real (kind=RKIND), dimension(:), pointer :: &
           thickness,        & ! ice thickness
           bedTopography       ! bed topography (negative below sea level)

      real (kind=RKIND), dimension(:), pointer :: &
           calvingThickness    ! thickness of ice that calves (computed in this subroutine)
                               ! typically the entire ice thickness, but will be a fraction of the thickness if calving_timescale > dt

      integer :: iCell, iCellOnCell, iCellNeighbor

      integer :: err_tmp

      real (kind=RKIND), parameter :: scyr = 31536000.0_RKIND  ! seconds per 365-day year; diagnostic only
                                                               !TODO - put this in a constants module?
      err = 0

      ! get config options
      dminfo => domain % dminfo
      call mpas_pool_get_config(liConfigs, 'config_calving', config_calving)
      call mpas_pool_get_config(liConfigs, 'config_calving_thickness_threshold', config_calving_thickness_threshold)
      call mpas_pool_get_config(liConfigs, 'config_calving_topographic_threshold', config_calving_topographic_threshold)
      call mpas_pool_get_config(liConfigs, 'config_calving_timescale', config_calving_timescale)
      call mpas_pool_get_config(liConfigs, 'config_print_calving_info', config_print_calving_info)
      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)

      ! based on the calving timescale, set the fraction of ice that calves
      if (config_calving_timescale > 0.0_RKIND) then
         calvingFraction = max(deltat/config_calving_timescale, 1.0_RKIND)
      else
         calvingFraction = 1.0_RKIND   ! calve the entire thickness in eligible columns
      endif

      ! block loop
      block => domain % blocklist
      do while (associated(block))

         ! get pools
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)  ! required for cellMask computation
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

         ! get dimensions
         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

         ! get required fields from the mesh pool
         call mpas_pool_get_array(meshPool, 'nCellsOnCell', nCellsOnCell)
         call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
         call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)  ! diagnostic only

         ! get required fields from the geometry pool
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
         call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)

         if (config_print_calving_info) then
            write(stderrUnit,*) 'Do ice calving, option =', trim(config_calving)
            write(stderrUnit,*) 'Calving timscale (yr) =', config_calving_timescale / scyr
         endif

         ! calculate masks - so we know where the ice is floating and/or dynamically active

         call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
         err = ior(err, err_tmp)

         ! initialize
         calvingThickness = 0.0_RKIND

         ! compute calving based on the calving_config option

         if (trim(config_calving) == 'none') then

            ! do nothing

         elseif (trim(config_calving) == 'thickness_threshold') then

            ! calve ice thinner than the threshold thickness

            !Note: This is not as simple as identifying floating ice thinner than the thickness threshold.
            !      The problem with that approach is that any ice advected in front of the calving front 
            !       would be instantly removed, making it impossible for the calving front to advance.
            !      Instead, we define an inactive margin containing cells that are inactive but border active cells.
            !      Cells on the inactive margin are protected from calving, but thin floating ice 
            !       beyond the inactive margin can calve.
            !
            !      Specifically, the rules are as follows:
            !      - Mark cells as ocean if not land and not active ice.
            !      - Mark cells as inactive margin if not active ice, but with an active ice neighbor.
            !      - Calve ice in ocean cells that are not on the inactive margin.

            if (config_print_calving_info) then
               write(stderrUnit,*) 'Calving thickness threshold (m) =', config_calving_thickness_threshold
            endif

            ! get scratch fields for calving
            ! 'true' flag means to allocate the field for a single block

            call mpas_pool_get_field(scratchPool, 'iceCellMask',  inactiveMarginMaskField)
            call mpas_allocate_scratch_field(inactiveMarginMaskField, .true.)
            inactiveMarginMask => inactiveMarginMaskField % array

            call mpas_pool_get_field(scratchPool, 'iceCellMask2', oceanMaskField)
            call mpas_allocate_scratch_field(oceanMaskField, .true.)
            oceanMask => oceanMaskField % array

            ! Identify cells that are inactive but border dynamically active cells
            !WHL - This might not work as intended if there are cells on the margin that are thick but have Dirichlet BC
            !      and thus are classified as inactive.
            do iCell = 1, nCells
               inactiveMarginMask(iCell) = 0
               if (.not. li_mask_is_dynamic_ice(cellMask(iCell))) then   ! either thin (inactive) ice or no ice
                  do iCellOnCell = 1, nCellsOnCell(iCell)
                     iCellNeighbor = cellsOnCell(iCellOnCell,iCell)
                     if (li_mask_is_dynamic_margin(cellMask(iCellNeighbor))) then  ! neighbor cell is on the dynamic ice margin
                        inactiveMarginMask(iCell) = 1
                        exit
                     endif
                  enddo  ! iCellOnCell
               endif   ! not dynamic ice
            enddo   ! iCell

            ! Identify ocean cells (not land and not dynamic ice; may include inactive floating ice)

            where (bedTopography < config_sea_level .and. .not.li_mask_is_dynamic_ice(cellMask))
               oceanMask = 1
            elsewhere
               oceanMask = 0
            endwhere

            ! Calve ice in ocean cells that are not on the inactive margin

            where (oceanMask == 1 .and. inactiveMarginMask == 0)
               calvingThickness = thickness * calvingFraction
               thickness = thickness - calvingThickness
            endwhere

         else   ! other calving options (floating, topographic_threshold)

            ! get scratch fields for calving

            call mpas_pool_get_field(scratchPool, 'iceCellMask',  calvingLawMaskField)
            call mpas_allocate_scratch_field(calvingLawMaskField, .true.)
            calvingLawMask => calvingLawMaskField % array
            
            if (trim(config_calving) == 'floating') then

               ! calve floating ice
               ! Note: The floating_ice mask includes all floating ice, both inactive and active

               where (li_mask_is_floating_ice(cellMask))
                  calvingLawMask = 1
               elsewhere
                  calvingLawMask = 0
               endwhere

            elseif (trim(config_calving) == 'topographic_threshold') then
            
               ! calve ice where the bed topography lies below a threshold depth

               if (config_print_calving_info) then
                  write(stderrUnit,*) 'Calving topographic threshold (m) =', config_calving_topographic_threshold
               endif
            
               where (bedTopography < config_calving_topographic_threshold + config_sea_level)
                  calvingLawMask = 1
               elsewhere
                  calvingLawMask = 0
               endwhere

               !Note: Could add other calving options here (e.g., damage-based calving)

            endif   ! floating or topographic_threshold

            ! Calve ice where specified by the calving law.
            ! If the calving timescale > 0, then will calve only a fraction of the thickness in each calving cell.

            where (calvingLawMask == 1)

               calvingThickness = thickness * calvingFraction
               thickness = thickness - calvingThickness

            endwhere
               
         endif   ! config_calving

         ! Optionally, print a list of cells with calving

         if (config_print_calving_info) then
            write(stderrUnit,*) ' '
            write(stderrUnit,*) 'Global cell ID, calving thickness:'
            do iCell = 1, nCells
               if (calvingThickness(iCell) > 0.0_RKIND) then
                  write(stderrUnit,*) indexToCellID(iCell), calvingThickness(iCell)
               endif
            enddo
         endif

         block => block % next
      enddo

      ! === error check
      if (err > 0) then
          write (stderrUnit,*) "An error has occurred in do_calving."
      endif

   !--------------------------------------------------------------------
    end subroutine do_calving

!***********************************************************************
!***********************************************************************
! Private subroutines:
!***********************************************************************
!***********************************************************************

! No private subroutines so far
   
end module li_calving


