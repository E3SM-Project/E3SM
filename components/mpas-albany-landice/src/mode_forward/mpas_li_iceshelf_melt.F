! Copyright (c) 2013-2018,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_iceshelf_melt
!
!> \brief MPAS land ice vertical temperature/enthalpy solver
!> \author William Lipscomb
!> \date   October 2015
!> \details
!>  This module contains solvers for the vertical temperature
!>  and/or enthalpy profile.
!
!-----------------------------------------------------------------------

module li_iceshelf_melt

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timer
   use mpas_abort
   use mpas_log

   use li_setup
   use li_mask
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

   public :: li_basal_melt_floating_ice, li_face_melt_grounded_ice

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
!  !  routine li_face_melt_grounded_ice
!
!> \brief MPAS land ice solver for face melt of grounded ice
!> \author Trevor Hillebrand
!> \date   February 2021
!> \details
!>  This routine computes face melting for floating ice.
!>  The following options are supported:
!>  (1) Do nothing (config_front_mass_bal_grounded = 'none')
!>  (2) Prescribed uniform face melt speed (config_front_mass_bal_grounded = 'uniform')
!>  (3) Face melt speed as in ISMIP6 (config_front_mass_bal_grounded = 'ismip6')

!-----------------------------------------------------------------------

   subroutine li_face_melt_grounded_ice(domain, err)
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

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

      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: geometryPool, velocityPool, meshPool
      character(len=StrKIND), pointer :: &
           config_front_mass_bal_grounded, & ! option for submarine mass balance at glacier front
           config_basal_mass_bal_float      ! option for mass balance under ice shelf
      real (kind=RKIND), dimension(:), pointer :: &
              faceMeltSpeed, &
              faceMeltingThickness, &
              thickness
      integer, dimension(:), pointer :: cellMask
      integer :: err_tmp
      logical :: applyToFloating, applyToGrounded, applyToGroundingLine

      err = 0
      err_tmp = 0

      call mpas_pool_get_config(liConfigs, 'config_front_mass_bal_grounded', config_front_mass_bal_grounded)
      call mpas_pool_get_config(liConfigs, 'config_basal_mass_bal_float', config_basal_mass_bal_float)

      ! Determine where to apply melt. This avoids redundant melting at the grounding line
      ! if both face-melting and ice shelf melting are active.
      applyToFloating = .false. ! As of 06/2021 there are no face-melt routines for floating ice
      if (trim(config_basal_mass_bal_float) .eq. 'none') then
         applyToGrounded = .false.
         applyToGroundingLine = .true.
      else
         applyToGrounded = .true.
         applyToGroundingLine = .false.
      endif


      if ( trim(config_front_mass_bal_grounded) == 'ismip6' &
          .or. trim(config_front_mass_bal_grounded) == 'uniform' ) then
         call grounded_face_melt_ismip6(domain, applyToGrounded, &
                    applyToFloating, applyToGroundingLine, err_tmp)
         err = ior(err, err_tmp)
         
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
            call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_array(geometryPool, 'faceMeltingThickness', faceMeltingThickness)
            call mpas_pool_get_array(geometryPool, 'thickness', thickness)

            ! Update halos on calvingThickness or faceMeltingThickness before applying it.
            ! Testing seemed to indicate this is not necessary, but I don't understand
            ! why not, so leaving it.
            ! NOTE: THIS WILL NOT WORK ON MULTIPLE BLOCKS PER PROCESSOR
            call mpas_timer_start("halo updates")
            call mpas_dmpar_field_halo_exch(domain, 'faceMeltingThickness')
            call mpas_timer_stop("halo updates")

            ! Apply facemelt: open the Ark
            thickness(:) = thickness(:) - faceMeltingThickness(:)

            ! update mask
            call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
            err = ior(err, err_tmp)

            block => block % next
         enddo   ! associated(block) 
   
      elseif ( trim(config_front_mass_bal_grounded) == 'none' ) then
         ! Zero entire field

         ! block loop
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
            call mpas_pool_get_array(geometryPool, 'faceMeltSpeed', faceMeltSpeed)
            call mpas_pool_get_array(geometryPool, 'faceMeltingThickness', faceMeltingThickness)

            faceMeltSpeed(:) = 0.0_RKIND
            faceMeltingThickness(:) = 0.0_RKIND

            block => block % next
         enddo   ! associated(block)
      else
         call mpas_log_write('Unknown option selected for config_front_mass_bal_grounded:' // &
                             trim(config_front_mass_bal_grounded), MPAS_LOG_ERR)
         err = ior(err, 1)
      endif

   end subroutine li_face_melt_grounded_ice
!-----------------------------------------------------------------------
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||r
!
!  !  routine li_basal_melt_floating_ice
!
!> \brief MPAS land ice solver for basal melt of floating ice
!> \author William Lipscomb
!> \date   November 2015
!> \details
!>  This routine computes basal melting for floating ice.
!>  The following options are supported:
!>  (1) Do nothing (config_basal_mass_bal_float = 'none')
!>  (2) Read melt rate from a file (config_basal_mass_bal_float = 'file')
!>  (2) Prescribed constant basal melt rate (config_basal_mass_bal_float = 'constant')
!>  (3) Basal melt rate as in MISMIP+ (config_basal_mass_bal_float = 'mismip')

!-----------------------------------------------------------------------

   subroutine li_basal_melt_floating_ice(domain, err)
      use li_calving, only: li_flood_fill
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

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

      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: velocityPool   ! needed for mask subroutine
      type (mpas_pool_type), pointer :: scratchPool

      type (field1dInteger), pointer :: seedMaskField
      type (field1dInteger), pointer :: growMaskField
      integer, dimension(:), pointer :: seedMask, growMask !masks to pass to flood-fill routine

      integer, pointer :: &
           nCellsSolve, nCells         ! number of locally owned cells & total cells on this proc

      logical, pointer :: &
           config_print_thermal_info   ! if true, print debug info

      character(len=StrKIND), pointer :: &
           config_basal_mass_bal_float ! option for basal mass balance of floating ice

      character(len=StrKIND), pointer :: &
           config_front_mass_bal_grounded ! option for submarine mass balance
                                          ! at grounded glacier front

      real(kind=RKIND), pointer :: &
           config_thermal_thickness, & ! minimum thickness (m) for temperature calculations
           config_sea_level,         & ! sea level (m) relative to z = 0
           config_bmlt_float_flux,   & ! constant heat flux (W/m^2) applied to the base of floating ice; positive upward
           config_bmlt_float_xlimit    ! x value (m) defining region where bmlt_float_flux is applied;
                                       !   melt only where abs(x) > xlimit
      real (kind=RKIND), pointer :: config_ice_density !< ice density

      integer, dimension(:), pointer :: &
           cellMask, &                    ! bit mask describing whether ice is floating, dynamically active, etc.
           nEdgesOnCell, &
           edgeMask

      integer, dimension(:,:), pointer :: edgesOnCell, cellsOnCell

      type (field1dInteger), pointer :: thermalCellMaskField

      integer, dimension(:), pointer :: &
           thermalCellMask             ! mask for thermal calculations
                                       ! = 1 where thickness > config_thermal_thickness, elsewhere = 0

      real (kind=RKIND), dimension(:), pointer :: &
           xCell                       ! x coordinate for each cell (m)

      real (kind=RKIND), dimension(:), pointer :: &
           floatingBasalMassBal,     & ! basal mass balance for floating ice
           thickness,                & ! ice thickness (m)
           lowerSurface,             & ! lower surface elevation (m)
           bedTopography               ! bed topography (m; negative below sea level)

      real(kind=RKIND), pointer :: daysSinceStart

      integer :: iCell, jCell, iEdge, iNeighbor, err_tmp

      ! Local variables for some melt methods

      real(kind=RKIND) :: &
           bmlt_float_rate                          ! constant basal melt rate (m/s)
                                                    ! = config_bmlt_float_flux / (config_ice_density*latent_heat_ice)
      real(kind=RKIND) :: hCavity               ! depth of ice cavity beneath floating ice (m)
      real(kind=RKIND) :: zDraft                ! draft of floating ice (m below sea level)

      ! basal melting parameters for MISMIP+ experiment
      ! Note: These could be made user-configurable, but are hardwired for now because there are no plans
      !       to run MISMIP+ with different values
      real(kind=RKIND), parameter :: &
           bmlt_float_omega = 0.20_RKIND / scyr,  & ! time scale for basal melting (s^-1)
                                                    ! MISMIP+ default value = 0.2 yr^-1
           bmlt_float_h0 = 75._RKIND,             & ! scale for sub-shelf cavity thickness (m)
                                                    ! MISMIP+ default value = 75 m
           bmlt_float_z0 = -100._RKIND              ! scale for ice draft (m)
                                                    ! MISMIP+ default value = -100 m

      real(kind=RKIND) :: GLdepth, CFdepth ! ice shelf draft stats needed by the temperature_profile method:
                                           ! characteristic grounding line depth, calving front depth


      err = 0
      err_tmp = 0

      call mpas_pool_get_config(liConfigs, 'config_basal_mass_bal_float', config_basal_mass_bal_float)
      call mpas_pool_get_config(liConfigs, 'config_front_mass_bal_grounded', config_front_mass_bal_grounded)
      call mpas_pool_get_config(liConfigs, 'config_print_thermal_info', config_print_thermal_info)

      if (config_print_thermal_info) then
         call mpas_log_write('Solving for basal melting of floating ice, config_basal_mass_bal_float = ' // &
              trim(config_basal_mass_bal_float) )
      endif

      ! Simple options handled here
      if (trim(config_basal_mass_bal_float) == 'none') then

         ! Zero entire field

         ! block loop
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
            call mpas_pool_get_array(geometryPool, 'floatingBasalMassBal', floatingBasalMassBal)

            floatingBasalMassBal = 0.0_RKIND

            block => block % next
         enddo   ! associated(block)

         return ! skip the rest of this routine

      elseif (trim(config_basal_mass_bal_float) == 'file') then

         return   ! already set; nothing to do

      endif

      ! get rest of config variables
      call mpas_pool_get_config(liConfigs, 'config_thermal_thickness', config_thermal_thickness)
      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
      call mpas_pool_get_config(liConfigs, 'config_bmlt_float_flux', config_bmlt_float_flux)
      call mpas_pool_get_config(liConfigs, 'config_bmlt_float_xlimit', config_bmlt_float_xlimit)
      call mpas_pool_get_config(liConfigs, 'config_ice_density', config_ice_density)

      ! The temp-profile method needs a global reduction before it can operate, so do that here
      ! before the main option block
      if (trim(config_basal_mass_bal_float) == 'temperature_profile') then
         call calc_iceshelf_draft_info(domain, GLdepth, CFdepth, err_tmp)
         err = ior(err, err_tmp)
         if (config_print_thermal_info) then
            call mpas_log_write('GLdepth=$r, CFdepth=$r', realArgs=(/GLdepth, CFdepth/))
         endif
      endif

      if (trim(config_basal_mass_bal_float) == 'ismip6') then
         call iceshelf_melt_ismip6(domain, err_tmp)
         err = ior(err, err_tmp)
      else
         ! block loop
         block => domain % blocklist
         do while (associated(block))

            ! get pools
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
            call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
            call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

            ! get dimensions
            call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

            ! get fields from the mesh pool
            call mpas_pool_get_array(meshPool, 'xCell', xCell)
            call mpas_pool_get_array(meshPool, 'daysSinceStart',daysSinceStart)

            ! get fields from the geometry pool
            call mpas_pool_get_array(geometryPool, 'thickness', thickness)
            call mpas_pool_get_array(geometryPool, 'lowerSurface', lowerSurface)
            call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
            call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
            call mpas_pool_get_array(geometryPool, 'floatingBasalMassBal', floatingBasalMassBal)

            ! get fields from the scratch pool
            call mpas_pool_get_field(scratchPool, 'iceCellMask', thermalCellMaskField)
            call mpas_allocate_scratch_field(thermalCellMaskField, .true.)
            thermalCellMask => thermalCellMaskField % array

            ! calculate masks - so we know where the ice is floating
            call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
            err = ior(err, err_tmp)

            ! calculate a mask to identify ice that is thick enough to be thermally active
            do iCell = 1, nCellsSolve
               if (thickness(iCell) > config_thermal_thickness) then
                  thermalCellMask(iCell) = 1
               else
                  thermalCellMask(iCell) = 0
               endif
            enddo

            ! -----------------
            ! Compute basal melting for floating ice.
            ! -----------------

            ! initialize to zero melt
            floatingBasalMassBal(:) = 0.0_RKIND
            if (trim(config_basal_mass_bal_float) == 'none') then
               ! Do nothing, handled above

            elseif (trim(config_basal_mass_bal_float) == 'file') then
               ! Do nothing, handled above

            elseif (trim(config_basal_mass_bal_float) == 'constant') then

               ! set melt rate to a constant value for floating ice
               ! allow basal melt in ice-free ocean cells, in case ice is advected to those cells by the transport scheme

               bmlt_float_rate = config_bmlt_float_flux / (config_ice_density*latent_heat_ice)  ! convert W/m^2 to m/s

               floatingBasalMassBal(:) = 0.0_RKIND

               do iCell = 1, nCellsSolve
                  if ( (li_mask_is_floating_ice(cellMask(iCell))) .or.  &
                       (bedTopography(iCell) < config_sea_level .and. thermalCellMask(iCell) == 0) ) then
                     ! ice is present and floating, or ice-free ocean

                     ! Provided xCell > bmlt_float_xlimit, prescribe a uniform basal melt rate.
                     ! The default is 0.0, but for MISMIP+ the prescribed value of xlimit is 480 km.
                     if (abs(xCell(iCell)) >= config_bmlt_float_xlimit) then  ! basal melting is allowed
                        floatingBasalMassBal(iCell) = -bmlt_float_rate
                     endif

                  endif   ! ice is present and floating

               enddo   ! iCell

               ! change units from m/s to kg/m2/s
               floatingBasalMassBal(:) = floatingBasalMassBal(:) * config_ice_density

            elseif (trim(config_basal_mass_bal_float) == 'mismip') then

               ! compute melt rate (m/s) based on bed depth and cavity thickness
               ! The MISMIP+ formula is as follows:
               !
               ! bmlt_float = omega * tanh(H_c/H_0) * max(z_0 - z_d, 0)
               !
               ! where H_c = lsrf - topg is the cavity thickness
               !       z_d = lsrf - eus is the ice draft
               !       omega = a time scale = 0.2 yr^{-1} by default
               !       H_0 = 75 m by default
               !       z_0 = -100 m by default

               ! allow basal melt in ice-free ocean cells, in case ice is advected to those cells by the transport scheme

               floatingBasalMassBal(:) = 0.0_RKIND

               do iCell = 1, nCellsSolve

                  if ( (li_mask_is_floating_ice(cellMask(iCell))) .or.  &
                       (bedTopography(iCell) < config_sea_level .and. thermalCellMask(iCell) == 0) ) then
                       ! ice is present and floating, or ice-free ocean

                     hCavity = lowerSurface(iCell) - bedTopography(iCell)
                     zDraft = lowerSurface(iCell) - config_sea_level
                     floatingBasalMassBal(iCell) = -bmlt_float_omega * tanh(hCavity/bmlt_float_h0) * max(bmlt_float_z0 - &
                        zDraft, 0.0_RKIND)

                  endif   ! ice is present and floating
               enddo   ! iCell

               ! change units from m/s to kg/m2/s
               floatingBasalMassBal(:) = floatingBasalMassBal(:) * config_ice_density

            elseif (trim(config_basal_mass_bal_float) == 'seroussi') then

               call basal_melt_thwaites_seroussi(floatingBasalMassBal, daysSinceStart, lowerSurface, cellMask, &
                       config_sea_level, config_ice_density, nCellsSolve, err_tmp)
               err = ior(err, err_tmp)

            elseif (trim(config_basal_mass_bal_float) == 'temperature_profile') then

               call iceshelf_melt_param_from_temp_profile(nCellsSolve, lowerSurface, cellMask, &
                       config_sea_level, config_ice_density, config_print_thermal_info, &
                       GLdepth, CFdepth, daysSinceStart, floatingBasalMassBal, err_tmp)
               err = ior(err, err_tmp)

            else

               call mpas_log_write('Unknown option selected for config_basal_mass_bal_float:' // &
                  trim(config_basal_mass_bal_float), MPAS_LOG_ERR)
               err = ior(err, 1)

            endif   ! config_basal_mass_bal_float


            ! clean up
            call mpas_deallocate_scratch_field(thermalCellMaskField, .true.)

            block => block % next
         enddo   ! associated(block)

      endif

      ! Allocate scratch fields for flood-fill
      ! Note: This only supports one block per processor
      call mpas_pool_get_subpool(domain % blocklist % structs, 'scratch', scratchPool)
      call mpas_pool_get_field(scratchPool, 'seedMask', seedMaskField)
      call mpas_allocate_scratch_field(seedMaskField, single_block_in = .true.)
      seedMask => seedMaskField % array
      seedMask(:) = 0
      call mpas_pool_get_field(scratchPool, 'growMask', growMaskField)
      call mpas_allocate_scratch_field(growMaskField, single_block_in = .true.)
      growMask => growMaskField % array
      growMask(:) = 0

      ! Prevent sub-shelf melt in subglacial lakes
      call mpas_pool_get_subpool(domain % blocklist % structs, 'geometry', geometryPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_array(geometryPool, 'floatingBasalMassBal', floatingBasalMassBal)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
      do iCell = 1, nCells
         if ( (.not. li_mask_is_ice(cellMask(iCell))) .and. (bedTopography(iCell) < config_sea_level) ) then
              seedMask(iCell) = 1
         else if ( li_mask_is_floating_ice(cellMask(iCell)) ) then
              growMask(iCell) = 1
         end if
      end do
      call li_flood_fill(seedMask, growMask, domain)

      floatingBasalMassBal(:) = floatingBasalMassBal(:) * seedMask(:)

      ! deallocate scratch fields used for flood fill
      call mpas_deallocate_scratch_field(seedMaskField, single_block_in=.true.)
      call mpas_deallocate_scratch_field(growMaskField, single_block_in=.true.)

      if (trim(config_front_mass_bal_grounded) .ne. 'none') then
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
            call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
            call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
            call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
            call mpas_pool_get_array(geometryPool, 'floatingBasalMassBal', floatingBasalMassBal)
            call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
            call mpas_pool_get_array(geometryPool, 'edgeMask', edgeMask)

            ! If config_front_mass_bal_grounded is not none, only apply ice shelf melt to active cells
            ! and stranded non-dynamic cells. i.e., if a floating non-dynamic
            ! cell has a grounded neighbor, do not apply ice shelf melt to it.
            do iCell = 1, nCellsSolve
               if ( li_mask_is_floating_ice(cellMask(iCell)) .and. (.not. li_mask_is_dynamic_ice(cellMask(iCell))) ) then
                  do iNeighbor = 1, nEdgesOnCell(iCell)
                     jCell = cellsOnCell(iNeighbor, iCell)
                     if ( li_mask_is_grounded_ice(cellMask(jCell)) ) then
                         floatingBasalMassBal(iCell) = 0.0_RKIND
                         exit ! No need to look over other neighbors
                      endif
                  enddo
                endif
            enddo

            block => block % next
         enddo   ! associated(block)
      endif

    end subroutine li_basal_melt_floating_ice
!-----------------------------------------------------------------------


!***********************************************************************
!***********************************************************************
! Private subroutines:
!***********************************************************************
!***********************************************************************



!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  !  routine basal_melt_thwaites_seroussi
!
!> \brief Calculate ice shelf melt rate from depth param.
!> \author William Lipscomb
!> \date   November 2015
!> \details
!>  Melt rate parameterization from:
!> Seroussi, H., Y. Nakayama, E. Larour, D. Menemenlis, M. Morlighem, E. Rignot, and A. Khazendar (2017),
!> Continued retreat of Thwaites Glacier, West Antarctica, controlled by bed topography and ocean circulation,
!> Geophys. Res. Lett., 1-9, doi:10.1002/2017GL072910.
!> for Thwaites Glacier.
!> Specifically, this is a linear fit of melt with shelf draft from the Supplemental Information, Figure S1.
!> The linear relation is modified by a:
!>   * depth above which there is no melt (Antarctic Surface Water saturation)
!>   * a maximum melt rate (Circumpolar Deep Water saturation)
!>   * a depth below which melt stops increasing (minimum sill height)

!-----------------------------------------------------------------------

    subroutine basal_melt_thwaites_seroussi(floatingBasalMassBal, daysSinceStart, lowerSurface, cellMask, &
                    config_sea_level, config_ice_density, nCellsSolve, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real(kind=RKIND), pointer, intent(in) :: daysSinceStart
      real (kind=RKIND), dimension(:), pointer, intent(in) :: lowerSurface  !< lower surface elevation (m)
      integer, dimension(:), pointer, intent(in) :: &
           cellMask                    !< bit mask describing whether ice is floating, dynamically active, etc.
      real(kind=RKIND), pointer, intent(in) :: &
           config_sea_level !< sea level (m) relative to z = 0
      real (kind=RKIND), pointer, intent(in) :: config_ice_density !< ice density
      integer, pointer :: &
           nCellsSolve                 !< number of locally owned cells

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      real (kind=RKIND), dimension(:), pointer, intent(out) :: &
           floatingBasalMassBal      !< basal mass balance for floating ice
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      real (kind=RKIND) :: slopeSer  ! slope of relation between depth and melt rate
      real (kind=RKIND) :: interceptSer  ! depth at which melting goes to 0
      real (kind=RKIND) :: maxMeltSer  ! maximum allowable melt rate
      real (kind=RKIND) :: sillDepth  ! depth below which melt rate no longer increases
      real (kind=RKIND), pointer :: config_basal_mass_bal_seroussi_amplitude
      real (kind=RKIND), pointer :: config_basal_mass_bal_seroussi_period
      real (kind=RKIND), pointer :: config_basal_mass_bal_seroussi_phase
      real(kind=RKIND) :: hCavity               ! depth of ice cavity beneath floating ice (m)
      real(kind=RKIND) :: zDraft                ! draft of floating ice (m below sea level)
      integer :: iCell


      err = 0

      call mpas_pool_get_config(liConfigs, 'config_basal_mass_bal_seroussi_amplitude', &
         config_basal_mass_bal_seroussi_amplitude) ! meters
      call mpas_pool_get_config(liConfigs, 'config_basal_mass_bal_seroussi_period', &
         config_basal_mass_bal_seroussi_period) ! years
      call mpas_pool_get_config(liConfigs, 'config_basal_mass_bal_seroussi_phase', &
         config_basal_mass_bal_seroussi_phase) ! cycles

      slopeSer = 0.088_RKIND ! slope of relation between depth and melt rate (melt (m/yr) per depth (m))
      interceptSer = -100.0_RKIND  ! depth (m) at which melting goes to 0 (negative meaning below sea level)
      maxMeltSer = 50.0_RKIND ! maximum allowable melt rate (m/yr) (positive meaning melting)
      sillDepth = -650.0_RKIND ! depth below which melt stops increasing (m) (negative meaning below sea level)

      if (config_basal_mass_bal_seroussi_period <= 0.0_RKIND) then
         call mpas_log_write("Value for config_basal_mass_bal_seroussi_period has to be a positive real value.", MPAS_LOG_ERR)
         err = ior(err, 1)
      endif

      ! Modify intercept height for variability parameters
      interceptSer = interceptSer + config_basal_mass_bal_seroussi_amplitude * &
              sin( (2.0_RKIND * pii / config_basal_mass_bal_seroussi_period) * (daysSinceStart/365.0_RKIND) &
                  + 2.0_RKIND * pii * config_basal_mass_bal_seroussi_phase)

      ! Initialize before computing
      floatingBasalMassBal(:) = 0.0_RKIND

      do iCell = 1, nCellsSolve

         ! Shut off melt at an arbitrary shallow depth to discourage ice from disappearing.
         if ( (li_mask_is_floating_ice(cellMask(iCell))) .and. (lowerSurface(iCell) < -10.0_RKIND) ) then
            ! ice is present and floating

            zDraft = lowerSurface(iCell) - config_sea_level
            ! Coefficients for m/yr melt rate (in units of Seroussi figure but without negative meaning melting)
            floatingBasalMassBal(iCell) = max(-1.0_RKIND * maxMeltSer, min(0.0_RKIND, slopeSer * &
               (max(zDraft, sillDepth) - interceptSer)))

         endif   ! ice is present
      enddo   ! iCell

      ! change units from m/yr to kg/m2/s
      floatingBasalMassBal(:) = floatingBasalMassBal(:) * config_ice_density / scyr


    end subroutine basal_melt_thwaites_seroussi
!-----------------------------------------------------------------------


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  !  routine iceshelf_melt_param_from_temp_profile
!
!> \brief Calculate ice shelf melt rate using parameterization based on ocean temperature profile
!> \author Matthew Hoffman, based on parameterization developed by Xylar Asay-Davis
!> \date   October 2018
!> \details  Melt rate parameterization by Xylar Asay-Davis used in Thwaites variability paper
!> in prep. by Hoffman et al.
!> Parameterization takes a ocean temperature profile and sill depth as input.
!> The temp. profile is adjusted for cavity temperature behind the sill,
!> assuming constant temp. below the sill depth.
!> Then a local temp profile in the boundary layer is estimated assuming exponential
!> decay from the cavity profile to a plume profile at infinite distance downstream
!> from the grounding line, accounting for the pressure-melt temperature.
!> The melt rate is calculated as proportional to the  product of the boundary layer
!> temperature profile and the mean cavity temperature.
!-----------------------------------------------------------------------

    subroutine iceshelf_melt_param_from_temp_profile(nCellsSolve, lowerSurface, cellMask, &
                    config_sea_level, config_ice_density, print_debug, zGL, zCavityMax, daysSinceStart, floatingBasalMassBal, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      integer, intent(in) :: nCellsSolve !< number of locally owned cells
      real (kind=RKIND), dimension(:), intent(in) :: lowerSurface !< ice shelf draft
      integer, dimension(:), intent(in) :: cellMask !< bitmask on cells
      real (kind=RKIND), intent(in) :: config_sea_level, config_ice_density
      logical, intent(in) :: print_debug
      real(kind=RKIND), intent(in) :: zGL, &  !< mean grounding line elevation (m),
                                              !< it's a scalar now, but eventually may want smoothed spatial field
                                      zCavityMax !< depth below which mean thermal forcing is averaged
      real(kind=RKIND), pointer, intent(in) :: daysSinceStart

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      real (kind=RKIND), dimension(:), intent(out) :: floatingBasalMassBal !< BMB for ice shelves
      integer, intent(out) :: err


      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      ! parameters for param.
      real(kind=RKIND), pointer :: slopeDraft !< characteristic slope of ice shelf draft
      real(kind=RKIND), pointer :: plumeThickness !< characteristic plume thickness (m)
      real(kind=RKIND), parameter :: E0 = 0.036_RKIND !< entrainment parameter (Jenkins 1991)
      real(kind=RKIND), parameter :: Stanton = 5.9e-4_RKIND !< Stanton number (Jenkins et al. 2010)

      ! All depths are written as elevation relative to sea level (negative values)
      ! variables describing regional profile
      real(kind=RKIND), pointer :: zThermoclineTopBase
         !< depth at top of thermocline (temp assumed constant above this), before variability adjustment
      real(kind=RKIND), pointer :: zThermoclineBotBase
         !< depth at bottom of thermocline (temp assumed constant below this), before variability adjustment
      real(kind=RKIND) :: zThermoclineTop
         !< depth at top of thermocline (temp assumed constant above this), after variability adjustment
      real(kind=RKIND) :: zThermoclineBot
         !< depth at bottom of thermocline (temp assumed constant below this), after variability adjustment
      real(kind=RKIND), pointer :: TThermoclineTop !< temperature at top of thermocline (constant shallow water mass temp.)
      real(kind=RKIND), pointer :: TThermoclineBot !< temperature at bottom of thermocline (constant deep water mass temp.)

      real(kind=RKIND), pointer :: variabilityAmplitude !< amplitude in optional sinusoidal variab. applied to thermocline depths
      real(kind=RKIND), pointer :: variabilityPeriod !< period in optional sinusoidal variability applied to thermocline depths
      real(kind=RKIND), pointer :: variabilityPhase !< phase in optional sinusoidal variability applied to thermocline depths
      real(kind=RKIND) :: depthOffset !< depth offset in thermocline depths due to variability

      ! variables that are general to any depth
      real(kind=RKIND), pointer :: zSill !< sill depth
      real(kind=RKIND) :: TSill !< temperature at the sill depth
      real(kind=RKIND) :: TGL !< temperature at grounding line depth
      real(kind=RKIND) :: TFmean !< mean thermal forcing within the cavity
      real(kind=RKIND), pointer :: kappa !< constant for non-local, quadratic melt param.
      ! variables specific to a given depth
      real(kind=RKIND) :: draft !< elevation of ice shelf draft (m)
      real(kind=RKIND) :: Tregional !< regional temperature at given depth (comes from imposed profile)
      real(kind=RKIND) :: Tcavity !< temperature in cavity at given depth
      real(kind=RKIND) :: Tfreeze !< freezing temperature at given depth
      real(kind=RKIND) :: Tinfinity !< limit of exponential decay of T at given depth
      real(kind=RKIND) :: Tplume !< temperature in plume at a given depth
      real(kind=RKIND) :: zeta !< length scale of decay of plume temperature

      ! other
      integer :: nLevels
      integer :: iCell, z
      real(kind=RKIND) :: zr

      err = 0

      ! get namelist options
      call mpas_pool_get_config(liConfigs, 'config_temperature_profile_melt_scale_factor', kappa)
      call mpas_pool_get_config(liConfigs, 'config_temperature_profile_sill_elevation', zSill)
      call mpas_pool_get_config(liConfigs, 'config_temperature_profile_plume_thickness', plumeThickness)
      call mpas_pool_get_config(liConfigs, 'config_temperature_profile_draft_slope', slopeDraft)
      ! Get temp. profile parameters
      call mpas_pool_get_config(liConfigs, 'config_temperature_profile_thermocline_upper_depth', zThermoclineTopBase)
      call mpas_pool_get_config(liConfigs, 'config_temperature_profile_thermocline_upper_temp', TThermoclineTop)
      call mpas_pool_get_config(liConfigs, 'config_temperature_profile_thermocline_lower_depth', zThermoclineBotBase)
      call mpas_pool_get_config(liConfigs, 'config_temperature_profile_thermocline_lower_temp', TThermoclineBot)
      ! Get variability parameters
      call mpas_pool_get_config(liConfigs, 'config_temperature_profile_variability_amplitude', variabilityAmplitude)
      call mpas_pool_get_config(liConfigs, 'config_temperature_profile_variability_period', variabilityPeriod)
      call mpas_pool_get_config(liConfigs, 'config_temperature_profile_variability_phase', variabilityPhase)


      ! Account for variability in thermocline depth
      if (variabilityPeriod <= 0.0_RKIND) then
         call mpas_log_write("Value for config_temperature_profile_variability_period must be a positive real value.", &
                 MPAS_LOG_ERR)
         err = ior(err, 1)
      endif
      depthOffset = variabilityAmplitude * sin( (2.0_RKIND * pii / variabilityPeriod) * (daysSinceStart/365.0_RKIND) &
                    + 2.0_RKIND * pii * variabilityPhase)
      zThermoclineTop = zThermoclineTopBase + depthOffset
      zThermoclineBot = zThermoclineBotBase + depthOffset

      ! Calculate variables independent of depth
      zeta = slopeDraft * plumeThickness / (E0 * slopeDraft + Stanton)

      ! calculate temperature at sill
      if (zSill >= zThermoclineTop) then
         TSill = TThermoclineTop
      elseif (zSill <= zThermoclineBot) then
         TSill = TThermoclineBot
      else
         TSill = (TThermoclineTop - TThermoclineBot) * &
                 (zSill - zThermoclineBot) / (zThermoclineTop - zThermoclineBot) + TThermoclineBot
      endif
      if (print_debug) then
         call mpas_log_write("kappa=$r, plumeThickness=$r, slopeDraft=$r, zeta=$r", realArgs=(/kappa, plumeThickness, &
                 slopeDraft, zeta/))
         call mpas_log_write("TSill=$r at zSill=$r", realArgs=(/TSill, zSill/))
      endif

      ! calculate temperature at GL (allowing for no sill - sill deeper than GL depth)
      if (zGL < zSill) then
         TGL = TSill
      elseif (zGL >= zThermoclineTop) then
         TGL = TThermoclineTop
      elseif (zGL <= zThermoclineBot) then
         TGL = TThermoclineBot
      else
         TGL = (TThermoclineTop - TThermoclineBot) * &
                 (zGL - zThermoclineBot) / (zThermoclineTop - zThermoclineBot) + TThermoclineBot
      endif
      if (print_debug) then
         call mpas_log_write("TGL=$r at zGL=$r", realArgs=(/TGL, zGL/))
      endif

      ! Need to calculate a mean thermal forcing within the plume
      ! Ideally this would be a function of the ice shelf hypsometry,
      ! but for simplicity calculate assuming a linear (ramp) ice shelf shape.
      ! Average between the GL depth and an arbitrary ice shelf thickness at calving front.
      ! (We think that is better than averaging all the way to surface, even if that thin
      !  shelf thickness does not match that of the modeled ice shelf.)
      TFmean = 0.0_RKIND
      nLevels = 0
      do z = ceiling(zGL), floor(zCavityMax), 1 ! Discretize in 1 m increments
         zr = real(z, kind=RKIND) ! convert once to real data type
         nLevels = nLevels + 1
         ! calculate Tcavity for this depth (allowing for no sill - sill deeper than GL depth)
         if (zr < zSill) then
            Tcavity = TSill
         elseif (zr >= zThermoclineTop) then
            Tcavity = TThermoclineTop
         elseif (zr <= zThermoclineBot) then
            Tcavity = TThermoclineBot
         else
            Tcavity = (TThermoclineTop - TThermoclineBot) * &
                 (zr - zThermoclineBot) / (zThermoclineTop - zThermoclineBot) + TThermoclineBot
         endif
         ! calculate Tfreeze for this depth
         Tfreeze = ocean_freezing_temperature(zr)
         ! calculate Tinfinity for this depth
         Tinfinity = (E0 * slopeDraft * Tcavity + Stanton * Tfreeze) / (E0 * slopeDraft + Stanton)
         ! calculate Tplume for this depth
         Tplume = Tinfinity + (TGL - Tinfinity) * exp(min(-1.0_RKIND * (zr - zGL) / zeta, 0.0_RKIND))
         ! now add to sum of TFmean
         TFmean = TFmean + (Tplume - Tfreeze)
      enddo
      ! now find the average
      if (nLevels > 0) then
         TFmean = TFmean / real(nLevels, kind=RKIND)
      else
         call mpas_log_write("iceshelf_melt_param_from_temp_profile found 0 levels over which to average mean thermal forcing!",&
               MPAS_LOG_ERR)
         err = ior(err, 1)
      endif
      if (print_debug) then
         call mpas_log_write("TFmean=$r", realArgs=(/TFmean/))
      endif

      ! Now calculate melt rate at each grid cell
      do iCell = 1, nCellsSolve
         if (li_mask_is_floating_ice(cellMask(iCell))) then
            if (lowerSurface(iCell) > -20.0_RKIND) then
               ! Shut off melt at an arbitrary shallow depth to prevent ice from melting out entirely
               floatingBasalMassBal(iCell) = 0.0_RKIND
            else
               ! Calculate melt rate here
               draft = lowerSurface(iCell) - config_sea_level
               ! calculate Tcavity for this depth (allowing for no sill - sill deeper than GL depth)
               if (draft < zSill) then
                  Tcavity = TSill
               elseif (draft >= zThermoclineTop) then
                  Tcavity = TThermoclineTop
               elseif (draft <= zThermoclineBot) then
                  Tcavity = TThermoclineBot
               else
                  Tcavity = (TThermoclineTop - TThermoclineBot) * &
                       (draft - zThermoclineBot) / (zThermoclineTop - zThermoclineBot) + TThermoclineBot
               endif
               ! calculate Tfreeze for this depth
               Tfreeze = ocean_freezing_temperature(draft)
               ! calculate Tinfinity for this depth
               Tinfinity = (E0 * slopeDraft * Tcavity + Stanton * Tfreeze) / (E0 * slopeDraft + Stanton)
               ! calculate Tplume for this depth
               Tplume = Tinfinity + (TGL - Tinfinity) * exp(min(-1.0_RKIND * (draft - zGL) / zeta, 0.0_RKIND))
               ! now add to sum of TFmean
               floatingBasalMassBal(iCell) = kappa * TFmean * (Tplume - Tfreeze) ! in m/yr, +=melting
               !floatingBasalMassBal(iCell) = Tplume
            endif
         else ! not floating
            floatingBasalMassBal(iCell) = 0.0_RKIND
         endif
      enddo
      floatingBasalMassBal(:) = floatingBasalMassBal(:) * (-1.0_RKIND) * config_ice_density / scyr ! convert to BMB

    end subroutine iceshelf_melt_param_from_temp_profile
!-----------------------------------------------------------------------


    function ocean_freezing_temperature(depth) result(temp)
      real(kind=RKIND), intent(in) :: depth
      real(kind=RKIND), parameter :: Sref = 34.4_RKIND  !< Reference salinity (made up, typical for MISOMIP1)
      real(kind=RKIND) :: temp

      ! function to calculate freezing temperature from Jenkins (1991)
      temp = 0.0901_RKIND - 0.0575_RKIND * SRef + 7.61e-4_RKIND * depth

    end function ocean_freezing_temperature
!-----------------------------------------------------------------------


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  !  routine calc_iceshelf_draft_info
!
!> \brief Calculate characteristic grounding line and calving front depths underneath ice shelves
!> \author Matthew Hoffman
!> \date   October 2018
!> \details Grounding line and calving front depths are used by the
!> temperature_profile iceshelf melt method, but must be calculated before the
!> method is called because they require a global reduce.
!> Different options for how to calculated these depths could be implemented.
!> For the GL depth, an approximation of the mean depth of deepest specified-fraction of the
!> grounding line is used.  For the CF depth, the mean of the calving front
!> is used.
!-----------------------------------------------------------------------

    subroutine calc_iceshelf_draft_info(domain, GLdepth, CFdepth, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: &
         domain          !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      real (kind=RKIND), intent(out) :: GLdepth !< mean grounding line depth
      real (kind=RKIND), intent(out) :: CFdepth !< mean calving front depth
      integer, intent(out) :: err

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: geometryPool

      integer, pointer :: &
           nCellsSolve                 ! number of locally owned cells

      logical, pointer :: &
           config_print_thermal_info   ! if true, print debug info

      real(kind=RKIND), pointer :: &
           config_sea_level ! sea level (m) relative to z = 0

      integer, dimension(:), pointer :: &
           cellMask                    ! bit mask describing whether ice is floating, dynamically active, etc.

      real (kind=RKIND), dimension(:), pointer :: lowerSurface ! lower surface elevation (m)

      real (kind=RKIND), dimension(:), pointer :: areaCell ! cell area (m^2)
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell

      real(kind=RKIND) :: sumCFdepth, sumCFarea
      integer :: iCell, jCell
      real(kind=RKIND), dimension(2) :: localCFSums, globalCFSums ! Vars for CF global reduction

      real(kind=RKIND), parameter :: binWidth = 10.0_RKIND !< width of depth bins in meters
      real(kind=RKIND), parameter :: binMinimum = -4000.0_RKIND !< depth of bottom of deepest bin in meters
      integer, parameter :: nBinsMax = 500 !< Maximum possible number of bins - needs to be set to cover the range
                                        !< defined by the previous two variables and sea level (which is run-time configurable)
      integer :: nBins !< number of bins - calculated below
      real(kind=RKIND), dimension(nBinsMax) :: binDepth !< depth in m at bottom of each bin
      real(kind=RKIND), dimension(nBinsMax) :: binWeightSum !< sum of cell weights in this bin
      real(kind=RKIND), dimension(nBinsMax) :: binWeightedDepthSum !< weighted-depth sum over cells in this bin
      real(kind=RKIND), dimension(2*nBinsMax) :: localGLSums, globalGLSums ! Vars for GL global reduction
      real(kind=RKIND) :: cumBinWeightedDepthSum, cumBinWeightSum !< cumulative sums over bins for final calculation
      real(kind=RKIND) :: totalBinWeight !< sum of weights in all bins
      real(kind=RKIND), pointer :: weightThreshold !< fraction of deepest GL depths to average overn

      integer :: indexBin


      err = 0

      call mpas_pool_get_config(liConfigs, 'config_temperature_profile_GL_depth_fraction', weightThreshold)

      ! need to first determine bin info (could happen on init, but cheap to recompute each timestep)
      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
      nBins = ceiling((config_sea_level - binMinimum) / binWidth)  ! round up
      do indexBin = 1, nBins
         binDepth = binMinimum + real(indexBin-1, kind=RKIND) * binWidth ! calc. lower limit of each bin
      enddo

      ! Calculate average GL and CF depths on this processor, looping over blocks.
      ! Note that these are somewhat approximate calculations.
      ! Note that the GL in cellMask includes grounded cells next to open ocean,
      ! which we do NOT want included, so instead we here identify the GL depth at
      ! cells that are floating and have a grounded neighbor.
      ! (Also, the GL cell mask identifies the last grounded cell, not the first floating cell,
      ! which also isn't what we want.)
      ! Also, ideally the weighted average for GL and CF depths would be weighted
      ! based on GL and CF lengths.  That is substantially more work for little
      ! expected improvement (and we would have to make a decision about what elevation to
      ! assign to edges), so I instead use cell area weighted averages.
      ! The GL and CF depths do not need to be highly accurate, so these approximations
      ! are assumed to be fine.
      binWeightSum(:) = 0.0_RKIND
      binWeightedDepthSum(:) = 0.0_RKIND
      sumCFdepth = 0.0_RKIND
      sumCFarea = 0.0_RKIND
      ! init bins for GL calcs
      ! block loop
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_array(geometryPool, 'lowerSurface', lowerSurface)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
         call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

         do iCell = 1, nCellsSolve
            ! First calc. GL stat
            if (li_mask_is_floating_ice(cellMask(iCell)) .and. li_mask_is_dynamic_ice(cellMask(iCell))) then
               ! See if a neighbor is grounded ice
               do jCell = 1, nEdgesOnCell(iCell)
                  if (li_mask_is_grounded_ice(cellMask(cellsOnCell(jCell, iCell)))) then ! we found the first floating cell
                     ! check for reasonable values first
                     if (lowerSurface(iCell) < binMinimum) then
                        call mpas_log_write("calc_iceshelf_draft_info found a grounding line depth below the minimum search " // &
                              "depth of $r.  Adjust binMinimum and related parameters and recompile.", &
                              MPAS_LOG_ERR, realArgs = (/binMinimum/))
                        err = ior(err, 1)
                        return
                     elseif (lowerSurface(iCell) > config_sea_level) then
                        call mpas_log_write("calc_iceshelf_draft_info found a grounding line depth above sea level."&
                                // "  Adjust hard-coded parameters and recompile.", MPAS_LOG_ERR)
                        err = ior(err, 1)
                        return
                     endif
                     ! Find which bin this location falls in
                     indexBin = floor((lowerSurface(iCell) - binMinimum) / binWidth) + 1


                     ! Use the sqrt of cell area as the weighting function (it should be a linear measure, not an area measure)
                     binWeightedDepthSum(indexBin) = binWeightedDepthSum(indexBin) + lowerSurface(iCell) * sqrt(areaCell(iCell))
                     binWeightSum(indexBin) = binWeightSum(indexBin) + sqrt(areaCell(iCell))
                     !call mpas_log_write("GL cell: $i, $r, $r", intArgs=(/iCell/), &
                     !realArgs=(/lowerSurface(iCell), areaCell(iCell)/))
                     exit  ! skip the rest of the loop around neighbors
                  endif
               enddo
            endif
            ! Now calc. CF stat
            if (li_mask_is_floating_ice(cellMask(iCell)) .and. li_mask_is_dynamic_margin(cellMask(iCell))) then
               ! accumulate area-weighted CF depth and area
               sumCFdepth = sumCFdepth + lowerSurface(iCell) * areaCell(iCell)
               sumCFarea = sumCFarea + areaCell(iCell)
            endif
         enddo

         block => block % next
      enddo   ! associated(block)

      ! --- Global reduction for GL stat ---
      ! pack communication array
      localGLSums(1:nBins) = binWeightedDepthSum(1:nBins)
      localGLSums(nBins+1:2*nBins) = binWeightSum(1:nBins)
      ! global reduction
      call mpas_dmpar_sum_real_array(domain % dminfo, 2*nBins, localGLSums(:), globalGLSums(:))
      ! unpack communication array
      binWeightedDepthSum(1:nBins) = globalGLSums(1:nBins)
      binWeightSum(1:nBins) = globalGLSums(nBins+1:2*nBins)
      ! sum up through bins until threshold is reached
      totalBinWeight = sum(binWeightSum(1:nBins))
      if (totalBinweight > 0.0_RKIND) then
         cumBinWeightedDepthSum = 0.0_RKIND
         cumBinWeightSum = 0.0_RKIND
         do indexBin = 1, nBins
            cumBinWeightedDepthSum = cumBinWeightedDepthSum + binWeightedDepthSum(indexBin)
            cumBinWeightSum = cumBinWeightSum + binWeightSum(indexBin)
            ! check if we've reached threshold
            if (cumBinWeightSum / totalBinweight > weightThreshold) then
               ! We reached the desired threshold of deepest samples.  Make final calculation and quit.
               GLdepth = cumBinWeightedDepthSum / cumBinWeightSum
               exit  ! skip the rest of the loop over bins
            endif
         enddo
      else
              ! Most likely situation: no GL!  So do nothing.
      endif

      ! --- Global reduction for CF stat ---
      localCFSums = (/sumCFdepth, sumCFarea/)
      ! sum up over all procs
      call mpas_dmpar_sum_real_array(domain % dminfo, 2, localCFSums(:), globalCFSums(:))
      if (globalCFSums(2) > 0.0_RKIND) then
         CFdepth = globalCFSums(1) / globalCFSums(2)
      else
         CFdepth = 0.0_RKIND
      endif

    end subroutine calc_iceshelf_draft_info
!-----------------------------------------------------------------------


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  !  routine iceshelf_melt_ismip6
!
!> \brief Calculate ice shelf melt rate using parameterization from ISMIP6
!> \author Matthew Hoffman
!> \date   July 2019
!> \details
!> http://www.climate-cryosphere.org/wiki/index.php?title=ISMIP6-Projections-Antarctica
!> It is a quadratic method of the form:
!> melt = gamma * physicalconstant * (TF(x,z,y_draft) + deltaT_basin) * |<TF> + deltatT_basin|
!> where TF is local thermal forcing at each grid cell and <TF> is the mean
!> thermal forcing for all of the ice shelf grid cells in an entire basin
!-----------------------------------------------------------------------

    subroutine iceshelf_melt_ismip6(domain, err)

      use li_constants, only: oceanFreezingTempDepthDependence

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: &
         domain          !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err


      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      integer :: iCell
      integer, pointer :: nOceanLayers
      integer :: ksup, kk, kinf
      integer, pointer :: nCellsSolve
      real(kind=RKIND), pointer :: rhoi, rhosw
      real (kind=RKIND), pointer :: invalid_value_TF
      real(kind=RKIND) :: cste
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: geometryPool, meshPool
      real (kind=RKIND), dimension(:,:), pointer :: TFocean
      real (kind=RKIND), dimension(:), pointer :: zOcean
      real (kind=RKIND), dimension(:), pointer :: TFdraft
      real (kind=RKIND), dimension(:), pointer :: deltaT
      real (kind=RKIND), dimension(:), pointer :: ismip6shelfMelt_offset
      integer, dimension(:), pointer :: basinNumber
      real (kind=RKIND), dimension(:), pointer :: lowerSurface
      real (kind=RKIND), dimension(:), pointer :: areaCell
      real (kind=RKIND), dimension(:), allocatable :: mean_TF, IS_area
      integer, parameter :: maxNumBasins = 32
      integer, dimension(:), pointer :: cellMask, indexToCellID
      real(kind=RKIND), dimension(maxNumBasins*2) :: localSums, globalSums
      real (kind=RKIND), pointer :: gamma0
      real (kind=RKIND), dimension(:), pointer :: floatingBasalMassBal
      real (kind=RKIND) :: coef
      integer :: i

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_ice_density', rhoi)
      call mpas_pool_get_config(liConfigs, 'config_ocean_density', rhosw)
      call mpas_pool_get_config(liConfigs, 'config_invalid_value_TF', invalid_value_TF)
      cste = (rhosw*cp_seawater/(rhoi*latent_heat_ice))**2  ! in K^(-2)

      allocate(mean_TF(maxNumBasins))
      allocate(IS_area(maxNumBasins))

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(meshPool, 'nISMIP6OceanLayers', nOceanLayers)
         call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)

         ! Get 3d thermal forcing (had to be read-in)
         call mpas_pool_get_array(geometryPool, 'ismip6shelfMelt_3dThermalForcing', TFocean)

         ! Get z-coords for ocean data
         call mpas_pool_get_array(geometryPool, 'ismip6shelfMelt_zOcean', zOcean)

         ! Get basin numbers
         call mpas_pool_get_array(geometryPool, 'ismip6shelfMelt_basin', basinNumber)

! NOTE: Would need to do a global max for this check to work
!         if (maxval(basinNumber)+1 > maxNumBasins) then
!            call mpas_log_write("Field ismip6shelfMelt_basin contains a value greater than the hard-coded value" // &
!               "for maxNumBasins in iceshelf_melt_ismip6.  Please adjust maxNumBasins and recompile or correct input data.", &
!               MPAS_LOG_ERR)
!         end if

         ! Get shelf draft and mask
         call mpas_pool_get_array(geometryPool, 'lowerSurface', lowerSurface)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

         ! Get TFdraft field that we will calculate
         call mpas_pool_get_array(geometryPool, 'ismip6shelfMelt_TFdraft', TFdraft)

         mean_TF(:) = 0.d0
         IS_area(:) = 0.d0

         ! Check zOcean for valid values
         if (minval(zOcean) >= 0.0_RKIND) then
            call mpas_log_write("Invalid value for zOcean. It should have negative values but min value >= 0.0 was found", &
                                MPAS_LOG_ERR)
            err = ior(err, 1)
         endif
         if (maxval(zOcean) > 0.0_RKIND) then
            call mpas_log_write("Invalid value for zOcean. It should have non-positive values but max value greater than 0.0 was found", &
                                MPAS_LOG_ERR)
            err = ior(err, 1)
         endif

         do iCell = 1, nCellsSolve

            if ( li_mask_is_floating_ice(cellMask(iCell)) ) then
               ! 1 -  Linear interpolation of the thermal forcing on the ice draft depth :
               ksup=0
               do kk=1,nOceanLayers
                  if ( zOcean(kk) >= lowerSurface(iCell) )  ksup = kk
               enddo
               kinf = ksup + 1

               !call mpas_log_write("kinf=$i, zOcean(kinf)=$r, TFocean=$r",realArgs=(/zOcean(kinf),TFocean(kinf,iCell)/), &
               !        intArgs=(/kinf/))
               !call mpas_log_write("ksup=$i, zOcean(ksup)=$r, TFocean=$r",realArgs=(/zOcean(ksup),TFocean(ksup,iCell)/), &
               !        intArgs=(/ksup/))

               ! check if any invalid TFocean value is used for calculating TF at the draft
               if (ksup >= 1) then
                  if (TFocean(ksup, iCell) == invalid_value_TF) then
                     call mpas_log_write("iceshelf_melt_ismip6: Invalid value for TFocean. " // &
                        "ksup=$i TFocean(ksup, iCell)=$r indexToCellID=$i", MPAS_LOG_ERR, &
                        intArgs=(/ksup,  indexToCellID(iCell)/), &
                        realArgs=(/TFocean(ksup, iCell) /) )
                     err = ior(err, 1)
                  endif
               endif
               if (kinf <= nOceanLayers) then
                  if (TFocean(kinf, iCell) == invalid_value_TF) then
                     call mpas_log_write("iceshelf_melt_ismip6: Invalid value for TFocean. " // &
                        "kinf=$i TFocean(kinf,iCell)=$r indexToCellID=$i", MPAS_LOG_ERR, &
                        intArgs=(/kinf, indexToCellID(iCell)/), &
                        realArgs=(/TFocean(kinf, iCell) /) )
                     err = ior(err, 1)
                  endif
               endif

               if (ksup == 0) then
                  ! For depths shallower than shallowest layer center, use shallowest layer
                  TFdraft(iCell) = TFocean(1, iCell)
               elseif (ksup == nOceanLayers) then
                  ! for depths below the deepest layer center, use deepest layer corrected for Tfreeze
                  TFdraft(iCell) = TFocean(nOceanLayers, iCell) - &
                     (zOcean(nOceanLayers) - lowerSurface(iCell)) * oceanFreezingTempDepthDependence
               else
                  ! for depths between the first and last layer centers, linearly interpolate
                  TFdraft(iCell) = ( (zOcean(ksup)-lowerSurface(iCell)) * TFocean(kinf, iCell) &
                     + (lowerSurface(iCell)-zOcean(kinf)) * TFocean(ksup, iCell) ) / (zOcean(ksup)-zOcean(kinf))
               endif

               ! 2 -  Mean Thermal forcing in individual basins (NB: fortran norm while basins start at zero):
               mean_TF(basinNumber(iCell)+1) = mean_TF(basinNumber(iCell)+1) + areaCell(iCell) * TFdraft(iCell)
               IS_area(basinNumber(iCell)+1) = IS_area(basinNumber(iCell)+1) + areaCell(iCell)

            else
               TFdraft(iCell) = 0.d0
            endif
         enddo

         block => block % next
      enddo   ! associated(block)

      ! global sum mean_TF and IS_area
      localSums(1:maxNumBasins) = mean_TF(:)
      localSums(maxNumBasins+1:2*maxNumBasins) = IS_area(:)
      call mpas_dmpar_sum_real_array(domain % dminfo, 2*maxNumBasins, localSums(:), globalSums(:))
      mean_TF(:) = globalSums(1:maxNumBasins)
      IS_area(:) = globalSums(maxNumBasins+1:2*maxNumBasins)
      do i = 1, maxNumBasins
         if (IS_area(i) > 0.0_RKIND) then
            mean_TF(i) = mean_TF(i) / IS_area(i) ! actual mean TF per basin
            call mpas_log_write("basin=$i, mean_TF=$r", intArgs=(/i-1/), realArgs=(/mean_TF(i)/))
         else ! avoid divide by zero for invalid indices
            mean_TF(i) = 0.0_RKIND
         endif
      enddo


      ! Now calculate melt rate
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

         call mpas_pool_get_array(geometryPool, 'ismip6shelfMelt_gamma0', gamma0)
         call mpas_pool_get_array(geometryPool, 'ismip6shelfMelt_basin', basinNumber)
         call mpas_pool_get_array(geometryPool, 'ismip6shelfMelt_TFdraft', TFdraft)
         call mpas_pool_get_array(geometryPool, 'ismip6shelfMelt_deltaT', deltaT)
         call mpas_pool_get_array(geometryPool, 'ismip6shelfMelt_offset', ismip6shelfMelt_offset)
         call mpas_pool_get_array(geometryPool, 'floatingBasalMassBal', floatingBasalMassBal)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)

         call mpas_log_write("iceshelf_melt_ismip6: gamma0=$r", realArgs=(/gamma0/))
         floatingBasalMassBal(:) = 0.0_RKIND
         coef = gamma0 * cste / scyr * rhoi
         do iCell = 1, nCellsSolve
            if ( li_mask_is_floating_ice(cellMask(iCell)) ) then
               floatingBasalMassBal(iCell) = -1.0_RKIND * coef * (TFdraft(iCell) + deltaT(iCell)) * &
                  abs(mean_TF(basinNumber(iCell)+1) + deltaT(iCell)) + ismip6shelfMelt_offset(iCell)
            endif
         enddo

         block => block % next
      enddo   ! associated(block)

      deallocate(mean_TF)
      deallocate(IS_area)

    end subroutine iceshelf_melt_ismip6
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  !  routine grounded_face_melt_ismip6
!
!> \brief Calculate grounded glacier submarine melt rate using parameterization from ISMIP6
!> \author Trevor Hillebrand
!> \date   April 2020
!> \details
!> http://www.climate-cryosphere.org/wiki/index.php?title=ISMIP6-Projections-Greenland
!> https://www.the-cryosphere.net/14/985/2020/tc-14-985-2020.pdf
!> https://agupubs-onlinelibrary-wiley-com.lanl.idm.oclc.org/doi/full/10.1002/2016GL068784
!> Melt parameterization of the form:
!> melt_rate = (A * h * q_sg^alpha + B) * TF^beta
!> where TF is local thermal forcing at each grid cell, h is water depth at ice front
!> q_sg is subglacial water glux in m/day. Default values A = 3e-4, alpha = 0.39, B = 0.15
!> beta=1.18
!-----------------------------------------------------------------------


    subroutine grounded_face_melt_ismip6(domain, applyToGrounded, &
                    applyToFloating, applyToGroundingLine, err)

      use li_calving
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: &
         domain          !< Input/Output: domain object
      logical, intent(in) :: applyToFloating, applyToGrounded, applyToGroundingLine
      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      integer :: iCell, iNeighbor, iEdge, nEmptyNeighbors
      real (kind=RKIND), pointer :: rhoi
      integer, pointer :: nCells, nISMIP6OceanLayers
      integer, dimension(:,:), pointer :: cellsOnCell, edgesOnCell
      real (kind=RKIND) :: waterDepth
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: geometryPool, meshPool, velocityPool
      real (kind=RKIND), pointer :: seaLevel
      real (kind=RKIND), pointer :: faceMeltingCFLdt
      real (kind=RKIND), pointer :: dtFaceMeltingCFLratio
      real (kind=RKIND), dimension(:), pointer :: thickness, faceMeltingThickness
      real (kind=RKIND), dimension(:), pointer :: bedTopography
      real (kind=RKIND), dimension(:), pointer :: lowerSurface
      real (kind=RKIND), dimension(:), pointer :: TFocean, ismip6Runoff
      real (kind=RKIND), dimension(:,:), pointer :: ismip6shelfMelt_3dThermalForcing
      real (kind=RKIND), dimension(:), pointer :: ismip6shelfMelt_deltaT
      real (kind=RKIND), dimension(:), pointer :: zOcean
      real (kind=RKIND), dimension(:), pointer :: areaCell
      integer, dimension(:), pointer :: cellMask, edgeMask, nEdgesOnCell
      real (kind=RKIND), pointer :: aSubglacial ! param A
      real (kind=RKIND), pointer :: alphaSubglacial ! param alpha
      real (kind=RKIND), pointer :: B ! param B
      real (kind=RKIND), pointer :: betaTF ! param beta
      real (kind=RKIND), pointer :: addTFocean ! adjust ocean thermal forcing
      real (kind=RKIND), pointer :: config_uniform_face_melt_rate
      character(len=StrKIND), pointer :: config_front_mass_bal_grounded
      logical, pointer :: config_use_3d_thermal_forcing_for_face_melt
      real (kind=RKIND), parameter :: secPerDay = 86400.0_RKIND
      real (kind=RKIND), parameter :: rhow = 1000.0_RKIND ! subglacial runoff density
      real (kind=RKIND), dimension(:), pointer :: faceMeltSpeed
      real (kind=RKIND), dimension(:), allocatable :: faceMeltSpeedVertAvg !vertically averaged faceMeltSpeed
                                                                !to pass to  li_apply_front_ablation_velocity,
                                                                !because faceMeltSpeed is only below water-line
      integer :: err_tmp, kk

      err = 0
      call mpas_log_write('Starting face melt routine')

      ! Get sea level, bedTopography, ice density
      call mpas_pool_get_config(liConfigs, 'config_ice_density', rhoi)
      call mpas_pool_get_config(liConfigs, 'config_sea_level', seaLevel)

      ! Get melt parameters
      call mpas_pool_get_config(liConfigs, 'config_beta_ocean_thermal_forcing', betaTF)
      call mpas_pool_get_config(liConfigs, 'config_alpha_subglacial_discharge', alphaSubglacial)
      call mpas_pool_get_config(liConfigs, 'config_subglacial_discharge_coefficient', aSubglacial)
      call mpas_pool_get_config(liConfigs, 'config_subglacial_discharge_intercept', B)
      call mpas_pool_get_config(liConfigs, 'config_add_ocean_thermal_forcing', addTFocean)
      call mpas_pool_get_config(liConfigs, 'config_uniform_face_melt_rate', config_uniform_face_melt_rate)
      call mpas_pool_get_config(liConfigs, 'config_front_mass_bal_grounded', config_front_mass_bal_grounded)
      call mpas_pool_get_config(liConfigs, 'config_use_3d_thermal_forcing_for_face_melt', config_use_3d_thermal_forcing_for_face_melt)

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_array(geometryPool, 'ismip6_2dThermalForcing', TFocean)
         call mpas_pool_get_array(geometryPool, 'ismip6Runoff', ismip6Runoff)

         ! Get mesh and geometry arrays
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_array(geometryPool, 'edgeMask', edgeMask)
         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
         call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
         call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(geometryPool, 'lowerSurface', lowerSurface)
         call mpas_pool_get_array(geometryPool, 'faceMeltingThickness', faceMeltingThickness)
         call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
         call mpas_pool_get_array(geometryPool, 'faceMeltSpeed', faceMeltSpeed)
         call mpas_pool_get_array(geometryPool, 'faceMeltingCFLdt', faceMeltingCFLdt)
         call mpas_pool_get_array(geometryPool, 'dtFaceMeltingCFLratio', dtFaceMeltingCFLratio)

         if ( config_use_3d_thermal_forcing_for_face_melt ) then
            call mpas_pool_get_dimension(meshPool, 'nISMIP6OceanLayers', nISMIP6OceanLayers)
            call mpas_pool_get_array(geometryPool, 'ismip6shelfMelt_3dThermalForcing', ismip6shelfMelt_3dThermalForcing)
            call mpas_pool_get_array(geometryPool, 'ismip6shelfMelt_zOcean', zOcean)
            call mpas_pool_get_array(geometryPool, 'ismip6shelfMelt_deltaT', ismip6shelfMelt_deltaT)

            ! Calculate thermal forcing based on bedTopography
            TFocean(:) = 0.0_RKIND
            do iCell = 1, nCells
               if (bedTopography(iCell) < 0.0_RKIND) then
                     kk = 1
                     ! Find deepest zOcean above ocean floor
                     do while ( (zOcean(kk) >= bedTopography(iCell)) .and. (kk < nISMIP6OceanLayers) )
                        kk = kk + 1
                     enddo
                     ! If bed is shallower than first layer, use TF from the first layer.
                     ! If bed is deeper than the bottomg ocean layer, use TF from the bottom layer.
                     if ( (kk == 1) .or. ( (kk == nISMIP6OceanLayers) .and. (zOcean(kk) >= bedTopography(iCell)) ) ) then
                        TFocean(iCell) = ismip6shelfMelt_3dThermalForcing(kk, iCell) + ismip6shelfMelt_deltaT(iCell)
                     ! For all other bed depths, interpolate linearly between layers above and below bed depth.
                     else
                        TFocean(iCell) = ( (zOcean(kk-1) - bedTopography(iCell)) * ismip6shelfMelt_3dThermalForcing(kk, iCell) + &
                                           (bedTopography(iCell) - zOcean(kk)) * ismip6shelfMelt_3dThermalForcing(kk-1, iCell) ) / &
                                           (zOcean(kk-1) - zOcean(kk)) + ismip6shelfMelt_deltaT(iCell)
                     endif
               endif
            end do
         endif

         faceMeltSpeed(:) = 0.0_RKIND
         allocate(faceMeltSpeedVertAvg(nCells+1))
         faceMeltSpeedVertAvg(:) = 0.0_RKIND

       ! Calculate face melting for each cell
         do iCell = 1, nCells

            if ( bedTopography(iCell) < 0.0_RKIND ) then
                waterDepth = seaLevel - bedTopography(iCell)
            else
                waterDepth = 0.0_RKIND
            endif

            if (trim(config_front_mass_bal_grounded) == 'ismip6') then
                ! Calculate ice front melt rate at each cell
               faceMeltSpeed(iCell) = (aSubglacial * waterDepth * & ! m s^-1
                         (ismip6Runoff(iCell) / rhow * secPerDay)**alphaSubglacial &
                         + B) *  max(0.0_RKIND, TFocean(iCell) + addTFocean)**betaTF / secPerDay
            elseif (trim(config_front_mass_bal_grounded) == 'uniform') then
                   faceMeltSpeed(iCell) = config_uniform_face_melt_rate
            endif

         enddo

         where ( (thickness > 0.0_RKIND) .and. (lowerSurface < 0.0_RKIND) )
               faceMeltSpeedVertAvg(:) = faceMeltSpeed(:) * abs(lowerSurface(:)) / thickness(:)
         end where

         call mpas_log_write("calling li_apply_front_ablation_velocity from grounded_face_melt_ismip6")
         ! Distribute melt among neighbors
         call li_apply_front_ablation_velocity(meshPool, geometryPool, velocityPool, &
                                            faceMeltingThickness, faceMeltSpeedVertAvg, applyToGrounded, &
                                            applyToFloating, applyToGroundingLine, domain, & 
                                            maxDt=faceMeltingCFLdt, CFLratio=dtFaceMeltingCFLratio, err=err_tmp)
         err = ior(err, err_tmp)


         block => block % next
      enddo    ! associated(block)

    deallocate(faceMeltSpeedVertAvg)

    end subroutine grounded_face_melt_ismip6
!-----------------------------------------------------------------------


  !***********************************************************************
  end module li_iceshelf_melt

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


