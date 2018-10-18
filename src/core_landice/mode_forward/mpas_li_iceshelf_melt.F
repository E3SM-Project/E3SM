! Copyright (c) 2013-2018,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
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

   public :: li_basal_melt_floating_ice

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

      integer, pointer :: &
           nCellsSolve                 ! number of locally owned cells

      logical, pointer :: &
           config_print_thermal_info   ! if true, print debug info

      character(len=StrKIND), pointer :: &
           config_basal_mass_bal_float ! option for basal mass balance of floating ice

      real(kind=RKIND), pointer :: &
           config_thermal_thickness, & ! minimum thickness (m) for temperature calculations
           config_sea_level,         & ! sea level (m) relative to z = 0
           config_bmlt_float_flux,   & ! constant heat flux (W/m^2) applied to the base of floating ice; positive upward
           config_bmlt_float_xlimit    ! x value (m) defining region where bmlt_float_flux is applied;
                                       !   melt only where abs(x) > xlimit

      integer, dimension(:), pointer :: &
           cellMask                    ! bit mask describing whether ice is floating, dynamically active, etc.

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

      integer :: iCell, err_tmp

      err = 0
      err_tmp = 0

      call mpas_pool_get_config(liConfigs, 'config_basal_mass_bal_float', config_basal_mass_bal_float)

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

      elseif (trim(config_basal_mass_bal_float) == 'file') then

         return   ! already set; nothing to do

      endif

      ! get rest of config variables
      call mpas_pool_get_config(liConfigs, 'config_print_thermal_info', config_print_thermal_info)
      call mpas_pool_get_config(liConfigs, 'config_thermal_thickness', config_thermal_thickness)
      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
      call mpas_pool_get_config(liConfigs, 'config_bmlt_float_flux', config_bmlt_float_flux)
      call mpas_pool_get_config(liConfigs, 'config_bmlt_float_xlimit', config_bmlt_float_xlimit)

      if (config_print_thermal_info) then
         call mpas_log_write('Solving for basal melting of floating ice, config_basal_mass_bal_float = ' // &
              trim(config_basal_mass_bal_float) )
      endif

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

         ! Compute basal melting for floating ice.

         call basal_melt_floating_ice(&
              config_basal_mass_bal_float,  &
              nCellsSolve,                  &
              xCell,                        &
              daysSinceStart,               &
              thermalCellMask,              &
              li_mask_is_floating_ice_int(cellMask), &
              lowerSurface,                 &
              bedTopography,                &
              config_sea_level,             &
              config_bmlt_float_flux,       &
              config_bmlt_float_xlimit,     &
              floatingBasalMassBal,         &
              err_tmp)
         err = ior(err, err_tmp)

         !WHL - debug
         ! write(stdoutUnit,*) 'Computed basal melt for floating ice'
         ! write(stdoutUnit,*) ' '
         ! write(stdoutUnit,*) 'iCell, thickness, basal mbal:'
         ! do iCell = 1, nCellsSolve
         !    if (li_mask_is_floating_ice(cellMask(iCell))) then
         !       write(stdoutUnit,*) iCell, thickness(iCell), floatingBasalMassBal(iCell)*scyr/config_ice_density
         !    endif
         ! enddo
         ! write(stdoutUnit,*) 'Done with basal melt for floating ice'

         ! clean up
         call mpas_deallocate_scratch_field(thermalCellMaskField, .true.)

         block => block % next
      enddo   ! associated(block)


    end subroutine li_basal_melt_floating_ice


!***********************************************************************
!***********************************************************************
! Private subroutines:
!***********************************************************************
!***********************************************************************



!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  !  routine basal_melt_floating_ice
!
!> \brief MPAS melt rate at base of floating ice
!> \author William Lipscomb
!> \date   November 2015
!> \details
!>  This routine computes the melt rate at the base of floating ice.
!-----------------------------------------------------------------------

    subroutine basal_melt_floating_ice(&
         config_basal_mass_bal_float,   &
         nCellsSolve,                   &
         xCell,                         &
         daysSinceStart,                &
         iceMask,                       &
         floatingMask,                  &
         lowerSurface,                  &
         bedTopography,                 &
         config_sea_level,              &
         config_bmlt_float_flux,        &
         config_bmlt_float_xlimit,      &
         floatingBasalMassBal,          &
         err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      character(len=strKIND), intent(in) :: &
           config_basal_mass_bal_float  !< Input: option for computing basal mass balance for floating ice

      integer, intent(in) :: &
           nCellsSolve              !< Input: number of locally owned cells

      real(kind=RKIND), dimension(:), intent(in) :: &
           xCell                    !< Input: x coordinate for each cell

      real(kind=RKIND), pointer, intent(in) :: daysSinceStart

      real(kind=RKIND), dimension(:), intent(in) :: &
           lowerSurface,          & !< Input: lower surface elevation (m)
           bedTopography            !< Input: elevation of bed topography (m)

      real(kind=RKIND), intent(in) :: &
           config_sea_level         !< Input: sea level relative to z = 0 (m)

      integer, dimension(:), intent(in) :: &
           iceMask,               & !< Input: = 1 where ice exists (thickness > config_thermal_thickness), else = 0
           floatingMask             !< Input: = 1 where ice is floating, else = 0

      ! inputs for constant melt rate as in MISMIP+ Ice2 experiments
      real(kind=RKIND), intent(in) ::  &
           config_bmlt_float_flux,   & !< Input: constant heat flux (W/m^2) applied to the base of floating ice; positive upward
                                       !<        MISMIP+ default value = 975.17 W/m^2 (gives melt rate of 100 m/yr)
           config_bmlt_float_xlimit    !< Input: x value (m) defining region where bmlt_float_flux is applied;
                                       !<        melt only where abs(x) > xlimit
                                       !<        MISMIP+ default value = 480 km

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real(kind=RKIND), dimension(:), intent(out):: &
           floatingBasalMassBal     !< Output: basal mass balance for floating ice

      integer, intent(out) :: err

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      real(kind=RKIND) :: &
           bmlt_float_rate                          ! constant basal melt rate (m/s)
                                                    ! = config_bmlt_float_flux / (config_ice_density*latent_heat_ice)
      integer :: iCell

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

      ! basal melt parameters for Seroussi param.
      real (kind=RKIND) :: slopeSer  ! slope of relation between depth and melt rate
      real (kind=RKIND) :: interceptSer  ! depth at which melting goes to 0
      real (kind=RKIND) :: maxMeltSer  ! maximum allowable melt rate
      real (kind=RKIND) :: sillDepth  ! depth below which melt rate no longer increases
      real (kind=RKIND), pointer :: config_basal_mass_bal_seroussi_amplitude
      real (kind=RKIND), pointer :: config_basal_mass_bal_seroussi_period
      real (kind=RKIND), pointer :: config_basal_mass_bal_seroussi_phase

      real (kind=RKIND), pointer :: config_ice_density !< ice density


      ! Compute melt rate for floating ice
      err = 0

      call mpas_pool_get_config(liConfigs, 'config_ice_density', config_ice_density)
      ! initialize to zero melt
      floatingBasalMassBal(:) = 0.0_RKIND
      if (trim(config_basal_mass_bal_float) == 'none') then
         ! Do nothing, handled in calling routine

      elseif (trim(config_basal_mass_bal_float) == 'file') then
         ! Do nothing, handled in calling routine

      elseif (trim(config_basal_mass_bal_float) == 'constant') then

         ! set melt rate to a constant value for floating ice
         ! allow basal melt in ice-free ocean cells, in case ice is advected to those cells by the transport scheme

         bmlt_float_rate = config_bmlt_float_flux / (config_ice_density*latent_heat_ice)  ! convert W/m^2 to m/s

         floatingBasalMassBal(:) = 0.0_RKIND

         do iCell = 1, nCellsSolve
            if ( floatingMask(iCell) == 1 .or.  &
                 (bedTopography(iCell) < config_sea_level .and. iceMask(iCell) == 0) ) then
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

            if ( floatingMask(iCell) == 1 .or.  &
                 (bedTopography(iCell) < config_sea_level .and. iceMask(iCell) == 0) ) then
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

         ! Melt rate parameterization from:
         ! Seroussi, H., Y. Nakayama, E. Larour, D. Menemenlis, M. Morlighem, E. Rignot, and A. Khazendar (2017),
         ! Continued retreat of Thwaites Glacier, West Antarctica, controlled by bed topography and ocean circulation,
         ! Geophys. Res. Lett., 1-9, doi:10.1002/2017GL072910.
         ! for Thwaites Glacier.
         ! Specifically, this is a linear fit of melt with shelf draft from the Supplemental Information, Figure S1.
         ! The linear relation is modified by a:
         !   * depth above which there is no melt (Antarctic Surface Water saturation)
         !   * a maximum melt rate (Circumpolar Deep Water saturation)
         !   * a depth below which melt stops increasing (minimum sill height)

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
            if ( (floatingMask(iCell) == 1) .and. (lowerSurface(iCell) < -10.0_RKIND) ) then
               ! ice is present and floating

               zDraft = lowerSurface(iCell) - config_sea_level
               ! Coefficients for m/yr melt rate (in units of Seroussi figure but without negative meaning melting)
               floatingBasalMassBal(iCell) = max(-1.0_RKIND * maxMeltSer, min(0.0_RKIND, slopeSer * &
                  (max(zDraft, sillDepth) - interceptSer)))

            endif   ! ice is present
         enddo   ! iCell

         ! change units from m/yr to kg/m2/s
         floatingBasalMassBal(:) = floatingBasalMassBal(:) * config_ice_density / scyr

      else

         call mpas_log_write('Unknown option selected for config_basal_mass_bal_float:' // &
            trim(config_basal_mass_bal_float), MPAS_LOG_ERR)
         err = ior(err, 1)

      endif   ! config_basal_mass_bal_float


    end subroutine basal_melt_floating_ice



    !***********************************************************************

  end module li_iceshelf_melt

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


