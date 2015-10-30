! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_thermal
!
!> \brief MPAS land ice vertical temperature/enthalpy solver
!> \author William Lipscomb
!> \date   October 2015
!> \details
!>  This module contains solvers for the vertical temperature
!>  and/or enthalpy profile.
!
!-----------------------------------------------------------------------

module li_thermal

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar
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

   public :: li_thermal_init, li_thermal_solver, li_heat_dissipation_sia

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   real (kind=RKIND), save :: rhoi    ! ice density (kg m^{-3}), copied from config_ice_density
   real (kind=RKIND), save :: rhoo    ! ocean density (kg m^{-3}), copied from config_ocean_density

   !TODO - dups is from CISM.  Choose a better name?
   real (kind=RKIND), dimension(:,:), allocatable :: dups   ! vertical grid quantities

   ! max and min allowed temperatures (Kelvin)
   ! Note: kelvin_to_celsius = 273.15 (perhaps it should be called celsius_to_kelvin?)

   real (kind=RKIND), parameter ::   &
        maxtemp_threshold =  100._RKIND + kelvin_to_celsius, &
        mintemp_threshold = -100._RKIND + kelvin_to_celsius

!***********************************************************************
   contains
!***********************************************************************


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  !  routine li_thermal_init
!
!> \brief MPAS land ice initialize vertical temperature
!> \author William Lipscomb
!> \date   October 2015
!> \details
!>  This routine initializes the vertical temperature profile in each column
!>  and computes some quantities required by the thermal solver.
!-----------------------------------------------------------------------

   subroutine li_thermal_init(domain, err)

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
      type (mpas_pool_type), pointer :: thermalPool 

      ! config options

      logical, pointer :: &
           config_print_thermal_info,  &
           config_do_restart

      character(len=StrKIND), pointer :: &
           config_thermal_solver,    & ! solver option ('temperature' or 'enthalpy')
           config_thermal_init         ! initialization option ('linear' or 'file')

      real (kind=RKIND), pointer ::  &
           config_ice_density,       & ! ice density
           config_ocean_density        ! ocean density

      integer, pointer :: &
           index_temperature,  &
           index_waterfrac

      integer, pointer :: &
           nCellsSolve,              & ! number of locally owned cells
           nVertLevels                 ! number of vertical layers

      real (kind=RKIND), dimension(:), pointer :: &
           layerCenterSigma,         & ! sigma coordinate at midpoint of each layer
           layerInterfaceSigma         ! sigma coordinate at layer interfaces (including top and bottom)

      real (kind=RKIND), dimension(:), pointer :: &
           thickness                   ! ice thickness

      real (kind=RKIND), dimension(:,:,:), pointer :: tracers

      real (kind=RKIND), dimension(:,:), pointer :: &
           temperature,              & ! interior ice temperature (K)
           waterfrac,                & ! interior water fraction (unitless)
           enthalpy                    ! interior ice enthalpy (J m^{-3})
      
      real (kind=RKIND), dimension(:), pointer :: &
           surfaceTemperature,       & ! surface ice temperature (K)
           basalTemperature,         & ! basal ice temperature (K)
           surfaceAirTemperature       ! surface air temperature (K)

      ! Note: The following fields are needed for halo updates
      ! TODO - Are halo updates needed at initialization?
      type (field1DReal), pointer :: &
           surfaceTemperatureField,  &
           basalTemperatureField

      type (field2DReal), pointer :: &
           temperatureField,  &
           waterfracField

      real (kind=RKIND), dimension(:), allocatable :: &
           pmptemp                     ! pressure melting point temp in ice interior

      real (kind=RKIND) :: &
           pmptemp_bed                 ! pressure melting point temp at bed

      real (kind=RKIND), parameter :: &
           pmpt_offset = 2.0_RKIND     ! offset of initial Tbed from pressure melting point temperature (K)
                                       ! Note: pmtp_offset is positive for T < Tpmp

      integer :: k, iLayer

      integer :: iCell

      integer :: err_tmp

      !WHL - debug - for circular shelf test case
      integer, parameter :: ncellsPerRow = 40
      integer, parameter :: nRows = 46
      integer :: i, iRow

      err = 0

      ! get config options
      call mpas_pool_get_config(liConfigs, 'config_thermal_solver', config_thermal_solver)
      call mpas_pool_get_config(liConfigs, 'config_thermal_init', config_thermal_init)
      call mpas_pool_get_config(liConfigs, 'config_print_thermal_info', config_print_thermal_info)
      call mpas_pool_get_config(liConfigs, 'config_do_restart', config_do_restart)
            
      ! set some physical constants
      ! (to avoid calling mpas_pool_get_config repeatedly in this module)
      call mpas_pool_get_config(liConfigs, 'config_ice_density', config_ice_density)
      call mpas_pool_get_config(liConfigs, 'config_ocean_density', config_ocean_density)
      rhoi = config_ice_density
      rhoo = config_ocean_density

      ! block loop
      block => domain % blocklist
      do while (associated(block))

         write(stderrUnit,*) 'Get pools'

         ! get pools
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'thermal', thermalPool)

         ! get dimensions
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(thermalPool, 'index_temperature', index_temperature)
         call mpas_pool_get_dimension(thermalPool, 'index_waterfrac', index_waterfrac)

         ! get fields from the mesh pool
         call mpas_pool_get_array(meshPool, 'layerCenterSigma', layerCenterSigma)
         call mpas_pool_get_array(meshPool, 'layerInterfaceSigma', layerInterfaceSigma)

         ! get fields from the geometry pool
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)

         ! get fields from the thermal pool
         call mpas_pool_get_array(thermalPool, 'tracers', tracers)
!         call mpas_pool_get_array(thermalPool, 'temperature', temperature)
!         call mpas_pool_get_array(thermalPool, 'waterfrac',   waterfrac)
         call mpas_pool_get_array(thermalPool, 'surfaceAirTemperature', surfaceAirTemperature)
         call mpas_pool_get_array(thermalPool, 'surfaceTemperature', surfaceTemperature)
         call mpas_pool_get_array(thermalPool, 'basalTemperature',   basalTemperature)

         ! temporary init for surfaceAirTemperature
         !TODO - Read surfaceAirTemperature from a file or create a simple field
         surfaceAirTemperature(:) = kelvin_to_celsius   ! 273.15

         !TODO - Is there a better way to access the temperature and waterfrac arrays?
         temperature => tracers(index_temperature,:,:)
         waterfrac => tracers(index_waterfrac,:,:)

         if (config_print_thermal_info) then

            write(stderrUnit,*) 'Initialize thermal solver, config_thermal_init =', trim(config_thermal_init)

            !WHL - temporary debugging code - for circular shelf test case
!            write(stderrUnit,*) ' '
!            write(stderrUnit,*) 'Surface ice temperature before init'
!            do iRow = nRows, 1, -1
!               if (mod(iRow,2) == 0) then  ! indent for even-numbered rows
!                  write(stderrUnit,'(a3)',advance='no') '    '
!               endif
!               do i = nCellsPerRow/2 - 2, nCellsPerRow
!                  iCell = (iRow-1)*nCellsPerRow + i
!                  write(stderrUnit,'(f8.2)',advance='no') surfaceTemperature(iCell)
!               enddo
!               write(stderrUnit,*) ' '
!            enddo

         endif  ! config_print_thermal_info

         ! Precompute some grid quantities used in the vertical temperature solve
         ! (Commented-out lines are from CISM)

         allocate(dups(nVertLevels,2))
         dups(:,:) = 0.0_RKIND
      
         k = 1
         !      dups(k,1) = 1.0_RKIND/((sigma(k+1) - sigma(k)) * (stagsigma(k) - sigma(k)) )
         dups(k,1) = 1.0_RKIND/((layerInterfaceSigma(k+1) - layerInterfaceSigma(k)) * (layerCenterSigma(k) - layerInterfaceSigma(k)) )
      
         do k = 2, nVertLevels
            !         dups(k,1) = 1.0_RKIND/((sigma(k+1) - sigma(k)) * (sigma(k) - sigma(k-1)) )
            dups(k,1) = 1.0_RKIND/((layerInterfaceSigma(k+1) - layerInterfaceSigma(k)) * (layerCenterSigma(k) - layerCenterSigma(k-1)) )
         enddo
         
         do k = 1, nVertLevels-1
            !         dups(k,2) = 1.0_RKIND/((sigma(k+1) - sigma(k)) * (stagsigma(k+1) - stagsigma(k)) )
            dups(k,2) = 1.0_RKIND/((layerInterfaceSigma(k+1) - layerInterfaceSigma(k)) * (layerCenterSigma(k+1) - layerCenterSigma(k)) )
         end do
         
         k = nVertLevels
         !      dups(k,2) = 1.0_RKIND/((sigma(k+1) - sigma(k)) * (sigma(k+1) - stagsigma(k)) )
         dups(k,2) = 1.0_RKIND/((layerInterfaceSigma(k+1) - layerInterfaceSigma(k)) * (layerInterfaceSigma(k+1) - layerCenterSigma(k)) )
         
         if (config_print_thermal_info) then
            write(stderrUnit,*) 'dups coefficients:'
            do k = 1, nVertLevels
               write(stderrUnit,*) k, dups(k,1), dups(k,2)
            enddo
         endif
         
         ! Initialize vertical temperature profile.
         ! Three possibilities:
         ! (1) Set up a linear temperature profile, with T = artm at the surface and T <= Tpmp
         !     at the bed (config_thermal_init = 'linear'). 
         !     A parameter (pmpt_offset) controls how far below Tpmp the initial bed temp is set.
         ! (2) Read ice temperature from an initial input file (config_thermal_init = 'file').
         ! (3) Read ice temperature from a restart file.
         !
         ! The default is (1).
         ! If restarting, we always do (3).
         ! If (2) or (3), then the temperature should already have been read in, and there is
         !  nothing to do here (except possibly to set waterfrac).

         if (config_do_restart) then

            ! nothing to do; temperature was read from the restart file
            !TODO - Make sure waterfrac is also read, if needed
            if (config_print_thermal_info) then
               write(stderrUnit,*) 'Initialized ice temperature from the restart file'
            endif

         elseif (trim(config_thermal_init) == 'file') then

            ! Temperature was read from the input file
            if (config_print_thermal_info) then
               write(stderrUnit,*) 'Initialized ice temperature from the input file'
            endif

            ! If using the enthalpy solver, initialize waterfrac here
            !TODO - Allow waterfrac to be read from the input file?

            if (trim(config_thermal_solver) == 'enthalpy') then
               waterfrac(:,:) = 0.0_RKIND
            endif

         elseif (trim(config_thermal_init) == 'linear') then

            ! set up a linear temperature profile in each column
            ! T = surfaceAirTemperature at the ice surface, and T <= Tpmp at the bed

            allocate(pmptemp(nVertLevels))

            ! initialize T = 273.15 K = 0 C everywhere

            temperature(:,:) = kelvin_to_celsius    ! = 273.15                                              

            do iCell = 1, nCellsSolve

               ! set surface temperature to the air temperature (or 273.15, whichever is less)
               surfaceTemperature(iCell) = min(surfaceAirTemperature(iCell), kelvin_to_celsius)

               ! compute the pressure melting point temperature in the column and at the bed
               !TODO - Change pmp subroutine to return temperature in Kelvin

               call pressure_melting_point_column(layerCenterSigma(:), thickness(iCell), pmptemp(:))
               call pressure_melting_point(thickness(iCell), pmptemp_bed)

               pmptemp(:) = pmptemp(:) + kelvin_to_celsius
               pmptemp_bed = pmptemp_bed + kelvin_to_celsius

               ! set the basal temperature to slightly below the pressure melting point temperature
               basalTemperature(iCell) = pmptemp_bed - pmpt_offset

               ! set the interior temperatures
               ! make sure T <= Tpmp - pmpt_offset in column interior

               temperature(:,iCell) = surfaceTemperature(iCell) +  &
                                     (basalTemperature(iCell) - surfaceTemperature(iCell)) * layerCenterSigma(:)
               
               temperature(:,iCell) = min(temperature(:,iCell), pmptemp(:) - pmpt_offset)

            enddo  ! iCell

            ! initialize the water fraction to zero
            waterfrac(:,:) = 0.0_RKIND

            if (config_print_thermal_info) then
               write(stderrUnit,*) 'Initialized a linear temperature profile in each column'
            endif

         endif    ! restart file, input file, or linear

         ! clean up
         if (allocated(pmptemp)) deallocate(pmptemp)

         block => block % next
      enddo

      !TODO - Add a debug check for bad values
      !       E.g., make sure the temperature read from a file is in Kelvin and not Celsius

      ! halo updates

      call mpas_pool_get_field(thermalPool, 'surfaceTemperature', surfaceTemperatureField)
      call mpas_dmpar_exch_halo_field(surfaceTemperatureField)

      call mpas_pool_get_field(thermalPool, 'basalTemperature', basalTemperatureField)
      call mpas_dmpar_exch_halo_field(basalTemperatureField)

      !TODO - Halo updates for components of the tracer array
!!      call mpas_pool_get_field(thermalPool, 'temperature', temperatureField)
!!      call mpas_dmpar_exch_halo_field(temperatureField)

      if (trim(config_thermal_solver) == 'enthalpy') then
!!         call mpas_pool_get_field(thermalPool, 'waterfrac', waterfracField)
!!         call mpas_dmpar_exch_halo_field(waterfracField)
      endif

      ! === error check
      if (err > 0) then
          write (stderrUnit,*) "An error has occurred in li_thermal_init."
      endif

   !--------------------------------------------------------------------
    end subroutine li_thermal_init

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  !  routine li_thermal_solver
!
!> \brief MPAS land ice solver for vertical temperature/enthalpy
!> \author William Lipscomb
!> \date   October 2015
!> \details
!>  This routine is the driver for the vertical temperature/enthalpy
!>  calculation in each ice column. The following options are supported:
!>  (1) Do nothing (config_thermal_solver = 'none')
!>  (2) Standard prognostic temperature solve (config_thermal_solver = 'temperature')
!>  (3) Prognostic solve for enthalpy (config_thermal_solver = 'enthalpy') 

!-----------------------------------------------------------------------

   subroutine li_thermal_solver(domain, deltat, err)

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

      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: meshPool 
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: thermalPool 
      type (mpas_pool_type), pointer :: velocityPool   ! needed for mask subroutine

      integer, pointer :: &
           nCellsSolve,              & ! number of locally owned cells
           nVertLevels,              & ! number of vertical layers
           nVertInterfaces             ! number of vertical interfaces (including top and bottom)

      logical, pointer :: &
           config_print_thermal_info   ! if true, print debug info

      character(len=StrKIND), pointer :: &
           config_thermal_solver       ! option for thermal solver

      real(kind=RKIND), pointer :: &
           config_thermal_thickness    ! minimum thickness for temperature calculations

      integer, pointer :: &
           config_stats_cell_ID        ! global ID for diagnostic cell

      integer, dimension(:), pointer :: &
           cellMask,                 & ! bit mask describing whether ice is floating, dynamically active, etc.
           thermalCellMask,          & ! mask for thermal calculations
                                       ! = 1 where thickness > config_thermal_thickness, elsewhere = 0 
           indexToCellID               ! list of global cell IDs
      
      real (kind=RKIND), dimension(:), pointer :: &
           layerCenterSigma,         & ! sigma coordinate at midpoint of each layer
           layerInterfaceSigma         ! sigma coordinate at layer interfaces (including top and bottom)

      real (kind=RKIND), dimension(:), pointer :: &
           surfaceTemperature,       & ! surface ice temperature (K)
           basalTemperature,         & ! basal ice temperature (K)
           surfaceAirTemperature,    & ! surface air temperature (K)
           basalHeatFlux,            & ! basal heat flux into the ice (W m^{-2}, positive upward)
           basalFrictionFlux,        & ! basal frictional flux into the ice (W m^{-2})
           surfaceConductiveFlux,    & ! conductive heat flux at the upper surface (W m^{-2}, positive downward) 
           basalConductiveFlux,      & ! conductive heat flux at the lower surface (W m^{-2}, positive downward) 
           basalMassBal,             & ! basal mass balance (kg m^{-2} s^{-1}); positive for freeze-on, negative for melting
           basalWaterThickness,      & ! basal water thickness (m)
           thickness,                & ! ice thickness (m)
           bedTopography               ! bed topography (m; negative below sea level)

      real (kind=RKIND), dimension(:,:), pointer :: &
           temperature,              & ! interior ice temperature (K)
           waterfrac,                & ! interior water fraction (unitless)
           enthalpy,                 & ! interior ice enthalpy (J m^{-3})
           heatDissipation             ! interior heat dissipation (deg/s)   ! TODO - These are CISM units.  Change?

      ! Note: The following fields are needed for halo updates

      type (field1DReal), pointer :: &
           surfaceTemperatureField,  &
           basalTemperatureField

      type (field2DReal), pointer :: &
           temperatureField,  &
           waterfracField

      real(kind=RKIND), dimension(:), allocatable :: &
           subdiagonal, diagonal, superdiagonal,   &  ! tridiagonal matrix elements
           rhs                                        ! matrix right-hand side

      !TODO - Change name of alpha_enth?
      real(kind=RKIND), dimension(:), allocatable :: &
           alpha_enth   ! diffusivity at interfaces (m2/s) for enthalpy solver
                        ! = coni / (rhoi*cp_ice) for cold ice  !TODO - Change name of coni?

      !TODO - Get rid of these temporary arrays?
      real(kind=RKIND), dimension(:), allocatable :: &
           temp, enth       ! like temperature/enthalpy, but including surface and bed

      real(kind=RKIND) :: &
           depth,                    & ! depth within ice column
           dTtop, dTbot,             & ! temperature differences
           denth_top, denth_bot,     & ! enthalpy differences
           columnHeatDissipation,    & ! integrated heat dissipation in column
           maxtemp, mintemp,         & ! max and min temperatures in column
           initialEnergy,            & ! initial energy in ice column (J m^{-2})
           finalEnergy,              & ! final energy in ice column (J m^{-2})
           deltaEnergy                 ! change in energy

      integer :: iCell, err_tmp

      logical :: verboseColumn

      integer :: k

      !WHL - debug - for circular shelf test case
      integer, parameter :: ncellsPerRow = 40
      integer, parameter :: nRows = 46
      integer :: i, iRow

      err = 0

      ! block loop
      block => domain % blocklist
      do while (associated(block))

         ! get pools
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'thermal', thermalPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)

         ! get dimensions
         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(meshPool, 'nVertInterfaces', nVertInterfaces)

         ! get fields from the mesh pool
         call mpas_pool_get_array(meshPool, 'layerCenterSigma', layerCenterSigma)
         call mpas_pool_get_array(meshPool, 'layerInterfaceSigma', layerInterfaceSigma)
         call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)  ! diagnostic only

         ! get fields from the geometry pool
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_array(geometryPool, 'basalMassBal', basalMassBal)
         call mpas_pool_get_array(geometryPool, 'basalWaterThickness', basalWaterThickness)

         ! get fields from the thermal pool
         call mpas_pool_get_array(thermalPool, 'surfaceTemperature', surfaceTemperature)
         call mpas_pool_get_array(thermalPool, 'basalTemperature',   basalTemperature)
         call mpas_pool_get_array(thermalPool, 'temperature', temperature)
         call mpas_pool_get_array(thermalPool, 'waterfrac',   waterfrac)
         call mpas_pool_get_array(thermalPool, 'enthalpy',    enthalpy)
         call mpas_pool_get_array(thermalPool, 'surfaceAirTemperature', surfaceAirTemperature)
         call mpas_pool_get_array(thermalPool, 'surfaceConductiveFlux', surfaceConductiveFlux)
         call mpas_pool_get_array(thermalPool, 'basalConductiveFlux', basalConductiveFlux)
         call mpas_pool_get_array(thermalPool, 'basalHeatFlux', basalHeatFlux)
         call mpas_pool_get_array(thermalPool, 'basalFrictionFlux', basalFrictionFlux)
         call mpas_pool_get_array(thermalPool, 'heatDissipation', heatDissipation)

         ! get fields from the scratch pool
         call mpas_pool_get_array(thermalPool, 'iceCellMask', thermalCellMask)

         ! get config variables
         call mpas_pool_get_config(liConfigs, 'config_print_thermal_info', config_print_thermal_info)
         call mpas_pool_get_config(liConfigs, 'config_thermal_solver', config_thermal_solver)
         call mpas_pool_get_config(liConfigs, 'config_thermal_thickness', config_thermal_thickness)
         call mpas_pool_get_config(liConfigs, 'config_stats_cell_ID', config_stats_cell_ID)

         if (config_print_thermal_info) then
            write(stderrUnit,*) 'Solving for temperature, config_thermal_solver = ', config_thermal_solver

            !WHL - debug - for circular shelf test case
!            write(stderrUnit,*) 'Surface ice temperature before thermal calc'
!            do iRow = nRows, 1, -1
!               if (mod(iRow,2) == 0) then  ! indent for even-numbered rows
!                  write(stderrUnit,'(a3)',advance='no') '    '
!               endif
!               do i = nCellsPerRow/2 - 2, nCellsPerRow
!                  iCell = (iRow-1)*nCellsPerRow + i
!                  write(stderrUnit,'(f8.2)',advance='no') surfaceTemperature(iCell)
!               enddo
!               write(stderrUnit,*) ' '
!            enddo

         endif

         ! calculate masks - so we know where the ice is floating
         call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
         err = ior(err, err_tmp)

         select case(config_thermal_solver)

         case ('none')

            ! Do nothing

         case ('temperature', 'enthalpy')

            ! Convert temperature from Kelvin to Celsius to avoid repeated use of kelvin_to_celsius below
            ! (Convert back at the end.)
            temperature(:,:) = temperature(:,:) - kelvin_to_celsius

            ! allocate some vertical arrays
            allocate(subdiagonal(nVertInterfaces+1))  ! temperature/enthalpy in each layer, plus surface and basal temperature
            allocate(diagonal(nVertInterfaces+1))    
            allocate(superdiagonal(nVertInterfaces+1))    
            allocate(rhs(nVertInterfaces+1))
            allocate(alpha_enth(nVertInterfaces))

            !TODO - Get rid of these temporary arrasy?
            allocate(temp(0:nVertInterfaces))
            allocate(enth(0:nVertInterfaces))

            ! loop over locally owned cells
            do iCell = 1, nCellsSolve
               
               if (config_print_thermal_info .and. indexToCellID(iCell) == config_stats_cell_ID) then
                  verboseColumn = .true.
               else
                  verboseColumn = .false.
               endif

               if (thickness(iCell) > config_thermal_thickness) then

                  ! set thermal mask
                  thermalCellMask(iCell) = 1

                  ! Set surface temperature (Celsius)
                  
                  surfaceTemperature(iCell) = min(0.0_RKIND, surfaceAirTemperature(iCell) - kelvin_to_celsius)

                  ! For floating ice, set the basal temperature to the freezing temperature of seawater
                  ! Values based on Ocean Water Freezing Point Calculator with S = 35 PSU 
                  if (li_mask_is_floating_ice(cellMask(iCell))) then
                     depth = thickness(iCell) * rhoi/rhoo
                     basalTemperature(iCell) = oceanFreezingTempSurface + oceanFreezingTempDepthDependence * depth  ! Celsius
                  endif
                  
                  if (trim(config_thermal_solver) == 'enthalpy') then

                     ! Given temperature and waterfrac in ice interior, compute enthalpy

                     call temperature_to_enthalpy(&
                          layerCenterSigma,      &
                          thickness(iCell),      &
                          temperature(:,iCell),  &
                          waterfrac(:,iCell),    &
                          enthalpy(:,iCell))

                     if (verboseColumn) then
                        write(stderrUnit,*) ' '
                        write(stderrUnit,*) 'Before prognostic enthalpy, iCell =', indexToCellID(iCell)
                        write(stderrUnit,*) 'thickness =', thickness(iCell)
                        write(stderrUnit,*) 'Temperature (C), waterfrac, enthalpy/(rhoi*cp_ice):'
                        write(stderrUnit,*) surfaceTemperature(iCell)
                        do k = 1, nVertLevels
                           write(stderrUnit,*) k, temperature(k,iCell), waterfrac(k,iCell), enthalpy(k,iCell)/(rhoi*cp_ice)
                        enddo
                        write(stderrUnit,*) basalTemperature(iCell)
                     endif

                     ! compute initial internal energy in column (for energy conservation check)
                     initialEnergy = 0.0_RKIND
                     do k = 1, nVertLevels
                        initialEnergy = initialEnergy + enthalpy(k,iCell) * (layerInterfaceSigma(k+1) - layerInterfaceSigma(k)) * thickness(iCell)
                     enddo

                     ! compute matrix elements using enthalpy gradient method

                     temp(0) = surfaceTemperature(iCell)
                     temp(1:nVertLevels) = temperature(:,iCell)
                     temp(nVertInterfaces) = basalTemperature(iCell)

                     enth(0) = surfaceTemperature(iCell) * rhoi*cp_ice
                     enth(1:nVertLevels) = enthalpy(:,iCell)
                     enth(nVertInterfaces) = basalTemperature(iCell) * rhoi*cp_ice

                     call enthalpy_matrix_elements(&
                          deltat,               &
                          nVertInterfaces,      &   ! CISM passes upn
                          layerCenterSigma,     &   ! CISM passes stagsigma
                          subdiagonal,          &
                          diagonal,             &
                          superdiagonal,        &
                          rhs,                  &
                          dups,                 &
                          li_mask_is_floating_ice_int(cellMask(iCell)), &
                          thickness(iCell),               &
                          temp(0:nVertInterfaces),        &
                          waterfrac(:,iCell),             & 
                          enth(0:nVertInterfaces),        & 
                          heatDissipation(:,iCell),       &
!!                          basalHeatFlux(iCell),           &
                          -basalHeatFlux(iCell),           &    ! CISM subroutine assumes positive down, so flip sign
                          basalFrictionFlux(iCell),       &
                          alpha_enth,                &
                          verboseColumn)

                     if (verboseColumn) then
                        write(stderrUnit,*) ' '
                        write(stderrUnit,*) 'After matrix elements, iCell =', indexToCellID(iCell)
                        write(stderrUnit,*) 'k, subd, diag, supd, rhs/(rhoi*ci):'
                        do k = 1, nVertInterfaces+1
                           write(stderrUnit,*) k-1, subdiagonal(k), diagonal(k), superdiagonal(k), rhs(k)/(rhoi*cp_ice)
                        enddo
                     endif

                     ! solve the tridiagonal system
                     ! Note: Enthalpy is indexed from 0 to nVertInterfaces, with indices 1 to nVertInterfaces-1 colocated
                     ! with layerCenterSigma values of the same index.
                     ! However, the matrix elements are indexed 1 to nVertInterfaces+1, with the first row
                     ! corresponding to the surface enthalpy, enthalpy(0).

                     call tridiag_solver(&
                          subdiagonal,      &
                          diagonal,         &
                          superdiagonal,    &
                          enth(0:nVertInterfaces),  &
                          rhs)

                     ! Compute conductive fluxes = (alpha/H * denth/dsigma) at upper and lower surfaces; positive down.
                     ! Here alpha = coni / (rhoi*shci) for cold ice, with a smaller value for temperate ice.
                     ! Assume implicit backward Euler time step.
                     ! Note: These fluxes should be computed before calling glissade_enth2temp (which might reset the bed enthalpy).

                     ! commented-out code is from CISM
!!                     denth_top = enthalpy(1,ew,ns) - enthalpy(0,ew,ns)
!!                     denth_bot = enthalpy(upn,ew,ns) - enthalpy(upn-1,ew,ns)
                
!!                     ucondflx(ew,ns) = -alpha_enth(1)  /thck(ew,ns) * denth_top/( stagsigma(1))
!!                     lcondflx(ew,ns) = -alpha_enth(upn)/thck(ew,ns) * denth_bot/(1.d0 - stagsigma(upn-1))

                     denth_top = enth(1) - enth(0)
                     denth_bot = enth(nVertInterfaces) - enth(nVertLevels)
                
                     surfaceConductiveFlux(iCell) = -alpha_enth(1)/thickness(iCell) * denth_top/layerCenterSigma(1)
                     basalConductiveFlux(iCell) = -alpha_enth(nVertInterfaces)/thickness(iCell) * denth_bot/(1.0_RKIND - layerCenterSigma(nVertLevels))
                                              
                     ! copy the enthalpy back into the full array
                     enthalpy(:,iCell) = enth(1:nVertLevels)

                     ! convert enthalpy back to temperature and waterfrac

                     call enthalpy_to_temperature(&
                          layerCenterSigma,      &
                          thickness(iCell),      &
                          enthalpy(:,iCell),     &
                          temperature(:,iCell),  &
                          waterfrac(:,iCell))

                     if (verboseColumn) then
                        write(stderrUnit,*) ' '
                        write(stderrUnit,*) 'After prognostic enthalpy, iCell =', indexToCellID(iCell)
                        write(stderrUnit,*) 'thickness =', thickness(iCell)
                        write(stderrUnit,*) 'Temp, waterfrac, enthalpy/(rhoi*cp_ice):'
                        write(stderrUnit,*) surfaceTemperature(iCell)
                        do k = 1, nVertLevels
                           write(stderrUnit,*) k, temperature(k,iCell), waterfrac(k,iCell), enthalpy(k,iCell)/(rhoi*cp_ice)
                        enddo
                        write(stderrUnit,*) basalTemperature(iCell)
                     endif

                     ! compute final internal energy in column (for energy conservation check)
                     finalEnergy = 0.0_RKIND
                     do k = 1, nVertLevels
                        finalEnergy = finalEnergy + enthalpy(k,iCell) * (layerInterfaceSigma(k+1) - layerInterfaceSigma(k)) * thickness(iCell)
                     enddo

                  else    ! temperature solver

                     if (verboseColumn) then
                        write(stderrUnit,*) ' '
                        write(stderrUnit,*) 'Before prognostic temperature, iCell =', indexToCellID(iCell)
                        write(stderrUnit,*) 'thickness =', thickness(iCell)
                        write(stderrUnit,*) surfaceTemperature(iCell)
                        do k = 1, nVertLevels
                           write(stderrUnit,*) k, temperature(k,iCell)
                        enddo
                        write(stderrUnit,*) basalTemperature(iCell)
                     endif

                     ! compute initial internal energy in column (for energy conservation check)
                     initialEnergy = 0.0_RKIND
                     do k = 1, nVertLevels
                        initialEnergy = initialEnergy + temperature(k,iCell) * (layerInterfaceSigma(k+1) - layerInterfaceSigma(k)) * thickness(iCell)
                     enddo

                     ! compute matrix elements using enthalpy gradient method
                     !TODO - Get rid of temporary arrays?

                     temp(0) = surfaceTemperature(iCell)
                     temp(1:nVertLevels) = temperature(:,iCell)
                     temp(nVertInterfaces) = basalTemperature(iCell)

                     call temperature_matrix_elements(&
                          deltat,               &
                          nVertInterfaces,      &   ! CISM passes upn
                          layerCenterSigma,     &   ! CISM passes stagsigma
                          subdiagonal,          &
                          diagonal,             &
                          superdiagonal,        &
                          rhs,                  &
                          dups,                 &
                          li_mask_is_floating_ice_int(cellMask(iCell)), &
                          thickness(iCell),               &
                          temp(0:nVertInterfaces),        &
                          heatDissipation(:,iCell),       &
!!                          basalHeatFlux(iCell),           &
                          -basalHeatFlux(iCell),           &    ! CISM subroutine assumes positive down, so flip sign
                          basalFrictionFlux(iCell))

                     if (verboseColumn) then
                        write(stderrUnit,*) ' '
                        write(stderrUnit,*) 'After matrix elements, iCell =', indexToCellID(iCell)
                        write(stderrUnit,*) 'k, subd, diag, supd, rhs:'
                        do k = 1, nVertInterfaces+1
                           write(stderrUnit,*) k-1, subdiagonal(k), diagonal(k), superdiagonal(k), rhs(k)
                        enddo
                     endif

                     ! solve the tridiagonal system
                     ! Note: Temperature is indexed from 0 to nVertInterfaces, with indices 1 to nVertInterfaces-1 colocated
                     ! with layerCenterSigma values of the same index.
                     ! However, the matrix elements are indexed 1 to nVertInterfaces+1, with the first row
                     ! corresponding to the surface temperature, temp(0).

                     call tridiag_solver(&
                          subdiagonal,      &
                          diagonal,         &
                          superdiagonal,    &
                          temp(0:nVertInterfaces),  &
                          rhs)

                     ! Compute conductive flux = (k/H * dT/dsigma) at upper and lower surfaces; positive down
                     ! Assume implicit backward Euler time step.
                     
                     dTtop = temp(1) - temp(0)
                     dTbot = temp(nVertInterfaces) - temp(nVertInterfaces-1)

                     surfaceConductiveFlux(iCell) = (-iceConductivity/thickness(iCell) ) * dTtop / layerCenterSigma(1)
                     basalConductiveFlux(iCell) = (-iceConductivity/thickness(iCell) ) * dTbot / (1.0_RKIND - layerCenterSigma(nVertLevels))

                     ! copy the temperature back into the full array
                     temperature(:,iCell) = temp(1:nVertLevels)

                     if (verboseColumn) then
                        write(stderrUnit,*) ' '
                        write(stderrUnit,*) 'After prognostic temperature, iCell =', indexToCellID(iCell)
                        write(stderrUnit,*) 'thickness =', thickness(iCell)
                        write(stderrUnit,*) surfaceTemperature(iCell)
                        do k = 1, nVertLevels
                           write(stderrUnit,*) k, temperature(k,iCell)
                        enddo
                        write(stderrUnit,*) basalTemperature(iCell)
                     endif

                     ! compute final internal energy in column (for energy conservation check)
                     finalEnergy = 0.0_RKIND
                     do k = 1, nVertLevels
                        finalEnergy = finalEnergy + temperature(k,iCell) * (layerInterfaceSigma(k+1) - layerInterfaceSigma(k)) * thickness(iCell)
                     enddo

                  endif   ! temperature or enthalpy solver

                  ! Convert temperature from Celsius back to Kelvin
                  temperature(:,iCell) = temperature(:,iCell) + kelvin_to_celsius
                  surfaceTemperature(iCell) = surfaceTemperature(iCell) + kelvin_to_celsius
                  basalTemperature(iCell) = basalTemperature(iCell) + kelvin_to_celsius

                  ! Compute total dissipation rate in column (W/m^2)                                                                                                                  
                  columnHeatDissipation = 0.0_RKIND
                  do k = 1, nVertLevels
                     columnHeatDissipation = columnHeatDissipation  &
                                           + heatDissipation(k,iCell) * (layerCenterSigma(k+1) - layerCenterSigma(k))
                  enddo
                  columnHeatDissipation = columnHeatDissipation * thickness(iCell)*rhoi*cp_ice

                  ! Verify that the net input of energy into the column is equal to the change in internal energy.

                  deltaEnergy = (surfaceConductiveFlux(iCell) - basalConductiveFlux(iCell) + columnHeatDissipation) * deltat

                  !TODO - Confirm that this is a reasonable error threshold
                  if (abs((finalEnergy - initialEnergy - deltaEnergy) / deltat) > 1.0e-8_RKIND) then

                     if (verboseColumn) then
                        print*, 'Ice thickness:', thickness(iCell)
                        print*, 'config_thermal_thickness:', config_thermal_thickness
                        print*, ' '
                        print*, 'Interior fluxes:'
                        print*, 'sfc conductive flx (positive up)=', -surfaceConductiveFlux(iCell)
                        print*, 'bed conductive flx (positive up)=', -basalConductiveFlux(iCell)
                        print*, 'column heat dissipation =', columnHeatDissipation
                        print*, 'Net flux =', deltaEnergy/deltat
                        print*, ' '
                        print*, 'deltaEnergy =', deltaEnergy
                        print*, 'initialEnergy =', initialEnergy
                        print*, 'finalEnergy =', finalEnergy
                        print*, ' '
                        print*, 'Energy imbalance =', finalEnergy - initialEnergy - deltaEnergy
                        print*, ' '
                        print*, 'Basal fluxes:'
                        print*, 'frictional =', basalFrictionFlux(iCell)
                        print*, 'geothermal =', basalHeatFlux(iCell)
                        print*, 'flux for bottom melting =', basalFrictionFlux(iCell) + basalHeatFlux(iCell) + basalConductiveFlux(iCell)
                     endif   ! verboseColumn

                     write(stderrUnit,*) 'li_thermal, energy conservation error: iCell, imbalance (W/m2):', &
                                         indexToCellID(iCell), (finalEnergy - initialEnergy - deltaEnergy)/deltat
                     err = 1
                     
                  endif  ! energy conservation error

               else    ! thickness <= config_thermal_thickness

                  ! set thermal mask
                  thermalCellMask(iCell) = 0

                  ! Set temperature of thin ice to 0 C = 273.15 K
                  !TODO - For cells that have just crossed this thickness threshold, energy is not conserved here.
                  !       Keep track of the energy difference?

                  surfaceTemperature(iCell) = kelvin_to_celsius
                  basalTemperature(iCell) = kelvin_to_celsius
                  temperature(:,iCell) = kelvin_to_celsius
                  waterfrac(:,iCell) = 0.0_RKIND

               endif   ! thickness > config_thermal_thickness
  
            enddo   ! iCell

            ! Calculate basal melt rate
            ! For the standard temperature scheme, temperatures above the pressure melting point
            !  are reset to Tpmp, with excess heat contributing to basal melt.
            ! For the enthalpy scheme, internal meltwater in excess of the prescribed maximum
            !  fraction (0.01 by default) is drained to the bed.  

            call basal_melting(&
                 config_thermal_solver,         &     
                 deltat,                        &
                 nCellsSolve,              &
                 nVertInterfaces,          &
                 layerInterfaceSigma,      &
                 layerCenterSigma,         &
                 thermalCellMask,           &
                 li_mask_is_floating_ice_int(cellMask), &
                 thickness,                &
                 temperature,              &
                 basalTemperature,         &
                 waterfrac,                &
                 enthalpy,                 &
                 basalFrictionFlux,        &
!!                 basalHeatFlux,            &
                 -basalHeatFlux,            &  !TODO - Switch sign convention in subroutine
                 basalConductiveFlux,      &
                 basalWaterThickness,      &
                 basalMassBal)

            ! Before exiting, check for temperatures that are physically unrealistic.
            ! Thresholds are set at the top of this module.

            do iCell = 1, nCellsSolve

               maxtemp = maxval(temperature(:,iCell))
               mintemp = minval(temperature(:,iCell))
          
               if (maxtemp > maxtemp_threshold) then
                  write(stderrUnit,*) 'maxtemp > 0: iCell, maxtemp =', iCell, maxtemp
                  write(stderrUnit,*) 'thickness =', thickness(iCell)
                  write(stderrUnit,*) 'temperature:'
                  do k = 1, nVertLevels
                     write(stderrUnit,*) k, temperature(k,iCell)
                  enddo
                  call mpas_dmpar_global_abort("An error has occurred in li_core_finalize. Aborting...")
               endif

               if (mintemp < mintemp_threshold) then
                  write(stderrUnit,*) 'mintemp < mintemp_threshold: iCell, mintemp =', iCell, mintemp
                  write(stderrUnit,*) 'thickness =', thickness(iCell)
                  write(stderrUnit,*) 'temperature:'
                  do k = 1, nVertLevels
                     write(stderrUnit,*) k, temperature(k,iCell)
                  enddo
                  call mpas_dmpar_global_abort("An error has occurred in li_core_finalize. Aborting...")
               endif
          
            enddo   ! iCell

            ! clean up
            if (allocated(subdiagonal)) deallocate(subdiagonal)
            if (allocated(diagonal)) deallocate(diagonal)
            if (allocated(superdiagonal)) deallocate(superdiagonal)
            if (allocated(rhs)) deallocate(rhs)
            if (allocated(alpha_enth)) deallocate(alpha_enth)
            if (allocated(temp)) deallocate(temp)
            if (allocated(enth)) deallocate(enth)

         end select

         block => block % next
      enddo   ! associated(block)

      ! halo updates

      call mpas_pool_get_field(thermalPool, 'surfaceTemperature', surfaceTemperatureField)
      call mpas_dmpar_exch_halo_field(surfaceTemperatureField)

      call mpas_pool_get_field(thermalPool, 'basalTemperature', basalTemperatureField)
      call mpas_dmpar_exch_halo_field(basalTemperatureField)

      !TODO - Halo updates for components of the tracer array
!!      call mpas_pool_get_field(thermalPool, 'temperature', temperatureField)
!!      call mpas_dmpar_exch_halo_field(temperatureField)

      if (trim(config_thermal_solver) == 'enthalpy') then
!!         call mpas_pool_get_field(thermalPool, 'waterfrac', waterfracField)
!!         call mpas_dmpar_exch_halo_field(waterfracField)
      endif

      ! === error check
      if (err > 0) then
          write (stderrUnit,*) 'An error has occurred in li_thermal_solver'
      endif


    end subroutine li_thermal_solver

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  !  routine li_heat_dissipation_sia
!
!> \brief MPAS land ice heat dissipation for SIA velocity solver
!> \author William Lipscomb
!> \date   October 2015
!> \details
!>  This routine computes heat dissipation in the ice interior for the
!>  SIA velocity solver.
!-----------------------------------------------------------------------

  !TODO - Move this subroutine to another module?
  subroutine li_heat_dissipation_sia(domain, err)

     ! Compute the dissipation source term associated with strain heating,
     ! based on the shallow-ice approximation.

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
     type (mpas_pool_type), pointer :: velocityPool
     type (mpas_pool_type), pointer :: thermalPool 
     type (mpas_pool_type), pointer :: scratchPool     
    
     integer, pointer :: &
          nCells,                & ! number of cells
          nEdges,                & ! number of edges
          nVertLevels              ! number of vertical layers

     integer, dimension(:), pointer ::  &
          nEdgesOnCell             ! number of edges on each cell

     integer, dimension(:,:), pointer ::  &
          cellsOnEdge,           & ! indices for 2 cells on each edge
          edgesOnCell              ! indices for edges on each cell

     real(kind=RKIND), dimension(:), pointer :: &
          areaCell                 ! area of each cell

     real(kind=RKIND), dimension(:), pointer :: &
          dcEdge,                & ! distance between neighboring cells across edge
          dvEdge                   ! distance between eighboring vertices along edge

     real(kind=RKIND), dimension(:), pointer :: &
          layerCenterSigma         ! vertical coordinate at center of each layer
     
     real(kind=RKIND), dimension(:), pointer :: &
          thickness                ! ice thickness in cells

     real(kind=RKIND), dimension(:), pointer :: &
          slopeEdge                ! surface slope at edges

     real(kind=RKIND), dimension(:,:), pointer :: &
          flowParamA               ! flow factor in each layer of each cell, Pa^(-n) s^(-1)

     real(kind=RKIND), dimension(:,:), pointer :: &
          heatDissipation          ! interior heat dissipation in each layer of each cell (deg/s)
                                   ! output from this subroutine

     type (field2dReal), pointer :: heatDissipationEdgeField

     real (kind=RKIND), dimension(:,:), pointer :: &
          heatDissipationEdge      ! heat dissipation on edges
 
     real (kind=RKIND) :: &
          thicknessEdge,         & ! thickness averaged to edge
          weightEdge               ! edge weight for averaging to cell center

     real (kind=RKIND), dimension(:), allocatable ::  &
          flowParamAEdge           ! flow parameter averaged to edge

     real (kind=RKIND), pointer :: n    ! flow law exponent

     integer :: iCell, iCell1, iCell2, iEdge, iEdgeOnCell

     ! Here are notes from the Glimmer calculation of heat dissipation:
     !
     !     "Two methods of doing this calculation: 
     !      1. find dissipation at u-pts and then average
     !      2. find dissipation at H-pts by averaging quantities from u-pts
     !     (2) works best for eismint divide (symmetry) but (1) may be better for full expts"
     !
     ! Glimmer uses (2).
     ! Here we use the C-grid variant of (1); we find the dissipation at edges and then average to cell centers.
     ! The heating rate phi, defined on an edge, is given by
     !
     !      phi = 2 * A(T) * (sigma * rhoi * g * H * |grad(s)|)^(n+1)
     ! 
     ! where A(T) is the flow factor, sigma is the vertical coordinate of the layer,
     ! H is the ice thickness averaged to the edge, and grad(s) is the surface elevation gradient.
     !
     ! phi has units of W m^{-3}.
     ! The heat dissipation in deg/s is given by phi / (rhoi * cp_ice).

     block => domain % blocklist
     do while (associated(block))

        ! get pools
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
        call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
        call mpas_pool_get_subpool(block % structs, 'thermal', thermalPool)
        call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
        
        ! get fields from the mesh pool
        call mpas_pool_get_array(meshPool, 'nCells', nCells)
        call mpas_pool_get_array(meshPool, 'nEdges', nEdges)
        call mpas_pool_get_array(meshPool, 'nVertLevels', nVertLevels)
        call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
        call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
        call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
        call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
        call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
        call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
                
        ! get fields from the geometry pool
        !TODO - Make sure slopeEdge is up to date.
        call mpas_pool_get_array(geometryPool, 'thickness', thickness)
        call mpas_pool_get_array(geometryPool, 'slopeEdge', slopeEdge)
        
        ! get fields from the velocity pool
        call mpas_pool_get_array(geometryPool, 'flowParamA', flowParamA)
        
        ! get fields from the thermal pool
        call mpas_pool_get_array(geometryPool, 'heatDissipation', heatDissipation)
        
        ! get scratch fields
        call mpas_pool_get_field(scratchPool, 'workLevelEdge', heatDissipationEdgeField)
        call mpas_allocate_scratch_field(heatDissipationEdgeField, .true.)
        heatDissipationEdge => heatDissipationEdgeField % array
        
        ! get config parameters
        call mpas_pool_get_config(liConfigs, 'config_flowLawExponent', n)
        
        allocate(flowParamAEdge(nVertLevels))
        
        ! compute the heat dissipation on edges
        
        do iEdge = 1, nEdges
           
           ! identify the cells on this edge
           iCell1 = cellsOnEdge(1,iEdge)
           iCell2 = cellsOnEdge(2,iEdge)
           
           if (iCell1 >= 1 .and. iCell1 <= nCells .and. iCell2 >= 1 .and. iCell2 <= nCells) then  ! both cells exist

              ! average the thickness and flow parameter to the edge
              thicknessEdge = 0.5_RKIND * (thickness(iCell1) + thickness(iCell2))
              flowParamAEdge(:) = 0.5_RKIND * (flowParamA(:,iCell1) + flowParamA(:,iCell2))

              ! compute the dissipation at each level
              ! Note: n = config_flowLawExponent
              !TODO - Verify that this equation gives the right answer
              heatDissipationEdge(:,iEdge) = 2.0_RKIND * flowParamAEdge(:) * &
                   (layerCenterSigma(:) * rhoi * gravity * thicknessEdge * abs(slopeEdge)) ** (n+1.0_RKIND)

           else  ! one neighbor cell does not exist
                 !TODO = Confirm that the dissipation is not needed at such edges

              heatDissipationEdge(:,iEdge) = 0.0_RKIND

           endif
           
         enddo  ! iEdge
        

        ! average the heat dissipation to cell centers
        
        do iCell = 1, nCells

           heatDissipation(:,iCell) = 0.0_RKIND

           do iEdgeOnCell = 1, nEdgesOnCell(iCell)
           
              iEdge = edgesOnCell(iEdgeOnCell,iCell)

              !TODO - Is this the preferred way of getting the edge weights?
              weightEdge = 0.25_RKIND*dcEdge(iEdge)*dvEdge(iEdge) / areaCell(iCell)

              heatDissipation(:,iCell) = heatDissipation(:,iCell) + weightEdge * heatDissipationEdge(:,iEdge)

           enddo  ! iEdgeOnCell

        enddo  ! iCell

        ! convert units from W/m^3 to deg/s
        !TODO - Confirm that deg/s are the desired units.  Might want to go with W/m^3?

        heatDissipation(:,:) = heatDissipation(:,:) / (rhoi * cp_ice)

     enddo  ! associated(block) 

   end subroutine li_heat_dissipation_sia


!***********************************************************************
!***********************************************************************
! Private subroutines:
!***********************************************************************
!***********************************************************************

!TODO - Add subroutine headers
!       Clean up subroutines and code in MPAS style/
!       Change variable names to agree with the driver subroutine

    subroutine temperature_matrix_elements(dttem,                        &   ! deltat
                                           upn,          stagsigma,      &   ! upn = nVertInterfaces; stagsigma = layerCenterSigma
                                           subd,         diag,           &
                                           supd,         rhsd,           &
                                           dups,         floating_mask,  &
                                           thck,         temp,           &
                                           dissip,                       &
                                           bheatflx,     bfricflx)

      ! solve for tridiagonal entries of sparse matrix

      ! Note: Matrix elements (subd, supd, diag, rhsd) are indexed from 1 to upn+1,
      !       whereas temperature is indexed from 0 to upn.
      !       The first row of the matrix is the equation for temperature(0),
      !       the last row is the equation for temperature(upn), and so on.
      
      real(kind=RKIND), intent(in) :: dttem       ! time step (s)
      integer, intent(in) :: upn          ! number of layer interfaces
      real(kind=RKIND), dimension(upn-1), intent(in) :: stagsigma    ! sigma coordinate at temp nodes
      real(kind=RKIND), dimension(:), intent(out) :: subd, diag, supd, rhsd
      real(kind=RKIND), dimension(:,:), intent(in) :: dups           ! vertical grid quantities
      integer, intent(in) :: floating_mask
      real(kind=RKIND), intent(in) ::  thck       ! ice thickness (m)
      real(kind=RKIND), dimension(0:upn), intent(in) ::  temp     ! ice temperature (deg C)
      real(kind=RKIND), dimension(upn-1), intent(in) :: dissip     ! interior heat dissipation (deg/s)
      real(kind=RKIND), intent(in) :: bheatflx    ! geothermal flux (W m-2), positive down
      real(kind=RKIND), intent(in) :: bfricflx    ! basal friction heat flux (W m-2), >= 0

      ! local variables

      real(kind=RKIND) :: pmptemp_bed  ! pressure melting temp at bed
      real(kind=RKIND) :: fact
      real(kind=RKIND) :: dsigbot      ! bottom layer thicknes in sigma coords

      ! Compute subdiagonal, diagonal, and superdiagonal matrix elements

      ! upper boundary: set to surface air temperature

      supd(1) = 0.0_RKIND
      subd(1) = 0.0_RKIND
      diag(1) = 1.0_RKIND
      rhsd(1) = temp(0)
      
      ! ice interior, layers 1:upn-1  (matrix elements 2:upn)

      fact = dttem * iceConductivity / (rhoi*cp_ice) / thck**2
      subd(2:upn) = -fact * dups(1:upn-1,1)
      supd(2:upn) = -fact * dups(1:upn-1,2)
      diag(2:upn) = 1.0_RKIND - subd(2:upn) - supd(2:upn)
      rhsd(2:upn) = temp(1:upn-1) + dissip(1:upn-1)*dttem

      ! basal boundary:
      ! for grounded ice, a heat flux is applied
      ! for floating ice, the basal temperature is held constant

      !Note: If T(upn) < T_pmp, then require dT/dsigma = H/k * (G + taub*ubas)
      !       That is, net heat flux at lower boundary must equal zero.
      !      If T(upn) >= Tpmp, then set T(upn) = Tpmp
      
      if (floating_mask == 1) then

         supd(upn+1) = 0.0_RKIND
         subd(upn+1) = 0.0_RKIND
         diag(upn+1) = 1.0_RKIND
         rhsd(upn+1) = temp(upn) 
         
      else    ! grounded ice

         call pressure_melting_point(thck, pmptemp_bed)

         if (abs(temp(upn) - pmptemp_bed) < 0.001_RKIND) then

            ! hold basal temperature at pressure melting point

            supd(upn+1) = 0.0_RKIND
            subd(upn+1) = 0.0_RKIND
            diag(upn+1) = 1.0_RKIND
            rhsd(upn+1) = pmptemp_bed

         else   ! frozen at bed
                ! maintain balance of heat sources and sinks
                ! (conductive flux, geothermal flux, and basal friction)

            ! Note: bheatflx is generally <= 0, since defined as positive down.

            ! calculate dsigma for the bottom layer between the basal boundary and the temp. point above
            dsigbot = 1.0_RKIND - stagsigma(upn-1)

            ! backward Euler flux basal boundary condition
            subd(upn+1) = -1.0_RKIND
            supd(upn+1) =  0.0_RKIND
            diag(upn+1) =  1.0_RKIND
            rhsd(upn+1) = (bfricflx - bheatflx) * dsigbot*thck / iceConductivity
            
         endif   ! melting or frozen
         
      end if     ! floating or grounded

    end subroutine temperature_matrix_elements


    subroutine enthalpy_matrix_elements(dttem,                       &   ! deltat
                                        upn,       stagsigma,        &   ! upn = nVertInterfaces; stagsigma = layerCenterSigma
                                        subd,      diag,             &
                                        supd,      rhsd,             &
                                        dups,      floating_mask,    &
                                        thck,                        &
                                        temp,      waterfrac,        &
                                        enthalpy,  dissip,           &
                                        bheatflx,  bfricflx,         &
                                        alpha_enth,                  &
                                        verbose_column_in)

      ! solve for tridiagonal entries of sparse matrix

      ! Note: Matrix elements (subd, supd, diag, rhsd) are indexed from 1 to upn+1,
      ! whereas temperature/enthalpy is indexed from 0 to upn.
      ! The first row of the matrix is the equation for enthalpy(0),
      ! the last row is the equation for enthalpy(upn), and so on.

      !I/O variables
      real(kind=RKIND), intent(in) :: dttem       ! time step (s)
      integer, intent(in) :: upn                  ! number of vertical interfaces
      real(kind=RKIND), dimension(upn-1), intent(in) :: stagsigma    ! sigma coordinate at temp/enthalpy nodes
      real(kind=RKIND), dimension(:,:), intent(in) :: dups           ! vertical grid quantities
      real(kind=RKIND), dimension(:), intent(out) :: subd, diag, supd, rhsd   ! matrix elements
      integer, intent(in) :: floating_mask
      real(kind=RKIND), intent(in) :: thck                         ! ice thickness (m)
      real(kind=RKIND), dimension(0:upn), intent(in) :: temp       ! temperature (deg C)            !TODO - Do units matter? K or C?
      real(kind=RKIND), dimension(upn-1), intent(in) :: waterfrac  ! water fraction (unitless)
      real(kind=RKIND), dimension(0:upn), intent(in) :: enthalpy   ! specific enthalpy (J/m^3)
      real(kind=RKIND), dimension(upn-1), intent(in) :: dissip     ! interior heat dissipation (deg/s)
      real(kind=RKIND), intent(in) :: bheatflx   ! geothermal flux (W m-2), positive down           !TODO - Flip sign to positive up
      real(kind=RKIND), intent(in) :: bfricflx   ! basal friction heat flux (W m-2), >= 0
      real(kind=RKIND), dimension(:), intent(out) :: alpha_enth  ! half-node diffusivity (m^2/s) for enthalpy
                                                                 ! located halfway between temperature points
      
      logical, intent(in), optional :: verbose_column_in   ! if true, print debug statements for this column
    
      ! local variables
      real(kind=RKIND) :: dsigbot  ! bottom layer thicknes in sigma coords.
      real(kind=RKIND) :: alphai   ! cold ice diffusivity
      real(kind=RKIND) :: alpha0   ! temperate ice diffusivity
      real(kind=RKIND) :: fact     ! coefficient in tridiag matrix
      real(kind=RKIND), dimension(1:upn-1) :: pmptemp    ! pressure melting point temp in interior (deg C)
      real(kind=RKIND) :: pmptemp_bed                    ! pressure melting point temp at bed (deg C)
      real(kind=RKIND), dimension(0:upn) :: enth_T       ! temperature part of specific enthalpy (J/m^3)
      real(kind=RKIND) :: denth       ! enthalpy difference between adjacent layers
      real(kind=RKIND) :: denth_T     ! difference in temperature component of enthalpy between adjacent layers
      real(kind=RKIND) :: alpha_fact  ! factor for averaging diffusivity, 0 <= fact <= 1
      logical :: verbose_column       ! if true, print debug statements for this column
      integer  :: k
      
      logical, parameter :: &
           alpha_harmonic_avg  = .false.  ! if true, take harmonic average of alpha in adjacent layers
                                          ! if false, take arithmetic average

      if (present(verbose_column_in)) then
         verbose_column = verbose_column_in
      else
         verbose_column = .false.
      endif

      ! define diffusivities alpha_i and alpha_0
      !!    alphai = coni / rhoi / cp_ice
      alphai = iceConductivity / rhoi / cp_ice
      alpha0 = alphai / 100.0_RKIND
        
      ! find pmptemp for this column (interior nodes and boundary)
      call pressure_melting_point_column(stagsigma(1:upn-1), thck, pmptemp(1:upn-1))
      call pressure_melting_point(thck, pmptemp_bed)

      !WHL - debug                                                                                                                       
      if (verbose_column) then
         print*, ' '
         print*, 'Computing enthalpy matrix elements'
         print*, 'k, temp, wfrac, enthalpy/(rhoi*ci), pmpt:'
         k = 0
         print*, k, temp(k), 0.0_RKIND, enthalpy(k)/(rhoi*cp_ice)
         do k = 1, upn-1
            print*, k, temp(k), waterfrac(k), &
                 enthalpy(k)/(rhoi*cp_ice), pmptemp(k)
         enddo
         k = upn
         print*, k, temp(k), 0.0_RKIND, enthalpy(k)/(rhoi*cp_ice), pmptemp_bed
      endif
      
      !--------------------------------------------------------------------
      !WHL - Commenting out the following and replacing it with a new way of computing alpha.
      !      The commented-out code can result in sudden large changes in alpha that
      !       lead to oscillations in the thickness, temperature and velocity fields.
      !      These oscillations have a period of ~1 yr or more, spatial scale of
      !       many grid cells, and amplitude of ~10 m in thickness, 1 deg in temperature,
      !       and 2 m/s in velocity.
      
      ! create a column vector of size (0:upn) of diffusivity based on 
      ! previous timestep's temp.  Boundary nodes need a value so half-node
      ! diffusivity can be calculated at interior nodes (1:upn-1)
      
!    do k = 0,upn
!       if (temp(k) < pmptemp(k)) then
!          alpha(k) = alphai
!       else
!          alpha(k) = alpha0
!       endif
!    end do
    
    ! Find half-node diffusivity using harmonic average between nodes.
    ! The vector will be size (1:upn) - the first value is the half-node
    ! between nodes 0 and 1, the last value is the half-node between
    ! nodes upn-1 and upn. 
    
!    do k = 1,upn
!       alpha_enth(k) = 2.0_RKIND / ((1.0_RKIND/alpha(k-1)) + (1.0_RKIND/alpha(k)))
!    end do

      ! end of commented-out method
      !--------------------------------------------------------------------
       
      !--------------------------------------------------------------------
      !WHL - Trying a different approach to the diffusivity at layer interfaces.
      ! Let d(enth)/dz = the gradient of enthalpy
      ! Can write
      !    d(enth)/dz = d(enth_T)/dz + d(enth_w)/dz,
      ! where
      !    enth_T = (1-phi_w) * rhoi*ci*T
      !    enth_w =    phi_w  * rhoo*(L + ci*Tpmp)
      !
      ! Now let f = d(enth_T)/z / d(enth)/dz
      !   (f -> 0 if f is computed to be negative)
      ! For cold ice, f = 1 and alpha = alphai
      ! For temperate ice, f ~ 0 and alpha = alpha0
      ! At the interface between cold and temperate ice,
      !  f ~ 0 if the temperate ice has large phi_w, but
      !  f ~ 1 if the temperate ice has close to zero phi_w.
      ! Two ways to average:
      ! (1) arithmetic average:  alpha = f*alphai + (1-f)*alpha0
      ! (2) harmonic average:    alpha = 1 / (f/alphai + (1-f)/alpha0).
      ! Both methods have the same asymptotic values at f = 0 or 1,
      !  but the arithmetic average gives greater diffusivity for
      !  intermediate values.
      !
      ! Still to be determined which is more accurate.
      ! The harmonic average allows large temperature gradients between the 
      !  bottom layer and the next layer up; the arithmetic average gives
      !  smoother gradients.
      !--------------------------------------------------------------------
      !
      ! At each temperature point, compute the temperature part of the enthalpy.
      ! enth_T = enth for cold ice, enth_T < enth for temperate ice
      
      do k = 0, upn
         enth_T(k) = (1.0_RKIND - waterfrac(k)) * rhoi*cp_ice*temp(k)
      enddo
      
      !WHL - debug
      if (verbose_column) then
         print*, ' '
         print*, 'k, denth_T/(rhoi*cp_ice), denth/(rhoi*cp_ice), alpha_fact, alpha_enth(up):'
      endif

      ! Compute factors relating the temperature gradient to the total enthalpy gradient.
      ! Use these factors to average the diffusivity between adjacent temperature points.
      do k = 1,upn
         denth   = enthalpy(k) - enthalpy(k-1)
         denth_T = enth_T(k) - enth_T(k-1)   ! = denth in cold ice, < denth in temperate ice
         if (abs(denth) > 1.e-20_RKIND * rhoo*latent_heat_ice) then
            alpha_fact = max(0.0_RKIND, denth_T/denth)
            alpha_fact = min(1.0_RKIND, alpha_fact)
         else
            alpha_fact = 0.0_RKIND
         endif
         
         if (alpha_harmonic_avg) then  ! take a harmonic average
                                       ! This gives slower cooling of temperate layers and allows
                                       !  large temperature gradients between cold and temperate layers
            alpha_enth(k) = 1.0_RKIND / ((alpha_fact/alphai) + (1.0_RKIND-alpha_fact)/alpha0)
         else   ! take an arithmetic average
                ! This gives faster cooling of temperate layers and smaller gradients 
            alpha_enth(k) = alpha_fact*alphai + (1.0_RKIND-alpha_fact)*alpha0
         endif

         !WHL - debug
         if (verbose_column) then
            print*, k, denth_T/(rhoi*cp_ice), denth/(rhoi*cp_ice), alpha_fact, alpha_enth(k)
         endif
         
      end do

      ! Compute subdiagonal, diagonal, and superdiagonal matrix elements
      ! Assume backward Euler time stepping
    
      ! upper boundary: set to surface air temperature
      supd(1) = 0.0_RKIND
      subd(1) = 0.0_RKIND
      diag(1) = 1.0_RKIND
      rhsd(1) = min(0.0_RKIND,temp(0)) * rhoi*cp_ice
      
      ! ice interior, layers 1:upn-1  (matrix elements 2:upn)

      fact = dttem / thck**2

      subd(2:upn) = -fact * alpha_enth(1:upn-1) * dups(1:upn-1,1)                                
      supd(2:upn) = -fact * alpha_enth(2:upn) * dups(1:upn-1,2)                                
      diag(2:upn) = 1.0_RKIND - subd(2:upn) - supd(2:upn)                                
      rhsd(2:upn) = enthalpy(1:upn-1) + dissip(1:upn-1)*dttem * rhoi * cp_ice
      
      ! BDM I'm assuming that dissip has units of phi/rhoi/cp_ice.
      ! For an enthalpy calc, we want just phi, hence dissip * rhoi * cp_ice
        
      ! basal boundary:
      ! for grounded ice, a heat flux is applied
      ! for floating ice, the basal temperature is held constant

      !NOTE: This lower BC is different from the one in glide_temp.
      !      If T(upn) < T_pmp, then require dT/dsigma = H/k * (G + taub*ubas)
      !       That is, net heat flux at lower boundary must equal zero.
      !      If T(upn) >= Tpmp, then set T(upn) = Tpmp
      
      if (floating_mask == 1) then

         supd(upn+1) = 0.0_RKIND 
         subd(upn+1) = 0.0_RKIND
         diag(upn+1) = 1.0_RKIND
         rhsd(upn+1) = enthalpy(upn)
         
      else    ! grounded ice
         
         !WHL - debug
         if (verbose_column) then
            k = upn-1
            print*, 'temp(upn-1), pmptemp(upn-1):', temp(k), pmptemp(k)
            k = upn
            print*, 'temp(upn), pmptemp(upn):', temp(k), pmptemp_bed
         endif
         
         ! Positive-Thickness Basal Temperate Boundary Layer

         !WHL - Not sure whether this condition is ideal.
         !      It implies that the enthalpy at the bed (upn) = enthalpy in layer (upn-1). 
         if (abs(temp(upn-1) - pmptemp(upn-1)) < 0.001_RKIND) then   
       
            subd(upn+1) = -1.0_RKIND
            supd(upn+1) =  0.0_RKIND 
            diag(upn+1) = 1.0_RKIND 
            rhsd(upn+1) = 0.0_RKIND
            
            !WHL - debug
            if (verbose_column) then
               print*, 'basal BC: branch 1 (finite-thck BL)'
            endif

            !Zero-Thickness Basal Temperate Boundary Layer
         elseif (abs(temp(upn) -  pmptemp_bed) < 0.001_RKIND) then  ! melting
          
            ! hold basal temperature at pressure melting point
            supd(upn+1) = 0.0_RKIND
            subd(upn+1) = 0.0_RKIND
            diag(upn+1) = 1.0_RKIND
            rhsd(upn+1) = pmptemp_bed * rhoi * cp_ice
            
            !WHL - debug
            if (verbose_column) then
               print*, 'basal BC: branch 2 (zero-thck BL)'
            endif
            
         else  
          
            !WHL - debug
            if (verbose_column) then
               print*, 'basal BC: branch 3 (cold ice)'
            endif
            
            ! frozen at bed
            ! maintain balance of heat sources and sinks
            ! (conductive flux, geothermal flux, and basal friction)
            
            ! Note: The heat source due to basal sliding (bfricflx) is computed in subroutine calcbfric.
            ! Also note that bheatflx is generally <= 0, since defined as positive down.
            
            ! calculate dsigma for the bottom layer between the basal boundary and the temp. point above
            dsigbot = (1.0_RKIND - stagsigma(upn-1))                                                                  
            
            ! =====Backward Euler flux basal boundary condition=====
            ! MJH: If Crank-Nicolson is desired for the b.c., it is necessary to
            ! ensure that the i.c. temperature for the boundary satisfies the
            ! b.c. - otherwise oscillations will occur because the C-N b.c. only
            ! specifies the basal flux averaged over two consecutive time steps.
            subd(upn+1) = -1.0_RKIND
            supd(upn+1) =  0.0_RKIND 
            diag(upn+1) =  1.0_RKIND 
            rhsd(upn+1) = (bfricflx - bheatflx) * dsigbot*thck * rhoi*cp_ice/iceConductivity
            ! BDM temp approach should work out to be dT/dsigma, so enthalpy approach
            ! should just need dT/dsigma * rhoi * cp_ice for correct units
            
         endif   ! melting or frozen
         
      end if     ! floating or grounded
      
    end subroutine enthalpy_matrix_elements
    
    subroutine basal_melting(&
                             config_thermal_solver,     &
                             dttem,         &
                             nCellsSolve,                  &
                             upn,                          &
                             sigma,         stagsigma,     &
                             ice_mask,      floating_mask, &
                             thck,                         &
                             temp,          basalTemp,     &
                             waterfrac,     enthalpy,      &
                             bfricflx,      bheatflx,      &
                             lcondflx,                     &
                             bwat,          bmlt)

    ! Compute the rate of basal melting.
    ! The basal melting computed here is applied to the ice thickness
    !  by glissade_transport_driver, conserving mass and energy.
    !
    ! For the standard prognostic temperature scheme, any internal temperatures 
    !  above the pressure melting point are reset to Tpmp.  Excess energy 
    !  is applied toward melting with immediate drainage to the bed. 
    ! For the enthalpy scheme, any meltwater in excess of the maximum allowed
    !  meltwater fraction (0.01 by default) is drained to the bed.
    !
    !TODO - Deal with basal melting for floating ice

    !-----------------------------------------------------------------
    ! Input/output arguments
    !-----------------------------------------------------------------

    character(len=StrKIND), intent(in) :: &
         config_thermal_solver                ! thermal solver option (temperature or enthalpy)

    real(kind=RKIND), intent(in) :: dttem     ! time step (s)
    integer, intent(in) :: nCellsSolve        ! number of locally owned cells
    integer, intent(in) :: upn                ! number of vertical interfaces
                                              !TODO - change to nVertInterfaces
    real(kind=RKIND), dimension(upn),    intent(in) :: sigma           ! vertical sigma coordinate
    real(kind=RKIND), dimension(upn-1),  intent(in) :: stagsigma       ! staggered vertical coordinate for temperature

!!    real(kind=RKIND), dimension(0:,:,:), intent(inout) :: temp         ! temperature (deg C)
    real(kind=RKIND), dimension(:,:), intent(inout) :: temp              ! temperature (deg C) in each layer
    real(kind=RKIND), dimension(:), intent(inout) :: basalTemp         ! basal temperature (deg C)

    real(kind=RKIND), dimension(:,:),  intent(inout) :: waterfrac    ! water fraction

    real(kind=RKIND), dimension(:,:), intent(in) :: enthalpy        ! enthalpy
!!    real(kind=RKIND), dimension(:), intent(in) :: basalEnthalpy        ! basal enthalpy

    real(kind=RKIND), dimension(:),   intent(in) :: &
         thck,                 & ! ice thickness (m)
         bfricflx,             & ! basal frictional heating flux (W m-2), >= 0
         bheatflx,             & ! geothermal heating flux (W m-2), positive down   !TODO - Switch sign convention
         lcondflx,             & ! heat conducted from ice interior to bed (W m-2), positive down
         bwat                    ! depth of basal water (m)

    integer, dimension(:), intent(in) ::  &
         ice_mask,           &! = 1 where ice exists (thickness > config_thermal_thickness), else = 0
         floating_mask        ! = 1 where ice is floating, else = 0

    !TODO - Fix sign of bmlt before passing out
    real(kind=RKIND), dimension(:), intent(out):: bmlt    ! melt rate (m/s)
                                                    ! > 0 for melting, < 0 for freeze-on

    !-----------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------

    integer :: k, iCell
    real(kind=RKIND), dimension(upn-1)  :: pmptemp   ! pressure melting point temp in ice interior
    real(kind=RKIND) :: pmptemp_bed     ! pressure melting point temp at bed
    real(kind=RKIND) :: bflx          ! heat flux available for basal melting (W/m^2)
    real(kind=RKIND) :: layer_thck    ! layer thickness (m)
    real(kind=RKIND) :: melt_energy   ! energy available for internal melting (J/m^2)
    real(kind=RKIND) :: internal_melt_rate   ! internal melt rate, transferred to bed (m/s)
    real(kind=RKIND) :: melt_fact     ! factor for bmlt calculation
    real(kind=RKIND) :: hmlt          ! melt thickness associated with excess meltwater

    real(kind=RKIND), parameter :: &
         max_waterfrac = 0.01_RKIND   ! maximum allowed water fraction; excess drains to bed

    real(kind=RKIND), parameter :: &
         eps11 = 1.0e-11_RKIND        ! small number

    bmlt(:) = 0.0_RKIND
    melt_fact = 1.0_RKIND / (latent_heat_ice * rhoi) 
    
    do iCell = 1, nCellsSolve
       
       if (ice_mask(iCell) == 1 .and. floating_mask(iCell) == 0) then  ! ice is present and grounded

          !TODO - Switch sign conventions for bmlt and bheatflux
          ! Compute basal melting
          ! Note: bmlt > 0 for melting, < 0 for freeze-on   
          !       bfricflx >= 0 by definition
          !       bheatflx is positive down, so usually bheatflx < 0 (with negative values contributing to melt)
          !       lcondflx is positive down, so lcondflx < 0 for heat is flowing from the bed toward the surface
          !
          !       This equation allows for freeze-on (bmlt < 0) if the conductive term 
          !        (lcondflx, positive down) is carrying enough heat away from the boundary.  
          !       But freeze-on requires a local water supply, bwat > 0.
          !       When bwat = 0, we reset the bed temperature to a value slightly below the melting point.
          !
          !TODO - For the enthalpy scheme, deal with the rare case that the bottom layer melts completely
          !       and overlying layers with a different enthalpy also melt.
          
          bflx = bfricflx(iCell) + lcondflx(iCell) - bheatflx(iCell)  ! W/m^2

          if (abs(bflx) < eps11) then  ! bflx might be slightly different from zero because of rounding errors; if so, then zero out
             bflx = 0.0_RKIND
          endif

          if (trim(config_thermal_solver) == 'enthalpy') then
!!             bmlt(ew,ns) = bflx / (lhci*rhoi - enthalpy(upn,ew,ns))
             !Note: basalEnthalpy = rhoi*cp_ice*basalTemp
             bmlt(iCell) = bflx / (latent_heat_ice*rhoi - rhoi*cp_ice*basalTemp(iCell))   
          else   ! temperature solver
             bmlt(iCell) = bflx * melt_fact   ! m/s   (melt_fact = 1/(rhoi*latent_heat_ice)
          endif
          
          ! Add internal melting

          if (trim(config_thermal_solver) == 'enthalpy') then

             ! Add internal melting associated with waterfrac > max_waterfrac (1%)

             !TODO - Any correction for rhoi/rhow here?  Or melting ice that is already partly melted?
             do k = 1, upn-1
                if (waterfrac(k,iCell) > max_waterfrac) then
                   ! compute melt rate associated with excess water
                   hmlt = (waterfrac(k,iCell) - max_waterfrac) * thck(iCell) * (sigma(k+1) - sigma(k))  ! m
                   internal_melt_rate = hmlt / dttem          ! m/s
                   ! transfer meltwater to the bed
                   bmlt(iCell) = bmlt(iCell) + internal_melt_rate      ! m/s
                   ! reset waterfrac to max value
                   waterfrac(k,iCell) = max_waterfrac
                endif
             enddo

          else   ! temperature solver

             ! Add internal melting associated with T > Tpmp

             call pressure_melting_point_column(&
                  stagsigma,     &
                  thck(iCell),   &
                  pmptemp)

             do k = 1, upn-1
                if (temp(k,iCell) > pmptemp(k)) then
                   ! compute excess energy available for melting
                   layer_thck = thck(iCell) * (sigma(k+1) - sigma(k))  ! m
                   melt_energy = rhoi * cp_ice * (temp(k,iCell) - pmptemp(k)) * layer_thck    ! J/m^2
                   internal_melt_rate = melt_energy / (rhoi * latent_heat_ice * dttem)  ! m/s
                   ! transfer internal melting to the bed
                   bmlt(iCell) = bmlt(iCell) + internal_melt_rate  ! m/s
                   ! reset T to Tpmp
                   temp(k,iCell) = pmptemp(k)
                endif
             enddo

          endif    ! config_thermal_solver

          ! Cap basal temp at pressure melting point, if necessary

          call pressure_melting_point(&
               thck(iCell), &
               pmptemp_bed)

          temp(upn,iCell) = min (temp(upn,iCell), pmptemp_bed)
          
          ! If freeze-on was computed above (bmlt < 0) and Tbed = Tpmp but no basal water is present, then set T(upn) < Tpmp.
          ! Note: In the matrix element subroutines, we solve for Tbed (instead of holding it at Tpmp) when Tbed < -0.001.
          !       With an offset here of 0.01, we will solve for T_bed at the next timestep.
          ! Note: I don't think energy conservation is violated here, because no energy is associated with
          !       the infinitesimally thin layer at the bed.
          
          if (bmlt(iCell) < 0.0_RKIND .and. bwat(iCell) == 0.0_RKIND .and. temp(upn,iCell) >= pmptemp_bed) then
             temp(upn,iCell) = pmptemp_bed - 0.01_RKIND
          endif
          
       endif   ! ice is present and grounded
       
    enddo   ! iCell

  end subroutine basal_melting

   
    subroutine temperature_to_enthalpy(&
         layerCenterSigma,  &
         thickness,         &
         temperature,       &
         waterfrac,         &
         enthalpy)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(in) :: &
           layerCenterSigma      !< Input: sigma coordinate at midpoint of each layer

      real (kind=RKIND), intent(in) ::  &
           thickness             !< Input: ice thickness
 
      real (kind=RKIND), dimension(:), intent(in) :: &
           temperature,        & !< Input: interior ice temperature
           waterfrac             !< Input: interior water fraction

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

      real (kind=RKIND), dimension(:), intent(out) :: &
           enthalpy              !< Output:  interior ice enthalpy

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(size(layerCenterSigma)) :: pmpTemperature

      integer :: k, nVertLevels

      nVertLevels = size(layerCenterSigma)

      ! Find pressure melting point temperature in column

      call pressure_melting_point_column(&
           layerCenterSigma,  &
           thickness,         &
           pmpTemperature)

      ! Solve for enthalpy

      do k = 1, nVertLevels
         enthalpy(k) = (1.0_RKIND - waterfrac(k)) * rhoi * cp_ice * temperature(k)   &
                      + waterfrac(k) * rhoo * (cp_ice * pmpTemperature(k) + latent_heat_ice)
      end do

    end subroutine temperature_to_enthalpy


    subroutine enthalpy_to_temperature(&
         layerCenterSigma,  &
         thickness,         &
         enthalpy,          &
         temperature,       &
         waterfrac)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(in) :: &
           layerCenterSigma      !< Input: sigma coordinate at midpoint of each layer

      real (kind=RKIND), intent(in) ::  &
           thickness             !< Input: ice thickness
 
      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(0:), intent(inout) :: &
           enthalpy              !< Input/output:  interior ice enthalpy

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(0:), intent(out) :: &
           temperature           !< Output: interior ice temperature

      real (kind=RKIND), dimension(:), intent(out) :: &
           waterfrac             !< Output: interior water fraction

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(size(layerCenterSigma)) :: pmpTemperature
      real (kind=RKIND) :: pmpTemperatureBed
      real (kind=RKIND), dimension(0:size(layerCenterSigma)+1) :: pmpEnthalpy

      integer :: k, nVertLevels

      nVertLevels = size(layerCenterSigma)

      ! Commented-out code below is from CISM
!    real(dp), dimension(size(stagsigma))                 :: pmptemp     ! (1:upn-1)
!    real(dp)                                             :: pmptemp_bed
!    real(dp), dimension(0:size(stagsigma)+1)             :: pmpenthalpy ! (0:upn)
!    integer :: up, upn
        
!    upn = size(stagsigma) + 1
        
      ! Find pressure melting point temperature in ice interior
      call pressure_melting_point_column(&
           layerCenterSigma,  &
           thickness,         &
           pmpTemperature)

      ! find pressure melting point temperature at bed
      call pressure_melting_point(&
           thickness, &
           pmpTemperatureBed)

!    upn = size(stagsigma) + 1
        
      ! find pressure melting point enthalpy in the column
      pmpEnthalpy(0) = 0.0_RKIND
      pmpEnthalpy(1:nVertLevels) = pmpTemperature(1:nVertLevels) * rhoi*cp_ice
      pmpEnthalpy(nVertLevels+1) = pmpTemperatureBed * rhoi*cp_ice

!    call glissade_pressure_melting_point_column(thck, stagsigma(1:upn-1), pmptemp(1:upn-1))
!    call glissade_pressure_melting_point(thck, pmptemp_bed)
!    pmpenthalpy(0) = 0.d0
!    pmpenthalpy(1:upn-1) = pmptemp(1:upn-1) * rhoi*shci
!    pmpenthalpy(upn) = pmptemp_bed * rhoi*shci

      ! solve for temperature and waterfrac
      
      ! upper surface
      if (enthalpy(0) >= pmpEnthalpy(0)) then   ! temperature ice
         temperature(0) = 0.0_RKIND
         ! Reset enthalpy to be consistent with the surface temperature.
         ! This is consistent with energy conservation because the top surface
         !  is infinitesimally thin.
         enthalpy(0) = pmpEnthalpy(0)
      else   ! cold ice
         temperature(0) = enthalpy(0) / (rhoi*cp_ice)
      endif

!    if (enthalpy(0) >= pmpenthalpy(0)) then ! temperate ice
!       temp(0) = 0.d0                      ! temperate ice
       ! Reset enthalpy to be consistent with the surface temperature.
       ! This is consistent with energy conservation because the top surface
       !  is infinitesimally thin.
!       enthalpy(0) = pmpenthalpy(0)
!    else
!       temp(0) = enthalpy(0) / (rhoi*shci) ! cold ice
!    endif

      ! interior
      do k = 1, nVertLevels
         if (enthalpy(k) >= pmpEnthalpy(k)) then   ! temperate ice
            temperature(k) = pmpTemperature(k)
            waterfrac(k) = (enthalpy(k) - pmpenthalpy(k)) /        &
                          ((rhoo-rhoi)*cp_ice*pmpTemperature(k) + rhoo*latent_heat_ice)
         else   ! cold ice
            temperature(k) = enthalpy(k) / (rhoi*cp_ice)
            waterfrac(k) = 0.0_RKIND
         endif
      enddo   ! k

!    do up = 1, upn-1
!       if (enthalpy(up) >= pmpenthalpy(up)) then ! temperate ice
!          temp(up) = pmptemp(up)
!          waterfrac(up) = (enthalpy(up)-pmpenthalpy(up)) /                 &
!                          ((rhow-rhoi) * shci * pmptemp(up) + rhow * lhci)
!       else ! cold ice
!          temp(up) = enthalpy(up) / (rhoi*shci)
!          waterfrac(up) = 0.0d0
!       endif
!    end do
        
      ! bed
      k = nVertLevels + 1
      if (enthalpy(k) >= pmpEnthalpy(k)) then   ! temperature ice
         temperature(k) = 0.0_RKIND
         ! Reset enthalpy to be consistent with the surface temperature.
         ! This is consistent with energy conservation because the top surface
         !  is infinitesimally thin.
         enthalpy(k) = pmpEnthalpy(k)
      else   ! cold ice
         temperature(k) = enthalpy(k) / (rhoi*cp_ice)
      endif

!    if (enthalpy(upn) >= pmpenthalpy(upn)) then  ! temperate ice
!       temp(upn) = pmptemp_bed
       ! Reset enthalpy to be consistent with the bed temperature.
       ! This is consistent with energy conservation because the basal surface
       !  is infinitesimally thin.
!       enthalpy(upn) = pmpenthalpy(upn)
!    else
!       temp(upn) = enthalpy(upn) / (rhoi*shci)   ! cold ice
!    endif

    end subroutine enthalpy_to_temperature


    subroutine pressure_melting_point_column(&
         layerCenterSigma,  &
         thickness,         &
         pmpTemperature)
      
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(in) :: &
           layerCenterSigma      !< Input: sigma coordinate at midpoint of each layer

      real (kind=RKIND), intent(in) ::  &
           thickness             !< Input: ice thickness
 
      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(out) :: &
           pmpTemperature          !< Output: pressure melting point temperature

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      pmpTemperature(:) = - iceMeltingPointPressureDependence * rhoi * gravity * thickness * layerCenterSigma(:) 

    end subroutine pressure_melting_point_column


    subroutine pressure_melting_point(&
         depth,           &
         pmpTemperature)
      
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(in) ::  &
           depth                !< Input: depth in column
 
      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(out) :: &
           pmpTemperature          !< Output: pressure melting point temperature

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      pmpTemperature = - iceMeltingPointPressureDependence * rhoi * gravity * depth

    end subroutine pressure_melting_point


    !TODO - Move the tridiag solver to a utility module?
    subroutine tridiag_solver(a,b,c,x,y)

      real(kind=RKIND), dimension(:), intent(in)  :: a !< Input: Lower diagonal; a(1) is ignored
      real(kind=RKIND), dimension(:), intent(in)  :: b !< Input: Main diagonal
      real(kind=RKIND), dimension(:), intent(in)  :: c !< Input: Upper diagonal; c(n) is ignored
      real(kind=RKIND), dimension(:), intent(in)  :: y !< Input: Right-hand side
      real(kind=RKIND), dimension(:), intent(out) :: x !< Output: Unknown vector
      
      real(kind=RKIND),dimension(size(a)) :: aa
      real(kind=RKIND),dimension(size(a)) :: bb

      integer :: n,i

      n = size(a)
      
      aa(1) = c(1)/b(1)
      bb(1) = y(1)/b(1)
      
      do i = 2,n
         aa(i) = c(i)/(b(i)-a(i)*aa(i-1))
         bb(i) = (y(i)-a(i)*bb(i-1)) / (b(i)-a(i)*aa(i-1))
      end do
      
      x(n) = bb(n)

      do i = n-1,1,-1
         x(i) = bb(i) - aa(i)*x(i+1)
      end do
      
    end subroutine tridiag_solver

    !***********************************************************************

  end module li_thermal

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||   



