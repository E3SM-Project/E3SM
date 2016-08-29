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
   use mpas_timer
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

   public :: li_thermal_init, li_thermal_solver, &
             li_init_linear_temperature_in_column,  &
             li_heat_dissipation_sia, li_basal_friction,  &
             li_temperature_to_enthalpy, li_enthalpy_to_temperature, &
             li_compute_pressure_melting_point_fields


   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   real (kind=RKIND), save :: rhoi    ! ice density (kg m^{-3}), copied from config_ice_density
   real (kind=RKIND), save :: rhoo    ! ocean density (kg m^{-3}), copied from config_ocean_density

   real (kind=RKIND), dimension(:,:), allocatable :: dsigmaTerm   ! vertical grid quantities

   ! max and min allowed temperatures (Kelvin)
   ! Note: kelvin_to_celsius = 273.15 (perhaps it should be called celsius_to_kelvin?)

   real (kind=RKIND), parameter ::   &
        maxtempThreshold =  100._RKIND + kelvin_to_celsius,   & 
        mintempThreshold = -100._RKIND + kelvin_to_celsius     

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
           config_thermal_solver,                 & ! solver option ('temperature' or 'enthalpy')
           config_temperature_init,               & ! temperature initialization option ('linear' or 'file')
           config_surface_air_temperature_source, & ! surface air temperature initialization option ('constant' or 'file') 
           config_basal_heat_flux_source            ! basal heat flux initialization option ('constant' or 'file') 

      real (kind=RKIND), pointer ::  &
           config_ice_density,            & ! ice density
           config_ocean_density             ! ocean density

      real (kind=RKIND), pointer ::  &
           config_surface_air_temperature_value, & ! constant value of surface air temperature
           config_basal_heat_flux_value            ! constant value of basal heat flux
  
      integer, pointer :: &
           config_stats_cell_ID             ! cell ID for diagnostic output

      integer, pointer :: &
           nCellsSolve,              & ! number of locally owned cells
           nVertLevels                 ! number of vertical layers

      integer, dimension(:), pointer :: &
           indexToCellID               ! global ID for each local cell

      real (kind=RKIND), dimension(:), pointer :: &
           layerCenterSigma,         & ! sigma coordinate at midpoint of each layer
           layerInterfaceSigma         ! sigma coordinate at layer interfaces (including top and bottom)

      real (kind=RKIND), dimension(:), pointer :: &
           thickness                   ! ice thickness

      real (kind=RKIND), dimension(:,:), pointer :: &
           temperature,              & ! interior ice temperature (K)
           waterfrac,                & ! interior water fraction (unitless)
           enthalpy                    ! interior ice enthalpy (J m^{-3})

      real (kind=RKIND), dimension(:), pointer :: &
           surfaceAirTemperature,    & ! surface air temperature (K)
           surfaceTemperature,       & ! surface ice temperature (K)
           basalTemperature            ! basal ice temperature (K)

      real (kind=RKIND), dimension(:), pointer :: &
           basalHeatFlux               ! basal heat flux into the ice (W m^{-2}, positive upward)

      integer :: k, iLayer

      integer :: iCell

      integer :: err_tmp

      real(kind=RKIND) :: temperatureValue

      !WHL - debug - for test-case diagnostics
      logical, parameter :: circular_shelf_test = .false.
      logical, parameter :: dome_test = .false.

      integer :: ncellsPerRow 
      integer :: nRows
      integer :: i, iRow
      real(kind=RKIND), parameter :: surfaceAirTemperatureCelsius = -15.0_RKIND


      call mpas_timer_start("thermal init")

      if (circular_shelf_test) then
         nCellsPerRow = 40
         nRows = 46
      elseif (dome_test) then
         nCellsPerRow = 30
         nRows = 34
      endif

      err = 0

      ! get config options
      call mpas_pool_get_config(liConfigs, 'config_thermal_solver', config_thermal_solver)
      call mpas_pool_get_config(liConfigs, 'config_temperature_init', config_temperature_init)
      call mpas_pool_get_config(liConfigs, 'config_surface_air_temperature_source', config_surface_air_temperature_source)
      call mpas_pool_get_config(liConfigs, 'config_surface_air_temperature_value', config_surface_air_temperature_value)
      call mpas_pool_get_config(liConfigs, 'config_basal_heat_flux_source', config_basal_heat_flux_source)
      call mpas_pool_get_config(liConfigs, 'config_basal_heat_flux_value', config_basal_heat_flux_value)
      call mpas_pool_get_config(liConfigs, 'config_print_thermal_info', config_print_thermal_info)
      call mpas_pool_get_config(liConfigs, 'config_do_restart', config_do_restart)
      call mpas_pool_get_config(liConfigs, 'config_stats_cell_ID', config_stats_cell_ID)
            
      ! set some physical constants
      ! (to avoid calling mpas_pool_get_config repeatedly in this module)
      call mpas_pool_get_config(liConfigs, 'config_ice_density', config_ice_density)
      call mpas_pool_get_config(liConfigs, 'config_ocean_density', config_ocean_density)
      rhoi = config_ice_density
      rhoo = config_ocean_density

      ! block loop
      block => domain % blocklist
      do while (associated(block))

         ! get pools
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'thermal', thermalPool)

         ! get dimensions
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

         ! get arrays from the mesh pool
         call mpas_pool_get_array(meshPool, 'layerCenterSigma', layerCenterSigma)
         call mpas_pool_get_array(meshPool, 'layerInterfaceSigma', layerInterfaceSigma)
         call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)

         ! get arrays from the geometry pool
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)

         ! get arrays from the thermal pool
         call mpas_pool_get_array(thermalPool, 'temperature', temperature)
         call mpas_pool_get_array(thermalPool, 'waterfrac', waterfrac)
         call mpas_pool_get_array(thermalPool, 'enthalpy', enthalpy)
         call mpas_pool_get_array(thermalPool, 'surfaceAirTemperature', surfaceAirTemperature)
         call mpas_pool_get_array(thermalPool, 'surfaceTemperature', surfaceTemperature)
         call mpas_pool_get_array(thermalPool, 'basalTemperature',   basalTemperature)
         call mpas_pool_get_array(thermalPool, 'basalHeatFlux', basalHeatFlux)

         if (config_print_thermal_info) then
            write(stdoutUnit,*) 'Initialize thermal solver, config_temperature_init = ', trim(config_temperature_init)
         endif

         ! init for surfaceAirTemperature
         !TODO - Allow for reading surfaceAirTemperature from an external file
         if (trim(config_surface_air_temperature_source) == 'constant') then
            surfaceAirTemperature(:) = config_surface_air_temperature_value
            if (config_print_thermal_info) write(stdoutUnit,*) 'Initialize sfc air temp:', config_surface_air_temperature_value
         endif

         ! init for basalHeatFlux
         !TODO - Allow for reading basalHeatFlux from an external file
         if (trim(config_basal_heat_flux_source) == 'constant') then
            basalHeatFlux(:) = config_basal_heat_flux_value
            if (config_print_thermal_info) write(stdoutUnit,*) 'Initialize basal heat flux:', config_basal_heat_flux_value
         endif

         ! Precompute some grid quantities used in the vertical temperature solve

         allocate(dsigmaTerm(nVertLevels,2))
         dsigmaTerm(:,:) = 0.0_RKIND
      
         k = 1
         dsigmaTerm(k,1) = 1.0_RKIND/( (layerInterfaceSigma(k+1) - layerInterfaceSigma(k)) * (layerCenterSigma(k) - layerInterfaceSigma(k)) )
      
         do k = 2, nVertLevels
            dsigmaTerm(k,1) = 1.0_RKIND/( (layerInterfaceSigma(k+1) - layerInterfaceSigma(k)) * (layerCenterSigma(k) - layerCenterSigma(k-1)) )
         enddo
         
         do k = 1, nVertLevels-1
            dsigmaTerm(k,2) = 1.0_RKIND/( (layerInterfaceSigma(k+1) - layerInterfaceSigma(k)) * (layerCenterSigma(k+1) - layerCenterSigma(k)) )
         end do
         
         k = nVertLevels
         dsigmaTerm(k,2) = 1.0_RKIND/( (layerInterfaceSigma(k+1) - layerInterfaceSigma(k)) * (layerInterfaceSigma(k+1) - layerCenterSigma(k)) )
         
         if (config_print_thermal_info) then
            write(stdoutUnit,*) 'dsigmaTerm coefficients:'
            do k = 1, nVertLevels
               write(stdoutUnit,*) k, dsigmaTerm(k,1), dsigmaTerm(k,2)
            enddo
         endif
         
         !WHL - debug - Temporary initialization for test cases

         if (config_print_thermal_info .and. (circular_shelf_test .or. dome_test)) then

            write(stdoutUnit,*) ' '

            if (circular_shelf_test) then
               
               write(stdoutUnit,*) 'Circular shelf test:'

            elseif (dome_test) then

               write(stdoutUnit,*) 'Dome test:'
!                  write(stdoutUnit,*) ' '
!                  write(stdoutUnit,*) 'Cells with ice: iCell, thickness'
!                  do iCell = 1, nCellsSolve
!                     if (thickness(iCell) > 0.0_RKIND) then
!                        write(stdoutUnit,*) iCell, thickness(iCell)
!                     endif
!                  enddo
                  
            endif

            !WHL - debug dome - Adjust thickness slightly so the central cell has the same thickness as in CISM
!!            thickness(:) = thickness(:) * 707.1068115234375000_RKIND / 707.10678118654755_RKIND
            
            write(stdoutUnit,*) ' '
            write(stdoutUnit,*) 'Initial thickness'
            do iRow = nRows, 1, -1
               if (mod(iRow,2) == 0) then  ! indent for even-numbered rows
                  write(stdoutUnit,'(a3)',advance='no') '    '
               endif
               do i = nCellsPerRow/2 - 2, nCellsPerRow
                  iCell = (iRow-1)*nCellsPerRow + i
                  write(stdoutUnit,'(f8.2)',advance='no') thickness(iCell)
               enddo
               write(stdoutUnit,*) ' '
            enddo
            
            write(stdoutUnit,*) ' '
            write(stdoutUnit,*) 'Initial surfaceAirTemperature'
            do iRow = nRows, 1, -1
               if (mod(iRow,2) == 0) then  ! indent for even-numbered rows
                  write(stdoutUnit,'(a3)',advance='no') '    '
               endif
               do i = nCellsPerRow/2 - 2, nCellsPerRow
                  iCell = (iRow-1)*nCellsPerRow + i
                  write(stdoutUnit,'(f8.2)',advance='no') surfaceAirTemperature(iCell)
               enddo
               write(stdoutUnit,*) ' '
            enddo
            
         endif  ! circ shelf or dome


         ! Initialize vertical temperature profile.
         ! Three possibilities:
         ! (1) Set up a linear temperature profile, with T = artm at the surface and T <= Tpmp
         !     at the bed (config_temperature_init = 'linear'). 
         !     A parameter (pmpt_offset) controls how far below Tpmp the initial bed temp is set.
         ! (2) Read ice temperature from an initial input file (config_temperature_init = 'file').
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
               write(stdoutUnit,*) 'Initialized ice temperature from the restart file'
            endif

         elseif (trim(config_temperature_init) == 'file') then

            ! Temperature was read from the input file
            if (config_print_thermal_info) then
               write(stdoutUnit,*) 'Initialized ice temperature from the input file'
            endif

            ! initialize waterfrac, in case we are using the enthalpy solver
            !TODO - Allow waterfrac to be read from the input file?
            waterfrac(:,:) = 0.0_RKIND

         elseif (trim(config_temperature_init) == 'sfc_air_temperature') then

            do iCell = 1, nCellsSolve
               ! set the column temperature to the surface air temperature or 273.15, whichever is less
               temperatureValue = min(surfaceAirTemperature(iCell), kelvin_to_celsius)
               surfaceTemperature(iCell) = temperatureValue
               temperature(:,iCell)      = temperatureValue
               basalTemperature(iCell)   = temperatureValue
            enddo

            waterfrac(:,:) = 0.0_RKIND

            if (config_print_thermal_info) then
               write(stdoutUnit,*) 'Initialized ice column temperature to the surface air temperature'
            endif

         elseif (trim(config_temperature_init) == 'linear') then

            ! set up a linear temperature profile in each column
            ! T = surfaceAirTemperature at the ice surface, and T <= Tpmp at the bed

            ! initialize T = 273.15 K = 0 C everywhere

            temperature(:,:) = kelvin_to_celsius    ! = 273.15                                              

            do iCell = 1, nCellsSolve

               call li_init_linear_temperature_in_column(&
                    nVertLevels,                   &
                    layerCenterSigma,              &
                    thickness(iCell),              &
                    surfaceAirTemperature(iCell),  &
                    temperature(:,iCell),          &
                    waterfrac(:,iCell),            &
                    surfaceTemperature(iCell),     &
                    basalTemperature(iCell))

            enddo  ! iCell

            if (config_print_thermal_info) then
               write(stdoutUnit,*) 'Initialized a linear temperature profile in each column'

               !WHL - debug
               if (circular_shelf_test .or. dome_test) then

                  write(stdoutUnit,*) ' '
                  write(stdoutUnit,*) 'Initial layer 1 temperature'
                  iLayer = 1
                  do iRow = nRows, 1, -1
                     if (mod(iRow,2) == 0) then  ! indent for even-numbered rows
                        write(stdoutUnit,'(a3)',advance='no') '    '
                     endif
                     do i = nCellsPerRow/2 - 2, nCellsPerRow
                        iCell = (iRow-1)*nCellsPerRow + i
                        write(stdoutUnit,'(f8.2)',advance='no') temperature(iLayer,iCell)
                     enddo
                     write(stdoutUnit,*) ' '
                  enddo

                  write(stdoutUnit,*) ' '
                  write(stdoutUnit,*) 'Initial basal temperature'
                  do iRow = nRows, 1, -1
                     if (mod(iRow,2) == 0) then  ! indent for even-numbered rows
                        write(stdoutUnit,'(a3)',advance='no') '    '
                     endif
                     do i = nCellsPerRow/2 - 2, nCellsPerRow
                        iCell = (iRow-1)*nCellsPerRow + i
                        write(stdoutUnit,'(f8.2)',advance='no') basalTemperature(iCell)
                     enddo
                     write(stdoutUnit,*) ' '
                  enddo

                  do iCell = 1, nCellsSolve
                     if (indexToCellID(iCell) == config_stats_cell_ID) then
                        write(stdoutUnit,*) ' '
                        write(stdoutUnit,*) 'Diagnostic cell =', iCell
                        write(stdoutUnit,*) 'Thickness =', thickness(iCell)
                        write(stdoutUnit,*) ' '
                        write(stdoutUnit,*) 'Temperature profile (C):'
                        write(stdoutUnit,*) 'Sfc:', surfaceTemperature(iCell) - kelvin_to_celsius
                        do iLayer = 1, nVertLevels
                           write(stdoutUnit,*) iLayer, layerCenterSigma(iLayer), temperature(iLayer,iCell) - kelvin_to_celsius
                        enddo
                        write(stdoutUnit,*) 'Bed:', basalTemperature(iCell) - kelvin_to_celsius
                     endif
                  enddo

               endif  ! shelf or dome test

            endif  ! config_print_thermal_info

         endif    ! restart file, input file, or linear

         block => block % next
      enddo

      !TODO - Add a debug check for bad values
      !       E.g., make sure the temperature read from a file is in Kelvin and not Celsius

      ! halo updates
      call mpas_timer_start("halo updates")

      call mpas_dmpar_field_halo_exch(domain, 'surfaceTemperature')
      call mpas_dmpar_field_halo_exch(domain, 'basalTemperature')
      call mpas_dmpar_field_halo_exch(domain, 'temperature')

      if (trim(config_thermal_solver) == 'enthalpy') then
         ! prognostic variables are temperature and waterfrac, so need to update waterfrac too
         call mpas_dmpar_field_halo_exch(domain, 'waterfrac')
      endif

      call mpas_timer_stop("halo updates")

      ! === error check
      if (err > 0) then
          write (stderrUnit,*) "An error has occurred in li_thermal_init."
      endif

      call mpas_timer_stop("thermal init")

   !--------------------------------------------------------------------
    end subroutine li_thermal_init

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  !  li_init_linear_temperature_in_column
!
!> \brief MPAS initialize a linear temperature profile in a column
!> \author William Lipscomb
!> \date   February 2016
!> \details
!>  This routine initializes a linear temperature profile in an ice column,
!>  with the surface ice temperature set to the surface air temperature,
!>  and the basal ice temperature set to the pressure melting point minus
!>  an offset.
!-----------------------------------------------------------------------

    subroutine li_init_linear_temperature_in_column(&
         nVertLevels,            &
         layerCenterSigma,       &
         thickness,              &
         surfaceAirTemperature,  &
         temperature,            &
         waterfrac,              &
         surfaceTemperature,     &
         basalTemperature)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      integer, intent(in) :: &
           nVertLevels                  !< Input: number of ice layers

      real(kind=RKIND), dimension(nVertLevels), intent(in) :: &
           layerCenterSigma             !< Input: vertical sigma coordinate at layer midpoints

      real(kind=RKIND), intent(in) ::  &
           thickness,    &              !< Input: ice thickness
           surfaceAirTemperature        !< Input: surface air temperature

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real(kind=RKIND), dimension(nVertLevels), intent(out) ::  &
           temperature,     &           !< Output: interior ice temperature at midpoint of each layer
           waterfrac                    !< Output: interior water fraction at midpoint of each layer

      real(kind=RKIND), intent(out) ::  &
           surfaceTemperature,     &    !< Output: surface ice temperature
           basalTemperature             !< Output: basal ice temperature

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(nVertLevels) :: &
           pmptemp                     ! pressure melting point temp in ice interior

      real (kind=RKIND) :: &
           pmptemp_bed                 ! pressure melting point temp at bed

      integer :: nVert

      !TODO - Make pmpt_offset a namelist parameter?
      real (kind=RKIND), parameter :: &
           pmpt_offset = 2.0_RKIND     ! offset of initial Tbed from pressure melting point temperature (K)
                                       ! Note: pmtp_offset is positive for T < Tpmp


      ! set the surface ice temperature to the surface air temperature (or 273.15, whichever is less)
      surfaceTemperature = min(surfaceAirTemperature, kelvin_to_celsius)

      ! compute the pressure melting point temperature in the column and at the bed
      ! (in degrees C)

      call pressure_melting_point_column(&
           layerCenterSigma,   &
           thickness,          &
           pmptemp)

      call pressure_melting_point(&
           thickness,          &
           pmptemp_bed)

      ! convert to Kelvin
      pmptemp(:) = pmptemp(:) + kelvin_to_celsius
      pmptemp_bed = pmptemp_bed + kelvin_to_celsius

      ! set the basal temperature to slightly below the pressure melting point temperature
      basalTemperature = pmptemp_bed - pmpt_offset

      ! set the interior temperatures
      ! make sure T <= Tpmp - pmpt_offset in column interior

      temperature(:) = surfaceTemperature +  &
           (basalTemperature - surfaceTemperature) * layerCenterSigma(:)

      temperature(:) = min(temperature(:), pmptemp(:) - pmpt_offset)

      ! set waterfrac = 0
      waterfrac(:) = 0.0_RKIND

    end subroutine li_init_linear_temperature_in_column

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

   subroutine li_thermal_solver(domain, err)

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
      type (mpas_pool_type), pointer :: velocityPool   ! needed for mask subroutine
      type (mpas_pool_type), pointer :: scratchPool 

      integer, pointer :: &
           nCellsSolve,              & ! number of locally owned cells
           nVertLevels                 ! number of vertical layers

      logical, pointer :: &
           config_print_thermal_info   ! if true, print debug info

      character(len=StrKIND), pointer :: &
           config_thermal_solver,    & ! option for thermal solver
           config_basal_mass_bal_float ! option for basal mass balance of floating ice

      real(kind=RKIND), pointer :: &
           config_thermal_thickness, & ! minimum thickness (m) for temperature calculations
           config_sea_level,         & ! sea level (m) relative to z = 0
           config_bmlt_float_flux,   & ! constant heat flux (W/m^2) applied to the base of floating ice; positive upward
           config_bmlt_float_xlimit    ! x value (m) defining region where bmlt_float_flux is applied; melt only where abs(x) > xlimit

      integer, pointer :: &
           config_stats_cell_ID        ! global ID for diagnostic cell

      integer, dimension(:), pointer :: &
           cellMask,                 & ! bit mask describing whether ice is floating, dynamically active, etc.
           indexToCellID               ! list of global cell IDs
      
      type (field1dInteger), pointer :: thermalCellMaskField
      integer, dimension(:), pointer :: &
           thermalCellMask             ! mask for thermal calculations
                                       ! = 1 where thickness > config_thermal_thickness, elsewhere = 0 

      real (kind=RKIND), pointer :: deltat  !< time step in seconds

      real (kind=RKIND), dimension(:), pointer :: &
           xCell                       ! x coordinate for each cell (m)

      real (kind=RKIND), dimension(:), pointer :: &
           layerCenterSigma,         & ! sigma coordinate at midpoint of each layer
           layerInterfaceSigma         ! sigma coordinate at layer interfaces (including top and bottom)

      real (kind=RKIND), dimension(:), pointer :: &
           surfaceTemperature,       & ! surface ice temperature (K)
           basalTemperature,         & ! basal ice temperature (K)
           surfaceAirTemperature,    & ! surface air temperature (K)
           basalHeatFlux,            & ! basal heat flux into the ice (W m^{-2}, positive upward)
           basalFrictionFlux,        & ! basal frictional flux into the ice (W m^{-2})
           surfaceConductiveFlux,    & ! conductive heat flux at the upper surface (W m^{-2}, positive down) 
           basalConductiveFlux,      & ! conductive heat flux at the lower surface (W m^{-2}, positive down) 
           basalMassBal,             & ! basal mass balance (kg m^{-2} s^{-1}); positive for freeze-on, negative for melting
           groundedBasalMassBal,     & ! basal mass balance for grounded ice
           floatingBasalMassBal,     & ! basal mass balance for floating ice
           basalWaterThickness,      & ! basal water thickness (m)
           thickness,                & ! ice thickness (m)
           lowerSurface,             & ! lower surface elevation (m)
           bedTopography               ! bed topography (m; negative below sea level)

      real (kind=RKIND), dimension(:,:), pointer :: &
           temperature,              & ! interior ice temperature (K)
           waterfrac,                & ! interior water fraction (unitless)
           enthalpy,                 & ! interior ice enthalpy (J m^{-3})
           heatDissipation             ! interior heat dissipation (deg/s)

      real(kind=RKIND), dimension(:), allocatable :: &
           subdiagonal, diagonal, superdiagonal,   &  ! tridiagonal matrix elements
           rhs                                        ! matrix right-hand side

      real(kind=RKIND), dimension(:), allocatable :: &
           diffusivity   ! diffusivity at interfaces (m^2/s) for enthalpy solver
                         ! = iceConductivity / (rhoi*cp_ice) for cold ice

      real(kind=RKIND), dimension(:), allocatable :: &
           solution     ! solution of tridiagonal matrix problem

      real(kind=RKIND) :: &
           surfaceEnthalpy,          & ! surface ice enthalpy
           basalEnthalpy,            & ! basal ice enthalpy
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

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_thermal_solver', config_thermal_solver)
      if (trim(config_thermal_solver) == 'none') then
         return ! nothing to do
      endif

      ! get rest of config variables
      call mpas_pool_get_config(liConfigs, 'config_print_thermal_info', config_print_thermal_info)
      call mpas_pool_get_config(liConfigs, 'config_thermal_thickness', config_thermal_thickness)
      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
      call mpas_pool_get_config(liConfigs, 'config_basal_mass_bal_float', config_basal_mass_bal_float)
      call mpas_pool_get_config(liConfigs, 'config_stats_cell_ID', config_stats_cell_ID)


      if (config_print_thermal_info) then
         write(stdoutUnit,*) 'Solving for temperature, config_thermal_solver = ', trim(config_thermal_solver)
      endif


      ! check that deltat is valid
      ! Note: With an adaptive time step, deltat is set to 0 at the start of the time step.
      !       Here we make sure the thermal solver is not called while deltat = 0.

      ! Get dt - same on all blocks so just grab first one
      block => domain % blocklist
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_array(meshPool, 'deltat', deltat)

      if (deltat <= 0.0_RKIND) then
         write(stderrUnit,*) 'ERROR: li_thermal_solver was called with invalid deltat =', deltat
         call mpas_dmpar_global_abort("MPAS-landice: An error has occurred in li_thermal. Aborting...")
      endif



      ! block loop
      block => domain % blocklist
      do while (associated(block))

         ! get pools
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'thermal', thermalPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

         ! get dimensions
         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

         ! get fields from the mesh pool
         call mpas_pool_get_array(meshPool, 'xCell', xCell)
         call mpas_pool_get_array(meshPool, 'layerCenterSigma', layerCenterSigma)
         call mpas_pool_get_array(meshPool, 'layerInterfaceSigma', layerInterfaceSigma)
         call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)  ! diagnostic only

         ! get fields from the geometry pool
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(geometryPool, 'lowerSurface', lowerSurface)
         call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_array(geometryPool, 'basalMassBal', basalMassBal)
         call mpas_pool_get_array(geometryPool, 'groundedBasalMassBal', groundedBasalMassBal)
         call mpas_pool_get_array(geometryPool, 'floatingBasalMassBal', floatingBasalMassBal)
         call mpas_pool_get_array(geometryPool, 'basalWaterThickness', basalWaterThickness)

         ! get fields from the thermal pool
         call mpas_pool_get_array(thermalPool, 'temperature', temperature)
         call mpas_pool_get_array(thermalPool, 'waterfrac', waterfrac)
         call mpas_pool_get_array(thermalPool, 'enthalpy', enthalpy)
         call mpas_pool_get_array(thermalPool, 'surfaceTemperature', surfaceTemperature)
         call mpas_pool_get_array(thermalPool, 'basalTemperature',   basalTemperature)
         call mpas_pool_get_array(thermalPool, 'surfaceAirTemperature', surfaceAirTemperature)
         call mpas_pool_get_array(thermalPool, 'surfaceConductiveFlux', surfaceConductiveFlux)
         call mpas_pool_get_array(thermalPool, 'basalConductiveFlux', basalConductiveFlux)
         call mpas_pool_get_array(thermalPool, 'basalHeatFlux', basalHeatFlux)
         call mpas_pool_get_array(thermalPool, 'basalFrictionFlux', basalFrictionFlux)
         call mpas_pool_get_array(thermalPool, 'heatDissipation', heatDissipation)

         ! get fields from the scratch pool
         call mpas_pool_get_field(scratchPool, 'iceCellMask', thermalCellMaskField)
         call mpas_allocate_scratch_field(thermalCellMaskField, .true.)
         thermalCellMask => thermalCellMaskField % array

         ! get config variables
         call mpas_pool_get_config(liConfigs, 'config_print_thermal_info', config_print_thermal_info)
         call mpas_pool_get_config(liConfigs, 'config_thermal_solver', config_thermal_solver)
         call mpas_pool_get_config(liConfigs, 'config_thermal_thickness', config_thermal_thickness)
         call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
         call mpas_pool_get_config(liConfigs, 'config_basal_mass_bal_float', config_basal_mass_bal_float)
         call mpas_pool_get_config(liConfigs, 'config_bmlt_float_flux', config_bmlt_float_flux)
         call mpas_pool_get_config(liConfigs, 'config_bmlt_float_xlimit', config_bmlt_float_xlimit)
         call mpas_pool_get_config(liConfigs, 'config_stats_cell_ID', config_stats_cell_ID)

         if (config_print_thermal_info) then
            write(stdoutUnit,*) 'Solving for temperature, config_thermal_solver = ', trim(config_thermal_solver)
         endif

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

         !WHL - debug dome - Zero out heat dissipation for easier comparison with dome problem
!!         heatDissipation(:,:) = 0.0d0

         select case(config_thermal_solver)

         case ('temperature', 'enthalpy')

            ! Convert temperature from Kelvin to Celsius to avoid repeated use of kelvin_to_celsius below
            ! (Convert back at the end.)
            temperature(:,:) = temperature(:,:) - kelvin_to_celsius
            surfaceAirTemperature(:) = surfaceAirTemperature(:) - kelvin_to_celsius
            surfaceTemperature(:) = surfaceTemperature(:) - kelvin_to_celsius
            basalTemperature(:) = basalTemperature(:) - kelvin_to_celsius

            ! allocate some vertical arrays
            allocate(subdiagonal(nVertLevels+2))  ! temperature/enthalpy in each layer, plus surface and basal temperature
            allocate(diagonal(nVertLevels+2))
            allocate(superdiagonal(nVertLevels+2))
            allocate(rhs(nVertLevels+2))
            allocate(solution(nVertLevels+2))
            allocate(diffusivity(nVertLevels+1))

            if (config_print_thermal_info) then
               write(stdoutUnit,*) ' '
               do iCell = 1, nCellsSolve
                  if (indexToCellID(iCell) == config_stats_cell_ID) then
                     write(stdoutUnit,*) 'thickness =', thickness(iCell)
                     write(stdoutUnit,*) 'surfaceAirTemperature =', surfaceAirTemperature(iCell)
                     write(stdoutUnit,*) ' '
                     write(stdoutUnit,*) 'Initial column temperatures, iCell =', iCell
                     write(stdoutUnit,*) 0, surfaceTemperature(iCell)
                     do k = 1, nVertLevels
                        write(stdoutUnit,*) k, temperature(k,iCell)
                     enddo
                     write(stdoutUnit,*) nVertLevels+1, basalTemperature(iCell)
                     write(stdoutUnit,*) 'Loop over cells'
                  endif
               enddo
            endif

            ! loop over locally owned cells
            do iCell = 1, nCellsSolve
               
               if (config_print_thermal_info .and. indexToCellID(iCell) == config_stats_cell_ID) then
                  verboseColumn = .true.
               else
                  verboseColumn = .false.
               endif

               if (thermalCellMask(iCell) == 1) then  ! thermally active ice is present

                  ! Set surface temperature (Celsius)
                  
                  surfaceTemperature(iCell) = min(0.0_RKIND, surfaceAirTemperature(iCell))

                  ! For floating ice, set the basal temperature to the freezing temperature of seawater.
                  ! Values based on Ocean Water Freezing Point Calculator with S = 35 PSU 
                  if (li_mask_is_floating_ice(cellMask(iCell))) then
                     depth = thickness(iCell) * rhoi/rhoo
                     basalTemperature(iCell) = oceanFreezingTempSurface + oceanFreezingTempDepthDependence * depth  ! Celsius
                  endif
                  
                  if (trim(config_thermal_solver) == 'enthalpy') then

                     ! Given temperature and waterfrac in ice interior, compute enthalpy

                     !WHL - debug - nonzero waterfrac for testing
!!                     if (verboseColumn) then
!!                        waterfrac(nVertLevels-1:nVertLevels,iCell) = 0.01d0
!!                     endif

                     call li_temperature_to_enthalpy(&
                          layerCenterSigma,      &
                          thickness(iCell),      &
                          temperature(:,iCell),  &
                          waterfrac(:,iCell),    &
                          enthalpy(:,iCell))

                     surfaceEnthalpy = surfaceTemperature(iCell) * rhoi*cp_ice
                     basalEnthalpy = basalTemperature(iCell) * rhoi*cp_ice

                     if (verboseColumn) then
                        write(stdoutUnit,*) ' '
                        write(stdoutUnit,*) 'Before prognostic enthalpy, iCell =', indexToCellID(iCell)
                        write(stdoutUnit,*) 'thickness =', thickness(iCell)
                        write(stdoutUnit,*) 'Temperature (C), waterfrac, enthalpy/(rhoi*cp_ice):'
                        write(stdoutUnit,*) 0, surfaceEnthalpy/(rhoi*cp_ice)
                        do k = 1, nVertLevels
                           write(stdoutUnit,*) k, temperature(k,iCell), waterfrac(k,iCell), enthalpy(k,iCell)/(rhoi*cp_ice)
                        enddo
                        write(stdoutUnit,*) nVertLevels+1, basalEnthalpy/(rhoi*cp_ice)
                     endif

                     ! Compute initial internal energy in column (for energy conservation check)
                     initialEnergy = 0.0_RKIND
                     do k = 1, nVertLevels
                        initialEnergy = initialEnergy + enthalpy(k,iCell) * (layerInterfaceSigma(k+1) - layerInterfaceSigma(k))
                     enddo
                     initialEnergy = initialEnergy * thickness(iCell)

                     ! Compute matrix elements using enthalpy gradient method

                     call enthalpy_matrix_elements(&
                          deltat,                         &
                          nVertLevels,                    &
                          layerCenterSigma,               &
                          dsigmaTerm,                     &
                          li_mask_is_floating_ice_int(cellMask(iCell)), &
                          thickness(iCell),               &
                          temperature(:,iCell),           &
                          surfaceTemperature(iCell),      &
                          basalTemperature(iCell),        &
                          waterfrac(:,iCell),             &
                          enthalpy(:,iCell),              &
                          heatDissipation(:,iCell),       &
                          basalHeatFlux(iCell),           &
                          basalFrictionFlux(iCell),       &
                          diffusivity,                    &
                          subdiagonal,                    &
                          diagonal,                       &
                          superdiagonal,                  &
                          rhs)

                     if (verboseColumn) then
                        write(stdoutUnit,*) ' '
                        write(stdoutUnit,*) 'After matrix elements, iCell =', indexToCellID(iCell)
                        write(stdoutUnit,*) 'k, subd, diag, supd, rhs/(rhoi*ci):'
                        do k = 1, nVertLevels+2
                           write(stdoutUnit,*) k-1, subdiagonal(k), diagonal(k), superdiagonal(k), rhs(k)/(rhoi*cp_ice)
                        enddo
                     endif

                     ! solve the tridiagonal system
                     ! Note: Temperature is indexed from 1 to nVertLevels, whereas the matrix elements
                     !        are indexed from 1 to nVertLevels+2.
                     !       Matrix row 1 corresponds to surface temperature, and matrix row nVertLevels+2
                     !        corresponds to the basal temperature.

                     call tridiag_solver(&
                          subdiagonal,   &
                          diagonal,      &
                          superdiagonal, &
                          solution,      &
                          rhs)

                     ! Copy the solution into the enthalpy variables
                     surfaceEnthalpy   = solution(1)
                     enthalpy(:,iCell) = solution(2:nVertLevels+1)
                     basalEnthalpy     = solution(nVertLevels+2)

                     ! Compute conductive fluxes = (diffusivity/thickness * denth/dsigma) at upper and lower surfaces; positive down.
                     ! Here diffusivity = iceConductivity / (rhoi*cp_ice) for cold ice, with a smaller value for temperate ice.
                     ! Assume implicit backward Euler time step.
                     ! Note: These fluxes should be computed before calling glissade_enth2temp (which might change basalEnthalpy).

                     denth_top = enthalpy(1,iCell) - surfaceEnthalpy
                     denth_bot = basalEnthalpy - enthalpy(nVertLevels,iCell)
                
                     surfaceConductiveFlux(iCell) = -diffusivity(1)/thickness(iCell) * denth_top/layerCenterSigma(1)
                     basalConductiveFlux(iCell) = -diffusivity(nVertLevels+1)/thickness(iCell) * denth_bot/(1.0_RKIND - layerCenterSigma(nVertLevels))
                                              
                     ! convert enthalpy in ice interior back to temperature and waterfrac

                     call li_enthalpy_to_temperature(&
                          layerCenterSigma,          &
                          thickness(iCell),          &
                          enthalpy(:,iCell),         &
                          temperature(:,iCell),      &
                          waterfrac(:,iCell),        &
                          surfaceEnthalpy,           &
                          surfaceTemperature(iCell), &
                          basalEnthalpy,             &
                          basalTemperature(iCell))

                     if (verboseColumn) then
                        write(stdoutUnit,*) ' '
                        write(stdoutUnit,*) 'After prognostic enthalpy, iCell =', indexToCellID(iCell)
                        write(stdoutUnit,*) 'thickness =', thickness(iCell)
                        write(stdoutUnit,*) 'Temperature, waterfrac, enthalpy/(rhoi*cp_ice):'
                        write(stdoutUnit,*) 0, surfaceEnthalpy/(rhoi*cp_ice)
                        do k = 1, nVertLevels
                           write(stdoutUnit,*) k, temperature(k,iCell), waterfrac(k,iCell), enthalpy(k,iCell)/(rhoi*cp_ice)
                        enddo
                        write(stdoutUnit,*) nVertLevels+1, basalEnthalpy/(rhoi*cp_ice)
                     endif

                     ! Compute final internal energy in column (for energy conservation check)
                     finalEnergy = 0.0_RKIND
                     do k = 1, nVertLevels
                        finalEnergy = finalEnergy + enthalpy(k,iCell) * (layerInterfaceSigma(k+1) - layerInterfaceSigma(k))
                     enddo
                     finalEnergy = finalEnergy * thickness(iCell)

                  else    ! temperature solver

                     if (verboseColumn) then
                        write(stdoutUnit,*) ' '
                        write(stdoutUnit,*) 'Before prognostic temperature, iCell =', indexToCellID(iCell)
                        write(stdoutUnit,*) 'thickness =', thickness(iCell)
                        write(stdoutUnit,*) 0, surfaceTemperature(iCell)
                        do k = 1, nVertLevels
                           write(stdoutUnit,*) k, temperature(k,iCell)
                        enddo
                        write(stdoutUnit,*) nVertLevels+1, basalTemperature(iCell)
                     endif

                     ! Compute initial internal energy in column (for energy conservation check)
                     initialEnergy = 0.0_RKIND
                     do k = 1, nVertLevels
                        initialEnergy = initialEnergy + temperature(k,iCell) * (layerInterfaceSigma(k+1) - layerInterfaceSigma(k))
                     enddo
                     initialEnergy = initialEnergy * thickness(iCell) * rhoi*cp_ice

                     ! Compute matrix elements

                     call temperature_matrix_elements(&
                          deltat,                         &
                          nVertLevels,                    &
                          layerCenterSigma,               &
                          dsigmaTerm,                     &
                          li_mask_is_floating_ice_int(cellMask(iCell)), &
                          thickness(iCell),               &
                          temperature(:,iCell),           &
                          surfaceTemperature(iCell),      &
                          basalTemperature(iCell),        &
                          heatDissipation(:,iCell),       &
                          basalHeatFlux(iCell),           &
                          basalFrictionFlux(iCell),       &
                          subdiagonal,                    &
                          diagonal,                       &
                          superdiagonal,                  &
                          rhs)

                     if (verboseColumn) then
                        write(stdoutUnit,*) ' '
                        write(stdoutUnit,*) 'deltat =', deltat
                        write(stdoutUnit,*) 'After matrix elements, iCell =', indexToCellID(iCell)
                        write(stdoutUnit,*) 'k, subd, diag, supd, rhs:'
                        do k = 1, nVertLevels+2
                           write(stdoutUnit,*) k-1, subdiagonal(k), diagonal(k), superdiagonal(k), rhs(k)
                        enddo
                     endif

                     ! Solve the tridiagonal system
                     ! Note: Temperature is indexed from 1 to nVertLevels, whereas the matrix elements
                     !        are indexed from 1 to nVertLevels+2.
                     !       Matrix row 1 corresponds to surface temperature, and matrix row nVertLevels+2
                     !        corresponds to the basal temperature.

                     call tridiag_solver(&
                          subdiagonal,   &
                          diagonal,      &
                          superdiagonal, &
                          solution,      &
                          rhs)

                     ! Copy the solution into the temperature variables
                     surfaceTemperature(iCell) = solution(1)
                     temperature(:,iCell)      = solution(2:nVertLevels+1)
                     basalTemperature(iCell)   = solution(nVertLevels+2)

                     ! Compute conductive flux = (k/H * dT/dsigma) at upper and lower surfaces; positive down
                     ! Assume implicit backward Euler time step.
                     
                     dTtop = temperature(1,iCell) - surfaceTemperature(iCell)
                     dTbot = basalTemperature(iCell) - temperature(nVertLevels,iCell)

                     surfaceConductiveFlux(iCell) = (-iceConductivity/thickness(iCell) ) * dTtop / layerCenterSigma(1)
                     basalConductiveFlux(iCell) = (-iceConductivity/thickness(iCell) ) * dTbot / (1.0_RKIND - layerCenterSigma(nVertLevels))

                     if (verboseColumn) then
                        write(stdoutUnit,*) ' '
                        write(stdoutUnit,*) 'After prognostic temperature, iCell =', indexToCellID(iCell)
                        write(stdoutUnit,*) 'thickness =', thickness(iCell)
                        write(stdoutUnit,*) 0, surfaceTemperature(iCell)
                        do k = 1, nVertLevels
                           write(stdoutUnit,*) k, temperature(k,iCell)
                        enddo
                        write(stdoutUnit,*) nVertLevels+1, basalTemperature(iCell)
                     endif

                     ! Compute final internal energy in column (for energy conservation check)
                     finalEnergy = 0.0_RKIND
                     do k = 1, nVertLevels
                        finalEnergy = finalEnergy + temperature(k,iCell) * (layerInterfaceSigma(k+1) - layerInterfaceSigma(k))
                     enddo
                     finalEnergy = finalEnergy * thickness(iCell) * rhoi*cp_ice

                  endif   ! temperature or enthalpy solver

                  ! Compute total dissipation rate in column (W/m^2)                                                                                                                  
                  columnHeatDissipation = 0.0_RKIND
                  do k = 1, nVertLevels
                     columnHeatDissipation = columnHeatDissipation  &
                                           + heatDissipation(k,iCell) * (layerInterfaceSigma(k+1) - layerInterfaceSigma(k))
                  enddo
                  columnHeatDissipation = columnHeatDissipation * thickness(iCell)*rhoi*cp_ice

                  if (verboseColumn) then
                     write(stdoutUnit,*) ' '
                     write(stdoutUnit,*) 'heatDissipation:', heatDissipation(:,iCell)
                     write(stdoutUnit,*) 'columnHeatDissipation:', columnHeatDissipation
                  endif

                  ! Verify that the net input of energy into the column is equal to the change in internal energy.

                  deltaEnergy = (surfaceConductiveFlux(iCell) - basalConductiveFlux(iCell) + columnHeatDissipation) * deltat

                  !TODO - Confirm that this is a reasonable error threshold
                  if (abs((finalEnergy - initialEnergy - deltaEnergy) / deltat) > 1.0e-8_RKIND) then

                     if (verboseColumn) then
                        write(stdoutUnit,*) 'Ice thickness:', thickness(iCell)
                        write(stdoutUnit,*) 'config_thermal_thickness:', config_thermal_thickness
                        write(stdoutUnit,*) ' '
                        write(stdoutUnit,*) 'Interior fluxes:'
                        write(stdoutUnit,*) 'sfc conductive flx (W/m^2, positive down)=', surfaceConductiveFlux(iCell)
                        write(stdoutUnit,*) 'bed conductive flx (W/m^2, positive down)=', basalConductiveFlux(iCell)
                        write(stdoutUnit,*) 'column heat dissipation (W/m^2) =', columnHeatDissipation
                        write(stdoutUnit,*) 'Net flux (W/m^2) =', deltaEnergy/deltat
                        write(stdoutUnit,*) ' '
                        write(stdoutUnit,*) 'deltaEnergy (J/m^2) =', deltaEnergy
                        write(stdoutUnit,*) 'initialEnergy (J/m^2) =', initialEnergy
                        write(stdoutUnit,*) 'finalEnergy (J/m^2) =', finalEnergy
                        write(stdoutUnit,*) ' '
                        write(stdoutUnit,*) 'Energy imbalance (J/m^2)=', finalEnergy - initialEnergy - deltaEnergy
                        write(stdoutUnit,*) ' '
                        write(stdoutUnit,*) 'Basal fluxes:'
                        write(stdoutUnit,*) 'frictional =', basalFrictionFlux(iCell)
                        write(stdoutUnit,*) 'geothermal =', basalHeatFlux(iCell)
                        write(stdoutUnit,*) 'flux for bottom melting =', basalFrictionFlux(iCell) + basalHeatFlux(iCell) + basalConductiveFlux(iCell)
                     endif   ! verboseColumn

                     write(stderrUnit,*) 'li_thermal, energy conservation error: iCell, imbalance (W/m2):', &
                                         indexToCellID(iCell), (finalEnergy - initialEnergy - deltaEnergy)/deltat
                     err = 1
                     
                  endif  ! energy conservation error

               else    ! thermalCellMask = 0; ice is not thermally active

                  ! Set temperature of thin ice to 0 C
                  !TODO - For cells that have just crossed the config_thermal_thickness threshold, energy is not conserved here.
                  !       Keep track of the energy difference?

                  surfaceTemperature(iCell) = 0.0_RKIND
                  basalTemperature(iCell) = 0.0_RKIND
                  temperature(:,iCell) = 0.0_RKIND
                  waterfrac(:,iCell) = 0.0_RKIND
                  enthalpy(:,iCell) = 0.0_RKIND

               endif   ! thickness > config_thermal_thickness
  
            enddo   ! iCell

            ! Compute basal melt rate for grounded ice.
            ! Note:
            ! * This subroutine needs to be called only if the temperature/enthalpy is prognostic.
            ! * It assumes Celsius units, so it should be called before converting back to Kelvin.
            ! * It includes internal melting.
            !   For the standard temperature scheme, temperatures above the pressure melting point
            !    are reset to Tpmp, with excess heat contributing to basal melt.
            !   For the enthalpy scheme, internal meltwater in excess of the prescribed maximum
            !    fraction (0.01 by default) is drained to the bed.

            call basal_melt_grounded_ice(&
                 config_thermal_solver,        &     
                 deltat,                       &
                 nCellsSolve,                  &
                 nVertLevels,                  &
                 layerInterfaceSigma,          &
                 layerCenterSigma,             &
                 thermalCellMask,              &
                 li_mask_is_floating_ice_int(cellMask), &
                 thickness,                    &
                 temperature,                  &
                 basalTemperature,             &
                 waterfrac,                    &
                 enthalpy,                     &
                 basalFrictionFlux,            &
                 basalHeatFlux,                &
                 basalConductiveFlux,          &
                 basalWaterThickness,          &
                 groundedBasalMassBal)

            ! Convert temperatures from Celsius back to Kelvin
            temperature(:,:) = temperature(:,:) + kelvin_to_celsius
            surfaceAirTemperature(:) = surfaceAirTemperature(:) + kelvin_to_celsius
            surfaceTemperature(:) = surfaceTemperature(:) + kelvin_to_celsius
            basalTemperature(:) = basalTemperature(:) + kelvin_to_celsius

            ! Check for temperatures that are physically unrealistic.
            ! Thresholds are set at the top of this module.

            do iCell = 1, nCellsSolve

               maxtemp = maxval(temperature(:,iCell))
               mintemp = minval(temperature(:,iCell))
          
               if (maxtemp > maxtempThreshold) then
                  write(stderrUnit,*) 'maxtemp > maxtempThreshold: iCell, maxtemp =', iCell, maxtemp
                  write(stderrUnit,*) 'thickness =', thickness(iCell)
                  write(stderrUnit,*) 'temperature:'
                  do k = 1, nVertLevels
                     write(stderrUnit,*) k, temperature(k,iCell)
                  enddo
                  call mpas_dmpar_global_abort("MPAS-landice: An error has occurred in li_thermal. Aborting...")
               endif

               if (mintemp < mintempThreshold) then
                  write(stderrUnit,*) 'mintemp < mintempThreshold: iCell, mintemp =', iCell, mintemp
                  write(stderrUnit,*) 'thickness =', thickness(iCell)
                  write(stderrUnit,*) 'temperature:'
                  do k = 1, nVertLevels
                     write(stderrUnit,*) k, temperature(k,iCell)
                  enddo
                  call mpas_dmpar_global_abort("MPAS-landice: An error has occurred in li_thermal. Aborting...")
               endif
          
            enddo   ! iCell

         end select   ! config_thermal_solver

         ! Compute basal melting for floating ice.
         ! Note: This subroutine is called whether or not the temperature/enthalpy is prognostic.
         !       Some tests (e.g., MISMIP+) have basal melting of floating ice without prognostic temperature.

         call basal_melt_floating_ice(&
              config_basal_mass_bal_float,  &
              nCellsSolve,                  &
              xCell,                        &
              thermalCellMask,              &
              li_mask_is_floating_ice_int(cellMask), &
              lowerSurface,                 &
              bedTopography,                &
              config_sea_level,             &
              config_bmlt_float_flux,       &
              config_bmlt_float_xlimit,     &
              floatingBasalMassBal)

         do iCell = 1, nCellsSolve

            ! It is possible that internal melting was computed above for floating ice and assigned
            !  to the groundedBasalMassBal array.  If so, then transfer it to floatingBasalMassBal.

            if (thermalCellMask(iCell) == 1 .and. li_mask_is_floating_ice(cellMask(iCell)) .and. &
                 groundedBasalMassBal(iCell) /= 0.0_RKIND) then
               floatingBasalMassBal(iCell) = floatingBasalMassBal(iCell) + groundedBasalMassBal(iCell)
               groundedBasalMassBal(iCell) = 0.0_RKIND
            endif

            ! Compute total basal mass balance
            basalMassBal(iCell) = groundedBasalMassBal(iCell) + floatingBasalMassBal(iCell)

            if (config_print_thermal_info .and. indexToCellID(iCell) == config_stats_cell_ID) then
               write(stdoutUnit,*) 'iCell, basal mass balance (m/yr) (grounded, floating, total):', &
                    iCell, groundedBasalMassBal(iCell)*scyr/rhoi, floatingBasalMassBal(iCell)*scyr/rhoi, basalMassBal(iCell)*scyr/rhoi
            endif

         enddo

         ! clean up
         call mpas_deallocate_scratch_field(thermalCellMaskField, .true.)
         if (allocated(subdiagonal)) deallocate(subdiagonal)
         if (allocated(diagonal)) deallocate(diagonal)
         if (allocated(superdiagonal)) deallocate(superdiagonal)
         if (allocated(rhs)) deallocate(rhs)
         if (allocated(solution)) deallocate(solution)
         if (allocated(diffusivity)) deallocate(diffusivity)

         block => block % next
      enddo   ! associated(block)

      ! halo updates
      ! Note: This subroutine does not change the ice thickness, so no thickness update is needed.

      call mpas_timer_start("halo updates")

      call mpas_dmpar_field_halo_exch(domain, 'surfaceTemperature')
      call mpas_dmpar_field_halo_exch(domain, 'basalTemperature')
      call mpas_dmpar_field_halo_exch(domain, 'temperature')

      if (trim(config_thermal_solver) == 'enthalpy') then
         ! prognostic variables are temperature and waterfrac, so need to update waterfrac too
         call mpas_dmpar_field_halo_exch(domain, 'waterfrac')
      endif

      call mpas_timer_stop("halo updates")

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
!>  This routine computes heat dissipation associated with strain heating
!>  in the ice interior, based on the shallow-ice approximation.
!-----------------------------------------------------------------------

  subroutine li_heat_dissipation_sia(domain, err)

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
    
     logical, pointer :: &
          config_print_thermal_info   ! if true, print debug info

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
 
     integer, dimension(:), pointer :: &
          indexToCellID            ! list of global cell IDs

     integer, pointer :: &
          config_stats_cell_ID     ! global ID for diagnostic cell

     real (kind=RKIND), pointer :: &
          config_flowLawExponent   ! flow law exponent

     real (kind=RKIND) :: &
          thicknessEdge,         & ! thickness averaged to edge
          weightEdge               ! edge weight for averaging to cell center

     real (kind=RKIND), dimension(:), allocatable ::  &
          flowParamAEdge           ! flow parameter averaged to edge

     integer :: iCell, iCell1, iCell2, iEdge, iEdgeOnCell, iLayer

     ! Here are notes from the Glimmer calculation of heat dissipation:
     !
     !     "Two methods of doing this calculation: 
     !      1. find dissipation at u-pts and then average
     !      2. find dissipation at H-pts by averaging quantities from u-pts
     !     (2) works best for eismint divide (symmetry) but (1) may be better for full expts"
     !
     ! Glimmer uses (2).
     !
     ! For MPAS LI we use the C-grid variant of (1); we find the dissipation at edges and then average to cell centers.
     ! The heating rate phi, defined on an edge, is given by
     !
     !      phi = 2 * A(T) * (sigma * rhoi * g * H * |grad(s)|)^(n+1)
     ! 
     ! where A(T) is the flow factor, sigma is the vertical coordinate of the layer,
     ! H is the ice thickness averaged to the edge, and grad(s) is the surface elevation gradient.
     !
     ! phi has units of W m^{-3}.
     ! The heat dissipation in deg/s is given by phi / (rhoi * cp_ice).

     err = 0

     ! block loop
     block => domain % blocklist
     do while (associated(block))

        ! get pools
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
        call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
        call mpas_pool_get_subpool(block % structs, 'thermal', thermalPool)
        call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
        
        ! get fields from the mesh pool
        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
        call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
        call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

        call mpas_pool_get_array(meshPool, 'layerCenterSigma', layerCenterSigma)
        call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
        call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
        call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
        call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
        call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
        call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
        call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)  ! diagnostic only
                
        ! get fields from the geometry pool
        call mpas_pool_get_array(geometryPool, 'thickness', thickness)
        call mpas_pool_get_array(geometryPool, 'slopeEdge', slopeEdge)
        
        ! get fields from the velocity pool
        call mpas_pool_get_array(velocityPool, 'flowParamA', flowParamA)
        
        ! get fields from the thermal pool
        call mpas_pool_get_array(thermalPool, 'heatDissipation', heatDissipation)
        
        ! get config parameters
        call mpas_pool_get_config(liConfigs, 'config_flowLawExponent', config_flowLawExponent)
        call mpas_pool_get_config(liConfigs, 'config_stats_cell_ID', config_stats_cell_ID)
        call mpas_pool_get_config(liConfigs, 'config_print_thermal_info', config_print_thermal_info)

        ! get scratch fields
        call mpas_pool_get_field(scratchPool, 'workLevelEdge', heatDissipationEdgeField)
        call mpas_allocate_scratch_field(heatDissipationEdgeField, .true.)
        heatDissipationEdge => heatDissipationEdgeField % array
        
        allocate(flowParamAEdge(nVertLevels))
        
        if (config_print_thermal_info) then
           write(stdoutUnit,*) 'Compute SIA heat dissipation'
        endif

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
              !TODO - Verify that this equation gives the right answer
              heatDissipationEdge(:,iEdge) = 2.0_RKIND * flowParamAEdge(:) * &
                   (layerCenterSigma(:) * rhoi * gravity * thicknessEdge * abs(slopeEdge(iEdge))) ** (config_flowLawExponent + 1.0_RKIND)

              if (config_print_thermal_info .and. indexToCellID(iCell1) == config_stats_cell_ID) then
!                 write(stdoutUnit,*) ' '
!                 write(stdoutUnit,*) 'Heat dissipation, iEdge, iCell1 =', iEdge, iCell1
!                 write(stdoutUnit,*) 'thicknessEdge:', thicknessEdge
!                 write(stdoutUnit,*) 'flowParamAEdge:', flowParamAEdge(:)
!                 write(stdoutUnit,*) 'slopeEdge:', slopeEdge(iEdge)
!                 write(stdoutUnit,*) 'heatDissipationEdge:', heatDissipationEdge(:,iEdge)
              endif

           else  ! one neighbor cell does not exist
                 !TODO = Confirm that the dissipation is not needed at such edges

              heatDissipationEdge(:,iEdge) = 0.0_RKIND

           endif
           
         enddo  ! iEdge
        
        ! average the heat dissipation to cell centers
        
        do iCell = 1, nCells

           heatDissipation(:,iCell) = 0.0_RKIND

           if (config_print_thermal_info .and. indexToCellID(iCell) == config_stats_cell_ID) then
              write(stdoutUnit,*) ' '
              write(stdoutUnit,*) 'Heat dissipation, iCell (local, global ID) =', iCell, indexToCellID(iCell)
              write(stdoutUnit,*) 'iEdgeOnCell, iEdge, weightEdge, heatDissipationEdge (layer 1)'
           endif

           do iEdgeOnCell = 1, nEdgesOnCell(iCell)
           
              iEdge = edgesOnCell(iEdgeOnCell,iCell)
              weightEdge = 0.25_RKIND*dcEdge(iEdge)*dvEdge(iEdge) / areaCell(iCell)

              heatDissipation(:,iCell) = heatDissipation(:,iCell) + weightEdge * heatDissipationEdge(:,iEdge)

              if (config_print_thermal_info .and. indexToCellID(iCell) == config_stats_cell_ID) then
                 write(stdoutUnit,*) iEdgeOnCell, iEdge, weightEdge, heatDissipationEdge(1,iEdge), heatDissipation(1,iCell)
              endif

           enddo  ! iEdgeOnCell

           if (config_print_thermal_info .and. indexToCellID(iCell) == config_stats_cell_ID) then
              write(stdoutUnit,*) 'heatDissipation (deg/s) for cell', iCell
              do iLayer = 1, nVertLevels
                 write(stdoutUnit,*) iLayer, heatDissipation(iLayer,iCell) / (rhoi*cp_ice)
              enddo
           endif

        enddo  ! iCell

        ! convert units from W/m^3 to deg/s

        heatDissipation(:,:) = heatDissipation(:,:) / (rhoi * cp_ice)

        ! clean up
        deallocate(flowParamAEdge)
        call mpas_deallocate_scratch_field(heatDissipationEdgeField, .true.)

        block => block % next
     enddo  ! associated(block) 

   end subroutine li_heat_dissipation_sia

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  !  routine li_basal_friction
!
!> \brief MPAS heat flux due to basal friction for SIA dynamics
!> \author William Lipscomb
!> \date   October 2015
!> \details
!>  This routine computes the heat flux due to basal friction at the
!>  base of the ice, based on the shallow-ice approximation.
!-----------------------------------------------------------------------

   subroutine li_basal_friction(domain, err)

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
     type (mpas_pool_type), pointer :: velocityPool
     type (mpas_pool_type), pointer :: thermalPool 
    
     logical, pointer :: &
          config_print_thermal_info  ! if true, print debug info

     integer, pointer :: &
          config_stats_cell_ID       ! global ID for diagnostic cell

     integer, dimension(:), pointer :: &
          indexToCellID              ! list of global cell IDs

     integer, pointer :: &
          nCellsSolve                ! number of cells

     real(kind=RKIND), dimension(:), pointer :: &
          basalSpeed,              & ! basal ice speed, reconstructed at cell centers (m s^{-1})
          betaSolve,               & ! basal traction parameter (Pa m^{-1} s); use betaSolve to treat floating ice correctly
          basalFrictionFlux          ! heat flux due to basal friction (W^{m-2}), computed in this subroutine

     integer :: iCell

     !----------------------------------------------------------------
     ! Compute the heat flux due to basal friction, given the basal speed
     ! and basal friction parameter (beta) fields.
     !
     ! Assume a sliding law of the form 
     !
     !     tau_b = beta * basalSpeed
     ! 
     ! The frictional heat flux (W/m^2) is given by (e.g., Cuffey & Paterson, p. 418)
     !
     !     basalFrictionFlux = tau_b * basalSpeed = beta * basalSpeed**2
     !
     ! Note: Currently (Nov. 2015), beta is used by the HO solver only,
     !       and the SIA solver assumes no slip (basalSpeed = 0).
     ! TODO: Modify the SIA solver to allow sliding.
     !----------------------------------------------------------------

     err = 0

     ! block loop
     block => domain % blocklist
     do while (associated(block))

        ! get pools
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
        call mpas_pool_get_subpool(block % structs, 'thermal', thermalPool)
        
        ! get dimensions
        call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

        ! get fields from the mesh pool
        call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)  ! diagnostic only
                
        ! get fields from the velocity pool
        call mpas_pool_get_array(velocityPool, 'basalSpeed', basalSpeed)
        call mpas_pool_get_array(velocityPool, 'betaSolve', betaSolve)
        
        ! get fields from the thermal pool
        call mpas_pool_get_array(thermalPool, 'basalFrictionFlux', basalFrictionFlux)
        
        ! get config parameters
        call mpas_pool_get_config(liConfigs, 'config_stats_cell_ID', config_stats_cell_ID)
        call mpas_pool_get_config(liConfigs, 'config_print_thermal_info', config_print_thermal_info)

        if (config_print_thermal_info) write(stdoutUnit,*) 'Compute basal friction flux'

        ! Compute basal frictional heating for each cell
        basalFrictionFlux(:) = betaSolve(:) * basalSpeed(:)**2

        if (config_print_thermal_info) then
           do iCell = 1, nCellsSolve
              if (indexToCellID(iCell) == config_stats_cell_ID) then
                 write(stdoutUnit,*) 'iCell, betaSolve, basalSpeed, basalFrictionFlux:', &
                   iCell, betaSolve(iCell), basalSpeed(iCell), basalFrictionFlux(iCell)
              endif
           enddo
        endif

        block => block % next
     enddo  ! associated(block) 

  end subroutine li_basal_friction

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  !  routine li_temperature_to_enthalpy
!
!> \brief MPAS convert temperature to enthalpy
!> \author William Lipscomb
!> \date   October 2015
!> \details
!>  This routine computes the enthalpy in each layer of an ice column,
!>  given the temperature and water fraction.
!-----------------------------------------------------------------------
   
    subroutine li_temperature_to_enthalpy(&
         layerCenterSigma,  &
         thickness,         &
         temperature,       &
         waterfrac,         &
         enthalpy)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(in) :: &
           layerCenterSigma      !< Input: sigma coordinate at midpoint of each layer

      real (kind=RKIND), intent(in) ::  &
           thickness             !< Input: ice thickness
 
      real (kind=RKIND), dimension(:), intent(in) :: &
           temperature,        & !< Input: interior ice temperature
           waterfrac             !< Input: interior water fraction

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(out) :: &
           enthalpy              !< Output:  interior ice enthalpy

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(size(layerCenterSigma)) :: &
           pmpTemperature         ! pressure melting point temperature

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
                      + waterfrac(k) * rho_water * (cp_ice * pmpTemperature(k) + latent_heat_ice)
      end do

    end subroutine li_temperature_to_enthalpy

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  !  routine li_enthalpy_to_temperature
!
!> \brief MPAS convert enthalpy to temperature
!> \author William Lipscomb
!> \date   October 2015
!> \details
!>  This routine computes the temperature and water fraction in each layer
!>  of an ice column, given the enthalpy.
!-----------------------------------------------------------------------

    subroutine li_enthalpy_to_temperature(&
         layerCenterSigma,   &
         thickness,          &
         enthalpy,           &
         temperature,        &
         waterfrac,          &
         surfaceEnthalpy,    &
         surfaceTemperature, &
         basalEnthalpy,      &
         basalTemperature)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(in) :: &
           layerCenterSigma      !< Input: sigma coordinate at midpoint of each layer

      real (kind=RKIND), intent(in) ::  &
           thickness             !< Input: ice thickness
 
      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(inout) :: &
           enthalpy              !< Input/output: interior ice enthalpy

      real (kind=RKIND), intent(inout) , optional :: &
           surfaceEnthalpy,    & !< Input/output: surface ice enthalpy
           basalEnthalpy         !< Input/output: basal ice enthalpy

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(out) :: &
           temperature           !< Output: interior ice temperature

      real (kind=RKIND), dimension(:), intent(out) :: &
           waterfrac             !< Output: interior water fraction

      real (kind=RKIND), intent(out), optional :: &
           surfaceTemperature, & !< Output: surface ice temperature
           basalTemperature      !< Output: basal ice temperature

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(size(layerCenterSigma)) :: pmpTemperature
      real (kind=RKIND), dimension(size(layerCenterSigma)) :: pmpEnthalpy

      real (kind=RKIND) :: pmpTemperatureBed
      real (kind=RKIND) :: pmpEnthalpyBed

      integer :: k, nVertLevels

      nVertLevels = size(layerCenterSigma)

      ! Find pressure melting point temperature in ice interior
      call pressure_melting_point_column(&
           layerCenterSigma,  &
           thickness,         &
           pmpTemperature)

      ! find pressure melting point temperature at bed
      call pressure_melting_point(&
           thickness, &
           pmpTemperatureBed)

      ! Find pressure melting point enthalpy
      pmpEnthalpy(:) = pmpTemperature(:) * rhoi*cp_ice
      pmpEnthalpyBed = pmpTemperatureBed * rhoi*cp_ice

      ! Solve for temperature and waterfrac in ice interior

      ! ice interior

      do k = 1, nVertLevels
         if (enthalpy(k) >= pmpEnthalpy(k)) then   ! temperate ice
            temperature(k) = pmpTemperature(k)
            waterfrac(k) = (enthalpy(k) - pmpenthalpy(k)) /        &
                          ((rho_water-rhoi)*cp_ice*pmpTemperature(k) + rho_water*latent_heat_ice)
         else   ! cold ice
            temperature(k) = enthalpy(k) / (rhoi*cp_ice)
            waterfrac(k) = 0.0_RKIND
         endif
      enddo   ! k

      ! surface temperature (optional)

      if (present(surfaceEnthalpy) .and. present(surfaceTemperature)) then

         if (surfaceEnthalpy >= 0.0_RKIND) then   ! temperate ice
            surfaceTemperature = 0.0_RKIND
            ! Reset surfaceEnthalpy to agree with the surface temperature.
            ! This is consistent with energy conservation because the top surface
            !  is infinitesimally thin.
            surfaceEnthalpy = 0.0_RKIND
         else   ! cold ice
            surfaceTemperature = surfaceEnthalpy / (rhoi*cp_ice)
         endif

      endif

      ! basal temperature (optional)
      
      if (present(basalEnthalpy) .and. present(basalTemperature)) then

         k = nVertLevels + 1
         if (basalEnthalpy >= pmpEnthalpyBed) then   ! temperate ice
            basalTemperature = pmpTemperatureBed
            ! Reset basalEnthalpy to agree with the surface temperature.
            ! This is consistent with energy conservation because the lower surface
            !  is infinitesimally thin.
            basalEnthalpy = pmpEnthalpyBed
         else   ! cold ice
            basalTemperature = basalEnthalpy / (rhoi*cp_ice)
         endif
         
      endif

    end subroutine li_enthalpy_to_temperature

!***********************************************************************
!***********************************************************************
! Private subroutines:
!***********************************************************************
!***********************************************************************

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  !  routine temperature_matrix_elements
!
!> \brief MPAS matrix elements for vertical temperature solver
!> \author William Lipscomb
!> \date   October 2015
!> \details
!>  This routine computes matrix elements for the tridiagonal solve
!>  of the vertical temperature profile.
!-----------------------------------------------------------------------

    subroutine temperature_matrix_elements(&
         deltat,                &
         nVertLevels,           &
         layerCenterSigma,      &
         dsigmaTerm,            &
         floatingMask,          &
         thickness,             &
         temperature,           &
         surfaceTemperature,    &
         basalTemperature,      &
         heatDissipation,       &
         basalHeatFlux,         &
         basalFrictionFlux,     &
         subd,                  &
         diag,                  &
         supd,                  &
         rhs)

      ! Note: Matrix elements (subd, supd, diag, rhs) are indexed from 1 to nVertLevels+2,
      !       whereas temperature is indexed from 1 to nVertLevels.
      !       The first row of the matrix is the equation for surfaceTemperature, and
      !       the last row is the equation for basalTemperature.
      
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real(kind=RKIND), intent(in) :: &
           deltat                   !< Input: time step (s)

      integer, intent(in) :: &
           nVertLevels              !< Input: number of vertical layers

      real(kind=RKIND), dimension(nVertLevels), intent(in) :: &
           layerCentersigma         !< Input: sigma coordinate at layer midpoints

      real(kind=RKIND), dimension(:,:), intent(in) :: &
           dsigmaTerm               !< Input: vertical grid quantities

      integer, intent(in) :: &
           floatingMask             !< Input: = 1 where ice is floating, else = 0

      real(kind=RKIND), intent(in) ::  &
           thickness                !< Input: ice thickness (m)

      real(kind=RKIND), dimension(nVertLevels), intent(in) :: &
           temperature              !< Input: ice temperature (deg C)

      real(kind=RKIND), intent(in) :: &
           surfaceTemperature,    & !< Input: surface ice temperature (deg C)
           basalTemperature         !< Input: basal ice temperature (deg C)

      real(kind=RKIND), dimension(nVertLevels), intent(in) :: &
           heatDissipation          !< Input: interior heat dissipation (deg/s)

      real(kind=RKIND), intent(in) :: &
           basalHeatFlux            !< Input: geothermal heat flux (W m-2), positive up

      real(kind=RKIND), intent(in) :: &
           basalFrictionFlux        !< Input: basal friction heat flux (W m-2), >= 0

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real(kind=RKIND), dimension(:), intent(out) :: &
           subd, diag, supd, rhs    !< Output: tridiagonal matrix elements

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      real(kind=RKIND) :: pmpTemperatureBed  ! pressure melting temperature at bed
      real(kind=RKIND) :: factor             ! term in tridiagonal matrix
      real(kind=RKIND) :: dsigmaBot          ! bottom layer thickness in sigma coords

      ! Compute subdiagonal, diagonal, and superdiagonal matrix elements,
      ! along with the right-hand-side vector

      ! upper boundary: set to surfaceTemperature

      subd(1) = 0.0_RKIND
      diag(1) = 1.0_RKIND
      supd(1) = 0.0_RKIND
      rhs(1)  = surfaceTemperature
      
      ! ice interior, layers 1:nVertLevels  (matrix elements 2:nVertLevels+1)

      factor = deltat * iceConductivity / (rhoi*cp_ice) / thickness**2

      subd(2:nVertLevels+1) = -factor * dsigmaTerm(1:nVertLevels,1)
      supd(2:nVertLevels+1) = -factor * dsigmaTerm(1:nVertLevels,2)
      diag(2:nVertLevels+1) = 1.0_RKIND - subd(2:nVertLevels+1) - supd(2:nVertLevels+1)
      rhs(2:nVertLevels+1)  = temperature(1:nVertLevels) + heatDissipation(1:nVertLevels)*deltat

      ! basal boundary:
      ! for grounded ice, a heat flux is applied
      ! for floating ice, the basal temperature is held constant

      !Note: If T(upn) < T_pmp, then require dT/dsigma = H/k * (G + taub*ubas)
      !       That is, net heat flux at lower boundary must equal zero.
      !      If T(upn) >= Tpmp, then set T(upn) = Tpmp
      
      if (floatingMask == 1) then

         subd(nVertLevels+2) = 0.0_RKIND
         diag(nVertLevels+2) = 1.0_RKIND
         supd(nVertLevels+2) = 0.0_RKIND
         rhs(nVertLevels+2)  = basalTemperature
         
      else    ! grounded ice

         call pressure_melting_point(thickness, pmpTemperatureBed)

         if (abs(basalTemperature - pmpTemperatureBed) < 0.001_RKIND) then

            ! hold basal temperature at pressure melting point

            subd(nVertLevels+2) = 0.0_RKIND
            diag(nVertLevels+2) = 1.0_RKIND
            supd(nVertLevels+2) = 0.0_RKIND
            rhs(nVertLevels+2)  = pmpTemperatureBed

         else   ! frozen at bed
                ! maintain balance of heat sources and sinks
                ! (conductive flux, geothermal flux, and basal friction)

            ! Note: basalHeatFlux is generally >= 0, since defined as positive up

            ! calculate dsigma for the bottom layer between the basal boundary and the temperature point above
            dsigmaBot = 1.0_RKIND - layerCentersigma(nVertLevels)

            ! backward Euler flux basal boundary condition
            subd(nVertLevels+2) = -1.0_RKIND
            diag(nVertLevels+2) =  1.0_RKIND
            supd(nVertLevels+2) =  0.0_RKIND
            rhs(nVertLevels+2)  = (basalFrictionFlux + basalHeatFlux) * dsigmaBot*thickness / iceConductivity
            
         endif   ! melting or frozen
         
      end if     ! floating or grounded

    end subroutine temperature_matrix_elements

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  !  routine enthalpy_matrix_elements
!
!> \brief MPAS matrix elements for vertical enthalpy solver
!> \author William Lipscomb
!> \date   October 2015
!> \details
!>  This routine computes matrix elements for the tridiagonal solve
!>  of the vertical enthalpy profile.
!-----------------------------------------------------------------------

    subroutine enthalpy_matrix_elements(&
         deltat,                &
         nVertLevels,           &
         layerCentersigma,      &
         dsigmaTerm,            &
         floatingMask,          &
         thickness,             &
         temperature,           &
         surfaceTemperature,    &
         basalTemperature,      &
         waterfrac,             &
         enthalpy,              &
         heatDissipation,       &
         basalHeatFlux,         &
         basalFrictionFlux,     &
         diffusivity,           &
         subd,                  &
         diag,                  &
         supd,                  &
         rhs)

      ! Note: Matrix elements (subd, supd, diag, rhs) are indexed from 1 to nVertLevels+2,
      ! whereas temperature and enthalpy are indexed from 1 to nVertLevels.
      ! The first row of the matrix is the equation for surface enthalpy, and
      ! the last row is the equation for basal enthalpy.

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real(kind=RKIND), intent(in) :: &
           deltat                   !< Input: time step (s)

      integer, intent(in) :: &
           nVertLevels              !< Input: number of vertical layers

      real(kind=RKIND), dimension(nVertLevels), intent(in) :: &
           layerCentersigma         !< Input: sigma coordinate at layer midpoints

      real(kind=RKIND), dimension(:,:), intent(in) :: &
           dsigmaTerm               !< Input: vertical grid quantities

      integer, intent(in) :: &
           floatingMask             !< Input: = 1 where ice is floating, else = 0

      real(kind=RKIND), intent(in) :: &
           thickness                !< Input: ice thickness (m)

      real(kind=RKIND), dimension(nVertLevels), intent(in) :: &
           temperature,           & !< Input: ice temperature (deg C)
           waterfrac,             & !< Input: water fraction (unitless)
           enthalpy                 !< Input: specific enthalpy (J/m^3)

      real(kind=RKIND), intent(in) :: &
           surfaceTemperature,    & !< Input: surface ice temperature (deg C)
           basalTemperature         !< Input: basal ice temperature (deg C)

      real(kind=RKIND), dimension(nVertLevels), intent(in) :: &
           heatDissipation          !< Input: interior heat dissipation (deg/s)

      real(kind=RKIND), intent(in) :: &
           basalHeatFlux            !< Input: geothermal heat flux (W m-2), positive up

      real(kind=RKIND), intent(in) :: &
           basalFrictionFlux        !< Input: basal friction heat flux (W m-2), >= 0

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real(kind=RKIND), dimension(:), intent(out) :: &
           diffusivity              !< Output: half-node diffusivity (m^2/s) for enthalpy
                                    !          located at layer interfaces
      
      real(kind=RKIND), dimension(:), intent(out) :: &
           subd, diag, supd, rhs    !< Output: matrix elements

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      real(kind=RKIND), dimension(nVertLevels) :: &
           pmpTemperature           ! pressure melting point temperature in interior (deg C)

      real(kind=RKIND) :: &
           pmpTemperatureBed        ! pressure melting point temperature at bed (deg C)

      real(kind=RKIND), dimension(0:nVertLevels+1) :: &
           enthalpyTemp             ! temperature part of specific enthalpy (J/m^3)

      real(kind=RKIND) :: surfaceEnthalpy       ! surface enthalpy = surfaceTemperature*rhoi*cp_ice
      real(kind=RKIND) :: basalEnthalpy         ! basal enthalpy = basalTemperature*rhoi*cp_ice
      real(kind=RKIND) :: dsigmaBot             ! bottom layer thickness in sigma coords
      real(kind=RKIND) :: diffusivityCold       ! cold ice diffusivity
      real(kind=RKIND) :: diffusivityTemperate  ! temperate ice diffusivity
      real(kind=RKIND) :: factor                ! term in tridiagonal matrix
      real(kind=RKIND) :: dEnthalpy             ! enthalpy difference between adjacent layers
      real(kind=RKIND) :: dEnthalpyTemp         ! difference in temperature component of enthalpy between adjacent layers
      real(kind=RKIND) :: avgFactor             ! factor for averaging diffusivity, 0 <= avgFactor <= 1

      integer  :: k
      
      logical, parameter :: &
           harmonic_avg  = .false.  ! if true, take the harmonic average of diffusivity in adjacent layers
                                    ! if false, take the arithmetic average

      diffusivityCold = iceConductivity / (rhoi*cp_ice)
      diffusivityTemperate = diffusivityCold / 100.0_RKIND
        
      ! set enthalpy at surface and bed
      surfaceEnthalpy = surfaceTemperature * rhoi*cp_ice
      basalEnthalpy = basalTemperature * rhoi*cp_ice

      ! find pmpTemperature for this column (interior nodes and boundary)

      call pressure_melting_point_column(&
           layerCenterSigma, &
           thickness,        &
           pmpTemperature)

      call pressure_melting_point(&
           thickness,        &
           pmpTemperatureBed)

      !--------------------------------------------------------------------
      ! Compute the enthalpy diffusivity at layer interfaces.
      ! Let d(enth)/dz = the gradient of enthalpy
      ! Can write
      !    d(enth)/dz = d(enth_T)/dz + d(enth_w)/dz,
      ! where
      !    enth_T = (1-phi_w) * rhoi*ci*T
      !    enth_w =    phi_w  * rhow*(L + ci*Tpmp)
      !    phi_w = water fraction
      ! 
      ! Now let f = d(enth_T)/dz / d(enth)/dz
      !   (f -> 0 if f is computed to be negative)
      ! For cold ice, f = 1 and diffusivity = diffusivityCold
      ! For temperate ice, f ~ 0 and diffusivity = diffusivityTemperate
      ! At the interface between cold and temperate ice,
      !  f ~ 0 if the temperate ice has large phi_w, but
      !  f ~ 1 if the temperate ice has close to zero phi_w.
      ! Two ways to average:
      ! (1) arithmetic average:  diffusivity = f*diffusivityCold + (1-f)*diffusivityTemperate
      ! (2) harmonic average:    diffusivity = 1 / (f/diffusivityCold + (1-f)/diffusivityTemperate).
      ! Both methods have the same asymptotic values at f = 0 or 1,
      !  but the arithmetic average gives greater diffusivity for
      !  intermediate values.
      !
      ! Still to be determined which is more accurate.
      ! The harmonic average allows large temperature gradients between the 
      !  bottom layer and the next layer up; the arithmetic average gives
      !  smoother gradients.
      ! Currently (as of Oct. 2015), the arithmetic average is the default.
      !--------------------------------------------------------------------
      !
      ! At each temperature point, compute the temperature part of the enthalpy.
      ! Note: enthalpyTemp = enthalpy for cold ice, enthalpyTemp < enthalpy for temperate ice
      
      enthalpyTemp(0) = surfaceEnthalpy
      do k = 1, nVertLevels
         enthalpyTemp(k) = (1.0_RKIND - waterfrac(k)) * rhoi*cp_ice*temperature(k)
      enddo
      enthalpyTemp(nVertLevels+1) = basalEnthalpy

      ! Compute factors relating the temperature gradient to the total enthalpy gradient.
      ! Use these factors to average the diffusivity between adjacent temperature points.

      do k = 1, nVertLevels+1

         if (k == 1) then
            dEnthalpy = enthalpy(k) - surfaceEnthalpy
         elseif (k == nVertLevels+1) then
            dEnthalpy = basalEnthalpy - enthalpy(k-1)
         else
            dEnthalpy = enthalpy(k) - enthalpy(k-1)
         endif

         dEnthalpyTemp = enthalpyTemp(k) - enthalpyTemp(k-1)   ! = dEnthalpy in cold ice, < dEnthalpy in temperate ice

         if (abs(dEnthalpy) > 1.e-20_RKIND * rho_water * latent_heat_ice) then
            avgFactor = max(0.0_RKIND, dEnthalpyTemp/dEnthalpy)
            avgFactor = min(1.0_RKIND, avgFactor)
         else
            avgFactor = 0.0_RKIND
         endif
         
         if (harmonic_avg) then  ! take a harmonic average
                                 ! This gives slower cooling of temperate layers and allows
                                 !  large temperature gradients between cold and temperate layers
            diffusivity(k) = 1.0_RKIND / ((avgFactor/diffusivityCold) + (1.0_RKIND - avgFactor)/diffusivityTemperate)
         else   ! take an arithmetic average
                ! This gives faster cooling of temperate layers and smaller gradients 
            diffusivity(k) = avgFactor*diffusivityCold + (1.0_RKIND - avgFactor)*diffusivityTemperate
         endif

      end do

      ! Compute subdiagonal, diagonal, and superdiagonal matrix elements
      ! Assume backward Euler time stepping
    
      ! upper boundary: set to surfaceEnthalpy = surface temperature*rhoi*cp_ice
      subd(1) = 0.0_RKIND
      diag(1) = 1.0_RKIND
      supd(1) = 0.0_RKIND
      rhs(1)  = surfaceEnthalpy
      
      ! ice interior, layers 1:nVertLevels  (matrix elements 2:nVertLevels+1)

      factor = deltat / thickness**2

      subd(2:nVertLevels+1) = -factor * diffusivity(1:nVertLevels)   * dsigmaTerm(1:nVertLevels,1)                                
      supd(2:nVertLevels+1) = -factor * diffusivity(2:nVertLevels+1) * dsigmaTerm(1:nVertLevels,2)                                
      diag(2:nVertLevels+1) = 1.0_RKIND - subd(2:nVertLevels+1) - supd(2:nVertLevels+1)                              
      rhs(2:nVertLevels+1)  = enthalpy(1:nVertLevels) + heatDissipation(1:nVertLevels) * deltat * rhoi * cp_ice
      
      ! Note: heatDissipation has units of phi/rhoi/cp_ice, where phi has units of W m^-3..
      ! For an enthalpy calculation, we want just phi, hence heatDissipation * rhoi * cp_ice
        
      ! Basal boundary:
      ! For grounded ice, a heat flux is applied.
      ! For floating ice, the basal temperature is held constant.
      !      If basalTemperature < T_pmp, then require dT/dsigma = H/k * (G + taub*ubas)
      !       That is, net heat flux at lower boundary must equal zero.
      !      If basalTemperature >= Tpmp, then set basalTemperature = Tpmp
      
      if (floatingMask == 1) then

         subd(nVertLevels+2) = 0.0_RKIND
         diag(nVertLevels+2) = 1.0_RKIND
         supd(nVertLevels+2) = 0.0_RKIND 
         rhs(nVertLevels+2)  = basalEnthalpy
         
      else    ! grounded ice
         
         if (abs(temperature(nVertLevels) - pmpTemperature(nVertLevels)) < 0.001_RKIND) then

            ! Positive-thickness basal temperate boundary Layer
            !WHL - Not sure whether this condition is ideal. It implies that basalEnthalpy = enthalpy(nVertLevels).

            subd(nVertLevels+2) = -1.0_RKIND
            diag(nVertLevels+2) =  1.0_RKIND
            supd(nVertLevels+2) =  0.0_RKIND 
            rhs(nVertLevels+2)  =  0.0_RKIND
            
         elseif (abs(basalTemperature -  pmpTemperatureBed) < 0.001_RKIND) then  ! melting
          
            ! Zero-thickness basal temperate boundary layer
            ! Hold basal temperature at pressure melting point
            subd(nVertLevels+2) = 0.0_RKIND
            diag(nVertLevels+2) = 1.0_RKIND
            supd(nVertLevels+2) = 0.0_RKIND
            rhs(nVertLevels+2)  = pmpTemperatureBed * rhoi*cp_ice
            
         else  ! cold ice at bed
          
            ! Frozen at bed: Maintain balance of heat sources and sinks
            ! (conductive flux, geothermal flux, and basal friction)
            
            ! Note: basalHeatFlux is generally >= 0, since defined as positive up
            
            ! calculate dsigma for the bottom layer between the basal boundary and the temperature point above
            dsigmaBot = (1.0_RKIND - layerCentersigma(nVertLevels))
            
            subd(nVertLevels+2) = -1.0_RKIND
            diag(nVertLevels+2) =  1.0_RKIND
            supd(nVertLevels+2) =  0.0_RKIND 
            rhs(nVertLevels+2)  = (basalFrictionFlux + basalHeatFlux) * dsigmaBot * thickness *rhoi * cp_ice / iceConductivity

         endif   ! melting or frozen
         
      end if     ! floating or grounded
      
    end subroutine enthalpy_matrix_elements
    
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  !  routine basal_melt_grounded_ice
!
!> \brief MPAS melt rate at base of grounded ice
!> \author William Lipscomb
!> \date   October 2015
!> \details
!>  This routine computes the melt rate at the base of grounded ice.
!>  Internal melting is included in this melt rate.
!-----------------------------------------------------------------------

    subroutine basal_melt_grounded_ice(&
         config_thermal_solver,       &
         deltat,                      &
         nCellsSolve,                 &
         nVertLevels,                 &
         layerInterfaceSigma,         &
         layerCenterSigma,            &
         iceMask,                     &
         floatingMask,                &
         thickness,                   &
         temperature,                 &
         basalTemperature,            &
         waterfrac,                   &
         enthalpy,                    &
         basalFrictionFlux,           &
         basalHeatFlux,               &
         basalConductiveFlux,         &
         basalWaterThickness,         &
         groundedBasalMassBal)
      
      !-----------------------------------------------------------------
      ! Note: For the temperature solver, any internal temperatures above
      !  the pressure melting point are reset to Tpmp.  Excess energy 
      !  is applied toward melting with immediate drainage to the bed. 
      ! For the enthalpy solver, any meltwater in excess of the maximum allowed
      !  meltwater fraction (0.01 by default) is drained to the bed.
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      character(len=StrKIND), intent(in) :: &
           config_thermal_solver    !< Input: thermal solver option (temperature or enthalpy)

      real(kind=RKIND), intent(in) :: &
           deltat                   !< Input: time step (s)

      integer, intent(in) :: &
           nCellsSolve              !< Input: number of locally owned cells

      integer, intent(in) :: &
           nVertLevels              !< Input: number of vertical layers

      real(kind=RKIND), dimension(:), intent(in) :: &
           layerInterfaceSigma      !< Input: vertical coordinate at layer interfaces

      real(kind=RKIND), dimension(:), intent(in) :: &
           layerCenterSigma         !< Input: vertical coordinate at layer centers

      real(kind=RKIND), dimension(:,:), intent(in) :: &
           enthalpy                 !< Input: enthalpy

      real(kind=RKIND), dimension(:), intent(in) :: &
           thickness,             & !< Input: ice thickness (m)
           basalFrictionFlux,     & !< Input: basal frictional heating flux (W m-2), >= 0
           basalHeatFlux,         & !< Input: geothermal heating flux (W m-2), positive up
           basalConductiveFlux,   & !< Input: heat conducted from ice interior to bed (W m-2), positive down
           basalWaterThickness      !< Input: thickness of basal water layer (m)
      
      integer, dimension(:), intent(in) :: &
           iceMask,               & !< Input: = 1 where ice exists (thickness > config_thermal_thickness), else = 0
           floatingMask             !< Input: = 1 where ice is floating, else = 0

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      real(kind=RKIND), dimension(:,:), intent(inout) :: &
           temperature              !< Input/output: temperature (deg C)

      real(kind=RKIND), dimension(:), intent(inout) :: &
           basalTemperature         !< Input/output: basal temperature (deg C)

      real(kind=RKIND), dimension(:,:), intent(inout) :: &
           waterfrac                !< Input/output: water fraction
      
      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real(kind=RKIND), dimension(:), intent(out):: &
           groundedBasalMassBal     !< Output: basal mass balance for grounded ice (kg/m^2/s): < 0 for melting, > 0 for freeze-on
                                   
      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: k, iCell

      real(kind=RKIND), dimension(nVertLevels) :: &
           pmpTemperature           ! pressure melting point temperature in ice interior

      real(kind=RKIND) :: pmpTemperatureBed     ! pressure melting point temperature at bed
      real(kind=RKIND) :: netBasalFlux          ! heat flux available for basal melting (W/m^2)
      real(kind=RKIND) :: layerThickness        ! layer thickness (m)
      real(kind=RKIND) :: meltEnergy            ! energy available for internal melting (J/m^2)
      real(kind=RKIND) :: internalMeltRate      ! internal melt rate, transferred to bed (m/s)
      real(kind=RKIND) :: excessWater           ! thickness of excess meltwater (m)

      real(kind=RKIND), parameter :: &
           maxWaterfrac = 0.01_RKIND       ! maximum allowed water fraction; excess drains to bed
      
      real(kind=RKIND), parameter :: &
           eps11 = 1.0e-11_RKIND           ! small number
      

      ! Compute melt rate for grounded ice

      groundedBasalMassBal(:) = 0.0_RKIND

      do iCell = 1, nCellsSolve

         if (iceMask(iCell) == 1 .and. floatingMask(iCell) == 0) then  ! ice is present and grounded

            ! Compute basal mass balance
            ! Note: basalMassBal < 0 for melting, > 0 for freeze-on   
            !       basalFrictionFlux >= 0 by definition
            !       basalHeatFlux is positive up, so generally basalHeatFlux >= 0
            !       basalConductiveFlux is positive down, so basalConductiveFlux < 0 for heat flowing from the bed toward the surface
            !
            !       This equation allows for freeze-on (basalMassBal > 0) if the conductive term 
            !        (basalConductiveFlux, positive down) is carrying enough heat away from the boundary.  
            !       But freeze-on requires a local water supply, basalWaterThickness > 0.
            !       When basalWaterThickness = 0, we reset the bed temperature to a value slightly below the melting point.
            !
            !TODO - For the enthalpy scheme, deal with the rare case that the bottom layer melts completely,
            !       and overlying layers with a different enthalpy also melt.
            
            netBasalFlux = basalFrictionFlux(iCell) + basalHeatFlux(iCell) + basalConductiveFlux(iCell)  ! W/m^2
            
            if (abs(netBasalFlux) < eps11) then  ! netBasalFlux might be slightly different from zero because of rounding errors; if so, then zero out
               netBasalFlux = 0.0_RKIND
            endif

            if (trim(config_thermal_solver) == 'enthalpy') then
               groundedBasalMassBal(iCell) = -netBasalFlux / (latent_heat_ice * rhoi - enthalpy(nVertLevels,iCell))
            else   ! temperature solver
               groundedBasalMassBal(iCell) = -netBasalFlux / (latent_heat_ice * rhoi)   ! m/s
            endif

         endif   ! ice is present and grounded

         ! Add internal melting
         ! Note: It is possible to have internal melting for floating ice.
         !       If so, this melting will be switched later from groundedBasalMassBal to floatingBasalMassBal

         if (iceMask(iCell) == 1) then   ! ice is present

            if (trim(config_thermal_solver) == 'enthalpy') then

               ! Add internal melting associated with waterfrac > maxWaterfrac (1%)
               !TODO - Add correction for rhoi/rhow here?  Or melting ice that is already partly melted?

               do k = 1, nVertLevels
                  if (waterfrac(k,iCell) > maxWaterfrac) then

                     ! compute melt rate associated with excess water
                     excessWater = (waterfrac(k,iCell) - maxWaterfrac) * thickness(iCell) * (layerInterfaceSigma(k+1) - layerInterfaceSigma(k))  ! m
                     internalMeltRate = excessWater / deltat

                     ! transfer meltwater to the bed
                     ! Note: It is possible to have internal melting for floating ice.
                     !       If so, then this melting will later be switched from groundedBasalMassBall to floatingBasalMassBal.
                     groundedBasalMassBal(iCell) = groundedBasalMassBal(iCell) - internalMeltRate      ! m/s

                     ! reset waterfrac to max value
                     waterfrac(k,iCell) = maxWaterfrac

                  endif   ! waterfrac > maxWaterfrac
               enddo   ! k
               
            elseif (config_thermal_solver == 'temperature') then
               
               ! Add internal melting associated with T > Tpmp

               call pressure_melting_point_column(&
                    layerCentersigma,     &
                    thickness(iCell),   &
                    pmpTemperature)

               do k = 1, nVertLevels
                  if (temperature(k,iCell) > pmpTemperature(k)) then

                     ! compute excess energy available for melting
                     layerThickness = thickness(iCell) * (layerInterfaceSigma(k+1) - layerInterfaceSigma(k))     ! m
                     meltEnergy = rhoi*cp_ice * (temperature(k,iCell) - pmpTemperature(k)) * layerThickness    ! J/m^2
                     internalMeltRate = meltEnergy / (rhoi * latent_heat_ice * deltat)  ! m/s

                     ! transfer meltwater to the bed
                     ! Note: It is possible to have internal melting for floating ice.
                     !       If so, then this melting will later be switched from groundedBasalMassBall to floatingBasalMassBal.
                     groundedBasalMassBal(iCell) = groundedBasalMassBal(iCell) - internalMeltRate      ! m/s

                     ! reset T to Tpmp
                     temperature(k,iCell) = pmpTemperature(k)

                  endif  ! temperature > pmpTemperature
               enddo   ! k

            endif   ! config_thermal_solver

         endif   ! ice is present

         ! Cap basal temperature at pressure melting point

         if (iceMask(iCell) == 1 .and. floatingMask(iCell) == 0) then  ! ice is present and grounded

            call pressure_melting_point(&
                 thickness(iCell), &
                 pmpTemperatureBed)
            
            basalTemperature(iCell) = min (basalTemperature(iCell), pmpTemperatureBed)
            
            ! If freeze-on was computed above (basalMassBal > 0) and Tbed = Tpmp but no basal water is present, 
            !  then set basalTemperature < Tpmp.
            ! Note: In the matrix element subroutines, we solve for Tbed (instead of holding it at Tpmp) when Tbed < -0.001.
            !       With an offset here of 0.01, we will solve for T_bed at the next timestep.
            ! Note: I don't think energy conservation is violated here, because no energy is associated with
            !       the infinitesimally thin layer at the bed.
            
            if (groundedBasalMassBal(iCell) > 0.0_RKIND .and. basalWaterThickness(iCell) == 0.0_RKIND .and. &
                 basalTemperature(iCell) >= pmpTemperatureBed) then
               basalTemperature(iCell) = pmpTemperatureBed - 0.01_RKIND
            endif
            
         endif   ! ice is present and grounded

         ! change units from m/s to kg/m2/s
         groundedBasalMassBal(iCell) = groundedBasalMassBal(iCell) * rhoi

      enddo   ! iCell

    end subroutine basal_melt_grounded_ice

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
         iceMask,                       &
         floatingMask,                  &
         lowerSurface,                  &
         bedTopography,                 &
         config_sea_level,              &
         config_bmlt_float_flux,        &
         config_bmlt_float_xlimit,      &
         floatingBasalMassBal)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      character(len=strKIND), intent(in) :: &
           config_basal_mass_bal_float  !< Input: option for computing basal mass balance for floating ice

      integer, intent(in) :: &
           nCellsSolve              !< Input: number of locally owned cells

      real(kind=RKIND), dimension(:), intent(in) :: &
           xCell                    !< Input: x coordinate for each cell

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
                                       !         MISMIP+ default value = 975.17 W/m^2 (gives melt rate of 100 m/yr)
           config_bmlt_float_xlimit    !< Input: x value (m) defining region where bmlt_float_flux is applied; melt only where abs(x) > xlimit
                                       !         MISMIP+ default value = 480 km

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real(kind=RKIND), dimension(:), intent(out):: &
           floatingBasalMassBal     !< Output: basal mass balance for floating ice 

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      real(kind=RKIND) :: &
           bmlt_float_rate                          ! constant basal melt rate (m/s)
                                                    ! = config_bmlt_float_flux / (rhoi*latent_heat_ice)
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

      ! Compute melt rate for floating ice

      ! initialize to zero melt
      floatingBasalMassBal(:) = 0.0_RKIND

      if (trim(config_basal_mass_bal_float) == 'none') then

         ! nothing to do; bmlt_float already set to zero

      elseif (trim(config_basal_mass_bal_float) == 'file') then

         ! Do nothing in this routine -- use whatever value was input as data

      elseif (trim(config_basal_mass_bal_float) == 'constant') then

         ! set melt rate to a constant value for floating ice 
         ! allow basal melt in ice-free ocean cells, in case ice is advected to those cells by the transport scheme

         bmlt_float_rate = config_bmlt_float_flux / (rhoi*latent_heat_ice)  ! convert W/m^2 to m/s

         floatingBasalMassBal(:) = 0.0_RKIND

         do iCell = 1, nCellsSolve
            if ( floatingMask(iCell) == 1 .or.  &
                 (bedTopography(iCell) < config_sea_level .and. iceMask(iCell) == 0) ) then   ! ice is present and floating, or ice-free ocean

               ! Provided xCell > bmlt_float_xlimit, prescribe a uniform basal melt rate.
               ! The default is 0.0, but for MISMIP+ the prescribed value of xlimit is 480 km.
               if (abs(xCell(iCell)) >= config_bmlt_float_xlimit) then  ! basal melting is allowed
                  floatingBasalMassBal(iCell) = -bmlt_float_rate
               endif

            endif   ! ice is present and floating
         enddo   ! iCell

         ! change units from m/s to kg/m2/s
         floatingBasalMassBal(:) = floatingBasalMassBal(:) * rhoi

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
                 (bedTopography(iCell) < config_sea_level .and. iceMask(iCell) == 0) ) then   ! ice is present and floating, or ice-free ocean

               hCavity = lowerSurface(iCell) - bedTopography(iCell)
               zDraft = lowerSurface(iCell) - config_sea_level
               floatingBasalMassBal(iCell) = -bmlt_float_omega * tanh(hCavity/bmlt_float_h0) * max(bmlt_float_z0 - zDraft, 0.0_RKIND)

            endif   ! ice is present and floating
         enddo   ! iCell

         ! change units from m/s to kg/m2/s
         floatingBasalMassBal(:) = floatingBasalMassBal(:) * rhoi

      else

         write(stderrUnit,*) 'ERROR: Unknown option selected for config_basal_mass_bal_float:', trim(config_basal_mass_bal_float)
         ! TODO: Actually throw an error and abort the model!

      endif   ! config_basal_mass_bal_float

      
    end subroutine basal_melt_floating_ice

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  !  routine pressure_melting_point_column
!
!> \brief MPAS compute pressure melting point temperature
!> \author William Lipscomb
!> \date   October 2015
!> \details
!>  This routine computes the pressure melting point temperature in each layer
!>  of an ice column, given the thickness.
!-----------------------------------------------------------------------

    subroutine pressure_melting_point_column(&
         layerCenterSigma,  &
         thickness,         &
         pmpTemperature)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(in) :: &
           layerCenterSigma      !< Input: sigma coordinate at midpoint of each layer

      real (kind=RKIND), intent(in) ::  &
           thickness             !< Input: ice thickness

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(out) :: &
           pmpTemperature          !< Output: pressure melting point temperature (deg C)

      ! Compute the pressure melting point temperature in the column

      pmpTemperature(:) = - iceMeltingPointPressureDependence * rhoi * gravity * thickness * layerCenterSigma(:)

    end subroutine pressure_melting_point_column

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  !  routine pressure_melting_point
!
!> \brief MPAS compute pressure melting point temperature
!> \author William Lipscomb
!> \date   October 2015
!> \details
!>  This routine computes the pressure melting point temperature at a
!>  given depth in an ice column.
!-----------------------------------------------------------------------

    subroutine pressure_melting_point(&
         depth,           &
         pmpTemperature)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(in) ::  &
           depth                !< Input: depth in column

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(out) :: &
           pmpTemperature          !< Output: pressure melting point temperature (deg C)

      ! Compute the pressure melting point temperature

      pmpTemperature = - iceMeltingPointPressureDependence * rhoi * gravity * depth

    end subroutine pressure_melting_point


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  !  routine li_compute_pressure_melting_point_fields
!
!> \brief MPAS compute pressure melting point temperature for entire domain
!> \author Matthew Hoffman
!> \date   February 2016
!> \details
!>  This routine computes the pressure melting point temperature for the
!>  entire interior ice domain and the basal boundary.
!>  It makes use of existing routines for PMP calculations.
!>  This is created for diagnostic output purposes, and is not used by the
!>  thermal solver.
!-----------------------------------------------------------------------

    subroutine li_compute_pressure_melting_point_fields(&
         nCells,                 &
         thickness,              &
         layerCenterSigma,       &
         pmpTemperature,         &
         basalPmpTemperature)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      integer, intent(in) :: nCells !< Input: number of cells

      real (kind=RKIND), dimension(:), intent(in) ::  &
           thickness                !< Input: thickness field

      real (kind=RKIND), dimension(:), intent(in) ::  &
           layerCenterSigma         !< Input: layer sigma levels

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      real (kind=RKIND), dimension(:,:), intent(out) :: &
           pmpTemperature          !< Output: pressure melting point temperature (K) for interior ice

      real (kind=RKIND), dimension(:), intent(out) :: &
           basalPmpTemperature     !< Output: pressure melting point temperature (K) at ice-bed interface

      ! Local variables
      integer :: iCell

      do iCell = 1, nCells
         ! Calculate internal ice PMP
         call pressure_melting_point_column(layerCenterSigma, thickness(iCell), pmpTemperature(:,iCell))
         ! Calculate basal PMP
         call pressure_melting_point(thickness(iCell), basalPmpTemperature(iCell))
      enddo

      ! Convert from Celsius to Kelvin
      pmpTemperature = pmpTemperature + kelvin_to_celsius
      basalPmpTemperature = basalPmpTemperature + kelvin_to_celsius

    end subroutine li_compute_pressure_melting_point_fields


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  !  routine tridiag_solver
!
!> \brief MPAS solve tridiagonal matrix
!> \author William Lipscomb
!> \date   October 2015
!> \details
!>  This routine solves a tridiagonal matrix equation, given the matrix 
!>  coefficients and right-hand side.
!-----------------------------------------------------------------------

    !TODO - Move the tridiag solver to a utility module?
    subroutine tridiag_solver(a,b,c,x,y)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real(kind=RKIND), dimension(:), intent(in)  :: a !< Input: Lower diagonal; a(1) is ignored
      real(kind=RKIND), dimension(:), intent(in)  :: b !< Input: Main diagonal
      real(kind=RKIND), dimension(:), intent(in)  :: c !< Input: Upper diagonal; c(n) is ignored
      real(kind=RKIND), dimension(:), intent(in)  :: y !< Input: Right-hand side

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      real(kind=RKIND), dimension(:), intent(out) :: x !< Output: Unknown vector
      
      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      real(kind=RKIND), dimension(size(a)) :: aa
      real(kind=RKIND), dimension(size(a)) :: bb

      integer :: n, i

      n = size(a)
      
      aa(1) = c(1) / b(1)
      bb(1) = y(1) / b(1)
      
      do i = 2, n
         aa(i) = c(i) / (b(i)-a(i)*aa(i-1))
         bb(i) = (y(i)-a(i)*bb(i-1)) / (b(i)-a(i)*aa(i-1))
      end do
      
      x(n) = bb(n)

      do i = n-1, 1, -1
         x(i) = bb(i) - aa(i)*x(i+1)
      end do
      
    end subroutine tridiag_solver

    !***********************************************************************

  end module li_thermal

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||   



