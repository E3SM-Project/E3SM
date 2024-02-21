!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_prescribed
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 7th November 2018
!> \details
!>
!
!-----------------------------------------------------------------------

module seaice_prescribed

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_forcing
  use seaice_forcing, only: seaiceForcingGroups

  implicit none

  private
  save

  public :: &
       seaice_init_prescribed_ice, &
       seaice_run_prescribed_ice

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_init_prescribed_ice
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 7th November 2018
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_init_prescribed_ice(domain)

    type (domain_type), intent(inout) :: &
         domain !< Input/Output:

    logical, pointer :: &
         config_use_prescribed_ice, &
         config_use_prescribed_ice_forcing, &
         config_use_dynamics, &
         config_use_column_itd_thermodynamics, &
         config_do_restart

    character(len=strKIND) :: &
         forcingIntervalMonthly, &
         forcingReferenceTimeMonthly

    call mpas_pool_get_config(domain % configs, "config_use_prescribed_ice", config_use_prescribed_ice)
    if (config_use_prescribed_ice) then

       ! check for compatable options
       call mpas_pool_get_config(domain % configs, "config_use_dynamics", config_use_dynamics)
       if (config_use_dynamics) then
          call mpas_log_write("Prescribed ice mode: Dynamics must be off (config_use_dynamics = false)", MPAS_LOG_CRIT)
       endif

       call mpas_pool_get_config(domain % configs, "config_use_column_itd_thermodynamics", config_use_column_itd_thermodynamics)
       if (config_use_column_itd_thermodynamics) then
          call mpas_log_write("Prescribed ice mode: ITD thermodynamics must be off (config_use_column_itd_thermodynamics = false)", MPAS_LOG_CRIT)
       endif

       call mpas_pool_get_config(domain % configs, "config_use_prescribed_ice_forcing", config_use_prescribed_ice_forcing)
       if (config_use_prescribed_ice_forcing) then

          call MPAS_pool_get_config(domain % configs, "config_do_restart", config_do_restart)

          ! create the monthly forcing group
          call MPAS_forcing_init_group(&
               seaiceForcingGroups, &
               "seaice_prescribed_ice_coverage", &
               domain, &
               '0000-01-01_00:00:00', &
               '0000-01-01_00:00:00', &
               '0001-00-00_00:00:00', &
               config_do_restart)

          forcingIntervalMonthly = "00-01-00_00:00:00"
          forcingReferenceTimeMonthly = "0001-01-15_00:00:00"

          call MPAS_forcing_init_field(&
               domain % streamManager, &
               seaiceForcingGroups, &
               "seaice_prescribed_ice_coverage", &
               "iceCoverage", &
               "prescribedIceForcing", &
               "prescribed_ice", &
               "iceCoverage", &
               "linear", &
               forcingReferenceTimeMonthly, &
               forcingIntervalMonthly)

          call MPAS_forcing_init_field_data(&
               seaiceForcingGroups, &
               "seaice_prescribed_ice_coverage", &
               domain % streamManager, &
               config_do_restart, &
               .false.)

       endif ! config_use_prescribed_ice_forcing

    endif ! config_use_prescribed_ice

  end subroutine seaice_init_prescribed_ice

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_run_prescribed_ice
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 7th November 2018
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_run_prescribed_ice(domain)

    use seaice_constants, only: &
         seaicePuny

    use ice_colpkg, only: &
         colpkg_enthalpy_snow, &
         colpkg_enthalpy_ice, &
         colpkg_salinity_profile

    use seaice_column, only: &     ! colpkg
         seaice_column_reinitialize_fluxes, &
         seaice_column_aggregate

    use icepack_intfc, only: &
         icepack_enthalpy_snow, &
         icepack_salinity_profile

    use seaice_icepack, only: &
         seaice_icepack_enthalpy_ice, &
         seaice_icepack_reinitialize_fluxes, &
         seaice_icepack_aggregate

    type (domain_type), intent(inout) :: &
         domain !< Input/Output:

    logical, pointer :: &
         config_use_prescribed_ice, &
         config_use_prescribed_ice_forcing

    character(len=strKIND), pointer :: &
         config_column_physics_type, &
         config_thermodynamics_type

    type(block_type), pointer :: &
         blockPtr

    type(mpas_pool_type), pointer :: &
         prescribedIcePool, &
         tracersPool, &
         initialPool, &
         oceanCouplingPool, &
         meshPool, &
         velocitySolverPool

    real(kind=RKIND), dimension(:), pointer :: &
         iceCoverage, &
         categoryThicknessLimits, &
         seaFreezingTemperature, &
         latCell, &
         freezingMeltingPotential, &
         uVelocity, &
         vVelocity, &
         oceanStressCellU, &
         oceanStressCellV

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         surfaceTemperature, &
         iceEnthalpy, &
         iceSalinity, &
         snowEnthalpy

    real(kind=RKIND) :: &
         iceThickness, &
         snowThickness, &
         depth, &
         iceTemperature, &
         temperatureGradient

    real(kind=RKIND), pointer :: &
         config_dt

    integer, pointer :: &
         nCells, &
         nCategories, &
         nIceLayers, &
         nSnowLayers

    integer :: &
         iCell, &
         iCategory, &
         iIceLayer, &
         iSnowLayer

    call MPAS_pool_get_config(domain % configs, "config_column_physics_type", config_column_physics_type)            
    call mpas_pool_get_config(domain % configs, "config_use_prescribed_ice", config_use_prescribed_ice)
    call mpas_pool_get_config(domain % configs, "config_thermodynamics_type", config_thermodynamics_type)

    if (config_use_prescribed_ice) then

       ! get ice coverage
       call mpas_pool_get_config(domain % configs, "config_use_prescribed_ice_forcing", config_use_prescribed_ice_forcing)
       if (config_use_prescribed_ice_forcing) then

          call mpas_pool_get_config(domain % configs, 'config_dt', config_dt)

          call MPAS_forcing_get_forcing(&
               seaiceForcingGroups, &
               "seaice_prescribed_ice_coverage", &
               domain % streamManager, &
               config_dt)

       endif ! config_use_prescribed_ice_forcing

       blockPtr => domain % blocklist
       do while (associated(blockPtr))

          call mpas_pool_get_subpool(blockPtr % structs, "prescribed_ice", prescribedIcePool)
          call mpas_pool_get_subpool(blockPtr % structs, "tracers", tracersPool)
          call mpas_pool_get_subpool(blockPtr % structs, "initial", initialPool)
          call mpas_pool_get_subpool(blockPtr % structs, "ocean_coupling", oceanCouplingPool)
          call mpas_pool_get_subpool(blockPtr % structs, "mesh", meshPool)

          call mpas_pool_get_array(prescribedIcePool, "iceCoverage", iceCoverage)

          call mpas_pool_get_array(tracersPool, "iceAreaCategory", iceAreaCategory, 1)
          call mpas_pool_get_array(tracersPool, "iceVolumeCategory", iceVolumeCategory, 1)
          call mpas_pool_get_array(tracersPool, "snowVolumeCategory", snowVolumeCategory, 1)
          call mpas_pool_get_array(tracersPool, "surfaceTemperature", surfaceTemperature, 1)
          call mpas_pool_get_array(tracersPool, "iceEnthalpy", iceEnthalpy, 1)
          call mpas_pool_get_array(tracersPool, "iceSalinity", iceSalinity, 1)
          call mpas_pool_get_array(tracersPool, "snowEnthalpy", snowEnthalpy, 1)

          call mpas_pool_get_array(initialPool, "categoryThicknessLimits", categoryThicknessLimits)

          call mpas_pool_get_array(oceanCouplingPool, "seaFreezingTemperature", seaFreezingTemperature)

          call mpas_pool_get_array(meshPool, "latCell", latCell)

          call mpas_pool_get_dimension(blockPtr % dimensions, "nCells", nCells)
          call mpas_pool_get_dimension(blockPtr % dimensions, "nCategories", nCategories)
          call mpas_pool_get_dimension(blockPtr % dimensions, "nIceLayers", nIceLayers)
          call mpas_pool_get_dimension(blockPtr % dimensions, "nSnowLayers", nSnowLayers)

          do iCell = 1, nCells

             ! limit ice coverage
             if (iceCoverage(iCell) < 1e-4_RKIND) iceCoverage(iCell) = 0.0_RKIND
             if (iceCoverage(iCell) > 1.0_RKIND)  iceCoverage(iCell) = 1.0_RKIND

             if (iceCoverage(iCell) >= 1e-4_RKIND) then

                ! set thickness based on hemisphere
                if (latCell(iCell) >= 0.0) then
                   iceThickness = 2.0_RKIND
                else if (latCell(iCell) < 0.0) then
                   iceThickness = 1.0_RKIND
                endif

                ! loop over categories
                do iCategory = 1, nCategories

                   ! find relevant thickness category
                   if (iceThickness >= categoryThicknessLimits(iCategory) .and. & !!!! index limits
                       iceThickness <  categoryThicknessLimits(iCategory+1)) then

                      if (iceAreaCategory(1,iCategory,iCell) > 0.0_RKIND) then
                         snowThickness = snowVolumeCategory(1,iCategory,iCell) / iceAreaCategory(1,iCategory,iCell)
                      else
                         snowThickness = 0.0_RKIND
                      endif ! prognostic ice present

                      ! reset ice state
                      iceAreaCategory(1,iCategory,iCell)    = iceCoverage(iCell)
                      iceVolumeCategory(1,iCategory,iCell)  = iceThickness  * iceAreaCategory(1,iCategory,iCell)
                      snowVolumeCategory(1,iCategory,iCell) = snowThickness * iceAreaCategory(1,iCategory,iCell)

                      ! set thermodynamic tracers
                      if (abs(iceEnthalpy(1,iCategory,iCell)) < seaicePuny) then !!!!! check indices

                         if (iceAreaCategory(1,iCategory,iCell) < seaicePuny) then
                            surfaceTemperature(1,iCategory,iCell) = seaFreezingTemperature(iCell)
                         endif

                         temperatureGradient = seaFreezingTemperature(iCell) - surfaceTemperature(1,iCategory,iCell)

                         if (trim(config_column_physics_type) == "icepack") then

                            ! ice quantities
                            do iIceLayer = 1, nIceLayers

                               depth = (real(iIceLayer,kind=RKIND) - 0.5_RKIND) / real(nIceLayers,kind=RKIND)
                               iceTemperature = surfaceTemperature(1,iCategory,iCell) + temperatureGradient * depth
                               iceSalinity(iIceLayer,iCategory,iCell) = icepack_salinity_profile(depth)
                               iceEnthalpy(iIceLayer,iCategory,iCell) = seaice_icepack_enthalpy_ice(iceTemperature,iceSalinity(iIceLayer,iCategory,iCell),config_thermodynamics_type)

                            enddo ! iIceLayer

                            ! snow quantities
                            do iSnowLayer = 1, nSnowLayers
                               snowEnthalpy(iSnowLayer,iCategory,iCell) = icepack_enthalpy_snow(surfaceTemperature(1,iCategory,iCell))
                            enddo ! iSnowLayer

                         else if (trim(config_column_physics_type) == "column_package") then

                            ! ice quantities
                            do iIceLayer = 1, nIceLayers

                               depth = (real(iIceLayer,kind=RKIND) - 0.5_RKIND) / real(nIceLayers,kind=RKIND)
                               iceTemperature = surfaceTemperature(1,iCategory,iCell) + temperatureGradient * depth
                               iceSalinity(iIceLayer,iCategory,iCell) = colpkg_salinity_profile(depth)
                               iceEnthalpy(iIceLayer,iCategory,iCell) = colpkg_enthalpy_ice(iceTemperature,iceSalinity(iIceLayer,iCategory,iCell))

                            enddo ! iIceLayer

                            ! snow quantities
                            do iSnowLayer = 1, nSnowLayers
                               snowEnthalpy(iSnowLayer,iCategory,iCell) = colpkg_enthalpy_snow(surfaceTemperature(1,iCategory,iCell))
                            enddo ! iSnowLayer

                         endif ! config_column_physics_type

                      endif

                   else

                      surfaceTemperature(1,iCategory,iCell) = seaFreezingTemperature(iCell)
                      iceAreaCategory(1,iCategory,iCell)    = 0.0_RKIND
                      iceVolumeCategory(1,iCategory,iCell)  = 0.0_RKIND
                      snowVolumeCategory(1,iCategory,iCell) = 0.0_RKIND
                      iceSalinity(:,iCategory,iCell)        = 0.0_RKIND
                      iceEnthalpy(:,iCategory,iCell)        = 0.0_RKIND
                      snowEnthalpy(:,iCategory,iCell)       = 0.0_RKIND

                   endif ! in category

                enddo ! iCategory

             else

                ! low ice coverage
                surfaceTemperature(1,:,iCell) = seaFreezingTemperature(iCell)
                iceAreaCategory(1,:,iCell)    = 0.0_RKIND
                iceVolumeCategory(1,:,iCell)  = 0.0_RKIND
                snowVolumeCategory(1,:,iCell) = 0.0_RKIND
                iceSalinity(:,:,iCell)        = 0.0_RKIND
                iceEnthalpy(:,:,iCell)        = 0.0_RKIND
                snowEnthalpy(:,:,iCell)       = 0.0_RKIND

             endif ! ice coverage

          enddo ! iCell

          blockPtr => blockPtr % next
       enddo

       ! aggregate tracers
       if (trim(config_column_physics_type) == "icepack") then
          call seaice_icepack_aggregate(domain)
       else if (trim(config_column_physics_type) == "column_package") then
          call seaice_column_aggregate(domain)
       endif ! config_column_physics_type

       ! set non-computed fluxes, ice velocities, ice-ocn stresses to zero
       blockPtr => domain % blocklist
       do while (associated(blockPtr))

          call mpas_pool_get_subpool(blockPtr % structs, "velocity_solver", velocitySolverPool)
          call mpas_pool_get_subpool(blockPtr % structs, "ocean_coupling", oceanCouplingPool)

          call mpas_pool_get_array(oceanCouplingPool, "freezingMeltingPotential", freezingMeltingPotential)
          call mpas_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
          call mpas_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)
          call mpas_pool_get_array(velocitySolverPool, "oceanStressCellU", oceanStressCellU)
          call mpas_pool_get_array(velocitySolverPool, "oceanStressCellV", oceanStressCellV)

          freezingMeltingPotential(:) = 0.0_RKIND
          uVelocity(:) = 0.0_RKIND
          vVelocity(:) = 0.0_RKIND
          oceanStressCellU(:) = 0.0_RKIND
          oceanStressCellV(:) = 0.0_RKIND

          blockPtr => blockPtr % next
       enddo

       ! reinitialize fluxes
       if (trim(config_column_physics_type) == "icepack") then
          call seaice_icepack_reinitialize_fluxes(domain)
       else if (trim(config_column_physics_type) == "column_package") then
          call seaice_column_reinitialize_fluxes(domain)
       endif ! config_column_physics_type

    endif ! prescribed ice mode

  end subroutine seaice_run_prescribed_ice

  !-----------------------------------------------------------------------

end module seaice_prescribed
