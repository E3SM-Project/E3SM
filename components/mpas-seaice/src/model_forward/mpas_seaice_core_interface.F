module seaice_core_interface

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_io_units
   use seaice_core
   use seaice_analysis_driver
   use mpas_log, only: mpas_log_write

   implicit none
   public

   contains

   !***********************************************************************
   !
   !  routine seaice_setup_core
   !
   !> \brief   Seaice core setup routine
   !> \author  Doug Jacobsen
   !> \date    03/18/2015
   !> \details
   !>  This routine is intended to setup the necessary variables within a core_type
   !>  for the ocean core.
   !
   !-----------------------------------------------------------------------
   subroutine seaice_setup_core(core)!{{{
      type (core_type), pointer :: core

      core % core_init => seaice_core_init
      core % core_run => seaice_core_run
      core % core_finalize => seaice_core_finalize
      core % define_packages => seaice_define_packages
      core % setup_packages => seaice_setup_packages
      core % setup_decompositions => seaice_setup_decompositions
      core % setup_clock => seaice_setup_clock
      core % setup_log => seaice_setup_log
      core % get_mesh_stream => seaice_get_mesh_stream
      core % setup_immutable_streams => seaice_setup_immutable_streams
      core % setup_derived_dimensions => seaice_setup_derived_dimensions
      core % setup_decomposed_dimensions => seaice_setup_decomposed_dimensions
      core % setup_block => seaice_setup_block
      core % setup_namelist => seaice_setup_namelists

      core % Conventions = 'MPAS'
      core % source = 'MPAS'
#include "../inc/core_variables.inc"

   end subroutine seaice_setup_core!}}}

   !***********************************************************************
   !
   !  routine seaice_setup_domain
   !
   !> \brief   Seaice domain setup routine
   !> \author  Doug Jacobsen
   !> \date    03/18/2015
   !> \details
   !>  This routine is intended to setup the necessary variables within a domain_type
   !>  for the ocean core.
   !
   !-----------------------------------------------------------------------
   subroutine seaice_setup_domain(domain)!{{{
      type (domain_type), pointer :: domain

#include "../inc/domain_variables.inc"

   end subroutine seaice_setup_domain!}}}

   !***********************************************************************
   !
   !  function seaice_setup_packages
   !
   !> \brief   Pacakge setup routine
   !> \author  Doug Jacobsen
   !> \date    03/12/2015
   !> \details
   !>  This function is intended to correctly configure the packages for this MPAS
   !>   core. It can use any Fortran logic to properly configure packages, and it
   !>   can also make use of any namelist options. All variables in the model are
   !>   *not* allocated until after this routine is called.
   !
   !-----------------------------------------------------------------------
   function seaice_setup_packages(configPool, packagePool, iocontext) result(ierr)!{{{

      use mpas_derived_types

      implicit none

      type (mpas_pool_type), intent(inout) :: configPool
      type (mpas_pool_type), intent(inout) :: packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext

      type (mpas_pool_iterator_type) :: pkgItr
      logical, pointer :: packageActive

      integer :: ierr

      ierr = 0

      ! dynamics
      call setup_packages_dynamics(configPool, packagePool, ierr)

      ! column physics
      call setup_packages_column_physics(configPool, packagePool, ierr)

      ! set up analysis member packages
      call seaice_analysis_setup_packages(configPool, packagePool, iocontext, ierr)

      ! icebergs
      call setup_packages_bergs(configPool, packagePool, ierr)

      ! testing system test
      call setup_packages_other(configPool, packagePool, ierr)

      call mpas_log_write( '')
      call mpas_log_write( '  **** Summary of sea ice packages ****')
      call mpas_pool_begin_iteration(packagePool)
      do while ( mpas_pool_get_next_member(packagePool, pkgItr) )

         if ( pkgItr % memberType == MPAS_POOL_PACKAGE ) then
            call mpas_pool_get_package(packagePool, pkgItr % memberName, packageActive)
            if ( packageActive ) then
               call mpas_log_write( '      ' // trim(pkgItr % memberName) // ' = ON')
            else
               call mpas_log_write( '      ' // trim(pkgItr % memberName) // ' = OFF')
            end if
         end if
      end do
      call mpas_log_write( '  ***********************************')
      call mpas_log_write( '')

   end function seaice_setup_packages!}}}

   !***********************************************************************
   !
   !  routine setup_packages_dynamics
   !
   !> \brief
   !> \author  Adrian K. Turner
   !> \date    9th Feburary 2015
   !> \details
   !>
   !
   !-----------------------------------------------------------------------

   subroutine setup_packages_dynamics(configPool, packagePool, ierr)!{{{

     type (mpas_pool_type), intent(in) :: configPool
     type (mpas_pool_type), intent(in) :: packagePool
     integer, intent(out) :: ierr

     logical, pointer :: &
          config_use_velocity_solver

     character(len=strKIND), pointer :: &
          config_strain_scheme, &
          config_stress_divergence_scheme

     logical, pointer :: &
          pkgWeakActive, &
          pkgVariationalActive, &
          pkgWeakVariationalActive

     !pkgWeak
     !pkgVariational

     call MPAS_pool_get_config(configPool, "config_use_velocity_solver", config_use_velocity_solver)
     call MPAS_pool_get_config(configPool, "config_strain_scheme", config_strain_scheme)
     call MPAS_pool_get_config(configPool, "config_stress_divergence_scheme", config_stress_divergence_scheme)

     call MPAS_pool_get_package(packagePool, "pkgWeakActive", pkgWeakActive)
     call MPAS_pool_get_package(packagePool, "pkgVariationalActive", pkgVariationalActive)
     call MPAS_pool_get_package(packagePool, "pkgWeakVariationalActive", pkgWeakVariationalActive)

     ! is the velocity solver on
     if (config_use_velocity_solver) then

        ! stress divergence scheme
        if (trim(config_strain_scheme)            == "weak" .or. &
            trim(config_stress_divergence_scheme) == "weak") then

           ! weak stress divergence scheme
           pkgWeakActive = .true.

        endif

        if (trim(config_strain_scheme)            == "variational" .or. &
            trim(config_stress_divergence_scheme) == "variational") then

           ! variational stress divergence scheme
           pkgVariationalActive = .true.

        endif

        if (trim(config_strain_scheme)            == "weak" .and. &
            trim(config_stress_divergence_scheme) == "variational") then

           ! weak strain / variational stress divergence
           pkgWeakVariationalActive = .true.

        endif

     endif ! config_use_velocity_solver

   end subroutine setup_packages_dynamics!}}}

   !***********************************************************************
   !
   !  routine setup_packages_dynamics
   !
   !> \brief
   !> \author  Adrian K. Turner
   !> \date    9th Feburary 2015
   !> \details
   !>
   !
   !-----------------------------------------------------------------------

   subroutine setup_packages_column_physics(configPool, packagePool, ierr)!{{{

     type (mpas_pool_type), intent(in) :: configPool
     type (mpas_pool_type), intent(in) :: packagePool
     integer, intent(out) :: ierr

     ! column physics package packages
     logical, pointer :: &
          config_use_column_physics, &
          config_use_column_biogeochemistry

     logical, pointer :: &
          pkgColumnPackageActive, &
          pkgColumnBiogeochemistryActive

     ! column tracers
     logical, pointer :: &
          config_use_ice_age, &
          config_use_first_year_ice, &
          config_use_level_ice, &
          config_use_cesm_meltponds, &
          config_use_level_meltponds, &
          config_use_topo_meltponds, &
          config_use_aerosols, &
          config_use_brine, &
          config_use_vertical_zsalinity, &
          config_use_vertical_biochemistry, &
          config_use_vertical_tracers, &
          config_use_skeletal_biochemistry, &
          config_use_nitrate, &
          config_use_carbon, &
          config_use_chlorophyll, &
          config_use_ammonium, &
          config_use_silicate, &
          config_use_DMS, &
          config_use_nonreactive, &
          config_use_humics, &
          config_use_DON, &
          config_use_iron, &
          config_use_zaerosols, &
          config_use_snow_grain_radius, &
          config_use_effective_snow_density

     logical, pointer :: &
          pkgColumnTracerIceAgeActive, &
          pkgColumnTracerFirstYearIceActive, &
          pkgColumnTracerLevelIceActive, &
          pkgColumnTracerPondsActive, &
          pkgColumnTracerLidThicknessActive, &
          pkgColumnTracerAerosolsActive, &
          pkgTracerBrineActive, &
          pkgTracerMobileFractionActive, &
          pkgTracerSkeletalAlgaeActive, &
          pkgTracerSkeletalNitrateActive, &
          pkgTracerSkeletalCarbonActive, &
          pkgTracerSkeletalAmmoniumActive, &
          pkgTracerSkeletalSilicateActive, &
          pkgTracerSkeletalDMSActive, &
          pkgTracerSkeletalNonreactiveActive, &
          pkgTracerSkeletalHumicsActive, &
          pkgTracerSkeletalDONActive, &
          pkgTracerSkeletalIronActive, &
          pkgTracerVerticalAlgaeActive, &
          pkgTracerVerticalNitrateActive, &
          pkgTracerVerticalCarbonActive, &
          pkgTracerVerticalAmmoniumActive, &
          pkgTracerVerticalSilicateActive, &
          pkgTracerVerticalDMSActive, &
          pkgTracerVerticalNonreactiveActive, &
          pkgTracerVerticalHumicsActive, &
          pkgTracerVerticalDONActive, &
          pkgTracerVerticalIronActive, &
          pkgTracerZAerosolsActive, &
          pkgTracerZSalinityActive, &
          pkgColumnTracerEffectiveSnowDensityActive, &
          pkgColumnTracerSnowGrainRadiusActive

     ! other packages
     logical, pointer :: &
          config_use_form_drag

     logical, pointer :: &
          pkgColumnFormDragActive

     logical :: &
          use_meltponds

     ierr = 0

     !-----------------------------------------------------------------------
     ! column main routines
     !-----------------------------------------------------------------------

     !pkgColumnPackage
     !pkgColumnBiogeochemistry

     call MPAS_pool_get_config(configPool, "config_use_column_physics", config_use_column_physics)
     call MPAS_pool_get_config(configPool, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)
     call MPAS_pool_get_config(configPool, "config_use_zaerosols", config_use_zaerosols)

     call MPAS_pool_get_package(packagePool, "pkgColumnPackageActive", pkgColumnPackageActive)
     call MPAS_pool_get_package(packagePool, "pkgColumnBiogeochemistryActive", pkgColumnBiogeochemistryActive)

     pkgColumnPackageActive         = config_use_column_physics
     pkgColumnBiogeochemistryActive = ((config_use_column_biogeochemistry .or. config_use_zaerosols) &
                                    .and. config_use_column_physics)

     !pkgColumnPackageActive         = .true.
     !pkgColumnBiogeochemistryActive = .true.

     !-----------------------------------------------------------------------
     ! tracers
     !-----------------------------------------------------------------------

     !pkgColumnTracerIceAge
     !pkgColumnTracerFirstYearIce
     !pkgColumnTracerLevelIce
     !pkgColumnTracerPonds
     !pkgColumnTracerLidThickness
     !pkgColumnTracerAerosols
     !pkgTracerBrine
     !pkgTracerMobileFraction
     !pkgTracerSkeletalAlgae
     !pkgTracerSkeletalNitrate
     !pkgTracerSkeletalCarbon
     !pkgTracerSkeletalAmmonium
     !pkgTracerSkeletalSilicate
     !pkgTracerSkeletalDMS
     !pkgTracerSkeletalNonreactive
     !pkgTracerSkeletalHumics
     !pkgTracerSkeletalDON
     !pkgTracerSkeletalIron
     !pkgTracerVerticalAlgae
     !pkgTracerVerticalNitrate
     !pkgTracerVerticalCarbon
     !pkgTracerVerticalAmmonium
     !pkgTracerVerticalSilicate
     !pkgTracerVerticalDMS
     !pkgTracerVerticalNonreactive
     !pkgTracerVerticalHumics
     !pkgTracerVerticalDON
     !pkgTracerVerticalIron
     !pkgTracerZAerosols
     !pkgTracerZSalinity
     !pkgColumnTracerEffectiveSnowDensity
     !pkgColumnTracerSnowGrainRadius

     call MPAS_pool_get_config(configPool, "config_use_ice_age", config_use_ice_age)
     call MPAS_pool_get_config(configPool, "config_use_first_year_ice", config_use_first_year_ice)
     call MPAS_pool_get_config(configPool, "config_use_level_ice", config_use_level_ice)
     call MPAS_pool_get_config(configPool, "config_use_cesm_meltponds", config_use_cesm_meltponds)
     call MPAS_pool_get_config(configPool, "config_use_level_meltponds", config_use_level_meltponds)
     call MPAS_pool_get_config(configPool, "config_use_topo_meltponds", config_use_topo_meltponds)
     call MPAS_pool_get_config(configPool, "config_use_aerosols", config_use_aerosols)
     call MPAS_pool_get_config(configPool, "config_use_brine", config_use_brine)
     call MPAS_pool_get_config(configPool, "config_use_vertical_zsalinity", config_use_vertical_zsalinity)
     call MPAS_pool_get_config(configPool, "config_use_vertical_tracers", config_use_vertical_tracers)
     call MPAS_pool_get_config(configPool, "config_use_skeletal_biochemistry", config_use_skeletal_biochemistry)
     call MPAS_pool_get_config(configPool, "config_use_vertical_biochemistry", config_use_vertical_biochemistry)
     call MPAS_pool_get_config(configPool, "config_use_nitrate", config_use_nitrate)
     call MPAS_pool_get_config(configPool, "config_use_carbon", config_use_carbon)
     call MPAS_pool_get_config(configPool, "config_use_chlorophyll", config_use_chlorophyll)
     call MPAS_pool_get_config(configPool, "config_use_ammonium", config_use_ammonium)
     call MPAS_pool_get_config(configPool, "config_use_silicate", config_use_silicate)
     call MPAS_pool_get_config(configPool, "config_use_DMS", config_use_DMS)
     call MPAS_pool_get_config(configPool, "config_use_nonreactive", config_use_nonreactive)
     call MPAS_pool_get_config(configPool, "config_use_humics", config_use_humics)
     call MPAS_pool_get_config(configPool, "config_use_DON", config_use_DON)
     call MPAS_pool_get_config(configPool, "config_use_iron", config_use_iron)
     call MPAS_pool_get_config(configPool, "config_use_effective_snow_density", config_use_effective_snow_density)
     call MPAS_pool_get_config(configPool, "config_use_snow_grain_radius", config_use_snow_grain_radius)

     call MPAS_pool_get_package(packagePool, "pkgColumnTracerIceAgeActive", pkgColumnTracerIceAgeActive)
     call MPAS_pool_get_package(packagePool, "pkgColumnTracerFirstYearIceActive", pkgColumnTracerFirstYearIceActive)
     call MPAS_pool_get_package(packagePool, "pkgColumnTracerLevelIceActive", pkgColumnTracerLevelIceActive)
     call MPAS_pool_get_package(packagePool, "pkgColumnTracerPondsActive", pkgColumnTracerPondsActive)
     call MPAS_pool_get_package(packagePool, "pkgColumnTracerLidThicknessActive", pkgColumnTracerLidThicknessActive)
     call MPAS_pool_get_package(packagePool, "pkgColumnTracerAerosolsActive", pkgColumnTracerAerosolsActive)
     call MPAS_pool_get_package(packagePool, "pkgTracerBrineActive", pkgTracerBrineActive)
     call MPAS_pool_get_package(packagePool, "pkgTracerMobileFractionActive", pkgTracerMobileFractionActive)
     call MPAS_pool_get_package(packagePool, "pkgTracerSkeletalAlgaeActive", pkgTracerSkeletalAlgaeActive)
     call MPAS_pool_get_package(packagePool, "pkgTracerSkeletalNitrateActive", pkgTracerSkeletalNitrateActive)
     call MPAS_pool_get_package(packagePool, "pkgTracerSkeletalCarbonActive", pkgTracerSkeletalCarbonActive)
     call MPAS_pool_get_package(packagePool, "pkgTracerSkeletalAmmoniumActive", pkgTracerSkeletalAmmoniumActive)
     call MPAS_pool_get_package(packagePool, "pkgTracerSkeletalSilicateActive", pkgTracerSkeletalSilicateActive)
     call MPAS_pool_get_package(packagePool, "pkgTracerSkeletalDMSActive", pkgTracerSkeletalDMSActive)
     call MPAS_pool_get_package(packagePool, "pkgTracerSkeletalNonreactiveActive", pkgTracerSkeletalNonreactiveActive)
     call MPAS_pool_get_package(packagePool, "pkgTracerSkeletalHumicsActive", pkgTracerSkeletalHumicsActive)
     call MPAS_pool_get_package(packagePool, "pkgTracerSkeletalDONActive", pkgTracerSkeletalDONActive)
     call MPAS_pool_get_package(packagePool, "pkgTracerSkeletalIronActive", pkgTracerSkeletalIronActive)
     call MPAS_pool_get_package(packagePool, "pkgTracerVerticalAlgaeActive", pkgTracerVerticalAlgaeActive)
     call MPAS_pool_get_package(packagePool, "pkgTracerVerticalNitrateActive", pkgTracerVerticalNitrateActive)
     call MPAS_pool_get_package(packagePool, "pkgTracerVerticalCarbonActive", pkgTracerVerticalCarbonActive)
     call MPAS_pool_get_package(packagePool, "pkgTracerVerticalAmmoniumActive", pkgTracerVerticalAmmoniumActive)
     call MPAS_pool_get_package(packagePool, "pkgTracerVerticalSilicateActive", pkgTracerVerticalSilicateActive)
     call MPAS_pool_get_package(packagePool, "pkgTracerVerticalDMSActive", pkgTracerVerticalDMSActive)
     call MPAS_pool_get_package(packagePool, "pkgTracerVerticalNonreactiveActive", pkgTracerVerticalNonreactiveActive)
     call MPAS_pool_get_package(packagePool, "pkgTracerVerticalHumicsActive", pkgTracerVerticalHumicsActive)
     call MPAS_pool_get_package(packagePool, "pkgTracerVerticalDONActive", pkgTracerVerticalDONActive)
     call MPAS_pool_get_package(packagePool, "pkgTracerVerticalIronActive", pkgTracerVerticalIronActive)
     call MPAS_pool_get_package(packagePool, "pkgTracerZAerosolsActive", pkgTracerZAerosolsActive)
     call MPAS_pool_get_package(packagePool, "pkgTracerZSalinityActive", pkgTracerZSalinityActive)
     call MPAS_pool_get_package(packagePool, "pkgColumnTracerEffectiveSnowDensityActive", pkgColumnTracerEffectiveSnowDensityActive)
     call MPAS_pool_get_package(packagePool, "pkgColumnTracerSnowGrainRadiusActive", pkgColumnTracerSnowGrainRadiusActive)

     use_meltponds = (config_use_cesm_meltponds .or. config_use_level_meltponds .or. config_use_topo_meltponds)

     pkgColumnTracerIceAgeActive       = config_use_ice_age
     pkgColumnTracerFirstYearIceActive = config_use_first_year_ice
     pkgColumnTracerLevelIceActive     = config_use_level_ice
     pkgColumnTracerPondsActive        = use_meltponds
     pkgColumnTracerLidThicknessActive = (config_use_level_meltponds .or. config_use_topo_meltponds)
     pkgColumnTracerAerosolsActive     = config_use_aerosols

     pkgTracerBrineActive                = config_use_brine
     pkgTracerMobileFractionActive       = config_use_vertical_tracers
     pkgTracerSkeletalAlgaeActive        = config_use_skeletal_biochemistry
     pkgTracerSkeletalNitrateActive      = (config_use_skeletal_biochemistry .and. config_use_nitrate)
     pkgTracerSkeletalCarbonActive       = (config_use_skeletal_biochemistry .and. config_use_carbon)
     pkgTracerSkeletalAmmoniumActive     = (config_use_skeletal_biochemistry .and. config_use_ammonium)
     pkgTracerSkeletalSilicateActive     = (config_use_skeletal_biochemistry .and. config_use_silicate)
     pkgTracerSkeletalDMSActive          = (config_use_skeletal_biochemistry .and. config_use_DMS)
     pkgTracerSkeletalNonreactiveActive  = (config_use_skeletal_biochemistry .and. config_use_nonreactive)
     pkgTracerSkeletalHumicsActive       = (config_use_skeletal_biochemistry .and. config_use_humics)
     pkgTracerSkeletalDONActive          = (config_use_skeletal_biochemistry .and. config_use_DON)
     pkgTracerSkeletalIronActive         = (config_use_skeletal_biochemistry .and. config_use_iron)
     pkgTracerVerticalAlgaeActive        = (config_use_vertical_tracers .and. config_use_vertical_biochemistry)
     pkgTracerVerticalNitrateActive      = (config_use_vertical_tracers .and. config_use_nitrate)
     pkgTracerVerticalCarbonActive       = (config_use_vertical_tracers .and. config_use_carbon)
     pkgTracerVerticalAmmoniumActive     = (config_use_vertical_tracers .and. config_use_ammonium)
     pkgTracerVerticalSilicateActive     = (config_use_vertical_tracers .and. config_use_silicate)
     pkgTracerVerticalDMSActive          = (config_use_vertical_tracers .and. config_use_DMS)
     pkgTracerVerticalNonreactiveActive  = (config_use_vertical_tracers .and. config_use_nonreactive)
     pkgTracerVerticalHumicsActive       = (config_use_vertical_tracers .and. config_use_humics)
     pkgTracerVerticalDONActive          = (config_use_vertical_tracers .and. config_use_DON)
     pkgTracerVerticalIronActive         = (config_use_vertical_tracers .and. config_use_iron)
     pkgTracerZAerosolsActive            = config_use_zaerosols
     pkgTracerZSalinityActive            = config_use_vertical_zsalinity

     pkgColumnTracerEffectiveSnowDensityActive = config_use_effective_snow_density
     pkgColumnTracerSnowGrainRadiusActive = config_use_snow_grain_radius

     if (.not. config_use_column_physics) then
        pkgColumnTracerIceAgeActive       = .false.
        pkgColumnTracerFirstYearIceActive = .false.
        pkgColumnTracerLevelIceActive     = .false.
        pkgColumnTracerPondsActive        = .false.
        pkgColumnTracerLidThicknessActive = .false.
        pkgColumnTracerAerosolsActive     = .false.
        pkgTracerBrineActive                = .false.
        pkgTracerMobileFractionActive       = .false.
        pkgTracerSkeletalAlgaeActive        = .false.
        pkgTracerSkeletalNitrateActive      = .false.
        pkgTracerSkeletalCarbonActive       = .false.
        pkgTracerSkeletalAmmoniumActive     = .false.
        pkgTracerSkeletalSilicateActive     = .false.
        pkgTracerSkeletalDMSActive          = .false.
        pkgTracerSkeletalNonreactiveActive  = .false.
        pkgTracerSkeletalHumicsActive       = .false.
        pkgTracerSkeletalDONActive          = .false.
        pkgTracerSkeletalIronActive         = .false.
        pkgTracerVerticalAlgaeActive        = .false.
        pkgTracerVerticalNitrateActive      = .false.
        pkgTracerVerticalCarbonActive       = .false.
        pkgTracerVerticalAmmoniumActive     = .false.
        pkgTracerVerticalSilicateActive     = .false.
        pkgTracerVerticalDMSActive          = .false.
        pkgTracerVerticalNonreactiveActive  = .false.
        pkgTracerVerticalHumicsActive       = .false.
        pkgTracerVerticalDONActive          = .false.
        pkgTracerVerticalIronActive         = .false.
        pkgTracerZAerosolsActive            = .false.
        pkgTracerZSalinityActive            = .false.
        pkgColumnTracerEffectiveSnowDensityActive = .false.
        pkgColumnTracerSnowGrainRadiusActive      = .false.
     endif

     if (.not. config_use_column_biogeochemistry .and. config_use_column_physics) then
        pkgTracerSkeletalAlgaeActive        = .false.
        pkgTracerSkeletalNitrateActive      = .false.
        pkgTracerSkeletalCarbonActive       = .false.
        pkgTracerSkeletalAmmoniumActive     = .false.
        pkgTracerSkeletalSilicateActive     = .false.
        pkgTracerSkeletalDMSActive          = .false.
        pkgTracerSkeletalNonreactiveActive  = .false.
        pkgTracerSkeletalHumicsActive       = .false.
        pkgTracerSkeletalDONActive          = .false.
        pkgTracerSkeletalIronActive         = .false.
        pkgTracerVerticalAlgaeActive        = .false.
        pkgTracerVerticalNitrateActive      = .false.
        pkgTracerVerticalCarbonActive       = .false.
        pkgTracerVerticalAmmoniumActive     = .false.
        pkgTracerVerticalSilicateActive     = .false.
        pkgTracerVerticalDMSActive          = .false.
        pkgTracerVerticalNonreactiveActive  = .false.
        pkgTracerVerticalHumicsActive       = .false.
        pkgTracerVerticalDONActive          = .false.
        pkgTracerVerticalIronActive         = .false.
        pkgTracerZSalinityActive            = .false.
     endif

     !pkgColumnTracerIceAgeActive       = .true.
     !pkgColumnTracerFirstYearIceActive = .true.
     !pkgColumnTracerLevelIceActive     = .true.
     !pkgColumnTracerPondsActive        = .true.
     !pkgColumnTracerLidThicknessActive = .true.
     !pkgColumnTracerAerosolsActive     = .true.
     !pkgColumnTracerArrayStandInActive = .true.

     !-----------------------------------------------------------------------
     ! other column packages
     !-----------------------------------------------------------------------

     !pkgColumnFormDrag

     ! form drag
     call MPAS_pool_get_config(configPool, "config_use_form_drag", config_use_form_drag)
     call MPAS_pool_get_package(packagePool, "pkgColumnFormDragActive", pkgColumnFormDragActive)
     pkgColumnFormDragActive = (config_use_column_physics .and. config_use_form_drag)

     !pkgColumnFormDragActive = .true.

   end subroutine setup_packages_column_physics!}}}

   !***********************************************************************
   !
   !  routine setup_packages_bergs
   !
   !> \brief   Setup icebergs package
   !> \author  Darin Comeau, Xylar Asay-Davis
   !> \date    August 2024
   !> \details This routine is intended to set the icebergs packages PkgBergs
   !> and pkgScaledDIB as active/deactive based on the namelist options.
   !
   !-----------------------------------------------------------------------

   subroutine setup_packages_bergs(configPool, packagePool, ierr)!{{{

     type (mpas_pool_type), intent(in) :: configPool
     type (mpas_pool_type), intent(in) :: packagePool
     integer, intent(out) :: ierr

          ! icebergs package
     logical, pointer :: &
          config_use_data_icebergs, &
          config_scale_dib_by_removed_ice_runoff


     logical, pointer :: &
          pkgBergsActive, &
          pkgScaledDIBActive

     ierr = 0

     !-----------------------------------------------------------------------
     ! iceberg routines
     !-----------------------------------------------------------------------

     call MPAS_pool_get_config(configPool, "config_use_data_icebergs", config_use_data_icebergs)
     call MPAS_pool_get_package(packagePool, "pkgBergsActive", pkgBergsActive)
     pkgBergsActive = config_use_data_icebergs

     call MPAS_pool_get_config(configPool, "config_scale_dib_by_removed_ice_runoff", config_scale_dib_by_removed_ice_runoff)
     call MPAS_pool_get_package(packagePool, "pkgScaledDIBActive", pkgScaledDIBActive)
     pkgScaledDIBActive = config_scale_dib_by_removed_ice_runoff

   end subroutine setup_packages_bergs!}}}

   !***********************************************************************
   !
   !  routine setup_packages_other
   !
   !> \brief
   !> \author  Adrian K. Turner
   !> \date    9th Feburary 2015
   !> \details
   !>
   !
   !-----------------------------------------------------------------------

   subroutine setup_packages_other(configPool, packagePool, ierr)!{{{

     type (mpas_pool_type), intent(in) :: configPool
     type (mpas_pool_type), intent(in) :: packagePool
     integer, intent(out) :: ierr

     logical, pointer :: &
          config_use_forcing, &
          config_use_data_icebergs, &
          config_use_zaerosols, &
          config_use_iron_solubility_file, &
          config_use_atm_dust_file, &
          config_testing_system_test, &
          config_use_snicar_ad, &
          config_use_prescribed_ice, &
          config_use_special_boundaries_velocity, &
          config_use_special_boundaries_velocity_masks, &
          config_use_special_boundaries_tracers

     logical, pointer :: &
          pkgForcingActive, &
          pkgTestingSystemTestActive, &
          pkgSnicarActive, &
          pkgPrescribedActive, &
          pkgSpecialBoundariesActive

     logical :: &
          config_use_dust_file

     ierr = 0

     ! pkgForcing

     call MPAS_pool_get_config(configPool, "config_use_forcing", config_use_forcing)
     call MPAS_pool_get_config(configPool, "config_use_data_icebergs", config_use_data_icebergs)
     call MPAS_pool_get_config(configPool, "config_use_zaerosols", config_use_zaerosols)
     call MPAS_pool_get_config(configPool, "config_use_atm_dust_file", config_use_atm_dust_file)
     call MPAS_pool_get_config(configPool, "config_use_iron_solubility_file", config_use_iron_solubility_file)

     call MPAS_pool_get_package(packagePool, "pkgForcingActive", pkgForcingActive)

     config_use_dust_file = .false.
     if (config_use_zaerosols .and. (config_use_iron_solubility_file .or. &
         config_use_atm_dust_file) ) config_use_dust_file = .true.

     ! see if we are using the forcing system
     if (config_use_forcing .or. config_use_data_icebergs .or. &
         config_use_dust_file) then

        pkgForcingActive = .true.

     endif

     !pkgTestingSystemTest

     call MPAS_pool_get_config(configPool, "config_testing_system_test", config_testing_system_test)

     call MPAS_pool_get_package(packagePool, "pkgTestingSystemTestActive", pkgTestingSystemTestActive)

     ! see if we are testing the testing system
     if (config_testing_system_test) then
        pkgTestingSystemTestActive = .true.
     endif

     ! pkgSnicar
     call MPAS_pool_get_config(configPool, "config_use_snicar_ad", config_use_snicar_ad)
     call MPAS_pool_get_package(packagePool, "pkgSnicarActive", pkgSnicarActive)
     ! see if we are using the snicar_ad system
     pkgSnicarActive = .true.

     ! pkgPrescribed
     call MPAS_pool_get_config(configPool, "config_use_prescribed_ice", config_use_prescribed_ice)
     call MPAS_pool_get_package(packagePool, "pkgPrescribedActive", pkgPrescribedActive)
     if (config_use_prescribed_ice) then
        pkgPrescribedActive = .true.
     endif

     ! pkgSpecialBoundaries
     call MPAS_pool_get_config(configPool, "config_use_special_boundaries_velocity", config_use_special_boundaries_velocity)
     call MPAS_pool_get_config(configPool, "config_use_special_boundaries_velocity_masks", config_use_special_boundaries_velocity_masks)
     call MPAS_pool_get_config(configPool, "config_use_special_boundaries_tracers", config_use_special_boundaries_tracers)
     call MPAS_pool_get_package(packagePool, "pkgSpecialBoundariesActive", pkgSpecialBoundariesActive)
     if (config_use_special_boundaries_velocity .or. &
         config_use_special_boundaries_velocity_masks .or. &
         config_use_special_boundaries_tracers) then
        pkgSpecialBoundariesActive = .true.
     endif

   end subroutine setup_packages_other!}}}

   !***********************************************************************
   !
   !  routine seaice_setup_decompositions
   !
   !> \brief   Decomposition setup routine
   !> \author  Doug Jacobsen
   !> \date    04/08/2015
   !> \details
   !>  This routine is intended to create the decomposition list within a
   !>  domain type, and register any decompositons the core wants within it.
   !
   !-----------------------------------------------------------------------
   function seaice_setup_decompositions(decompList) result(ierr)!{{{

      use mpas_derived_types
      use mpas_decomp

      implicit none

      type (mpas_decomp_list), pointer :: decompList

      integer :: ierr
      procedure (mpas_decomp_function), pointer :: decompFunc

      ierr = 0

      nullify(decompList)

      call mpas_decomp_create_decomp_list(decompList)

      decompFunc => mpas_uniform_decomp

      if ( associated(decompList) ) then
         call mpas_log_write(' It is associated (decompList that is)')
      end if

      call mpas_decomp_register_method(decompList, 'uniform', decompFunc, iErr)

      if ( iErr == MPAS_DECOMP_NOERR ) then
         iErr = 0
      end if

   end function seaice_setup_decompositions!}}}


   !***********************************************************************
   !
   !  function seaice_setup_clock
   !
   !> \brief   Pacakge setup routine
   !> \author  Michael Duda
   !> \date    6 August 2014
   !> \details
   !>  The purpose of this function is to allow the core to set up a simulation
   !>  clock that will be used by the I/O subsystem for timing reads and writes
   !>  of I/O streams.
   !>  This function is called from the superstructure after the framework
   !>  has been initialized but before any fields have been allocated and
   !>  initial fields have been read from input files. However, all namelist
   !>  options are available.
   !
   !-----------------------------------------------------------------------
   function seaice_setup_clock(core_clock, configs) result(ierr)!{{{

      use mpas_derived_types

      implicit none

      type (MPAS_Clock_type), intent(inout) :: core_clock
      type (mpas_pool_type), intent(inout) :: configs
      integer :: ierr

      call seaice_simulation_clock_init(core_clock, configs, ierr)

   end function seaice_setup_clock!}}}

   !***********************************************************************
   !
   !  function seaice_setup_log
   !
   !> \brief   Log setup routine
   !> \author  Matt Hoffman, Adrian K. Turner
   !> \date    21 March 2017
   !> \details
   !>  The purpose of this routine is to set up the logging manager
   !>  and allow the core to specify details of the configuration.
   !
   !-----------------------------------------------------------------------
   function seaice_setup_log(logInfo, domain, unitNumbers) result(iErr)!{{{

      use mpas_derived_types
      use mpas_log

      implicit none

      type (mpas_log_type), intent(inout), pointer :: logInfo    !< logging information object to set up
      type (domain_type), intent(in), pointer :: domain          !< domain object to provide info for setting up log manager
      integer, dimension(2), intent(in), optional :: unitNumbers !< Fortran unit numbers to use for output and error logs
      integer :: iErr

      ! Local variables
      integer :: local_err

      iErr = 0

      ! Initialize log manager
      call mpas_log_init(logInfo, domain, unitNumbers=unitNumbers, err=local_err)
      iErr = ior(iErr, local_err)

      ! Set core specific options here
      ! (At present, there are not any.  There could eventually be choices about the file naming conventions
      !  or other settings controlling behavior.)

      ! After core has had a chance to modify log defaults, open the output log
      call mpas_log_open(err=local_err)
      iErr = ior(iErr, local_err)

   end function seaice_setup_log!}}}

   !***********************************************************************
   !
   !  function seaice_get_mesh_stream
   !
   !> \brief   Returns the name of the stream containing mesh information
   !> \author  Michael Duda
   !> \date    8 August 2014
   !> \details
   !>  This function returns the name of the I/O stream containing dimensions,
   !>  attributes, and mesh fields needed by the framework bootstrapping
   !>  routine. At the time this routine is called, only namelist options
   !>  are available.
   !
   !-----------------------------------------------------------------------
   function seaice_get_mesh_stream(configs, stream) result(ierr)!{{{

      use mpas_derived_types
      use mpas_pool_routines

      implicit none

      type (mpas_pool_type), intent(inout) :: configs
      character(len=StrKIND), intent(out) :: stream
      integer :: ierr

      logical, pointer :: config_do_restart
      character(len=strKIND), pointer :: config_initial_condition_type

      ierr = 0

      call mpas_pool_get_config(configs, 'config_initial_condition_type', config_initial_condition_type)
      call mpas_pool_get_config(configs, 'config_do_restart', config_do_restart)

      if (.not. associated(config_do_restart)) then
         write(stream,'(a)') ''
         ierr = 1
      else
         if (config_do_restart) then
            write(stream,'(a)') 'restart'
         else
            if (trim(config_initial_condition_type) == "restart") then
               write(stream,'(a)') 'restart_ic'
            else
               write(stream,'(a)') 'input'
            end if
         endif
      endif

   end function seaice_get_mesh_stream!}}}


   !***********************************************************************
   !
   !  function seaice_setup_block
   !
   !> \brief   Ocean block setup function
   !> \author  Doug Jacobsen
   !> \date    03/18/2015
   !> \details
   !>  This function is a wrapper function to properly setup a block to be a
   !>  ocean core block.
   !
   !-----------------------------------------------------------------------
   function seaice_setup_block(block) result(iErr)!{{{
      use mpas_derived_types
      type (block_type), pointer :: block
      integer :: iErr

      iErr = 0
      call seaice_generate_structs(block, block % structs, block % dimensions, block % packages)
   end function seaice_setup_block!}}}

#include "../inc/setup_immutable_streams.inc"

#include "../inc/block_dimension_routines.inc"

#include "../inc/define_packages.inc"

#include "../inc/structs_and_variables.inc"

#include "../inc/namelist_call.inc"

#include "../inc/namelist_defines.inc"

end module seaice_core_interface
