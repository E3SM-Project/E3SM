! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_atm_diagnostics_manager

    use mpas_io_units, only : stderrUnit

    private

    public :: mpas_atm_diag_setup, &
              mpas_atm_diag_update, &
              mpas_atm_diag_compute, &
              mpas_atm_diag_reset, &
              mpas_atm_diag_cleanup

    contains


    !-----------------------------------------------------------------------
    !  routine MPAS_atm_diag_setup
    !
    !> \brief Initialize the diagnostics manager and all diagnostics
    !> \author Michael Duda
    !> \date   6 September 2016
    !> \details
    !>  Initialize the diagnostics manager and all diagnostics.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_atm_diag_setup(stream_mgr, configs, structs, clock)

        use mpas_atm_diagnostics_utils, only : mpas_atm_diag_utils_init
        use mpas_derived_types, only : MPAS_streamManager_type, MPAS_pool_type, MPAS_clock_type
        use diagnostic_template, only : diagnostic_template_setup
        use isobaric_diagnostics, only : isobaric_diagnostics_setup
        use convective_diagnostics, only : convective_diagnostics_setup
        use pv_diagnostics, only : pv_diagnostics_setup

        implicit none

        type (MPAS_streamManager_type), target, intent(inout) :: stream_mgr
        type (MPAS_pool_type), pointer :: configs
        type (MPAS_pool_type), pointer :: structs
        type (MPAS_clock_type), pointer :: clock


        !
        ! Prepare the diagnostics utilities module for later use by diagnostics
        !
        call mpas_atm_diag_utils_init(stream_mgr)

        call diagnostic_template_setup(configs, structs, clock)
        call isobaric_diagnostics_setup(structs, clock)
        call convective_diagnostics_setup(structs, clock)
        call pv_diagnostics_setup(structs, clock)

    end subroutine mpas_atm_diag_setup


    !-----------------------------------------------------------------------
    !  routine MPAS_atm_diag_update
    !
    !> \brief Handle diagnostics accumulation at the end of each timestep
    !> \author Michael Duda
    !> \date   6 September 2016
    !> \details
    !>  MPAS_atm_diag_update.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_atm_diag_update()

        use diagnostic_template, only : diagnostic_template_update
        use convective_diagnostics, only : convective_diagnostics_update

        implicit none


        call diagnostic_template_update()
        call convective_diagnostics_update()

    end subroutine mpas_atm_diag_update


    !-----------------------------------------------------------------------
    !  routine MPAS_atm_diag_compute
    !
    !> \brief Compute diagnostics before they are written to output streams
    !> \author Michael Duda
    !> \date   6 September 2016
    !> \details
    !>  MPAS_atm_diag_compute.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_atm_diag_compute()

        use diagnostic_template, only : diagnostic_template_compute
        use isobaric_diagnostics, only : isobaric_diagnostics_compute
        use convective_diagnostics, only : convective_diagnostics_compute
        use pv_diagnostics, only : pv_diagnostics_compute

        implicit none


        call diagnostic_template_compute()
        call isobaric_diagnostics_compute()
        call convective_diagnostics_compute()
        call pv_diagnostics_compute()

    end subroutine mpas_atm_diag_compute


    !-----------------------------------------------------------------------
    !  routine MPAS_atm_diag_reset
    !
    !> \brief Resets a diagnostic after it has been computed and written
    !> \author Michael Duda
    !> \date   6 September 2016
    !> \details
    !>  MPAS_atm_diag_reset.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_atm_diag_reset()

        use diagnostic_template, only : diagnostic_template_reset
        use convective_diagnostics, only : convective_diagnostics_reset

        implicit none


        call diagnostic_template_reset()
        call convective_diagnostics_reset()

    end subroutine mpas_atm_diag_reset


    !-----------------------------------------------------------------------
    !  routine MPAS_atm_diag_cleanup
    !
    !> \brief Finalizes diagnostics manager and all diagnostics
    !> \author Michael Duda
    !> \date   6 September 2016
    !> \details
    !>  Finalizes diagnostics manager and all diagnostics.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_atm_diag_cleanup()

        use mpas_atm_diagnostics_utils, only : mpas_atm_diag_utils_finalize
        use diagnostic_template, only : diagnostic_template_cleanup

        implicit none


        call diagnostic_template_cleanup()

        !
        ! Take care of any needed cleanup in the diagnostics utility module
        !
        call mpas_atm_diag_utils_finalize()

    end subroutine mpas_atm_diag_cleanup

end module mpas_atm_diagnostics_manager
