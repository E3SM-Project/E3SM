! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module diagnostic_template

    use mpas_derived_types, only : MPAS_pool_type, MPAS_clock_type

    type (MPAS_pool_type), pointer :: mesh
    type (MPAS_clock_type), pointer :: clock

    public :: diagnostic_template_setup, &
              diagnostic_template_update, &
              diagnostic_template_compute, &
              diagnostic_template_reset, &
              diagnostic_template_cleanup

    private


    contains


    !-----------------------------------------------------------------------
    !  routine diagnostic_template_setup
    !
    !> \brief Initialize the diagnostic
    !> \author 
    !> \date   
    !> \details
    !>  Initialize the diagnostic
    !
    !-----------------------------------------------------------------------
    subroutine diagnostic_template_setup(configs, all_pools, simulation_clock)

        use mpas_derived_types, only : MPAS_pool_type, MPAS_clock_type
        use mpas_pool_routines, only : mpas_pool_get_subpool

        implicit none

        type (MPAS_pool_type), pointer :: configs
        type (MPAS_pool_type), pointer :: all_pools
        type (MPAS_clock_type), pointer :: simulation_clock

        ! Perform initialization, memory allocation, etc.

        ! Also, save pointers to any pools that will be used by this diagnostic
        ! E.g.,
        call mpas_pool_get_subpool(all_pools, 'mesh', mesh)

        clock => simulation_clock
   
    end subroutine diagnostic_template_setup


    !-----------------------------------------------------------------------
    !  routine diagnostic_template_update
    !
    !> \brief Handle diagnostic calculation needed after each timestep
    !> \author 
    !> \date   
    !> \details
    !>  Handle diagnostic calculation needed after each timestep
    !
    !-----------------------------------------------------------------------
    subroutine diagnostic_template_update()

        implicit none

        ! Called at the end of every timestep
        ! Update extrema, accumulations, etc.
   
    end subroutine diagnostic_template_update


    !-----------------------------------------------------------------------
    !  routine diagnostic_template_compute
    !
    !> \brief Compute diagnostic before model output is written
    !> \author 
    !> \date   
    !> \details
    !>  Compute diagnostic before model output is written
    !
    !-----------------------------------------------------------------------
    subroutine diagnostic_template_compute()

        implicit none

        ! Called immediately before diagnostics will be written
        ! Compute the diagnostic
   
    end subroutine diagnostic_template_compute


    !-----------------------------------------------------------------------
    !  routine diagnostic_template_reset
    !
    !> \brief Reset diagnostic after it has been written
    !> \author 
    !> \date   
    !> \details
    !>  Reset diagnostic after it has been written
    !
    !-----------------------------------------------------------------------
    subroutine diagnostic_template_reset()

        implicit none

        ! Called immediately after diagnostics have been written
        ! Reset counters, accumulations, etc.
   
    end subroutine diagnostic_template_reset


    !-----------------------------------------------------------------------
    !  routine diagnostic_template_cleanup
    !
    !> \brief Finalizes diagnostic
    !> \author Michael Duda
    !> \date   6 September 2016
    !> \details
    !>  Finalizes diagnostic
    !
    !-----------------------------------------------------------------------
    subroutine diagnostic_template_cleanup()

        implicit none

        ! Deallocate scratch arrays, etc.
   
    end subroutine diagnostic_template_cleanup

end module diagnostic_template
