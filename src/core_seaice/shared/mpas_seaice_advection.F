!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_advection
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 19th March 2015
!> \details
!>
!
!-----------------------------------------------------------------------

module seaice_advection

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timer

  implicit none

  private
  save

  public :: &
       seaice_init_advection, &
       seaice_run_advection

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_init_advection
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 19th March 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_init_advection(domain)

    use seaice_advection_upwind, only: &
         seaice_init_advection_upwind

    use seaice_advection_incremental_remap, only: &
         seaice_init_advection_incremental_remap

    type (domain_type), intent(inout) :: &
         domain !< Input/Output:

    logical, pointer :: &
         config_use_advection

    character(len=strKIND), pointer :: &
         config_advection_type

    call MPAS_pool_get_config(domain % blocklist % configs, "config_use_advection", config_use_advection)
    call MPAS_pool_get_config(domain % blocklist % configs, "config_advection_type", config_advection_type)

    if (config_use_advection) then

       if (trim(config_advection_type) == "upwind") then

          call seaice_init_advection_upwind(domain)

       else if (trim(config_advection_type) == "incremental_remap") then

          call seaice_init_advection_incremental_remap(domain)

       endif

    endif

  end subroutine seaice_init_advection

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_run_advection
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 19th March 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_run_advection(domain, clock, ierr)

    use seaice_advection_upwind, only: &
         seaice_run_advection_upwind

    use seaice_advection_incremental_remap, only: &
         seaice_run_advection_incremental_remap

    type (domain_type), intent(inout) :: &
         domain !< Input/Output:

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    integer, intent(inout) :: &
         ierr  !< Input/Output:

    logical, pointer :: &
         config_use_advection

    character(len=strKIND), pointer :: &
         config_advection_type

    call MPAS_pool_get_config(domain % blocklist % configs, "config_use_advection", config_use_advection)
    call MPAS_pool_get_config(domain % blocklist % configs, "config_advection_type", config_advection_type)

    if (config_use_advection) then

       if (trim(config_advection_type) == "upwind") then

          call mpas_timer_start("advection upwind")
          call seaice_run_advection_upwind(domain, clock)
          call mpas_timer_stop("advection upwind")

       else if (trim(config_advection_type) == "incremental_remap") then

          !TODO WHL - Convert ice/snow volume to thickness here
          !           (Later, change volume to thickness throughout code?)

          call mpas_timer_start("advection incr remap")
          call seaice_run_advection_incremental_remap(domain, clock, ierr)
          call mpas_timer_stop("advection incr remap")

          !TODO WHL - Convert ice/snow thickness to volume here

       endif

    endif

  end subroutine seaice_run_advection

!-----------------------------------------------------------------------

end module seaice_advection
