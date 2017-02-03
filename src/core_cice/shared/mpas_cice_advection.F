!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_advection
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 19th March 2015
!> \details
!>
!
!-----------------------------------------------------------------------

module cice_advection

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timer

  implicit none

  private
  save

  public :: &
       cice_init_advection, &
       cice_run_advection

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_init_advection
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 19th March 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_init_advection(domain)

    use cice_advection_upwind, only: &
         cice_init_advection_upwind

    use cice_advection_incremental_remap, only: &
         cice_init_advection_incremental_remap

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

          call cice_init_advection_upwind(domain)

       else if (trim(config_advection_type) == "incremental_remap") then

          call cice_init_advection_incremental_remap(domain)

       endif

    endif

  end subroutine cice_init_advection

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_run_advection
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 19th March 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_run_advection(domain, clock)

    use cice_advection_upwind, only: &
         cice_run_advection_upwind

    use cice_advection_incremental_remap, only: &
         cice_run_advection_incremental_remap

    type (domain_type), intent(inout) :: &
         domain !< Input/Output:

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    logical, pointer :: &
         config_use_advection

    character(len=strKIND), pointer :: &
         config_advection_type

    call MPAS_pool_get_config(domain % blocklist % configs, "config_use_advection", config_use_advection)
    call MPAS_pool_get_config(domain % blocklist % configs, "config_advection_type", config_advection_type)

    if (config_use_advection) then

       if (trim(config_advection_type) == "upwind") then

          call mpas_timer_start("advection upwind")
          call cice_run_advection_upwind(domain, clock)
          call mpas_timer_stop("advection upwind")

       else if (trim(config_advection_type) == "incremental_remap") then

          !TODO WHL - Convert ice/snow volume to thickness here
          !           (Later, change volume to thickness throughout code?)

          call mpas_timer_start("advection incr remap")
          call cice_run_advection_incremental_remap(domain, clock)
          call mpas_timer_stop("advection incr remap")

          !TODO WHL - Convert ice/snow thickness to volume here

       endif

    endif

  end subroutine cice_run_advection

!-----------------------------------------------------------------------

end module cice_advection
