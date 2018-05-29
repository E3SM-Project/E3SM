!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_unit_test
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

module seaice_unit_test

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_log, only: mpas_log_write

  private
  save

  public :: &
       seaice_perform_unit_test

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_perform_unit_test
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_perform_unit_test(&
       domain)!{{{

    use seaice_velocity_solver_unit_tests, only: &
         seaice_strain_rate_operator_unit_test, &
         seaice_stress_divergence_operator_unit_test, &
         seaice_constitutive_relationship_unit_test

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    character(len=strKIND), pointer :: &
         config_unit_test_type, &
         config_unit_test_subtype

    call MPAS_pool_get_config(domain % configs, "config_unit_test_type", config_unit_test_type)
    call MPAS_pool_get_config(domain % configs, "config_unit_test_subtype", config_unit_test_subtype)

    call mpas_log_write("Unit test: "//trim(config_unit_test_type)//" "//trim(config_unit_test_subtype))

    select case (trim(config_unit_test_type))
    case ("strain rate operator")
       call seaice_strain_rate_operator_unit_test(domain, trim(config_unit_test_subtype))
    case ("stress divergence operator")
       call seaice_stress_divergence_operator_unit_test(domain, trim(config_unit_test_subtype))
    case ("constitutive relationship")
       call seaice_constitutive_relationship_unit_test(domain)
    case default
       call mpas_log_write("seaice_perform_unit_test: config_unit_test_type unknown: "//trim(config_unit_test_type))
    end select

  end subroutine seaice_perform_unit_test!}}}

!-----------------------------------------------------------------------

end module seaice_unit_test
