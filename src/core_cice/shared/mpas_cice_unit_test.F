!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_unit_test
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

module cice_unit_test

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_io_units

  private
  save

  public :: &
       cice_perform_unit_test

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_perform_unit_test
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine cice_perform_unit_test(&
       domain)!{{{

    use cice_velocity_solver_unit_tests, only: &
         cice_strain_rate_operator_unit_test, &
         cice_stress_divergence_operator_unit_test, &
         cice_constitutive_relationship_unit_test

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    character(len=strKIND), pointer :: &
         config_unit_test_type, &
         config_unit_test_subtype

    call MPAS_pool_get_config(domain % configs, "config_unit_test_type", config_unit_test_type)
    call MPAS_pool_get_config(domain % configs, "config_unit_test_subtype", config_unit_test_subtype)

    write(stderrUnit,*) "Unit test: ", trim(config_unit_test_type), " ", trim(config_unit_test_subtype)

    select case (trim(config_unit_test_type))
    case ("strain rate operator")
       call cice_strain_rate_operator_unit_test(domain, trim(config_unit_test_subtype))
    case ("stress divergence operator")
       call cice_stress_divergence_operator_unit_test(domain, trim(config_unit_test_subtype))
    case ("constitutive relationship")
       call cice_constitutive_relationship_unit_test(domain)
    case default
       write(*,*) "Unit test unknown!"
       write(*,*) "config_unit_test_type:", trim(config_unit_test_type)
    end select

  end subroutine cice_perform_unit_test!}}}

!-----------------------------------------------------------------------

end module cice_unit_test
