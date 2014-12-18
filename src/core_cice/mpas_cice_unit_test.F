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

  use mpas_grid_types

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
       block, &
       configs)!{{{

    use cice_velocity_solver_unit_tests, only: &
         cice_strain_rate_operator_unit_test, &
         cice_stress_divergence_operator_unit_test, &
         cice_constitutive_relationship_unit_test

    use cice_therm_zerolayer, only: &
         cice_zero_layer_thermo_semtner_unit_test

    type(block_type), intent(inout) :: &
         block !< Input/Output: 

    type(MPAS_pool_type), pointer, intent(in) :: &
         configs !< Input: 

    character(len=strKIND), pointer :: &
         config_unit_test_type, &
         config_unit_test_subtype

    call MPAS_pool_get_config(configs, "config_unit_test_type", config_unit_test_type)
    call MPAS_pool_get_config(configs, "config_unit_test_subtype", config_unit_test_subtype)

    write(stderrUnit,*) "Unit test: ", trim(config_unit_test_type), " ", trim(config_unit_test_subtype)

    select case (trim(config_unit_test_type))
    case ("strain rate operator")
       call cice_strain_rate_operator_unit_test(block, trim(config_unit_test_subtype))
    case ("stress divergence operator")
       call cice_stress_divergence_operator_unit_test(block, trim(config_unit_test_subtype))
    case ("constitutive relationship")
       call cice_constitutive_relationship_unit_test(block)
    case ("semtner 1d thermo test")
       call cice_zero_layer_thermo_semtner_unit_test()
    case default
       write(*,*) "Unit test unknown!"
       write(*,*) "config_unit_test_type:", trim(config_unit_test_type)
    end select

  end subroutine cice_perform_unit_test!}}}

!-----------------------------------------------------------------------

end module cice_unit_test

