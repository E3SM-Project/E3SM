module cice_unit_test

  use mpas_grid_types

  private
  public :: &
       perform_unit_test

contains

  !--------------------------------------------------------------------------

  subroutine perform_unit_test(block, configs)

    use cice_dynamics_unit_tests, only: &
         strain_rate_operator_unit_test, &
         stress_divergence_operator_unit_test, &
         constitutive_relationship_unit_test
    use cice_therm_zerolayer, only: &
         zero_layer_thermo_semtner_unit_test

    type(block_type), intent(inout) :: block

    type(MPAS_pool_type), pointer, intent(in) :: &
         configs

    character(len=strKIND), pointer :: &
         config_unit_test_type, &
         config_unit_test_subtype

    call MPAS_pool_get_config(configs, "config_unit_test_type", config_unit_test_type)
    call MPAS_pool_get_config(configs, "config_unit_test_subtype", config_unit_test_subtype)

    select case (trim(config_unit_test_type))
    case ("strain rate operator")
       call strain_rate_operator_unit_test(block, trim(config_unit_test_subtype))
    case ("stress divergence operator")
       call stress_divergence_operator_unit_test(block, trim(config_unit_test_subtype))
    case ("constitutive relationship")
       call constitutive_relationship_unit_test(block)
    case ("semtner 1d thermo test")
       call zero_layer_thermo_semtner_unit_test()
    case default
       write(*,*) "Unit test unknown!"
       write(*,*) "config_unit_test_type:", trim(config_unit_test_type)
    end select

    stop

  end subroutine perform_unit_test

  !--------------------------------------------------------------------------

end module cice_unit_test

