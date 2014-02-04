module cice_unit_test

  use mpas_grid_types

  private
  public :: &
       perform_unit_test

contains

  !--------------------------------------------------------------------------

  subroutine perform_unit_test(block, unit_test_type, unit_test_subtype)

    use cice_dynamics_unit_tests, only: &
         strain_rate_operator_unit_test, &
         stress_divergence_operator_unit_test

    type(block_type), intent(inout) :: block

    character(len=*), intent(in) :: &
         unit_test_type, &
         unit_test_subtype

    select case (trim(unit_test_type))
    case ("strain rate operator")
       call strain_rate_operator_unit_test(block, trim(unit_test_subtype))
    case ("stress divergence operator")
       call stress_divergence_operator_unit_test(block, trim(unit_test_subtype))
    case default
       write(*,*) "Unit test unknown!"
       write(*,*) "config_unit_test_type:", trim(unit_test_type)
    end select

    stop

  end subroutine perform_unit_test

  !--------------------------------------------------------------------------

end module cice_unit_test

