!$Id: extrapolation.F90 5623 2012-01-17 17:55:26Z connork@uwm.edu $
module extrapolation

  implicit none

  public :: lin_ext_zm_bottom, lin_ext_zt_bottom

  private ! Default scope

  contains
!===============================================================================
  pure function lin_ext_zm_bottom( var_zmp2, var_zmp1,  & 
                                   zmp2, zmp1, zm )  & 
  result( var_zm )

    ! Description:
    !   This function computes the value of a momentum-level variable at a bottom
    !   grid level by using a linear extension of the values of the variable at
    !   the two levels immediately above the level where the result value is
    !   needed.

    ! References:
    !   None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none


    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      var_zmp2,    & ! Momentum level variable at level (k+2)  [units vary]
      var_zmp1,    & ! Momentum level variable at level (k+1)  [units vary]
      zmp2,        & ! Altitude at momentum level (k+2)        [m]
      zmp1,        & ! Altitude at momentum level (k+1)        [m]
      zm             ! Altitude at momentum level (k)          [m]

    ! Return Variable
    real( kind = core_rknd ) :: var_zm   ! Momentum level variable at level (k)    [units vary]

    ! ---- Begin Code -----

    var_zm = ( ( var_zmp2 - var_zmp1 ) / ( zmp2 - zmp1 ) ) & 
             * ( zm - zmp1 ) + var_zmp1

    return
  end function lin_ext_zm_bottom

!===============================================================================
  pure function lin_ext_zt_bottom( var_ztp2, var_ztp1,  & 
                                   ztp2, ztp1, zt )  & 
  result( var_zt )

    ! Description:
    !   This function computes the value of a thermodynamic-level variable at a
    !   bottom grid level by using a linear extension of the values of the
    !   variable at the two levels immediately above the level where the result
    !   value is needed.
    !
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      var_ztp2,    & ! Thermodynamic level variable at level (k+2)  [units vary]
      var_ztp1,    & ! Thermodynamic level variable at level (k+1)  [units vary]
      ztp2,        & ! Altitude at thermodynamic level (k+2)        [m]
      ztp1,        & ! Altitude at thermodynamic level (k+1)        [m]
      zt             ! Altitude at thermodynamic level (k)          [m]

    ! Return Variable
    real( kind = core_rknd ) :: var_zt   ! Thermodynamic level variable at level (k)    [units vary]

    ! ---- Begin Code -----

    var_zt = ( ( var_ztp2 - var_ztp1 ) / ( ztp2 - ztp1 ) ) & 
             * ( zt - ztp1 ) + var_ztp1

    return
  end function lin_ext_zt_bottom

end module extrapolation
