! $Id: T_in_K_module.F90 5623 2012-01-17 17:55:26Z connork@uwm.edu $ 

module T_in_K_module

  implicit none

  private ! Default scope

  public :: thlm2T_in_K, T_in_K2thlm

  contains

!-------------------------------------------------------------------------------
  elemental function thlm2T_in_K( thlm, exner, rcm )  & 
    result( T_in_K )

! Description:
!   Calculates absolute temperature from liquid water potential
!   temperature.  (Does not include ice.)

! References: 
!   Cotton and Anthes (1989), "Storm and Cloud Dynamics", Eqn. (2.51). 
!-------------------------------------------------------------------------------
    use constants_clubb, only: & 
      ! Variable(s) 
      Cp,  & ! Dry air specific heat at constant p [J/kg/K]
      Lv     ! Latent heat of vaporization         [J/kg]

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input 
    real( kind = core_rknd ), intent(in) :: & 
      thlm,   & ! Liquid potential temperature  [K]
      exner,  & ! Exner function                [-]
      rcm       ! Liquid water mixing ratio     [kg/kg]

    real( kind = core_rknd ) :: & 
      T_in_K ! Result temperature [K]

    ! ---- Begin Code ----

    T_in_K = thlm * exner + Lv * rcm / Cp

    return
  end function thlm2T_in_K
!-------------------------------------------------------------------------------
  elemental function T_in_K2thlm( T_in_K, exner, rcm )  & 
    result( thlm )

! Description:
!   Calculates liquid water potential temperature from absolute temperature 

! References: 
!   None
!-------------------------------------------------------------------------------
    use constants_clubb, only: & 
      ! Variable(s) 
      Cp,  & ! Dry air specific heat at constant p [J/kg/K]
      Lv     ! Latent heat of vaporization         [J/kg]

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input 
    real( kind = core_rknd ), intent(in) :: & 
      T_in_K, &! Result temperature [K]
      exner,  & ! Exner function                [-]
      rcm       ! Liquid water mixing ratio     [kg/kg]

    real( kind = core_rknd ) :: & 
      thlm    ! Liquid potential temperature  [K]

    ! ---- Begin Code ----

    thlm = ( T_in_K - Lv/Cp * rcm ) / exner 

    return
  end function T_in_K2thlm
!-------------------------------------------------------------------------------

end module T_in_K_module
