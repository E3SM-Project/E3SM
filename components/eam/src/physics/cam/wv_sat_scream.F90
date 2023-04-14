! Include bit-for-bit math macros.
#include "bfb_math.inc"

module wv_sat_scream

  !------------------------------------------------------------------------------------
  !This module has functions which scream uses for computing saturation vapor pressure
  !of ice and liquid.
  !Note: I just copied them over from micro_p3.F90 without any cleanup as this code
  ! not be maintained on long term basis
  !------------------------------------------------------------------------------------

  ! get real kind from utils
  use physics_utils,  only: rtype
  use micro_p3_utils, only: T_zerodegc
#ifdef SCREAM_CONFIG_IS_CMAKE
  use physics_share_f2c, only: cxx_pow, cxx_sqrt, cxx_cbrt, cxx_gamma, cxx_log, &
                                 cxx_log10, cxx_exp, cxx_tanh
#endif


  implicit none
  private

  public:: qv_sat, MurphyKoop_svp

contains

  !===========================================================================================
  real(rtype) function qv_sat(t_atm,p_atm,i_wrt)

    !------------------------------------------------------------------------------------
    ! Calls polysvp1 to obtain the saturation vapor pressure, and then computes
    ! and returns the saturation mixing ratio, with respect to either liquid or ice,
    ! depending on value of 'i_wrt'
    !------------------------------------------------------------------------------------

    use micro_p3_utils, only: ep_2
    implicit none

    !Calling parameters:
    real(rtype), intent(in)    :: t_atm  !temperature [K]
    real(rtype), intent(in)    :: p_atm  !pressure    [Pa]
    integer, intent(in) :: i_wrt  !index, 0 = w.r.t. liquid, 1 = w.r.t. ice

    !Local variables:
    real(rtype)            :: e_pres         !saturation vapor pressure [Pa]

    !e_pres = polysvp1(t_atm,i_wrt)
    e_pres = MurphyKoop_svp(t_atm,i_wrt)
    qv_sat = ep_2*e_pres/max(1.e-3_rtype,(p_atm-e_pres))

    return

  end function qv_sat
  !===========================================================================================

  !==========================================================================================!

  real(rtype) function MurphyKoop_svp(t, i_type)

    use scream_abortutils, only : endscreamrun

    implicit none

    !-------------------------------------------------------------------
    ! Compute saturation vapor pressure (returned in units of pa)
    ! Inputs:
    ! "t", units [K]
    !  i_type refers to saturation with respect to liquid (0) or ice (1)
    ! Author: Balwinder Singh (algorithm from CAM routine)
    !--------------------------------------------------------------------

    !Murphy & Koop (2005)
    real(rtype), intent(in) :: t
    integer, intent(in)     :: i_type

    !local vars
    character(len=1000) :: err_msg
    real(rtype)         :: logt, tmp

    !parameters for ice saturation eqn
    real(rtype), parameter :: ic(4)  =(/9.550426_rtype, 5723.265_rtype, 3.53068_rtype, &
         0.00728332_rtype/)

    !parameters for liq saturation eqn
    real(rtype), parameter :: lq(10) = (/54.842763_rtype, 6763.22_rtype, 4.210_rtype, &
         0.000367_rtype, 0.0415_rtype, 218.8_rtype, 53.878_rtype, 1331.22_rtype,       &
         9.44523_rtype, 0.014025_rtype /)

    !------------------
    !Check if temprature is within legitimate range
    call check_temp(t, "MurphyKoop_svp")

    logt = bfb_log(t)

    if (i_type .eq. 1 .and. t .lt. T_zerodegc) then

       !(good down to 110 K)
       MurphyKoop_svp = bfb_exp(ic(1) - (ic(2) / t) + (ic(3) * logt) - (ic(4) * t))

    elseif (i_type .eq. 0 .or. t .ge. T_zerodegc) then

       ! (good for 123 < T < 332 K)
       !For some reason, we cannot add line breaks if we use "bfb_exp", storing experssion in "tmp"
       tmp = lq(1) - (lq(2) / t) - (lq(3) * logt) + (lq(4) * t) + &
            (bfb_tanh(lq(5) * (t - lq(6))) * (lq(7) - (lq(8) / t) - &
            (lq(9) * logt) + lq(10) * t))
       MurphyKoop_svp = bfb_exp(tmp)
    else

       write(err_msg,*)'Error: Either MurphyKoop_svp i_type is not 0 or 1 or t=NaN. itype= ', &
            i_type,' and temperature t=',t,' in file: ',__FILE__,' at line:',__LINE__
       call endscreamrun(err_msg)
    endif

    return
  end function MurphyKoop_svp

  !_rtype
  real(rtype) function polysvp1(t,i_type)

    !-------------------------------------------
    !  COMPUTE SATURATION VAPOR PRESSURE
    !  POLYSVP1 RETURNED IN UNITS OF PA.
    !  T IS INPUT IN UNITS OF K.
    !  i_type REFERS TO SATURATION WITH RESPECT TO LIQUID (0) OR ICE (1)
    !-------------------------------------------

    use scream_abortutils, only : endscreamrun

    use debug_info, only: report_error_info
    implicit none

    real(rtype), intent(in) :: t
    integer, intent(in)     :: i_type

    ! REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM FLATAU ET AL. 1992, TABLE 4 (RIGHT-HAND COLUMN)

    !local variables
    character(len=1000) :: err_msg

    ! ice
    real(rtype) a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i
    data a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i /&
         6.11147274_rtype,     0.503160820_rtype,     0.188439774e-1_rtype, &
         0.420895665e-3_rtype, 0.615021634e-5_rtype,  0.602588177e-7_rtype, &
         0.385852041e-9_rtype, 0.146898966e-11_rtype, 0.252751365e-14_rtype/

    ! liquid
    real(rtype) a0,a1,a2,a3,a4,a5,a6,a7,a8

    ! V1.7
    data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
         6.11239921_rtype,      0.443987641_rtype,     0.142986287e-1_rtype, &
         0.264847430e-3_rtype,  0.302950461e-5_rtype,  0.206739458e-7_rtype, &
         0.640689451e-10_rtype,-0.952447341e-13_rtype,-0.976195544e-15_rtype/
    real(rtype) dt

    !-------------------------------------------

    !------------------
    !Check if temprature is within legitimate range
    call check_temp(t, "polysvp1")


    if (i_type.eq.1 .and. t.lt.T_zerodegc) then
       ! ICE

       !       Flatau formulation:
       dt       = max(-80._rtype,t-273.15_rtype)
       polysvp1 = a0i + dt*(a1i+dt*(a2i+dt*(a3i+dt*(a4i+dt*(a5i+dt*(a6i+dt*(a7i+       &
            a8i*dt)))))))
       polysvp1 = polysvp1*100._rtype

       !       Goff-Gratch formulation:
       !        POLYSVP1 = 10.**(-9.09718*(273.16/T-1.)-3.56654*                 &
       !          log10(273.16/T)+0.876793*(1.-T/273.16)+                        &
       !          log10(6.1071))*100.


    elseif (i_type.eq.0 .or. t.ge.T_zerodegc) then
       ! LIQUID

       !       Flatau formulation:
       dt       = max(-80._rtype,t-273.15_rtype)
       polysvp1 = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
       polysvp1 = polysvp1*100._rtype

       !       Goff-Gratch formulation:
       !        POLYSVP1 = 10.**(-7.90298*(373.16/T-1.)+                         &
       !             5.02808*log10(373.16/T)-                                    &
       !             1.3816E-7*(10**(11.344*(1.-T/373.16))-1.)+                  &
       !             8.1328E-3*(10**(-3.49149*(373.16/T-1.))-1.)+                &
       !             log10(1013.246))*100.

       !PMC added error checking
    else

       call report_error_info('Something went wrong', 'polysvp1')
       write(err_msg,*)'** polysvp1 i_type must be 0 or 1 but is: ', &
            i_type,' temperature is:',t,' in file: ',__FILE__, &
            ' at line:',__LINE__
       call endscreamrun(err_msg)
    endif

    return

  end function polysvp1

  subroutine check_temp(t, subname)
    !Check if temprature values are in legit range
    use scream_abortutils, only : endscreamrun
    use ieee_arithmetic,   only : ieee_is_finite, ieee_is_nan

    implicit none

    real(rtype),      intent(in) :: t
    character(len=*), intent(in) :: subname

    character(len=1000) :: err_msg

    if(t <= 0.0_rtype) then
       write(err_msg,*)'Error: Called from:',trim(adjustl(subname)),'; Temperature is:',t,' which is <= 0._r8 in file:',__FILE__, &
            ' at line:',__LINE__
       call endscreamrun(err_msg)
    elseif(.not. ieee_is_finite(t)) then
       write(err_msg,*)'Error: Called from:',trim(adjustl(subname)),'; Temperature is:',t,' which is not finite in file:', &
            __FILE__,' at line:',__LINE__
       call endscreamrun(err_msg)
    elseif(ieee_is_nan(t)) then
       write(err_msg,*)'Error: Called from:',trim(adjustl(subname)),'; Temperature is:',t,' which is NaN in file:',__FILE__, &
            'at line:',__LINE__
       call endscreamrun(err_msg)
    endif

    return
  end subroutine check_temp

  !------------------------------------------------------------------------------------------!

end module wv_sat_scream
