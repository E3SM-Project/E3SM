module QSatMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Computes saturation mixing ratio and the change in saturation
  !
  ! !PUBLIC TYPES:
  use shr_kind_mod , only: r8 => shr_kind_r8
  implicit none
  save
  private 
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: QSat
  public :: rhoSat
  !-----------------------------------------------------------------------

    ! For water vapor (temperature range 0C-100C)
    real(r8), parameter :: a0 =  6.11213476_r8
    real(r8), parameter :: a1 =  0.444007856_r8
    real(r8), parameter :: a2 =  0.143064234e-01_r8
    real(r8), parameter :: a3 =  0.264461437e-03_r8
    real(r8), parameter :: a4 =  0.305903558e-05_r8
    real(r8), parameter :: a5 =  0.196237241e-07_r8
    real(r8), parameter :: a6 =  0.892344772e-10_r8
    real(r8), parameter :: a7 = -0.373208410e-12_r8
    real(r8), parameter :: a8 =  0.209339997e-15_r8
    ! For derivative:water vapor
    real(r8), parameter :: b0 =  0.444017302_r8
    real(r8), parameter :: b1 =  0.286064092e-01_r8
    real(r8), parameter :: b2 =  0.794683137e-03_r8
    real(r8), parameter :: b3 =  0.121211669e-04_r8
    real(r8), parameter :: b4 =  0.103354611e-06_r8
    real(r8), parameter :: b5 =  0.404125005e-09_r8
    real(r8), parameter :: b6 = -0.788037859e-12_r8
    real(r8), parameter :: b7 = -0.114596802e-13_r8
    real(r8), parameter :: b8 =  0.381294516e-16_r8
    ! For ice (temperature range -75C-0C)
    real(r8), parameter :: c0 =  6.11123516_r8
    real(r8), parameter :: c1 =  0.503109514_r8
    real(r8), parameter :: c2 =  0.188369801e-01_r8
    real(r8), parameter :: c3 =  0.420547422e-03_r8
    real(r8), parameter :: c4 =  0.614396778e-05_r8
    real(r8), parameter :: c5 =  0.602780717e-07_r8
    real(r8), parameter :: c6 =  0.387940929e-09_r8
    real(r8), parameter :: c7 =  0.149436277e-11_r8
    real(r8), parameter :: c8 =  0.262655803e-14_r8
    ! For derivative:ice
    real(r8), parameter :: d0 =  0.503277922_r8
    real(r8), parameter :: d1 =  0.377289173e-01_r8
    real(r8), parameter :: d2 =  0.126801703e-02_r8
    real(r8), parameter :: d3 =  0.249468427e-04_r8
    real(r8), parameter :: d4 =  0.313703411e-06_r8
    real(r8), parameter :: d5 =  0.257180651e-08_r8
    real(r8), parameter :: d6 =  0.133268878e-10_r8
    real(r8), parameter :: d7 =  0.394116744e-13_r8
    real(r8), parameter :: d8 =  0.498070196e-16_r8  
contains



  !-----------------------------------------------------------------------
  subroutine QSat (T, p, es, esdT, qs, qsdT)
    !
    ! !DESCRIPTION:
    ! Computes saturation mixing ratio and the change in saturation
    ! mixing ratio with respect to temperature.
    ! Reference:  Polynomial approximations from:
    !             Piotr J. Flatau, et al.,1992:  Polynomial fits to saturation
    !             vapor pressure.  Journal of Applied Meteorology, 31, 1507-1513.
    !
    ! !USES:
    use shr_kind_mod , only: r8 => shr_kind_r8
    use shr_const_mod, only: SHR_CONST_TKFRZ
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: T        ! temperature (K)
    real(r8), intent(in)  :: p        ! surface atmospheric pressure (pa)
    real(r8), intent(out) :: es       ! vapor pressure (pa)
    real(r8), intent(out) :: esdT     ! d(es)/d(T)
    real(r8), intent(out) :: qs       ! humidity (kg/kg)
    real(r8), intent(out) :: qsdT     ! d(qs)/d(T)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: T_limit
    real(r8) :: td,vp,vp1,vp2
    !-----------------------------------------------------------------------

    T_limit = T - SHR_CONST_TKFRZ
    if (T_limit > 100.0_r8) T_limit=100.0_r8
    if (T_limit < -75.0_r8) T_limit=-75.0_r8

    td       = T_limit
    if (td >= 0.0_r8) then
       es   = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 &
            + td*(a5 + td*(a6 + td*(a7 + td*a8)))))))
            
      esdT = b0 + td*(b1 + td*(b2 + td*(b3 + td*(b4 &
            + td*(b5 + td*(b6 + td*(b7 + td*b8)))))))

    else
       es   = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 &
            + td*(c5 + td*(c6 + td*(c7 + td*c8)))))))

       esdT = d0 + td*(d1 + td*(d2 + td*(d3 + td*(d4 &
            + td*(d5 + td*(d6 + td*(d7 + td*d8)))))))

    endif

    es    = es    * 100._r8            ! pa

    esdT  = esdT  * 100._r8            ! pa/K

    
    vp    = 1.0_r8   / (p - 0.378_r8*es)
    vp1   = 0.622_r8 * vp
    vp2   = vp1   * vp

    qs    = es    * vp1             ! kg/kg
    qsdT  = esdT  * vp2 * p         ! 1 / K


  end subroutine QSat


  
!-------------------------------------------------------------------------------
  subroutine rhoSat(T, rho, rhodT)
    ! compute the saturated vapor pressure density and its derivative against the temperature
    ! jyt
    use elm_varcon,    only: rwat
    use shr_const_mod, only: SHR_CONST_TKFRZ

    implicit none
    real(r8)           , intent(in)  :: T
    real(r8)           , intent(out) :: rho
    real(r8), optional , intent(out) :: rhodT


    !------------------

    real(r8) :: T_limit
    real(r8) :: td, es, esdT

    T_limit = T - SHR_CONST_TKFRZ
    if (T_limit > 100.0_r8) T_limit=100.0_r8
    if (T_limit < -75.0_r8) T_limit=-75.0_r8

    td       = T_limit
    if (td >= 0.0_r8) then
       es   = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 &
               + td*(a5 + td*(a6 + td*(a7 + td*a8)))))))
       esdT = b0 + td*(b1 + td*(b2 + td*(b3 + td*(b4 &
               + td*(b5 + td*(b6 + td*(b7 + td*b8)))))))
    else
       es   = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 &
               + td*(c5 + td*(c6 + td*(c7 + td*c8)))))))
       esdT = d0 + td*(d1 + td*(d2 + td*(d3 + td*(d4 &
               + td*(d5 + td*(d6 + td*(d7 + td*d8)))))))
    endif

    es    = es    * 100._r8            ! pa
    rho   = es/(rwat*T)                !kg  m^-3

    if(present(rhodT)) rhodT= esdT/(rwat*T)-rho/T         !kg  m^-3 K^-1

  end subroutine rhoSat
end module QSatMod
