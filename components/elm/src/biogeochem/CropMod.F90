module CropMod
  !-----------------------------------------------------------------------
  ! !MODULE: CropMod
  !
  ! !DESCRIPTION:
  ! Module holding routines used in crop model 
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use elm_varcon          , only : tfrz
  use VegetationType      , only : veg_pp
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: calculate_eto
  public :: plant_month
  real(r8) :: dt                            ! radiation time step delta t (seconds)
  real(r8) :: fracday                       ! dtime as a fraction of day

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine calculate_eto(T, rn, g, p, rh, u, es, dt, eto)

    ! !DESCRIPTION
    ! calculates the reference evapotranspiration for precipication seasonality
    ! to determine crop planting date
    ! method from penman-monteith

    ! The basic equation looks like this:
    ! ETo = (0.408*delta*(Rn-G) + gamma*(900/(T+273))*u2*(es-ea))/
    !         (delta + gamma*(1+0.34*u2))

    ! !USES:
    use elm_varcon       , only : tfrz

    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: T        ! temperature (K)
    real(r8), intent(in)  :: rn       ! net radiation (W/m^2)
    real(r8), intent(in)  :: g        ! soil heat flux (W/m^2)
    real(r8), intent(in)  :: p        ! surface atmospheric pressure (pa)
    real(r8), intent(in)  :: rh       ! relative humidity
    real(r8), intent(in)  :: u        ! wind speed (m/s)
    real(r8), intent(in)  :: es       ! saturated vapor pressure (Pa)
    real(r8), intent(in)  :: dt       ! seconds in timestep
    real(r8), intent(inout) :: eto    ! reference evapotranspiration (mm)

    ! !Local Variables:
    real(r8) :: t_c, rn_c, g_c, p_c, delta, ea, gamma1, es_c
    real(r8) :: c1 = 0.408_r8       ! constant in ETo eqn
    real(r8) :: c2                  ! from ETo eqn
    real(r8) :: c3 = 0.34_r8        ! constant in ETo eqn
    real(r8) :: c4 = 900._r8        ! constant in ETo eqn
    real(r8) :: c5 = 237.3_r8       ! constant in ETo eqn
    real(r8) :: conv1 = 1000000._r8 ! convert from W m-2 to MJ m-2
    real(r8) :: conv2 = 1000._r8    ! convert to kPa
    real(r8) :: esc1 = 0.6108_r8    ! constant in es eqn
    real(r8) :: esc2 = 17.27_r8     ! constant in es eqn
    real(r8) :: dc1 = 4098._r8      ! constant in delta eqn
    real(r8) :: gc1 = 0.000665_r8   ! constant in phychrometric eqn

    c2 = c4/(86400._r8/dt)  ! convert to /timestep
    t_c   = T - tfrz      ! convert to degrees C
    rn_c  = rn*dt/conv1   ! convert from W m-2 to MJ m-2
    g_c   = g*dt/conv1    ! convert from W m-2 to MJ m-2
    p_c   = p/conv2       ! convert to kPa

    ! calculate saturation vapor pressure (Tetens, 1930)
    es_c  = esc1*exp((esc2*t_c)/(t_c+c5))

    ! vapor pressure
    ea   = (rh/100)*(es_c) 

    ! calculate the slope of the vapor pressure curve
    delta = (dc1 * (esc1 * exp((esc2 * t_c)/(t_c + c5)))) / &
            (t_c + c5)**2

    ! calculate phychrometric constant
    gamma1 = gc1 * p_c

    ! calculate the potential evapotranspiration
    eto = (c1 * delta * (rn_c - g_c) + gamma1 * (c2/(t_c + c5)) * u * (es_c - ea)) / &
          (delta + gamma1 * (1 + c3 * u))

  end subroutine calculate_eto

  !-----------------------------------------------------------------------
  subroutine plant_month(p, cvt, cvp, temp, p2e, minplantjday, plantmonth)

    ! !DESCRIPTION
    ! calculates the plant month based on the type of seasonality

    ! !USES:
    use pftvarcon        , only : planttemp
    use CropType          , only : tcvp, tcvt, cst
    ! !ARGUMENTS:
    implicit none
    integer,  intent(in)  :: p
    real(r8), intent(in)  :: cvt      ! CV of temperature (K)
    real(r8), intent(in)  :: cvp      ! CV of precipitation (mm)
    real(r8), intent(in)  :: temp(:)     ! average monthly temperture
    real(r8), intent(in)  :: p2e(:)      ! 4-month sum of P:ETo
    integer,  intent(in)  :: minplantjday
    real(r8), intent(inout):: plantmonth(1)

    ! !Local Variables:
    integer :: i
    real(r8) :: r, t_m
    integer,  dimension(12) ::  month = (/10,11,12,1,2,3,4,5,6,7,8,9/)
    integer,  dimension(12) ::  jday  = (/1,32,60,91,121,152,182,213,244,274,305,335/)
    real(r8), dimension(12) ::  indx  = (/1,2,3,4,5,6,7,8,9,10,11,12/)
    !-----------------------------------------------------------------------
    associate(                                    &
         ivt        => veg_pp%itype                  & ! Input:  [integer  (:) ] pft vegetation type
         )

    i = 1
    t_m = indx(1)
    do while ((temp(i) .lt. planttemp(ivt(p)) .or. jday(i) .lt. minplantjday) .and. i .lt. 12)
        i = i+1
        t_m = indx(i) ! the reason this is set after the index count update is to make
                      ! sure the plant month is the first month temp > planttemp
    end do
    if (cvp > tcvp) then
       if (cvt >= tcvt) then    ! Both temperature and precip seasonality
          if (minval(temp) .lt. cst)then ! cold season exists
              plantmonth = t_m
          else                            ! no cold season
              plantmonth = month(maxloc(p2e))
          end if
       else                     ! precipitation seasonality
         plantmonth = month(maxloc(p2e))
       end if
    else if (cvt > tcvt) then   ! temperature seasonality
       plantmonth = t_m
    else                        ! no seasonality at all - tropics
       plantmonth = 1
    end if
  end associate

  end subroutine plant_month

end module CropMod

