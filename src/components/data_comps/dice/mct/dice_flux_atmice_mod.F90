module dice_flux_atmice_mod

  use shr_kind_mod    , only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_const_mod   ! shared constants
  use shr_sys_mod     ! shared system routines

  implicit none

  real(R8) :: loc_zvir   = shr_const_zvir
  real(R8) :: loc_cpdair = shr_const_cpdair
  real(R8) :: loc_cpvir  = shr_const_cpvir
  real(R8) :: loc_karman = shr_const_karman
  real(R8) :: loc_g      = shr_const_g
  real(R8) :: loc_latvap = shr_const_latvap
  real(R8) :: loc_latice = shr_const_latice
  real(R8) :: loc_stebol = shr_const_stebol

  integer,parameter :: dbug = 0 ! internal debug level

!===============================================================================
contains
!===============================================================================

  subroutine dice_flux_atmice(             &
       mask  ,zbot  ,ubot  ,vbot  ,thbot,  &
       qbot  ,rbot  ,tbot  ,ts    ,sen,    &
       lat   ,lwup  ,evap  ,taux  ,tauy,   &
       tref  ,qref  ,logunit               )

    !-------------------------------------------------------------------------------
    ! PURPOSE:
    !   using atm & ice state variables, compute atm/ice fluxes
    !   and diagnostic 10m air temperature and humidity
    !
    ! NOTE:
    !   o all fluxes are positive downward
    !   o net heat flux = net sw + lw up + lw down + sen + lat
    !   o here, tstar = <WT>/U*, and qstar = <WQ>/U*.
    !   o wind speeds should all be above a minimum speed (eg. 1.0 m/s)
    !
    ! ASSUME:
    !   o The saturation humidity of air at T(K): qsat(T)  (kg/m^3)
    !-------------------------------------------------------------------------------

    !--- input arguments --------------------------------
    integer(IN),intent(in)  :: mask (:)    ! 0 <=> cell NOT in model domain
    real(R8)   ,intent(in)  :: zbot (:)    ! atm level height  (m)
    real(R8)   ,intent(in)  :: ubot (:)    ! atm u wind     (m/s)
    real(R8)   ,intent(in)  :: vbot (:)    ! atm v wind     (m/s)
    real(R8)   ,intent(in)  :: thbot(:)    ! atm potential T   (K)
    real(R8)   ,intent(in)  :: qbot (:)    ! atm specific humidity (kg/kg)
    real(R8)   ,intent(in)  :: rbot (:)    ! atm air density   (kg/m^3)
    real(R8)   ,intent(in)  :: tbot (:)    ! atm T       (K)
    real(R8)   ,intent(in)  :: ts   (:)    ! surface temperature
    integer(IN),intent(in)  :: logunit     ! logging unit number

    !--- output arguments -------------------------------
    real(R8)   ,intent(out) :: sen  (:)    ! sensible      heat flux  (W/m^2)
    real(R8)   ,intent(out) :: lat  (:)    ! latent        heat flux  (W/m^2)
    real(R8)   ,intent(out) :: lwup (:)    ! long-wave upward heat flux  (W/m^2)
    real(R8)   ,intent(out) :: evap (:)    ! evaporative water flux ((kg/s)/m^2)
    real(R8)   ,intent(out) :: taux (:)    ! x surface stress (N)
    real(R8)   ,intent(out) :: tauy (:)    ! y surface stress (N)
    real(R8)   ,intent(out) :: tref (:)    ! 2m reference height temperature
    real(R8)   ,intent(out) :: qref (:)    ! 2m reference height humidity

    !--- local constants --------------------------------
    real(R8),parameter :: umin   =  1.0_R8            ! minimum wind speed (m/s)
    real(R8),parameter :: zref   = 10.0_R8            ! ref height           ~ m
    real(R8),parameter :: ztref  =  2.0_R8            ! ref height for air T ~ m
    real(R8),parameter :: spval  = shr_const_spval    ! special value
    real(R8),parameter :: zzsice = 0.0005_R8          ! ice surface roughness

    !--- local variables --------------------------------
    integer(IN) :: lsize  ! array dimensions
    integer(IN) :: n      ! array indicies
    real(R8)    :: vmag   ! surface wind magnitude   (m/s)
    real(R8)    :: thvbot ! virtual temperature      (K)
    real(R8)    :: ssq    ! sea surface humidity     (kg/kg)
    real(R8)    :: dssqdt ! derivative of ssq wrt Ts (kg/kg/K)
    real(R8)    :: delt   ! potential T difference   (K)
    real(R8)    :: delq   ! humidity difference      (kg/kg)
    real(R8)    :: stable ! stability factor
    real(R8)    :: rdn    ! sqrt of neutral exchange coefficient (momentum)
    real(R8)    :: rhn    ! sqrt of neutral exchange coefficient (heat)
    real(R8)    :: ren    ! sqrt of neutral exchange coefficient (water)
    real(R8)    :: rd     ! sqrt of exchange coefficient (momentum)
    real(R8)    :: rh     ! sqrt of exchange coefficient (heat)
    real(R8)    :: re     ! sqrt of exchange coefficient (water)
    real(R8)    :: ustar  ! ustar
    real(R8)    :: qstar  ! qstar
    real(R8)    :: tstar  ! tstar
    real(R8)    :: hol    ! H (at zbot) over L
    real(R8)    :: xsq    ! temporary variable
    real(R8)    :: xqq    ! temporary variable
    real(R8)    :: psimh  ! stability function at zbot (momentum)
    real(R8)    :: psixh  ! stability function at zbot (heat and water)
    real(R8)    :: alz    ! ln(zbot/z10)
    real(R8)    :: ltheat ! latent heat for surface
    real(R8)    :: tau    ! stress at zbot
    real(R8)    :: cp     ! specific heat of moist air

    real(R8)    :: bn     ! exchange coef funct for interpolation
    real(R8)    :: bh     ! exchange coef funct for interpolation
    real(R8)    :: fac    ! interpolation factor
    real(R8)    :: ln0    ! log factor for interpolation
    real(R8)    :: ln3    ! log factor for interpolation

    !--- local functions --------------------------------
    real(R8)   :: Tk      ! temperature (K)
    real(R8)   :: qsat    ! the saturation humidity of air (kg/m^3)
    real(R8)   :: dqsatdt ! derivative of qsat wrt surface temperature
    real(R8)   :: xd      ! dummy argument
    real(R8)   :: psimhu  ! unstable part of psimh
    real(R8)   :: psixhu  ! unstable part of psimx

    qsat(Tk)    = 627572.4_R8 / exp(5107.4_R8/Tk)
    dqsatdt(Tk) = (5107.4_R8 / Tk**2) * 627572.4_R8 / exp(5107.4_R8/Tk)
    psimhu(xd)  = log((1.0_R8+xd*(2.0_R8+xd))*(1.0_R8+xd*xd)/8.0_R8) - 2.0_R8*atan(xd) + 1.571_R8
    psixhu(xd)  =  2.0_R8 * log((1.0_R8 + xd*xd)/2.0_R8)

    !--- formats ----------------------------------------
    character(*),parameter ::    F01 = "('(dice_flux_atmIce) ',a, i7,2x,d21.14)"
    character(*),parameter :: subName =  "(dice_flux_atmIce) "
    !-------------------------------------------------------------------------------

    lsize = size(tbot)

    do n = 1,lsize

       if (mask(n) == 0) then
          sen  (n) = spval
          lat  (n) = spval
          lwup (n) = spval
          evap (n) = spval
          taux (n) = spval
          tauy (n) = spval
          tref (n) = spval
          qref (n) = spval
       else
          !--- define some needed variables ---
          vmag   = max(umin, sqrt(ubot(n)**2+vbot(n)**2))
          thvbot = thbot(n)*(1.0_R8 + loc_zvir * qbot(n)) ! virtual pot temp (K)
          ssq   =  qsat  (ts(n)) / rbot(n)           ! sea surf hum (kg/kg)
          dssqdt = dqsatdt(ts(n)) / rbot(n)           ! deriv of ssq wrt Ts
          delt   = thbot(n) - ts(n)                   ! pot temp diff (K)
          delq   = qbot(n) - ssq                        ! spec hum dif (kg/kg)
          alz    = log(zbot(n)/zref)
          cp     = loc_cpdair*(1.0_R8 + loc_cpvir*ssq)
          ltheat = loc_latvap + loc_latice

          !----------------------------------------------------------
          ! first estimate of Z/L and ustar, tstar and qstar
          !----------------------------------------------------------

          !--- neutral coefficients, z/L = 0.0 ---
          rdn = loc_karman/log(zref/zzsice)
          rhn = rdn
          ren = rdn

          !--- ustar,tstar,qstar ----
          ustar = rdn * vmag
          tstar = rhn * delt
          qstar = ren * delq

          !--- compute stability & evaluate all stability functions ---
          hol    = loc_karman * loc_g * zbot(n) &
               &     * (tstar/thvbot+qstar/(1.0_R8/loc_zvir+qbot(n))) / ustar**2
          hol    = sign( min(abs(hol),10.0_R8), hol )
          stable = 0.5_R8 + sign(0.5_R8 , hol)
          xsq    = max(sqrt(abs(1.0_R8 - 16.0_R8*hol)) , 1.0_R8)
          xqq    = sqrt(xsq)
          psimh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psimhu(xqq)
          psixh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)

          !--- shift all coeffs to measurement height and stability ---
          rd = rdn / (1.0_R8+rdn/loc_karman*(alz-psimh))
          rh = rhn / (1.0_R8+rhn/loc_karman*(alz-psixh))
          re = ren / (1.0_R8+ren/loc_karman*(alz-psixh))

          !--- update ustar, tstar, qstar w/ updated, shifted coeffs --
          ustar = rd * vmag
          tstar = rh * delt
          qstar = re * delq

          !----------------------------------------------------------
          ! iterate to converge on Z/L, ustar, tstar and qstar
          !----------------------------------------------------------

          !--- compute stability & evaluate all stability functions ---
          hol    = loc_karman * loc_g * zbot(n) &
               &      * (tstar/thvbot+qstar/(1.0_R8/loc_zvir+qbot(n))) / ustar**2
          hol    = sign( min(abs(hol),10.0_R8), hol )
          stable = 0.5_R8 + sign(0.5_R8 , hol)
          xsq    = max(sqrt(abs(1.0_R8 - 16.0_R8*hol)) , 1.0_R8)
          xqq    = sqrt(xsq)
          psimh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psimhu(xqq)
          psixh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)

          !--- shift all coeffs to measurement height and stability ---
          rd = rdn / (1.0_R8+rdn/loc_karman*(alz-psimh))
          rh = rhn / (1.0_R8+rhn/loc_karman*(alz-psixh))
          re = ren / (1.0_R8+ren/loc_karman*(alz-psixh))

          !--- update ustar, tstar, qstar w/ updated, shifted coeffs --
          ustar = rd * vmag
          tstar = rh * delt
          qstar = re * delq

          !----------------------------------------------------------
          ! compute the fluxes
          !----------------------------------------------------------

          tau = rbot(n) * ustar * ustar

          !--- momentum flux ---
          taux(n) = tau * ubot(n) / vmag
          tauy(n) = tau * vbot(n) / vmag

          !--- heat flux ---
          sen (n) =   cp * tau * tstar / ustar
          lat (n) =  ltheat * tau * qstar / ustar
          lwup(n) = -loc_stebol * ts(n)**4

          !--- water flux ---
          evap(n) = lat(n)/ltheat

          !----------------------------------------------------------
          ! compute diagnostic: 2m reference height temperature
          !----------------------------------------------------------

          !--- Compute function of exchange coefficients. Assume that
          !--- cn = rdn*rdn, cm=rd*rd and ch=rh*rd, and therefore
          !--- 1/sqrt(cn(n))=1/rdn and sqrt(cm(n))/ch(n)=1/rh
          bn = loc_karman/rdn
          bh = loc_karman/rh

          !--- Interpolation factor for stable and unstable cases
          ln0 = log(1.0_R8 + (ztref/zbot(n))*(exp(bn) - 1.0_R8))
          ln3 = log(1.0_R8 + (ztref/zbot(n))*(exp(bn - bh) - 1.0_R8))
          fac = (ln0 - ztref/zbot(n)*(bn - bh))/bh * stable &
               &   + (ln0 - ln3)/bh * (1.0_R8-stable)
          fac = min(max(fac,0.0_R8),1.0_R8)

          !--- actual interpolation
          tref(n) = ts(n) + (tbot(n) - ts(n))*fac
          qref(n) = qbot(n) - delq*fac

       endif
    enddo

    if (dbug > 0) then
       do n = 1,lsize
          if (mask(n) /= 0) then
             write(logunit, F01)'n,mask  = ',n,mask(n)
             write(logunit, F01)'n,zbot  = ',n,zbot(n)
             write(logunit, F01)'n,ubot  = ',n,ubot(n)
             write(logunit, F01)'n,vbot  = ',n,vbot(n)
             write(logunit, F01)'n,thbot = ',n,thbot(n)
             write(logunit, F01)'n,qbot  = ',n,qbot(n)
             write(logunit, F01)'n,tbot  = ',n,tbot(n)
             write(logunit, F01)'n,ts    = ',n,ts(n)
             write(logunit, F01)'n,lat   = ',n,lat(n)
             write(logunit, F01)'n,sen   = ',n,sen(n)
             write(logunit, F01)'n,taux  = ',n,taux(n)
             write(logunit, F01)'n,taux  = ',n,tauy(n)
             write(logunit, F01)'n,lwup  = ',n,lwup(n)
             write(logunit, F01)'n,evap  = ',n,evap(n)
             write(logunit, F01)'n,tref  = ',n,tref(n)
             write(logunit, F01)'n,qref  = ',n,qref(n)
          end if
       end do
    end if

  end subroutine dice_flux_atmIce

end module dice_flux_atmice_mod
