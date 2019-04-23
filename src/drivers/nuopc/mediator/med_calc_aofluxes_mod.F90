module med_calc_aofluxes_mod

  !-------------------------------------------------------------------------------
  ! PURPOSE:
  !   computes atm/ocn surface fluxes
  !
  ! NOTES:
  !   o all fluxes are positive downward
  !   o net heat flux = net sw + lw up + lw down + sen + lat
  !   o here, tstar = <WT>/U*, and qstar = <WQ>/U*.
  !   o wind speeds should all be above a minimum speed (eg. 1.0 m/s)
  !
  ! ASSUMPTIONS:
  !   o Neutral 10m drag coeff: cdn = .0027/U10 + .000142 + .0000764 U10
  !   o Neutral 10m stanton number: ctn = .0327 sqrt(cdn), unstable
  !                                 ctn = .0180 sqrt(cdn), stable
  !   o Neutral 10m dalton number:  cen = .0346 sqrt(cdn)
  !   o The saturation humidity of air at T(K): qsat(T)  (kg/m^3)
  !-------------------------------------------------------------------------------
  
  use ESMF                  , only : ESMF_FAILURE, ESMF_SUCCESS, ESMF_LogWrite
  use med_constants_mod     , only : R8
  use med_internalstate_mod , only : logunit
  use water_isotopes        , only : wiso_flxoce ! calculate water isotope fluxes.
  use shr_const_mod  

  implicit none
  private ! default private

  ! public member functions:
  public :: flux_adjust_constants ! adjust constant values used in flux calculations.
  public :: flux_atmOcn           ! computes atm/ocn fluxes

  ! The follow variables are not declared as parameters so that they can be
  ! adjusted to support aquaplanet and potentially other simple model modes.
  ! The shr_flux_adjust_constants subroutine is called to set the desired
  ! values.  The default values are from shr_const_mod.  Currently they are
  ! only used by the shr_flux_atmocn and shr_flux_atmice routines.

  real(R8) :: loc_zvir   = shr_const_zvir
  real(R8) :: loc_cpdair = shr_const_cpdair
  real(R8) :: loc_cpvir  = shr_const_cpvir
  real(R8) :: loc_karman = shr_const_karman
  real(R8) :: loc_g      = shr_const_g
  real(R8) :: loc_latvap = shr_const_latvap
  real(R8) :: loc_latice = shr_const_latice
  real(R8) :: loc_stebol = shr_const_stebol

  ! These control convergence of the iterative flux calculation
  real(r8) :: flux_con_tol = 0.0_R8
  integer  :: flux_con_max_iter = 2

  ! cold air outbreak parameters  (Mahrt & Sun 1995,MWR)
  logical             :: use_coldair_outbreak_mod = .false.
  real(R8),parameter  :: alpha = 1.4_R8
  real(R8),parameter  :: maxscl =2._R8  ! maximum wind scaling for flux
  real(R8),parameter  :: td0 = -10._R8  ! start t-ts for scaling

  character(len=*), parameter :: sourcefile = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine flux_adjust_constants( zvir, cpair, cpvir, karman, gravit, &
       latvap, latice, stebol, flux_convergence_tolerance, &
       flux_convergence_max_iteration, coldair_outbreak_mod)

    ! Adjust local constants.  Used to support simple models.

    real(R8)    , optional, intent(in) :: zvir
    real(R8)    , optional, intent(in) :: cpair
    real(R8)    , optional, intent(in) :: cpvir
    real(R8)    , optional, intent(in) :: karman
    real(R8)    , optional, intent(in) :: gravit
    real(R8)    , optional, intent(in) :: latvap
    real(R8)    , optional, intent(in) :: latice
    real(R8)    , optional, intent(in) :: stebol
    real(r8)    , optional, intent(in) :: flux_convergence_tolerance
    integer     , optional, intent(in) :: flux_convergence_max_iteration
    logical     , optional, intent(in) :: coldair_outbreak_mod
    !----------------------------------------------------------------------------

    if (present(zvir))   loc_zvir   = zvir
    if (present(cpair))  loc_cpdair = cpair
    if (present(cpvir))  loc_cpvir  = cpvir
    if (present(karman)) loc_karman = karman
    if (present(gravit)) loc_g      = gravit
    if (present(latvap)) loc_latvap = latvap
    if (present(latice)) loc_latice = latice
    if (present(stebol)) loc_stebol = stebol
    if (present(flux_convergence_tolerance     )) flux_con_tol = flux_convergence_tolerance
    if (present(flux_convergence_max_iteration )) flux_con_max_iter = flux_convergence_max_iteration
    if (present(coldair_outbreak_mod           )) use_coldair_outbreak_mod = coldair_outbreak_mod

  end subroutine flux_adjust_constants

  !===============================================================================

  subroutine flux_atmOcn(nMax  , zbot    , ubot    , vbot  , thbot ,  prec_gust, gust_fac, &
                         qbot  , s16O    , sHDO    , s18O  , rbot  ,                       &
                         tbot  , us      , vs      ,                                       &
                         ts    , mask    , sen     , lat   , lwup  ,                       &
                         r16O  , rhdo    , r18O,                                           &
                         evap  , evap_16O, evap_HDO, evap_18O,                             &
                         taux  , tauy    , tref    , qref  ,                               &
                         duu10n, ustar_sv, re_sv   , ssq_sv,                               &
                         missval, rc)

    ! Calculate atm/ocn fluxes

    ! input/output variables
    integer    ,intent(in) ::       nMax       ! data vector length
    integer    ,intent(in) :: mask (nMax)      ! ocn domain mask       0 <=> out of domain
    real(R8)   ,intent(in) :: zbot (nMax)      ! atm level height      (m)
    real(R8)   ,intent(in) :: ubot (nMax)      ! atm u wind            (m/s)
    real(R8)   ,intent(in) :: vbot (nMax)      ! atm v wind            (m/s)
    real(R8)   ,intent(in) :: thbot(nMax)      ! atm potential T       (K)
    real(R8)   ,intent(in) :: qbot (nMax)      ! atm specific humidity (kg/kg)
    real(R8)   ,intent(in) :: s16O (nMax)      ! atm H216O tracer conc. (kg/kg)
    real(R8)   ,intent(in) :: sHDO (nMax)      ! atm HDO tracer conc.  (kg/kg)
    real(R8)   ,intent(in) :: s18O (nMax)      ! atm H218O tracer conc. (kg/kg)
    real(R8)   ,intent(in) :: r16O (nMax)      ! ocn H216O tracer ratio/Rstd
    real(R8)   ,intent(in) :: rHDO (nMax)      ! ocn HDO tracer ratio/Rstd
    real(R8)   ,intent(in) :: r18O (nMax)      ! ocn H218O tracer ratio/Rstd
    real(R8)   ,intent(in) :: rbot (nMax)      ! atm air density       (kg/m^3)
    real(R8)   ,intent(in) :: tbot (nMax)      ! atm T                 (K)
    real(R8)   ,intent(in) :: us   (nMax)      ! ocn u-velocity        (m/s)
    real(R8)   ,intent(in) :: vs   (nMax)      ! ocn v-velocity        (m/s)
    real(R8)   ,intent(in) :: ts   (nMax)      ! ocn temperature       (K)
    real(R8)   ,intent(in) :: prec_gust (nMax) ! atm precip for convective gustiness (kg/m^3)
    real(R8)   ,intent(in) :: gust_fac         ! wind gustiness factor

    !--- output arguments -------------------------------
    real(R8),intent(out)  ::  sen  (nMax)     ! heat flux: sensible    (W/m^2)
    real(R8),intent(out)  ::  lat  (nMax)     ! heat flux: latent      (W/m^2)
    real(R8),intent(out)  ::  lwup (nMax)     ! heat flux: lw upward   (W/m^2)
    real(R8),intent(out)  ::  evap (nMax)     ! water flux: evap  ((kg/s)/m^2)
    real(R8),intent(out)  ::  evap_16O (nMax) ! water flux: evap ((kg/s/m^2)
    real(R8),intent(out)  ::  evap_HDO (nMax) ! water flux: evap ((kg/s)/m^2)
    real(R8),intent(out)  ::  evap_18O (nMax) ! water flux: evap ((kg/s/m^2)
    real(R8),intent(out)  ::  taux (nMax)     ! surface stress, zonal      (N)
    real(R8),intent(out)  ::  tauy (nMax)     ! surface stress, maridional (N)
    real(R8),intent(out)  ::  tref (nMax)     ! diag:  2m ref height T     (K)
    real(R8),intent(out)  ::  qref (nMax)     ! diag:  2m ref humidity (kg/kg)
    real(R8),intent(out)  :: duu10n(nMax)     ! diag: 10m wind speed squared (m/s)^2
    integer ,intent(out)  :: rc

    real(R8),intent(out),optional :: ustar_sv(nMax) ! diag: ustar
    real(R8),intent(out),optional :: re_sv   (nMax) ! diag: sqrt of exchange coefficient (water)
    real(R8),intent(out),optional :: ssq_sv  (nMax) ! diag: sea surface humidity  (kg/kg)
    real(R8),intent(in) ,optional :: missval        ! masked value

    !--- local constants --------------------------------
    real(R8),parameter :: umin  =  0.5_R8 ! minimum wind speed       (m/s)
    real(R8),parameter :: zref  = 10.0_R8 ! reference height           (m)
    real(R8),parameter :: ztref =  2.0_R8 ! reference height for air T (m)

    !--- local variables --------------------------------
    integer     :: n      ! vector loop index
    integer     :: iter
    real(R8)    :: vmag   ! surface wind magnitude   (m/s)
    real(R8)    :: vmag_old   ! surface wind magnitude without gustiness (m/s)
    real(R8)    :: ssq    ! sea surface humidity     (kg/kg)
    real(R8)    :: delt   ! potential T difference   (K)
    real(R8)    :: delq   ! humidity difference      (kg/kg)
    real(R8)    :: stable ! stability factor
    real(R8)    :: rdn    ! sqrt of neutral exchange coeff (momentum)
    real(R8)    :: rhn    ! sqrt of neutral exchange coeff (heat)
    real(R8)    :: ren    ! sqrt of neutral exchange coeff (water)
    real(R8)    :: rd     ! sqrt of exchange coefficient (momentum)
    real(R8)    :: rh     ! sqrt of exchange coefficient (heat)
    real(R8)    :: re     ! sqrt of exchange coefficient (water)
    real(R8)    :: ustar  ! ustar
    real(r8)    :: ustar_prev
    real(R8)    :: qstar  ! qstar
    real(R8)    :: tstar  ! tstar
    real(R8)    :: hol    ! H (at zbot) over L
    real(R8)    :: xsq    ! ?
    real(R8)    :: xqq    ! ?
    real(R8)    :: psimh  ! stability function at zbot (momentum)
    real(R8)    :: psixh  ! stability function at zbot (heat and water)
    real(R8)    :: psix2  ! stability function at ztref reference height
    real(R8)    :: alz    ! ln(zbot/zref)
    real(R8)    :: al2    ! ln(zref/ztref)
    real(R8)    :: u10n   ! 10m neutral wind
    real(R8)    :: tau    ! stress at zbot
    real(R8)    :: cp     ! specific heat of moist air
    real(R8)    :: fac    ! vertical interpolation factor
    real(R8)    :: spval  ! local missing value

    !--- local functions --------------------------------
    real(R8)    :: qsat   ! function: the saturation humididty of air (kg/m^3)
    real(R8)    :: cdn    ! function: neutral drag coeff at 10m
    real(R8)    :: psimhu ! function: unstable part of psimh
    real(R8)    :: psixhu ! function: unstable part of psimx
    real(R8)    :: ugust  ! function: gustiness as a function of convective rainfall
    real(R8)    :: Umps   ! dummy arg ~ wind velocity (m/s)
    real(R8)    :: Tk     ! dummy arg ~ temperature (K)
    real(R8)    :: xd     ! dummy arg ~ ?
    real(R8)    :: gprec  ! dummy arg ~ ?
    !--- for cold air outbreak calc --------------------------------
    real(R8)    :: tdiff(nMax)               ! tbot - ts
    real(R8)    :: vscl

    qsat(Tk)   = 640380.0_R8 / exp(5107.4_R8/Tk)
    cdn(Umps)  =   0.0027_R8 / Umps + 0.000142_R8 + 0.0000764_R8 * Umps
    psimhu(xd) = log((1.0_R8+xd*(2.0_R8+xd))*(1.0_R8+xd*xd)/8.0_R8) - 2.0_R8*atan(xd) + 1.571_R8
    psixhu(xd) = 2.0_R8 * log((1.0_R8 + xd*xd)/2.0_R8)

    ! Convective gustiness appropriate for input precipitation.
    ! Following Redelsperger et al. (2000, J. Clim)
    ! Ug = log(1.0+6.69R-0.476R^2)
    ! Coefficients X by 8640 for mm/s (from cam) -> cm/day (for above forumla)
    ugust(gprec) = gust_fac*log(1._R8+57801.6_R8*gprec-3.55332096e7_R8*(gprec**2.0_R8))

    !--- formats ----------------------------------------
    character(*),parameter :: subName = '(shr_flux_atmOcn) '
    character(*),parameter ::   F00 = "('(shr_flux_atmOcn) ',4a)"
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (present(missval)) then
       spval = missval
    else
       spval = shr_const_spval
    endif
    u10n = spval
    rh = spval
    psixh = spval
    hol=spval

    !--- for cold air outbreak calc --------------------------------
    tdiff= tbot - ts

    al2 = log(zref/ztref)
    DO n=1,nMax
       if (mask(n) /= 0) then

          !--- compute some needed quantities ---

          ! old version
          !vmag   = max(umin, sqrt( (ubot(n)-us(n))**2 + (vbot(n)-vs(n))**2) )

          !--- vmag+ugust (convective gustiness) Limit to a max precip 6 cm/day = 0.00069444 m/s.
          !--- reverts to original formula if gust_fac=0

          !PMA saves vmag_old for taux tauy computation

          vmag_old    = max(umin, sqrt( (ubot(n)-us(n))**2 + (vbot(n)-vs(n))**2) )

          if (gust_fac .gt. 1.e-12_R8) then
             vmag       = max(umin, sqrt( (ubot(n)-us(n))**2 + (vbot(n)-vs(n))**2) + ugust(min(prec_gust(n),6.94444e-4_R8)))
          else
             vmag       = vmag_old
          endif
          if (use_coldair_outbreak_mod) then
             ! Cold Air Outbreak Modification:
             ! Increase windspeed for negative tbot-ts
             ! based on Mahrt & Sun 1995,MWR

             if (tdiff(n).lt.td0) then
                vscl=min((1._R8+alpha*(abs(tdiff(n)-td0)**0.5_R8/abs(vmag))),maxscl)
                vmag=vmag*vscl
                vmag_old=vmag_old*vscl
             endif
          endif
          ssq    = 0.98_R8 * qsat(ts(n)) / rbot(n)   ! sea surf hum (kg/kg)
          delt   = thbot(n) - ts(n)                  ! pot temp diff (K)
          delq   = qbot(n) - ssq                     ! spec hum dif (kg/kg)
          alz    = log(zbot(n)/zref)
          cp     = loc_cpdair*(1.0_R8 + loc_cpvir*ssq)

          !------------------------------------------------------------
          ! first estimate of Z/L and ustar, tstar and qstar
          !------------------------------------------------------------
          !--- neutral coefficients, z/L = 0.0 ---
          stable = 0.5_R8 + sign(0.5_R8 , delt)
          rdn    = sqrt(cdn(vmag))
          rhn    = (1.0_R8-stable) * 0.0327_R8 + stable * 0.018_R8
          ren    = 0.0346_R8

          !--- ustar, tstar, qstar ---
          ustar = rdn * vmag
          tstar = rhn * delt
          qstar = ren * delq
          ustar_prev = ustar*2.0_R8
          iter = 0
          do while( abs((ustar - ustar_prev)/ustar) > flux_con_tol .and. iter < flux_con_max_iter)
             iter = iter + 1
             ustar_prev = ustar
             !--- compute stability & evaluate all stability functions ---
             hol  = loc_karman*loc_g*zbot(n)*  &
                  (tstar/thbot(n)+qstar/(1.0_R8/loc_zvir+qbot(n)))/ustar**2
             hol  = sign( min(abs(hol),10.0_R8), hol )
             stable = 0.5_R8 + sign(0.5_R8 , hol)
             xsq    = max(sqrt(abs(1.0_R8 - 16.0_R8*hol)) , 1.0_R8)
             xqq    = sqrt(xsq)
             psimh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psimhu(xqq)
             psixh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)

             !--- shift wind speed using old coefficient ---
             rd   = rdn / (1.0_R8 + rdn/loc_karman*(alz-psimh))
             u10n = vmag * rd / rdn

             !--- update transfer coeffs at 10m and neutral stability ---
             rdn = sqrt(cdn(u10n))
             ren = 0.0346_R8
             rhn = (1.0_R8-stable)*0.0327_R8 + stable * 0.018_R8

             !--- shift all coeffs to measurement height and stability ---
             rd = rdn / (1.0_R8 + rdn/loc_karman*(alz-psimh))
             rh = rhn / (1.0_R8 + rhn/loc_karman*(alz-psixh))
             re = ren / (1.0_R8 + ren/loc_karman*(alz-psixh))

             !--- update ustar, tstar, qstar using updated, shifted coeffs --
             ustar = rd * vmag
             tstar = rh * delt
             qstar = re * delq
          enddo
          if (iter < 1) then
             write(logunit,*) ustar,ustar_prev,flux_con_tol,flux_con_max_iter
             call ESMF_LogWrite('No iterations performed - ERROR in med_calc_aofluxe')
             rc=ESMF_Failure
             return
          end if

          !------------------------------------------------------------
          ! compute the fluxes
          !------------------------------------------------------------

          tau = rbot(n) * ustar * ustar

          !--- momentum flux ---
          taux(n) = tau * (ubot(n)-us(n)) / vmag_old !PMA uses vmag_old for taux
          tauy(n) = tau * (vbot(n)-vs(n)) / vmag_old !    tauy c20170620

          !--- heat flux ---
          sen (n) =          cp * tau * tstar / ustar
          lat (n) =  loc_latvap * tau * qstar / ustar
          lwup(n) = -loc_stebol * ts(n)**4

          !--- water flux ---
          evap(n) = lat(n)/loc_latvap

          !---water isotope flux ---
          call wiso_flxoce(2,rbot(n),zbot(n),s16O(n),ts(n),r16O(n),ustar,re,ssq,evap_16O(n), &
               qbot(n),evap(n))
          call wiso_flxoce(3,rbot(n),zbot(n),sHDO(n),ts(n),rHDO(n),ustar,re,ssq, evap_HDO(n),&
               qbot(n),evap(n))
          call wiso_flxoce(4,rbot(n),zbot(n),s18O(n),ts(n),r18O(n),ustar,re,ssq, evap_18O(n), &
               qbot(n),evap(n))

          !------------------------------------------------------------
          ! compute diagnositcs: 2m ref T & Q, 10m wind speed squared
          !------------------------------------------------------------
          hol = hol*ztref/zbot(n)
          xsq = max( 1.0_R8, sqrt(abs(1.0_R8-16.0_R8*hol)) )
          xqq = sqrt(xsq)
          psix2   = -5.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)
          fac     = (rh/loc_karman) * (alz + al2 - psixh + psix2 )
          tref(n) = thbot(n) - delt*fac
          tref(n) = tref(n) - 0.01_R8*ztref   ! pot temp to temp correction
          fac     = (re/loc_karman) * (alz + al2 - psixh + psix2 )
          qref(n) =  qbot(n) - delq*fac

          duu10n(n) = u10n*u10n ! 10m wind speed squared

          !------------------------------------------------------------
          ! optional diagnostics, needed for water tracer fluxes (dcn)
          !------------------------------------------------------------
          if (present(ustar_sv)) ustar_sv(n) = ustar
          if (present(re_sv   )) re_sv(n)    = re
          if (present(ssq_sv  )) ssq_sv(n)   = ssq

       else
          !------------------------------------------------------------
          ! no valid data here -- out of domain
          !------------------------------------------------------------
          sen   (n) = spval  ! sensible         heat flux  (W/m^2)
          lat   (n) = spval  ! latent           heat flux  (W/m^2)
          lwup  (n) = spval  ! long-wave upward heat flux  (W/m^2)
          evap  (n) = spval  ! evaporative water flux ((kg/s)/m^2)
          evap_16O (n) = spval !water tracer flux (kg/s)/m^2)
          evap_HDO (n) = spval !HDO tracer flux  (kg/s)/m^2)
          evap_18O (n) = spval !H218O tracer flux (kg/s)/m^2)
          taux  (n) = spval  ! x surface stress (N)
          tauy  (n) = spval  ! y surface stress (N)
          tref  (n) = spval  !  2m reference height temperature (K)
          qref  (n) = spval  !  2m reference height humidity (kg/kg)
          duu10n(n) = spval  ! 10m wind speed squared (m/s)^2

          if (present(ustar_sv)) ustar_sv(n) = spval
          if (present(re_sv   )) re_sv   (n) = spval
          if (present(ssq_sv  )) ssq_sv  (n) = spval
       endif
    end DO

  end subroutine flux_atmOcn

end module med_calc_aofluxes_mod
