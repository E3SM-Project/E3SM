! !MODULE: flux_mod -- CCSM shared flux calculations.
!
! !DESCRIPTION:
!
!     CCSM shared flux calculations.
!
! !REVISION HISTORY:
!     2006-Nov-07 - B. Kauffman - first version, code taken/migrated from cpl6
!
! !INTERFACE: ------------------------------------------------------------------

module shr_flux_mod

! !USES:

   use shr_kind_mod, only : R8=>SHR_KIND_R8, IN=>SHR_KIND_IN    ! shared kinds
   use shr_const_mod   ! shared constants
   use shr_sys_mod     ! shared system routines
   use shr_log_mod, only: s_loglev  => shr_log_Level
   use shr_log_mod, only: s_logunit => shr_log_Unit
   use shr_log_mod, only: errMsg => shr_log_errMsg

   implicit none

   private ! default private

! !PUBLIC TYPES:

  ! none

! !PUBLIC MEMBER FUNCTIONS:

   public :: shr_flux_atmOcn      ! computes atm/ocn fluxes
   public :: shr_flux_atmOcn_diurnal   ! computes atm/ocn fluxes with diurnal cycle
   public :: shr_flux_atmOcn_UA   ! computes atm/ocn fluxes using University of
                                  ! Arizona algorithm (Zeng et al., 1998)
   public :: shr_flux_atmIce      ! computes atm/ice fluxes
   public :: shr_flux_MOstability ! boundary layer stability scales/functions
   public :: shr_flux_adjust_constants ! adjust constant values used in flux calculations.

! !PUBLIC DATA MEMBERS:

  integer(IN),parameter,public :: shr_flux_MOwScales   = 1 ! w scales  option
  integer(IN),parameter,public :: shr_flux_MOfunctions = 2 ! functions option
  real   (R8),parameter,public :: shr_flux_MOgammaM = 3.59_R8
  real   (R8),parameter,public :: shr_flux_MOgammaS = 7.86_R8

!EOP

   !--- rename kinds for local readability only ---

   integer,parameter :: debug = 0 ! internal debug level

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
   real(R8) :: loc_tkfrz  = shr_const_tkfrz

   ! These control convergence of the iterative flux calculation
   ! (For Large and Pond scheme only; not UA or COARE).
   real(r8) :: flux_con_tol = 0.0_R8
   integer(IN) :: flux_con_max_iter = 2

   character(len=*), parameter :: sourcefile = &
     __FILE__

   !--- cold air outbreak parameters  (Mahrt & Sun 1995,MWR) -------------
   logical :: use_coldair_outbreak_mod = .false.
   real(R8),parameter    :: alpha = 1.4_R8
   real(R8),parameter    :: maxscl =2._R8  ! maximum wind scaling for flux
   real(R8),parameter    :: td0 = -10._R8   ! start t-ts for scaling

!===============================================================================
contains
!===============================================================================
!===============================================================================
subroutine shr_flux_adjust_constants( &
   zvir, cpair, cpvir, karman, gravit, &
   latvap, latice, stebol, flux_convergence_tolerance, &
   flux_convergence_max_iteration, &
   coldair_outbreak_mod)

   ! Adjust local constants.  Used to support simple models.

   real(R8), optional, intent(in) :: zvir
   real(R8), optional, intent(in) :: cpair
   real(R8), optional, intent(in) :: cpvir
   real(R8), optional, intent(in) :: karman
   real(R8), optional, intent(in) :: gravit
   real(R8), optional, intent(in) :: latvap
   real(R8), optional, intent(in) :: latice
   real(R8), optional, intent(in) :: stebol
   real(r8), optional, intent(in)  :: flux_convergence_tolerance
   integer(in), optional, intent(in) :: flux_convergence_max_iteration
   logical, optional, intent(in) :: coldair_outbreak_mod
   !----------------------------------------------------------------------------

   if (present(zvir))   loc_zvir   = zvir
   if (present(cpair))  loc_cpdair = cpair
   if (present(cpvir))  loc_cpvir  = cpvir
   if (present(karman)) loc_karman = karman
   if (present(gravit)) loc_g      = gravit
   if (present(latvap)) loc_latvap = latvap
   if (present(latice)) loc_latice = latice
   if (present(stebol)) loc_stebol = stebol
   if (present(flux_convergence_tolerance)) flux_con_tol = flux_convergence_tolerance
   if (present(flux_convergence_max_iteration)) flux_con_max_iter = flux_convergence_max_iteration
   if(present(coldair_outbreak_mod)) use_coldair_outbreak_mod = coldair_outbreak_mod
end subroutine shr_flux_adjust_constants
!===============================================================================
! !BOP =========================================================================
!
! !IROUTINE: shr_flux_atmOcn -- internal atm/ocn flux calculation
!
! !DESCRIPTION:
!
!     Internal atm/ocn flux calculation
!
! !REVISION HISTORY:
!     2002-Jun-10 - B. Kauffman - code migrated from cpl5 to cpl6
!     2003-Apr-02 - B. Kauffman - taux & tauy now utilize ocn velocity
!     2003-Apr-02 - B. Kauffman - tref,qref,duu10n mods as per Bill Large
!     2006-Nov-07 - B. Kauffman - code migrated from cpl6 to share
!
!     2011-Mar-13 - J. Nusbaumer - Water Isotope ocean flux added.
!     2019-May-16 - Jack Reeves Eyre (UA) and Kai Zhang (PNNL) - Added COARE/Fairall surface flux scheme option (ocn_surface_flux_scheme .eq. 1) based on code from Thomas Toniazzo (Bjerknes Centre, Bergen) ”
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE shr_flux_atmOcn(nMax  ,zbot  ,ubot  ,vbot  ,thbot ,   &
           &               qbot  ,s16O  ,sHDO  ,s18O  ,rbot  ,   &
           &               tbot  ,us    ,vs    ,   &
           &               ts    ,mask  , seq_flux_atmocn_minwind, &
           &               sen   ,lat   ,lwup  ,   &
           &               r16O, rhdo, r18O, &
           &               evap  ,evap_16O, evap_HDO, evap_18O, &
           &               taux  ,tauy  ,tref  ,qref  ,   &
           &               ocn_surface_flux_scheme, &
           &               duu10n,  ustar_sv   ,re_sv ,ssq_sv,   &
           &               missval    )

! !USES:

   use water_isotopes, only: wiso_flxoce !subroutine used to calculate water isotope fluxes.

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   !--- input arguments --------------------------------
   integer(IN),intent(in) ::       nMax  ! data vector length
   integer(IN),intent(in) :: mask (nMax) ! ocn domain mask       0 <=> out of domain
   integer(IN),intent(in) :: ocn_surface_flux_scheme
   real(R8)   ,intent(in) :: zbot (nMax) ! atm level height      (m)
   real(R8)   ,intent(in) :: ubot (nMax) ! atm u wind            (m/s)
   real(R8)   ,intent(in) :: vbot (nMax) ! atm v wind            (m/s)
   real(R8)   ,intent(in) :: thbot(nMax) ! atm potential T       (K)
   real(R8)   ,intent(in) :: qbot (nMax) ! atm specific humidity (kg/kg)
   real(R8)   ,intent(in) :: s16O (nMax) ! atm H216O tracer conc. (kg/kg)
   real(R8)   ,intent(in) :: sHDO (nMax) ! atm HDO tracer conc.  (kg/kg)
   real(R8)   ,intent(in) :: s18O (nMax) ! atm H218O tracer conc. (kg/kg)
   real(R8)   ,intent(in) :: r16O (nMax) ! ocn H216O tracer ratio/Rstd
   real(R8)   ,intent(in) :: rHDO (nMax) ! ocn HDO tracer ratio/Rstd
   real(R8)   ,intent(in) :: r18O (nMax) ! ocn H218O tracer ratio/Rstd
   real(R8)   ,intent(in) :: rbot (nMax) ! atm air density       (kg/m^3)
   real(R8)   ,intent(in) :: tbot (nMax) ! atm T                 (K)
   real(R8)   ,intent(in) :: us   (nMax) ! ocn u-velocity        (m/s)
   real(R8)   ,intent(in) :: vs   (nMax) ! ocn v-velocity        (m/s)
   real(R8)   ,intent(in) :: ts   (nMax) ! ocn temperature       (K)
   real(R8)   ,intent(in) :: seq_flux_atmocn_minwind        ! minimum wind speed for atmocn      (m/s)

   !--- output arguments -------------------------------
   real(R8),intent(out)  ::  sen  (nMax) ! heat flux: sensible    (W/m^2)
   real(R8),intent(out)  ::  lat  (nMax) ! heat flux: latent      (W/m^2)
   real(R8),intent(out)  ::  lwup (nMax) ! heat flux: lw upward   (W/m^2)
   real(R8),intent(out)  ::  evap (nMax) ! water flux: evap  ((kg/s)/m^2)
   real(R8),intent(out)  ::  evap_16O (nMax) ! water flux: evap ((kg/s/m^2)
   real(R8),intent(out)  ::  evap_HDO (nMax) ! water flux: evap ((kg/s)/m^2)
   real(R8),intent(out)  ::  evap_18O (nMax) ! water flux: evap ((kg/s/m^2)
   real(R8),intent(out)  ::  taux (nMax) ! surface stress, zonal      (N)
   real(R8),intent(out)  ::  tauy (nMax) ! surface stress, maridional (N)
   real(R8),intent(out)  ::  tref (nMax) ! diag:  2m ref height T     (K)
   real(R8),intent(out)  ::  qref (nMax) ! diag:  2m ref humidity (kg/kg)
   real(R8),intent(out)  :: duu10n(nMax) ! diag: 10m wind speed squared (m/s)^2

   real(R8),intent(out),optional :: ustar_sv(nMax) ! diag: ustar
   real(R8),intent(out),optional :: re_sv   (nMax) ! diag: sqrt of exchange coefficient (water)
   real(R8),intent(out),optional :: ssq_sv  (nMax) ! diag: sea surface humidity  (kg/kg)

   real(R8),intent(in) ,optional :: missval        ! masked value

! !EOP

   !--- local constants --------------------------------
   real(R8),parameter :: zref  = 10.0_R8 ! reference height           (m)
   real(R8),parameter :: ztref =  2.0_R8 ! reference height for air T (m)
!!++ Large only
   !real(R8),parameter :: cexcd  = 0.0346_R8 ! ratio Ch(water)/CD
   !real(R8),parameter :: chxcds = 0.018_R8  ! ratio Ch(heat)/CD for stable case
   !real(R8),parameter :: chxcdu = 0.0327_R8 ! ratio Ch(heat)/CD for unstable case
!!++ COARE only
   real(R8),parameter :: zpbl =700.0_R8 ! PBL depth [m] for gustiness parametriz.

   !--- local variables --------------------------------
   integer(IN) :: n      ! vector loop index
   integer(IN) :: iter
   real(R8)    :: vmag   ! surface wind magnitude   (m/s)
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
   real(r8)     :: ustar_prev
   real(R8)    :: qstar  ! qstar
   real(R8)    :: tstar  ! tstar
   real(R8)    :: hol    ! H (at zbot) over L
   real(R8)    :: xsq    ! ?
   real(R8)    :: xqq    ! ?
 !!++ Large only
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
!!++ COARE only
   real(R8)    :: zo,zot,zoq      ! roughness lengths
   real(R8)    :: hsb,hlb         ! sens & lat heat flxs at zbot
   real(R8) :: trf,qrf,urf,vrf ! reference-height quantities


    !--- local functions --------------------------------
   real(R8)    :: qsat   ! function: the saturation humididty of air (kg/m^3)
!!++ Large only (formula v*=[c4/U10+c5+c6*U10]*U10 in Large et al. 1994)
   real(R8)    :: cdn    ! function: neutral drag coeff at 10m
!!++ Large only (stability functions)
   real(R8)    :: psimhu ! function: unstable part of psimh
   real(R8)    :: psixhu ! function: unstable part of psimx
   real(R8)    :: Umps   ! dummy arg ~ wind velocity (m/s)
   real(R8)    :: Tk     ! dummy arg ~ temperature (K)
   real(R8)    :: xd     ! dummy arg ~ ?
   !--- for cold air outbreak calc --------------------------------
   real(R8)    :: tdiff(nMax)               ! tbot - ts
   real(R8)    :: vscl


   qsat(Tk)   = 640380.0_R8 / exp(5107.4_R8/Tk)
   cdn(Umps)  =   0.0027_R8 / Umps + 0.000142_R8 + 0.0000764_R8 * Umps
   psimhu(xd) = log((1.0_R8+xd*(2.0_R8+xd))*(1.0_R8+xd*xd)/8.0_R8) - 2.0_R8*atan(xd) + 1.571_R8
   psixhu(xd) = 2.0_R8 * log((1.0_R8 + xd*xd)/2.0_R8)

   !--- formats ----------------------------------------
   character(*),parameter :: subName = '(shr_flux_atmOcn) '
   character(*),parameter ::   F00 = "('(shr_flux_atmOcn) ',4a)"

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
!  Large:
!   o Neutral 10m drag coeff: cdn = .0027/U10 + .000142 + .0000764 U10
!   o Neutral 10m stanton number: ctn = .0327 sqrt(cdn), unstable
!                                 ctn = .0180 sqrt(cdn), stable
!   o Neutral 10m dalton number:  cen = .0346 sqrt(cdn)
!   o The saturation humidity of air at T(K): qsat(T)  (kg/m^3)
!  COARE:
!   o use COAREv3.0 function (tht 22/11/2013)
!-------------------------------------------------------------------------------

   if (debug > 0 .and. s_loglev > 0) write(s_logunit,F00) "enter"

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

!!.................................................................
!! ocn_surface_flux_scheme = 0 : Default CESM1.2
!!                         = 1 : COARE algorithm
!!                         = 2 : UA algorithm (separate subroutine)
!!.................................................................

   ! Default flux scheme.
   if (ocn_surface_flux_scheme .eq. 0) then

   al2 = log(zref/ztref)
   DO n=1,nMax
     if (mask(n) /= 0) then

        !--- compute some needed quantities ---
        vmag   = max(seq_flux_atmocn_minwind, sqrt( (ubot(n)-us(n))**2 + (vbot(n)-vs(n))**2) )
        if (use_coldair_outbreak_mod) then
            ! Cold Air Outbreak Modification:
            ! Increase windspeed for negative tbot-ts
            ! based on Mahrt & Sun 1995,MWR

            if (tdiff(n).lt.td0) then
               vscl=min((1._R8+alpha*(abs(tdiff(n)-td0)**0.5_R8/abs(vmag))),maxscl)
               vmag=vmag*vscl
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
                 !(1.0_R8-stable) * chxcdu + stable * chxcds
        ren    = 0.0346_R8 !cexcd

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
           ren = 0.0346_R8 !cexcd
           rhn = (1.0_R8-stable)*0.0327_R8 + stable * 0.018_R8
                 !(1.0_R8-stable) * chxcdu + stable * chxcds

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
           write(s_logunit,*) ustar,ustar_prev,flux_con_tol,flux_con_max_iter
           call shr_sys_abort('No iterations performed ' // errMsg(sourcefile, __LINE__))
        end if
        !------------------------------------------------------------
        ! compute the fluxes
        !------------------------------------------------------------

        tau = rbot(n) * ustar * ustar

        !--- momentum flux ---
        taux(n) = tau * (ubot(n)-us(n)) / vmag
        tauy(n) = tau * (vbot(n)-vs(n)) / vmag

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
   ENDDO

   else if (ocn_surface_flux_scheme .eq. 1) then
    !!.................................
    !! use COARE algorithm
    !!.................................


    DO n=1,nMax
     if (mask(n) /= 0) then

        !--- compute some needed quantities ---
        vmag    = max(seq_flux_atmocn_minwind, sqrt( (ubot(n)-us(n))**2 + (vbot(n)-vs(n))**2) )

         if (use_coldair_outbreak_mod) then
            ! Cold Air Outbreak Modification:
            ! Increase windspeed for negative tbot-ts
            ! based on Mahrt & Sun 1995,MWR

            if (tdiff(n).lt.td0) then
               vscl=min((1._R8+alpha*(abs(tdiff(n)-td0)**0.5_R8/abs(vmag))),maxscl)
               vmag=vmag*vscl
            endif
         endif
        ssq    = 0.98_R8 * qsat(ts(n)) / rbot(n)   ! sea surf hum (kg/kg)

        call cor30a(ubot(n),vbot(n),tbot(n),qbot(n),rbot(n) &  ! in atm params
                 & ,us(n),vs(n),ts(n),ssq                   &  ! in surf params
                 & ,zpbl,zbot(n),zbot(n),zref,ztref,ztref   &  ! in heights
                 & ,tau,hsb,hlb                             &  ! out: fluxes
                 & ,zo,zot,zoq,hol,ustar,tstar,qstar        &  ! out: ss scales
                 & ,rd,rh,re                                &  ! out: exch. coeffs
                 & ,trf,qrf,urf,vrf) ! out: reference-height params

! for the sake of maintaining same defs
        hol=zbot(n)/hol
        rd=sqrt(rd)
        rh=sqrt(rh)
        re=sqrt(re)

        !--- momentum flux ---
        taux(n) = tau * (ubot(n)-us(n)) / vmag
        tauy(n) = tau * (vbot(n)-vs(n)) / vmag

        !--- heat flux ---
        sen (n) =  hsb
        lat (n) =  hlb
        lwup(n) = -shr_const_stebol * ts(n)**4

        !--- water flux ---
        evap(n) = lat(n)/shr_const_latvap

        !---water isotope flux ---
        call wiso_flxoce(2,rbot(n),zbot(n),s16O(n),ts(n),r16O(n),ustar,re,ssq, evap_16O(n), &
                         qbot(n),evap(n))
        call wiso_flxoce(3,rbot(n),zbot(n),sHDO(n),ts(n),rHDO(n),ustar,re,ssq, evap_HDO(n),&
                         qbot(n),evap(n))
        call wiso_flxoce(4,rbot(n),zbot(n),s18O(n),ts(n),r18O(n),ustar,re,ssq, evap_18O(n), &
                         qbot(n),evap(n))

        !------------------------------------------------------------
        ! compute diagnositcs: 2m ref T & Q, 10m wind speed squared
        !------------------------------------------------------------
        tref(n) = trf
        qref(n) = qrf
        duu10n(n) = urf**2+vrf**2

        !------------------------------------------------------------
        ! optional diagnostics, needed for water tracer fluxes (dcn)
        !------------------------------------------------------------
        if (present(ustar_sv)) ustar_sv(n) = ustar
        if (present(re_sv   )) re_sv(n)    = re
        if (present(ssq_sv )) ssq_sv(n) = ssq

     else
        !------------------------------------------------------------
        ! no valid data here -- out of domain
        !------------------------------------------------------------
        sen      (n) = spval  ! sensible         heat flux  (W/m^2)
        lat      (n) = spval  ! latent           heat flux  (W/m^2)
        lwup     (n) = spval  ! long-wave upward heat flux  (W/m^2)
        evap     (n) = spval  ! evaporative water flux ((kg/s)/m^2)
        evap_16O (n) = spval  ! water tracer flux (kg/s)/m^2)
        evap_HDO (n) = spval  ! HDO tracer flux  (kg/s)/m^2)
        evap_18O (n) = spval  ! H218O tracer flux (kg/s)/m^2)
        taux     (n) = spval  ! x surface stress (N)
        tauy     (n) = spval  ! y surface stress (N)
        tref     (n) = spval  !  2m reference height temperature (K)
        qref     (n) = spval  !  2m reference height humidity (kg/kg)
        duu10n   (n) = spval  ! 10m wind speed squared (m/s)^2

        if (present(ustar_sv)) ustar_sv(n) = spval
        if (present(re_sv   )) re_sv   (n) = spval
        if (present(ssq_sv  )) ssq_sv  (n) = spval
     endif
    ENDDO

   else

      call shr_sys_abort(subName//" subroutine shr_flux_atmOcn requires ocn_surface_flux_scheme = 0 or 1")

   endif  !! ocn_surface_flux_scheme


END subroutine shr_flux_atmOcn


!===============================================================================
! !BOP =========================================================================
!
! !IROUTINE: shr_flux_atmOcn_UA -- internal atm/ocn flux calculation
!
! !DESCRIPTION:
!
!     Internal atm/ocn flux calculation
!     using University of Arizona method.
!
!     Reference:
!         Zeng, X., M. Zhao, and R.E. Dickinson, 1998: Intercomparison of Bulk
!             Aerodynamic Algorithms for the Computation of Sea Surface Fluxes
!             Using TOGA COARE and TAO Data. J. Climate, 11, 2628–2644,
!             https://doi.org/10.1175/1520-0442(1998)011<2628%3AIOBAAF>2.0.CO%3B2
!
!     Equation numbers are from this paper.
!
! !REVISION HISTORY:
!     2017-Aug-28 - J. Reeves Eyre - code re-written for E3SM
!     2018-Oct-30 - J. Reeves Eyre - bug fix and add
!                   convective gustiness.
!     2019-May-08 - J. Reeves Eyre - remove convective gustiness
!                   and add cold air outbreak modification.
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE shr_flux_atmOcn_UA(   &
           &               nMax  ,zbot  ,ubot  ,vbot  ,thbot ,  &
           &               qbot  ,s16O  ,sHDO  ,s18O  ,rbot  ,   &
           &               tbot  , pslv ,us    , vs   ,   &
           &               ts    ,mask  ,sen   ,lat   ,lwup  ,   &
           &               r16O, rhdo, r18O, &
           &               evap  ,evap_16O, evap_HDO, evap_18O, &
           &               taux  ,tauy  ,tref  ,qref  ,   &
           &               duu10n,  ustar_sv   ,re_sv ,ssq_sv,   &
           &               missval    )


! !USES:
   use water_isotopes, only: wiso_flxoce !subroutine used to calculate water isotope fluxes.

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   !--- input arguments --------------------------------
   integer(IN),intent(in) ::       nMax  ! data vector length
   integer(IN),intent(in) :: mask (nMax) ! ocn domain mask       0 <=> out of domain
   real(R8)   ,intent(in) :: zbot (nMax) ! atm level height      (m)
   real(R8)   ,intent(in) :: ubot (nMax) ! atm u wind            (m/s)
   real(R8)   ,intent(in) :: vbot (nMax) ! atm v wind            (m/s)
   real(R8)   ,intent(in) :: thbot(nMax) ! atm potential T       (K)
   real(R8)   ,intent(in) :: qbot (nMax) ! atm specific humidity (kg/kg)
   real(R8)   ,intent(in) :: s16O (nMax) ! atm H216O tracer conc. (kg/kg)
   real(R8)   ,intent(in) :: sHDO (nMax) ! atm HDO tracer conc.  (kg/kg)
   real(R8)   ,intent(in) :: s18O (nMax) ! atm H218O tracer conc. (kg/kg)
   real(R8)   ,intent(in) :: r16O (nMax) ! ocn H216O tracer ratio/Rstd
   real(R8)   ,intent(in) :: rHDO (nMax) ! ocn HDO tracer ratio/Rstd
   real(R8)   ,intent(in) :: r18O (nMax) ! ocn H218O tracer ratio/Rstd
   real(R8)   ,intent(in) :: rbot (nMax) ! atm air density       (kg/m^3)
   real(R8)   ,intent(in) :: tbot (nMax) ! atm T                 (K)
   real(R8)   ,intent(in) :: pslv (nMax) ! sea level pressure    (Pa)
   real(R8)   ,intent(in) :: us   (nMax) ! ocn u-velocity        (m/s)
   real(R8)   ,intent(in) :: vs   (nMax) ! ocn v-velocity        (m/s)
   real(R8)   ,intent(in) :: ts   (nMax) ! ocn temperature       (K)

   !--- output arguments -------------------------------
   real(R8),intent(out)  ::  sen  (nMax) ! heat flux: sensible    (W/m^2)
   real(R8),intent(out)  ::  lat  (nMax) ! heat flux: latent      (W/m^2)
   real(R8),intent(out)  ::  lwup (nMax) ! heat flux: lw upward   (W/m^2)
   real(R8),intent(out)  ::  evap (nMax) ! water flux: evap  ((kg/s)/m^2)
   real(R8),intent(out)  ::  evap_16O (nMax) ! water flux: evap ((kg/s/m^2)
   real(R8),intent(out)  ::  evap_HDO (nMax) ! water flux: evap ((kg/s)/m^2)
   real(R8),intent(out)  ::  evap_18O (nMax) ! water flux: evap ((kg/s/m^2)
   real(R8),intent(out)  ::  taux (nMax) ! surface stress, zonal      (N)
   real(R8),intent(out)  ::  tauy (nMax) ! surface stress, maridional (N)
   real(R8),intent(out)  ::  tref (nMax) ! diag:  2m ref height T     (K)
   real(R8),intent(out)  ::  qref (nMax) ! diag:  2m ref humidity (kg/kg)
   real(R8),intent(out)  :: duu10n(nMax) ! diag: 10m wind speed squared (m/s)^2

   real(R8),intent(out),optional :: ustar_sv(nMax) ! diag: ustar
   real(R8),intent(out),optional :: re_sv   (nMax) ! diag: sqrt of exchange coefficient (water)
   real(R8),intent(out),optional :: ssq_sv  (nMax) ! diag: sea surface humidity  (kg/kg)

   real(R8),intent(in) ,optional :: missval        ! masked value

! !EOP

   !--- local constants --------------------------------
   real(R8),parameter :: zetam = -1.574_R8 ! Very unstable zeta cutoff for momentum (-)
   real(R8),parameter :: zetat = -0.465_R8 ! Very unstable zeta cutoff for T/q (-)
   real(R8),parameter :: umin  = 0.1_R8    ! minimum wind speed       (m/s)
   real(R8),parameter :: zref  = 10.0_R8   ! reference height           (m)
   real(R8),parameter :: ztref = 2.0_R8    ! reference height for air T (m)
   real(R8),parameter :: beta = 1.0_R8     ! constant used in W* calculation (-)
   real(R8),parameter :: zpbl = 1000.0_R8  ! PBL height used in W* calculation (m)
   real(R8),parameter :: gamma = 0.0098_R8         ! Dry adiabatic lapse rate (K/m)
   real(R8),parameter :: onethird = 1.0_R8/3.0_R8  ! Used repeatedly.

   !--- local variables --------------------------------
   integer(IN) :: n          ! vector loop index
   integer(IN) :: i          ! iteration loop index
   real(R8)    :: vmag_abs   ! surface wind magnitude   (m s-1)
   real(R8)    :: vmag_rel   ! surface wind magnitude relative to
                             ! surface current   (m s-1)
   real(R8)    :: vmag       ! surface wind magnitude with large
                             ! eddy correction and minimum value (m s-1)
                             ! (This can change on each iteration.)
   real(R8)    :: thv        ! virtual temperature      (K)
   real(R8)    :: ssq        ! sea surface humidity     (kg/kg)
   real(R8)    :: delth      ! potential T difference   (K)
   real(R8)    :: delthv     ! virtual potential T difference   (K)
   real(R8)    :: delq       ! humidity difference      (kg/kg)
   real(R8)    :: ustar      ! friction velocity (m s-1)
   real(R8)    :: qstar      ! humidity scaling parameter (kg/kg)
   real(R8)    :: tstar      ! temperature scaling parameter (K)
   real(R8)    :: thvstar    ! virtual temperature scaling parameter (K)
   real(R8)    :: wstar      ! convective velocity scale (m s-1)
   real(R8)    :: zeta       ! dimensionless height (z / Obukhov length)
   real(R8)    :: obu        ! Obukhov length (m)
   real(R8)    :: tau        ! magnitude of wind stress (N m-2)
   real(R8)    :: cp         ! specific heat of moist air (J kg-1 K-1)
   real(R8)    :: xlv        ! Latent heat of vaporization (J kg-1)
   real(R8)    :: visa       ! Kinematic viscosity of dry air (m2 s-1)
   real(R8)    :: tbot_oC    ! Temperature used in visa (deg C)
   real(R8)    :: rb         ! Bulk Richardson number (-)
   real(R8)    :: zo         ! Roughness length for momentum (m)
   real(R8)    :: zoq        ! Roughness length for moisture (m)
   real(R8)    :: zot        ! Roughness length for heat (m)
   real(R8)    :: u10        ! 10-metre wind speed (m s-1)
   real(R8)    :: re         ! Moisture exchange coefficient for compatibility
                             ! with default algorithm.
   real(R8)    :: spval      ! local missing value
   real(R8)    :: loc_epsilon  ! Ratio of gas constants (-)

   !--- for cold air outbreak calc --------------------------------
   real(R8)    :: tdiff(nMax)  ! tbot - ts
   real(R8)    :: vscl

   !--- formats ----------------------------------------
   character(*),parameter :: subName = '(shr_flux_atmOcn) '
   character(*),parameter ::   F00 = "('(shr_flux_atmOcn) ',4a)"

   !-----
   ! Straight from original subroutine.
   if (debug > 0 .and. s_loglev > 0) write(s_logunit,F00) "enter"

   if (present(missval)) then
      spval = missval
   else
      spval = shr_const_spval
   endif
   !-----

   ! Evaluate loc_epsilon.
   loc_epsilon = 1.0_R8 / (1.0_R8 + loc_zvir)

  !--- for cold air outbreak calc --------------------------------
   tdiff = tbot - ts

   ! Loop over grid points.
   DO n=1,nMax
     if (mask(n) /= 0) then

     !-----Calculate some required near surface variables.---------
        vmag_abs = sqrt( ubot(n)**2 + vbot(n)**2 )
        vmag_rel = sqrt( (ubot(n)-us(n))**2 + (vbot(n)-vs(n))**2 )

        ! For Cold Air Outbreak Modification (based on Mahrt & Sun 1995,MWR):
        if (use_coldair_outbreak_mod) then
            ! Increase windspeed for negative tbot-ts
            if (tdiff(n).lt.td0) then
               vscl=min((1._R8+alpha*(abs(tdiff(n)-td0)**0.5_R8/abs(vmag_rel))),maxscl)
               vmag_rel=vmag_rel*vscl
            endif
         endif

        delth = thbot(n) - ts(n)                     ! Pot. temp. difference with surface (K)
                                                     ! Note this is equivalent to Zeng et al
                                                     ! (1998) version = delt + 0.0098*zbot
        thv = thbot(n)*(1.0_R8+0.61_R8*qbot(n))      ! Virtual potential temperature (K)
        ! EQN (17):
        !ssq = 0.98_R8 * qsat_ua(ts(n),ps, &          ! Surface specific humidity (kg kg-1)
        !                        loc_epsilon)
        ssq = 0.98_R8 * qsat_ua(ts(n),pslv(n), &     ! Surface specific humidity (kg kg-1)
                                loc_epsilon)
        delq = qbot(n) - ssq                         ! Difference to surface (kg kg-1)
        delthv = delth*(1.0_R8+0.61_R8*qbot(n)) + &  ! Difference of virtual potential
               & 0.61_R8*thbot(n)*delq               ! temperature with surface (K)

        xlv = 1.0e+6_R8 * &                           ! Latent heat of vaporization (J kg-1)
             & (2.501_R8 - 0.00237_R8 * (ts(n) - loc_tkfrz))
        tbot_oC = tbot(n) - loc_tkfrz
        visa = 1.326e-5_R8 * (1.0_R8 + &             ! Kinematic viscosity of dry
             & 6.542e-3_R8*tbot_oC + &               ! air (m2 s-1) from Andreas (1989)
             & 8.301e-6_R8*tbot_oC*tbot_oC - &       ! CRREL Rep. 89-11
             & 4.84e-9_R8*tbot_oC*tbot_oC*tbot_oC)
        cp = loc_cpdair*(1.0_R8 + loc_cpvir*ssq)     ! specific heat of moist air (J kg-1 K-1)

     !-----Initial values of u* and convective velocity.-----------
        ustar = 0.06_R8
        wstar = 0.5_R8
        ! Update wind speed if unstable regime.
        if (delthv.lt.0.0_R8) then
            ! EQN (19)
            vmag = sqrt( vmag_rel**2 + beta*beta*wstar*wstar )
        else
           ! EQN (18)
            vmag = max(umin,vmag_rel)
        endif

     !-----Iterate to compute new u* and z0.-----------------------
        do i = 1,5
           ! EQN (24)
           zo = 0.013_R8*ustar*ustar/loc_g + 0.11_R8*visa/ustar
           ! EQN (9) assuming neutral
           ustar = loc_karman*vmag/log(zbot(n)/zo)
        enddo

     !-----Assess stability.---------------------------------------
        rb = loc_g*zbot(n)*delthv / (thv*vmag*vmag)    ! bulk Richardson number

        if(rb.ge.0.0_R8) then
            ! Neutral or stable: EQNs (4), (9), (13) and definition of rb.
            zeta = rb*log(zbot(n)/zo) / &
                 & (1.0_R8 - 5.0_R8*min(rb,0.19_R8))
        else
            ! Unstable: EQNs (4), (8), (12) and definition of rb.
            zeta = rb*log(zbot(n)/zo)
        endif

        obu = zbot(n)/zeta                             ! Obukhov length
        obu = sign(max(zbot(n)/10.0_R8, abs(obu)), obu)

     !-----Main iterations (2-10 iterations would be fine).-------
        do i=1,10

            ! Update roughness lengths.
            call rough_ua(zo,zot,zoq,ustar,visa)

            ! Wind variables.
            zeta = zbot(n) / obu
            if (zeta.lt.zetam) then
                ! Very unstable regime
                ! EQN (7) with extra z0 term.
                ustar = loc_karman * vmag / (log(zetam*obu/zo) - &
                     & psi_ua(1_IN, zetam) + &
                     & psi_ua(1_IN, zo/obu) + &
                     & 1.14_R8 * ((-zeta)**onethird - (-zetam)**onethird) )
            else if (zeta.lt.0.0_R8) then
                ! Unstable regime
                ! EQN (8) with extra z0 term.
                ustar = loc_karman * vmag / (log(zbot(n)/zo) - &
                    & psi_ua(1_IN,zeta) + psi_ua(1_IN,zo/obu) )
            else if (zeta.le.1.0_R8) then
                ! Stable regime
                ! EQN (9) with extra z0 term.
                ustar = loc_karman * vmag / (log(zbot(n)/zo) + &
                    & 5.0_R8*zeta - 5.0_R8*zo/obu)
            else
                ! Very stable regime
                ! EQN (10) with extra z0 term.
                ustar = loc_karman * vmag / (log(obu/zo) + 5.0_R8 - &
                     &  5.0_R8*zo/obu + &
                     &  (5.0_R8*log(zeta) + zeta - 1.0_R8) )
            endif

            ! Temperature variables.
            if(zeta.lt.zetat) then
                ! Very unstable regime
                ! EQN (11) with extra z0 term.
                tstar = loc_karman * delth / (log(zetat*obu/zot) - &
                      & psi_ua(2_IN, zetat) + &
                      & psi_ua(2_IN, zot/obu) + &
                      & 0.8_R8*((-zetat)**(-onethird) - (-zeta)**(-onethird)) )
            else if (zeta.lt.0.0_R8) then
                ! Unstable regime
                ! EQN (12) with extra z0 term.
                tstar = loc_karman * delth / &
                      & (log(zbot(n)/zot) - psi_ua(2_IN,zeta) + psi_ua(2_IN,zot/obu))
            else if (zeta.le.1.0_R8) then
                ! Stable regime
                ! EQN (13) with extra z0 term.
                tstar = loc_karman * delth / (log(zbot(n)/zot) + &
                    &   5.0_R8*zeta - 5.0_R8*zot/obu)
            else
                ! Very stable regime
                ! EQN (14) with extra z0 term.
                tstar = loc_karman * delth / (log(obu/zot) + &
                    &   5.0_R8 - 5.0_R8*zot/obu  + &
                    &   (5.0_R8*log(zeta) + zeta - 1.0_R8) )
            endif

            ! Humidity variables.
            ! This is done with re to give variable to save out like
            ! in old algorithm.
            if (zeta.lt.zetat) then
                ! Very unstable regime
                ! EQN (11) with extra z0 term.
                re = loc_karman / (log(zetat*obu/zoq) - psi_ua(2_IN,zetat) + &
                   & psi_ua(2_IN,zoq/obu) + &
                   & 0.8_R8*((-zetat)**(-onethird) - (-zeta)**(-onethird)) )
            else if (zeta.lt.0.0_R8) then
                ! Unstable regime
                ! EQN (12) with extra z0 term.
                re = loc_karman / &
                   & (log(zbot(n)/zoq) - psi_ua(2_IN,zeta) + psi_ua(2_IN,zoq/obu))
            else if (zeta.le.1.0_R8) then
                ! Stable regime
                ! EQN (13) with extra z0 term.
                re = loc_karman / &
                   & (log(zbot(n)/zoq) + 5.0_R8*zeta - 5.0_R8*zoq/obu)
            else
                ! Very stable regime
                ! EQN (14) with extra z0 term.
                re = loc_karman / &
                   & (log(obu/zoq) + 5.0_R8 - 5.0_R8*zoq/obu + &
                   & (5.0_R8*log(zeta) + zeta - 1.0_R8) )
             endif
             qstar = re * delq

            ! Update Obukhov length.
            thvstar = tstar*(1.0_R8 + 0.61_R8*qbot(n)) + 0.61_R8*thbot(n)*qstar
            ! EQN (4)
            obu = ustar*ustar * thv / (loc_karman*loc_g*thvstar)
            obu = sign( max(zbot(n)/10.0_R8, abs(obu)) ,obu)

            ! Update wind speed if in unstable regime.
            if (delthv.lt.0.0_R8) then
                ! EQN (20)
                wstar = beta * (-loc_g*ustar*thvstar*zpbl/thv)**onethird
                ! EQN (19)
                vmag = sqrt(vmag_rel**2 + wstar*wstar)
             else
                ! EQN (18)
                vmag = max(umin,vmag_rel)
            endif

        enddo ! End of iterations for ustar, tstar, qstar etc.


     !-----Calculate fluxes and wind stress.---------------------

        !--- momentum flux ---
        ! This should ensure zero wind stress when (relative) wind speed is zero,
        ! components are consistent with total, and we don't ever divide by zero.
        ! EQN (21)
        tau = rbot(n) * ustar * ustar
        taux(n) = tau * (ubot(n)-us(n)) / max(umin, vmag_rel)
        tauy(n) = tau * (vbot(n)-vs(n)) / max(umin, vmag_rel)

        !--- heat flux ---
        ! EQNs (22) and (23)
        sen (n) =  cp * rbot(n) * tstar * ustar
        lat (n) = xlv * rbot(n) * qstar * ustar
        lwup(n) = -loc_stebol * ts(n)**4

        !--- water flux ---
        evap(n) = lat(n)/xlv

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

        zeta = zbot(n) / obu
        if (zeta.lt.zetat) then
            if (zeta.lt.zetam) then
               ! Very unstable regime for U.
               ! EQN (7)
               u10 = vmag_abs + (ustar/loc_karman) * &
                   & 1.14_R8 * ((-zref/obu)**onethird - (-zeta)**onethird)
            else
               ! Unstable regime for U.
               ! EQN (8)
               u10 = vmag_abs + (ustar/loc_karman) * &
                   & (log(zref/zbot(n)) - (psi_ua(1_IN,zref/obu) - psi_ua(1_IN,zeta)) )
            endif
            ! Very unstable regime for T and q.
            ! EQN (11)
            tref(n) = thbot(n) + (tstar/loc_karman) * &
                    & 0.8_R8 * ((-zeta)**(-onethird) - (-ztref/obu)**(-onethird))
            qref(n) = qbot(n) + (qstar/loc_karman) * &
                    & 0.8_R8 * ((-zeta)**(-onethird) - (-ztref/obu)**(-onethird))

        else if (zeta.lt.0.0_R8) then
            ! Unstable regime.
            ! EQN (8)
            u10 = vmag_abs + (ustar/loc_karman) * &
                & (log(zref/zbot(n)) - (psi_ua(1_IN,zref/obu) - psi_ua(1_IN,zeta)) )
            ! EQN (12)
            tref(n) = thbot(n) + (tstar/loc_karman) * &
                    & (log(ztref/zbot(n)) - (psi_ua(2_IN,ztref/obu) - psi_ua(2_IN,zeta)) )
            qref(n) = qbot(n) + (qstar/loc_karman) * &
                    & (log(ztref/zbot(n)) - (psi_ua(2_IN,ztref/obu) - psi_ua(2_IN,zeta)) )
        else if (zeta.le.1.0_R8) then
            ! Stable regime.
            ! EQN (9)
            u10 = vmag_abs + (ustar/loc_karman) * &
                & (log(zref/zbot(n)) + 5.0_R8*zref/obu - 5.0_R8*zeta)
            ! EQN (13)
            tref(n) = thbot(n) + (tstar/loc_karman) * &
                    & (log(ztref/zbot(n)) + 5.0_R8*ztref/obu - 5.0_R8*zeta)
            qref(n) = qbot(n) + (qstar/loc_karman) * &
                 & (log(ztref/zbot(n)) + 5.0_R8*ztref/obu - 5.0_R8*zeta)
         else
            ! Very stable regime.
            ! EQN (10)
            u10 = vmag_abs + (ustar/loc_karman) * &
                & (5.0_R8*log(zref/zbot(n)) + zref/obu - zeta)
            ! EQN (14)
            tref(n) = thbot(n) + (tstar/loc_karman) * &
                    & (5.0_R8*log(ztref/zbot(n)) + ztref/obu - zeta)
            qref(n) = qbot(n) + (qstar/loc_karman) * &
                    & (5.0_R8*log(ztref/zbot(n)) + ztref/obu - zeta)

        endif

        tref(n) = tref(n) - gamma*ztref   ! pot. temp to temp correction
        duu10n(n) = u10*u10 ! 10m wind speed squared

        !------------------------------------------------------------
        ! optional diagnostics, needed for water tracer fluxes (dcn)
        !------------------------------------------------------------
        if (present(ustar_sv)) ustar_sv(n) = ustar
        if (present(ssq_sv  )) ssq_sv(n)   = ssq
        if (present(re_sv   )) re_sv(n)    = re


     else

        !------------------------------------------------------------
        ! no valid data here -- out of ocean domain
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
        ! Optional diagnostics too:
        if (present(ustar_sv)) ustar_sv(n) = spval
        if (present(re_sv   )) re_sv   (n) = spval
        if (present(ssq_sv  )) ssq_sv  (n) = spval

     endif

   ENDDO ! loop over grid points

END subroutine shr_flux_atmOcn_UA

!===============================================================================
! Functions/subroutines used by UA surface flux scheme.
!===============================================================================

   ! Stability function for rb < 0

real(R8) function psi_ua(k,zeta)

       implicit none

       !-----Input variables.----------
       integer(IN), intent(in) :: k       ! Indicates whether this is for momentum (k=1)
                                          ! or for heat/moisture (k=2)
       real(R8), intent(in) :: zeta       ! Dimensionless height (=z/L)

       !-----Local variables.----------
       real(R8) :: chik                   ! Function of zeta.

       ! EQN (16)
       chik = (1.0_R8 - 16.0_R8*zeta)**0.25_R8

       if(k.eq.1) then
          ! EQN (15) for momentum
          psi_ua = 2.0_R8 * log((1.0_R8 + chik)*0.5_R8) + &
                 &      log((1.0_R8 + chik*chik)*0.5_R8) - &
                 & 2.0_R8 * atan(chik) + 2.0_R8 * atan(1.0_R8)
       else
          ! EQN (15) for heat/moisture
          psi_ua = 2.0_R8 * log((1.0_R8 + chik*chik)*0.5_R8)
       endif

end function psi_ua

!===============================================================================
   ! Uses Tetens' formula for saturation vapor pressure from
   ! Buck(1981) JAM 20, 1527-1532

real(R8) function qsat_ua(t,p,loc_epsilon)

       implicit none

       !-----Input variables.----------
       real(R8), intent(in) :: t           ! temperature (K)
       real(R8), intent(in) :: p           ! pressure (Pa)
       real(R8), intent(in) :: loc_epsilon ! Ratio of gas constants (-)

       !-----Local variables.----------
       real(R8) :: esat                    ! saturated vapor pressure (hPa)

       ! Calculate saturated vapor pressure in hPa.
       esat = (1.0007_R8 + 0.00000346_R8 * (p/100.0_R8)) * 6.1121_R8 * &
            & exp(17.502_R8 * (t - loc_tkfrz) / (240.97_R8 + (t - loc_tkfrz)))

       ! Convert to specific humidity (kg kg-1).
       qsat_ua = loc_epsilon * esat / ((p/100.0_R8) - (1.0_R8 - loc_epsilon)*esat)

end function qsat_ua

!===============================================================================
   !Calculate roughness lengths: zo, zot, zoq.

subroutine rough_ua(zo,zot,zoq,ustar,visa)

       implicit none

       !-----Input variables.----------
       real(R8), intent(in) :: ustar      ! friction velocity (m s-1)
       real(R8), intent(in) :: visa       ! kinematic viscosity of dry air (m2 s-1)

       !-----Output variables.---------
       real(R8), intent(out) :: zo        ! roughness length for momentum (m)
       real(R8), intent(out) :: zot       ! roughness length for heat (m)
       real(R8), intent(out) :: zoq       ! roughness length for water vapor (m)

       !-----Local variables.----------
       real(R8) :: re_rough               ! Rougness Reynold's number (-)
       real(R8) :: xq                     ! Logarithm of roughness length ratios (moisture)
       real(R8) :: xt                     ! Logarithm of roughness length ratios (heat)

       zo = 0.013_R8*ustar*ustar/loc_g + 0.11_R8*visa/ustar      ! EQN (24)
       re_rough = ustar*zo/visa                                  ! By definition.
       xq = 2.67_R8*re_rough**0.25_R8 - 2.57_R8                  ! EQN (25)
       xt = xq                                                   ! EQN (26)
       zoq = zo/exp(xq)                                          ! By definition of xq
       zot = zo/exp(xt)                                          ! By definition of xt

end subroutine rough_ua



real(R8) elemental function cuberoot(a)
  real(R8), intent(in) :: a
  real(R8), parameter :: one_third = 1._R8/3._R8
  cuberoot = sign(abs(a)**one_third, a)
end function cuberoot

!===============================================================================
! !BOP =========================================================================
!
! !IROUTINE: shr_flux_atmOcn_diurnal -- internal atm/ocn flux calculation
!
! !DESCRIPTION:
!
!     Internal atm/ocn flux calculation
!
! !REVISION HISTORY:
!     2002-Jun-10 - B. Kauffman - code migrated from cpl5 to cpl6
!     2003-Apr-02 - B. Kauffman - taux & tauy now utilize ocn velocity
!     2003-Apr-02 - B. Kauffman - tref,qref,duu10n mods as per Bill Large
!     2006-Nov-07 - B. Kauffman - code migrated from cpl6 to share
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE shr_flux_atmOcn_diurnal &
                          (nMax  ,zbot  ,ubot  ,vbot  ,thbot ,             &
                           qbot  ,s16O  ,sHDO  ,s18O  ,rbot  ,             &
                           tbot  ,us    ,vs    ,                           &
                           ts    ,mask  , seq_flux_atmocn_minwind,         &
                           sen   ,lat   ,lwup  ,                           &
                           r16O  ,rhdo  ,r18O  ,evap  ,evap_16O,           &
                           evap_HDO     ,evap_18O,                         &
                           taux  ,tauy  ,tref  ,qref  ,                    &
                           uGust, lwdn , swdn , swup, prec   ,             &
                           swpen, ocnsal, ocn_prognostic, flux_diurnal,    &
                           ocn_surface_flux_scheme,                                    &
                           latt, long , warm , salt , speed, regime,       &
                           warmMax, windMax, qSolAvg, windAvg,             &
                           warmMaxInc, windMaxInc, qSolInc, windInc, nInc, &
                           tBulk, tSkin, tSkin_day, tSkin_night,           &
                           cSkin, cSkin_night, secs ,dt,                   &
                           duu10n,  ustar_sv   ,re_sv ,ssq_sv,             &
                           missval, cold_start    )
! !USES:

   use water_isotopes, only: wiso_flxoce !subroutine used to calculate water isotope fluxes.

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   !--- input arguments --------------------------------
   integer(IN),intent(in) ::       nMax  ! data vector length
   integer(IN),intent(in) :: mask (nMax) ! ocn domain mask       0 <=> out of domain
   real(R8)   ,intent(in) :: zbot (nMax) ! atm level height      (m)
   real(R8)   ,intent(in) :: ubot (nMax) ! atm u wind            (m/s)
   real(R8)   ,intent(in) :: vbot (nMax) ! atm v wind            (m/s)
   real(R8)   ,intent(in) :: thbot(nMax) ! atm potential T       (K)
   real(R8)   ,intent(in) :: qbot (nMax) ! atm specific humidity (kg/kg)
   real(R8)   ,intent(in) :: s16O (nMax) ! atm H216O tracer conc. (kg/kg)
   real(R8)   ,intent(in) :: sHDO (nMax) ! atm HDO tracer conc.  (kg/kg)
   real(R8)   ,intent(in) :: s18O (nMax) ! atm H218O tracer conc. (kg/kg)
   real(R8)   ,intent(in) :: r16O (nMax) ! ocn H216O tracer ratio/Rstd
   real(R8)   ,intent(in) :: rHDO (nMax) ! ocn HDO tracer ratio/Rstd
   real(R8)   ,intent(in) :: r18O (nMax) ! ocn H218O tracer ratio/Rstd
   real(R8)   ,intent(in) :: rbot (nMax) ! atm air density       (kg/m^3)
   real(R8)   ,intent(in) :: tbot (nMax) ! atm T                 (K)
   real(R8)   ,intent(in) :: us   (nMax) ! ocn u-velocity        (m/s)
   real(R8)   ,intent(in) :: vs   (nMax) ! ocn v-velocity        (m/s)
   real(R8)   ,intent(in) :: ts   (nMax) ! ocn temperature       (K)

   !--- new    arguments -------------------------------
   real(R8),intent(inout) :: swpen (nMax)       ! NEW
   real(R8),intent(inout) :: ocnsal(nMax)       ! NEW (kg/kg)
   logical ,intent(in)    :: ocn_prognostic     ! NEW
   logical ,intent(in)    :: flux_diurnal       ! NEW logical for diurnal on/off
   integer(IN) ,intent(in)    :: ocn_surface_flux_scheme

   real(R8),intent(in)    :: uGust (nMax)      ! NEW not used
   real(R8),intent(in)    :: lwdn  (nMax)       ! NEW
   real(R8),intent(in)    :: swdn  (nMax)       ! NEW
   real(R8),intent(in)    :: swup  (nMax)       ! NEW
   real(R8),intent(in)    :: prec  (nMax)       ! NEW
   real(R8),intent(in)    :: latt  (nMax)       ! NEW
   real(R8),intent(in)    :: long  (nMax)       ! NEW
   real(R8),intent(inout) :: warm  (nMax)       ! NEW
   real(R8),intent(inout) :: salt  (nMax)       ! NEW
   real(R8),intent(inout) :: speed (nMax)       ! NEW
   real(R8),intent(inout) :: regime(nMax)       ! NEW
   real(R8),intent(out)   :: warmMax(nMax)      ! NEW
   real(R8),intent(out)   :: windMax(nMax)      ! NEW
   real(R8),intent(inout) :: qSolAvg(nMax)      ! NEW
   real(R8),intent(inout) :: windAvg(nMax)      ! NEW
   real(R8),intent(inout) :: warmMaxInc(nMax)   ! NEW
   real(R8),intent(inout) :: windMaxInc(nMax)   ! NEW
   real(R8),intent(inout) :: qSolInc(nMax)      ! NEW
   real(R8),intent(inout) :: windInc(nMax)      ! NEW
   real(R8),intent(inout) :: nInc(nMax)         ! NEW

   real(R8),intent(out)   :: tBulk (nMax)       ! NEW
   real(R8),intent(out)   :: tSkin (nMax)       ! NEW
   real(R8),intent(out)   :: tSkin_day (nMax)   ! NEW
   real(R8),intent(out)   :: tSkin_night (nMax) ! NEW
   real(R8),intent(out)   :: cSkin (nMax)       ! NEW
   real(R8),intent(out)   :: cSkin_night (nMax) ! NEW
   integer(IN),intent(in) :: secs               ! NEW  elsapsed seconds in day (GMT)
   integer(IN),intent(in) :: dt                 ! NEW
   logical ,intent(in)    :: cold_start         ! cold start flag
   real(R8),intent(in)    :: seq_flux_atmocn_minwind   ! minimum wind speed for atmocn      (m/s)

   real(R8),intent(in) ,optional :: missval     ! masked value

   !--- output arguments -------------------------------
   real(R8),intent(out)  ::  sen  (nMax) ! heat flux: sensible    (W/m^2)
   real(R8),intent(out)  ::  lat  (nMax) ! heat flux: latent      (W/m^2)
   real(R8),intent(out)  ::  lwup (nMax) ! heat flux: lw upward   (W/m^2)
   real(R8),intent(out)  ::  evap (nMax) ! water flux: evap  ((kg/s)/m^2)
   real(R8),intent(out)  ::  evap_16O (nMax) ! water flux: evap ((kg/s/m^2)
   real(R8),intent(out)  ::  evap_HDO (nMax) ! water flux: evap ((kg/s)/m^2)
   real(R8),intent(out)  ::  evap_18O (nMax) ! water flux: evap ((kg/s/m^2)
   real(R8),intent(out)  ::  taux (nMax) ! surface stress, zonal      (N)
   real(R8),intent(out)  ::  tauy (nMax) ! surface stress, maridional (N)
   real(R8),intent(out)  ::  tref (nMax) ! diag:  2m ref height T     (K)
   real(R8),intent(out)  ::  qref (nMax) ! diag:  2m ref humidity (kg/kg)
   real(R8),intent(out)  :: duu10n(nMax) ! diag: 10m wind speed squared (m/s)^2

   real(R8),intent(out),optional :: ustar_sv(nMax) ! diag: ustar
   real(R8),intent(out),optional :: re_sv   (nMax) ! diag: sqrt of exchange coefficient (water)
   real(R8),intent(out),optional :: ssq_sv  (nMax) ! diag: sea surface humidity  (kg/kg)

! !EOP


   !--- local constants --------------------------------
   real(R8),parameter :: zref  = 10.0_R8 ! reference height           (m)
   real(R8),parameter :: ztref =  2.0_R8 ! reference height for air T (m)

   real(R8),parameter :: lambdaC  = 6.0_R8
   real(R8),parameter :: lambdaL  = 0.0_R8
   real(R8),parameter :: doLMax   = 1.0_R8
   real(R8),parameter :: pwr      = 0.2_R8
   real(R8),parameter :: Rizero   = 1.0_R8
   real(R8),parameter :: NUzero   = 40.0e-4_R8
   real(R8),parameter :: Prandtl  = 1.0_R8
   real(R8),parameter :: kappa0   = 0.2e-4_R8

   real(R8),parameter :: F0       = 0.5_R8
   real(R8),parameter :: F1       = 0.15_R8
   real(R8),parameter :: R1       = 10.0_R8

   real(R8),parameter :: Ricr     = 0.30_R8
   real(R8),parameter :: tiny     = 1.0e-12_R8
   real(R8),parameter :: tiny2    = 1.0e-6_R8
   real(R8),parameter :: pi       = SHR_CONST_PI

!!++ COARE only
   real(R8),parameter :: zpbl =700.0_R8 ! PBL depth [m] for gustiness parametriz.

   !--- local variables --------------------------------
   integer(IN) :: n       ! vector loop index
   integer(IN) :: iter       ! iteration loop index
   integer(IN) :: lsecs   ! local seconds elapsed
   integer(IN) :: lonsecs ! incrememnt due to lon offset
   real(R8)    :: vmag    ! surface wind magnitude   (m/s)
   real(R8)    :: ssq     ! sea surface humidity     (kg/kg)
   real(R8)    :: delt    ! potential T difference   (K)
   real(R8)    :: delq    ! humidity difference      (kg/kg)
   real(R8)    :: stable  ! stability factor
   real(R8)    :: rdn     ! sqrt of neutral exchange coeff (momentum)
   real(R8)    :: rhn     ! sqrt of neutral exchange coeff (heat)
   real(R8)    :: ren     ! sqrt of neutral exchange coeff (water)
   real(R8)    :: rd      ! sqrt of exchange coefficient (momentum)
   real(R8)    :: rh      ! sqrt of exchange coefficient (heat)
   real(R8)    :: re      ! sqrt of exchange coefficient (water)
   real(R8)    :: ustar   ! ustar
   real(R8)    :: ustar_prev   ! ustar
   real(R8)    :: qstar   ! qstar
   real(R8)    :: tstar   ! tstar
   real(R8)    :: hol     ! H (at zbot) over L
   real(R8)    :: xsq     ! ?
   real(R8)    :: xqq     ! ?
   real(R8)    :: psimh   ! stability function at zbot (momentum)
   real(R8)    :: psixh   ! stability function at zbot (heat and water)
   real(R8)    :: psix2   ! stability function at ztref reference height
   real(R8)    :: alz     ! ln(zbot/zref)
   real(R8)    :: al2     ! ln(zref/ztref)
   real(R8)    :: u10n    ! 10m neutral wind
   real(R8)    :: tau     ! stress at zbot
   real(R8)    :: cp      ! specific heat of moist air
   real(R8)    :: fac     ! vertical interpolation factor
   real(R8)    :: DTiter  !
   real(R8)    :: DSiter  !
   real(R8)    :: DViter  !

   real(R8)    :: Dcool   !
   real(R8)    :: Qdel    ! net cool skin heating
   real(R8)    :: Hd      ! net heating above -z=d
   real(R8)    :: Hb      ! net kinematic heating above -z = delta
   real(R8)    :: lambdaV !
   real(R8)    :: Fd      ! net fresh water forcing above -z=d
   real(R8)    :: ustarw  ! surface wind forcing of layer above -z=d

   real(R8)    :: Qsol   ! solar heat flux (W/m2)
   real(R8)    :: Qnsol  ! non-solar heat flux (W/m2)

   real(R8)    :: SSS  ! sea surface salinity
   real(R8)    :: alphaT  !
   real(R8)    :: betaS  !

   real(R8)    :: doL     ! ocean forcing stablity parameter
   real(R8)    :: Rid     ! Richardson number at depth d
   real(R8)    :: Ribulk  ! Bulk  Richardson number at depth d
   real(R8)    :: FofRi   ! Richardon number dependent diffusivity
   real(R8)    :: Smult   ! multiplicative term based on regime
   real(R8)    :: Sfact   ! multiplicative term based on regime
   real(R8)    :: Kdiff   ! diffusive term based on regime
   real(R8)    :: Kvisc   ! viscosity term based on regime
   real(R8)    :: rhocn   !
   real(R8)    :: rcpocn  !
   real(R8)    :: Nreset  ! value for multiplicative reset factor
   logical     :: lmidnight
   logical     :: ltwopm
   logical     :: ltwoam
   logical     :: lfullday
   integer     :: nsum
   real(R8)    :: pexp   ! eqn 19
   real(R8)    :: AMP    ! eqn 18
   real(R8)    :: dif3
   real(R8)    :: phid
   real(R8)    :: spval

!!++ COARE only
   real(R8)    :: zo,zot,zoq      ! roughness lengths
   real(R8)    :: hsb,hlb         ! sens & lat heat flxs at zbot
   real(R8)    :: trf,qrf,urf,vrf ! reference-height quantities

   !--- local functions --------------------------------
   real(R8)    :: qsat   ! function: the saturation humididty of air (kg/m^3)
   real(R8)    :: cdn    ! function: neutral drag coeff at 10m
   real(R8)    :: psimhu ! function: unstable part of psimh
   real(R8)    :: psixhu ! function: unstable part of psimx
   real(R8)    :: Umps   ! dummy arg ~ wind velocity (m/s)
   real(R8)    :: Tk     ! dummy arg ~ temperature (K)
   real(R8)    :: xd     ! dummy arg ~ ?
   real(R8)    :: molvisc ! molecular viscosity
   real(R8)    :: molPr   ! molecular Prandtl number

   !--- for cold air outbreak calc --------------------------------
   real(R8)    :: tdiff(nMax)               ! tbot - ts
   real(R8)    :: vscl

   qsat(Tk)   = 640380.0_R8 / exp(5107.4_R8/Tk)
   cdn(Umps)  =   0.0027_R8 / Umps + 0.000142_R8 + 0.0000764_R8 * Umps
   psimhu(xd) = log((1.0_R8+xd*(2.0_R8+xd))*(1.0_R8+xd*xd)/8.0_R8) - 2.0_R8*atan(xd) + 1.571_R8
   psixhu(xd) = 2.0_R8 * log((1.0_R8 + xd*xd)/2.0_R8)
   molvisc(Tk)  = 1.623e-6_R8 * exp((-1.0_R8*(Tk-273.15_R8))/45.2_R8)
   molPr(Tk)    = 11.64_R8 * exp((-1.0_R8*(Tk-273.15_R8))/40.7_R8)

   !--- formats ----------------------------------------
   character(*),parameter :: subName = '(shr_flux_atmOcn_diurnal) '
   character(*),parameter ::   F00 = "('(shr_flux_atmOcn_diurnal) ',4a)"

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

   if (debug > 0 .and. s_loglev > 0) write(s_logunit,F00) "enter"

   ! this is especially for flux_diurnal calculations
   if (.not. flux_diurnal) then
      write(s_logunit,F00) "ERROR: flux_diurnal must be true"
      call shr_sys_abort(subName//"flux diurnal must be true")
   endif
   spval = shr_const_spval
   rh = spval
   dviter = spval
   dtiter = spval
   dsiter = spval
   al2 = log(zref/ztref)
  !--- for cold air outbreak calc --------------------------------
   tdiff= tbot - ts

   ! equations 18 and 19
   AMP = 1.0_R8/F0-1.0_R8
   pexp = log( (1.0_R8/F1-F0) / (1.0_R8-F0) ) / log(R1)

   if (.not. ocn_prognostic) then
      ! Set swpen and ocean salinity from following analytic expressions
      swpen(:) = 0.67_R8*(exp((-1._R8*shr_const_zsrflyr)/1.0_R8)) + &
           0.33_R8*exp((-1._R8*shr_const_zsrflyr)/17.0_R8)
      ocnsal(:) = shr_const_ocn_ref_sal/1000.0_R8
   else
      ! use swpen and ocnsal from input argument
   endif

   if (cold_start) then
      !         if (s_loglev > 0) then
      write(s_logunit,F00) "Initialize diurnal cycle fields"
      !         end if
      warm       (:) = 0.0_R8
      salt       (:) = 0.0_R8
      speed      (:) = 0.0_R8
      regime     (:) = 0.0_R8
      qSolAvg    (:) = 0.0_R8
      windAvg    (:) = 0.0_R8
      warmMax    (:) = 0.0_R8
      windMax    (:) = 0.0_R8
      warmMaxInc (:) = 0.0_R8
      windMaxInc (:) = 0.0_R8
      qSolInc    (:) = 0.0_R8
      windInc    (:) = 0.0_R8
      nInc       (:) = 0.0_R8
      tSkin_day  (:) = ts(:)
      tSkin_night(:) = ts(:)
      cSkin_night(:) = 0.0_R8
   endif

   DO n=1,nMax

      if (mask(n) /= 0) then

         !--- compute some initial and useful flux quantities ---

         vmag     = max(seq_flux_atmocn_minwind, sqrt( (ubot(n)-us(n))**2 + (vbot(n)-vs(n))**2) )
         if (use_coldair_outbreak_mod) then
            ! Cold Air Outbreak Modification:
            ! Increase windspeed for negative tbot-ts
            ! based on Mahrt & Sun 1995,MWR

            if (tdiff(n).lt.td0) then
               vscl=min((1._R8+alpha*(abs(tdiff(n)-td0)**0.5_R8/abs(vmag))),maxscl)
               vmag=vmag*vscl
            endif
         endif
         alz      = log(zbot(n)/zref)
         hol      = 0.0
         psimh    = 0.0
         psixh    = 0.0
         rdn      = sqrt(cdn(vmag))

         tBulk(n) = ts(n)+warm(n)    ! first guess for tBulk from read in ts,warm
         tSkin(n) = tBulk(n)
         Qsol     = swdn(n) + swup(n)
         SSS      = 1000.0_R8*ocnsal(n)+salt(n)
         lambdaV = lambdaC

         alphaT   = 0.000297_R8*(1.0_R8+0.0256_R8*(ts(n)-298.15_R8)+0.003_R8*(SSS - 35.0_R8))
         betaS    = 0.000756_R8*(1.0_R8-0.0016_R8*(ts(n)-298.15_R8))
         rhocn    = 1023.342_R8*(1.0_R8-0.000297_R8*(ts(n)-298.15_R8)+0.000756_R8 * (SSS - 35.0_R8))
         rcpocn   = rhocn * 3990.0_R8*(1.0_R8-0.0012_R8*(SSS - 35.0_R8))

         Rid =  shr_const_g * (alphaT*warm(n) - betaS*salt(n)) *pwr*shr_const_zsrflyr  / &
              ( pwr*MAX(tiny,speed(n)) )**2

         Ribulk = 0.0

         !----------------------------------------------------------
         ! convert elapsed time from GMT to local &
         ! check elapsed time. reset warm if near lsecs = reset_sec
         !----------------------------------------------------------
         Nreset = 1.0_R8

         lonsecs   = ceiling(long(n)/360.0_R8*86400.0)
         lsecs     = mod(secs + lonsecs,86400)

         lmidnight = (lsecs >= 0     .and. lsecs < dt)        ! 0 = midnight
         ltwopm    = (lsecs >= 48600 .and. lsecs < 48600+dt)  ! 48600 = 1:30pm
         ltwoam    = (lsecs >= 5400  .and. lsecs < 5400 +dt)  ! 5400 = 1:30am
         lfullday  = (lsecs > 86400-dt .and. lsecs <= 86400)
         nsum = nint(nInc(n))

         if ( lmidnight ) then
            Regime(n)  = 1.0_R8               !  RESET DIURNAL
            warm(n)    = 0.0_R8
            salt(n)    = 0.0_R8
            speed(n)   = 0.0_R8
         endif

         ssq    = 0.98_R8 * qsat(tBulk(n)) / rbot(n)   ! sea surf hum (kg/kg)
         delt   = thbot(n) - tBulk(n)                  ! pot temp diff (K)
         delq   = qbot(n) - ssq                     ! spec hum dif (kg/kg)
         cp     = shr_const_cpdair*(1.0_R8 + shr_const_cpvir*ssq)

!!.................................................................
!! ocn_surface_flux_scheme = 0 : Default E3SMv1
!!                         = 1 : COARE algorithm
!!.................................................................
         if (ocn_surface_flux_scheme .eq. 0) then! use Large algorithm
            stable = 0.5_R8 + sign(0.5_R8 , delt)


            !--- shift wind speed using old coefficient  and stability function

            rd   = rdn / (1.0_R8 + rdn/shr_const_karman*(alz-psimh))
            u10n = vmag * rd / rdn

            !--- initial neutral  transfer coeffs at 10m
            rdn    = sqrt(cdn(u10n))
            rhn    = (1.0_R8-stable) * 0.0327_R8 + stable * 0.018_R8
            ren    = 0.0346_R8

            !--- initial ustar, tstar, qstar ---
            ustar = rdn * vmag
            tstar = rhn * delt
            qstar = ren * delq

         else if (ocn_surface_flux_scheme .eq. 1) then! use COARE algorithm

              call cor30a(ubot(n),vbot(n),tbot(n),qbot(n),rbot(n) &  ! in atm params
                       & ,us(n),vs(n),tBulk(n),ssq                &  ! in surf params (NB ts -> tBulk)
                       & ,zpbl,zbot(n),zbot(n),zref,ztref,ztref   &  ! in heights
                       & ,tau,hsb,hlb                             &  ! out: fluxes
                       & ,zo,zot,zoq,hol,ustar,tstar,qstar        &  ! out: ss scales
                       & ,rd,rh,re                                &  ! out: exch. coeffs
                       & ,trf,qrf,urf,vrf)			       ! out: reference-height params
             ! for the sake of maintaining same defs
              hol=zbot(n)/hol
              rd=sqrt(rd)
              rh=sqrt(rh)
              re=sqrt(re)

         ELSE  ! N.B.: *no* valid ocn_surface_flux_scheme=2 option if diurnal=.true.

            call shr_sys_abort(subName//" shr_flux_atmOcn_diurnal requires ocn_surface_flux_scheme = 0 or 1")
         ENDIF


        ustar_prev = ustar * 2.0_R8
        iter = 0
        ! --- iterate ---
        ! Originally this code did three iterations while the non-diurnal version did two
        ! So in the new loop this is <= flux_con_max_iter instead of < so that the same defaults
        ! will give the same answers in both cases.
        do while( abs((ustar - ustar_prev)/ustar) > flux_con_tol .and. iter <= flux_con_max_iter)
            iter = iter + 1
            ustar_prev = ustar
            !------------------------------------------------------------
            ! iterate to converge on FLUXES  Z/L, ustar, tstar and qstar
            ! and on Rid  in the DIURNAL CYCLE
            !------------------------------------------------------------
            Smult = 0.0_R8
            Sfact = 0.0_R8
            Kdiff = 0.0_R8
            Kvisc = 0.0_R8
            dif3 = 0.0_R8

            ustarw  = ustar*sqrt(max(tiny,rbot(n)/rhocn))
            Qnsol   = lwdn(n) - shr_const_stebol*(tSkin(n))**4 + &
                 rbot(n)*ustar*(cp*tstar + shr_const_latvap*qstar)
            Hd      = (Qnsol   + Qsol*(1.0_R8-swpen(n)) ) / rcpocn
            Fd      = (prec(n) + rbot(n)*ustar*qstar ) * SSS / rhocn

            !--- COOL SKIN EFFECT ---
            Dcool  = lambdaV*molvisc(tBulk(n)) / ustarw
            Qdel   = Qnsol + Qsol * &
                 (0.137_R8 + 11.0_R8*Dcool - 6.6e-5/Dcool *(1.0_R8 - exp((-1.0_R8*Dcool)/8.0e-4)))
            Hb = (Qdel/rcpocn)+(Fd*betaS/alphaT)
            Hb = min(Hb , 0.0_R8)

!            lambdaV = lambdaC*(1.0_R8 + ( (0.0_R8-Hb)*16.0_R8*molvisc(tBulk(n))* &
!                 shr_const_g*alphaT*molPr(tBulk(n))**2/ustarw**4)**0.75)**(-1._R8/3._R8)
            lambdaV = 6.5_R8
            cSkin(n) =  MIN(0.0_R8, lambdaV * molPr(tBulk(n)) * Qdel / ustarw / rcpocn )

            !--- REGIME ---
            doL = shr_const_zsrflyr*shr_const_karman*shr_const_g* &
                 (alphaT*Hd + betaS*Fd ) / ustarw**3
            Rid = MAX(0.0_R8,Rid)
            Smult = dt * (pwr+1.0_R8) / (shr_const_zsrflyr*pwr)
            Sfact = dt * (pwr+1.0_R8) / (shr_const_zsrflyr)**2
            FofRi = 1.0_R8/(1.0_R8 + AMP*(Rid/Rizero)**pexp)

            if ( (doL.gt.0.0_R8) .and. (Qsol.gt.0.0)  ) then
               phid  = MIN(1.0_R8 + 5.0_R8 * doL, 5.0_R8 + doL)
               FofRi = 1.0_R8/(1.0_R8 + AMP*(Rid/Rizero)**pexp)
               dif3 = (kappa0 + NUzero *FofRi)

               if ((doL.le.lambdaL).and.(NINT(regime(n)).le.2)) then
                  regime(n) = 2.0_R8
                  Kdiff =  shr_const_karman * ustarw * shr_const_zsrflyr / phid
                  Kvisc = Kdiff * (1.0_R8 - doL/lambdaL)**2 + &
                       dif3 * (doL/lambdaL)**2 * (3.0_R8 - 2.0_R8 * doL/lambdaL)
                  Kdiff = Kvisc
               else
                  regime(n) = 3.0_R8
                  Kdiff =          kappa0 + NUzero * FofRi
                  Kvisc = Prandtl* kappa0 + NUzero * FofRi
               endif
            else
               if (regime(n).eq.1.0_R8) then
                  Smult      = 0.0_R8
               else
                  if (Ribulk .gt. Ricr) then
                     regime(n) = 3.0_R8
                     Kdiff =          kappa0 + NUzero * FofRi
                     Kvisc = Prandtl* kappa0 + NUzero * FofRi
                  else
                     regime(n) = 4.0_R8
                     Kdiff = shr_const_karman*ustarw*shr_const_zsrflyr *cuberoot(1.0_R8-7.0_R8*doL)
                     Kvisc = Kdiff
                  endif
               endif

            endif

            !--- IMPLICIT INTEGRATION ---

            DTiter = (warm(n)  +(Smult*Hd))               /(1.+ Sfact*Kdiff)
            DSiter = (salt(n)  -(Smult*Fd))               /(1.+ Sfact*Kdiff)
            DViter = (speed(n) +(Smult*ustarw*ustarw))    /(1.+ Sfact*Kvisc)
            DTiter = MAX( 0.0_R8, DTiter)
            DViter = MAX( 0.0_R8, DViter)

            Rid =(shr_const_g*(alphaT*DTiter-betaS*DSiter)*pwr*shr_const_zsrflyr)  / &
                 (pwr*MAX(tiny,DViter))**2
            Ribulk = Rid * pwr
            Ribulk = 0.0_R8
            tBulk(n) = ts(n) + DTiter
            tSkin(n) = tBulk(n) + cskin(n)

            !--need to update ssq,delt,delq as function of tBulk ----

            ssq    = 0.98_R8 * qsat(tBulk(n)) / rbot(n)   ! sea surf hum (kg/kg)
            delt   = thbot(n) - tBulk(n)                  ! pot temp diff (K)
            delq   = qbot(n) - ssq                        ! spec hum dif (kg/kg)

            !--- UPDATE FLUX ITERATION ---

!!.................................................................
!! ocn_surface_flux_scheme = 0 : Default CESM1.2
!!                         = 1 : COARE algorithm
!!.................................................................
         if (ocn_surface_flux_scheme .eq. 0) then! use Large algorithm

            !--- compute stability & evaluate all stability functions ---
            hol  = shr_const_karman*shr_const_g*zbot(n)*  &
                   (tstar/thbot(n)+qstar/(1.0_R8/shr_const_zvir+qbot(n)))/ustar**2
            hol  = sign( min(abs(hol),10.0_R8), hol )
            stable = 0.5_R8 + sign(0.5_R8 , hol)
            xsq    = max(sqrt(abs(1.0_R8 - 16.0_R8*hol)) , 1.0_R8)
            xqq    = sqrt(xsq)
            psimh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psimhu(xqq)
            psixh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)

            !--- shift wind speed using old coefficient  and stability function  ---
            rd   = rdn / (1.0_R8 + rdn/shr_const_karman*(alz-psimh))
            u10n = vmag * rd / rdn

            !--- update neutral  transfer coeffs at 10m
            rdn    = sqrt(cdn(u10n))
            rhn    = (1.0_R8-stable) * 0.0327_R8 + stable * 0.018_R8
            ren    = 0.0346_R8

            !--- shift all coeffs to measurement height and stability ---
            rd = rdn / (1.0_R8 + rdn/shr_const_karman*(alz-psimh))
            rh = rhn / (1.0_R8 + rhn/shr_const_karman*(alz-psixh))
            re = ren / (1.0_R8 + ren/shr_const_karman*(alz-psixh))

            ustar = rd * vmag
            tstar = rh * delt
            qstar = re * delq

            !--- heat flux ---

            tau     = rbot(n) * ustar * ustar
            sen (n) =                cp * tau * tstar / ustar
            lat (n) = shr_const_latvap * tau * qstar / ustar

         else if (ocn_surface_flux_scheme .eq. 1) then! use COARE algorithm

            call cor30a(ubot(n),vbot(n),tbot(n),qbot(n),rbot(n) &  ! in atm params
                     & ,us(n),vs(n),tBulk(n),ssq                &  ! in surf params (NB ts -> tBulk)
                     & ,zpbl,zbot(n),zbot(n),zref,ztref,ztref   &  ! in heights
                     & ,tau,hsb,hlb                             &  ! out: fluxes
                     & ,zo,zot,zoq,hol,ustar,tstar,qstar        &  ! out: ss scales
                     & ,rd,rh,re                                &  ! out: exch. coeffs
                     & ,trf,qrf,urf,vrf)			       ! out: reference-height params
            ! for the sake of maintaining same defs
            hol=zbot(n)/hol
            rd=sqrt(rd)
            rh=sqrt(rh)
            re=sqrt(re)

            !--- heat flux ---

            sen (n) =  hsb
            lat (n) =  hlb

         else ! N.B.: NO ocn_surface_flux_scheme=2 option
               call shr_sys_abort(subName//", flux_diurnal requires ocn_surface_flux_scheme = 0 or 1")
         endif

         ENDDO   ! end iteration loop
         if (iter < 1) then
            call shr_sys_abort('No iterations performed ' // errMsg(sourcefile, __LINE__))
         end if
         !--- COMPUTE FLUXES TO ATMOSPHERE AND OCEAN ---

         ! Now calculated further up in subroutine.
         !tau = rbot(n) * ustar * ustar
         !sen (n) =                cp * tau * tstar / ustar
         !lat (n) =  shr_const_latvap * tau * qstar / ustar

         !--- momentum flux ---
         taux(n) = tau * (ubot(n)-us(n)) / vmag
         tauy(n) = tau * (vbot(n)-vs(n)) / vmag

         !--- LW radiation ---
         lwup(n) = -shr_const_stebol * Tskin(n)**4

         !--- water flux ---
         evap(n) = lat(n)/shr_const_latvap

         !---water isotope flux ---
!!ZZZ bugfix to be done
         call wiso_flxoce(2,rbot(n),zbot(n),s16O(n),ts(n),r16O(n),ustar,re,ssq, evap_16O(n),&
                          qbot(n),evap(n))
         call wiso_flxoce(3,rbot(n),zbot(n),sHDO(n),ts(n),rHDO(n),ustar,re,ssq, evap_HDO(n),&
                          qbot(n),evap(n))
         call wiso_flxoce(4,rbot(n),zbot(n),s18O(n),ts(n),r18O(n),ustar,re,ssq, evap_18O(n),&
                          qbot(n),evap(n))

         !------------------------------------------------------------
         ! compute diagnostics: 2m ref T & Q, 10m wind speed squared
         !------------------------------------------------------------

      if (ocn_surface_flux_scheme .eq. 0) then ! use Large algorithm

         hol = hol*ztref/zbot(n)
         xsq = max( 1.0_R8, sqrt(abs(1.0_R8-16.0_R8*hol)) )
         xqq = sqrt(xsq)
         psix2   = -5.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)
         fac     = (rh/shr_const_karman) * (alz + al2 - psixh + psix2 )
         tref(n) = thbot(n) - delt*fac
         tref(n) = tref(n) - 0.01_R8*ztref   ! pot temp to temp correction
         fac     = (re/shr_const_karman) * (alz + al2 - psixh + psix2 )
         qref(n) =  qbot(n) - delq*fac

         duu10n(n) = u10n*u10n ! 10m wind speed squared

      else if (ocn_surface_flux_scheme .eq. 1) then! use COARE algorithm

         tref(n) = trf
         qref(n) = qrf
         duu10n(n) = urf**2+vrf**2
         u10n = sqrt(duu10n(n))
      endif

         if (flux_diurnal) then

            !------------------------------------------------------------
            ! update new prognostic variables
            !------------------------------------------------------------

            warm  (n) = DTiter
            salt  (n) = DSiter
            speed (n) = DViter

            if (ltwopm) then
               tSkin_day(n) = tSkin(n)
               warmmax(n) = max(DTiter,0.0_R8)
            endif

            if (ltwoam) then
               tSkin_night(n) = tSkin(n)
               cSkin_night(n) = cSkin(n)
            endif

            if ((lmidnight).and.(lfullday)) then
               qSolAvg(n) = qSolInc(n)/real(nsum+1,R8)
               windAvg(n) = windInc(n)/real(nsum+1,R8)
               ! warmMax(n) = max(DTiter,warmMaxInc(n))
               windMax(n) = max(u10n,windMaxInc(n))

               nsum = 0

               qSolInc(n) = Qsol
               windInc(n) = u10n

               ! warmMaxInc(n) = 0.0_R8
               windMaxInc(n) = 0.0_R8

!               tSkin_night(n) = tSkin(n)
!               cSkin_night(n) = cSkin(n)

            else

               if ((lmidnight).and.(.not.(lfullday))) then

                  nsum = 0

                  qSolInc(n) = Qsol
                  windInc(n) = u10n

                  ! warmMaxInc(n) = 0.0_R8
                  windMaxInc(n) = 0.0_R8

               else

                  nsum = nsum + 1

                  ! warmMaxInc (n) = max(DTiter,warmMaxInc(n))
                  windMaxInc (n) = max(u10n, windMaxInc(n))
                  ! windMaxInc (n) = max(Qsol, windMaxInc(n))
                  qSolInc    (n) = qSolInc(n)+Qsol
                  windInc    (n) = windInc(n)+u10n

               endif
            endif

            nInc(n) = real(nsum,R8) ! set nInc to incremented or reset nsum


            if (present(ustar_sv)) ustar_sv(n) = ustar
            if (present(re_sv   )) re_sv   (n) = re
            if (present(ssq_sv  )) ssq_sv  (n) = ssq

         else              ! mask = 0

            !------------------------------------------------------------
            ! no valid data here -- out of domain
            !------------------------------------------------------------
            warm       (n) = spval ! NEW
            salt       (n) = spval ! NEW
            speed      (n) = spval ! NEW
            regime     (n) = spval ! NEW
            tBulk      (n) = spval ! NEW
            tSkin      (n) = spval ! NEW
            tSkin_night(n) = spval ! NEW
            tSkin_day  (n) = spval ! NEW
            cSkin      (n) = spval ! NEW
            cSkin_night(n) = spval ! NEW
            warmMax    (n) = spval ! NEW
            windMax    (n) = spval ! NEW
            qSolAvg    (n) = spval ! NEW
            windAvg    (n) = spval ! NEW
            warmMaxInc (n) = spval ! NEW
            windMaxInc (n) = spval ! NEW
            qSolInc    (n) = spval ! NEW
            windInc    (n) = spval ! NEW
            nInc       (n) = 0.0_R8 ! NEW

            sen   (n)    = spval  ! sensible         heat flux  (W/m^2)
            lat   (n)    = spval  ! latent           heat flux  (W/m^2)
            lwup  (n)    = spval  ! long-wave upward heat flux  (W/m^2)
            evap  (n)    = spval  ! evaporative water flux ((kg/s)/m^2)
            evap_16O (n) = spval  ! water tracer flux (kg/s)/m^2)
            evap_HDO (n) = spval  ! HDO tracer flux  (kg/s)/m^2)
            evap_18O (n) = spval  ! H218O tracer flux (kg/s)/m^2)
            taux  (n)    = spval  ! x surface stress (N)
            tauy  (n)    = spval  ! y surface stress (N)
            tref  (n)    = spval  ! 2m reference height temperature (K)
            qref  (n)    = spval  ! 2m reference height humidity (kg/kg)
            duu10n(n)    = spval  ! 10m wind speed squared (m/s)^2

            if (present(ustar_sv)) ustar_sv(n) = spval
            if (present(re_sv   )) re_sv   (n) = spval
            if (present(ssq_sv  )) ssq_sv  (n) = spval

         endif   ! mask

      endif ! flux diurnal logic

   ENDDO ! end n loop

END subroutine shr_flux_atmOcn_diurnal

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_flux_atmIce -- computes atm/ice fluxes
!
! !DESCRIPTION:
!    Computes atm/ice fluxes
!
! !REVISION HISTORY:
!    2006-Jun-12 - B. Kauffman, first version, adapted from dice6 code
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_flux_atmIce(mask  ,zbot  ,ubot  ,vbot  ,thbot  &
               &          ,qbot  ,rbot  ,tbot  ,ts    ,sen    &
               &          ,lat   ,lwup  ,evap  ,taux  ,tauy   &
               &          ,tref  ,qref                        )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

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

   !--- output arguments -------------------------------
   real(R8)   ,intent(out) :: sen  (:)    ! sensible      heat flux  (W/m^2)
   real(R8)   ,intent(out) :: lat  (:)    ! latent        heat flux  (W/m^2)
   real(R8)   ,intent(out) :: lwup (:)    ! long-wave upward heat flux  (W/m^2)
   real(R8)   ,intent(out) :: evap (:)    ! evaporative water flux ((kg/s)/m^2)
   real(R8)   ,intent(out) :: taux (:)    ! x surface stress (N)
   real(R8)   ,intent(out) :: tauy (:)    ! y surface stress (N)
   real(R8)   ,intent(out) :: tref (:)    ! 2m reference height temperature
   real(R8)   ,intent(out) :: qref (:)    ! 2m reference height humidity

!EOP

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
   character(*),parameter :: subName =  "(shr_flux_atmIce) "

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

end subroutine shr_flux_atmIce

!===============================================================================
! !BOP =========================================================================
!
! !IROUTINE: shr_flux_MOstability -- Monin-Obukhov BL stability functions
!
! !DESCRIPTION:
!
!    Monin-Obukhov boundary layer stability functions, two options:
!    turbulent velocity scales or gradient and integral functions
!    via option = shr_flux_MOwScales or shr_flux_MOfunctions
!
! !REVISION HISTORY:
!    2007-Sep-19 - B. Kauffman, Bill Large - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_flux_MOstability(option,arg1,arg2,arg3,arg4,arg5)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(IN),intent(in)           :: option ! shr_flux_MOwScales or MOfunctions
   real(R8)   ,intent(in)           :: arg1   ! scales: uStar (in)  funct: zeta (in)
   real(R8)   ,intent(inout)        :: arg2   ! scales: zkB   (in)  funct: phim (out)
   real(R8)   ,intent(out)          :: arg3   ! scales: phim  (out) funct: phis (out)
   real(R8)   ,intent(out)          :: arg4   ! scales: phis  (out) funct: psim (out)
   real(R8)   ,intent(out),optional :: arg5   ! scales:    (unused) funct: psis (out)

! !EOP

   !----- local variables -----
   real(R8)           :: zeta  ! z/L
   real(R8)           :: uStar ! friction velocity
   real(R8)           :: zkB   ! (height)*(von Karman)*(surface bouyancy flux)
   real(R8)           :: phim  ! momentum    gradient function or scale
   real(R8)           :: phis  ! temperature gradient function or scale
   real(R8)           :: psim  ! momentum    integral function or scale
   real(R8)           :: psis  ! temperature integral function or scale
   real(R8)           :: temp  ! temporary-variable/partial calculation

   !----- local variables, stable case -----
   real(R8),parameter :: uStarMin = 0.001_R8 ! lower bound on uStar
   real(R8),parameter :: a = 1.000_R8  ! constant from Holtslag & de Bruin, equation 12
   real(R8),parameter :: b = 0.667_R8  ! constant from Holtslag & de Bruin, equation 12
   real(R8),parameter :: c = 5.000_R8  ! constant from Holtslag & de Bruin, equation 12
   real(R8),parameter :: d = 0.350_R8  ! constant from Holtslag & de Bruin, equation 12

   !----- local variables, unstable case -----
   real(R8),parameter :: a2 = 3.0_R8   ! constant from Wilson, equation 10

   !----- formats -----
   character(*),parameter :: subName = '(shr_flux_MOstability) '
   character(*),parameter ::   F00 = "('(shr_flux_MOstability) ',4a)"
   character(*),parameter ::   F01 = "('(shr_flux_MOstability) ',a,i5)"

!-------------------------------------------------------------------------------
! Notes::
!   o this could be two routines, but are one to help keep them aligned
!   o the stable calculation is taken from...
!     A.A.M. HoltSlag and H.A.R. de Bruin, 1988:
!     "Applied Modeling of the Nighttime Surface Energy Balance over Land",
!     Journal of Applied Meteorology, Vol. 27, No. 6, June 1988, 659-704
!   o the unstable calculation is taken from...
!     D. Keith Wilson, 2001: "An Alternative Function for the Wind and
!     Temperature Gradients in Unstable Surface Layers",
!     Boundary-Layer Meteorology, 99 (2001), 151-158
!-------------------------------------------------------------------------------

   !----- check for consistancy between option and arguments ------------------
   if (debug > 1 .and. s_loglev > 0) then
      if (debug > 2) write(s_logunit,F01) "enter, option = ",option
      if ( option == shr_flux_MOwScales .and. present(arg5) ) then
         write(s_logunit,F01) "ERROR: option1 must have four arguments"
         call shr_sys_abort(subName//"option inconsistant with arguments")
      else if ( option == shr_flux_MOfunctions .and. .not. present(arg5) ) then
         write(s_logunit,F01) "ERROR: option2 must have five arguments"
         call shr_sys_abort(subName//"option inconsistant with arguments")
      else
         write(s_logunit,F01) "invalid option = ",option
         call shr_sys_abort(subName//"invalid option")
      end if
   end if

   !------ velocity scales option ----------------------------------------------
   if (option == shr_flux_MOwScales) then

      !--- input ---
      uStar = arg1
      zkB   = arg2

      if (zkB >= 0.0_R8) then ! ----- stable -----
         zeta = zkB/(max(uStar,uStarMin)**3)
         temp = exp(-d*zeta)
         phim = uStar/(1.0_R8 + zeta*(a + b*(1.0_R8 + c - d*zeta)*temp))
         phis = phim
      else                    ! ----- unstable -----
         temp = (zkB*zkB)**(1.0_R8/a2)   ! note: zkB < 0, zkB*zkB > 0
         phim = sqrt(uStar**2 + shr_flux_MOgammaM*temp)
         phis = sqrt(uStar**2 + shr_flux_MOgammaS*temp)
      end if

      !--- output ---
      arg3 = phim
      arg4 = phis
   !  arg5 = <unused>

   !------ stability function option -------------------------------------------
   else if (option == shr_flux_MOfunctions) then

      !--- input ---
      zeta  = arg1

      if (zeta >= 0.0_R8) then ! ----- stable -----
         temp = exp(-d*zeta)
         phim =        1.0_R8 + zeta*(a + b*(1.0_R8 + c - d*zeta)*temp)
         phis = phim
         psim = -a*zeta - b*(zeta - c/d)*temp - b*c/d
         psis = psim
      else                    ! ----- unstable ----
         temp = (zeta*zeta)**(1.0_R8/a2)   ! note: zeta < 0, zeta*zeta > 0
         phim = 1.0_R8/sqrt(1.0_R8 + shr_flux_MOgammaM*temp)
         phis = 1.0_R8/sqrt(1.0_R8 + shr_flux_MOgammaS*temp)
         psim = a2*log(0.5_R8 + 0.5_R8/phim)
         psis = a2*log(0.5_R8 + 0.5_R8/phis)
      end if

      !--- output ---
      arg2 = phim
      arg3 = phis
      arg4 = psim
      arg5 = psis
   !----------------------------------------------------------------------------
   else
      write(s_logunit,F01) "invalid option = ",option
      call shr_sys_abort(subName//"invalid option")
   endif

end subroutine shr_flux_MOstability

!===============================================================================
!===============================================================================

!===============================================================================
! !DESCRIPTION:
!
!   COARE v3.0 parametrisation
!
! !REVISION HISTORY:
!   2013-Nov-22: Thomas Toniazzo's adaptation of Chris Fairall's code,
!    downloaded from
!    ftp://ftp1.esrl.noaa.gov/users/cfairall/wcrp_wgsf/computer_programs/cor3_0/
!     * no wave, standard coare 2.6 charnock
!     * skin parametrisation also off (would require radiative fluxes and
!      rainrate in input)
!     * added diagnostics, comments and references
!===============================================================================
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cor30a(ubt,vbt,tbt,qbt,rbt        &    ! in atm params
               & ,uss,vss,tss,qss            &    ! in surf params
               & ,zbl,zbu,zbt,zrfu,zrfq,zrft &    ! in heights
               & ,tau,hsb,hlb                &    ! out: fluxes
               & ,zo,zot,zoq,L,usr,tsr,qsr   &    ! out: ss scales
               & ,Cd,Ch,Ce                   &    ! out: exch. coeffs
               & ,trf,qrf,urf,vrf)                ! out: reference-height params

! !USES:

IMPLICIT NONE

! !INPUT/OUTPUT PARAMETERS:

real(R8),intent(in) :: ubt,vbt,tbt,qbt,rbt,uss,vss,tss,qss
real(R8),intent(in) :: zbl,zbu,zbt,zrfu,zrfq,zrft
real(R8),intent(out):: tau,hsb,hlb,zo,zot,zoq,L,usr,tsr,qsr,Cd,Ch,Ce &
                    & ,trf,qrf,urf,vrf
! !EOP

real(R8) ua,va,ta,q,rb,us,vs,ts,qs,zi,zu,zt,zq,zru,zrq,zrt ! internal vars

real(R8):: cpa,rgas,grav,pi,von,beta ! phys. params
real(R8):: le,rhoa,cpv               ! derived phys. params
real(R8):: t,visa,du,dq,dt           ! params of problem

real(R8):: u10,zo10,zot10,cd10,ch10,ct10,ct,cc,ribu,zetu,l10,charn ! init vars
real(R8):: zet,rr,bf,ug,ut     ! loop iter vars
real(R8):: cdn_10,chn_10,cen_10  ! aux. output vars

integer(IN):: i,nits ! iter loop counters

integer(IN):: jcool                  ! aux. cool-skin vars
real(R8):: dter,wetc,dqer

ua=ubt  !wind components (m/s) at height zu (m)
va=vbt
ta=tbt  !bulk air temperature (K), height zt
Q =qbt  !bulk air spec hum (kg/kg), height zq
rb=rbt  ! air density
us=uss  !surface current components (m/s)
vs=vss
ts=tss  !bulk water temperature (K) if jcool=1, interface water T if jcool=0
qs=qss  !bulk water spec hum (kg/kg) if jcool=1 etc
zi=zbl  !PBL depth (m)
zu=zbu  !wind speed measurement height (m)
zt=zbt  !air T measurement height (m)
zq=zbt  !air q measurement height (m)
zru=zrfu ! reference height for st.diagn.U
zrq=zrfq ! reference height for st.diagn.T,q
zrt=zrft ! reference height for st.diagn.T,q

!**** constants
    Beta= 1.2_R8
    von = 0.4_R8
    pi  = 3.141593_R8
    grav= SHR_CONST_G
    Rgas= SHR_CONST_RGAS
    cpa = SHR_CONST_CPDAIR

!*** physical parameters
    Le  = SHR_CONST_LATVAP -.00237e6_R8*(ts-273.16_R8)
!   cpv = shr_const_cpdair*(1.0_R8 + shr_const_cpvir*Qs) ! form in NCAR code
    cpv = cpa*(1.0_R8+0.84_R8*Q)
!   rhoa= P/(Rgas*ta*(1+0.61*Q)) ! if input were pressure
    rhoa= rb

! parametrisation for air kinematic viscosity (Andreas 1989,p.31)
    t   = ta-273.16_R8
    visa= 1.326e-5_R8*(1.0_R8+6.542e-3_R8*t+8.301e-6_R8*t*t-4.84e-9_R8*t*t*t)

    du  = sqrt((ua-us)**2+(va-vs)**2)
    dt  = ts-ta -.0098_R8*zt
    dq  = Qs-Q

!*** don't use cool-skin params for now, but assign values to Ter and Qer
    jcool=0_IN
    dter=0.3_R8
    wetc=0.622_R8*Le*Qs/(Rgas*ts**2)
    dqer=wetc*dter

!***************** Begin bulk-model calculations ***************

!*************** first guess
    ug=0.5_R8

    ut   = sqrt(du*du+ug*ug)
    u10  = ut*log(10.0_R8/1.0e-4_R8)/log(zu/1.0e-4_R8)
    usr  = .035_R8*u10
    zo10 = 0.011_R8*usr*usr/grav+0.11_R8*visa/usr
    Cd10 = (von/log(10.0_R8/zo10))**2
    Ch10 = 0.00115_R8
    Ct10 = Ch10/sqrt(Cd10)
    zot10= 10.0_R8/exp(von/Ct10)
    Cd   =(von/log(zu/zo10))**2
    Ct   = von/log(zt/zot10)
    CC   = von*Ct/Cd

! Bulk Richardson number
    Ribu=-grav*zu/ta*((dt-dter*jcool)+.61_R8*ta*dq)/ut**2
! initial guess for stability parameter...
    if (Ribu .LT. 0.0_R8) then
    ! pbl-height dependent
        zetu=CC*Ribu/( 1.0_R8 - (.004_R8*Beta**3*zi/zu) * Ribu )
    else
        zetu=CC*Ribu*(1.0_R8 + 27.0_R8/9.0_R8*Ribu/CC)
    endif
! ...and MO length
    L10=zu/zetu

    if (zetu .GT. 50.0_R8) then
        nits=1_IN
    else
        nits=3_IN
    endif

    usr =  ut*von/(log(zu/zo10)-psiuo(zu/L10))
    tsr = (dt-dter*jcool)*von/(log(zt/zot10)-psit_30(zt/L10))
    qsr = (dq-dqer*jcool)*von/(log(zq/zot10)-psit_30(zq/L10))

! parametrisation for Charney parameter (section 3c of Fairall et al. 2003)
    charn=0.011_R8
    if (ut .GT. 10.0_R8) then
      charn=0.011_R8+(ut-10.0_R8)/(18.0_R8-10.0_R8)*(0.018_R8-0.011_R8)
    endif
    if (ut .GT. 18.0_R8) then
      charn=0.018_R8
    endif

!***************  iteration loop ************
    do i=1, nits

     ! stability parameter
     zet=-von*grav*zu/ta*(tsr*(1.0_R8+0.61_R8*Q)+.61_R8*ta*qsr)/(usr*usr)/(1.0_R8+0.61_R8*Q)

     ! momentum roughness length...
     zo = charn*usr*usr/grav+0.11_R8*visa/usr
     ! ...& MO length
     L  = zu/zet

     ! tracer roughness length
     rr = zo*usr/visa
     zoq= min(1.15e-4_R8,5.5e-5_R8/rr**.6_R8)
     zot= zoq ! N.B. same for vapour and heat

     ! new surface-layer scales
     usr =  ut            *von/(log(zu/zo )-psiuo(zu/L))
     tsr = (dt-dter*jcool)*von/(log(zt/zot)-psit_30(zt/L))
     qsr = (dq-dqer*jcool)*von/(log(zq/zoq)-psit_30(zq/L))

     ! gustiness parametrisation
     Bf=-grav/ta*usr*(tsr+.61_R8*ta*qsr)
     if (Bf .GT. 0.0_R8) then
       ug=Beta*(Bf*zi)**.333_R8
     else
       ug=.2_R8
     endif
     ut=sqrt(du*du+ug*ug)

    enddo
!***************     end loop    ************


   !******** fluxes @ measurement heights zu,zt,zq ********
   tau= rhoa*usr*usr*du/ut                !stress magnitude
   hsb=-rhoa*cpa*usr*tsr                  !heat downwards
   hlb=-rhoa*Le*usr*qsr                   !wv downwards

   !****** transfer coeffs relative to ut @meas. hts ******
   Cd= tau/rhoa/ut/max(.1_R8,du)
   if (tsr.ne.0._r8) then
    Ch= usr/ut*tsr/(dt-dter*jcool)
   else
    Ch= usr/ut* von/(log(zt/zot)-psit_30(zt/L))
   endif
   if (qsr.ne.0.0_R8) then
    Ce= usr/ut*qsr/(dq-dqer*jcool)
   else
    Ce= usr/ut* von/(log(zq/zoq)-psit_30(zq/L))
   endif

   !**********  10-m neutral coeff relative to ut *********
   Cdn_10=von*von/log(10.0_R8/zo)/log(10.0_R8/zo)
   Chn_10=von*von/log(10.0_R8/zo)/log(10.0_R8/zot)
   Cen_10=von*von/log(10.0_R8/zo)/log(10.0_R8/zoq)

   !**********  reference-height values for u,q,T *********
   urf=us+(ua-us)*(log(zru/zo)-psiuo(zru/L))/(log(zu/zo)-psiuo(zu/L))
   vrf=vs+(va-vs)*(log(zru/zo)-psiuo(zru/L))/(log(zu/zo)-psiuo(zu/L))
   qrf=qs-dq*(log(zrq/zoq)-psit_30(zrq/L))/(log(zq/zoq)-psit_30(zq/L))
   trf=ts-dt*(log(zrt/zot)-psit_30(zrt/L))/(log(zt/zot)-psit_30(zt/L))
   trf=trf+.0098_R8*zrt

end subroutine cor30a


!===============================================================================
! !BOP =========================================================================
!
! !IROUTINE: PSIUo
!
! !DESCRIPTION:
!
!   momentum stability functions adopted in COARE v3.0 parametrisation.
!   Chris Fairall's code (see cor30a)
!
! !REVISION HISTORY:
!   22/11/2013: Thomas Toniazzo: comments added
!
! !INTERFACE: ------------------------------------------------------------------
real (R8) function psiuo(zet)
! !INPUT/OUTPUT PARAMETERS:
real(R8),intent(in)  :: zet
! !EOP
real(R8) ::c,x,psik,psic,f
!-----------------------------------------------------------------
! N.B.: z0/L always neglected compared to z/L and to 1
!-----------------------------------------------------------------
    if(zet>0.0_R8)then
! Beljaars & Holtslag (1991)
     c=min(50._R8,.35_R8*zet)
     psiuo=-((1.0_R8+1.0_R8*zet)**1.0_R8+.667_R8*(zet-14.28_R8)/exp(c)+8.525_R8)
    else
! Dyer & Hicks (1974) for weak instability
     x=(1.0_R8-15.0_R8*zet)**.25_R8                   ! 15 instead of 16
     psik=2.0_R8*log((1.0_R8+x)/2.0_R8)+log((1.0_R8+x*x)/2.0_R8)-2.0_R8*atan(x)+2.0_R8*atan(1.0_R8)
! Fairall et al. (1996) for strong instability (Eq.(13))
     x=(1.0_R8-10.15_R8*zet)**.3333_R8
     psic= 1.5_R8*log((1.0_R8+x+x*x)/3.0_R8)-sqrt(3.0_R8)*atan((1.0_R8+2.0_R8*x)/sqrt(3.0_R8)) &
         & +4.0_R8*atan(1.0_R8)/sqrt(3.0_R8)
     f=zet*zet/(1.0_R8+zet*zet)
     psiuo=(1.0_R8-f)*psik+f*psic
    endif
END FUNCTION psiuo



!===============================================================================
! !BOP =========================================================================
!
! !IROUTINE: PSIT_30
!
! !DESCRIPTION:
!
!   momentum stability functions adopted in COARE v3.0 parametrisation.
!   Chris Fairall's code (see cor30a)
!
! !REVISION HISTORY:
!   22/11/2013: Thomas Toniazzo: comments added
!
! !INTERFACE: ------------------------------------------------------------------
real (R8) function psit_30(zet)
! !INPUT/OUTPUT PARAMETERS:
real(R8),intent(in)  :: zet
! !EOP
real(R8) ::c,x,psik,psic,f
!-----------------------------------------------------------------
! N.B.: z0/L always neglected compared to z/L and to 1
!-----------------------------------------------------------------
    if(zet>0.0_R8)then
! Beljaars & Holtslag (1991)
     c=min(50._R8,.35_R8*zet)
     psit_30=-((1.0_R8+2.0_R8/3.0_R8*zet)**1.5_R8+.667_R8*(zet-14.28_R8)/exp(c)+8.525_R8)
    else
! Dyer & Hicks (1974) for weak instability
     x=(1.0_R8-15.0_R8*zet)**.5_R8                    ! 15 instead of 16
     psik=2.0_R8*log((1.0_R8+x)/2.0_R8)
! Fairall et al. (1996) for strong instability
     x=(1.0_R8-(34.15_R8*zet))**.3333_R8
     psic= 1.5_R8*log((1.0_R8+x+x*x)/3.0_R8)-sqrt(3.0_R8)*atan((1.0_R8+2.0_R8*x)/sqrt(3.0_R8)) &
         & +4.0_R8*atan(1.0_R8)/sqrt(3.0_R8)
     f=zet*zet/(1.0_R8+zet*zet)
     psit_30=(1.0_R8-f)*psik+f*psic
   endif
end FUNCTION psit_30



!!! !DESCRIPTION:
!!!     set docoare flag
!!!     \newline
!!!     call shr\_flux\_setDopole(flag)
!!!
!!! !REVISION HISTORY:
!!!     2009-Jun-22 - T. Craig - first version
!!!
!!! !INTERFACE: ------------------------------------------------------------------
!!
!!subroutine shr_flux_docoare(iflag)
!!
!!  implicit none
!!
!!! !INPUT/OUTPUT PARAMETERS:
!!
!!  integer, intent(in) :: iflag
!!
!!
!!  flux_scheme = iflag
!!
!!
!!end subroutine shr_flux_docoare


end module shr_flux_mod
