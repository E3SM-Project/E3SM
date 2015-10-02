!===============================================================================
! SVN $Id: shr_flux_mod.F90 25015 2010-10-06 21:06:16Z olson $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/branches/cam3_restore_cam5_0_15/shr/shr_flux_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: flux_mod -- CCSM shared flux calculations.
!
! !DESCRIPTION:
!
!     CCSM shared flux calculations.
!     
! !REVISION HISTORY:
!     2006-Nov-07 - B. Kauffman - first version, code taken/migrated from cpl6
!     2010-Dec-21 - J. Olson - revised for CAM aqua-planet use
!
! !INTERFACE: ------------------------------------------------------------------

module shr_flux_mod

! !USES:

   use shr_kind_mod    ! shared kinds
   use physconst       ! shared constants (CAM defined)
   use shr_sys_mod     ! shared system routines
   use shr_log_mod, only: s_loglev  => shr_log_Level
   use shr_log_mod, only: s_logunit => shr_log_Unit
   use phys_control,only: cam_physpkg_is

   implicit none

   private ! default private

! !PUBLIC TYPES:

  ! none

! !PUBLIC MEMBER FUNCTIONS:

   public :: shr_flux_atmOcn      ! computes atm/ocn fluxes
   public :: shr_flux_atmIce      ! computes atm/ice fluxes

!EOP

   !--- rename kinds for local readability only ---
   integer,parameter :: R8 = SHR_KIND_R8  ! 8 byte real
   integer,parameter :: IN = SHR_KIND_IN  ! native/default integer

   integer,parameter :: debug = 0 ! internal debug level

!===============================================================================
contains
!===============================================================================
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
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE shr_flux_atmOcn(nMax  ,zbot  ,ubot  ,vbot  ,thbot ,   & 
           &               qbot  ,rbot  ,tbot  ,us    ,vs    ,   &
           &               ts    ,mask  ,sen   ,lat   ,lwup  ,   &
           &               evap  ,taux  ,tauy  ,tref  ,qref  ,   &
           &               duu10n,  ustar_sv   ,re_sv ,ssq_sv,   &
           &               missval    )

! !USES:

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
   real(R8)   ,intent(in) :: rbot (nMax) ! atm air density       (kg/m^3)
   real(R8)   ,intent(in) :: tbot (nMax) ! atm T                 (K) 
   real(R8)   ,intent(in) :: us   (nMax) ! ocn u-velocity        (m/s)
   real(R8)   ,intent(in) :: vs   (nMax) ! ocn v-velocity        (m/s)
   real(R8)   ,intent(in) :: ts   (nMax) ! ocn temperature       (K)

   !--- output arguments -------------------------------
   real(R8),intent(out)  ::  sen  (nMax) ! heat flux: sensible    (W/m^2)
   real(R8),intent(out)  ::  lat  (nMax) ! heat flux: latent      (W/m^2)
   real(R8),intent(out)  ::  lwup (nMax) ! heat flux: lw upward   (W/m^2)
   real(R8),intent(out)  ::  evap (nMax) ! water flux: evap  ((kg/s)/m^2)
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

   !--- local variables --------------------------------
   integer(IN) :: n      ! vector loop index
   real(R8)    :: umin   ! minimum wind speed       (m/s)
   real(R8)    :: vmag   ! surface wind magnitude   (m/s)
   real(R8)    :: thvbot ! virtual temperature      (K)
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
   real(R8)    :: bn     ! exchange coef funct for interpolation
   real(R8)    :: bh     ! exchange coef funct for interpolation
   real(R8)    :: fac    ! vertical interpolation factor

   !--- local functions --------------------------------
   real(R8)    :: qsat   ! function: the saturation humididty of air (kg/m^3)
   real(R8)    :: cdn    ! function: neutral drag coeff at 10m
   real(R8)    :: psimhu ! function: unstable part of psimh
   real(R8)    :: psixhu ! function: unstable part of psimx
   real(R8)    :: Umps   ! dummy arg ~ wind velocity (m/s)
   real(R8)    :: Tk     ! dummy arg ~ temperature (K)
   real(R8)    :: xd     ! dummy arg ~ ?
 
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
!   o Neutral 10m drag coeff: cdn = .0027/U10 + .000142 + .0000764 U10
!   o Neutral 10m stanton number: ctn = .0327 sqrt(cdn), unstable
!                                 ctn = .0180 sqrt(cdn), stable
!   o Neutral 10m dalton number:  cen = .0346 sqrt(cdn)
!   o The saturation humidity of air at T(K): qsat(T)  (kg/m^3)
!-------------------------------------------------------------------------------

   if (debug > 0 .and. s_loglev > 0) write(s_logunit,F00) "enter"

   if( cam_physpkg_is('cam3') ) then
      umin  =  1.0_R8
   else
      umin  =  0.5_R8
   endif

   al2 = log(zref/ztref)

   DO n=1,nMax
     if (mask(n) /= 0) then
    
        !--- compute some needed quantities ---
        vmag   = max(umin, sqrt( (ubot(n)-us(n))**2 + (vbot(n)-vs(n))**2) )
        thvbot = thbot(n) * (1.0_R8 + zvir * qbot(n)) ! virtual temp (K)
        ssq    = 0.98_R8 * qsat(ts(n)) / rbot(n)   ! sea surf hum (kg/kg)
        delt   = thbot(n) - ts(n)                  ! pot temp diff (K)
        delq   = qbot(n) - ssq                     ! spec hum dif (kg/kg)
        alz    = log(zbot(n)/zref) 
        cp     = cpair*(1.0_R8 + cpvir*ssq) 
   
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
   
        !--- compute stability & evaluate all stability functions ---
        if( cam_physpkg_is('cam3') ) then
           hol  = karman*gravit*zbot(n)*  &
                  (tstar/thvbot+qstar/(1.0_R8/zvir+qbot(n)))/ustar**2
        else
           hol  = karman*gravit*zbot(n)*  &
                  (tstar/thbot(n)+qstar/(1.0_R8/zvir+qbot(n)))/ustar**2
        endif
        hol  = sign( min(abs(hol),10.0_R8), hol )
        stable = 0.5_R8 + sign(0.5_R8 , hol)
        xsq    = max(sqrt(abs(1.0_R8 - 16.0_R8*hol)) , 1.0_R8)
        xqq    = sqrt(xsq)
        psimh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psimhu(xqq)
        psixh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)
   
        !--- shift wind speed using old coefficient ---
        rd   = rdn / (1.0_R8 + rdn/karman*(alz-psimh))
        u10n = vmag * rd / rdn 
   
        !--- update transfer coeffs at 10m and neutral stability ---
        rdn = sqrt(cdn(u10n))
        ren = 0.0346_R8
        rhn = (1.0_R8-stable)*0.0327_R8 + stable * 0.018_R8 
    
        !--- shift all coeffs to measurement height and stability ---
        rd = rdn / (1.0_R8 + rdn/karman*(alz-psimh)) 
        rh = rhn / (1.0_R8 + rhn/karman*(alz-psixh)) 
        re = ren / (1.0_R8 + ren/karman*(alz-psixh)) 
   
        !--- update ustar, tstar, qstar using updated, shifted coeffs --
        ustar = rd * vmag 
        tstar = rh * delt 
        qstar = re * delq 
    
        !------------------------------------------------------------
        ! iterate to converge on Z/L, ustar, tstar and qstar
        !------------------------------------------------------------
    
        !--- compute stability & evaluate all stability functions ---
        if( cam_physpkg_is('cam3') ) then
           hol  = karman*gravit*zbot(n)* &
                  (tstar/thvbot+qstar/(1.0_R8/zvir+qbot(n)))/ustar**2
        else
           hol  = karman*gravit*zbot(n)* &
                  (tstar/thbot(n)+qstar/(1.0_R8/zvir+qbot(n)))/ustar**2
        endif
        hol  = sign( min(abs(hol),10.0_R8), hol )
        stable = 0.5_R8 + sign(0.5_R8 , hol)
        xsq    = max(sqrt(abs(1.0_R8 - 16.0_R8*hol)) , 1.0_R8)
        xqq    = sqrt(xsq)
        psimh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psimhu(xqq)
        psixh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)
    
        !--- shift wind speed using old coeffs ---
        rd   = rdn / (1.0_R8 + rdn/karman*(alz-psimh))
        u10n = vmag * rd/rdn 
    
        !--- update transfer coeffs at 10m and neutral stability ---
        rdn = sqrt(cdn(u10n))
        ren = 0.0346_R8
        rhn = (1.0_R8 - stable)*0.0327_R8 + stable * 0.018_R8 
   
        !--- shift all coeffs to measurement height and stability ---
        rd = rdn / (1.0_R8 + rdn/karman*(alz-psimh)) 
        rh = rhn / (1.0_R8 + rhn/karman*(alz-psixh)) 
        re = ren / (1.0_R8 + ren/karman*(alz-psixh)) 
    
        !--- update ustar, tstar, qstar using updated, shifted coeffs ---
        ustar = rd * vmag 
        tstar = rh * delt 
        qstar = re * delq 
    
        !------------------------------------------------------------
        ! compute the fluxes
        !------------------------------------------------------------
    
        tau = rbot(n) * ustar * ustar 
       
        !--- momentum flux ---
        taux(n) = tau * (ubot(n)-us(n)) / vmag 
        tauy(n) = tau * (vbot(n)-vs(n)) / vmag 
        
        !--- heat flux ---
        sen (n) =                cp * tau * tstar / ustar 
        lat (n) =  latvap * tau * qstar / ustar
        lwup(n) = -stebol * ts(n)**4 
      
        !--- water flux ---
        evap(n) = lat(n)/latvap 
    
        !------------------------------------------------------------
        ! compute diagnositcs: 2m ref T & Q, 10m wind speed squared
        !------------------------------------------------------------
        hol = hol*ztref/zbot(n)
        xsq = max( 1.0_R8, sqrt(abs(1.0_R8-16.0_R8*hol)) )
        xqq = sqrt(xsq)
        psix2   = -5.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)
        fac     = (rh/karman) * (alz + al2 - psixh + psix2 )
        tref(n) = thbot(n) - delt*fac 
        tref(n) = tref(n) - 0.01_R8*ztref   ! pot temp to temp correction
        fac     = (re/karman) * (alz + al2 - psixh + psix2 )
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

END subroutine shr_flux_atmOcn

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
        thvbot = thbot(n)*(1.0_R8 + zvir * qbot(n)) ! virtual pot temp (K)
         ssq   =  qsat  (ts(n)) / rbot(n)           ! sea surf hum (kg/kg)
        dssqdt = dqsatdt(ts(n)) / rbot(n)           ! deriv of ssq wrt Ts 
        delt   = thbot(n) - ts(n)                   ! pot temp diff (K)
        delq   = qbot(n) - ssq                        ! spec hum dif (kg/kg)
        alz    = log(zbot(n)/zref) 
        cp     = cpair*(1.0_R8 + cpvir*ssq) 
        ltheat = latvap + latice

        !----------------------------------------------------------
        ! first estimate of Z/L and ustar, tstar and qstar
        !----------------------------------------------------------

        !--- neutral coefficients, z/L = 0.0 ---
        rdn = karman/log(zref/zzsice)
        rhn = rdn
        ren = rdn

        !--- ustar,tstar,qstar ----
        ustar = rdn * vmag
        tstar = rhn * delt  
        qstar = ren * delq  

        !--- compute stability & evaluate all stability functions ---
        hol    = karman * gravit * zbot(n) &
        &     * (tstar/thvbot+qstar/(1.0_R8/zvir+qbot(n))) / ustar**2
        hol    = sign( min(abs(hol),10.0_R8), hol )
        stable = 0.5_R8 + sign(0.5_R8 , hol)
        xsq    = max(sqrt(abs(1.0_R8 - 16.0_R8*hol)) , 1.0_R8)
        xqq    = sqrt(xsq)
        psimh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psimhu(xqq)
        psixh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)

        !--- shift all coeffs to measurement height and stability ---
        rd = rdn / (1.0_R8+rdn/karman*(alz-psimh))
        rh = rhn / (1.0_R8+rhn/karman*(alz-psixh))
        re = ren / (1.0_R8+ren/karman*(alz-psixh))

        !--- update ustar, tstar, qstar w/ updated, shifted coeffs --
        ustar = rd * vmag 
        tstar = rh * delt 
        qstar = re * delq 

        !----------------------------------------------------------
        ! iterate to converge on Z/L, ustar, tstar and qstar
        !----------------------------------------------------------

        !--- compute stability & evaluate all stability functions ---
        hol    = karman * gravit * zbot(n) &
        &      * (tstar/thvbot+qstar/(1.0_R8/zvir+qbot(n))) / ustar**2
        hol    = sign( min(abs(hol),10.0_R8), hol )
        stable = 0.5_R8 + sign(0.5_R8 , hol)
        xsq    = max(sqrt(abs(1.0_R8 - 16.0_R8*hol)) , 1.0_R8)
        xqq    = sqrt(xsq)
        psimh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psimhu(xqq)
        psixh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)

        !--- shift all coeffs to measurement height and stability ---
        rd = rdn / (1.0_R8+rdn/karman*(alz-psimh)) 
        rh = rhn / (1.0_R8+rhn/karman*(alz-psixh)) 
        re = ren / (1.0_R8+ren/karman*(alz-psixh)) 

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
        lwup(n) = -stebol * ts(n)**4 
     
        !--- water flux ---
        evap(n) = lat(n)/ltheat 

        !----------------------------------------------------------
        ! compute diagnostic: 2m reference height temperature
        !----------------------------------------------------------

        !--- Compute function of exchange coefficients. Assume that 
        !--- cn = rdn*rdn, cm=rd*rd and ch=rh*rd, and therefore 
        !--- 1/sqrt(cn(n))=1/rdn and sqrt(cm(n))/ch(n)=1/rh 
        bn = karman/rdn
        bh = karman/rh

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
!===============================================================================

end module shr_flux_mod
