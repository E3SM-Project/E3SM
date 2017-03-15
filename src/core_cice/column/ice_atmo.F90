!  SVN:$Id: ice_atmo.F90 1173 2017-03-02 03:57:43Z njeffery $
!=======================================================================

! Atmospheric boundary interface (stability based flux calculations)

! author: Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb
! 2004: Block structure added by William Lipscomb
! 2006: Converted to free source form (F90) by Elizabeth Hunke
! 2013: Form drag routine added (neutral_drag_coeffs) by David Schroeder 
! 2014: Adjusted form drag and added high frequency coupling by Andrew Roberts

      module ice_atmo

      use ice_kinds_mod
      use ice_constants_colpkg, only: c0, c1, c2, c4, c5, c8, c10, &
           c16, c20, p001, p01, p2, p4, p5, p75, puny, &
           cp_wv, cp_air, iceruf, zref, qqqice, TTTice, qqqocn, TTTocn, &
           Lsub, Lvap, vonkar, Tffresh, zvir, gravit, &
           pih, dragio, rhoi, rhos, rhow

      implicit none
      save

      private
      public :: atmo_boundary_layer, atmo_boundary_const, neutral_drag_coeffs

!=======================================================================

      contains

!=======================================================================

! Compute coefficients for atm/ice fluxes, stress, and reference
! temperature and humidity. NOTE: \\
! (1) All fluxes are positive downward,  \\
! (2) Here, tstar = (WT)/U*, and qstar = (WQ)/U*,  \\
! (3a) wind speeds should all be above a minimum speed (eg. 1.0 m/s). \\
!
! ASSUME:
!  The saturation humidity of air at T(K): qsat(T)  (kg/m**3)
!
! Code originally based on CSM1

      subroutine atmo_boundary_layer (sfctype,            &
                                      calc_strair, formdrag, &
                                      highfreq, natmiter, &
                                      Tsf,      potT,     &
                                      uatm,     vatm,     &  
                                      wind,     zlvl,     &  
                                      Qa,       rhoa,     &
                                      strx,     stry,     &   
                                      Tref,     Qref,     &
                                      delt,     delq,     &
                                      lhcoef,   shcoef,   &
                                      Cdn_atm,            &
                                      Cdn_atm_ratio_n,    &
                                      uvel,     vvel,     &
                                      Uref                )     

      character (len=3), intent(in) :: &
         sfctype      ! ice or ocean

      logical (kind=log_kind), intent(in) :: &
         calc_strair, &  ! if true, calculate wind stress components
         formdrag,    &  ! if true, calculate form drag
         highfreq        ! if true, use high frequency coupling

      integer (kind=int_kind), intent(in) :: &
         natmiter        ! number of iterations for boundary layer calculations

      real (kind=dbl_kind), intent(in) :: &
         Tsf      , & ! surface temperature of ice or ocean
         potT     , & ! air potential temperature  (K)
         uatm     , & ! x-direction wind speed (m/s)
         vatm     , & ! y-direction wind speed (m/s)
         wind     , & ! wind speed (m/s)
         zlvl     , & ! atm level height (m)
         Qa       , & ! specific humidity (kg/kg)
         rhoa         ! air density (kg/m^3)

      real (kind=dbl_kind), intent(inout) :: &
         Cdn_atm      ! neutral drag coefficient
 
      real (kind=dbl_kind), intent(out) :: &
         Cdn_atm_ratio_n ! ratio drag coeff / neutral drag coeff

      real (kind=dbl_kind), &
         intent(inout) :: &
         strx     , & ! x surface stress (N)
         stry         ! y surface stress (N)

      real (kind=dbl_kind), intent(inout) :: &
         Tref     , & ! reference height temperature  (K)
         Qref     , & ! reference height specific humidity (kg/kg)
         delt     , & ! potential T difference   (K)
         delq     , & ! humidity difference      (kg/kg)
         shcoef   , & ! transfer coefficient for sensible heat
         lhcoef       ! transfer coefficient for latent heat

      real (kind=dbl_kind), intent(in) :: &
         uvel     , & ! x-direction ice speed (m/s)
         vvel         ! y-direction ice speed (m/s)

      real (kind=dbl_kind), intent(out) :: &
         Uref         ! reference height wind speed (m/s)

      ! local variables

      integer (kind=int_kind) :: &
         k         ! iteration index

      real (kind=dbl_kind) :: &
         TsfK  , & ! surface temperature in Kelvin (K)
         xqq   , & ! temporary variable
         psimh , & ! stability function at zlvl   (momentum)
         tau   , & ! stress at zlvl
         fac   , & ! interpolation factor
         al2   , & ! ln(z10   /zTrf)
         psix2 , & ! stability function at zTrf   (heat and water)
         psimhs, & ! stable profile
         ssq   , & ! sat surface humidity     (kg/kg)
         qqq   , & ! for qsat, dqsfcdt
         TTT   , & ! for qsat, dqsfcdt
         qsat  , & ! the saturation humidity of air (kg/m^3)
         Lheat , & ! Lvap or Lsub, depending on surface type
         umin      ! minimum wind speed (m/s)

      real (kind=dbl_kind) :: &
         ustar , & ! ustar (m/s)
         tstar , & ! tstar
         qstar , & ! qstar
         rdn   , & ! sqrt of neutral exchange coefficient (momentum)
         rhn   , & ! sqrt of neutral exchange coefficient (heat)
         ren   , & ! sqrt of neutral exchange coefficient (water)
         rd    , & ! sqrt of exchange coefficient (momentum)
         re    , & ! sqrt of exchange coefficient (water)
         rh    , & ! sqrt of exchange coefficient (heat)
         vmag  , & ! surface wind magnitude   (m/s)
         alz   , & ! ln(zlvl  /z10)
         thva  , & ! virtual temperature      (K)
         cp    , & ! specific heat of moist air
         hol   , & ! H (at zlvl  ) over L
         stable, & ! stability factor
         psixh     ! stability function at zlvl   (heat and water)

      real (kind=dbl_kind), parameter :: &
         cpvir = cp_wv/cp_air-c1, & ! defined as cp_wv/cp_air - 1.
         zTrf  = c2                 ! reference height for air temp (m)

      ! local functions
      real (kind=dbl_kind) :: &
         xd    , & ! dummy argument
         psimhu, & ! unstable part of psimh
         psixhu    ! unstable part of psimx

      !------------------------------------------------------------
      ! Define functions
      !------------------------------------------------------------

      psimhu(xd) = log((c1+xd*(c2+xd))*(c1+xd*xd)/c8) &
                 - c2*atan(xd) + pih
!ech                 - c2*atan(xd) + 1.571_dbl_kind

      psixhu(xd) =  c2 * log((c1 + xd*xd)/c2)

      al2 = log(zref/zTrf)

      !------------------------------------------------------------
      ! Initialize
      !------------------------------------------------------------

      if (highfreq) then       
       umin  = p5 ! minumum allowable wind-ice speed difference of 0.5 m/s
      else
       umin  = c1 ! minumum allowable wind speed of 1m/s
      endif

      Tref = c0
      Qref = c0
      Uref = c0
      delt = c0
      delq = c0
      shcoef = c0
      lhcoef = c0

      !------------------------------------------------------------
      ! Compute turbulent flux coefficients, wind stress, and
      ! reference temperature and humidity.
      !------------------------------------------------------------

      !------------------------------------------------------------
      ! define variables that depend on surface type
      !------------------------------------------------------------

      if (sfctype(1:3)=='ice') then

         qqq   = qqqice          ! for qsat
         TTT   = TTTice          ! for qsat
         Lheat = Lsub            ! ice to vapor

            if (highfreq) then
               vmag = max(umin, sqrt( (uatm-uvel)**2 + &
                                      (vatm-vvel)**2) )
            else
               vmag = max(umin, wind)
            endif

            if (formdrag .and. Cdn_atm > puny) then 
               rdn = sqrt(Cdn_atm)               
            else
               rdn  = vonkar/log(zref/iceruf) ! neutral coefficient
               Cdn_atm = rdn * rdn
            endif

      elseif (sfctype(1:3)=='ocn') then

         qqq   = qqqocn
         TTT   = TTTocn
         Lheat = Lvap           ! liquid to vapor
         vmag = max(umin, wind)
         rdn  = sqrt(0.0027_dbl_kind/vmag &
              + .000142_dbl_kind + .0000764_dbl_kind*vmag)

      endif   ! sfctype

      !------------------------------------------------------------
      ! define some more needed variables
      !------------------------------------------------------------

      TsfK = Tsf + Tffresh     ! surface temp (K)
      qsat = qqq * exp(-TTT/TsfK)   ! saturation humidity (kg/m^3)
      ssq  = qsat / rhoa       ! sat surf hum (kg/kg)
      
      thva = potT * (c1 + zvir * Qa) ! virtual pot temp (K)
      delt = potT - TsfK       ! pot temp diff (K)
      delq = Qa - ssq          ! spec hum dif (kg/kg)
      alz  = log(zlvl/zref)
      cp   = cp_air*(c1 + cpvir*ssq)
      
      !------------------------------------------------------------
      ! first estimate of Z/L and ustar, tstar and qstar
      !------------------------------------------------------------

      ! neutral coefficients, z/L = 0.0
      rhn = rdn
      ren = rdn
      
      ! ustar,tstar,qstar
      ustar = rdn * vmag
      tstar = rhn * delt
      qstar = ren * delq

      !------------------------------------------------------------
      ! iterate to converge on Z/L, ustar, tstar and qstar
      !------------------------------------------------------------

      do k = 1, natmiter

         ! compute stability & evaluate all stability functions
         hol = vonkar * gravit * zlvl &
                 * (tstar/thva &
                 + qstar/(c1/zvir+Qa)) &
                 / ustar**2
         hol    = sign( min(abs(hol),c10), hol )
         stable = p5 + sign(p5 , hol)
         xqq    = max(sqrt(abs(c1 - c16*hol)) , c1)
         xqq    = sqrt(xqq)

         ! Jordan et al 1999
         psimhs = -(0.7_dbl_kind*hol &
                + 0.75_dbl_kind*(hol-14.3_dbl_kind) &
                * exp(-0.35_dbl_kind*hol) + 10.7_dbl_kind)
         psimh  = psimhs*stable &
                + (c1 - stable)*psimhu(xqq)
         psixh  = psimhs*stable &
                + (c1 - stable)*psixhu(xqq)

         ! shift all coeffs to measurement height and stability
         rd = rdn / (c1+rdn/vonkar*(alz-psimh))
         rh = rhn / (c1+rhn/vonkar*(alz-psixh))
         re = ren / (c1+ren/vonkar*(alz-psixh))
         
         ! update ustar, tstar, qstar using updated, shifted coeffs
         ustar = rd * vmag
         tstar = rh * delt
         qstar = re * delq

      enddo                     ! end iteration

      if (calc_strair) then

         ! initialize
         strx = c0
         stry = c0

         if (highfreq .and. sfctype(1:3)=='ice') then

            !------------------------------------------------------------
            ! momentum flux for RASM
            !------------------------------------------------------------
            ! tau = rhoa * rd * rd
            ! strx = tau * |Uatm-U| * (uatm-u)
            ! stry = tau * |Uatm-U| * (vatm-v)
            !------------------------------------------------------------

            tau = rhoa * rd * rd ! not the stress at zlvl

            ! high frequency momentum coupling following Roberts et al. (2014)
            strx = tau * sqrt((uatm-uvel)**2 + (vatm-vvel)**2) * (uatm-uvel)
            stry = tau * sqrt((uatm-uvel)**2 + (vatm-vvel)**2) * (vatm-vvel)

         else

            !------------------------------------------------------------
            ! momentum flux
            !------------------------------------------------------------
            ! tau = rhoa * ustar * ustar
            ! strx = tau * uatm / vmag
            ! stry = tau * vatm / vmag
            !------------------------------------------------------------

            tau = rhoa * ustar * rd ! not the stress at zlvl
            strx = tau * uatm
            stry = tau * vatm

         endif

         Cdn_atm_ratio_n = rd * rd / rdn / rdn

      endif                     ! calc_strair

      !------------------------------------------------------------
      ! coefficients for turbulent flux calculation
      !------------------------------------------------------------
      ! add windless coefficient for sensible heat flux
      ! as in Jordan et al (JGR, 1999)
      !------------------------------------------------------------

      shcoef = rhoa * ustar * cp * rh + c1
      lhcoef = rhoa * ustar * Lheat  * re

      !------------------------------------------------------------
      ! Compute diagnostics: 2m ref T, Q, U
      !------------------------------------------------------------

      hol   = hol*zTrf/zlvl
      xqq   = max( c1, sqrt(abs(c1-c16*hol)) )
      xqq   = sqrt(xqq)
      psix2 = -c5*hol*stable + (c1-stable)*psixhu(xqq)
      fac   = (rh/vonkar) &
            * (alz + al2 - psixh + psix2)
      Tref  = potT - delt*fac
      Tref  = Tref - p01*zTrf ! pot temp to temp correction
      fac   = (re/vonkar) &
            * (alz + al2 - psixh + psix2)
      Qref  = Qa - delq*fac

      if (highfreq .and. sfctype(1:3)=='ice') then
         Uref = sqrt((uatm-uvel)**2 + (vatm-vvel)**2) * rd / rdn
      else
         Uref = vmag * rd / rdn
      endif

      end subroutine atmo_boundary_layer

!=======================================================================

! Compute coefficients for atm/ice fluxes, stress
! NOTE: \\
! (1) all fluxes are positive downward,  \\
! (2) reference temperature and humidity are NOT computed

      subroutine atmo_boundary_const (sfctype,  calc_strair, &
                                      uatm,     vatm,     &  
                                      wind,     rhoa,     &
                                      strx,     stry,     &   
                                      Tsf,      potT,     &
                                      Qa,                 &
                                      delt,     delq,     &
                                      lhcoef,   shcoef,   &
                                      Cdn_atm)  

      character (len=3), intent(in) :: &
         sfctype      ! ice or ocean

      logical (kind=log_kind), intent(in) :: &
         calc_strair  ! if true, calculate wind stress components

      real (kind=dbl_kind), intent(in) :: &
         Tsf      , & ! surface temperature of ice or ocean
         potT     , & ! air potential temperature  (K)
         Qa       , & ! specific humidity (kg/kg)
         uatm     , & ! x-direction wind speed (m/s)
         vatm     , & ! y-direction wind speed (m/s)
         wind     , & ! wind speed (m/s)
         rhoa         ! air density (kg/m^3)

      real (kind=dbl_kind), intent(in) :: &
         Cdn_atm      ! neutral drag coefficient

      real (kind=dbl_kind), intent(inout):: &
         strx     , & ! x surface stress (N)
         stry         ! y surface stress (N)

      real (kind=dbl_kind), intent(out):: &
         delt     , & ! potential T difference   (K)
         delq     , & ! humidity difference      (kg/kg)
         shcoef   , & ! transfer coefficient for sensible heat
         lhcoef       ! transfer coefficient for latent heat

       ! local variables

      real (kind=dbl_kind) :: &
         TsfK, & ! surface temperature in Kelvin (K)
         qsat, & ! the saturation humidity of air (kg/m^3)
         ssq , & ! sat surface humidity     (kg/kg)
         tau, &  ! stress at zlvl
         Lheat   ! Lvap or Lsub, depending on surface type

      !------------------------------------------------------------
      ! Initialize
      !------------------------------------------------------------

      delt = c0
      delq = c0
      shcoef = c0
      lhcoef = c0
      
      if (calc_strair) then

         strx = c0
         stry = c0

      !------------------------------------------------------------
      ! momentum flux
      !------------------------------------------------------------
         tau = rhoa * 0.0012_dbl_kind * wind
!AOMIP         tau = rhoa * (1.10_dbl_kind + c4*p01*wind) &
!AOMIP                         * wind * p001
         strx = tau * uatm
         stry = tau * vatm

      endif                     ! calc_strair

      !------------------------------------------------------------
      ! define variables that depend on surface type
      !------------------------------------------------------------

      if (sfctype(1:3)=='ice') then
         Lheat = Lsub           ! ice to vapor
      elseif (sfctype(1:3)=='ocn') then
         Lheat = Lvap           ! liquid to vapor
      endif   ! sfctype

      !------------------------------------------------------------
      ! potential temperature and specific humidity differences
      !------------------------------------------------------------

      TsfK     = Tsf + Tffresh    ! surface temp (K)
      qsat     = qqqocn * exp(-TTTocn/TsfK) ! sat humidity (kg/m^3)
      ssq      = qsat / rhoa      ! sat surf hum (kg/kg)
      
      delt= potT - TsfK      ! pot temp diff (K)
      delq= Qa - ssq         ! spec hum dif (kg/kg)
      
      !------------------------------------------------------------
      ! coefficients for turbulent flux calculation
      !------------------------------------------------------------

      shcoef = (1.20e-3_dbl_kind)*cp_air*rhoa*wind
      lhcoef = (1.50e-3_dbl_kind)*Lheat *rhoa*wind

      end subroutine atmo_boundary_const

!=======================================================================

! Neutral drag coefficients for ocean and atmosphere also compute the 
! intermediate necessary variables ridge height, distance, floe size	 
! based upon Tsamados et al. (2014), JPO, DOI: 10.1175/JPO-D-13-0215.1.
! Places where the code varies from the paper are commented.
!
! authors: Michel Tsamados, CPOM
!          David Schroeder, CPOM
!
! changes: Andrew Roberts, NPS (RASM/CESM coupling and documentation)

      subroutine neutral_drag_coeffs (apnd,     hpnd,     &
                                      ipnd,               &
                                      alvl,     vlvl,     &
                                      aice,     vice,     &
                                      vsno,     aicen,    &
                                      vicen,    vsnon,    &
                                      Cdn_ocn,  Cdn_ocn_skin,    &
                                      Cdn_ocn_floe, Cdn_ocn_keel,&
                                      Cdn_atm,  Cdn_atm_skin,    &
                                      Cdn_atm_floe, Cdn_atm_pond,&
                                      Cdn_atm_rdg, hfreebd,      &
                                      hdraft,   hridge,          &
                                      distrdg,  hkeel,           &
                                      dkeel,    lfloe,           &
                                      dfloe,    ncat)

        use ice_colpkg_tracers, only: &
             tr_pond, tr_pond_lvl, tr_pond_topo

      integer (kind=int_kind), intent(in) :: &
         ncat

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         apnd     ,& ! melt pond fraction of sea ice 
         hpnd     ,& ! mean melt pond depth over sea ice 
         ipnd     ,& ! mean ice pond depth over sea ice in cat n
         alvl     ,& ! level ice area fraction (of grid cell ?)
         vlvl        ! level ice mean thickness 
         
      real (kind=dbl_kind), intent(in) :: &
         aice     , & ! concentration of ice
         vice     , & ! volume per unit area of ice
         vsno         ! volume per unit area of snow 
         
      real (kind=dbl_kind), dimension (:), intent(in) :: &
         aicen    , & ! concentration of ice
         vicen    , & ! volume per unit area of ice (m)
         vsnon        ! volume per unit area of snow (m)   
     
      real (kind=dbl_kind), &
         intent(out) :: &
         hfreebd      , & ! freeboard (m)
         hdraft       , & ! draught of ice + snow column (Stoessel1993)
         hridge       , & ! ridge height
         distrdg      , & ! distance between ridges
         hkeel        , & ! keel depth
         dkeel        , & ! distance between keels
         lfloe        , & ! floe length (m)
         dfloe        , & ! distance between floes
         Cdn_ocn      , & ! ocean-ice neutral drag coefficient 
         Cdn_ocn_skin , & ! drag coefficient due to skin drag 
         Cdn_ocn_floe , & ! drag coefficient due to floe edges 
         Cdn_ocn_keel , & ! drag coefficient due to keels 
         Cdn_atm      , & ! ice-atmosphere drag coefficient 
         Cdn_atm_skin , & ! drag coefficient due to skin drag 
         Cdn_atm_floe , & ! drag coefficient due to floe edges 
         Cdn_atm_pond , & ! drag coefficient due to ponds 
         Cdn_atm_rdg      ! drag coefficient due to ridges 

      real (kind=dbl_kind), parameter :: & 
                                      ! [,] = range of values that can be tested 
         csw       = 0.002_dbl_kind ,&! ice-ocn drag coefficient [0.0005,0.005]
         csa       = 0.0005_dbl_kind,&! ice-air drag coefficient [0.0001,0.001] 
         dragia    = 0.0012_dbl_kind,&! ice-air drag coefficient [0.0005,0.002] 
         mrdg      = c20            ,&! screening effect see Lu2011 [5,50]
         mrdgo     = c10            ,&! screening effect see Lu2011 [5,50]
         beta      = p5             ,&! power exponent appearing in astar and 
                                      ! L=Lmin(A*/(A*-A))**beta [0,1]
         Lmin      = c8             ,&! min length of floe (m) [5,100]
         Lmax      = 300._dbl_kind  ,&! max length of floe (m) [30,3000]
         Lmoy      = 300._dbl_kind  ,&! average length of floe (m) [30,1000]
         cfa       = p2             ,&! Eq. 12 ratio of local from drag over 
                                      ! geometrical parameter [0,1] 
         cfw       = p2             ,&! Eq. 15 ratio of local from drag over 
                                      ! geometrical parameter [0,1]
         cpa       = p2             ,&! Eq. 16 ratio of local form drag over 
                                      ! geometrical parameter [0,1]
         cra       = p2             ,&! Eq. 10 local form drag coefficient [0,1]
         crw       = p2             ,&! Eq. 11 local form drag coefficient [0,1]
         sl        = 22._dbl_kind   ,&! Sheltering parameter Lupkes2012 [10,30]
         lpmin     = 2.26_dbl_kind  ,&! min pond length (m) see Eq. 17 [1,10]
         lpmax     = 24.63_dbl_kind ,&! max pond length (m) see Eq. 17 [10,100]
         tanar     = p4             ,&! 0.25 sail slope = 14 deg [0.4,1]
         tanak     = p4             ,&! 0.58 keel slope = 30 deg [0.4,1]
         invsqrte  = 0.6065_dbl_kind,&!
         phir      = 0.8_dbl_kind   ,&! porosity of ridges [0.4,1]
         phik      = 0.8_dbl_kind   ,&! porosity of keels  [0.4,1]
         hkoverhr  = c4             ,&! hkeel/hridge ratio [4,8]
         dkoverdr  = c1             ,&! dkeel/distrdg ratio [1,5]
         sHGB      = 0.18_dbl_kind  ,&! Lupkes2012 Eq. 28, Hanssen1988, 
                                      ! Steele1989 suggest instead 0.18
         alpha2    = c0             ,&! weight functions for area of 
         beta2     = p75              ! ridged ice [0,1]

       integer (kind=int_kind) :: &
         n            ! category index       

      real (kind=dbl_kind) :: &
         astar,     & ! new constant for form drag
         ctecaf,    & ! constante
         ctecwf,    & ! constante
         sca,       & ! wind attenuation function
         scw,       & ! ocean attenuation function	     
         lp,        & ! pond length (m)
         ctecar,    &
         ctecwk,    &
         ai, aii,   & ! ice area and its inverse
         tmp1         ! temporary

      real (kind=dbl_kind) :: &
         apond    , & ! melt pond fraction of grid cell
         vpond    , & ! mean melt pond depth over grid cell
         ipond    , & ! mean melt pond ice depth over grid cell
         ardg     , & ! ridged ice area fraction of grid cell
         vrdg         ! ridged ice mean thickness  

      real (kind=dbl_kind), parameter :: &
         ocnruf   = 0.000327_dbl_kind, & ! ocean surface roughness (m)
         ocnrufi  = c1/ocnruf, & ! inverse ocean roughness
         icerufi  = c1/iceruf    ! inverse ice roughness

      real (kind=dbl_kind), parameter :: &
         camax    = 0.02_dbl_kind , & ! Maximum for atmospheric drag
         cwmax    = 0.06_dbl_kind     ! Maximum for ocean drag

      astar = c1/(c1-(Lmin/Lmax)**(c1/beta))


      !-----------------------------------------------------------------
      ! Initialize across entire grid
      !-----------------------------------------------------------------

      hfreebd=c0
      hdraft =c0       
      hridge =c0       
      distrdg=c0    
      hkeel  =c0 
      dkeel  =c0 
      lfloe  =c0
      dfloe  =c0 
      Cdn_ocn=dragio
      Cdn_ocn_skin=c0 
      Cdn_ocn_floe=c0 
      Cdn_ocn_keel=c0 
      Cdn_atm = (vonkar/log(zref/iceruf)) * (vonkar/log(zref/iceruf))
      Cdn_atm_skin=c0 
      Cdn_atm_floe=c0 
      Cdn_atm_pond=c0 
      Cdn_atm_rdg =c0

      if (aice > p001) then 
      
         Cdn_atm_skin = csa
         Cdn_ocn_skin = csw

         ai  = aice
         aii = c1/ai
         
      !------------------------------------------------------------
      ! Compute average quantities
      !------------------------------------------------------------

         ! ponds
         apond = c0
         vpond = c0
         ipond = c0
         if (tr_pond) then
            do n = 1,ncat
               ! area of pond per unit area of grid cell 
               apond = apond+apnd(n)*aicen(n)  
               ! volume of pond per unit area of grid cell
               vpond = vpond+apnd(n)*hpnd(n)*aicen(n)
            enddo
         endif
         if (tr_pond_lvl .and. tr_pond_topo) then
            do n = 1,ncat
               ! volume of lid per unit area of grid cell
               ipond = ipond+apnd(n)*ipnd(n)*aicen(n)
            enddo
         endif
         
         ! draft and freeboard (see Eq. 27)
         hdraft = (rhoi*vice+rhos*vsno)*aii/rhow ! without ponds
         hfreebd = (vice+vsno)*aii-hdraft
         
         ! Do not allow draft larger than ice thickness (see Eq. 28)  
         if (hdraft >= vice*aii) then
            ! replace excess snow with ice so hi~=hdraft
            hfreebd = (hdraft*ai*(c1-rhoi/rhow) + &
                      (vsno-(vice-hdraft*ai)*rhoi/rhos) * &
                      (c1-rhos/rhow))*aii ! Stoessel1993  
         endif
         
         ! floe size parameterization see Eq. 13
         lfloe = Lmin * (astar / (astar - ai))**beta
         
         ! distance between floes parameterization see Eq. 14
         dfloe = lfloe * (c1/sqrt(ai) - c1)
         
         ! Relate ridge height and distance between ridges to 
         ! ridged ice area fraction and ridged ice mean thickness
         ! Assumes total volume of ridged ice is split into ridges and keels.
         ! Then assume total ridges volume over total area of ridges = 
         ! volume of one average ridge / area of one average ridge
         ! Same for keels.
         
         ardg=c0
         vrdg=c0
         do n=1,ncat
            ! ridged ice area fraction over grid cell 
            ardg=ardg+(c1-alvl(n))*aicen(n)
            ! total ridged ice volume per unit grid cell area
            vrdg=vrdg+(c1-vlvl(n))*vicen(n)
         enddo
         
         ! hridge, hkeel, distrdg and dkeel estimates from CICE for 
         ! simple triangular geometry
         if (ardg > p001) then 
            ! see Eq. 25 and Eq. 26
            hridge = vrdg/ardg*c2 &
                   * (alpha2+beta2*hkoverhr/dkoverdr*tanar/tanak) &
                   / (phir*c1+phik*tanar/tanak*hkoverhr**c2/dkoverdr)
            distrdg = c2*hridge*ai/ardg &
                    * (alpha2/tanar+beta2/tanak*hkoverhr/dkoverdr)
            hkeel = hkoverhr * hridge
            dkeel = dkoverdr * distrdg
            
          ! Use the height of ridges relative to the mean freeboard of
          ! the pack.  Therefore skin drag and ridge drag differ in
          ! this code as compared to  Tsamados et al. (2014) equations
          ! 10 and 18, which reference both to sea level. 
          tmp1 = max(c0,hridge - hfreebd)

      !------------------------------------------------------------
      ! Skin drag (atmo)
      !------------------------------------------------------------	  

          Cdn_atm_skin = csa*(c1 - mrdg*tmp1/distrdg)
          Cdn_atm_skin = max(min(Cdn_atm_skin,camax),c0)

      !------------------------------------------------------------
      ! Ridge effect (atmo)
      !------------------------------------------------------------

          if (tmp1 > puny) then
            sca = c1 - exp(-sHGB*distrdg/tmp1) ! see Eq. 9
            ctecar = cra*p5
            Cdn_atm_rdg = ctecar*tmp1/distrdg*sca* &
                       (log(tmp1*icerufi)/log(zref*icerufi))**c2
            Cdn_atm_rdg = min(Cdn_atm_rdg,camax)
          endif

          ! Use the depth of keels relative to the mean draft of
          ! the pack.  Therefore skin drag and keel drag differ in
          ! this code as compared to  Tsamados et al. (2014) equations
          ! 11 and 19, which reference both to  sea level. In some
          ! circumstances, hkeel can be less than hdraft because hkoverhr
          ! is constant, and max(c0,...) temporarily addresses this.
          tmp1 = max(c0,hkeel - hdraft)

      !------------------------------------------------------------
      ! Skin drag bottom ice (ocean)
      !------------------------------------------------------------	  
  
          Cdn_ocn_skin = csw * (c1 - mrdgo*tmp1/dkeel)
          Cdn_ocn_skin = max(min(Cdn_ocn_skin,cwmax), c0)
  
      !------------------------------------------------------------
      ! Keel effect (ocean)
      !------------------------------------------------------------

          if (tmp1 > puny) then
            scw = c1 - exp(-sHGB*dkeel/tmp1) 
            ctecwk = crw*p5
            Cdn_ocn_keel = ctecwk*tmp1/dkeel*scw* &
                        (log(tmp1*icerufi)/log(zref*icerufi))**c2  
            Cdn_ocn_keel = max(min(Cdn_ocn_keel,cwmax),c0)
          endif
  
         endif ! ardg > 0.001

      !------------------------------------------------------------
      ! Floe edge drag effect (atmo)
      !------------------------------------------------------------

        if (hfreebd > puny) then
          sca = c1 - exp(-sl*beta*(c1-ai))
          ctecaf = cfa*p5*(log(hfreebd*ocnrufi)/log(zref*ocnrufi))**c2*sca
          Cdn_atm_floe = ctecaf * hfreebd / lfloe
          Cdn_atm_floe = max(min(Cdn_atm_floe,camax),c0)
        endif

      !------------------------------------------------------------
      ! Pond edge effect (atmo)
      !------------------------------------------------------------

        if (hfreebd > puny) then
          sca = (apond)**(c1/(zref*beta))
          lp  = lpmin*(1-apond)+lpmax*apond
          Cdn_atm_pond = cpa*p5*sca*apond*hfreebd/lp &
                   * (log(hfreebd*ocnrufi)/log(zref*ocnrufi))**c2
          Cdn_atm_pond = min(Cdn_atm_pond,camax)
        endif
  
      !------------------------------------------------------------
      ! Floe edge drag effect (ocean)
      !------------------------------------------------------------

        if (hdraft > puny) then
          scw = c1 - exp(-sl*beta*(c1-ai))
          ctecwf = cfw*p5*(log(hdraft*ocnrufi)/log(zref*ocnrufi))**c2*scw
          Cdn_ocn_floe = ctecwf * hdraft / lfloe
          Cdn_ocn_floe = max(min(Cdn_ocn_floe,cwmax),c0)
        endif

      !------------------------------------------------------------
      ! Total drag coefficient (atmo)
      !------------------------------------------------------------

         Cdn_atm = Cdn_atm_skin + Cdn_atm_floe + Cdn_atm_pond + Cdn_atm_rdg
         Cdn_atm = min(Cdn_atm,camax)

      !------------------------------------------------------------
      ! Total drag coefficient (ocean)
      !------------------------------------------------------------

         Cdn_ocn = Cdn_ocn_skin + Cdn_ocn_floe + Cdn_ocn_keel
         Cdn_ocn = min(Cdn_ocn,cwmax)  

      endif

      end subroutine neutral_drag_coeffs

!=======================================================================

      end module ice_atmo

!=======================================================================
