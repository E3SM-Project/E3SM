!=======================================================================
!BOP
!
! !MODULE: ice_atmo - atm-ice interface: stability based flux calculations
!
! !DESCRIPTION:
!
! Atmospheric boundary interface (stability based flux calculations)
!
! !REVISION HISTORY:
!  SVN:$Id: ice_atmo.F90 49 2007-01-11 22:07:00Z eclare $
!
! author: Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb
! 2004: Block structure added by William Lipscomb
! 2006: Converted to free source form (F90) by Elizabeth Hunke
!
! !INTERFACE:
!
      module ice_atmo
!
! !USES:
!
      use ice_kinds_mod
      use ice_constants
!
!EOP
!
      implicit none
      save

      character (len=char_len) :: &
         atmbndy ! atmo boundary method, 'default' ('ccsm3') or 'constant'

      logical (kind=log_kind) :: &
         calc_strair ! if true, calculate wind stress components

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: atmo_boundary_layer - compute coefficients for atm-ice fluxes, 
!                                  stress and Tref/Qref
!
! !INTERFACE:
!
      subroutine atmo_boundary_layer (nx_block, ny_block, &
                                      sfctype,  icells,   &
                                      indxi,    indxj,    & 
                                      Tsf,      potT,     &
                                      uatm,     vatm,     &  
                                      uvel,     vvel,     &  
                                      wind,     zlvl,     &  
                                      Qa,       rhoa,     &
                                      strx,     stry,     &   
                                      Uref,               &
                                      Tref,     Qref,     &
                                      delt,     delq,     &
                                      lhcoef,   shcoef)

! !DESCRIPTION:
!
! Compute coefficients for atm/ice fluxes, stress, and reference
! temperature and humidity. NOTE: \\
! (1) all fluxes are positive downward,  \\
! (2) here, tstar = (WT)/U*, and qstar = (WQ)/U*,  \\
! (3a) wind speeds should all be above a minimum speed (eg. 1.0 m/s). \\
! (3b) (wind-ice) speeds should all be above a minimum speed (eg. 1.0 m/s). \\
!
! ASSUME:
!  The saturation humidity of air at T(K): qsat(T)  (kg/m**3)
!
! Code originally based on CSM1
!
! !REVISION HISTORY: same as module
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells                ! number of cells that require atmo fluxes

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj    ! compressed i and j indices

      character (len=3), intent(in) :: &
         sfctype      ! ice or ocean

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         Tsf      , & ! surface temperature of ice or ocean
         potT     , & ! air potential temperature  (K)
         uatm     , & ! x-direction wind speed (m/s)
         vatm     , & ! y-direction wind speed (m/s)
         uvel     , & ! x-direction ice speed (m/s)
         vvel     , & ! y-direction ice speed (m/s)
         wind     , & ! wind speed (m/s)
         zlvl     , & ! atm level height (m)
         Qa       , & ! specific humidity (kg/kg)
         rhoa         ! air density (kg/m^3)

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         strx     , & ! x surface stress (N)
         stry         ! y surface stress (N)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out) :: &
         Uref     , & ! reference height wind speed (m/s)
         Tref     , & ! reference height temperature  (K)
         Qref     , & ! reference height specific humidity (kg/kg)
         delt     , & ! potential T difference   (K)
         delq     , & ! humidity difference      (kg/kg)
         shcoef   , & ! transfer coefficient for sensible heat
         lhcoef       ! transfer coefficient for latent heat
!
!EOP
!
       integer (kind=int_kind) :: &
         k     , & ! iteration index
         i, j  , & ! horizontal indices
         ij        ! combined ij index

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
         Lheat     ! Lvap or Lsub, depending on surface type

      real (kind=dbl_kind), dimension (icells) :: &
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
         zTrf  = c2             , & ! reference height for air temp (m)
         umin  = c1                 ! minimum wind speed (m/s)

      ! local functions
      real (kind=dbl_kind) :: &
         xd    , & ! dummy argument
         psimhu, & ! unstable part of psimh
         psixhu    ! unstable part of psimx

      !------------------------------------------------------------
      ! Define functions
      !------------------------------------------------------------

      psimhu(xd)  = log((c1+xd*(c2+xd))*(c1+xd*xd)/c8) &
                  - c2*atan(xd) + pih
!ech                  - c2*atan(xd) + 1.571_dbl_kind

      psixhu(xd)  =  c2 * log((c1 + xd*xd)/c2)

      al2 = log(zref/zTrf)

      !------------------------------------------------------------
      ! Initialize
      !------------------------------------------------------------

      do j = 1, ny_block
      do i = 1, nx_block
         Uref(i,j) = c0
         Tref(i,j) = c0
         Qref(i,j) = c0
         delt(i,j) = c0
         delq(i,j) = c0
         shcoef(i,j) = c0
         lhcoef(i,j) = c0
      enddo
      enddo

      !------------------------------------------------------------
      ! Compute turbulent flux coefficients, wind stress, and
      ! reference temperature and humidity.
      !------------------------------------------------------------

      !------------------------------------------------------------
      ! define variables that depend on surface type
      !------------------------------------------------------------

      if (sfctype(1:3)=='ice') then

         qqq  = qqqice          ! for qsat
         TTT  = TTTice          ! for qsat
         Lheat = Lsub           ! ice to vapor
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            vmag(ij) = max(umin, wind(i,j))
!---------- (3b) option by Andrew Roberts
!            vmag(ij)   = max(umin, sqrt( (uatm(i,j)-uvel(i,j))**2 + (vatm(i,j)-vvel(i,j))**2) )
!---------- (3b) option end
            rdn(ij)  = vonkar/log(zref/iceruf) ! neutral coefficient
         enddo   ! ij

      elseif (sfctype(1:3)=='ocn') then

         qqq  = qqqocn
         TTT  = TTTocn
         Lheat = Lvap           ! liquid to vapor
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
            vmag(ij) = max(umin, wind(i,j))
!---------- (3b) option by Andrew Roberts
!            vmag(ij)   = max(umin, sqrt( (uatm(i,j)-uvel(i,j))**2 + (vatm(i,j)-vvel(i,j))**2) )
!---------- (3b) option end
            rdn(ij)  = sqrt(0.0027_dbl_kind/vmag(ij) &
                    + .000142_dbl_kind + .0000764_dbl_kind*vmag(ij))
         enddo   ! ij

      endif   ! sfctype

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

      !------------------------------------------------------------
      ! define some more needed variables
      !------------------------------------------------------------

         TsfK       = Tsf(i,j) + Tffresh     ! surface temp (K)
         qsat       = qqq * exp(-TTT/TsfK)   ! saturation humidity (kg/m^3)
         ssq        = qsat / rhoa(i,j)       ! sat surf hum (kg/kg)

         thva(ij)   = potT(i,j) * (c1 + zvir * Qa(i,j)) ! virtual pot temp (K)
         delt(i,j)  = potT(i,j) - TsfK       ! pot temp diff (K)
         delq(i,j)  = Qa(i,j) - ssq          ! spec hum dif (kg/kg)
         alz(ij)    = log(zlvl(i,j)/zref)
         cp(ij)     = cp_air*(c1 + cpvir*ssq)

      !------------------------------------------------------------
      ! first estimate of Z/L and ustar, tstar and qstar
      !------------------------------------------------------------

         ! neutral coefficients, z/L = 0.0
         rhn(ij) = rdn(ij)
         ren(ij) = rdn(ij)

         ! ustar,tstar,qstar
         ustar(ij) = rdn(ij) * vmag(ij)
         tstar(ij) = rhn(ij) * delt(i,j)
         qstar(ij) = ren(ij) * delq(i,j)

      enddo                     ! ij

      !------------------------------------------------------------
      ! iterate to converge on Z/L, ustar, tstar and qstar
      !------------------------------------------------------------

      do k=1,5

         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

        ! compute stability & evaluate all stability functions
            hol(ij) = vonkar * gravit * zlvl(i,j) &
                   * (tstar(ij)/thva(ij) &
                    + qstar(ij)/(c1/zvir+Qa(i,j))) &
                   / ustar(ij)**2
            hol(ij)    = sign( min(abs(hol(ij)),c10), hol(ij) )
            stable(ij) = p5 + sign(p5 , hol(ij))
            xqq    = max(sqrt(abs(c1 - c16*hol(ij))) , c1)
            xqq    = sqrt(xqq)

            ! Jordan et al 1999
            psimhs = -(0.7_dbl_kind*hol(ij) &
                   + 0.75_dbl_kind*(hol(ij)-14.3_dbl_kind) &
                   * exp(-0.35_dbl_kind*hol(ij)) + 10.7_dbl_kind)
            psimh  = psimhs*stable(ij) &
                    + (c1 - stable(ij))*psimhu(xqq)
            psixh(ij)  = psimhs*stable(ij) &
                    + (c1 - stable(ij))*psixhu(xqq)

        ! shift all coeffs to measurement height and stability
            rd(ij) = rdn(ij) / (c1+rdn(ij)/vonkar*(alz(ij)-psimh))
            rh(ij) = rhn(ij) / (c1+rhn(ij)/vonkar*(alz(ij)-psixh(ij)))
            re(ij) = ren(ij) / (c1+ren(ij)/vonkar*(alz(ij)-psixh(ij)))

        ! update ustar, tstar, qstar using updated, shifted coeffs
            ustar(ij) = rd(ij) * vmag(ij)
            tstar(ij) = rh(ij) * delt(i,j)
            qstar(ij) = re(ij) * delq(i,j)

         enddo                  ! ij
      enddo                     ! end iteration

      if (calc_strair) then

      ! initialize
      do j = 1, ny_block
      do i = 1, nx_block
         strx(i,j) = c0
         stry(i,j) = c0
      enddo
      enddo

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

      !------------------------------------------------------------
      ! momentum flux
      !------------------------------------------------------------
      ! tau = rhoa(i,j) * ustar * ustar
      ! strx = tau * uatm(i,j) / vmag
      ! stry = tau * vatm(i,j) / vmag
      !------------------------------------------------------------

         tau = rhoa(i,j) * ustar(ij) * rd(ij) ! not the stress at zlvl(i,j)
         strx(i,j) = tau * (uatm(i,j)-uvel(i,j))
         stry(i,j) = tau * (vatm(i,j)-vvel(i,j))

      enddo                     ! ij

      endif                     ! calc_strair

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

      !------------------------------------------------------------
      ! coefficients for turbulent flux calculation
      !------------------------------------------------------------
      ! add windless coefficient for sensible heat flux
      ! as in Jordan et al (JGR, 1999)
      !------------------------------------------------------------

         shcoef(i,j) = rhoa(i,j) * ustar(ij) * cp(ij) * rh(ij) + c1
         lhcoef(i,j) = rhoa(i,j) * ustar(ij) * Lheat  * re(ij)

      !------------------------------------------------------------
      ! Compute diagnostics: 2m ref T & Q and 10m wind speed.
      !------------------------------------------------------------
         hol(ij)  = hol(ij)*zTrf/zlvl(i,j)
         xqq      = max( c1, sqrt(abs(c1-c16*hol(ij))) )
         xqq      = sqrt(xqq)
         psix2    = -c5*hol(ij)*stable(ij) + (c1-stable(ij))*psixhu(xqq)
         fac      = (rh(ij)/vonkar) &
                  * (alz(ij) + al2 - psixh(ij) + psix2)
         Tref(i,j)= potT(i,j) - delt(i,j)*fac
         Tref(i,j)= Tref(i,j) - p01*zTrf ! pot temp to temp correction
         fac      = (re(ij)/vonkar) &
                  * (alz(ij) + al2 - psixh(ij) + psix2)
         Qref(i,j)= Qa(i,j) - delq(i,j)*fac

         Uref(i,j)= vmag(ij) * rd(ij) / rdn(ij)
      enddo                     ! ij

      end subroutine atmo_boundary_layer

!=======================================================================
!BOP
!
! !IROUTINE: atmo_boundary_const - compute coeeficients for atm-ice fluxes
!
!
! !INTERFACE:
!
      subroutine atmo_boundary_const (nx_block, ny_block, &
                                      sfctype,  icells,   &
                                      indxi,    indxj,    & 
                                      uatm,     vatm,     &  
                                      wind,     rhoa,     &
                                      strx,     stry,     &   
                                      lhcoef,   shcoef)

! !DESCRIPTION:
!
! Compute coefficients for atm/ice fluxes, stress
! NOTE: \\
! (1) all fluxes are positive downward,  \\
! (2) reference temperature and humidity are NOT computed
!
! !REVISION HISTORY: same as module
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells                ! number of cells that require atmo fluxes

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj    ! compressed i and j indices

      character (len=3), intent(in) :: &
         sfctype      ! ice or ocean

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         uatm     , & ! x-direction wind speed (m/s)
         vatm     , & ! y-direction wind speed (m/s)
         wind     , & ! wind speed (m/s)
         rhoa         ! air density (kg/m^3)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout):: &
         strx     , & ! x surface stress (N)
         stry         ! y surface stress (N)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out):: &
         shcoef   , & ! transfer coefficient for sensible heat
         lhcoef       ! transfer coefficient for latent heat
!
!EOP
!
       integer (kind=int_kind) :: &
         i, j, & ! horizontal indices
         ij      ! combined ij index

      real (kind=dbl_kind) :: &
         tau, &  ! stress at zlvl
         Lheat   ! Lvap or Lsub, depending on surface type

      !------------------------------------------------------------
      ! Initialize
      !------------------------------------------------------------

      do j = 1, ny_block
      do i = 1, nx_block
         shcoef(i,j) = c0
         lhcoef(i,j) = c0
      enddo
      enddo

      if (calc_strair) then

      do j = 1, ny_block
      do i = 1, nx_block
         strx(i,j) = c0
         stry(i,j) = c0
      enddo
      enddo

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

      !------------------------------------------------------------
      ! momentum flux
      !------------------------------------------------------------

         tau = rhoa(i,j) * 0.0012_dbl_kind * wind(i,j)
!AOMIP         tau = rhoa(i,j) * (1.10_dbl_kind + c4*p01*wind(i,j)) &
!AOMIP                         * wind(i,j) * p001
         strx(i,j) = tau * uatm(i,j)
         stry(i,j) = tau * vatm(i,j)

      enddo                     ! ij

      endif                     ! calc_strair

      !------------------------------------------------------------
      ! define variables that depend on surface type
      !------------------------------------------------------------

      if (sfctype(1:3)=='ice') then
         Lheat = Lsub           ! ice to vapor
      elseif (sfctype(1:3)=='ocn') then
         Lheat = Lvap           ! liquid to vapor
      endif   ! sfctype

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

      !------------------------------------------------------------
      ! coefficients for turbulent flux calculation
      !------------------------------------------------------------

         shcoef(i,j) = (1.20e-3_dbl_kind)*cp_air*rhoa(i,j)*wind(i,j)
         lhcoef(i,j) = (1.50e-3_dbl_kind)*Lheat *rhoa(i,j)*wind(i,j)

      enddo                     ! ij

      end subroutine atmo_boundary_const

!=======================================================================

      end module ice_atmo

!=======================================================================
