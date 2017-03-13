!  SVN:$Id: ice_meltpond_lvl.F90 1177 2017-03-08 18:17:21Z eclare $
!=======================================================================

! Level-ice meltpond parameterization
!
! This meltpond parameterization was developed for use with the delta-
! Eddington radiation scheme, and only affects the radiation budget in
! the model.  That is, although the pond volume is tracked, that liquid
! water is not used elsewhere in the model for mass budgets or other
! physical processes.
!
! authors Elizabeth Hunke (LANL)
!         David Hebert (NRL Stennis)
!         Olivier Lecomte (Univ. Louvain)

      module ice_meltpond_lvl

      use ice_kinds_mod
      use ice_constants_colpkg, only: c0, c1, c2, c10, p01, p5, puny, &
          viscosity_dyn, rhoi, rhos, rhow, Timelt, Tffresh, Lfresh, &
          gravit, depressT, rhofresh, kice

      implicit none

      private
      public :: compute_ponds_lvl

!=======================================================================

      contains

!=======================================================================

      subroutine compute_ponds_lvl(dt,     nilyr,        &
                                   ktherm,               &
                                   hi_min, dpscale,      &
                                   frzpnd, pndaspect,    &
                                   rfrac,  meltt, melts, &
                                   frain,  Tair,  fsurfn,&
                                   dhs,    ffrac,        &
                                   aicen,  vicen, vsnon, &
                                   qicen,  sicen,        &
                                   Tsfcn,  alvl,         &
                                   apnd,   hpnd,  ipnd)

      integer (kind=int_kind), intent(in) :: &
         nilyr, &    ! number of ice layers
         ktherm      ! type of thermodynamics (0 0-layer, 1 BL99, 2 mushy)

      real (kind=dbl_kind), intent(in) :: &
         dt,       & ! time step (s)  
         hi_min,   & ! minimum ice thickness allowed for thermo (m)
         dpscale,  & ! alter e-folding time scale for flushing 
         pndaspect   ! ratio of pond depth to pond fraction

      character (len=char_len), intent(in) :: &
         frzpnd      ! pond refreezing parameterization

      real (kind=dbl_kind), &
         intent(in) :: &
         Tsfcn, &    ! surface temperature (C)
         alvl,  &    ! fraction of level ice
         rfrac, &    ! water fraction retained for melt ponds
         meltt, &    ! top melt rate (m/s)
         melts, &    ! snow melt rate (m/s)
         frain, &    ! rainfall rate (kg/m2/s)
         Tair,  &    ! air temperature (K)
         fsurfn,&    ! atm-ice surface heat flux  (W/m2)
         aicen, &    ! ice area fraction
         vicen, &    ! ice volume (m)
         vsnon       ! snow volume (m)

      real (kind=dbl_kind), &
         intent(inout) :: &
         apnd, hpnd, ipnd

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         qicen, &  ! ice layer enthalpy (J m-3)
         sicen     ! salinity (ppt)   

      real (kind=dbl_kind), &
         intent(in) :: &
         dhs       ! depth difference for snow on sea ice and pond ice

      real (kind=dbl_kind), &
         intent(out) :: &
         ffrac     ! fraction of fsurfn over pond used to melt ipond

      ! local temporary variables

      real (kind=dbl_kind) :: &
         volpn     ! pond volume per unit area (m)

      real (kind=dbl_kind), dimension (nilyr) :: &
         Tmlt      ! melting temperature (C)

      real (kind=dbl_kind) :: &
         hi                     , & ! ice thickness (m)
         hs                     , & ! snow depth (m)
         dTs                    , & ! surface temperature diff for freeze-up (C)
         Tp                     , & ! pond freezing temperature (C)
         Ts                     , & ! surface air temperature (C)
         apondn                 , & ! local pond area 
         hpondn                 , & ! local pond depth (m)
         dvn                    , & ! change in pond volume (m)
         hlid, alid             , & ! refrozen lid thickness, area
         dhlid                  , & ! change in refrozen lid thickness
         bdt                    , & ! 2 kice dT dt / (rhoi Lfresh)
         alvl_tmp               , & ! level ice fraction of ice area
         draft, deltah, pressure_head, perm, drain ! for permeability

      real (kind=dbl_kind), parameter :: &
         Td       = c2          , & ! temperature difference for freeze-up (C)
         rexp     = p01             ! pond contraction scaling

      !-----------------------------------------------------------------
      ! Initialize 
      !-----------------------------------------------------------------

      volpn = hpnd * aicen * alvl * apnd
      ffrac = c0

      !-----------------------------------------------------------------
      ! Identify grid cells where ponds can be
      !-----------------------------------------------------------------

      if (aicen*alvl > puny**2) then
         
         hi = vicen/aicen
         hs = vsnon/aicen
         alvl_tmp = alvl

         if (hi < hi_min) then

            !--------------------------------------------------------------
            ! Remove ponds on thin ice
            !--------------------------------------------------------------
            apondn = c0
            hpondn = c0
            volpn  = c0
            hlid = c0

         else

            !-----------------------------------------------------------
            ! initialize pond area as fraction of ice
            !-----------------------------------------------------------
            apondn = apnd*alvl_tmp

            !-----------------------------------------------------------
            ! update pond volume
            !-----------------------------------------------------------
            ! add melt water
            dvn = rfrac/rhofresh*(meltt*rhoi &
                +                 melts*rhos &
                +                 frain*  dt)*aicen

            ! shrink pond volume under freezing conditions
            if (trim(frzpnd) == 'cesm') then
               Tp = Timelt - Td
               dTs = max(Tp - Tsfcn,c0)
               dvn = dvn - volpn * (c1 - exp(rexp*dTs/Tp))

            else 
               ! trim(frzpnd) == 'hlid' Stefan approximation
               ! assumes pond is fresh (freezing temperature = 0 C)
               ! and ice grows from existing pond ice
               hlid = ipnd
               if (dvn == c0) then ! freeze pond
                  Ts = Tair - Tffresh
                  if (Ts < c0) then
                     ! if (Ts < -c2) then ! as in meltpond_cesm
                     bdt = -c2*Ts*kice*dt/(rhoi*Lfresh)
                     dhlid = p5*sqrt(bdt)                  ! open water freezing
                     if (hlid > dhlid) dhlid = p5*bdt/hlid ! existing ice
                     dhlid = min(dhlid, hpnd*rhofresh/rhoi)
                     hlid = hlid + dhlid
                  else
                     dhlid = c0 ! to account for surface inversions
                  endif
               else ! convert refrozen pond ice back to water
                  dhlid = max(fsurfn*dt / (rhoi*Lfresh), c0) ! > 0
                  dhlid = -min(dhlid, hlid) ! < 0
                  hlid = max(hlid + dhlid, c0)
                  if (hs - dhs < puny) then ! pond ice is snow-free
                     ffrac = c1 ! fraction of fsurfn over pond used to melt ipond
                     if (fsurfn > puny) &
                          ffrac = min(-dhlid*rhoi*Lfresh/(dt*fsurfn), c1)
                  endif
               endif
               alid = apondn * aicen
               dvn = dvn - dhlid*alid*rhoi/rhofresh
            endif

            volpn = volpn + dvn

            !-----------------------------------------------------------
            ! update pond area and depth
            !-----------------------------------------------------------
            if (volpn <= c0) then
               volpn = c0
               apondn = c0
            endif

            if (apondn*aicen > puny) then ! existing ponds
               apondn = max(c0, min(alvl_tmp, &
                    apondn + 0.5*dvn/(pndaspect*apondn*aicen)))
               hpondn = c0
               if (apondn > puny) &
                    hpondn = volpn/(apondn*aicen)

            elseif (alvl_tmp*aicen > c10*puny) then ! new ponds
               apondn = min (sqrt(volpn/(pndaspect*aicen)), alvl_tmp)
               hpondn = pndaspect * apondn

            else           ! melt water runs off deformed ice      
               apondn = c0
               hpondn = c0
            endif
            apondn = max(apondn, c0)

            ! limit pond depth to maintain nonnegative freeboard
            hpondn = min(hpondn, ((rhow-rhoi)*hi - rhos*hs)/rhofresh)

            ! fraction of grid cell covered by ponds
            apondn = apondn * aicen

            volpn = hpondn*apondn
            if (volpn <= c0) then
               volpn = c0
               apondn = c0
               hpondn = c0
               hlid = c0
            endif

            !-----------------------------------------------------------
            ! drainage due to permeability (flushing)
            ! setting dpscale = 0 turns this off
            ! NOTE this uses the initial salinity and melting T profiles
            !-----------------------------------------------------------

            if (ktherm /= 2 .and. hpondn > c0 .and. dpscale > puny) then
               draft = (rhos*hs + rhoi*hi)/rhow + hpondn
               deltah = hpondn + hi - draft
               pressure_head = gravit * rhow * max(deltah, c0)
               Tmlt(:) = -sicen(:) * depressT
               call brine_permeability(nilyr, qicen, &
                    vicen, sicen, Tmlt, perm)
               drain = perm*pressure_head*dt / (viscosity_dyn*hi) * dpscale
               deltah = min(drain, hpondn)
               dvn = -deltah*apondn
               volpn = volpn + dvn
               apondn = max(c0, min(apondn &
                    + 0.5*dvn/(pndaspect*apondn), alvl_tmp*aicen))
               hpondn = c0
               if (apondn > puny) hpondn = volpn/apondn
            endif

         endif

         !-----------------------------------------------------------
         ! Reload tracer array
         !-----------------------------------------------------------

         hpnd = hpondn
         apnd = apondn / (aicen*alvl_tmp)
         if (trim(frzpnd) == 'hlid') ipnd = hlid

      endif

      end subroutine compute_ponds_lvl

!=======================================================================

! determine the liquid fraction of brine in the ice and the permeability

      subroutine brine_permeability(nilyr, qicen, vicen, salin, Tmlt, perm)

      use ice_therm_shared, only: calculate_Tin_from_qin

      integer (kind=int_kind), intent(in) :: &
         nilyr     ! number of ice layers

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         qicen, &  ! enthalpy for each ice layer (J m-3)
         salin, &  ! salinity (ppt)   
         Tmlt      ! melting temperature (C)
    
      real (kind=dbl_kind), intent(in) :: &
         vicen     ! ice volume (m)
    
      real (kind=dbl_kind), intent(out) :: &
         perm      ! permeability (m^2)

      ! local variables

      real (kind=dbl_kind) ::   &
         Sbr       ! brine salinity

      real (kind=dbl_kind), dimension(nilyr) ::   &
         Tin, &    ! ice temperature (C)
         phi       ! liquid fraction

      integer (kind=int_kind) :: k
    
      !-----------------------------------------------------------------
      ! Compute ice temperatures from enthalpies using quadratic formula
      !-----------------------------------------------------------------

      do k = 1,nilyr
         Tin(k) = calculate_Tin_from_qin(qicen(k),Tmlt(k))
      enddo

      !-----------------------------------------------------------------
      ! brine salinity and liquid fraction
      !-----------------------------------------------------------------

      do k = 1,nilyr
         Sbr = c1/(1.e-3_dbl_kind - depressT/Tin(k)) ! Notz thesis eq 3.6
         phi(k) = salin(k)/Sbr ! liquid fraction
         if (phi(k) < 0.05) phi(k) = c0 ! impermeable
      enddo

      !-----------------------------------------------------------------
      ! permeability
      !-----------------------------------------------------------------

      perm = 3.0e-8_dbl_kind * (minval(phi))**3
    
      end subroutine brine_permeability
  
!=======================================================================

      end module ice_meltpond_lvl

!=======================================================================
