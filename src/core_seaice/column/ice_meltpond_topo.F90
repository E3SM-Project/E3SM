!  SVN:$Id: ice_meltpond_topo.F90 1112 2016-03-24 22:49:56Z eclare $
!=======================================================================

! Melt pond evolution based on the ice topography as inferred from
! the ice thickness distribution.  This code is based on (but differs
! from) that described in
!
! Flocco, D. and D. L. Feltham, 2007.  A continuum model of melt pond 
! evolution on Arctic sea ice.  J. Geophys. Res. 112, C08016, doi: 
! 10.1029/2006JC003836.
!
! Flocco, D., D. L. Feltham and A. K. Turner, 2010.  Incorporation of a
! physically based melt pond scheme into the sea ice component of a
! climate model.  J. Geophys. Res. 115, C08012, doi: 10.1029/2009JC005568.
!
! authors Daniela Flocco (UCL)
!         Adrian Turner (UCL)
! 2010 ECH added module based on original code from Daniela Flocco, UCL
! 2012 DSCHR modifications

      module ice_meltpond_topo

      use ice_kinds_mod
      use ice_constants_colpkg, only: c0, c1, c2, p01, p1, p15, p4, p6, &
          puny, viscosity_dyn, rhoi, rhos, rhow, Timelt, Lfresh, &
          gravit, depressT, kice, ice_ref_salinity

      implicit none

      private
      public :: compute_ponds_topo

!=======================================================================

      contains

!=======================================================================

      subroutine compute_ponds_topo(dt,    ncat, nilyr, &
                                    ktherm, heat_capacity, &
                                    aice,  aicen,       &
                                    vice,  vicen,       &
                                    vsno,  vsnon,       &
                                    potT,  meltt,       &
                                    fsurf, fpond,       &
                                    Tsfcn, Tf,          &
                                    qicen, sicen,       &
                                    apnd,  hpnd, ipnd,  &
                                    l_stop,stop_label)

      integer (kind=int_kind), intent(in) :: &
         ncat , &   ! number of thickness categories
         nilyr, &   ! number of ice layers
         ktherm     ! type of thermodynamics (0 0-layer, 1 BL99, 2 mushy)

      logical (kind=log_kind), intent(in) :: &
         heat_capacity   ! if true, ice has nonzero heat capacity
                         ! if false, use zero-layer thermodynamics

      real (kind=dbl_kind), intent(in) :: &
         dt ! time step (s)

      real (kind=dbl_kind), intent(in) :: &
         aice, &    ! total ice area fraction
         vsno, &    ! total snow volume (m)
         Tf   ! ocean freezing temperature [= ice bottom temperature] (degC) 

      real (kind=dbl_kind), intent(inout) :: &
         vice, &    ! total ice volume (m)
         fpond      ! fresh water flux to ponds (m)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         aicen, &   ! ice area fraction, per category
         vsnon      ! snow volume, per category (m)

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         vicen      ! ice volume, per category (m)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         Tsfcn

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         qicen, &
         sicen

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         apnd, &
         hpnd, &
         ipnd

      real (kind=dbl_kind), intent(in) :: &
         potT,  &   ! air potential temperature
         meltt, &   ! total surface meltwater flux
         fsurf      ! thermodynamic heat flux at ice/snow surface (W/m^2)

      logical (kind=log_kind), intent(out) :: &
         l_stop          ! if true, abort model

      character (len=char_len), intent(out) :: &
         stop_label

      ! local variables

      real (kind=dbl_kind), dimension (ncat) :: &
         volpn, & ! pond volume per unit area, per category (m)
         vuin     ! water-equivalent volume of ice lid on melt pond ('upper ice', m) 

      real (kind=dbl_kind), dimension (ncat) :: &
         apondn,& ! pond area fraction, per category
         hpondn   ! pond depth, per category (m)

      real (kind=dbl_kind) :: &
         volp       ! total volume of pond, per unit area of pond (m)

      real (kind=dbl_kind) :: &
         hi,    & ! ice thickness (m)
         dHui,  & ! change in thickness of ice lid (m)
         omega,	& ! conduction
         dTice, & ! temperature difference across ice lid (C)
         dvice, & ! change in ice volume (m)
         Tavg,  & ! mean surface temperature across categories (C)
         Tp,    & ! pond freezing temperature (C)
         dvn      ! change in melt pond volume for fresh water budget

      integer (kind=int_kind) :: n ! loop indices

      real (kind=dbl_kind), parameter :: &
         hicemin = p1           , & ! minimum ice thickness with ponds (m) 
         Td      = p15          , & ! temperature difference for freeze-up (C)
         rhoi_L  = Lfresh * rhoi, & ! (J/m^3)
         min_volp = 1.e-4_dbl_kind  ! minimum pond volume (m)

      !---------------------------------------------------------------
      ! initialize
      !---------------------------------------------------------------
   
      volp = c0

      do n = 1, ncat
         ! load tracers
         volp = volp + hpnd(n) &
              * apnd(n) * aicen(n)
         vuin (n) = ipnd(n) &
                  * apnd(n) * aicen(n)

         hpondn(n) = c0     ! pond depth, per category
         apondn(n) = c0     ! pond area,  per category
      enddo

      ! The freezing temperature for meltponds is assumed slightly below 0C,
      ! as if meltponds had a little salt in them.  The salt budget is not
      ! altered for meltponds, but if it were then an actual pond freezing 
      ! temperature could be computed.

      Tp = Timelt - Td

      !-----------------------------------------------------------------
      ! Identify grid cells with ponds
      !-----------------------------------------------------------------

      hi = c0
      if (aice > puny) hi = vice/aice
      if ( aice > p01 .and. hi > hicemin .and. &
           volp > min_volp*aice) then

         !--------------------------------------------------------------
         ! calculate pond area and depth
         !--------------------------------------------------------------
         call pond_area(dt,         ncat,     nilyr,    &
                        ktherm,     heat_capacity,      &
                        aice,       vice,     vsno,     &
                        aicen,      vicen,    vsnon,    &
                        qicen,      sicen,              &
                        volpn,      volp,               &
                        Tsfcn,      Tf,                 & 
                        apondn,     hpondn,    dvn,     &
                        l_stop,     stop_label)

         fpond = fpond - dvn
         
         ! mean surface temperature
         Tavg = c0
         do n = 1, ncat
            Tavg = Tavg + Tsfcn(n)*aicen(n)
         enddo
         Tavg = Tavg / aice
         
         do n = 1, ncat-1
            
            if (vuin(n) > puny) then
               
         !----------------------------------------------------------------
         ! melting: floating upper ice layer melts in whole or part
         !----------------------------------------------------------------
               ! Use Tsfc for each category
               if (Tsfcn(n) > Tp) then

                  dvice = min(meltt*apondn(n), vuin(n))
                  if (dvice > puny) then
                     vuin (n) = vuin (n) - dvice
                     volpn(n) = volpn(n) + dvice
                     volp     = volp     + dvice
                     fpond    = fpond    + dvice
                     
                     if (vuin(n) < puny .and. volpn(n) > puny) then
                        ! ice lid melted and category is pond covered
                        volpn(n) = volpn(n) + vuin(n)
                        fpond    = fpond    + vuin(n)
                        vuin(n)  = c0
                     endif
                     hpondn(n) = volpn(n) / apondn(n)
                  endif
                  
         !----------------------------------------------------------------
         ! freezing: existing upper ice layer grows
         !----------------------------------------------------------------

               else if (volpn(n) > puny) then ! Tsfcn(i,j,n) <= Tp

                  ! differential growth of base of surface floating ice layer
                  dTice = max(-Tsfcn(n)-Td, c0) ! > 0   
                  omega = kice*DTice/rhoi_L
                  dHui = sqrt(c2*omega*dt + (vuin(n)/aicen(n))**2) &
                                           - vuin(n)/aicen(n)

                  dvice = min(dHui*apondn(n), volpn(n))   
                  if (dvice > puny) then
                     vuin (n) = vuin (n) + dvice
                     volpn(n) = volpn(n) - dvice
                     volp     = volp     - dvice
                     fpond    = fpond    - dvice
                     hpondn(n) = volpn(n) / apondn(n)
                  endif
                  
               endif ! Tsfcn(i,j,n)

         !----------------------------------------------------------------
         ! freezing: upper ice layer begins to form
         ! note: albedo does not change
         !----------------------------------------------------------------
            else ! vuin < puny
                    
               ! thickness of newly formed ice
               ! the surface temperature of a meltpond is the same as that
               ! of the ice underneath (0C), and the thermodynamic surface 
               ! flux is the same
               dHui = max(-fsurf*dt/rhoi_L, c0)
               dvice = min(dHui*apondn(n), volpn(n))  
               if (dvice > puny) then
                  vuin (n) = dvice
                  volpn(n) = volpn(n) - dvice
                  volp     = volp     - dvice
                  fpond    = fpond    - dvice
                  hpondn(n)= volpn(n) / apondn(n)
               endif
               
            endif  ! vuin
            
         enddo ! ncat

      else  ! remove ponds on thin ice
         fpond = fpond - volp
         volpn(:) = c0
         vuin (:) = c0
         volp = c0         
      endif

      !---------------------------------------------------------------
      ! remove ice lid if there is no liquid pond
      ! vuin may be nonzero on category ncat due to dynamics
      !---------------------------------------------------------------

      do n = 1, ncat
         if (aicen(n) > puny .and. volpn(n) < puny &
                             .and. vuin (n) > puny) then
            vuin(n) = c0
         endif

         ! reload tracers
         if (apondn(n) > puny) then
            ipnd(n) = vuin(n) / apondn(n)
         else
            vuin(n) = c0
            ipnd(n) = c0
         endif
         if (aicen(n) > puny) then
            apnd(n) = apondn(n) / aicen(n)
            hpnd(n) = hpondn(n)
         else
            apnd(n) = c0
            hpnd(n) = c0
            ipnd(n) = c0
         endif
      enddo       ! n

 end subroutine compute_ponds_topo

!=======================================================================

! Computes melt pond area, pond depth and melting rates

      subroutine pond_area(dt,    ncat,  nilyr,&
                           ktherm, heat_capacity, &
                           aice,  vice,  vsno, &
                           aicen, vicen, vsnon,& 
                           qicen, sicen,       &
                           volpn, volp,        &
                           Tsfcn,  Tf,         &
                           apondn,hpondn,dvolp,&
                           l_stop,stop_label)

      integer (kind=int_kind), intent(in) :: &
         ncat , & ! number of thickness categories
         nilyr, & ! number of ice layers
         ktherm   ! type of thermodynamics (0 0-layer, 1 BL99, 2 mushy)

      logical (kind=log_kind), intent(in) :: &
         heat_capacity   ! if true, ice has nonzero heat capacity
                         ! if false, use zero-layer thermodynamics

      real (kind=dbl_kind), intent(in) :: &
         dt, aice, vice, vsno, Tf

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         aicen, vicen, vsnon, Tsfcn

      real (kind=dbl_kind), dimension(:,:), intent(in) :: &
         qicen, &
         sicen

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         volpn

      real (kind=dbl_kind), intent(inout) :: &
         volp, dvolp

      real (kind=dbl_kind), dimension(:), intent(out) :: &
         apondn, hpondn

      logical (kind=log_kind), intent(out) :: &
         l_stop          ! if true, abort model

      character (len=char_len), intent(out) :: &
         stop_label

      ! local variables

      integer (kind=int_kind) :: &
         n, ns,   &
	 m_index, &
         permflag

      real (kind=dbl_kind), dimension(ncat) :: &
         hicen, &
         hsnon, &
         asnon, &
         alfan, &
         betan, &
         cum_max_vol, &
         reduced_aicen        

      real (kind=dbl_kind), dimension(0:ncat) :: &
         cum_max_vol_tmp

      real (kind=dbl_kind) :: &
         hpond, &
         drain, &
         floe_weight, &
         pressure_head, &
         hsl_rel, &
         deltah, &
         perm, &
         apond

 !-----------|
 !           |
 !           |-----------|
 !___________|___________|______________________________________sea-level
 !           |           |
 !           |           |---^--------|
 !           |           |   |        |
 !           |           |   |        |-----------|              |-------
 !           |           |   |alfan(n)|           |              |
 !           |           |   |        |           |--------------|
 !           |           |   |        |           |              |
 !---------------------------v-------------------------------------------
 !           |           |   ^        |           |              |
 !           |           |   |        |           |--------------|
 !           |           |   |betan(n)|           |              |
 !           |           |   |        |-----------|              |-------
 !           |           |   |        |
 !           |           |---v------- |
 !           |           |
 !           |-----------|
 !           |
 !-----------|
    
      !-------------------------------------------------------------------
      ! initialize
      !-------------------------------------------------------------------

      do n = 1, ncat

         apondn(n) = c0
         hpondn(n) = c0

         if (aicen(n) < puny)  then
            hicen(n) =  c0 
            hsnon(n) = c0
            reduced_aicen(n) = c0
            asnon(n) = c0
         else
            hicen(n) = vicen(n) / aicen(n)
            hsnon(n) = vsnon(n) / aicen(n)
            reduced_aicen(n) = c1 ! n=ncat
            if (n < ncat) reduced_aicen(n) = aicen(n) &
                * max(0.2_dbl_kind,(-0.024_dbl_kind*hicen(n) + 0.832_dbl_kind))
            asnon(n) = reduced_aicen(n) 
         endif

! This choice for alfa and beta ignores hydrostatic equilibium of categories.
! Hydrostatic equilibium of the entire ITD is accounted for below, assuming
! a surface topography implied by alfa=0.6 and beta=0.4, and rigidity across all
! categories.  alfa and beta partition the ITD - they are areas not thicknesses!
! Multiplying by hicen, alfan and betan (below) are thus volumes per unit area.
! Here, alfa = 60% of the ice area (and since hice is constant in a category, 
! alfan = 60% of the ice volume) in each category lies above the reference line, 
! and 40% below. Note: p6 is an arbitrary choice, but alfa+beta=1 is required.

         alfan(n) = p6 * hicen(n)
         betan(n) = p4 * hicen(n)
       
         cum_max_vol(n)     = c0
         cum_max_vol_tmp(n) = c0
    
      enddo ! ncat

      cum_max_vol_tmp(0) = c0
      drain = c0
      dvolp = c0
    
      !--------------------------------------------------------------------------
      ! the maximum amount of water that can be contained up to each ice category
      !--------------------------------------------------------------------------
    
      do n = 1, ncat-1 ! last category can not hold any volume

         if (alfan(n+1) >= alfan(n) .and. alfan(n+1) > c0) then

            ! total volume in level including snow
            cum_max_vol_tmp(n) = cum_max_vol_tmp(n-1) + &
               (alfan(n+1) - alfan(n)) * sum(reduced_aicen(1:n)) 


            ! subtract snow solid volumes from lower categories in current level
            do ns = 1, n 
               cum_max_vol_tmp(n) = cum_max_vol_tmp(n) &
                  - rhos/rhow  * &    ! fraction of snow that is occupied by solid
                    asnon(ns)  * &    ! area of snow from that category
                    max(min(hsnon(ns)+alfan(ns)-alfan(n), alfan(n+1)-alfan(n)), c0)  
                                      ! thickness of snow from ns layer in n layer
            enddo

         else ! assume higher categories unoccupied
            cum_max_vol_tmp(n) = cum_max_vol_tmp(n-1)
         endif
         if (cum_max_vol_tmp(n) < c0) then
            l_stop = .true.
            stop_label =  'topo ponds: negative melt pond volume'
            return
         endif
      enddo
      cum_max_vol_tmp(ncat) = cum_max_vol_tmp(ncat-1)  ! last category holds no volume
      cum_max_vol  (1:ncat) = cum_max_vol_tmp(1:ncat)
    
      !----------------------------------------------------------------
      ! is there more meltwater than can be held in the floe?
      !----------------------------------------------------------------
      if (volp >= cum_max_vol(ncat)) then
         drain = volp - cum_max_vol(ncat) + puny
         volp = volp - drain
         dvolp = drain
         if (volp < puny) then
            dvolp = dvolp + volp
            volp = c0
         endif
      endif
    
      ! height and area corresponding to the remaining volume

      call calc_hpond(ncat, reduced_aicen, asnon, hsnon, &
                      alfan, volp, cum_max_vol, hpond, m_index)
    
      do n=1, m_index
         hpondn(n) = max((hpond - alfan(n) + alfan(1)), c0)
         apondn(n) = reduced_aicen(n) 
      enddo
      apond = sum(apondn(1:m_index))
    
      !------------------------------------------------------------------------
      ! drainage due to ice permeability - Darcy's law
      !------------------------------------------------------------------------
    
      ! sea water level 
      floe_weight = (vsno*rhos + rhoi*vice + rhow*volp) / aice
      hsl_rel = floe_weight / rhow &
              - ((sum(betan(:)*aicen(:))/aice) + alfan(1))
    
      deltah = hpond - hsl_rel
      pressure_head = gravit * rhow * max(deltah, c0)

      ! drain if ice is permeable    
      permflag = 0
      if (ktherm /= 2 .and. pressure_head > c0) then
      do n = 1, ncat-1
         if (hicen(n) > c0) then
            call permeability_phi(heat_capacity, nilyr, &
                                  qicen(:,n), sicen(:,n), Tsfcn(n), Tf, &
                                  vicen(n),   perm,       l_stop,   stop_label)
            if (l_stop) return
            if (perm > c0) permflag = 1
            drain = perm*apondn(n)*pressure_head*dt / (viscosity_dyn*hicen(n))
            dvolp = dvolp + min(drain, volp)
            volp = max(volp - drain, c0)
            if (volp < puny) then
               dvolp = dvolp + volp
               volp = c0
            endif
         endif
      enddo
 
      ! adjust melt pond dimensions
      if (permflag > 0) then
         ! recompute pond depth    
         call calc_hpond(ncat, reduced_aicen, asnon, hsnon, &
                         alfan, volp, cum_max_vol, hpond, m_index)
         do n=1, m_index
            hpondn(n) = hpond - alfan(n) + alfan(1)
            apondn(n) = reduced_aicen(n) 
         enddo
         apond = sum(apondn(1:m_index))
      endif
      endif ! pressure_head

      !------------------------------------------------------------------------
      ! total melt pond volume in category does not include snow volume
      ! snow in melt ponds is not melted
      !------------------------------------------------------------------------

      ! Calculate pond volume for lower categories
      do n=1,m_index-1
         volpn(n) = apondn(n) * hpondn(n) &
                  - (rhos/rhow) * asnon(n) * min(hsnon(n), hpondn(n))
      enddo

      ! Calculate pond volume for highest category = remaining pond volume
      if (m_index == 1) volpn(m_index) = volp
      if (m_index > 1) then
        if (volp > sum(volpn(1:m_index-1))) then
          volpn(m_index) = volp - sum(volpn(1:m_index-1))
        else
          volpn(m_index) = c0
          hpondn(m_index) = c0
          apondn(m_index) = c0
          ! If remaining pond volume is negative reduce pond volume of 
          ! lower category
          if (volp+puny < sum(volpn(1:m_index-1))) & 
            volpn(m_index-1) = volpn(m_index-1) - sum(volpn(1:m_index-1)) + &
                               volp
        endif
      endif

      do n=1,m_index
         if (apondn(n) > puny) then
             hpondn(n) = volpn(n) / apondn(n)
         else
            dvolp = dvolp + volpn(n)
            hpondn(n) = c0
            volpn(n) = c0
            apondn(n) = c0
         end if
      enddo
      do n = m_index+1, ncat
         hpondn(n) = c0
         apondn(n) = c0
         volpn (n) = c0
      enddo

      end subroutine pond_area
  
!=======================================================================
  
  subroutine calc_hpond(ncat, aicen, asnon, hsnon, &
                        alfan, volp, cum_max_vol, hpond, m_index)
    
      integer (kind=int_kind), intent(in) :: &
         ncat       ! number of thickness categories

    real (kind=dbl_kind), dimension(:), intent(in) :: &
         aicen, &
         asnon, &
         hsnon, &
         alfan, &
         cum_max_vol
    
    real (kind=dbl_kind), intent(in) :: &
         volp
    
    real (kind=dbl_kind), intent(out) :: &
         hpond
    
    integer (kind=int_kind), intent(out) :: &
         m_index
    
    integer :: n, ns
    
    real (kind=dbl_kind), dimension(0:ncat+1) :: &
         hitl, &
         aicetl
    
    real (kind=dbl_kind) :: &
         rem_vol, &
         area, &
         vol, &
         tmp
    
    !----------------------------------------------------------------
    ! hpond is zero if volp is zero - have we fully drained? 
    !----------------------------------------------------------------
    
    if (volp < puny) then
       hpond = c0
       m_index = 0
    else
       
       !----------------------------------------------------------------
       ! Calculate the category where water fills up to 
       !----------------------------------------------------------------
       
       !----------|
       !          |
       !          |
       !          |----------|                                     -- --
       !__________|__________|_________________________________________ ^
       !          |          |             rem_vol     ^                | Semi-filled
       !          |          |----------|-- -- -- - ---|-- ---- -- -- --v layer
       !          |          |          |              |
       !          |          |          |              |hpond
       !          |          |          |----------|   |     |-------
       !          |          |          |          |   |     |
       !          |          |          |          |---v-----|
       !          |          | m_index  |          |         |
       !-------------------------------------------------------------
       
       m_index = 0  ! 1:m_index categories have water in them
       do n = 1, ncat
          if (volp <= cum_max_vol(n)) then
             m_index = n
             if (n == 1) then
                rem_vol = volp
             else
                rem_vol = volp - cum_max_vol(n-1)
             endif
             exit ! to break out of the loop
          endif
       enddo
       m_index = min(ncat-1, m_index)
       
       !----------------------------------------------------------------
       ! semi-filled layer may have m_index different snows in it
       !----------------------------------------------------------------
       
       !-----------------------------------------------------------  ^
       !                                                             |  alfan(m_index+1)
       !                                                             |
       !hitl(3)-->                             |----------|          |
       !hitl(2)-->                |------------| * * * * *|          |
       !hitl(1)-->     |----------|* * * * * * |* * * * * |          |
       !hitl(0)-->-------------------------------------------------  |  ^
       !                various snows from lower categories          |  |alfa(m_index)
       
       ! hitl - heights of the snow layers from thinner and current categories
       ! aicetl - area of each snow depth in this layer
       
       hitl(:) = c0
       aicetl(:) = c0
       do n = 1, m_index
          hitl(n)   = max(min(hsnon(n) + alfan(n) - alfan(m_index), &
                                 alfan(m_index+1) - alfan(m_index)), c0)
          aicetl(n) = asnon(n)
          
          aicetl(0) = aicetl(0) + (aicen(n) - asnon(n))
       enddo
       hitl(m_index+1) = alfan(m_index+1) - alfan(m_index)
       aicetl(m_index+1) = c0
       
       !----------------------------------------------------------------
       ! reorder array according to hitl 
       ! snow heights not necessarily in height order
       !----------------------------------------------------------------
       
       do ns = 1, m_index+1
          do n = 0, m_index - ns + 1
             if (hitl(n) > hitl(n+1)) then ! swap order
                tmp = hitl(n)
                hitl(n) = hitl(n+1)
                hitl(n+1) = tmp
                tmp = aicetl(n)
                aicetl(n) = aicetl(n+1)
                aicetl(n+1) = tmp
             endif
          enddo
       enddo
       
       !----------------------------------------------------------------
       ! divide semi-filled layer into set of sublayers each vertically homogenous
       !----------------------------------------------------------------
       
       !hitl(3)----------------------------------------------------------------
       !                                                       | * * * * * * * *  
       !                                                       |* * * * * * * * * 
       !hitl(2)----------------------------------------------------------------
       !                                    | * * * * * * * *  | * * * * * * * *  
       !                                    |* * * * * * * * * |* * * * * * * * * 
       !hitl(1)----------------------------------------------------------------
       !                 | * * * * * * * *  | * * * * * * * *  | * * * * * * * *  
       !                 |* * * * * * * * * |* * * * * * * * * |* * * * * * * * * 
       !hitl(0)----------------------------------------------------------------
       !    aicetl(0)         aicetl(1)           aicetl(2)          aicetl(3)            
       
       ! move up over layers incrementing volume
       do n = 1, m_index+1
          
          area = sum(aicetl(:)) - &                 ! total area of sub-layer
               (rhos/rhow) * sum(aicetl(n:ncat+1)) ! area of sub-layer occupied by snow
          
          vol = (hitl(n) - hitl(n-1)) * area      ! thickness of sub-layer times area
          
          if (vol >= rem_vol) then  ! have reached the sub-layer with the depth within
             hpond = rem_vol / area + hitl(n-1) + alfan(m_index) - alfan(1)
             exit
          else  ! still in sub-layer below the sub-layer with the depth
             rem_vol = rem_vol - vol
          endif
          
       enddo
       
    endif
    
  end subroutine calc_hpond
  
!=======================================================================

! determine the liquid fraction of brine in the ice and the permeability

      subroutine permeability_phi(heat_capacity, nilyr, &
                                  qicen, sicen, Tsfcn, Tf, &
                                  vicen, perm,  l_stop, stop_label)

      use ice_therm_shared, only: calculate_Tin_from_qin

      logical (kind=log_kind), intent(in) :: &
         heat_capacity   ! if true, ice has nonzero heat capacity
                         ! if false, use zero-layer thermodynamics

      integer (kind=int_kind), intent(in) :: &
         nilyr       ! number of ice layers

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         qicen, &  ! energy of melting for each ice layer (J/m2)
         sicen     ! salinity (ppt)   
    
      real (kind=dbl_kind), intent(in) :: &
         vicen, &  ! ice volume
         Tsfcn, &  ! sea ice surface skin temperature (degC)     
         Tf     ! ocean freezing temperature [= ice bottom temperature] (degC) 
    
      real (kind=dbl_kind), intent(out) :: &
         perm      ! permeability

      logical (kind=log_kind), intent(out) :: &
         l_stop          ! if true, abort model

      character (len=char_len), intent(out) :: &
         stop_label

      ! local variables

      real (kind=dbl_kind) ::   &
         Tmlt, &   ! melting temperature 
         Sbr       ! brine salinity

      real (kind=dbl_kind), dimension(nilyr) ::   &
         Tin, &    ! ice temperature
         phi       ! liquid fraction

      integer (kind=int_kind) :: k
    
      !-----------------------------------------------------------------
      ! Compute ice temperatures from enthalpies using quadratic formula
      ! NOTE this assumes Tmlt = Si * depressT
      !-----------------------------------------------------------------

      if (heat_capacity) then
        do k = 1,nilyr
           Tmlt = -sicen(k) * depressT
           Tin(k) = calculate_Tin_from_qin(qicen(k),Tmlt)
        enddo
      else
        Tin(1) = (Tsfcn + Tf) / c2
      endif  

      !-----------------------------------------------------------------
      ! brine salinity and liquid fraction
      !-----------------------------------------------------------------

      if (maxval(Tin) <= -c2) then

         ! Assur 1958
         do k = 1,nilyr
            Sbr = - 1.2_dbl_kind                 &
                  -21.8_dbl_kind     * Tin(k)    &
                  - 0.919_dbl_kind   * Tin(k)**2 &
                  - 0.01878_dbl_kind * Tin(k)**3
            if (heat_capacity) then
              phi(k) = sicen(k)/Sbr ! liquid fraction
            else
              phi(k) = ice_ref_salinity / Sbr ! liquid fraction
            endif
         enddo ! k
       
      else

         ! Notz 2005 thesis eq. 3.2
         do k = 1,nilyr
            Sbr = -17.6_dbl_kind    * Tin(k)    &
                  - 0.389_dbl_kind  * Tin(k)**2 &
                  - 0.00362_dbl_kind* Tin(k)**3
            if (Sbr == c0) then
               l_stop = .true.
               stop_label =  'topo ponds: zero brine salinity in permeability'
               return
            endif
            if (heat_capacity) then
              phi(k) = sicen(k) / Sbr         ! liquid fraction
            else
              phi(k) = ice_ref_salinity / Sbr ! liquid fraction
            endif

         enddo

      endif

      !-----------------------------------------------------------------
      ! permeability
      !-----------------------------------------------------------------

      perm = 3.0e-08_dbl_kind * (minval(phi))**3
    
      end subroutine permeability_phi

!=======================================================================

      end module ice_meltpond_topo

!=======================================================================
