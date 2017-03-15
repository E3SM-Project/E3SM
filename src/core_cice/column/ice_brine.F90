!  SVN:$Id: ice_brine.F90 1008 2015-06-20 23:55:12Z eclare $
!=======================================================================
!
! Computes ice microstructural information for use in biogeochemistry
!
! authors: Nicole Jeffery, LANL
!
      module ice_brine

      use ice_kinds_mod
      use ice_constants_colpkg
      use ice_colpkg_tracers, only: ntrcr, nt_qice, nt_sice, nt_bgc_S 
      use ice_zbgc_shared
      use ice_warnings, only: add_warning

      implicit none

      private
      public :: preflushing_changes, compute_microS_mushy, &
                update_hbrine, compute_microS, calculate_drho
 
      real (kind=dbl_kind), parameter :: &   
         maxhbr  = 1.25_dbl_kind  , & ! brine overflows if hbr > maxhbr*hin
         viscos  = 2.1e-6_dbl_kind, & ! kinematic viscosity (m^2/s) 
         ! Brine salinity as a cubic function of temperature
         a1      = -21.4_dbl_kind , & ! (psu/C)  
         a2      = -0.886_dbl_kind, & ! (psu/C^2)
         a3      = -0.012_dbl_kind, & ! (psu/C^3)
         ! Brine density as a quadratic of brine salinity
         b1      = 1000.0_dbl_kind, & ! (kg/m^3)  
         b2      = 0.8_dbl_kind       ! (kg/m^3/ppt)

!=======================================================================

      contains

!=======================================================================
! Computes the top and bottom brine boundary changes for flushing
! works for zsalinity and tr_salinity
!
! NOTE: In this subroutine, trcrn(nt_fbri) is the volume fraction of ice with 
! dynamic salinity or the height ratio = hbr/vicen*aicen, where hbr is the 
! height of the brine surface relative to the bottom of the ice.  This volume fraction
! may be > 1 in which case there is brine above the ice surface (meltponds). 

      subroutine preflushing_changes (n_cat,                &
                                      aicen,    vicen,    vsnon,      &
                                      meltb,    meltt,    congel,     &
                                      snoice,   hice_old, dhice,      &
                                      fbri,     dhbr_top, dhbr_bot,   &
                                      hbr_old,  hin,hsn,  firstice,   &
                                      l_stop,   stop_label)
 
      integer (kind=int_kind), intent(in) :: &
         n_cat           ! category
                              
      real (kind=dbl_kind), intent(in) :: &
         aicen       , & ! concentration of ice
         vicen       , & ! volume per unit area of ice          (m)
         vsnon       , & ! volume per unit area of snow         (m)
         meltb       , & ! bottom ice melt                      (m)
         meltt       , & ! top ice melt                         (m)
         congel      , & ! bottom ice growth                    (m)
         snoice          ! top ice growth from flooding         (m)
  
      real (kind=dbl_kind), intent(out) :: &
         hbr_old          ! old brine height (m)

      real (kind=dbl_kind), intent(inout) :: &
         hin          , & ! ice thickness (m) 
         hsn          , & ! snow thickness (m) 
         dhice            ! change due to sublimation (<0)/condensation (>0) (m)

      real (kind=dbl_kind), intent(inout) :: &
         fbri         , & ! trcrn(nt_fbri)
         dhbr_top     , & ! brine change in top for diagnostics (m)
         dhbr_bot     , & ! brine change in bottom for diagnostics (m)
         hice_old         ! old ice thickness (m)

      logical (kind=log_kind), intent(in) :: &
         firstice         ! if true, initialized values should be used     

      logical (kind=log_kind), intent(out) :: &  
         l_stop            ! if true, abort the model

      character (char_len), intent(out) :: stop_label

      ! local variables

      real (kind=dbl_kind) :: &
         hin_old          ! ice thickness before current melt/growth (m)

      character(len=char_len_long) :: &
         warning ! warning message

      !-----------------------------------------------------------------
      ! initialize
      !-----------------------------------------------------------------

      l_stop = .false.
      if (fbri <= c0) then
         write(warning, *) 'fbri, hice_old', fbri, hice_old
         call add_warning(warning)
         write(warning, *) 'vicen, aicen', vicen, aicen
         call add_warning(warning)         
         l_stop = .true.
         stop_label = 'ice_brine preflushing: fbri <= c0'
      endif

      hin = vicen / aicen
      hsn = vsnon / aicen
      hin_old = max(c0, hin + meltb  + meltt - congel - snoice)
      dhice = hin_old - hice_old   ! change due to subl/cond
      dhbr_top = meltt - snoice - dhice 
      dhbr_bot = congel - meltb

      if ((hice_old < puny) .OR. (hin_old < puny) .OR. firstice) then
         hice_old   = hin
         dhbr_top   = c0
         dhbr_bot   = c0
         dhice       = c0
         fbri       = c1 
      endif

      hbr_old = fbri * hice_old

      end subroutine preflushing_changes

!=======================================================================

! Computes ice microstructural properties for updating hbrine
!
! NOTE: This subroutine uses thermosaline_vertical output to compute
! average ice permeability and the surface ice porosity

      subroutine compute_microS_mushy (n_cat,    nilyr,      nblyr,      &
                                       bgrid,    cgrid,      igrid,      &
                                       trcrn,    hice_old,   hbr_old,    &
                                       sss,      sst,        bTin,       &
                                       iTin,     bphin,                  &
                                       kperm,    bphi_min,   phi_snow,   &
                                       bSin,     brine_sal,  brine_rho,  &
                                       iphin,    ibrine_rho, ibrine_sal, &
                                       sice_rho, iDin,       l_stop,     &
                                       stop_label)  

      use ice_therm_mushy, only: permeability
      use ice_mushy_physics, only: temperature_mush, liquid_fraction
      use ice_colpkg_shared, only: l_sk, min_salin

      integer (kind=int_kind), intent(in) :: &
         n_cat       , & ! ice category
         nilyr       , & ! number of ice layers
         nblyr           ! number of bio layers
                              
      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         bgrid              ! biology nondimensional vertical grid points

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid              ! biology vertical interface points
 
      real (kind=dbl_kind), dimension (nilyr+1), intent(in) :: &
         cgrid              ! CICE vertical coordinate   

      real (kind=dbl_kind), &
         intent(in) :: &
         hice_old    , & ! previous timestep ice height (m)
         phi_snow    , & ! porosity of snow
         sss         , & ! ocean salinity (ppt)
         sst             ! ocean temperature (C)
       
      real (kind=dbl_kind), dimension(ntrcr), &
         intent(in) :: &
         trcrn           

      real (kind=dbl_kind), intent(out) :: & 
         kperm       , & ! average ice permeability (m^2)
         bphi_min        ! surface porosity

      real (kind=dbl_kind), intent(inout) :: &
         hbr_old           ! previous timestep brine height (m)

      real (kind=dbl_kind), dimension (nblyr+1), &
         intent(inout)  :: &
         iDin           ! tracer diffusivity/h^2 (1/s) includes gravity drainage/molecular

      real (kind=dbl_kind), dimension (nblyr+1), &
         intent(inout)  :: &
         iphin       , & ! porosity on the igrid 
         ibrine_rho  , & ! brine rho on interface  
         ibrine_sal  , & ! brine sal on interface   
         iTin            ! Temperature on the igrid (oC)

      real (kind=dbl_kind), dimension (nblyr+2), &
         intent(inout)  :: &
         bSin        , &    ! bulk salinity (ppt) on bgrid
         brine_sal   , & ! equilibrium brine salinity (ppt) 
         brine_rho       ! internal brine density (kg/m^3) 

      real (kind=dbl_kind), dimension (nblyr+2), intent(inout) :: &
         bTin        , & ! Temperature on bgrid
         bphin           ! porosity on bgrid

      real (kind=dbl_kind), intent(inout) :: &
         sice_rho        ! average ice density  

      logical (kind=log_kind), intent(inout) :: &
         l_stop            ! if true, print diagnostics and abort on return
        
      character (char_len), intent(inout) :: stop_label

      ! local variables

      real (kind=dbl_kind), dimension (nilyr) :: &
         cSin        , & ! bulk salinity (ppt)
         cqin            ! enthalpy ()

      real (kind=dbl_kind), dimension (nblyr+2) :: &
         zTin        , & ! Temperature of ice layers on bgrid (C) 
         zSin        , & ! Salinity of ice layers on bgrid (C) 
         bqin            ! enthalpy on the bgrid ()

      real (kind=dbl_kind), dimension (nblyr+1) :: &
         ikin            ! permeability (m^2) 

      integer (kind=int_kind) :: &
         k               ! vertical biology layer index 
      
      real (kind=dbl_kind) :: &
         surface_S   , & ! salinity of ice above hin > hbr 
         hinc_old    , & ! mean ice thickness before current melt/growth (m)
         hbrc_old        ! mean brine thickness before current melt/growth (m)
      
      real (kind=dbl_kind), dimension (ntrcr+2) :: & ! nblyr+2)
         trtmp_s     , & ! temporary, remapped tracers   
         trtmp_q         ! temporary, remapped tracers   
      
      real (kind=dbl_kind), dimension(nblyr+1) :: &   
         drho            ! brine density difference (kg/m^3)
     
      real(kind=dbl_kind), parameter :: &
         Smin = p01
    
      !-----------------------------------------------------------------
      ! Define ice salinity and temperature on bgrid
      !-----------------------------------------------------------------

      trtmp_s(:) = c0
      trtmp_q(:) = c0
      iDin(:) = c0

      do k = 1, nilyr
         cSin(k) = trcrn(nt_sice+k-1)
         cqin(k) = trcrn(nt_qice+k-1)
      enddo

      ! map Sin and qin (cice) profiles to bgc grid 
      surface_S = min_salin
      hbr_old   = min(hbr_old, maxhbr*hice_old)
      hinc_old  = hice_old
      hbrc_old  = hbr_old 

      call remap_zbgc(ntrcr,            nilyr,          &
                      nt_sice,                          &
                      trcrn,            trtmp_s,        &
                      0,                nblyr,          &
                      hinc_old,         hinc_old,       &
                      cgrid(2:nilyr+1),                 &
                      bgrid(2:nblyr+1), surface_S,      &
                      l_stop,           stop_label)
      if (l_stop) return
     
      call remap_zbgc(ntrcr,            nilyr,          &
                      nt_qice,                          &
                      trcrn,            trtmp_q,        &
                      0,                nblyr,          &
                      hinc_old,         hinc_old,       &
                      cgrid(2:nilyr+1),                 &
                      bgrid(2:nblyr+1), surface_S,      &
                      l_stop,           stop_label)
      if (l_stop) return

      do k = 1, nblyr
         bqin (k+1) = min(c0,   trtmp_q(nt_qice+k-1))
         bSin (k+1) = max(Smin, trtmp_s(nt_sice+k-1))
         bTin (k+1) = temperature_mush(bqin(k+1), bSin(k+1))
         bphin(k+1) = liquid_fraction (bTin(k+1), bSin(k+1))
      enddo    ! k

      bSin (1)       = bSin(2)
      bTin (1)       = bTin(2)
      bphin(1)       = bphin(2)
      bphin(nblyr+2) = c1
      bSin (nblyr+2) = sss
      bTin (nblyr+2) = sst
      bphin(nblyr+2) = c1

      !-----------------------------------------------------------------
      ! Define ice multiphase structure
      !-----------------------------------------------------------------
     
      call prepare_hbrine (nblyr, &
                           bSin,          bTin,          iTin,       &
                           brine_sal,     brine_rho,                 &
                           ibrine_sal,    ibrine_rho,                &
                           sice_rho,                                 &
                           bphin,         iphin,                     &
                           kperm,         bphi_min,      phi_snow,   &
                           igrid,         sss)

      call calculate_drho(nblyr, igrid, bgrid,             &
                          brine_rho,    ibrine_rho, drho)   

      do k= 2, nblyr+1
         ikin(k) = k_o*iphin(k)**exp_h 
         iDin(k) = iphin(k)*Dm/hbr_old**2  
         if (hbr_old .GE. Ra_c) &
            iDin(k) = iDin(k) &
                    + l_sk*ikin(k)*gravit/viscos_dynamic*drho(k)/hbr_old**2  
      enddo    ! k

      end subroutine compute_microS_mushy

!=======================================================================

      subroutine prepare_hbrine (nblyr, &
                                 bSin,       bTin,      iTin, &
                                 brine_sal,  brine_rho,       &
                                 ibrine_sal, ibrine_rho,      &
                                 sice_rho,   bphin,     iphin,& 
                                 kperm,      bphi_min,  phi_snow, &
                                 i_grid,     sss)

      use ice_colpkg_shared, only: rhosi
      use ice_therm_shared, only: calculate_Tin_from_qin

      integer (kind=int_kind), intent(in) :: &
         nblyr           ! number of bio layers

      real (kind=dbl_kind), dimension (:), &
         intent(in) :: &
         bSin      , & ! salinity of ice layers on bio grid (ppt)
         bTin      , & ! temperature of ice layers on bio grid for history (C)
         i_grid        ! biology grid interface points

      real (kind=dbl_kind), dimension (:), &
         intent(inout) :: &
         brine_sal  , & ! equilibrium brine salinity (ppt)  
         brine_rho  , & ! internal brine density (kg/m^3)
         ibrine_rho , & ! brine density on interface (kg/m^3)
         ibrine_sal , & ! brine salinity on interface (ppt)
         iphin      , & ! porosity on interface
         iTin       , & ! Temperature on interface 
         bphin          ! porosity of layers

      real (kind=dbl_kind), intent(in) :: &
         phi_snow,    & ! porosity of snow
         sss            ! sea surface salinity (ppt)

      real (kind=dbl_kind), intent(out) :: &
         kperm      , & ! harmonic average permeability (m^2)
         bphi_min       ! minimum porosity

      real (kind=dbl_kind), intent(inout) :: &
         sice_rho       ! avg sea ice density 

      ! local variables

      real (kind=dbl_kind), dimension(nblyr+1) :: &
          kin       !  permeability  
    
      real (kind=dbl_kind) :: &
          k_min, ktemp, &
          igrp, igrm, rigr  ! grid finite differences

      integer (kind=int_kind) :: &
           k   ! layer index

      !-----------------------------------------------------------------
      !  calculate equilibrium brine density and gradients 
      !-----------------------------------------------------------------
    
      sice_rho = c0
      
      do k = 1, nblyr+1
       
         if (k == 1) then 
            igrm = 0
         else
            igrm = i_grid(k) - i_grid(k-1)
         endif

         brine_sal(k) = a1*bTin(k)    &
                      + a2*bTin(k)**2 &
                      + a3*bTin(k)**3
         brine_rho(k) = b1 + b2*brine_sal(k)
         bphin    (k) = max(puny, bSin(k)*rhosi &
                      / (brine_sal(k)*brine_rho(k)))
         bphin    (k) = min(c1, bphin(k))
         kin      (k) = k_o*bphin(k)**exp_h 
         sice_rho     = sice_rho + (rhoi*(c1-bphin(k)) &
                      + brine_rho(k)*bphin(k))*igrm
      enddo    ! k         

      brine_sal (nblyr+2) = sss
      brine_rho (nblyr+2) = rhow
      bphin     (nblyr+2) = c1
      ibrine_sal(1)       = brine_sal (2)
      ibrine_sal(nblyr+1) = brine_sal (nblyr+2)
      ibrine_rho(1)       = brine_rho (2)
      ibrine_rho(nblyr+1) = brine_rho (nblyr+2)
      iTin      (1)       = bTin(2)
      iTin      (nblyr+1) = bTin(nblyr+1)
      iphin     (1)       = bphin     (2)
      iphin     (nblyr+1) = bphin     (nblyr+1)
      k_min               = MINVAL(kin(2:nblyr+1))
      kperm               = c0  ! initialize
      ktemp               = c0
      bphi_min            = bphin     (1)
!     bphi_min            = max(bphin(1),bSin(2)*rhosi/bphin(2) &  
!                        / (brine_sal(1)*brine_rho(1))*phi_snow)

      do k = 2, nblyr
         if (k_min > c0) then
            ktemp = ktemp + c1/kin(k)
            kperm = k_min
         endif

         igrp = i_grid(k+1) - i_grid(k  )
         igrm = i_grid(k  ) - i_grid(k-1)
         rigr = c1 / (i_grid(k+1)-i_grid(k-1))

         ibrine_sal(k) = (brine_sal(k+1)*igrp + brine_sal(k)*igrm) * rigr
         ibrine_rho(k) = (brine_rho(k+1)*igrp + brine_rho(k)*igrm) * rigr
         iTin      (k) = (bTin     (k+1)*igrp + bTin     (k)*igrm) * rigr
         iphin     (k) = max(puny, &
                         (bphin    (k+1)*igrp + bphin    (k)*igrm) * rigr)
         iphin     (k) = min(c1, iphin (k))
      enddo    ! k         

      if (k_min > c0) then
         ktemp = ktemp + c1/kin(nblyr+1) 
         kperm = real(nblyr,kind=dbl_kind)/ktemp
      endif

      end subroutine prepare_hbrine

!=======================================================================

! Changes include brine height increases from ice and snow surface melt, 
! congelation growth, and upward pressure driven flow from snow loading.
!  
! Decreases arise from downward flushing and bottom melt.  
!
! NOTE: In this subroutine, trcrn(nt_fbri) is  the volume fraction of ice 
! with dynamic salinity or the height ratio == hbr/vicen*aicen, where 
! hbr is the height of the brine surface relative to the bottom of the 
! ice.  This volume fraction may be > 1 in which case there is brine 
! above the ice surface (ponds).

      subroutine update_hbrine (meltb,      meltt,       &
                                melts,      dt,          &
                                hin,        hsn,         &
                                hin_old,    hbr,         &
                                hbr_old,    phi_snow,    &
                                fbri,       snoice,      &
                                dhS_top,    dhS_bottom,  &
                                dh_top_chl, dh_bot_chl,  &
                                kperm,      bphi_min,    &
                                darcy_V, darcy_V_chl,    &
                                bphin,      aice0,       &
                                dh_direct)

      use ice_colpkg_shared, only: rhosi

      real (kind=dbl_kind), intent(in) :: &
         dt             ! timestep
           
      real (kind=dbl_kind), intent(in):: &
         meltb,       & ! bottom melt over dt (m)
         meltt,       & ! true top melt over dt (m)
         melts,       & ! true snow melt over dt (m)
         hin,         & ! ice thickness (m)
         hsn,         & ! snow thickness (m)
         hin_old,     & ! past timestep ice thickness (m)
         hbr_old,     & ! previous timestep hbr
         phi_snow,    & ! porosity of snow
         kperm,       & ! avg ice permeability 
         bphin,       & ! upper brine porosity  
         snoice,      & ! snoice change (m)
         dhS_bottom,  & ! change in bottom hbr initially before darcy flow
         aice0          ! open water area fraction

      real (kind=dbl_kind), intent(inout):: &
         darcy_V    , & ! Darcy velocity: m/s
         darcy_V_chl, & ! Darcy velocity: m/s for bgc 
         dhS_top    , & ! change in top hbr before darcy flow
         dh_bot_chl , & ! change in bottom for algae  
         dh_top_chl , & ! change in bottom for algae  
         hbr        , & ! thickness of brine (m) 
         fbri       , & ! brine height ratio tracer (hbr/hin) 
         bphi_min       ! surface porosity   

      real (kind=dbl_kind), intent(out):: &
         dh_direct      ! surface flooding or runoff (m)
 
      ! local variables

      real (kind=dbl_kind) :: &  
         hbrmin    , & ! thinS or hin 
         dhbr_hin   , & ! hbr-hin
         hbrocn     , & ! brine height above sea level (m) hbr-h_ocn
         dhbr       , & ! change in brine surface
         h_ocn      , & ! new ocean surface from ice bottom (m)
         darcy_coeff, & ! magnitude of the Darcy velocity/hbrocn (1/s)
         hbrocn_new , & ! hbrocn after flushing
         dhflood    , & ! surface flooding by ocean
         dhrunoff       ! direct runoff to ocean

      real (kind=dbl_kind), parameter :: &
         dh_min = p001  ! brine remains within dh_min of sea level
                        ! when ice thickness is less than thinS
 
         hbrocn      = c0
         darcy_V     = c0
         darcy_V_chl = c0  
         hbrocn_new  = c0
         h_ocn = rhosi/rhow*hin + rhos/rhow*hsn 
         
         if (hbr_old > thinS .AND. hin_old > thinS .AND. hin > thinS ) then
            hbrmin = thinS
            dhS_top = -max(c0, min(hin_old-hbr_old, meltt)) * rhoi/rhow 
            dhS_top = dhS_top - max(c0, melts) * rhos/rhow
            dh_top_chl = dhS_top
            dhbr    = dhS_bottom - dhS_top  
            hbr     = max(puny, hbr_old+dhbr) 
            hbrocn  = hbr - h_ocn
            darcy_coeff = max(c0, kperm*gravit/(viscos*hbr_old))

            if (hbrocn > c0 .AND. hbr > thinS ) then 
               bphi_min   = bphin  
               dhrunoff  = -dhS_top*aice0
               hbrocn    = max(c0,hbrocn - dhrunoff)
               hbrocn_new = hbrocn*exp(-darcy_coeff/bphi_min*dt)
               hbr = max(hbrmin, h_ocn + hbrocn_new)
               hbrocn_new = hbr-h_ocn
               darcy_V = -SIGN((hbrocn-hbrocn_new)/dt*bphi_min, hbrocn)
               darcy_V_chl= darcy_V 
               dhS_top    = dhS_top - darcy_V*dt/bphi_min + dhrunoff
               dh_top_chl = dh_top_chl - darcy_V_chl*dt/bphi_min + dhrunoff
               dh_direct  = dhrunoff
            elseif (hbrocn < c0 .AND. hbr > thinS) then
               hbrocn_new = hbrocn*exp(-darcy_coeff/bphi_min*dt)
               dhflood  = max(c0,hbrocn_new - hbrocn)*aice0               
               hbr = max(hbrmin, h_ocn + hbrocn_new)
               darcy_V    = -SIGN((hbrocn-hbrocn_new + dhflood)/dt*bphi_min, hbrocn)
               darcy_V_chl= darcy_V 
               dhS_top    = dhS_top - darcy_V*dt/bphi_min - dhflood
               dh_top_chl = dh_top_chl - darcy_V_chl*dt/bphi_min - dhflood
               dh_direct  = -dhflood
            endif
            
            dh_bot_chl = dhS_bottom 
  
         else    ! very thin brine height 
            hbrmin  = min(thinS, hin)
            hbr = max(hbrmin, hbr_old+dhS_bottom-dhS_top)
            dhbr_hin = hbr - h_ocn
            if (abs(dhbr_hin) > dh_min) &
               hbr = max(hbrmin, h_ocn + SIGN(dh_min,dhbr_hin)) 
            dhS_top = hbr_old - hbr + dhS_bottom
            dh_top_chl = dhS_top
            dh_bot_chl = dhS_bottom
         endif 
        
         fbri = hbr/hin 

      end subroutine update_hbrine

!=======================================================================
!
! Computes ice microstructural properties for zbgc and zsalinity 
!
! NOTE: In this subroutine, trcrn(nt_fbri) is  the volume fraction of ice with 
! dynamic salinity or the height ratio == hbr/vicen*aicen, where hbr is the 
! height of the brine surface relative to the bottom of the ice.  
! This volume fraction
! may be > 1 in which case there is brine above the ice surface (meltponds). 
! 
      subroutine compute_microS   (n_cat,      nilyr,    nblyr,      &
                                   bgrid,      cgrid,    igrid,      &
                                   trcrn,      hice_old,             &
                                   hbr_old,    sss,      sst,        &
                                   bTin,       iTin,     bphin,      &
                                   kperm,      bphi_min, phi_snow,   &
                                   Rayleigh_criteria,    firstice,   &
                                   bSin,                 brine_sal,  &
                                   brine_rho,  iphin,    ibrine_rho, &
                                   ibrine_sal, sice_rho, sloss,      &
                                   salinz,     l_stop,   stop_label)
 
      use ice_therm_shared, only: calculate_Tin_from_qin
      use ice_colpkg_tracers, only: nt_fbri, nt_Tsfc
      use ice_colpkg_shared, only: min_salin, rhosi, salt_loss

      integer (kind=int_kind), intent(in) :: &
         n_cat       , & ! ice category
         nilyr       , & ! number of ice layers
         nblyr           ! number of bio layers
                              
      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         bgrid              ! biology nondimensional vertical grid points

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid              ! biology vertical interface points
 
      real (kind=dbl_kind), dimension (nilyr+1), intent(in) :: &
         cgrid              ! CICE vertical coordinate   

      real (kind=dbl_kind), intent(in) :: &
         hice_old      , & ! previous timestep ice height (m)
         phi_snow      , & ! porosity of snow
         sss           , & ! ocean salinity (ppt)
         sst               ! ocean temperature (oC)
 
      real (kind=dbl_kind), dimension(ntrcr), intent(inout) :: &
         trcrn           

      real (kind=dbl_kind), intent(inout) :: &
         hbr_old       , & ! old brine height
         sice_rho          ! average ice density

      real (kind=dbl_kind), dimension (nblyr+2), intent(inout) :: &
         bTin          , & ! Temperature of ice layers on bio grid for history file (^oC) 
         bphin         , & ! Porosity of layers
         brine_sal     , & ! equilibrium brine salinity (ppt) 
         brine_rho         ! Internal brine density (kg/m^3) 

      real (kind=dbl_kind), dimension (nblyr+2), intent(out) :: &
         bSin

      real (kind=dbl_kind), dimension (nblyr+1), intent(out) :: &
         iTin              ! Temperature on the interface grid

      real (kind=dbl_kind), dimension (nilyr), &
         intent(in) :: &
         salinz            ! initial salinity profile for new ice (on cice grid)

      real (kind=dbl_kind), intent(out) :: & 
         bphi_min      , & ! surface porosity
         kperm         , & ! average ice permeability (m^2)
         sloss             ! (g/m^2) salt from brine runoff for hbr > maxhbr*hin

      logical (kind=log_kind), intent(inout) :: &
         Rayleigh_criteria ! .true. if ice exceeded a minimum thickness hin >= Ra_c 
        
      logical (kind=log_kind), intent(in) :: &
         firstice          ! .true. if ice is newly formed

      real (kind=dbl_kind), dimension (nblyr+1), intent(inout)  :: &
         iphin         , & ! porosity on the igrid 
         ibrine_rho    , & ! brine rho on interface  
         ibrine_sal        ! brine sal on interface   

      logical (kind=log_kind), intent(out) :: &  
         l_stop            ! if true, abort the model

      character (char_len), intent(out) :: stop_label

      ! local variables
 
      integer (kind=int_kind) :: &
         k                 ! vertical biology layer index 
      
      real (kind=dbl_kind) :: &
         surface_S     , & ! salinity of ice above hin > hbr 
         hinc_old      , & ! ice thickness (cell quantity) before current melt/growth (m)
         hbrc_old      , & ! brine thickness(cell quantity) before current melt/growth (m)
         h_o               ! freeboard height (m)

      logical (kind=log_kind) :: &
         Rayleigh          ! .true. if ice exceeded a minimum thickness hin >= Ra_c 

      real (kind=dbl_kind), dimension (ntrcr+2) :: &
         trtmp0        , & ! temporary, remapped tracers  
         trtmp             ! temporary, remapped tracers   
     
      real (kind=dbl_kind) :: &
         Tmlts             ! melting temperature

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      l_stop = .false.

      sloss = c0  
      bTin(:) = c0
      bSin(:) = c0

      trtmp(:) = c0    
      surface_S = min_salin   

      hinc_old = hice_old

      !-----------------------------------------------------------------
      ! Rayleigh condition for salinity and bgc:
      ! Implemented as a minimum thickness criteria for category 1 ice only.
      ! When hin >= Ra_c (m),  pressure flow is allowed. 
      ! Turn off by putting Ra_c = 0 in ice_in namelist.
      !-----------------------------------------------------------------

      Rayleigh = .true.
      if (n_cat == 1 .AND. hbr_old < Ra_c) then
         Rayleigh = Rayleigh_criteria ! only category 1 ice can be false 
      endif
                     
      !-----------------------------------------------------------------
      ! Define ice salinity on Sin
      !-----------------------------------------------------------------

         if (firstice) then

            do k = 1, nilyr
               trcrn(nt_sice+k-1) =  sss*salt_loss 
            enddo

            call remap_zbgc(ntrcr,            nilyr,     &
                            nt_sice,                     &
                            trcrn,            trtmp,     &
                            0,                nblyr,     &
                            hinc_old,         hinc_old,  &
                            cgrid(2:nilyr+1),            &
                            bgrid(2:nblyr+1), surface_S, &
                            l_stop,           stop_label)
            if (l_stop) return

            do k = 1, nblyr    
               trcrn(nt_bgc_S+k-1) = max(min_salin,trtmp(nt_sice+k-1)) 
               bSin(k+1) = max(min_salin,trcrn(nt_bgc_S+k-1)) 
               if (trcrn(nt_bgc_S+k-1) < min_salin-puny) l_stop = .true.
            enddo ! k  

            bSin(1) = bSin(2) 
            bSin(nblyr+2) =  sss 

         elseif (hbr_old > maxhbr*hice_old) then

            call remap_zbgc(ntrcr,            nblyr,    &
                            nt_bgc_S,                   &
                            trcrn,            trtmp,    &
                            0,                nblyr,    &
                            hbr_old,                    &
                            maxhbr*hinc_old,            &
                            bgrid(2:nblyr+1),           &
                            bgrid(2:nblyr+1), surface_S,&
                            l_stop,           stop_label)
            if (l_stop) return
      
            do k = 1, nblyr    
               bSin(k+1) = max(min_salin,trtmp(nt_bgc_S+k-1)) 
               sloss = sloss + rhosi*(hbr_old*trcrn(nt_bgc_S+k-1) &
                     - maxhbr*hice_old*bSin(k+1))*(igrid(k+1)-igrid(k))
               trcrn(nt_bgc_S+k-1) = bSin(k+1)                                           
               if (trcrn(nt_bgc_S+k-1) < min_salin-puny) l_stop = .true.
            enddo ! k

            bSin(1) = bSin(2)
            bSin(nblyr+2) =  sss
            hbr_old = maxhbr*hinc_old

         else ! old, thin ice

            do k = 1, nblyr   
               trcrn(nt_bgc_S+k-1) = max(min_salin,trcrn(nt_bgc_S+k-1)) 
               bSin (k+1)     = trcrn(nt_bgc_S+k-1)
            enddo ! k

            bSin (1)       = bSin(2)
            bSin (nblyr+2) = sss

         endif         ! ice type
      
      !-----------------------------------------------------------------
      ! sea ice temperature for bio grid
      !-----------------------------------------------------------------
      
      do k = 1, nilyr
         Tmlts               = -trcrn(nt_sice+k-1)*depressT
         trtmp0(nt_qice+k-1) = calculate_Tin_from_qin(trcrn(nt_qice+k-1),Tmlts)
      enddo   ! k

      trtmp(:) = c0                
      
      ! CICE to Bio: remap temperatures
      call remap_zbgc (ntrcr,                       &
                       nilyr,            nt_qice,   &
                       trtmp0(1:ntrcr),  trtmp,     &
                       0,                nblyr,     &
                       hinc_old,         hbr_old,   &
                       cgrid(2:nilyr+1),            & 
                       bgrid(2:nblyr+1), surface_S, &
                       l_stop,           stop_label)
      if (l_stop) return

      do k = 1, nblyr
         Tmlts          = -bSin(k+1) * depressT
         bTin (k+1)     = min(Tmlts,trtmp(nt_qice+k-1))
      enddo !k

      Tmlts          = -min_salin* depressT 
      bTin (1)       = min(Tmlts,(bTin(2) + trcrn(nt_Tsfc))*p5) 
      Tmlts          = -bSin(nblyr+2)* depressT  
      bTin (nblyr+2) = sst

      !-----------------------------------------------------------------
      ! Define ice multiphase structure
      !-----------------------------------------------------------------
     
      call prepare_hbrine (nblyr, &
                           bSin,          bTin,          iTin,       &
                           brine_sal,     brine_rho,                 &
                           ibrine_sal,    ibrine_rho,                &
                           sice_rho,                                 &
                           bphin,         iphin,                     &
                           kperm,         bphi_min,      phi_snow,   &
                           igrid,         sss)

      if (l_stop) then
         stop_label = 'CICE ice_brine:zsalin < min_salin'
      endif

      end subroutine compute_microS

!==========================================================================================
!
! Find density difference about interface grid points
! for gravity drainage parameterization

      subroutine calculate_drho (nblyr,     i_grid,     b_grid, &
                                 brine_rho, ibrine_rho, drho)

      integer (kind=int_kind), intent(in) :: &
         nblyr          ! number of bio layers

      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         b_grid         ! biology nondimensional grid layer points 

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         i_grid         ! biology grid interface points 

      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         brine_rho     ! Internal brine density (kg/m^3)

      real (kind=dbl_kind), dimension (nblyr + 1), intent(in) :: &
         ibrine_rho    ! Internal brine density (kg/m^3)

      real (kind=dbl_kind), dimension (nblyr+1), intent(out) :: & 
         drho          ! brine difference about grid point (kg/m^3)

      ! local variables

      integer (kind=int_kind) :: &
         k, m, mm ! indices

      integer (kind=int_kind) :: &
         mstop, mstart

      real (kind=dbl_kind), dimension (nblyr + 1) :: &  !on the zbgc vertical grid
         rho_a   , &  ! average brine density  above grid point (kg/m^3)
         rho_2a       ! average brine density  above and below grid points (kg/m^3)

      real (kind=dbl_kind), dimension (nblyr + 1) :: &  !on the zbgc vertical grid
         rho_b   , &  ! brine density  above grid point (kg/m^3)
         rho_2b       ! brine density  above and below grid points (kg/m^3)

       rho_a (:) = c0
       rho_2a(:) = c0
       rho_b (:) = c0
       rho_2b(:) = c0
       drho  (:) = c0 ! surface is snow or atmosphere 

       do k = 1, nblyr+1   ! i_grid values

         !----------------------------------------------
         ! h_avg(k) = i_grid(k)                        
         ! Calculate rho_a(k), ie  average rho above i_grid(k) 
         ! first part is good
         !----------------------------------------------

         if (k == 2) then
            rho_a(2) = (brine_rho(2)*b_grid(2) &
                     + (ibrine_rho(2) + brine_rho(2)) &
                     * p5*(i_grid(2)-b_grid(2)) )/i_grid(2)
            rho_b(2) = brine_rho(2)

         elseif (k > 2 .AND. k < nblyr+1) then 
            rho_a(k) = (rho_a(k-1)*i_grid(k-1)   + (ibrine_rho(k-1) + brine_rho(k)) &
                     * p5*(b_grid(k)-i_grid(k-1)) + (ibrine_rho(k  ) + brine_rho(k)) &
                     * p5*(i_grid(k)-b_grid(k)))/i_grid(k)
            rho_b(k) = brine_rho(k)
         else
            rho_a(nblyr+1) = (rho_a(nblyr)*i_grid(nblyr) + (ibrine_rho(nblyr) + &
                        brine_rho(nblyr+1))*p5*(b_grid(nblyr+1)-i_grid(nblyr)) + &
                        brine_rho(nblyr+1)*(i_grid(nblyr+1)-b_grid(nblyr+1)))/i_grid(nblyr+1)
            rho_a(1) = brine_rho(2)   !for k == 1 use grid point value
            rho_b(nblyr+1) = brine_rho(nblyr+1)
            rho_b(1) =  brine_rho(2)
         endif

     enddo     !k

     !----------------------------------------------
     ! Calculate average above and below k rho_2a
     !----------------------------------------------

     do k = 1, nblyr+1   !i_grid values
        if (k == 1) then
           rho_2a(1) = (rho_a(1)*b_grid(2) + p5*(brine_rho(2) + ibrine_rho(2)) &
                     * (i_grid(2)-b_grid(2)))/i_grid(2)
           rho_2b(1) = brine_rho(2)
        else
           mstop = 2*(k-1) + 1
           if (mstop < nblyr+1) then
              rho_2a(k) = rho_a(mstop)
              mstart = 2
              mstop = 1
           else
              mstart = nblyr+2
              mstop = nblyr+3
           endif                     

           do mm = mstart,mstop
              rho_2a(k) =(rho_a(nblyr+1) + rhow*(c2*i_grid(k)-c1))*p5/i_grid(k)
           enddo
           rho_2b(k) = brine_rho(k+1)
        endif
        drho(k) = max(rho_b(k) - rho_2b(k),max(c0,c2*(rho_a(k)-rho_2a(k)), &
              c2*(brine_rho(k)-brine_rho(k+1))/real(nblyr,kind=dbl_kind)))
     enddo

     end subroutine calculate_drho

!=======================================================================

      end module ice_brine

!=======================================================================
