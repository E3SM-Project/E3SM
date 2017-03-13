!  SVN:$Id: ice_zbgc_shared.F90 1177 2017-03-08 18:17:21Z eclare $
!=======================================================================
!
! Biogeochemistry variables
!
! authors: Nicole Jeffery, LANL
!          Scott Elliot,   LANL
!          Elizabeth C. Hunke, LANL
!
      module ice_zbgc_shared

      use ice_kinds_mod
      use ice_constants_colpkg, only: p01, p1, p5, c0, c1
      use ice_colpkg_shared, only: max_nbtrcr, max_algae, max_doc, &
                                   max_dic, max_aero, max_don, max_fe

      implicit none 

      private
      public :: calculate_qin_from_Sin, remap_zbgc, &
                zap_small_bgc, regrid_stationary

      ! bio parameters for algal_dyn
 
       real (kind=dbl_kind),  dimension(max_algae), public :: &
         R_Si2N     , & ! algal Sil to N (mole/mole) 
         R_S2N      , & ! algal S to N (mole/mole)
         ! Marchetti et al 2006, 3 umol Fe/mol C for iron limited Pseudo-nitzschia
         R_Fe2C     , & ! algal Fe to carbon (umol/mmol)
         R_Fe2N         ! algal Fe to N (umol/mmol)

      real (kind=dbl_kind), dimension(max_don), public :: & 
         R_Fe2DON       ! Fe to N of DON (nmol/umol)

      real (kind=dbl_kind), dimension(max_doc), public :: &  
         R_Fe2DOC       ! Fe to C of DOC (nmol/umol)

      real (kind=dbl_kind), parameter, public :: &
         R_gC2molC  = 12.01_dbl_kind ! mg/mmol C

      ! scavenging coefficient for tracers in snow
      ! bottom to last 6 are from Flanner et al., 2007
      ! very last one is for humic material
      real (kind=dbl_kind), parameter, dimension(max_nbtrcr),  public :: &
         kscavz    = (/ 0.03_dbl_kind, 0.03_dbl_kind, 0.03_dbl_kind, &
                        0.03_dbl_kind, 0.03_dbl_kind, 0.03_dbl_kind, &
                        0.03_dbl_kind, 0.03_dbl_kind, 0.03_dbl_kind, &
                        0.03_dbl_kind, 0.03_dbl_kind, 0.03_dbl_kind, &
                        0.03_dbl_kind, 0.03_dbl_kind, 0.03_dbl_kind, &
                        0.03_dbl_kind, 0.03_dbl_kind, 0.03_dbl_kind, &
                        0.03_dbl_kind, 0.03_dbl_kind, 0.03_dbl_kind, &
                        0.03_dbl_kind, &
                        0.03_dbl_kind, 0.20_dbl_kind, 0.02_dbl_kind, &
                        0.02_dbl_kind, 0.01_dbl_kind, 0.01_dbl_kind, &
                        0.03_dbl_kind /)

      !-----------------------------------------------------------------
      ! skeletal layer biogeochemistry
      !-----------------------------------------------------------------

      real (kind=dbl_kind), parameter, public :: &
         phi_sk     = 0.30_dbl_kind     ! skeletal layer porosity

      !-----------------------------------------------------------------
      ! general biogeochemistry
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension(max_nbtrcr), public :: &
         zbgc_frac_init,&! initializes mobile fraction
         bgc_tracer_type ! described tracer in mobile or stationary phases      
                         ! < 0 is purely mobile (eg. nitrate)
                         ! > 0 has timescales for transitions between 
                         ! phases based on whether the ice is melting or growing

      real (kind=dbl_kind), dimension(max_nbtrcr), public :: & 
         zbgc_init_frac, &   ! fraction of ocean tracer  concentration in new ice
         tau_ret,        &   ! retention timescale  (s), mobile to stationary phase
         tau_rel             ! release timescale    (s), stationary to mobile phase         

      !-----------------------------------------------------------------
      ! From algal_dyn in ice_algae.F90 but not in namelist
      !-----------------------------------------------------------------

      real (kind=dbl_kind), dimension(max_algae), public :: &
         chlabs           , & ! chla absorption 1/m/(mg/m^3)
         alpha2max_low    , & ! light limitation (1/(W/m^2))
         beta2max         , & ! light inhibition (1/(W/m^2))
         mu_max           , & ! maximum growth rate (1/d)
         grow_Tdep        , & ! T dependence of growth (1/C)
         fr_graze         , & ! fraction of algae grazed
         mort_pre         , & ! mortality (1/day)
         mort_Tdep        , & ! T dependence of mortality (1/C)
         k_exude          , & ! algal carbon  exudation rate (1/d)
         K_Nit            , & ! nitrate half saturation (mmol/m^3) 
         K_Am             , & ! ammonium half saturation (mmol/m^3) 
         K_Sil            , & ! silicon half saturation (mmol/m^3)
         K_Fe                 ! iron half saturation  or micromol/m^3
            
      real (kind=dbl_kind), dimension(max_DON), public :: &
         f_don            , & ! fraction of spilled grazing to DON
         kn_bac           , & ! Bacterial degredation of DON (1/d)
         f_don_Am             ! fraction of remineralized DON to Am

      real (kind=dbl_kind), dimension(max_DOC), public :: &
         f_doc            , & ! fraction of mort_N that goes to each doc pool
         f_exude          , & ! fraction of exuded carbon to each DOC pool
         k_bac                ! Bacterial degredation of DOC (1/d)    

      !-----------------------------------------------------------------
      ! brine
      !-----------------------------------------------------------------

      integer (kind=int_kind), parameter, public :: &
         exp_h     = 3              ! power law for hierarchical model  

      real (kind=dbl_kind), parameter, public :: & 
         k_o       = 3.e-8_dbl_kind, & ! permeability scaling factor (m^2)
         thinS     = 0.05_dbl_kind     ! minimum ice thickness for brine

      real (kind=dbl_kind), public :: & 
         flood_frac     ! fraction of ocean/meltwater that floods  !*****

      real (kind=dbl_kind), parameter, public :: & 
         bphimin = 0.03_dbl_kind      ! minimum porosity for zbgc only

!-----------------------------------------------------------------------
! Parameters for zsalinity
!-----------------------------------------------------------------------

      real (kind=dbl_kind), parameter, public :: & 
         viscos_dynamic = 2.2_dbl_kind   , & ! 1.8e-3_dbl_kind (pure water at 0^oC) (kg/m/s)
         Dm             = 1.0e-9_dbl_kind, & ! molecular diffusion (m^2/s)
         Ra_c           = 0.05_dbl_kind      ! critical Rayleigh number for bottom convection

!=======================================================================

      contains

!=======================================================================
! 
! Compute the internal ice enthalpy using new salinity and Tin
!

      function calculate_qin_from_Sin (Tin, Tmltk) &
               result(qin)
            
      use ice_constants_colpkg, only: c1, rhoi, cp_ocn, cp_ice, Lfresh  

      real (kind=dbl_kind), intent(in) :: &
         Tin                ,&  ! internal temperature
         Tmltk                  ! melting temperature at one level

      ! local variables

      real (kind=dbl_kind) :: &
         qin                    ! melting temperature at one level   

      qin =-rhoi*(cp_ice*(Tmltk-Tin) + Lfresh*(c1-Tmltk/Tin) - cp_ocn*Tmltk)

      end function calculate_qin_from_Sin

!=======================================================================
!
! Remaps tracer fields in a given category from one set of layers to another.
! Grids can be very different and  so can  vertical spaces.  

      subroutine remap_zbgc(ntrcr,    nlyrn,    &
                            it,                 &
                            trcrn,    trtmp,    &
                            nr0,      nbyrn,    &
                            hice,     hinS,     &
                            ice_grid, bio_grid, &
                            S_min,    l_stop,   &
                            stop_label)

      integer (kind=int_kind), intent(in) :: &
         ntrcr         , & ! number of tracers in use
         it            , & ! tracer index in top layer
         nr0           , & ! receiver category
         nlyrn         , & ! number of ice layers
         nbyrn             ! number of biology layers

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         trcrn             ! ice tracers

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         trtmp             ! temporary, remapped ice tracers

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         ice_grid          ! CICE grid  cgrid(2:nilyr+1)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         bio_grid          ! CICE grid  grid(2:nbyrn+1)

      real(kind=dbl_kind), intent(in) :: &
         hice          , & ! CICE ice thickness
         hinS          , & ! brine height 
         S_min             ! for salinity on CICE grid        

      logical (kind=log_kind), intent(inout) :: &
         l_stop            ! if true, print diagnostics and abort on return
        
      character (char_len), intent(inout) :: stop_label

      ! local variables

      integer (kind=int_kind) :: &
           kd, kr, kdr , & ! more indices
           kdi         , & ! more indices
           n_nd        , & ! number of layers in donor
           n_nr, n_plus    ! number of layers in receiver

      real (kind=dbl_kind), dimension (nbyrn+3+nlyrn) :: &
           trdr        , & ! combined tracer 
           trgrid          ! combined grid 

      real (kind=dbl_kind), dimension (nbyrn+nlyrn+3) :: &
           tracer      , & ! temporary, ice tracers values
           dgrid       , & ! temporary, donor grid dimensional
           rgrid           ! temporary, receiver grid dimensional

      if ((hinS < c0) .OR. (hice < c0)) then
         l_stop = .true.
         stop_label = 'ice: remap_layers_bgc error'
         return
      endif
         
      if (nr0 == 0) then ! cice to bio

         n_nd            = nlyrn
         n_nr            = nbyrn
         n_plus          = 2
         dgrid (1)       = min(-hice+hinS, -hinS+hice, c0)            
         dgrid (nlyrn+2) = min(hinS, hice) 
	 tracer(1)       = trcrn(it)
	 tracer(nlyrn+2) = trcrn(it+nlyrn-1)
         rgrid (nbyrn+2) = min(hinS, hice)
         if (hice > hinS) then
            rgrid(1) = c0 
	    do kr = 1,n_nr
	       rgrid(kr+1) = bio_grid(kr)*hinS
            enddo
            do kd = 1,n_nd
	       dgrid(kd+1) = (ice_grid(kd)-c1)*hice+hinS
               tracer(kd+1) = trcrn(it+kd-1)
            enddo
         else
            rgrid(1) = -hinS + hice 
            do kr = 1,n_nr
	       rgrid(kr+1) = (bio_grid(kr)-c1)*hinS + hice
            enddo
            do kd = 1,n_nd
	       dgrid(kd+1) = ice_grid(kd)*hice
               tracer(kd+1) = trcrn(it+kd-1)
            enddo
         endif
              
      else               ! bio to cice

         n_nd = nbyrn
         n_nr = nlyrn
         if (hice > hinS) then
            n_plus          = 3
            tracer(1)       = S_min
            tracer(2)       = S_min
            dgrid (1)       = -hice+hinS
            dgrid (2)       = p5*(hinS-hice)
            dgrid (nbyrn+3) = hinS
            tracer(nbyrn+3) = trcrn(it+nbyrn-1)
            rgrid (1)       = -hice + hinS
            rgrid (nlyrn+2) = hinS 
            do kd = 1,n_nd
               dgrid(kd+2) = bio_grid(kd)*hinS
               tracer(kd+2) = trcrn(it+kd-1)
            enddo
            do kr = 1,n_nr
               rgrid(kr+1) = (ice_grid(kr)-c1)*hice+ hinS
            enddo
         else
            n_plus          = 2
            tracer(1)       = trcrn(it)
            tracer(nbyrn+2) = trcrn(it+nbyrn-1)
            dgrid (1)       = hice-hinS
            dgrid (nbyrn+2) = hice
            rgrid (nlyrn+2) = hice
            rgrid (1)       = c0
            do kd = 1,n_nd
              dgrid(kd+1) = (bio_grid(kd)-c1)*hinS + hice
              tracer(kd+1) = trcrn(it+kd-1)
            enddo
            do kr = 1,n_nr
              rgrid(kr+1) = ice_grid(kr)*hice
            enddo
         endif

      endif

      kdr = 0  ! combined indices
      kdi = 1  

      do kr = 1, n_nr
         do kd = kdi, n_nd+n_plus
            if (dgrid(kd) < rgrid(kr+1)) then
               kdr = kdr+1
               trgrid(kdr) = dgrid(kd)
               trdr  (kdr) = tracer(kd)
            elseif (dgrid(kd) > rgrid(kr+1)) then
               kdr = kdr + 1
               kdi = kd
               trgrid(kdr) = rgrid(kr+1)
               trtmp (it+kr-1)  = trdr(kdr-1) &
                           + (rgrid(kr+1) - trgrid(kdr-1)) &
                           * (tracer(kd) - trdr(kdr-1)) &
                           / (dgrid(kd) - trgrid(kdr-1))
               trdr(kdr) = trtmp(it+kr-1) 
               EXIT
            else
               kdr = kdr+1
               kdi = kd+1
               trgrid(kdr) = rgrid(kr+1)
               trtmp (it+kr-1)  = tracer(kd)              
               trdr  (kdr) = tracer(kd)
               EXIT
            endif
         enddo
      enddo

      end subroutine remap_zbgc

!=======================================================================

! remove tracer for very small fractional areas

      subroutine zap_small_bgc (zlevels,  dflux_bio, &
                                dt, zvol, btrcr)

      integer (kind=int_kind), intent(in) :: &
         zlevels    ! number of vertical levels in ice

      real (kind=dbl_kind), intent(in) :: &
         dt         ! time step (s)

      real (kind=dbl_kind), intent(inout) :: &
         dflux_bio  ! zapped bio tracer flux from biology (mmol/m^2/s)

      real (kind=dbl_kind), dimension (zlevels), intent(in) :: &
         btrcr  , & ! zapped bio tracer flux from biology (mmol/m^2/s)
         zvol       ! ice volume (m)

      ! local variables

      integer (kind=int_kind) :: &
         k          ! layer index

      do k = 1, zlevels
         dflux_bio = dflux_bio + btrcr(k)*zvol(k)/dt
      enddo
          
      end subroutine zap_small_bgc

!=======================================================================
!
! authors     Nicole Jeffery, LANL

      subroutine regrid_stationary (C_stationary, hbri_old, &
                                    hbri,         dt,       &
                                    ntrcr,        nblyr,    &
                                    top_conc,     igrid,    &
                                    flux_bio,               &
                                    l_stop,       stop_label, &
                                    melt_b)
      
      use ice_constants_colpkg, only: c0, c1, p5

      integer (kind=int_kind), intent(in) :: &
         ntrcr,         & ! number of tracers
         nblyr            ! number of bio layers

      real (kind=dbl_kind), intent(inout) ::  &
         flux_bio         ! ocean tracer flux (mmol/m^2/s) positive into ocean
 
      real (kind=dbl_kind), dimension (nblyr+1), intent(inout) ::  &     
         C_stationary     ! stationary bulk concentration*h (mmol/m^2)

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid            ! CICE bio grid 
         
      real(kind=dbl_kind),  intent(in) :: &
         dt           , & ! time step
         top_conc     , & ! c0 or frazil concentration
         hbri_old     , & ! previous timestep brine height
         hbri             ! brine height 

      logical (kind=log_kind), intent(inout) :: &
         l_stop            ! if true, print diagnostics and abort on return
        
      character (char_len), intent(inout) :: stop_label

      real(kind=dbl_kind), intent(in), optional :: &
         melt_b           ! bottom melt

      !  local variables

      integer (kind=int_kind) :: k, n, nt

      real (kind=dbl_kind), dimension (ntrcr+2) :: &
         trtmp0,   &    ! temporary, remapped tracers
         trtmp

      real (kind=dbl_kind):: &
         meltb,    &    !
         htemp,    &    ! ice thickness after melt (m)
         zspace,   &    ! bio grid spacing
         sum_old,  &    ! total tracer before melt loss
         sum_new        ! total tracer after melt

      ! initialize

      zspace = c1/(real(nblyr,kind=dbl_kind))
      trtmp0(:) = c0
      trtmp(:) = c0
      meltb = c0
      nt = 1
      if (present(melt_b)) then
         meltb = melt_b
      endif
      if (hbri_old > c0) then
         do k = 1, nblyr+1
            trtmp0(nblyr+2-k) = C_stationary(k)/hbri_old  ! reverse order
         enddo   ! k
      endif
      htemp = max(c0,max(hbri,hbri_old - meltb))

      if (meltb > c0) then

      !-----------------------------------------------------------------
      ! Regrid C_stationary to remove bottom melt
      !-----------------------------------------------------------------
         if (hbri_old > htemp) then
            call remap_zbgc (ntrcr,  nblyr+1,            &
                             nt,                         &
                             trtmp0(1:ntrcr),            &
                             trtmp,                      &
                             0,                nblyr+1,  &
                             hbri_old,         htemp,    &
                             igrid(1:nblyr+1),           &
                             igrid(1:nblyr+1), top_conc, &
                             l_stop,           stop_label)
            if (l_stop) return
    
            sum_new = (trtmp(1)+trtmp(nblyr+1))*htemp*p5*zspace
            sum_old = (C_stationary(1) + C_stationary(nblyr+1))*p5*zspace
            trtmp0(:) = c0
            trtmp0(nblyr+1) = trtmp(nt)
            trtmp0(1) = trtmp(nt + nblyr)
            do k = 2,nblyr
               sum_old = sum_old + C_stationary(k)*zspace
               trtmp0(nblyr+2-k) = trtmp(nt + k-1)
               sum_new = sum_new + trtmp0(nblyr+2-k)*htemp*zspace
            enddo       !k
            flux_bio = flux_bio + max(c0,(sum_old - sum_new)/dt)
            do k = 1, nblyr+1
               C_stationary(k) = trtmp0(k)*hbri
            enddo   ! k
         endif

      elseif (hbri > hbri_old) then

      !-----------------------------------------------------------------
      ! Regrid C_stationary to migrate if there is bottom growth
      !-----------------------------------------------------------------
         call remap_zbgc    (ntrcr,            nblyr+1,  &
                             nt,                         &
                             trtmp0(1:ntrcr),            &
                             trtmp,                      &
                             0,                nblyr+1,  &
                             hbri_old,         hbri,     &
                             igrid(1:nblyr+1),           &
                             igrid(1:nblyr+1), top_conc, &
                             l_stop,           stop_label)
         if (l_stop) return

         trtmp0(:) = c0
         do k = 1,nblyr+1
            trtmp0(nblyr+2-k) = trtmp(nt+k-1)
         enddo    ! k
         flux_bio = flux_bio  - (hbri-hbri_old)*top_conc/dt
         do k = 1, nblyr+1
            C_stationary(k) = trtmp0(k)*hbri
          enddo   ! k
      endif 

      end subroutine regrid_stationary

!=======================================================================

      end module ice_zbgc_shared

!=======================================================================
