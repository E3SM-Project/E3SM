!  SVN:$Id: ice_colpkg.F90 1175 2017-03-02 19:53:26Z akt $
!=========================================================================
!
! flags and interface routines for the column package
!
! authors: Elizabeth C. Hunke, LANL

      module ice_colpkg

      use ice_kinds_mod
      use ice_colpkg_shared ! namelist and other parameters
      use ice_warnings, only: add_warning

      implicit none

      private

      ! initialization
      public :: &
           colpkg_init_itd, &
           colpkg_init_itd_hist, &
           colpkg_init_thermo, &
           colpkg_init_orbit, &
           colpkg_init_trcr, &
           colpkg_init_bgc, &
           colpkg_init_zbgc, &
           colpkg_init_hbrine, &
           colpkg_init_zsalinity, &
           colpkg_init_ocean_conc, &
           colpkg_init_OceanConcArray, &
           colpkg_init_bgc_trcr, &
           colpkg_init_parameters, &
           colpkg_init_tracer_flags, &
           colpkg_init_tracer_indices, &
           colpkg_init_tracer_numbers

      ! time stepping
      public :: &
           colpkg_step_therm1, &
           colpkg_biogeochemistry, &
           colpkg_step_therm2, &
           colpkg_prep_radiation, &
           colpkg_step_radiation, &
           colpkg_step_ridge

      ! other column routines
      public :: &
           colpkg_aggregate, &
           colpkg_ice_strength, &
           colpkg_atm_boundary, &
           colpkg_ocn_mixed_layer

      ! temperature inquiry functions
      public :: &
           colpkg_ice_temperature, &
           colpkg_snow_temperature, &
           colpkg_liquidus_temperature, &
           colpkg_sea_freezing_temperature, &
           colpkg_enthalpy_snow

      ! warning messages
      public :: &
           colpkg_clear_warnings, &
           colpkg_get_warnings, &
           colpkg_print_warnings



!=======================================================================

      contains

!=======================================================================
!     Initialization routines
!=======================================================================

! Initialize area fraction and thickness boundaries for the itd model
!
! authors: William H. Lipscomb and Elizabeth C. Hunke, LANL
!          C. M. Bitz, UW

      subroutine colpkg_init_itd(ncat, hin_max, l_stop, stop_label)

      use ice_colpkg_shared, only: kcatbound, kitd
      use ice_therm_shared, only: hi_min
      use ice_constants_colpkg, only: p01, p1, c0, c1, c2, c3, c15, c25, c100

      integer (kind=int_kind), intent(in) :: &
           ncat ! number of thickness categories

      real (kind=dbl_kind), intent(out) :: &
           hin_max(0:ncat)  ! category limits (m)

      logical (kind=log_kind), intent(inout) :: &
         l_stop          ! if true, print diagnostics and abort model

      character (len=*), intent(out) :: &
         stop_label   ! abort error message

      ! local variables

      integer (kind=int_kind) :: &
           n    ! thickness category index

      real (kind=dbl_kind) :: &
           cc1, cc2, cc3, & ! parameters for kcatbound = 0
           x1           , &
           rn           , & ! real(n)
           rncat        , & ! real(ncat)
           d1           , & ! parameters for kcatbound = 1 (m)
           d2           , & !
           b1           , & ! parameters for kcatbound = 3
           b2           , & !
           b3

      real (kind=dbl_kind), dimension(5) :: wmo5 ! data for wmo itd
      real (kind=dbl_kind), dimension(6) :: wmo6 ! data for wmo itd
      real (kind=dbl_kind), dimension(7) :: wmo7 ! data for wmo itd

      l_stop = .false.

      rncat = real(ncat, kind=dbl_kind)
      d1 = 3.0_dbl_kind / rncat
      d2 = 0.5_dbl_kind / rncat
      b1 = p1         ! asymptotic category width (m)
      b2 = c3         ! thickness for which participation function is small (m)
      b3 = max(rncat*(rncat-1), c2*b2/b1)

      hi_min = p01    ! minimum ice thickness allowed (m) for thermo
                      ! note hi_min is reset to 0.1 for kitd=0, below

      !-----------------------------------------------------------------
      ! Choose category boundaries based on one of four options.
      !
      ! The first formula (kcatbound = 0) was used in Lipscomb (2001) 
      !  and in CICE versions 3.0 and 3.1.
      !
      ! The second formula is more user-friendly in the sense that it
      !  is easy to obtain round numbers for category boundaries:
      !
      !    H(n) = n * [d1 + d2*(n-1)] 
      ! 
      ! Default values are d1 = 300/ncat, d2 = 50/ncat.
      ! For ncat = 5, boundaries in cm are 60, 140, 240, 360, which are 
      !  close to the standard values given by the first formula.
      ! For ncat = 10, boundaries in cm are 30, 70, 120, 180, 250, 330,
      !  420, 520, 630.    
      !
      ! The third option provides support for World Meteorological
      !  Organization classification based on thickness.  The full
      !  WMO thickness distribution is used if ncat = 7;  if ncat=5 
      !  or ncat = 6, some of the thinner categories are combined.
      ! For ncat = 5,  boundaries are         30, 70, 120, 200, >200 cm.
      ! For ncat = 6,  boundaries are     15, 30, 70, 120, 200, >200 cm.
      ! For ncat = 7,  boundaries are 10, 15, 30, 70, 120, 200, >200 cm.
      !
      ! The fourth formula asymptotes to a particular category width as
      ! the number of categories increases, given by the parameter b1.
      ! The parameter b3 is computed so that the category boundaries
      ! are even numbers.
      !
      !    H(n) = b1 * [n + b3*n*(n+1)/(2*N*(N-1))] for N=ncat
      !
      ! kcatbound=-1 is available only for 1-category runs, with
      ! boundaries 0 and 100 m.
      !-----------------------------------------------------------------

      if (kcatbound == -1) then ! single category
         hin_max(0) = c0
         hin_max(1) = c100

      elseif (kcatbound == 0) then   ! original scheme

         if (kitd == 1) then
            ! linear remapping itd category limits
            cc1 = c3/rncat
            cc2 = c15*cc1
            cc3 = c3

            hin_max(0) = c0     ! minimum ice thickness, m
         else
            ! delta function itd category limits
#ifndef CCSMCOUPLED
            hi_min = p1    ! minimum ice thickness allowed (m) for thermo
#endif
            cc1 = max(1.1_dbl_kind/rncat,hi_min)
            cc2 = c25*cc1
            cc3 = 2.25_dbl_kind

            ! hin_max(0) should not be zero
            ! use some caution in making it less than 0.10
            hin_max(0) = hi_min ! minimum ice thickness, m
         endif                  ! kitd

         do n = 1, ncat
            x1 = real(n-1,kind=dbl_kind) / rncat
            hin_max(n) = hin_max(n-1) &
                       + cc1 + cc2*(c1 + tanh(cc3*(x1-c1)))
         enddo

      elseif (kcatbound == 1) then  ! new scheme

         hin_max(0) = c0
         do n = 1, ncat
            rn = real(n, kind=dbl_kind)
            hin_max(n) = rn * (d1 + (rn-c1)*d2)
         enddo

      elseif (kcatbound == 2) then  ! WMO standard

        if (ncat == 5) then
         ! thinnest 3 categories combined
         data wmo5 / 0.30_dbl_kind, 0.70_dbl_kind, &
                    1.20_dbl_kind, 2.00_dbl_kind,  &
                    999._dbl_kind  /
         hin_max(0) = c0
         do n = 1, ncat
            hin_max(n) = wmo5(n)
         enddo
       elseif (ncat == 6) then
         ! thinnest 2 categories combined
         data wmo6 / 0.15_dbl_kind, &
                    0.30_dbl_kind, 0.70_dbl_kind,  &
                    1.20_dbl_kind, 2.00_dbl_kind,  &
                    999._dbl_kind /
!echmod wmo6a
!         data wmo6 /0.30_dbl_kind, 0.70_dbl_kind,  &
!                    1.20_dbl_kind, 2.00_dbl_kind,  &
!                    4.56729_dbl_kind, &
!                    999._dbl_kind /

         hin_max(0) = c0
         do n = 1, ncat
            hin_max(n) = wmo6(n)
         enddo
       elseif (ncat == 7) then
         ! all thickness categories 
         data wmo7 / 0.10_dbl_kind, 0.15_dbl_kind, &
                    0.30_dbl_kind, 0.70_dbl_kind,  &
                    1.20_dbl_kind, 2.00_dbl_kind,  &
                    999._dbl_kind  /
         hin_max(0) = c0
         do n = 1, ncat
            hin_max(n) = wmo7(n)
         enddo
       else
         stop_label = 'kcatbound=2 (WMO) must have ncat=5, 6 or 7'
         l_stop = .true. 
         return
       endif

      elseif (kcatbound == 3) then  ! asymptotic scheme

         hin_max(0) = c0
         do n = 1, ncat
            rn = real(n, kind=dbl_kind)
            hin_max(n) = b1 * (rn + b3*rn*(rn+c1)/(c2*rncat*(rncat-c1)))
         enddo

      endif ! kcatbound

      end subroutine colpkg_init_itd

!=======================================================================

! Initialize area fraction and thickness boundaries for the itd model
!
! authors: William H. Lipscomb and Elizabeth C. Hunke, LANL
!          C. M. Bitz, UW

      subroutine colpkg_init_itd_hist (ncat, hin_max, c_hi_range)

      use ice_colpkg_shared, only: kcatbound, kitd
      use ice_constants_colpkg, only: p01, p1, c2, c3, c15, c25, c100

      integer (kind=int_kind), intent(in) :: &
           ncat ! number of thickness categories

      real (kind=dbl_kind), intent(in) :: &
           hin_max(0:ncat)  ! category limits (m)

      character (len=35), intent(out) :: &
           c_hi_range(ncat) ! string for history output

      ! local variables

      integer (kind=int_kind) :: &
           n    ! thickness category index

      character(len=8) :: c_hinmax1,c_hinmax2
      character(len=2) :: c_nc

      character(len=char_len_long) :: &
           warning ! warning message

         write(warning,*) ' '
         call add_warning(warning)
         write(warning,*) 'hin_max(n-1) < Cat n < hin_max(n)'
         call add_warning(warning)
         do n = 1, ncat
            write(warning,*) hin_max(n-1),' < Cat ',n, ' < ',hin_max(n)
            call add_warning(warning)
            ! Write integer n to character string
            write (c_nc, '(i2)') n    

            ! Write hin_max to character string
            write (c_hinmax1, '(f6.3)') hin_max(n-1)
            write (c_hinmax2, '(f6.3)') hin_max(n)

            ! Save character string to write to history file
            c_hi_range(n)=c_hinmax1//'m < hi Cat '//c_nc//' < '//c_hinmax2//'m'
         enddo

         write(warning,*) ' '
         call add_warning(warning)

      end subroutine colpkg_init_itd_hist

!=======================================================================
!
! Initialize the vertical profile of ice salinity and melting temperature.
!
! authors: C. M. Bitz, UW
!          William H. Lipscomb, LANL

      subroutine colpkg_init_thermo(nilyr, sprofile)

      use ice_colpkg_shared, only: saltmax, ktherm, heat_capacity, &
          min_salin
      use ice_constants_colpkg, only: p5, c0, c1, c2, pi
      use ice_therm_shared, only: l_brine

      integer (kind=int_kind), intent(in) :: &
         nilyr                            ! number of ice layers

      real (kind=dbl_kind), dimension(:), intent(out) :: &
         sprofile                         ! vertical salinity profile

      real (kind=dbl_kind), parameter :: &
         nsal    = 0.407_dbl_kind, &
         msal    = 0.573_dbl_kind

      integer (kind=int_kind) :: k        ! ice layer index
      real (kind=dbl_kind)    :: zn       ! normalized ice thickness

      !-----------------------------------------------------------------
      ! Determine l_brine based on saltmax.
      ! Thermodynamic solver will not converge if l_brine is true and
      !  saltmax is close to zero.
      ! Set l_brine to false for zero layer thermodynamics
      !-----------------------------------------------------------------

      heat_capacity = .true.      
      if (ktherm == 0) heat_capacity = .false. ! 0-layer thermodynamics

      l_brine = .false.
      if (saltmax > min_salin .and. heat_capacity) l_brine = .true.

      !-----------------------------------------------------------------
      ! Prescibe vertical profile of salinity and melting temperature.
      ! Note this profile is only used for BL99 thermodynamics.
      !-----------------------------------------------------------------

      if (l_brine) then
         do k = 1, nilyr
            zn = (real(k,kind=dbl_kind)-p5) /  &
                  real(nilyr,kind=dbl_kind)
            sprofile(k)=(saltmax/c2)*(c1-cos(pi*zn**(nsal/(msal+zn))))
            sprofile(k) = max(sprofile(k), min_salin)
         enddo ! k
         sprofile(nilyr+1) = saltmax

      else ! .not. l_brine
         do k = 1, nilyr+1
            sprofile(k) = c0
         enddo
      endif ! l_brine

      end subroutine colpkg_init_thermo

!=======================================================================

! Compute orbital parameters for the specified date.
!
! author:  Bruce P. Briegleb, NCAR 

      subroutine colpkg_init_orbit(l_stop, stop_label)

      use ice_constants_colpkg, only: iyear_AD, eccen, obliqr, lambm0, &
         mvelpp, obliq, mvelp, decln, eccf, log_print

#ifndef CCSMCOUPLED
      use ice_orbital, only: shr_orb_params
#endif

      logical (kind=log_kind), intent(out) :: &
         l_stop          ! if true, abort the model

      character (len=*), intent(out) :: stop_label

      l_stop = .false.      ! initialized for CCSMCOUPLED
      stop_label = ''       ! initialized for CCSMCOUPLED
      iyear_AD  = 1950
      log_print = .false.   ! if true, write out orbital parameters

#ifndef CCSMCOUPLED
      call shr_orb_params( iyear_AD, eccen , obliq , mvelp    , &
                           obliqr  , lambm0, mvelpp, log_print, &
                           l_stop, stop_label)
#endif

      end subroutine colpkg_init_orbit
 
!=======================================================================

      subroutine colpkg_init_trcr(Tair,     Tf,       &
                                  Sprofile, Tprofile, &
                                  Tsfc,               &
                                  nilyr,    nslyr,    &
                                  qin,      qsn)

      use ice_colpkg_shared, only: calc_Tsfc
      use ice_constants_colpkg, only: Tsmelt, Tffresh, p5, cp_ice, cp_ocn, &
          Lfresh, rhoi, rhos, c0, c1
      use ice_mushy_physics, only: enthalpy_mush

      integer (kind=int_kind), intent(in) :: &
         nilyr, &    ! number of ice layers
         nslyr       ! number of snow layers

      real (kind=dbl_kind), intent(in) :: &
         Tair, &     ! air temperature (C)
         Tf          ! freezing temperature (C)

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         Sprofile, & ! vertical salinity profile (ppt)
         Tprofile    ! vertical temperature profile (C)

      real (kind=dbl_kind), intent(out) :: &
         Tsfc        ! surface temperature (C)

      real (kind=dbl_kind), dimension(:), intent(out) :: &
         qin, &      ! ice enthalpy profile (J/m3)
         qsn         ! snow enthalpy profile (J/m3)

      ! local variables

      integer (kind=int_kind) :: k

      real (kind=dbl_kind) :: &
         slope, Ti

            ! surface temperature
            Tsfc = Tf ! default
            if (calc_Tsfc) Tsfc = min(Tsmelt, Tair - Tffresh) ! deg C

            if (heat_capacity) then

               ! ice enthalpy
               do k = 1, nilyr
                  ! assume linear temp profile and compute enthalpy
                  slope = Tf - Tsfc
                  Ti = Tsfc + slope*(real(k,kind=dbl_kind)-p5) &
                                    /real(nilyr,kind=dbl_kind)
                  if (ktherm == 2) then
                     qin(k) = enthalpy_mush(Ti, Sprofile(k))
                  else
                     qin(k) = -(rhoi * (cp_ice*(Tprofile(k)-Ti) &
                         + Lfresh*(c1-Tprofile(k)/Ti) - cp_ocn*Tprofile(k)))
                  endif
               enddo               ! nilyr

               ! snow enthalpy
               do k = 1, nslyr
                  Ti = min(c0, Tsfc)
                  qsn(k) = -rhos*(Lfresh - cp_ice*Ti)
               enddo               ! nslyr

            else  ! one layer with zero heat capacity

               ! ice energy
               qin(1) = -rhoi * Lfresh 

               ! snow energy
               qsn(1) = -rhos * Lfresh 

            endif               ! heat_capacity

      end subroutine colpkg_init_trcr

!=======================================================================

      subroutine colpkg_init_bgc(dt, ncat, nblyr, nilyr, ntrcr_o, cgrid, igrid, &
         restart_bgc, ntrcr, nbtrcr, sicen, trcrn, &
         sss, nit, amm, sil, dmsp, dms, algalN, &
         doc, don, dic, fed, fep, zaeros, hum,  &
         ocean_bio_all, &
         max_algae, max_doc, max_dic, max_don,  max_fe, max_nbtrcr, max_aero, &
         l_stop, stop_label)

      use ice_constants_colpkg, only: c0, c1, c2, p1, p15, p5
      use ice_zbgc_shared, only: R_S2N, zbgc_frac_init, zbgc_init_frac, remap_zbgc

      ! column package includes
      use ice_colpkg_tracers, only: nt_fbri, nt_bgc_S, nt_sice, nt_zbgc_frac, &
         bio_index_o,  bio_index  
      use ice_colpkg_shared, only: solve_zsal, ktherm, hs_ssl,  &
         skl_bgc, scale_bgc, grid_o_t,  fe_data_type, &
         R_C2N, R_chl2N

      real (kind=dbl_kind), intent(in) :: &
         dt        ! time step

      integer (kind=int_kind), intent(in) :: &
         ncat  , & ! number of thickness categories
         nilyr , & ! number of ice layers
         nblyr , & ! number of bio layers
         ntrcr_o, & ! number of tracers not including bgc
         ntrcr , & ! number of tracers in use
         nbtrcr, & ! number of bio tracers in use
         max_algae, &
         max_doc, &
         max_dic, &
         max_don, &
         max_fe, &
         max_nbtrcr, &
         max_aero
 
      real (kind=dbl_kind), dimension (nblyr+1), intent(inout) :: &
         igrid     ! biology vertical interface points
 
      real (kind=dbl_kind), dimension (nilyr+1), intent(inout) :: &
         cgrid     ! CICE vertical coordinate   

      logical (kind=log_kind), intent(in) :: & 
         restart_bgc ! if .true., read bgc restart file

      real (kind=dbl_kind), dimension(nilyr, ncat), intent(in) :: &
         sicen     ! salinity on the cice grid

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcrn     ! subset of tracer array (only bgc) 

      real (kind=dbl_kind), intent(in) :: &
         sss       ! sea surface salinity (ppt)

      real (kind=dbl_kind), intent(inout) :: &
         nit   , & ! ocean nitrate (mmol/m^3)          
         amm   , & ! ammonia/um (mmol/m^3)
         sil   , & ! silicate (mmol/m^3)
         dmsp  , & ! dmsp (mmol/m^3)
         dms   , & ! dms (mmol/m^3)
         hum       ! hum (mmol/m^3)

      real (kind=dbl_kind), dimension (max_algae), intent(inout) :: &
         algalN    ! ocean algal nitrogen (mmol/m^3) (diatoms, pico, phaeocystis)

      real (kind=dbl_kind), dimension (max_doc), intent(inout) :: &
         doc       ! ocean doc (mmol/m^3)  (proteins, EPS, lipid)

      real (kind=dbl_kind), dimension (max_don), intent(inout) :: &
         don       ! ocean don (mmol/m^3) 

      real (kind=dbl_kind), dimension (max_dic), intent(inout) :: &
         dic       ! ocean dic (mmol/m^3) 

      real (kind=dbl_kind), dimension (max_fe), intent(inout) :: &
         fed, fep  ! ocean disolved and particulate fe (nM) 

      real (kind=dbl_kind), dimension (max_aero), intent(inout) :: &
         zaeros    ! ocean aerosols (mmol/m^3) 

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         ocean_bio_all   ! fixed order, all values even for tracers false

      logical (kind=log_kind), intent(inout) :: &
         l_stop            ! if true, print diagnostics and abort on return

      character (len=*), intent(inout) :: stop_label

      ! local variables

      integer (kind=int_kind) :: &
         k     , & ! vertical index 
         n     , & ! category index 
         mm    , & ! bio tracer index
         ki    , & ! loop index
         ks    , & ! 
         ntrcr_bgc

      real (kind=dbl_kind), dimension (ntrcr+2) :: & 
         trtmp     ! temporary, remapped tracers   
      
      real (kind=dbl_kind), dimension (nblyr+1) :: &
         zspace    ! vertical grid spacing

      real (kind=dbl_kind) :: & 
         dvssl , & ! volume of snow surface layer (m)
         dvint , & ! volume of snow interior      (m)
         nit_dum, &!
         sil_dum

      zspace(:)       = c1/real(nblyr,kind=dbl_kind)
      zspace(1)       = p5*zspace(1)
      zspace(nblyr+1) = p5*zspace(nblyr+1)
      ntrcr_bgc       = ntrcr-ntrcr_o

      call colpkg_init_OceanConcArray(max_nbtrcr,                &
                                 max_algae, max_don,  max_doc,   &
                                 max_dic,   max_aero, max_fe,    &
                                 nit,       amm,      sil,       &
                                 dmsp,      dms,      algalN,    &
                                 doc,       don,      dic,       &  
                                 fed,       fep,      zaeros,    &
                                 ocean_bio_all,       hum)

      if (.not. restart_bgc) then  ! not restarting

      !-----------------------------------------------------------------------------   
      !     Skeletal Layer Model
      !  All bgc tracers are Bulk quantities in units of mmol or mg per m^3
      !  The skeletal layer model assumes a constant 
      !  layer depth (sk_l) and porosity (phi_sk)
      !-----------------------------------------------------------------------------   
         if (skl_bgc) then
       
            do  n = 1,ncat
            do mm = 1,nbtrcr
               ! bulk concentration (mmol or mg per m^3, or 10^-3 mmol/m^3)
               trcrn(bio_index(mm)-ntrcr_o, n) = ocean_bio_all(bio_index_o(mm))
            enddo       ! nbtrcr
            enddo       ! n 

      !-----------------------------------------------------------------------------   
      !    zbgc Model
      !  All bgc tracers are Bulk quantities in units of mmol or mg per m^3
      !  The vertical layer model uses prognosed porosity and layer depth
      !-----------------------------------------------------------------------------   

         else   ! not skl_bgc

            if (scale_bgc .and. solve_zsal) then ! bulk concentration (mmol or mg per m^3)
               do n = 1,ncat
               do mm = 1,nbtrcr
                  do k = 2, nblyr
                     trcrn(bio_index(mm)+k-1-ntrcr_o,n) = &
                          (p5*(trcrn(nt_bgc_S+k-1-ntrcr_o,n)+ trcrn(nt_bgc_S+k-2-ntrcr_o,n)) &
                         / sss*ocean_bio_all(bio_index_o(mm))) 
                  enddo  !k
                  trcrn(nt_zbgc_frac-1+mm-ntrcr_o,n) = zbgc_frac_init(mm)
                  trcrn(bio_index(mm)-ntrcr_o,n) = (trcrn(nt_bgc_S-ntrcr_o,n) &
                                         / sss*ocean_bio_all(bio_index_o(mm))) 
                  trcrn(bio_index(mm)+nblyr-ntrcr_o,n) = (trcrn(nt_bgc_S+nblyr-1-ntrcr_o,n) &
                                               / sss*ocean_bio_all(bio_index_o(mm)))
                  trcrn(bio_index(mm)+nblyr+1-ntrcr_o:bio_index(mm)+nblyr+2-ntrcr_o,n) = c0 ! snow
               enddo ! mm
               enddo ! n 
    
            elseif (scale_bgc .and. ktherm == 2) then
               trtmp(:) = c0
               do n = 1,ncat     
                  call remap_zbgc(nilyr,            nilyr,    &
                                  1,                          &
                                  sicen(:,n),       trtmp,    &
                                  0,                nblyr+1,  &
                                  c1,               c1,       &
                                  cgrid(2:nilyr+1),           &
                                  igrid(1:nblyr+1),           &
                                  sicen(1,n),                 &
                                  l_stop,           stop_label)
                  if (l_stop) return

                  do mm = 1,nbtrcr
                  do k = 1, nblyr + 1            
                     trcrn(bio_index(mm)+k-1-ntrcr_o,n) =   &
                          (trtmp(k)/sss*ocean_bio_all(bio_index_o(mm)))
                     trcrn(bio_index(mm)+nblyr+1-ntrcr_o:bio_index(mm)+nblyr+2-ntrcr_o,n) = c0 ! snow
                  enddo  ! k
                  enddo  ! mm
               enddo     ! n 

            elseif (nbtrcr > 0 .and. nt_fbri > 0) then ! not scale_bgc         
     
               do n = 1,ncat
               do mm = 1,nbtrcr
               do k = 1, nblyr+1
                  trcrn(bio_index(mm)+k-1-ntrcr_o,n) = ocean_bio_all(bio_index_o(mm)) &
                                             * zbgc_init_frac(mm) 
                  trcrn(bio_index(mm)+nblyr+1-ntrcr_o:bio_index(mm)+nblyr+2-ntrcr_o,n) = c0 ! snow
               enddo    ! k
               trcrn(nt_zbgc_frac-1+mm-ntrcr_o,n) = zbgc_frac_init(mm)
               enddo    ! mm
               enddo    ! n 
              
            endif  ! scale_bgc
         endif     ! skl_bgc
      endif        ! restart

      end subroutine colpkg_init_bgc

!=======================================================================

      subroutine colpkg_init_zbgc (nblyr, nilyr, nslyr, &
                 n_algae, n_zaero, n_doc, n_dic, n_don, n_fed, n_fep, &
                 trcr_base, trcr_depend, n_trcr_strata, nt_strata, nbtrcr_sw, &
                 tr_brine, nt_fbri, ntrcr, nbtrcr, nt_bgc_Nit, nt_bgc_Am, &
                 nt_bgc_Sil, nt_bgc_DMS, nt_bgc_PON, nt_bgc_S, nt_bgc_N, &
                 nt_bgc_C, nt_bgc_chl, nt_bgc_DOC, nt_bgc_DON, nt_bgc_DIC, & 
                 nt_zaero, nt_bgc_DMSPp, nt_bgc_DMSPd, nt_bgc_Fed, nt_bgc_Fep, &
                 nt_zbgc_frac, tr_bgc_Nit, tr_bgc_Am, tr_bgc_Sil, tr_bgc_DMS, &
                 tr_bgc_PON, tr_bgc_S, tr_bgc_N, tr_bgc_C, tr_bgc_chl, &
                 tr_bgc_DON, tr_bgc_Fe, tr_zaero, nlt_zaero_sw, nlt_chl_sw, &
                 nlt_bgc_N, nlt_bgc_Nit, nlt_bgc_Am, nlt_bgc_Sil, &
                 nlt_bgc_DMS, nlt_bgc_DMSPp, nlt_bgc_DMSPd, &
                 nlt_bgc_C, nlt_bgc_chl, nlt_bgc_DIC, nlt_bgc_DOC, &
                 nlt_bgc_PON, nlt_bgc_DON, nlt_bgc_Fed, nlt_bgc_Fep, &
                 nlt_zaero, &
                 nt_bgc_hum, nlt_bgc_hum, tr_bgc_hum, solve_zsal, &
                 skl_bgc, z_tracers, dEdd_algae, solve_zbgc, &
                 frazil_scav, initbio_frac, bio_index_o, bio_index, ntrcr_o, &
                 max_algae, max_doc, max_dic, max_don, max_fe, &
                 ratio_Si2N_diatoms, ratio_Si2N_sp, ratio_Si2N_phaeo, &
                 ratio_S2N_diatoms, ratio_S2N_sp, ratio_S2N_phaeo, &
                 ratio_Fe2C_diatoms, ratio_Fe2C_sp, ratio_Fe2C_phaeo, &
                 ratio_Fe2N_diatoms, ratio_Fe2N_sp, ratio_Fe2N_phaeo, &
                 ratio_Fe2DON, ratio_Fe2DOC_s,  ratio_Fe2DOC_l, & 
                 chlabs_diatoms, chlabs_sp, chlabs_phaeo, &    
                 alpha2max_low_diatoms, alpha2max_low_sp, alpha2max_low_phaeo, &  
                 beta2max_diatoms, beta2max_sp, beta2max_phaeo, &    
                 mu_max_diatoms, mu_max_sp, mu_max_phaeo, &      
                 grow_Tdep_diatoms, grow_Tdep_sp, grow_Tdep_phaeo, &      
                 fr_graze_diatoms, fr_graze_sp, fr_graze_phaeo, &    
                 mort_pre_diatoms, mort_pre_sp, mort_pre_phaeo, &        
                 mort_Tdep_diatoms, mort_Tdep_sp, mort_Tdep_phaeo, &
                 k_exude_diatoms, k_exude_sp, k_exude_phaeo, &   
                 K_Nit_diatoms, K_Nit_sp, K_Nit_phaeo, &     
                 K_Am_diatoms, K_Am_sp, K_Am_phaeo, &     
                 K_Sil_diatoms, K_Sil_sp, K_Sil_phaeo, &     
                 K_Fe_diatoms, K_Fe_sp, K_Fe_phaeo, &  
                 f_don_protein, kn_bac_protein, & 
                 f_don_Am_protein ,f_doc_s, f_doc_l, &
                 f_exude_s, f_exude_l, k_bac_s,  k_bac_l, &
                 algaltype_diatoms, algaltype_sp, algaltype_phaeo, &
                 doctype_s, doctype_l, dontype_protein, &
                 fedtype_1, feptype_1, zaerotype_bc1, zaerotype_bc2, &
                 zaerotype_dust1, zaerotype_dust2, zaerotype_dust3, &
                 zaerotype_dust4, &
                 ratio_C2N_diatoms, ratio_C2N_sp, ratio_C2N_phaeo, &
                 ratio_chl2N_diatoms, ratio_chl2N_sp, ratio_chl2N_phaeo, &
                 F_abs_chl_diatoms, F_abs_chl_sp, F_abs_chl_phaeo, &
                 ratio_C2N_proteins, &
                 nitratetype, ammoniumtype, dmspptype, dmspdtype, &
                 silicatetype, humtype, tau_min, tau_max)
                    
      use ice_constants_colpkg, only: c1, p5, c0, c2

      use ice_colpkg_shared, only: &
         algaltype, doctype, dictype, dontype, fedtype, feptype, zaerotype, &
         R_C2N, R_chl2N, F_abs_chl, R_C2N_DON

      use ice_zbgc_shared, only: zbgc_init_frac, &
         bgc_tracer_type, zbgc_frac_init, &
         tau_ret, tau_rel, R_Si2N, R_S2N, R_Fe2C, &
         R_Fe2N, R_Fe2DON, R_Fe2DOC, &
         chlabs, alpha2max_low, beta2max, &
         mu_max, grow_Tdep, fr_graze, &
         mort_pre, mort_Tdep, k_exude, &
         K_Nit, K_Am, K_Sil, K_Fe, &
         f_don, kn_bac, f_don_Am, &
         f_doc, f_exude, k_bac
         

      integer (kind=int_kind), intent(in) :: &
         nblyr     , & ! number of bio/brine layers per category 
         nilyr     , & ! number of ice layers per category
         nslyr     , & ! number of snow layers per category
         n_zaero   , & ! number of z aerosols in use 
         n_algae   , & ! number of algae in use 
         n_doc     , & ! number of DOC pools in use
         n_dic     , & ! number of DIC pools in use
         n_don     , & ! number of DON pools in use
         n_fed     , & ! number of Fe  pools in use dissolved Fe
         n_fep     , & ! number of Fe  pools in use particulate Fe
         max_algae , &
         max_doc   , &
         max_dic   , &
         max_don   , &
         max_fe

      integer (kind=int_kind), intent(inout) :: &
         ntrcr_o,     & ! number of non-bio tracers in use
         ntrcr,       & ! number of tracers
         nbtrcr,      & ! number of bgc tracers in use
         nbtrcr_sw      ! size of shorwave tracer vector

      integer (kind=int_kind), dimension (:), intent(inout) :: &
         trcr_depend   ! = 0 for ice area tracers
                       ! = 1 for ice volume tracers
                       ! = 2 for snow volume tracers

      integer (kind=int_kind), dimension (:), intent(inout) :: &
         n_trcr_strata ! number of underlying tracer layers

      integer (kind=int_kind), dimension (:,:), intent(inout) :: &
         nt_strata     ! indices of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcr_base     ! = 0 or 1 depending on tracer dependency
                       ! argument 2:  (1) aice, (2) vice, (3) vsno

      logical (kind=log_kind), intent(in) :: &
         tr_brine,       & ! if .true., brine height differs from ice thickness
         tr_bgc_S,       & ! if .true., use zsalinity
         tr_zaero,       & ! if .true., black carbon is tracers  (n_zaero)
         tr_bgc_Nit,     & ! if .true. Nitrate tracer in ice 
         tr_bgc_N,       & ! if .true., algal nitrogen tracers  (n_algae)
         tr_bgc_DON,     & ! if .true., DON pools are tracers  (n_don)
         tr_bgc_C,       & ! if .true., algal carbon tracers + DOC and DIC 
         tr_bgc_chl,     & ! if .true., algal chlorophyll tracers 
         tr_bgc_Am,      & ! if .true., ammonia/um as nutrient tracer 
         tr_bgc_Sil,     & ! if .true., silicon as nutrient tracer 
         tr_bgc_DMS,     & ! if .true., DMS as  tracer 
         tr_bgc_Fe,      & ! if .true., Fe as  tracer 
         tr_bgc_PON,     & ! if .true., PON as tracer 
         tr_bgc_hum,     & ! if .true., humic material as tracer
         solve_zsal,     & ! if true, update salinity profile from solve_S_dt
         z_tracers,      & ! if .true., bgc or aerosol tracers are vertically resolved
         solve_zbgc,     & ! if .true., solve vertical biochemistry portion of code
         dEdd_algae,     & ! if .true., algal absorption of Shortwave is computed in the 
         skl_bgc           ! if true, solve skeletal biochemistry

       integer (kind=int_kind), intent(out) :: &
         nt_fbri,      & ! volume fraction of ice with dynamic salt (hinS/vicen*aicen)
         nt_bgc_Nit,   & ! nutrients  
         nt_bgc_Am,    & ! 
         nt_bgc_Sil,   & !
         nt_bgc_DMSPp, & ! trace gases (skeletal layer)
         nt_bgc_DMSPd, & ! 
         nt_bgc_DMS,   & ! 
         nt_bgc_PON,   & ! zooplankton and detritus 
         nt_bgc_hum,   & ! humic material 
                         ! bio layer indicess
         nlt_bgc_Nit,  & ! nutrients  
         nlt_bgc_Am,   & ! 
         nlt_bgc_Sil,  & !
         nlt_bgc_DMSPp,& ! trace gases (skeletal layer)
         nlt_bgc_DMSPd,& ! 
         nlt_bgc_DMS,  & ! 
         nlt_bgc_PON,  & ! zooplankton and detritus 
         nlt_bgc_hum,  & ! humic material 
         nlt_chl_sw,   & ! points to total chla in trcrn_sw
         nt_zbgc_frac, & ! fraction of tracer in the mobile phase
         nt_bgc_S        ! Bulk salinity in fraction ice with dynamic salinity (Bio grid)

      integer (kind=int_kind), dimension(:), intent(out) :: &  
         nt_bgc_N , & ! diatoms, phaeocystis, pico/small   
         nt_bgc_C , & ! diatoms, phaeocystis, pico/small   
         nt_bgc_chl,& ! diatoms, phaeocystis, pico/small  
         nlt_bgc_N ,& ! diatoms, phaeocystis, pico/small   
         nlt_bgc_C ,& ! diatoms, phaeocystis, pico/small   
         nlt_bgc_chl   ! diatoms, phaeocystis, pico/small 

      integer (kind=int_kind), dimension(:), intent(out) :: &  
         nt_bgc_DOC,   & !  dissolved organic carbon  
         nlt_bgc_DOC     !  dissolved organic carbon

      integer (kind=int_kind), dimension(:), intent(out) :: & 
         nt_bgc_DON,   & !  dissolved organic nitrogen
         nlt_bgc_DON     !  dissolved organic nitrogen

      integer (kind=int_kind), dimension(:), intent(out) :: &  
         nt_bgc_DIC,    & !  dissolved inorganic carbon
         nlt_bgc_DIC      !  dissolved inorganic carbon

      integer (kind=int_kind), dimension(:), intent(out) :: & 
         nt_bgc_Fed,     & !  dissolved iron
         nt_bgc_Fep,     & !  particulate iron
         nlt_bgc_Fed,    & !  dissolved iron
         nlt_bgc_Fep       !  particulate iron

      integer (kind=int_kind), dimension(:), intent(out) :: &  
         nt_zaero,    & !  black carbon and other aerosols 
         nlt_zaero,   & !  black carbon and other aerosols
         nlt_zaero_sw   
    
      integer (kind=int_kind), dimension(:), intent(out) :: &   
         bio_index_o , & ! nlt  to appropriate value in ocean data array
         bio_index       ! nlt to nt

      real (kind=dbl_kind), intent(in) :: &
         initbio_frac, & ! fraction of ocean tracer concentration used to initialize tracer 
         frazil_scav     ! multiple of ocean tracer concentration due to frazil scavenging

      real (kind=dbl_kind), intent(in) :: &
        ratio_Si2N_diatoms, &   ! algal Si to N (mol/mol)
        ratio_Si2N_sp     , &
        ratio_Si2N_phaeo  , &
        ratio_S2N_diatoms , &   ! algal S  to N (mol/mol)
        ratio_S2N_sp      , &
        ratio_S2N_phaeo   , &
        ratio_Fe2C_diatoms, &   ! algal Fe to C  (umol/mol)
        ratio_Fe2C_sp     , &
        ratio_Fe2C_phaeo  , &
        ratio_Fe2N_diatoms, &   ! algal Fe to N  (umol/mol)
        ratio_Fe2N_sp     , &
        ratio_Fe2N_phaeo  , &
        ratio_Fe2DON      , &   ! Fe to N of DON (nmol/umol)
        ratio_Fe2DOC_s    , &   ! Fe to C of DOC (nmol/umol) saccharids
        ratio_Fe2DOC_l    , &   ! Fe to C of DOC (nmol/umol) lipids 
        tau_min           , &   ! rapid mobile to stationary exchanges (s) = 1.5 hours
        tau_max           , &   ! long time mobile to stationary exchanges (s) = 2 days
        chlabs_diatoms   , & ! chl absorption (1/m/(mg/m^3))
        chlabs_sp        , & !
        chlabs_phaeo     , & !
        alpha2max_low_diatoms , & ! light limitation (1/(W/m^2))  
        alpha2max_low_sp      , & 
        alpha2max_low_phaeo   , & 
        beta2max_diatoms , & ! light inhibition (1/(W/m^2))  
        beta2max_sp      , & 
        beta2max_phaeo   , & 
        mu_max_diatoms   , & ! maximum growth rate (1/day)       
        mu_max_sp        , & 
        mu_max_phaeo     , & 
        grow_Tdep_diatoms, & ! Temperature dependence of growth (1/C)
        grow_Tdep_sp     , & 
        grow_Tdep_phaeo  , & 
        fr_graze_diatoms , & ! Fraction grazed
        fr_graze_sp      , & 
        fr_graze_phaeo   , & 
        mort_pre_diatoms , & ! Mortality (1/day)
        mort_pre_sp      , & 
        mort_pre_phaeo   , & 
        mort_Tdep_diatoms, & ! T dependence of mortality (1/C)
        mort_Tdep_sp     , &  
        mort_Tdep_phaeo  , &  
        k_exude_diatoms  , & ! algal exudation (1/d)
        k_exude_sp       , &  
        k_exude_phaeo    , &  
        K_Nit_diatoms    , & ! nitrate half saturation (mmol/m^3)
        K_Nit_sp        , &  
        K_Nit_phaeo      , &  
        K_Am_diatoms     , & ! ammonium half saturation (mmol/m^3)
        K_Am_sp         , &   
        K_Am_phaeo       , &   
        K_Sil_diatoms    , & ! silicate half saturation (mmol/m^3)
        K_Sil_sp        , &   
        K_Sil_phaeo      , &   
        K_Fe_diatoms     , & ! iron half saturation (nM)
        K_Fe_sp         , &   
        K_Fe_phaeo       , &    
        f_don_protein    , & ! fraction of spilled grazing to proteins          
        kn_bac_protein   , & ! Bacterial degredation of DON (1/d)               
        f_don_Am_protein , & ! fraction of remineralized DON to ammonium        
        f_doc_s         , & ! fraction of mortality to DOC 
        f_doc_l         , &   
        f_exude_s        , & ! fraction of exudation to DOC
        f_exude_l        , & 
        k_bac_s         , & ! Bacterial degredation of DOC (1/d)
        k_bac_l         , & 
        algaltype_diatoms  , & ! mobility type
        algaltype_sp       , & !
        algaltype_phaeo    , & !
        nitratetype        , & !
        ammoniumtype       , & !
        silicatetype       , & !
        dmspptype         , & !
        dmspdtype         , & !
        humtype           , & !
        doctype_s         , & !
        doctype_l         , & !
        dontype_protein    , & !
        fedtype_1         , & !
        feptype_1         , & !
        zaerotype_bc1      , & !
        zaerotype_bc2      , & !
        zaerotype_dust1    , & !
        zaerotype_dust2    , & !
        zaerotype_dust3    , & !
        zaerotype_dust4    , & !
        ratio_C2N_diatoms  , & ! algal C to N ratio (mol/mol)
        ratio_C2N_sp       , & !
        ratio_C2N_phaeo    , & !
        ratio_chl2N_diatoms, & ! algal chlorophyll to N ratio (mg/mmol)
        ratio_chl2N_sp     , & !
        ratio_chl2N_phaeo  , & !
        F_abs_chl_diatoms  , & ! scales absorbed radiation for dEdd
        F_abs_chl_sp       , & !
        F_abs_chl_phaeo    , & !
        ratio_C2N_proteins     ! ratio of C to N in proteins (mol/mol)   

      ! local variables

      integer (kind=int_kind) :: &
        k, mm    , & ! loop index  
        ntd      , & ! for tracer dependency calculation
        nk       , & !
        nt_depend

      ntrcr_o = ntrcr
      nt_fbri = 0
      if (tr_brine) then
          nt_fbri = ntrcr + 1   ! ice volume fraction with salt
          ntrcr = ntrcr + 1
          trcr_depend(nt_fbri)   = 1   ! volume-weighted
          trcr_base  (nt_fbri,1) = c0  ! volume-weighted
          trcr_base  (nt_fbri,2) = c1  ! volume-weighted
          trcr_base  (nt_fbri,3) = c0  ! volume-weighted
          n_trcr_strata(nt_fbri) = 0
          nt_strata  (nt_fbri,1) = 0
          nt_strata  (nt_fbri,2) = 0
      endif

      ntd = 0                    ! if nt_fbri /= 0 then use fbri dependency
      if (nt_fbri == 0) ntd = -1 ! otherwise make tracers depend on ice volume

      if (solve_zsal) then       ! .true. only if tr_brine = .true.
          nt_bgc_S = ntrcr + 1
          ntrcr = ntrcr + nblyr
          do k = 1,nblyr
             trcr_depend(nt_bgc_S + k - 1) = 2 + nt_fbri + ntd
             trcr_base  (nt_bgc_S,1) = c0  ! default: ice area
             trcr_base  (nt_bgc_S,2) = c1 
             trcr_base  (nt_bgc_S,3) = c0  
             n_trcr_strata(nt_bgc_S) = 1
             nt_strata(nt_bgc_S,1) = nt_fbri
             nt_strata(nt_bgc_S,2) = 0
          enddo
      endif 

      !-----------------------------------------------------------------
      ! biogeochemistry
      !-----------------------------------------------------------------

      nbtrcr = 0
      nbtrcr_sw = 0

      ! vectors of size max_algae
      nlt_bgc_N(:) = 0
      nlt_bgc_C(:) = 0
      nlt_bgc_chl(:) = 0
      nt_bgc_N(:) = 0
      nt_bgc_C(:) = 0
      nt_bgc_chl(:) = 0

      ! vectors of size max_dic
      nlt_bgc_DIC(:) = 0
      nt_bgc_DIC(:) = 0

      ! vectors of size max_doc
      nlt_bgc_DOC(:) = 0
      nt_bgc_DOC(:) = 0

      ! vectors of size max_don
      nlt_bgc_DON(:) = 0
      nt_bgc_DON(:) = 0

      ! vectors of size max_fe 
      nlt_bgc_Fed(:) = 0
      nlt_bgc_Fep(:) = 0
      nt_bgc_Fed(:) = 0
      nt_bgc_Fep(:) = 0

      ! vectors of size max_aero
      nlt_zaero(:) = 0
      nlt_zaero_sw(:) = 0
      nt_zaero(:) = 0

      nlt_bgc_Nit    = 0
      nlt_bgc_Am     = 0
      nlt_bgc_Sil    = 0
      nlt_bgc_DMSPp  = 0
      nlt_bgc_DMSPd  = 0
      nlt_bgc_DMS    = 0
      nlt_bgc_PON    = 0
      nlt_bgc_hum    = 0
      nlt_chl_sw     = 0
      bio_index(:)   = 0
      bio_index_o(:) = 0

      nt_bgc_Nit    = 0
      nt_bgc_Am     = 0
      nt_bgc_Sil    = 0
      nt_bgc_DMSPp  = 0
      nt_bgc_DMSPd  = 0
      nt_bgc_DMS    = 0
      nt_bgc_PON    = 0
      nt_bgc_hum    = 0

      !-----------------------------------------------------------------
      ! Define array parameters
      !-----------------------------------------------------------------
      R_Si2N(1) = ratio_Si2N_diatoms
      R_Si2N(2) = ratio_Si2N_sp
      R_Si2N(3) = ratio_Si2N_phaeo

      R_S2N(1) = ratio_S2N_diatoms
      R_S2N(2) = ratio_S2N_sp
      R_S2N(3) = ratio_S2N_phaeo

      R_Fe2C(1) = ratio_Fe2C_diatoms
      R_Fe2C(2) = ratio_Fe2C_sp
      R_Fe2C(3) = ratio_Fe2C_phaeo

      R_Fe2N(1) = ratio_Fe2N_diatoms
      R_Fe2N(2) = ratio_Fe2N_sp
      R_Fe2N(3) = ratio_Fe2N_phaeo

      R_C2N(1) = ratio_C2N_diatoms
      R_C2N(2) = ratio_C2N_sp
      R_C2N(3) = ratio_C2N_phaeo

      R_chl2N(1) = ratio_chl2N_diatoms
      R_chl2N(2) = ratio_chl2N_sp
      R_chl2N(3) = ratio_chl2N_phaeo

      F_abs_chl(1) = F_abs_chl_diatoms
      F_abs_chl(2) = F_abs_chl_sp
      F_abs_chl(3) = F_abs_chl_phaeo

      R_Fe2DON(1) = ratio_Fe2DON
      R_C2N(1) = ratio_C2N_proteins
     
      R_Fe2DOC(1) = ratio_Fe2DOC_s
      R_Fe2DOC(2) = ratio_Fe2DOC_l
      R_Fe2DOC(3) = c0

      chlabs(1) = chlabs_diatoms
      chlabs(2) = chlabs_sp
      chlabs(3) = chlabs_phaeo

      alpha2max_low(1) = alpha2max_low_diatoms
      alpha2max_low(2) = alpha2max_low_sp
      alpha2max_low(3) = alpha2max_low_phaeo

      beta2max(1) = beta2max_diatoms
      beta2max(2) = beta2max_sp
      beta2max(3) = beta2max_phaeo

      mu_max(1) = mu_max_diatoms
      mu_max(2) = mu_max_sp
      mu_max(3) = mu_max_phaeo

      grow_Tdep(1) = grow_Tdep_diatoms
      grow_Tdep(2) = grow_Tdep_sp
      grow_Tdep(3) = grow_Tdep_phaeo

      fr_graze(1) = fr_graze_diatoms
      fr_graze(2) = fr_graze_sp
      fr_graze(3) = fr_graze_phaeo

      mort_pre(1) = mort_pre_diatoms
      mort_pre(2) = mort_pre_sp
      mort_pre(3) = mort_pre_phaeo

      mort_Tdep(1) = mort_Tdep_diatoms
      mort_Tdep(2) = mort_Tdep_sp
      mort_Tdep(3) = mort_Tdep_phaeo

      k_exude(1) = k_exude_diatoms
      k_exude(2) = k_exude_sp
      k_exude(3) = k_exude_phaeo

      K_Nit(1) = K_Nit_diatoms
      K_Nit(2) = K_Nit_sp
      K_Nit(3) = K_Nit_phaeo

      K_Am(1) = K_Am_diatoms
      K_Am(2) = K_Am_sp
      K_Am(3) = K_Am_phaeo

      K_Sil(1) = K_Sil_diatoms
      K_Sil(2) = K_Sil_sp
      K_Sil(3) = K_Sil_phaeo

      K_Fe(1) = K_Fe_diatoms
      K_Fe(2) = K_Fe_sp
      K_Fe(3) = K_Fe_phaeo

      f_don(1) = f_don_protein
      kn_bac(1) = kn_bac_protein
      f_don_Am(1) = f_don_Am_protein

      f_exude(1) = f_exude_s
      f_exude(2) = f_exude_l
      k_bac(1) = k_bac_s
      k_bac(2) = k_bac_l
      
      algaltype(1) = algaltype_diatoms
      algaltype(2) = algaltype_sp
      algaltype(3) = algaltype_phaeo

      doctype(1) = doctype_s
      doctype(2) = doctype_l
 
      dontype(1) = dontype_protein

      fedtype(1) = fedtype_1
      feptype(1) = feptype_1

      zaerotype(1) = zaerotype_bc1
      zaerotype(2) = zaerotype_bc2
      zaerotype(3) = zaerotype_dust1
      zaerotype(4) = zaerotype_dust2
      zaerotype(5) = zaerotype_dust3
      zaerotype(6) = zaerotype_dust4     

      if (skl_bgc) then

         nk = 1
         nt_depend = 0

         if (dEdd_algae) then 
           nlt_chl_sw = 1
           nbtrcr_sw = nilyr+nslyr+2  ! only the bottom layer 
                                                 ! will be nonzero    
         endif  
         
      elseif (z_tracers) then ! defined on nblyr+1 in ice 
                              ! and 2 snow layers (snow surface + interior)

         nk = nblyr + 1
         nt_depend = 2 + nt_fbri + ntd 

         if (tr_bgc_N) then
            if (dEdd_algae) then
               nlt_chl_sw = 1
               nbtrcr_sw =  nilyr+nslyr+2
            endif
         endif ! tr_bgc_N

      endif ! skl_bgc or z_tracers

      if (skl_bgc .or. z_tracers) then

      !-----------------------------------------------------------------
      ! assign tracer indices and dependencies
      ! bgc_tracer_type: < 0  purely mobile , >= 0 stationary 
      !------------------------------------------------------------------

      if (tr_bgc_N) then
         do mm = 1, n_algae      
            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_N(mm),    nlt_bgc_N(mm), &
                                      algaltype(mm),   nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_N(mm)) = mm
         enddo   ! mm
      endif ! tr_bgc_N

      if (tr_bgc_Nit) then
            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_Nit,      nlt_bgc_Nit,   &
                                      nitratetype,     nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_Nit) = max_algae + 1
      endif ! tr_bgc_Nit
         
      if (tr_bgc_C) then
       !
       ! Algal C is not yet distinct from algal N
       ! * Reqires exudation and/or changing C:N ratios
       ! for implementation
       !
       !  do mm = 1,n_algae      
       !     call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
       !                               nt_bgc_C(mm),    nlt_bgc_C(mm), &
       !                               algaltype(mm),   nt_depend,     &
       !                               ntrcr,           nbtrcr,        &
       !                               bgc_tracer_type, trcr_depend,   &
       !                               trcr_base,       n_trcr_strata, &
       !                               nt_strata,       bio_index)
       !     bio_index_o(nlt_bgc_C(mm)) = max_algae + 1 + mm
       !  enddo   ! mm

         do mm = 1, n_doc
            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_DOC(mm),  nlt_bgc_DOC(mm), &
                                      doctype(mm),     nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DOC(mm)) = max_algae + 1 + mm
         enddo   ! mm
         do mm = 1, n_dic
            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_DIC(mm),  nlt_bgc_DIC(mm), &
                                      dictype(mm),     nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DIC(mm)) = max_algae + max_doc + 1 + mm
         enddo   ! mm
      endif      ! tr_bgc_C

      if (tr_bgc_chl) then
         do mm = 1, n_algae
            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_chl(mm),  nlt_bgc_chl(mm), &
                                      algaltype(mm),   nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_chl(mm)) = max_algae + 1 + max_doc + max_dic + mm
         enddo   ! mm
      endif      ! tr_bgc_chl

      if (tr_bgc_Am) then
            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_Am,       nlt_bgc_Am,    &
                                      ammoniumtype,    nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_Am) = 2*max_algae + max_doc + max_dic + 2
      endif    
      if (tr_bgc_Sil) then
            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_Sil,      nlt_bgc_Sil,   &
                                      silicatetype,    nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_Sil) = 2*max_algae + max_doc + max_dic + 3
      endif    
      if (tr_bgc_DMS) then   ! all together
            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_DMSPp,    nlt_bgc_DMSPp, &
                                      dmspptype,       nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DMSPp) = 2*max_algae + max_doc + max_dic + 4

            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_DMSPd,    nlt_bgc_DMSPd, &
                                      dmspdtype,       nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DMSPd) = 2*max_algae + max_doc + max_dic + 5

            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_DMS,      nlt_bgc_DMS,   &
                                      dmspdtype,       nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DMS) = 2*max_algae + max_doc + max_dic + 6
      endif    
      if (tr_bgc_PON) then
            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_PON,      nlt_bgc_PON, &
                                      nitratetype,     nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_PON) =  2*max_algae + max_doc + max_dic + 7
      endif
      if (tr_bgc_DON) then
         do mm = 1, n_don
            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_DON(mm),  nlt_bgc_DON(mm), &
                                      dontype(mm),     nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_DON(mm)) = 2*max_algae + max_doc + max_dic + 7 + mm
         enddo   ! mm
      endif      ! tr_bgc_DON
      if (tr_bgc_Fe) then
         do mm = 1, n_fed
            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_Fed(mm),  nlt_bgc_Fed(mm), &
                                      fedtype(mm),     nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_Fed(mm)) = 2*max_algae + max_doc + max_dic &
                                         + max_don + 7 + mm
         enddo   ! mm
         do mm = 1, n_fep
            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_Fep(mm),  nlt_bgc_Fep(mm), &
                                      feptype(mm),     nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_Fep(mm)) = 2*max_algae + max_doc + max_dic &
                                         + max_don + max_fe + 7 + mm
         enddo   ! mm
      endif      ! tr_bgc_Fe 
  
      if (tr_bgc_hum) then
            call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc_hum,      nlt_bgc_hum,   &
                                      humtype,         nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)
            bio_index_o(nlt_bgc_hum) =   2*max_algae + max_doc + 8 + max_dic &
                                         + max_don + 2*max_fe + max_aero 
      endif
      endif  ! skl_bgc or z_tracers

      if (z_tracers) then ! defined on nblyr+1 in ice 
                              ! and 2 snow layers (snow surface + interior)

         nk = nblyr + 1
         nt_depend = 2 + nt_fbri + ntd 

         ! z layer aerosols
         if (tr_zaero) then
            do mm = 1, n_zaero
               if (dEdd_algae) then
                  nlt_zaero_sw(mm) = nbtrcr_sw + 1
                  nbtrcr_sw = nbtrcr_sw + nilyr + nslyr+2
               endif
               call colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                         nt_zaero(mm),    nlt_zaero(mm), &
                                         zaerotype(mm),   nt_depend,     &
                                         ntrcr,           nbtrcr,        &
                                         bgc_tracer_type, trcr_depend,   &
                                         trcr_base,       n_trcr_strata, &
                                         nt_strata,       bio_index)
               bio_index_o(nlt_zaero(mm)) = 2*max_algae + max_doc + max_dic &
                                          + max_don + 2*max_fe + 7 + mm
            enddo   ! mm
         endif      ! tr_zaero

         nt_zbgc_frac = 0
         if (nbtrcr > 0) then
            nt_zbgc_frac = ntrcr + 1
            ntrcr = ntrcr + nbtrcr
            do k = 1,nbtrcr 
               zbgc_frac_init(k) = c1   
               trcr_depend(nt_zbgc_frac+k-1) =  2+nt_fbri 
               trcr_base(nt_zbgc_frac+ k - 1,1)  = c0
               trcr_base(nt_zbgc_frac+ k - 1,2)  = c1
               trcr_base(nt_zbgc_frac+ k - 1,3)  = c0
               n_trcr_strata(nt_zbgc_frac+ k - 1)= 1  
               nt_strata(nt_zbgc_frac+ k - 1,1)  = nt_fbri
               nt_strata(nt_zbgc_frac+ k - 1,2)  = 0 
               tau_ret(k) = c1
               tau_rel(k) = c1
               if (bgc_tracer_type(k) >=  c0 .and. bgc_tracer_type(k) < p5) then
                  tau_ret(k) = tau_min
                  tau_rel(k) = tau_max
                  zbgc_frac_init(k) = c1
               elseif (bgc_tracer_type(k) >= p5 .and. bgc_tracer_type(k) < c1) then
                  tau_ret(k) = tau_min
                  tau_rel(k) = tau_min
                  zbgc_frac_init(k) = c1
               elseif (bgc_tracer_type(k) >= c1 .and. bgc_tracer_type(k) < c2) then
                  tau_ret(k) = tau_max
                  tau_rel(k) = tau_min
                  zbgc_frac_init(k) = c1
               elseif (bgc_tracer_type(k) >= c2 ) then
                  tau_ret(k) = tau_max
                  tau_rel(k) = tau_max
                  zbgc_frac_init(k) = c1
               endif
            enddo
         endif

      endif ! z_tracers

      do k = 1, nbtrcr
         zbgc_init_frac(k) = frazil_scav
         if (bgc_tracer_type(k) < c0)  zbgc_init_frac(k) = initbio_frac
      enddo  

      end subroutine colpkg_init_zbgc

!=======================================================================

      subroutine colpkg_init_bgc_trcr(nk,              nt_fbri,       &
                                      nt_bgc,          nlt_bgc,       &
                                      bgctype,         nt_depend,     &
                                      ntrcr,           nbtrcr,        &
                                      bgc_tracer_type, trcr_depend,   &
                                      trcr_base,       n_trcr_strata, &
                                      nt_strata,       bio_index)

      use ice_constants_colpkg, only: c0, c1

      integer (kind=int_kind), intent(in) :: &
         nk           , & ! counter
         nt_depend    , & ! tracer dependency index
         nt_fbri

      integer (kind=int_kind), intent(inout) :: &
         ntrcr        , & ! number of tracers
         nbtrcr       , & ! number of bio tracers
         nt_bgc       , & ! tracer index
         nlt_bgc          ! bio tracer index

      integer (kind=int_kind), dimension(:), intent(inout) :: &
         trcr_depend  , & ! tracer dependencies
         n_trcr_strata, & ! number of underlying tracer layers
         bio_index        !

      integer (kind=int_kind), dimension(:,:), intent(inout) :: &
         nt_strata        ! indices of underlying tracer layers

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         trcr_base        ! = 0 or 1 depending on tracer dependency
                          ! argument 2:  (1) aice, (2) vice, (3) vsno

      real (kind=dbl_kind), intent(in) :: &
         bgctype          ! bio tracer transport type (mobile vs stationary)

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         bgc_tracer_type  ! bio tracer transport type array

      ! local variables

      integer (kind=int_kind) :: &
         k         , & ! loop index
         n_strata  , & ! temporary values
         nt_strata1, & ! 
         nt_strata2

      real (kind=dbl_kind) :: &
         trcr_base1, & ! temporary values
         trcr_base2, &
         trcr_base3

         nt_bgc = ntrcr + 1 
         nbtrcr = nbtrcr + 1
         nlt_bgc = nbtrcr
         bgc_tracer_type(nbtrcr) = bgctype
         
         if (nk > 1) then 
            ! include vertical bgc in snow
            do k = nk, nk+1
               ntrcr = ntrcr + 1
               trcr_depend  (nt_bgc + k  ) = 2 ! snow volume
               trcr_base    (nt_bgc + k,1) = c0
               trcr_base    (nt_bgc + k,2) = c0
               trcr_base    (nt_bgc + k,3) = c1
               n_trcr_strata(nt_bgc + k  ) = 0
               nt_strata    (nt_bgc + k,1) = 0
               nt_strata    (nt_bgc + k,2) = 0
            enddo

            trcr_base1 = c0      
            trcr_base2 = c1     
            trcr_base3 = c0
            n_strata = 1    
            nt_strata1 = nt_fbri
            nt_strata2 = 0
         else  ! nk = 1
            trcr_base1 = c1
            trcr_base2 = c0
            trcr_base3 = c0
            n_strata = 0
            nt_strata1 = 0
            nt_strata2 = 0
         endif ! nk

         do k = 1, nk     !in ice
            ntrcr = ntrcr + 1
            trcr_depend  (nt_bgc + k - 1  ) = nt_depend
            trcr_base    (nt_bgc + k - 1,1) = trcr_base1
            trcr_base    (nt_bgc + k - 1,2) = trcr_base2
            trcr_base    (nt_bgc + k - 1,3) = trcr_base3
            n_trcr_strata(nt_bgc + k - 1  ) = n_strata
            nt_strata    (nt_bgc + k - 1,1) = nt_strata1
            nt_strata    (nt_bgc + k - 1,2) = nt_strata2
         enddo

         bio_index (nlt_bgc) = nt_bgc

      end subroutine colpkg_init_bgc_trcr

!=======================================================================
!     Temperature functions
!=======================================================================

      function colpkg_liquidus_temperature(Sin) result(Tmlt)

        use ice_colpkg_shared, only: ktherm
        use ice_constants_colpkg, only: depressT
        use ice_mushy_physics, only: liquidus_temperature_mush

        real(dbl_kind), intent(in) :: Sin
        real(dbl_kind) :: Tmlt

        if (ktherm == 2) then

           Tmlt = liquidus_temperature_mush(Sin)

        else

           Tmlt = -depressT * Sin

        endif

      end function colpkg_liquidus_temperature

!=======================================================================

      function colpkg_sea_freezing_temperature(sss) result(Tf)

        use ice_colpkg_shared, only: tfrz_option
        use ice_constants_colpkg, only: depressT

        real(dbl_kind), intent(in) :: sss
        real(dbl_kind) :: Tf

        if (trim(tfrz_option) == 'mushy') then

           Tf = colpkg_liquidus_temperature(sss) ! deg C
           
        elseif (trim(tfrz_option) == 'linear_salt') then

           Tf = -depressT * sss ! deg C

        else

           Tf = -1.8_dbl_kind

        endif

      end function colpkg_sea_freezing_temperature

!=======================================================================

      function colpkg_ice_temperature(qin, Sin) result(Tin)

        use ice_colpkg_shared, only: ktherm
        use ice_constants_colpkg, only: depressT
        use ice_mushy_physics, only: temperature_mush
        use ice_therm_shared, only: calculate_Tin_from_qin

        real(kind=dbl_kind), intent(in) :: qin, Sin
        real(kind=dbl_kind) :: Tin

        real(kind=dbl_kind) :: Tmlts

        if (ktherm == 2) then

           Tin = temperature_mush(qin, Sin)

        else

           Tmlts = -depressT * Sin
           Tin = calculate_Tin_from_qin(qin,Tmlts)

        endif

      end function colpkg_ice_temperature

!=======================================================================

      function colpkg_snow_temperature(qin) result(Tsn)

        use ice_colpkg_shared, only: ktherm
        use ice_mushy_physics, only: temperature_snow
        use ice_constants_colpkg, only: Lfresh, rhos, cp_ice

        real(kind=dbl_kind), intent(in) :: qin
        real(kind=dbl_kind) :: Tsn

        if (ktherm == 2) then

           Tsn = temperature_snow(qin)

        else

           Tsn = (Lfresh + qin/rhos)/cp_ice

        endif

      end function colpkg_snow_temperature

!=======================================================================

      function colpkg_enthalpy_snow(zTsn) result(qsn)

        use ice_mushy_physics, only: enthalpy_snow

        real(kind=dbl_kind), intent(in) :: zTsn
        real(kind=dbl_kind) :: qsn

        qsn = enthalpy_snow(zTsn)

      end function colpkg_enthalpy_snow

!=======================================================================
!     Time-stepping routines
!=======================================================================

! Driver for thermodynamic changes not needed for coupling:
! transport in thickness space, lateral growth and melting.
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL

      subroutine colpkg_step_therm1(dt, ncat, nilyr, nslyr, n_aero, &
                                    aicen_init  ,               &
                                    vicen_init  , vsnon_init  , &
                                    aice        , aicen       , &
                                    vice        , vicen       , &
                                    vsno        , vsnon       , &
                                    uvel        , vvel        , &
                                    Tsfc        , zqsn        , &
                                    zqin        , zSin        , &
                                    alvl        , vlvl        , &
                                    apnd        , hpnd        , &
                                    ipnd        ,               &
                                    iage        , FY          , &
                                    aerosno     , aeroice     , &
                                    uatm        , vatm        , &
                                    wind        , zlvl        , &
                                    Qa          , rhoa        , &
                                    Tair        , Tref        , &
                                    Qref        , Uref        , &
                                    Cdn_atm_ratio,              &
                                    Cdn_ocn     , Cdn_ocn_skin, &
                                    Cdn_ocn_floe, Cdn_ocn_keel, &
                                    Cdn_atm     , Cdn_atm_skin, &
                                    Cdn_atm_floe, Cdn_atm_pond, &
                                    Cdn_atm_rdg , hfreebd     , &
                                    hdraft      , hridge      , &
                                    distrdg     , hkeel       , &
                                    dkeel       , lfloe       , &
                                    dfloe       ,               &
                                    strax       , stray       , &
                                    strairxT    , strairyT    , &
                                    potT        , sst         , &
                                    sss         , Tf          , &
                                    strocnxT    , strocnyT    , &
                                    fbot        ,               &
                                    frzmlt      , rside       , &
                                    fsnow       , frain       , &
                                    fpond       ,               &
                                    fsurf       , fsurfn      , &
                                    fcondtop    , fcondtopn   , &
                                    fswsfcn     , fswintn     , &
                                    fswthrun    , fswabs      , &
                                    flwout      ,               &
                                    Sswabsn     , Iswabsn     , &
                                    flw         , coszen      , & 
                                    fsens       , fsensn      , &
                                    flat        , flatn       , &
                                    evap        ,               &
                                    fresh       , fsalt       , &
                                    fhocn       , fswthru     , &
                                    flatn_f     , fsensn_f    , &
                                    fsurfn_f    , fcondtopn_f , &
                                    faero_atm   , faero_ocn   , &
                                    dhsn        , ffracn      , &
                                    meltt       , melttn      , &
                                    meltb       , meltbn      , &
                                    meltl       ,               &
                                    melts       , meltsn      , &
                                    congel      , congeln     , &
                                    snoice      , snoicen     , &
                                    dsnown      , frazil      , &
                                    lmask_n     , lmask_s     , &
                                    mlt_onset   , frz_onset   , &
                                    yday        , l_stop      , &
                                    stop_label  , prescribed_ice)

      use ice_aerosol, only: update_aerosol
      use ice_atmo, only: neutral_drag_coeffs
      use ice_age, only: increment_age
      use ice_constants_colpkg, only: rhofresh, rhoi, rhos, c0, c1, puny
      use ice_firstyear, only: update_FYarea
      use ice_flux_colpkg, only: set_sfcflux, merge_fluxes
      use ice_meltpond_cesm, only: compute_ponds_cesm
      use ice_meltpond_lvl, only: compute_ponds_lvl
      use ice_meltpond_topo, only: compute_ponds_topo
      use ice_therm_shared, only: hi_min
      use ice_therm_vertical, only: frzmlt_bottom_lateral, thermo_vertical
      use ice_colpkg_tracers, only: tr_iage, tr_FY, tr_aero, tr_pond, &
          tr_pond_cesm, tr_pond_lvl, tr_pond_topo

      integer (kind=int_kind), intent(in) :: &
         ncat    , & ! number of thickness categories
         nilyr   , & ! number of ice layers
         nslyr   , & ! number of snow layers
         n_aero      ! number of aerosol tracers in use

      real (kind=dbl_kind), intent(in) :: &
         dt          , & ! time step
         uvel        , & ! x-component of velocity (m/s)
         vvel        , & ! y-component of velocity (m/s)
         strax       , & ! wind stress components (N/m^2)
         stray       , & ! 
         yday            ! day of year

      logical (kind=log_kind), intent(in) :: &
         lmask_n     , & ! northern hemisphere mask
         lmask_s         ! southern hemisphere mask

      logical (kind=log_kind), intent(in), optional :: &
         prescribed_ice  ! if .true., use prescribed ice instead of computed

      real (kind=dbl_kind), intent(inout) :: &
         aice        , & ! sea ice concentration
         vice        , & ! volume per unit area of ice          (m)
         vsno        , & ! volume per unit area of snow         (m)
         zlvl        , & ! atm level height (m)
         uatm        , & ! wind velocity components (m/s)
         vatm        , &
         wind        , & ! wind speed (m/s)
         potT        , & ! air potential temperature  (K)
         Tair        , & ! air temperature  (K)
         Qa          , & ! specific humidity (kg/kg)
         rhoa        , & ! air density (kg/m^3)
         frain       , & ! rainfall rate (kg/m^2 s)
         fsnow       , & ! snowfall rate (kg/m^2 s)
         fpond       , & ! fresh water flux to ponds (kg/m^2/s)
         fresh       , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt       , & ! salt flux to ocean (kg/m^2/s)
         fhocn       , & ! net heat flux to ocean (W/m^2)
         fswthru     , & ! shortwave penetrating to ocean (W/m^2)
         fsurf       , & ! net surface heat flux (excluding fcondtop)(W/m^2)
         fcondtop    , & ! top surface conductive flux        (W/m^2)
         fsens       , & ! sensible heat flux (W/m^2)
         flat        , & ! latent heat flux   (W/m^2)
         fswabs      , & ! shortwave flux absorbed in ice and ocean (W/m^2)
         coszen      , & ! cosine solar zenith angle, < 0 for sun below horizon 
         flw         , & ! incoming longwave radiation (W/m^2)
         flwout      , & ! outgoing longwave radiation (W/m^2)
         evap        , & ! evaporative water flux (kg/m^2/s)
         congel      , & ! basal ice growth         (m/step-->cm/day)
         frazil      , & ! frazil ice growth        (m/step-->cm/day)
         snoice      , & ! snow-ice formation       (m/step-->cm/day)
         Tref        , & ! 2m atm reference temperature (K)
         Qref        , & ! 2m atm reference spec humidity (kg/kg)
         Uref        , & ! 10m atm reference wind speed (m/s)
         Cdn_atm     , & ! atm drag coefficient
         Cdn_ocn     , & ! ocn drag coefficient
         hfreebd     , & ! freeboard (m)
         hdraft      , & ! draft of ice + snow column (Stoessel1993)
         hridge      , & ! ridge height
         distrdg     , & ! distance between ridges
         hkeel       , & ! keel depth
         dkeel       , & ! distance between keels
         lfloe       , & ! floe length
         dfloe       , & ! distance between floes
         Cdn_atm_skin, & ! neutral skin drag coefficient
         Cdn_atm_floe, & ! neutral floe edge drag coefficient
         Cdn_atm_pond, & ! neutral pond edge drag coefficient
         Cdn_atm_rdg , & ! neutral ridge drag coefficient
         Cdn_ocn_skin, & ! skin drag coefficient
         Cdn_ocn_floe, & ! floe edge drag coefficient
         Cdn_ocn_keel, & ! keel drag coefficient
         Cdn_atm_ratio,& ! ratio drag atm / neutral drag atm
         strairxT    , & ! stress on ice by air, x-direction
         strairyT    , & ! stress on ice by air, y-direction
         strocnxT    , & ! ice-ocean stress, x-direction
         strocnyT    , & ! ice-ocean stress, y-direction
         fbot        , & ! ice-ocean heat flux at bottom surface (W/m^2)
         frzmlt      , & ! freezing/melting potential (W/m^2)
         rside       , & ! fraction of ice that melts laterally
         sst         , & ! sea surface temperature (C)
         Tf          , & ! freezing temperature (C)
         sss         , & ! sea surface salinity (ppt)
         meltt       , & ! top ice melt             (m/step-->cm/day)
         melts       , & ! snow melt                (m/step-->cm/day)
         meltb       , & ! basal ice melt           (m/step-->cm/day)
         meltl       , & ! lateral ice melt         (m/step-->cm/day)
         mlt_onset   , & ! day of year that sfc melting begins
         frz_onset       ! day of year that freezing begins (congel or frazil)

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         aicen_init  , & ! fractional area of ice
         vicen_init  , & ! volume per unit area of ice (m)
         vsnon_init  , & ! volume per unit area of snow (m)
         aicen       , & ! concentration of ice
         vicen       , & ! volume per unit area of ice          (m)
         vsnon       , & ! volume per unit area of snow         (m)
         Tsfc        , & ! ice/snow surface temperature, Tsfcn
         alvl        , & ! level ice area fraction
         vlvl        , & ! level ice volume fraction
         apnd        , & ! melt pond area fraction
         hpnd        , & ! melt pond depth (m)
         ipnd        , & ! melt pond refrozen lid thickness (m)
         iage        , & ! volume-weighted ice age
         FY          , & ! area-weighted first-year ice area
         fsurfn      , & ! net flux to top surface, excluding fcondtop
         fcondtopn   , & ! downward cond flux at top surface (W m-2)
         flatn       , & ! latent heat flux (W m-2)
         fsensn      , & ! sensible heat flux (W m-2)
         fsurfn_f    , & ! net flux to top surface, excluding fcondtop
         fcondtopn_f , & ! downward cond flux at top surface (W m-2)
         flatn_f     , & ! latent heat flux (W m-2)
         fsensn_f    , & ! sensible heat flux (W m-2)
         fswsfcn     , & ! SW absorbed at ice/snow surface (W m-2)
         fswthrun    , & ! SW through ice to ocean            (W/m^2)
         fswintn     , & ! SW absorbed in ice interior, below surface (W m-2)
         faero_atm   , & ! aerosol deposition rate (kg/m^2 s)
         faero_ocn   , & ! aerosol flux to ocean  (kg/m^2/s)
         dhsn        , & ! depth difference for snow on sea ice and pond ice
         ffracn      , & ! fraction of fsurfn used to melt ipond
         meltsn      , & ! snow melt                       (m)
         melttn      , & ! top ice melt                    (m)
         meltbn      , & ! bottom ice melt                 (m)
         congeln     , & ! congelation ice growth          (m)
         snoicen     , & ! snow-ice growth                 (m)
         dsnown          ! change in snow thickness (m/step-->cm/day)

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         zqsn        , & ! snow layer enthalpy (J m-3)
         zqin        , & ! ice layer enthalpy (J m-3)
         zSin        , & ! internal ice layer salinities
         Sswabsn     , & ! SW radiation absorbed in snow layers (W m-2)
         Iswabsn         ! SW radiation absorbed in ice layers (W m-2)

      real (kind=dbl_kind), dimension(:,:,:), intent(inout) :: &
         aerosno    , &  ! snow aerosol tracer (kg/m^2)
         aeroice         ! ice aerosol tracer (kg/m^2)

      logical (kind=log_kind), intent(out) :: &
         l_stop          ! if true, abort model

      character (len=*), intent(out) :: &
         stop_label      ! abort error message

      ! local variables

      integer (kind=int_kind) :: &
         n               ! category index

      real (kind=dbl_kind) :: &
         worka, workb    ! temporary variables

      ! 2D coupler variables (computed for each category, then aggregated)
      real (kind=dbl_kind) :: &
         fswabsn     , & ! shortwave absorbed by ice          (W/m^2)
         flwoutn     , & ! upward LW at surface               (W/m^2)
         evapn       , & ! flux of vapor, atmos to ice   (kg m-2 s-1)
         freshn      , & ! flux of water, ice to ocean     (kg/m^2/s)
         fsaltn      , & ! flux of salt, ice to ocean      (kg/m^2/s)
         fhocnn      , & ! fbot corrected for leftover energy (W/m^2)
         strairxn    , & ! air/ice zonal  stress,             (N/m^2)
         strairyn    , & ! air/ice meridional stress,         (N/m^2)
         Cdn_atm_ratio_n, & ! drag coefficient ratio
         Trefn       , & ! air tmp reference level                (K)
         Urefn       , & ! air speed reference level            (m/s)
         Qrefn       , & ! air sp hum reference level         (kg/kg)
         Tbot        , & ! ice bottom surface temperature (deg C)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef      , & ! transfer coefficient for latent heat
         rfrac           ! water fraction retained for melt ponds

      real (kind=dbl_kind) :: &
         raice       , & ! 1/aice
         pond            ! water retained in ponds (m)

      !-----------------------------------------------------------------
      ! Adjust frzmlt to account for ice-ocean heat fluxes since last
      !  call to coupler.
      ! Compute lateral and bottom heat fluxes.
      !-----------------------------------------------------------------

      call frzmlt_bottom_lateral (dt,        ncat,      &
                                  nilyr,     nslyr,     &
                                  aice,      frzmlt,    &
                                  vicen,     vsnon,     &
                                  zqin,      zqsn,      &
                                  sst,       Tf,        &
                                  ustar_min,            &
                                  fbot_xfer_type,       &
                                  strocnxT,  strocnyT,  &
                                  Tbot,      fbot,      &
                                  rside,     Cdn_ocn)
      
      !-----------------------------------------------------------------
      ! Update the neutral drag coefficients to account for form drag
      ! Oceanic and atmospheric drag coefficients
      !-----------------------------------------------------------------

      if (formdrag) then
         call neutral_drag_coeffs (apnd         , &
                                   hpnd        , ipnd         , &
                                   alvl        , vlvl         , &
                                   aice        , vice,          &
                                   vsno        , aicen        , &
                                   vicen       , vsnon        , &
                                   Cdn_ocn     , Cdn_ocn_skin, &
                                   Cdn_ocn_floe, Cdn_ocn_keel, &
                                   Cdn_atm     , Cdn_atm_skin, &
                                   Cdn_atm_floe, Cdn_atm_pond, &
                                   Cdn_atm_rdg , hfreebd     , &
                                   hdraft      , hridge      , &
                                   distrdg     , hkeel       , &
                                   dkeel       , lfloe       , &
                                   dfloe       , ncat)
      endif

      do n = 1, ncat

         meltsn (n) = c0
         melttn (n) = c0
         meltbn (n) = c0
         congeln(n) = c0
         snoicen(n) = c0
         dsnown (n) = c0

         Trefn  = c0
         Qrefn  = c0
         Urefn  = c0
         lhcoef = c0
         shcoef = c0
         worka  = c0
         workb  = c0

         if (aicen_init(n) > puny) then

            if (calc_Tsfc .or. calc_strair) then 

      !-----------------------------------------------------------------
      ! Atmosphere boundary layer calculation; compute coefficients
      ! for sensible and latent heat fluxes.
      !
      ! NOTE: The wind stress is computed here for later use if 
      !       calc_strair = .true.   Otherwise, the wind stress
      !       components are set to the data values.
      !-----------------------------------------------------------------

               call colpkg_atm_boundary( 'ice',                  &
                                        Tsfc(n),  potT,          &
                                        uatm,     vatm,          &
                                        wind,     zlvl,          &
                                        Qa,       rhoa,          &
                                        strairxn, strairyn,      &
                                        Trefn,    Qrefn,         &
                                        worka,    workb,         &
                                        lhcoef,   shcoef,        &
                                        Cdn_atm,                 &
                                        Cdn_atm_ratio_n,         &
                                        uvel,     vvel,          &
                                        Uref=Urefn)

            endif   ! calc_Tsfc or calc_strair

            if (.not.(calc_strair)) then
#ifndef CICE_IN_NEMO
               ! Set to data values (on T points)
               strairxn = strax
               strairyn = stray
#else
               ! NEMO wind stress is supplied on u grid, multipied 
               ! by ice concentration and set directly in evp, so
               ! strairxT/yT = 0. Zero u-components here for safety.
               strairxn = c0
               strairyn = c0
#endif
            endif

      !-----------------------------------------------------------------
      ! Update ice age
      ! This is further adjusted for freezing in the thermodynamics.
      ! Melting does not alter the ice age.
      !-----------------------------------------------------------------

            if (tr_iage) call increment_age (dt, iage(n))
            if (tr_FY)   call update_FYarea (dt,               &
                                             lmask_n, lmask_s, &
                                             yday,    FY(n))

      !-----------------------------------------------------------------
      ! Vertical thermodynamics: Heat conduction, growth and melting.
      !----------------------------------------------------------------- 

            if (.not.(calc_Tsfc)) then

               ! If not calculating surface temperature and fluxes, set 
               ! surface fluxes (flatn, fsurfn, and fcondtopn) to be used 
               ! in thickness_changes
 
               ! hadgem routine sets fluxes to default values in ice-only mode
               call set_sfcflux(aicen      (n),                 &
                                flatn_f    (n), fsensn_f   (n), &
                                fcondtopn_f(n),                 &
                                fsurfn_f   (n),                 &
                                flatn      (n), fsensn     (n), &
                                fsurfn     (n),                 &
                                fcondtopn  (n))
            endif

            call thermo_vertical(nilyr,        nslyr,        &
                                 dt,           aicen    (n), &
                                 vicen    (n), vsnon    (n), &
                                 Tsfc     (n), zSin   (:,n), &
                                 zqin   (:,n), zqsn   (:,n), &
                                 apnd     (n), hpnd     (n), &
                                 iage     (n), tr_pond_topo, &
                                 flw,          potT,         &
                                 Qa,           rhoa,         &
                                 fsnow,        fpond,        &
                                 fbot,         Tbot,         &
                                 sss,                        &
                                 lhcoef,       shcoef,       &
                                 fswsfcn  (n), fswintn  (n), &
                                 Sswabsn(:,n), Iswabsn(:,n), &
                                 fsurfn   (n), fcondtopn(n), &
                                 fsensn   (n), flatn    (n), &
                                 flwoutn,      evapn,        &
                                 freshn,       fsaltn,       &
                                 fhocnn,                     &
                                 melttn   (n), meltsn   (n), &
                                 meltbn   (n),               &
                                 congeln  (n), snoicen  (n), &
                                 mlt_onset,    frz_onset,    &
                                 yday,         dsnown   (n), &
                                 l_stop,       stop_label,   &
                                 prescribed_ice)
               
            if (l_stop) then
               stop_label = 'ice: Vertical thermo error: '//trim(stop_label)
               return
            endif
               
      !-----------------------------------------------------------------
      ! Total absorbed shortwave radiation
      !-----------------------------------------------------------------

            fswabsn = fswsfcn(n) + fswintn(n) + fswthrun(n)

      !-----------------------------------------------------------------
      ! Aerosol update
      !-----------------------------------------------------------------

            if (tr_aero) then
               call update_aerosol (dt,                             &
                                    nilyr, nslyr, n_aero,           &
                                    melttn     (n), meltsn     (n), &
                                    meltbn     (n), congeln    (n), &
                                    snoicen    (n), fsnow,          &
                                    aerosno(:,:,n), aeroice(:,:,n), &
                                    aicen_init (n), vicen_init (n), &
                                    vsnon_init (n),                 &
                                    vicen      (n), vsnon      (n), &
                                    aicen      (n),                 &
                                    faero_atm     ,  faero_ocn)
            endif

         endif   ! aicen_init

      !-----------------------------------------------------------------
      ! Melt ponds
      ! If using tr_pond_cesm, the full calculation is performed here.
      ! If using tr_pond_topo, the rest of the calculation is done after
      ! the surface fluxes are merged, below.
      !-----------------------------------------------------------------

         !call ice_timer_start(timer_ponds)
         if (tr_pond) then
               
            if (tr_pond_cesm) then
               rfrac = rfracmin + (rfracmax-rfracmin) * aicen(n) 
               call compute_ponds_cesm(dt,        hi_min,    &
                                       pndaspect, rfrac,     &
                                       melttn(n), meltsn(n), &
                                       frain,                &
                                       aicen (n), vicen (n), &
                                       vsnon (n), Tsfc  (n), &
                                       apnd  (n), hpnd  (n))
                  
            elseif (tr_pond_lvl) then
               rfrac = rfracmin + (rfracmax-rfracmin) * aicen(n)
               call compute_ponds_lvl(dt,        nilyr,     &
                                      ktherm,               &
                                      hi_min,               &
                                      dpscale,   frzpnd,    &
                                      pndaspect, rfrac,     &
                                      melttn(n), meltsn(n), &
                                      frain,     Tair,      &
                                      fsurfn(n),            &
                                      dhsn  (n), ffracn(n), &
                                      aicen (n), vicen (n), &
                                      vsnon (n),            &
                                      zqin(:,n), zSin(:,n), &
                                      Tsfc  (n), alvl  (n), &
                                      apnd  (n), hpnd  (n), &
                                      ipnd  (n))
                  
            elseif (tr_pond_topo) then
               if (aicen_init(n) > puny) then
                     
                  ! collect liquid water in ponds
                  ! assume salt still runs off
                  rfrac = rfracmin + (rfracmax-rfracmin) * aicen(n)
                  pond = rfrac/rhofresh * (melttn(n)*rhoi &
                       +                   meltsn(n)*rhos &
                       +                   frain *dt)

                  ! if pond does not exist, create new pond over full ice area
                  ! otherwise increase pond depth without changing pond area
                  if (apnd(n) < puny) then
                     hpnd(n) = c0
                     apnd(n) = c1
                  endif
                  hpnd(n) = (pond + hpnd(n)*apnd(n)) / apnd(n)
                  fpond = fpond + pond * aicen(n) ! m
               endif ! aicen_init
            endif

         endif ! tr_pond
         !call ice_timer_stop(timer_ponds)

      !-----------------------------------------------------------------
      ! Increment area-weighted fluxes.
      !-----------------------------------------------------------------

         if (aicen_init(n) > puny) &
            call merge_fluxes (aicen_init(n),            &
                               flw,        coszen,       & 
                               strairxn,   strairyn,     &
                               Cdn_atm_ratio_n,          &
                               fsurfn(n),  fcondtopn(n), &
                               fsensn(n),  flatn(n),     &
                               fswabsn,    flwoutn,      &
                               evapn,                    &
                               Trefn,      Qrefn,        &
                               freshn,     fsaltn,       &
                               fhocnn,     fswthrun(n),  &
                               strairxT,   strairyT,     &
                               Cdn_atm_ratio,            &
                               fsurf,      fcondtop,     &
                               fsens,      flat,         &
                               fswabs,     flwout,       &
                               evap,                     &
                               Tref,       Qref,         &
                               fresh,      fsalt,        &
                               fhocn,      fswthru,      &
                               melttn (n), meltsn(n),    &
                               meltbn (n), congeln(n),   &
                               snoicen(n),               &
                               meltt,      melts,        &
                               meltb,      congel,       &
                               snoice,                   &
                               Uref,       Urefn)

      enddo                  ! ncat

      !-----------------------------------------------------------------
      ! Calculate ponds from the topographic scheme
      !-----------------------------------------------------------------
      !call ice_timer_start(timer_ponds)
      if (tr_pond_topo) then
         call compute_ponds_topo(dt,       ncat,      nilyr,     &
                                 ktherm,   heat_capacity,        &
                                 aice,     aicen,                &
                                 vice,     vicen,                &
                                 vsno,     vsnon,                &
                                 potT,     meltt,                &
                                 fsurf,    fpond,                &
                                 Tsfc,     Tf,                   &
                                 zqin,     zSin,                 &
                                 apnd,     hpnd,      ipnd,      &
                                 l_stop,   stop_label)
      endif
      !call ice_timer_stop(timer_ponds)

      end subroutine colpkg_step_therm1

!=======================================================================
! Driver for thermodynamic changes not needed for coupling:
! transport in thickness space, lateral growth and melting.
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL

      subroutine colpkg_step_therm2 (dt, ncat, n_aero, nltrcr,           &
                                     nilyr,        nslyr,         &
                                     hin_max,      nblyr,         &
                                     aicen,                       &
                                     vicen,        vsnon,         &
                                     aicen_init,   vicen_init,    &
                                     trcrn,                       &
                                     aice0,        aice,          &
                                     trcr_depend,                 &
                                     trcr_base,    n_trcr_strata, &
                                     nt_strata,                   &
                                     Tf,           sss,           &
                                     salinz,                      &
                                     rside,        meltl,         &
                                     frzmlt,       frazil,        &
                                     frain,        fpond,         &
                                     fresh,        fsalt,         &
                                     fhocn,        update_ocn_f,  &
                                     bgrid,        cgrid,         &
                                     igrid,        faero_ocn,     &
                                     first_ice,    fzsal,         &
                                     flux_bio,     ocean_bio,     &
                                     l_stop,       stop_label,    &
                                     frazil_diag,                 &
                                     frz_onset,    yday)

      use ice_constants_colpkg, only: puny
      use ice_itd, only: aggregate_area, reduce_area, cleanup_itd
      use ice_therm_itd, only: linear_itd, add_new_ice, lateral_melt
      use ice_colpkg_tracers, only: ntrcr, nbtrcr, tr_aero, tr_pond_topo

      integer (kind=int_kind), intent(in) :: &
         ncat     , & ! number of thickness categories
         nltrcr   , & ! number of zbgc tracers
         nblyr    , & ! number of bio layers
         nilyr    , & ! number of ice layers
         nslyr    , & ! number of snow layers
         n_aero       ! number of aerosol tracers

      logical (kind=log_kind), intent(in) :: &
         update_ocn_f     ! if true, update fresh water and salt fluxes

      real (kind=dbl_kind), dimension(0:ncat), intent(inout) :: &
         hin_max      ! category boundaries (m)

      real (kind=dbl_kind), intent(in) :: &
         dt       , & ! time step
         Tf       , & ! freezing temperature (C)
         sss      , & ! sea surface salinity (ppt)
         rside    , & ! fraction of ice that melts laterally
         frzmlt       ! freezing/melting potential (W/m^2)

      integer (kind=int_kind), dimension (:), intent(in) :: &
         trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  ! number of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcr_base      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nt_strata      ! indices of underlying tracer layers

      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         bgrid              ! biology nondimensional vertical grid points

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid              ! biology vertical interface points
 
      real (kind=dbl_kind), dimension (nilyr+1), intent(in) :: &
         cgrid              ! CICE vertical coordinate   

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         salinz   , & ! initial salinity profile
         ocean_bio    ! ocean concentration of biological tracer

      real (kind=dbl_kind), intent(inout) :: &
         aice     , & ! sea ice concentration
         aice0    , & ! concentration of open water
         frain    , & ! rainfall rate (kg/m^2 s)
         fpond    , & ! fresh water flux to ponds (kg/m^2/s)
         fresh    , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt    , & ! salt flux to ocean (kg/m^2/s)
         fhocn    , & ! net heat flux to ocean (W/m^2)
         fzsal    , & ! salt flux to ocean from zsalinity (kg/m^2/s)
         meltl    , & ! lateral ice melt         (m/step-->cm/day)
         frazil   , & ! frazil ice growth        (m/step-->cm/day)
         frazil_diag  ! frazil ice growth diagnostic (m/step-->cm/day)

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         aicen_init,& ! initial concentration of ice
         vicen_init,& ! initial volume per unit area of ice          (m)
         aicen    , & ! concentration of ice
         vicen    , & ! volume per unit area of ice          (m)
         vsnon    , & ! volume per unit area of snow         (m)
         faero_ocn, & ! aerosol flux to ocean  (kg/m^2/s)
         flux_bio     ! all bio fluxes to ocean

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         trcrn        ! tracers
 
      logical (kind=log_kind), dimension(:), intent(inout) :: &
         first_ice      ! true until ice forms

      logical (kind=log_kind), intent(out) :: &
         l_stop         ! if true, abort model

      character (len=*), intent(out) :: stop_label

      real (kind=dbl_kind), intent(inout), optional :: &
         frz_onset    ! day of year that freezing begins (congel or frazil)

      real (kind=dbl_kind), intent(in), optional :: &
         yday         ! day of year

      l_stop = .false.

      !-----------------------------------------------------------------
      ! Let rain drain through to the ocean.
      !-----------------------------------------------------------------

      fresh  = fresh + frain * aice

      !-----------------------------------------------------------------
      ! Given thermodynamic growth rates, transport ice between
      ! thickness categories.
      !-----------------------------------------------------------------

!      call ice_timer_start(timer_catconv)    ! category conversions

      !-----------------------------------------------------------------
      ! Compute fractional ice area in each grid cell.
      !-----------------------------------------------------------------

      call aggregate_area (ncat, aicen, aice, aice0)

      if (kitd == 1) then

      !-----------------------------------------------------------------
      ! Identify grid cells with ice.
      !-----------------------------------------------------------------

         if (aice > puny) then

            call linear_itd (ncat,     hin_max,        &
                             nilyr,    nslyr,          &
                             ntrcr,    trcr_depend,    &
                             trcr_base,        & 
                             n_trcr_strata,   &
                             nt_strata,                &
                             aicen_init,            &
                             vicen_init,            &
                             aicen,                 &
                             trcrn,           & 
                             vicen,                 &
                             vsnon,                 &
                             aice      ,         &
                             aice0     ,         &
                             fpond,       l_stop,      &
                             stop_label)

            if (l_stop) return

         endif ! aice > puny

      endif  ! kitd = 1

!      call ice_timer_stop(timer_catconv)    ! category conversions

      !-----------------------------------------------------------------
      ! Add frazil ice growing in leads.
      !-----------------------------------------------------------------

      ! identify ice-ocean cells

         call add_new_ice (ncat,          nilyr,        &
                           nblyr,                       &
                           n_aero,        dt,           &
                           ntrcr,         nltrcr,       &
                           hin_max,       ktherm,       &
                           aicen,         trcrn,        &
                           vicen,         vsnon(1),     &
                           aice0,         aice,         &
                           frzmlt,        frazil,       &
                           frz_onset,     yday,         &
                           update_ocn_f,                &
                           fresh,         fsalt,        &
                           Tf,            sss,          &
                           salinz,        phi_init,     &
                           dSin0_frazil,  bgrid,        &
                           cgrid,         igrid,        &
                           nbtrcr,        flux_bio,     &
                           ocean_bio,     fzsal,        &
                           frazil_diag,                 &
                           l_stop,        stop_label)

         if (l_stop) return

      !-----------------------------------------------------------------
      ! Melt ice laterally.
      !-----------------------------------------------------------------

      call lateral_melt (dt,        ncat,          &
                         nilyr,     nslyr,         &
                         n_aero,    fpond,         &
                         fresh,     fsalt,         &
                         fhocn,     faero_ocn,     &
                         rside,     meltl,         &
                         aicen,     vicen,         &
                         vsnon,     trcrn,         &
                         fzsal,     flux_bio,      &
                         nbtrcr,    nblyr)

      !-----------------------------------------------------------------
      ! For the special case of a single category, adjust the area and
      ! volume (assuming that half the volume change decreases the
      ! thickness, and the other half decreases the area).  
      !-----------------------------------------------------------------

!echmod: test this
      if (ncat==1) &
          call reduce_area (hin_max   (0),                &
                            aicen     (1), vicen     (1), &
                            aicen_init(1), vicen_init(1))

      !-----------------------------------------------------------------
      ! ITD cleanup: Rebin thickness categories if necessary, and remove
      !  categories with very small areas.
      !-----------------------------------------------------------------

      call cleanup_itd (dt,                   ntrcr,            &
                        nilyr,                nslyr,            &
                        ncat,                 hin_max,          &
                        aicen,                trcrn(1:ntrcr,:), &
                        vicen,                vsnon,            &
                        aice0,                aice,             & 
                        n_aero,                                 &
                        nbtrcr,               nblyr,            &
                        l_stop,               stop_label,       &
                        tr_aero,                                &
                        tr_pond_topo,         heat_capacity,    &
                        first_ice,                              &
                        trcr_depend,          trcr_base,        &
                        n_trcr_strata,        nt_strata,        &
                        fpond,                fresh,            &
                        fsalt,                fhocn,            &
                        faero_ocn,            fzsal,            &
                        flux_bio)   

      end subroutine colpkg_step_therm2

!=======================================================================
!
! Scales radiation fields computed on the previous time step.
!
! authors: Elizabeth Hunke, LANL

      subroutine colpkg_prep_radiation (ncat, nilyr, nslyr,    &
                                        aice,        aicen,    &
                                        swvdr,       swvdf,    &
                                        swidr,       swidf,    &
                                        alvdr_ai,    alvdf_ai, &
                                        alidr_ai,    alidf_ai, &
                                        scale_factor,          &
                                        fswsfcn,     fswintn,  &
                                        fswthrun,    fswpenln, &
                                        Sswabsn,     Iswabsn)

      use ice_constants_colpkg, only: c0, c1, puny

      integer (kind=int_kind), intent(in) :: &
         ncat    , & ! number of ice thickness categories
         nilyr   , & ! number of ice layers
         nslyr       ! number of snow layers

      real (kind=dbl_kind), intent(in) :: &
         aice        , & ! ice area fraction
         swvdr       , & ! sw down, visible, direct  (W/m^2)
         swvdf       , & ! sw down, visible, diffuse (W/m^2)
         swidr       , & ! sw down, near IR, direct  (W/m^2)
         swidf       , & ! sw down, near IR, diffuse (W/m^2)
         ! grid-box-mean albedos aggregated over categories (if calc_Tsfc)
         alvdr_ai    , & ! visible, direct   (fraction)
         alidr_ai    , & ! near-ir, direct   (fraction)
         alvdf_ai    , & ! visible, diffuse  (fraction)
         alidf_ai        ! near-ir, diffuse  (fraction)

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         aicen           ! ice area fraction in each category

      real (kind=dbl_kind), intent(inout) :: &
         scale_factor    ! shortwave scaling factor, ratio new:old

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         fswsfcn     , & ! SW absorbed at ice/snow surface (W m-2)
         fswintn     , & ! SW absorbed in ice interior, below surface (W m-2)
         fswthrun        ! SW through ice to ocean (W/m^2)

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         fswpenln    , & ! visible SW entering ice layers (W m-2)
         Iswabsn     , & ! SW radiation absorbed in ice layers (W m-2)
         Sswabsn         ! SW radiation absorbed in snow layers (W m-2)

      ! local variables

      integer (kind=int_kind) :: &
         k           , & ! vertical index       
         n               ! thickness category index

      real (kind=dbl_kind) :: netsw 

      !-----------------------------------------------------------------
      ! Compute netsw scaling factor (new netsw / old netsw)
      !-----------------------------------------------------------------

         if (aice > c0 .and. scale_factor > puny) then
            netsw = swvdr*(c1 - alvdr_ai) &
                  + swvdf*(c1 - alvdf_ai) &
                  + swidr*(c1 - alidr_ai) &
                  + swidf*(c1 - alidf_ai)
            scale_factor = netsw / scale_factor
         else
            scale_factor = c1
         endif

         do n = 1, ncat

            if (aicen(n) > puny) then

      !-----------------------------------------------------------------
      ! Scale absorbed solar radiation for change in net shortwave
      !-----------------------------------------------------------------

               fswsfcn(n)  = scale_factor*fswsfcn (n)
               fswintn(n)  = scale_factor*fswintn (n)
               fswthrun(n) = scale_factor*fswthrun(n)
               do k = 1,nilyr+1
                  fswpenln(k,n) = scale_factor*fswpenln(k,n)
               enddo       !k
               do k=1,nslyr
                  Sswabsn(k,n) = scale_factor*Sswabsn(k,n)
               enddo
               do k=1,nilyr
                  Iswabsn(k,n) = scale_factor*Iswabsn(k,n)
               enddo

            endif
         enddo                  ! ncat

      end subroutine colpkg_prep_radiation

!=======================================================================
!
! Computes radiation fields
!
! authors: William H. Lipscomb, LANL
!          David Bailey, NCAR
!          Elizabeth C. Hunke, LANL

      subroutine colpkg_step_radiation (dt,       ncat,      & 
                                        n_algae,  tr_zaero,  &
                                        nblyr,    ntrcr,     &
                                        nbtrcr,   nbtrcr_sw, &
                                        nilyr,    nslyr,     &
                                        n_aero,   n_zaero,   &
                                        dEdd_algae,          &
                                        nlt_chl_sw,          &
                                        nlt_zaero_sw,        &
                                        swgrid,   igrid,     &
                                        fbri,                &
                                        aicen,    vicen,     &
                                        vsnon,    Tsfcn,     &
                                        alvln,    apndn,     &
                                        hpndn,    ipndn,     &
                                        aeron,               &
                                        zbion,               &
                                        trcrn,               &
                                        TLAT,     TLON,      &
                                        calendar_type,       &
                                        days_per_year,       &
                                        nextsw_cday,         &
                                        yday,     sec,       &
                                        kaer_tab, waer_tab,  &
                                        gaer_tab,            &
                                        kaer_bc_tab,         &
                                        waer_bc_tab,         &
                                        gaer_bc_tab,         &
                                        bcenh,               &
                                        modal_aero,          &
                                        swvdr,    swvdf,     &
                                        swidr,    swidf,     &
                                        coszen,   fsnow,     &
                                        alvdrn,   alvdfn,    &
                                        alidrn,   alidfn,    &
                                        fswsfcn,  fswintn,   &
                                        fswthrun, fswpenln,  &
                                        Sswabsn,  Iswabsn,   &
                                        albicen,  albsnon,   &
                                        albpndn,  apeffn,    &
                                        snowfracn,           &
                                        dhsn,     ffracn,    &
                                        l_print_point, &
                                        initonly)

      use ice_constants_colpkg, only: c0, puny
      use ice_shortwave, only: run_dEdd, shortwave_ccsm3, compute_shortwave_trcr
      use ice_colpkg_tracers, only: tr_pond_cesm, tr_pond_lvl, tr_pond_topo, &
                                    tr_bgc_N, tr_aero
      use ice_colpkg_shared, only:  z_tracers, skl_bgc

      integer (kind=int_kind), intent(in) :: &
         ncat      , & ! number of ice thickness categories
         nilyr     , & ! number of ice layers
         nslyr     , & ! number of snow layers
         n_aero    , & ! number of aerosols
         n_zaero   , & ! number of zaerosols 
         nlt_chl_sw, & ! index for chla
         nblyr     , &
         ntrcr     , &
         nbtrcr    , &
         nbtrcr_sw , &
         n_algae

      integer (kind=int_kind), dimension(:), intent(in) :: &
        nlt_zaero_sw   ! index for zaerosols

      real (kind=dbl_kind), intent(in) :: &
         dt        , & ! time step (s)
         swvdr     , & ! sw down, visible, direct  (W/m^2)
         swvdf     , & ! sw down, visible, diffuse (W/m^2)
         swidr     , & ! sw down, near IR, direct  (W/m^2)
         swidf     , & ! sw down, near IR, diffuse (W/m^2)
         fsnow     , & ! snowfall rate (kg/m^2 s)
         TLAT, TLON    ! latitude and longitude (radian)

      character (len=char_len), intent(in) :: &
         calendar_type       ! differentiates Gregorian from other calendars

      integer (kind=int_kind), intent(in) :: &
         days_per_year, &    ! number of days in one year
         sec                 ! elapsed seconds into date

      real (kind=dbl_kind), intent(in) :: &
         nextsw_cday     , & ! julian day of next shortwave calculation
         yday                ! day of the year

      real (kind=dbl_kind), intent(inout) :: &
         coszen        ! cosine solar zenith angle, < 0 for sun below horizon 

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         igrid              ! biology vertical interface points
 
      real (kind=dbl_kind), dimension (:), intent(in) :: &
         swgrid                ! grid for ice tracers used in dEdd scheme
        
      real (kind=dbl_kind), dimension(:,:), intent(in) :: & 
         kaer_tab, & ! aerosol mass extinction cross section (m2/kg)
         waer_tab, & ! aerosol single scatter albedo (fraction)
         gaer_tab    ! aerosol asymmetry parameter (cos(theta))

      real (kind=dbl_kind), dimension(:,:), intent(in) :: & 
         kaer_bc_tab, & ! aerosol mass extinction cross section (m2/kg)
         waer_bc_tab, & ! aerosol single scatter albedo (fraction)
         gaer_bc_tab    ! aerosol asymmetry parameter (cos(theta))

      real (kind=dbl_kind), dimension(:,:,:), intent(in) :: & 
         bcenh 

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         aicen     , & ! ice area fraction in each category
         vicen     , & ! ice volume in each category (m)
         vsnon     , & ! snow volume in each category (m)
         Tsfcn     , & ! surface temperature (deg C)
         alvln     , & ! level-ice area fraction
         apndn     , & ! pond area fraction
         hpndn     , & ! pond depth (m)
         ipndn     , & ! pond refrozen lid thickness (m)
         fbri           ! brine fraction 

      real(kind=dbl_kind), dimension(:,:), intent(in) :: &
         aeron     , & ! aerosols (kg/m^3)
         trcrn         ! tracers

      real(kind=dbl_kind), dimension(:,:), intent(inout) :: &
         zbion         ! zaerosols (kg/m^3) and chla (mg/m^3)

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         alvdrn    , & ! visible, direct  albedo (fraction)
         alidrn    , & ! near-ir, direct   (fraction)
         alvdfn    , & ! visible, diffuse  (fraction)
         alidfn    , & ! near-ir, diffuse  (fraction)
         fswsfcn   , & ! SW absorbed at ice/snow surface (W m-2)
         fswintn   , & ! SW absorbed in ice interior, below surface (W m-2)
         fswthrun  , & ! SW through ice to ocean (W/m^2)
         snowfracn , & ! snow fraction on each category
         dhsn      , & ! depth difference for snow on sea ice and pond ice
         ffracn    , & ! fraction of fsurfn used to melt ipond
                       ! albedo components for history
         albicen   , & ! bare ice 
         albsnon   , & ! snow 
         albpndn   , & ! pond 
         apeffn        ! effective pond area used for radiation calculation

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         fswpenln  , & ! visible SW entering ice layers (W m-2)
         Iswabsn   , & ! SW radiation absorbed in ice layers (W m-2)
         Sswabsn       ! SW radiation absorbed in snow layers (W m-2)

      logical (kind=log_kind), intent(in) :: &
         l_print_point, & ! flag for printing diagnostics
         dEdd_algae   , & ! .true. use prognostic chla in dEdd
         modal_aero   , & ! .true. use modal aerosol optical treatment
         tr_zaero

      logical (kind=log_kind), optional :: &
         initonly         ! flag to indicate init only, default is false

      ! local variables

      integer (kind=int_kind) :: &
         n                  ! thickness category index

      logical (kind=log_kind) :: &
         l_stop      ,&  ! if true, abort the model
         linitonly       ! local flag for initonly

      character (char_len) :: stop_label

      real(kind=dbl_kind) :: &
        hin,         & ! Ice thickness (m)
        hbri           ! brine thickness (m)

        hin = c0
        hbri = c0
        linitonly = .false.
        if (present(initonly)) then
           linitonly = initonly
        endif

         ! Initialize
         do n = 1, ncat
            alvdrn  (n) = c0
            alidrn  (n) = c0
            alvdfn  (n) = c0
            alidfn  (n) = c0
            fswsfcn (n) = c0
            fswintn (n) = c0
            fswthrun(n) = c0
         enddo   ! ncat
         fswpenln (:,:) = c0
         Iswabsn  (:,:) = c0
         Sswabsn  (:,:) = c0
         zbion(:,:) = c0

         ! Interpolate z-shortwave tracers to shortwave grid
         if (dEdd_algae) then
         do n = 1, ncat
              if (aicen(n) .gt. puny) then
                 hin = vicen(n)/aicen(n)
                 hbri= fbri(n)*hin
                 call compute_shortwave_trcr(n_algae, nslyr,  &
                                     trcrn(1:ntrcr,n),        &
                                     zbion(1:nbtrcr_sw,n),    &
                                     swgrid,       hin,       &
                                     hbri,         ntrcr,     &
                                     nilyr,        nblyr,     &
                                     igrid,                   &
                                     nbtrcr_sw,    n_zaero,   &
                                     skl_bgc,      z_tracers, &
                                     l_stop,       stop_label)
              endif
         enddo
         endif

         if (calc_Tsfc) then
         if (trim(shortwave) == 'dEdd') then ! delta Eddington
            
            call run_dEdd(dt,           tr_aero,        &
                          tr_pond_cesm,                 &
                          tr_pond_lvl,                  &
                          tr_pond_topo,                 &
                          ncat,         n_aero,         &
                          n_zaero,      dEdd_algae,     &
                          nlt_chl_sw,   nlt_zaero_sw,   &
                          tr_bgc_N,     tr_zaero,       &
                          nilyr,        nslyr,          &
                          aicen,        vicen,          &
                          vsnon,        Tsfcn,          &
                          alvln,        apndn,          &
                          hpndn,        ipndn,          &
                          aeron,        kalg,           &
                          zbion,                        &
                          heat_capacity,                &
                          TLAT,         TLON,           &
                          calendar_type,days_per_year,  &
                          nextsw_cday,  yday,           &
                          sec,          R_ice,          &
                          R_pnd,        R_snw,          &
                          dT_mlt,       rsnw_mlt,       &
                          hs0,          hs1,            &
                          hp1,          pndaspect,      &
                          kaer_tab,     waer_tab,       &
                          gaer_tab,                     &
                          kaer_bc_tab,                  &
                          waer_bc_tab,                  &
                          gaer_bc_tab,                  &
                          bcenh,                        &
                          modal_aero,                   &
                          swvdr,        swvdf,          &
                          swidr,        swidf,          &
                          coszen,       fsnow,          &
                          alvdrn,       alvdfn,         &
                          alidrn,       alidfn,         &
                          fswsfcn,      fswintn,        &
                          fswthrun,     fswpenln,       &
                          Sswabsn,      Iswabsn,        &
                          albicen,      albsnon,        &
                          albpndn,      apeffn,         &
                          snowfracn,                    &
                          dhsn,         ffracn,         &
                          l_print_point,                &
                          linitonly)
 
         else  ! .not. dEdd

            call shortwave_ccsm3(aicen,      vicen,      &
                                 vsnon,                  &
                                 Tsfcn,                  &
                                 swvdr,      swvdf,      &
                                 swidr,      swidf,      &
                                 heat_capacity,          &
                                 albedo_type,            &
                                 albicev,    albicei,    &
                                 albsnowv,   albsnowi,   &
                                 ahmax,                  &
                                 alvdrn,     alidrn,     &
                                 alvdfn,     alidfn,     &
                                 fswsfcn,    fswintn,    &
                                 fswthrun,               &
                                 fswpenln,               &
                                 Iswabsn,                &
                                 Sswabsn,                &
                                 albicen,    albsnon,    &
                                 coszen,     ncat)

         endif   ! shortwave

      else    ! .not. calc_Tsfc

      ! Calculate effective pond area for HadGEM

         if (tr_pond_topo) then
            do n = 1, ncat
               apeffn(n) = c0 
               if (aicen(n) > puny) then
               ! Lid effective if thicker than hp1
                 if (apndn(n)*aicen(n) > puny .and. ipndn(n) < hp1) then
                    apeffn(n) = apndn(n)
                 else
                    apeffn(n) = c0
                 endif
                 if (apndn(n) < puny) apeffn(n) = c0
               endif
            enddo  ! ncat
 
         endif ! tr_pond_topo

         ! Initialize for safety
         do n = 1, ncat
            alvdrn(n) = c0
            alidrn(n) = c0
            alvdfn(n) = c0
            alidfn(n) = c0
            fswsfcn(n) = c0
            fswintn(n) = c0
            fswthrun(n) = c0
         enddo   ! ncat
         Iswabsn(:,:) = c0
         Sswabsn(:,:) = c0

      endif    ! calc_Tsfc

      end subroutine colpkg_step_radiation

!=======================================================================
!
! Computes sea ice mechanical deformation
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL

      subroutine colpkg_step_ridge (dt,           ndtd,          &
                                    nilyr,        nslyr,         &
                                    nblyr,                       &
                                    ncat,         hin_max,       &
                                    rdg_conv,     rdg_shear,     &
                                    aicen,                       &
                                    trcrn,                       &
                                    vicen,        vsnon,         &
                                    aice0,        trcr_depend,   &
                                    trcr_base,    n_trcr_strata, &
                                    nt_strata,                   &
                                    dardg1dt,     dardg2dt,      &
                                    dvirdgdt,     opening,       &
                                    fpond,                       &
                                    fresh,        fhocn,         &
                                    n_aero,                      &
                                    faero_ocn,                   &
                                    aparticn,     krdgn,         &
                                    aredistn,     vredistn,      &
                                    dardg1ndt,    dardg2ndt,     &
                                    dvirdgndt,                   &
                                    araftn,       vraftn,        &
                                    aice,         fsalt,         &
                                    first_ice,    fzsal,         &
                                    flux_bio,                    &
                                    l_stop,       stop_label)

      use ice_mechred, only: ridge_ice
      use ice_itd, only: cleanup_itd
      use ice_colpkg_tracers, only: tr_pond_topo, tr_aero, tr_brine, ntrcr, nbtrcr

      real (kind=dbl_kind), intent(in) :: &
         dt           ! time step

      integer (kind=int_kind), intent(in) :: &
         ncat  , & ! number of thickness categories
         ndtd  , & ! number of dynamics supercycles
         nblyr , & ! number of bio layers
         nilyr , & ! number of ice layers
         nslyr , & ! number of snow layers
         n_aero    ! number of aerosol tracers

      real (kind=dbl_kind), dimension(0:ncat), intent(inout) :: &
         hin_max   ! category limits (m)

      integer (kind=int_kind), dimension (:), intent(in) :: &
         trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  ! number of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcr_base      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nt_strata      ! indices of underlying tracer layers

      real (kind=dbl_kind), intent(inout) :: &
         aice     , & ! sea ice concentration
         aice0    , & ! concentration of open water
         rdg_conv , & ! convergence term for ridging (1/s)
         rdg_shear, & ! shear term for ridging (1/s)
         dardg1dt , & ! rate of area loss by ridging ice (1/s)
         dardg2dt , & ! rate of area gain by new ridges (1/s)
         dvirdgdt , & ! rate of ice volume ridged (m/s)
         opening  , & ! rate of opening due to divergence/shear (1/s)
         fpond    , & ! fresh water flux to ponds (kg/m^2/s)
         fresh    , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt    , & ! salt flux to ocean (kg/m^2/s)
         fhocn    , & ! net heat flux to ocean (W/m^2)
         fzsal        ! zsalinity flux to ocean(kg/m^2/s)

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         aicen    , & ! concentration of ice
         vicen    , & ! volume per unit area of ice          (m)
         vsnon    , & ! volume per unit area of snow         (m)
         dardg1ndt, & ! rate of area loss by ridging ice (1/s)
         dardg2ndt, & ! rate of area gain by new ridges (1/s)
         dvirdgndt, & ! rate of ice volume ridged (m/s)
         aparticn , & ! participation function
         krdgn    , & ! mean ridge thickness/thickness of ridging ice
         araftn   , & ! rafting ice area
         vraftn   , & ! rafting ice volume 
         aredistn , & ! redistribution function: fraction of new ridge area
         vredistn , & ! redistribution function: fraction of new ridge volume
         faero_ocn, & ! aerosol flux to ocean  (kg/m^2/s)
         flux_bio     ! all bio fluxes to ocean

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         trcrn        ! tracers

      !logical (kind=log_kind), intent(in) :: &
         !tr_pond_topo,& ! if .true., use explicit topography-based ponds
         !tr_aero     ,& ! if .true., use aerosol tracers
         !tr_brine    !,& ! if .true., brine height differs from ice thickness
         !heat_capacity  ! if true, ice has nonzero heat capacity

      logical (kind=log_kind), dimension(:), intent(inout) :: &
         first_ice    ! true until ice forms

      logical (kind=log_kind), intent(out) :: &
         l_stop       ! if true, abort the model

      character (len=*), intent(out) :: &
         stop_label   ! diagnostic information for abort

      ! local variables

      real (kind=dbl_kind) :: &
         dtt      , & ! thermo time step
         atmp     , & ! temporary ice area
         atmp0        ! temporary open water area

      l_stop = .false.

      !-----------------------------------------------------------------
      ! Identify ice-ocean cells.
      ! Note:  We can not limit the loop here using aice>puny because
      !        aice has not yet been updated since the transport (and
      !        it may be out of whack, which the ridging helps fix).-ECH
      !-----------------------------------------------------------------
           
         call ridge_ice (dt,           ndtd,           &
                         ncat,         n_aero,         &
                         nilyr,        nslyr,          &
                         ntrcr,        hin_max,        &
                         rdg_conv,     rdg_shear,      &
                         aicen,                        &
                         trcrn,                        &
                         vicen,        vsnon,          &
                         aice0,                        &
                         trcr_depend,                  &
                         trcr_base,                    &
                         n_trcr_strata,                &
                         nt_strata,                    &
                         l_stop,                       &
                         stop_label,                   &
                         krdg_partic, krdg_redist, &
                         mu_rdg,                   &
                         dardg1dt,     dardg2dt,       &
                         dvirdgdt,     opening,        &
                         fpond,                        &
                         fresh,        fhocn,          &
                         tr_brine,     faero_ocn,      &
                         aparticn,     krdgn,          &
                         aredistn,     vredistn,       &
                         dardg1ndt,    dardg2ndt,      &
                         dvirdgndt,                    &
                         araftn,       vraftn)        

         if (l_stop) return

      !-----------------------------------------------------------------
      ! ITD cleanup: Rebin thickness categories if necessary, and remove
      !  categories with very small areas.
      !-----------------------------------------------------------------

      dtt = dt * ndtd  ! for proper averaging over thermo timestep
      call cleanup_itd (dtt,                  ntrcr,            &
                        nilyr,                nslyr,            &
                        ncat,                 hin_max,          &
                        aicen,                trcrn,            &
                        vicen,                vsnon,            &
                        aice0,                aice,             &          
                        n_aero,                                 &
                        nbtrcr,               nblyr,            &
                        l_stop,               stop_label,       &
                        tr_aero,                                &
                        tr_pond_topo,         heat_capacity,    &  
                        first_ice,                              &                
                        trcr_depend,          trcr_base,        &
                        n_trcr_strata,        nt_strata,        &
                        fpond,                fresh,            &
                        fsalt,                fhocn,            &
                        faero_ocn,            fzsal,            &
                        flux_bio)

      if (l_stop) then
         stop_label = 'ice: ITD cleanup error in colpkg_step_ridge'
      endif

      end subroutine colpkg_step_ridge

!=======================================================================

! Aggregate ice state variables over thickness categories.
!
! authors: C. M. Bitz, UW
!          W. H. Lipscomb, LANL

      subroutine colpkg_aggregate (ncat,               &
                                   aicen,    trcrn,    &
                                   vicen,    vsnon,    &
                                   aice,     trcr,     &
                                   vice,     vsno,     &
                                   aice0,              &
                                   ntrcr,              &
                                   trcr_depend,        &
                                   trcr_base,          & 
                                   n_trcr_strata,      &
                                   nt_strata)

      use ice_constants_colpkg, only: c0, c1
      use ice_colpkg_tracers, only: colpkg_compute_tracers

      integer (kind=int_kind), intent(in) :: &
         ncat  , & ! number of thickness categories
         ntrcr     ! number of tracers in use

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (:,:), &
         intent(inout) :: &
         trcrn     ! ice tracers

      integer (kind=int_kind), dimension (:), intent(in) :: &
         trcr_depend, & ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon
         n_trcr_strata  ! number of underlying tracer layers

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcr_base      ! = 0 or 1 depending on tracer dependency
                        ! argument 2:  (1) aice, (2) vice, (3) vsno

      integer (kind=int_kind), dimension (:,:), intent(in) :: &
         nt_strata      ! indices of underlying tracer layers

      real (kind=dbl_kind), intent(out) :: &
         aice  , & ! concentration of ice
         vice  , & ! volume per unit area of ice          (m)
         vsno  , & ! volume per unit area of snow         (m)
         aice0     ! concentration of open water

      real (kind=dbl_kind), dimension (:),  &
         intent(out) :: &
         trcr      ! ice tracers

      ! local variables

      integer (kind=int_kind) :: &
         n, it, itl, & ! loop indices
         ntr           ! tracer index

      real (kind=dbl_kind), dimension (:), allocatable :: &
         atrcr     ! sum of aicen*trcrn or vicen*trcrn or vsnon*trcrn

      real (kind=dbl_kind) :: &
         atrcrn    ! category value

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      aice0 = c1
      aice  = c0
      vice  = c0
      vsno  = c0

      allocate (atrcr(ntrcr))

      !-----------------------------------------------------------------
      ! Aggregate
      !-----------------------------------------------------------------

      atrcr(:) = c0

      do n = 1, ncat

            aice = aice + aicen(n)
            vice = vice + vicen(n)
            vsno = vsno + vsnon(n)

         do it = 1, ntrcr
            atrcrn = trcrn(it,n)*(trcr_base(it,1) * aicen(n) &
                                + trcr_base(it,2) * vicen(n) &
                                + trcr_base(it,3) * vsnon(n))
            if (n_trcr_strata(it) > 0) then  ! additional tracer layers
               do itl = 1, n_trcr_strata(it)
                  ntr = nt_strata(it,itl)
                  atrcrn = atrcrn * trcrn(ntr,n)
               enddo
            endif
            atrcr(it) = atrcr(it) + atrcrn
         enddo                  ! ntrcr
      enddo                     ! ncat

      ! Open water fraction
      aice0 = max (c1 - aice, c0)

      ! Tracers
      call colpkg_compute_tracers (ntrcr,     trcr_depend,   &
                                   atrcr,     aice,          &
                                   vice ,     vsno,          &
                                   trcr_base, n_trcr_strata, &
                                   nt_strata, trcr)   

      deallocate (atrcr)

      end subroutine colpkg_aggregate

!=======================================================================

! Compute the strength of the ice pack, defined as the energy (J m-2)
! dissipated per unit area removed from the ice pack under compression,
! and assumed proportional to the change in potential energy caused
! by ridging.
!
! See Rothrock (1975) and Hibler (1980).
!
! For simpler strength parameterization, see this reference:
! Hibler, W. D. III, 1979: A dynamic-thermodynamic sea ice model,
!  J. Phys. Oceanog., 9, 817-846.
!
! authors: William H. Lipscomb, LANL
!          Elizabeth C. Hunke, LANL

      subroutine colpkg_ice_strength (ncat,               &
                                      aice,     vice,     &
                                      aice0,    aicen,    &
                                      vicen,    &
                                      strength)

      use ice_constants_colpkg, only: p333, c0, c1, c2, Cf, Cp, Pstar, Cstar, &
          rhoi, puny
      use ice_mechred, only: asum_ridging, ridge_itd

      integer (kind=int_kind), intent(in) :: & 
         ncat       ! number of thickness categories

      real (kind=dbl_kind), intent(in) :: &
         aice   , & ! concentration of ice
         vice   , & ! volume per unit area of ice  (m)
         aice0      ! concentration of open water

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         aicen  , & ! concentration of ice
         vicen      ! volume per unit area of ice  (m)

      real (kind=dbl_kind), intent(inout) :: &
         strength   ! ice strength (N/m)

      ! local variables

      real (kind=dbl_kind) :: &
         asum   , & ! sum of ice and open water area
         aksum      ! ratio of area removed to area ridged

      real (kind=dbl_kind), dimension (0:ncat) :: &
         apartic    ! participation function; fraction of ridging
                    ! and closing associated w/ category n

      real (kind=dbl_kind), dimension (ncat) :: &
         hrmin  , & ! minimum ridge thickness
         hrmax  , & ! maximum ridge thickness (krdg_redist = 0)
         hrexp  , & ! ridge e-folding thickness (krdg_redist = 1) 
         krdg       ! mean ridge thickness/thickness of ridging ice

      integer (kind=int_kind) :: &
         n          ! thickness category index

      real (kind=dbl_kind) :: &
         hi     , & ! ice thickness (m)
         h2rdg  , & ! mean value of h^2 for new ridge
         dh2rdg     ! change in mean value of h^2 per unit area
                    ! consumed by ridging 

      if (kstrength == 1) then  ! Rothrock '75 formulation

      !-----------------------------------------------------------------
      ! Compute thickness distribution of ridging and ridged ice.
      !-----------------------------------------------------------------

         call asum_ridging (ncat, aicen, aice0, asum)

         call ridge_itd (ncat,     aice0,      &
                         aicen,    vicen,      &
                         krdg_partic, krdg_redist, &
                         mu_rdg,                   &
                         aksum,    apartic,    &
                         hrmin,    hrmax,      &
                         hrexp,    krdg)   

      !-----------------------------------------------------------------
      ! Compute ice strength based on change in potential energy,
      ! as in Rothrock (1975)
      !-----------------------------------------------------------------

         if (krdg_redist==0) then ! Hibler 1980 formulation

            do n = 1, ncat
               if (aicen(n) > puny .and. apartic(n) > c0)then
                  hi = vicen(n) / aicen(n)
                  h2rdg = p333 * (hrmax(n)**3 - hrmin(n)**3)  &
                               / (hrmax(n) - hrmin(n)) 
                  dh2rdg = -hi*hi + h2rdg/krdg(n)
                  strength = strength + apartic(n) * dh2rdg
               endif         ! aicen > puny
            enddo               ! n

         elseif (krdg_redist==1) then ! exponential formulation

            do n = 1, ncat
               if (aicen(n) > puny .and. apartic(n) > c0) then
                  hi = vicen(n) / aicen(n)
                  h2rdg =    hrmin(n)*hrmin(n) &
                        + c2*hrmin(n)*hrexp(n) &
                        + c2*hrexp(n)*hrexp(n)
                  dh2rdg = -hi*hi + h2rdg/krdg(n)
                  strength = strength + apartic(n) * dh2rdg
               endif
            enddo               ! n

         endif                  ! krdg_redist

         strength = Cf * Cp * strength / aksum
                       ! Cp = (g/2)*(rhow-rhoi)*(rhoi/rhow)
                       ! Cf accounts for frictional dissipation

      else                      ! kstrength /= 1:  Hibler (1979) form

      !-----------------------------------------------------------------
      ! Compute ice strength as in Hibler (1979)
      !-----------------------------------------------------------------

         strength = Pstar*vice*exp(-Cstar*(c1-aice))

      endif                     ! kstrength

      end subroutine colpkg_ice_strength

!=======================================================================

      subroutine colpkg_atm_boundary(sfctype,                    &
                                     Tsf,         potT,          &
                                     uatm,        vatm,          &
                                     wind,        zlvl,          &
                                     Qa,          rhoa,          &
                                     strx,        stry,          &
                                     Tref,        Qref,          &
                                     delt,        delq,          &
                                     lhcoef,      shcoef,        &
                                     Cdn_atm,                    &
                                     Cdn_atm_ratio_n,            &
                                     uvel,        vvel,          &
                                     Uref)

      use ice_atmo, only: atmo_boundary_const, atmo_boundary_layer
      use ice_constants_colpkg, only: c0

      character (len=3), intent(in) :: &
         sfctype      ! ice or ocean

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
         Cdn_atm  , &    ! neutral drag coefficient
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

      real (kind=dbl_kind), optional, intent(in) :: &
         uvel     , & ! x-direction ice speed (m/s)
         vvel         ! y-direction ice speed (m/s)

      real (kind=dbl_kind), optional, intent(out) :: &
         Uref         ! reference height wind speed (m/s)

      real (kind=dbl_kind) :: &
         worku, workv, workr

      worku = c0
      workv = c0
      workr = c0
      if (present(uvel)) then
         worku = uvel
      endif
      if (present(uvel)) then
         worku = uvel
      endif

               if (trim(atmbndy) == 'constant') then
                  call atmo_boundary_const (sfctype,  calc_strair, &
                                            uatm,     vatm,     &
                                            wind,     rhoa,     &
                                            strx,     stry,     &
                                            Tsf,      potT,     &
                                            Qa,                 &
                                            delt,     delq,     &
                                            lhcoef,   shcoef,   &
                                            Cdn_atm)
               else ! default
                  call atmo_boundary_layer (sfctype,                 &
                                            calc_strair, formdrag,   &
                                            highfreq, natmiter,      &
                                            Tsf,      potT,          &
                                            uatm,     vatm,          &
                                            wind,     zlvl,          &
                                            Qa,       rhoa,          &
                                            strx,     stry,          &
                                            Tref,     Qref,          &
                                            delt,     delq,          &
                                            lhcoef,   shcoef,        &
                                            Cdn_atm,                 &
                                            Cdn_atm_ratio_n,         &
                                            worku,    workv,         &
                                            workr)
               endif ! atmbndy

      if (present(Uref)) then
         Uref = workr
      endif

      end subroutine colpkg_atm_boundary

!=======================================================================
! Compute the mixed layer heat balance and update the SST.
! Compute the energy available to freeze or melt ice.
! NOTE: SST changes due to fluxes through the ice are computed in
!       ice_therm_vertical.

      subroutine colpkg_ocn_mixed_layer (alvdr_ocn, swvdr,      &
                                         alidr_ocn, swidr,      &
                                         alvdf_ocn, swvdf,      &
                                         alidf_ocn, swidf,      &
                                         sst,       flwout_ocn, &
                                         fsens_ocn, shcoef,     &
                                         flat_ocn,  lhcoef,     &
                                         evap_ocn,  flw,        &
                                         delt,      delq,       &
                                         aice,      fhocn,      &
                                         fswthru,   hmix,       &
                                         Tf,        qdp,        &
                                         frzmlt,    dt)

      use ice_constants_colpkg, only: c0, c1, c1000, &
          cp_ocn, Tffresh, stefan_boltzmann, Lvap, cprho

      real (kind=dbl_kind), intent(in) :: &
         alvdr_ocn , & ! visible, direct   (fraction)
         alidr_ocn , & ! near-ir, direct   (fraction)
         alvdf_ocn , & ! visible, diffuse  (fraction)
         alidf_ocn , & ! near-ir, diffuse  (fraction)
         swvdr     , & ! sw down, visible, direct  (W/m^2)
         swvdf     , & ! sw down, visible, diffuse (W/m^2)
         swidr     , & ! sw down, near IR, direct  (W/m^2)
         swidf     , & ! sw down, near IR, diffuse (W/m^2)
         flw       , & ! incoming longwave radiation (W/m^2)
         Tf        , & ! freezing temperature (C)
         hmix      , & ! mixed layer depth (m)
         delt      , & ! potential temperature difference   (K)
         delq      , & ! specific humidity difference   (kg/kg)
         shcoef    , & ! transfer coefficient for sensible heat
         lhcoef    , & ! transfer coefficient for latent heat
         fhocn     , & ! net heat flux to ocean (W/m^2)
         fswthru   , & ! shortwave penetrating to ocean (W/m^2)
         aice      , & ! ice area fraction
         dt            ! time step (s)

      real (kind=dbl_kind), intent(inout) :: &
         flwout_ocn, & ! outgoing longwave radiation (W/m^2)
         fsens_ocn , & ! sensible heat flux (W/m^2)
         flat_ocn  , & ! latent heat flux   (W/m^2)
         evap_ocn  , & ! evaporative water flux (kg/m^2/s)
         qdp       , & ! deep ocean heat flux (W/m^2), negative upward
         sst       , & ! sea surface temperature (C)
         frzmlt        ! freezing/melting potential (W/m^2)

      ! local variables

      real (kind=dbl_kind), parameter :: &
         frzmlt_max = c1000   ! max magnitude of frzmlt (W/m^2)

      real (kind=dbl_kind) :: &
         TsfK , & ! surface temperature (K)
         swabs    ! surface absorbed shortwave heat flux (W/m^2)

      ! shortwave radiative flux
      swabs = (c1-alvdr_ocn) * swvdr + (c1-alidr_ocn) * swidr &
            + (c1-alvdf_ocn) * swvdf + (c1-alidf_ocn) * swidf 

      ! ocean surface temperature in Kelvin
      TsfK = sst + Tffresh

      ! longwave radiative flux
      flwout_ocn = -stefan_boltzmann * TsfK**4

      ! downward latent and sensible heat fluxes
      fsens_ocn =  shcoef * delt
      flat_ocn  =  lhcoef * delq
      evap_ocn  = -flat_ocn / Lvap

      ! Compute sst change due to exchange with atm/ice above
      sst = sst + dt * ( &
            (fsens_ocn + flat_ocn + flwout_ocn + flw + swabs) * (c1-aice) &
          + fhocn + fswthru)         &  ! these are *aice
          / (cprho*hmix)

      ! adjust qdp if cooling of mixed layer would occur when sst <= Tf
      if (sst <= Tf .and. qdp > c0) qdp = c0

      ! computed T change due to exchange with deep layers:
      sst = sst - qdp*dt/(cprho*hmix)

      ! compute potential to freeze or melt ice
      frzmlt = (Tf-sst)*cprho*hmix/dt
      frzmlt = min(max(frzmlt,-frzmlt_max),frzmlt_max)

      ! if sst is below freezing, reset sst to Tf
      if (sst <= Tf) sst = Tf

      end subroutine colpkg_ocn_mixed_layer

!=======================================================================
! subroutine to set the column package internal parameters

      subroutine colpkg_init_parameters(&
           ktherm_in, &
           conduct_in, &
           fbot_xfer_type_in, &
           calc_Tsfc_in, &
           ustar_min_in, &
           a_rapid_mode_in, &
           Rac_rapid_mode_in, &
           aspect_rapid_mode_in, &
           dSdt_slow_mode_in, &
           phi_c_slow_mode_in, &
           phi_i_mushy_in, &
           shortwave_in, &
           albedo_type_in, &
           albicev_in, &
           albicei_in, &
           albsnowv_in, &
           albsnowi_in, &
           ahmax_in, &
           R_ice_in, &
           R_pnd_in, &
           R_snw_in, &
           dT_mlt_in, &
           rsnw_mlt_in, &
           kalg_in, &
           kstrength_in, &
           krdg_partic_in, &
           krdg_redist_in, &
           mu_rdg_in, &
           Cf_in, &
           atmbndy_in, &
           calc_strair_in, &
           formdrag_in, &
           highfreq_in, &
           natmiter_in, &
           oceanmixed_ice_in, &
           tfrz_option_in, &
           kitd_in, &
           kcatbound_in, &
           hs0_in, &
           frzpnd_in, &
           dpscale_in, &
           rfracmin_in, &
           rfracmax_in, &
           pndaspect_in, &
           hs1_in, &
           hp1_in, &
         ! bgc_data_dir_in, &
         ! sil_data_type_in, &
         ! nit_data_type_in, &
         ! fe_data_type_in, &
           bgc_flux_type_in, &
           z_tracers_in, &
           scale_bgc_in, &
           solve_zbgc_in, &
           dEdd_algae_in, &
           modal_aero_in, &
           skl_bgc_in, &
           solve_zsal_in, &
           grid_o_in, &
           l_sk_in, &
           grid_o_t_in, &
           initbio_frac_in, &
           frazil_scav_in, &
           grid_oS_in, &
           l_skS_in, &
           phi_snow_in, &
           ratio_Si2N_diatoms_in, &
           ratio_Si2N_sp_in, &
           ratio_Si2N_phaeo_in, &
           ratio_S2N_diatoms_in, &
           ratio_S2N_sp_in, &      
           ratio_S2N_phaeo_in, &   
           ratio_Fe2C_diatoms_in, & 
           ratio_Fe2C_sp_in, &     
           ratio_Fe2C_phaeo_in, &  
           ratio_Fe2N_diatoms_in, & 
           ratio_Fe2N_sp_in, &     
           ratio_Fe2N_phaeo_in, &  
           ratio_Fe2DON_in, &       
           ratio_Fe2DOC_s_in, &     
           ratio_Fe2DOC_l_in, &     
           fr_resp_in, &            
           tau_min_in, &            
           tau_max_in, &            
           algal_vel_in, &          
           R_dFe2dust_in, &         
           dustFe_sol_in, &         
           chlabs_diatoms_in, &    
           chlabs_sp_in, &         
           chlabs_phaeo_in, &      
           alpha2max_low_diatoms_in, &  
           alpha2max_low_sp_in, &       
           alpha2max_low_phaeo_in, &    
           beta2max_diatoms_in, & 
           beta2max_sp_in, &       
           beta2max_phaeo_in, &    
           mu_max_diatoms_in, &   
           mu_max_sp_in, &         
           mu_max_phaeo_in, &      
           grow_Tdep_diatoms_in, &
           grow_Tdep_sp_in, &      
           grow_Tdep_phaeo_in, &   
           fr_graze_diatoms_in, & 
           fr_graze_sp_in, &       
           fr_graze_phaeo_in, &    
           mort_pre_diatoms_in, & 
           mort_pre_sp_in, &       
           mort_pre_phaeo_in, &    
           mort_Tdep_diatoms_in, &
           mort_Tdep_sp_in, &       
           mort_Tdep_phaeo_in, &    
           k_exude_diatoms_in, &  
           k_exude_sp_in, &         
           k_exude_phaeo_in, &      
           K_Nit_diatoms_in, &    
           K_Nit_sp_in, &           
           K_Nit_phaeo_in, &        
           K_Am_diatoms_in, &     
           K_Am_sp_in, &             
           K_Am_phaeo_in, &          
           K_Sil_diatoms_in, &    
           K_Sil_sp_in, &            
           K_Sil_phaeo_in, &         
           K_Fe_diatoms_in, &     
           K_Fe_sp_in, &             
           K_Fe_phaeo_in, &           
           f_don_protein_in, &    
           kn_bac_protein_in, &   
           f_don_Am_protein_in, & 
           f_doc_s_in, &            
           f_doc_l_in, &               
           f_exude_s_in, &          
           f_exude_l_in, &           
           k_bac_s_in, &            
           k_bac_l_in, &             
           T_max_in, &              
           fsal_in, &               
           op_dep_min_in, &         
           fr_graze_s_in, &         
           fr_graze_e_in, &         
           fr_mort2min_in, &        
           fr_dFe_in, &             
           k_nitrif_in, &           
           t_iron_conv_in, &        
           max_loss_in, &           
           max_dfe_doc1_in, &       
           fr_resp_s_in, &          
           y_sk_DMS_in, &           
           t_sk_conv_in, &          
           t_sk_ox_in, &             
           algaltype_diatoms_in, &   
           algaltype_sp_in, &       
           algaltype_phaeo_in, &    
           nitratetype_in, &        
           ammoniumtype_in, &       
           silicatetype_in, &       
           dmspptype_in, &          
           dmspdtype_in, &          
           humtype_in, &            
           doctype_s_in, &          
           doctype_l_in, &          
           dontype_protein_in, &     
           fedtype_1_in, &           
           feptype_1_in, &           
           zaerotype_bc1_in, &       
           zaerotype_bc2_in, &       
           zaerotype_dust1_in, &     
           zaerotype_dust2_in, &     
           zaerotype_dust3_in, &     
           zaerotype_dust4_in, &     
           ratio_C2N_diatoms_in, &   
           ratio_C2N_sp_in, &        
           ratio_C2N_phaeo_in, &     
           ratio_chl2N_diatoms_in, & 
           ratio_chl2N_sp_in, &      
           ratio_chl2N_phaeo_in, &   
           F_abs_chl_diatoms_in, &   
           F_abs_chl_sp_in, &        
           F_abs_chl_phaeo_in, &
           ratio_C2N_proteins_in)
           !restore_bgc_in)

        use ice_colpkg_shared, only: &
             ktherm, &
             conduct, &
             fbot_xfer_type, &
             calc_Tsfc, &
             ustar_min, &
             a_rapid_mode, &
             Rac_rapid_mode, &
             aspect_rapid_mode, &
             dSdt_slow_mode, &
             phi_c_slow_mode, &
             phi_i_mushy, &
             shortwave, &
             albedo_type, &
             albicev, &
             albicei, &
             albsnowv, &
             albsnowi, &
             ahmax, &
             R_ice, &
             R_pnd, &
             R_snw, &
             dT_mlt, &
             rsnw_mlt, &
             kalg, &
             kstrength, &
             krdg_partic, &
             krdg_redist, &
             mu_rdg, &
             Cf, &
             atmbndy, &
             calc_strair, &
             formdrag, &
             highfreq, &
             natmiter, &
             oceanmixed_ice, &
             tfrz_option, &
             kitd, &
             kcatbound, &
             hs0, &
             frzpnd, &
             dpscale, &
             rfracmin, &
             rfracmax, &
             pndaspect, &
             hs1, &
             hp1, &
           ! bgc_data_dir, &
           ! sil_data_type, &
           ! nit_data_type, &
           ! fe_data_type, &
             bgc_flux_type, &
             z_tracers, &
             scale_bgc, &
             solve_zbgc, &
             dEdd_algae, &
             modal_aero, &
             skl_bgc, &
             solve_zsal, &
             grid_o, &
             l_sk, &
             grid_o_t, &
             initbio_frac, &
             frazil_scav, &
             grid_oS, &
             l_skS, &
             phi_snow, &
             ratio_Si2N_diatoms, & 
             ratio_Si2N_sp     , &
             ratio_Si2N_phaeo  , &
             ratio_S2N_diatoms , & 
             ratio_S2N_sp      , &
             ratio_S2N_phaeo   , &
             ratio_Fe2C_diatoms, & 
             ratio_Fe2C_sp     , &
             ratio_Fe2C_phaeo  , &
             ratio_Fe2N_diatoms, & 
             ratio_Fe2N_sp     , &
             ratio_Fe2N_phaeo  , &
             ratio_Fe2DON      , & 
             ratio_Fe2DOC_s    , & 
             ratio_Fe2DOC_l    , & 
             fr_resp           , & 
             tau_min           , & 
             tau_max           , & 
             algal_vel         , & 
             R_dFe2dust        , & 
             dustFe_sol        , & 
             chlabs_diatoms    , &
             chlabs_sp         , &
             chlabs_phaeo      , &
             alpha2max_low_diatoms , & 
             alpha2max_low_sp      , & 
             alpha2max_low_phaeo   , & 
             beta2max_diatoms , &
             beta2max_sp      , & 
             beta2max_phaeo   , & 
             mu_max_diatoms   , &
             mu_max_sp        , & 
             mu_max_phaeo     , & 
             grow_Tdep_diatoms, &
             grow_Tdep_sp     , & 
             grow_Tdep_phaeo  , & 
             fr_graze_diatoms , &
             fr_graze_sp      , & 
             fr_graze_phaeo   , & 
             mort_pre_diatoms , &
             mort_pre_sp      , & 
             mort_pre_phaeo   , & 
             mort_Tdep_diatoms, &
             mort_Tdep_sp     , &  
             mort_Tdep_phaeo  , &  
             k_exude_diatoms  , &
             k_exude_sp       , &  
             k_exude_phaeo    , &  
             K_Nit_diatoms    , &
             K_Nit_sp         , &  
             K_Nit_phaeo      , &  
             K_Am_diatoms     , &
             K_Am_sp          , &   
             K_Am_phaeo       , &   
             K_Sil_diatoms    , &
             K_Sil_sp         , &   
             K_Sil_phaeo      , &   
             K_Fe_diatoms     , &
             K_Fe_sp          , &   
             K_Fe_phaeo       , &    
             f_don_protein    , &
             kn_bac_protein   , &
             f_don_Am_protein , &
             f_doc_s            , &
             f_doc_l            , &   
             f_exude_s          , &
             f_exude_l          , & 
             k_bac_s            , &
             k_bac_l            , & 
             T_max              , &
             fsal               , &
             op_dep_min         , &
             fr_graze_s         , &
             fr_graze_e         , &
             fr_mort2min        , &
             fr_dFe             , &
             k_nitrif           , &
             t_iron_conv        , &
             max_loss           , &
             max_dfe_doc1       , &
             fr_resp_s          , &
             y_sk_DMS           , &
             t_sk_conv          , &
             t_sk_ox            , & 
             algaltype_diatoms  , & 
             algaltype_sp       , &
             algaltype_phaeo    , &
             nitratetype        , &
             ammoniumtype       , &
             silicatetype       , &
             dmspptype          , &
             dmspdtype          , &
             humtype            , &
             doctype_s          , &
             doctype_l          , &
             dontype_protein    , & 
             fedtype_1          , & 
             feptype_1          , & 
             zaerotype_bc1      , & 
             zaerotype_bc2      , & 
             zaerotype_dust1    , & 
             zaerotype_dust2    , & 
             zaerotype_dust3    , & 
             zaerotype_dust4    , & 
             ratio_C2N_diatoms  , & 
             ratio_C2N_sp       , & 
             ratio_C2N_phaeo    , & 
             ratio_chl2N_diatoms, & 
             ratio_chl2N_sp     , & 
             ratio_chl2N_phaeo  , & 
             F_abs_chl_diatoms  , & 
             F_abs_chl_sp       , & 
             F_abs_chl_phaeo    , & 
             ratio_C2N_proteins
            !restore_bgc

!-----------------------------------------------------------------------
! Parameters for thermodynamics
!-----------------------------------------------------------------------

        integer (kind=int_kind), intent(in) :: &
             ktherm_in          ! type of thermodynamics
                                ! 0 = 0-layer approximation
                                ! 1 = Bitz and Lipscomb 1999
                                ! 2 = mushy layer theory

        character (char_len), intent(in) :: &
             conduct_in, &      ! 'MU71' or 'bubbly'
             fbot_xfer_type_in  ! transfer coefficient type for ice-ocean heat flux
        
        logical (kind=log_kind), intent(in) :: &
             calc_Tsfc_in       ! if true, calculate surface temperature
                                ! if false, Tsfc is computed elsewhere and
                                ! atmos-ice fluxes are provided to CICE

        real (kind=dbl_kind), intent(in) :: &
             ustar_min_in       ! minimum friction velocity for ice-ocean heat flux
 
        ! mushy thermo
        real(kind=dbl_kind), intent(in) :: &
             a_rapid_mode_in      , & ! channel radius for rapid drainage mode (m)
             Rac_rapid_mode_in    , & ! critical Rayleigh number for rapid drainage mode
             aspect_rapid_mode_in , & ! aspect ratio for rapid drainage mode (larger=wider)
             dSdt_slow_mode_in    , & ! slow mode drainage strength (m s-1 K-1)
             phi_c_slow_mode_in   , & ! liquid fraction porosity cutoff for slow mode
             phi_i_mushy_in           ! liquid fraction of congelation ice
        
!-----------------------------------------------------------------------
! Parameters for radiation
!-----------------------------------------------------------------------

        character (len=char_len), intent(in) :: &
             shortwave_in, & ! shortwave method, 'default' ('ccsm3') or 'dEdd'
             albedo_type_in  ! albedo parameterization, 'default' ('ccsm3') or 'constant'
                             ! shortwave='dEdd' overrides this parameter

        ! baseline albedos for ccsm3 shortwave, set in namelist
        real (kind=dbl_kind), intent(in) :: &
             albicev_in  , & ! visible ice albedo for h > ahmax
             albicei_in  , & ! near-ir ice albedo for h > ahmax
             albsnowv_in , & ! cold snow albedo, visible
             albsnowi_in , & ! cold snow albedo, near IR
             ahmax_in        ! thickness above which ice albedo is constant (m)
        
        ! dEdd tuning parameters, set in namelist
        real (kind=dbl_kind), intent(in) :: &
             R_ice_in    , & ! sea ice tuning parameter; +1 > 1sig increase in albedo
             R_pnd_in    , & ! ponded ice tuning parameter; +1 > 1sig increase in albedo
             R_snw_in    , & ! snow tuning parameter; +1 > ~.01 change in broadband albedo
             dT_mlt_in   , & ! change in temp for non-melt to melt snow grain 
                             ! radius change (C)
             rsnw_mlt_in , & ! maximum melting snow grain radius (10^-6 m)
             kalg_in         ! algae absorption coefficient for 0.5 m thick layer

!-----------------------------------------------------------------------
! Parameters for ridging and strength
!-----------------------------------------------------------------------

        integer (kind=int_kind), intent(in) :: & ! defined in namelist 
             kstrength_in  , & ! 0 for simple Hibler (1979) formulation 
                               ! 1 for Rothrock (1975) pressure formulation 
             krdg_partic_in, & ! 0 for Thorndike et al. (1975) formulation 
                               ! 1 for exponential participation function 
             krdg_redist_in    ! 0 for Hibler (1980) formulation 
                               ! 1 for exponential redistribution function 
 
        real (kind=dbl_kind), intent(in) :: &  
             mu_rdg_in, &      ! gives e-folding scale of ridged ice (m^.5) 
                               ! (krdg_redist = 1) 
             Cf_in             ! ratio of ridging work to PE change in ridging (kstrength = 1)

!-----------------------------------------------------------------------
! Parameters for atmosphere
!-----------------------------------------------------------------------

        character (len=char_len), intent(in) :: &
             atmbndy_in ! atmo boundary method, 'default' ('ccsm3') or 'constant'
        
        logical (kind=log_kind), intent(in) :: &
             calc_strair_in, &  ! if true, calculate wind stress components
             formdrag_in,    &  ! if true, calculate form drag
             highfreq_in        ! if true, use high frequency coupling
        
        integer (kind=int_kind), intent(in) :: &
             natmiter_in        ! number of iterations for boundary layer calculations
        
!-----------------------------------------------------------------------
! Parameters for ocean
!-----------------------------------------------------------------------

        logical (kind=log_kind), intent(in) :: &
             oceanmixed_ice_in           ! if true, use ocean mixed layer
        
        character(len=char_len), intent(in) :: &
             tfrz_option_in              ! form of ocean freezing temperature
                                         ! 'minus1p8' = -1.8 C
                                         ! 'linear_salt' = -depressT * sss
                                         ! 'mushy' conforms with ktherm=2

!-----------------------------------------------------------------------
! Parameters for the ice thickness distribution
!-----------------------------------------------------------------------

        integer (kind=int_kind), intent(in) :: &
             kitd_in        , & ! type of itd conversions
                                !   0 = delta function
                                !   1 = linear remap
             kcatbound_in       !   0 = old category boundary formula
                                !   1 = new formula giving round numbers
                                !   2 = WMO standard
                                !   3 = asymptotic formula

!-----------------------------------------------------------------------
! Parameters for biogeochemistry
!-----------------------------------------------------------------------

   !  character(char_len_long), intent(in) :: & 
   !     bgc_data_dir_in   ! directory for biogeochemistry data

     character(char_len), intent(in) :: &     
        bgc_flux_type_in    ! type of ocean-ice piston velocity 
                            ! 'constant', 'Jin2006'      
    !    sil_data_type_in  , & ! 'default', 'clim'
    !    nit_data_type_in  , & ! 'default', 'clim'   
    !    fe_data_type_in   , & ! 'default', 'clim'      

      logical (kind=log_kind), intent(in) :: &
         z_tracers_in,      & ! if .true., bgc or aerosol tracers are vertically resolved
         scale_bgc_in,      & ! if .true., initialize bgc tracers proportionally with salinity
         solve_zbgc_in,     & ! if .true., solve vertical biochemistry portion of code
         dEdd_algae_in,     & ! if .true., algal absorptionof Shortwave is computed in the
         modal_aero_in        ! if .true., use modal aerosol formulation in shortwave
        
      logical (kind=log_kind), intent(in) :: & 
         skl_bgc_in,        &   ! if true, solve skeletal biochemistry
         solve_zsal_in          ! if true, update salinity profile from solve_S_dt

      real (kind=dbl_kind), intent(in) :: & 
         grid_o_in      , & ! for bottom flux        
         l_sk_in        , & ! characteristic diffusive scale (zsalinity) (m)
         grid_o_t_in    , & ! top grid point length scale 
         initbio_frac_in, & ! fraction of ocean tracer concentration used to initialize tracer 
         frazil_scav_in , & ! multiple of ocean tracer concentration due to frazil scavenging
         phi_snow_in        ! snow porosity at the ice/snow interface 

      real (kind=dbl_kind), intent(in) :: & 
         grid_oS_in     , & ! for bottom flux (zsalinity)
         l_skS_in           ! 0.02 characteristic skeletal layer thickness (m) (zsalinity)
      real (kind=dbl_kind), intent(in) :: &
         ratio_Si2N_diatoms_in, &   ! algal Si to N (mol/mol)
         ratio_Si2N_sp_in     , &
         ratio_Si2N_phaeo_in  , &
         ratio_S2N_diatoms_in , &   ! algal S  to N (mol/mol)
         ratio_S2N_sp_in      , &
         ratio_S2N_phaeo_in   , &
         ratio_Fe2C_diatoms_in, &   ! algal Fe to C  (umol/mol)
         ratio_Fe2C_sp_in     , &
         ratio_Fe2C_phaeo_in  , &
         ratio_Fe2N_diatoms_in, &   ! algal Fe to N  (umol/mol)
         ratio_Fe2N_sp_in     , &
         ratio_Fe2N_phaeo_in  , &
         ratio_Fe2DON_in      , &   ! Fe to N of DON (nmol/umol)
         ratio_Fe2DOC_s_in    , &   ! Fe to C of DOC (nmol/umol) saccharids
         ratio_Fe2DOC_l_in    , &   ! Fe to C of DOC (nmol/umol) lipids
         fr_resp_in           , &   ! fraction of algal growth lost due to respiration
         tau_min_in           , &   ! rapid mobile to stationary exchanges (s) = 1.5 hours
         tau_max_in           , &   ! long time mobile to stationary exchanges (s) = 2 days
         algal_vel_in         , &   ! 0.5 cm/d(m/s) Lavoie 2005  1.5 cm/day
         R_dFe2dust_in        , &   !  g/g (3.5% content) Tagliabue 2009
         dustFe_sol_in        , &   ! solubility fraction
         chlabs_diatoms_in   , & ! chl absorption (1/m/(mg/m^3))
         chlabs_sp_in        , & !
         chlabs_phaeo_in     , & !
         alpha2max_low_diatoms_in , & ! light limitation (1/(W/m^2))  
         alpha2max_low_sp_in      , & 
         alpha2max_low_phaeo_in   , & 
         beta2max_diatoms_in , & ! light inhibition (1/(W/m^2))  
         beta2max_sp_in      , & 
         beta2max_phaeo_in   , & 
         mu_max_diatoms_in   , & ! maximum growth rate (1/day)       
         mu_max_sp_in        , & 
         mu_max_phaeo_in     , & 
         grow_Tdep_diatoms_in, & ! Temperature dependence of growth (1/C)
         grow_Tdep_sp_in     , & 
         grow_Tdep_phaeo_in  , & 
         fr_graze_diatoms_in , & ! Fraction grazed
         fr_graze_sp_in      , & 
         fr_graze_phaeo_in   , & 
         mort_pre_diatoms_in , & ! Mortality (1/day)
         mort_pre_sp_in      , & 
         mort_pre_phaeo_in   , & 
         mort_Tdep_diatoms_in, & ! T dependence of mortality (1/C) 
         mort_Tdep_sp_in     , &  
         mort_Tdep_phaeo_in  , &  
         k_exude_diatoms_in  , & ! algal exudation (1/d)
         k_exude_sp_in       , &  
         k_exude_phaeo_in    , &  
         K_Nit_diatoms_in    , & ! nitrate half saturation (mmol/m^3)
         K_Nit_sp_in         , &  
         K_Nit_phaeo_in      , &  
         K_Am_diatoms_in     , & ! ammonium half saturation (mmol/m^3)
         K_Am_sp_in          , &   
         K_Am_phaeo_in       , &   
         K_Sil_diatoms_in    , & ! silicate half saturation (mmol/m^3)
         K_Sil_sp_in         , &   
         K_Sil_phaeo_in      , &   
         K_Fe_diatoms_in     , & ! iron half saturation (nM)
         K_Fe_sp_in          , &   
         K_Fe_phaeo_in       , &    
         f_don_protein_in    , & ! fraction of spilled grazing to proteins            
         kn_bac_protein_in   , & ! Bacterial degredation of DON (1/d)                  
         f_don_Am_protein_in , & ! fraction of remineralized DON to ammonium          
         f_doc_s_in          , & ! fraction of mortality to DOC 
         f_doc_l_in          , &   
         f_exude_s_in        , & ! fraction of exudation to DOC
         f_exude_l_in        , & 
         k_bac_s_in          , & ! Bacterial degredation of DOC (1/d)
         k_bac_l_in          , & 
         T_max_in            , & ! maximum temperature (C)
         fsal_in             , & ! Salinity limitation (ppt)
         op_dep_min_in       , & ! Light attenuates for optical depths exceeding min
         fr_graze_s_in       , & ! fraction of grazing spilled or slopped
         fr_graze_e_in       , & ! fraction of assimilation excreted 
         fr_mort2min_in      , & ! fractionation of mortality to Am
         fr_dFe_in           , & ! fraction of remineralized nitrogen 
                                    ! (in units of algal iron)
         k_nitrif_in         , & ! nitrification rate (1/day)            
         t_iron_conv_in      , & ! desorption loss pFe to dFe (day)
         max_loss_in         , & ! restrict uptake to % of remaining value 
         max_dfe_doc1_in     , & ! max ratio of dFe to saccharides in the ice 
                                    ! (nM Fe/muM C)    
         fr_resp_s_in        , & ! DMSPd fraction of respiration loss as DMSPd
         y_sk_DMS_in         , & ! fraction conversion given high yield
         t_sk_conv_in        , & ! Stefels conversion time (d)
         t_sk_ox_in          , &   ! DMS oxidation time (d)
         algaltype_diatoms_in  , & ! mobility type
         algaltype_sp_in       , & !
         algaltype_phaeo_in    , & !
         nitratetype_in        , & !
         ammoniumtype_in       , & !
         silicatetype_in       , & !
         dmspptype_in          , & !
         dmspdtype_in          , & !
         humtype_in            , & !
         doctype_s_in          , & !
         doctype_l_in          , & !
         dontype_protein_in    , & !
         fedtype_1_in          , & !
         feptype_1_in          , & !
         zaerotype_bc1_in      , & !
         zaerotype_bc2_in      , & !
         zaerotype_dust1_in    , & !
         zaerotype_dust2_in    , & !
         zaerotype_dust3_in    , & !
         zaerotype_dust4_in    , & !
         ratio_C2N_diatoms_in  , & ! algal C to N ratio (mol/mol)
         ratio_C2N_sp_in       , & !
         ratio_C2N_phaeo_in    , & !
         ratio_chl2N_diatoms_in, & ! algal chlorophyll to N ratio (mg/mmol)
         ratio_chl2N_sp_in     , & !
         ratio_chl2N_phaeo_in  , & !
         F_abs_chl_diatoms_in  , & ! scales absorbed radiation for dEdd
         F_abs_chl_sp_in       , & !
         F_abs_chl_phaeo_in    , & !
         ratio_C2N_proteins_in     ! ratio of C to N in proteins (mol/mol)       

     !logical (kind=log_kind), intent(in) :: & 
     !   restore_bgc_in      ! if true, restore nitrate

!-----------------------------------------------------------------------
! Parameters for melt ponds
!-----------------------------------------------------------------------

        real (kind=dbl_kind), intent(in) :: &
             hs0_in             ! snow depth for transition to bare sea ice (m)
        
        ! level-ice ponds
        character (len=char_len), intent(in) :: &
             frzpnd_in          ! pond refreezing parameterization
        
        real (kind=dbl_kind), intent(in) :: &
             dpscale_in, &      ! alter e-folding time scale for flushing 
             rfracmin_in, &     ! minimum retained fraction of meltwater
             rfracmax_in, &     ! maximum retained fraction of meltwater
             pndaspect_in, &    ! ratio of pond depth to pond fraction
             hs1_in             ! tapering parameter for snow on pond ice
        
        ! topo ponds
        real (kind=dbl_kind), intent(in) :: &
             hp1_in             ! critical parameter for pond ice thickness
        
        ktherm = ktherm_in
        conduct = conduct_in
        fbot_xfer_type = fbot_xfer_type_in
        calc_Tsfc = calc_Tsfc_in
        ustar_min = ustar_min_in
        a_rapid_mode = a_rapid_mode_in
        Rac_rapid_mode = Rac_rapid_mode_in
        aspect_rapid_mode = aspect_rapid_mode_in
        dSdt_slow_mode = dSdt_slow_mode_in
        phi_c_slow_mode = phi_c_slow_mode_in
        phi_i_mushy = phi_i_mushy_in
        shortwave = shortwave_in
        albedo_type = albedo_type_in
        albicev = albicev_in
        albicei = albicei_in
        albsnowv = albsnowv_in
        albsnowi = albsnowi_in
        ahmax = ahmax_in
        R_ice = R_ice_in
        R_pnd = R_pnd_in
        R_snw = R_snw_in
        dT_mlt = dT_mlt_in
        rsnw_mlt = rsnw_mlt_in
        kalg = kalg_in
        kstrength = kstrength_in
        krdg_partic = krdg_partic_in
        krdg_redist = krdg_redist_in
        mu_rdg = mu_rdg_in
        Cf = Cf_in
        atmbndy = atmbndy_in
        calc_strair = calc_strair_in
        formdrag = formdrag_in
        highfreq = highfreq_in
        natmiter = natmiter_in
        oceanmixed_ice = oceanmixed_ice_in
        tfrz_option = tfrz_option_in
        kitd = kitd_in
        kcatbound = kcatbound_in
        hs0 = hs0_in
        frzpnd = frzpnd_in
        dpscale = dpscale_in
        rfracmin = rfracmin_in
        rfracmax = rfracmax_in
        pndaspect = pndaspect_in
        hs1 = hs1_in
        hp1 = hp1_in
     !  bgc_data_dir = bgc_data_dir_in
     !  sil_data_type= sil_data_type_in
     !  nit_data_type = nit_data_type_in
     !  fe_data_type = fe_data_type_in
        bgc_flux_type = bgc_flux_type_in
        z_tracers = z_tracers_in
        scale_bgc = scale_bgc_in
        solve_zbgc = solve_zbgc_in
        dEdd_algae = dEdd_algae_in
        skl_bgc = skl_bgc_in
        grid_o = grid_o_in
        l_sk = l_sk_in
        grid_o_t = grid_o_t_in
        initbio_frac = initbio_frac_in
        frazil_scav = frazil_scav_in
        grid_oS = grid_oS_in
        l_skS = l_skS_in
        phi_snow = phi_snow_in
     !  restore_bgc = restore_bgc_in
        ratio_Si2N_diatoms= ratio_Si2N_diatoms_in 
        ratio_Si2N_sp     = ratio_Si2N_sp_in
        ratio_Si2N_phaeo  = ratio_Si2N_phaeo_in
        ratio_S2N_diatoms = ratio_S2N_diatoms_in
        ratio_S2N_sp      = ratio_S2N_sp_in
        ratio_S2N_phaeo   = ratio_S2N_phaeo_in
        ratio_Fe2C_diatoms= ratio_Fe2C_diatoms_in 
        ratio_Fe2C_sp     = ratio_Fe2C_sp_in
        ratio_Fe2C_phaeo  = ratio_Fe2C_phaeo_in
        ratio_Fe2N_diatoms= ratio_Fe2N_diatoms_in 
        ratio_Fe2N_sp     = ratio_Fe2N_sp_in
        ratio_Fe2N_phaeo  = ratio_Fe2N_phaeo_in
        ratio_Fe2DON      = ratio_Fe2DON_in
        ratio_Fe2DOC_s    = ratio_Fe2DOC_s_in
        ratio_Fe2DOC_l    = ratio_Fe2DOC_l_in
        fr_resp           = fr_resp_in
        tau_min           = tau_min_in
        tau_max           = tau_max_in
        algal_vel         = algal_vel_in
        R_dFe2dust        = R_dFe2dust_in
        dustFe_sol        = dustFe_sol_in
        chlabs_diatoms    = chlabs_diatoms_in
        chlabs_sp         = chlabs_sp_in
        chlabs_phaeo      = chlabs_phaeo_in
        alpha2max_low_diatoms = alpha2max_low_diatoms_in
        alpha2max_low_sp      = alpha2max_low_sp_in
        alpha2max_low_phaeo   = alpha2max_low_phaeo_in
        beta2max_diatoms = beta2max_diatoms_in
        beta2max_sp      = beta2max_sp_in
        beta2max_phaeo   = beta2max_phaeo_in
        mu_max_diatoms   = mu_max_diatoms_in
        mu_max_sp        = mu_max_sp_in
        mu_max_phaeo     = mu_max_phaeo_in
        grow_Tdep_diatoms= grow_Tdep_diatoms_in
        grow_Tdep_sp     = grow_Tdep_sp_in
        grow_Tdep_phaeo  = grow_Tdep_phaeo_in
        fr_graze_diatoms = fr_graze_diatoms_in
        fr_graze_sp      = fr_graze_sp_in
        fr_graze_phaeo   = fr_graze_phaeo_in
        mort_pre_diatoms = mort_pre_diatoms_in
        mort_pre_sp      = mort_pre_sp_in
        mort_pre_phaeo   = mort_pre_phaeo_in
        mort_Tdep_diatoms= mort_Tdep_diatoms_in
        mort_Tdep_sp     = mort_Tdep_sp_in
        mort_Tdep_phaeo  = mort_Tdep_phaeo_in
        k_exude_diatoms  = k_exude_diatoms_in
        k_exude_sp       = k_exude_sp_in
        k_exude_phaeo    = k_exude_phaeo_in
        K_Nit_diatoms    = K_Nit_diatoms_in
        K_Nit_sp         = K_Nit_sp_in
        K_Nit_phaeo      = K_Nit_phaeo_in
        K_Am_diatoms     = K_Am_diatoms_in
        K_Am_sp          = K_Am_sp_in
        K_Am_phaeo       = K_Am_phaeo_in
        K_Sil_diatoms    = K_Sil_diatoms_in
        K_Sil_sp         = K_Sil_sp_in
        K_Sil_phaeo      = K_Sil_phaeo_in
        K_Fe_diatoms     = K_Fe_diatoms_in
        K_Fe_sp          = K_Fe_sp_in
        K_Fe_phaeo       = K_Fe_phaeo_in
        f_don_protein    = f_don_protein_in
        kn_bac_protein   = kn_bac_protein_in
        f_don_Am_protein = f_don_Am_protein_in
        f_doc_s          = f_doc_s_in
        f_doc_l          = f_doc_l_in
        f_exude_s        = f_exude_s_in
        f_exude_l        = f_exude_l_in
        k_bac_s          = k_bac_s_in
        k_bac_l          = k_bac_l_in
        T_max            = T_max_in
        fsal             = fsal_in
        op_dep_min       = op_dep_min_in
        fr_graze_s       = fr_graze_s_in
        fr_graze_e       = fr_graze_e_in
        fr_mort2min      = fr_mort2min_in
        fr_dFe           = fr_dFe_in
        k_nitrif         = k_nitrif_in
        t_iron_conv      = t_iron_conv_in
        max_loss         = max_loss_in
        max_dfe_doc1     = max_dfe_doc1_in
        fr_resp_s        = fr_resp_s_in
        y_sk_DMS         = y_sk_DMS_in
        t_sk_conv        = t_sk_conv_in
        t_sk_ox          = t_sk_ox_in
        algaltype_diatoms  = algaltype_diatoms_in
        algaltype_sp       = algaltype_sp_in
        algaltype_phaeo    = algaltype_phaeo_in
        nitratetype        = nitratetype_in
        ammoniumtype       = ammoniumtype_in
        silicatetype       = silicatetype_in
        dmspptype          = dmspptype_in
        dmspdtype          = dmspdtype_in
        humtype            = humtype_in
        doctype_s          = doctype_s_in
        doctype_l          = doctype_l_in
        dontype_protein    = dontype_protein_in
        fedtype_1          = fedtype_1_in
        feptype_1          = feptype_1_in
        zaerotype_bc1      = zaerotype_bc1_in
        zaerotype_bc2      = zaerotype_bc2_in
        zaerotype_dust1    = zaerotype_dust1_in
        zaerotype_dust2    = zaerotype_dust2_in
        zaerotype_dust3    = zaerotype_dust3_in
        zaerotype_dust4    = zaerotype_dust4_in
        ratio_C2N_diatoms  = ratio_C2N_diatoms_in
        ratio_C2N_sp       = ratio_C2N_sp_in
        ratio_C2N_phaeo    = ratio_C2N_phaeo_in
        ratio_chl2N_diatoms= ratio_chl2N_diatoms_in
        ratio_chl2N_sp     = ratio_chl2N_sp_in
        ratio_chl2N_phaeo  = ratio_chl2N_phaeo_in
        F_abs_chl_diatoms  = F_abs_chl_diatoms_in
        F_abs_chl_sp       = F_abs_chl_sp_in
        F_abs_chl_phaeo    = F_abs_chl_phaeo_in
        ratio_C2N_proteins = ratio_C2N_proteins_in

      end subroutine colpkg_init_parameters

!=======================================================================
! set tracer active flags

      subroutine colpkg_init_tracer_flags(&
           tr_iage_in      , & ! if .true., use age tracer
           tr_FY_in        , & ! if .true., use first-year area tracer
           tr_lvl_in       , & ! if .true., use level ice tracer
           tr_pond_in      , & ! if .true., use melt pond tracer
           tr_pond_cesm_in , & ! if .true., use cesm pond tracer
           tr_pond_lvl_in  , & ! if .true., use level-ice pond tracer
           tr_pond_topo_in , & ! if .true., use explicit topography-based ponds
           tr_aero_in      , & ! if .true., use aerosol tracers
           tr_brine_in     , & ! if .true., brine height differs from ice thickness
           tr_bgc_S_in     , & ! if .true., use zsalinity
           tr_zaero_in     , & ! if .true., black carbon is tracers  (n_zaero)
           tr_bgc_Nit_in   , & ! if .true., Nitrate tracer in ice 
           tr_bgc_N_in     , & ! if .true., algal nitrogen tracers  (n_algae)
           tr_bgc_DON_in   , & ! if .true., DON pools are tracers  (n_don)
           tr_bgc_C_in     , & ! if .true., algal carbon tracers + DOC and DIC 
           tr_bgc_chl_in   , & ! if .true., algal chlorophyll tracers 
           tr_bgc_Am_in    , & ! if .true., ammonia/um as nutrient tracer 
           tr_bgc_Sil_in   , & ! if .true., silicon as nutrient tracer 
           tr_bgc_DMS_in   , & ! if .true., DMS as product tracer 
           tr_bgc_Fe_in    , & ! if .true., Fe as product tracer 
           tr_bgc_hum_in   , & ! if .true., hum as tracer 
           tr_bgc_PON_in)      ! if .true., PON as product tracer 


        use ice_colpkg_tracers, only: &
             tr_iage      , & ! if .true., use age tracer
             tr_FY        , & ! if .true., use first-year area tracer
             tr_lvl       , & ! if .true., use level ice tracer
             tr_pond      , & ! if .true., use melt pond tracer
             tr_pond_cesm , & ! if .true., use cesm pond tracer
             tr_pond_lvl  , & ! if .true., use level-ice pond tracer
             tr_pond_topo , & ! if .true., use explicit topography-based ponds
             tr_aero      , & ! if .true., use aerosol tracers
             tr_brine     , & ! if .true., brine height differs from ice thickness
             tr_bgc_S     , & ! if .true., use zsalinity
             tr_zaero     , & ! if .true., black carbon is tracers  (n_zaero)
             tr_bgc_Nit   , & ! if .true., Nitrate tracer in ice 
             tr_bgc_N     , & ! if .true., algal nitrogen tracers  (n_algae)
             tr_bgc_DON   , & ! if .true., DON pools are tracers  (n_don)
             tr_bgc_C     , & ! if .true., algal carbon tracers + DOC and DIC 
             tr_bgc_chl   , & ! if .true., algal chlorophyll tracers 
             tr_bgc_Am    , & ! if .true., ammonia/um as nutrient tracer 
             tr_bgc_Sil   , & ! if .true., silicon as nutrient tracer 
             tr_bgc_DMS   , & ! if .true., DMS as product tracer 
             tr_bgc_Fe    , & ! if .true., Fe as product tracer 
             tr_bgc_hum   , & ! if .true., hum as product tracer 
             tr_bgc_PON       ! if .true., PON as product tracer 


        logical, intent(in) :: &
             tr_iage_in      , & ! if .true., use age tracer
             tr_FY_in        , & ! if .true., use first-year area tracer
             tr_lvl_in       , & ! if .true., use level ice tracer
             tr_pond_in      , & ! if .true., use melt pond tracer
             tr_pond_cesm_in , & ! if .true., use cesm pond tracer
             tr_pond_lvl_in  , & ! if .true., use level-ice pond tracer
             tr_pond_topo_in , & ! if .true., use explicit topography-based ponds
             tr_aero_in      , & ! if .true., use aerosol tracers
             tr_brine_in     , & ! if .true., brine height differs from ice thickness
             tr_bgc_S_in     , & ! if .true., use zsalinity
             tr_zaero_in     , & ! if .true., black carbon is tracers  (n_zaero)
             tr_bgc_Nit_in   , & ! if .true., Nitrate tracer in ice 
             tr_bgc_N_in     , & ! if .true., algal nitrogen tracers  (n_algae)
             tr_bgc_DON_in   , & ! if .true., DON pools are tracers  (n_don)
             tr_bgc_C_in     , & ! if .true., algal carbon tracers + DOC and DIC 
             tr_bgc_chl_in   , & ! if .true., algal chlorophyll tracers 
             tr_bgc_Am_in    , & ! if .true., ammonia/um as nutrient tracer 
             tr_bgc_Sil_in   , & ! if .true., silicon as nutrient tracer 
             tr_bgc_DMS_in   , & ! if .true., DMS as product tracer 
             tr_bgc_Fe_in    , & ! if .true., Fe as product tracer 
             tr_bgc_hum_in   , & ! if .true., hum as product tracer 
             tr_bgc_PON_in       ! if .true., PON as product tracer 

        tr_iage      = tr_iage_in
        tr_FY        = tr_FY_in
        tr_lvl       = tr_lvl_in
        tr_pond      = tr_pond_in
        tr_pond_cesm = tr_pond_cesm_in
        tr_pond_lvl  = tr_pond_lvl_in
        tr_pond_topo = tr_pond_topo_in
        tr_aero      = tr_aero_in
        tr_brine     = tr_brine_in
        tr_bgc_S     = tr_bgc_S_in
        tr_zaero     = tr_zaero_in 
        tr_bgc_Nit   = tr_bgc_Nit_in
        tr_bgc_N     = tr_bgc_N_in 
        tr_bgc_DON   = tr_bgc_DON_in
        tr_bgc_C     = tr_bgc_C_in 
        tr_bgc_chl   = tr_bgc_chl_in
        tr_bgc_Am    = tr_bgc_Am_in
        tr_bgc_Sil   = tr_bgc_Sil_in
        tr_bgc_DMS   = tr_bgc_DMS_in
        tr_bgc_Fe    = tr_bgc_Fe_in 
        tr_bgc_hum   = tr_bgc_hum_in
        tr_bgc_PON   = tr_bgc_PON_in 

      end subroutine colpkg_init_tracer_flags

!=======================================================================

      subroutine colpkg_init_tracer_indices(&
           nt_Tsfc_in, & ! ice/snow temperature
           nt_qice_in, & ! volume-weighted ice enthalpy (in layers)
           nt_qsno_in, & ! volume-weighted snow enthalpy (in layers)
           nt_sice_in, & ! volume-weighted ice bulk salinity (CICE grid layers)
           nt_fbri_in, & ! volume fraction of ice with dynamic salt (hinS/vicen*aicen)
           nt_iage_in, & ! volume-weighted ice age
           nt_FY_in, & ! area-weighted first-year ice area
           nt_alvl_in, & ! level ice area fraction
           nt_vlvl_in, & ! level ice volume fraction
           nt_apnd_in, & ! melt pond area fraction
           nt_hpnd_in, & ! melt pond depth
           nt_ipnd_in, & ! melt pond refrozen lid thickness
           nt_aero_in, & ! starting index for aerosols in ice 
           nt_zaero_in,   & !  black carbon and other aerosols
           nt_bgc_N_in ,  & ! diatoms, phaeocystis, pico/small   
           nt_bgc_C_in ,  & ! diatoms, phaeocystis, pico/small   
           nt_bgc_chl_in, & ! diatoms, phaeocystis, pico/small 
           nt_bgc_DOC_in, & !  dissolved organic carbon
           nt_bgc_DON_in, & !  dissolved organic nitrogen
           nt_bgc_DIC_in, & !  dissolved inorganic carbon
           nt_bgc_Fed_in, & !  dissolved iron
           nt_bgc_Fep_in, & !  particulate iron
           nt_bgc_Nit_in, & ! nutrients  
           nt_bgc_Am_in,  & ! 
           nt_bgc_Sil_in, & !
           nt_bgc_DMSPp_in,&! trace gases (skeletal layer)
           nt_bgc_DMSPd_in,&! 
           nt_bgc_DMS_in, & ! 
           nt_bgc_hum_in, & ! 
           nt_bgc_PON_in, & ! zooplankton and detritus  
           nlt_zaero_in,  & !  black carbon and other aerosols
           nlt_bgc_N_in , & ! diatoms, phaeocystis, pico/small   
           nlt_bgc_C_in , & ! diatoms, phaeocystis, pico/small   
           nlt_bgc_chl_in,& ! diatoms, phaeocystis, pico/small 
           nlt_bgc_DOC_in,& !  dissolved organic carbon
           nlt_bgc_DON_in,& !  dissolved organic nitrogen
           nlt_bgc_DIC_in,& !  dissolved inorganic carbon
           nlt_bgc_Fed_in,& !  dissolved iron
           nlt_bgc_Fep_in,& !  particulate iron
           nlt_bgc_Nit_in,& ! nutrients  
           nlt_bgc_Am_in, & ! 
           nlt_bgc_Sil_in,& !
           nlt_bgc_DMSPp_in,&! trace gases (skeletal layer)
           nlt_bgc_DMSPd_in,&! 
           nlt_bgc_DMS_in,& ! 
           nlt_bgc_hum_in,& ! 
           nlt_bgc_PON_in,& ! zooplankton and detritus  
           nt_zbgc_frac_in,&! fraction of tracer in the mobile phase
           nt_bgc_S_in,   & ! Bulk salinity in fraction ice with dynamic salinity (Bio grid))
           nlt_chl_sw_in, & ! points to total chla in trcrn_sw
           nlt_zaero_sw_in,&! black carbon and dust in trcrn_sw
                            ! Index Dimensions: 
           n_algae, n_algalC, & !
           n_algalchl, n_DOC, & !
           n_DON,n_DIC,n_dFe, & !
           n_pFe, n_aerosols, & !
           bio_index_o_in,    & ! nlt index to fixed data array
           bio_index_in,      & ! nlt index to nt index
           nbtrcr)

        use ice_colpkg_tracers, only: &
             nt_Tsfc, & ! ice/snow temperature
             nt_qice, & ! volume-weighted ice enthalpy (in layers)
             nt_qsno, & ! volume-weighted snow enthalpy (in layers)
             nt_sice, & ! volume-weighted ice bulk salinity (CICE grid layers)
             nt_fbri, & ! volume fraction of ice with dynamic salt (hinS/vicen*aicen)
             nt_iage, & ! volume-weighted ice age
             nt_FY, & ! area-weighted first-year ice area
             nt_alvl, & ! level ice area fraction
             nt_vlvl, & ! level ice volume fraction
             nt_apnd, & ! melt pond area fraction
             nt_hpnd, & ! melt pond depth
             nt_ipnd, & ! melt pond refrozen lid thickness
             nt_aero, & ! starting index for aerosols in ice
             nt_zaero,   & !  black carbon and other aerosols
             nt_bgc_N ,  & ! diatoms, phaeocystis, pico/small   
             nt_bgc_C ,  & ! diatoms, phaeocystis, pico/small   
             nt_bgc_chl, & ! diatoms, phaeocystis, pico/small 
             nt_bgc_DOC, & !  dissolved organic carbon
             nt_bgc_DON, & !  dissolved organic nitrogen
             nt_bgc_DIC, & !  dissolved inorganic carbon
             nt_bgc_Fed, & !  dissolved iron
             nt_bgc_Fep, & !  particulate iron
             nt_bgc_Nit, & ! nutrients  
             nt_bgc_Am,  & ! 
             nt_bgc_Sil, & !
             nt_bgc_DMSPp,&! trace gases (skeletal layer)
             nt_bgc_DMSPd,&! 
             nt_bgc_DMS, & ! 
             nt_bgc_hum, & ! 
             nt_bgc_PON, & ! zooplankton and detritus 
             nlt_zaero,  & !  black carbon and other aerosols
             nlt_bgc_N , & ! diatoms, phaeocystis, pico/small   
             nlt_bgc_C , & ! diatoms, phaeocystis, pico/small   
             nlt_bgc_chl,& ! diatoms, phaeocystis, pico/small 
             nlt_bgc_DOC,& !  dissolved organic carbon
             nlt_bgc_DON,& !  dissolved organic nitrogen
             nlt_bgc_DIC,& !  dissolved inorganic carbon
             nlt_bgc_Fed,& !  dissolved iron
             nlt_bgc_Fep,& !  particulate iron
             nlt_bgc_Nit,& ! nutrients  
             nlt_bgc_Am, & ! 
             nlt_bgc_Sil,& !
             nlt_bgc_DMSPp,&! trace gases (skeletal layer)
             nlt_bgc_DMSPd,&! 
             nlt_bgc_DMS,& ! 
             nlt_bgc_hum,& ! 
             nlt_bgc_PON,& ! zooplankton and detritus   
             nt_zbgc_frac,&! fraction of tracer in the mobile phase
             nt_bgc_S,   & ! Bulk salinity in fraction ice with dynamic salinity (Bio grid))
             nlt_chl_sw, & ! points to total chla in trcrn_sw
             nlt_zaero_sw,&! black carbon and dust in trcrn_sw
             bio_index_o,& !
             bio_index 
        
        integer, intent(in) :: &
             nt_Tsfc_in, & ! ice/snow temperature
             nt_qice_in, & ! volume-weighted ice enthalpy (in layers)
             nt_qsno_in, & ! volume-weighted snow enthalpy (in layers)
             nt_sice_in, & ! volume-weighted ice bulk salinity (CICE grid layers)
             nt_fbri_in, & ! volume fraction of ice with dynamic salt (hinS/vicen*aicen)
             nt_iage_in, & ! volume-weighted ice age
             nt_FY_in, & ! area-weighted first-year ice area
             nt_alvl_in, & ! level ice area fraction
             nt_vlvl_in, & ! level ice volume fraction
             nt_apnd_in, & ! melt pond area fraction
             nt_hpnd_in, & ! melt pond depth
             nt_ipnd_in, & ! melt pond refrozen lid thickness
             nt_aero_in, & ! starting index for aerosols in ice
             nt_bgc_Nit_in, & ! nutrients  
             nt_bgc_Am_in,  & ! 
             nt_bgc_Sil_in, & !
             nt_bgc_DMSPp_in,&! trace gases (skeletal layer)
             nt_bgc_DMSPd_in,&! 
             nt_bgc_DMS_in, & ! 
             nt_bgc_hum_in, & ! 
             nt_bgc_PON_in, & ! zooplankton and detritus   
             nlt_bgc_Nit_in,& ! nutrients  
             nlt_bgc_Am_in, & ! 
             nlt_bgc_Sil_in,& !
             nlt_bgc_DMSPp_in,&! trace gases (skeletal layer)
             nlt_bgc_DMSPd_in,&! 
             nlt_bgc_DMS_in,& ! 
             nlt_bgc_hum_in,& ! 
             nlt_bgc_PON_in,& ! zooplankton and detritus  
             nt_zbgc_frac_in,&! fraction of tracer in the mobile phase
             nt_bgc_S_in,   & ! Bulk salinity in fraction ice with dynamic salinity (Bio grid))
             nlt_chl_sw_in    ! points to total chla in trcrn_sw

       integer, intent(in) :: &
             n_algae,    & !  Dimensions
             n_algalC,   & !
             n_algalchl, & !
             n_DOC,      & !
             n_DON,      & !
             n_DIC,      & !
             n_dFe,      & !
             n_pFe,      & ! 
             n_aerosols, & !
             nbtrcr

        integer (kind=int_kind), dimension(:), intent(in) :: &
             bio_index_o_in, & 
             bio_index_in  

        integer (kind=int_kind), dimension(:), intent(in) :: &
             nt_bgc_N_in ,  & ! diatoms, phaeocystis, pico/small   
             nt_bgc_C_in ,  & ! diatoms, phaeocystis, pico/small   
             nt_bgc_chl_in, & ! diatoms, phaeocystis, pico/small 
             nlt_bgc_N_in , & ! diatoms, phaeocystis, pico/small   
             nlt_bgc_C_in , & ! diatoms, phaeocystis, pico/small   
             nlt_bgc_chl_in   ! diatoms, phaeocystis, pico/small 

        integer (kind=int_kind), dimension(:), intent(in) :: &
             nt_bgc_DOC_in, & !  dissolved organic carbon
             nlt_bgc_DOC_in   !  dissolved organic carbon

        integer (kind=int_kind), dimension(:), intent(in) :: &
             nt_bgc_DON_in, & !  dissolved organic nitrogen
             nlt_bgc_DON_in   !  dissolved organic nitrogen

        integer (kind=int_kind), dimension(:), intent(in) :: &
             nt_bgc_DIC_in, & ! dissolved inorganic carbon
             nlt_bgc_DIC_in   !  dissolved inorganic carbon

        integer (kind=int_kind), dimension(:), intent(in) :: &
             nt_bgc_Fed_in, & !  dissolved iron
             nt_bgc_Fep_in, & !  particulate iron
             nlt_bgc_Fed_in,& !  dissolved iron
             nlt_bgc_Fep_in   !  particulate iron

        integer (kind=int_kind), dimension(:), intent(in) :: &
             nt_zaero_in,   & !  black carbon and other aerosols
             nlt_zaero_in,  & !  black carbon and other aerosols
             nlt_zaero_sw_in  ! black carbon and dust in trcrn_sw

        ! local
        integer (kind=int_kind) :: k

        nt_Tsfc = nt_Tsfc_in
        nt_qice = nt_qice_in
        nt_qsno = nt_qsno_in
        nt_sice = nt_sice_in
        nt_fbri = nt_fbri_in
        nt_iage = nt_iage_in
        nt_FY = nt_FY_in
        nt_alvl = nt_alvl_in
        nt_vlvl = nt_vlvl_in
        nt_apnd = nt_apnd_in
        nt_hpnd = nt_hpnd_in
        nt_ipnd = nt_ipnd_in
        nt_aero = nt_aero_in
        nt_bgc_Nit = nt_bgc_Nit_in
        nt_bgc_Am  = nt_bgc_Am_in
        nt_bgc_Sil = nt_bgc_Sil_in
        nt_bgc_DMSPp=nt_bgc_DMSPp_in
        nt_bgc_DMSPd=nt_bgc_DMSPd_in
        nt_bgc_DMS = nt_bgc_DMS_in
        nt_bgc_hum = nt_bgc_hum_in
        nt_bgc_PON = nt_bgc_PON_in
        nlt_bgc_Nit = nlt_bgc_Nit_in
        nlt_bgc_Am  = nlt_bgc_Am_in
        nlt_bgc_Sil = nlt_bgc_Sil_in
        nlt_bgc_DMSPp=nlt_bgc_DMSPp_in
        nlt_bgc_DMSPd=nlt_bgc_DMSPd_in
        nlt_bgc_DMS = nlt_bgc_DMS_in
        nlt_bgc_hum = nlt_bgc_hum_in
        nlt_bgc_PON = nlt_bgc_PON_in
        nlt_chl_sw  = nlt_chl_sw_in
        nt_zbgc_frac=nt_zbgc_frac_in
        nt_bgc_S   = nt_bgc_S_in

        nt_bgc_N(:)    = 0
        nt_bgc_C(:)    = 0
        nt_bgc_chl(:)  = 0
        nlt_bgc_N(:)   = 0
        nlt_bgc_C(:)   = 0
        nlt_bgc_chl(:) = 0
        nt_bgc_DOC(:)  = 0
        nlt_bgc_DOC(:) = 0
        nt_bgc_DIC(:)  = 0
        nlt_bgc_DIC(:) = 0
        nt_bgc_DON(:)  = 0
        nlt_bgc_DON(:) = 0
        nt_bgc_Fed(:)  = 0
        nt_bgc_Fep(:)  = 0
        nlt_bgc_Fed(:) = 0
        nlt_bgc_Fep(:) = 0
        nt_zaero(:)    = 0
        nlt_zaero(:)   = 0
        nlt_zaero_sw(:)= 0
        bio_index(:)   = 0
        bio_index_o(:) = 0

        do k = 1, nbtrcr
           bio_index_o(k)= bio_index_o_in(k)
           bio_index(k)  = bio_index_in(k)
        enddo
        do k = 1, n_algae
           nt_bgc_N(k) = nt_bgc_N_in(k) 
           nlt_bgc_N(k)= nlt_bgc_N_in(k) 
        enddo
        do k = 1, n_algalC
           nt_bgc_C(k) = nt_bgc_C_in(k) 
           nlt_bgc_C(k)= nlt_bgc_C_in(k) 
        enddo
        do k = 1, n_algalchl
           nt_bgc_chl(k) = nt_bgc_chl_in(k) 
           nlt_bgc_chl(k)= nlt_bgc_chl_in(k) 
        enddo
        do k = 1, n_DOC
           nt_bgc_DOC(k) = nt_bgc_DOC_in(k) 
           nlt_bgc_DOC(k)= nlt_bgc_DOC_in(k) 
        enddo
        do k = 1, n_DON
           nt_bgc_DON(k) = nt_bgc_DON_in(k) 
           nlt_bgc_DON(k)= nlt_bgc_DON_in(k) 
        enddo
        do k = 1, n_DIC
           nt_bgc_DIC(k) = nt_bgc_DIC_in(k) 
           nlt_bgc_DIC(k)= nlt_bgc_DIC_in(k) 
        enddo
        do k = 1, n_dFe  
           nt_bgc_Fed(k) = nt_bgc_Fed_in(k) 
           nlt_bgc_Fed(k)= nlt_bgc_Fed_in(k) 
        enddo
        do k = 1, n_pFe  
           nt_bgc_Fep(k) = nt_bgc_Fep_in(k) 
           nlt_bgc_Fep(k)= nlt_bgc_Fep_in(k) 
        enddo
        do k = 1, n_aerosols
           nt_zaero(k)    = nt_zaero_in(k)   
           nlt_zaero(k)   = nlt_zaero_in(k)   
           nlt_zaero_sw(k)= nlt_zaero_sw_in(k)   
        enddo

      end subroutine colpkg_init_tracer_indices

!=======================================================================
! set the number of column tracers

      subroutine colpkg_init_tracer_numbers(&
         ntrcr_in, nbtrcr_in, nbtrcr_sw_in)

      use ice_colpkg_tracers, only: &
         ntrcr, nbtrcr, nbtrcr_sw

      integer (kind=int_kind), intent(in) :: &
         ntrcr_in  , &! number of tracers in use
         nbtrcr_in , &! number of bio tracers in use
         nbtrcr_sw_in ! number of shortwave bio tracers in use
        
         ntrcr     = ntrcr_in
         nbtrcr    = nbtrcr_in
         nbtrcr_sw = nbtrcr_sw_in

      end subroutine colpkg_init_tracer_numbers

!=======================================================================

      subroutine colpkg_biogeochemistry(dt, &
                           ntrcr, nbtrcr,  &
                           upNO, upNH, iDi, iki, zfswin, &
                           zsal_tot, darcy_V, grow_net,  &
                           PP_net, hbri,dhbr_bot, dhbr_top, Zoo,&
                           fbio_snoice, fbio_atmice, ocean_bio, &
                           first_ice, fswpenln, bphi, bTiz, ice_bio_net,  &
                           snow_bio_net, fswthrun, Rayleigh_criteria, &
                           sice_rho, fzsal, fzsal_g, &
                           bgrid, igrid, icgrid, cgrid,  &
                           nblyr, nilyr, nslyr, n_algae, n_zaero, ncat, &
                           n_doc, n_dic,  n_don, n_fed, n_fep,  &
                           meltbn, melttn, congeln, snoicen, &
                           sst, sss, fsnow, meltsn, hmix, salinz, &
                           hin_old, flux_bio, flux_bio_atm, &
                           aicen_init, vicen_init, aicen, vicen, vsnon, &
                           aice0, trcrn, vsnon_init, skl_bgc, &
                           max_algae, max_nbtrcr, &
                           l_stop, stop_label)

      use ice_algae, only: zbio, sklbio
      use ice_brine, only: preflushing_changes, compute_microS_mushy, &
                           update_hbrine, compute_microS 
      use ice_colpkg_shared, only: solve_zsal, z_tracers, phi_snow
      use ice_colpkg_tracers, only: nt_fbri, tr_brine, &
          nt_bgc_S, nt_qice, nt_sice, nt_zbgc_frac, bio_index 
      use ice_constants_colpkg, only: c0, c1, puny
      use ice_zsalinity, only: zsalinity
      use ice_zbgc_shared, only:  zbgc_frac_init

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      integer (kind=int_kind), intent(in) :: &
         ncat, &
         nilyr, &
         nslyr, &
         nblyr, &
         ntrcr, &
         nbtrcr, &
         n_algae, n_zaero, &
         n_doc, n_dic,  n_don, n_fed, n_fep, &
         max_algae, max_nbtrcr

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         bgrid         , &  ! biology nondimensional vertical grid points
         igrid         , &  ! biology vertical interface points
         cgrid         , &  ! CICE vertical coordinate   
         icgrid        , &  ! interface grid for CICE (shortwave variable)
         ocean_bio     , &  ! contains all the ocean bgc tracer concentrations
         fbio_snoice   , &  ! fluxes from snow to ice
         fbio_atmice   , &  ! fluxes from atm to ice
         dhbr_top      , &  ! brine top change
         dhbr_bot      , &  ! brine bottom change
         darcy_V       , &  ! darcy velocity positive up (m/s)
         hin_old       , &  ! old ice thickness
         sice_rho      , &  ! avg sea ice density  (kg/m^3) 
         ice_bio_net   , &  ! depth integrated tracer (mmol/m^2) 
         snow_bio_net  , &  ! depth integrated snow tracer (mmol/m^2) 
         flux_bio     ! all bio fluxes to ocean

      logical (kind=log_kind), dimension (:), intent(inout) :: &
         first_ice      ! distinguishes ice that disappears (e.g. melts)
                        ! and reappears (e.g. transport) in a grid cell
                        ! during a single time step from ice that was
                        ! there the entire time step (true until ice forms)

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         Zoo            , & ! N losses accumulated in timestep (ie. zooplankton/bacteria)
                            ! mmol/m^3
         bphi           , & ! porosity of layers    
         bTiz           , & ! layer temperatures interpolated on bio grid (C)
         zfswin         , & ! Shortwave flux into layers interpolated on bio grid  (W/m^2)
         iDi            , & ! igrid Diffusivity (m^2/s)    
         iki            , & ! Ice permeability (m^2)  
         trcrn     ! tracers   

      real (kind=dbl_kind), intent(inout) :: &
         grow_net       , & ! Specific growth rate (/s) per grid cell
         PP_net         , & ! Total production (mg C/m^2/s) per grid cell
         hbri           , & ! brine height, area-averaged for comparison with hi (m)
         zsal_tot       , & ! Total ice salinity in per grid cell (g/m^2) 
         fzsal          , & ! Total flux  of salt to ocean at time step for conservation
         fzsal_g        , & ! Total gravity drainage flux
         upNO           , & ! nitrate uptake rate (mmol/m^2/d) times aice
         upNH         ! ammonium uptake rate (mmol/m^2/d) times aice

      logical (kind=log_kind), intent(inout) :: &
         Rayleigh_criteria    ! .true. means Ra_c was reached  

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         fswpenln        ! visible SW entering ice layers (W m-2)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         fswthrun    , & ! SW through ice to ocean            (W/m^2)
         meltsn      , & ! snow melt in category n (m)
         melttn      , & ! top melt in category n (m)
         meltbn      , & ! bottom melt in category n (m)
         congeln     , & ! congelation ice formation in category n (m)
         snoicen     , & ! snow-ice formation in category n (m)
         salinz      , & ! initial salinity  profile (ppt) 
         flux_bio_atm, & ! all bio fluxes to ice from atmosphere  
         aicen_init  , & ! initial ice concentration, for linear ITD
         vicen_init  , & ! initial ice volume (m), for linear ITD
         vsnon_init  , & ! initial snow volume (m), for aerosol 
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), intent(in) :: &
         aice0   , & ! open water area fraction
         sss     , & ! sea surface salinity (ppt)
         sst     , & ! sea surface temperature (C)
         hmix    , & ! mixed layer depth (m)
         fsnow       ! snowfall rate (kg/m^2 s)

      logical (kind=log_kind), intent(in) :: &
         skl_bgc       ! if true, solve skeletal biochemistry

      logical (kind=log_kind), intent(inout) :: &  
         l_stop          ! if true, abort the model

      character (len=*), intent(inout) :: stop_label

      ! local variables

      integer (kind=int_kind) :: &
         k              , & ! vertical index
         n, mm              ! thickness category index

      real (kind=dbl_kind) :: &
         hin         , & ! new ice thickness
         hsn         , & ! snow thickness  (m)
         hbr_old     , & ! old brine thickness before growh/melt
         dhice       , & ! change due to sublimation/condensation (m)
         kavg        , & ! average ice permeability (m^2)
         bphi_o      , & ! surface ice porosity 
         hbrin       , & ! brine height
         dh_direct       ! surface flooding or runoff

      real (kind=dbl_kind), dimension (nblyr+2) :: &
      ! Defined on Bio Grid points
         bSin        , & ! salinity on the bio grid  (ppt)
         brine_sal   , & ! brine salinity (ppt)
         brine_rho       ! brine_density (kg/m^3)

      real (kind=dbl_kind), dimension (nblyr+1) :: &
      ! Defined on Bio Grid interfaces
         iphin       , & ! porosity 
         ibrine_sal  , & ! brine salinity  (ppt)
         ibrine_rho  , & ! brine_density (kg/m^3)
         iTin            ! Temperature on the interface grid (oC)

      real (kind=dbl_kind) :: & 
         sloss            ! brine flux contribution from surface runoff (g/m^2)

      ! for bgc sk
      real (kind=dbl_kind) :: & 
         dh_bot_chl  , & ! Chlorophyll may or may not flush
         dh_top_chl  , & ! Chlorophyll may or may not flush
         darcy_V_chl     

      l_stop = .false.

      do n = 1, ncat

      !-----------------------------------------------------------------
      ! initialize
      !-----------------------------------------------------------------
         hin_old(n) = c0
         if (aicen_init(n) > puny) then 
            hin_old(n) = vicen_init(n) &
                                / aicen_init(n)
         else
            first_ice(n) = .true.
            if (tr_brine) trcrn(nt_fbri,n) = c1
            do mm = 1,nbtrcr
               trcrn(nt_zbgc_frac-1+mm,n) = zbgc_frac_init(mm)
            enddo
            if (n == 1) Rayleigh_criteria = .false.
            if (solve_zsal) trcrn(nt_bgc_S:nt_bgc_S+nblyr-1,n) = c0
         endif

         if (aicen(n) > puny) then
          
            dh_top_chl = c0
            dh_bot_chl = c0
            darcy_V_chl= c0
            bSin(:)    = c0
            hsn        = c0
            hin        = c0
            hbrin      = c0
            kavg       = c0
            bphi_o     = c0
            sloss      = c0
  
      !-----------------------------------------------------------------
      ! brine dynamics
      !-----------------------------------------------------------------

            dhbr_top(n) = c0
            dhbr_bot(n) = c0

            if (tr_brine) then 
               if (trcrn(nt_fbri,n) .le. c0) trcrn(nt_fbri,n) = c1

               dhice = c0
               call preflushing_changes  (n,  aicen  (n),   &
                                 vicen   (n), vsnon  (n),   &
                                 meltbn  (n), melttn (n),   &
                                 congeln (n), snoicen(n),   &
                                 hin_old (n), dhice,        & 
                                 trcrn(nt_fbri,n),          &
                                 dhbr_top(n), dhbr_bot(n),  &
                                 hbr_old,     hin,          &
                                 hsn,         first_ice(n), &
                                 l_stop,      stop_label)

               if (l_stop) return

               if (solve_zsal)  then  

                  call compute_microS (n,         nilyr,       nblyr,             &
                                bgrid,            cgrid,       igrid,             &
                                trcrn(1:ntrcr,n), hin_old(n),  hbr_old,           &
                                sss,              sst,         bTiz(:,n),         &
                                iTin,             bphi(:,n),   kavg,              &
                                bphi_o,           phi_snow,    Rayleigh_criteria, &
                                first_ice(n),     bSin,        brine_sal,         &
                                brine_rho,        iphin,       ibrine_rho,        &
                                ibrine_sal,       sice_rho(n), sloss,             &
                                salinz(1:nilyr),  l_stop,      stop_label)

                  if (l_stop) return
               else     

                 ! Requires the average ice permeability = kavg(:)
                 ! and the surface ice porosity = zphi_o(:)
                 ! computed in "compute_microS" or from "thermosaline_vertical"

                  iDi(:,n) = c0

                  call compute_microS_mushy (n,   nilyr,         nblyr,       &
                                   bgrid,         cgrid,         igrid,       &
                                   trcrn(:,n),    hin_old(n),    hbr_old,     &
                                   sss,           sst,           bTiz(:,n),   & 
                                   iTin(:),       bphi(:,n),     kavg,        &
                                   bphi_o,        phi_snow,      bSin(:),     &
                                   brine_sal(:),  brine_rho(:),  iphin(:),    &
                                   ibrine_rho(:), ibrine_sal(:), sice_rho(n), &
                                   iDi(:,n),      l_stop,        stop_label)

               endif ! solve_zsal  

               call update_hbrine (meltbn  (n), melttn(n),   &
                                   meltsn  (n), dt,          &
                                   hin,         hsn,         &
                                   hin_old (n), hbrin,       &

                                   hbr_old,     phi_snow,    &
                                   trcrn(nt_fbri,n),         &
                                   snoicen(n),               &
                                   dhbr_top(n), dhbr_bot(n), &
                                   dh_top_chl,  dh_bot_chl,  & 
                                   kavg,        bphi_o,      &
                                   darcy_V (n), darcy_V_chl, &  
                                   bphi(2,n),   aice0,       &
                                   dh_direct)
               
               hbri = hbri + hbrin * aicen(n)  

               if (solve_zsal) then 

                  call zsalinity (n,             dt,                  &
                                  nilyr,         bgrid,               & 
                                  cgrid,         igrid,               &
                                  trcrn(nt_bgc_S:nt_bgc_S+nblyr-1,n), &
                                  trcrn(nt_qice:nt_qice+nilyr-1,n),   &
                                  trcrn(nt_sice:nt_sice+nilyr-1,n),   &
                                  ntrcr,         trcrn(nt_fbri,n),    &
                                  bSin,          bTiz(:,n),           &
                                  bphi(:,n),     iphin,               &
                                  iki(:,n),      hbr_old,             &
                                  hbrin,         hin,                 &
                                  hin_old(n),    iDi(:,n),            &
                                  darcy_V(n),    brine_sal,           & 
                                  brine_rho,     ibrine_sal,          & 
                                  ibrine_rho,    dh_direct,           &
                                  Rayleigh_criteria,                  &
                                  first_ice(n),  sss,                 &
                                  sst,           dhbr_top(n),         &
                                  dhbr_bot(n),                        &
                                  l_stop,        stop_label,          &
                                  fzsal,         fzsal_g,             &
                                  bphi_o,        nblyr,               & 
                                  vicen(n),      aicen_init(n),       &
                                  zsal_tot) 

                  if (l_stop) return

               endif  ! solve_zsal

            endif ! tr_brine

      !-----------------------------------------------------------------
      ! biogeochemistry
      !-----------------------------------------------------------------

            if (z_tracers) then 
 
               call zbio (dt,                    nblyr,                  &
                          nslyr,                 nilyr,                  &
                          melttn(n),                                     &
                          meltsn(n),             meltbn  (n),            &
                          congeln(n),            snoicen(n),             & 
                          nbtrcr,                fsnow,                  &
                          ntrcr,                 trcrn(1:ntrcr,n),       &
                          bio_index(1:nbtrcr),   aicen_init(n),          &
                          vicen_init(n),         vsnon_init(n),          &
                          vicen(n),              vsnon(n),               &
                          aicen(n),              flux_bio_atm(1:nbtrcr), &
                          n,                     n_algae,                &
                          n_doc,                 n_dic,                  &
                          n_don,                                         &
                          n_fed,                 n_fep,                  &
                          n_zaero,               first_ice(n),           &
                          hin_old(n),            ocean_bio(1:nbtrcr),    &
                          bphi(:,n),             iphin,                  &     
                          iDi(:,n),              sss,                    &
                          fswpenln(:,n),                                 &
                          dhbr_top(n),           dhbr_bot(n),            &
                          dh_top_chl,            dh_bot_chl,             &
                          zfswin(:,n),                                   &
                          hbrin,                 hbr_old,                &
                          darcy_V(n),            darcy_V_chl,            &
                          bgrid,                 cgrid,                  &
                          igrid,                 icgrid,                 &
                          bphi_o,                                        &
                          dhice,                 iTin,                   &
                          Zoo(:,n),                                      &
                          flux_bio(1:nbtrcr),    dh_direct,              &
                          upNO,                  upNH,                   &
                          fbio_snoice,           fbio_atmice,            &
                          PP_net,                ice_bio_net (1:nbtrcr), &
                          snow_bio_net(1:nbtrcr),grow_net,               &
                          l_stop,                stop_label)
            
               if (l_stop) return
     
            elseif (skl_bgc) then

               call sklbio (dt,                      ntrcr,               &
                            nilyr,                                        &
                            nbtrcr,                  n_algae,             &
                            n_zaero,                 n_doc,               &
                            n_dic,                   n_don,               &
                            n_fed,                   n_fep,               &
                            flux_bio (1:nbtrcr),     ocean_bio(1:nbtrcr), &
                            hmix,                    aicen    (n),        &
                            meltbn   (n),            congeln  (n),        &
                            fswthrun (n),            first_ice(n),        &
                            trcrn    (1:ntrcr,n),    hin,                 &
                            PP_net,                  upNO,                &
                            upNH,                    grow_net,            &
                            l_stop,                  stop_label)

               if (l_stop) return

            endif  ! skl_bgc

            first_ice(n) = .false.

         endif             ! aicen > puny
      enddo                ! ncat

      end subroutine colpkg_biogeochemistry

!=======================================================================

!  Initialize brine height tracer

      subroutine colpkg_init_hbrine(bgrid, igrid, cgrid, &
          icgrid, swgrid, nblyr, nilyr, phi_snow)

      use ice_constants_colpkg, only: c1, c1p5, c2, p5, c0, rhoi, rhos
 
      integer (kind=int_kind), intent(in) :: &
         nilyr, & ! number of ice layers
         nblyr    ! number of bio layers

      real (kind=dbl_kind), intent(inout) :: &
         phi_snow           !porosity at the ice-snow interface

      real (kind=dbl_kind), dimension (nblyr+2), intent(out) :: &
         bgrid              ! biology nondimensional vertical grid points

      real (kind=dbl_kind), dimension (nblyr+1), intent(out) :: &
         igrid              ! biology vertical interface points
 
      real (kind=dbl_kind), dimension (nilyr+1), intent(out) :: &
         cgrid            , &  ! CICE vertical coordinate   
         icgrid           , &  ! interface grid for CICE (shortwave variable)
         swgrid                ! grid for ice tracers used in dEdd scheme

      integer (kind=int_kind) :: &
         k           , & ! vertical index
         n               ! thickness category index

      real (kind=dbl_kind) :: & 
         zspace            ! grid spacing for CICE vertical grid


      if (phi_snow .le. c0) phi_snow = c1-rhos/rhoi

      !-----------------------------------------------------------------
      ! Calculate bio gridn: 0 to 1 corresponds to ice top to bottom 
      !-----------------------------------------------------------------

      bgrid(:)       = c0 ! zsalinity grid points         
      bgrid(nblyr+2) = c1 ! bottom value
      igrid(:)       = c0 ! bgc interface grid points   
      igrid(1)       = c0 ! ice top
      igrid(nblyr+1) = c1 ! ice bottom
      
      zspace = c1/max(c1,(real(nblyr,kind=dbl_kind)))
      do k = 2, nblyr+1
         bgrid(k) = zspace*(real(k,kind=dbl_kind) - c1p5)
      enddo
      
      do k = 2, nblyr
         igrid(k) = p5*(bgrid(k+1)+bgrid(k))
      enddo

      !-----------------------------------------------------------------
      ! Calculate CICE cgrid for interpolation ice top (0) to bottom (1) 
      !-----------------------------------------------------------------
       
      cgrid(1) = c0                           ! CICE vertical grid top point
      zspace = c1/(real(nilyr,kind=dbl_kind)) ! CICE grid spacing
    
      do k = 2, nilyr+1
         cgrid(k) = zspace * (real(k,kind=dbl_kind) - c1p5) 
      enddo 

      !-----------------------------------------------------------------
      ! Calculate CICE icgrid for ishortwave interpolation top(0) , bottom (1)
      !-----------------------------------------------------------------
       
      icgrid(1) = c0                        
      zspace = c1/(real(nilyr,kind=dbl_kind)) ! CICE grid spacing
    
      do k = 2, nilyr+1
         icgrid(k) = zspace * (real(k,kind=dbl_kind)-c1)
      enddo 

      !------------------------------------------------------------------------
      ! Calculate CICE swgrid for dEdd ice: top of ice (0) , bottom of ice (1)
      ! Does not include snow
      ! see ice_shortwave.F90
      ! swgrid represents the layer index of the delta-eddington ice layer index
      !------------------------------------------------------------------------ 
      zspace = c1/(real(nilyr,kind=dbl_kind)) ! CICE grid spacing
      swgrid(1) = min(c1/60.0_dbl_kind, zspace/c2)      
      swgrid(2) = zspace/c2                   !+ swgrid(1)
      do k = 3, nilyr+1
         swgrid(k) = zspace * (real(k,kind=dbl_kind)-c1p5)
      enddo 

      end subroutine colpkg_init_hbrine

!=======================================================================

!  Initialize ocean concentration

      subroutine colpkg_init_ocean_conc (amm, dmsp, dms, algalN, doc, dic, don, &
             fed, fep, hum, nit, sil, zaeros, max_dic, max_don, max_fe, max_aero,&
             CToN, CToN_DON)

      use ice_constants_colpkg, only: c1,  c2, p5, c0, p1
      use ice_colpkg_shared, only: R_C2N, R_C2N_DON
 
      integer (kind=int_kind), intent(in) :: &
        max_dic, &
        max_don, &
        max_fe, &
        max_aero

      real (kind=dbl_kind), intent(out):: &
       amm      , & ! ammonium
       dmsp     , & ! DMSPp
       dms      , & ! DMS
       hum      , & ! humic material
       nit      , & ! nitrate
       sil          ! silicate

      real (kind=dbl_kind), dimension(:), intent(out):: &
       algalN   , & ! algae
       doc      , & ! DOC
       dic      , & ! DIC
       don      , & ! DON
       fed      , & ! Dissolved Iron
       fep      , & ! Particulate Iron
       zaeros       ! BC and dust

      real (kind=dbl_kind), dimension(:), intent(inout), optional :: &
       CToN     , & ! carbon to nitrogen ratio for algae
       CToN_DON     ! nitrogen to carbon ratio for proteins

      integer (kind=int_kind) :: &
        k 

       if (present(CToN)) then
         CToN(1) = R_C2N(1)
         CToN(2) = R_C2N(2)    
         CToN(3) = R_C2N(3)     
       endif

       if (present(CToN_DON)) then
         CToN_DON(1) = R_C2N_DON(1)
       endif

       amm  = c1 ! ISPOL < 1 mmol/m^3 
       dmsp = p1  
       dms  = p1    
       algalN(1) = c1  !0.0026_dbl_kind ! ISPOL, Lannuzel 2013(pennate) 
       algalN(2) = 0.0057_dbl_kind ! ISPOL, Lannuzel 2013(small plankton)
       algalN(3) = 0.0027_dbl_kind ! ISPOL, Lannuzel 2013(Phaeocystis)
                                     ! 0.024_dbl_kind ! 5% of 1 mgchl/m^3 
       doc(1) = 16.2_dbl_kind ! 18% saccharides
       doc(2) = 9.0_dbl_kind  ! lipids
       doc(3) = c1 ! 
       do k = 1, max_dic
            dic(k) = c1
       enddo  
       do k = 1, max_don
            don(k) = 12.9_dbl_kind              
            ! 64.3_dbl_kind ! 72% Total DOC~90 mmolC/m^3  ISPOL with N:C of 0.2
       enddo  
       !ki = 1
       !if (trim(fe_data_type) == 'clim') ki = 2
       do k = 1, max_fe ! ki, max_fe
            fed(k) = 0.4_dbl_kind ! c1 (nM) Lannuzel2007 DFe, 
                                  ! range 0.14-2.6 (nM) van der Merwe 2011
                                  ! Tagliabue 2012 (0.4 nM)
            fep(k) = c2 ! (nM) van der Merwe 2011
                        ! (0.6 to 2.9 nM ocean)
       enddo 
       hum  = c1        ! mmol C/m^3
       nit  = 12.0_dbl_kind
       sil  = 25.0_dbl_kind
       do k = 1, max_aero
         zaeros(k) = c0
       enddo
 

      end subroutine colpkg_init_ocean_conc

!=======================================================================

!  Initialize zSalinity

      subroutine colpkg_init_zsalinity(nblyr,ntrcr_o, restart_zsal,  Rayleigh_criteria, &
               Rayleigh_real, trcrn, nt_bgc_S, ncat, sss)

      use ice_constants_colpkg, only: c1,  c2, p5, c0, p1
      use ice_colpkg_shared, only: dts_b, salt_loss
 
      integer (kind=int_kind), intent(in) :: &
       nblyr, & ! number of biolayers
       ntrcr_o, & ! number of non bio tracers
       ncat , & ! number of categories
       nt_bgc_S ! zsalinity index

      logical (kind=log_kind), intent(in) :: &
       restart_zsal

      logical (kind=log_kind), intent(inout) :: &
       Rayleigh_criteria

      real (kind=dbl_kind), intent(inout):: &
       Rayleigh_real

      real (kind=dbl_kind), intent(in):: &
       sss

      real (kind=dbl_kind), dimension(:,:), intent(inout):: &
       trcrn ! bgc subset of trcrn

      integer (kind=int_kind) :: &
        k , n
      
      if (nblyr .LE. 7) then
          dts_b = 300.0_dbl_kind
      else
          dts_b = 50.0_dbl_kind 
      endif

      if (.not. restart_zsal) then
         Rayleigh_criteria = .false.    ! no ice initial condition 
         Rayleigh_real     = c0
         do n = 1,ncat
             do k = 1,nblyr
                trcrn(nt_bgc_S+k-1-ntrcr_o,n) = sss*salt_loss
             enddo   ! k
         enddo      ! n
      endif

      end subroutine colpkg_init_zsalinity

!=======================================================================

! basic initialization for ocean_bio_all

      subroutine colpkg_init_OceanConcArray(max_nbtrcr, &
          max_algae, max_don, max_doc, max_dic, max_aero, max_fe, &
          nit, amm, sil, dmsp, dms, algalN, &
          doc, don, dic, fed, fep, zaeros, ocean_bio_all, hum)

      use ice_constants_colpkg, only: c0
      use ice_colpkg_shared, only:  R_C2N, R_chl2N
      use ice_zbgc_shared, only: R_S2N

      integer (kind=int_kind), intent(in) :: &
         max_algae   , & ! maximum number of algal types 
         max_dic     , & ! maximum number of dissolved inorganic carbon types 
         max_doc     , & ! maximum number of dissolved organic carbon types
         max_don     , & ! maximum number of dissolved organic nitrogen types
         max_fe      , & ! maximum number of iron types
         max_aero    , & ! maximum number of aerosols 
         max_nbtrcr      ! maximum number of bio tracers

      real (kind=dbl_kind), intent(in) :: &
         nit         , & ! ocean nitrate (mmol/m^3)          
         amm         , & ! ammonia/um (mmol/m^3)
         sil         , & ! silicate (mmol/m^3)
         dmsp        , & ! dmsp (mmol/m^3)
         dms         , & ! dms (mmol/m^3)
         hum             ! humic material (mmol/m^3)

      real (kind=dbl_kind), dimension (max_algae), intent(in) :: &
         algalN          ! ocean algal nitrogen (mmol/m^3) (diatoms, phaeo, pico)

      real (kind=dbl_kind), dimension (max_doc), intent(in) :: &
         doc             ! ocean doc (mmol/m^3)  (proteins, EPS, lipid)

      real (kind=dbl_kind), dimension (max_don), intent(in) :: &
         don             ! ocean don (mmol/m^3) 

      real (kind=dbl_kind), dimension (max_dic), intent(in) :: &
         dic             ! ocean dic (mmol/m^3) 

      real (kind=dbl_kind), dimension (max_fe), intent(in) :: &
         fed, fep        ! ocean disolved and particulate fe (nM) 

      real (kind=dbl_kind), dimension (max_aero), intent(in) :: &
         zaeros          ! ocean aerosols (mmol/m^3) 

      real (kind=dbl_kind), dimension (max_nbtrcr), intent(inout) :: &
         ocean_bio_all   ! fixed order, all values even for tracers false

      ! local variables

      integer (kind=int_kind) :: &
         k, ks           ! tracer indices

      ocean_bio_all(:) = c0

      do k = 1, max_algae           
         ocean_bio_all(k)      = algalN(k)           ! N
         ks = max_algae + max_doc + max_dic + 1
         ocean_bio_all(ks + k) = R_chl2N(k)*algalN(k)!chl
      enddo   

      ks = max_algae + 1
      do k = 1, max_doc
         ocean_bio_all(ks + k) = doc(k)              ! doc
      enddo  
      ks = ks + max_doc
      do k = 1, max_dic
         ocean_bio_all(ks + k) = dic(k)              ! dic
      enddo 

      ks = 2*max_algae + max_doc + max_dic + 7
      do k = 1, max_don
         ocean_bio_all(ks + k) = don(k)              ! don
      enddo  

      ks = max_algae + 1
      ocean_bio_all(ks) = nit                        ! nit

      ks = 2*max_algae + max_doc + 2 + max_dic
      ocean_bio_all(ks) = amm                        ! Am
      ks = ks + 1
      ocean_bio_all(ks) = sil                        ! Sil
      ks = ks + 1
      ocean_bio_all(ks) =  R_S2N(1)*algalN(1) &      ! DMSPp
                        +  R_S2N(2)*algalN(2) &
                        +  R_S2N(3)*algalN(3) 
      ks = ks + 1
      ocean_bio_all(ks) = dmsp                       ! DMSPd
      ks = ks + 1
      ocean_bio_all(ks) = dms                        ! DMS
      ks = ks + 1
      ocean_bio_all(ks) = nit                        ! PON
      ks = 2*max_algae + max_doc + 7 + max_dic + max_don
      do k = 1, max_fe
         ocean_bio_all(ks + k) = fed(k)              ! fed
      enddo  
      ks = ks + max_fe
      do k = 1, max_fe
         ocean_bio_all(ks + k) = fep(k)              ! fep
      enddo  
      ks = ks + max_fe
      do k = 1, max_aero
         ocean_bio_all(ks+k) = zaeros(k)             ! zaero
      enddo
      ks = ks + max_aero + 1 
      ocean_bio_all(ks)  = hum                       ! humics

      end subroutine colpkg_init_OceanConcArray

!=======================================================================
! Warning messages
!=======================================================================

      subroutine colpkg_clear_warnings()

        use ice_warnings, only: reset_warnings

        call reset_warnings()

      end subroutine colpkg_clear_warnings

!=======================================================================
      
      subroutine colpkg_get_warnings(warningsOut)

        use ice_warnings, only: &
             get_number_warnings, &
             get_warning

        character(len=*), dimension(:), allocatable, intent(out) :: &
             warningsOut

        integer :: &
             iWarning

        if (allocated(warningsOut)) deallocate(warningsOut)
        allocate(warningsOut(get_number_warnings()))

        do iWarning = 1, get_number_warnings()
           warningsOut(iWarning) = trim(get_warning(iWarning))
        enddo

      end subroutine colpkg_get_warnings

!=======================================================================

      subroutine colpkg_print_warnings(nu_diag)

        use ice_warnings, only: &
             get_number_warnings, &
             get_warning

        integer, intent(in) :: nu_diag

        integer :: &
             iWarning

        do iWarning = 1, get_number_warnings()
           write(nu_diag,*) trim(get_warning(iWarning))
        enddo

      end subroutine colpkg_print_warnings

!=======================================================================

      end module ice_colpkg

!=======================================================================
