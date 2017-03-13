!  SVN:$Id: ice_zbgc.F90 1178 2017-03-08 19:24:07Z eclare $
!=======================================================================
!
! Biogeochemistry driver
!
! authors: Nicole Jeffery, LANL
!          Scott Elliot,   LANL
!          Elizabeth C. Hunke, LANL
!
      module ice_zbgc

      use ice_kinds_mod
      use ice_zbgc_shared ! everything

      implicit none 

      private
      public :: add_new_ice_bgc, lateral_melt_bgc, &
                merge_bgc_fluxes, merge_bgc_fluxes_skl

!=======================================================================

      contains

!=======================================================================

! Adjust biogeochemical tracers when new frazil ice forms

      subroutine add_new_ice_bgc (dt,         nblyr,                &
                                  ncat,       nilyr,      nltrcr,   &
                                  bgrid,      cgrid,      igrid,    &
                                  aicen_init, vicen_init, vi0_init, &
                                  aicen,      vicen,      vsnon1,   &
                                  vi0new,                           &
                                  ntrcr,      trcrn,      nbtrcr,   &
                                  sss,        ocean_bio,  flux_bio, &
                                  hsurp,      l_stop,   &
                                  stop_label, l_conservation_check)

      use ice_constants_colpkg, only: c0, c1, puny, depressT
      use ice_itd, only: column_sum, &
                         column_conservation_check
      use ice_colpkg_tracers, only: tr_brine, nt_fbri, nt_sice, nt_qice, nt_Tsfc
      use ice_colpkg_shared, only: solve_zsal
      use ice_therm_shared, only: calculate_Tin_from_qin

      integer (kind=int_kind), intent(in) :: &
         nblyr   , & ! number of bio layers
         ncat     , & ! number of thickness categories
         nilyr    , & ! number of ice layers
         nltrcr, & ! number of zbgc tracers
         nbtrcr  , & ! number of biology tracers
         ntrcr       ! number of tracers in use

      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         bgrid              ! biology nondimensional vertical grid points

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid              ! biology vertical interface points
 
      real (kind=dbl_kind), dimension (nilyr+1), intent(in) :: &
         cgrid              ! CICE vertical coordinate   

      real (kind=dbl_kind), intent(in) :: &
         dt              ! time step (s)

      real (kind=dbl_kind), dimension (:), &
         intent(in) :: &
         aicen_init  , & ! initial concentration of ice
         vicen_init  , & ! intiial volume per unit area of ice  (m)
         aicen       , & ! concentration of ice
         vicen           ! volume per unit area of ice          (m)

      real (kind=dbl_kind), intent(in) :: &
         vsnon1          ! category 1 snow volume per unit area (m)

      real (kind=dbl_kind), dimension (:,:), &
         intent(inout) :: &
         trcrn           ! ice tracers

      real (kind=dbl_kind), intent(in) :: &
         sss              !sea surface salinity (ppt)

      real (kind=dbl_kind), intent(in) :: &
         vi0_init    , & ! volume of new ice added to cat 1 (intial)
         vi0new          ! volume of new ice added to cat 1

      real (kind=dbl_kind), intent(in) :: &
         hsurp           ! thickness of new ice added to each cat

      real (kind=dbl_kind), dimension (:), &
         intent(inout) :: &
         flux_bio   ! tracer flux to ocean from biology (mmol/m^2/s) 
        
      real (kind=dbl_kind), dimension (:), &
         intent(in) :: &
         ocean_bio       ! ocean concentration of biological tracer

      logical (kind=log_kind), intent(in) :: &
         l_conservation_check

      logical (kind=log_kind), intent(inout) :: &
         l_stop     
        
      character (char_len), intent(inout) :: stop_label

! local

      integer (kind=int_kind) :: &
         location    , & ! 1 (add frazil to bottom), 0 (add frazil throughout)
         n           , & ! ice category index
         k               ! ice layer index

      real (kind=dbl_kind) :: &
         vbri1       , & ! starting volume of existing brine
         vbri_init   , & ! brine volume summed over categories
         vbri_final      ! brine volume summed over categories

      real (kind=dbl_kind) :: &
         vsurp       , & ! volume of new ice added to each cat
         vtmp            ! total volume of new and old ice
        
      real (kind=dbl_kind), dimension (ncat) :: &
         vbrin           ! trcrn(nt_fbri,n)*vicen(n) 
       
      real (kind=dbl_kind) :: &
         vice_new        ! vicen_init + vsurp

      real (kind=dbl_kind) :: &
         Tmlts       ! melting temperature (oC)

      character (len=char_len) :: &
         fieldid         ! field identifier

      !-----------------------------------------------------------------     
      ! brine
      !-----------------------------------------------------------------
      vbrin(:) = c0
      do n = 1, ncat
         vbrin(n) = vicen_init(n)
         if (tr_brine) vbrin(n) =  trcrn(nt_fbri,n)*vicen_init(n)
      enddo
     
      call column_sum (ncat,  vbrin,  vbri_init)

      vbri_init = vbri_init + vi0_init
      do k = 1, nbtrcr  
         flux_bio(k) = flux_bio(k) &
                            - vi0_init/dt*ocean_bio(k)*zbgc_init_frac(k)
      enddo
      !-----------------------------------------------------------------
      ! Distribute bgc in new ice volume among all ice categories by 
      ! increasing ice thickness, leaving ice area unchanged.
      !-----------------------------------------------------------------

         ! Diffuse_bio handles concentration changes from ice growth/melt
         ! ice area does not change
         ! add salt to the bottom , location = 1 

       vsurp = c0
       vtmp = c0

      do n = 1,ncat
 
      if (hsurp > c0) then

         vtmp = vbrin(n)
         vsurp = hsurp * aicen_init(n) 
         vbrin(n) = vbrin(n) + vsurp
         vice_new = vicen_init(n) + vsurp
         if (tr_brine .and. vicen(n) > c0) then
            trcrn(nt_fbri,n) = vbrin(n)/vicen(n)
         elseif (tr_brine .and. vicen(n) <= c0) then
            trcrn(nt_fbri,n) = c1
         endif

         if (nltrcr > 0) then 
            location = 1  
            call adjust_tracer_profile(nbtrcr,   dt, ntrcr, &
                                       aicen_init(n),       &
                                       vbrin(n),            &
                                       vice_new,            &
                                       trcrn(:,n),          &
                                       vtmp,                &
                                       vsurp,        sss,   &
                                       nilyr,        nblyr, &
                                       solve_zsal,   bgrid, & 
                                       cgrid,               &
                                       ocean_bio,    igrid, &
                                       location,            &
                                       l_stop,     stop_label)
            if (l_stop) return
         endif       ! nltrcr       
      endif          ! hsurp > 0
      enddo          ! n

      !-----------------------------------------------------------------
      ! Combine bgc in new ice grown in open water with category 1 ice.
      !-----------------------------------------------------------------
       
      if (vi0new > c0) then

         vbri1    = vbrin(1) 
         vbrin(1) = vbrin(1) + vi0new
         if (tr_brine .and. vicen(1) > c0) then
            trcrn(nt_fbri,1) = vbrin(1)/vicen(1)
         elseif (tr_brine .and. vicen(1) <= c0) then
            trcrn(nt_fbri,1) = c1
         endif
       
      ! Diffuse_bio handles concentration changes from ice growth/melt
      ! ice area changes
      ! add salt throughout, location = 0

         if (nltrcr > 0) then 
            location = 0  
            call adjust_tracer_profile(nbtrcr,  dt,     ntrcr,  &
                                       aicen(1),                &
                                       vbrin(1),                &
                                       vicen(1),                &
                                       trcrn(:,1),              &
                                       vbri1,                   &
                                       vi0new,          sss,    &
                                       nilyr,           nblyr,  &
                                       solve_zsal,      bgrid,  &
                                       cgrid,                   &
                                       ocean_bio,       igrid,  &
                                       location,                &
                                       l_stop,     stop_label)
            if (l_stop) return

            if (solve_zsal .and. vsnon1 .le. c0) then
               Tmlts = -trcrn(nt_sice,1)*depressT
               trcrn(nt_Tsfc,1) =  calculate_Tin_from_qin(trcrn(nt_qice,1),Tmlts)
            endif        ! solve_zsal 
         endif           ! nltrcr > 0
      endif              ! vi0new > 0

      if (tr_brine .and. l_conservation_check) then
         call column_sum (ncat,   vbrin,  vbri_final)

         fieldid = 'vbrin, add_new_ice_bgc'
         call column_conservation_check (fieldid,                  &
                                         vbri_init, vbri_final,    &
                                         puny,      l_stop)

         if (l_stop) then
            stop_label = 'add_new_ice_bgc: Column conservation error'
            return
         endif
      endif   ! l_conservation_check

      end subroutine add_new_ice_bgc

!=======================================================================

! When sea ice melts laterally, flux bgc to ocean

      subroutine lateral_melt_bgc (dt,                 &
                                   ncat,     nblyr,    &
                                   rside,    vicen,    &
                                   trcrn,    fzsal,    &
                                   flux_bio, nbltrcr)

      use ice_colpkg_tracers, only: nt_fbri, nt_bgc_S, bio_index
      use ice_colpkg_shared, only: solve_zsal, rhosi
      use ice_constants_colpkg, only: c1, p001

      integer (kind=int_kind), intent(in) :: &
         ncat  , & ! number of thickness categories
         nblyr , & ! number of bio layers
         nbltrcr   ! number of biology tracers

      real (kind=dbl_kind), intent(in) :: &
         dt        ! time step (s)

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         vicen     ! volume per unit area of ice          (m)

      real (kind=dbl_kind), dimension (:,:), intent(in) :: &
         trcrn     ! tracer array

      real (kind=dbl_kind), intent(in) :: &
         rside     ! fraction of ice that melts laterally

      real (kind=dbl_kind), intent(inout) :: &
         fzsal     ! salt flux from layer Salinity (kg/m^2/s)
  
      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         flux_bio  ! biology tracer flux from layer bgc (mmol/m^2/s)

      ! local variables

      integer (kind=int_kind) :: &
         k     , & ! layer index
         m     , & !
         n         ! category index

      real (kind=dbl_kind) :: &
         zspace    ! bio grid spacing

      zspace = c1/(real(nblyr,kind=dbl_kind))

      if (solve_zsal) then
         do n = 1, ncat
         do k = 1,nblyr
            fzsal = fzsal + rhosi*trcrn(nt_fbri,n) &
                  * vicen(n)*p001*zspace*trcrn(nt_bgc_S+k-1,n) &
                  * rside/dt
         enddo
         enddo
      endif

      do m = 1, nbltrcr
         do n = 1, ncat
         do k = 1, nblyr+1
            flux_bio(m) = flux_bio(m) + trcrn(nt_fbri,n) &
                        * vicen(n)*zspace*trcrn(bio_index(m)+k-1,n) &
                        * rside/dt
         enddo
         enddo
      enddo

      end subroutine lateral_melt_bgc 

!=======================================================================
!
! Add new ice tracers to the ice bottom and adjust the vertical profile 
!
! author: Nicole Jeffery, LANL

      subroutine adjust_tracer_profile (nbtrcr, dt, ntrcr, &
                                        aicen,      vbrin, &
                                        vicen,      trcrn, &
                                        vtmp,              &
                                        vsurp,      sss,   &
                                        nilyr,      nblyr, &
                                        solve_zsal, bgrid, &
                                        cgrid,      ocean_bio, &
                                        igrid,      location, &
                                        l_stop,     stop_label)

      use ice_constants_colpkg, only: c1, c0
      use ice_colpkg_tracers, only: nt_sice, nt_bgc_S, bio_index 
      use ice_colpkg_shared, only: min_salin, salt_loss

      integer (kind=int_kind), intent(in) :: &
         location          , & ! 1 (add frazil to bottom), 0 (add frazil throughout)
         ntrcr             , & ! number of tracers in use
         nilyr             , & ! number of ice layers
         nbtrcr            , & ! number of biology tracers
         nblyr                 ! number of biology layers

      real (kind=dbl_kind), intent(in) :: &
         dt              ! timestep (s)

      real (kind=dbl_kind), intent(in) :: &
         aicen   , & ! concentration of ice
         vicen   , & ! volume of ice
         sss     , & ! ocean salinity (ppt)
       ! hsurp   , & ! flags new ice added to each cat
         vsurp   , & ! volume of new ice added to each cat
         vtmp        ! total volume of new and old ice

      real (kind=dbl_kind), dimension (nbtrcr), intent(in) :: &
         ocean_bio

      real (kind=dbl_kind), intent(in) :: &
         vbrin       ! fbri*volume per unit area of ice  (m)
       
      logical (kind=log_kind), intent(in) :: &
         solve_zsal 

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid       ! zbio grid

      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         bgrid       ! zsal grid

      real (kind=dbl_kind), dimension (nilyr+1), intent(in) :: &
         cgrid       ! CICE grid

      real (kind=dbl_kind), dimension (ntrcr), &
         intent(inout) :: &
         trcrn       ! ice tracers
      
      logical (kind=log_kind), intent(inout) :: &
         l_stop            ! if true, print diagnostics and abort on return
        
      character (char_len), intent(inout) :: stop_label

      ! local variables

      real (kind=dbl_kind), dimension (ntrcr+2) :: &
         trtmp0, &      ! temporary, remapped tracers
         trtmp          ! temporary, remapped tracers
        
      real (kind=dbl_kind) :: &
         hin     , & ! ice height
         hinS_new, & ! brine height
         temp_S         

      integer (kind=int_kind) :: &
         k, m 

      real (kind=dbl_kind), dimension (nblyr+1) ::  &     
         C_stationary      ! stationary bulk concentration*h (mmol/m^2)

      real (kind=dbl_kind), dimension (nblyr) ::  &     
         S_stationary      ! stationary bulk concentration*h (ppt*m)

      real(kind=dbl_kind) :: &
         top_conc     , & ! salinity or bgc ocean concentration of frazil
         fluxb        , & ! needed for regrid (set to zero here)
         hbri_old     , & ! previous timestep brine height
         hbri             ! brine height 

      trtmp0(:) = c0
      trtmp(:) = c0
      fluxb = c0

      if (location == 1 .and. vbrin > c0) then  ! add frazil to bottom

         hbri     = vbrin
         hbri_old = vtmp
         if (solve_zsal) then
            top_conc = sss * salt_loss
            do k = 1, nblyr 
               S_stationary(k) = trcrn(nt_bgc_S+k-1)* hbri_old
            enddo
            call regrid_stationary (S_stationary, hbri_old, &
                                    hbri,         dt,       &
                                    ntrcr,                  &
                                    nblyr-1,      top_conc, &
                                    bgrid(2:nblyr+1), fluxb,&
                                    l_stop,       stop_label)
            if (l_stop) return
            do k = 1, nblyr 
               trcrn(nt_bgc_S+k-1) =  S_stationary(k)/hbri
               trtmp0(nt_sice+k-1) = trcrn(nt_bgc_S+k-1)
            enddo
         endif  ! solve_zsal

         do m = 1, nbtrcr
            top_conc = ocean_bio(m)*zbgc_init_frac(m)
            do k = 1, nblyr+1 
               C_stationary(k) = trcrn(bio_index(m) + k-1)* hbri_old
            enddo !k
            call regrid_stationary (C_stationary, hbri_old, &
                                    hbri,         dt,       &
                                    ntrcr,                  &
                                    nblyr,        top_conc, &
                                    igrid,        fluxb,    &
                                    l_stop,       stop_label)
            if (l_stop) return
            do k = 1, nblyr+1 
               trcrn(bio_index(m) + k-1) =  C_stationary(k)/hbri
            enddo !k                  
         enddo !m

         if (solve_zsal) then
            if (aicen > c0) then
               hinS_new  = vbrin/aicen
               hin       = vicen/aicen
            else
               hinS_new  = c0
               hin       = c0
            endif                   ! aicen
            temp_S    = min_salin   ! bio to cice
            call remap_zbgc(ntrcr,           nilyr,    &
                            nt_sice,                   &
                            trtmp0(1:ntrcr), trtmp,    &
                            1,               nblyr,    &
                            hin,             hinS_new, &
                            cgrid(2:nilyr+1),          &
                            bgrid(2:nblyr+1), temp_S,  &
                            l_stop,           stop_label)
            do k = 1, nilyr
               trcrn(nt_sice+k-1) = trtmp(nt_sice+k-1)   
            enddo        ! k
         endif           ! solve_zsal

      elseif (vbrin > c0) then   ! add frazil throughout  location == 0 .and.

         do k = 1, nblyr+1
            if (solve_zsal .and. k < nblyr + 1) then
               trcrn(nt_bgc_S+k-1) = (trcrn(nt_bgc_S+k-1) * vtmp &
                                          + sss*salt_loss * vsurp) / vbrin
               trtmp0(nt_sice+k-1) = trcrn(nt_bgc_S+k-1)
            endif                    ! solve_zsal
            do m = 1, nbtrcr
               trcrn(bio_index(m) + k-1) = (trcrn(bio_index(m) + k-1) * vtmp &
                         + ocean_bio(m)*zbgc_init_frac(m) * vsurp) / vbrin
            enddo
         enddo

         if (solve_zsal) then
            if (aicen > c0) then
               hinS_new  = vbrin/aicen
               hin       = vicen/aicen
            else
               hinS_new  = c0
               hin       = c0
            endif              !aicen
            temp_S    = min_salin   ! bio to cice
            call remap_zbgc(ntrcr,        nilyr,    &
                         nt_sice,                   &
                         trtmp0(1:ntrcr), trtmp,    &
                         1,               nblyr,    &
                         hin,             hinS_new, &
                         cgrid(2:nilyr+1),          &        
                         bgrid(2:nblyr+1),temp_S,   &
                         l_stop,           stop_label)
            do k = 1, nilyr
               trcrn(nt_sice+k-1) = trtmp(nt_sice+k-1)   
            enddo        !k
         endif   ! solve_zsal

      endif     ! location

      end subroutine adjust_tracer_profile

!=======================================================================
!
! Aggregate flux information from all ice thickness categories
! for z layer biogeochemistry
!
      subroutine merge_bgc_fluxes (dt,       nblyr,      &
                               bio_index,    n_algae,    &
                               nbtrcr,       aicen,      &    
                               vicen,        vsnon,      &
                               ntrcr,        iphin,      &
                               trcrn,      &
                               flux_bion,    flux_bio,   &
                               upNOn,        upNHn,      &
                               upNO,         upNH,       &
                               zbgc_snown,   zbgc_atmn,  &
                               zbgc_snow,    zbgc_atm,   &
                               PP_net,       ice_bio_net,&
                               snow_bio_net, grow_alg,   &
                               grow_net)
 
      use ice_constants_colpkg, only: c1, c0, p5, secday, puny
      use ice_colpkg_shared, only: solve_zbgc, max_nbtrcr, hs_ssl, R_C2N, &
                             fr_resp
      use ice_colpkg_tracers, only: nt_bgc_N, nt_fbri

      real (kind=dbl_kind), intent(in) :: &          
         dt             ! timestep (s)

      integer (kind=int_kind), intent(in) :: &
         nblyr, &
         n_algae, &     !
         ntrcr, &       ! number of tracers
         nbtrcr         ! number of biology tracer tracers

      integer (kind=int_kind), dimension(:), intent(in) :: &
         bio_index      ! relates bio indices, ie.  nlt_bgc_N to nt_bgc_N 

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         trcrn     , &  ! input tracer fields
         iphin          ! porosity

      real (kind=dbl_kind), intent(in):: &          
         aicen      , & ! concentration of ice
         vicen      , & ! volume of ice (m)
         vsnon          ! volume of snow(m)

      ! single category rates
      real (kind=dbl_kind), dimension(:), intent(in):: &
         zbgc_snown , & ! bio flux from snow to ice per cat (mmol/m^3*m) 
         zbgc_atmn  , & ! bio flux from atm to ice per cat (mmol/m^3*m)
         flux_bion

      ! single category rates
      real (kind=dbl_kind), dimension(:,:), intent(in):: &
         upNOn      , & ! nitrate uptake rate per cat (mmol/m^3/s)
         upNHn      , & ! ammonium uptake rate per cat (mmol/m^3/s)   
         grow_alg       ! algal growth rate per cat (mmolN/m^3/s)

      ! cumulative fluxes
      real (kind=dbl_kind), dimension(:), intent(inout):: &     
         flux_bio   , & ! 
         zbgc_snow  , & ! bio flux from snow to ice per cat (mmol/m^2/s) 
         zbgc_atm   , & ! bio flux from atm to ice per cat (mmol/m^2/s)
         ice_bio_net, & ! integrated ice tracers mmol or mg/m^2)
         snow_bio_net   ! integrated snow tracers mmol or mg/m^2)

      ! cumulative variables and rates
      real (kind=dbl_kind), intent(inout):: & 
         PP_net     , & ! net PP (mg C/m^2/d)  times aice
         grow_net   , & ! net specific growth (m/d) times vice
         upNO       , & ! tot nitrate uptake rate (mmol/m^2/d) times aice 
         upNH           ! tot ammonium uptake rate (mmol/m^2/d) times aice

      ! local variables

      real (kind=dbl_kind) :: &
         tmp        , & ! temporary
         dvssl      , & ! volume of snow surface layer (m)
         dvint          ! volume of snow interior      (m)

      integer (kind=int_kind) :: &
         k, mm         ! tracer indice

      real (kind=dbl_kind), dimension (nblyr+1) :: & 
         zspace

      !-----------------------------------------------------------------
      ! Column summation
      !-----------------------------------------------------------------
      zspace(:) = c1/real(nblyr,kind=dbl_kind)
      zspace(1) = p5/real(nblyr,kind=dbl_kind)
      zspace(nblyr+1) =  p5/real(nblyr,kind=dbl_kind)

      do mm = 1, nbtrcr
         do k = 1, nblyr+1
            ice_bio_net(mm) = ice_bio_net(mm) &
                            + trcrn(bio_index(mm)+k-1) &
                            * trcrn(nt_fbri) &
                            * vicen*zspace(k)
         enddo    ! k
      
      !-----------------------------------------------------------------
      ! Merge fluxes
      !-----------------------------------------------------------------
         dvssl  = min(p5*vsnon, hs_ssl*aicen) ! snow surface layer
         dvint  = vsnon - dvssl               ! snow interior
         snow_bio_net(mm) = snow_bio_net(mm) &
                          + trcrn(bio_index(mm)+nblyr+1)*dvssl &
                          + trcrn(bio_index(mm)+nblyr+2)*dvint
         flux_bio    (mm) = flux_bio (mm) + flux_bion (mm)*aicen
         zbgc_snow   (mm) = zbgc_snow(mm) + zbgc_snown(mm)*aicen/dt
         zbgc_atm    (mm) = zbgc_atm (mm) + zbgc_atmn (mm)*aicen/dt
      enddo     ! mm

      if (solve_zbgc) then
         do mm = 1, n_algae
            do k = 1, nblyr+1
               tmp      = iphin(k)*trcrn(nt_fbri)*vicen*zspace(k)*secday 
               PP_net   = PP_net   + grow_alg(k,mm)*tmp &
                        * (c1-fr_resp)* R_C2N(mm)*R_gC2molC 
               grow_net = grow_net + grow_alg(k,mm)*tmp &
                        / (trcrn(nt_bgc_N(mm)+k-1)+puny)
               upNO     = upNO     + upNOn   (k,mm)*tmp 
               upNH     = upNH     + upNHn   (k,mm)*tmp
            enddo   ! k
         enddo      ! mm
      endif

      end subroutine merge_bgc_fluxes

!=======================================================================

! Aggregate flux information from all ice thickness categories
! for skeletal layer biogeochemistry
!
! author: Elizabeth C. Hunke and William H. Lipscomb, LANL

      subroutine merge_bgc_fluxes_skl (ntrcr,           &
                               nbtrcr,    n_algae,         &
                               aicen,     trcrn,           &
                               flux_bion, flux_bio,        &
                               PP_net,    upNOn,           &
                               upNHn,     upNO,            &
                               upNH,      grow_net,        &
                               grow_alg)

      use ice_constants_colpkg, only: c1, secday, puny
      use ice_colpkg_tracers, only: nt_bgc_N
      use ice_colpkg_shared, only: sk_l, R_C2N, fr_resp

      integer (kind=int_kind), intent(in) :: &
         ntrcr   , & ! number of cells with aicen > puny
         nbtrcr  , & ! number of bgc tracers
         n_algae     ! number of autotrophs

      ! single category fluxes
      real (kind=dbl_kind), intent(in):: &          
         aicen       ! category ice area fraction

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         trcrn       ! Bulk tracer concentration (mmol N or mg/m^3)
     
      real (kind=dbl_kind), dimension(:), intent(in):: &
         flux_bion   ! all bio fluxes to ocean, on categories

      real (kind=dbl_kind), dimension(:), intent(inout):: &
         flux_bio    ! all bio fluxes to ocean, aggregated

      real (kind=dbl_kind), dimension(:), intent(in):: & 
         grow_alg, & ! algal growth rate (mmol/m^3/s) 
         upNOn   , & ! nitrate uptake rate per cat (mmol/m^3/s)
         upNHn       ! ammonium uptake rate per cat (mmol/m^3/s)   

      ! history output
      real (kind=dbl_kind), intent(inout):: & 
         PP_net  , & ! Bulk net PP (mg C/m^2/s)
         grow_net, & ! net specific growth (/s)
         upNO    , & ! tot nitrate uptake rate (mmol/m^2/s) 
         upNH        ! tot ammonium uptake rate (mmol/m^2/s)
      
      ! local variables

      integer (kind=int_kind) :: &
         k, mm       ! tracer indices

      real (kind=dbl_kind) :: &
         tmp         ! temporary
    
      !-----------------------------------------------------------------
      ! Merge fluxes
      !-----------------------------------------------------------------

      do k = 1,nbtrcr
         flux_bio (k) = flux_bio(k) + flux_bion(k)*aicen
      enddo

      do mm = 1, n_algae
         tmp = phi_sk * sk_l * aicen * secday 
         PP_net   = PP_net   &
                  + grow_alg(mm) * tmp &
                  * R_C2N(mm) * R_gC2molC * (c1-fr_resp) 
         grow_net = grow_net &
                  + grow_alg(mm) * tmp &
                  / (trcrn(nt_bgc_N(mm))+puny)
         upNO     = upNO  + upNOn(mm) * tmp
         upNH     = upNH  + upNHn(mm) * tmp
      enddo

      end subroutine merge_bgc_fluxes_skl

!=======================================================================

      end module ice_zbgc

!=======================================================================
