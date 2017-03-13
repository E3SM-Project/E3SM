!=======================================================================
!
! Vertical salinity (trcrn(nt_bgc_S)) is solved on the bio grid (bgrid and igrid)
! with domain defined by the dynamic brine height (trcrn(nt_fbri) * vicen/aicen).
! The CICE Bitz and Lipscomb thermodynamics is solved on the cgrid with height
! vicen/aicen.
! Gravity drainage is parameterized as nonlinear advection
! Flushing is incorporated in the boundary changes and a darcy flow. 
! (see Jeffery et al., JGR, 2011).  
!
! authors: Nicole Jeffery, LANL
!          Elizabeth C. Hunke, LANL
!
      module ice_zsalinity

      use ice_kinds_mod
      use ice_constants_colpkg
      use ice_zbgc_shared
      use ice_warnings, only: add_warning

      implicit none

      private
      public :: zsalinity

      real (kind=dbl_kind), parameter :: & 
         max_salin = 200.0_dbl_kind, & !(ppt) maximum bulk salinity
         lapidus_g = 0.3_dbl_kind  , & ! constant for artificial 
                                       ! viscosity/diffusion during growth
         lapidus_m = 0.007_dbl_kind    ! constant for artificial diffusion during melt

!=======================================================================

      contains

!=======================================================================

      subroutine zsalinity (n_cat,              dt,           &
                            nilyr,              bgrid,        &
                            cgrid,              igrid,        &
                            trcrn_S,            trcrn_q,      &
                            trcrn_Si,           ntrcr,        &
                            fbri,                             &
                            bSin,               bTin,         &
                            bphin,              iphin,        &
                            ikin,               hbr_old,      &
                            hbrin,              hin,          &
                            hin_old,            iDin,         &
                            darcy_V,            brine_sal,    &
                            brine_rho,          ibrine_sal,   &
                            ibrine_rho,                       &
                            Rayleigh_criteria,                &
                            first_ice,          sss,          &
                            sst,                dh_top,       &
                            dh_bot,                           &
                            l_stop,             stop_label,   &
                            fzsal,                            &
                            fzsal_g,            bphi_min,     &
                            nblyr,              vicen,        &
                            aicen,              zsal_tot) 
             
      use ice_constants_colpkg, only: c0, c1, puny
 
      integer (kind=int_kind), intent(in) :: &
         nilyr         , & ! number of ice layers
         nblyr         , & ! number of bio layers
         ntrcr         , & ! number of tracers
         n_cat             ! category number 
                    
      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         bgrid              ! biology nondimensional vertical grid points

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid              ! biology vertical interface points
 
      real (kind=dbl_kind), dimension (nilyr+1), intent(in) :: &
         cgrid              ! CICE vertical coordinate   

      real (kind=dbl_kind), intent(in) :: &
         sss           , & ! ocean salinity (ppt)
         sst           , & ! ocean temperature (oC)
         hin_old       , & ! old ice thickness (m)
         dh_top        , & ! brine change in top and bottom for diagnostics (m)
         dh_bot        , & ! minimum porosity
         darcy_V       , & ! darcy velocity (m/s)
         dt            , & ! time step
         fbri          , & ! ratio of brine height to ice thickness
         hbr_old       , & ! old brine height  (m)
         hin           , & ! new ice thickness (m)
         hbrin         , & ! new brine height  (m)
         vicen         , & ! ice volume (m)
         aicen         , & ! ice area (m)
         bphi_min   

      real (kind=dbl_kind), intent(inout) :: &
         zsal_tot      , & ! tot salinity (psu*rhosi*total vol ice)
         fzsal         , & ! total flux of salt out of ice over timestep(kg/m^2/s)
         fzsal_g           ! gravity drainage flux of salt over timestep(kg/m^2/s)

      real (kind=dbl_kind), dimension (nblyr+2), intent(inout) :: &
         bTin          , & ! Ice Temperature ^oC (on bio grid)
         bphin             ! Ice porosity (on bio grid)

      real (kind=dbl_kind), dimension (nblyr+2), intent(inout) :: &
         bSin          , & ! Ice salinity ppt (on bio  grid)
         brine_sal     , & ! brine salinity (ppt)
         brine_rho         ! brine density  (kg/m^3)

      real (kind=dbl_kind), dimension (nblyr), &
         intent(inout) :: &
         trcrn_S           ! salinity tracer ppt (on bio grid)

      real (kind=dbl_kind), dimension (nilyr), &
         intent(inout) :: &
         trcrn_q       , & ! enthalpy tracer 
         trcrn_Si          ! salinity on CICE grid

      logical (kind=log_kind), intent(inout) :: &
         Rayleigh_criteria ! .true. if minimun ice thickness (Ra_c)  was reached 
      
      logical (kind=log_kind), intent(in) :: &
         first_ice         ! for first category ice only .true. 
                           !initialized values should be used 

      real (kind=dbl_kind), dimension (nblyr+1), intent(out) :: &
         iDin          , & ! Diffusivity on the igrid   (1/s)
         ikin              ! permeability on the igrid 

      real (kind=dbl_kind), dimension (nblyr+1), intent(inout) :: &
         iphin         , & ! porosity on the igrid 
         ibrine_rho    , & ! brine rho on interface  
         ibrine_sal        ! brine sal on interface   
         
      logical (kind=log_kind), intent(inout) :: &
         l_stop            ! if true, print diagnostics and abort on return

      character (char_len) :: stop_label

      ! local variables

      integer (kind=int_kind) :: &
         k              , & ! vertical index
         n, mm              ! thickness category index

      real (kind=dbl_kind) :: &
         fzsaln        , & ! category flux of salt out of ice over timestep(kg/m^2/s)
         fzsaln_g      , & ! category gravity drainage flux of salt over timestep(kg/m^2/s)
         zsal_totn         ! total salt content

      call solve_zsalinity        (nilyr, nblyr, n_cat, dt,          &
                                   bgrid,      cgrid,    igrid,      &
                                   trcrn_S,            trcrn_q,      &
                                   trcrn_Si,           ntrcr,        &
                                   bSin,               bTin,         &
                                   bphin,              iphin,        &
                                   ikin,               hbr_old,      &
                                   hbrin,              hin,          &
                                   hin_old,            iDin,         &
                                   darcy_V,            brine_sal,    &
                                   brine_rho,          ibrine_sal,   &
                                   ibrine_rho,                       &
                                   Rayleigh_criteria,                &
                                   first_ice,          sss,          &
                                   sst,                dh_top,       &
                                   dh_bot,                           &
                                   l_stop,             stop_label,   &
                                   fzsaln,                           &
                                   fzsaln_g,           bphi_min)

         zsal_totn = c0

         call column_sum_zsal    (zsal_totn,          nblyr,      &
                                  vicen,              trcrn_S,    &
                                  fbri)

         call merge_zsal_fluxes (aicen, &
                                 zsal_totn,          zsal_tot,       &
                                 fzsal,              fzsaln,         &
                                 fzsal_g,            fzsaln_g)            

      end subroutine zsalinity

!=======================================================================
!
! update vertical salinity 
! 
      subroutine solve_zsalinity  (nilyr,              nblyr,        &
                                   n_cat,              dt,           &
                                   bgrid,    cgrid,    igrid,        &
                                   trcrn_S,            trcrn_q,      &
                                   trcrn_Si,           ntrcr,        &
                                   bSin,               bTin,         &
                                   bphin,              iphin,        &
                                   ikin,               hbr_old,      &
                                   hbrin,              hin,          &
                                   hin_old,            iDin,         &
                                   darcy_V,            brine_sal,    &
                                   brine_rho,          ibrine_sal,   &
                                   ibrine_rho,                       &
                                   Rayleigh_criteria,                &
                                   first_ice,          sss,          &
                                   sst,                dh_top,       &
                                   dh_bot,                           &  
                                   l_stop,             stop_label,   &
                                   fzsaln,                           &
                                   fzsaln_g,           bphi_min)

      use ice_colpkg_tracers, only: nt_sice
      use ice_colpkg_shared, only: solve_zsal, min_salin, dts_b, rhosi
      use ice_therm_shared, only: calculate_Tin_from_qin

      integer (kind=int_kind), intent(in) :: &
         nilyr,          & ! number of ice layers
         nblyr,          & ! number of bio layers
         ntrcr,          & ! number of tracers
         n_cat             ! category number 
                    
      real (kind=dbl_kind), intent(in) :: &
         dt                ! time step

      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         bgrid              ! biology nondimensional vertical grid points

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid              ! biology vertical interface points
 
      real (kind=dbl_kind), dimension (nilyr+1), intent(in) :: &
         cgrid              ! CICE vertical coordinate   

      real (kind=dbl_kind), intent(in) :: &
         sss           , & ! ocean salinity (ppt)
         sst           , & ! ocean temperature (oC)
         hin_old       , & ! old ice thickness (m)
         dh_top        , & ! brine change in top and bottom for diagnostics (m)
         dh_bot        , &
         darcy_V    

      real (kind=dbl_kind), intent(in) :: &
         hbr_old       , & ! old brine height  (m)
         hin           , & ! new ice thickness (m)
         hbrin         , & ! new brine height  (m)
         bphi_min   
 
      real (kind=dbl_kind), intent(out) :: &
         fzsaln        , & ! total flux of salt out of ice over timestep(kg/m^2/s)
         fzsaln_g          ! gravity drainage flux of salt  over timestep(kg/m^2/s)

      real (kind=dbl_kind), dimension (nblyr+2), intent(inout) :: &
         bTin          , & ! Ice Temperature ^oC (on bio grid)
         bphin             ! Ice porosity (on bio grid)

      real (kind=dbl_kind), dimension (nblyr+2), intent(inout) :: &
         bSin          , & ! Ice salinity ppt (on bio  grid)
         brine_sal     , & ! brine salinity (ppt)
         brine_rho         ! brine density  (kg/m^3)

      real (kind=dbl_kind), dimension (nblyr), &
         intent(inout) :: &
         trcrn_S           ! salinity tracer ppt (on bio grid)

      real (kind=dbl_kind), dimension (nilyr), &
         intent(inout) :: &
         trcrn_q       , & ! enthalpy tracer 
         trcrn_Si          ! salinity on CICE grid

      logical (kind=log_kind), intent(inout) :: &
         Rayleigh_criteria ! .true. if minimun ice thickness (Ra_c)  was reached 
      
      logical (kind=log_kind), intent(in) :: &
         first_ice         ! for first category ice only .true. 
                           !initialized values should be used 

      real (kind=dbl_kind), dimension (nblyr+1), intent(out) :: &
         iDin          , & ! Diffusivity on the igrid   (1/s)
         ikin              ! permeability on the igrid 

      real (kind=dbl_kind), dimension (nblyr+1), intent(inout) :: &
         iphin         , & ! porosity on the igrid 
         ibrine_rho    , & ! brine rho on interface  
         ibrine_sal        ! brine sal on interface   
         
      logical (kind=log_kind), intent(inout) :: &
         l_stop            ! if true, print diagnostics and abort on return

      character (char_len) :: stop_label

      ! local variables

      integer (kind=int_kind) :: &
         k, m, nint        ! vertical biology layer index 

      real (kind=dbl_kind) :: &
         surface_S         ! salinity of ice above hin > hbrin
      
      real (kind=dbl_kind), dimension(2) :: &
         S_bot           

      real (kind=dbl_kind) :: &
         Tmlts         , & ! melting temperature
         dts               ! local timestep (s)

      logical (kind=log_kind) :: &
         Rayleigh
     
      real (kind=dbl_kind):: &
         Ttemp             ! initial temp profile on the CICE grid

      real (kind=dbl_kind), dimension (ntrcr+2) :: &
         trtmp0        , & ! temporary, remapped tracers     !need extra 
         trtmp             ! temporary, remapped tracers     !

      logical (kind=log_kind) :: &
         cflag

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      dts = dts_b
      nint = max(1,INT(dt/dts))  
      dts = dt/nint

      l_stop = .false.
     
      !----------------------------------------------------------------
      ! Update boundary conditions
      !----------------------------------------------------------------
         
      surface_S = min_salin

      Rayleigh = .true.
      if (n_cat == 1 .AND. hbr_old < Ra_c) then
         Rayleigh = Rayleigh_criteria ! only category 1 ice can be false 
      endif

      if (dh_bot + darcy_V*dt > c0) then  
            
         bSin     (nblyr+2) = sss
         bTin     (nblyr+2) = sst
         brine_sal(nblyr+2) = sss 
         brine_rho(nblyr+2) = rhow
         bphin    (nblyr+2) = c1 
         S_bot    (1)       = c0
         S_bot    (2)       = c1  
   
      ! bottom melt
      else  
         bSin (nblyr+2) = bSin(nblyr+1)  
         Tmlts          =  -bSin(nblyr+2)* depressT 
         bTin (nblyr+2) = bTin(nblyr+1)
         bphin(nblyr+2) = iphin(nblyr+1)
         S_bot(1)       = c1
         S_bot(2)       = c0
      endif

      if (abs(dh_top) > puny .AND. abs(darcy_V) > puny) then 
         bSin(1) = max(min_salin,-(brine_rho(2)*brine_sal(2)/rhosi &
                 * darcy_V*dt - (dh_top +darcy_V*dt/bphi_min)*min_salin)/dh_top)
         brine_sal(1) = brine_sal(2)
         brine_rho(1) = brine_rho(2)
         bphin(1) = bphi_min
      else
         bSin(1) = min_salin
      endif
            
      !-----------------------------------------------------------------
      ! Solve for S using CICE T.  If solve_zsal = .true., then couple back
      ! to the thermodynamics
      !-----------------------------------------------------------------

      call solve_S_dt (cflag,         nblyr,         &
                       nint         , dts          , &
                       bSin         , bTin         , &
                       bphin        , iphin        , &
                       igrid        , bgrid        , &
                       ikin         ,                &
                       hbr_old      , hbrin        , &
                       hin          , hin_old      , &
                       iDin         , darcy_V      , &
                       brine_sal    , Rayleigh     , &
                       first_ice    , sss          , &
                       dt           , dh_top       , &
                       dh_bot       , brine_rho    , &
                       ibrine_sal   , ibrine_rho   , &
                       fzsaln       , fzsaln_g     , &
                       S_bot        , l_stop       , &
                       stop_label) 

      if (l_stop) return
  
      if (n_cat == 1)   Rayleigh_criteria = Rayleigh

      trtmp0(:) = c0
      trtmp (:) = c0
       
      do k = 1,nblyr                  ! back to bulk quantity 
            trcrn_S(k) =   bSin(k+1) 
            trtmp0(nt_sice+k-1) = trcrn_S(k)
      enddo           ! k

      call remap_zbgc   (ntrcr,     nilyr, &
                         nt_sice,          &
                         trtmp0(1:ntrcr),  &
                         trtmp,            &
                         1,         nblyr, &
                         hin,       hbrin, &
                         cgrid(2:nilyr+1), &
                         bgrid(2:nblyr+1), &
                         surface_S, l_stop,&
                         stop_label)
               
      do k = 1, nilyr
            Tmlts = -trcrn_Si(k)*depressT
            Ttemp = min(-(min_salin+puny)*depressT, &
                       calculate_Tin_from_qin(trcrn_q(k),Tmlts)) 
            trcrn_Si(k) = min(-Ttemp/depressT, max(min_salin, &
                                 trtmp(nt_sice+k-1)))
            Tmlts = - trcrn_Si(k)*depressT 
           ! if (cflag)  trcrn_q(k) = calculate_qin_from_Sin(Ttemp,Tmlts)  
      enddo ! k

      end subroutine solve_zsalinity

!=======================================================================
!
!  solves salt continuity explicitly using 
!  Lax-Wendroff-type scheme (MacCormack)
!  (Mendez-Nunez and Carroll,  Monthly Weather Review, 1993)
!
! authors     Nicole Jeffery, LANL
!
      subroutine solve_S_dt          (cflag,  nblyr, nint,          &
                                      dts,    bSin,  bTin,          &
                                      bphin,  iphin, igrid,         &
                                      bgrid,  ikin,  hbri_old,      &
                                      hbrin,  hice,  hice_old,      &
                                      iDin,          darcy_V,       &
                                      brine_sal,     Rayleigh,      &
                                      first_ice,     sss,           &
                                      dt,            dht,           &
                                      dhb,           brine_rho,     &
                                      ibrine_sal,    ibrine_rho,    &
                                      fzsaln,        fzsaln_g,      &
                                      S_bot,         l_stop,        &
                                      stop_label)    

      use ice_brine, only: calculate_drho
      use ice_colpkg_shared, only: l_skS, grid_oS, l_sk, min_salin, rhosi, salt_loss

      integer (kind=int_kind), intent(in) :: &
         nblyr            , & ! number of bio layers
         nint                 ! number of interations

      logical (kind=log_kind), intent(out) :: &
         cflag                ! thin or not

      real (kind=dbl_kind), intent(in) :: &
         dt               , & ! timestep (s)
         dts              , & ! local timestep (s)
         sss              , & ! sea surface salinity
         dht              , & ! change in the ice top  (positive for melting)
         dhb              , & ! change in the ice bottom (positive for freezing)
         hice_old         , & ! old ice thickness (m)
         hbri_old         , & ! brine thickness (m) 
         hbrin            , & ! new brine thickness (m)
         hice             , & ! ice thickness (m
         darcy_V              ! Darcy velocity due to a pressure head (m/s) or melt      

      real (kind=dbl_kind), intent(out) :: &
         fzsaln           , & ! salt flux +ive to  ocean (kg/m^2/s)
         fzsaln_g             ! gravity drainage salt flux +ive to ocean (kg/m^2/s)

      logical (kind=log_kind), intent(inout) :: &
         Rayleigh             ! if .true. convection is allowed; if .false. not yet

      logical (kind=log_kind), intent(in) :: &
         first_ice

      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         brine_sal        , & ! Internal brine salinity (ppt)
         brine_rho        , & ! Internal brine density (kg/m^3)
         bgrid            , & ! biology nondimensional grid layer points 
         bTin                 ! Temperature of ice layers on bio grid for history (C) 

      real (kind=dbl_kind), dimension (nblyr+2), intent(inout) :: &
         bphin            , & ! Porosity of layers
         bSin                 ! Bulk Salinity (ppt) contains previous timestep
                              ! and ocean ss

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         ibrine_rho       , & ! brine rho on interface 
         ibrine_sal       , & ! brine sal on interface 
         igrid                ! biology grid interface points 

      real (kind=dbl_kind), dimension (nblyr+1), intent(inout) :: &
         iphin                ! Porosity of layers on interface

      real (kind=dbl_kind), dimension (nblyr+1), intent(out) :: &
         iDin             , & ! Diffusivity on the igrid (1/s) with minimum bphi condition
         ikin                 ! permeability on interface

      real (kind=dbl_kind), dimension (2), intent(in) :: &
         S_bot
       
      logical (kind=log_kind), intent(out) :: &
         l_stop               ! if true, print diagnostics and abort on return

      character (char_len) :: stop_label

      ! local variables

      integer (kind=int_kind) :: &
         k, m , mm            ! vertical biology layer index 

      real (kind=dbl_kind), dimension (nblyr+1) :: &
         iDin_p           , & ! Diffusivity on the igrid (1/s)/bphi^3 
         dSbdx            , & ! gradient of brine rho on grid
         drho             , & ! brine difference rho_a-rho_b  (kg/m^3)
         Ci_s             , & !
         Ui_s             , & ! interface function
         Vi_s             , & ! for conservation check
         ivel

      real (kind=dbl_kind), dimension (nblyr+2) :: &
         Din_p            , & ! Diffusivity on the igrid (1/s)/bphi^3  
         Sintemp          , & ! initial salinity
         pre_sin          , & ! estimate of  salinity of layers
         pre_sinb         , & ! estimate of  salinity of layers
         bgrid_temp       , & ! biology nondimensional grid layer points 
                              ! with boundary values 
         Q_s, C_s         , & ! Functions in continuity equation
         V_s, U_s, F_s    

      real (kind=dbl_kind) :: &
         dh               , & ! (m) change in hbrine over dts
         dbgrid           , & ! ratio of grid space to spacing across boundary 
                              ! i.e. 1/nilyr/(dbgrid(2)-dbgrid(1))
         lapidus          , & ! artificial viscosity:  use lapidus_g for growth
         Ssum_old,Ssum_new, & ! depth integrated salt before and after timestep
         fluxcorr,          & ! flux correction to prevent S < min_salin
         Ssum_corr,         & ! numerical boundary flux correction
         fluxb            , & ! bottom, top and total boundary flux (g/kg/m^2)
         fluxg            , & ! bottom, top and total gravity drainage flux (g/kg/m^2)
         fluxm            , & ! bottom, top and total molecular diffusion flux (g/kg/m^2)
         sum_old,sum_new  , & ! integrated salinity at t and t+dt
         dh_dt, dS_dt     , &
         Ssum_tmp

      real (kind=dbl_kind), dimension (nblyr) :: &
         vel              , & ! advective velocity times dt (m)
         lapidus_diff     , & ! lapidus term and 
         flux_corr        , &
         lapA             , &
         lapB    

         logical (kind=log_kind) :: &   
         write_flag       , &    ! set to true at each timestep        
         test_conservation       ! test that salt change is balanced by fluxes 

      !-----------------------------------------------------------------
      !  Initialize
      !-----------------------------------------------------------------

      l_stop = .false.
      cflag = .false.
      write_flag = .true.
      test_conservation = .false. 
      iDin_p(:) = c0   
      Din_p(:) = c0 
      lapA(:) = c1
      lapB(:) = c1
      lapA(nblyr) = c0
      lapB(1) = c0
      V_s(:) = c0
      U_s(:) = c0
      Q_s(:) = c0
      C_s(:) = c0
      Ci_s(:) = c0
      F_s(:) = c0
      Ui_s(:) = c0
      Vi_s(:) = c0
      ivel(:) = c0
      vel(:) = c0
      dh = c0
      dbgrid = c2
      fzsaln = c0
      fzsaln_g = c0

      !-----------------------------------------------------------------
      ! Find brine density gradient for gravity drainage parameterization
      !-----------------------------------------------------------------

      call calculate_drho(nblyr, igrid, bgrid,&
                          brine_rho, ibrine_rho, drho)   

      !-----------------------------------------------------------------
      ! Calculate bphi diffusivity on the grid points
      ! rhosi = 919-974 kg/m^2  set in bio_in
      ! rhow = 1026.0 density of sea water: uses kinematic viscosity (m^2/s) in Q18
      ! dynamic viscosity  divided by density = kinematic. 
      !-----------------------------------------------------------------

      do k = 2, nblyr+1
         iDin_p(k) = k_o*gravit*l_skS/viscos_dynamic*drho(k)/(hbri_old**2)
          Din_p(k) = (iDin_p(k)*(igrid(k)-bgrid(k)) &
                   + iDin_p(k-1)*(bgrid(k)-igrid(k-1)))/(igrid(k)-igrid(k-1))
      enddo                !k

      !-----------------------------------------------------------------
      ! Critical Ra_c value is only for the onset of convection in thinS 
      ! ice and not throughout, therefore I need a flag to indicate the 
      ! Ra_c was reached for a particular ice block.
      ! Using a thickness minimum (Ra_c) for simplicity.
      !-----------------------------------------------------------------

      bgrid_temp(:) = bgrid(:)
      Din_p(nblyr+2) =  iDin_p(nblyr+1)
      if (.NOT. Rayleigh .AND. hbrin < Ra_c) then
         Din_p(:) =  c0  
         iDin_p(:) = c0  
      else
         Rayleigh = .true.
      endif

      if (hbri_old > thinS .AND. hbrin > thinS .and. &
          hice_old > thinS .AND. .NOT. first_ice) then

         cflag = .true.

         bgrid_temp(1) = c0
         bgrid_temp(nblyr+2) = c1
         dbgrid = igrid(2)/(bgrid_temp(2)-bgrid_temp(1))

         !-----------------------------------
         ! surface boundary terms
         !-----------------------------------
                             
         lapidus = lapidus_g/real(nblyr,kind=dbl_kind)**2
         ivel(1) = dht/hbri_old
         U_s (1) = ivel(1)/dt*dts 
         Ui_s(1) = U_s(1) 
         Ci_s(1) = c0
         F_s (1) = brine_rho(2)*brine_sal(2)/rhosi*darcy_V*dts/hbri_old/bSin(1)

         !-----------------------------------
         ! bottom boundary terms
         !-----------------------------------

         ivel(nblyr+1) =  dhb/hbri_old           
         Ui_s(nblyr+1) = ivel(nblyr+1)/dt*dts  
         dSbdx(nblyr) = (ibrine_sal(nblyr+1)*ibrine_rho(nblyr+1) &
                      -  ibrine_sal(nblyr)*ibrine_rho(nblyr)) &
                      / (igrid(nblyr+1)-igrid(nblyr))	              
         C_s(nblyr+1) = Dm/brine_sal(nblyr+1)/brine_rho(nblyr+1)*dts/hbri_old**2 &
                      * (ibrine_sal(nblyr+1)*ibrine_rho(nblyr+1) &
                      -  ibrine_sal(nblyr)*ibrine_rho(nblyr)) &
                      / (igrid(nblyr+1)-igrid(nblyr))
         F_s(nblyr+1) = darcy_V*dts/hbri_old/bphin(nblyr+1)
         F_s(nblyr+2) = darcy_V*dts/hbri_old/bphin(nblyr+2)  
         vel(nblyr) =(bgrid(nblyr+1)*(dhb) -(bgrid(nblyr+1) - c1)*dht)/hbri_old
         U_s(nblyr+1) = vel(nblyr)/dt*dts  
         V_s(nblyr+1) = Din_p(nblyr+1)/rhosi &
                      * (rhosi/brine_sal(nblyr+1)/brine_rho(nblyr+1))**exp_h &
                      * dts*dSbdx(nblyr) 
         dSbdx(nblyr+1) = (brine_sal(nblyr+2)*brine_rho(nblyr+2) &
                        -  brine_sal(nblyr+1)*brine_rho(nblyr+1)) &
                        / (bgrid(nblyr+2)-bgrid(nblyr+1)+ grid_oS/hbri_old )  
         C_s( nblyr+2) = Dm/brine_sal(nblyr+2)/brine_rho(nblyr+2)*dts/hbri_old**2 &
                       * (brine_sal(nblyr+2)*brine_rho(nblyr+2) & 
                       -  brine_sal(nblyr+1)*brine_rho(nblyr+1)) &
                       / (bgrid(nblyr+2)-bgrid(nblyr+1) + grid_oS/hbri_old ) 
         U_s(nblyr+2) = ivel(nblyr+1)/dt*dts 
         V_s(nblyr+2) = Din_p(nblyr+2)/rhosi &
                      * (bphin(nblyr+1)/bSin(nblyr+2))**exp_h &
                      * dts*dSbdx(nblyr+1)
         Ci_s(nblyr+1) = C_s(nblyr+2)
         Vi_s(nblyr+1) = V_s(nblyr+2) 
         dh = (dhb-dht)/dt*dts

         do k = 2, nblyr  
            ivel(k) =  (igrid(k)*dhb - (igrid(k)-c1)*dht)/hbri_old
            Ui_s(k) = ivel(k)/dt*dts   
            Vi_s(k) = iDin_p(k)/rhosi &
                    *(rhosi/ibrine_rho(k)/ibrine_sal(k))**exp_h*dts &
                    * (brine_sal(k+1)*brine_rho(k+1) &
                    -  brine_sal(k)*brine_rho(k)) &
                    / (bgrid(k+1)-bgrid(k)) 
            dSbdx(k-1) = (ibrine_sal(k)*ibrine_rho(k) &
                       -  ibrine_sal(k-1)*ibrine_rho(k-1))/(igrid(k)-igrid(k-1))
            F_s(k) = darcy_V*dts/hbri_old/bphin(k)
            C_s(k) = Dm/brine_sal(k)/brine_rho(k)*dts/hbri_old**2 &
                   * (ibrine_sal(k)*ibrine_rho(k) &
                   -  ibrine_sal(k-1)*ibrine_rho(k-1))/(igrid(k)-igrid(k-1))
            Ci_s(k) = Dm/ibrine_sal(k)/ibrine_rho(k)*dts/hbri_old**2 &
                    * (brine_sal(k+1)*brine_rho(k+1) &
                    -  brine_sal(k)*brine_rho(k))/(bgrid(k+1)-bgrid(k))
            vel(k-1) = (bgrid(k)*(dhb) - (bgrid(k) - c1)* dht)/hbri_old
            U_s(k) = vel(k-1)/dt*dts 
            V_s(k) = Din_p(k)/rhosi &
                   * (rhosi/brine_sal(k)/brine_rho(k))**exp_h*dts*dSbdx(k-1) 
            C_s(2) = c0
            V_s(2) = c0
         enddo !k

      !-----------------------------------------------------------------
      ! Solve
      !-----------------------------------------------------------------

         do m = 1, nint     

            Sintemp(:) = bSin(:)
            pre_sin(:) = bSin(:)  
            pre_sinb(:) = bSin(:)
            Ssum_old = bSin(nblyr+1)*(igrid(nblyr+1)-igrid(nblyr))

            ! forward-difference 

            do k = 2, nblyr
               Ssum_old = Ssum_old + bSin(k)*(igrid(k)-igrid(k-1))

               pre_sin(k) =bSin(k) + (Ui_s(k)*(bSin(k+1) - bSin(k)) + &
                 V_s(k+1)*bSin(k+1)**3 - V_s(k)*bSin(k)**3 + &
                 (C_s(k+1)+F_s(k+1))*bSin(k+1)-&
                 (C_s(k)+F_s(k))*bSin(k))/(bgrid_temp(k+1)-bgrid_temp(k)) 
            enddo    !k

            pre_sin(nblyr+1) = bSin(nblyr+1) + (Ui_s(nblyr+1)*(bSin(nblyr+2) - &
                 bSin(nblyr+1)) +  V_s(nblyr+2)*bSin(nblyr+2)**3 - &
                 V_s(nblyr+1)*bSin(nblyr+1)**3+ (C_s(nblyr+2)+F_s(nblyr+2))*&
                 bSin(nblyr+2)-(C_s(nblyr+1)+F_s(nblyr+1))*bSin(nblyr+1) )/&
                 (bgrid_temp(nblyr+2)- bgrid_temp(nblyr+1))
              
            ! backward-difference 

            pre_sinb(2) = p5*(bSin(2) + pre_sin(2) +  (Ui_s(1)&
                  *(pre_sin(2) -pre_sin(1)) + &
                  V_s(2)*pre_sin(2)**3 - &
                  V_s(1)*pre_sin(1)**3 + (C_s(2)+F_s(2))*pre_sin(2)-&
                  (C_s(1)+F_s(1))*pre_sin(1) )/(bgrid_temp(2)-bgrid_temp(1)) )
              
            do k = nblyr+1, 3, -1  !nblyr+1
               pre_sinb(k) =p5*(bSin(k) + pre_sin(k) +  (Ui_s(k-1)&
                  *(pre_sin(k) - pre_sin(k-1)) + &
                  V_s(k)*pre_sin(k)**3 - &
                  V_s(k-1)*pre_sin(k-1)**3 + (C_s(k)+F_s(k))*pre_sin(k)-&
                  (C_s(k-1)+F_s(k-1))*pre_sin(k-1))/(bgrid_temp(k)-bgrid_temp(k-1)) )

               Q_s(k) = V_s(k)*pre_sin(k)**2 + U_s(k) + C_s(k) + F_s(k) 
            enddo   !k

            Q_s(2) = V_s(2)*pre_sin(2)**2 + U_s(2) + C_s(2) + F_s(2)
            Q_s(1) = V_s(1)*pre_sin(2)**2 + Ui_s(1) + C_s(1)+ F_s(1)
            Q_s(nblyr+2) = V_s(nblyr+2)*pre_sin(nblyr+1)**2 + & 
            Ui_s(nblyr+1) + C_s(nblyr+2) +  F_s(nblyr+2)

      !-----------------------------------------------------------------
      ! Add artificial viscosity   [Lapidus,1967] [Lohner et al, 1985]
      ! * more for melting ice
      !-----------------------------------------------------------------

            lapidus_diff(:) = c0
            flux_corr(:) = c0
            Ssum_new = c0
            Ssum_corr = c0
            fluxcorr = c0
            fluxg = c0
            fluxb = c0
            fluxm = c0

            do k = 2, nblyr+1  

               lapidus_diff(k-1) =    lapidus/& ! lapidus/real(nblyr,kind=dbl_kind)**2/&
                  (igrid(k)-igrid(k-1))* &
                  ( lapA(k-1)*ABS(Q_s(k+1)-Q_s(k))*(abs(pre_sinb(k+1))-abs(pre_sinb(k)))/&
                  (bgrid_temp(k+1)-bgrid_temp(k) )**2 - &
                  lapB(k-1)*ABS(Q_s(k)-Q_s(k-1))*(abs(pre_sinb(k))-abs(pre_sinb(k-1)))/&
                  (bgrid_temp(k)-bgrid_temp(k-1))**2)
                          
               bSin(k) = pre_sinb(k) + lapidus_diff(k-1)

               if (bSin(k) < min_salin) then
                  flux_corr(k-1) = min_salin - bSin(k) !  flux into the ice
                  bSin(k) = min_salin 
               elseif (bSin(k) > -bTin(k)/depressT) then
                  flux_corr(k-1) = bSin(k)+bTin(k)/depressT !  flux into the ice
                  bSin(k) = -bTin(k)/depressT
               elseif (bSin(k) > max_salin) then
                  l_stop = .true.
                  stop_label = 'bSin(k) > max_salin'
               endif
            
               if (k == nblyr+1) bSin(nblyr+2) = S_bot(1)*bSin(nblyr+1) &
                                                  + S_bot(2)*bSin(nblyr+2) 

               Ssum_new = Ssum_new + bSin(k)*(igrid(k)-igrid(k-1))
               fluxcorr = fluxcorr + (flux_corr(k-1) &
                        + lapidus_diff(k-1))*(igrid(k)-igrid(k-1))

            enddo   !k

            Ssum_tmp = Ssum_old

            call calc_salt_fluxes (nint, m, nblyr, igrid, &
                 Ui_s,   dh,dbgrid,hbri_old,Sintemp,    &
                 pre_sin,   fluxb,fluxg,fluxm,V_s,    &
                 C_s,   F_s,   Ssum_corr,fzsaln_g,fzsaln, &
                 Ssum_tmp,fluxcorr,dts, Ssum_new)

            if (test_conservation) then
               call check_conserve_salt(nint, m, dt, dts,&
                                Ssum_tmp, Ssum_new, Ssum_corr,&
                                fluxcorr, fluxb, fluxg, fluxm, &
                                hbrin, hbri_old, l_stop)
               stop_label = 'check_conserve_salt fails'
               if (l_stop) return
            endif  ! test_conservation

         enddo !m

      else    !  add/melt ice only 

         sum_old = c0
         sum_new = c0
         dh_dt = hbrin-hbri_old
         dS_dt = c0
         if (dh_dt > c0) then 
            dS_dt = sss*dh_dt*salt_loss
            do k = 2, nblyr+1 
               bSin(k) = max(min_salin,(bSin(k)*hbri_old + dS_dt)/hbrin)
            enddo  !k
         else
            do k = 2, nblyr+1 
               sum_old = sum_old + bSin(k)*hbri_old*(igrid(k)-igrid(k-1))
               bSin(k) = max(min_salin,bSin(k) * (c1-abs(dh_dt)/hbri_old))
               sum_new = sum_new + bSin(k)*hbrin*(igrid(k)-igrid(k-1))
            enddo  !k
         endif
         fzsaln = rhosi*(sum_old-sum_new + dS_dt)*p001/dt   !kg/m^2/s
         fzsaln_g = c0

      endif  ! (hbri_old > thinS .AND. hbrin > thinS &
             ! .and.  hice_old > thinS .AND. .NOT. first_ice) 

      !-----------------------------------------------------------------
      ! Move this to bgc calculation if using tr_salinity
      ! Calculate bphin, iphin, ikin, iDin and iDin_N 
      !-----------------------------------------------------------------

      iDin(:) = c0
      iphin(:)  = c1
      ikin(:) = c0    

      do k = 1, nblyr+1
         if (k < nblyr+1) bphin(k+1) = min(c1,max(puny, &
                          bSin(k+1)*rhosi/(brine_sal(k+1)*brine_rho(k+1)))) 
         if (k == 1) then     
            bphin(k) = min(c1,max(puny, bSin(k)*rhosi/(brine_sal(k)*brine_rho(k))))  
            iphin(k) = bphin(2)
         elseif (k == nblyr+1) then
            iphin(nblyr+1) = bphin(nblyr+1)
         else
            iphin(k) = min(c1, max(c0,(bphin(k+1) - bphin(k))/(bgrid(k+1) &
                           - bgrid(k))*(igrid(k)-bgrid(k)) + bphin(k)))
         endif 
         ikin(k) = k_o*iphin(k)**exp_h 
      enddo    !k

      if (cflag) then
               
         do k = 2, nblyr+1
            iDin(k) =  iphin(k)*Dm/hbri_old**2  
            if (Rayleigh .AND. hbrin .GE. Ra_c) &
            iDin(k) = iDin(k) + l_sk*ikin(k)*gravit/viscos_dynamic &
                         * drho(k)/hbri_old**2 
         enddo       !k
      else   !  .not. cflag
         do k = 2, nblyr+1
            iDin(k) = iphin(k)*Dm/hbri_old**2
         enddo       !k
      endif

      end subroutine solve_S_dt

!=======================================================================
!
! Calculate salt fluxes
! 
      subroutine calc_salt_fluxes (mmax, mint, nblyr, igrid, &
                                   Ui_s,dh,dbgrid,hbri_old,Sintemp,pre_sin,&
                                   fluxb,fluxg,fluxm,V_s,&
                                   C_s,F_s,Ssum_corr,fzsaln_g,fzsaln,Ssum_old, &
                                   fluxcorr,dts, Ssum_new)

      use ice_colpkg_shared, only: rhosi

      integer(kind=int_kind), intent(in) :: &
         nblyr,          & ! number of bio layers
         mint ,          & ! current iteration
         mmax              ! total number of iterations     

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: &
         igrid              ! biology vertical interface points
 
      real (kind=dbl_kind), intent(in) :: &
         dts      , & !  halodynamic timesteps (s)
        ! hbrin    , & ! new brine height after all iterations (m)
         dh       , &  ! (m) change in hbrine over dts
         dbgrid   , & ! ratio of grid space to spacing across boundary 
                      ! ie. 1/nilyr/(dbgrid(2)-dbgrid(1))
         fluxcorr , & ! flux correction to ensure S >= min_salin
         hbri_old     ! initial brine height (m)

      real (kind=dbl_kind), dimension (nblyr+1), intent(in) :: & 
         Ui_s         ! interface function

      real (kind=dbl_kind), dimension (nblyr+2), intent(in) :: &
         Sintemp  , & ! initial salinity
         pre_sin  , & ! estimate of salinity of layers
         C_s      , & ! Functions in continuity equation
         F_s      , &
         V_s

      real (kind=dbl_kind), intent(in) :: &
         Ssum_old , & ! initial integrated salt content (ppt)/h  
         Ssum_new     ! next integrated salt content(ppt)/h

      real (kind=dbl_kind), intent(inout) :: &
         fzsaln   , & ! total salt flux  out of ice over timestep(kg/m^2/s)
         fzsaln_g , & ! gravity drainage flux of salt  over timestep(kg/m^2/s)
         Ssum_corr, & ! boundary flux correction due to numerics
         fluxb    , & ! total boundary salt flux into the ice (+ into ice)
         fluxg    , & ! total gravity drainage salt flux into the ice (+ into ice)
         fluxm        ! total molecular diffusive salt flux into the ice (+ into ice)

      ! local variables

      real (kind=dbl_kind) :: &
         Ssum_corr_flux  , & ! numerical boundary flux correction
         fluxb_b, fluxb_t, & ! bottom, top and total boundary flux (g/kg/m^2)
         fluxg_b, fluxg_t, & ! bottom, top and total gravity drainage flux (g/kg/m^2)
         fluxm_b, fluxm_t    ! bottom, top and total molecular diffusion flux (g/kg/m^2)

      real (kind=dbl_kind) :: hin_old, hin_next, dhtmp !, dh

      dhtmp = c1-dh/hbri_old
      hin_next = hbri_old +  real(mint,kind=dbl_kind)*dh
      hin_old  = hbri_old + (real(mint,kind=dbl_kind)-c1)*dh

      !-----------------------------------------------------------------
      ! boundary fluxes (positive into the ice)
      !---------------------------------------------
      !  without higher order numerics corrections 
      ! fluxb = (Ui_s(nblyr+1) + F_s(nblyr+2))*Sintemp(nblyr+2)  - (Ui_s(1) + F_s(1))*Sintemp(1)
      !-----------------------------------------------------------------

      fluxb_b = p5* Ui_s(nblyr+1) * (dhtmp*Sintemp(nblyr+2)*dbgrid &
                                  +        pre_sin(nblyr+1) &
                                  +  dhtmp*Sintemp(nblyr+1)*(c1-dbgrid)) &
              + p5*(F_s(nblyr+2)  *  dhtmp*Sintemp(nblyr+2)*dbgrid &
                  + F_s(nblyr+1)  *       (pre_sin(nblyr+1) &
                                  +  dhtmp*Sintemp(nblyr+1)*(c1-dbgrid)))

      fluxb_t = -p5*Ui_s(1)*(pre_sin(1)*dbgrid +  &
                              dhtmp*Sintemp(2) -  &
                              (dbgrid-c1)*pre_sin(2)) - &
                              p5*(dbgrid*F_s(1)*pre_sin(1) + &
                              F_s(2)*(dhtmp*Sintemp(2) &
                              +(c1-dbgrid)*pre_sin(2)))

      fluxb = fluxb_b + fluxb_t 

      !-----------------------------------------------------------------
      ! gravity drainage fluxes (positive into the ice)
      ! without higher order numerics corrections 
      ! fluxg =  V_s(nblyr+2)*Sintemp(nblyr+1)**3
      !-----------------------------------------------------------------

      fluxg_b =  p5*(dhtmp* dbgrid* &
                               V_s(nblyr+2)*Sintemp(nblyr+1)**3  +  &
                               V_s(nblyr+1)*(pre_sin(nblyr+1)**3 - &
                               dhtmp*(dbgrid - c1)* &
                               Sintemp(nblyr+1)**3))
                
      fluxg_t =  -p5*(dbgrid*V_s(1)*pre_sin(1)**3 + &
                               V_s(2)*(dhtmp*Sintemp(2)**3- &
                               (dbgrid-c1)*pre_sin(2)**3))
                
      fluxg =  fluxg_b + fluxg_t
      
      !-----------------------------------------------------------------
      ! diffusion fluxes (positive into the ice)
      ! without higher order numerics corrections 
      ! fluxm = C_s(nblyr+2)*Sintemp(nblyr+2)
      !-----------------------------------------------------------------

      fluxm_b = p5*(dhtmp*C_s(nblyr+2)* Sintemp(nblyr+2)*dbgrid &
                       +  C_s(nblyr+1)*(pre_sin(nblyr+1) &
                  + dhtmp * Sintemp(nblyr+1)*(c1-dbgrid)))

      fluxm_t = -p5 * (C_s(1) *  pre_sin(1)*dbgrid &
              +        C_s(2) * (pre_sin(2)*(c1-dbgrid) + dhtmp*Sintemp(2)))
           
      fluxm =  fluxm_b + fluxm_t
              
      Ssum_corr      = (-dh/hbri_old + p5*(dh/hbri_old)**2)*Ssum_old 
      Ssum_corr_flux = dh*Ssum_old/hin_next + Ssum_corr
      Ssum_corr      = Ssum_corr_flux
          
      fzsaln_g = fzsaln_g - hin_next *  fluxg_b &
                          * rhosi*p001/dts
       
      !approximate fluxes
      !fzsaln   = fzsaln   - hin_next * (fluxg &
      !                    + fluxb + fluxm + fluxcorr + Ssum_corr_flux) & 
      !                    * rhosi*p001/dts  

      fzsaln   = fzsaln + (Ssum_old*hin_old - Ssum_new*hin_next) &
                          * rhosi*p001/dts  ! positive into the ocean

      end subroutine calc_salt_fluxes

!=======================================================================
! 
! Test salt conservation:   flux conservative form d(hSin)/dt = -dF(x,Sin)/dx 
!  
      subroutine check_conserve_salt (mmax, mint,     dt,       dts,        &
                                      Ssum_old, Ssum_new, Ssum_corr,        & 
                                      fluxcorr, fluxb,    fluxg,     fluxm, &
                                      hbrin,    hbri_old, l_stop)

      use ice_colpkg_shared, only: rhosi

      integer(kind=int_kind), intent(in) :: &
         mint      , & ! current iteration
         mmax          ! maximum number of iterations

      real (kind=dbl_kind), intent(in) :: &
         dt, dts        , &  ! thermodynamic and halodynamic timesteps (s)
         hbrin          , &  ! (m) final brine height 
         hbri_old       , &  ! (m) initial brine height
         Ssum_old       , &  ! initial integrated salt content
         Ssum_new       , &  ! final integrated salt content
         fluxcorr       , &  ! flux correction to ensure S >= min_salin
         Ssum_corr      , &  ! boundary flux correction due to numerics
         fluxb          , &  ! total boundary salt flux into the ice (+ into ice)
         fluxg          , &  ! total gravity drainage salt flux into the ice (+ into ice)
         fluxm               ! 

      logical (kind=log_kind), intent(inout) :: &   
         l_stop              ! if false, conservation satisfied within error

     ! local variables

     real (kind=dbl_kind):: &
         diff2     , & !
         dsum_flux , & ! salt change in kg/m^2/s
         flux_tot  , & ! fluxg + fluxb
         order     , & !
         dh

     real (kind=dbl_kind), parameter :: &
         accuracy = 1.0e-7_dbl_kind ! g/kg/m^2/s difference between boundary fluxes 

     character(len=char_len_long) :: &
         warning ! warning message
     
         dh = (hbrin-hbri_old)/real(mmax,kind=dbl_kind)

         flux_tot = (fluxb + fluxg + fluxm + fluxcorr + Ssum_corr)*&
                    (hbri_old + (real(mint,kind=dbl_kind))*dh)/dt
         dsum_flux =(Ssum_new*(hbri_old + (real(mint,kind=dbl_kind))*dh) - &
                     Ssum_old*(hbri_old + (real(mint,kind=dbl_kind)-c1)* &
                     dh) )/dt
         order = abs(dh/min(hbri_old,hbrin))
         if (abs(dsum_flux) > accuracy) then
           diff2 = abs(dsum_flux - flux_tot)
           if (diff2 >  puny .AND. diff2 > order ) then 
              l_stop = .true.
              write(warning,*) 'Poor salt conservation: check_conserve_salt'
              call add_warning(warning)
              write(warning,*) 'mint:', mint
              call add_warning(warning)
              write(warning,*) 'Ssum_corr',Ssum_corr
              call add_warning(warning)
              write(warning,*) 'fluxb,fluxg,fluxm,flux_tot,fluxcorr:'
              call add_warning(warning)
              write(warning,*)  fluxb,fluxg,fluxm,flux_tot,fluxcorr
              call add_warning(warning)
              write(warning,*) 'fluxg,',fluxg
              call add_warning(warning)
              write(warning,*) 'dsum_flux,',dsum_flux
              call add_warning(warning)
              write(warning,*) 'Ssum_new,Ssum_old,hbri_old,dh:'
              call add_warning(warning)
              write(warning,*)  Ssum_new,Ssum_old,hbri_old,dh
              call add_warning(warning)
              write(warning,*) 'diff2,order,puny',diff2,order,puny
              call add_warning(warning)
           endif
         endif

     end subroutine check_conserve_salt

!=======================================================================
!
! Aggregate flux information from all ice thickness categories
!
      subroutine merge_zsal_fluxes(aicenS,               &
                                   zsal_totn,  zsal_tot, & 
                                   fzsal,      fzsaln,   &
                                   fzsal_g,    fzsaln_g)

      ! single category fluxes
      real (kind=dbl_kind), intent(in):: &          
          aicenS   , & ! concentration of ice
          fzsaln  , &  ! salt flux                       (kg/m**2/s)
          fzsaln_g     ! gravity drainage salt flux      (kg/m**2/s)

      real (kind=dbl_kind), intent(in):: &
          zsal_totn   ! tot salinity in category (psu*volume*rhosi)

      real (kind=dbl_kind), intent(inout):: &          
          zsal_tot, & ! tot salinity (psu*rhosi*total vol ice)
          fzsal   , & ! salt flux                       (kg/m**2/s)
          fzsal_g     ! gravity drainage salt flux      (kg/m**2/s)

      !-----------------------------------------------------------------
      ! Merge fluxes
      !-----------------------------------------------------------------

      zsal_tot  = zsal_tot + zsal_totn  ! already *aicenS

      ! ocean tot and gravity drainage salt fluxes
      fzsal    = fzsal   + fzsaln   * aicenS
      fzsal_g  = fzsal_g + fzsaln_g * aicenS

      end subroutine merge_zsal_fluxes

!==========================================================================
!
! For each grid cell, sum field over all ice layers.  "Net" refers to  the column
! integration while "avg"  is normalized by the ice thickness

      subroutine column_sum_zsal (zsal_totn, nblyr,   &
                                  vicenS, trcrn_S, fbri)

      use ice_colpkg_shared, only: rhosi

      integer (kind=int_kind), intent(in) :: &
         nblyr         ! number of layers

      real (kind=dbl_kind),  intent(in):: &          
         vicenS    , & ! volume of ice (m)
         fbri          ! brine height to ice thickness ratio

      real (kind=dbl_kind), dimension (nblyr), intent(in) :: &
         trcrn_S       ! input field

      real (kind=dbl_kind), intent(inout) :: &
         zsal_totn    ! avg salinity (psu*rhosi*vol of ice)

      ! local variables

      integer (kind=int_kind) :: &
           k           ! layer index

      do k = 1, nblyr
         zsal_totn = zsal_totn &
                   + rhosi * trcrn_S(k) &
                           * fbri  &
                           * vicenS/real(nblyr,kind=dbl_kind)
      enddo ! k

      end subroutine column_sum_zsal

!=======================================================================

      end module ice_zsalinity

!=======================================================================
