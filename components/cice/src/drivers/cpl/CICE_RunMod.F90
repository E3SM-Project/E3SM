!=======================================================================
!
!BOP
!
! !MODULE: CICE_RunMod
!
! !DESCRIPTION:
!
!  Contains CICE component driver routines particular to the specific
!  driver.
!
! !REVISION HISTORY:
!  SVN:$Id: CICE_RunMod.F90 138 2008-07-08 20:39:37Z eclare $
!
!  authors Elizabeth C. Hunke, LANL
!          Philip W. Jones, LANL
!          William H. Lipscomb, LANL
!
! 2008 ECH: created module by moving subroutines from drivers/cice4/
!
! !INTERFACE:
!
      module CICE_RunMod
!
! !USES:
!
      use ice_age
      use ice_aerosol
      use ice_atmo
      use ice_calendar
      use ice_communicate
      use ice_diagnostics
      use ice_domain
      use ice_dyn_evp
      use ice_fileunits
      use ice_flux
      use ice_forcing
      use ice_FY
      use ice_grid
      use ice_history
      use ice_restart
      use ice_itd
      use ice_kinds_mod
      use ice_mechred
      use ice_meltpond
      use ice_ocean
      use ice_orbital
      use ice_shortwave
      use ice_state
      use ice_therm_itd
      use ice_therm_vertical
      use ice_timers
      use ice_transport_driver
      use ice_transport_remap

      implicit none
      private
      save

! !PUBLIC MEMBER FUNCTIONS:

      public :: step_therm1, coupling_prep
!
!EOP
!
!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: step_therm1 - step pre-coupler thermodynamics
!
! !DESCRIPTION:
!
! Wrapper to driver for updating ice and snow internal temperatures and
! computing thermodynamic growth rates and coupler fluxes. Needed for 
! introducing OpenMP threading more simply into code.
!
! !REVISION HISTORY:
!
! authors: Mariana Vertenstein, NCAR
!
! !INTERFACE:

      subroutine step_therm1(dt)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step
!
!EOP
!
      integer (kind=int_kind) :: &
         iblk    ! block index

      call ice_timer_start(timer_column)  ! column physics
      call ice_timer_start(timer_thermo)  ! thermodynamics

      call init_history_therm    ! initialize thermo history variables
      call init_flux_ocn         ! initialize ocean fluxes sent to coupler

      if (oceanmixed_ice) &
           call ocean_mixed_layer (dt)   ! ocean surface fluxes and sst

!      call ice_timer_start(timer_tmp)  ! temporary timer

      !$OMP PARALLEL DO PRIVATE(iblk)
      do iblk = 1, nblocks
         call step_therm1_iblk(dt, iblk)
      end do
      !$OMP END PARALLEL DO

!      call ice_timer_stop(timer_tmp)  ! temporary timer
      call ice_timer_stop(timer_thermo) ! thermodynamics
      call ice_timer_stop(timer_column) ! column physics

      end subroutine step_therm1

!=======================================================================
!BOP
!
! !ROUTINE: step_therm1_iblk - step pre-coupler thermodynamics
!
! !DESCRIPTION:
!
! Driver for updating ice and snow internal temperatures and
! computing thermodynamic growth rates and coupler fluxes.
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!
! !INTERFACE:

      subroutine step_therm1_iblk (dt, iblk)
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      integer (kind=int_kind), intent(in) :: &
         iblk ! block index

!
!EOP
!
      integer (kind=int_kind) :: &
         i, j, ij    , & ! horizontal indices
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         n           , & ! thickness category index
         il1, il2    , & ! ice layer indices for eice
         sl1, sl2        ! snow layer indices for esno

      integer (kind=int_kind) :: &  
         icells          ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         indxi, indxj    ! indirect indices for cells with aicen > puny

! 2D coupler variables (computed for each category, then aggregated)

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         fsensn      , & ! surface downward sensible heat     (W/m^2)
         fswabsn     , & ! shortwave absorbed by ice          (W/m^2)
         flwoutn     , & ! upward LW at surface               (W/m^2)
         evapn       , & ! flux of vapor, atmos to ice   (kg m-2 s-1)
         freshn      , & ! flux of water, ice to ocean     (kg/m^2/s)
         fsaltn      , & ! flux of salt, ice to ocean      (kg/m^2/s)
         fhocnn      , & ! fbot corrected for leftover energy (W/m^2)
         strairxn    , & ! air/ice zonal  strss,              (N/m^2)
         strairyn    , & ! air/ice merdnl strss,              (N/m^2)
         Urefn       , & ! wind speed reference level         (m/s)
         Trefn       , & ! air tmp reference level            (K)
         Qrefn           ! air sp hum reference level         (kg/kg)

      ! other local variables
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         Tbot        , & ! ice bottom surface temperature (deg C)
         fbot        , & ! ice-ocean heat flux at bottom surface (W/m^2)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef          ! transfer coefficient for latent heat

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         melttn      , & ! top melt in category n (m)
         meltbn      , & ! bottom melt in category n (m)
         meltsn      , & ! snow melt in category n (m)
         congeln     , & ! congelation ice formation in category n (m)
         snoicen         ! snow-ice formation in category n (m)

      ! Local variables to keep track of mass budget for aerosols MH
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         vsno_old   

      type (block) :: &
         this_block      ! block information for current block

      logical (kind=log_kind) :: &
         l_stop          ! if true, abort the model

      integer (kind=int_kind) :: &
         istop, jstop    ! indices of grid cell where model aborts 

      ! local (single block)
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         worka, workb

      l_stop = .false.

      worka(:,:) = c0
      workb(:,:) = c0

      this_block = get_block(blocks_ice(iblk),iblk)         
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi

      !-----------------------------------------------------------------
      ! Save the ice area passed to the coupler (so that history fields
      !  can be made consistent with coupler fields).
      ! Save the initial ice area and volume in each category.
      !-----------------------------------------------------------------

      do j = 1, ny_block
      do i = 1, nx_block
         aice_init (i,j,  iblk) = aice (i,j,  iblk)
      enddo
      enddo

      do n = 1, ncat
      do j = 1, ny_block
      do i = 1, nx_block
         aicen_init(i,j,n,iblk) = aicen(i,j,n,iblk)
         vicen_init(i,j,n,iblk) = vicen(i,j,n,iblk)
      enddo
      enddo
      enddo

      if (mod(time-dt_thm,dt_dyn) == c0) then
         strairxT_accum(:,:,iblk) = c0
         strairyT_accum(:,:,iblk) = c0
      endif

      !-----------------------------------------------------------------
      ! Adjust frzmlt to account for ice-ocean heat fluxes since last
      !  call to coupler.
      ! Compute lateral and bottom heat fluxes.
      !-----------------------------------------------------------------

      call frzmlt_bottom_lateral                                      &
                             (nx_block,           ny_block,           &
                              ilo, ihi,           jlo, jhi,           &
                              dt,                                     &
                              aice  (:,:,  iblk), frzmlt(:,:,  iblk), &
                              eicen (:,:,:,iblk), esnon (:,:,:,iblk), &
                              sst   (:,:,  iblk), Tf    (:,:,  iblk), &
                              strocnxT(:,:,iblk), strocnyT(:,:,iblk), &
                              Tbot,               fbot,               &
                              rside (:,:,  iblk) )


      do n = 1, ncat

      !-----------------------------------------------------------------
      ! Identify cells with nonzero ice area
      !-----------------------------------------------------------------
           
         icells = 0
         do j = jlo, jhi
         do i = ilo, ihi
            if (aicen(i,j,n,iblk) > puny) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif
         enddo               ! i
         enddo               ! j

         if (calc_Tsfc .or. calc_strair .and. icells > 0) then

      !-----------------------------------------------------------------
      ! Atmosphere boundary layer calculation; compute coefficients
      ! for sensible and latent heat fluxes.
      !
      ! NOTE: The wind stress is computed here for later use if 
      !       calc_strair = .true.   Otherwise, the wind stress
      !       components are set to the data values.
      !-----------------------------------------------------------------

            if (trim(atmbndy) == 'constant') then
               call atmo_boundary_const(nx_block,      ny_block,        &
                                        'ice',          icells,         &
                                        indxi,          indxj,          &
                                        uatm(:,:,iblk), vatm(:,:,iblk), &
                                        wind(:,:,iblk), rhoa(:,:,iblk), &
                                        strairxn,       strairyn,       &
                                        lhcoef,         shcoef)
            else ! default
               call atmo_boundary_layer(nx_block,       ny_block,       &
                                        'ice',          icells,         &
                                        indxi,          indxj,          &
                                        trcrn(:,:,nt_Tsfc,n,iblk),      &
                                        potT(:,:,iblk),                 &
                                        uatm(:,:,iblk), vatm(:,:,iblk), &
                                        uvel(:,:,iblk), vvel(:,:,iblk), &
                                        wind(:,:,iblk), zlvl(:,:,iblk), &
                                        Qa  (:,:,iblk), rhoa(:,:,iblk), &
                                        strairxn,       strairyn,       &
                                        Urefn,                          &
                                        Trefn,          Qrefn,          &
                                        worka,          workb,          &
                                        lhcoef,         shcoef)
            endif ! atmbndy

         else

            ! Initialize for safety
            Urefn (:,:)  = c0
            Trefn (:,:)  = c0
            Qrefn (:,:)  = c0
            lhcoef(:,:)  = c0
            shcoef(:,:)  = c0

         endif   ! calc_Tsfc or calc_strair

         if (.not.(calc_strair)) then
            strairxn(:,:) = strax(:,:,iblk)
            strairyn(:,:) = stray(:,:,iblk)
         endif

      !-----------------------------------------------------------------
      ! Update ice age
      ! This is further adjusted for freezing in the thermodynamics.
      ! Melting does not alter the ice age.
      !-----------------------------------------------------------------

         if (tr_iage) then
            call increment_age (nx_block, ny_block,      &
                                dt, icells,              &
                                indxi, indxj,            &
                                trcrn(:,:,nt_iage,n,iblk))
         endif
         if (tr_FY) then
            call update_FYarea (nx_block, ny_block,      &
                                dt, icells,              &
                                indxi, indxj,            &
                                lmask_n(:,:,iblk),       &
                                lmask_s(:,:,iblk),       &
                                trcrn(:,:,nt_FY,n,iblk))
         endif

      !-----------------------------------------------------------------
      ! Vertical thermodynamics: Heat conduction, growth and melting.
      !----------------------------------------------------------------- 

         il1 = ilyr1(n)
         il2 = ilyrn(n)
         sl1 = slyr1(n)
         sl2 = slyrn(n)

         vsno_old   = vsnon(:,:,n,iblk) !MH

         call thermo_vertical                                       &
                         (nx_block,            ny_block,            &
                          dt,                  icells,              &
                          indxi,               indxj,               &
                          aicen(:,:,n,iblk),                        &
                          trcrn(:,:,:,n,iblk),                      &
                          vicen(:,:,n,iblk),   vsnon(:,:,n,iblk),   &
                          eicen  (:,:,il1:il2,iblk),                &
                          esnon  (:,:,sl1:sl2,iblk),                &
                          flw    (:,:,iblk),   potT (:,:,iblk),     &
                          Qa     (:,:,iblk),   rhoa (:,:,iblk),     &
                          fsnow  (:,:,iblk),                        &
                          fbot,                Tbot,                &
                          lhcoef,              shcoef,              &
                          fswsfcn(:,:,n,iblk), fswintn(:,:,n,iblk), &
                          fswthrun(:,:,n,iblk),                     &
                          Sswabsn(:,:,sl1:sl2,iblk),                &
                          Iswabsn(:,:,il1:il2,iblk),                &
                          fsurfn(:,:,n,iblk),  fcondtopn(:,:,n,iblk),&
                          fsensn,              flatn(:,:,n,iblk),   &
                          fswabsn,             flwoutn,             &
                          evapn,               freshn,              &
                          fsaltn,              fhocnn,              &
                          melttn,              meltsn,              &
                          meltbn,                                   &
                          congeln,             snoicen,             &
                          mlt_onset(:,:,iblk), frz_onset(:,:,iblk), &
                          yday,                l_stop,              &
                          istop,               jstop)

         if (l_stop) then
            write (nu_diag,*) 'istep1, my_task, iblk =', &
                               istep1, my_task, iblk
            write (nu_diag,*) 'Global block:', this_block%block_id
            if (istop > 0 .and. jstop > 0) &
                 write(nu_diag,*) 'Global i and j:', &
                                  this_block%i_glob(istop), &
                                  this_block%j_glob(jstop) 
                 write(nu_diag,*) 'Lat, Lon:', &
                                 TLAT(istop,jstop,iblk)*rad_to_deg, &
                                 TLON(istop,jstop,iblk)*rad_to_deg
            call abort_ice ('ice: Vertical thermo error')
         endif

      !-----------------------------------------------------------------
      ! Aerosol update  MH
      !-----------------------------------------------------------------
         if (tr_aero) then

            if (icells > 0) then

               call update_aerosol (nx_block, ny_block,                  &
                                    dt, icells,                          &
                                    indxi, indxj,                        &
                                    melttn, meltsn,                      &
                                    meltbn, congeln, snoicen,            &
                                    fsnow(:,:,iblk),                     &
                                    trcrn(:,:,:,n,iblk),                 &
                                    aicen_init(:,:,n,iblk),              &
                                    vicen_init(:,:,n,iblk),              &
                                    vsno_old,                            &
                                    vicen(:,:,n,iblk),                   &
                                    vsnon(:,:,n,iblk),                   &
                                    aicen(:,:,n,iblk),faero(:,:,:,iblk), &
                                    fsoot(:,:,:,iblk))

            endif

         endif

      !-----------------------------------------------------------------
      ! Melt ponds
      !-----------------------------------------------------------------

         if (tr_pond) then

            call compute_ponds(nx_block, ny_block,                      &
                               ilo, ihi, jlo, jhi,                      &
                               melttn, meltsn, frain(:,:,iblk),         &
                               aicen (:,:,n,iblk), vicen (:,:,n,iblk),  &
                               vsnon (:,:,n,iblk), trcrn (:,:,:,n,iblk),&
                               apondn(:,:,n,iblk), hpondn(:,:,n,iblk))

         endif

      !-----------------------------------------------------------------
      ! Increment area-weighted fluxes.
      !-----------------------------------------------------------------

         call merge_fluxes (nx_block,           ny_block,             &
                            icells,                                   &
                            indxi,              indxj,                &
                            aicen_init(:,:,n,iblk),                   &
                            flw(:,:,iblk),      coszen(:,:,iblk),     &
                            strairxn,           strairyn,             &
                            fsurfn(:,:,n,iblk), fcondtopn(:,:,n,iblk),&
                            fsensn,             flatn(:,:,n,iblk),    &
                            fswabsn,            flwoutn,              &
#if (defined AEROFRC) || (defined PONDFRC) || (defined CCSM3FRC)
                            fswsfcn(:,:,n,iblk),fswintn(:,:,n,iblk),  &
#endif
#ifdef AEROFRC
                            dfswabsn_noaero(:,:,n,iblk),               &
                            dfswsfcn_noaero(:,:,n,iblk),               &
                            dfswintn_noaero(:,:,n,iblk),               &
                            dfswthrun_noaero(:,:,n,iblk),              &
#endif
#ifdef CCSM3FRC
                            dfswabsn_ccsm3(:,:,n,iblk),               &
                            dfswsfcn_ccsm3(:,:,n,iblk),               &
                            dfswintn_ccsm3(:,:,n,iblk),               &
                            dfswthrun_ccsm3(:,:,n,iblk),              &
#endif
#ifdef PONDFRC
                            dfswabsn_nopond(:,:,n,iblk),               &
                            dfswsfcn_nopond(:,:,n,iblk),               &
                            dfswintn_nopond(:,:,n,iblk),               &
                            dfswthrun_nopond(:,:,n,iblk),              &
#endif
                            evapn,                                    &
                            Urefn,                                    &
                            Trefn,              Qrefn,                &
                            freshn,             fsaltn,               &
                            fhocnn,             fswthrun(:,:,n,iblk), &
                            strairxT(:,:,iblk), strairyT  (:,:,iblk), &
                            fsurf   (:,:,iblk), fcondtop  (:,:,iblk), &
                            fsens   (:,:,iblk), flat      (:,:,iblk), &
                            fswabs  (:,:,iblk), flwout    (:,:,iblk), &
#if (defined AEROFRC) || (defined PONDFRC) || (defined CCSM3FRC)
                            fswsfc  (:,:,iblk), fswint    (:,:,iblk), &
#endif
#ifdef AEROFRC
                            dfswabs_noaero(:,:,iblk),dfswsfc_noaero(:,:,iblk), &
                            dfswint_noaero(:,:,iblk),dfswthru_noaero(:,:,iblk),&
#endif
#ifdef CCSM3FRC
                            dfswabs_ccsm3(:,:,iblk),dfswsfc_ccsm3(:,:,iblk), &
                            dfswint_ccsm3(:,:,iblk),dfswthru_ccsm3(:,:,iblk),&
#endif
#ifdef PONDFRC
                            dfswabs_nopond(:,:,iblk),dfswsfc_nopond(:,:,iblk), &
                            dfswint_nopond(:,:,iblk),dfswthru_nopond(:,:,iblk),&
#endif
                            evap    (:,:,iblk),                       &
                            Uref    (:,:,iblk),                       &
                            Tref    (:,:,iblk), Qref      (:,:,iblk), &
                            fresh   (:,:,iblk), fsalt   (:,:,iblk),   &
                            fhocn   (:,:,iblk), fswthru (:,:,iblk),   &
                            melttn, meltsn, meltbn, congeln, snoicen, &
                            meltt   (:,:,iblk),  melts   (:,:,iblk),  &
                            meltb   (:,:,iblk),                       &
                            congel  (:,:,iblk),  snoice  (:,:,iblk))

      enddo                  ! ncat


! Accumulate stresses when super-cycling the dynamics. Otherwise just
! use the stresses as computed.

      if (dt_thm < dt_dyn) then
         strairxT_accum(:,:,iblk) = strairxT_accum(:,:,iblk) &
                                  + strairxT(:,:,iblk) * dt_thm / dt_dyn
         strairyT_accum(:,:,iblk) = strairyT_accum(:,:,iblk) &
                                  + strairyT(:,:,iblk) * dt_thm / dt_dyn
      else
         strairxT_accum(:,:,iblk) = strairxT(:,:,iblk)
         strairyT_accum(:,:,iblk) = strairyT(:,:,iblk)
      endif

      !-----------------------------------------------------------------
      ! Update mixed layer with heat and radiation from ice.
      !-----------------------------------------------------------------

      if (oceanmixed_ice) then
         do j = jlo, jhi
         do i = ilo, ihi
            if (hmix(i,j,iblk) > puny) then
               sst(i,j,iblk) = sst(i,j,iblk) &
                    + (fhocn(i,j,iblk) + fswthru(i,j,iblk))*dt &
                    / (cprho*hmix(i,j,iblk))
            endif
         enddo
         enddo
      endif

      end subroutine step_therm1_iblk
      
!=======================================================================
!BOP
!
! !ROUTINE: coupling_prep
!
! !DESCRIPTION:
!
! Prepare for coupling
!
! !REVISION HISTORY:
!
! authors: Elizabeth C. Hunke, LANL
!
! !INTERFACE:

      subroutine coupling_prep
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      type (block) :: &
         this_block      ! block information for current block

      integer (kind=int_kind) :: & 
         iblk        , & ! block index 
         n           , & ! thickness category index
         i,j         , & ! horizontal indices
         ilo,ihi,jlo,jhi ! beginning and end of physical domain

      real (kind=dbl_kind) :: cszn ! counter for history averaging

      call ice_timer_start(timer_column)

      !-----------------------------------------------------------------
      ! Update mixed layer with heat and radiation from ice.
      !-----------------------------------------------------------------

      if (oceanmixed_ice) &
         call ocean_mixed_layer (dt)   ! ocean surface fluxes and sst

      do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! Aggregate albedos
      !-----------------------------------------------------------------

         do j = 1, ny_block
         do i = 1, nx_block
            alvdf(i,j,iblk) = c0
            alidf(i,j,iblk) = c0
            alvdr(i,j,iblk) = c0
            alidr(i,j,iblk) = c0

            albice(i,j,iblk) = c0
            albsno(i,j,iblk) = c0
            albpnd(i,j,iblk) = c0
#ifdef AEROFRC
            dalvdf_noaero(i,j,iblk) = c0
            dalidf_noaero(i,j,iblk) = c0
            dalvdr_noaero(i,j,iblk) = c0
            dalidr_noaero(i,j,iblk) = c0

            dalbice_noaero(i,j,iblk) = c0
            dalbsno_noaero(i,j,iblk) = c0
            dalbpnd_noaero(i,j,iblk) = c0
#endif
#ifdef CCSM3FRC
            dalvdf_ccsm3(i,j,iblk) = c0
            dalidf_ccsm3(i,j,iblk) = c0
            dalvdr_ccsm3(i,j,iblk) = c0
            dalidr_ccsm3(i,j,iblk) = c0

            dalbice_ccsm3(i,j,iblk) = c0
            dalbsno_ccsm3(i,j,iblk) = c0
#endif
#ifdef PONDFRC
            dalvdf_nopond(i,j,iblk) = c0
            dalidf_nopond(i,j,iblk) = c0
            dalvdr_nopond(i,j,iblk) = c0
            dalidr_nopond(i,j,iblk) = c0

            dalbice_nopond(i,j,iblk) = c0
            dalbsno_nopond(i,j,iblk) = c0
            dalbpnd_nopond(i,j,iblk) = c0
#endif

            ! for history averaging
            cszn = c0
            if (coszen(i,j,iblk) > puny) cszn = c1
            do n = 1, nstreams
               albcnt(i,j,iblk,n) = albcnt(i,j,iblk,n) + cszn
            enddo
         enddo
         enddo
         do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            alvdf(i,j,iblk) = alvdf(i,j,iblk) &
               + alvdfn(i,j,n,iblk)*aicen(i,j,n,iblk)
            alidf(i,j,iblk) = alidf(i,j,iblk) &
               + alidfn(i,j,n,iblk)*aicen(i,j,n,iblk)
            alvdr(i,j,iblk) = alvdr(i,j,iblk) &
               + alvdrn(i,j,n,iblk)*aicen(i,j,n,iblk)
            alidr(i,j,iblk) = alidr(i,j,iblk) &
               + alidrn(i,j,n,iblk)*aicen(i,j,n,iblk)

            if (coszen(i,j,iblk) > puny) then ! sun above horizon
               albice(i,j,iblk) = albice(i,j,iblk) &
                  + albicen(i,j,n,iblk)*aicen(i,j,n,iblk)
               albsno(i,j,iblk) = albsno(i,j,iblk) &
                  + albsnon(i,j,n,iblk)*aicen(i,j,n,iblk)
               albpnd(i,j,iblk) = albpnd(i,j,iblk) &
                  + albpndn(i,j,n,iblk)*aicen(i,j,n,iblk)
            endif
#ifdef AEROFRC
            dalvdf_noaero(i,j,iblk) = dalvdf_noaero(i,j,iblk) &
               + dalvdfn_noaero(i,j,n,iblk)*aicen(i,j,n,iblk)
            dalidf_noaero(i,j,iblk) = dalidf_noaero(i,j,iblk) &
               + dalidfn_noaero(i,j,n,iblk)*aicen(i,j,n,iblk)
            dalvdr_noaero(i,j,iblk) = dalvdr_noaero(i,j,iblk) &
               + dalvdrn_noaero(i,j,n,iblk)*aicen(i,j,n,iblk)
            dalidr_noaero(i,j,iblk) = dalidr_noaero(i,j,iblk) &
               + dalidrn_noaero(i,j,n,iblk)*aicen(i,j,n,iblk)

            if (coszen(i,j,iblk) > puny) then ! sun above horizon
               dalbice_noaero(i,j,iblk) = dalbice_noaero(i,j,iblk) &
                  + dalbicen_noaero(i,j,n,iblk)*aicen(i,j,n,iblk)
               dalbsno_noaero(i,j,iblk) = dalbsno_noaero(i,j,iblk) &
                  + dalbsnon_noaero(i,j,n,iblk)*aicen(i,j,n,iblk)
               dalbpnd_noaero(i,j,iblk) = dalbpnd_noaero(i,j,iblk) &
                  + dalbpndn_noaero(i,j,n,iblk)*aicen(i,j,n,iblk)
            endif
#endif
#ifdef CCSM3FRC
            dalvdf_ccsm3(i,j,iblk) = dalvdf_ccsm3(i,j,iblk) &
               + dalvdfn_ccsm3(i,j,n,iblk)*aicen(i,j,n,iblk)
            dalidf_ccsm3(i,j,iblk) = dalidf_ccsm3(i,j,iblk) &
               + dalidfn_ccsm3(i,j,n,iblk)*aicen(i,j,n,iblk)
            dalvdr_ccsm3(i,j,iblk) = dalvdr_ccsm3(i,j,iblk) &
               + dalvdrn_ccsm3(i,j,n,iblk)*aicen(i,j,n,iblk)
            dalidr_ccsm3(i,j,iblk) = dalidr_ccsm3(i,j,iblk) &
               + dalidrn_ccsm3(i,j,n,iblk)*aicen(i,j,n,iblk)

            if (coszen(i,j,iblk) > puny) then ! sun above horizon
               dalbice_ccsm3(i,j,iblk) = dalbice_ccsm3(i,j,iblk) &
                  + dalbicen_ccsm3(i,j,n,iblk)*aicen(i,j,n,iblk)
               dalbsno_ccsm3(i,j,iblk) = dalbsno_ccsm3(i,j,iblk) &
                  + dalbsnon_ccsm3(i,j,n,iblk)*aicen(i,j,n,iblk)
            endif
#endif
#ifdef PONDFRC
            dalvdf_nopond(i,j,iblk) = dalvdf_nopond(i,j,iblk) &
               + dalvdfn_nopond(i,j,n,iblk)*aicen(i,j,n,iblk)
            dalidf_nopond(i,j,iblk) = dalidf_nopond(i,j,iblk) &
               + dalidfn_nopond(i,j,n,iblk)*aicen(i,j,n,iblk)
            dalvdr_nopond(i,j,iblk) = dalvdr_nopond(i,j,iblk) &
               + dalvdrn_nopond(i,j,n,iblk)*aicen(i,j,n,iblk)
            dalidr_nopond(i,j,iblk) = dalidr_nopond(i,j,iblk) &
               + dalidrn_nopond(i,j,n,iblk)*aicen(i,j,n,iblk)

            if (coszen(i,j,iblk) > puny) then ! sun above horizon
               dalbice_nopond(i,j,iblk) = dalbice_nopond(i,j,iblk) &
                  + dalbicen_nopond(i,j,n,iblk)*aicen(i,j,n,iblk)
               dalbsno_nopond(i,j,iblk) = dalbsno_nopond(i,j,iblk) &
                  + dalbsnon_nopond(i,j,n,iblk)*aicen(i,j,n,iblk)
               dalbpnd_nopond(i,j,iblk) = dalbpnd_nopond(i,j,iblk) &
                  + dalbpndn_nopond(i,j,n,iblk)*aicen(i,j,n,iblk)
            endif
#endif
         enddo
         enddo
         enddo

      !----------------------------------------------------------------
      ! Store grid box mean albedos and fluxes before scaling by aice
      !----------------------------------------------------------------

         do j = 1, ny_block
         do i = 1, nx_block
            alvdf_gbm  (i,j,iblk) = alvdf  (i,j,iblk)
            alidf_gbm  (i,j,iblk) = alidf  (i,j,iblk)
            alvdr_gbm  (i,j,iblk) = alvdr  (i,j,iblk)
            alidr_gbm  (i,j,iblk) = alidr  (i,j,iblk)
            fresh_gbm  (i,j,iblk) = fresh  (i,j,iblk)
            fsalt_gbm  (i,j,iblk) = fsalt  (i,j,iblk)
            fhocn_gbm  (i,j,iblk) = fhocn  (i,j,iblk)
            fswthru_gbm(i,j,iblk) = fswthru(i,j,iblk)

!echmod
      !-----------------------------------------------------------------
      ! Save net shortwave for scaling factor in scale_factor
      !-----------------------------------------------------------------
            if (istep > 0) then

            scale_factor(i,j,iblk) = &
                       swvdr(i,j,iblk)*(c1 - alvdr_gbm(i,j,iblk)) &
                     + swvdf(i,j,iblk)*(c1 - alvdf_gbm(i,j,iblk)) &
                     + swidr(i,j,iblk)*(c1 - alidr_gbm(i,j,iblk)) &
                     + swidf(i,j,iblk)*(c1 - alidf_gbm(i,j,iblk))

            endif
!echmod

         enddo
         enddo

         call scale_fluxes (nx_block,            ny_block,           &
                            tmask    (:,:,iblk),                     &
                            aice     (:,:,iblk), Tf      (:,:,iblk), &
                            Tair     (:,:,iblk), Qa      (:,:,iblk), &
                            wind     (:,:,iblk),                     &
                            strairxT (:,:,iblk), strairyT(:,:,iblk), &
                            fsens    (:,:,iblk), flat    (:,:,iblk), &
                            fswabs   (:,:,iblk), flwout  (:,:,iblk), &
                            evap     (:,:,iblk),                     &
                            Uref     (:,:,iblk),                     &
                            Tref     (:,:,iblk), Qref    (:,:,iblk), &
                            fresh    (:,:,iblk), fsalt   (:,:,iblk), &
                            fhocn    (:,:,iblk), fswthru (:,:,iblk), &
                            fsoot    (:,:,:,iblk),                   &
                            alvdr    (:,:,iblk), alidr   (:,:,iblk), &
                            alvdf    (:,:,iblk), alidf   (:,:,iblk))

      enddo

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (scale_factor,     halo_info, &
                           field_loc_center, field_type_scalar)
      call ice_timer_stop(timer_bound)

      call ice_timer_stop(timer_column)

      end subroutine coupling_prep

!=======================================================================

      end module CICE_RunMod

!=======================================================================
