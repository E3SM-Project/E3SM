!=======================================================================
!BOP
!
! !MODULE: ice_shortwave - reflected, transmitted, and absorbed solar rad
!
! !DESCRIPTION:
!
! The albedo and absorbed/transmitted flux parameterizations for
! snow over ice, bare ice and ponded ice.
!
! Presently, two methods are included:
!   (1) CCSM3 
!   (2) Delta-Eddington 
! as two distinct routines.
! Either can be called from the ice driver.
!
! The Delta-Eddington method is described here:
!
! Briegleb, B. P., and B. Light (2007): A Delta-Eddington Multiple 
!    Scattering Parameterization for Solar Radiation in the Sea Ice 
!    Component of the Community Climate System Model, NCAR Technical 
!    Note  NCAR/TN-472+STR  February 2007
!
! !REVISION HISTORY:
!  SVN:$Id: ice_shortwave.F90 144 2008-08-12 21:37:19Z eclare $
!
! name: originally ice_albedo
!
! authors:  Bruce P. Briegleb, NCAR
!           Elizabeth C. Hunke and William H. Lipscomb, LANL
! 2005, WHL: Moved absorbed_solar from ice_therm_vertical to this 
!            module and changed name from ice_albedo
! 2006, WHL: Added Delta Eddington routines from Bruce Briegleb
! 2006, ECH: Changed data statements in Delta Eddington routines (no 
!            longer hardwired)
!            Converted to free source form (F90)
! 2007, BPB: Completely updated Delta-Eddington code, so that:
!            (1) multiple snow layers enabled (i.e. nslyr > 1)
!            (2) included SSL for snow surface absorption
!            (3) added Sswabs for internal snow layer absorption
!            (4) variable sea ice layers allowed (i.e. not hardwired)
!            (5) updated all inherent optical properties
!            (6) included algae absorption for sea ice lowest layer
!            (7) very complete internal documentation included
! 2007, ECH: Improved efficiency
! 2008, BPB: Added aerosols to Delta Eddington code
!
! !INTERFACE:
!
      module ice_shortwave
!
! !USES:
!
      use ice_kinds_mod
      use ice_domain_size
      use ice_constants
      use ice_blocks
      use ice_diagnostics
      use ice_communicate, only: my_task
      use perf_mod, only: t_startf, t_stopf, t_barrierf
!
!EOP
!
      implicit none
      save

      character (len=char_len) :: &
         shortwave, & ! shortwave method, 'default' ('ccsm3') or 'dEdd'
         albedo_type  ! albedo parameterization, 'default' ('ccsm3') or 'constant'
                      ! shortwave='dEdd' overrides this parameter

      ! baseline albedos for ccsm3 shortwave, set in namelist
      real (kind=dbl_kind) :: &
         albicev , & ! visible ice albedo for h > ahmax
         albicei , & ! near-ir ice albedo for h > ahmax
         albsnowv, & ! cold snow albedo, visible
         albsnowi    ! cold snow albedo, near IR

      ! category albedos
      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,ncat,max_blocks) :: &
#ifdef AEROFRC
         dalvdrn_noaero, & ! visible direct albedo (diag)    (fraction)
         dalidrn_noaero, & ! near-ir direct albedo (diag)    (fraction)
         dalvdfn_noaero, & ! visible diffuse albedo (diag)   (fraction)
         dalidfn_noaero, & ! near-ir diffuse albedo (diag)   (fraction)
#endif
#ifdef CCSM3FRC
         dalvdrn_ccsm3, & ! visible direct albedo (diag)    (fraction)
         dalidrn_ccsm3, & ! near-ir direct albedo (diag)    (fraction)
         dalvdfn_ccsm3, & ! visible diffuse albedo (diag)   (fraction)
         dalidfn_ccsm3, & ! near-ir diffuse albedo (diag)   (fraction)
#endif
#ifdef PONDFRC
         dalvdrn_nopond, & ! visible direct albedo (diag)    (fraction)
         dalidrn_nopond, & ! near-ir direct albedo (diag)    (fraction)
         dalvdfn_nopond, & ! visible diffuse albedo (diag)   (fraction)
         dalidfn_nopond, & ! near-ir diffuse albedo (diag)   (fraction)
#endif
         alvdrn       , & ! visible direct albedo           (fraction)
         alidrn       , & ! near-ir direct albedo           (fraction)
         alvdfn       , & ! visible diffuse albedo          (fraction)
         alidfn           ! near-ir diffuse albedo          (fraction)

      ! albedo components for history
      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,ncat,max_blocks) :: &
#ifdef AEROFRC
         dalbicen_noaero, & ! bare ice (diag)
         dalbsnon_noaero, & ! snow (diag)
         dalbpndn_noaero, & ! pond (diag)
#endif
#ifdef CCSM3FRC
         dalbicen_ccsm3, & ! bare ice (diag)
         dalbsnon_ccsm3, & ! snow (diag)
#endif
#ifdef PONDFRC
         dalbicen_nopond, & ! bare ice (diag)
         dalbsnon_nopond, & ! snow (diag)
         dalbpndn_nopond, & ! pond (diag)
#endif
         albicen  , & ! bare ice 
         albsnon  , & ! snow 
         albpndn      ! pond 

      ! shortwave components
      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,ntilyr,max_blocks) :: &
#ifdef AEROFRC
         dIswabsn_noaero, & ! SW radiation absorbed in ice layers (diag) (W m-2)
#endif
#ifdef CCSM3FRC
         dIswabsn_ccsm3, & ! SW radiation absorbed in ice layers (diag) (W m-2)
#endif
#ifdef PONDFRC
         dIswabsn_nopond, & ! SW radiation absorbed in ice layers (diag) (W m-2)
#endif
         Iswabsn           ! SW radiation absorbed in ice layers (W m-2)

      real (kind=dbl_kind), &
         dimension (nx_block,ny_block,ntslyr,max_blocks) :: &
#ifdef AEROFRC
         dSswabsn_noaero, & ! SW radiation absorbed in snow layers (diag) (W m-2)
#endif
#ifdef CCSM3FRC
         dSswabsn_ccsm3, & ! SW radiation absorbed in snow layers (diag) (W m-2)
#endif
#ifdef PONDFRC
         dSswabsn_nopond, & ! SW radiation absorbed in snow layers (diag) (W m-2)
#endif
         Sswabsn           ! SW radiation absorbed in snow layers (W m-2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat,max_blocks) :: &
#ifdef AEROFRC
         dfswabsn_noaero , & ! SW absorbed in ice/snow (diag) (W m-2)
         dfswsfcn_noaero , & ! SW absorbed at ice/snow surface (diag) (W m-2)
         dfswthrun_noaero, & ! SW through ice to ocean (diag)     (W/m^2)
         dfswintn_noaero , & ! SW absorbed in ice interior, below surface (W m-2)
#endif
#ifdef CCSM3FRC
         dfswabsn_ccsm3 , & ! SW absorbed in ice/snow (diag) (W m-2)
         dfswsfcn_ccsm3 , & ! SW absorbed at ice/snow surface (diag) (W m-2)
         dfswthrun_ccsm3, & ! SW through ice to ocean (diag)     (W/m^2)
         dfswintn_ccsm3 , & ! SW absorbed in ice interior, below surface (W m-2)
#endif
#ifdef PONDFRC
         dfswabsn_nopond , & ! SW absorbed in ice/snow (diag) (W m-2)
         dfswsfcn_nopond , & ! SW absorbed at ice/snow surface (diag) (W m-2)
         dfswthrun_nopond, & ! SW through ice to ocean (diag)     (W/m^2)
         dfswintn_nopond , & ! SW absorbed in ice interior, below surface (W m-2)
#endif
         fswsfcn        , & ! SW absorbed at ice/snow surface (W m-2)
         fswthrun       , & ! SW through ice to ocean            (W/m^2)
         fswintn            ! SW absorbed in ice interior, below surface (W m-2)

      real (kind=dbl_kind) :: &
         rnilyr      , & ! real(nilyr)
         rnslyr          ! real(nslyr)

      ! melt pond tuning parameters, set in namelist
      real (kind=dbl_kind) :: &
         R_ice , & ! sea ice tuning parameter; +1 > 1sig increase in albedo
         R_pnd , & ! ponded ice tuning parameter; +1 > 1sig increase in albedo
         R_snw     ! snow tuning parameter; +1 > ~.01 change in broadband albedo

      real (kind=dbl_kind) :: &
         dT_mlt_in          , &  ! temperature at which melt begins (tuning)
         rsnw_melt_in            ! maximum snow grain radius (tuning)

      ! for delta Eddington
      real (kind=dbl_kind) :: &
         exp_min              ! minimum exponential value
      logical :: tflag = .false.   ! timer flag
      character(len=1) :: tstr    ! timer string
#ifdef MODAL_AER
      !mgf++
      real (kind=dbl_kind), dimension (3,10) :: kaer_bc_tab ! BC mass extinction cross section (m2/kg)
      real (kind=dbl_kind), dimension (3,10) :: waer_bc_tab ! BC single scatter albedo (fraction)
      real (kind=dbl_kind), dimension (3,10) :: gaer_bc_tab ! BC aerosol asymmetry parameter (cos(theta))
      real (kind=dbl_kind), dimension (3,10,8) :: bcenh
      !mgf--
#endif

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: init_shortwave
!
! !DESCRIPTION:
!
!  Initialize shortwave
! 
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine init_shortwave
!
! !USES:
!
      use ice_domain_size
      use ice_blocks
      use ice_calendar
      use ice_domain
      use ice_flux
      use ice_grid
      use ice_itd
      use ice_orbital
      use ice_state
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
!     local temporary variables

      integer (kind=int_kind) :: &
         icells          ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         indxi, indxj    ! indirect indices for cells with aicen > puny

      integer (kind=int_kind) :: &
         i, j, ij    , & ! horizontal indices
         iblk        , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         n           , & ! thickness category index
         il1, il2    , & ! ice layer indices for eice
         sl1, sl2        ! snow layer indices for esno

      real (kind=dbl_kind) :: cszn ! counter for history averaging

      type (block) :: &
         this_block      ! block information for current block

      do iblk=1,nblocks
      do j = 1, ny_block
      do i = 1, nx_block
         alvdr_gbm(i,j,iblk) = c0
         alidr_gbm(i,j,iblk) = c0
         alvdf_gbm(i,j,iblk) = c0
         alidf_gbm(i,j,iblk) = c0
      enddo
      enddo
      enddo

      if (trim(shortwave) == 'dEdd') then ! delta Eddington

#ifndef CCSMCOUPLED
         ! These come from the driver in the coupled model.
         call init_orbit       ! initialize orbital parameters
#endif
         call init_dEdd        ! initialize delta Eddington
 
      else                     ! basic (ccsm3) shortwave

         coszen(:,:,:) = p5 ! sun above the horizon

         do iblk=1,nblocks
            this_block = get_block(blocks_ice(iblk),iblk)         
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            do n = 1, ncat

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

               il1 = ilyr1(n)
               il2 = ilyrn(n)
               sl1 = slyr1(n)
               sl2 = slyrn(n)

               Sswabsn(:,:,sl1:sl2,iblk) = c0

               call shortwave_ccsm3(nx_block,     ny_block,           &
                              icells,                                 &
                              indxi,             indxj,               &
                              aicen(:,:,n,iblk), vicen(:,:,n,iblk),   &
                              vsnon(:,:,n,iblk),                      &
                              trcrn(:,:,nt_Tsfc,n,iblk),              &
                              swvdr(:,:,  iblk), swvdf(:,:,  iblk),   &
                              swidr(:,:,  iblk), swidf(:,:,  iblk),   &
                              alvdrn(:,:,n,iblk),alidrn(:,:,n,iblk),  &
                              alvdfn(:,:,n,iblk),alidfn(:,:,n,iblk),  &
                              fswsfcn(:,:,n,iblk),fswintn(:,:,n,iblk),&
                              fswthrun(:,:,n,iblk),                   &
                              Iswabsn(:,:,il1:il2,iblk),              &
                              albicen(:,:,n,iblk),albsnon(:,:,n,iblk))


            enddo  ! ncat
         enddo     ! nblocks

      endif

      !-----------------------------------------------------------------
      ! Aggregate albedos 
      !-----------------------------------------------------------------

      do iblk=1,nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do n = 1, ncat

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

            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

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

         enddo  ! ncat

      !----------------------------------------------------------------
      ! Store grid box mean albedos and fluxes before scaling by aice
      !----------------------------------------------------------------

         do j = 1, ny_block
         do i = 1, nx_block
            alvdf_gbm  (i,j,iblk) = alvdf  (i,j,iblk)
            alidf_gbm  (i,j,iblk) = alidf  (i,j,iblk)
            alvdr_gbm  (i,j,iblk) = alvdr  (i,j,iblk)
            alidr_gbm  (i,j,iblk) = alidr  (i,j,iblk)

            ! for history averaging
            cszn = c0
            if (coszen(i,j,iblk) > puny) cszn = c1
            do n = 1, nstreams
               albcnt(i,j,iblk,n) = albcnt(i,j,iblk,n) + cszn
            enddo
         enddo
         enddo

      enddo     ! nblocks

      end subroutine init_shortwave

!=======================================================================
!BOP
!
! !IROUTINE: shortwave_ccsm3 - driver for CCSM3 shortwave radiation
!
! !INTERFACE:
!
      subroutine shortwave_ccsm3 (nx_block, ny_block, &
                                  icells,             &
                                  indxi,    indxj,    &
                                  aicen,    vicen,    &
                                  vsnon,    Tsfcn,    &
                                  swvdr,    swvdf,    &
                                  swidr,    swidf,    &
                                  alvdrn,   alidrn,   &
                                  alvdfn,   alidfn,   &
                                  fswsfc,   fswint,   &
                                  fswthru,  Iswabs,   &
                                  albin,    albsn)
!
! !DESCRIPTION:
!
! Driver for basic solar radiation from CCSM3.  Albedos and absorbed solar.
!
! !REVISION HISTORY:
!
! authors:  same as module
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells              ! number of ice-covered grid cells

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi    , & ! indices for ice-covered cells
         indxj

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         aicen    , & ! concentration of ice per category
         vicen    , & ! volume of ice per category
         vsnon    , & ! volume of ice per category
         Tsfcn    , & ! surface temperature
         swvdr    , & ! sw down, visible, direct  (W/m^2)
         swvdf    , & ! sw down, visible, diffuse (W/m^2)
         swidr    , & ! sw down, near IR, direct  (W/m^2)
         swidf        ! sw down, near IR, diffuse (W/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(out) :: &
         alvdrn   , & ! visible, direct, avg   (fraction)
         alidrn   , & ! near-ir, direct, avg   (fraction)
         alvdfn   , & ! visible, diffuse, avg  (fraction)
         alidfn   , & ! near-ir, diffuse, avg  (fraction)
         fswsfc   , & ! SW absorbed at ice/snow surface (W m-2)
         fswint   , & ! SW absorbed in ice interior, below surface (W m-2)
         fswthru  , & ! SW through ice to ocean (W m-2)
         albin    , & ! bare ice albedo
         albsn        ! snow albedo

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr), &
           intent(out) :: &
         Iswabs       ! SW absorbed in particular layer (W m-2)
!
!EOP
!
      ! ice and snow albedo for each category

      real (kind=dbl_kind), dimension (nx_block,ny_block):: &
         alvdrni, & ! visible, direct, ice    (fraction)
         alidrni, & ! near-ir, direct, ice    (fraction)
         alvdfni, & ! visible, diffuse, ice   (fraction)
         alidfni, & ! near-ir, diffuse, ice   (fraction)
         alvdrns, & ! visible, direct, snow   (fraction)
         alidrns, & ! near-ir, direct, snow   (fraction)
         alvdfns, & ! visible, diffuse, snow  (fraction)
         alidfns    ! near-ir, diffuse, snow  (fraction)

      !-----------------------------------------------------------------
      ! Compute albedos for ice and snow.
      !-----------------------------------------------------------------

      if (trim(albedo_type) == 'constant') then
         call constant_albedos (nx_block,   ny_block, &
                                icells,               &
                                indxi,      indxj,    &
                                aicen,                &
                                vsnon,      Tsfcn,    &
                                alvdrni,    alidrni,  &
                                alvdfni,    alidfni,  &
                                alvdrns,    alidrns,  &
                                alvdfns,    alidfns,  &
                                alvdrn,     alidrn,   &
                                alvdfn,     alidfn,   &
                                albin,      albsn)
      else ! default
         call compute_albedos (nx_block,   ny_block, &
                               icells,               &
                               indxi,      indxj,    &
                               aicen,      vicen,    &
                               vsnon,      Tsfcn,    &
                               alvdrni,    alidrni,  &
                               alvdfni,    alidfni,  &
                               alvdrns,    alidrns,  &
                               alvdfns,    alidfns,  &
                               alvdrn,     alidrn,   &
                               alvdfn,     alidfn,   &
                               albin,      albsn)
      endif

      !-----------------------------------------------------------------
      ! Compute solar radiation absorbed in ice and penetrating to ocean.
      !-----------------------------------------------------------------

      call absorbed_solar  (nx_block,   ny_block, &
                            icells,               &
                            indxi,      indxj,    &
                            aicen,                &
                            vicen,      vsnon,    &
                            swvdr,      swvdf,    &
                            swidr,      swidf,    &
                            alvdrni,    alvdfni,  &
                            alidrni,    alidfni,  &
                            alvdrns,    alvdfns,  &
                            alidrns,    alidfns,  &
                            fswsfc,     fswint,   &
                            fswthru,    Iswabs)

      end subroutine shortwave_ccsm3

!=======================================================================
!BOP
!
! !IROUTINE: compute_albedos - compute albedos for each thickness ategory
!
! !INTERFACE:
!
      subroutine compute_albedos (nx_block, ny_block, &
                                  icells,             &
                                  indxi,    indxj,    &
                                  aicen,    vicen,    &
                                  vsnon,    Tsfcn,    &
                                  alvdrni,  alidrni,  &
                                  alvdfni,  alidfni,  &
                                  alvdrns,  alidrns,  &
                                  alvdfns,  alidfns,  &
                                  alvdrn,   alidrn,   &
                                  alvdfn,   alidfn,   &
                                  albin,    albsn)
!
! !DESCRIPTION:
!
! Compute albedos for each thickness category
!
! !REVISION HISTORY:
!
! authors:  same as module
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells              ! number of ice-covered grid cells

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi   , & ! compressed indices for ice-covered cells
         indxj

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         aicen   , & ! concentration of ice per category
         vicen   , & ! volume of ice per category
         vsnon   , & ! volume of ice per category
         Tsfcn       ! surface temperature

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(out) :: &
         alvdrni  , & ! visible, direct, ice   (fraction)
         alidrni  , & ! near-ir, direct, ice   (fraction)
         alvdfni  , & ! visible, diffuse, ice  (fraction)
         alidfni  , & ! near-ir, diffuse, ice  (fraction)
         alvdrns  , & ! visible, direct, snow  (fraction)
         alidrns  , & ! near-ir, direct, snow  (fraction)
         alvdfns  , & ! visible, diffuse, snow (fraction)
         alidfns  , & ! near-ir, diffuse, snow (fraction)
         alvdrn   , & ! visible, direct, avg   (fraction)
         alidrn   , & ! near-ir, direct, avg   (fraction)
         alvdfn   , & ! visible, diffuse, avg  (fraction)
         alidfn   , & ! near-ir, diffuse, avg  (fraction)
         albin    , & ! bare ice 
         albsn        ! snow 
!
!EOP
!
      real (kind=dbl_kind), parameter :: &
         ahmax     = p5          , & ! thickness above which ice albedo 
                                     ! is constant (m)
         dT_mlt    = c1          , & ! change in temp to give dalb_mlt 
                                     ! albedo change
         dalb_mlt  = -0.075_dbl_kind, & ! albedo change per dT_mlt change
                                     ! in temp for ice
         dalb_mltv = -p1         , & ! albedo vis change per dT_mlt change
                                     ! in temp for snow
         dalb_mlti = -p15            ! albedo nir change per dT_mlt change
                                     ! in temp for snow

      integer (kind=int_kind) :: &
         i, j, n

      real (kind=dbl_kind) :: &
         hi  , & ! ice thickness  (m)
         hs  , & ! snow thickness  (m)
         albo, & ! effective ocean albedo, function of ice thickness
         fh  , & ! piecewise linear function of thickness
         fT  , & ! piecewise linear function of surface temperature
         dTs , & ! difference of Tsfc and Timelt
         fhtan,& ! factor used in albedo dependence on ice thickness
         asnow   ! fractional area of snow cover

      integer (kind=int_kind) :: &
         ij      ! horizontal index, combines i and j loops

      fhtan = atan(ahmax*c4)

      do j = 1, ny_block
      do i = 1, nx_block
         alvdrni(i,j) = albocn
         alidrni(i,j) = albocn
         alvdfni(i,j) = albocn
         alidfni(i,j) = albocn

         alvdrns(i,j) = albocn
         alidrns(i,j) = albocn
         alvdfns(i,j) = albocn
         alidfns(i,j) = albocn

         alvdrn(i,j) = albocn
         alidrn(i,j) = albocn
         alvdfn(i,j) = albocn
         alidfn(i,j) = albocn

         albin(i,j) = c0
         albsn(i,j) = c0
      enddo         
      enddo         

      !-----------------------------------------------------------------
      ! Compute albedo for each thickness category.
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         hi = vicen(i,j) / aicen(i,j)
         hs = vsnon(i,j) / aicen(i,j)            

         ! bare ice, thickness dependence
         fh = min(atan(hi*c4)/fhtan,c1)
         albo = albocn*(c1-fh)
         alvdfni(i,j) = albicev*fh + albo
         alidfni(i,j) = albicei*fh + albo

         ! bare ice, temperature dependence
         dTs = Timelt - Tsfcn(i,j)
         fT = min(dTs/dT_mlt-c1,c0)
         alvdfni(i,j) = alvdfni(i,j) - dalb_mlt*fT
         alidfni(i,j) = alidfni(i,j) - dalb_mlt*fT

         ! avoid negative albedos for thin, bare, melting ice
         alvdfni(i,j) = max (alvdfni(i,j), albocn)
         alidfni(i,j) = max (alidfni(i,j), albocn)

         if (hs > puny) then

            alvdfns(i,j) = albsnowv
            alidfns(i,j) = albsnowi

            ! snow on ice, temperature dependence
            alvdfns(i,j) = alvdfns(i,j) - dalb_mltv*fT
            alidfns(i,j) = alidfns(i,j) - dalb_mlti*fT

         endif                  ! hs > puny

         ! direct albedos (same as diffuse for now)
         alvdrni(i,j) = alvdfni(i,j)
         alidrni(i,j) = alidfni(i,j)
         alvdrns(i,j) = alvdfns(i,j)
         alidrns(i,j) = alidfns(i,j)


         ! fractional area of snow cover
         if (hs > puny) then
            asnow = hs / (hs + snowpatch)
         else
            asnow = c0
         endif

         ! combine ice and snow albedos (for coupler)
         alvdfn(i,j) = alvdfni(i,j)*(c1-asnow) + &
                       alvdfns(i,j)*asnow
         alidfn(i,j) = alidfni(i,j)*(c1-asnow) + &
                       alidfns(i,j)*asnow
         alvdrn(i,j) = alvdrni(i,j)*(c1-asnow) + &
                       alvdrns(i,j)*asnow
         alidrn(i,j) = alidrni(i,j)*(c1-asnow) + &
                       alidrns(i,j)*asnow

         ! save ice and snow albedos (for history)
         albin(i,j) = awtvdr*alvdrni(i,j) + awtidr*alidrni(i,j) &
                    + awtvdf*alvdfni(i,j) + awtidf*alidfni(i,j) 
         albsn(i,j) = awtvdr*alvdrns(i,j) + awtidr*alidrns(i,j) &
                    + awtvdf*alvdfns(i,j) + awtidf*alidfns(i,j) 

      enddo                     ! ij

      end subroutine compute_albedos

!=======================================================================
!BOP
!
! !IROUTINE: constant_albedos - set albedos for each thickness ategory
!
! !INTERFACE:
!
      subroutine constant_albedos (nx_block, ny_block, &
                                  icells,             &
                                  indxi,    indxj,    &
                                  aicen,              &
                                  vsnon,    Tsfcn,    &
                                  alvdrni,  alidrni,  &
                                  alvdfni,  alidfni,  &
                                  alvdrns,  alidrns,  &
                                  alvdfns,  alidfns,  &
                                  alvdrn,   alidrn,   &
                                  alvdfn,   alidfn,   &
                                  albin,    albsn)
!
! !DESCRIPTION:
!
! Compute albedos for each thickness category
!
! !REVISION HISTORY:
!
! authors:  same as module
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells              ! number of ice-covered grid cells

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi   , & ! compressed indices for ice-covered cells
         indxj

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         aicen   , & ! concentration of ice per category
         vsnon   , & ! volume of ice per category
         Tsfcn       ! surface temperature

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(out) :: &
         alvdrni  , & ! visible, direct, ice   (fraction)
         alidrni  , & ! near-ir, direct, ice   (fraction)
         alvdfni  , & ! visible, diffuse, ice  (fraction)
         alidfni  , & ! near-ir, diffuse, ice  (fraction)
         alvdrns  , & ! visible, direct, snow  (fraction)
         alidrns  , & ! near-ir, direct, snow  (fraction)
         alvdfns  , & ! visible, diffuse, snow (fraction)
         alidfns  , & ! near-ir, diffuse, snow (fraction)
         alvdrn   , & ! visible, direct, avg   (fraction)
         alidrn   , & ! near-ir, direct, avg   (fraction)
         alvdfn   , & ! visible, diffuse, avg  (fraction)
         alidfn   , & ! near-ir, diffuse, avg  (fraction)
         albin    , & ! bare ice 
         albsn        ! snow 
!
!EOP
!
      real (kind=dbl_kind), parameter :: &
         warmice  = 0.68_dbl_kind, &
         coldice  = 0.70_dbl_kind, &
         warmsnow = 0.77_dbl_kind, &
         coldsnow = 0.81_dbl_kind

      integer (kind=int_kind) :: &
         i, j, n

      real (kind=dbl_kind) :: &
         hs      ! snow thickness  (m)

      integer (kind=int_kind) :: &
         ij      ! horizontal index, combines i and j loops

      do j = 1, ny_block
      do i = 1, nx_block
         alvdrn(i,j) = albocn
         alidrn(i,j) = albocn
         alvdfn(i,j) = albocn
         alidfn(i,j) = albocn

         albin(i,j) = c0
         albsn(i,j) = c0
      enddo
      enddo

      !-----------------------------------------------------------------
      ! Compute albedo for each thickness category.
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         hs = vsnon(i,j) / aicen(i,j)

         if (hs > puny) then
            ! snow, temperature dependence
            if (Tsfcn(i,j) >= -c2*puny) then
               alvdfn(i,j) = warmsnow
               alidfn(i,j) = warmsnow
            else
               alvdfn(i,j) = coldsnow
               alidfn(i,j) = coldsnow
            endif
         else      ! hs < puny
            ! bare ice, temperature dependence
            if (Tsfcn(i,j) >= -c2*puny) then
               alvdfn(i,j) = warmice
               alidfn(i,j) = warmice
            else
               alvdfn(i,j) = coldice
               alidfn(i,j) = coldice
            endif
         endif                  ! hs > puny

         ! direct albedos (same as diffuse for now)
         alvdrn (i,j) = alvdfn(i,j)
         alidrn (i,j) = alidfn(i,j)

         alvdrni(i,j) = alvdrn(i,j)
         alidrni(i,j) = alidrn(i,j)
         alvdrns(i,j) = alvdrn(i,j)
         alidrns(i,j) = alidrn(i,j)
         alvdfni(i,j) = alvdfn(i,j)
         alidfni(i,j) = alidfn(i,j)
         alvdfns(i,j) = alvdfn(i,j)
         alidfns(i,j) = alidfn(i,j)

         ! save ice and snow albedos (for history)
         albin(i,j) = awtvdr*alvdrni(i,j) + awtidr*alidrni(i,j) &
                    + awtvdf*alvdfni(i,j) + awtidf*alidfni(i,j) 
         albsn(i,j) = awtvdr*alvdrns(i,j) + awtidr*alidrns(i,j) &
                    + awtvdf*alvdfns(i,j) + awtidf*alidfns(i,j) 

      enddo                     ! ij

      end subroutine constant_albedos

!=======================================================================
!BOP
!
! !ROUTINE: absorbed_solar - shortwave radiation absorbed by ice, ocean
!
! !DESCRIPTION:
!
! Compute solar radiation absorbed in ice and penetrating to ocean
!
! !REVISION HISTORY:
!
! authors William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! !INTERFACE:
!
      subroutine absorbed_solar (nx_block, ny_block, &
                                 icells,             &
                                 indxi,    indxj,    &
                                 aicen,              &
                                 vicen,    vsnon,    &
                                 swvdr,    swvdf,    &
                                 swidr,    swidf,    &
                                 alvdrni,  alvdfni,  &
                                 alidrni,  alidfni,  &
                                 alvdrns,  alvdfns,  &
                                 alidrns,  alidfns,  &
                                 fswsfc,   fswint,   &
                                 fswthru,  Iswabs)
!
! !USES:
!
	use ice_therm_vertical, only: heat_capacity
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells          ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         aicen       , & ! fractional ice area
         vicen       , & ! ice volume
         vsnon       , & ! snow volume
         swvdr       , & ! sw down, visible, direct  (W/m^2)
         swvdf       , & ! sw down, visible, diffuse (W/m^2)
         swidr       , & ! sw down, near IR, direct  (W/m^2)
         swidf       , & ! sw down, near IR, diffuse (W/m^2)
         alvdrni     , & ! visible, direct albedo,ice
         alidrni     , & ! near-ir, direct albedo,ice
         alvdfni     , & ! visible, diffuse albedo,ice
         alidfni     , & ! near-ir, diffuse albedo,ice
         alvdrns     , & ! visible, direct albedo, snow
         alidrns     , & ! near-ir, direct albedo, snow
         alvdfns     , & ! visible, diffuse albedo, snow
         alidfns         ! near-ir, diffuse albedo, snow

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(out):: &
         fswsfc      , & ! SW absorbed at ice/snow surface (W m-2)
         fswint      , & ! SW absorbed in ice interior, below surface (W m-2)
         fswthru         ! SW through ice to ocean (W m-2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr), &
         intent(out) :: &
         Iswabs          ! SW absorbed in particular layer (W m-2)
!
!EOP
!
      real (kind=dbl_kind), parameter :: &
         i0vis = 0.70_dbl_kind  ! fraction of penetrating solar rad (visible)

      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij          , & ! horizontal index, combines i and j loops
         k               ! ice layer index

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         fswpen      , & ! SW penetrating beneath surface (W m-2)
         trantop     , & ! transmitted frac of penetrating SW at layer top
         tranbot         ! transmitted frac of penetrating SW at layer bot

      real (kind=dbl_kind) :: &
         swabs       , & ! net SW down at surface (W m-2)
         swabsv      , & ! swabs in vis (wvlngth < 700nm)  (W/m^2)
         swabsi      , & ! swabs in nir (wvlngth > 700nm)  (W/m^2)
         fswpenvdr   , & ! penetrating SW, vis direct
         fswpenvdf   , & ! penetrating SW, vis diffuse
         hi          , & ! ice thickness (m)
         hs          , & ! snow thickness (m)
         hilyr       , & ! ice layer thickness
         asnow           ! fractional area of snow cover

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      rnilyr = real(nilyr,kind=dbl_kind)

      do j = 1, ny_block
      do i = 1, nx_block
         fswsfc (i,j) = c0
         fswint (i,j) = c0
         fswthru(i,j) = c0
         fswpen (i,j) = c0
         trantop(i,j) = c0
         tranbot(i,j) = c0
      enddo
      enddo
      Iswabs (:,:,:) = c0

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         hs  = vsnon(i,j) / aicen(i,j)

      !-----------------------------------------------------------------
      ! Fractional snow cover
      !-----------------------------------------------------------------
         if (hs > puny) then
            asnow = hs / (hs + snowpatch)
         else
            asnow = c0
         endif

      !-----------------------------------------------------------------
      ! Shortwave flux absorbed at surface, absorbed internally,
      !  and penetrating to mixed layer.
      ! This parameterization assumes that all IR is absorbed at the
      !  surface; only visible is absorbed in the ice interior or
      !  transmitted to the ocean.
      !-----------------------------------------------------------------

         swabsv  = swvdr(i,j) * ( (c1-alvdrni(i,j))*(c1-asnow) &
                                + (c1-alvdrns(i,j))*asnow ) &
                 + swvdf(i,j) * ( (c1-alvdfni(i,j))*(c1-asnow) &
                                + (c1-alvdfns(i,j))*asnow )

         swabsi  = swidr(i,j) * ( (c1-alidrni(i,j))*(c1-asnow) &
                                + (c1-alidrns(i,j))*asnow ) &
                 + swidf(i,j) * ( (c1-alidfni(i,j))*(c1-asnow) &
                                + (c1-alidfns(i,j))*asnow )

         swabs   = swabsv + swabsi

         fswpenvdr = swvdr(i,j) * (c1-alvdrni(i,j)) * (c1-asnow) * i0vis
         fswpenvdf = swvdf(i,j) * (c1-alvdfni(i,j)) * (c1-asnow) * i0vis

          ! no penetrating radiation in near IR
!         fswpenidr = swidr(i,j) * (c1-alidrni(i,j)) * (c1-asnow) * i0nir
!         fswpenidf = swidf(i,j) * (c1-alidfni(i,j)) * (c1-asnow) * i0nir  

         fswpen(i,j) = fswpenvdr + fswpenvdf
                      
         fswsfc(i,j) = swabs - fswpen(i,j)

         trantop(i,j) = c1  ! transmittance at top of ice

      enddo                     ! ij

      !-----------------------------------------------------------------
      ! penetrating SW absorbed in each ice layer
      !-----------------------------------------------------------------

      do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            hi  = vicen(i,j) / aicen(i,j)
            hilyr = hi / rnilyr

            tranbot(i,j) = exp (-kappav * hilyr * real(k,kind=dbl_kind))
            Iswabs(i,j,k) = fswpen(i,j) * (trantop(i,j)-tranbot(i,j))

            ! bottom of layer k = top of layer k+1
            trantop(i,j) = tranbot(i,j)

         enddo                  ! ij
      enddo                     ! nilyr

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         ! SW penetrating thru ice into ocean
         fswthru(i,j) = fswpen(i,j) * tranbot(i,j)

         ! SW absorbed in ice interior
         fswint(i,j)  = fswpen(i,j) - fswthru(i,j)

      enddo                     ! ij

      !----------------------------------------------------------------
      ! if zero-layer model (no heat capacity), no SW is absorbed in ice
      ! interior, so add to surface absorption
      !----------------------------------------------------------------

      if (.not. heat_capacity) then

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
        do ij = 1, icells
           i = indxi(ij)
           j = indxj(ij)

           ! SW absorbed at snow/ice surface
           fswsfc(i,j) = fswsfc(i,j) + fswint(i,j)

           ! SW absorbed in ice interior (nilyr = 1)
           fswint(i,j)   = c0
           Iswabs(i,j,1) = c0

        enddo                     ! ij

      endif                       ! heat_capacity

      end subroutine absorbed_solar

! End ccsm3 shortwave method
!=======================================================================
! Begin Delta-Eddington shortwave method
!
!BOP
!
! !IROUTINE: init_dEdd - initialize Delta-Eddington parameters
!
! !INTERFACE:
!
      subroutine init_dEdd
!
! !DESCRIPTION:
!
! Compute initial data for Delta-Eddington method, specifically, 
! the approximate exponential look-up table.
!
! !REVISION HISTORY:
!
! author:  Bruce P. Briegleb, NCAR
!
! !USES:
!
      use ice_domain_size
      use ice_blocks
      use ice_calendar
      use ice_domain
      use ice_flux
      use ice_grid
      use ice_itd
      use ice_meltpond
      use ice_orbital
      use ice_state
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
!     local temporary variables

      integer (kind=int_kind) :: &
         icells          ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         indxi, indxj    ! indirect indices for cells with aicen > puny

      ! other local variables
      ! snow variables for Delta-Eddington shortwave
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         fsn             ! snow horizontal fraction
      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr) :: &
         rhosnwn     , & ! snow density (kg/m3)
         rsnwn           ! snow grain radius (micrometers)

      ! pond variables for Delta-Eddington shortwave
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         fpn         , & ! pond fraction
         hpn             ! pond depth (m)

      integer (kind=int_kind) :: &
         i, j, ij    , & ! horizontal indices
         iblk        , & ! block index
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         n           , & ! thickness category index
         il1, il2    , & ! ice layer indices for eice
         sl1, sl2        ! snow layer indices for esno

      type (block) :: &
         this_block      ! block information for current block

      exp_min = exp(-c10)

      do iblk=1,nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         ! identify ice-ocean cells
         icells = 0
         do j = 1, ny_block
         do i = 1, nx_block
            if (tmask(i,j,iblk)) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
            endif
         enddo               ! i
         enddo               ! j

         call compute_coszen (nx_block,         ny_block,       &
                              icells,                           &
                              indxi,            indxj,          &
                              tlat  (:,:,iblk), tlon(:,:,iblk), &
                              coszen(:,:,iblk), dt)

         do n = 1, ncat

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

            il1 = ilyr1(n)
            il2 = ilyrn(n)
            sl1 = slyr1(n)
            sl2 = slyrn(n)

      ! note that rhoswn, rsnw, fp, hp and Sswabs ARE NOT dimensioned with ncat
      ! BPB 19 Dec 2006

            ! set snow properties
            call shortwave_dEdd_set_snow(nx_block, ny_block,           &
                              icells,                                  &
                              indxi,               indxj,              &
                              aicen(:,:,n,iblk),   vsnon(:,:,n,iblk),  &
                              trcrn(:,:,1,n,iblk), fsn,                &
                              rhosnwn,             rsnwn)


            if (.not. tr_pond) then

            ! set pond properties
               call shortwave_dEdd_set_pond(nx_block, ny_block,            &
                                 icells,                                   &
                                 indxi,               indxj,               &
                                 aicen(:,:,n,iblk),   trcrn(:,:,1,n,iblk), &
                                 fsn,                 fpn,                 &
                                 hpn)

            else

               fpn(:,:) = apondn(:,:,n,iblk)
               hpn(:,:) = hpondn(:,:,n,iblk)

            endif

#ifdef AEROFRC
            if (tr_aero) then

            tr_aero = .false.

            call shortwave_dEdd(nx_block,        ny_block,            &
                              icells,                                 &
                              indxi,             indxj,               &
                              coszen(:,:, iblk),                      &
                              aicen(:,:,n,iblk), vicen(:,:,n,iblk),   &
                              vsnon(:,:,n,iblk), fsn,                 &
                              rhosnwn,           rsnwn,               &
                              fpn,               hpn,                 &
                              trcrn(:,:,:,n,iblk),tarea(:,:,iblk),    &
                              swvdr(:,:,  iblk), swvdf(:,:,  iblk),   &
                              swidr(:,:,  iblk), swidf(:,:,  iblk),   &
                              dalvdrn_noaero(:,:,n,iblk),              &
                              dalvdfn_noaero(:,:,n,iblk),              &
                              dalidrn_noaero(:,:,n,iblk),              &
                              dalidfn_noaero(:,:,n,iblk),              &
                              dfswsfcn_noaero(:,:,n,iblk),             &
                              dfswintn_noaero(:,:,n,iblk),             &
                              dfswthrun_noaero(:,:,n,iblk),            &
                              dSswabsn_noaero(:,:,sl1:sl2,iblk),       &
                              dIswabsn_noaero(:,:,il1:il2,iblk),       &
                              dalbicen_noaero(:,:,n,iblk),             &
                              dalbsnon_noaero(:,:,n,iblk),             &
                              dalbpndn_noaero(:,:,n,iblk))

            tr_aero = .true.

            endif
#endif
#ifdef CCSM3FRC
            call shortwave_ccsm3(nx_block,     ny_block,              &
                              icells,                                 &
                              indxi,             indxj,               &
                              aicen(:,:,n,iblk), vicen(:,:,n,iblk),   &
                              vsnon(:,:,n,iblk),                      &
                              trcrn(:,:,nt_Tsfc,n,iblk),              &
                              swvdr(:,:,  iblk), swvdf(:,:,  iblk),   &
                              swidr(:,:,  iblk), swidf(:,:,  iblk),   &
                              dalvdrn_ccsm3(:,:,n,iblk),               &
                              dalidrn_ccsm3(:,:,n,iblk),               &
                              dalvdfn_ccsm3(:,:,n,iblk),               &
                              dalidfn_ccsm3(:,:,n,iblk),               &
                              dfswsfcn_ccsm3(:,:,n,iblk),              &
                              dfswintn_ccsm3(:,:,n,iblk),              &
                              dfswthrun_ccsm3(:,:,n,iblk),             &
                              dIswabsn_ccsm3(:,:,il1:il2,iblk),        &
                              dalbicen_ccsm3(:,:,n,iblk),              &
                              dalbsnon_ccsm3(:,:,n,iblk))
#endif
#ifdef PONDFRC
            if (tr_pond) then

            fpn(:,:) = c0
            hpn(:,:) = c0

            call shortwave_dEdd(nx_block,        ny_block,            &
                              icells,                                 &
                              indxi,             indxj,               &
                              coszen(:,:, iblk),                      &
                              aicen(:,:,n,iblk), vicen(:,:,n,iblk),   &
                              vsnon(:,:,n,iblk), fsn,                 &
                              rhosnwn,           rsnwn,               &
                              fpn,               hpn,                 &
                              trcrn(:,:,:,n,iblk),tarea(:,:,iblk),    &
                              swvdr(:,:,  iblk), swvdf(:,:,  iblk),   &
                              swidr(:,:,  iblk), swidf(:,:,  iblk),   &
                              dalvdrn_nopond(:,:,n,iblk),              &
                              dalvdfn_nopond(:,:,n,iblk),              &
                              dalidrn_nopond(:,:,n,iblk),              &
                              dalidfn_nopond(:,:,n,iblk),              &
                              dfswsfcn_nopond(:,:,n,iblk),             &
                              dfswintn_nopond(:,:,n,iblk),             &
                              dfswthrun_nopond(:,:,n,iblk),            &
                              dSswabsn_nopond(:,:,sl1:sl2,iblk),       &
                              dIswabsn_nopond(:,:,il1:il2,iblk),       &
                              dalbicen_nopond(:,:,n,iblk),             &
                              dalbsnon_nopond(:,:,n,iblk),             &
                              dalbpndn_nopond(:,:,n,iblk))

            fpn(:,:) = apondn(:,:,n,iblk)
            hpn(:,:) = hpondn(:,:,n,iblk)

            endif
#endif

            call shortwave_dEdd(nx_block,        ny_block,            &
                              icells,                                 &
                              indxi,             indxj,               &
                              coszen(:,:, iblk),                      &
                              aicen(:,:,n,iblk), vicen(:,:,n,iblk),   &
                              vsnon(:,:,n,iblk), fsn,                 &
                              rhosnwn,           rsnwn,               &
                              fpn,               hpn,                 &
                              trcrn(:,:,:,n,iblk),tarea(:,:,iblk),    &
                              swvdr(:,:,  iblk), swvdf(:,:,  iblk),   &
                              swidr(:,:,  iblk), swidf(:,:,  iblk),   &
                              alvdrn(:,:,n,iblk),alvdfn(:,:,n,iblk),  &
                              alidrn(:,:,n,iblk),alidfn(:,:,n,iblk),  &
                              fswsfcn(:,:,n,iblk),fswintn(:,:,n,iblk),&
                              fswthrun(:,:,n,iblk), &
                              Sswabsn(:,:,sl1:sl2,iblk), &
                              Iswabsn(:,:,il1:il2,iblk), &
                              albicen(:,:,n,iblk),albsnon(:,:,n,iblk),&
                              albpndn(:,:,n,iblk))
#ifdef AEROFRC
            dalvdrn_noaero(:,:,n,iblk) = dalvdrn_noaero(:,:,n,iblk)-alvdrn(:,:,n,iblk)
            dalvdfn_noaero(:,:,n,iblk) = dalvdfn_noaero(:,:,n,iblk)-alvdfn(:,:,n,iblk)
            dalidrn_noaero(:,:,n,iblk) = dalidrn_noaero(:,:,n,iblk)-alidrn(:,:,n,iblk)
            dalidfn_noaero(:,:,n,iblk) = dalidfn_noaero(:,:,n,iblk)-alidfn(:,:,n,iblk)
            dfswsfcn_noaero(:,:,n,iblk) = dfswsfcn_noaero(:,:,n,iblk)-fswsfcn(:,:,n,iblk)
            dfswintn_noaero(:,:,n,iblk) = dfswintn_noaero(:,:,n,iblk)-fswintn(:,:,n,iblk)
            dfswthrun_noaero(:,:,n,iblk) = dfswthrun_noaero(:,:,n,iblk)-fswthrun(:,:,n,iblk)
            dfswabsn_noaero(:,:,n,iblk) = dfswsfcn_noaero(:,:,n,iblk)+dfswintn_noaero(:,:,n,iblk)+dfswthrun_noaero(:,:,n,iblk)
            dalbicen_noaero(:,:,n,iblk) = dalbicen_noaero(:,:,n,iblk)-albicen(:,:,n,iblk)
            dalbsnon_noaero(:,:,n,iblk) = dalbsnon_noaero(:,:,n,iblk)-albsnon(:,:,n,iblk)
            dalbpndn_noaero(:,:,n,iblk) = dalbpndn_noaero(:,:,n,iblk)-albpndn(:,:,n,iblk)
            dSswabsn_noaero(:,:,sl1:sl2,iblk) = dSswabsn_noaero(:,:,sl1:sl2,iblk)-Sswabsn(:,:,sl1:sl2,iblk)
            dIswabsn_noaero(:,:,il1:il2,iblk) = dIswabsn_noaero(:,:,il1:il2,iblk)-Iswabsn(:,:,il1:il2,iblk)
#endif
#ifdef CCSM3FRC
            dalvdrn_ccsm3(:,:,n,iblk) = dalvdrn_ccsm3(:,:,n,iblk)-alvdrn(:,:,n,iblk)
            dalvdfn_ccsm3(:,:,n,iblk) = dalvdfn_ccsm3(:,:,n,iblk)-alvdfn(:,:,n,iblk)
            dalidrn_ccsm3(:,:,n,iblk) = dalidrn_ccsm3(:,:,n,iblk)-alidrn(:,:,n,iblk)
            dalidfn_ccsm3(:,:,n,iblk) = dalidfn_ccsm3(:,:,n,iblk)-alidfn(:,:,n,iblk)
            dfswsfcn_ccsm3(:,:,n,iblk) = dfswsfcn_ccsm3(:,:,n,iblk)-fswsfcn(:,:,n,iblk)
            dfswintn_ccsm3(:,:,n,iblk) = dfswintn_ccsm3(:,:,n,iblk)-fswintn(:,:,n,iblk)
            dfswthrun_ccsm3(:,:,n,iblk) = dfswthrun_ccsm3(:,:,n,iblk)-fswthrun(:,:,n,iblk)
            dfswabsn_ccsm3(:,:,n,iblk) = dfswsfcn_ccsm3(:,:,n,iblk)+dfswintn_ccsm3(:,:,n,iblk)+dfswthrun_ccsm3(:,:,n,iblk)
            dalbicen_ccsm3(:,:,n,iblk) = dalbicen_ccsm3(:,:,n,iblk)-albicen(:,:,n,iblk)
            dalbsnon_ccsm3(:,:,n,iblk) = dalbsnon_ccsm3(:,:,n,iblk)-albsnon(:,:,n,iblk)
            dIswabsn_ccsm3(:,:,il1:il2,iblk) = dIswabsn_ccsm3(:,:,il1:il2,iblk)-Iswabsn(:,:,il1:il2,iblk)
#endif
#ifdef PONDFRC
            dalvdrn_nopond(:,:,n,iblk) = dalvdrn_nopond(:,:,n,iblk)-alvdrn(:,:,n,iblk)
            dalvdfn_nopond(:,:,n,iblk) = dalvdfn_nopond(:,:,n,iblk)-alvdfn(:,:,n,iblk)
            dalidrn_nopond(:,:,n,iblk) = dalidrn_nopond(:,:,n,iblk)-alidrn(:,:,n,iblk)
            dalidfn_nopond(:,:,n,iblk) = dalidfn_nopond(:,:,n,iblk)-alidfn(:,:,n,iblk)
            dfswsfcn_nopond(:,:,n,iblk) = dfswsfcn_nopond(:,:,n,iblk)-fswsfcn(:,:,n,iblk)
            dfswintn_nopond(:,:,n,iblk) = dfswintn_nopond(:,:,n,iblk)-fswintn(:,:,n,iblk)
            dfswthrun_nopond(:,:,n,iblk) = dfswthrun_nopond(:,:,n,iblk)-fswthrun(:,:,n,iblk)
            dfswabsn_nopond(:,:,n,iblk) = dfswsfcn_nopond(:,:,n,iblk)+dfswintn_nopond(:,:,n,iblk)+dfswthrun_nopond(:,:,n,iblk)
            dalbicen_nopond(:,:,n,iblk) = dalbicen_nopond(:,:,n,iblk)-albicen(:,:,n,iblk)
            dalbsnon_nopond(:,:,n,iblk) = dalbsnon_nopond(:,:,n,iblk)-albsnon(:,:,n,iblk)
            dalbpndn_nopond(:,:,n,iblk) = dalbpndn_nopond(:,:,n,iblk)-albpndn(:,:,n,iblk)
            dSswabsn_nopond(:,:,sl1:sl2,iblk) = dSswabsn_nopond(:,:,sl1:sl2,iblk)-Sswabsn(:,:,sl1:sl2,iblk)
            dIswabsn_nopond(:,:,il1:il2,iblk) = dIswabsn_nopond(:,:,il1:il2,iblk)-Iswabsn(:,:,il1:il2,iblk)
#endif
         enddo  ! ncat
      enddo     ! nblocks
 
      end subroutine init_dEdd
 
!=======================================================================
!BOP
!
! !IROUTINE: shortwave_dEdd - driver for Delta-Eddington shortwave
!
! !INTERFACE:
!
      subroutine shortwave_dEdd  (nx_block, ny_block,    &
                                  icells,   indxi,       &
                                  indxj,    coszen,      &
                                  aice,     vice,        &
                                  vsno,     fs,          & 
                                  rhosnw,   rsnw,        &
                                  fp,       hp,          &
                                  trcr,     tarea,       &
                                  swvdr,    swvdf,       &
                                  swidr,    swidf,       &
                                  alvdr,    alvdf,       &
                                  alidr,    alidf,       &
                                  fswsfc,   fswint,      &
                                  fswthru,  Sswabs,      &
                                  Iswabs,   albice,      &
                                  albsno,   albpnd)
!
! !DESCRIPTION:
!
!   Compute snow/bare ice/ponded ice shortwave albedos, absorbed and transmitted 
!   flux using the Delta-Eddington solar radiation method as described in:
!
!   "A Delta-Eddington Multiple Scattering Parameterization for Solar Radiation
!        in the Sea Ice Component of the Community Climate System Model"
!            B.P.Briegleb and B.Light   NCAR/TN-472+STR  February 2007
!
!   Compute shortwave albedos and fluxes for three surface types: 
!   snow over ice, bare ice and ponded ice. 
!   
!   Albedos and fluxes are output for later use by thermodynamic routines. 
!   Invokes three calls to compute_dEdd, which sets inherent optical properties 
!   appropriate for the surface type. Within compute_dEdd, a call to solution_dEdd 
!   evaluates the Delta-Eddington solution. The final albedos and fluxes are then
!   evaluated in compute_dEdd. Albedos and fluxes are transferred to output in 
!   this routine.
!
!   NOTE regarding albedo diagnostics:  This method yields zero albedo values
!   if there is no incoming solar and thus the albedo diagnostics are masked
!   out when the sun is below the horizon.  To estimate albedo from the history 
!   output (post-processing), compute ice albedo using
!   (1 - albedo)*swdn = swabs. -ECH
!
! !REVISION HISTORY:
!
! author:  Bruce P. Briegleb, NCAR 
! update:  8 February 2007
! update:  September 2008 added aerosols
!
! !USES:
!
      use ice_calendar
      use ice_state, only: nt_aero, tr_aero
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), &
         intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells                ! number of ice-covered grid cells

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi   , & ! compressed indices for ice-covered cells
         indxj

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         coszen  , & ! cosine of solar zenith angle 
         aice    , & ! concentration of ice 
         vice    , & ! volume of ice 
         vsno    , & ! volume of snow 
         fs          ! horizontal coverage of snow

      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr), &
         intent(in) :: &
         rhosnw  , & ! density in snow layer (kg/m3)
         rsnw        ! grain radius in snow layer (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr), &
         intent(in) :: &
         trcr        ! aerosol tracers

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         tarea       ! t-grid cell area in m2

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         fp      , & ! pond fractional coverage (0 to 1) 
         hp      , & ! pond depth (m) 
         swvdr   , & ! sw down, visible, direct  (W/m^2)
         swvdf   , & ! sw down, visible, diffuse (W/m^2)
         swidr   , & ! sw down, near IR, direct  (W/m^2)
         swidf       ! sw down, near IR, diffuse (W/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(out) :: &
         alvdr   , & ! visible, direct, albedo (fraction) 
         alvdf   , & ! visible, diffuse, albedo (fraction) 
         alidr   , & ! near-ir, direct, albedo (fraction) 
         alidf   , & ! near-ir, diffuse, albedo (fraction) 
         fswsfc  , & ! SW absorbed at snow/bare ice/pondedi ice surface (W m-2)
         fswint  , & ! SW interior absorption (below surface, above ocean,W m-2)
         fswthru     ! SW through snow/bare ice/ponded ice into ocean (W m-2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr), &
         intent(out) :: &
         Sswabs      ! SW absorbed in snow layer (W m-2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr), &
         intent(out) :: &
         Iswabs      ! SW absorbed in ice layer (W m-2)

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(out) :: &
         albice  , & ! bare ice albedo, for history  
         albsno  , & ! snow albedo, for history  
         albpnd      ! pond albedo, for history  
!      
!EOP
! 
! !LOCAL PARAMETERS:
!
      real (kind=dbl_kind),dimension (nx_block,ny_block) :: &
         fnidr        ! fraction of direct to total down surface flux in nir

      real (kind=dbl_kind), dimension(nx_block,ny_block) :: &
         hs       , & ! snow thickness (all snow layers, m)
         hi       , & ! ice thickness (all sea ice layers, m)
         fi           ! snow/bare ice fractional coverage (0 to 1)

      real (kind=dbl_kind), dimension (nx_block,ny_block,4*n_aeromx) :: &
         aero_mp      ! aerosol mass path in kg/m2

      integer (kind=int_kind), dimension(nx_block,ny_block) :: &
         srftyp       ! surface type over ice: (0=air, 1=snow, 2=pond)
 
      integer (kind=int_kind) :: &
         i        , & ! longitude index
         j        , & ! latitude index
         ij       , & ! horizontal index, combines i and j loops
         k        , & ! level index
         na       , & ! aerosol index
         icells_DE    ! number of cells in Delta-Eddington calculation
 
      integer (kind=int_kind), dimension (nx_block*ny_block) :: &
         indxi_DE , & ! compressed indices for Delta-Eddington cells
         indxj_DE  

      real (kind=dbl_kind) :: &
         hpmin    , & ! minimum allowed melt pond depth
         hsmax    , & ! maximum snow depth below which Sswabs adjustment
         hs_ssl   , & ! assumed snow surface scattering layer for Sswabs adj
         frcadj       ! fractional Sswabs adjustment
      data hpmin  / .005_dbl_kind /
      data hs_ssl / .040_dbl_kind /

      ! for printing points
      integer (kind=int_kind) :: &
         n            ! point number for prints

      logical (kind=log_kind) :: &
         dbug         ! true/false flag

      real (kind=dbl_kind) :: &
               swdn  , & ! swvdr(i,j)+swvdf(i,j)+swidr(i,j)+swidf(i,j)
               swab  , & ! fswsfc(i,j)+fswint(i,j)+fswthru(i,j)
               swalb     ! (1.-swab/(swdn+.0001))

      ! for history
      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         avdrl   , & ! visible, direct, albedo (fraction) 
         avdfl   , & ! visible, diffuse, albedo (fraction) 
         aidrl   , & ! near-ir, direct, albedo (fraction) 
         aidfl       ! near-ir, diffuse, albedo (fraction) 

!-----------------------------------------------------------------------

      do j = 1, ny_block
      do i = 1, nx_block
         ! zero storage albedos and fluxes for accumulation over surface types:
         hs(i,j)       = c0
         hi(i,j)       = c0
         fi(i,j)       = c0
         srftyp(i,j)   =  0
         alvdr(i,j)    = c0
         alvdf(i,j)    = c0
         alidr(i,j)    = c0
         alidf(i,j)    = c0
         avdrl(i,j)    = c0
         avdfl(i,j)    = c0
         aidrl(i,j)    = c0
         aidfl(i,j)    = c0
         fswsfc(i,j)   = c0
         fswint(i,j)   = c0
         fswthru(i,j)  = c0
      ! compute fraction of nir down direct to total over all points:
         fnidr(i,j) = c0
         if( swidr(i,j) + swidf(i,j) > puny ) then
            fnidr(i,j) = swidr(i,j)/(swidr(i,j)+swidf(i,j))
         endif
         albice(i,j)    = c0
         albsno(i,j)    = c0
         albpnd(i,j)    = c0
      enddo
      enddo
      Sswabs(:,:,:) = c0
      Iswabs(:,:,:) = c0

      ! compute aerosol mass path

      aero_mp(:,:,:) = c0
      if( tr_aero ) then
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         ! assume 4 layers for each aerosol, a snow SSL, snow below SSL,
         ! sea ice SSL, and sea ice below SSL, in that order.
         do na=1,4*n_aero,4
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               ! sea ice points with sun above horizon
               if (aice(i,j) > puny .and. coszen(i,j) > puny) then
                  aero_mp(i,j,na  ) = trcr(i,j,nt_aero-1+na  )*vsno(i,j)
                  aero_mp(i,j,na+1) = trcr(i,j,nt_aero-1+na+1)*vsno(i,j)
                  aero_mp(i,j,na+2) = trcr(i,j,nt_aero-1+na+2)*vice(i,j)
                  aero_mp(i,j,na+3) = trcr(i,j,nt_aero-1+na+3)*vice(i,j)
               endif                  ! aice > 0 and coszen > 0
            enddo                     ! ij
         enddo      ! na
      endif      ! if aerosols

      ! compute shortwave radiation accounting for snow/ice (both snow over 
      ! ice and bare ice) and ponded ice (if any):
 
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      ! find bare ice points 
      icells_DE = 0
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         ! sea ice points with sun above horizon
         if (aice(i,j) > puny .and. coszen(i,j) > puny) then
            ! evaluate sea ice thickness and fraction
            hi(i,j)  = vice(i,j) / aice(i,j)
            fi(i,j)  = c1 - fs(i,j) - fp(i,j)
            ! bare sea ice points
            if(fi(i,j) > c0) then
              icells_DE = icells_DE + 1
              indxi_DE(icells_DE) = i
              indxj_DE(icells_DE) = j 
              ! bare ice
              srftyp(i,j) = 0
            endif               ! fi > 0
         endif                  ! aice > 0 and coszen > 0
      enddo                     ! ij

      ! calculate bare sea ice
      if (tflag) call t_startf('cice_swdedd_computedEdd1')
      tstr = '1'
      call compute_dEdd                                    &
            (nx_block,ny_block,                            &
             icells_DE, indxi_DE, indxj_DE, fnidr, coszen, &
             swvdr,     swvdf,    swidr,    swidf, srftyp, &
             hs,        rhosnw,   rsnw,     hi,    hp,     &
             fi,       aero_mp,   avdrl,    avdfl,         &
                                  aidrl,    aidfl,         &
                                  fswsfc,   fswint,        &
                                  fswthru,  Sswabs(:,:,:), &
                                  Iswabs)
      if (tflag) call t_stopf('cice_swdedd_computedEdd1')

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells_DE
         i = indxi_DE(ij)
         j = indxj_DE(ij)
         alvdr(i,j)   = alvdr(i,j)   + avdrl(i,j) *fi(i,j)
         alvdf(i,j)   = alvdf(i,j)   + avdfl(i,j) *fi(i,j)
         alidr(i,j)   = alidr(i,j)   + aidrl(i,j) *fi(i,j)
         alidf(i,j)   = alidf(i,j)   + aidfl(i,j) *fi(i,j)
         ! for history
         albice(i,j) = albice(i,j) &
                     + awtvdr*avdrl(i,j) + awtidr*aidrl(i,j) &
                     + awtvdf*avdfl(i,j) + awtidf*aidfl(i,j) 
      enddo

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      ! find snow-covered ice points
      icells_DE = 0
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         ! sea ice points with sun above horizon
         if (aice(i,j) > puny .and. coszen(i,j) > puny) then
            ! evaluate snow thickness
            hs(i,j)  = vsno(i,j) / aice(i,j)
            ! snow-covered sea ice points
            if(fs(i,j) > c0) then
              icells_DE = icells_DE + 1
              indxi_DE(icells_DE) = i
              indxj_DE(icells_DE) = j 
              ! snow-covered ice
              srftyp(i,j) = 1
            endif               ! fs > 0
         endif                  ! aice > 0 and coszen > 0
      enddo                     ! ij

      ! calculate snow covered sea ice
      if (tflag) call t_startf('cice_swdedd_computedEdd2')
      tstr = '2'
      call compute_dEdd                                    &
            (nx_block,ny_block,                            &
             icells_DE, indxi_DE, indxj_DE, fnidr, coszen, &
             swvdr,     swvdf,    swidr,    swidf, srftyp, &
             hs,        rhosnw,   rsnw,     hi,    hp,     &
             fs,       aero_mp,   avdrl,    avdfl,         &
                                  aidrl,    aidfl,         &
                                  fswsfc,   fswint,        &
                                  fswthru,  Sswabs(:,:,:), &
                                  Iswabs)
      if (tflag) call t_stopf('cice_swdedd_computedEdd2')

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells_DE
         i = indxi_DE(ij)
         j = indxj_DE(ij)
         alvdr(i,j)   = alvdr(i,j)   + avdrl(i,j) *fs(i,j)
         alvdf(i,j)   = alvdf(i,j)   + avdfl(i,j) *fs(i,j)
         alidr(i,j)   = alidr(i,j)   + aidrl(i,j) *fs(i,j)
         alidf(i,j)   = alidf(i,j)   + aidfl(i,j) *fs(i,j)
         ! for history
         albsno(i,j) = albsno(i,j) &
                     + awtvdr*avdrl(i,j) + awtidr*aidrl(i,j) &
                     + awtvdf*avdfl(i,j) + awtidf*aidfl(i,j) 
      enddo

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      ! find ponded points
      icells_DE = 0
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         hi(i,j) = c0
         ! sea ice points with sun above horizon
         if (aice(i,j) > puny .and. coszen(i,j) > puny) then
            hi(i,j)  = vice(i,j) / aice(i,j)
            ! if non-zero pond fraction and sufficient pond depth
            if( fp(i,j) > puny .and. hp(i,j) > hpmin ) then
               icells_DE = icells_DE + 1
               indxi_DE(icells_DE) = i
               indxj_DE(icells_DE) = j
               ! ponded ice
               srftyp(i,j)   = 2
            endif               
         endif                  ! aice > puny, coszen > puny
      enddo                     ! ij

      ! calculate ponded ice
      if (tflag) call t_startf('cice_swdedd_computedEdd3')
      tstr = '3'
      call compute_dEdd                                    &
            (nx_block,ny_block,                            &
             icells_DE, indxi_DE, indxj_DE, fnidr, coszen, &
             swvdr,     swvdf,    swidr,    swidf, srftyp, &
             hs,        rhosnw,   rsnw,     hi,    hp,     &
             fp,       aero_mp,   avdrl,    avdfl,         &
                                  aidrl,    aidfl,         &
                                  fswsfc,   fswint,        &
                                  fswthru,  Sswabs(:,:,:), &
                                  Iswabs)
      if (tflag) call t_stopf('cice_swdedd_computedEdd3')

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells_DE
         i = indxi_DE(ij)
         j = indxj_DE(ij)
         alvdr(i,j)   = alvdr(i,j)   + avdrl(i,j) *fp(i,j)
         alvdf(i,j)   = alvdf(i,j)   + avdfl(i,j) *fp(i,j)
         alidr(i,j)   = alidr(i,j)   + aidrl(i,j) *fp(i,j)
         alidf(i,j)   = alidf(i,j)   + aidfl(i,j) *fp(i,j)
         ! for history
         albpnd(i,j) = albpnd(i,j) &
                     + awtvdr*avdrl(i,j) + awtidr*aidrl(i,j) &
                     + awtvdf*avdfl(i,j) + awtidf*aidfl(i,j) 
      enddo

      dbug = .false.
      if (dbug .and. print_points) then
         do n = 1, npnt
            if (my_task == pmloc(n)) then
               i = piloc(n)
               j = pjloc(n)
               if( coszen(i,j) > .01_dbl_kind ) then
               write(nu_diag,*) ' my_task = ',my_task &
                               ,' printing point = ',n &
                               ,' i and j = ',i,j
               write(nu_diag,*) ' coszen = ', &
                                  coszen(i,j)
               write(nu_diag,*) ' swvdr  swvdf = ', &
                                  swvdr(i,j),swvdf(i,j)
               write(nu_diag,*) ' swidr  swidf = ', &
                                  swidr(i,j),swidf(i,j)
               write(nu_diag,*) ' aice = ', &
                                  aice(i,j)
               write(nu_diag,*) ' hs = ', &
                                  hs(i,j)
               write(nu_diag,*) ' hp = ', &
                                  hp(i,j)
               write(nu_diag,*) ' fs = ', &
                                  fs(i,j)
               write(nu_diag,*) ' fi = ', &
                                  fi(i,j)
               write(nu_diag,*) ' fp = ', &
                                  fp(i,j)
               write(nu_diag,*) ' hi = ', &
                                  hi(i,j)
               write(nu_diag,*) ' srftyp = ', &
                                  srftyp(i,j)
               write(nu_diag,*) ' alvdr  alvdf = ', &
                                  alvdr(i,j),alvdf(i,j)
               write(nu_diag,*) ' alidr  alidf = ', &
                                  alidr(i,j),alidf(i,j)
               write(nu_diag,*) ' fswsfc fswint fswthru = ', &
                                  fswsfc(i,j),fswint(i,j),fswthru(i,j)
               swdn  = swvdr(i,j)+swvdf(i,j)+swidr(i,j)+swidf(i,j)
               swab  = fswsfc(i,j)+fswint(i,j)+fswthru(i,j)
               swalb = (1.-swab/(swdn+.0001))
               write(nu_diag,*) ' swdn swab swalb = ',swdn,swab,swalb
               do k = 1, nslyr               
                 write(nu_diag,*) ' snow layer k    = ', k, &
                                  ' rhosnw = ', &
                                    rhosnw(i,j,k), &
                                  ' rsnw = ', &
                                    rsnw(i,j,k)
               enddo
               do k = 1, nslyr               
                 write(nu_diag,*) ' snow layer k    = ', k, &
                                  ' Sswabs(k)       = ', Sswabs(i,j,k)
               enddo
               do k = 1, nilyr               
                 write(nu_diag,*) ' sea ice layer k = ', k, &
                                  ' Iswabs(k)       = ', Iswabs(i,j,k)
               enddo
               endif  ! coszen(i,j) > .01
            endif     ! my_task
         enddo        ! n for printing points
      endif           ! if print_points

      end subroutine shortwave_dEdd

!=======================================================================
!BOP
!
! !IROUTINE: compute_dEdd - evaluate Delta-Edd IOPs and compute solution
!
! !INTERFACE:
!
      subroutine compute_dEdd                              &
            (nx_block,ny_block,                            &
             icells_DE, indxi_DE, indxj_DE, fnidr, coszen, &
             swvdr,     swvdf,    swidr,    swidf, srftyp, &
             hs,        rhosnw,   rsnw,     hi,    hp,     &
             fi,       aero_mp,   alvdr,    alvdf,       &
                                  alidr,    alidf,       &
                                  fswsfc,   fswint,      &
                                  fswthru,  Sswabs,      &
                                  Iswabs)
!
! !DESCRIPTION:
!
! Evaluate snow/ice/ponded ice inherent optical properties (IOPs), and 
! then calculate the multiple scattering solution by calling solution_dEdd.
!
! !REVISION HISTORY:
!
! author:  Bruce P. Briegleb, NCAR 
! update:  8 February 2007
! update:  September 2008 added aerosols
!
! !USES:
!
      use ice_therm_vertical, only: heat_capacity
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), &
         intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells_DE             ! number of sea ice grid cells for surface type
 
      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi_DE, & ! compressed indices for sea ice cells for surface type
         indxj_DE
 
      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         fnidr   , & ! fraction of direct to total down flux in nir
         coszen  , & ! cosine solar zenith angle
         swvdr   , & ! shortwave down at surface, visible, direct  (W/m^2)
         swvdf   , & ! shortwave down at surface, visible, diffuse (W/m^2)
         swidr   , & ! shortwave down at surface, near IR, direct  (W/m^2)
         swidf       ! shortwave down at surface, near IR, diffuse (W/m^2)
 
      integer (kind=int_kind), dimension(nx_block,ny_block), &
         intent(in) :: &
         srftyp      ! surface type over ice: (0=air, 1=snow, 2=pond)
 
      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(in) :: &
         hs          ! snow thickness (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr), &
         intent(in) :: &
         rhosnw  , & ! snow density in snow layer (kg/m3)
         rsnw        ! snow grain radius in snow layer (m)

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(in) :: &
         hi      , & ! ice thickness (m)
         hp      , & ! pond depth (m)
         fi          ! snow/bare ice fractional coverage (0 to 1)
 
      real (kind=dbl_kind), dimension (nx_block,ny_block,4*n_aeromx), &
         intent(in) :: &
         aero_mp     ! aerosol mass path in kg/m2

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         alvdr   , & ! visible, direct, albedo (fraction) 
         alvdf   , & ! visible, diffuse, albedo (fraction) 
         alidr   , & ! near-ir, direct, albedo (fraction) 
         alidf   , & ! near-ir, diffuse, albedo (fraction) 
         fswsfc  , & ! SW absorbed at snow/bare ice/pondedi ice surface (W m-2)
         fswint  , & ! SW interior absorption (below surface, above ocean,W m-2)
         fswthru     ! SW through snow/bare ice/ponded ice into ocean (W m-2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr), &
         intent(inout) :: &
         Sswabs      ! SW absorbed in snow layer (W m-2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr), &
         intent(inout) :: &
         Iswabs      ! SW absorbed in ice layer (W m-2)
! 
!EOP
!-----------------------------------------------------------------------
!
! Set up optical property profiles, based on snow, sea ice and ponded 
! ice IOPs from:
!
! Briegleb, B. P., and B. Light (2007): A Delta-Eddington Multiple 
!    Scattering Parameterization for Solar Radiation in the Sea Ice 
!    Component of the Community Climate System Model, NCAR Technical 
!    Note  NCAR/TN-472+STR  February 2007
!
! Computes column Delta-Eddington radiation solution for specific
! surface type: either snow over sea ice, bare sea ice, or ponded sea ice.
!
! Divides solar spectrum into 3 intervals: 0.2-0.7, 0.7-1.19, and
! 1.19-5.0 micro-meters. The latter two are added (using an assumed
! partition of incident shortwave in the 0.7-5.0 micro-meter band between
! the 0.7-1.19 and 1.19-5.0 micro-meter band) to give the final output 
! of 0.2-0.7 visible and 0.7-5.0 near-infrared albedos and fluxes.
!
! Specifies vertical layer optical properties based on input snow depth,
! density and grain radius, along with ice and pond depths, then computes
! layer by layer Delta-Eddington reflectivity, transmissivity and combines
! layers (done by calling routine solution_dEdd). Finally, surface albedos
! and internal fluxes/flux divergences are evaluated.
!
!  Description of the level and layer index conventions. This is
!  for the standard case of one snow layer and four sea ice layers.
!
!  Please read the following; otherwise, there is 99.9% chance you
!  will be confused about indices at some point in time........ :)
!
!  CICE4.0 snow treatment has one snow layer above the sea ice. This 
!  snow layer has finite heat capacity, so that surface absorption must
!  be distinguished from internal. The Delta-Eddington solar radiation
!  thus adds extra surface scattering layers to both snow and sea ice.
!  Note that in the following, we assume a fixed vertical layer structure
!  for the radiation calculation. In other words, we always have the 
!  structure shown below for one snow and four sea ice layers, but for 
!  ponded ice the pond fills "snow" layer 1 over the sea ice, and for 
!  bare sea ice the top layers over sea ice are treated as transparent air.
!
!  SSL = surface scattering layer for either snow or sea ice
!  DL  = drained layer for sea ice immediately under sea ice SSL
!  INT = interior layers for sea ice below the drained layer.
!
!  Notice that the radiation level starts with 0 at the top. Thus,
!  the total number radiation layers is klev+1, where klev is the
!  sum of nslyr, the number of CCSM snow layers, and nilyr, the
!  number of CCSM sea ice layers, plus the sea ice SSL:
!  klev = 1 + nslyr + nilyr
!
!  For the standard case illustrated below, nslyr=1, nilyr=4,
!  and klev=6, with the number of layer interfaces klevp=klev+1.
!  Layer interfaces are the surfaces on which reflectivities,
!  transmissivities and fluxes are evaluated.
!
!  CCSM3 Sea Ice Model            Delta-Eddington Solar Radiation
!                                     Layers and Interfaces
!                             Layer Index             Interface Index
!    ---------------------            ---------------------  0
!                                  0  \\\   snow SSL    \\\
!       snow layer 1                  ---------------------  1
!                                  1    rest of snow layer
!    +++++++++++++++++++++            +++++++++++++++++++++  2
!                                  2  \\\ sea ice SSL   \\\
!      sea ice layer 1                ---------------------  3
!                                  3      sea ice  DL
!    ---------------------            ---------------------  4
!
!      sea ice layer 2             4      sea ice INT
!
!    ---------------------            ---------------------  5
!
!      sea ice layer 3             5      sea ice INT
!
!    ---------------------            ---------------------  6
!
!      sea ice layer 4             6      sea ice INT
!
!    ---------------------            ---------------------  7
!
! When snow lies over sea ice, the radiation absorbed in the
! snow SSL is used for surface heating, and that in the rest
! of the snow layer for its internal heating. For sea ice in
! this case, all of the radiant heat absorbed in both the
! sea ice SSL and the DL are used for sea ice layer 1 heating.
!
! When pond lies over sea ice, and for bare sea ice, all of the
! radiant heat absorbed within and above the sea ice SSL is used
! for surface heating, and that absorbed in the sea ice DL is
! used for sea ice layer 1 heating.
!
! Basically, vertical profiles of the layer extinction optical depth (tau), 
! single scattering albedo (w0) and asymmetry parameter (g) are required over
! the klev+1 layers, where klev+1 = 2 + nslyr + nilyr. All of the surface type
! information and snow/ice iop properties are evaulated in this routine, so
! the tau,w0,g profiles can be passed to solution_dEdd for multiple scattering
! evaluation. Snow, bare ice and ponded ice iops are contained in data arrays
! in this routine.
!
!-----------------------------------------------------------------------
!
! !LOCAL PARAMETERS
!
      integer (kind=int_kind) :: &
         i       , & ! longitude index
         j       , & ! latitude index
         k       , & ! level index
         ij      , & ! horizontal index, combines i and j loops
         ns      , & ! spectral index
         nr      , & ! index for grain radius tables
         ksa     , & ! index for snow internal absorption
         ki      , & ! index for sea ice internal absorption
         km      , & ! k starting index for snow, sea ice internal absorption
         kp      , & ! k+1 or k+2 index for snow, sea ice internal absorption
         ksrf    , & ! level index for surface absorption
         ksnow   , & ! level index for snow density and grain size
         kii        ! level starting index for sea ice (nslyr+1)

      integer (kind=int_kind), parameter :: & 
         klev    = nslyr + nilyr + 1   , & ! number of radiation layers - 1
         klevp   = klev  + 1               ! number of radiation interfaces - 1
                                           ! (0 layer is included also)
 
      integer (kind=int_kind), parameter :: & 
         nspint  = 3     , & ! number of solar spectral intervals
         nmbrad  = 32        ! number of snow grain radii in tables
 
      real (kind=dbl_kind), dimension(icells_DE) :: & 
         avdr    , & ! visible albedo, direct   (fraction)
         avdf    , & ! visible albedo, diffuse  (fraction)
         aidr    , & ! near-ir albedo, direct   (fraction)
         aidf        ! near-ir albedo, diffuse  (fraction)
 
      real (kind=dbl_kind), dimension(icells_DE) :: & 
         fsfc    , & ! shortwave absorbed at snow/bare ice/ponded ice surface (W m-2)
         fint    , & ! shortwave absorbed in interior (below surface but above ocean, W m-2)
         fthru       ! shortwave through snow/bare ice/ponded ice to ocean (W/m^2)

      real (kind=dbl_kind), dimension(icells_DE,nslyr) :: & 
         Sabs        ! shortwave absorbed in snow layer (W m-2)

      real (kind=dbl_kind), dimension(icells_DE,nilyr) :: & 
         Iabs        ! shortwave absorbed in ice layer (W m-2)

      real (kind=dbl_kind), dimension (icells_DE,nspint) :: &
         wghtns              ! spectral weights
 
      real (kind=dbl_kind), parameter :: & 
         cp67    = 0.67_dbl_kind   , & ! nir band weight parameter
         cp33    = 0.33_dbl_kind   , & ! nir band weight parameter
         cp78    = 0.78_dbl_kind   , & ! nir band weight parameter
         cp22    = 0.22_dbl_kind   , & ! nir band weight parameter
         cp01    = 0.01_dbl_kind       ! for ocean visible albedo
 
      real (kind=dbl_kind), dimension (0:klev,icells_DE) :: &
         tau     , & ! layer extinction optical depth
         w0      , & ! layer single scattering albedo
         g           ! layer asymmetry parameter
 
      ! following arrays are defined at model interfaces; 0 is the top of the
      ! layer above the sea ice; klevp is the sea ice/ocean interface.
      real (kind=dbl_kind), dimension (0:klevp,icells_DE) :: &
         trndir  , & ! solar beam down transmission from top
         trntdr  , & ! total transmission to direct beam for layers above
         trndif  , & ! diffuse transmission to diffuse beam for layers above
         rupdir  , & ! reflectivity to direct radiation for layers below
         rupdif  , & ! reflectivity to diffuse radiation for layers below
         rdndif      ! reflectivity to diffuse radiation for layers above
 
      real (kind=dbl_kind) :: &
         refk        ! interface k multiple scattering term
 
      real (kind=dbl_kind), dimension (0:klevp,icells_DE) :: &
         fdirup  , & ! up   flux at model interface due to direct beam at top surface
         fdirdn  , & ! down flux at model interface due to direct beam at top surface
         fdifup  , & ! up   flux at model interface due to diffuse beam at top surface
         fdifdn      ! down flux at model interface due to diffuse beam at top surface
 
      ! inherent optical property (iop) arrays for snow
      real (kind=dbl_kind), dimension (nspint) :: &
         Qs      , & ! Snow extinction efficiency
         ks      , & ! Snow extinction coefficient (/m)
         ws      , & ! Snow single scattering albedo
         gs          ! Snow asymmetry parameter

      real (kind=dbl_kind), dimension (nmbrad) :: &
         rsnw_tab    ! snow grain radius for each table entry (micro-meters)

      real (kind=dbl_kind), dimension (nspint,nmbrad) :: &
         Qs_tab  , & ! extinction efficiency for each snow grain radius
         ws_tab  , & ! single scatter albedo for each snow grain radius
         gs_tab      ! assymetry parameter   for each snow grain radius
      real (kind=dbl_kind) :: &
         delr    , & ! snow grain radius interpolation parameter
         rhoi    , & ! pure ice density (kg/m3)
         fr      , & ! snow grain adjustment factor 
         fr_max  , & ! snow grain adjustment factor max
         fr_min      ! snow grain adjustment factor min

      ! inherent optical property (iop) arrays for ice and ponded ice
      ! mn = specified mean (or base) value
      real (kind=dbl_kind), dimension (nspint) :: &
         ki_ssl_mn    , & ! Surface-scattering-layer ice extinction coefficient (/m)
         wi_ssl_mn    , & ! Surface-scattering-layer ice single scattering albedo
         gi_ssl_mn    , & ! Surface-scattering-layer ice asymmetry parameter
         ki_dl_mn     , & ! Drained-layer ice extinction coefficient (/m)
         wi_dl_mn     , & ! Drained-layer ice single scattering albedo
         gi_dl_mn     , & ! Drained-layer ice asymmetry parameter
         ki_int_mn    , & ! Interior-layer ice extinction coefficient (/m)
         wi_int_mn    , & ! Interior-layer ice single scattering albedo
         gi_int_mn    , & ! Interior-layer ice asymmetry parameter
         ki_p_ssl_mn  , & ! Ice under pond surface-scattering-layer extinction coefficient (/m)
         wi_p_ssl_mn  , & ! Ice under pond surface-scattering-layer single scattering albedo
         gi_p_ssl_mn  , & ! Ice under pond surface-scattering-layer asymmetry parameter
         ki_p_int_mn  , & ! Ice under pond interior extinction coefficient (/m)
         wi_p_int_mn  , & ! Ice under pond interior single scattering albedo
         gi_p_int_mn      ! Ice under pond interior asymmetry parameter

      ! actual used ice and ponded ice IOPs, allowing for tuning 
      ! modifications of the above "_mn" value
      real (kind=dbl_kind), dimension (nspint) :: &
         ki_ssl       , & ! Surface-scattering-layer ice extinction coefficient (/m)
         wi_ssl       , & ! Surface-scattering-layer ice single scattering albedo
         gi_ssl       , & ! Surface-scattering-layer ice asymmetry parameter
         ki_dl        , & ! Drained-layer ice extinction coefficient (/m)
         wi_dl        , & ! Drained-layer ice single scattering albedo
         gi_dl        , & ! Drained-layer ice asymmetry parameter
         ki_int       , & ! Interior-layer ice extinction coefficient (/m)
         wi_int       , & ! Interior-layer ice single scattering albedo
         gi_int       , & ! Interior-layer ice asymmetry parameter
         ki_p_ssl     , & ! Ice under pond srf scat layer extinction coefficient (/m)
         wi_p_ssl     , & ! Ice under pond srf scat layer single scattering albedo
         gi_p_ssl     , & ! Ice under pond srf scat layer asymmetry parameter
         ki_p_int     , & ! Ice under pond extinction coefficient (/m)
         wi_p_int     , & ! Ice under pond single scattering albedo
         gi_p_int         ! Ice under pond asymmetry parameter

      real (kind=dbl_kind) :: &
         hi_ssl       , & ! sea ice surface scattering layer thickness (m)
         hs_ssl       , & ! snow surface scattering layer thickness (m)
         dz           , & ! snow, sea ice or pond water layer thickness
         dz_ssl       , & ! snow or sea ice surface scattering layer thickness
         fs               ! scaling factor to reduce (nilyr<4) or increase (nilyr>4) DL
                          ! extinction coefficient to maintain DL optical depth constant
                          ! with changing number of sea ice layers, to approximately 
                          ! conserve computed albedo for constant physical depth of sea
                          ! ice when the number of sea ice layers vary
      real (kind=dbl_kind) :: &
         kalg         , & ! algae absorption coefficient for 0.5 m thick layer
         sig          , & ! scattering coefficient for tuning
         kabs         , & ! absorption coefficient for tuning
         sigp             ! modified scattering coefficient for tuning

      ! inherent optical property (iop) arrays for pond water and underlying ocean
      real (kind=dbl_kind), dimension (nspint) :: &
         kw           , & ! Pond water extinction coefficient (/m)
         ww           , & ! Pond water single scattering albedo
         gw               ! Pond water asymmetry parameter
      real (kind=dbl_kind), dimension (icells_DE) :: &
         albodr       , & ! spectral ocean albedo to direct rad
         albodf           ! spectral ocean albedo to diffuse rad
      
      ! tuning parameters
      real (kind=dbl_kind) :: &
         fp_ice       , & ! ice fraction of scat coeff for + stn dev in alb
         fm_ice       , & ! ice fraction of scat coeff for - stn dev in alb
         fp_pnd       , & ! ponded ice fraction of scat coeff for + stn dev in alb
         fm_pnd           ! ponded ice fraction of scat coeff for - stn dev in alb
 
      ! for melt pond transition to bare sea ice for small pond depths 
      real (kind=dbl_kind) :: &
         hpmin        , & ! minimum allowed melt pond depth (m)
         hp0          , & ! melt pond depth below which iops are weighted bare ice + pond (m)
         sig_i        , & ! ice scattering coefficient (/m)
         sig_p        , & ! pond scattering coefficient (/m)
         kext             ! weighted extinction coefficient (/m)

      ! aerosol optical properties from Mark Flanner, 26 June 2008
      ! order assumed: hydrophobic black carbon, hydrophilic black carbon,
      ! four dust aerosols by particle size range:
      ! dust1(.05-0.5 micron), dust2(0.5-1.25 micron),
      ! dust3(1.25-2.5 micron), dust4(2.5-5.0 micron)
      ! spectral bands same as snow/sea ice: (0.3-0.7 micron, 0.7-1.19 micron
      ! and 1.19-5.0 micron in wavelength)
      integer (kind=int_kind), parameter :: &
         nmbaer  = 6                ! number of aerosols

      integer (kind=int_kind) :: &
         nmbaer_actual,          &  ! actual number of aerosols used
         na                         ! aerosol index

      real (kind=dbl_kind) :: &
         kaer_tab(nspint,nmbaer), & ! aerosol mass extinction cross section (m2/kg)
         waer_tab(nspint,nmbaer), & ! aerosol single scatter albedo (fraction)
         gaer_tab(nspint,nmbaer), & ! aerosol asymmetry parameter (cos(theta))
         taer                   , & ! total aerosol extinction optical depth
         waer                   , & ! total aerosol single scatter albedo
         gaer                       ! total aerosol asymmetry parameter
#ifdef MODAL_AER
      !mgf++
      integer(kind=int_kind) :: idx_bcint_icerds
      integer(kind=int_kind) :: idx_bcint_nclrds
      integer(kind=int_kind) :: idx_bcext_nclrds
      real(kind=dbl_kind):: tmp_gs
      real(kind=dbl_kind):: tmp1
      !mgf--
#endif

      ! snow grain radii (micro-meters) for table
      data rsnw_tab/ &
          5._dbl_kind,    7._dbl_kind,   10._dbl_kind,   15._dbl_kind, &
         20._dbl_kind,   30._dbl_kind,   40._dbl_kind,   50._dbl_kind, &
         65._dbl_kind,   80._dbl_kind,  100._dbl_kind,  120._dbl_kind, &
        140._dbl_kind,  170._dbl_kind,  200._dbl_kind,  240._dbl_kind, &
        290._dbl_kind,  350._dbl_kind,  420._dbl_kind,  500._dbl_kind, &
        570._dbl_kind,  660._dbl_kind,  760._dbl_kind,  870._dbl_kind, &
       1000._dbl_kind, 1100._dbl_kind, 1250._dbl_kind, 1400._dbl_kind, &
       1600._dbl_kind, 1800._dbl_kind, 2000._dbl_kind, 2500._dbl_kind/

      ! snow extinction efficiency (unitless)
      data Qs_tab/ &
          2.131798_dbl_kind,  2.187756_dbl_kind,  2.267358_dbl_kind, &
          2.104499_dbl_kind,  2.148345_dbl_kind,  2.236078_dbl_kind, &
          2.081580_dbl_kind,  2.116885_dbl_kind,  2.175067_dbl_kind, &
          2.062595_dbl_kind,  2.088937_dbl_kind,  2.130242_dbl_kind, &
          2.051403_dbl_kind,  2.072422_dbl_kind,  2.106610_dbl_kind, &
          2.039223_dbl_kind,  2.055389_dbl_kind,  2.080586_dbl_kind, &
          2.032383_dbl_kind,  2.045751_dbl_kind,  2.066394_dbl_kind, &
          2.027920_dbl_kind,  2.039388_dbl_kind,  2.057224_dbl_kind, &
          2.023444_dbl_kind,  2.033137_dbl_kind,  2.048055_dbl_kind, &
          2.020412_dbl_kind,  2.028840_dbl_kind,  2.041874_dbl_kind, &
          2.017608_dbl_kind,  2.024863_dbl_kind,  2.036046_dbl_kind, &
          2.015592_dbl_kind,  2.022021_dbl_kind,  2.031954_dbl_kind, &
          2.014083_dbl_kind,  2.019887_dbl_kind,  2.028853_dbl_kind, &
          2.012368_dbl_kind,  2.017471_dbl_kind,  2.025353_dbl_kind, &
          2.011092_dbl_kind,  2.015675_dbl_kind,  2.022759_dbl_kind, &
          2.009837_dbl_kind,  2.013897_dbl_kind,  2.020168_dbl_kind, &
          2.008668_dbl_kind,  2.012252_dbl_kind,  2.017781_dbl_kind, &
          2.007627_dbl_kind,  2.010813_dbl_kind,  2.015678_dbl_kind, &
          2.006764_dbl_kind,  2.009577_dbl_kind,  2.013880_dbl_kind, &
          2.006037_dbl_kind,  2.008520_dbl_kind,  2.012382_dbl_kind, &
          2.005528_dbl_kind,  2.007807_dbl_kind,  2.011307_dbl_kind, &
          2.005025_dbl_kind,  2.007079_dbl_kind,  2.010280_dbl_kind, &
          2.004562_dbl_kind,  2.006440_dbl_kind,  2.009333_dbl_kind, &
          2.004155_dbl_kind,  2.005898_dbl_kind,  2.008523_dbl_kind, &
          2.003794_dbl_kind,  2.005379_dbl_kind,  2.007795_dbl_kind, &
          2.003555_dbl_kind,  2.005041_dbl_kind,  2.007329_dbl_kind, &
          2.003264_dbl_kind,  2.004624_dbl_kind,  2.006729_dbl_kind, &
          2.003037_dbl_kind,  2.004291_dbl_kind,  2.006230_dbl_kind, &
          2.002776_dbl_kind,  2.003929_dbl_kind,  2.005700_dbl_kind, &
          2.002590_dbl_kind,  2.003627_dbl_kind,  2.005276_dbl_kind, &
          2.002395_dbl_kind,  2.003391_dbl_kind,  2.004904_dbl_kind, &
          2.002071_dbl_kind,  2.002922_dbl_kind,  2.004241_dbl_kind/

      ! snow single scattering albedo (unitless)
      data ws_tab/ &
         0.9999994_dbl_kind,  0.9999673_dbl_kind,  0.9954589_dbl_kind, &
         0.9999992_dbl_kind,  0.9999547_dbl_kind,  0.9938576_dbl_kind, &
         0.9999990_dbl_kind,  0.9999382_dbl_kind,  0.9917989_dbl_kind, &
         0.9999985_dbl_kind,  0.9999123_dbl_kind,  0.9889724_dbl_kind, &
         0.9999979_dbl_kind,  0.9998844_dbl_kind,  0.9866190_dbl_kind, &
         0.9999970_dbl_kind,  0.9998317_dbl_kind,  0.9823021_dbl_kind, &
         0.9999960_dbl_kind,  0.9997800_dbl_kind,  0.9785269_dbl_kind, &
         0.9999951_dbl_kind,  0.9997288_dbl_kind,  0.9751601_dbl_kind, &
         0.9999936_dbl_kind,  0.9996531_dbl_kind,  0.9706974_dbl_kind, &
         0.9999922_dbl_kind,  0.9995783_dbl_kind,  0.9667577_dbl_kind, &
         0.9999903_dbl_kind,  0.9994798_dbl_kind,  0.9621007_dbl_kind, &
         0.9999885_dbl_kind,  0.9993825_dbl_kind,  0.9579541_dbl_kind, &
         0.9999866_dbl_kind,  0.9992862_dbl_kind,  0.9541924_dbl_kind, &
         0.9999838_dbl_kind,  0.9991434_dbl_kind,  0.9490959_dbl_kind, &
         0.9999810_dbl_kind,  0.9990025_dbl_kind,  0.9444940_dbl_kind, &
         0.9999772_dbl_kind,  0.9988171_dbl_kind,  0.9389141_dbl_kind, &
         0.9999726_dbl_kind,  0.9985890_dbl_kind,  0.9325819_dbl_kind, &
         0.9999670_dbl_kind,  0.9983199_dbl_kind,  0.9256405_dbl_kind, &
         0.9999605_dbl_kind,  0.9980117_dbl_kind,  0.9181533_dbl_kind, &
         0.9999530_dbl_kind,  0.9976663_dbl_kind,  0.9101540_dbl_kind, &
         0.9999465_dbl_kind,  0.9973693_dbl_kind,  0.9035031_dbl_kind, &
         0.9999382_dbl_kind,  0.9969939_dbl_kind,  0.8953134_dbl_kind, &
         0.9999289_dbl_kind,  0.9965848_dbl_kind,  0.8865789_dbl_kind, &
         0.9999188_dbl_kind,  0.9961434_dbl_kind,  0.8773350_dbl_kind, &
         0.9999068_dbl_kind,  0.9956323_dbl_kind,  0.8668233_dbl_kind, &
         0.9998975_dbl_kind,  0.9952464_dbl_kind,  0.8589990_dbl_kind, &
         0.9998837_dbl_kind,  0.9946782_dbl_kind,  0.8476493_dbl_kind, &
         0.9998699_dbl_kind,  0.9941218_dbl_kind,  0.8367318_dbl_kind, &
         0.9998515_dbl_kind,  0.9933966_dbl_kind,  0.8227881_dbl_kind, &
         0.9998332_dbl_kind,  0.9926888_dbl_kind,  0.8095131_dbl_kind, &
         0.9998148_dbl_kind,  0.9919968_dbl_kind,  0.7968620_dbl_kind, &
         0.9997691_dbl_kind,  0.9903277_dbl_kind,  0.7677887_dbl_kind/

      ! snow asymmetry parameter (unitless)
      data gs_tab / &
          0.859913_dbl_kind,  0.848003_dbl_kind,  0.824415_dbl_kind, &
          0.867130_dbl_kind,  0.858150_dbl_kind,  0.848445_dbl_kind, &
          0.873381_dbl_kind,  0.867221_dbl_kind,  0.861714_dbl_kind, &
          0.878368_dbl_kind,  0.874879_dbl_kind,  0.874036_dbl_kind, &
          0.881462_dbl_kind,  0.879661_dbl_kind,  0.881299_dbl_kind, &
          0.884361_dbl_kind,  0.883903_dbl_kind,  0.890184_dbl_kind, &
          0.885937_dbl_kind,  0.886256_dbl_kind,  0.895393_dbl_kind, &
          0.886931_dbl_kind,  0.887769_dbl_kind,  0.899072_dbl_kind, &
          0.887894_dbl_kind,  0.889255_dbl_kind,  0.903285_dbl_kind, &
          0.888515_dbl_kind,  0.890236_dbl_kind,  0.906588_dbl_kind, &
          0.889073_dbl_kind,  0.891127_dbl_kind,  0.910152_dbl_kind, &
          0.889452_dbl_kind,  0.891750_dbl_kind,  0.913100_dbl_kind, &
          0.889730_dbl_kind,  0.892213_dbl_kind,  0.915621_dbl_kind, &
          0.890026_dbl_kind,  0.892723_dbl_kind,  0.918831_dbl_kind, &
          0.890238_dbl_kind,  0.893099_dbl_kind,  0.921540_dbl_kind, &
          0.890441_dbl_kind,  0.893474_dbl_kind,  0.924581_dbl_kind, &
          0.890618_dbl_kind,  0.893816_dbl_kind,  0.927701_dbl_kind, &
          0.890762_dbl_kind,  0.894123_dbl_kind,  0.930737_dbl_kind, &
          0.890881_dbl_kind,  0.894397_dbl_kind,  0.933568_dbl_kind, &
          0.890975_dbl_kind,  0.894645_dbl_kind,  0.936148_dbl_kind, &
          0.891035_dbl_kind,  0.894822_dbl_kind,  0.937989_dbl_kind, &
          0.891097_dbl_kind,  0.895020_dbl_kind,  0.939949_dbl_kind, &
          0.891147_dbl_kind,  0.895212_dbl_kind,  0.941727_dbl_kind, &
          0.891189_dbl_kind,  0.895399_dbl_kind,  0.943339_dbl_kind, &
          0.891225_dbl_kind,  0.895601_dbl_kind,  0.944915_dbl_kind, &
          0.891248_dbl_kind,  0.895745_dbl_kind,  0.945950_dbl_kind, &
          0.891277_dbl_kind,  0.895951_dbl_kind,  0.947288_dbl_kind, &
          0.891299_dbl_kind,  0.896142_dbl_kind,  0.948438_dbl_kind, &
          0.891323_dbl_kind,  0.896388_dbl_kind,  0.949762_dbl_kind, &
          0.891340_dbl_kind,  0.896623_dbl_kind,  0.950916_dbl_kind, &
          0.891356_dbl_kind,  0.896851_dbl_kind,  0.951945_dbl_kind, &
          0.891386_dbl_kind,  0.897399_dbl_kind,  0.954156_dbl_kind/

      ! ice surface scattering layer (ssl) iops (units of k = /m)
      data ki_ssl_mn / 1000.1_dbl_kind, 1003.7_dbl_kind, 7042._dbl_kind/
      data wi_ssl_mn / .9999_dbl_kind,  .9963_dbl_kind,  .9088_dbl_kind/
      data gi_ssl_mn /  .94_dbl_kind,     .94_dbl_kind,    .94_dbl_kind/

      ! ice drained layer (dl) iops (units of k = /m)
      data ki_dl_mn  / 100.2_dbl_kind, 107.7_dbl_kind,  1309._dbl_kind /
      data wi_dl_mn  / .9980_dbl_kind,  .9287_dbl_kind, .0305_dbl_kind /
      data gi_dl_mn  / .94_dbl_kind,     .94_dbl_kind,    .94_dbl_kind /

      ! ice interior layer (int) iops (units of k = /m)
      data ki_int_mn /  20.2_dbl_kind,  27.7_dbl_kind,  1445._dbl_kind /
      data wi_int_mn / .9901_dbl_kind, .7223_dbl_kind,  .0277_dbl_kind /
      data gi_int_mn / .94_dbl_kind,    .94_dbl_kind,     .94_dbl_kind /

      ! ponded ice surface scattering layer (ssl) iops (units of k = /m)
      data ki_p_ssl_mn / 70.2_dbl_kind,  77.7_dbl_kind,  1309._dbl_kind/
      data wi_p_ssl_mn / .9972_dbl_kind, .9009_dbl_kind, .0305_dbl_kind/
      data gi_p_ssl_mn / .94_dbl_kind,   .94_dbl_kind,   .94_dbl_kind  /

      ! ponded ice interior layer (int) iops (units of k = /m)
      data ki_p_int_mn /  20.2_dbl_kind,  27.7_dbl_kind, 1445._dbl_kind/
      data wi_p_int_mn / .9901_dbl_kind, .7223_dbl_kind, .0277_dbl_kind/
      data gi_p_int_mn / .94_dbl_kind,   .94_dbl_kind,   .94_dbl_kind  /

      ! pond water iops (units of k = /m)
      data kw   /    0.20_dbl_kind,   12.0_dbl_kind,   729._dbl_kind /
      data ww   /    0.00_dbl_kind,   0.00_dbl_kind,   0.00_dbl_kind /
      data gw   /    0.00_dbl_kind,   0.00_dbl_kind,   0.00_dbl_kind /

      ! snow data
      data hs_ssl / 0.040_dbl_kind / ! snow surface scattering layer thickness (m)
      data rhoi   /917.0_dbl_kind /  ! snow mass density (kg/m3)
      data fr_max / 1.00_dbl_kind /  ! snow grain adjustment factor max
      data fr_min / 0.80_dbl_kind /  ! snow grain adjustment factor min

      ! ice data
      data hi_ssl / 0.050_dbl_kind / ! sea ice surface scattering layer thickness (m)
!     data kalg   / 0.60_dbl_kind /  ! for 0.5 m path of 75 mg Chl a / m2
!     turn off algae absorption for now - DAB
      data kalg   / 0.00_dbl_kind /  ! for 0.5 m path of 75 mg Chl a / m2

      ! ice and pond scat coeff fractional change for +- one-sigma in albedo
      data fp_ice / 0.15_dbl_kind /
      data fm_ice / 0.15_dbl_kind /
      data fp_pnd / 2.00_dbl_kind /
      data fm_pnd / 0.50_dbl_kind /

      ! ice to pond parameters
      data hpmin  / .005_dbl_kind / ! minimum allowable pond depth (m)
      data hp0    / .200_dbl_kind / ! pond depth below which transition to bare sea ice

      ! aerosol optical properties   -> band  |
      !                                       v aerosol
      ! for combined dust category, let's use category 4 properties
      data kaer_tab/ &
          11580.61872,   5535.41835,   2793.79690, &
          25798.96479,  11536.03871,   4688.24207, &
            196.49772,    204.14078,    214.42287, &
           2665.85867,   2256.71027,    820.36024, &
            840.78295,   1028.24656,   1163.03298, &
            387.51211,    414.68808,    450.29814/
      data waer_tab/ &
              0.29003,      0.17349,      0.06613, &
              0.51731,      0.41609,      0.21324, &
              0.84467,      0.94216,      0.95666, &
              0.97764,      0.99402,      0.98552, &
              0.94146,      0.98527,      0.99093, &
              0.90034,      0.96543,      0.97678/
      data gaer_tab/ &
              0.35445,      0.19838,      0.08857, &
              0.52581,      0.32384,      0.14970, &
              0.83162,      0.78306,      0.74375, &
              0.68861,      0.70836,      0.54171, &
              0.70239,      0.66115,      0.71983, &
              0.78734,      0.73580,      0.64411/
!     data kaer_tab/ &
!         11580.61872,   5535.41835,   2793.79690, &
!         25798.96479,  11536.03871,   4688.24207, &
!          2665.85867,   2256.71027,    820.36024, &
!           840.78295,   1028.24656,   1163.03298, &
!           387.51211,    414.68808,    450.29814, &
!           196.49772,    204.14078,    214.42287  /
!     data waer_tab/ &
!             0.29003,      0.17349,      0.06613, &
!             0.51731,      0.41609,      0.21324, &
!             0.97764,      0.99402,      0.98552, &
!             0.94146,      0.98527,      0.99093, &
!             0.90034,      0.96543,      0.97678, &
!             0.84467,      0.94216,      0.95666  /
!     data gaer_tab/ &
!             0.35445,      0.19838,      0.08857, &
!             0.52581,      0.32384,      0.14970, &
!             0.68861,      0.70836,      0.54171, &
!             0.70239,      0.66115,      0.71983, &
!             0.78734,      0.73580,      0.64411, &
!             0.83162,      0.78306,      0.74375  /

!-----------------------------------------------------------------------
! Initialize and tune bare ice/ponded ice iops
 
      rnilyr = real(nilyr,kind=dbl_kind)
      rnslyr = real(nslyr,kind=dbl_kind)

      ! initialize albedos and fluxes to 0
      do ij = 1, icells_DE
         avdr(ij)   = c0
         avdf(ij)   = c0
         aidr(ij)   = c0
         aidf(ij)   = c0
         fsfc(ij)   = c0
         fint(ij)   = c0
         fthru(ij)  = c0
      enddo                ! ij
      Sabs(:,:) = c0
      Iabs(:,:) = c0
 
      ! spectral weights; weights 2 (0.7-1.19 micro-meters) and 3 (1.19-5.0 micro-meters) 
      ! are chosen based on 1D calculations using ratio of direct to total near-infrared
      ! solar (0.7-5.0 micro-meter) which indicates clear/cloudy conditions: more cloud,
      ! the less 1.19-5.0 relative to the 0.7-1.19 micro-meter due to cloud absorption.
      do ij = 1, icells_DE
        i = indxi_DE(ij)
        j = indxj_DE(ij)
        wghtns(ij,1) = c1
        wghtns(ij,2) = cp67 + (cp78-cp67)*(c1-fnidr(i,j))
        wghtns(ij,3) = cp33 + (cp22-cp33)*(c1-fnidr(i,j))
      enddo

      ! adjust sea ice iops with tuning parameters; tune only the
      ! scattering coefficient by factors of R_ice, R_pnd, where
      ! R values of +1 correspond approximately to +1 sigma changes in albedo, and
      ! R values of -1 correspond approximately to -1 sigma changes in albedo
      ! Note: the albedo change becomes non-linear for R values > +1 or < -1
      if( R_ice >= c0 ) then
        do ns = 1, nspint
          sigp       = ki_ssl_mn(ns)*wi_ssl_mn(ns)*(c1+fp_ice*R_ice)
          ki_ssl(ns) = sigp+ki_ssl_mn(ns)*(c1-wi_ssl_mn(ns))
          wi_ssl(ns) = sigp/ki_ssl(ns)
          gi_ssl(ns) = gi_ssl_mn(ns)

          sigp       = ki_dl_mn(ns)*wi_dl_mn(ns)*(c1+fp_ice*R_ice)
          ki_dl(ns)  = sigp+ki_dl_mn(ns)*(c1-wi_dl_mn(ns))
          wi_dl(ns)  = sigp/ki_dl(ns)
          gi_dl(ns)  = gi_dl_mn(ns)

          sigp       = ki_int_mn(ns)*wi_int_mn(ns)*(c1+fp_ice*R_ice)
          ki_int(ns) = sigp+ki_int_mn(ns)*(c1-wi_int_mn(ns))
          wi_int(ns) = sigp/ki_int(ns)
          gi_int(ns) = gi_int_mn(ns)
        enddo
      else !if( R_ice < c0 ) then
        do ns = 1, nspint
          sigp       = ki_ssl_mn(ns)*wi_ssl_mn(ns)*(c1+fm_ice*R_ice)
          sigp       = max(sigp, c0)
          ki_ssl(ns) = sigp+ki_ssl_mn(ns)*(c1-wi_ssl_mn(ns))
          wi_ssl(ns) = sigp/ki_ssl(ns)
          gi_ssl(ns) = gi_ssl_mn(ns)

          sigp       = ki_dl_mn(ns)*wi_dl_mn(ns)*(c1+fm_ice*R_ice)
          sigp       = max(sigp, c0)
          ki_dl(ns)  = sigp+ki_dl_mn(ns)*(c1-wi_dl_mn(ns))
          wi_dl(ns)  = sigp/ki_dl(ns)
          gi_dl(ns)  = gi_dl_mn(ns)

          sigp       = ki_int_mn(ns)*wi_int_mn(ns)*(c1+fm_ice*R_ice)
          sigp       = max(sigp, c0)
          ki_int(ns) = sigp+ki_int_mn(ns)*(c1-wi_int_mn(ns))
          wi_int(ns) = sigp/ki_int(ns)
          gi_int(ns) = gi_int_mn(ns)
        enddo
      endif          ! adjust ice iops

      ! adjust ponded ice iops with tuning parameters
      if( R_pnd >= c0 ) then
        do ns = 1, nspint
          sigp         = ki_p_ssl_mn(ns)*wi_p_ssl_mn(ns)*(c1+fp_pnd*R_pnd)
          ki_p_ssl(ns) = sigp+ki_p_ssl_mn(ns)*(c1-wi_p_ssl_mn(ns))
          wi_p_ssl(ns) = sigp/ki_p_ssl(ns)
          gi_p_ssl(ns) = gi_p_ssl_mn(ns)

          sigp         = ki_p_int_mn(ns)*wi_p_int_mn(ns)*(c1+fp_pnd*R_pnd)
          ki_p_int(ns) = sigp+ki_p_int_mn(ns)*(c1-wi_p_int_mn(ns))
          wi_p_int(ns) = sigp/ki_p_int(ns)
          gi_p_int(ns) = gi_p_int_mn(ns)
        enddo
      else !if( R_pnd < c0 ) then
        do ns = 1, nspint
          sigp         = ki_p_ssl_mn(ns)*wi_p_ssl_mn(ns)*(c1+fm_pnd*R_pnd)
          sigp         = max(sigp, c0)
          ki_p_ssl(ns) = sigp+ki_p_ssl_mn(ns)*(c1-wi_p_ssl_mn(ns))
          wi_p_ssl(ns) = sigp/ki_p_ssl(ns)
          gi_p_ssl(ns) = gi_p_ssl_mn(ns)

          sigp         = ki_p_int_mn(ns)*wi_p_int_mn(ns)*(c1+fm_pnd*R_pnd)
          sigp         = max(sigp, c0)
          ki_p_int(ns) = sigp+ki_p_int_mn(ns)*(c1-wi_p_int_mn(ns))
          wi_p_int(ns) = sigp/ki_p_int(ns)
          gi_p_int(ns) = gi_p_int_mn(ns)
        enddo
      endif            ! adjust ponded ice iops

!-----------------------------------------------------------------------
 
      ! begin spectral loop
      do ns = 1, nspint
 
        if (tflag) call t_startf('cice_computededd1_'//tstr)
        ! set optical properties of air/snow/pond overlying sea ice
        do ij = 1, icells_DE
          i = indxi_DE(ij)
          j = indxj_DE(ij)
          ! air
          if( srftyp(i,j) == 0 ) then
            do k=0,nslyr 
              tau(k,ij) = c0
              w0(k,ij)  = c0
              g(k,ij)   = c0
            enddo
          ! snow
          else if( srftyp(i,j) == 1 ) then
            dz_ssl = hs_ssl
            dz     = hs(i,j)/rnslyr
            ! for small enough snow thickness, ssl thickness half of snow layer
            dz_ssl = min(dz_ssl, dz/c2)
            ! find snow grain adjustment factor, dependent upon clear/overcast sky
            ! estimate. comparisons with SNICAR show better agreement with DE when
            ! this factor is included (clear sky near 1 and overcast near 0.8 give
            ! best agreement).
            fr     = fr_max*fnidr(i,j) + fr_min*(c1-fnidr(i,j))
            ! interpolate snow iops using input snow grain radius,
            ! snow density and tabular data
            ksnow = 1
            do k=0,nslyr
              ! use top rsnw, rhosnw for snow ssl and rest of top layer
              if( k > 1 ) ksnow = k
              ! find snow iops using input snow density and snow grain radius:
              if( fr*rsnw(i,j,ksnow) < rsnw_tab(1) ) then
                Qs(ns) = Qs_tab(ns,1)
                ws(ns) = ws_tab(ns,1)
                gs(ns) = gs_tab(ns,1)
              else if( fr*rsnw(i,j,ksnow) >= rsnw_tab(nmbrad) ) then
                Qs(ns) = Qs_tab(ns,nmbrad)
                ws(ns) = ws_tab(ns,nmbrad)
                gs(ns) = gs_tab(ns,nmbrad)
              else
                ! linear interpolation in rsnw
                do nr=2,nmbrad
                  if( rsnw_tab(nr-1) <= fr*rsnw(i,j,ksnow) .and. &
                      fr*rsnw(i,j,ksnow) < rsnw_tab(nr)) then
                        delr = (fr*rsnw(i,j,ksnow) - rsnw_tab(nr-1)) / &
                               (rsnw_tab(nr) - rsnw_tab(nr-1))
                        Qs(ns) = Qs_tab(ns,nr-1)*(c1-delr) + &
                                 Qs_tab(ns,nr)*delr
                        ws(ns) = ws_tab(ns,nr-1)*(c1-delr) + &
                                 ws_tab(ns,nr)*delr
                        gs(ns) = gs_tab(ns,nr-1)*(c1-delr) + &
                                 gs_tab(ns,nr)*delr
                  endif
                enddo       ! nr
              endif
              ks(ns) = Qs(ns)*((rhosnw(i,j,ksnow)/rhoi)*3._dbl_kind / &
                       (4._dbl_kind*fr*rsnw(i,j,ksnow)*1.0e-6_dbl_kind))
              if( k == 0 ) then
                tau(k,ij) = ks(ns)*dz_ssl
              else if( k == 1 ) then
                tau(k,ij) = ks(ns)*(dz-dz_ssl)
              else !if( k >= 2 ) then
                tau(k,ij) = ks(ns)*dz
              endif
              w0(k,ij)  = ws(ns)
              g(k,ij)   = gs(ns)
              ! aerosol in snow
              nmbaer_actual = min(n_aero,nmbaer)

#ifdef MODAL_AER
              !mgf++
              ! snow grain size (um)
              tmp_gs = fr*rsnw(i,j,ksnow)
                
              ! get grain size index:
              ! works for 25 < snw_rds < 1625 um:
              if (tmp_gs < 125) then
                 tmp1 = tmp_gs/50
                 idx_bcint_icerds = nint(tmp1)
              elseif (tmp_gs < 175) then
                 idx_bcint_icerds = 2
              else
                 tmp1 = (tmp_gs/250)+2
                 idx_bcint_icerds = nint(tmp1)
              endif

              ! Set index corresponding to BC effective radius.  Here,
              ! asssume constant BC effective radius of 100nm
              ! (corresponding to index 2)
              idx_bcint_nclrds = 2
              idx_bcext_nclrds = 2

              ! check bounds:
              if (idx_bcint_icerds < 1)  idx_bcint_icerds = 1
              if (idx_bcint_icerds > 8)  idx_bcint_icerds = 8
              if (idx_bcint_nclrds < 1)  idx_bcint_nclrds = 1
              if (idx_bcint_nclrds > 10) idx_bcint_nclrds = 10
              if (idx_bcext_nclrds < 1)  idx_bcext_nclrds = 1
              if (idx_bcext_nclrds > 10) idx_bcext_nclrds = 10

              ! print ice radius index:
              !write(6,*) "MGFICE2: ice index= ", idx_bcint_icerds
              !write(6,*) "MGFICE2: ext_cff_mss_bc = ", kaer_bc_tab(1,1), kaer_bc_tab(2,1), kaer_bc_tab(1,2), kaer_bc_tab(3,1),kaer_bc_tab(3,10) 
              !write(6,*) "MGFICE2: ss_alb_bc = ", waer_bc_tab(1,1), waer_bc_tab(2,1), waer_bc_tab(1,2), waer_bc_tab(3,1), waer_bc_tab(3,10)
              !write(6,*) "MGFICE2: asm_prm_bc = ", gaer_bc_tab(1,1), gaer_bc_tab(2,1), gaer_bc_tab(1,2), gaer_bc_tab(3,1), gaer_bc_tab(3,10)
              !write(6,*) "MGFICE2: bcenh = ", bcenh(1,1,1), bcenh(1,2,1), bcenh(1,1,2), bcenh(2,1,1), bcenh(3,10,1), bcenh(3,1,8), bcenh(3,10,8)
              !mgf--
#endif
              if( k == 0 ) then  ! snow SSL
                taer = c0
                waer = c0
                gaer = c0
                do na=1,4*nmbaer_actual,4
#ifdef MODAL_AER
                   ! mgf++
                   if (na==1) then
                      !  interstitial BC
                      taer = taer + &
                           aero_mp(i,j,na)*kaer_bc_tab(ns,idx_bcext_nclrds)
                      waer = waer + &
                           aero_mp(i,j,na)*kaer_bc_tab(ns,idx_bcext_nclrds)* &
                           waer_bc_tab(ns,idx_bcext_nclrds)
                      gaer = gaer + &
                           aero_mp(i,j,na)*kaer_bc_tab(ns,idx_bcext_nclrds)* &
                           waer_bc_tab(ns,idx_bcext_nclrds)*gaer_bc_tab(ns,idx_bcext_nclrds)
                   elseif (na==5) then
                      ! within-ice BC
                      taer = taer + &
                           aero_mp(i,j,na)*kaer_bc_tab(ns,idx_bcint_nclrds)*bcenh(ns,idx_bcint_nclrds,idx_bcint_icerds)
                      waer = waer + &
                           aero_mp(i,j,na)*kaer_bc_tab(ns,idx_bcint_nclrds)* &
                           waer_bc_tab(ns,idx_bcint_nclrds)
                      gaer = gaer + &
                           aero_mp(i,j,na)*kaer_bc_tab(ns,idx_bcint_nclrds)* &
                           waer_bc_tab(ns,idx_bcint_nclrds)*gaer_bc_tab(ns,idx_bcint_nclrds)
                      
                   else
                      ! other species (dust)
                      taer = taer + &
                           aero_mp(i,j,na)*kaer_tab(ns,(1+(na-1)/4))
                      waer = waer + &
                           aero_mp(i,j,na)*kaer_tab(ns,(1+(na-1)/4))* &
                           waer_tab(ns,(1+(na-1)/4))
                      gaer = gaer + &
                           aero_mp(i,j,na)*kaer_tab(ns,(1+(na-1)/4))* &
                           waer_tab(ns,(1+(na-1)/4))*gaer_tab(ns,(1+(na-1)/4))
                   endif
                   !mgf--
#else
                  taer = taer + &
                       aero_mp(i,j,na)*kaer_tab(ns,(1+(na-1)/4))
                  waer = waer + &
                       aero_mp(i,j,na)*kaer_tab(ns,(1+(na-1)/4))* &
                         waer_tab(ns,(1+(na-1)/4))
                  gaer = gaer + &
                       aero_mp(i,j,na)*kaer_tab(ns,(1+(na-1)/4))* &
                         waer_tab(ns,(1+(na-1)/4))*gaer_tab(ns,(1+(na-1)/4))
#endif

                enddo       ! na
                gaer = gaer/(waer+puny)
                waer = waer/(taer+puny)
              else if ( k > 0 ) then  ! snow below SSL
                taer = c0
                waer = c0
                gaer = c0
                do na=1,4*nmbaer_actual,4
#ifdef MODAL_AER
                   !mgf++
                   if (na==1) then
                      ! interstitial BC
                      taer = taer + &
                           (aero_mp(i,j,na+1)/rnslyr)*kaer_bc_tab(ns,idx_bcext_nclrds)
                      waer = waer + &
                           (aero_mp(i,j,na+1)/rnslyr)*kaer_bc_tab(ns,idx_bcext_nclrds)* &
                           waer_bc_tab(ns,idx_bcext_nclrds)
                      gaer = gaer + &
                           (aero_mp(i,j,na+1)/rnslyr)*kaer_bc_tab(ns,idx_bcext_nclrds)* &
                           waer_bc_tab(ns,idx_bcext_nclrds)*gaer_bc_tab(ns,idx_bcext_nclrds)
                   elseif (na==5) then
                      ! within-ice BC
                      taer = taer + &
                           (aero_mp(i,j,na+1)/rnslyr)*kaer_bc_tab(ns,idx_bcint_nclrds)*bcenh(ns,idx_bcint_nclrds,idx_bcint_icerds)
                      waer = waer + &
                           (aero_mp(i,j,na+1)/rnslyr)*kaer_bc_tab(ns,idx_bcint_nclrds)* &
                           waer_bc_tab(ns,idx_bcint_nclrds)
                      gaer = gaer + &
                           (aero_mp(i,j,na+1)/rnslyr)*kaer_bc_tab(ns,idx_bcint_nclrds)* &
                           waer_bc_tab(ns,idx_bcint_nclrds)*gaer_bc_tab(ns,idx_bcint_nclrds)
                      
                   else
                      ! other species (dust)
                      taer = taer + &
                           (aero_mp(i,j,na+1)/rnslyr)*kaer_tab(ns,(1+(na-1)/4))
                      waer = waer + &
                           (aero_mp(i,j,na+1)/rnslyr)*kaer_tab(ns,(1+(na-1)/4))* &
                           waer_tab(ns,(1+(na-1)/4))
                      gaer = gaer + &
                           (aero_mp(i,j,na+1)/rnslyr)*kaer_tab(ns,(1+(na-1)/4))* &
                           waer_tab(ns,(1+(na-1)/4))*gaer_tab(ns,(1+(na-1)/4))
                   endif
                   !mgf--

#else
                  taer = taer + &
                       (aero_mp(i,j,na+1)/rnslyr)*kaer_tab(ns,(1+(na-1)/4))
                  waer = waer + &
                       (aero_mp(i,j,na+1)/rnslyr)*kaer_tab(ns,(1+(na-1)/4))* &
                         waer_tab(ns,(1+(na-1)/4))
                  gaer = gaer + &
                       (aero_mp(i,j,na+1)/rnslyr)*kaer_tab(ns,(1+(na-1)/4))* &
                         waer_tab(ns,(1+(na-1)/4))*gaer_tab(ns,(1+(na-1)/4))
#endif
                enddo       ! na
                gaer = gaer/(waer+puny)
                waer = waer/(taer+puny)
              endif
              g(k,ij)   = (g(k,ij)*w0(k,ij)*tau(k,ij) + gaer*waer*taer) / &
                                  (w0(k,ij)*tau(k,ij) + waer*taer)
              w0(k,ij)  = (w0(k,ij)*tau(k,ij) + waer*taer) / &
                                   (tau(k,ij) + taer)
              tau(k,ij) = tau(k,ij) + taer
            enddo       ! k
          ! pond
          else !if( srftyp(i,j) == 2 ) then
            ! pond water layers evenly spaced
            dz = hp(i,j)/(rnslyr+c1)
            do k=0,nslyr
              tau(k,ij) = kw(ns)*dz
              w0(k,ij)  = ww(ns)
              g(k,ij)   = gw(ns)
              ! no aerosol in pond
            enddo       ! k
          endif        ! srftyp
        enddo         ! ij ... optical properties above sea ice set
        if (tflag) call t_stopf('cice_computededd1_'//tstr)

        if (tflag) call t_startf('cice_computededd2_'//tstr)
        ! set optical properties of sea ice
        kii = nslyr + 1
        do ij = 1, icells_DE
          i = indxi_DE(ij)
          j = indxj_DE(ij)
          dz_ssl = hi_ssl
          dz     = hi(i,j)/rnilyr
          ! empirical reduction in sea ice ssl thickness for ice thinner than 1.5m;
          ! factor of 30 selected to give best albedo comparison with limited observations
          if( hi(i,j) < 1.5_dbl_kind ) dz_ssl = hi(i,j)/30._dbl_kind
          ! set sea ice ssl thickness to half top layer if sea ice thin enough
          dz_ssl = min(dz_ssl, dz/c2)
          ! bare or snow-covered sea ice layers
          if( srftyp(i,j) <= 1 ) then
#ifdef MODAL_AER
             !mgf++ 
             ! we are in sea-ice, so use largest snow grain size for BC optical property lookup
             idx_bcint_icerds = 8
             
             ! Set index corresponding to BC effective radius.  Here,
             ! asssume constant BC effective radius of 100nm
             ! (corresponding to index 2)
             idx_bcint_nclrds = 2
             idx_bcext_nclrds = 2
             !mgf--
#endif
              ! ssl
              k = kii
                tau(k,ij) = ki_ssl(ns)*dz_ssl
                w0(k,ij)  = wi_ssl(ns)
                g(k,ij)   = gi_ssl(ns)
              ! dl
              k = kii + 1
                ! scale dz for dl relative to 4 even-layer-thickness 1.5m case
                fs = rnilyr/c4
                tau(k,ij) = ki_dl(ns)*(dz-dz_ssl)*fs
                w0(k,ij)  = wi_dl(ns)
                g(k,ij)   = gi_dl(ns)
              ! int above lowest layer
              if (kii+2 <= klev-1) then
              do k = kii+2, klev-1
                tau(k,ij) = ki_int(ns)*dz
                w0(k,ij)  = wi_int(ns)
                g(k,ij)   = gi_int(ns)
              enddo
              endif
              ! lowest layer
              k = klev
                ! add algae to lowest sea ice layer, visible only:
                kabs       = ki_int(ns)*(c1-wi_int(ns))
                if( ns == 1 ) then
                  ! total layer absorption optical depth fixed at value
                  ! of kalg*0.50m, independent of actual layer thickness
                  kabs = kabs + kalg*(0.50_dbl_kind/dz)
                endif
                sig        = ki_int(ns)*wi_int(ns)
                tau(k,ij) = (kabs+sig)*dz
                w0(k,ij)  = (sig/(sig+kabs))
                g(k,ij)   = gi_int(ns)
              ! aerosol in sea ice
              nmbaer_actual = min(n_aero,nmbaer)
              do k = kii, klev
                if( k == kii ) then  ! sea ice SSL
                  taer = c0
                  waer = c0
                  gaer = c0
                  do na=1,4*nmbaer_actual,4
#ifdef MODAL_AER
                     !mgf++
                     if (na==1) then
                        ! interstitial BC
                        taer = taer + &
                             aero_mp(i,j,na+2)*kaer_bc_tab(ns,idx_bcext_nclrds)
                        waer = waer + &
                             aero_mp(i,j,na+2)*kaer_bc_tab(ns,idx_bcext_nclrds)* &
                             waer_bc_tab(ns,idx_bcext_nclrds)
                        gaer = gaer + &
                             aero_mp(i,j,na+2)*kaer_bc_tab(ns,idx_bcext_nclrds)* &
                             waer_bc_tab(ns,idx_bcext_nclrds)*gaer_bc_tab(ns,idx_bcext_nclrds)
                     elseif (na==5) then
                        ! within-ice BC
                        taer = taer + &
                             aero_mp(i,j,na+2)*kaer_bc_tab(ns,idx_bcint_nclrds)*bcenh(ns,idx_bcint_nclrds,idx_bcint_icerds)
                        waer = waer + &
                             aero_mp(i,j,na+2)*kaer_bc_tab(ns,idx_bcint_nclrds)* &
                             waer_bc_tab(ns,idx_bcint_nclrds)
                        gaer = gaer + &
                             aero_mp(i,j,na+2)*kaer_bc_tab(ns,idx_bcint_nclrds)* &
                             waer_bc_tab(ns,idx_bcint_nclrds)*gaer_bc_tab(ns,idx_bcint_nclrds)
                        
                     else
                        ! other species (dust)
                        taer = taer + &
                             aero_mp(i,j,na+2)*kaer_tab(ns,(1+(na-1)/4))
                        waer = waer + &
                             aero_mp(i,j,na+2)*kaer_tab(ns,(1+(na-1)/4))* &
                             waer_tab(ns,(1+(na-1)/4))
                        gaer = gaer + &
                             aero_mp(i,j,na+2)*kaer_tab(ns,(1+(na-1)/4))* &
                             waer_tab(ns,(1+(na-1)/4))*gaer_tab(ns,(1+(na-1)/4))
                     endif
#else
                    taer = taer + &
                         aero_mp(i,j,na+2)*kaer_tab(ns,(1+(na-1)/4))
                    waer = waer + &
                         aero_mp(i,j,na+2)*kaer_tab(ns,(1+(na-1)/4))* &
                           waer_tab(ns,(1+(na-1)/4))
                    gaer = gaer + &
                         aero_mp(i,j,na+2)*kaer_tab(ns,(1+(na-1)/4))* &
                           waer_tab(ns,(1+(na-1)/4))*gaer_tab(ns,(1+(na-1)/4))
#endif
                  enddo       ! na
                  gaer = gaer/(waer+puny)
                  waer = waer/(taer+puny)
                else if ( k > kii ) then  ! sea ice below SSL
                  taer = c0
                  waer = c0
                  gaer = c0
                  do na=1,4*nmbaer_actual,4
#ifdef MODAL_AER
                     !mgf++
                     if (na==1) then
                        ! interstitial BC
                        taer = taer + &
                             (aero_mp(i,j,na+3)/rnilyr)*kaer_bc_tab(ns,idx_bcext_nclrds)
                        waer = waer + &
                             (aero_mp(i,j,na+3)/rnilyr)*kaer_bc_tab(ns,idx_bcext_nclrds)* &
                             waer_bc_tab(ns,idx_bcext_nclrds)
                        gaer = gaer + &
                             (aero_mp(i,j,na+3)/rnilyr)*kaer_bc_tab(ns,idx_bcext_nclrds)* &
                             waer_bc_tab(ns,idx_bcext_nclrds)*gaer_bc_tab(ns,idx_bcext_nclrds)
                     elseif (na==5) then
                        ! within-ice BC
                        taer = taer + &
                             (aero_mp(i,j,na+3)/rnilyr)*kaer_bc_tab(ns,idx_bcint_nclrds)*bcenh(ns,idx_bcint_nclrds,idx_bcint_icerds)
                        waer = waer + &
                             (aero_mp(i,j,na+3)/rnilyr)*kaer_bc_tab(ns,idx_bcint_nclrds)* &
                             waer_bc_tab(ns,idx_bcint_nclrds)
                        gaer = gaer + &
                             (aero_mp(i,j,na+3)/rnilyr)*kaer_bc_tab(ns,idx_bcint_nclrds)* &
                             waer_bc_tab(ns,idx_bcint_nclrds)*gaer_bc_tab(ns,idx_bcint_nclrds)
                        
                     else
                        ! other species (dust)
                        taer = taer + &
                             (aero_mp(i,j,na+3)/rnilyr)*kaer_tab(ns,(1+(na-1)/4))
                        waer = waer + &
                             (aero_mp(i,j,na+3)/rnilyr)*kaer_tab(ns,(1+(na-1)/4))* &
                             waer_tab(ns,(1+(na-1)/4))
                        gaer = gaer + &
                             (aero_mp(i,j,na+3)/rnilyr)*kaer_tab(ns,(1+(na-1)/4))* &
                             waer_tab(ns,(1+(na-1)/4))*gaer_tab(ns,(1+(na-1)/4))
                     endif
                     !mgf--
                      
#else
                    taer = taer + &
                         (aero_mp(i,j,na+3)/rnilyr)*kaer_tab(ns,(1+(na-1)/4))
                    waer = waer + &
                         (aero_mp(i,j,na+3)/rnilyr)*kaer_tab(ns,(1+(na-1)/4))* &
                           waer_tab(ns,(1+(na-1)/4))
                    gaer = gaer + &
                         (aero_mp(i,j,na+3)/rnilyr)*kaer_tab(ns,(1+(na-1)/4))* &
                           waer_tab(ns,(1+(na-1)/4))*gaer_tab(ns,(1+(na-1)/4))
#endif
                  enddo       ! na
                  gaer = gaer/(waer+puny)
                  waer = waer/(taer+puny)
                endif
                g(k,ij)   = (g(k,ij)*w0(k,ij)*tau(k,ij) + gaer*waer*taer) / &
                                    (w0(k,ij)*tau(k,ij) + waer*taer)
                w0(k,ij)  = (w0(k,ij)*tau(k,ij) + waer*taer) / &
                                     (tau(k,ij) + taer)
                tau(k,ij) = tau(k,ij) + taer
              enddo ! k
          ! sea ice layers under ponds
          else !if( srftyp(i,j) == 2 ) then
              k = kii
                tau(k,ij) = ki_p_ssl(ns)*dz_ssl
                w0(k,ij)  = wi_p_ssl(ns)
                g(k,ij)   = gi_p_ssl(ns)
              k = kii + 1
                tau(k,ij) = ki_p_int(ns)*(dz-dz_ssl)
                w0(k,ij)  = wi_p_int(ns)
                g(k,ij)   = gi_p_int(ns)
              if (kii+2 <= klev) then
              do k = kii+2, klev
                tau(k,ij) = ki_p_int(ns)*dz
                w0(k,ij)  = wi_p_int(ns)
                g(k,ij)   = gi_p_int(ns)
              enddo       ! k
              endif
            ! adjust pond iops if pond depth within specified range
            if( hpmin <= hp(i,j) .and. hp(i,j) <= hp0 ) then
              k = kii
                  sig_i      = ki_ssl(ns)*wi_ssl(ns)
                  sig_p      = ki_p_ssl(ns)*wi_p_ssl(ns)
                  sig        = sig_i + (sig_p-sig_i)*(hp(i,j)/hp0)
                  kext       = sig + ki_p_ssl(ns)*(c1-wi_p_ssl(ns))
                  tau(k,ij) = kext*dz_ssl
                  w0(k,ij) = sig/kext
                  g(k,ij)  = gi_p_int(ns)
              k = kii + 1
                  ! scale dz for dl relative to 4 even-layer-thickness 1.5m case
                  fs = rnilyr/c4
                  sig_i      = ki_dl(ns)*wi_dl(ns)*fs
                  sig_p      = ki_p_int(ns)*wi_p_int(ns)
                  sig        = sig_i + (sig_p-sig_i)*(hp(i,j)/hp0)
                  kext       = sig + ki_p_int(ns)*(c1-wi_p_int(ns))
                  tau(k,ij) = kext*(dz-dz_ssl)
                  w0(k,ij) = sig/kext
                  g(k,ij)  = gi_p_int(ns)
              if (kii+2 <= klev) then
              do k = kii+2, klev
                  sig_i      = ki_int(ns)*wi_int(ns)
                  sig_p      = ki_p_int(ns)*wi_p_int(ns)
                  sig        = sig_i + (sig_p-sig_i)*(hp(i,j)/hp0)
                  kext       = sig + ki_p_int(ns)*(c1-wi_p_int(ns))
                  tau(k,ij) = kext*dz
                  w0(k,ij) = sig/kext
                  g(k,ij)  = gi_p_int(ns)
              enddo       ! k
              endif
            endif        ! small pond depth transition to bare sea ice
          endif         ! srftyp  
        enddo          ! ij ... optical properties of sea ice set
        if (tflag) call t_stopf('cice_computededd2_'//tstr)
 
        ! set reflectivities for ocean underlying sea ice
        if(ns == 1) then
          do ij = 1, icells_DE
            i = indxi_DE(ij)
            j = indxj_DE(ij)
            albodr(ij) = cp01
            albodf(ij) = cp01
          enddo       ! ij
        else !if(ns >= 2) then
          do ij = 1, icells_DE
            i = indxi_DE(ij)
            j = indxj_DE(ij)
            albodr(ij) = c0
            albodf(ij) = c0
          enddo       ! ij
        endif
 
        ! layer input properties now completely specified: tau, w0, g,
        ! albodr, albodf; now compute the Delta-Eddington solution 
        ! reflectivities and transmissivities for each layer; then,
        ! combine the layers going downwards accounting for multiple
        ! scattering between layers, and finally start from the 
        ! underlying ocean and combine successive layers upwards to
        ! the surface; see comments in solution_dEdd for more details.
 
        if (icells_DE > 0) then
           if (tflag) call t_startf('cice_computededdsoln_'//tstr)
           call solution_dEdd                                    &
               (nx_block, ny_block,                              &
                icells_DE, indxi_DE,  indxj_DE,  coszen, srftyp, &
                tau,       w0,        g,         albodr, albodf, &
                trndir,    trntdr,    trndif,    rupdir, rupdif, &
                rdndif)
           if (tflag) call t_stopf('cice_computededdsoln_'//tstr)
        endif

        ! the interface reflectivities and transmissivities required
        ! to evaluate interface fluxes are returned from solution_dEdd;
        ! now compute up and down fluxes for each interface, using the 
        ! combined layer properties at each interface:
        !
        !              layers       interface
        !
        !       ---------------------  k
        !                 k
        !       --------------------- 

        if (tflag) call t_startf('cice_computededd4_'//tstr)
        do ij = 1, icells_DE
          i = indxi_DE(ij)
          j = indxj_DE(ij)
          do k=0,klevp 
            ! interface scattering
            refk          = c1/(c1 - rdndif(k,ij)*rupdif(k,ij))
            ! dir tran ref from below times interface scattering, plus diff
            ! tran and ref from below times interface scattering
            fdirup(k,ij) = (trndir(k,ij)*rupdir(k,ij) + &
                            (trntdr(k,ij)-trndir(k,ij))  &
                            *rupdif(k,ij))*refk
            ! dir tran plus total diff trans times interface scattering plus
            ! dir tran with up dir ref and down dif ref times interface scattering 
            fdirdn(k,ij) = trndir(k,ij) + (trntdr(k,ij) &
                          - trndir(k,ij) + trndir(k,ij)  &
                          *rupdir(k,ij)*rdndif(k,ij))*refk
            ! diffuse tran ref from below times interface scattering
            fdifup(k,ij) = trndif(k,ij)*rupdif(k,ij)*refk
            ! diffuse tran times interface scattering
            fdifdn(k,ij) = trndif(k,ij)*refk
          enddo       ! k 
        enddo       ! ij
        if (tflag) call t_stopf('cice_computededd4_'//tstr)
 
        ! calculate final surface albedos and fluxes-
        ! all absorbed flux above ksrf is included in surface absorption

        if( ns == 1) then      ! visible

          if (tflag) call t_startf('cice_computededd5_'//tstr)
          do ij = 1, icells_DE
            i = indxi_DE(ij)
            j = indxj_DE(ij)
            avdr(ij)  = rupdir(0,ij)
            avdf(ij)  = rupdif(0,ij)

            ! use srftyp to determine interface index of surface absorption
            if( srftyp(i,j) == 1 ) then
              ! snow covered sea ice
              ksrf = 1
            else
              ! bare sea ice or ponded ice
              ksrf = nslyr + 2 
            endif
 
            fsfc(ij)  = fsfc(ij) + &
              ((fdirdn(0,ij)-fdirup(0,ij))*swvdr(i,j) + &
               (fdifdn(0,ij)-fdifup(0,ij))*swvdf(i,j)) - &
              ((fdirdn(ksrf,ij)-fdirup(ksrf,ij))*swvdr(i,j) + &
               (fdifdn(ksrf,ij)-fdifup(ksrf,ij))*swvdf(i,j))
 
            fint(ij)  = fint(ij) + &
              ((fdirdn(ksrf,ij)-fdirup(ksrf,ij))*swvdr(i,j) + &
               (fdifdn(ksrf,ij)-fdifup(ksrf,ij))*swvdf(i,j)) - &
              ((fdirdn(klevp,ij)-fdirup(klevp,ij))*swvdr(i,j) + &
               (fdifdn(klevp,ij)-fdifup(klevp,ij))*swvdf(i,j))
 
            fthru(ij)  = fthru(ij) + &
               (fdirdn(klevp,ij)-fdirup(klevp,ij))*swvdr(i,j) + &
               (fdifdn(klevp,ij)-fdifup(klevp,ij))*swvdf(i,j)

            ! if snow covered ice, set snow internal absorption; else, Sabs=0
            if( srftyp(i,j) == 1 ) then
              ksa = 0
              do k=1,nslyr
                ! skip snow SSL, since SSL absorption included in the surface
                ! absorption fsfc above
                km  = k
                kp  = km + 1
                ksa = ksa + 1
                Sabs(ij,ksa) = Sabs(ij,ksa) + &
                ((fdirdn(km,ij)-fdirup(km,ij))*swvdr(i,j) + &
                 (fdifdn(km,ij)-fdifup(km,ij))*swvdf(i,j)) - &
                ((fdirdn(kp,ij)-fdirup(kp,ij))*swvdr(i,j) + &
                 (fdifdn(kp,ij)-fdifup(kp,ij))*swvdf(i,j))
              enddo       ! k
            endif

            ! complex indexing to insure proper absorptions for sea ice
            ki = 0
            do k=nslyr+2,nslyr+1+nilyr
              ! for bare ice, DL absorption for sea ice layer 1
              km = k  
              kp = km + 1
              ! modify for top sea ice layer for snow over sea ice
              if( srftyp(i,j) == 1 ) then
                ! must add SSL and DL absorption for sea ice layer 1
                if( k == nslyr+2 ) then
                  km = k  - 1
                  kp = km + 2
                endif
              endif
              ki = ki + 1
              Iabs(ij,ki) = Iabs(ij,ki) + &
              ((fdirdn(km,ij)-fdirup(km,ij))*swvdr(i,j) + &
               (fdifdn(km,ij)-fdifup(km,ij))*swvdf(i,j)) - &
              ((fdirdn(kp,ij)-fdirup(kp,ij))*swvdr(i,j) + &
               (fdifdn(kp,ij)-fdifup(kp,ij))*swvdf(i,j))
            enddo       ! k
          enddo        ! ij
          if (tflag) call t_stopf('cice_computededd5_'//tstr)

        else !if(ns > 1) then  ! near IR

          if (tflag) call t_startf('cice_computededd6_'//tstr)
          do ij = 1, icells_DE
            i = indxi_DE(ij)
            j = indxj_DE(ij)

            ! let fr1 = alb_1*swd*wght1 and fr2 = alb_2*swd*wght2 be the ns=2,3
            ! reflected fluxes respectively, where alb_1, alb_2 are the band
            ! albedos, swd = nir incident shortwave flux, and wght1, wght2 are
            ! the 2,3 band weights. thus, the total reflected flux is:
            ! fr = fr1 + fr2 = alb_1*swd*wght1 + alb_2*swd*wght2  hence, the
            ! 2,3 nir band albedo is alb = fr/swd = alb_1*wght1 + alb_2*wght2

            aidr(ij)   = aidr(ij) + rupdir(0,ij)*wghtns(ij,ns)
            aidf(ij)   = aidf(ij) + rupdif(0,ij)*wghtns(ij,ns)

            ! use srftyp to determine interface index of surface absorption
            if( srftyp(i,j) == 1 ) then
              ! snow covered sea ice
              ksrf = 1
            else
              ! bare sea ice or ponded ice
              ksrf = nslyr + 2
            endif
 
            fsfc(ij)  = fsfc(ij) + &
            ( ((fdirdn(0,ij)-fdirup(0,ij))*swidr(i,j) + &
               (fdifdn(0,ij)-fdifup(0,ij))*swidf(i,j)) - &
              ((fdirdn(ksrf,ij)-fdirup(ksrf,ij))*swidr(i,j) + &
               (fdifdn(ksrf,ij)-fdifup(ksrf,ij))*swidf(i,j)) ) &
               *wghtns(ij,ns)
 
            fint(ij)  = fint(ij) + &
            ( ((fdirdn(ksrf,ij)-fdirup(ksrf,ij))*swidr(i,j) + &
               (fdifdn(ksrf,ij)-fdifup(ksrf,ij))*swidf(i,j)) - &
              ((fdirdn(klevp,ij)-fdirup(klevp,ij))*swidr(i,j) + &
               (fdifdn(klevp,ij)-fdifup(klevp,ij))*swidf(i,j)) ) &
               *wghtns(ij,ns)
 
            fthru(ij) = fthru(ij) + &
              ((fdirdn(klevp,ij)-fdirup(klevp,ij))*swidr(i,j) + &
               (fdifdn(klevp,ij)-fdifup(klevp,ij))*swidf(i,j)) &
               *wghtns(ij,ns)

            ! if snow covered ice, set snow internal absorption; else, Sabs=0
            if( srftyp(i,j) == 1 ) then
              ksa = 0
              do k=1,nslyr
                ! skip snow SSL, since SSL absorption included in the surface
                ! absorption fsfc above
                km  = k
                kp  = km + 1
                ksa = ksa + 1
                Sabs(ij,ksa) = Sabs(ij,ksa) + &
                ( ((fdirdn(km,ij)-fdirup(km,ij))*swidr(i,j) + &
                   (fdifdn(km,ij)-fdifup(km,ij))*swidf(i,j)) - &
                  ((fdirdn(kp,ij)-fdirup(kp,ij))*swidr(i,j) + &
                   (fdifdn(kp,ij)-fdifup(kp,ij))*swidf(i,j)) ) & 
                   *wghtns(ij,ns)
              enddo       ! k
            endif

            ! complex indexing to insure proper absorptions for sea ice
            ki = 0
            do k=nslyr+2,nslyr+1+nilyr
              ! for bare ice, DL absorption for sea ice layer 1
              km = k  
              kp = km + 1
              ! modify for top sea ice layer for snow over sea ice
              if( srftyp(i,j) == 1 ) then
                ! must add SSL and DL absorption for sea ice layer 1
                if( k == nslyr+2 ) then
                  km = k  - 1
                  kp = km + 2
                endif
              endif
              ki = ki + 1
              Iabs(ij,ki) = Iabs(ij,ki) + &
              ( ((fdirdn(km,ij)-fdirup(km,ij))*swidr(i,j) + &
                 (fdifdn(km,ij)-fdifup(km,ij))*swidf(i,j)) - &
                ((fdirdn(kp,ij)-fdirup(kp,ij))*swidr(i,j) + &
                 (fdifdn(kp,ij)-fdifup(kp,ij))*swidf(i,j)) ) &
                 *wghtns(ij,ns)
            enddo       ! k
          enddo        ! ij
          if (tflag) call t_stopf('cice_computededd6_'//tstr)

        endif        ! ns = 1, ns > 1

      enddo         ! end spectral loop  ns

      ! accumulate fluxes over bare sea ice
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells_DE
         i = indxi_DE(ij)
         j = indxj_DE(ij)
         alvdr(i,j)   = avdr(ij)
         alvdf(i,j)   = avdf(ij)
         alidr(i,j)   = aidr(ij)
         alidf(i,j)   = aidf(ij)
         fswsfc(i,j)  = fswsfc(i,j)  + fsfc(ij) *fi(i,j)
         fswint(i,j)  = fswint(i,j)  + fint(ij) *fi(i,j)
         fswthru(i,j) = fswthru(i,j) + fthru(ij)*fi(i,j)
      enddo                     ! ij

      do k = 1, nslyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells_DE
            i = indxi_DE(ij)
            j = indxj_DE(ij)
            Sswabs(i,j,k) = Sswabs(i,j,k) + Sabs(ij,k)*fi(i,j)
         enddo                  ! ij
      enddo                     ! k

      do k = 1, nilyr
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells_DE
            i = indxi_DE(ij)
            j = indxj_DE(ij)
            Iswabs(i,j,k) = Iswabs(i,j,k) + Iabs(ij,k)*fi(i,j)
         enddo                  ! ij
      enddo                     ! k

      !----------------------------------------------------------------
      ! if ice has zero heat capacity, no SW can be absorbed 
      ! in the ice/snow interior, so add to surface absorption.
      ! Note: nilyr = nslyr = 1 for this case
      !----------------------------------------------------------------

      if (.not. heat_capacity) then

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
        do ij = 1, icells_DE
           i = indxi_DE(ij)
           j = indxj_DE(ij)

           ! SW absorbed at snow/ice surface
           fswsfc(i,j) = fswsfc(i,j) + Iswabs(i,j,1) + Sswabs(i,j,1)

           ! SW absorbed in ice interior
           fswint(i,j)   = c0
           Iswabs(i,j,1) = c0
           Sswabs(i,j,1) = c0

        enddo                     ! ij

      endif                       ! heat_capacity

      end subroutine compute_dEdd

!=======================================================================
!BOP
!
! !IROUTINE: solution_dEdd - evaluate solution for Delta-Edddington solar
!
! !INTERFACE:
!
      subroutine solution_dEdd                                 &
            (nx_block, ny_block,                               &
             icells_DE,  indxi_DE,  indxj_DE,  coszen, srftyp, &
             tau,        w0,        g,         albodr, albodf, &
             trndir,     trntdr,    trndif,    rupdir, rupdif, &
             rdndif)
!
! !DESCRIPTION:
!
! Given input vertical profiles of optical properties, evaluate the 
! monochromatic Delta-Eddington solution.
!
! !REVISION HISTORY:
!
! author:  Bruce P. Briegleb, NCAR
! updated: 8 February 2007
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:

      integer (kind=int_kind), &
         intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells_DE             ! number of sea ice grid cells for surface type

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi_DE, & ! compressed indices for sea ice cells for surface type
         indxj_DE

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         coszen      ! cosine solar zenith angle

      integer (kind=int_kind), dimension(nx_block,ny_block), &
         intent(in) :: &
         srftyp      ! surface type over ice: (0=air, 1=snow, 2=pond)
 
      integer (kind=int_kind), parameter :: &
         klev    = nslyr + nilyr + 1 , &  ! number of radiation layers - 1
         klevp   = klev  + 1              ! number of radiation interfaces - 1

      real (kind=dbl_kind), dimension(0:klev,icells_DE), &
         intent(in) :: &
         tau     , & ! layer extinction optical depth
         w0      , & ! layer single scattering albedo
         g           ! layer asymmetry parameter
 
      real (kind=dbl_kind), dimension(icells_DE), &
         intent(in) :: &
         albodr  , & ! ocean albedo to direct rad
         albodf      ! ocean albedo to diffuse rad
 
      ! following arrays are defined at model interfaces; 0 is the top of the
      ! layer above the sea ice; klevp is the sea ice/ocean interface.
      real (kind=dbl_kind), dimension (0:klevp,icells_DE), &
         intent(out) :: &
         trndir  , & ! solar beam down transmission from top
         trntdr  , & ! total transmission to direct beam for layers above
         trndif  , & ! diffuse transmission to diffuse beam for layers above
         rupdir  , & ! reflectivity to direct radiation for layers below
         rupdif  , & ! reflectivity to diffuse radiation for layers below
         rdndif      ! reflectivity to diffuse radiation for layers above

!
!EOP
!-----------------------------------------------------------------------
!
! Delta-Eddington solution for snow/air/pond over sea ice
!
! Generic solution for a snow/air/pond input column of klev+1 layers,
! with srftyp determining at what interface fresnel refraction occurs.
!
! Computes layer reflectivities and transmissivities, from the top down
! to the lowest interface using the Delta-Eddington solutions for each
! layer; combines layers from top down to lowest interface, and from the
! lowest interface (underlying ocean) up to the top of the column.
!
! Note that layer diffuse reflectivity and transmissivity are computed
! by integrating the direct over several gaussian angles. This is
! because the diffuse reflectivity expression sometimes is negative,
! but the direct reflectivity is always well-behaved. We assume isotropic
! radiation in the upward and downward hemispheres for this integration.
!
! Assumes monochromatic (spectrally uniform) properties across a band
! for the input optical parameters.
!
! If total transmission of the direct beam to the interface above a particular 
! layer is less than trmin, then no further Delta-Eddington solutions are
! evaluated for layers below.
!
! The following describes how refraction is handled in the calculation.
!
! First, we assume that radiation is refracted when entering either
! sea ice at the base of the surface scattering layer, or water (i.e. melt
! pond); we assume that radiation does not refract when entering snow, nor 
! upon entering sea ice from a melt pond, nor upon entering the underlying 
! ocean from sea ice.
!
! To handle refraction, we define a "fresnel" layer, which physically
! is of neglible thickness and is non-absorbing, which can be combined to 
! any sea ice layer or top of melt pond. The fresnel layer accounts for 
! refraction of direct beam and associated reflection and transmission for
! solar radiation. A fresnel layer is combined with the top of a melt pond 
! or to the surface scattering layer of sea ice if no melt pond lies over it. 
!
! Some caution must be exercised for the fresnel layer, because any layer
! to which it is combined is no longer a homogeneous layer, as are all other
! individual layers. For all other layers for example, the direct and diffuse
! reflectivities/transmissivities (R/T) are the same for radiation above or
! below the layer. This is the meaning of homogeneous! But for the fresnel
! layer this is not so. Thus, the R/T for this layer must be distinguished
! for radiation above from that from radiation below. For generality, we
! treat all layers to be combined as inhomogeneous.
!
!-----------------------------------------------------------------------

! Local

      integer (kind=int_kind) :: & 
         kfrsnl      ! radiation interface index for fresnel layer
 
      ! following variables are defined for each layer; 0 refers to the top
      ! layer. In general we must distinguish directions above and below in 
      ! the diffuse reflectivity and transmissivity, as layers are not assumed
      ! to be homogeneous (apart from the single layer Delta-Edd solutions); 
      ! the direct is always from above.
      real (kind=dbl_kind), dimension (0:klev,icells_DE) :: &
         rdir    , & ! layer reflectivity to direct radiation
         rdif_a  , & ! layer reflectivity to diffuse radiation from above
         rdif_b  , & ! layer reflectivity to diffuse radiation from below
         tdir    , & ! layer transmission to direct radiation (solar beam + diffuse)
         tdif_a  , & ! layer transmission to diffuse radiation from above
         tdif_b  , & ! layer transmission to diffuse radiation from below
         trnlay      ! solar beam transm for layer (direct beam only)

      integer (kind=int_kind) :: & 
         i       , & ! longitude index
         j       , & ! latitude index
         ij      , & ! longitude/latitude index
         k           ! level index
 
      real (kind=dbl_kind), parameter :: &
         trmin = 0.001_dbl_kind   ! minimum total transmission allowed
      ! total transmission is that due to the direct beam; i.e. it includes
      ! both the directly transmitted solar beam and the diffuse downwards
      ! transmitted radiation resulting from scattering out of the direct beam 
      real (kind=dbl_kind) :: &
         tautot   , & ! layer optical depth
         wtot     , & ! layer single scattering albedo
         gtot     , & ! layer asymmetry parameter
         ftot     , & ! layer forward scattering fraction
         ts       , & ! layer scaled extinction optical depth
         ws       , & ! layer scaled single scattering albedo
         gs       , & ! layer scaled asymmetry parameter
         rintfc   , & ! reflection (multiple) at an interface
         refkp1   , & ! interface multiple scattering for k+1
         refkm1   , & ! interface multiple scattering for k-1
         tdrrdir  , & ! direct tran times layer direct ref 
         tdndif       ! total down diffuse = tot tran - direct tran
 
      ! perpendicular and parallel relative to plane of incidence and scattering
      real (kind=dbl_kind) :: &
         R1       , & ! perpendicular polarization reflection amplitude
         R2       , & ! parallel polarization reflection amplitude
         T1       , & ! perpendicular polarization transmission amplitude
         T2       , & ! parallel polarization transmission amplitude
         Rf_dir_a , & ! fresnel reflection to direct radiation
         Tf_dir_a , & ! fresnel transmission to direct radiation
         Rf_dif_a , & ! fresnel reflection to diff radiation from above
         Rf_dif_b , & ! fresnel reflection to diff radiation from below
         Tf_dif_a , & ! fresnel transmission to diff radiation from above
         Tf_dif_b     ! fresnel transmission to diff radiation from below
 
      ! refractive index for sea ice, water; pre-computed, band-independent,
      ! diffuse fresnel reflectivities
      real (kind=dbl_kind), parameter :: & 
         refindx = 1.310_dbl_kind  , & ! refractive index of sea ice (used for water also)
         cp063   = 0.063_dbl_kind  , & ! diffuse fresnel reflectivity from above
         cp455   = 0.455_dbl_kind      ! diffuse fresnel reflectivity from below
 
      real (kind=dbl_kind) :: &
         mu0      , & ! cosine solar zenith angle incident
         mu0n         ! cosine solar zenith angle in medium
 
      real (kind=dbl_kind) :: &
         alpha    , & ! term in direct reflectivity and transmissivity
         gamma    , & ! term in direct reflectivity and transmissivity
         el       , & ! term in alpha,gamma,n,u
         taus     , & ! scaled extinction optical depth
         omgs     , & ! scaled single particle scattering albedo
         asys     , & ! scaled asymmetry parameter
         u        , & ! term in diffuse reflectivity and transmissivity
         n        , & ! term in diffuse reflectivity and transmissivity
         lm       , & ! temporary for el
         mu       , & ! cosine solar zenith for either snow or water
         ne           ! temporary for n
 
      real (kind=dbl_kind) :: &
         w        , & ! dummy argument for statement function
         uu       , & ! dummy argument for statement function
         gg       , & ! dummy argument for statement function
         e        , & ! dummy argument for statement function
         f        , & ! dummy argument for statement function
         t        , & ! dummy argument for statement function
         et           ! dummy argument for statement function
 
      real (kind=dbl_kind) :: &
         alp      , & ! temporary for alpha
         gam      , & ! temporary for gamma
         ue       , & ! temporary for u
         arg      , & ! exponential argument
         extins   , & ! extinction
         amg      , & ! alp - gam
         apg          ! alp + gam
 
      integer (kind=int_kind), parameter :: &
         ngmax = 8    ! number of gaussian angles in hemisphere
 
      real (kind=dbl_kind), dimension (ngmax) :: &
         gauspt   , & ! gaussian angles (radians)
         gauswt       ! gaussian weights

      data gauspt/ &
                 .9894009_dbl_kind,  .9445750_dbl_kind, &
                 .8656312_dbl_kind,  .7554044_dbl_kind, &
                 .6178762_dbl_kind,  .4580168_dbl_kind, &
                 .2816036_dbl_kind,  .0950125_dbl_kind  /
      data gauswt/ &     
                 .0271525_dbl_kind,  .0622535_dbl_kind, &
                 .0951585_dbl_kind,  .1246290_dbl_kind, &
                 .1495960_dbl_kind,  .1691565_dbl_kind, &
                 .1826034_dbl_kind,  .1894506_dbl_kind  /
  
      integer (kind=int_kind) :: &
         ng           ! gaussian integration index
 
      real (kind=dbl_kind) :: &
         gwt      , & ! gaussian weight
         swt      , & ! sum of weights
         trn      , & ! layer transmission
         rdr      , & ! rdir for gaussian integration
         tdr      , & ! tdir for gaussian integration
         smr      , & ! accumulator for rdif gaussian integration
         smt          ! accumulator for tdif gaussian integration
 
      ! Delta-Eddington solution expressions
      alpha(w,uu,gg,e) = p75*w*uu*((c1 + gg*(c1-w))/(c1 - e*e*uu*uu))
      gamma(w,uu,gg,e) = p5*w*((c1 + c3*gg*(c1-w)*uu*uu) &
                        / (c1-e*e*uu*uu))
      n(uu,et)         = ((uu+c1)*(uu+c1)/et ) - ((uu-c1)*(uu-c1)*et)
      u(w,gg,e)        = c1p5*(c1 - w*gg)/e
      el(w,gg)         = sqrt(c3*(c1-w)*(c1 - w*gg))
      taus(w,f,t)      = (c1 - w*f)*t
      omgs(w,f)        = (c1 - f)*w/(c1 - w*f)
      asys(gg,f)       = (gg - f)/(c1 - f)
 
!-----------------------------------------------------------------------

      ! initialize all output to 0
      do ij = 1, icells_DE
      do k = 0, klevp 
         trndir(k,ij) = c0
         trntdr(k,ij) = c0
         trndif(k,ij) = c0
         rupdir(k,ij) = c0
         rupdif(k,ij) = c0
         rdndif(k,ij) = c0
      enddo
 
      ! initialize all layer apparent optical properties to 0
      do k = 0, klev
         rdir(k,ij)   = c0
         rdif_a(k,ij) = c0
         rdif_b(k,ij) = c0
         tdir(k,ij)   = c0
         tdif_a(k,ij) = c0
         tdif_b(k,ij) = c0
         trnlay(k,ij) = c0
      enddo

      ! initialize top interface of top layer 
        trndir(0,ij) =   c1
        trntdr(0,ij) =   c1
        trndif(0,ij) =   c1
        rdndif(0,ij) =   c0
      enddo       ! ij
 
      ! proceed down one layer at a time; if the total transmission to
      ! the interface just above a given layer is less than trmin, then no
      ! Delta-Eddington computation for that layer is done.

      if (tflag) call t_startf('cice_solEdd1x_'//tstr)
      do ij = 1, icells_DE 
        ! if (tflag) call t_startf('cice_solEdd1ij_'//tstr)
        i = indxi_DE(ij)
        j = indxj_DE(ij)
        ! if (tflag) call t_stopf('cice_solEdd1ij_'//tstr)

        ! begin main level loop
        ! if (tflag) call t_startf('cice_solEdd1kl_'//tstr)
        do k=0,klev
 
           ! initialize current layer properties to zero; only if total
           ! transmission to the top interface of the current layer exceeds the
           ! minimum, will these values be computed below:
           if ( k > 0 ) then  
              ! if (tflag) call t_startf('cice_solEdd1k0_'//tstr)
              ! Calculate the solar beam transmission, total transmission, and
              ! reflectivity for diffuse radiation from below at interface k, 
              ! the top of the current layer k:
              !
              !              layers       interface
              !         
              !       ---------------------  k-1 
              !                k-1
              !       ---------------------  k
              !                 k
              !       ---------------------  
              trndir(k,ij) = trndir(k-1,ij)*trnlay(k-1,ij)
              refkm1        = c1/(c1 - rdndif(k-1,ij)*rdif_a(k-1,ij))
              tdrrdir       = trndir(k-1,ij)*rdir(k-1,ij)
              tdndif        = trntdr(k-1,ij) - trndir(k-1,ij)
              trntdr(k,ij) = trndir(k-1,ij)*tdir(k-1,ij) + &
               (tdndif + tdrrdir*rdndif(k-1,ij))*refkm1*tdif_a(k-1,ij)
              rdndif(k,ij) = rdif_b(k-1,ij) + &
                (tdif_b(k-1,ij)*rdndif(k-1,ij)*refkm1*tdif_a(k-1,ij))
              trndif(k,ij) = trndif(k-1,ij)*refkm1*tdif_a(k-1,ij)
              ! if (tflag) call t_stopf('cice_solEdd1k0_'//tstr)
           endif       ! k > 0  
 
           ! compute next layer Delta-eddington solution only if total transmission
           ! of radiation to the interface just above the layer exceeds trmin.
      
           if (trntdr(k,ij) > trmin ) then
 
              ! calculation over layers with penetrating radiation
 
              ! if (tflag) call t_startf('cice_solEdd1c1_'//tstr)
              tautot  = tau(k,ij)
              wtot    = w0(k,ij)
              gtot    = g(k,ij)
              ftot    = gtot*gtot
 
              ts   = taus(wtot,ftot,tautot)
              ws   = omgs(wtot,ftot)
              gs   = asys(gtot,ftot)
              lm   = el(ws,gs)
              ue   = u(ws,gs,lm)
              ! if (tflag) call t_stopf('cice_solEdd1c1_'//tstr)

              ! if (tflag) call t_startf('cice_solEdd1c2_'//tstr)
              ! compute level of fresnel refraction
              if( srftyp(i,j) < 2 ) then
                ! if snow over sea ice or bare sea ice, fresnel level is
                ! at base of sea ice SSL (and top of the sea ice DL); the
                ! snow SSL counts for one, then the number of snow layers,
                ! then the sea ice SSL which also counts for one:
                kfrsnl = nslyr + 2 
              else
                ! if ponded sea ice, fresnel level is the top of the pond 
                kfrsnl = 0
              endif

              ! mu0 is cosine solar zenith angle above the fresnel level; make 
              ! sure mu0 is large enough for stable and meaningful radiation
              ! solution: .01 is like sun just touching horizon with its lower edge

              mu0  = max(coszen(i,j),p01)

              ! mu0n is cosine solar zenith angle used to compute the layer
              ! Delta-Eddington solution; it is initially computed to be the
              ! value below the fresnel level, i.e. the cosine solar zenith 
              ! angle below the fresnel level for the refracted solar beam:

              mu0n = sqrt(c1-((c1-mu0*mu0)/(refindx*refindx)))

              ! if level k is above fresnel level and the cell is non-pond, use the
              ! non-refracted beam instead

              if( srftyp(i,j) < 2 .and. k < kfrsnl ) mu0n = mu0
              ! if (tflag) call t_stopf('cice_solEdd1c2_'//tstr)
 
              ! if (tflag) call t_startf('cice_solEdd1c3_'//tstr)
              extins = max(exp_min, exp(-lm*ts))
              ne = n(ue,extins)
 
              ! first calculation of rdif, tdif using Delta-Eddington formulas

              rdif_a(k,ij) = (ue+c1)*(ue-c1)*(c1/extins - extins)/ne
              tdif_a(k,ij) = c4*ue/ne
 
              ! evaluate rdir,tdir for direct beam
              trnlay(k,ij) = max(exp_min, exp(-ts/mu0n))
              alp = alpha(ws,mu0n,gs,lm)
              gam = gamma(ws,mu0n,gs,lm)
              apg = alp + gam
              amg = alp - gam
              rdir(k,ij) = amg*(tdif_a(k,ij)*trnlay(k,ij) - c1) + &
                           apg*rdif_a(k,ij)
              tdir(k,ij) = apg*tdif_a(k,ij) + &
                          (amg*rdif_a(k,ij) - (apg-c1))*trnlay(k,ij)
              ! if (tflag) call t_stopf('cice_solEdd1c3_'//tstr)
 
              ! recalculate rdif,tdif using direct angular integration over rdir,tdir,
              ! since Delta-Eddington rdif formula is not well-behaved (it is usually
              ! biased low and can even be negative); use ngmax angles and gaussian
              ! integration for most accuracy:
              ! if (tflag) call t_startf('cice_solEdd1ng_'//tstr)
              swt = c0
              smr = c0
              smt = c0
              do ng=1,ngmax
                mu  = gauspt(ng)
                gwt = gauswt(ng)
                swt = swt + mu*gwt
                trn = max(exp_min, exp(-ts/mu))
                alp = alpha(ws,mu,gs,lm)
                gam = gamma(ws,mu,gs,lm)
                apg = alp + gam
                amg = alp - gam
                rdr = amg*(tdif_a(k,ij)*trn-c1) + &
                      apg*rdif_a(k,ij)
                tdr = apg*tdif_a(k,ij) + &
                     (amg*rdif_a(k,ij)-(apg-c1))*trn
                smr = smr + mu*rdr*gwt
                smt = smt + mu*tdr*gwt
              enddo      ! ng
              rdif_a(k,ij) = smr/swt
              tdif_a(k,ij) = smt/swt

              ! homogeneous layer
              rdif_b(k,ij) = rdif_a(k,ij)
              tdif_b(k,ij) = tdif_a(k,ij)
              ! if (tflag) call t_stopf('cice_solEdd1ng_'//tstr)
 
              ! add fresnel layer to top of desired layer if either 
              ! air or snow overlies ice; we ignore refraction in ice 
              ! if a melt pond overlies it:

              if( k == kfrsnl ) then
                 ! if (tflag) call t_startf('cice_solEdd1kf_'//tstr)
                 ! compute fresnel reflection and transmission amplitudes
                 ! for two polarizations: 1=perpendicular and 2=parallel to
                 ! the plane containing incident, reflected and refracted rays.
                 R1 = (mu0 - refindx*mu0n) / & 
                      (mu0 + refindx*mu0n)
                 R2 = (refindx*mu0 - mu0n) / &
                      (refindx*mu0 + mu0n)
                 T1 = c2*mu0 / &
                      (mu0 + refindx*mu0n)
                 T2 = c2*mu0 / &
                      (refindx*mu0 + mu0n)
 
                 ! unpolarized light for direct beam
                 Rf_dir_a = p5 * (R1*R1 + R2*R2)
                 Tf_dir_a = p5 * (T1*T1 + T2*T2)*refindx*mu0n/mu0
 
                 ! precalculated diffuse reflectivities and transmissivities
                 ! for incident radiation above and below fresnel layer, using
                 ! the direct albedos and accounting for complete internal
                 ! reflection from below; precalculated because high order
                 ! number of gaussian points (~256) is required for convergence:
 
                 ! above
                 Rf_dif_a = cp063
                 Tf_dif_a = c1 - Rf_dif_a
                 ! below
                 Rf_dif_b = cp455
                 Tf_dif_b = c1 - Rf_dif_b

                 ! the k = kfrsnl layer properties are updated to combined 
                 ! the fresnel (refractive) layer, always taken to be above
                 ! the present layer k (i.e. be the top interface):

                 rintfc   = c1 / (c1-Rf_dif_b*rdif_a(kfrsnl,ij))
                 tdir(kfrsnl,ij)   = Tf_dir_a*tdir(kfrsnl,ij) + &
                                     Tf_dir_a*rdir(kfrsnl,ij) * &
                                     Rf_dif_b*rintfc*tdif_a(kfrsnl,ij)
                 rdir(kfrsnl,ij)   = Rf_dir_a + &
                                     Tf_dir_a*rdir(kfrsnl,ij) * &
                                     rintfc*Tf_dif_b
                 rdif_a(kfrsnl,ij) = Rf_dif_a + &
                                     Tf_dif_a*rdif_a(kfrsnl,ij) * &
                                     rintfc*Tf_dif_b
                 rdif_b(kfrsnl,ij) = rdif_b(kfrsnl,ij) + &
                                     tdif_b(kfrsnl,ij)*Rf_dif_b * &
                                     rintfc*tdif_a(kfrsnl,ij)
                 tdif_a(kfrsnl,ij) = Tf_dif_a*rintfc*tdif_a(kfrsnl,ij)
                 tdif_b(kfrsnl,ij) = tdif_b(kfrsnl,ij)*rintfc*Tf_dif_b

                 ! update trnlay to include fresnel transmission
                 trnlay(kfrsnl,ij) = Tf_dir_a*trnlay(kfrsnl,ij)
                 ! if (tflag) call t_stopf('cice_solEdd1kf_'//tstr)
              endif      ! k = kfrsnl

          endif ! trntdr(k,ij) > trmin

        enddo       ! k   end main level loop

        ! if (tflag) call t_stopf('cice_solEdd1kl_'//tstr)
        ! if (tflag) call t_startf('cice_solEdd1rf_'//tstr)
        ! compute total direct beam transmission, total transmission, and
        ! reflectivity for diffuse radiation (from below) for all layers
        ! above the underlying ocean; note that we ignore refraction between 
        ! sea ice and underlying ocean:
        !
        !       For k = klevp
        !
        !              layers       interface
        !
        !       ---------------------  k-1 
        !                k-1
        !       ---------------------  k
        !       \\\\\\\ ocean \\\\\\\

        k = klevp
        trndir(k,ij) = trndir(k-1,ij)*trnlay(k-1,ij)
        refkm1        = c1/(c1 - rdndif(k-1,ij)*rdif_a(k-1,ij))
        tdrrdir       = trndir(k-1,ij)*rdir(k-1,ij)
        tdndif        = trntdr(k-1,ij) - trndir(k-1,ij)
        trntdr(k,ij) = trndir(k-1,ij)*tdir(k-1,ij) + &
          (tdndif + tdrrdir*rdndif(k-1,ij))*refkm1*tdif_a(k-1,ij)
        rdndif(k,ij) = rdif_b(k-1,ij) + &
          (tdif_b(k-1,ij)*rdndif(k-1,ij)*refkm1*tdif_a(k-1,ij))
        trndif(k,ij) = trndif(k-1,ij)*refkm1*tdif_a(k-1,ij)
 
        ! compute reflectivity to direct and diffuse radiation for layers 
        ! below by adding succesive layers starting from the underlying 
        ! ocean and working upwards:
        !
        !              layers       interface
        !
        !       ---------------------  k
        !                 k
        !       ---------------------  k+1
        !                k+1
        !       ---------------------

        rupdir(klevp,ij) = albodr(ij)
        rupdif(klevp,ij) = albodf(ij)
        ! if (tflag) call t_stopf('cice_solEdd1rf_'//tstr)

      enddo       ! ij  
      if (tflag) call t_stopf('cice_solEdd1x_'//tstr)

      if (tflag) call t_startf('cice_solEdd2x_'//tstr)
      do k=klev,0,-1
        do ij = 1, icells_DE
          i = indxi_DE(ij)
          j = indxj_DE(ij)
          ! interface scattering
          refkp1        = c1/( c1 - rdif_b(k,ij)*rupdif(k+1,ij))
          ! dir from top layer plus exp tran ref from lower layer, interface
          ! scattered and tran thru top layer from below, plus diff tran ref
          ! from lower layer with interface scattering tran thru top from below
          rupdir(k,ij) = rdir(k,ij) + &
                        ( trnlay(k,ij)*rupdir(k+1,ij) + &
                         (tdir(k,ij)-trnlay(k,ij))*rupdif(k+1,ij) ) * &
                          refkp1*tdif_b(k,ij)
          ! dif from top layer from above, plus dif tran upwards reflected and
          ! interface scattered which tran top from below
          rupdif(k,ij) = rdif_a(k,ij) + &
                          tdif_a(k,ij)*rupdif(k+1,ij)* &
                          refkp1*tdif_b(k,ij)
        enddo       ! ij  
      enddo       ! k
      if (tflag) call t_stopf('cice_solEdd2x_'//tstr)

      end subroutine solution_dEdd

!=======================================================================
!BOP
!
! !IROUTINE: shortwave_dEdd_set_snow - routine for Delta-Eddington shortwave
!            that sets snow density and snow grain radius.
!
! !INTERFACE:
!
      subroutine shortwave_dEdd_set_snow(nx_block, ny_block, &
                                         icells,             &
                                         indxi,    indxj,    &
                                         aice,     vsno,     &
                                         Tsfc,     fs,       &
                                         rhosnw,   rsnw)
!
! !DESCRIPTION:
!
!   Set snow horizontal coverage, density and grain radius diagnostically 
!   for the Delta-Eddington solar radiation method.
!
! !REVISION HISTORY:
!
! author:  Bruce P. Briegleb, NCAR 
!          8 February 2007
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), &
         intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells                ! number of ice-covered grid cells

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi   , & ! compressed indices for ice-covered cells
         indxj

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         aice   , & ! concentration of ice
         vsno   , & ! volume of snow
         Tsfc       ! surface temperature 

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(out) :: &
         fs         ! horizontal coverage of snow

      real (kind=dbl_kind), dimension (nx_block,ny_block,nslyr), &
         intent(out) :: &
         rhosnw , & ! density in snow layer (kg/m3)
         rsnw       ! grain radius in snow layer (micro-meters)
!      
!EOP
! 
! !LOCAL PARAMETERS:
!

      integer (kind=int_kind) :: &
         i        , & ! longitude index
         j        , & ! latitude index
         ij       , & ! horizontal index, combines i and j loops
         ks           ! snow vertical index

      real (kind=dbl_kind) :: &
         hs  , & ! snow depth (m)
         fT  , & ! piecewise linear function of surface temperature
         dTs , & ! difference of Tsfc and Timelt
         rsnw_nm ! actual used nonmelt snow grain radius (micro-meters)

      real (kind=dbl_kind), parameter :: &
!        Move these to ice_constants
!        hsmin  = .0001_dbl_kind, & ! minimum allowed snow depth (m) for DE
!        hs0    = .0300_dbl_kind, & ! snow depth for transition to bare sea ice
         ! units for the following are 1.e-6 m (micro-meters)
         rsnw_fresh    =  100._dbl_kind, & ! freshly-fallen snow grain radius 
         rsnw_nonmelt  =  500._dbl_kind, & ! nonmelt snow grain radius
         rsnw_sig      =  250._dbl_kind    ! assumed sigma for snow grain radius

       real (kind=dbl_kind) :: &
         dT_mlt            , & ! change in temp to give non-melt to melt change
                               ! in snow grain radius
         rsnw_melt             ! melting snow grain radius

!-----------------------------------------------------------------------

      dT_mlt    = dT_mlt_in
      rsnw_melt = rsnw_melt_in

      fs(:,:)       = c0
      do ks = 1, nslyr
      do j = 1, ny_block
      do i = 1, nx_block
         rhosnw(i,j,ks) = c0
         rsnw(i,j,ks)   = c0
      enddo
      enddo
      enddo

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         ! set snow horizontal fraction
         if( aice(i,j) > puny ) then
           hs = vsno(i,j) / aice(i,j)
           if( hs < hsmin ) then
             fs(i,j) = c0
           else if( hs <= hs0 ) then
             fs(i,j) = hs/hs0
           else
             fs(i,j) = c1
           endif
           ! bare ice, temperature dependence
           dTs = Timelt - Tsfc(i,j)
           fT  = -min(dTs/dT_mlt-c1,c0)
           ! tune nonmelt snow grain radius if desired: note that
           ! the sign is negative so that if R_snw is 1, then the
           ! snow grain radius is reduced and thus albedo increased.
           rsnw_nm = rsnw_nonmelt - R_snw*rsnw_sig
           rsnw_nm = max(rsnw_nm, rsnw_fresh)
           rsnw_nm = min(rsnw_nm, rsnw_melt) 
           do ks = 1, nslyr
             ! snow density ccsm3 constant value
             rhosnw(i,j,ks) = rhos
             ! snow grain radius between rsnw_nonmelt and rsnw_melt
             rsnw(i,j,ks) = rsnw_nm + (rsnw_melt-rsnw_nm)*fT
             rsnw(i,j,ks) = max(rsnw(i,j,ks), rsnw_fresh)
             rsnw(i,j,ks) = min(rsnw(i,j,ks), rsnw_melt) 
           enddo        ! ks
         endif         ! aice(i,j) > puny
      enddo          ! ij

      end subroutine shortwave_dEdd_set_snow

!=======================================================================
!BOP
!
! !IROUTINE: shortwave_dEdd_set_pond - routine for Delta-Eddington shortwave
!            that sets pond fraction and depth.
!
! !INTERFACE:
!
      subroutine shortwave_dEdd_set_pond(nx_block, ny_block, &
                                         icells,             &
                                         indxi,    indxj,    &
                                         aice,     Tsfc,     &
                                         fs,       fp,       &
                                         hp)
!
! !DESCRIPTION:
!
!   Set pond fraction and depth diagnostically for
!   the Delta-Eddington solar radiation method.
!
! !REVISION HISTORY:
!
! author:  Bruce P. Briegleb, NCAR 
!          8 February 2007
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), &
         intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells                ! number of ice-covered grid cells

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi   , & ! compressed indices for ice-covered cells
         indxj

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         aice   , & ! concentration of ice
         Tsfc   , & ! surface temperature
         fs         ! horizontal coverage of snow

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(out) :: &
         fp     , & ! pond fractional coverage (0 to 1)
         hp         ! pond depth (m)

!      
!EOP
! 
! !LOCAL PARAMETERS:
!

      integer (kind=int_kind) :: &
         i   , & ! longitude index
         j   , & ! latitude index
         ij      ! horizontal index, combines i and j loops

      real (kind=dbl_kind) :: &
         fT  , & ! piecewise linear function of surface temperature
         dTs     ! difference of Tsfc and Timelt

      real (kind=dbl_kind), parameter :: &
         dT_mlt    = c1   ! change in temp for pond fraction and depth

!-----------------------------------------------------------------------

      do j = 1, ny_block
      do i = 1, nx_block
         fp(i,j) = c0
         hp(i,j) = c0
      enddo
      enddo

      ! find pond fraction and depth for ice points
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         if (aice(i,j) > puny) then
           ! bare ice, temperature dependence
           dTs = Timelt - Tsfc(i,j)
           fT  = -min(dTs/dT_mlt-c1,c0)
           ! pond
           fp(i,j) = 0.3_dbl_kind*fT*(c1-fs(i,j))
           hp(i,j) = 0.3_dbl_kind*fT*(c1-fs(i,j))
         endif
      enddo          ! ij

      end subroutine shortwave_dEdd_set_pond

! End Delta-Eddington shortwave method
!=======================================================================

      end module ice_shortwave

!=======================================================================
