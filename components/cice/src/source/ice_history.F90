!=======================================================================
!
!BOP
!
! !MODULE: ice_history - ice model history files
!
! Output files: netCDF or binary data, Fortran unformatted dumps
!
! The following variables are currently hard-wired as snapshots 
!   (instantaneous rather than time-averages):
!   divu, shear, sig1, sig2, trsig, mlt_onset, frz_onset, hisnap, aisnap
!
! Options for histfreq: '1','h','d','m','y','x', where x means that
!   output stream will not be used (recommended for efficiency).  
! histfreq_n can be any nonnegative integer, where 0 means that the 
!   corresponding histfreq frequency will not be used.
! The flags (f_<field>) can be set to '1','h','d','m','y' or 'x', where
!   n means the field will not be written.  To output the same field at
!   more than one frequency, for instance monthy and daily, set 
!   f_<field> = 'md'.
!
!
! !REVISION HISTORY:
!  SVN:$Id: ice_history.F90 61 2007-04-25 17:50:16Z dbailey $
!
! authors Tony Craig and Bruce Briegleb, NCAR
!         Elizabeth C. Hunke and William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! 2004 WHL: Block structure added 
! 2006 ECH: Accepted some CCSM code into mainstream CICE
!           Added ice_present, aicen, vicen; removed aice1...10, vice1...1.
!           Added histfreq_n and histfreq='h' options, removed histfreq='w'
!           Converted to free source form (F90)
!           Added option for binary output instead of netCDF
! 2009 D Bailey and ECH: Generalized for multiple frequency output
!
! !INTERFACE:
!
      module ice_history
!
! !USES:
!
      use ice_kinds_mod
      use ice_broadcast
      use ice_communicate, only: my_task, master_task
      use ice_blocks
      use ice_grid
      use ice_read_write
      use ice_fileunits
      use ice_history_fields
      use ice_history_write
!
!EOP
!
      implicit none
      save
      
      character (len=16) :: vname_in     ! variable name
      character (len=55) :: vdesc_in     ! variable description
      character (len=55) :: vcomment_in  ! variable description

      !---------------------------------------------------------------
      ! primary info for the history file
      !---------------------------------------------------------------

      character (len=16), parameter :: &
         tstr  = 'TLON TLAT time', & ! vcoord for T cell quantities
         ustr  = 'ULON ULAT time', & ! vcoord for U cell quantities
         tcstr = 'area: tarea'   , & ! vcellmeas for T cell quantities
         ucstr = 'area: uarea'       ! vcellmeas for U cell quantities

      real (kind=dbl_kind) :: &
         avgct(max_nstrm)            ! average sample counter

      !---------------------------------------------------------------
      ! logical flags: write to output file if true
      !---------------------------------------------------------------

       logical (kind=log_kind) :: &
           f_tmask     = .true., f_blkmask    = .true., &
           f_tarea     = .true., f_uarea      = .true., &
           f_dxt       = .true., f_dyt        = .true., &
           f_dxu       = .true., f_dyu        = .true., &
           f_HTN       = .true., f_HTE        = .true., &
           f_ANGLE     = .true., f_ANGLET     = .true.

      character (len=max_nstrm) :: &
!          f_example   = 'mdxxx', &
           f_hi        = 'mxxxx', f_hs         = 'mxxxx', &
           f_fs        = 'mxxxx', &
           f_Tsfc      = 'mxxxx', f_aice       = 'mxxxx', &
           f_uvel      = 'mxxxx', f_vvel       = 'mxxxx', &
           f_transix   = 'mxxxx', f_transiy    = 'mxxxx', &
           f_qi        = 'mxxxx', f_qs         = 'mxxxx', &
           f_fswdn     = 'mxxxx', f_fswup      = 'mxxxx', &
           f_flwdn     = 'mxxxx', &
           f_snow      = 'mxxxx', f_snow_ai    = 'mxxxx', &
           f_rain      = 'mxxxx', f_rain_ai    = 'mxxxx', &
           f_faero_atm = 'mxxxx', f_faero_ocn  = 'mxxxx', &
           f_sst       = 'mxxxx', f_sss        = 'mxxxx', &
           f_uocn      = 'mxxxx', f_vocn       = 'mxxxx', &
           f_frzmlt    = 'mxxxx', &
           f_fswfac    = 'mxxxx', &
#if (defined AEROFRC) || (defined PONDFRC) || (defined CCSM3FRC)
           f_fswsfc_ai  = 'mxxxx', &
           f_fswint_ai  = 'mxxxx', &
#endif
#ifdef AEROFRC
           f_dfswabs_noaero = 'mxxxx', &
           f_dfswsfc_noaero = 'mxxxx', &
           f_dfswint_noaero = 'mxxxx', &
           f_dfswthru_noaero= 'mxxxx', &
           f_dalvdr_noaero  = 'mxxxx', f_dalidr_noaero     = 'mxxxx', &
           f_dalvdf_noaero  = 'mxxxx', f_dalidf_noaero     = 'mxxxx', &
           f_dalbice_noaero = 'mxxxx', f_dalbsno_noaero    = 'mxxxx', &
           f_dalbpnd_noaero = 'mxxxx',                               &
#endif
#ifdef CCSM3FRC
           f_dfswabs_ccsm3 = 'mxxxx', &
           f_dfswsfc_ccsm3 = 'mxxxx', &
           f_dfswint_ccsm3 = 'mxxxx', &
           f_dfswthru_ccsm3= 'mxxxx', &
           f_dalvdr_ccsm3  = 'mxxxx', f_dalidr_ccsm3     = 'mxxxx', &
           f_dalvdf_ccsm3  = 'mxxxx', f_dalidf_ccsm3     = 'mxxxx', &
           f_dalbice_ccsm3 = 'mxxxx', f_dalbsno_ccsm3    = 'mxxxx', &
#endif
#ifdef PONDFRC
           f_dfswabs_nopond = 'mxxxx', &
           f_dfswsfc_nopond = 'mxxxx', &
           f_dfswint_nopond = 'mxxxx', &
           f_dfswthru_nopond= 'mxxxx', &
           f_dalvdr_nopond  = 'mxxxx', f_dalidr_nopond     = 'mxxxx', &
           f_dalvdf_nopond  = 'mxxxx', f_dalidf_nopond     = 'mxxxx', &
           f_dalbice_nopond = 'mxxxx', f_dalbsno_nopond    = 'mxxxx', &
           f_dalbpnd_nopond = 'mxxxx',                               &
#endif
           f_fswabs    = 'mxxxx', f_fswabs_ai  = 'mxxxx', &
           f_alvdr     = 'mxxxx', f_alidr      = 'mxxxx', &
           f_alvdf     = 'mxxxx', f_alidf      = 'mxxxx', &
           f_albice    = 'mxxxx', f_albsno     = 'mxxxx', &
           f_albpnd    = 'mxxxx', f_coszen     = 'mxxxx', &
           f_flat      = 'mxxxx', f_flat_ai    = 'mxxxx', &
           f_fsens     = 'mxxxx', f_fsens_ai   = 'mxxxx', &
           f_flwup     = 'mxxxx', f_flwup_ai   = 'mxxxx', &
           f_evap      = 'mxxxx', f_evap_ai    = 'mxxxx', &
           f_Tair      = 'mxxxx', &
           f_Tref      = 'mxxxx', f_Qref       = 'mxxxx', &
           f_congel    = 'mxxxx', f_frazil     = 'mxxxx', &
           f_snoice    = 'mxxxx', &
           f_meltt     = 'mxxxx', f_melts      = 'mxxxx', &
           f_meltb     = 'mxxxx', f_meltl      = 'mxxxx', &
           f_fresh     = 'mxxxx', f_fresh_ai   = 'mxxxx', &
           f_fsalt     = 'mxxxx', f_fsalt_ai   = 'mxxxx', &
           f_fhocn     = 'mxxxx', f_fhocn_ai   = 'mxxxx', &
           f_fswthru   = 'mxxxx', f_fswthru_ai = 'mxxxx', &
           f_strairx   = 'mxxxx', f_strairy    = 'mxxxx', &
           f_strtltx   = 'mxxxx', f_strtlty    = 'mxxxx', &
           f_strcorx   = 'mxxxx', f_strcory    = 'mxxxx', &
           f_strocnx   = 'mxxxx', f_strocny    = 'mxxxx', &
           f_strintx   = 'mxxxx', f_strinty    = 'mxxxx', &
           f_strength  = 'mxxxx', f_opening    = 'mxxxx', &
           f_divu      = 'mxxxx', f_shear      = 'mxxxx', &
           f_sig1      = 'mxxxx', f_sig2       = 'mxxxx', &
           f_dvidtt    = 'mxxxx', f_dvidtd     = 'mxxxx', &
           f_daidtt    = 'mxxxx', f_daidtd     = 'mxxxx', &
           f_mlt_onset = 'mxxxx', f_frz_onset  = 'mxxxx', &
           f_dardg1dt  = 'mxxxx', f_dardg2dt   = 'mxxxx', &
           f_dvirdgdt  = 'mxxxx', f_iage       = 'mxxxx', &
           f_ardg      = 'mxxxx', f_vrdg       = 'mxxxx', &
           f_alvl      = 'mxxxx', f_vlvl       = 'mxxxx', &
           f_FY        = 'mxxxx',                         &
           f_aeron     = 'xxxxx', f_aero       = 'xxxxx', &
           f_apond     = 'xxxxx', f_apondn     = 'xxxxx', &
           f_hisnap    = 'mxxxx', f_aisnap     = 'mxxxx', &
           f_aicen     = 'mxxxx', f_vicen      = 'mxxxx', &
           f_trsig     = 'mxxxx', f_icepresent = 'mxxxx', &
           f_fsurf_ai  = 'mxxxx', f_fcondtop_ai= 'mxxxx', &
           f_fmeltt_ai = 'xxxxx',                     &
           f_fsurfn_ai = 'xxxxx',f_fcondtopn_ai= 'xxxxx', &
           f_fmelttn_ai= 'xxxxx', f_flatn_ai   = 'xxxxx'

      !---------------------------------------------------------------
      ! namelist variables 
      !---------------------------------------------------------------

      namelist / icefields_nml /     &
           f_tmask    , f_blkmask  , &
           f_tarea    , f_uarea    , &
           f_dxt      , f_dyt      , &
           f_dxu      , f_dyu      , &
           f_HTN      , f_HTE      , &
           f_ANGLE    , f_ANGLET   , &
           f_bounds   , &
!
!          f_example               , &
           f_hi,        f_hs       , &
           f_fs,                     &
           f_Tsfc,      f_aice     , &
           f_uvel,      f_vvel     , &
           f_transix,   f_transiy  , &
           f_qi,        f_qs       , &
           f_fswdn,     f_fswup    , &
           f_flwdn,                  &
           f_snow,      f_snow_ai  , &     
           f_rain,      f_rain_ai  , &
           f_faero_atm, f_faero_ocn, &
           f_sst,       f_sss      , &
           f_uocn,      f_vocn     , &
           f_frzmlt                , &
           f_fswfac                , &
           f_fswabs,    f_fswabs_ai, &
           f_alvdr,     f_alidr    , &
           f_alvdf,     f_alidf    , &
           f_albice,    f_albsno   , &
           f_albpnd,    f_coszen   , &
           f_flat,      f_flat_ai  , &
           f_fsens,     f_fsens_ai , &
           f_flwup,     f_flwup_ai , &
           f_evap,      f_evap_ai  , &
           f_Tair                  , &
           f_Tref,      f_Qref     , &
           f_congel,    f_frazil   , &
           f_snoice,                 & 
           f_meltt,     f_melts    , &
           f_meltb,     f_meltl    , &
           f_fresh,     f_fresh_ai , &  
           f_fsalt,     f_fsalt_ai , &
           f_fhocn,     f_fhocn_ai , &
           f_fswthru,   f_fswthru_ai,&
#if (defined AEROFRC) || (defined CCSM3FRC) || (defined PONDFRC)
           f_fswsfc_ai, f_fswint_ai ,&
#endif
#ifdef AEROFRC
           f_dfswabs_noaero, &
           f_dfswsfc_noaero, &
           f_dfswint_noaero, &
           f_dfswthru_noaero, &
           f_dalvdr_noaero , f_dalidr_noaero, &
           f_dalvdf_noaero , f_dalidf_noaero, &
           f_dalbice_noaero, f_dalbsno_noaero, &
           f_dalbpnd_noaero,                   &
#endif
#ifdef CCSM3FRC
           f_dfswabs_ccsm3, &
           f_dfswsfc_ccsm3, &
           f_dfswint_ccsm3, &
           f_dfswthru_ccsm3, &
           f_dalvdr_ccsm3 , f_dalidr_ccsm3, &
           f_dalvdf_ccsm3 , f_dalidf_ccsm3, &
           f_dalbice_ccsm3, f_dalbsno_ccsm3, &
#endif
#ifdef PONDFRC
           f_dfswabs_nopond, &
           f_dfswsfc_nopond, &
           f_dfswint_nopond, &
           f_dfswthru_nopond, &
           f_dalvdr_nopond , f_dalidr_nopond, &
           f_dalvdf_nopond , f_dalidf_nopond, &
           f_dalbice_nopond, f_dalbsno_nopond, &
           f_dalbpnd_nopond,                   &
#endif
           f_strairx,   f_strairy  , &
           f_strtltx,   f_strtlty  , &
           f_strcorx,   f_strcory  , &
           f_strocnx,   f_strocny  , &
           f_strintx,   f_strinty  , &
           f_strength,  f_opening  , &
           f_divu,      f_shear    , &
           f_sig1,      f_sig2     , &
           f_dvidtt,    f_dvidtd   , &
           f_daidtt,    f_daidtd   , &
           f_mlt_onset, f_frz_onset, &
           f_dardg1dt,  f_dardg2dt , &
           f_dvirdgdt,  f_iage     , &
           f_ardg,      f_vrdg     , &
           f_alvl,      f_vlvl     , &
           f_aeron,     f_aero     , &
           f_FY                    , &
           f_apond,     f_apondn   , &
           f_hisnap,    f_aisnap   , &
           f_aicen,     f_vicen    , &
           f_trsig,     f_icepresent,&
           f_fsurf_ai,  f_fcondtop_ai,&
           f_fmeltt_ai,              &
           f_fsurfn_ai,f_fcondtopn_ai,&
           f_fmelttn_ai,f_flatn_ai

!=======================================================================

      contains

!=======================================================================
!
!BOP
!
! !IROUTINE: init_hist - initialize history files
!
! !INTERFACE:
!
      subroutine init_hist (dt)
!
! !DESCRIPTION:
!
! Initialize history files
!
! !REVISION HISTORY:
!
! authors Tony Craig, NCAR
!         Elizabeth C. Hunke, LANL
!         C.M. Bitz, UW
!         Bruce P. Briegleb, NCAR
!         William H. Lipscomb, LANL
!
! !USES:
!
      use ice_constants
      use ice_calendar, only: yday, days_per_year, histfreq, &
                              histfreq_n, nstreams
      use ice_flux, only: mlt_onset, frz_onset, albcnt
      use ice_restart, only: restart
      use ice_state, only: tr_aero, tr_iage, tr_FY, tr_pond, tr_lvl
      use ice_exit
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step
!
!EOP
!
      integer (kind=int_kind) :: n, k, ns, ns1, ns2, lenf
      integer (kind=int_kind) :: hfreqn
      integer (kind=int_kind), dimension(max_nstrm) :: &
         ntmp
      integer (kind=int_kind) :: nml_error ! namelist i/o error flag

      character (len=3) :: nchar
      character (len=40) :: stmp

      !-----------------------------------------------------------------
      ! read namelist
      !-----------------------------------------------------------------

      call get_fileunit(nu_nml)
      if (my_task == master_task) then
         open (nu_nml, file=nml_filename, status='old',iostat=nml_error)
         if (nml_error /= 0) then
            nml_error = -1
         else
            nml_error =  1
         endif
         do while (nml_error > 0)
            read(nu_nml, nml=icefields_nml,iostat=nml_error)
         end do
         if (nml_error == 0) close(nu_nml)
      endif
      call release_fileunit(nu_nml)

      call broadcast_scalar(nml_error, master_task)
      if (nml_error /= 0) then
         close (nu_nml)
         call abort_ice('ice: error reading icefields_nml')
      endif

      ! histfreq options ('1','h','d','m','y')
      nstreams = 0
      do ns = 1, max_nstrm
         if (histfreq(ns) == '1' .or. histfreq(ns) == 'h' .or. &
             histfreq(ns) == 'd' .or. histfreq(ns) == 'm' .or. &
             histfreq(ns) == 'y') then
                nstreams = nstreams + 1
         else if (histfreq(ns) /= 'x') then
             call abort_ice('ice: histfreq contains illegal element')
         endif
      enddo
      if (nstreams == 0) write (nu_diag,*) 'WARNING: No history output'
      do ns1 = 1, nstreams
         do ns2 = 1, nstreams
            if (histfreq(ns1) == histfreq(ns2) .and. ns1/=ns2 &
               .and. my_task == master_task) then
               call abort_ice('ice: histfreq elements must be unique')
            endif
         enddo
      enddo

      if (.not. tr_iage) f_iage   = 'xxxxx'
      if (.not. tr_FY)   f_FY     = 'xxxxx'
      if (.not. tr_pond) then
          f_apond  = 'xxxxx'
          f_apondn = 'xxxxx'
      endif
      if (.not. tr_lvl) then
          f_ardg = 'xxxxx'
          f_vrdg = 'xxxxx'
          f_alvl = 'xxxxx'
          f_vlvl = 'xxxxx'
      endif
      if (.not. tr_aero) then
         f_faero_atm = 'xxxxx'
         f_faero_ocn = 'xxxxx'
         f_aero      = 'xxxxx' 
         f_aeron     = 'xxxxx'
      endif

      ! these must be output at the same frequency because of 
      ! cos(zenith angle) averaging
      if (f_albsno(1:1) /= 'x') f_albsno = f_albice
      if (f_albpnd(1:1) /= 'x') f_albpnd = f_albice
      if (f_coszen(1:1) /= 'x') f_coszen = f_albice

      ! to prevent array-out-of-bounds when aggregating
      if (f_fmeltt_ai(1:1) /= 'x') f_fmelttn_ai = f_fmeltt_ai

#ifndef ncdf
      f_bounds = .false.
#endif

      call broadcast_scalar (f_tmask, master_task)
      call broadcast_scalar (f_blkmask, master_task)
      call broadcast_scalar (f_tarea, master_task)
      call broadcast_scalar (f_uarea, master_task)
      call broadcast_scalar (f_dxt, master_task)
      call broadcast_scalar (f_dyt, master_task)
      call broadcast_scalar (f_dxu, master_task)
      call broadcast_scalar (f_dyu, master_task)
      call broadcast_scalar (f_HTN, master_task)
      call broadcast_scalar (f_HTE, master_task)
      call broadcast_scalar (f_ANGLE, master_task)
      call broadcast_scalar (f_ANGLET, master_task)
      call broadcast_scalar (f_bounds, master_task)

!     call broadcast_scalar (f_example, master_task)
      call broadcast_scalar (f_hi, master_task)
      call broadcast_scalar (f_hs, master_task)
      call broadcast_scalar (f_fs, master_task)
      call broadcast_scalar (f_Tsfc, master_task)
      call broadcast_scalar (f_aice, master_task)
      call broadcast_scalar (f_uvel, master_task)
      call broadcast_scalar (f_vvel, master_task)
      call broadcast_scalar (f_transix, master_task)
      call broadcast_scalar (f_transiy, master_task)
      call broadcast_scalar (f_fswdn, master_task)
      call broadcast_scalar (f_fswup, master_task)
      call broadcast_scalar (f_flwdn, master_task)
      call broadcast_scalar (f_snow, master_task)
      call broadcast_scalar (f_snow_ai, master_task)
      call broadcast_scalar (f_rain, master_task)
      call broadcast_scalar (f_rain_ai, master_task)
      call broadcast_scalar (f_faero_atm, master_task)
      call broadcast_scalar (f_faero_ocn, master_task)
      call broadcast_scalar (f_sst, master_task)
      call broadcast_scalar (f_sss, master_task)
      call broadcast_scalar (f_uocn, master_task)
      call broadcast_scalar (f_vocn, master_task)
      call broadcast_scalar (f_frzmlt, master_task)
      call broadcast_scalar (f_fswfac, master_task)
      call broadcast_scalar (f_fswabs, master_task)
      call broadcast_scalar (f_fswabs_ai, master_task)
#if (defined AEROFRC) || (defined PONDFRC) || (defined CCSM3FRC)
      call broadcast_scalar (f_fswsfc_ai, master_task)
      call broadcast_scalar (f_fswint_ai, master_task)
#endif
#ifdef AEROFRC
      call broadcast_scalar (f_dfswabs_noaero, master_task)
      call broadcast_scalar (f_dfswsfc_noaero, master_task)
      call broadcast_scalar (f_dfswint_noaero, master_task)
      call broadcast_scalar (f_dfswthru_noaero, master_task)
      call broadcast_scalar (f_dalvdr_noaero, master_task)
      call broadcast_scalar (f_dalidr_noaero, master_task)
      call broadcast_scalar (f_dalvdf_noaero, master_task)
      call broadcast_scalar (f_dalidf_noaero, master_task)
      call broadcast_scalar (f_dalbice_noaero, master_task)
      call broadcast_scalar (f_dalbsno_noaero, master_task)
      call broadcast_scalar (f_dalbpnd_noaero, master_task)
#endif
#ifdef CCSM3FRC
      call broadcast_scalar (f_dfswabs_ccsm3, master_task)
      call broadcast_scalar (f_dfswsfc_ccsm3, master_task)
      call broadcast_scalar (f_dfswint_ccsm3, master_task)
      call broadcast_scalar (f_dfswthru_ccsm3, master_task)
      call broadcast_scalar (f_dalvdr_ccsm3, master_task)
      call broadcast_scalar (f_dalidr_ccsm3, master_task)
      call broadcast_scalar (f_dalvdf_ccsm3, master_task)
      call broadcast_scalar (f_dalidf_ccsm3, master_task)
      call broadcast_scalar (f_dalbice_ccsm3, master_task)
      call broadcast_scalar (f_dalbsno_ccsm3, master_task)
#endif
#ifdef PONDFRC
      call broadcast_scalar (f_dfswabs_nopond, master_task)
      call broadcast_scalar (f_dfswsfc_nopond, master_task)
      call broadcast_scalar (f_dfswint_nopond, master_task)
      call broadcast_scalar (f_dfswthru_nopond, master_task)
      call broadcast_scalar (f_dalvdr_nopond, master_task)
      call broadcast_scalar (f_dalidr_nopond, master_task)
      call broadcast_scalar (f_dalvdf_nopond, master_task)
      call broadcast_scalar (f_dalidf_nopond, master_task)
      call broadcast_scalar (f_dalbice_nopond, master_task)
      call broadcast_scalar (f_dalbsno_nopond, master_task)
      call broadcast_scalar (f_dalbpnd_nopond, master_task)
#endif
      call broadcast_scalar (f_alvdr, master_task)
      call broadcast_scalar (f_alidr, master_task)
      call broadcast_scalar (f_alvdf, master_task)
      call broadcast_scalar (f_alidf, master_task)
      call broadcast_scalar (f_albice, master_task)
      call broadcast_scalar (f_albsno, master_task)
      call broadcast_scalar (f_albpnd, master_task)
      call broadcast_scalar (f_coszen, master_task)
      call broadcast_scalar (f_flat, master_task)
      call broadcast_scalar (f_flat_ai, master_task)
      call broadcast_scalar (f_fsens, master_task)
      call broadcast_scalar (f_fsens_ai, master_task)
      call broadcast_scalar (f_flwup, master_task)
      call broadcast_scalar (f_flwup_ai, master_task)
      call broadcast_scalar (f_evap, master_task)
      call broadcast_scalar (f_evap_ai, master_task)
      call broadcast_scalar (f_qi, master_task)
      call broadcast_scalar (f_qs, master_task)
      call broadcast_scalar (f_Tair, master_task)
      call broadcast_scalar (f_Tref, master_task)
      call broadcast_scalar (f_Qref, master_task)
      call broadcast_scalar (f_congel, master_task)
      call broadcast_scalar (f_frazil, master_task)
      call broadcast_scalar (f_snoice, master_task)
      call broadcast_scalar (f_meltt, master_task)
      call broadcast_scalar (f_meltb, master_task)
      call broadcast_scalar (f_meltl, master_task)
      call broadcast_scalar (f_melts, master_task)
      call broadcast_scalar (f_fresh, master_task)
      call broadcast_scalar (f_fresh_ai, master_task)
      call broadcast_scalar (f_fsalt, master_task)
      call broadcast_scalar (f_fsalt_ai, master_task)
      call broadcast_scalar (f_fhocn, master_task)
      call broadcast_scalar (f_fhocn_ai, master_task)
      call broadcast_scalar (f_fswthru, master_task)
      call broadcast_scalar (f_fswthru_ai, master_task)
      call broadcast_scalar (f_strairx, master_task)
      call broadcast_scalar (f_strairy, master_task)
      call broadcast_scalar (f_strtltx, master_task)
      call broadcast_scalar (f_strtlty, master_task)
      call broadcast_scalar (f_strcorx, master_task)
      call broadcast_scalar (f_strcory, master_task)
      call broadcast_scalar (f_strocnx, master_task)
      call broadcast_scalar (f_strocny, master_task)
      call broadcast_scalar (f_strintx, master_task)
      call broadcast_scalar (f_strinty, master_task)
      call broadcast_scalar (f_strength, master_task)
      call broadcast_scalar (f_opening, master_task)
      call broadcast_scalar (f_divu, master_task)
      call broadcast_scalar (f_shear, master_task)
      call broadcast_scalar (f_sig1, master_task)
      call broadcast_scalar (f_sig2, master_task)
      call broadcast_scalar (f_dvidtt, master_task)
      call broadcast_scalar (f_dvidtd, master_task)
      call broadcast_scalar (f_daidtt, master_task)
      call broadcast_scalar (f_daidtd, master_task)
      call broadcast_scalar (f_mlt_onset, master_task)
      call broadcast_scalar (f_frz_onset, master_task)
      call broadcast_scalar (f_dardg1dt, master_task)
      call broadcast_scalar (f_dardg2dt, master_task)
      call broadcast_scalar (f_dvirdgdt, master_task)
      call broadcast_scalar (f_aisnap, master_task)
      call broadcast_scalar (f_hisnap, master_task)
      call broadcast_scalar (f_aicen, master_task)
      call broadcast_scalar (f_vicen, master_task)
      call broadcast_scalar (f_trsig, master_task)
      call broadcast_scalar (f_icepresent, master_task)
      call broadcast_scalar (f_fsurf_ai, master_task)
      call broadcast_scalar (f_fcondtop_ai, master_task)
      call broadcast_scalar (f_fmeltt_ai, master_task)
      call broadcast_scalar (f_fsurfn_ai, master_task)
      call broadcast_scalar (f_fcondtopn_ai, master_task)
      call broadcast_scalar (f_fmelttn_ai, master_task)
      call broadcast_scalar (f_flatn_ai, master_task)

      call broadcast_scalar (f_aero, master_task)
      call broadcast_scalar (f_aeron, master_task)
      call broadcast_scalar (f_iage, master_task)
      call broadcast_scalar (f_FY, master_task)
      call broadcast_scalar (f_alvl, master_task)
      call broadcast_scalar (f_vlvl, master_task)
      call broadcast_scalar (f_apond, master_task)
      call broadcast_scalar (f_apondn, master_task)

      do ns1 = 1, nstreams

!!!!! begin example
!      if (f_example(1:1) /= 'x') &
!         call define_hist_field(n_example,"example","m",tstr, tcstr, & 
!            "example: mean ice thickness",                           &
!            "ice volume per unit grid cell area", c1, c0,            &
!            ns1, f_example)
!!!!! end example

      if (f_hi(1:1) /= 'x') &
         call define_hist_field(n_hi,"hi","m",tstr, tcstr,         & 
             "grid cell mean ice thickness",                       &
             "ice volume per unit grid cell area", c1, c0,         &
             ns1, f_hi)

      if (f_hs(1:1) /= 'x') &
         call define_hist_field(n_hs,"hs","m",tstr, tcstr,         &
             "grid cell mean snow thickness",                      &
             "snow volume per unit grid cell area", c1, c0,        &
             ns1, f_hs)

      if (f_fs(1:1) /= 'x') &
         call define_hist_field(n_fs,"fs"," ",tstr, tcstr,         &
             "grid cell mean snow fraction",                       &
             "none", c1, c0,                                       &
             ns1, f_fs)

      if (f_Tsfc(1:1) /= 'x') &
         call define_hist_field(n_Tsfc,"Tsfc","degC",tstr, tcstr,  &
             "snow/ice surface temperature",                       &
             "averaged with Tf if no ice is present", c1, c0,      &
             ns1, f_Tsfc)

      if (f_aice(1:1) /= 'x') &
         call define_hist_field(n_aice,"aice","%",tstr, tcstr,     &
             "ice area  (aggregate)",                              &
             "none", c100, c0,                                     &
             ns1, f_aice)

      if (f_qi(1:1) /= 'x') &
         call define_hist_field(n_qi,"qi","J",tstr, tcstr, &
             "internal ice heat content",                  &
             "none", c1, c0,                               &
             ns1, f_qi)
      
      if (f_qs(1:1) /= 'x') &
         call define_hist_field(n_qs,"qs","J",tstr, tcstr, &
             "internal snow heat content",                 &
             "none", c1, c0,                               &
             ns1, f_qs)
      
      if (f_uvel(1:1) /= 'x') &
         call define_hist_field(n_uvel,"uvel","cm/s",ustr, ucstr,  &
             "ice velocity (x)",                                   &
             "positive is x direction on U grid", m_to_cm, c0,     &
             ns1, f_uvel)

      if (f_vvel(1:1) /= 'x') &
         call define_hist_field(n_vvel,"vvel","cm/s",ustr, ucstr,  &
             "ice velocity (y)",                                   &
             "positive is y direction on U grid", m_to_cm, c0,     &
             ns1, f_vvel)

      if (f_transix(1:1) /= 'x') &
         call define_hist_field(n_transix,"transix","kg/s",tstr, tcstr,  &
             "ice mass transport (x) on East side",                      &
             "positive is x direction on U grid", c1, c0,     &
             ns1, f_transix)

      if (f_transiy(1:1) /= 'x') &
         call define_hist_field(n_transiy,"transiy","kg/s",tstr, tcstr,  &
             "ice mass transport (y) on North side",                     &
             "positive is y direction on U grid", c1, c0,     &
             ns1, f_transiy)

      if (f_fswdn(1:1) /= 'x') &
         call define_hist_field(n_fswdn,"fswdn","W/m^2",tstr, tcstr, &
             "down solar flux",                                      &
             "positive downward", c1, c0,                            &
             ns1, f_fswdn)

      if (f_fswup(1:1) /= 'x') &
         call define_hist_field(n_fswup,"fswup","W/m^2",tstr, tcstr, &
             "upward solar flux",                                    &
             "positive downward", c1, c0,                            &
             ns1, f_fswup)

      if (f_flwdn(1:1) /= 'x') &
         call define_hist_field(n_flwdn,"flwdn","W/m^2",tstr, tcstr, &
             "down longwave flux",                                   &
             "positive downward", c1, c0,                            &
             ns1, f_flwdn)

      if (f_snow(1:1) /= 'x') &
         call define_hist_field(n_snow,"snow","cm/day",tstr, tcstr, &
             "snowfall rate (cpl)",                                 &
             "none", mps_to_cmpdy/rhofresh, c0,                     &
             ns1, f_snow)

      if (f_snow_ai(1:1) /= 'x') &
         call define_hist_field(n_snow_ai,"snow_ai","cm/day",tstr, tcstr, &
             "snowfall rate",                                             &
             "weighted by ice", mps_to_cmpdy/rhofresh, c0,                &
             ns1, f_snow_ai)

      if (f_rain(1:1) /= 'x') &
         call define_hist_field(n_rain,"rain","cm/day",tstr, tcstr, &
             "rainfall rate (cpl)",                                 &
             "none", mps_to_cmpdy/rhofresh, c0,                     &
             ns1, f_rain)

      if (f_rain_ai(1:1) /= 'x') &
         call define_hist_field(n_rain_ai,"rain_ai","cm/day",tstr, tcstr, &
             "rainfall rate",                                             &
             "weighted by ice", mps_to_cmpdy/rhofresh, c0,                &
             ns1, f_rain_ai)

      if (f_sst(1:1) /= 'x') &
         call define_hist_field(n_sst,"sst","C",tstr, tcstr, &
             "sea surface temperature",                      &
             "none", c1, c0,                                 &
             ns1, f_sst)

      if (f_sss(1:1) /= 'x') &
         call define_hist_field(n_sss,"sss","ppt",tstr, tcstr, &
             "sea surface salinity",                           &
             "none", c1, c0,                                   &
             ns1, f_sss)

      if (f_uocn(1:1) /= 'x') &
         call define_hist_field(n_uocn,"uocn","cm/s",ustr, ucstr, &
             "ocean current (x)",                                 &
             "positive is x direction on U grid", m_to_cm, c0,    &
             ns1, f_uocn)

      if (f_vocn(1:1) /= 'x') &
         call define_hist_field(n_vocn,"vocn","cm/s",ustr, ucstr, &
             "ocean current (y)",                                 &
             "positive is y direction on U grid", m_to_cm, c0,    &
             ns1, f_vocn)

      if (f_fswfac(1:1) /= 'x') &
         call define_hist_field(n_fswfac,"fswfac","1",tstr, tcstr, &
             "shortwave scaling factor",                           &
             "ratio of netsw new:old", c1, c0,                     &
             ns1, f_fswfac)

      if (f_frzmlt(1:1) /= 'x') &
         call define_hist_field(n_frzmlt,"frzmlt","W/m^2",tstr, tcstr, &
             "freeze/melt potential",                                  &
             "if >0, new ice forms; if <0, ice melts", c1, c0,         &
             ns1, f_frzmlt)

      if (f_fswabs(1:1) /= 'x') &
         call define_hist_field(n_fswabs,"fswabs","W/m^2",tstr, tcstr, &
             "snow/ice/ocn absorbed solar flux (cpl)",                 &
             "positive downward", c1, c0,                              &
             ns1, f_fswabs)
      
      if (f_fswabs_ai(1:1) /= 'x') &
         call define_hist_field(n_fswabs_ai,"fswabs_ai","W/m^2",tstr, tcstr, &
             "snow/ice/ocn absorbed solar flux",                             &
             "weighted by ice area", c1, c0,                                 &
             ns1, f_fswabs_ai)

#if (defined AEROFRC) || (defined PONDFRC) || (defined CCSM3FRC)
      if (f_fswsfc_ai(1:1) /= 'x') &
         call define_hist_field(n_fswsfc_ai,"fswsfc_ai","W/m^2",tstr, tcstr, &
             "snow/ice/ocn surface absorbed solar flux",                     &
             "weighted by ice area", c1, c0,                                 &
             ns1, f_fswsfc_ai)

      if (f_fswint_ai(1:1) /= 'x') &
         call define_hist_field(n_fswint_ai,"fswint_ai","W/m^2",tstr, tcstr, &
             "snow/ice/ocn internal absorbed solar flux",                    &
             "weighted by ice area", c1, c0,                                 &
             ns1, f_fswint_ai)
#endif
#ifdef AEROFRC
      if (f_dfswabs_noaero(1:1) /= 'x') &
         call define_hist_field(n_dfswabs_noaero,"dfswabs_noaero",       &
             "W/m^2",tstr, tcstr,                                            &
             "snow/ice/ocn diagnostic absorbed solar flux",                  &
             "weighted by ice area", c1, c0,                                 &
             ns1, f_dfswabs_noaero)

      if (f_dfswsfc_noaero(1:1) /= 'x') &
         call define_hist_field(n_dfswsfc_noaero,"dfswsfc_noaero",       &
             "W/m^2",tstr, tcstr,                                            &
             "snow/ice/ocn diagnostic surface abs solar flux",               &
             "weighted by ice area", c1, c0,                                 &
             ns1, f_dfswsfc_noaero)

      if (f_dfswint_noaero(1:1) /= 'x') &
         call define_hist_field(n_dfswint_noaero,"dfswint_noaero",       &
             "W/m^2",tstr, tcstr,                                            &
             "snow/ice/ocn diagnostic internal abs solar flux",              &
             "weighted by ice area", c1, c0,                                 &
             ns1, f_dfswint_noaero)

      if (f_dfswthru_noaero(1:1) /= 'x') &
         call define_hist_field(n_dfswthru_noaero,"dfswthru_noaero",     &
             "W/m^2",tstr, tcstr,                                            &
             "snow/ice/ocn diagnostic penetrating solar flux",               &
             "weighted by ice area", c1, c0,                                 &
             ns1, f_dfswthru_noaero)

      if (f_dalvdr_noaero(1:1) /= 'x') &
         call define_hist_field(n_dalvdr_noaero,"dalvdr_noaero",   &
             "%",tstr, tcstr,                                    &
             "diagnostic visible direct albedo",                 &
             "none", c100, c0,                                   &
             ns1, f_dalvdr_noaero)
      
      if (f_dalidr_noaero(1:1) /= 'x') &
         call define_hist_field(n_dalidr_noaero,"dalidr_noaero",   &
             "%",tstr, tcstr,                                    &
             "diagnostic infrared direct albedo",                &
             "none", c100, c0,                                   &
             ns1, f_dalidr_noaero)
      
      if (f_dalvdf_noaero(1:1) /= 'x') &
         call define_hist_field(n_dalvdf_noaero,"dalvdf_noaero",   &
             "%",tstr, tcstr,                                    &
             "diagnostic visible diffuse albedo",                &
             "none", c100, c0,                                   &
             ns1, f_dalvdf_noaero)
      
      if (f_dalidf_noaero(1:1) /= 'x') &
         call define_hist_field(n_dalidf_noaero,"dalidf_noaero",   &
             "%",tstr, tcstr,                                    &
             "diagnostic infrared diffuse albedo",               &
             "none", c100, c0,                                   &
             ns1, f_dalidf_noaero)
      
      if (f_dalbice_noaero(1:1) /= 'x') &
         call define_hist_field(n_dalbice_noaero,"dalbice_noaero",   &
             "%",tstr, tcstr,                                      &
             "diagnostic bare ice albedo",                         &
             "averaged for coszen>0, weighted by aice", c100, c0,  &
             ns1, f_dalbice_noaero)
      
      if (f_dalbsno_noaero(1:1) /= 'x') &
         call define_hist_field(n_dalbsno_noaero,"dalbsno_noaero",   &
             "%",tstr, tcstr,                                      &
             "diagnostic snow albedo",                             &
             "averaged for coszen>0, weighted by aice", c100, c0,  &
             ns1, f_dalbsno_noaero)
      
      if (f_dalbpnd_noaero(1:1) /= 'x') &
         call define_hist_field(n_dalbpnd_noaero,"dalbpnd_noaero",   &
             "%",tstr, tcstr,                                      &
             "diagnostic pond albedo",                             &
             "averaged for coszen>0, weighted by aice", c100, c0,  &
             ns1, f_dalbpnd_noaero)
#endif
#ifdef CCSM3FRC
      if (f_dfswabs_ccsm3(1:1) /= 'x') &
         call define_hist_field(n_dfswabs_ccsm3,"dfswabs_ccsm3",       &
             "W/m^2",tstr, tcstr,                                            &
             "snow/ice/ocn diagnostic absorbed solar flux",                  &
             "weighted by ice area", c1, c0,                                 &
             ns1, f_dfswabs_ccsm3)

      if (f_dfswsfc_ccsm3(1:1) /= 'x') &
         call define_hist_field(n_dfswsfc_ccsm3,"dfswsfc_ccsm3",       &
             "W/m^2",tstr, tcstr,                                            &
             "snow/ice/ocn diagnostic surface abs solar flux",               &
             "weighted by ice area", c1, c0,                                 &
             ns1, f_dfswsfc_ccsm3)

      if (f_dfswint_ccsm3(1:1) /= 'x') &
         call define_hist_field(n_dfswint_ccsm3,"dfswint_ccsm3",       &
             "W/m^2",tstr, tcstr,                                            &
             "snow/ice/ocn diagnostic internal abs solar flux",              &
             "weighted by ice area", c1, c0,                                 &
             ns1, f_dfswint_ccsm3)

      if (f_dfswthru_ccsm3(1:1) /= 'x') &
         call define_hist_field(n_dfswthru_ccsm3,"dfswthru_ccsm3",     &
             "W/m^2",tstr, tcstr,                                            &
             "snow/ice/ocn diagnostic penetrating solar flux",               &
             "weighted by ice area", c1, c0,                                 &
             ns1, f_dfswthru_ccsm3)

      if (f_dalvdr_ccsm3(1:1) /= 'x') &
         call define_hist_field(n_dalvdr_ccsm3,"dalvdr_ccsm3",   &
             "%",tstr, tcstr,                                    &
             "diagnostic visible direct albedo",                 &
             "none", c100, c0,                                   &
             ns1, f_dalvdr_ccsm3)
      
      if (f_dalidr_ccsm3(1:1) /= 'x') &
         call define_hist_field(n_dalidr_ccsm3,"dalidr_ccsm3",   &
             "%",tstr, tcstr,                                    &
             "diagnostic infrared direct albedo",                &
             "none", c100, c0,                                   &
             ns1, f_dalidr_ccsm3)
      
      if (f_dalvdf_ccsm3(1:1) /= 'x') &
         call define_hist_field(n_dalvdf_ccsm3,"dalvdf_ccsm3",   &
             "%",tstr, tcstr,                                    &
             "diagnostic visible diffuse albedo",                &
             "none", c100, c0,                                   &
             ns1, f_dalvdf_ccsm3)
      
      if (f_dalidf_ccsm3(1:1) /= 'x') &
         call define_hist_field(n_dalidf_ccsm3,"dalidf_ccsm3",   &
             "%",tstr, tcstr,                                    &
             "diagnostic infrared diffuse albedo",               &
             "none", c100, c0,                                   &
             ns1, f_dalidf_ccsm3)
      
      if (f_dalbice_ccsm3(1:1) /= 'x') &
         call define_hist_field(n_dalbice_ccsm3,"dalbice_ccsm3",   &
             "%",tstr, tcstr,                                      &
             "diagnostic bare ice albedo",                         &
             "averaged for coszen>0, weighted by aice", c100, c0,  &
             ns1, f_dalbice_ccsm3)
      
      if (f_dalbsno_ccsm3(1:1) /= 'x') &
         call define_hist_field(n_dalbsno_ccsm3,"dalbsno_ccsm3",   &
             "%",tstr, tcstr,                                      &
             "diagnostic snow albedo",                             &
             "averaged for coszen>0, weighted by aice", c100, c0,  &
             ns1, f_dalbsno_ccsm3)
#endif
#ifdef PONDFRC
      if (f_dfswabs_nopond(1:1) /= 'x') &
         call define_hist_field(n_dfswabs_nopond,"dfswabs_nopond",       &
             "W/m^2",tstr, tcstr,                                            &
             "snow/ice/ocn diagnostic absorbed solar flux",                  &
             "weighted by ice area", c1, c0,                                 &
             ns1, f_dfswabs_nopond)

      if (f_dfswsfc_nopond(1:1) /= 'x') &
         call define_hist_field(n_dfswsfc_nopond,"dfswsfc_nopond",       &
             "W/m^2",tstr, tcstr,                                            &
             "snow/ice/ocn diagnostic surface abs solar flux",               &
             "weighted by ice area", c1, c0,                                 &
             ns1, f_dfswsfc_nopond)

      if (f_dfswint_nopond(1:1) /= 'x') &
         call define_hist_field(n_dfswint_nopond,"dfswint_nopond",       &
             "W/m^2",tstr, tcstr,                                            &
             "snow/ice/ocn diagnostic internal abs solar flux",              &
             "weighted by ice area", c1, c0,                                 &
             ns1, f_dfswint_nopond)

      if (f_dfswthru_nopond(1:1) /= 'x') &
         call define_hist_field(n_dfswthru_nopond,"dfswthru_nopond",     &
             "W/m^2",tstr, tcstr,                                            &
             "snow/ice/ocn diagnostic penetrating solar flux",               &
             "weighted by ice area", c1, c0,                                 &
             ns1, f_dfswthru_nopond)

      if (f_dalvdr_nopond(1:1) /= 'x') &
         call define_hist_field(n_dalvdr_nopond,"dalvdr_nopond",   &
             "%",tstr, tcstr,                                    &
             "diagnostic visible direct albedo",                 &
             "none", c100, c0,                                   &
             ns1, f_dalvdr_nopond)
      
      if (f_dalidr_nopond(1:1) /= 'x') &
         call define_hist_field(n_dalidr_nopond,"dalidr_nopond",   &
             "%",tstr, tcstr,                                    &
             "diagnostic infrared direct albedo",                &
             "none", c100, c0,                                   &
             ns1, f_dalidr_nopond)
      
      if (f_dalvdf_nopond(1:1) /= 'x') &
         call define_hist_field(n_dalvdf_nopond,"dalvdf_nopond",   &
             "%",tstr, tcstr,                                    &
             "diagnostic visible diffuse albedo",                &
             "none", c100, c0,                                   &
             ns1, f_dalvdf_nopond)
      
      if (f_dalidf_nopond(1:1) /= 'x') &
         call define_hist_field(n_dalidf_nopond,"dalidf_nopond",   &
             "%",tstr, tcstr,                                    &
             "diagnostic infrared diffuse albedo",               &
             "none", c100, c0,                                   &
             ns1, f_dalidf_nopond)
      
      if (f_dalbice_nopond(1:1) /= 'x') &
         call define_hist_field(n_dalbice_nopond,"dalbice_nopond",   &
             "%",tstr, tcstr,                                      &
             "diagnostic bare ice albedo",                         &
             "averaged for coszen>0, weighted by aice", c100, c0,  &
             ns1, f_dalbice_nopond)
      
      if (f_dalbsno_nopond(1:1) /= 'x') &
         call define_hist_field(n_dalbsno_nopond,"dalbsno_nopond",   &
             "%",tstr, tcstr,                                      &
             "diagnostic snow albedo",                             &
             "averaged for coszen>0, weighted by aice", c100, c0,  &
             ns1, f_dalbsno_nopond)
      
      if (f_dalbpnd_nopond(1:1) /= 'x') &
         call define_hist_field(n_dalbpnd_nopond,"dalbpnd_nopond",   &
             "%",tstr, tcstr,                                      &
             "diagnostic pond albedo",                             &
             "averaged for coszen>0, weighted by aice", c100, c0,  &
             ns1, f_dalbpnd_nopond)
#endif
      if (f_alvdr(1:1) /= 'x') &
         call define_hist_field(n_alvdr,"alvdr","%",tstr, tcstr, &
             "visible direct albedo",                            &
             "none", c100, c0,                                   &
             ns1, f_alvdr)
      
      if (f_alidr(1:1) /= 'x') &
         call define_hist_field(n_alidr,"alidr","%",tstr, tcstr, &
             "near IR direct albedo",                            &
             "none", c100, c0,                                   &
             ns1, f_alidr)
      
      if (f_alvdf(1:1) /= 'x') &
         call define_hist_field(n_alvdf,"alvdf","%",tstr, tcstr, &
             "visible diffuse albedo",                           &
             "none", c100, c0,                                   &
             ns1, f_alvdf)
      
      if (f_alidf(1:1) /= 'x') &
         call define_hist_field(n_alidf,"alidf","%",tstr, tcstr, &
             "near IR diffuse albedo",                           &
             "none", c100, c0,                                   &
             ns1, f_alidf)
      
      if (f_albice(1:1) /= 'x') &
         call define_hist_field(n_albice,"albice","%",tstr, tcstr, &
             "bare ice albedo",                                    &
             "averaged for coszen>0, weighted by aice", c100, c0,  &
             ns1, f_albice)
      
      if (f_albsno(1:1) /= 'x') &
         call define_hist_field(n_albsno,"albsno","%",tstr, tcstr, &
             "snow albedo",                                        &
             "averaged for coszen>0, weighted by aice", c100, c0,  &
             ns1, f_albsno)
      
      if (f_albpnd(1:1) /= 'x') &
         call define_hist_field(n_albpnd,"albpnd","%",tstr, tcstr, &
             "melt pond albedo",                                   &
             "averaged for coszen>0, weighted by aice", c100, c0,  &
             ns1, f_albpnd)
      
      if (f_coszen(1:1) /= 'x') &
         call define_hist_field(n_coszen,"coszen","radian",tstr, tcstr, &
             "cosine of the zenith angle",                              &
             "negative below horizon", c1, c0,                          &
             ns1, f_coszen)

      if (f_flat(1:1) /= 'x') &
         call define_hist_field(n_flat,"flat","W/m^2",tstr, tcstr, &
             "latent heat flux (cpl)",                             &
             "positive downward", c1, c0,                          &
             ns1, f_flat)
      
      if (f_flat_ai(1:1) /= 'x') &
         call define_hist_field(n_flat_ai,"flat_ai","W/m^2",tstr, tcstr, &
             "latent heat flux",                                         &
             "weighted by ice area", c1, c0,                             &
             ns1, f_flat_ai)
      
      if (f_fsens(1:1) /= 'x') &
         call define_hist_field(n_fsens,"fsens","W/m^2",tstr, tcstr, &
             "sensible heat flux (cpl)",                             &
             "positive downward", c1, c0,                            &
             ns1, f_fsens)
      
      if (f_fsens_ai(1:1) /= 'x') &
         call define_hist_field(n_fsens_ai,"fsens_ai","W/m^2",tstr, tcstr, &
             "sensible heat flux",                                         &
             "weighted by ice area", c1, c0,                               &
             ns1, f_fsens_ai)
      
      if (f_flwup(1:1) /= 'x') &
         call define_hist_field(n_flwup,"flwup","W/m^2",tstr, tcstr, &
             "upward longwave flux (cpl)",                           &
             "positive downward", c1, c0,                            &
             ns1, f_flwup)
      
      if (f_flwup_ai(1:1) /= 'x') &
         call define_hist_field(n_flwup_ai,"flwup_ai","W/m^2",tstr, tcstr, &
             "upward longwave flux",                                       &
             "weighted by ice area", c1, c0,                               &
             ns1, f_flwup_ai)
      
      if (f_evap(1:1) /= 'x') &
         call define_hist_field(n_evap,"evap","cm/day",tstr, tcstr, &
             "evaporative water flux (cpl)",                        &
             "none", mps_to_cmpdy/rhofresh, c0,                     &
             ns1, f_evap)
      
      if (f_evap_ai(1:1) /= 'x') &
         call define_hist_field(n_evap_ai,"evap_ai","cm/day",tstr, tcstr, &
             "evaporative water flux",                                    &
             "weighted by ice area", mps_to_cmpdy/rhofresh, c0,           &
             ns1, f_evap_ai)
      
      if (f_Tair(1:1) /= 'x') &
         call define_hist_field(n_Tair,"Tair","C",tstr, tcstr, &
            "air temperature",                                &
             "none", c1, -tffresh,                             &
             ns1, f_Tair)
      
      if (f_Tref(1:1) /= 'x') &
         call define_hist_field(n_Tref,"Tref","C",tstr, tcstr, &
             "2m reference temperature",                       &
             "none", c1, c0,                                   &
             ns1, f_Tref)
      
      if (f_Qref(1:1) /= 'x') &
         call define_hist_field(n_Qref,"Qref","g/kg",tstr, tcstr, &
             "2m reference specific humidity",                    &
             "none", kg_to_g, c0,                                 &
             ns1, f_Qref)
      
      if (f_congel(1:1) /= 'x') &
         call define_hist_field(n_congel,"congel","cm/day",tstr, tcstr, &
             "congelation ice growth",                                  &
             "none", mps_to_cmpdy/dt, c0,                               &
             ns1, f_congel)
      
      if (f_frazil(1:1) /= 'x') &
         call define_hist_field(n_frazil,"frazil","cm/day",tstr, tcstr, &
             "frazil ice growth",                                       &
             "none", mps_to_cmpdy/dt, c0,                               &
             ns1, f_frazil)
      
      if (f_snoice(1:1) /= 'x') &
         call define_hist_field(n_snoice,"snoice","cm/day",tstr, tcstr, &
             "snow-ice formation",                                      &
             "none", mps_to_cmpdy/dt, c0,                               &
             ns1, f_snoice)
      
      if (f_meltt(1:1) /= 'x') &
         call define_hist_field(n_meltt,"meltt","cm/day",tstr, tcstr, &
             "top ice melt",                                          &
             "none", mps_to_cmpdy/dt, c0,                             &
             ns1, f_meltt)
      
      if (f_meltb(1:1) /= 'x') &
         call define_hist_field(n_meltb,"meltb","cm/day",tstr, tcstr, &
             "basal ice melt",                                        &
             "none", mps_to_cmpdy/dt, c0,                             &
             ns1, f_meltb)
      
      if (f_meltl(1:1) /= 'x') &
         call define_hist_field(n_meltl,"meltl","cm/day",tstr, tcstr, &
             "lateral ice melt",                                      &
             "none", mps_to_cmpdy/dt, c0,                             &
             ns1, f_meltl)
      
      if (f_melts(1:1) /= 'x') &
         call define_hist_field(n_melts,"melts","cm/day",tstr, tcstr, &
             "snow melt",                                             &
             "none", mps_to_cmpdy/dt, c0,                             &
             ns1, f_melts)
      
      if (f_fresh(1:1) /= 'x') &
         call define_hist_field(n_fresh,"fresh","cm/day",tstr, tcstr,   &
             "freshwtr flx ice to ocn (cpl)",                           &
             "if positive, ocean gains fresh water",                    &
             mps_to_cmpdy/rhofresh, c0,                                 &
             ns1, f_fresh)
      
      if (f_fresh_ai(1:1) /= 'x') &
         call define_hist_field(n_fresh_ai,"fresh_ai","cm/day",tstr, tcstr, &
             "freshwtr flx ice to ocn",                                     &
             "weighted by ice area", mps_to_cmpdy/rhofresh, c0,             &
             ns1, f_fresh_ai)
      
      if (f_fsalt(1:1) /= 'x') &
         call define_hist_field(n_fsalt,"fsalt","kg/m^2/day",tstr, tcstr, &
             "salt flux ice to ocn (cpl)",                                &
             "if positive, ocean gains salt", secday, c0,                 &
             ns1, f_fsalt)
      
      if (f_fsalt_ai(1:1) /= 'x') &
         call define_hist_field(n_fsalt_ai,"fsalt_ai","kg/m^2/day",tstr, tcstr,&
             "salt flux ice to ocean",                                         &
             "weighted by ice area", secday, c0,                               &
             ns1, f_fsalt_ai)
      
      if (f_fhocn(1:1) /= 'x') &
         call define_hist_field(n_fhocn,"fhocn","W/m^2",tstr, tcstr, &
             "heat flux ice to ocn (cpl)",                           &
             "if positive, ocean gains heat", c1, c0,                &
             ns1, f_fhocn)
      
      if (f_fhocn_ai(1:1) /= 'x') &
         call define_hist_field(n_fhocn_ai,"fhocn_ai","W/m^2",tstr, tcstr, &
             "heat flux ice to ocean",                                     &
             "weighted by ice area", c1, c0,                               &
             ns1, f_fhocn_ai)
      
      if (f_fswthru(1:1) /= 'x') &
         call define_hist_field(n_fswthru,"fswthru","W/m^2",tstr, tcstr, &
             "SW thru ice to ocean (cpl)",                               &
             "if positive, ocean gains heat", c1, c0,                    &
             ns1, f_fswthru)
      
      if (f_fswthru_ai(1:1) /= 'x') &
         call define_hist_field(n_fswthru_ai,"fswthru_ai","W/m^2",tstr, tcstr,&
             "SW flux thru ice to ocean",                                     &
             "weighted by ice area", c1, c0,                                  &
             ns1, f_fswthru_ai)
      
      if (f_strairx(1:1) /= 'x') &
         call define_hist_field(n_strairx,"strairx","N/m^2",ustr, ucstr, &
             "atm/ice stress (x)",                                       &
             "positive is x direction on U grid", c1, c0,                &
             ns1, f_strairx)
      
      if (f_strairy(1:1) /= 'x') &
         call define_hist_field(n_strairy,"strairy","N/m^2",ustr, ucstr, &
             "atm/ice stress (y)",                                       &
             "positive is y direction on U grid", c1, c0,                &
             ns1, f_strairy)
      
      if (f_strtltx(1:1) /= 'x') &
         call define_hist_field(n_strtltx,"strtltx","N/m^2",ustr, ucstr, &
             "sea sfc tilt stress (x)",                                  &
             "positive is x direction on U grid", c1, c0,                &
             ns1, f_strtltx)
      
      if (f_strtlty(1:1) /= 'x') &
         call define_hist_field(n_strtlty,"strtlty","N/m^2",ustr, ucstr, &
             "sea sfc tilt stress (y)",                                  &
             "positive is y direction on U grid", c1, c0,                &
             ns1, f_strtlty)
      
      if (f_strcorx(1:1) /= 'x') &
         call define_hist_field(n_strcorx,"strcorx","N/m^2",ustr, ucstr, &
             "coriolis stress (x)",                                      &
             "positive is x direction on U grid", c1, c0,                &
             ns1, f_strcorx)
      
      if (f_strcory(1:1) /= 'x') &
         call define_hist_field(n_strcory,"strcory","N/m^2",ustr, ucstr, &
             "coriolis stress (y)",                                      &
             "positive is y direction on U grid", c1, c0,                &
             ns1, f_strcory)
      
      if (f_strocnx(1:1) /= 'x') &
         call define_hist_field(n_strocnx,"strocnx","N/m^2",ustr, ucstr, &
             "ocean/ice stress (x)",                                     &
             "positive is x direction on U grid", c1, c0,                &
             ns1, f_strocnx)
      
      if (f_strocny(1:1) /= 'x') &
         call define_hist_field(n_strocny,"strocny","N/m^2",ustr, ucstr, &
             "ocean/ice stress (y)",                                     &
             "positive is y direction on U grid", c1, c0,                &
             ns1, f_strocny)
      
      if (f_strintx(1:1) /= 'x') &
         call define_hist_field(n_strintx,"strintx","N/m^2",ustr, ucstr, &
             "internal ice stress (x)",                                  &
             "positive is x direction on U grid", c1, c0,                &
             ns1, f_strintx)
      
      if (f_strinty(1:1) /= 'x') &
         call define_hist_field(n_strinty,"strinty","N/m^2",ustr, ucstr, &
             "internal ice stress (y)",                                  &
             "positive is y direction on U grid", c1, c0,                &
             ns1, f_strinty)
      
      if (f_strength(1:1) /= 'x') &
         call define_hist_field(n_strength,"strength","N/m",tstr, tcstr, &
             "compressive ice strength",                                 &
             "none", c1, c0,                                             &
             ns1, f_strength)
      
      if (f_opening(1:1) /= 'x') &
         call define_hist_field(n_opening,"opening","%/day",tstr, tcstr, &
             "lead area opening rate",                                   &
             "none", secday*c100, c0,                                    &
             ns1, f_opening)
      
      if (f_divu(1:1) /= 'x') &
         call define_hist_field(n_divu,"divu","%/day",tstr, tcstr, &
             "strain rate (divergence)",                           &
             "none", secday*c100, c0,                              &
             ns1, f_divu)
      
      if (f_shear(1:1) /= 'x') &
         call define_hist_field(n_shear,"shear","%/day",tstr, tcstr, &
             "strain rate (shear)",                                  &
             "none", secday*c100, c0,                                &
             ns1, f_shear)
      
      if (f_sig1(1:1) /= 'x') &
         call define_hist_field(n_sig1,"sig1"," ",ustr, ucstr, &
             "norm. principal stress 1",                       &
             "sig1 is instantaneous", c1, c0,                  &
             ns1, f_sig1)
      
      if (f_sig2(1:1) /= 'x') &
         call define_hist_field(n_sig2,"sig2"," ",ustr, ucstr, &
             "norm. principal stress 2",                       &
             "sig2 is instantaneous", c1, c0,                  &
             ns1, f_sig2)
      
      if (f_dvidtt(1:1) /= 'x') &
         call define_hist_field(n_dvidtt,"dvidtt","cm/day",tstr, tcstr, &
             "volume tendency thermo",                                  &
             "none", mps_to_cmpdy, c0,                                  &
             ns1, f_dvidtt)
      
      if (f_dvidtd(1:1) /= 'x') &
         call define_hist_field(n_dvidtd,"dvidtd","cm/day",tstr, tcstr, &
             "volume tendency dynamics",                                &
             "none", mps_to_cmpdy, c0,                                  &
             ns1, f_dvidtd)
      
      if (f_daidtt(1:1) /= 'x') &
         call define_hist_field(n_daidtt,"daidtt","%/day",tstr, tcstr, &
             "area tendency thermo",                                   &
             "none", secday*c100, c0,                                  &
             ns1, f_daidtt)
      
      if (f_daidtd(1:1) /= 'x') &
         call define_hist_field(n_daidtd,"daidtd","%/day",tstr, tcstr, &
             "area tendency dynamics",                                 &
             "none", secday*c100, c0,                                  &
             ns1, f_daidtd)
      
      if (f_mlt_onset(1:1) /= 'x') &
         call define_hist_field(n_mlt_onset,"mlt_onset","day of year", &
             tstr, tcstr,"melt onset date",                            &
             "midyear restart gives erroneous dates", c1, c0,          &
             ns1, f_mlt_onset)
      
      if (f_frz_onset(1:1) /= 'x') &
         call define_hist_field(n_frz_onset,"frz_onset","day of year", &
             tstr, tcstr,"freeze onset date",                          &
             "midyear restart gives erroneous dates", c1, c0,          &
             ns1, f_frz_onset)
      
      if (f_dardg1dt(1:1) /= 'x') &
         call define_hist_field(n_dardg1dt,"dardg1dt","%/day",tstr, tcstr, &
             "ice area ridging rate",                                      &
             "none", secday*c100, c0,                                      &
             ns1, f_dardg1dt)
      
      if (f_dardg2dt(1:1) /= 'x') &
         call define_hist_field(n_dardg2dt,"dardg2dt","%/day",tstr, tcstr, &
             "ridge area formation rate",                                  &
             "none", secday*c100, c0,                                      &
             ns1, f_dardg2dt)
      
      if (f_dvirdgdt(1:1) /= 'x') &
         call define_hist_field(n_dvirdgdt,"dvirdgdt","cm/day",tstr, tcstr, &
             "ice volume ridging rate",                                     &
             "none", mps_to_cmpdy, c0,                                      &
             ns1, f_dvirdgdt)
      
      if (f_hisnap(1:1) /= 'x') &
         call define_hist_field(n_hisnap,"hisnap","m",tstr, tcstr, &
             "ice volume snapshot",                                &
             "none", c1, c0,                                       &
             ns1, f_hisnap)
      
      if (f_aisnap(1:1) /= 'x') &
         call define_hist_field(n_aisnap,"aisnap"," ",tstr, tcstr, &
             "ice area snapshot",                                  &
             "none", c1, c0,                                       &
             ns1, f_aisnap)
      
      if (f_trsig(1:1) /= 'x') &
         call define_hist_field(n_trsig,"trsig","N/m^2",tstr, tcstr, &
             "internal stress tensor trace",                         &
             "ice strength approximation", c1, c0,                   &
             ns1, f_trsig)
      
      if (f_icepresent(1:1) /= 'x') &
         call define_hist_field(n_icepresent,"ice_present","1",tstr, tcstr, &
             "fraction of time-avg interval that any ice is present",       &
             "ice extent flag", c1, c0,                                     &
             ns1, f_icepresent)
      
      if (f_fsurf_ai(1:1) /= 'x') &
         call define_hist_field(n_fsurf_ai,"fsurf_ai","W/m^2",tstr, tcstr, &
             "net surface heat flux",                                      &
             "positive downward, excludes conductive flux, weighted by ice area", &
             c1, c0, ns1, f_fsurf_ai)
      
      if (f_fcondtop_ai(1:1) /= 'x') &
         call define_hist_field(n_fcondtop_ai,"fcondtop_ai","W/m^2", &
             tstr, tcstr,"top surface conductive heat flux",         &
             "positive downwards, weighted by ice area", c1, c0,     &
             ns1, f_fcondtop_ai)
      
      if (f_fmeltt_ai(1:1) /= 'x') &
         call define_hist_field(n_fmeltt_ai,"fmeltt_ai","W/m^2",tstr, tcstr, &
             "net surface heat flux causing melt",                           &
             "always >= 0, weighted by ice area", c1, c0,                    &
             ns1, f_fmeltt_ai)
      
      ! Category variables
      if (f_aicen(1:1) /= 'x') then
         do n=1,ncat_hist
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'aicen', trim(nchar) ! aicen
            stmp = 'ice area, category '   ! aicen
            write(vdesc_in,'(a,2x,a)') trim(stmp), trim(nchar)
            call define_hist_field(n_aicen(n,:),vname_in,"%",tstr, tcstr, &
                vdesc_in, "Ice range:", c100, c0,                       &
                ns1, f_aicen)
         enddo
      endif

      if (f_vicen(1:1) /= 'x') then
         do n=1,ncat_hist
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'vicen', trim(nchar) ! vicen
            stmp = 'ice volume, category '   ! vicen
            write(vdesc_in,'(a,2x,a)') trim(stmp), trim(nchar)
            call define_hist_field(n_vicen(n,:),vname_in,"m",tstr, tcstr, &
                vdesc_in, "none", c1, c0,                               &
                ns1, f_vicen)
         enddo
      endif

      if (f_fsurfn_ai(1:1) /= 'x') then
         do n=1,ncat_hist
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'fsurfn_ai', trim(nchar) ! fsurfn_ai
            stmp = 'net surface heat flux, category '   ! fsurfn_ai
            write(vdesc_in,'(a,2x,a)') trim(stmp), trim(nchar)
            call define_hist_field(n_fsurfn_ai(n,:),vname_in,"W/m^2",    &
                tstr, tcstr, vdesc_in, "weighted by ice area", c1, c0, &
                ns1, f_fsurfn_ai)
         enddo
      endif

      if (f_fcondtopn_ai(1:1) /= 'x') then
         do n=1,ncat_hist
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'fcondtopn_ai', trim(nchar) ! fcondtopn_ai
            stmp = 'top sfc conductive heat flux, cat '   ! fcondtopn_ai
            write(vdesc_in,'(a,2x,a)') trim(stmp), trim(nchar)
            call define_hist_field(n_fcondtopn_ai(n,:),vname_in,"W/m^2", &
                tstr, tcstr, vdesc_in, "weighted by ice area", c1, c0, &
                ns1, f_fcondtopn_ai)
         enddo
      endif

      if (f_fmelttn_ai(1:1) /= 'x') then
         do n=1,ncat_hist
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'fmelttn_ai', trim(nchar) ! fmelttn_ai
            stmp = 'net sfc heat flux causing melt, cat '   ! fmelttn_ai
            write(vdesc_in,'(a,2x,a)') trim(stmp), trim(nchar)
            call define_hist_field(n_fmelttn_ai(n,:),vname_in,"W/m^2",   &
                tstr, tcstr, vdesc_in, "weighted by ice area", c1, c0, &
                ns1, f_fmelttn_ai)
         enddo
      endif

      if (f_flatn_ai(1:1) /= 'x') then
         do n=1,ncat_hist
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'flatn_ai', trim(nchar) ! flatn_ai
            stmp = 'latent heat flux, category '   ! flatn_ai
            write(vdesc_in,'(a,2x,a)') trim(stmp), trim(nchar)
            call define_hist_field(n_flatn_ai(n,:),vname_in,"W/m^2", &
                tstr, tcstr, vdesc_in, "weighted by ice area", c1, c0, &
                ns1, f_flatn_ai)
         enddo
      endif

      ! Tracers

      ! Ice Age
      if (f_iage(1:1) /= 'x') &
         call define_hist_field(n_iage,"iage","years",tstr, tcstr, &
             "sea ice age",                                        &
             "none", c1/(secday*days_per_year), c0,                &
              ns1, f_iage)

      ! FY Ice Concentration
      if (f_FY(1:1) /= 'x') &
         call define_hist_field(n_FY,"FYarea"," ",tstr, tcstr, &
             "first-year ice area",                            &
             "weighted by ice area", c1, c0,                   &
              ns1, f_FY)

      ! Aerosols
      if (f_aero(1:1) /= 'x') then
         do n=1,n_aero
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'aerosnossl', trim(nchar)
            call define_hist_field(n_aerosn1(n,:),vname_in,"kg/kg", &
                tstr, tcstr,"snow ssl aerosol mass","none", c1, c0, &
                ns1, f_aero)
            write(vname_in,'(a,a)') 'aerosnoint', trim(nchar)
            call define_hist_field(n_aerosn2(n,:),vname_in,"kg/kg", &
                tstr, tcstr,"snow int aerosol mass","none", c1, c0, &
                ns1, f_aero)
            write(vname_in,'(a,a)') 'aeroicessl', trim(nchar)
            call define_hist_field(n_aeroic1(n,:),vname_in,"kg/kg", &
                tstr, tcstr,"ice ssl aerosol mass","none", c1, c0,  &
                ns1, f_aero)
            write(vname_in,'(a,a)') 'aeroiceint', trim(nchar)
            call define_hist_field(n_aeroic2(n,:),vname_in,"kg/kg", &
                tstr, tcstr,"ice int aerosol mass","none", c1, c0,  &
                ns1, f_aero)
         enddo
      endif

      if (f_faero_atm(1:1) /= 'x') then
         do n=1,n_aero
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'faero_atm', trim(nchar)
            call define_hist_field(n_faero_atm(n,:),vname_in,"kg/m^2 s", &
                tstr, tcstr,"aerosol deposition rate","none", c1, c0,    &
                ns1, f_faero_atm)
         enddo
      endif

      if (f_faero_ocn(1:1) /= 'x') then
         do n=1,n_aero
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'faero_ocn', trim(nchar)
            call define_hist_field(n_faero_ocn(n,:),vname_in,"kg/m^2 s", &
                tstr, tcstr,"aerosol flux to ocean","none", c1, c0,      &
                ns1, f_faero_ocn)
         enddo
      endif

      ! Level and Ridged ice
      if (f_alvl(1:1) /= 'x') &
         call define_hist_field(n_alvl,"alvl","1",tstr, tcstr, &
             "level ice area fraction",                            &
             "none", c1, c0,                                       &
             ns1, f_alvl)
      if (f_vlvl(1:1) /= 'x') &
         call define_hist_field(n_vlvl,"vlvl","m",tstr, tcstr, &
             "level ice mean thickness",                           &
             "none", c1, c0,                                       &
             ns1, f_vlvl)
      if (f_ardg(1:1) /= 'x') &
         call define_hist_field(n_ardg,"ardg","1",tstr, tcstr, &
             "ridged ice area fraction",                           &
             "none", c1, c0,                                       &
             ns1, f_ardg)
      if (f_vrdg(1:1) /= 'x') &
         call define_hist_field(n_vrdg,"vrdg","m",tstr, tcstr, &
             "ridged ice mean thickness",                          &
             "none", c1, c0,                                       &
             ns1, f_vrdg)

      ! Melt ponds
      if (f_apond(1:1) /= 'x') &
         call define_hist_field(n_apond,"apond","%",tstr, tcstr, &
             "melt pond concentration",                          &
             "none", c100, c0,                                   &
             ns1, f_apond)

      if (f_apondn(1:1) /= 'x') then
         do n=1,ncat_hist
            write(nchar,'(i3.3)') n
            write(vname_in,'(a,a)') 'apond', trim(nchar) ! apondn
            stmp = 'melt pond concentration, category '   ! apondn
            write(vdesc_in,'(a,2x,a)') trim(stmp), trim(nchar)
            call define_hist_field(n_apondn(n,:),vname_in,"%", &
                tstr, tcstr, vdesc_in, "none", c100, c0,       &
                ns1, f_apondn)
         enddo
      endif

      enddo ! ns1

      allocate(aa(nx_block,ny_block,num_avail_hist_fields,max_blocks))

      !-----------------------------------------------------------------
      ! fill igrd array with namelist values
      !-----------------------------------------------------------------

      igrd=.true.

      igrd(n_tmask     ) = f_tmask
      igrd(n_blkmask   ) = f_blkmask
      igrd(n_tarea     ) = f_tarea
      igrd(n_uarea     ) = f_uarea
      igrd(n_dxt       ) = f_dxt
      igrd(n_dyt       ) = f_dyt
      igrd(n_dxu       ) = f_dxu
      igrd(n_dyu       ) = f_dyu
      igrd(n_HTN       ) = f_HTN
      igrd(n_HTE       ) = f_HTE
      igrd(n_ANGLE     ) = f_ANGLE
      igrd(n_ANGLET    ) = f_ANGLET

      ntmp(:) = 0
      if (my_task == master_task) then
        write(nu_diag,*) ' '
        write(nu_diag,*) 'The following variables will be ', &
                         'written to the history tape: '
        write(nu_diag,101) 'description','units','variable','frequency','x'
        do n=1,num_avail_hist_fields
           if (avail_hist_fields(n)%vhistfreq_n /= 0) &
           write(nu_diag,100) avail_hist_fields(n)%vdesc, &
              avail_hist_fields(n)%vunit, avail_hist_fields(n)%vname, &
              avail_hist_fields(n)%vhistfreq,avail_hist_fields(n)%vhistfreq_n
           do ns = 1, nstreams
              if (avail_hist_fields(n)%vhistfreq == histfreq(ns)) &
                 ntmp(ns)=ntmp(ns)+1
           enddo
        enddo ! num_avail_hist_fields
        write(nu_diag,*) ' '
      endif
  100 format (1x,a40,2x,a12,2x,a12,1x,a1,2x,i6)
  101 format (2x,a19,10x,a12,9x,a12,2x,a,3x,a1)

      call broadcast_array(ntmp, master_task)
      do ns = 1, nstreams
         if (ntmp(ns)==0) histfreq_n(ns) = 0
      enddo

      !-----------------------------------------------------------------
      ! initialize the history arrays
      !-----------------------------------------------------------------
      aa(:,:,:,:) = c0
      avgct(:) = c0
      albcnt(:,:,:,:) = c0

      if (restart .and. yday >= c2) then
! restarting midyear gives erroneous onset dates
         mlt_onset = 999._dbl_kind 
         frz_onset = 999._dbl_kind 
      else
         mlt_onset = c0
         frz_onset = c0
      endif

      end subroutine init_hist

!=======================================================================
!
!BOP
!
! !IROUTINE: ice_write_hist - write average ice quantities or snapshots
!
! !INTERFACE:
!
      subroutine ice_write_hist (dt)
!
! !DESCRIPTION:
!
! write average ice quantities or snapshots
!
! !REVISION HISTORY:
!
! author:   Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_blocks
      use ice_domain
      use ice_calendar, only: new_year, secday, yday, write_history, &
                              write_ic, time, histfreq, nstreams
      use ice_state
      use ice_constants
      use ice_flux
      use ice_dyn_evp
      use ice_itd, only: ilyr1, slyr1
      use ice_timers
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step
!EOP
!
      integer (kind=int_kind) :: &
           i,j,k,n,nct,ns   , &
           iblk             , & ! block index
           ilo,ihi,jlo,jhi  , & ! beginning and end of physical domain
           nstrm                ! nstreams (1 if writing initial condition)

      real (kind=dbl_kind) :: &
           ravgct            ,& ! 1/avgct
           ravgctz              ! 1/avgct

      type (block) :: &
         this_block             ! block information for current block

      real (kind=dbl_kind) :: &
           worka(nx_block,ny_block), &
           workb(nx_block,ny_block)

      real (kind=dbl_kind) :: &
           ai                ,& ! aice_init
           ain               ,& ! aicen_init
           hs                ,& ! snow depth
           qs                ,& ! snow heat content
           qi                ,& ! snow heat content
           uee, vnn          ,& ! velocity component on east and north edge
           hiee, hinn           ! ice volume on east and north edge

      !---------------------------------------------------------------
      ! increment step counter
      !---------------------------------------------------------------

      do ns=1,nstreams
         if (.not. hist_avg .or. histfreq(ns) == '1') then  ! write snapshots
           do k=1,num_avail_hist_fields
              if (avail_hist_fields(k)%vhistfreq == histfreq(ns)) &
                  aa(:,:,k,:) = c0
           enddo
           avgct(ns) = c1
         else                      ! write averages over time histfreq
           avgct(ns) = avgct(ns) + c1
           if (avgct(ns) == c1) time_beg(ns) = (time-dt)/int(secday)
         endif
      enddo

      !---------------------------------------------------------------
      ! increment field
      !---------------------------------------------------------------

     !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,worka,workb,i,j,k,n,nct,hs)
      do iblk = 1, nblocks

         workb(:,:) = aice_init(:,:,iblk)

         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         if (f_hi(1:1) /= 'x') &
            call accum_hist_field(n_hi,   iblk, vice(:,:,iblk))
         if (f_hs(1:1) /= 'x') &
            call accum_hist_field(n_hs,   iblk, vsno(:,:,iblk))
         if (f_Tsfc(1:1) /= 'x') &
            call accum_hist_field(n_Tsfc, iblk, trcr(:,:,nt_Tsfc,iblk))
         if (f_aice(1:1) /= 'x') &
            call accum_hist_field(n_aice, iblk, aice(:,:,iblk))
         if (f_uvel(1:1) /= 'x') &
            call accum_hist_field(n_uvel, iblk, uvel(:,:,iblk))
         if (f_vvel(1:1) /= 'x') &
            call accum_hist_field(n_vvel, iblk, vvel(:,:,iblk))
         if (f_transix(1:1) /= 'x') then
            worka(:,:) = c0
            do j = jlo, jhi
            do i = ilo, ihi
               uee = p5*(uvel(i,j,iblk)+uvel(i,j-1,iblk))
               hiee = p5*(vice(i,j,iblk)+vice(i+1,j,iblk))
               worka(i,j) = uee*hiee*HTE(i,j,iblk)*rhoi
            enddo
            enddo
            call accum_hist_field(n_transix, iblk, worka(:,:))
         endif
         if (f_transiy(1:1) /= 'x') then
            worka(:,:) = c0
            do j = jlo, jhi
            do i = ilo, ihi
               vnn = p5*(vvel(i,j,iblk)+vvel(i-1,j,iblk))
               hinn = p5*(vice(i,j,iblk)+vice(i,j+1,iblk))
               worka(i,j) = vnn*hinn*HTN(i,j,iblk)*rhoi
            enddo
            enddo
            call accum_hist_field(n_transiy, iblk, worka(:,:))
         endif

         if (f_fswdn(1:1) /= 'x') &
            call accum_hist_field(n_fswdn, iblk, fsw(:,:,iblk))
         if (f_fswup(1:1) /= 'x') &
            call accum_hist_field(n_fswup, iblk, &
                 (fsw(:,:,iblk)-fswabs(:,:,iblk)*workb(:,:)))
         if (f_flwdn(1:1) /= 'x') &
            call accum_hist_field(n_flwdn, iblk, flw(:,:,iblk))
         if (f_snow(1:1) /= 'x') &
            call accum_hist_field(n_snow,  iblk, fsnow(:,:,iblk))
         if (f_snow_ai(1:1) /= 'x') &
            call accum_hist_field(n_snow_ai, iblk, fsnow(:,:,iblk)*workb(:,:))
         if (f_rain(1:1) /= 'x') &
            call accum_hist_field(n_rain,  iblk, frain(:,:,iblk))
         if (f_rain_ai(1:1) /= 'x') &
            call accum_hist_field(n_rain_ai, iblk, frain(:,:,iblk)*workb(:,:))

         if (f_sst(1:1) /= 'x') &
            call accum_hist_field(n_sst,  iblk, sst(:,:,iblk))
         if (f_sss(1:1) /= 'x') &
            call accum_hist_field(n_sss,  iblk, sss(:,:,iblk))
         if (f_uocn(1:1) /= 'x') &
            call accum_hist_field(n_uocn, iblk, uocn(:,:,iblk))
         if (f_vocn(1:1) /= 'x') &
            call accum_hist_field(n_vocn, iblk, vocn(:,:,iblk))
         if (f_frzmlt(1:1) /= 'x') &
            call accum_hist_field(n_frzmlt,  iblk, frzmlt(:,:,iblk))

         if (f_fswfac(1:1) /= 'x') &
            call accum_hist_field(n_fswfac,  iblk, fswfac(:,:,iblk))
         if (f_fswabs(1:1) /= 'x') &
            call accum_hist_field(n_fswabs,  iblk, fswabs(:,:,iblk))
         if (f_fswabs_ai(1:1) /= 'x') &
            call accum_hist_field(n_fswabs_ai,iblk,fswabs(:,:,iblk)*workb(:,:))

#if (defined AEROFRC) || (defined PONDFRC) || (defined CCSM3FRC)
         if (f_fswsfc_ai(1:1) /= 'x') &
            call accum_hist_field(n_fswsfc_ai,iblk,fswsfc(:,:,iblk))
         if (f_fswint_ai(1:1) /= 'x') &
            call accum_hist_field(n_fswint_ai,iblk,fswint(:,:,iblk))
#endif
#ifdef AEROFRC
         if (f_dfswabs_noaero(1:1) /= 'x') &
            call accum_hist_field(n_dfswabs_noaero,iblk, &
                                    dfswabs_noaero(:,:,iblk))
         if (f_dfswsfc_noaero(1:1) /= 'x') &
            call accum_hist_field(n_dfswsfc_noaero,iblk, &
                                    dfswsfc_noaero(:,:,iblk))
         if (f_dfswint_noaero(1:1) /= 'x') &
            call accum_hist_field(n_dfswint_noaero,iblk, &
                                    dfswint_noaero(:,:,iblk))
         if (f_dfswthru_noaero(1:1) /= 'x') &
            call accum_hist_field(n_dfswthru_noaero,iblk, &
                                    dfswthru_noaero(:,:,iblk))
         if (f_dalvdr_noaero(1:1) /= 'x') &
            call accum_hist_field(n_dalvdr_noaero,iblk, &
                                    dalvdr_noaero(:,:,iblk))
         if (f_dalidr_noaero(1:1) /= 'x') &
            call accum_hist_field(n_dalidr_noaero,iblk, &
                                    dalidr_noaero(:,:,iblk))
         if (f_dalvdf_noaero(1:1) /= 'x') &
            call accum_hist_field(n_dalvdf_noaero,iblk, &
                                    dalvdf_noaero(:,:,iblk))
         if (f_dalidf_noaero(1:1) /= 'x') &
            call accum_hist_field(n_dalidf_noaero,iblk, &
                                    dalidf_noaero(:,:,iblk))
         if (f_dalbice_noaero(1:1) /= 'x') &
            call accum_hist_field(n_dalbice_noaero,iblk, &
                                    dalbice_noaero(:,:,iblk))
         if (f_dalbsno_noaero(1:1) /= 'x') &
            call accum_hist_field(n_dalbsno_noaero,iblk, &
                                    dalbsno_noaero(:,:,iblk))
         if (f_dalbpnd_noaero(1:1) /= 'x') &
            call accum_hist_field(n_dalbpnd_noaero,iblk, &
                                    dalbpnd_noaero(:,:,iblk))
#endif
#ifdef CCSM3FRC
         if (f_dfswabs_ccsm3(1:1) /= 'x') &
            call accum_hist_field(n_dfswabs_ccsm3,iblk, &
                                    dfswabs_ccsm3(:,:,iblk))
         if (f_dfswsfc_ccsm3(1:1) /= 'x') &
            call accum_hist_field(n_dfswsfc_ccsm3,iblk, &
                                    dfswsfc_ccsm3(:,:,iblk))
         if (f_dfswint_ccsm3(1:1) /= 'x') &
            call accum_hist_field(n_dfswint_ccsm3,iblk, &
                                    dfswint_ccsm3(:,:,iblk))
         if (f_dfswthru_ccsm3(1:1) /= 'x') &
            call accum_hist_field(n_dfswthru_ccsm3,iblk, &
                                    dfswthru_ccsm3(:,:,iblk))
         if (f_dalvdr_ccsm3(1:1) /= 'x') &
            call accum_hist_field(n_dalvdr_ccsm3,iblk, &
                                    dalvdr_ccsm3(:,:,iblk))
         if (f_dalidr_ccsm3(1:1) /= 'x') &
            call accum_hist_field(n_dalidr_ccsm3,iblk, &
                                    dalidr_ccsm3(:,:,iblk))
         if (f_dalvdf_ccsm3(1:1) /= 'x') &
            call accum_hist_field(n_dalvdf_ccsm3,iblk, &
                                    dalvdf_ccsm3(:,:,iblk))
         if (f_dalidf_ccsm3(1:1) /= 'x') &
            call accum_hist_field(n_dalidf_ccsm3,iblk, &
                                    dalidf_ccsm3(:,:,iblk))
         if (f_dalbice_ccsm3(1:1) /= 'x') &
            call accum_hist_field(n_dalbice_ccsm3,iblk, &
                                    dalbice_ccsm3(:,:,iblk))
         if (f_dalbsno_ccsm3(1:1) /= 'x') &
            call accum_hist_field(n_dalbsno_ccsm3,iblk, &
                                    dalbsno_ccsm3(:,:,iblk))
#endif
#ifdef PONDFRC
         if (f_dfswabs_nopond(1:1) /= 'x') &
            call accum_hist_field(n_dfswabs_nopond,iblk, &
                                    dfswabs_nopond(:,:,iblk))
         if (f_dfswsfc_nopond(1:1) /= 'x') &
            call accum_hist_field(n_dfswsfc_nopond,iblk, &
                                    dfswsfc_nopond(:,:,iblk))
         if (f_dfswint_nopond(1:1) /= 'x') &
            call accum_hist_field(n_dfswint_nopond,iblk, &
                                    dfswint_nopond(:,:,iblk))
         if (f_dfswthru_nopond(1:1) /= 'x') &
            call accum_hist_field(n_dfswthru_nopond,iblk, &
                                    dfswthru_nopond(:,:,iblk))
         if (f_dalvdr_nopond(1:1) /= 'x') &
            call accum_hist_field(n_dalvdr_nopond,iblk, &
                                    dalvdr_nopond(:,:,iblk))
         if (f_dalidr_nopond(1:1) /= 'x') &
            call accum_hist_field(n_dalidr_nopond,iblk, &
                                    dalidr_nopond(:,:,iblk))
         if (f_dalvdf_nopond(1:1) /= 'x') &
            call accum_hist_field(n_dalvdf_nopond,iblk, &
                                    dalvdf_nopond(:,:,iblk))
         if (f_dalidf_nopond(1:1) /= 'x') &
            call accum_hist_field(n_dalidf_nopond,iblk, &
                                    dalidf_nopond(:,:,iblk))
         if (f_dalbice_nopond(1:1) /= 'x') &
            call accum_hist_field(n_dalbice_nopond,iblk, &
                                    dalbice_nopond(:,:,iblk))
         if (f_dalbsno_nopond(1:1) /= 'x') &
            call accum_hist_field(n_dalbsno_nopond,iblk, &
                                    dalbsno_nopond(:,:,iblk))
         if (f_dalbpnd_nopond(1:1) /= 'x') &
            call accum_hist_field(n_dalbpnd_nopond,iblk, &
                                    dalbpnd_nopond(:,:,iblk))
#endif
         if (f_alvdr(1:1) /= 'x') &
            call accum_hist_field(n_alvdr,  iblk, alvdr(:,:,iblk)*workb(:,:))
         if (f_alidr(1:1) /= 'x') &
            call accum_hist_field(n_alidr,  iblk, alidr(:,:,iblk)*workb(:,:))
         if (f_alvdf(1:1) /= 'x') &
            call accum_hist_field(n_alvdf,  iblk, alvdf(:,:,iblk)*workb(:,:))
         if (f_alidf(1:1) /= 'x') &
            call accum_hist_field(n_alidf,  iblk, alidf(:,:,iblk)*workb(:,:))

         if (f_albice (1:1) /= 'x') &
             call accum_hist_field(n_albice, iblk, albice(:,:,iblk))
         if (f_albsno (1:1) /= 'x') &
             call accum_hist_field(n_albsno, iblk, albsno(:,:,iblk))
         if (f_albpnd (1:1) /= 'x') &
             call accum_hist_field(n_albpnd, iblk, albpnd(:,:,iblk))
         if (f_coszen (1:1) /= 'x') &
             call accum_hist_field(n_coszen, iblk, coszen(:,:,iblk))

         if (f_flat(1:1) /= 'x') &
            call accum_hist_field(n_flat,   iblk, flat(:,:,iblk))
         if (f_flat_ai(1:1) /= 'x') &
            call accum_hist_field(n_flat_ai,iblk, flat(:,:,iblk)*workb(:,:))
         if (f_fsens(1:1) /= 'x') &
            call accum_hist_field(n_fsens,   iblk, fsens(:,:,iblk))
         if (f_fsens_ai(1:1) /= 'x') &
            call accum_hist_field(n_fsens_ai,iblk, fsens(:,:,iblk)*workb(:,:))
         if (f_flwup(1:1) /= 'x') &
            call accum_hist_field(n_flwup,   iblk, flwout(:,:,iblk))
         if (f_flwup_ai(1:1) /= 'x') &
            call accum_hist_field(n_flwup_ai,iblk, flwout(:,:,iblk)*workb(:,:))
         if (f_evap(1:1) /= 'x') &
            call accum_hist_field(n_evap,   iblk, evap(:,:,iblk))
         if (f_evap_ai(1:1) /= 'x') &
            call accum_hist_field(n_evap_ai,iblk, evap(:,:,iblk)*workb(:,:))

         if (f_Tair(1:1) /= 'x') &
            call accum_hist_field(n_Tair,   iblk, Tair(:,:,iblk))
         if (f_Tref(1:1) /= 'x') &
            call accum_hist_field(n_Tref,   iblk, Tref(:,:,iblk))
         if (f_Qref(1:1) /= 'x') &
            call accum_hist_field(n_Qref,   iblk, Qref(:,:,iblk))
         if (f_congel(1:1) /= 'x') &
            call accum_hist_field(n_congel, iblk, congel(:,:,iblk))
         if (f_frazil(1:1) /= 'x') &
            call accum_hist_field(n_frazil, iblk, frazil(:,:,iblk))
         if (f_snoice(1:1) /= 'x') &
            call accum_hist_field(n_snoice, iblk, snoice(:,:,iblk))
         if (f_meltt(1:1) /= 'x') &
            call accum_hist_field(n_meltt, iblk, meltt(:,:,iblk))
         if (f_meltb(1:1) /= 'x') &
            call accum_hist_field(n_meltb, iblk, meltb(:,:,iblk))
         if (f_meltl(1:1) /= 'x') &
            call accum_hist_field(n_meltl, iblk, meltl(:,:,iblk))
         if (f_melts(1:1) /= 'x') &
            call accum_hist_field(n_melts, iblk, melts(:,:,iblk))

         if (f_fresh(1:1) /= 'x') &
            call accum_hist_field(n_fresh, iblk, fresh(:,:,iblk))
         if (f_fresh_ai(1:1) /= 'x') &
            call accum_hist_field(n_fresh_ai,iblk, fresh_gbm(:,:,iblk))
         if (f_fsalt(1:1) /= 'x') &
            call accum_hist_field(n_fsalt, iblk, fsalt(:,:,iblk))
         if (f_fsalt_ai(1:1) /= 'x') &
            call accum_hist_field(n_fsalt_ai,iblk, fsalt_gbm(:,:,iblk))
         if (f_fhocn(1:1) /= 'x') &
            call accum_hist_field(n_fhocn, iblk, fhocn(:,:,iblk))
         if (f_fhocn_ai(1:1) /= 'x') &
            call accum_hist_field(n_fhocn_ai,iblk, fhocn_gbm(:,:,iblk))
         if (f_fswthru(1:1) /= 'x') &
            call accum_hist_field(n_fswthru, iblk, fswthru(:,:,iblk))
         if (f_fswthru_ai(1:1) /= 'x') &
            call accum_hist_field(n_fswthru_ai,iblk, fswthru_gbm(:,:,iblk))
               
         if (f_strairx(1:1) /= 'x') &
            call accum_hist_field(n_strairx, iblk, strairx(:,:,iblk))
         if (f_strairy(1:1) /= 'x') &
            call accum_hist_field(n_strairy, iblk, strairy(:,:,iblk))
         if (f_strtltx(1:1) /= 'x') &
            call accum_hist_field(n_strtltx, iblk, strtltx(:,:,iblk))
         if (f_strtlty(1:1) /= 'x') &
            call accum_hist_field(n_strtlty, iblk, strtlty(:,:,iblk))
         if (f_strcorx(1:1) /= 'x') &
            call accum_hist_field(n_strcorx, iblk, fm(:,:,iblk)*vvel(:,:,iblk))
         if (f_strcory(1:1) /= 'x') &
            call accum_hist_field(n_strcory, iblk,-fm(:,:,iblk)*uvel(:,:,iblk))
         if (f_strocnx(1:1) /= 'x') &
            call accum_hist_field(n_strocnx, iblk, strocnx(:,:,iblk))
         if (f_strocny(1:1) /= 'x') &
            call accum_hist_field(n_strocny, iblk, strocny(:,:,iblk))
         if (f_strintx(1:1) /= 'x') &
            call accum_hist_field(n_strintx, iblk, strintx(:,:,iblk))
         if (f_strinty(1:1) /= 'x') &
            call accum_hist_field(n_strinty, iblk, strinty(:,:,iblk))
         if (f_strength(1:1) /= 'x') &
            call accum_hist_field(n_strength, iblk, strength(:,:,iblk))

! The following fields (divu, shear, sig1, and sig2) will be smeared
!  if averaged over more than a few days.
! Snapshots may be more useful (see below).

! Need divu and shear as monthly means for IPCC/CMIP.
         if (f_divu(1:1) /= 'x') &
            call accum_hist_field(n_divu, iblk, divu(:,:,iblk))
         if (f_shear(1:1) /= 'x') &
            call accum_hist_field(n_shear, iblk, shear(:,:,iblk))
!        if (f_sig1(1:1) /= 'x') &
!           call accum_hist_field(n_sig1, iblk, sig1(:,:,iblk))
!        if (f_sig2(1:1) /= 'x') &
!           call accum_hist_field(n_sig2, iblk, sig2(:,:,iblk))

         if (f_dvidtt(1:1) /= 'x') &
            call accum_hist_field(n_dvidtt, iblk, dvidtt(:,:,iblk))
         if (f_dvidtd(1:1) /= 'x') &
            call accum_hist_field(n_dvidtd, iblk, dvidtd(:,:,iblk))
         if (f_daidtt(1:1) /= 'x') &
            call accum_hist_field(n_daidtt, iblk, daidtt(:,:,iblk))
         if (f_daidtd(1:1) /= 'x') &
            call accum_hist_field(n_daidtd, iblk, daidtd(:,:,iblk))

         if (f_opening(1:1) /= 'x') &
            call accum_hist_field(n_opening, iblk, opening(:,:,iblk))
         if (f_dardg1dt(1:1) /= 'x') &
            call accum_hist_field(n_dardg1dt, iblk, dardg1dt(:,:,iblk))
         if (f_dardg2dt(1:1) /= 'x') &
            call accum_hist_field(n_dardg2dt, iblk, dardg2dt(:,:,iblk))
         if (f_dvirdgdt(1:1) /= 'x') &
            call accum_hist_field(n_dvirdgdt, iblk, dvirdgdt(:,:,iblk))

         if (f_alvl(1:1)/= 'x') &
             call accum_hist_field(n_alvl,   iblk, &
                                   aice(:,:,iblk) * trcr(:,:,nt_alvl,iblk))
         if (f_vlvl(1:1)/= 'x') &
             call accum_hist_field(n_vlvl,   iblk, &
                                   vice(:,:,iblk) * trcr(:,:,nt_vlvl,iblk))
         if (f_ardg(1:1)/= 'x') &
             call accum_hist_field(n_ardg,   iblk, &
                             aice(:,:,iblk) * (c1 - trcr(:,:,nt_alvl,iblk)))
         if (f_vrdg(1:1)/= 'x') &
             call accum_hist_field(n_vrdg,   iblk, &
                             vice(:,:,iblk) * (c1 - trcr(:,:,nt_vlvl,iblk)))

         if (f_fsurf_ai(1:1) /= 'x') &
            call accum_hist_field(n_fsurf_ai, iblk, fsurf(:,:,iblk)*workb(:,:))
         if (f_fcondtop_ai(1:1) /= 'x') &
            call accum_hist_field(n_fcondtop_ai, iblk, &
                                                 fcondtop(:,:,iblk)*workb(:,:))

         if (f_icepresent(1:1) /= 'x') then
           worka(:,:) = c0
           do j = jlo, jhi
           do i = ilo, ihi
              if (aice(i,j,iblk) > puny) worka(i,j) = c1
           enddo
           enddo
           call accum_hist_field(n_icepresent, iblk, worka(:,:))
         endif

         nct = min(ncat, ncat_hist)
         do n=1,nct
            workb(:,:) = aicen_init(:,:,n,iblk)
            if (f_aicen(1:1) /= 'x') &
               call accum_hist_field(n_aicen(n,:), iblk, aicen(:,:,n,iblk))
            if (f_vicen(1:1) /= 'x') &
               call accum_hist_field(n_vicen(n,:), iblk, vicen(:,:,n,iblk))
            if (f_apondn(1:1) /= 'x') &
               call accum_hist_field(n_apondn(n,:), iblk, apondn(:,:,n,iblk))
            if (f_fsurfn_ai(1:1) /= 'x') &
               call accum_hist_field(n_fsurfn_ai(n,:), iblk, &
               fsurfn(:,:,n,iblk)*workb(:,:))
            if (f_fcondtopn_ai(1:1) /= 'x') &
               call accum_hist_field(n_fcondtopn_ai(n,:), iblk, &
               fcondtopn(:,:,n,iblk)*workb(:,:))
            if (f_flatn_ai(1:1) /= 'x') &
               call accum_hist_field(n_flatn_ai(n,:), iblk, &
               flatn(:,:,n,iblk)*workb(:,:))
            ! Calculate surface heat flux that causes melt as this
            ! is what is calculated by the atmos in HadGEM3 so
            ! needed for checking purposes
            if (f_fmelttn_ai(1:1) /= 'x') &
               call accum_hist_field(n_fmelttn_ai(n,:), iblk, &
               max(fsurfn(:,:,n,iblk) - fcondtopn(:,:,n,iblk),c0)*workb(:,:))
        enddo                    ! n
        if (f_fs(1:1) /= 'x') then
           worka(:,:) = c0
           nct = min(ncat, ncat_hist)
           do n=1,nct
              workb(:,:) = aicen_init(:,:,n,iblk)
              do j = jlo, jhi
              do i = ilo, ihi
                 hs = c0
                 if (workb(i,j) > puny) &
                    hs = vsnon(i,j,n,iblk)/workb(i,j)
                 if (hs >= hsmin) &
                    worka(i,j) = worka(i,j) + workb(i,j)*min( hs/hs0, c1 )
              enddo
              enddo
           enddo
           call accum_hist_field(n_fs, iblk, worka(:,:))
        endif
        ! Compute the internal heat content (enthalpy) of the ice 
        if (f_qi(1:1) /= 'x') then
           worka(:,:) = c0
           nct = min(ncat, ncat_hist)
           do n=1,nct
              workb(:,:) = aicen_init(:,:,n,iblk)
              do k=1,nilyr
                 do j = jlo, jhi
                 do i = ilo, ihi
                    qi = c0
                    if (vicen(i,j,n,iblk) > puny) &
                       qi = eicen(i,j,ilyr1(n)+k-1,iblk)*tarea(i,j,iblk)
                    worka(i,j) = worka(i,j) + workb(i,j)*qi
                 enddo
                 enddo
              enddo
           enddo
           call accum_hist_field(n_qi, iblk, worka(:,:))
         endif
        ! Compute the internal heat content (enthalpy) of the snow
         if (f_qs(1:1) /= 'x') then
            worka(:,:) = c0
            nct = min(ncat, ncat_hist)
            do n=1,nct
              workb(:,:) = aicen_init(:,:,n,iblk)
               do k=1,nslyr
                  do j = jlo, jhi
                  do i = ilo, ihi
                     qs = c0
                     if (vsnon(i,j,n,iblk) > puny) &
                        qs = esnon(i,j,slyr1(n)+k-1,iblk)*tarea(i,j,iblk)
                     worka(i,j) = worka(i,j) + workb(i,j)*qs
                  enddo
                  enddo
               enddo
            enddo
            call accum_hist_field(n_qs, iblk, worka(:,:))
         endif

        ! Calculate aggregate melt pond area by summing category values
        if (f_apond(1:1) /= 'x') then
           do ns=1, nstreams
              if (n_apond(ns) /= 0) then
                 worka(:,:) = c0
                 do j = jlo, jhi
                 do i = ilo, ihi
                  if (tmask(i,j,iblk)) then
                    do n=1,nct
                       worka(i,j)  = worka(i,j) + aa(i,j,n_apondn(n,ns),iblk)
                    enddo            ! n
                  endif              ! tmask
                 enddo                ! i
                 enddo                ! j
                 aa(:,:,n_apond(ns),iblk) = worka(:,:)
              endif
           enddo
        endif
        ! Calculate aggregate surface melt flux by summing category values
        if (f_fmeltt_ai(1:1) /= 'x') then
           do ns=1, nstreams
              if (n_fmeltt_ai(ns) /= 0) then
                 worka(:,:) = c0
                 do j = jlo, jhi
                 do i = ilo, ihi
                  if (tmask(i,j,iblk)) then
                    do n=1,nct
                       worka(i,j) = worka(i,j) + aa(i,j,n_fmelttn_ai(n,ns),iblk)
                    enddo            ! n
                  endif              ! tmask
                 enddo                ! i
                 enddo                ! j
                 aa(:,:,n_fmeltt_ai(ns),iblk) = worka(:,:)
              endif
           enddo
        endif

        ! Aerosols

        if (f_faero_atm(1:1) /= 'x') then
           do n=1,n_aero
              call accum_hist_field(n_faero_atm(n,:), iblk, &
                                    faero(:,:,n,iblk))
           enddo
        endif

        if (f_faero_ocn(1:1) /= 'x') then
           do n=1,n_aero
              call accum_hist_field(n_faero_ocn(n,:), iblk, &
                                    fsoot(:,:,n,iblk))
           enddo
        endif

        if (f_aero(1:1) /= 'x') then
           do n=1,n_aero
              call accum_hist_field(n_aerosn1(n,:), iblk, &
                                 trcr(:,:,nt_aero  +4*(n-1),iblk)/rhos)
              call accum_hist_field(n_aerosn2(n,:), iblk, &
                                 trcr(:,:,nt_aero+1+4*(n-1),iblk)/rhos)
              call accum_hist_field(n_aeroic1(n,:), iblk, &
                                 trcr(:,:,nt_aero+2+4*(n-1),iblk)/rhoi)
              call accum_hist_field(n_aeroic2(n,:), iblk, &
                                 trcr(:,:,nt_aero+3+4*(n-1),iblk)/rhoi)
           enddo
        endif

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      !---------------------------------------------------------------
      ! Write output files at prescribed intervals
      !---------------------------------------------------------------

      nstrm = nstreams
      if (write_ic) nstrm = 1

      do ns = 1, nstrm

      if (write_history(ns) .or. write_ic) then

      !---------------------------------------------------------------
      ! Mask out land points and convert units 
      !---------------------------------------------------------------

        ravgct = c1/avgct(ns)
	!$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j,k)
        do iblk = 1, nblocks
           this_block = get_block(blocks_ice(iblk),iblk)         
           ilo = this_block%ilo
           ihi = this_block%ihi
           jlo = this_block%jlo
           jhi = this_block%jhi

           do k = 1, num_avail_hist_fields

              if (avail_hist_fields(k)%vhistfreq == histfreq(ns)) then

                 do j = jlo, jhi
                 do i = ilo, ihi
                    if (.not. tmask(i,j,iblk)) then ! mask out land points
                       aa(i,j,k,iblk) = spval
                    else                            ! convert units
                       aa(i,j,k,iblk) &
                          = avail_hist_fields(k)%cona*aa(i,j,k,iblk) &
                          * ravgct + avail_hist_fields(k)%conb
                    endif
                 enddo             ! i
                 enddo             ! j

                 ! back out albedo/zenith angle dependence
                 if (avail_hist_fields(k)%vname == 'albice') then
                    do j = jlo, jhi
                    do i = ilo, ihi
                       if (tmask(i,j,iblk)) then 
                          ravgctz = c0
                          if (albcnt(i,j,iblk,ns) > puny) &
                              ravgctz = c1/albcnt(i,j,iblk,ns)
                          if (n_albice(ns) /= 0) aa(i,j,n_albice(ns),iblk) = &
                             aa(i,j,n_albice(ns),iblk)*avgct(ns)*ravgctz
                          if (n_albsno(ns) /= 0) aa(i,j,n_albsno(ns),iblk) = &
                             aa(i,j,n_albsno(ns),iblk)*avgct(ns)*ravgctz
                          if (n_albpnd(ns) /= 0) aa(i,j,n_albpnd(ns),iblk) = &
                             aa(i,j,n_albpnd(ns),iblk)*avgct(ns)*ravgctz
#ifdef AEROFRC
                          if (n_dalbice_noaero(ns) /= 0) &
                             aa(i,j,n_dalbice_noaero(ns),iblk) = &
                             aa(i,j,n_dalbice_noaero(ns),iblk)*avgct(ns)*ravgctz
                          if (n_dalbsno_noaero(ns) /= 0) &
                             aa(i,j,n_dalbsno_noaero(ns),iblk) = &
                             aa(i,j,n_dalbsno_noaero(ns),iblk)*avgct(ns)*ravgctz
                          if (n_dalbpnd_noaero(ns) /= 0) &
                             aa(i,j,n_dalbpnd_noaero(ns),iblk) = &
                             aa(i,j,n_dalbpnd_noaero(ns),iblk)*avgct(ns)*ravgctz
#endif
#ifdef CCSM3FRC
                          if (n_dalbice_ccsm3(ns) /= 0) &
                             aa(i,j,n_dalbice_ccsm3(ns),iblk) = &
                             aa(i,j,n_dalbice_ccsm3(ns),iblk)*avgct(ns)*ravgctz
                          if (n_dalbsno_ccsm3(ns) /= 0) &
                             aa(i,j,n_dalbsno_ccsm3(ns),iblk) = &
                             aa(i,j,n_dalbsno_ccsm3(ns),iblk)*avgct(ns)*ravgctz
#endif
#ifdef PONDFRC
                          if (n_dalbice_nopond(ns) /= 0) &
                             aa(i,j,n_dalbice_nopond(ns),iblk) = &
                             aa(i,j,n_dalbice_nopond(ns),iblk)*avgct(ns)*ravgctz
                          if (n_dalbsno_nopond(ns) /= 0) &
                             aa(i,j,n_dalbsno_nopond(ns),iblk) = &
                             aa(i,j,n_dalbsno_nopond(ns),iblk)*avgct(ns)*ravgctz
                          if (n_dalbpnd_nopond(ns) /= 0) &
                             aa(i,j,n_dalbpnd_nopond(ns),iblk) = &
                             aa(i,j,n_dalbpnd_nopond(ns),iblk)*avgct(ns)*ravgctz
#endif
                       endif
                    enddo             ! i
                    enddo             ! j
                 endif

              endif
           enddo                ! k

      !---------------------------------------------------------------
      ! snapshots
      !---------------------------------------------------------------

          ! compute sig1 and sig2
        
           call principal_stress (nx_block,  ny_block,  &
                                  stressp_1 (:,:,iblk), &
                                  stressm_1 (:,:,iblk), &
                                  stress12_1(:,:,iblk), &
                                  prs_sig   (:,:,iblk), &
                                  sig1      (:,:,iblk), &
                                  sig2      (:,:,iblk))
 
           do j = jlo, jhi
           do i = ilo, ihi
              if (.not. tmask(i,j,iblk)) then ! mask out land points
!                if (n_divu     (ns) /= 0) aa(i,j,n_divu(ns),iblk)      = spval
!                if (n_shear    (ns) /= 0) aa(i,j,n_shear(ns),iblk)     = spval
                 if (n_sig1     (ns) /= 0) aa(i,j,n_sig1(ns),iblk )     = spval
                 if (n_sig2     (ns) /= 0) aa(i,j,n_sig2(ns),iblk )     = spval
                 if (n_mlt_onset(ns) /= 0) aa(i,j,n_mlt_onset(ns),iblk) = spval
                 if (n_frz_onset(ns) /= 0) aa(i,j,n_frz_onset(ns),iblk) = spval
                 if (n_hisnap   (ns) /= 0) aa(i,j,n_hisnap(ns),iblk)    = spval
                 if (n_aisnap   (ns) /= 0) aa(i,j,n_aisnap(ns),iblk)    = spval
                 if (n_trsig    (ns) /= 0) aa(i,j,n_trsig(ns),iblk )    = spval
                 if (n_iage     (ns) /= 0) aa(i,j,n_iage(ns),iblk )     = spval
                 if (n_FY       (ns) /= 0) aa(i,j,n_FY(ns),iblk )       = spval
              else
!                if (n_divu     (ns) /= 0) aa(i,j,n_divu(ns),iblk)      = &
!                      divu (i,j,iblk)*avail_hist_fields(n_divu(ns))%cona
!                if (n_shear    (ns) /= 0) aa(i,j,n_shear(ns),iblk)     = &
!                      shear(i,j,iblk)*avail_hist_fields(n_shear(ns))%cona
                 if (n_sig1     (ns) /= 0) aa(i,j,n_sig1(ns),iblk)      = &
                       sig1 (i,j,iblk)*avail_hist_fields(n_sig1(ns))%cona
                 if (n_sig2     (ns) /= 0) aa(i,j,n_sig2(ns),iblk)      = &
                       sig2 (i,j,iblk)*avail_hist_fields(n_sig2(ns))%cona
                 if (n_mlt_onset(ns) /= 0) aa(i,j,n_mlt_onset(ns),iblk) = &
                       mlt_onset(i,j,iblk)
                 if (n_frz_onset(ns) /= 0) aa(i,j,n_frz_onset(ns),iblk) = &
                       frz_onset(i,j,iblk)
                 if (n_hisnap   (ns) /= 0) aa(i,j,n_hisnap(ns),iblk)    = &
                       vice(i,j,iblk)
                 if (n_aisnap   (ns) /= 0) aa(i,j,n_aisnap(ns),iblk)    = &
                       aice(i,j,iblk)
                 if (n_trsig    (ns) /= 0) aa(i,j,n_trsig(ns),iblk )    = &
                                       p25*(stressp_1(i,j,iblk) &
                                          + stressp_2(i,j,iblk) &
                                          + stressp_3(i,j,iblk) &
                                          + stressp_4(i,j,iblk))
                 if (n_iage     (ns) /= 0) aa(i,j,n_iage(ns),iblk)      = &
                       trcr(i,j,nt_iage,iblk)*avail_hist_fields(n_iage(ns))%cona
                 if (n_FY       (ns) /= 0) aa(i,j,n_FY(ns),iblk)        = &
                       trcr(i,j,nt_FY,iblk)*avail_hist_fields(n_FY(ns))%cona
            endif
           enddo                ! i
           enddo                ! j

        enddo                   ! iblk
        !$OMP END PARALLEL DO

        time_end(ns) = time/int(secday)

      !---------------------------------------------------------------
      ! write file
      !---------------------------------------------------------------

        call ice_timer_start(timer_readwrite)  ! reading/writing

        if (trim(history_format) == 'nc') then
           call icecdf(ns)         ! netcdf output
        else
           call icebin(ns)         ! binary output
        endif

        call ice_timer_stop(timer_readwrite)  ! reading/writing

      !---------------------------------------------------------------
      ! reset to zero
      !------------------------------------------------------------
        if (write_ic) then
           aa(:,:,:,:) = c0
           avgct(:) = c0
           albcnt(:,:,:,:) = c0
           write_ic = .false.        ! write initial condition once at most
        else
           avgct(ns) = c0
           albcnt(:,:,:,ns) = c0
        endif

        do k=1,num_avail_hist_fields
           if (avail_hist_fields(k)%vhistfreq == histfreq(ns)) aa(:,:,k,:) = c0
        enddo

      endif  ! write_history or write_ic
      enddo  ! nstreams

      !$OMP PARALLEL DO PRIVATE(iblk,this_block,ilo,ihi,jlo,jhi,i,j,k)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         if (new_year) then

            do j=jlo,jhi
            do i=ilo,ihi
               ! reset NH Jan 1
               if (lmask_n(i,j,iblk)) mlt_onset(i,j,iblk) = c0
               ! reset SH Jan 1 
               if (lmask_s(i,j,iblk)) frz_onset(i,j,iblk) = c0
            enddo
            enddo
         endif                  ! new_year

         if ((yday >= 181._dbl_kind) .and. &
             (yday <  181._dbl_kind+dt/secday)) then

            do j=jlo,jhi
            do i=ilo,ihi

               ! reset SH Jul 1
               if (lmask_s(i,j,iblk)) mlt_onset(i,j,iblk) = c0

               ! reset NH Jul 1
               if (lmask_n(i,j,iblk)) frz_onset(i,j,iblk) = c0
            enddo
            enddo

         endif                  ! yday
      enddo                     ! iblk

      write_ic = .false.        ! write initial condition once at most

      end subroutine ice_write_hist

!=======================================================================

      subroutine define_hist_field(id, vname, vunit, vcoord, vcellmeas, &
                                   vdesc, vcomment, cona, conb, &
                                   ns1, vhistfreq)

!     !DESCRIPTION:
!     Initializes description of an available field and returns location
!     in the available fields array for use in later calls.
!
!     !REVISION HISTORY:
!     !REVISION HISTORY:
!     2009 Created by D. Bailey following POP

!     !USES:
      use ice_exit
      use ice_calendar, only: histfreq, histfreq_n, nstreams

!     !OUTPUT PARAMETERS:

      integer (int_kind), dimension(max_nstrm), intent(out) :: &
         id                ! location in avail_fields array for use in
                           ! later routines

!     !INPUT PARAMETERS

      character (len=*), intent(in) :: &
         vname      , & ! variable names
         vunit      , & ! variable units
         vcoord     , & ! variable coordinates
         vcellmeas  , & ! variables cell measures
         vdesc      , & ! variable descriptions
         vcomment       ! variable comments

      real (kind=dbl_kind), intent(in) :: &
         cona       , & ! multiplicative conversion factor
         conb           ! additive conversion factor

      character (len=*), intent(in) :: &
         vhistfreq      ! history frequency
 
      integer (kind=int_kind), intent(in) :: &
         ns1            ! stream index

      integer (kind=int_kind) :: &
         ns         , & ! loop index
         lenf           ! length of namelist string

      character (len=40) :: stmp

      lenf = len(trim(vhistfreq))
      if (ns1 == 1) id(:) = 0

      do ns = 1, nstreams
         if (vhistfreq(ns1:ns1) == histfreq(ns)) then

            num_avail_hist_fields = num_avail_hist_fields + 1

            if (num_avail_hist_fields > max_avail_hist_fields) &
               call abort_ice("Need to increase max_avail_hist_fields")

            id(ns) = num_avail_hist_fields

            stmp = vname
            if (lenf > 1 .and. ns1 > 1) &
               write(stmp,'(a,a1,a1)') trim(stmp),'_',vhistfreq(ns1:ns1)

            avail_hist_fields(id(ns))%vname = trim(stmp)
            avail_hist_fields(id(ns))%vunit = trim(vunit)
            avail_hist_fields(id(ns))%vcoord = trim(vcoord)
            avail_hist_fields(id(ns))%vcellmeas = trim(vcellmeas)
            avail_hist_fields(id(ns))%vdesc = trim(vdesc)
            avail_hist_fields(id(ns))%vcomment = trim(vcomment)
            avail_hist_fields(id(ns))%cona = cona
            avail_hist_fields(id(ns))%conb = conb
            avail_hist_fields(id(ns))%vhistfreq = vhistfreq(ns1:ns1)
            avail_hist_fields(id(ns))%vhistfreq_n = histfreq_n(ns)

         endif
      enddo

      end subroutine define_hist_field

!=======================================================================

      subroutine accum_hist_field(id, iblk, field_accum)

!     !DESCRIPTION:
!     Accumulates a history field
!
!     !REVISION HISTORY:
!     2009 Created by D. Bailey following POP

      use ice_domain
      use ice_grid, only: tmask
      use ice_calendar, only: nstreams

!     !OUTPUT PARAMETERS:

      integer (int_kind), dimension(max_nstrm), intent(in) :: &
         id                ! location in avail_fields array for use in
                           ! later routines

      integer (kind=int_kind), intent(in) :: iblk

      real (kind=dbl_kind), intent(in) :: field_accum(nx_block,ny_block)

      type (block) :: &
         this_block           ! block information for current block

      integer (kind=int_kind) :: i,j, ilo, ihi, jlo, jhi, ns, idns

      !---------------------------------------------------------------
      ! increment field
      !---------------------------------------------------------------

      do ns = 1, nstreams
         idns = id(ns)
         if (idns > 0) then

            this_block = get_block(blocks_ice(iblk),iblk)
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi

            do j = jlo, jhi
            do i = ilo, ihi
               if (tmask(i,j,iblk)) then
                  aa(i,j,idns, iblk) = aa(i,j,idns, iblk) + field_accum(i,j)
               endif
            enddo
            enddo

         endif
      enddo

      end subroutine accum_hist_field

!=======================================================================

      end module ice_history

!=======================================================================
