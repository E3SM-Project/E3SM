!=======================================================================
!
!BOP
!
! !MODULE: ice_history_fields - data types for ice_history

      module ice_history_fields

! !REVISION HISTORY:
! authors Mariana Vertenstein, NCAR
!
! !USES:
!
      use ice_kinds_mod
      use ice_domain_size	
!
!EOP
!
      implicit none
      save

      logical (kind=log_kind) :: &
         hist_avg  ! if true, write averaged data instead of snapshots

      character (len=char_len_long) :: &
         history_dir   , & ! directory name for history file
         incond_dir        ! directory for snapshot initial conditions

      character (len=char_len) :: &
         history_file  , & ! output file for history
         incond_file       ! output file for snapshot initial conditions

      character (len=3) :: &
         history_format    ! file format ('bin'=binary or 'nc'=netcdf)

      !---------------------------------------------------------------
      ! Instructions for adding a field: (search for 'example')
      !     Here:
      ! (1) Add to frequency flags (f_<field>)
      ! (2) Add to namelist (here and also in ice_in)
      ! (3) Add to index list
      !     In init_hist:
      ! (4) Add define_hist_field call with vname, vdesc, vunit,
      !     and vcomment, vcellmeas, and conversion factor if necessary.
      ! (5) Add flag to broadcast list
      ! (6) Add accum_hist_field call with appropriate variable
      !---------------------------------------------------------------

      type, public :: ice_hist_field
          character (len=16) :: vname     ! variable name
          character (len=16) :: vunit     ! variable units
          character (len=16) :: vcoord    ! variable coordinates
          character (len=16) :: vcellmeas ! variable cell measures
          character (len=55) :: vdesc     ! variable description
          character (len=55) :: vcomment  ! variable description
          real (kind=dbl_kind) :: cona    ! multiplicative conversion factor
          real (kind=dbl_kind) :: conb    ! additive conversion factor
          character (len=1) :: vhistfreq  ! frequency of history output
          integer (kind=int_kind) :: vhistfreq_n ! number of vhistfreq intervals
      end type

      integer (kind=int_kind), parameter :: &
         max_avail_hist_fields = 600      ! Max number of history fields

      integer (kind=int_kind) :: &
         num_avail_hist_fields   = 0      ! Current number of defined fields

      type (ice_hist_field), dimension(max_avail_hist_fields) :: &
         avail_hist_fields

      !---------------------------------------------------------------
      ! primary info for the history file
      !---------------------------------------------------------------

      integer (kind=int_kind), parameter :: &
         ncat_hist = ncat         ! number of ice categories written <= ncat

      integer (kind=int_kind), parameter :: &
         nvar = 12                  ! number of grid fields that can be written
                                    !   excluding grid vertices

      real (kind=real_kind) :: time_beg(max_nstrm), &
                               time_end(max_nstrm) ! bounds for averaging
      real (kind=real_kind) :: time_bounds(2)

      real (kind=dbl_kind), allocatable :: &
         aa(:,:,:,:)       ! field accumulations and averages
         
      logical (kind=log_kind) :: &
         igrd(nvar)        ! true if grid field is written to output file

      !---------------------------------------------------------------
      ! logical flags: write to output file if true
      !---------------------------------------------------------------

       logical (kind=log_kind) :: &
           f_bounds    = .true.

      !---------------------------------------------------------------
      ! field indices
      !---------------------------------------------------------------

       integer (kind=int_kind), parameter :: &
           n_tmask      = 1,  &
           n_tarea      = 2,  &
           n_uarea      = 3,  &
           n_dxt        = 4,  &
           n_dyt        = 5,  &
           n_dxu        = 6,  &
           n_dyu        = 7,  &
           n_HTN        = 8,  &
           n_HTE        = 9,  &
           n_ANGLE      = 10, &
           n_ANGLET     = 11, &
           n_blkmask    = 12, &

           n_lont_bnds  = 1, &
           n_latt_bnds  = 2, &
           n_lonu_bnds  = 3, &
           n_latu_bnds  = 4

       integer (kind=int_kind), dimension(max_nstrm) :: &
!          n_example    , &
           n_hi         , n_hs         , &
           n_fs         ,                &
           n_Tsfc       , n_aice       , &
           n_uvel       , n_vvel       , &
           n_transix    , n_transiy    , &
           n_fswdn      , n_fswup      , &
           n_flwdn      ,                &
           n_snow       , n_snow_ai    , &
           n_rain       , n_rain_ai    , &
           n_sst        , n_sss        , &
           n_uocn       , n_vocn       , &
           n_frzmlt     , n_fswfac     , &
#if (defined AEROFRC) || (defined PONDFRC) || (defined CCSM3FRC)
           n_fswsfc_ai , &
           n_fswint_ai , &
#endif
#ifdef AEROFRC
           n_dfswabs_noaero   , &
           n_dfswsfc_noaero  , &
           n_dfswint_noaero  , &
           n_dfswthru_noaero  , &
           n_dalvdr_noaero,   n_dalidr_noaero      , &
           n_dalvdf_noaero,   n_dalidf_noaero      , &
           n_dalbice_noaero,  n_dalbsno_noaero     , &
           n_dalbpnd_noaero,                        &
#endif
#ifdef CCSM3FRC
           n_dfswabs_ccsm3   , &
           n_dfswsfc_ccsm3  , &
           n_dfswint_ccsm3  , &
           n_dfswthru_ccsm3  , &
           n_dalvdr_ccsm3,   n_dalidr_ccsm3      , &
           n_dalvdf_ccsm3,   n_dalidf_ccsm3      , &
           n_dalbice_ccsm3,  n_dalbsno_ccsm3     , &
#endif
#ifdef PONDFRC
           n_dfswabs_nopond   , &
           n_dfswsfc_nopond  , &
           n_dfswint_nopond  , &
           n_dfswthru_nopond  , &
           n_dalvdr_nopond,   n_dalidr_nopond      , &
           n_dalvdf_nopond,   n_dalidf_nopond      , &
           n_dalbice_nopond,  n_dalbsno_nopond     , &
           n_dalbpnd_nopond,                        &
#endif
           n_fswabs     , n_fswabs_ai  , &
           n_alvdr      , n_alidr      , &
           n_alvdf      , n_alidf      , &
           n_albice     , n_albsno     , &
           n_albpnd     , n_coszen     , &
           n_flat       , n_flat_ai    , &
           n_fsens      , n_fsens_ai   , &
           n_flwup      , n_flwup_ai   , &
           n_evap       , n_evap_ai    , &
           n_qi         , n_qs         , &
           n_Tair       ,                &
           n_Tref       , n_Qref       , &
           n_congel     , n_frazil     , &
           n_snoice     ,                &
           n_meltt      , n_melts      , &
           n_meltb      , n_meltl      , &
           n_fresh      , n_fresh_ai   , &
           n_fsalt      , n_fsalt_ai   , &
           n_fhocn      , n_fhocn_ai   , &
           n_fswthru    , n_fswthru_ai , &
           n_strairx    , n_strairy    , &
           n_strtltx    , n_strtlty    , &
           n_strcorx    , n_strcory    , &
           n_strocnx    , n_strocny    , &
           n_strintx    , n_strinty    , &
           n_strength   , n_opening    , &
           n_divu       , n_shear      , &
           n_sig1       , n_sig2       , &
           n_dvidtt     , n_dvidtd     , &
           n_daidtt     , n_daidtd     , &
           n_mlt_onset  , n_frz_onset  , &
           n_dardg1dt   , n_dardg2dt   , &
           n_dvirdgdt   ,                &
           n_hisnap     , n_aisnap     , &
           n_trsig      , n_icepresent , &
           n_iage       , n_FY         , &
           n_apond      , &
           n_ardg       , n_vrdg       , &
           n_alvl       , n_vlvl       , &
           n_fsurf_ai   , n_fcondtop_ai, &
           n_fmeltt_ai

      ! Category dependent variables
      integer(kind=int_kind), dimension(ncat_hist,max_nstrm) :: &
           n_aicen       , &
           n_vicen       , &
           n_apondn      , &
           n_fsurfn_ai   , &
           n_fcondtopn_ai, &
           n_fmelttn_ai  , &
           n_flatn_ai

      integer(kind=int_kind), dimension(n_aeromx,max_nstrm) :: &
           n_faero_atm    , &
           n_faero_ocn    , &
           n_aerosn1      , &
           n_aerosn2      , &
           n_aeroic1      , &
           n_aeroic2

      integer(kind=int_kind), dimension(n_aeromx*ncat_hist,max_nstrm) :: &
           n_aerosn1n, &
           n_aerosn2n, &
           n_aeroic1n, &
           n_aeroic2n

!=======================================================================

      contains

!=======================================================================

      subroutine construct_filename(ncfile,suffix,ns)

      use ice_calendar, only: time, sec, idate, nyr, month, daymo,  &
                              mday, write_ic, histfreq, histfreq_n, &
                              year_init, new_year, new_month, new_day, &
                              dayyr, dt
      use ice_restart, only: lenstr

      integer (kind=int_kind), intent(in) :: ns
      character (char_len_long), intent(inout) :: ncfile
      character (len=2), intent(in) :: suffix
      character (len=1) :: cstream

      integer (kind=int_kind) :: iyear, imonth, iday, isec

        iyear = nyr + year_init - 1 ! set year_init=1 in ice_in to get iyear=nyr
        imonth = month
        iday = mday
!tcx        isec = sec - dt
        isec = sec

        ! construct filename
        if (write_ic) then
           write(ncfile,'(a,i4.4,a,i2.2,a,i2.2,a,i5.5,a,a)')  &
              incond_file(1:lenstr(incond_file)),iyear,'-', &
              imonth,'-',iday,'-',isec,'.',suffix
        else

         if (hist_avg) then
          if (histfreq(ns).eq.'1') then
           ! do nothing
          elseif (histfreq(ns).eq.'h'.or.histfreq(ns).eq.'H') then
           ! do nothing
          elseif (new_year) then
           iyear = iyear - 1
           imonth = 12
           iday = daymo(imonth)
          elseif (new_month) then
           imonth = month - 1
           iday = daymo(imonth)
          elseif (new_day) then
           iday = iday - 1
          endif
         endif

         cstream = ''
         if (ns > 1) write(cstream,'(i1.1)') ns-1

         if (histfreq(ns) == '1') then ! instantaneous, write every dt

          write(ncfile,'(a,a,i4.4,a,i2.2,a,i2.2,a,i5.5,a,a)')  &
            history_file(1:lenstr(history_file))//trim(cstream),'_inst.', &
             iyear,'-',imonth,'-',iday,'-',sec,'.',suffix

         elseif (hist_avg) then    ! write averaged data

          if (histfreq(ns).eq.'d'.or.histfreq(ns).eq.'D') then     ! daily
           write(ncfile,'(a,a,i4.4,a,i2.2,a,i2.2,a,a)')  &
            history_file(1:lenstr(history_file))//trim(cstream), &
             '.',iyear,'-',imonth,'-',iday,'.',suffix
          elseif (histfreq(ns).eq.'h'.or.histfreq(ns).eq.'H') then ! hourly
           write(ncfile,'(a,a,i2.2,a,i4.4,a,i2.2,a,i2.2,a,i5.5,a,a)')  &
            history_file(1:lenstr(history_file))//trim(cstream),'_',   &
             histfreq_n(ns),'h.',iyear,'-',imonth,'-',iday,'-',sec,'.',suffix
          elseif (histfreq(ns).eq.'m'.or.histfreq(ns).eq.'M') then ! monthly
           write(ncfile,'(a,a,i4.4,a,i2.2,a,a)')  &
            history_file(1:lenstr(history_file))//trim(cstream),'.', &
             iyear,'-',imonth,'.',suffix
          elseif (histfreq(ns).eq.'y'.or.histfreq(ns).eq.'Y') then ! yearly
           write(ncfile,'(a,a,i4.4,a,a)') &
            history_file(1:lenstr(history_file))//trim(cstream),'.', &
	      iyear,'.',suffix
          endif

         else                     ! instantaneous with histfreq > dt
           write(ncfile,'(a,a,i4.4,a,i2.2,a,i2.2,a,i5.5,a,a)')  &
            history_file(1:lenstr(history_file)),'_inst.', &
             iyear,'-',imonth,'-',iday,'-',sec,'.',suffix
         endif
        endif

      end subroutine construct_filename

      end module ice_history_fields
