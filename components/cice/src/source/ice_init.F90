!=======================================================================
!BOP
!
! !MODULE:   ice_init - parameter and variable initializations
!
! !DESCRIPTION:
!
! parameter and variable initializations
!
! !REVISION HISTORY:
!  SVN:$Id: ice_init.F90 143 2008-08-08 23:08:22Z eclare $
!
! authors Elizabeth C. Hunke and William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! 2004 WHL: Block structure added 
! 2006 ECH: Added namelist variables, warnings.
!           Replaced old default initial ice conditions with 3.14 version.
!           Converted to free source form (F90).
!
! !INTERFACE:
!
      module ice_init
!
! !USES:
!
      use ice_kinds_mod
      use ice_communicate, only: my_task, master_task, lprint_stats
      use ice_domain_size
      use ice_constants
!
!EOP
!
      implicit none
      save

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: input_data - namelist variables
!
! !INTERFACE:
!
      subroutine input_data
!
! !DESCRIPTION:
!
! Namelist variables, set to default values; may be altered
! at run time
!
! !REVISION HISTORY:
!
! author Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_broadcast
      use ice_diagnostics
      use ice_fileunits
      use ice_calendar, only: year_init, istep0, histfreq, histfreq_n, &
                              dumpfreq, dumpfreq_n, diagfreq, &
                              npt, dt, xndt_dyn, days_per_year, write_ic
      use ice_restart, only: &
          restart, restart_dir, restart_file, pointer_file, &
          runid, runtype, ice_ic, resttype, restart_format, lcdf64
      use ice_history_fields, only: hist_avg, &
                             history_format, history_dir, history_file, &
                             incond_dir, incond_file
      use ice_exit
      use ice_itd, only: kitd, kcatbound
      use ice_ocean, only: oceanmixed_ice
      use ice_flux, only: Tfrzpt, update_ocn_f
      use ice_forcing, only: &
          ycycle,          fyear_init,    dbug, &
          atm_data_type,   atm_data_dir,  precip_units, &
          atm_data_format, ocn_data_format, &
          sss_data_type,   sst_data_type, ocn_data_dir, &
          oceanmixed_file, restore_sst,   trestore 
      use ice_grid, only: grid_file, kmt_file, grid_type, grid_format, &
                          gridcpl_file
      use ice_mechred, only: kstrength, krdg_partic, krdg_redist
      use ice_dyn_evp, only: ndte, kdyn, evp_damping, yield_curve, &
          maskhalo_dyn, maskhalo_stress, splitcomm_dyn
#ifdef MODAL_AER
      !mgf++
      use ice_shortwave, only: albicev, albicei, albsnowv, albsnowi, &
                               shortwave, albedo_type, R_ice, R_pnd, &
                               R_snw, dT_mlt_in, rsnw_melt_in, &
                               kaer_bc_tab, waer_bc_tab, gaer_bc_tab, bcenh
      !mgf--
#else
      use ice_shortwave, only: albicev, albicei, albsnowv, albsnowi, &
                               shortwave, albedo_type, R_ice, R_pnd, &
                               R_snw, dT_mlt_in, rsnw_melt_in
#endif
      use ice_atmo, only: atmbndy, calc_strair
      use ice_transport_driver, only: advection
      use ice_state, only: nt_Tsfc, nt_iage, nt_FY, nt_volpn, nt_aero, &
                           tr_aero, tr_iage, tr_FY, tr_pond, &
                           nt_alvl, nt_vlvl, tr_lvl, maskhalo_bound
      use ice_transport_remap, only: maskhalo_remap
      use ice_aerosol, only: restart_aero
      use ice_age, only: restart_age
      use ice_FY, only: restart_FY
      use ice_lvl, only: restart_lvl
      use ice_meltpond, only: restart_pond
      use ice_therm_vertical, only: calc_Tsfc, heat_capacity
      use ice_restoring
#ifdef MODAL_AER
      !mgf++
      use netcdf 
      !mgf--
#endif
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
        nml_error, & ! namelist i/o error flag
        ntr      , & ! counter for number of tracers turned on
        ns

      character (len=6) :: chartmp

#ifdef CCSMCOUPLED
      logical :: exists
#endif

#ifdef MODAL_AER
      !mgf++
      integer (kind=int_kind) :: fid   ! file id
      integer (kind=int_kind) :: & 
           varid,           & ! netcdf id for field
           status             ! status output from netcdf routines
      !mgf--
#endif
      character (len=256) :: fsnowoptics = 'unknown_fsnowoptics_file'

      !-----------------------------------------------------------------
      ! Namelist variables.
      !-----------------------------------------------------------------

      namelist /setup_nml/ &
        days_per_year,  year_init,      istep0,          dt,            &
        npt,            xndt_dyn,                                       &
        runtype,        runid,                                          &
        ice_ic,         restart,                                        &
        restart_dir,    restart_file,   restart_format, lcdf64,         &
        pointer_file,   dumpfreq,       dumpfreq_n,                     &
        diagfreq,       diag_type,      diag_file,                      &
        print_global,   print_points,   latpnt,          lonpnt,        &
        dbug,           histfreq,       histfreq_n,      hist_avg,      &
        history_dir,    history_file,   history_format,                 &
        write_ic,       incond_dir,     incond_file, lprint_stats

      namelist /grid_nml/ &
        grid_format,    grid_type,       grid_file,     kmt_file,       &
        kcatbound, gridcpl_file

      namelist /ice_nml/ &
        maskhalo_dyn,   maskhalo_remap,  maskhalo_bound,splitcomm_dyn,  &
        maskhalo_stress, &
        kitd,           kdyn,            ndte,                          &
        evp_damping,    yield_curve,                                    &
        kstrength,      krdg_partic,     krdg_redist,   advection,      &
        heat_capacity,  shortwave,       albedo_type,                   &
        albicev,        albicei,         albsnowv,      albsnowi,       &
        R_ice,          R_pnd,           R_snw,                         &
        dT_mlt_in,      rsnw_melt_in,                                   &
        atmbndy,        fyear_init,      ycycle,        atm_data_format,&
        atm_data_type,  atm_data_dir,    calc_strair,   calc_Tsfc,      &
        precip_units,   Tfrzpt,          update_ocn_f,                  &
        oceanmixed_ice, ocn_data_format, sss_data_type, sst_data_type,  &
        ocn_data_dir,   oceanmixed_file, restore_sst,   trestore,       &
        restore_ice,    fsnowoptics    

      namelist /tracer_nml/    &
        tr_iage, restart_age,  &
        tr_FY,   restart_FY,   &
        tr_lvl,  restart_lvl,  &
        tr_pond, restart_pond, &
        tr_aero, restart_aero         !MH for soot

      !-----------------------------------------------------------------
      ! default values
      !-----------------------------------------------------------------

      days_per_year = 365    ! number of days in a year
      year_init = 0          ! initial year
      istep0 = 0             ! no. of steps taken in previous integrations,
                             ! real (dumped) or imagined (to set calendar)
#ifndef CCSMCOUPLED
      dt = 3600.0_dbl_kind   ! time step, s 
#endif
      npt = 99999            ! total number of time steps (dt) 
      diagfreq = 24          ! how often diag output is written
      print_points = .false. ! if true, print point data
      print_global = .true.  ! if true, print global diagnostic data
      lprint_stats = .false. ! if true, prints decomposition statistics
      diag_type = 'stdout'
      diag_file = 'ice_diag.d'
      histfreq(:) = 'd'      ! output frequency option for different streams  **** note this is only for testing - needs to get changed back'
      histfreq_n(:) = 1      ! output frequency of histfreq
      hist_avg = .true.      ! if true, write time-averages (not snapshots)
      history_dir  = ' '     ! write to executable dir for default
      history_file = 'iceh'  ! history file name prefix
      history_format = 'nc'  ! file format ('bin'=binary or 'nc'=netcdf)
      write_ic = .false.     ! write out initial condition
      incond_dir = history_dir ! write to history dir for default
      incond_file = 'iceh_ic'! file prefix
      dumpfreq='y'           ! restart frequency option
      dumpfreq_n = 1         ! restart frequency
      restart = .false.      ! if true, read restart files for initialization
      restart_dir  = ' '     ! write to executable dir for default
      restart_file = 'iced'  ! restart file name prefix
      restart_format = 'nc'  ! file format ('bin'=binary or 'nc'=netcdf)
      lcdf64       = .false. ! 64 bit offset for netCDF
#ifdef CCSMCOUPLED
      pointer_file = 'rpointer.ice' // trim(inst_suffix)
#else
      pointer_file = 'ice.restart_file'
#endif
      ice_ic       = 'default'      ! latitude and sst-dependent
      grid_format  = 'bin'          ! file format ('bin'=binary or 'nc'=netcdf)
      grid_type    = 'rectangular'  ! define rectangular grid internally
      grid_file    = 'unknown_grid_file'
      kmt_file     = 'unknown_kmt_file'
      gridcpl_file = 'unknown_gridcpl_file'

      kitd = 1           ! type of itd conversions (0 = delta, 1 = linear)
      kcatbound = 1      ! category boundary formula (0 = old, 1 = new)
      kdyn = 1           ! type of dynamics (1 = evp)
      xndt_dyn = c1      ! dynamic time steps per thermodynamic time step
      ndte = 120         ! subcycles per dynamics timestep:  ndte=dt_dyn/dte
      evp_damping = .false.  ! if true, use damping procedure in evp dynamics
      yield_curve = 'ellipse'
      splitcomm_dyn = .false.    ! allow overlap of computation and halo update in subcycling
      maskhalo_dyn = .true.     ! compute masked halo before subcycling for use in halo update
      maskhalo_stress = .true.  ! compute masked halo before tripole stress update
      maskhalo_remap = .true.   ! compute masked halo before large tracer update in remap
      maskhalo_bound = .true.   ! compute masked halo before state bound
      kstrength = 1          ! 1 = Rothrock 75 strength, 0 = Hibler 79
      krdg_partic = 1        ! 1 = new participation, 0 = Thorndike et al 75
      krdg_redist = 1        ! 1 = new redistribution, 0 = Hibler 80
      advection  = 'remap'   ! incremental remapping transport scheme
      shortwave = 'default'  ! or 'dEdd' (delta-Eddington)
      albedo_type = 'default'! or 'constant'
      heat_capacity = .true. ! nonzero heat capacity (F => 0-layer thermo)
      calc_Tsfc = .true.     ! calculate surface temperature
      Tfrzpt    = 'linear_S' ! ocean freezing temperature, 'constant'=-1.8C
      update_ocn_f = .false. ! include fresh water and salt fluxes for frazil
      R_ice     = 0.00_dbl_kind   ! tuning parameter for sea ice
      R_pnd     = 0.00_dbl_kind   ! tuning parameter for ponded sea ice
      R_snw     = 0.00_dbl_kind   ! tuning parameter for snow over sea ice
      dT_mlt_in = 1.50_dbl_kind   ! melt transition (tuning)
      rsnw_melt_in = 1500._dbl_kind ! max snow grain radius (tuning)
      albicev   = 0.78_dbl_kind   ! visible ice albedo for h > ahmax
      albicei   = 0.36_dbl_kind   ! near-ir ice albedo for h > ahmax
      albsnowv  = 0.98_dbl_kind   ! cold snow albedo, visible
      albsnowi  = 0.70_dbl_kind   ! cold snow albedo, near IR
      atmbndy   = 'default'  ! or 'constant'

      fyear_init = 1900           ! first year of forcing cycle
      ycycle = 1                  ! number of years in forcing cycle
      atm_data_format = 'bin'     ! file format ('bin'=binary or 'nc'=netcdf)
      atm_data_type   = 'default'
      atm_data_dir    = ' '
      calc_strair     = .true.    ! calculate wind stress
      precip_units    = 'mks'     ! 'mm_per_month' or
                                  ! 'mm_per_sec' = 'mks' = kg/m^2 s
      oceanmixed_ice  = .false.   ! if true, use internal ocean mixed layer
      ocn_data_format = 'bin'     ! file format ('bin'=binary or 'nc'=netcdf)
      sss_data_type   = 'default'
      sst_data_type   = 'default'
      ocn_data_dir    = ' '
      oceanmixed_file = 'unknown_oceanmixed_file' ! ocean forcing data
      restore_sst     = .false.   ! restore sst if true
      trestore        = 90        ! restoring timescale, days (0 instantaneous)
      restore_ice     = .false.   ! restore ice state on grid edges if true
      dbug      = .false.         ! true writes diagnostics for input forcing

      latpnt(1) =  90._dbl_kind   ! latitude of diagnostic point 1 (deg)
      lonpnt(1) =   0._dbl_kind   ! longitude of point 1 (deg)
      latpnt(2) = -65._dbl_kind   ! latitude of diagnostic point 2 (deg)
      lonpnt(2) = -45._dbl_kind   ! longitude of point 2 (deg)

#ifndef CCSMCOUPLED
      runid   = 'unknown'   ! run ID, only used in CCSM
      runtype = 'initial'   ! run type: 'initial', 'continue'
#endif

      !-----------------------------------------------------------------
      ! extra tracers (no longer namelist variables set in ice_domain_size)  
      !-----------------------------------------------------------------

      tr_iage      = .false. ! ice age
      restart_age  = .false. ! ice age restart
      filename_iage = 'none'
      tr_FY        = .false. ! FY ice area
      restart_FY   = .false. ! FY ice area restart
      filename_FY  = 'none'
      tr_lvl       = .false. ! FY ice area
      restart_lvl  = .false. ! FY ice area restart
      filename_lvl = 'none'
      tr_pond      = .false. ! explicit melt ponds
      restart_pond = .false. ! melt ponds restart
      filename_volpn = 'none'
      tr_aero      = .false. ! aerosols          MH
      restart_aero = .false. ! aerosol restart   MH
      filename_aero = 'none'

      resttype = 'old'

      !-----------------------------------------------------------------
      ! read from input file
      !-----------------------------------------------------------------

      call get_fileunit(nu_nml)

#ifdef CCSMCOUPLED
      nml_filename  = 'ice_in'//trim(inst_suffix)
#endif
      if (my_task == master_task) then
         open (nu_nml, file=nml_filename, status='old',iostat=nml_error)
         if (nml_error /= 0) then
            nml_error = -1
         else
            nml_error =  1
         endif
         do while (nml_error > 0)
            print*,'Reading setup_nml'
            read(nu_nml, nml=setup_nml,iostat=nml_error)
            print*,'Reading grid_nml'
            read(nu_nml, nml=grid_nml,iostat=nml_error)
            print*,'Reading ice_nml'
            read(nu_nml, nml=ice_nml,iostat=nml_error)
            print*,'Reading tracer_nml'
            read(nu_nml, nml=tracer_nml,iostat=nml_error)
         end do
         if (nml_error == 0) close(nu_nml)
      endif
      call broadcast_scalar(nml_error, master_task)
      if (nml_error /= 0) then
         call abort_ice('ice: error reading namelist')
      endif

      call release_fileunit(nu_nml)

      ! Need to broadcast ice_ic at this point.

      call broadcast_scalar(ice_ic,             master_task)

      !-----------------------------------------------------------------
      ! set up diagnostics output and resolve conflicts
      !-----------------------------------------------------------------

#ifdef CCSMCOUPLED
      ! Note in CCSMCOUPLED mode diag_file is not utilized and
      ! runid and runtype are obtained from the driver, not from the namelist

      if (my_task == master_task) then
         history_file  = trim(runid) // ".cice" // trim(inst_suffix) //".h"
         restart_file  = trim(runid) // ".cice" // trim(inst_suffix) //".r"
         incond_file   = trim(runid) // ".cice" // trim(inst_suffix) //".i."
         inquire(file='ice_modelio.nml'//trim(inst_suffix),exist=exists)
         if (exists) then
            nu_diag = shr_file_getUnit()
            call shr_file_setIO('ice_modelio.nml'//trim(inst_suffix),nu_diag)
         end if
      end if
      if (trim(ice_ic) /= 'default' .and. trim(ice_ic) /= 'none') then
         restart = .true.
      end if
#else
      if (trim(diag_type) == 'file') call get_fileunit(nu_diag)
#endif


      if (my_task == master_task) then
         if (trim(diag_type) == 'file') then
            write(ice_stdout,*) 'Diagnostic output will be in file ',diag_file
            open (nu_diag, file=diag_file, status='unknown')
         endif
         write(nu_diag,*) '--------------------------------'
         write(nu_diag,*) '  CICE model diagnostic output  '
         write(nu_diag,*) '--------------------------------'
         write(nu_diag,*) ' '
      endif

      if (runtype == 'continue') restart = .true.
      if (runtype /= 'continue' .and. (restart)) then
         if (ice_ic == 'none' .or. ice_ic == 'default') then
            if (my_task == master_task) then
            write(nu_diag,*) &
            'WARNING: runtype, restart, ice_ic are inconsistent:'
            write(nu_diag,*) runtype, restart, ice_ic
            write(nu_diag,*) &
            'WARNING: Need ice_ic = <filename>.'
            write(nu_diag,*) &
            'WARNING: Initializing using ice_ic conditions'
            endif
            restart = .false.
         endif
      endif
      if (runtype == 'initial' .and. .not.(restart)) then
         if (ice_ic /= 'none' .and. ice_ic /= 'default') then
            if (my_task == master_task) then
            write(nu_diag,*) &
            'WARNING: runtype, restart, ice_ic are inconsistent:'
            write(nu_diag,*) runtype, restart, ice_ic
            write(nu_diag,*) &
            'WARNING: Initializing with NO ICE: '
            endif
            ice_ic = 'none'
         endif
      endif

#ifndef ncdf
      ! netcdf is unavailable
      history_format  = 'bin'
      restart_format  = 'bin'
      grid_format     = 'bin'
      atm_data_format = 'bin'
      ocn_data_format = 'bin'
#endif

      if (days_per_year /= 365) shortwave = 'default' ! definite conflict

      chartmp = advection(1:6)
      if (chartmp /= 'upwind' .and. chartmp /= 'remap ') advection = 'remap'

      if (ncat == 1 .and. kitd == 1) then
         if (my_task == master_task) then
            write (nu_diag,*) 'Remapping the ITD is not allowed for ncat=1'
            write (nu_diag,*) 'Using the delta function ITD option instead'
         endif
         kitd = 0
      endif

      if (trim(atm_data_type) == 'monthly' .and. calc_strair) &
         calc_strair = .false.

      if (trim(atm_data_type) == 'hadgem' .and. & 
             trim(precip_units) /= 'mks') then
         if (my_task == master_task) &
         write (nu_diag,*) &
         'WARNING: HadGEM atmospheric data chosen with wrong precip_units'
         write (nu_diag,*) &
         'WARNING: Changing precip_units to mks (i.e. kg/m2 s).'
         precip_units='mks'
      endif

      call broadcast_scalar(days_per_year,      master_task)
      call broadcast_scalar(year_init,          master_task)
      call broadcast_scalar(istep0,             master_task)
      call broadcast_scalar(dt,                 master_task)
      call broadcast_scalar(npt,                master_task)
      call broadcast_scalar(diagfreq,           master_task)
      call broadcast_scalar(print_points,       master_task)
      call broadcast_scalar(print_global,       master_task)
      call broadcast_scalar(diag_type,          master_task)
      call broadcast_scalar(diag_file,          master_task)
      call broadcast_scalar(history_format,     master_task)
!     call broadcast_scalar(histfreq,           master_task)
!     call broadcast_scalar(histfreq_n,         master_task)
      do ns=1,max_nstrm
         call broadcast_scalar(histfreq(ns),         master_task)
         call broadcast_scalar(histfreq_n(ns),       master_task)
      enddo
!     call broadcast_array(histfreq_n(:),       master_task)
      call broadcast_scalar(hist_avg,           master_task)
      call broadcast_scalar(history_dir,        master_task)
      call broadcast_scalar(write_ic,           master_task)
      call broadcast_scalar(incond_dir,         master_task)
      call broadcast_scalar(incond_file,        master_task)
      call broadcast_scalar(dumpfreq,           master_task)
      call broadcast_scalar(dumpfreq_n,         master_task)
      call broadcast_scalar(restart_file,       master_task)
      call broadcast_scalar(restart,            master_task)
      call broadcast_scalar(restart_dir,        master_task)
      call broadcast_scalar(restart_format,     master_task)
      call broadcast_scalar(lcdf64,             master_task)
#ifndef CCSM_COUPLED
      call broadcast_scalar(pointer_file,       master_task)
#endif
      call broadcast_scalar(grid_format,        master_task)
      call broadcast_scalar(grid_type,          master_task)
      call broadcast_scalar(grid_file,          master_task)
      call broadcast_scalar(kmt_file,           master_task)
      call broadcast_scalar(kitd,               master_task)
      call broadcast_scalar(gridcpl_file,       master_task)
      call broadcast_scalar(kcatbound,          master_task)
      call broadcast_scalar(kdyn,               master_task)
      call broadcast_scalar(xndt_dyn,           master_task)
      call broadcast_scalar(ndte,               master_task)
      call broadcast_scalar(evp_damping,        master_task)
      call broadcast_scalar(yield_curve,        master_task)
      call broadcast_scalar(splitcomm_dyn,      master_task)
      call broadcast_scalar(maskhalo_dyn,       master_task)
      call broadcast_scalar(maskhalo_stress,    master_task)
      call broadcast_scalar(maskhalo_remap,     master_task)
      call broadcast_scalar(maskhalo_bound,     master_task)
      call broadcast_scalar(kstrength,          master_task)
      call broadcast_scalar(krdg_partic,        master_task)
      call broadcast_scalar(krdg_redist,        master_task)
      call broadcast_scalar(advection,          master_task)
      call broadcast_scalar(shortwave,          master_task)
      call broadcast_scalar(albedo_type,        master_task)
      call broadcast_scalar(heat_capacity,      master_task)
      call broadcast_scalar(R_ice,              master_task)
      call broadcast_scalar(R_pnd,              master_task)
      call broadcast_scalar(R_snw,              master_task)
      call broadcast_scalar(dT_mlt_in,          master_task)
      call broadcast_scalar(rsnw_melt_in,       master_task)
      call broadcast_scalar(albicev,            master_task)
      call broadcast_scalar(albicei,            master_task)
      call broadcast_scalar(albsnowv,           master_task)
      call broadcast_scalar(albsnowi,           master_task)
      call broadcast_scalar(atmbndy,            master_task)
      call broadcast_scalar(fyear_init,         master_task)
      call broadcast_scalar(ycycle,             master_task)
      call broadcast_scalar(atm_data_format,    master_task)
      call broadcast_scalar(atm_data_type,      master_task)
      call broadcast_scalar(atm_data_dir,       master_task)
      call broadcast_scalar(calc_strair,        master_task)
      call broadcast_scalar(calc_Tsfc,          master_task)
      call broadcast_scalar(Tfrzpt,             master_task)
      call broadcast_scalar(update_ocn_f,       master_task)
      call broadcast_scalar(precip_units,       master_task)
      call broadcast_scalar(oceanmixed_ice,     master_task)
      call broadcast_scalar(ocn_data_format,    master_task)
      call broadcast_scalar(sss_data_type,      master_task)
      call broadcast_scalar(sst_data_type,      master_task)
      call broadcast_scalar(ocn_data_dir,       master_task)
      call broadcast_scalar(oceanmixed_file,    master_task)
      call broadcast_scalar(restore_sst,        master_task)
      call broadcast_scalar(trestore,           master_task)
      call broadcast_scalar(restore_ice,        master_task)
      call broadcast_scalar(dbug,               master_task)
      call broadcast_array (latpnt(1:2),        master_task)
      call broadcast_array (lonpnt(1:2),        master_task)
      call broadcast_scalar(runid,              master_task)
      call broadcast_scalar(runtype,            master_task)
      if (dbug) & ! else only master_task writes to file
      call broadcast_scalar(nu_diag,            master_task)
      ! tracers
      call broadcast_scalar(tr_iage,            master_task)
      call broadcast_scalar(restart_age,        master_task)
      call broadcast_scalar(tr_FY,              master_task)
      call broadcast_scalar(restart_FY,         master_task)
      call broadcast_scalar(tr_lvl,             master_task)
      call broadcast_scalar(restart_lvl,        master_task)
      call broadcast_scalar(tr_pond,            master_task)
      call broadcast_scalar(restart_pond,       master_task)
      call broadcast_scalar(tr_aero,            master_task)
      call broadcast_scalar(restart_aero,       master_task)

      !-----------------------------------------------------------------
      ! spew
      !-----------------------------------------------------------------

      if (my_task == master_task) then

         write(nu_diag,*) ' Document ice_in namelist parameters:'
         write(nu_diag,*) ' ==================================== '
         write(nu_diag,*) ' '
         if (trim(runid) /= 'unknown') &
         write(nu_diag,*)    ' runid                     = ', &
                               trim(runid)
         write(nu_diag,1030) ' runtype                   = ', &
                               trim(runtype)
         write(nu_diag,1020) ' days_per_year             = ', days_per_year
         write(nu_diag,1020) ' year_init                 = ', year_init
         write(nu_diag,1020) ' istep0                    = ', istep0
         write(nu_diag,1000) ' dt                        = ', dt
         write(nu_diag,1020) ' npt                       = ', npt
         write(nu_diag,1020) ' diagfreq                  = ', diagfreq
         write(nu_diag,1010) ' print_global              = ', &
                               print_global
         write(nu_diag,1010) ' print_points              = ', &
                               print_points
         write(nu_diag,1050) ' histfreq                  = ', histfreq(:)
         write(nu_diag,1040) ' histfreq_n                = ', histfreq_n(:)
         write(nu_diag,1010) ' hist_avg                  = ', hist_avg
         if (hist_avg) then
            write (nu_diag,*) 'History data will be averaged over ', &
                               histfreq_n,' ',histfreq
         else
            write (nu_diag,*) 'History data will be snapshots'
         endif
         write(nu_diag,*)    ' history_dir               = ', &
                               trim(history_dir)
         write(nu_diag,*)    ' history_file              = ', &
                               trim(history_file)
         write(nu_diag,*)    ' history_format            = ', &
                               trim(history_format)
         if (write_ic) then
            write (nu_diag,*) 'Initial condition will be written in ', &
                               trim(incond_dir)
         endif
         write(nu_diag,1030) ' dumpfreq                  = ', &
                               trim(dumpfreq)
         write(nu_diag,1020) ' dumpfreq_n                = ', dumpfreq_n
         write(nu_diag,1010) ' restart                   = ', restart
         write(nu_diag,*)    ' restart_dir               = ', &
                               trim(restart_dir)
         write(nu_diag,*)    ' restart_file              = ', &
                               trim(restart_file)
         write(nu_diag,*)    ' restart_format            = ', &
                               trim(restart_format)
         write(nu_diag,*)    ' lcdf64                    = ', &
                               lcdf64
         write(nu_diag,*)    ' pointer_file              = ', &
                               trim(pointer_file)
         write(nu_diag,*   ) ' ice_ic                    = ', &
                               trim(ice_ic)
         write(nu_diag,*)    ' grid_type                 = ', &
                               trim(grid_type)
         if (trim(grid_type) /= 'rectangular' .or. &
             trim(grid_type) /= 'column') then
            write(nu_diag,*) ' grid_file                 = ', &
                               trim(grid_file)
            write(nu_diag,*) ' kmt_file                  = ', &
                               trim(kmt_file)
            write(nu_diag,*) ' gridcpl_file              = ', &
                               trim(gridcpl_file)

         endif

         write(nu_diag,1020) ' kitd                      = ', kitd
         write(nu_diag,1020) ' kcatbound                 = ', &
                               kcatbound
         write(nu_diag,1020) ' kdyn                      = ', kdyn
         write(nu_diag,1000) ' xndt_dyn                  = ', xndt_dyn
         write(nu_diag,1020) ' ndte                      = ', ndte
         write(nu_diag,1010) ' evp_damping               = ', &
                               evp_damping
         write(nu_diag,*)    ' yield_curve               = ', &
                               trim(yield_curve)
         write(nu_diag,1010) ' splitcomm_dyn             = ', splitcomm_dyn
         write(nu_diag,1010) ' maskhalo_dyn              = ', maskhalo_dyn
         write(nu_diag,1010) ' maskhalo_stress           = ', maskhalo_stress
         write(nu_diag,1010) ' maskhalo_remap            = ', maskhalo_remap
         write(nu_diag,1010) ' maskhalo_bound            = ', maskhalo_bound
         write(nu_diag,1020) ' kstrength                 = ', kstrength
         write(nu_diag,1020) ' krdg_partic               = ', &
                               krdg_partic
         write(nu_diag,1020) ' krdg_redist               = ', &
                               krdg_redist
         write(nu_diag,1030) ' advection                 = ', &
                               trim(advection)
         write(nu_diag,1030) ' shortwave                 = ', &
                               trim(shortwave)
         write(nu_diag,1030) ' albedo_type               = ', &
                               trim(albedo_type)
         write(nu_diag,1000) ' R_ice                     = ', R_ice
         write(nu_diag,1000) ' R_pnd                     = ', R_pnd
         write(nu_diag,1000) ' R_snw                     = ', R_snw
         write(nu_diag,1000) ' dT_mlt_in                 = ', dT_mlt_in
         write(nu_diag,1000) ' rsnw_melt_in             = ', rsnw_melt_in
         write(nu_diag,1000) ' albicev                   = ', albicev
         write(nu_diag,1000) ' albicei                   = ', albicei
         write(nu_diag,1000) ' albsnowv                  = ', albsnowv
         write(nu_diag,1000) ' albsnowi                  = ', albsnowi
         write(nu_diag,1010) ' heat_capacity             = ', & 
                               heat_capacity
         write(nu_diag,1030) ' atmbndy                   = ', &
                               trim(atmbndy)

         write(nu_diag,1020) ' fyear_init                = ', &
                               fyear_init
         write(nu_diag,1020) ' ycycle                    = ', ycycle
         write(nu_diag,*)    ' atm_data_type             = ', &
                               trim(atm_data_type)
         write(nu_diag,1010) ' calc_strair               = ', calc_strair
         write(nu_diag,1010) ' calc_Tsfc                 = ', calc_Tsfc
         write(nu_diag,*)    ' Tfrzpt                    = ', trim(Tfrzpt)
         write(nu_diag,1010) ' update_ocn_f              = ', update_ocn_f
         if (trim(atm_data_type) /= 'default') then
            write(nu_diag,*) ' atm_data_dir              = ', &
                               trim(atm_data_dir)
            write(nu_diag,*) ' precip_units              = ', &
                               trim(precip_units)
         endif 

         write(nu_diag,1010) ' oceanmixed_ice            = ', &
                               oceanmixed_ice
         write(nu_diag,*)    ' sss_data_type             = ', &
                               trim(sss_data_type)
         write(nu_diag,*)    ' sst_data_type             = ', &
                               trim(sst_data_type)
         if (trim(sss_data_type) /= 'default' .or. &
             trim(sst_data_type) /= 'default') then
            write(nu_diag,*) ' ocn_data_dir              = ', &
                               trim(ocn_data_dir)
         endif 
         if (trim(sss_data_type) == 'ncar' .or. &
             trim(sst_data_type) == 'ncar') then
            write(nu_diag,*) ' oceanmixed_file           = ', &
                               trim(oceanmixed_file)
         endif

#ifdef coupled
         if( oceanmixed_ice ) then
            write (nu_diag,*) 'WARNING WARNING WARNING WARNING '
            write (nu_diag,*) '*Coupled and oceanmixed flags are  *'
            write (nu_diag,*) '*BOTH ON.  Ocean data received from*'
            write (nu_diag,*) '*coupler will be altered by mixed  *'
            write (nu_diag,*) '*layer routine!                    *'
            write (nu_diag,*) ' '
         endif
#endif

         write (nu_diag,*) ' '
         write (nu_diag,'(a30,2f8.2)') 'Diagnostic point 1: lat, lon =', &
                            latpnt(1), lonpnt(1)
         write (nu_diag,'(a30,2f8.2)') 'Diagnostic point 2: lat, lon =', &
                            latpnt(2), lonpnt(2)

         ! tracers
         write(nu_diag,1010) ' tr_iage                   = ', tr_iage
         write(nu_diag,1010) ' restart_age               = ', restart_age
         write(nu_diag,1010) ' tr_FY                     = ', tr_FY
         write(nu_diag,1010) ' restart_FY                = ', restart_FY
         write(nu_diag,1010) ' tr_lvl                    = ', tr_lvl
         write(nu_diag,1010) ' restart_lvl               = ', restart_lvl
         write(nu_diag,1010) ' tr_pond                   = ', tr_pond
         write(nu_diag,1010) ' restart_pond              = ', restart_pond
         write(nu_diag,1010) ' tr_aero                   = ', tr_aero      !MH
         write(nu_diag,1010) ' restart_aero              = ', restart_aero !MH

         nt_Tsfc = 1           ! index tracers, starting with Tsfc = 1
         ntrcr = 1             ! count tracers, starting with Tsfc = 1

         if (tr_iage) then
            nt_iage = ntrcr + 1
            ntrcr = ntrcr + 1
         else
            nt_iage = max_ntrcr
         endif

         if (tr_FY) then
            nt_FY = ntrcr + 1
            ntrcr = ntrcr + 1
         else
            nt_FY = max_ntrcr
         endif

         if (tr_lvl) then
            nt_alvl = ntrcr + 1
            ntrcr = ntrcr + 1
            nt_vlvl = ntrcr + 1
            ntrcr = ntrcr + 1
         else
            nt_alvl = max_ntrcr
            nt_vlvl = max_ntrcr
         endif

         if (tr_pond) then
            nt_volpn = ntrcr + 1
            ntrcr = ntrcr + 1
         else
            nt_volpn = max_ntrcr
         endif

         if (tr_aero) then
            nt_aero = ntrcr + 1
            ntrcr = ntrcr + n_aero*4 !MH 2 for snow soot and 2 for ice
                                     !MH for multiple (n_aero) aerosols
         else
            nt_aero = max_ntrcr
         endif

         if (ntrcr > max_ntrcr) then
            write(nu_diag,*) 'max_ntrcr < number of namelist tracers'
            call abort_ice('max_ntrcr < number of namelist tracers')
         endif                               

 1000    format (a30,2x,f9.2)  ! a30 to align formatted, unformatted statements
 1010    format (a30,2x,l6)    ! logical
 1020    format (a30,2x,i6)    ! integer
 1030    format (a30,   a8)    ! character
 1040    format (a30,2x,6i6)   ! integer
 1050    format (a30,   6a8)   ! character

         write (nu_diag,*) ' '
         if (grid_type  /=  'displaced_pole' .and. &
             grid_type  /=  'tripole'        .and. &
             grid_type  /=  'column'         .and. &
             grid_type  /=  'rectangular'    .and. &
             grid_type  /=  'panarctic'      .and. &
             grid_type  /=  'latlon' ) then 
            call abort_ice('ice_init: unknown grid_type')
         endif

      endif                     ! my_task = master_task

      call broadcast_scalar(ntrcr,    master_task)
      call broadcast_scalar(nt_Tsfc,  master_task)
      call broadcast_scalar(nt_iage,  master_task)
      call broadcast_scalar(nt_FY,    master_task)
      call broadcast_scalar(nt_alvl,  master_task)
      call broadcast_scalar(nt_vlvl,  master_task)
      call broadcast_scalar(nt_volpn, master_task)
      call broadcast_scalar(nt_aero,  master_task)
      call broadcast_scalar(fsnowoptics, master_task)
#ifdef MODAL_AER
      !mgf++
      if (my_task == master_task) then

         status = nf90_open(fsnowoptics, NF90_NOWRITE, fid)
         if (status .ne. NF90_NOERR) call abort_ice( 'ice: Error opening file: '//fsnowoptics)            

         status = nf90_inq_varid(fid, 'ext_cff_mss_bc_mam_cice', varid)
         if (status .ne. NF90_NOERR) call abort_ice( 'ice: Error inquiring ext_cff_mss_bc_mam_cice variable in file: '//fsnowoptics)            

         status = nf90_get_var(fid, varid, kaer_bc_tab)
         if (status .ne. NF90_NOERR) call abort_ice( 'ice: Error reading ext_cff_mss_bc_mam_cice variable in file: '//fsnowoptics)            

         status = nf90_inq_varid(fid, 'ss_alb_bc_mam_cice', varid)
         if (status .ne. NF90_NOERR) call abort_ice( 'ice: Error inquiring ss_alb_bc_mam_cice  variable in file: '//fsnowoptics)
         status = nf90_get_var(fid, varid, waer_bc_tab)
         if (status .ne. NF90_NOERR) call abort_ice( 'ice: Error reading ss_alb_bc_mam_cice variable in file: '//fsnowoptics)
         
         status = nf90_inq_varid(fid, 'asm_prm_bc_mam_cice', varid)
         if (status .ne. NF90_NOERR) call abort_ice( 'ice: Error inquiring asm_prm_bc_mam_cice variable in file: '//fsnowoptics)
         status = nf90_get_var(fid, varid, gaer_bc_tab)
         if (status .ne. NF90_NOERR) call abort_ice( 'ice: Error reading  asm_prm_bc_mam_cice variable in file: '//fsnowoptics)
         
         status = nf90_inq_varid(fid, 'bcint_enh_mam_cice', varid)
         if (status .ne. NF90_NOERR) call abort_ice( 'ice: Error inquiring bcint_enh_mam_cice variable in file: '//fsnowoptics)
         status = nf90_get_var(fid, varid, bcenh)
         if (status .ne. NF90_NOERR) call abort_ice( 'ice: Error reading  bcint_enh_mam_cice variable in file: '//fsnowoptics)
         
         status = nf90_close(fid)

         write (nu_diag,*) "MGFICE ext_cff_mss_bc = ", kaer_bc_tab(1,1), kaer_bc_tab(2,1), kaer_bc_tab(1,2), kaer_bc_tab(3,1),kaer_bc_tab(3,10) 
         write (nu_diag,*) "MGFICE ss_alb_bc = ", waer_bc_tab(1,1), waer_bc_tab(2,1), waer_bc_tab(1,2), waer_bc_tab(3,1), waer_bc_tab(3,10)
         write (nu_diag,*) "MGFICE asm_prm_bc = ", gaer_bc_tab(1,1), gaer_bc_tab(2,1), gaer_bc_tab(1,2), gaer_bc_tab(3,1), gaer_bc_tab(3,10)
         write (nu_diag,*) "MGFICE bcenh = ", bcenh(1,1,1), bcenh(1,2,1), bcenh(1,1,2), bcenh(2,1,1), bcenh(3,10,1), bcenh(3,1,8), bcenh(3,10,8)

      endif

      call broadcast_array(kaer_bc_tab, master_task)
      call broadcast_array(waer_bc_tab, master_task)
      call broadcast_array(gaer_bc_tab, master_task)
      call broadcast_array(bcenh, master_task)
      !mgf--
#endif
      end subroutine input_data

!=======================================================================
!BOP
!
! !IROUTINE: init_state - initialize ice state variables
!
! !INTERFACE:
!
      subroutine init_state
!
! !DESCRIPTION:
!
! Initialize state for the itd model
!
! !REVISION HISTORY:
!
! authors: C. M. Bitz, UW
!          William H. Lipscomb, LANL
!
! !USES:
!
      use ice_blocks
      use ice_domain
      use ice_flux, only: sst, Tf, Tair
      use ice_grid
      use ice_state
      use ice_itd
      use ice_exit
      use ice_therm_vertical, only: heat_capacity
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         it, n       , & ! tracer index
         iblk            ! block index

      !-----------------------------------------------------------------
      ! Check number of layers in ice and snow.
      !-----------------------------------------------------------------

      if (my_task == master_task) then
 
         if (nilyr < 1) then
            write (nu_diag,*) 'nilyr =', nilyr
            write (nu_diag,*) 'Must have at least one ice layer'
            call abort_ice ('ice_init: Not enough ice layers')
         endif

         if (nslyr < 1) then
            write (nu_diag,*) 'nslyr =', nslyr
            write (nu_diag,*) 'Must have at least one snow layer'
            call abort_ice('ice_init: Not enough snow layers')
         endif

         if (.not.heat_capacity) then

            write (nu_diag,*) 'WARNING - Zero-layer thermodynamics'

            if (nilyr > 1) then
               write (nu_diag,*) 'nilyr =', nilyr
               write (nu_diag,*)        &
                    'Must have nilyr = 1 if heat_capacity = F'
               call abort_ice('ice_init: Too many ice layers')
            endif

            if (nslyr > 1) then
               write (nu_diag,*) 'nslyr =', nslyr
               write (nu_diag,*)        &
                    'Must have nslyr = 1 if heat_capacity = F'
               call abort_ice('ice_init: Too many snow layers')
            endif

         endif   ! heat_capacity = F

      endif      ! my_task

      !-----------------------------------------------------------------
      ! Set tracer types
      !-----------------------------------------------------------------

      trcr_depend(nt_Tsfc)  = 0   ! ice/snow surface temperature
      if (tr_iage) trcr_depend(nt_iage)  = 1   ! volume-weighted ice age
      if (tr_FY)   trcr_depend(nt_FY)    = 0   ! area-weighted FY conc
      if (tr_lvl)  trcr_depend(nt_alvl)  = 0   ! level ice area
      if (tr_lvl)  trcr_depend(nt_vlvl)  = 1   ! level ice volume
      if (tr_pond) trcr_depend(nt_volpn) = 0   ! melt pond volume
      if (tr_aero) then
        do n=1,n_aero
         trcr_depend(nt_aero+(n-1)*4  :nt_aero+(n-1)*4+1) = 2 ! snow volume-weighted MH
         trcr_depend(nt_aero+(n-1)*4+2:nt_aero+(n-1)*4+3) = 1 ! volume-weighted MH
        enddo
      endif

      !$OMP PARALLEL DO PRIVATE(iblk,it)
      do iblk = 1, nblocks

      !-----------------------------------------------------------------
      ! Set state variables
      !-----------------------------------------------------------------

         call set_state_var (nx_block,            ny_block,            &
                             tmask(:,:,    iblk),                      &
                             ULON (:,:,    iblk), ULAT (:,:,    iblk), &
                             Tair (:,:,    iblk), sst  (:,:,    iblk), &
                             Tf   (:,:,    iblk), trcr_depend,         &
                             aicen(:,:,  :,iblk), trcrn(:,:,:,:,iblk), &
                             vicen(:,:,  :,iblk), vsnon(:,:,  :,iblk), &
                             eicen(:,:,  :,iblk), esnon(:,:,  :,iblk))

      !-----------------------------------------------------------------
      ! compute aggregate ice state and open water area
      !-----------------------------------------------------------------

         aice(:,:,iblk) = c0
         vice(:,:,iblk) = c0
         vsno(:,:,iblk) = c0
         eice(:,:,iblk) = c0
         esno(:,:,iblk) = c0
         do it = 1, max_ntrcr
            trcr(:,:,it,iblk) = c0
         enddo

         call aggregate (nx_block, ny_block,  &
                         aicen(:,:,:,iblk),   &
                         trcrn(:,:,:,:,iblk), &
                         vicen(:,:,:,iblk),   &
                         vsnon(:,:,:,iblk),   &
                         eicen(:,:,:,iblk),   &
                         esnon(:,:,:,iblk),   &
                         aice (:,:,  iblk),   &
                         trcr (:,:,:,iblk),   &
                         vice (:,:,  iblk),   &
                         vsno (:,:,  iblk),   &
                         eice (:,:,  iblk),   &
                         esno (:,:,  iblk),   &
                         aice0(:,:,  iblk),   &
                         tmask(:,:,  iblk),   &
                         ntrcr,               &
                         trcr_depend)

         aice_init(:,:,iblk) = aice(:,:,iblk)

      enddo                     ! iblk
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      ! ghost cell updates
      !-----------------------------------------------------------------

      call bound_state (aicen, trcrn, &
                        vicen, vsnon, &
                        eicen, esnon)

      end subroutine init_state

!=======================================================================
!BOP
!
! !IROUTINE: set_state_var - initialize single-category state variables
!
! !INTERFACE:
!
      subroutine set_state_var (nx_block, ny_block, &
                                tmask,    ULON,  ULAT, &
                                Tair,     sst,  &
                                Tf,       trcr_depend, &
                                aicen,    trcrn, &
                                vicen,    vsnon, &
                                eicen,    esnon) 
!
! !DESCRIPTION:
!
! Initialize state in each ice thickness category
!
! !REVISION HISTORY:
!
! authors: C. M. Bitz
!          William H. Lipscomb, LANL
!
! !USES:
!
      use ice_state, only: nt_Tsfc
      use ice_therm_vertical, only: heat_capacity, calc_Tsfc, Tmlt
      use ice_itd, only: ilyr1, slyr1, hin_max
      use ice_grid, only: grid_type
      use ice_restart, only: ice_ic
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block  ! block dimensions

      logical (kind=log_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         tmask      ! true for ice/ocean cells

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         ULON   , & ! latitude of velocity pts (radians)
         ULAT       ! latitude of velocity pts (radians)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         Tair    , & ! air temperature  (K)
         Tf      , & ! freezing temperature (C) 
         sst         ! sea surface temperature (C) 

      integer (kind=int_kind), dimension (max_ntrcr), intent(inout) :: &
         trcr_depend ! = 0 for aicen tracers, 1 for vicen, 2 for vsnon

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
         intent(out) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_ntrcr,ncat), &
         intent(out) :: &
         trcrn     ! ice tracers
                   ! 1: surface temperature of ice/snow (C)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntilyr), &
         intent(out) :: &
         eicen     ! energy of melting for each ice layer  (J/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntslyr), &
         intent(out) :: &
         esnon     ! energy of melting for each ice layer  (J/m^2)
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij          , & ! horizontal index, combines i and j loops
         k           , & ! ice layer index
         n           , & ! thickness category index
         it          , & ! tracer index
         icells          ! number of cells initialized with ice

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind) :: &
         slope, Ti, sum, hbar, &
         ainit(ncat), &
         hinit(ncat)

      real (kind=dbl_kind), parameter :: &
         hsno_init = 0.20_dbl_kind   , & ! initial snow thickness (m)
         edge_init_nh =  70._dbl_kind, & ! initial ice edge, N.Hem. (deg) 
         edge_init_sh = -60._dbl_kind    ! initial ice edge, S.Hem. (deg)


      indxi(:) = 0
      indxj(:) = 0

      ! Initialize state variables.
      ! If restarting, these values are overwritten.

      do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            aicen(i,j,n) = c0
            vicen(i,j,n) = c0
            vsnon(i,j,n) = c0
            trcrn(i,j,nt_Tsfc,n) = Tf(i,j)  ! surface temperature
            if (max_ntrcr >= 2) then
               do it = 2, max_ntrcr
                  trcrn(i,j,it,n) = c0
               enddo
            endif
         enddo
         enddo
      enddo
      eicen(:,:,:) = c0
      esnon(:,:,:) = c0

      if (trim(ice_ic) == 'default') then

      !-----------------------------------------------------------------
      ! Place ice where ocean surface is cold.
      ! Note: If SST is not read from a file, then the ocean is assumed
      !       to be at its freezing point everywhere, and ice will
      !       extend to the prescribed edges.
      !-----------------------------------------------------------------

      ! initial category areas in cells with ice
         hbar = c3  ! initial ice thickness with greatest area
                    ! Note: the resulting average ice thickness 
                    ! tends to be less than hbar due to the
                    ! nonlinear distribution of ice thicknesses 
         sum = c0
         do n = 1, ncat
            if (n < ncat) then
               hinit(n) = p5*(hin_max(n-1) + hin_max(n)) ! m
            else                ! n=ncat
               hinit(n) = (hin_max(n-1) + c1) ! m
            endif
            ! parabola, max at h=hbar, zero at h=0, 2*hbar
            ainit(n) = max(c0, (c2*hbar*hinit(n) - hinit(n)**2))
            sum = sum + ainit(n)
         enddo
         do n = 1, ncat
            ainit(n) = ainit(n) / (sum + puny/ncat) ! normalize
         enddo

         if (trim(grid_type) == 'rectangular') then

         ! place ice on left side of domain
         icells = 0
         do j = 1, ny_block
         do i = 1, nx_block
            if (tmask(i,j)) then
               if (ULON(i,j) < -50./rad_to_deg) then
                  icells = icells + 1
                  indxi(icells) = i
                  indxj(icells) = j
               endif            ! ULON
            endif               ! tmask
         enddo                  ! i
         enddo                  ! j

         else

         ! place ice at high latitudes where ocean sfc is cold
         icells = 0
         do j = 1, ny_block
         do i = 1, nx_block
            if (tmask(i,j)) then
               ! place ice in high latitudes where ocean sfc is cold
               if ( (sst (i,j) <= Tf(i,j)+p2) .and. &
                    (ULAT(i,j) < edge_init_sh/rad_to_deg .or. &
                     ULAT(i,j) > edge_init_nh/rad_to_deg) ) then
                  icells = icells + 1
                  indxi(icells) = i
                  indxj(icells) = j
               endif            ! cold surface
            endif               ! tmask
         enddo                  ! i
         enddo                  ! j

         endif                  ! rectgrid

         do n = 1, ncat

            ! ice volume, snow volume
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               aicen(i,j,n) = ainit(n)
               vicen(i,j,n) = hinit(n) * ainit(n) ! m
               vsnon(i,j,n) =min(aicen(i,j,n)*hsno_init,p2*vicen(i,j,n))
            enddo               ! ij

            ! surface temperature
            if (calc_Tsfc) then
        
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  trcrn(i,j,nt_Tsfc,n) = min(Tsmelt, Tair(i,j) - Tffresh) !deg C
               enddo

            else    ! Tsfc is not calculated by the ice model

               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  trcrn(i,j,nt_Tsfc,n) = Tf(i,j)   ! not used
               enddo

            endif       ! calc_Tsfc

            ! other tracers (none at present)

            if (heat_capacity) then

               ! ice energy
               do k = 1, nilyr
                  do ij = 1, icells
                     i = indxi(ij)
                     j = indxj(ij)

                     ! assume linear temp profile and compute enthalpy
                     slope = Tf(i,j) - trcrn(i,j,nt_Tsfc,n)
                     Ti = trcrn(i,j,nt_Tsfc,n) &
                        + slope*(real(k,kind=dbl_kind)-p5) &
                                /real(nilyr,kind=dbl_kind)

                     eicen(i,j,ilyr1(n)+k-1) = &
                          -(rhoi * (cp_ice*(Tmlt(k)-Ti) &
                          + Lfresh*(c1-Tmlt(k)/Ti) - cp_ocn*Tmlt(k))) &
                          * vicen(i,j,n)/real(nilyr,kind=dbl_kind)
                  enddo            ! ij
               enddo               ! nilyr

               ! snow energy
               do k = 1, nslyr
                  do ij = 1, icells
                     i = indxi(ij)
                     j = indxj(ij)

                     Ti = min(c0, trcrn(i,j,nt_Tsfc,n))
                     esnon(i,j,slyr1(n)+k-1) = -rhos*(Lfresh - cp_ice*Ti) &
                                               *vsnon(i,j,n) &
                                               /real(nslyr,kind=dbl_kind)
                  enddo            ! ij
               enddo               ! nslyr

            else  ! one layer with zero heat capacity

               ! ice energy
               k = 1

               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  eicen(i,j,ilyr1(n)+k-1) = &
                          - rhoi * Lfresh * vicen(i,j,n)
               enddo            ! ij

               ! snow energy
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)
                  esnon(i,j,slyr1(n)+k-1) = & 
                          - rhos * Lfresh * vsnon(i,j,n)
               enddo            ! ij

            endif               ! heat_capacity
         enddo                  ! ncat
      endif                     ! ice_ic

      end subroutine set_state_var

!=======================================================================

      end module ice_init

!=======================================================================
