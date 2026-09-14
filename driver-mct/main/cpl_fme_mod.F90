module cpl_fme_mod

  !-----------------------------------------------------------------------------
  ! Coupler-native FME (Full Model Emulation) output.
  !
  ! Emits the merged, post-coupler forcing that each component actually
  ! integrates (x2o / x2i / xao) as horizontally-remapped lat-lon monthly
  ! NetCDF files.  This is the one quantity no component history tape can
  ! provide: it exists only on the mediator after the coupler merge.
  !
  ! This module mirrors the MPAS-Ocean FME analysis-member pattern
  ! (components/mpas-ocean/src/shared/mpas_ocn_fme_horiz_remap.F):
  !   init -> accumulate -> remap -> direct-PIO monthly write + sidecar restart
  ! but is component-agnostic and lives entirely on the coupler PEs.
  !
  ! Reuse, not reinvention:
  !   * Horizontal remap: shared share/util/shr_horiz_remap_mod.F90 (the same
  !     CRS-SpMV + masked renormalization used by EAM and the MPAS AMs).  No
  !     new map files -- the coupler ocean/ice bundles live on the same
  !     IcoswISC30E3r5 decomposition we already have a map for.
  !   * Accumulation: the merged ocean forcing is accumulated in SOURCE
  !     (ocean-cell) space as an mct_aVect, exactly like prep_ocn_mod's
  !     x2oacc_ox, and remapped only at flush.  For a linear masked remap with
  !     a static land mask this is bit-equivalent to accumulating in remap
  !     space, but (a) it reuses the proven x2oacc persistence machinery for
  !     restart-BFB trust, (b) it is immune to the map-change contamination of
  !     AGENTS.md gotcha #46, and (c) it is cheaper (one remap per field per
  !     output window instead of per coupling step).
  !   * Restart: a per-stream sidecar file persists the source-space
  !     accumulator + sample count + drift-free output schedule, written on
  !     the coupler restart alarm and read back at init.  Honors gotcha #11
  !     (advance schedule by exactly one interval -- no drift), #29 (append-
  !     mode reopen of the straddling monthly file with frame overwrite), and
  !     #45 (rewind/branch restart -> reset schedule + clobber-create).
  !   * Remap fill flow: apply_masked with the level-numlev+1 mask packing
  !     (gotcha #3), SHR_FILL_VALUE flow (gotcha #5), real(4)/PIO_REAL output
  !     (gotcha #1).
  !
  ! Purely additive: nothing here touches the EAM h0 tape or the MPAS FME AMs.
  !-----------------------------------------------------------------------------

  use shr_kind_mod,        only: r8 => shr_kind_r8, r4 => shr_kind_r4, &
                                 CL => shr_kind_cl, CS => shr_kind_cs, &
                                 CXX => shr_kind_cxx
  use shr_sys_mod,         only: shr_sys_abort, shr_sys_flush
  use shr_cal_mod,         only: shr_cal_date2ymd
  use shr_pio_mod,         only: shr_pio_getiosys, shr_pio_getiotype
  use shr_horiz_remap_mod, only: shr_horiz_remap_t, SHR_FILL_VALUE
  use mct_mod
  use seq_comm_mct,        only: CPLID, seq_comm_setptrs, seq_comm_iamin, logunit
  use seq_timemgr_mod,     only: seq_timemgr_EClockGetData
  use seq_infodata_mod,    only: seq_infodata_type, seq_infodata_GetData
  use seq_io_mod,          only: seq_io_date2yyyymmdd, seq_io_sec2hms
  use seq_flds_mod,        only: seq_flds_lookup
  use component_type_mod,  only: component_type, component_get_gsmap_cx, &
                                 component_get_x2c_cx, component_get_c2x_cx
  use prep_ocn_mod,        only: prep_ocn_get_x2oacc_ox
  use ESMF,                only: ESMF_Clock
  use pio,                 only: file_desc_t, var_desc_t, io_desc_t, &
                                 iosystem_desc_t, &
                                 pio_createfile, pio_openfile, pio_closefile, &
                                 pio_def_dim, pio_def_var, pio_enddef, &
                                 pio_put_att, pio_get_att, pio_put_var, &
                                 pio_get_var, pio_initdecomp, pio_freedecomp, &
                                 pio_write_darray, pio_read_darray, &
                                 pio_setframe, pio_inq_dimid, pio_inq_varid, &
                                 pio_inq_dimlen, &
                                 PIO_CLOBBER, PIO_WRITE, PIO_NOWRITE, &
                                 PIO_REAL, PIO_DOUBLE, PIO_INT, PIO_GLOBAL, &
                                 PIO_UNLIMITED, PIO_OFFSET_KIND, PIO_NOERR

  implicit none
  private
  save

  public :: cpl_fme_init          ! one-time setup (+ warm-restart read)
  public :: cpl_fme_accum         ! sample the merged forcing (per coupling phase)
  public :: cpl_fme_restart_write ! write accumulator sidecars on restart alarm
  public :: cpl_fme_final         ! cleanup

  !-----------------------------------------------------------------------------
  ! Coupling-sequence phase at which a stream's source bundle is validly populated
  ! (gotcha #51).  The merged IMPORT bundles must be sampled right after their
  ! own merge -- x2i after prep_ice_mrg, x2a after prep_atm_mrg -- NOT at the
  ! ocean hook, where x2i is still the previous step's value (merge runs later)
  ! and so reads a re-initialized ZERO on the first step of a warm-restart leg,
  ! polluting that window's mean (one bad sample / nAccum).  Everything else
  ! (x2o from the restart-safe x2oacc window-mean, xao recomputed each step,
  ! and the o2x/i2x/a2x component EXPORTS) is valid at the ocean hook.
  integer, parameter, public :: CPL_FME_PHASE_OCN = 1  ! after prep_ocn_accum_avg
  integer, parameter, public :: CPL_FME_PHASE_ICE = 2  ! after prep_ice_mrg
  integer, parameter, public :: CPL_FME_PHASE_ATM = 3  ! after prep_atm_mrg

  !-----------------------------------------------------------------------------
  ! Per-stream monthly lat-lon PIO output file (2D, time-averaged records).
  ! A simplified analog of ocn_fme_remap_file_t: 2D vars only, no scalars,
  ! no static masks, no instantaneous companions (the cpl FME tape is
  ! averaged-only by spec).
  !-----------------------------------------------------------------------------
  integer, parameter :: max_fme_vars = 128

  type :: cpl_fme_file_t
     type(file_desc_t)  :: pio_file
     logical            :: is_open = .false.
     integer            :: dim_lon = 0, dim_lat = 0, dim_time = 0, dim_nbnd = 0
     type(var_desc_t)   :: var_lon, var_lat
     type(var_desc_t)   :: var_time, var_tbnds, var_date, var_datesec
     integer            :: nvars = 0
     type(var_desc_t)   :: var_descs(max_fme_vars)
     integer            :: time_record = 0       ! frames already written
     character(len=CL)  :: current_filename = ''
     logical            :: reopened = .false.    ! append-mode this open?
  end type cpl_fme_file_t

  !-----------------------------------------------------------------------------
  ! A coupler FME stream: one merged-forcing bundle, one cadence, one map.
  !-----------------------------------------------------------------------------
  type :: cpl_fme_stream_t
     logical            :: enabled  = .false.
     character(len=CS)  :: name     = ''         ! e.g. 'x2o1D'
     character(len=CS)  :: bundle   = ''         ! 'x2o' | 'xao' | 'x2i'
     character(len=CL)  :: mapfile  = ''
     real(r8)           :: interval = 1.0_r8     ! output window, days
     integer            :: phase    = CPL_FME_PHASE_OCN  ! when to sample (gotcha #51)

     ! live source bundle + its decomposition (resolved at setup)
     type(mct_aVect), pointer :: src_bundle => null()
     type(mct_gsMap), pointer :: src_gsmap  => null()
     integer            :: src_lsize = 0
     integer            :: src_gsize = 0

     ! field selection (resolved against the live bundle at setup)
     integer            :: nflds = 0
     character(len=64)  :: fldnames(max_fme_vars) = ''
     integer            :: fld_kidx(max_fme_vars) = 0   ! index into source aVect
     logical            :: fld_is_state(max_fme_vars) = .false.  ! S* snapshot vs F* mean

     ! horizontal remap (shared module) + send-side packing
     type(shr_horiz_remap_t) :: map
     integer            :: n_send = 0
     integer, allocatable :: send_local_idx(:)         ! aVect col per send slot
     real(r8), allocatable :: send_buf(:)              ! packed (val,mask) buffer

     ! source-space accumulator (ocean-cell space), drift-free schedule
     type(mct_aVect)    :: accum                       ! running sum of window-means
     integer            :: nAccum = 0
     real(r8)           :: last_output_time = 0.0_r8   ! days since epoch
     real(r8)           :: next_output_time = 0.0_r8   ! days since epoch

     ! target-grid PIO decomposition (built once)
     type(io_desc_t)    :: iodesc_dst
     logical            :: iodesc_valid = .false.

     ! monthly output file
     type(cpl_fme_file_t) :: outfile
     logical            :: reopen_pending = .true.     ! gotcha #29/#44 latch
  end type cpl_fme_stream_t

  integer, parameter :: max_streams = 16
  type(cpl_fme_stream_t) :: streams(max_streams)
  integer :: nstreams = 0

  ! module-scope coupler context
  logical                        :: initialized   = .false.
  logical                        :: active        = .false.  ! any stream on?
  type(iosystem_desc_t), pointer :: cpl_iosys => null()
  integer :: cpl_iotype = 0
  integer :: cpl_mpicom = 0, cpl_iam = -1, cpl_npes = 0
  logical :: iamin_cpl = .false.

  ! restart / schedule state
  logical :: is_restart_run = .false.   ! gotcha #29
  logical :: rewind_restart = .false.   ! gotcha #45

  ! CF time epoch (= case RUN_STARTDATE) and provenance
  integer           :: epoch_ymd = 0, epoch_tod = 0
  character(len=CL) :: time_units = ''
  character(len=CL) :: calendar   = ''
  character(len=CL) :: case_name  = ''
  character(len=CL) :: model_version = ''

  ! source decompositions (coupler) + live merged-forcing bundle pointers.
  ! ocn and ice share the same physical mesh (so the same map file), but have
  ! distinct gsmaps.  Bundle pointers are resolved once at init; the targets
  ! are stable while the merged DATA is updated in place each coupling step.
  type(mct_gsMap), pointer :: ocn_gsmap => null()
  type(mct_gsMap), pointer :: ice_gsmap => null()
  type(mct_gsMap), pointer :: atm_gsmap => null()
  integer :: lsize_ocn = 0, gsize_ocn = 0
  integer :: lsize_ice = 0, gsize_ice = 0
  integer :: lsize_atm = 0, gsize_atm = 0
  type(mct_aVect), pointer :: b_x2o => null()   ! ocn merged window-mean (x2oacc)
  type(mct_aVect), pointer :: b_xao => null()   ! atm-ocn bulk fluxes (Faox_*)
  type(mct_aVect), pointer :: b_x2i => null()   ! ice merged import (x2i)
  type(mct_aVect), pointer :: b_o2x => null()   ! ocn export (Samudra state out)
  type(mct_aVect), pointer :: b_i2x => null()   ! ice export
  type(mct_aVect), pointer :: b_x2a => null()   ! atm merged import (ACE forcing)
  type(mct_aVect), pointer :: b_a2x => null()   ! atm export (ACE forcing out)

  real(r8), parameter :: SCHED_EPS  = 1.0e-9_r8   ! day-boundary tolerance
  real(r8), parameter :: FILL_DETECT = SHR_FILL_VALUE * 0.1_r8

  ! authoritative merged-forcing field lists (AGENTS.md PLAN section;
  ! seq_flds_mod.F90).  Absent fields are skipped+warned at setup, so these
  ! are safe across compsets.  Field names are length-16 (longest is 12).
  integer, parameter :: n_x2o_default = 17
  character(len=16), parameter :: x2o_default(n_x2o_default) = (/ &
       'Foxx_taux       ', 'Foxx_tauy       ', 'Foxx_sen        ', &
       'Foxx_lat        ', 'Foxx_lwup       ', 'Foxx_evap       ', &
       'Foxx_swnet      ', 'Faxa_lwdn       ', 'Faxa_rain       ', &
       'Faxa_snow       ', 'Fioi_melth      ', 'Fioi_meltw      ', &
       'Fioi_salt       ', 'Foxx_rofl       ', 'Foxx_rofi       ', &
       'Sa_pslv         ', 'Si_ifrac        ' /)

  integer, parameter :: n_xao_default = 13
  character(len=16), parameter :: xao_default(n_xao_default) = (/ &
       'Faox_taux       ', 'Faox_tauy       ', 'Faox_lat        ', &
       'Faox_sen        ', 'Faox_lwup       ', 'Faox_evap       ', &
       'Faox_swdn       ', 'Faox_swup       ', 'So_tref         ', &
       'So_qref         ', 'So_u10          ', 'So_ustar        ', &
       'So_duu10n       ' /)

  integer, parameter :: n_x2i_default = 21
  character(len=16), parameter :: x2i_default(n_x2i_default) = (/ &
       'Sa_z            ', 'Sa_u            ', 'Sa_v            ', &
       'Sa_tbot         ', 'Sa_ptem         ', 'Sa_shum         ', &
       'Sa_pbot         ', 'Sa_dens         ', 'Faxa_rain       ', &
       'Faxa_snow       ', 'Faxa_lwdn       ', 'Faxa_swndr      ', &
       'Faxa_swvdr      ', 'Faxa_swndf      ', 'Faxa_swvdf      ', &
       'So_t            ', 'So_s            ', 'So_u            ', &
       'So_v            ', 'So_dhdx         ', 'Fioo_q          ' /)

  ! o2x : ocean export (what the ocean sends up -- Samudra's prognostic state)
  integer, parameter :: n_o2x_default = 10
  character(len=16), parameter :: o2x_default(n_o2x_default) = (/ &
       'So_t            ', 'So_s            ', 'So_u            ', &
       'So_v            ', 'So_dhdx         ', 'So_dhdy         ', &
       'So_ssh          ', 'So_bldepth      ', 'So_fswpen       ', &
       'Fioo_q          ' /)

  ! i2x : sea-ice export (ice fraction + albedos + ice->atm/ocn fluxes)
  integer, parameter :: n_i2x_default = 20
  character(len=16), parameter :: i2x_default(n_i2x_default) = (/ &
       'Si_ifrac        ', 'Si_t            ', 'Si_avsdr        ', &
       'Si_anidr        ', 'Si_avsdf        ', 'Si_anidf        ', &
       'Si_tref         ', 'Si_qref         ', 'Si_u10          ', &
       'Si_snowh        ', 'Faii_taux       ', 'Faii_tauy       ', &
       'Faii_lat        ', 'Faii_sen        ', 'Faii_lwup       ', &
       'Faii_evap       ', 'Fioi_melth      ', 'Fioi_meltw      ', &
       'Fioi_salt       ', 'Fioi_swpen      ' /)

  ! x2a : atm import (the atmosphere's lower-boundary forcing -- ACE inputs).
  ! On the ATM grid (ne30pg2) -> needs the atm map file, not the ocean map.
  integer, parameter :: n_x2a_default = 19
  character(len=16), parameter :: x2a_default(n_x2a_default) = (/ &
       'Sf_lfrac        ', 'Sf_ifrac        ', 'Sf_ofrac        ', &
       'Sx_avsdr        ', 'Sx_anidr        ', 'Sx_avsdf        ', &
       'Sx_anidf        ', 'Sx_tref         ', 'Sx_qref         ', &
       'So_t            ', 'Sx_t            ', 'Sx_u10          ', &
       'Faxx_taux       ', 'Faxx_tauy       ', 'Faxx_lat        ', &
       'Faxx_sen        ', 'Faxx_lwup       ', 'Faxx_evap       ', &
       'So_ssq          ' /)

  ! a2x: atm export to the coupler (what the atm sends out each coupling step).
  ! States (Sa_*) are written as instantaneous snapshots at the record time;
  ! downwelling/precip flux densities (Faxa_*) are written as window means --
  ! see cpl_fme_accum/flush and the per-field cell_methods.  Absent fields are
  ! skipped at setup (perrWith='quiet').
  integer, parameter :: n_a2x_default = 18
  character(len=16), parameter :: a2x_default(n_a2x_default) = (/ &
       'Sa_z            ', 'Sa_u            ', 'Sa_v            ', &
       'Sa_tbot         ', 'Sa_ptem         ', 'Sa_shum         ', &
       'Sa_pbot         ', 'Sa_pslv         ', 'Sa_dens         ', &
       'Faxa_lwdn       ', 'Faxa_rainc      ', 'Faxa_rainl      ', &
       'Faxa_snowc      ', 'Faxa_snowl      ', 'Faxa_swndr      ', &
       'Faxa_swvdr      ', 'Faxa_swndf      ', 'Faxa_swvdf      ' /)

CONTAINS

  !=============================================================================
  subroutine cpl_fme_init(infodata, EClock, ocn, ocn_present, ice, ice_present, &
       atm, atm_present, read_restart, ocn_cpl_dt)
    !---------------------------------------------------------------------------
    ! One-time setup.  Reads the cpl_fme_inparm namelist group from drv_in,
    ! loads maps, builds remap comm patterns + decompositions, allocates the
    ! source-space accumulators, and (on a warm restart) restores accumulator
    ! state from the per-stream sidecars.  Inert (no streams) when the
    ! namelist group is absent or all streams are disabled -- purely additive.
    ! ocn_cpl_dt (ocean coupling step, seconds) is used only for a cadence
    ! sanity WARN -- see the OCN_NCPL check at the end of this routine.
    !---------------------------------------------------------------------------
    use prep_aoflux_mod, only: prep_aoflux_get_xao_ox

    type(seq_infodata_type), intent(inout) :: infodata
    type(ESMF_Clock),        intent(in)    :: EClock
    type(component_type),    intent(in)    :: ocn(:)
    logical,                 intent(in)    :: ocn_present
    type(component_type),    intent(in)    :: ice(:)
    logical,                 intent(in)    :: ice_present
    type(component_type),    intent(in)    :: atm(:)
    logical,                 intent(in)    :: atm_present
    logical,                 intent(in)    :: read_restart
    integer,                 intent(in)    :: ocn_cpl_dt

    integer :: s, ierr, base_dt
    type(mct_aVect), pointer :: avp(:)
    character(len=*), parameter :: subname = '(cpl_fme_init)'

    if (initialized) return
    initialized = .true.

    iamin_cpl = seq_comm_iamin(CPLID)
    if (.not. iamin_cpl) return

    ! all current streams ride the ocn/ice merged forcing; nothing without ocn
    if (.not. ocn_present) return

    ! coupler PIO + MPI context
    cpl_iosys  => shr_pio_getiosys(CPLID)
    cpl_iotype =  shr_pio_getiotype(CPLID)
    call seq_comm_setptrs(CPLID, iam=cpl_iam, mpicom=cpl_mpicom, npes=cpl_npes)

    is_restart_run = read_restart

    ! CF epoch (= RUN_STARTDATE) + calendar, identical to seq_hist aux files
    call seq_timemgr_EClockGetData(EClock, start_ymd=epoch_ymd, &
         start_tod=epoch_tod, calendar=calendar)
    time_units = 'days since ' // trim(seq_io_date2yyyymmdd(epoch_ymd)) &
         // ' ' // seq_io_sec2hms(epoch_tod)

    call seq_infodata_GetData(infodata, case_name=case_name, &
         model_version=model_version)

    ! ocean source decomposition + live bundle pointers (stable targets).
    ! x2o : merged window-mean (import); o2x : ocean export (state out).
    ocn_gsmap => component_get_gsmap_cx(ocn(1))
    gsize_ocn =  mct_gsMap_gsize(ocn_gsmap)
    lsize_ocn =  mct_gsMap_lsize(ocn_gsmap, cpl_mpicom)
    avp => prep_ocn_get_x2oacc_ox()   ; if (associated(avp)) b_x2o => avp(1)
    avp => prep_aoflux_get_xao_ox()   ; if (associated(avp)) b_xao => avp(1)
    b_o2x => component_get_c2x_cx(ocn(1))

    ! ice source decomposition + merged-import + export bundles (shares the
    ! ocean mesh, hence the same map file, but a distinct gsmap)
    if (ice_present) then
       ice_gsmap => component_get_gsmap_cx(ice(1))
       gsize_ice =  mct_gsMap_gsize(ice_gsmap)
       lsize_ice =  mct_gsMap_lsize(ice_gsmap, cpl_mpicom)
       b_x2i     => component_get_x2c_cx(ice(1))
       b_i2x     => component_get_c2x_cx(ice(1))
    end if

    ! atm source decomposition + merged-import bundle (atm grid -> atm map)
    if (atm_present) then
       atm_gsmap => component_get_gsmap_cx(atm(1))
       gsize_atm =  mct_gsMap_gsize(atm_gsmap)
       lsize_atm =  mct_gsMap_lsize(atm_gsmap, cpl_mpicom)
       b_x2a     => component_get_x2c_cx(atm(1))   ! merged import (into atm)
       b_a2x     => component_get_c2x_cx(atm(1))   ! atm export (out of atm)
    end if

    ! read namelist + populate streams(:) (default: all disabled)
    call cpl_fme_read_namelist()

    active = .false.
    do s = 1, nstreams
       if (.not. streams(s)%enabled) cycle
       call cpl_fme_stream_setup(streams(s), EClock, ierr)
       if (ierr /= 0) then
          if (cpl_iam == 0) write(logunit,*) subname, &
               ': WARNING disabling stream ', trim(streams(s)%name), &
               ' (setup ierr=', ierr, ')'
          streams(s)%enabled = .false.
          cycle
       end if
       active = .true.
    end do

    ! warm-restart: restore accumulator + schedule from sidecars
    if (is_restart_run .and. active) then
       do s = 1, nstreams
          if (streams(s)%enabled) call cpl_fme_read_sidecar(streams(s), EClock)
       end do
    end if

    if (cpl_iam == 0 .and. active) then
       write(logunit,*) subname, ': coupler-native FME output ACTIVE, ', &
            'nstreams=', count(streams(1:nstreams)%enabled)
    end if

    ! Cadence sanity (gotcha #52): the cpl-FME reductions assume the ocean
    ! couples every driver step (OCN_NCPL == ATM_NCPL).  If the ocean couples
    ! COARSER than the base step, the coupler's x2oacc accumulator pre-averages
    ! the ocean forcing over its window, so an x2o STATE snapshot becomes the
    ! last ocean-coupling-window mean (not a true instant), and the x2a/a2x
    ! window means undersample.  Output is still produced -- this is a loud
    ! heads-up, not a fatal error.
    if (active .and. cpl_iam == 0) then
       call seq_timemgr_EClockGetData(EClock, dtime=base_dt)
       if (ocn_cpl_dt > base_dt) then
          write(logunit,*) subname, ': WARNING OCN couples coarser than the ', &
               'base step (ocn_cpl_dt=', ocn_cpl_dt, 's > base dt=', base_dt, &
               's, i.e. OCN_NCPL < ATM_NCPL). cpl-FME assumes equal cadence: ', &
               'x2o STATE fields (Sa_pslv, Si_ifrac) snapshot the last ', &
               'ocean-coupling-window mean rather than a true instant, and ', &
               'x2a/a2x means undersample. See AGENTS.md gotcha #52.'
       end if
    end if

  end subroutine cpl_fme_init

  !=============================================================================
  subroutine cpl_fme_read_namelist()
    !---------------------------------------------------------------------------
    ! Read the &cpl_fme_inparm group from drv_in.  All PEs read the file
    ! (cheap; matches other driver namelist reads).  Missing group => no
    ! streams (inert).  Populates streams(:) for the ocean 1D and 5D cases;
    ! xao / x2i streams are added in a later phase.
    !---------------------------------------------------------------------------
    integer :: unitn, ios
    logical :: exists

    ! &cpl_fme_inparm namelist variables.  ocn/ice/xao/o2x/i2x share the ocean
    ! mesh map (cpl_fme_ocn_map); x2a uses the atm-grid map (cpl_fme_atm_map).
    ! One shared pair of windows; per-bundle enable flags.
    logical           :: cpl_fme_x2o_enable, cpl_fme_x2o_5d_enable
    logical           :: cpl_fme_xao_enable, cpl_fme_xao_5d_enable
    logical           :: cpl_fme_x2i_enable, cpl_fme_x2i_5d_enable
    logical           :: cpl_fme_o2x_enable, cpl_fme_o2x_5d_enable
    logical           :: cpl_fme_i2x_enable, cpl_fme_i2x_5d_enable
    logical           :: cpl_fme_x2a_enable, cpl_fme_a2x_enable
    character(len=CL) :: cpl_fme_ocn_map, cpl_fme_atm_map
    real(r8)          :: cpl_fme_ocn_interval, cpl_fme_ocn_5d_interval
    real(r8)          :: cpl_fme_atm_interval

    namelist /cpl_fme_inparm/ cpl_fme_x2o_enable, cpl_fme_x2o_5d_enable, &
         cpl_fme_xao_enable, cpl_fme_xao_5d_enable, &
         cpl_fme_x2i_enable, cpl_fme_x2i_5d_enable, &
         cpl_fme_o2x_enable, cpl_fme_o2x_5d_enable, &
         cpl_fme_i2x_enable, cpl_fme_i2x_5d_enable, &
         cpl_fme_x2a_enable, cpl_fme_a2x_enable, &
         cpl_fme_ocn_map, cpl_fme_atm_map, &
         cpl_fme_ocn_interval, cpl_fme_ocn_5d_interval, &
         cpl_fme_atm_interval

    ! defaults: everything off
    cpl_fme_x2o_enable      = .false. ; cpl_fme_x2o_5d_enable = .false.
    cpl_fme_xao_enable      = .false. ; cpl_fme_xao_5d_enable = .false.
    cpl_fme_x2i_enable      = .false. ; cpl_fme_x2i_5d_enable = .false.
    cpl_fme_o2x_enable      = .false. ; cpl_fme_o2x_5d_enable = .false.
    cpl_fme_i2x_enable      = .false. ; cpl_fme_i2x_5d_enable = .false.
    cpl_fme_x2a_enable      = .false. ; cpl_fme_a2x_enable    = .false.
    cpl_fme_ocn_map         = ''
    cpl_fme_atm_map         = ''
    cpl_fme_ocn_interval    = 1.0_r8
    cpl_fme_ocn_5d_interval = 5.0_r8
    cpl_fme_atm_interval    = 0.25_r8   ! 6 h (atm import/export cadence)

    inquire(file='drv_in', exist=exists)
    if (exists) then
       open(newunit=unitn, file='drv_in', status='old', iostat=ios)
       if (ios == 0) then
          read(unitn, nml=cpl_fme_inparm, iostat=ios)  ! ios/=0 if group absent
          close(unitn)
       end if
    end if

    nstreams = 0
    ! ocean merged forcing x2o (the prize): 1D + 5D
    call add_stream(cpl_fme_x2o_enable,    'x2o1D', 'x2o', cpl_fme_ocn_map, cpl_fme_ocn_interval)
    call add_stream(cpl_fme_x2o_5d_enable, 'x2o5D', 'x2o', cpl_fme_ocn_map, cpl_fme_ocn_5d_interval)
    ! atm-ocean bulk fluxes xao (Faox_*): 1D + 5D
    call add_stream(cpl_fme_xao_enable,    'xao1D', 'xao', cpl_fme_ocn_map, cpl_fme_ocn_interval)
    call add_stream(cpl_fme_xao_5d_enable, 'xao5D', 'xao', cpl_fme_ocn_map, cpl_fme_ocn_5d_interval)
    ! sea-ice merged import x2i: 1D + 5D
    call add_stream(cpl_fme_x2i_enable,    'x2i1D', 'x2i', cpl_fme_ocn_map, cpl_fme_ocn_interval)
    call add_stream(cpl_fme_x2i_5d_enable, 'x2i5D', 'x2i', cpl_fme_ocn_map, cpl_fme_ocn_5d_interval)
    ! ocean export o2x (Samudra state out): 1D + 5D
    call add_stream(cpl_fme_o2x_enable,    'o2x1D', 'o2x', cpl_fme_ocn_map, cpl_fme_ocn_interval)
    call add_stream(cpl_fme_o2x_5d_enable, 'o2x5D', 'o2x', cpl_fme_ocn_map, cpl_fme_ocn_5d_interval)
    ! sea-ice export i2x: 1D + 5D
    call add_stream(cpl_fme_i2x_enable,    'i2x1D', 'i2x', cpl_fme_ocn_map, cpl_fme_ocn_interval)
    call add_stream(cpl_fme_i2x_5d_enable, 'i2x5D', 'i2x', cpl_fme_ocn_map, cpl_fme_ocn_5d_interval)
    ! atm streams (ACE boundary forcing in + out), ATM grid -> atm map, at a
    ! configurable sub-daily cadence (cpl_fme_atm_interval, default 6 h):
    !   x2a = coupler merged import the atm integrates;
    !   a2x = atm export the coupler redistributes to ocn/ice.
    call add_stream(cpl_fme_x2a_enable, 'x2a', 'x2a', cpl_fme_atm_map, cpl_fme_atm_interval)
    call add_stream(cpl_fme_a2x_enable, 'a2x', 'a2x', cpl_fme_atm_map, cpl_fme_atm_interval)

  contains
    subroutine add_stream(en, nm, bn, mf, iv)
      logical,          intent(in) :: en
      character(len=*), intent(in) :: nm, bn, mf
      real(r8),         intent(in) :: iv
      nstreams = nstreams + 1
      streams(nstreams)%enabled  = en
      streams(nstreams)%name     = nm
      streams(nstreams)%bundle   = bn
      streams(nstreams)%mapfile  = mf
      streams(nstreams)%interval = iv
      ! Sample the merged IMPORT bundles right after their own merge (gotcha
      ! #51): x2i after prep_ice_mrg, x2a after prep_atm_mrg.  All other
      ! bundles are valid at the ocean hook.
      select case (trim(bn))
      case ('x2i') ; streams(nstreams)%phase = CPL_FME_PHASE_ICE
      case ('x2a') ; streams(nstreams)%phase = CPL_FME_PHASE_ATM
      case default ; streams(nstreams)%phase = CPL_FME_PHASE_OCN
      end select
    end subroutine add_stream

  end subroutine cpl_fme_read_namelist

  !=============================================================================
  subroutine cpl_fme_stream_setup(st, EClock, ierr)
    !---------------------------------------------------------------------------
    ! Load the map, build the source-gather comm pattern + send-side packing
    ! map, build the target-grid PIO decomposition, allocate the source-space
    ! accumulator, and initialize the drift-free output schedule.
    !---------------------------------------------------------------------------
    type(cpl_fme_stream_t), intent(inout) :: st
    type(ESMF_Clock),       intent(in)    :: EClock
    integer,                intent(out)   :: ierr

    integer :: i, s, j, kf, n_my, l
    integer, allocatable :: gcol_to_rank(:), send_gcol_list(:)
    integer, allocatable :: gcol_to_myidx(:)
    integer, pointer     :: dof(:) => null()
    integer(PIO_OFFSET_KIND), allocatable :: idof(:)
    integer :: gr, ilon, ilat, nfl
    character(len=16) :: flist(max_fme_vars)
    real(r8) :: curr_time
    character(len=*), parameter :: subname = '(cpl_fme_stream_setup)'

    ierr = 0

    ! --- bind to the live source bundle + its decomposition -----------------
    call cpl_fme_resolve_bundle(st, ierr)
    if (ierr /= 0) then
       if (cpl_iam == 0) write(logunit,*) subname, ': bundle ', trim(st%bundle), &
            ' unavailable for ', trim(st%name), ' (ierr=', ierr, ')'
       return
    end if

    if (len_trim(st%mapfile) == 0) then
       if (cpl_iam == 0) write(logunit,*) subname, ': no map file for ', &
            trim(st%name), ' -- disabling'
       ierr = 7
       return
    end if

    ! --- load map (shared module) -------------------------------------------
    call st%map%read_mapfile(trim(st%mapfile), cpl_mpicom, cpl_iam, &
         cpl_npes, cpl_iosys, ierr)
    if (ierr /= 0) return

    ! The masked renormalization path is only safe for non-negative maps; a
    ! negative-weight map (TempestRemap 2nd-order FV, patch, ...) silently
    ! corrupts partially-covered cells (gotcha #47).  Refuse it loudly.
    if (st%map%has_negative_weights) then
       if (cpl_iam == 0) write(logunit,*) subname, ': ERROR map ', &
            trim(st%mapfile), ' has NEGATIVE weights (min=', st%map%min_weight, &
            '); unsafe for masked FME remap. Use a non-negative map.'
       ierr = 1
       return
    end if

    if (st%map%n_a /= st%src_gsize) then
       if (cpl_iam == 0) write(logunit,*) subname, ': ERROR map n_a=', &
            st%map%n_a, ' /= source gsize=', st%src_gsize, ' for ', trim(st%name)
       ierr = 2
       return
    end if

    ! --- build comm pattern: gcol_to_rank from the source gsmap segments -----
    allocate(gcol_to_rank(st%map%n_a))
    gcol_to_rank(:) = -1
    do s = 1, st%src_gsmap%ngseg
       do j = 0, st%src_gsmap%length(s) - 1
          gcol_to_rank(st%src_gsmap%start(s) + j) = st%src_gsmap%pe_loc(s)
       end do
    end do
    call st%map%build_comm(gcol_to_rank, cpl_mpicom, cpl_iam, cpl_npes, &
         send_gcol_list, ierr)
    deallocate(gcol_to_rank)
    if (ierr /= 0) return

    ! --- send-side packing: global send gcol -> local aVect column ----------
    ! The aVect storage order matches mct_gsMap_OrderedPoints (the same
    ! assumption seq_io makes when it writes decomposed bundles).
    call mct_gsMap_OrderedPoints(st%src_gsmap, cpl_iam, dof)
    n_my = size(dof)
    allocate(gcol_to_myidx(st%map%n_a))
    gcol_to_myidx(:) = 0
    do l = 1, n_my
       gcol_to_myidx(dof(l)) = l
    end do
    st%n_send = st%map%n_send_total
    allocate(st%send_local_idx(max(1, st%n_send)))
    do i = 1, st%n_send
       st%send_local_idx(i) = gcol_to_myidx(send_gcol_list(i))
       if (st%send_local_idx(i) == 0) then
          if (cpl_iam == 0) write(logunit,*) subname, &
               ': ERROR send gcol not owned by this rank (map/mesh mismatch)'
          ierr = 3
          deallocate(gcol_to_myidx, dof, send_gcol_list)
          return
       end if
    end do
    deallocate(gcol_to_myidx, dof, send_gcol_list)
    allocate(st%send_buf(max(1, st%n_send * 2)))   ! numlev=1 -> (val,mask)

    ! --- resolve field list against the live source bundle ------------------
    call cpl_fme_field_list(st%bundle, flist, nfl)
    st%nflds = 0
    do i = 1, nfl
       kf = mct_aVect_indexRA(st%src_bundle, trim(flist(i)), perrWith='quiet')
       if (kf > 0) then
          st%nflds = st%nflds + 1
          st%fldnames(st%nflds) = trim(flist(i))
          st%fld_kidx(st%nflds) = kf
          ! MCT field-naming convention: 'S*' = state (instantaneous snapshot),
          ! 'F*' = flux (time-mean over the output window).  Everything else
          ! defaults to mean (treated as a flux) and is warned below.
          st%fld_is_state(st%nflds) = (flist(i)(1:1) == 'S')
          if (cpl_iam == 0 .and. flist(i)(1:1) /= 'S' .and. flist(i)(1:1) /= 'F') then
             write(logunit,*) subname, ': WARNING field ', trim(flist(i)), &
                  ' has no S/F prefix -- treating as flux (time: mean)'
          end if
       else if (cpl_iam == 0) then
          write(logunit,*) subname, ': note field ', trim(flist(i)), &
               ' absent from ', trim(st%bundle), ' -- skipped'
       end if
    end do
    if (st%nflds == 0) then
       ierr = 4
       return
    end if

    ! --- source-space accumulator aVect (one slot per selected field) -------
    call cpl_fme_build_accum_avect(st)

    ! --- target-grid PIO decomposition (n_b_local, 2D lat-lon, PIO_REAL) ----
    allocate(idof(max(1, st%map%n_b_local)))
    do i = 1, st%map%n_b_local
       gr   = st%map%row_start + i - 1
       ilon = mod(gr - 1, st%map%nlon) + 1
       ilat = (gr - 1) / st%map%nlon + 1
       idof(i) = int(ilon + st%map%nlon * (ilat - 1), PIO_OFFSET_KIND)
    end do
    call pio_initdecomp(cpl_iosys, PIO_REAL, &
         (/st%map%nlon, st%map%nlat/), idof, st%iodesc_dst)
    deallocate(idof)
    st%iodesc_valid = .true.

    ! --- drift-free output schedule (gotcha #11) ----------------------------
    call seq_timemgr_EClockGetData(EClock, curr_time=curr_time)
    st%last_output_time = curr_time
    st%next_output_time = curr_time + st%interval
    st%nAccum = 0

  end subroutine cpl_fme_stream_setup

  !=============================================================================
  subroutine cpl_fme_build_accum_avect(st)
    ! Build the source-space accumulator aVect holding exactly the selected
    ! fields, sized to the local ocean decomposition, zeroed.
    type(cpl_fme_stream_t), intent(inout) :: st
    character(len=CXX) :: rlist
    integer :: i

    rlist = trim(st%fldnames(1))
    do i = 2, st%nflds
       rlist = trim(rlist) // ':' // trim(st%fldnames(i))
    end do
    call mct_aVect_init(st%accum, rList=trim(rlist), lsize=st%src_lsize)
    call mct_aVect_zero(st%accum)
  end subroutine cpl_fme_build_accum_avect

  !=============================================================================
  subroutine cpl_fme_resolve_bundle(st, ierr)
    ! Bind a stream to its live source bundle + decomposition.
    !   'x2o' : ocean window-mean (x2oacc, in-place-averaged after
    !           prep_ocn_accum_avg -- PLAN key enabler #1), ocn gsmap.
    !   'xao' : atm-ocean bulk fluxes (Faox_*), ocn gsmap (no native
    !           accumulator -> we average snapshots ourselves).
    !   'x2i' : sea-ice merged import, ice gsmap (shares the ocean mesh; no
    !           native accumulator -> we average snapshots ourselves).
    type(cpl_fme_stream_t), intent(inout) :: st
    integer,                intent(out)   :: ierr
    ierr = 0
    select case (trim(st%bundle))
    case ('x2o')
       st%src_bundle => b_x2o ; st%src_gsmap => ocn_gsmap
       st%src_lsize = lsize_ocn ; st%src_gsize = gsize_ocn
    case ('xao')
       st%src_bundle => b_xao ; st%src_gsmap => ocn_gsmap
       st%src_lsize = lsize_ocn ; st%src_gsize = gsize_ocn
    case ('x2i')
       st%src_bundle => b_x2i ; st%src_gsmap => ice_gsmap
       st%src_lsize = lsize_ice ; st%src_gsize = gsize_ice
    case ('o2x')
       st%src_bundle => b_o2x ; st%src_gsmap => ocn_gsmap
       st%src_lsize = lsize_ocn ; st%src_gsize = gsize_ocn
    case ('i2x')
       st%src_bundle => b_i2x ; st%src_gsmap => ice_gsmap
       st%src_lsize = lsize_ice ; st%src_gsize = gsize_ice
    case ('x2a')
       st%src_bundle => b_x2a ; st%src_gsmap => atm_gsmap
       st%src_lsize = lsize_atm ; st%src_gsize = gsize_atm
    case ('a2x')
       st%src_bundle => b_a2x ; st%src_gsmap => atm_gsmap
       st%src_lsize = lsize_atm ; st%src_gsize = gsize_atm
    case default
       ierr = 5
       return
    end select
    if (.not. associated(st%src_bundle) .or. .not. associated(st%src_gsmap)) &
         ierr = 6   ! bundle/decomp not available in this configuration
  end subroutine cpl_fme_resolve_bundle

  !=============================================================================
  subroutine cpl_fme_field_list(bundle, names, n)
    ! Return the authoritative default field list for a bundle.
    character(len=*),  intent(in)  :: bundle
    character(len=16), intent(out) :: names(:)
    integer,           intent(out) :: n
    select case (trim(bundle))
    case ('x2o') ; n = n_x2o_default ; names(1:n) = x2o_default
    case ('xao') ; n = n_xao_default ; names(1:n) = xao_default
    case ('x2i') ; n = n_x2i_default ; names(1:n) = x2i_default
    case ('o2x') ; n = n_o2x_default ; names(1:n) = o2x_default
    case ('i2x') ; n = n_i2x_default ; names(1:n) = i2x_default
    case ('x2a') ; n = n_x2a_default ; names(1:n) = x2a_default
    case ('a2x') ; n = n_a2x_default ; names(1:n) = a2x_default
    case default ; n = 0
    end select
  end subroutine cpl_fme_field_list

  !=============================================================================
  subroutine cpl_fme_accum(EClock, phase)
    !---------------------------------------------------------------------------
    ! Sample the merged forcing into each enabled stream's source-space
    ! accumulator, then flush+remap+write any stream whose output window has
    ! elapsed.  Called once per coupling step at EACH of three phases
    ! (gotcha #51): CPL_FME_PHASE_OCN right after prep_ocn_accum_avg (x2o
    ! window-mean + xao/o2x/i2x/a2x), CPL_FME_PHASE_ICE right after
    ! prep_ice_mrg (x2i, freshly merged), CPL_FME_PHASE_ATM right after
    ! prep_atm_mrg (x2a, freshly merged).  Each stream is processed at exactly
    ! one phase (streams(s)%phase) so it samples its bundle when validly
    ! populated -- this is what makes x2i/x2a warm-restart BFB.  All three
    ! phases see the same curr_time (the driver clock advances once per step,
    ! before any prep), so the drift-free flush schedule stays consistent.
    !---------------------------------------------------------------------------
    type(ESMF_Clock), intent(in) :: EClock
    integer,          intent(in) :: phase

    integer :: s, kf, f
    real(r8) :: curr_time
    character(len=*), parameter :: subname = '(cpl_fme_accum)'

    if (.not. active) return
    if (.not. iamin_cpl) return

    call seq_timemgr_EClockGetData(EClock, curr_time=curr_time)

    do s = 1, nstreams
       if (.not. streams(s)%enabled) cycle
       if (streams(s)%phase /= phase) cycle

       ! Sample the latest values (source space) for the selected fields.
       ! State ('S*') fields are captured as an instantaneous snapshot -- the
       ! accumulator slot just holds the most recent sample (overwrite).  Flux
       ! ('F*') fields are summed and divided by nAccum at flush to form the
       ! time-mean over the output window (matching the coupler's own
       ! accumulate-then-average treatment of fluxes).
       do f = 1, streams(s)%nflds
          kf = streams(s)%fld_kidx(f)
          if (streams(s)%fld_is_state(f)) then
             streams(s)%accum%rAttr(f, :) = streams(s)%src_bundle%rAttr(kf, :)
          else
             streams(s)%accum%rAttr(f, :) = streams(s)%accum%rAttr(f, :) &
                  + streams(s)%src_bundle%rAttr(kf, :)
          end if
       end do
       streams(s)%nAccum = streams(s)%nAccum + 1

       ! flush when the (drift-free) window boundary is reached
       if (curr_time >= streams(s)%next_output_time - SCHED_EPS) then
          call cpl_fme_flush(streams(s), EClock, curr_time)
       end if
    end do

  end subroutine cpl_fme_accum

  !=============================================================================
  subroutine cpl_fme_flush(st, EClock, curr_time)
    !---------------------------------------------------------------------------
    ! Average the source-space accumulator, remap each field to the lat-lon
    ! grid, write a record to the monthly file (rotating at month boundary),
    ! then reset the accumulator and advance the schedule by exactly one
    ! interval (gotcha #11).
    !---------------------------------------------------------------------------
    type(cpl_fme_stream_t), intent(inout) :: st
    type(ESMF_Clock),       intent(in)    :: EClock
    real(r8),               intent(in)    :: curr_time

    integer :: f, i, ierr, curr_ymd, curr_tod
    real(r8) :: rinv, fred, tlo, thi
    real(r8), allocatable :: fld_out(:,:)
    character(len=*), parameter :: subname = '(cpl_fme_flush)'

    if (st%nAccum <= 0) return

    call seq_timemgr_EClockGetData(EClock, curr_ymd=curr_ymd, curr_tod=curr_tod)

    rinv = 1.0_r8 / real(st%nAccum, r8)
    tlo  = st%last_output_time
    thi  = curr_time

    ! rotate to the correct monthly file and (re)define vars as needed
    call cpl_fme_file_rotate(st, curr_ymd, ierr)

    ! write the time/bnds/date record for this frame
    call cpl_fme_file_write_time(st, tlo, thi, curr_ymd, curr_tod)

    allocate(fld_out(max(1, st%map%n_b_local), 1))
    do f = 1, st%nflds
       ! Reduce in source space, pack (value, mask=1) and remap (apply_masked
       ! -> land/no-coverage target cells get SHR_FILL_VALUE; gotchas #3/#5).
       ! State fields are the stored snapshot (no division); flux fields are
       ! divided by nAccum to form the window-mean.
       if (st%fld_is_state(f)) then
          fred = 1.0_r8        ! snapshot: accumulator already holds the value
       else
          fred = rinv          ! flux: divide the running sum by nAccum
       end if
       do i = 1, st%n_send
          st%send_buf((i-1)*2 + 1) = st%accum%rAttr(f, st%send_local_idx(i)) * fred
          st%send_buf((i-1)*2 + 2) = 1.0_r8
       end do
       call st%map%apply_masked(st%send_buf, 1, fld_out, cpl_mpicom, &
            cpl_npes, ierr)
       call cpl_fme_file_write_var(st, f, fld_out(:,1))
    end do
    deallocate(fld_out)

    ! reset accumulator + advance schedule by exactly one interval (no drift)
    call mct_aVect_zero(st%accum)
    st%nAccum = 0
    st%last_output_time = st%next_output_time
    st%next_output_time = st%next_output_time + st%interval

  end subroutine cpl_fme_flush

  !=============================================================================
  subroutine cpl_fme_file_rotate(st, curr_ymd, ierr)
    !---------------------------------------------------------------------------
    ! Open the monthly file matching curr_ymd, rotating (closing + reopening)
    ! at a calendar-month boundary.  Honors append-mode reopen on the first
    ! open of a warm restart (gotcha #29/#44) unless a rewind happened
    ! (gotcha #45), in which case it clobber-creates.
    !---------------------------------------------------------------------------
    type(cpl_fme_stream_t), intent(inout) :: st
    integer,                intent(in)    :: curr_ymd
    integer,                intent(out)   :: ierr

    character(len=CL) :: fname
    integer :: yy, mm, dd

    ierr = 0
    call shr_cal_date2ymd(curr_ymd, yy, mm, dd)
    write(fname, '(a,a,a,a,i4.4,a,i2.2,a)') trim(case_name), &
         '.cpl.fme.', trim(st%name), '.', yy, '-', mm, '.nc'

    if (st%outfile%is_open .and. &
         trim(fname) == trim(st%outfile%current_filename)) return

    if (st%outfile%is_open) call cpl_fme_file_close(st%outfile)

    call cpl_fme_file_open(st, trim(fname))

  end subroutine cpl_fme_file_rotate

  !=============================================================================
  subroutine cpl_fme_file_open(st, fname)
    !---------------------------------------------------------------------------
    ! Open (append) or create the monthly lat-lon file, define dims/coords/
    ! time vars + each data var, and seed the time_record on append.
    !---------------------------------------------------------------------------
    type(cpl_fme_stream_t), intent(inout), target :: st
    character(len=*),       intent(in)    :: fname

    type(cpl_fme_file_t), pointer :: fh
    logical :: file_exists, may_append
    integer :: rc, f, nt, it
    real(r8), allocatable :: existing_time(:)
    real(r8) :: now_days
    character(len=*), parameter :: subname = '(cpl_fme_file_open)'

    fh => st%outfile
    inquire(file=trim(fname), exist=file_exists)

    ! append only on the first open after a (non-rewind) warm restart
    may_append = file_exists .and. is_restart_run .and. &
         st%reopen_pending .and. .not. rewind_restart
    st%reopen_pending = .false.

    fh%current_filename = trim(fname)
    fh%nvars = 0

    if (may_append) then
       fh%reopened = .true.
       rc = pio_openfile(cpl_iosys, fh%pio_file, cpl_iotype, trim(fname), PIO_WRITE)
       call pio_chk(rc, subname//': openfile '//trim(fname))
       rc = pio_inq_dimid(fh%pio_file, 'lon',  fh%dim_lon)
       rc = pio_inq_dimid(fh%pio_file, 'lat',  fh%dim_lat)
       rc = pio_inq_dimid(fh%pio_file, 'time', fh%dim_time)
       rc = pio_inq_dimid(fh%pio_file, 'nbnd', fh%dim_nbnd)
       rc = pio_inq_varid(fh%pio_file, 'time',      fh%var_time)
       rc = pio_inq_varid(fh%pio_file, 'time_bnds', fh%var_tbnds)
       rc = pio_inq_varid(fh%pio_file, 'date',      fh%var_date)
       rc = pio_inq_varid(fh%pio_file, 'datesec',   fh%var_datesec)
       do f = 1, st%nflds
          rc = pio_inq_varid(fh%pio_file, trim(st%fldnames(f)), fh%var_descs(f))
       end do
       fh%nvars = st%nflds
       ! frame seeding (gotcha #29 pt 2): overwrite the first existing record
       ! whose time >= now (the leg-1 leftover straddling the restart); use
       ! STRICT-less-than count so that record is overwritten, not appended.
       rc = pio_inq_dimlen(fh%pio_file, fh%dim_time, nt)
       fh%time_record = nt
       if (nt > 0) then
          allocate(existing_time(nt))
          rc = pio_get_var(fh%pio_file, fh%var_time, existing_time)
          now_days = st%next_output_time
          fh%time_record = 0
          do it = 1, nt
             if (existing_time(it) < now_days - SCHED_EPS) &
                  fh%time_record = fh%time_record + 1
          end do
          deallocate(existing_time)
       end if
    else
       fh%reopened = .false.
       rc = pio_createfile(cpl_iosys, fh%pio_file, cpl_iotype, trim(fname), PIO_CLOBBER)
       call pio_chk(rc, subname//': createfile '//trim(fname))
       call cpl_fme_file_define(st)
       fh%time_record = 0
    end if

    fh%is_open = .true.

  end subroutine cpl_fme_file_open

  !=============================================================================
  subroutine cpl_fme_file_define(st)
    ! Define dims, lon/lat coords, CF time vars, global attrs and each data
    ! var on a freshly created file, then leave define mode.
    type(cpl_fme_stream_t), intent(inout), target :: st
    type(cpl_fme_file_t), pointer :: fh
    integer :: rc, f
    real(r8), parameter :: fillv = SHR_FILL_VALUE
    character(len=CL) :: lname, sname, cunit

    fh => st%outfile

    rc = pio_def_dim(fh%pio_file, 'lon',  st%map%nlon,   fh%dim_lon)
    rc = pio_def_dim(fh%pio_file, 'lat',  st%map%nlat,   fh%dim_lat)
    rc = pio_def_dim(fh%pio_file, 'time', PIO_UNLIMITED, fh%dim_time)
    rc = pio_def_dim(fh%pio_file, 'nbnd', 2,             fh%dim_nbnd)

    rc = pio_def_var(fh%pio_file, 'lon', PIO_DOUBLE, (/fh%dim_lon/), fh%var_lon)
    rc = pio_put_att(fh%pio_file, fh%var_lon, 'units', 'degrees_east')
    rc = pio_put_att(fh%pio_file, fh%var_lon, 'long_name', 'longitude')
    rc = pio_def_var(fh%pio_file, 'lat', PIO_DOUBLE, (/fh%dim_lat/), fh%var_lat)
    rc = pio_put_att(fh%pio_file, fh%var_lat, 'units', 'degrees_north')
    rc = pio_put_att(fh%pio_file, fh%var_lat, 'long_name', 'latitude')

    rc = pio_def_var(fh%pio_file, 'time', PIO_DOUBLE, (/fh%dim_time/), fh%var_time)
    rc = pio_put_att(fh%pio_file, fh%var_time, 'units', trim(time_units))
    rc = pio_put_att(fh%pio_file, fh%var_time, 'calendar', trim(calendar))
    rc = pio_put_att(fh%pio_file, fh%var_time, 'bounds', 'time_bnds')
    rc = pio_put_att(fh%pio_file, fh%var_time, 'axis', 'T')
    rc = pio_def_var(fh%pio_file, 'time_bnds', PIO_DOUBLE, &
         (/fh%dim_nbnd, fh%dim_time/), fh%var_tbnds)
    rc = pio_def_var(fh%pio_file, 'date', PIO_INT, (/fh%dim_time/), fh%var_date)
    rc = pio_put_att(fh%pio_file, fh%var_date, 'long_name', 'current date (YYYYMMDD)')
    rc = pio_def_var(fh%pio_file, 'datesec', PIO_INT, (/fh%dim_time/), fh%var_datesec)
    rc = pio_put_att(fh%pio_file, fh%var_datesec, 'long_name', 'current seconds of current date')

    ! global provenance attrs (gotcha #40)
    rc = pio_put_att(fh%pio_file, PIO_GLOBAL, 'Conventions', 'CF-1.8')
    rc = pio_put_att(fh%pio_file, PIO_GLOBAL, 'source', 'E3SM coupler-native FME output')
    rc = pio_put_att(fh%pio_file, PIO_GLOBAL, 'realm', 'cpl')
    rc = pio_put_att(fh%pio_file, PIO_GLOBAL, 'title', 'As-exchanged merged forcing ('//trim(st%bundle)//')')
    rc = pio_put_att(fh%pio_file, PIO_GLOBAL, 'time_reference_date', trim(seq_io_date2yyyymmdd(epoch_ymd)))
    rc = pio_put_att(fh%pio_file, PIO_GLOBAL, 'case', trim(case_name))
    rc = pio_put_att(fh%pio_file, PIO_GLOBAL, 'git_version', trim(model_version))
    ! document the per-field time-reduction convention so the files are
    ! self-describing (see also each var's cell_methods attribute)
    rc = pio_put_att(fh%pio_file, PIO_GLOBAL, 'time_reduction', &
         'State fields (name prefix "S") are instantaneous snapshots at the '// &
         'record time (cell_methods="time: point"); flux fields (name prefix '// &
         '"F") are means over the output window [time_bnds(1), time_bnds(2)] '// &
         '(cell_methods="time: mean").')

    do f = 1, st%nflds
       ! pull long_name / standard_name / units from the coupler field registry
       ! (same lookup the cpl history aux files use; falls back to the
       ! underscore-stripped shortname, then to 'unknown' if unregistered)
       call seq_flds_lookup(trim(st%fldnames(f)), longname=lname, &
            stdname=sname, units=cunit)
       rc = pio_def_var(fh%pio_file, trim(st%fldnames(f)), PIO_REAL, &
            (/fh%dim_lon, fh%dim_lat, fh%dim_time/), fh%var_descs(f))
       rc = pio_put_att(fh%pio_file, fh%var_descs(f), 'long_name', trim(lname))
       rc = pio_put_att(fh%pio_file, fh%var_descs(f), 'standard_name', trim(sname))
       rc = pio_put_att(fh%pio_file, fh%var_descs(f), 'units', trim(cunit))
       ! state ('S*') vars are instantaneous snapshots at the record time;
       ! flux ('F*') vars are means over [time_bnds(1), time_bnds(2)]
       if (st%fld_is_state(f)) then
          rc = pio_put_att(fh%pio_file, fh%var_descs(f), 'cell_methods', 'time: point')
       else
          rc = pio_put_att(fh%pio_file, fh%var_descs(f), 'cell_methods', 'time: mean')
       end if
       rc = pio_put_att(fh%pio_file, fh%var_descs(f), '_FillValue', real(fillv, r4))
       rc = pio_put_att(fh%pio_file, fh%var_descs(f), 'missing_value', real(fillv, r4))
    end do
    fh%nvars = st%nflds

    rc = pio_enddef(fh%pio_file)

    ! seed the static lon/lat coordinate variables once
    rc = pio_put_var(fh%pio_file, fh%var_lon, st%map%lon)
    rc = pio_put_var(fh%pio_file, fh%var_lat, st%map%lat)

  end subroutine cpl_fme_file_define

  !=============================================================================
  subroutine cpl_fme_file_write_time(st, tlo, thi, curr_ymd, curr_tod)
    ! Write one CF time record: numeric time = END of window (matches EAM,
    ! gotcha #39), time_bnds = [tlo, thi], date/datesec from the clock.
    type(cpl_fme_stream_t), intent(inout), target :: st
    real(r8),               intent(in)    :: tlo, thi
    integer,                intent(in)    :: curr_ymd, curr_tod
    type(cpl_fme_file_t), pointer :: fh
    integer :: rc
    real(r8) :: tbuf(1), bnds(2)
    integer  :: ibuf(1)

    fh => st%outfile
    fh%time_record = fh%time_record + 1

    tbuf(1) = thi
    rc = pio_put_var(fh%pio_file, fh%var_time, (/fh%time_record/), (/1/), tbuf)
    bnds(1) = tlo ; bnds(2) = thi
    rc = pio_put_var(fh%pio_file, fh%var_tbnds, (/1, fh%time_record/), (/2, 1/), bnds)
    ibuf(1) = curr_ymd
    rc = pio_put_var(fh%pio_file, fh%var_date, (/fh%time_record/), (/1/), ibuf)
    ibuf(1) = curr_tod
    rc = pio_put_var(fh%pio_file, fh%var_datesec, (/fh%time_record/), (/1/), ibuf)

  end subroutine cpl_fme_file_write_time

  !=============================================================================
  subroutine cpl_fme_file_write_var(st, f, fld)
    ! Write one remapped 2D field to the current frame.  PIO_REAL decomp =>
    ! convert to real(4) before writing (gotcha #1).
    type(cpl_fme_stream_t), intent(inout), target :: st
    integer,                intent(in)    :: f
    real(r8),               intent(in)    :: fld(:)
    type(cpl_fme_file_t), pointer :: fh
    integer :: rc
    real(r4), allocatable :: fld_r4(:)

    fh => st%outfile
    allocate(fld_r4(max(1, st%map%n_b_local)))
    fld_r4(1:st%map%n_b_local) = real(fld(1:st%map%n_b_local), r4)
    call pio_setframe(fh%pio_file, fh%var_descs(f), int(fh%time_record, PIO_OFFSET_KIND))
    call pio_write_darray(fh%pio_file, fh%var_descs(f), st%iodesc_dst, fld_r4, rc)
    deallocate(fld_r4)

  end subroutine cpl_fme_file_write_var

  !=============================================================================
  subroutine cpl_fme_file_close(fh)
    type(cpl_fme_file_t), intent(inout) :: fh
    if (fh%is_open) then
       call pio_closefile(fh%pio_file)
       fh%is_open = .false.
       fh%current_filename = ''
    end if
  end subroutine cpl_fme_file_close

  !=============================================================================
  subroutine cpl_fme_restart_write(EClock)
    !---------------------------------------------------------------------------
    ! Persist each enabled stream's source-space accumulator + sample count +
    ! drift-free schedule to a sidecar, on the coupler restart alarm.  The
    ! base x2oacc is already cpl-restart-safe (we do NOT re-accumulate it);
    ! only this 1D/5D super-accumulator layer is ours to persist.
    !---------------------------------------------------------------------------
    type(ESMF_Clock), intent(in) :: EClock
    integer :: s
    if (.not. active) return
    if (.not. iamin_cpl) return
    do s = 1, nstreams
       if (streams(s)%enabled) call cpl_fme_write_sidecar(streams(s))
    end do
  end subroutine cpl_fme_restart_write

  !=============================================================================
  subroutine cpl_fme_write_sidecar(st)
    ! Write <case>.cpl.fme.<stream>.accum.nc: per-field source-space sums on
    ! the ocean decomposition + nAccum + schedule + epoch.  Direct PIO with a
    ! 1D ocean-grid decomposition (built here, freed at close).
    type(cpl_fme_stream_t), intent(inout) :: st
    type(file_desc_t) :: pf
    type(io_desc_t)   :: iod
    type(var_desc_t)  :: vd(max_fme_vars)
    integer :: rc, f, dim_ncol
    integer, pointer  :: dofs(:) => null()
    real(r8), allocatable :: slab(:)
    character(len=CL) :: fname
    character(len=*), parameter :: subname = '(cpl_fme_write_sidecar)'

    write(fname, '(a,a,a,a)') trim(case_name), '.cpl.fme.', trim(st%name), '.accum.nc'

    rc = pio_createfile(cpl_iosys, pf, cpl_iotype, trim(fname), PIO_CLOBBER)
    call pio_chk(rc, subname//': create '//trim(fname))
    rc = pio_def_dim(pf, 'ncol', st%src_gsize, dim_ncol)
    do f = 1, st%nflds
       rc = pio_def_var(pf, 'accum_'//trim(st%fldnames(f)), PIO_DOUBLE, &
            (/dim_ncol/), vd(f))
    end do
    rc = pio_put_att(pf, PIO_GLOBAL, 'source', 'E3SM coupler FME accumulator sidecar')
    rc = pio_put_att(pf, PIO_GLOBAL, 'nAccum', st%nAccum)
    rc = pio_put_att(pf, PIO_GLOBAL, 'last_output_time', st%last_output_time)
    rc = pio_put_att(pf, PIO_GLOBAL, 'next_output_time', st%next_output_time)
    rc = pio_put_att(pf, PIO_GLOBAL, 'time_reference_date', trim(seq_io_date2yyyymmdd(epoch_ymd)))
    rc = pio_enddef(pf)

    call mct_gsMap_OrderedPoints(st%src_gsmap, cpl_iam, dofs)
    call pio_initdecomp(cpl_iosys, PIO_DOUBLE, (/st%src_gsize/), dofs, iod)
    allocate(slab(size(dofs)))
    do f = 1, st%nflds
       slab(:) = st%accum%rAttr(f, :)
       call pio_write_darray(pf, vd(f), iod, slab, rc)
    end do
    deallocate(slab, dofs)
    call pio_freedecomp(pf, iod)
    call pio_closefile(pf)

  end subroutine cpl_fme_write_sidecar

  !=============================================================================
  subroutine cpl_fme_read_sidecar(st, EClock)
    ! Restore accumulator + schedule from the sidecar on a warm restart.
    ! Detects rewind/branch restart (restored schedule in the FUTURE relative
    ! to MPAS_NOW): reset schedule to now, zero the accumulator, and latch
    ! rewind_restart so monthly files clobber-create (gotcha #45).
    type(cpl_fme_stream_t), intent(inout) :: st
    type(ESMF_Clock),       intent(in)    :: EClock
    type(file_desc_t) :: pf
    type(io_desc_t)   :: iod
    type(var_desc_t)  :: vd
    integer :: rc, f, na
    integer, pointer  :: dofs(:) => null()
    real(r8), allocatable :: slab(:)
    real(r8) :: lot, not_, curr_time
    logical :: exists
    character(len=CL) :: fname
    character(len=*), parameter :: subname = '(cpl_fme_read_sidecar)'

    write(fname, '(a,a,a,a)') trim(case_name), '.cpl.fme.', trim(st%name), '.accum.nc'
    inquire(file=trim(fname), exist=exists)
    if (.not. exists) then
       if (cpl_iam == 0) write(logunit,*) subname, ': WARN no sidecar ', &
            trim(fname), ' on warm restart -- cold-starting accumulator for ', &
            trim(st%name)
       return
    end if

    rc = pio_openfile(cpl_iosys, pf, cpl_iotype, trim(fname), PIO_NOWRITE)
    call pio_chk(rc, subname//': open '//trim(fname))
    rc = pio_get_att(pf, PIO_GLOBAL, 'nAccum', na)
    rc = pio_get_att(pf, PIO_GLOBAL, 'last_output_time', lot)
    rc = pio_get_att(pf, PIO_GLOBAL, 'next_output_time', not_)

    call mct_gsMap_OrderedPoints(st%src_gsmap, cpl_iam, dofs)
    call pio_initdecomp(cpl_iosys, PIO_DOUBLE, (/st%src_gsize/), dofs, iod)
    allocate(slab(size(dofs)))
    do f = 1, st%nflds
       rc = pio_inq_varid(pf, 'accum_'//trim(st%fldnames(f)), vd)
       call pio_read_darray(pf, vd, iod, slab, rc)
       st%accum%rAttr(f, :) = slab(:)
    end do
    deallocate(slab, dofs)
    call pio_freedecomp(pf, iod)
    call pio_closefile(pf)

    st%nAccum = na
    st%last_output_time = lot
    st%next_output_time = not_

    ! rewind detection (gotcha #45): strict '>' keeps an exact same-point warm
    ! restart byte-identical (BFB), but a rewind (schedule in the future) is
    ! reset and the re-crossed monthly files are clobber-created.
    call seq_timemgr_EClockGetData(EClock, curr_time=curr_time)
    if (st%last_output_time > curr_time + SCHED_EPS) then
       if (cpl_iam == 0) write(logunit,*) subname, &
            ': rewind restart detected for ', trim(st%name), &
            ' -- resetting schedule to now and clobber-creating files'
       st%last_output_time = curr_time
       st%next_output_time = curr_time + st%interval
       st%nAccum = 0
       call mct_aVect_zero(st%accum)
       rewind_restart = .true.
    end if

  end subroutine cpl_fme_read_sidecar

  !=============================================================================
  subroutine cpl_fme_final()
    integer :: s
    if (.not. iamin_cpl) return
    do s = 1, nstreams
       if (.not. streams(s)%enabled) cycle
       ! Close the file FIRST: pio_write_darray buffers field writes for the
       ! whole leg and flushes them inside pio_closefile, and that flush
       ! dereferences iodesc_dst.  Freeing the decomp before the close would
       ! leave the buffered flush pointing at a freed decomposition and PIO
       ! aborts (PIOc_write_darray_multi -> pio_err -> piodie).
       if (streams(s)%outfile%is_open) call cpl_fme_file_close(streams(s)%outfile)
       if (streams(s)%iodesc_valid) then
          call pio_freedecomp(cpl_iosys, streams(s)%iodesc_dst)
          streams(s)%iodesc_valid = .false.
       end if
    end do
  end subroutine cpl_fme_final

  !=============================================================================
  subroutine pio_chk(rc, msg)
    integer,          intent(in) :: rc
    character(len=*), intent(in) :: msg
    if (rc /= PIO_NOERR) then
       if (cpl_iam == 0) write(logunit,*) 'cpl_fme PIO error: ', trim(msg), ' rc=', rc
       call shr_sys_abort('cpl_fme PIO error: '//trim(msg))
    end if
  end subroutine pio_chk

end module cpl_fme_mod
