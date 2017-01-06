module cesm_comp_mod

!-------------------------------------------------------------------------------
!
! Purpose: Main program for NCAR CESM4/cpl7. Can have different
!          land, sea-ice, and ocean models plugged in at compile-time.
!          These models can be either: stub, dead, data, or active
!          components or some combination of the above.
!
!               stub -------- Do nothing.
!               dead -------- Send analytic data back.
!               data -------- Send data back interpolated from input files.
!               prognostic -- Prognostically simulate the given component.
!
! Method: Call appropriate initialization, run (time-stepping), and
!         finalization routines.
!
!-------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! share code & libs
   !----------------------------------------------------------------------------
   use shr_kind_mod,      only: r8 => SHR_KIND_R8
   use shr_kind_mod,      only: cs => SHR_KIND_CS
   use shr_kind_mod,      only: cl => SHR_KIND_CL
   use shr_sys_mod,       only: shr_sys_abort, shr_sys_flush
   use shr_const_mod,     only: shr_const_cday
   use shr_file_mod,      only: shr_file_setLogLevel, shr_file_setLogUnit
   use shr_file_mod,      only: shr_file_setIO, shr_file_getUnit
   use shr_scam_mod,      only: shr_scam_checkSurface
   use shr_map_mod,       only: shr_map_setDopole
   use shr_mpi_mod,       only: shr_mpi_min, shr_mpi_max
   use shr_mem_mod,       only: shr_mem_init, shr_mem_getusage
   use shr_cal_mod,       only: shr_cal_date2ymd, shr_cal_ymd2date, shr_cal_advdateInt
   use shr_orb_mod,       only: shr_orb_params
   use shr_frz_mod,       only: shr_frz_freezetemp_init
   use shr_reprosum_mod,  only: shr_reprosum_setopts
   use mct_mod            ! mct_ wrappers for mct lib
   use perf_mod
   use ESMF

   !----------------------------------------------------------------------------
   ! component model interfaces (init, run, final methods)
   !----------------------------------------------------------------------------

   use atm_comp_mct  , only: atm_init=>atm_init_mct, atm_run=>atm_run_mct, atm_final=>atm_final_mct
   use lnd_comp_mct  , only: lnd_init=>lnd_init_mct, lnd_run=>lnd_run_mct, lnd_final=>lnd_final_mct
   use ocn_comp_mct  , only: ocn_init=>ocn_init_mct, ocn_run=>ocn_run_mct, ocn_final=>ocn_final_mct
   use ice_comp_mct  , only: ice_init=>ice_init_mct, ice_run=>ice_run_mct, ice_final=>ice_final_mct
   use glc_comp_mct  , only: glc_init=>glc_init_mct, glc_run=>glc_run_mct, glc_final=>glc_final_mct
   use wav_comp_mct  , only: wav_init=>wav_init_mct, wav_run=>wav_run_mct, wav_final=>wav_final_mct
   use rof_comp_mct  , only: rof_init=>rof_init_mct, rof_run=>rof_run_mct, rof_final=>rof_final_mct
   use esp_comp_mct  , only: esp_init=>esp_init_mct, esp_run=>esp_run_mct, esp_final=>esp_final_mct

   !----------------------------------------------------------------------------
   ! cpl7 modules
   !----------------------------------------------------------------------------

   ! mpi comm data & routines, plus logunit and loglevel
    use seq_comm_mct, only: CPLID, GLOID, logunit, loglevel
    use seq_comm_mct, only: ATMID, LNDID, OCNID, ICEID, GLCID, ROFID, WAVID, ESPID
    use seq_comm_mct, only: ALLATMID,ALLLNDID,ALLOCNID,ALLICEID,ALLGLCID,ALLROFID,ALLWAVID,ALLESPID
    use seq_comm_mct, only: CPLALLATMID,CPLALLLNDID,CPLALLOCNID,CPLALLICEID
    use seq_comm_mct, only: CPLALLGLCID,CPLALLROFID,CPLALLWAVID,CPLALLESPID
    use seq_comm_mct, only: CPLATMID,CPLLNDID,CPLOCNID,CPLICEID,CPLGLCID,CPLROFID,CPLWAVID,CPLESPID
    use seq_comm_mct, only: num_inst_atm, num_inst_lnd, num_inst_rof
    use seq_comm_mct, only: num_inst_ocn, num_inst_ice, num_inst_glc
    use seq_comm_mct, only: num_inst_wav, num_inst_esp
    use seq_comm_mct, only: num_inst_xao, num_inst_frc, num_inst_phys
    use seq_comm_mct, only: num_inst_total, num_inst_max
    use seq_comm_mct, only: seq_comm_iamin, seq_comm_name, seq_comm_namelen
    use seq_comm_mct, only: seq_comm_init, seq_comm_setnthreads, seq_comm_getnthreads
    use seq_comm_mct, only: seq_comm_getinfo => seq_comm_setptrs
    use seq_comm_mct, only: seq_comm_petlist

   ! clock & alarm routines and variables
   use seq_timemgr_mod, only: seq_timemgr_type
   use seq_timemgr_mod, only: seq_timemgr_clockInit
   use seq_timemgr_mod, only: seq_timemgr_clockAdvance
   use seq_timemgr_mod, only: seq_timemgr_clockPrint
   use seq_timemgr_mod, only: seq_timemgr_EClockGetData
   use seq_timemgr_mod, only: seq_timemgr_alarmIsOn
   use seq_timemgr_mod, only: seq_timemgr_histavg_type
   use seq_timemgr_mod, only: seq_timemgr_type_never
   use seq_timemgr_mod, only: seq_timemgr_alarm_restart
   use seq_timemgr_mod, only: seq_timemgr_alarm_stop
   use seq_timemgr_mod, only: seq_timemgr_alarm_datestop
   use seq_timemgr_mod, only: seq_timemgr_alarm_history
   use seq_timemgr_mod, only: seq_timemgr_alarm_atmrun
   use seq_timemgr_mod, only: seq_timemgr_alarm_lndrun
   use seq_timemgr_mod, only: seq_timemgr_alarm_ocnrun
   use seq_timemgr_mod, only: seq_timemgr_alarm_icerun
   use seq_timemgr_mod, only: seq_timemgr_alarm_glcrun
   use seq_timemgr_mod, only: seq_timemgr_alarm_ocnnext
   use seq_timemgr_mod, only: seq_timemgr_alarm_tprof
   use seq_timemgr_mod, only: seq_timemgr_alarm_histavg
   use seq_timemgr_mod, only: seq_timemgr_alarm_rofrun
   use seq_timemgr_mod, only: seq_timemgr_alarm_wavrun
   use seq_timemgr_mod, only: seq_timemgr_alarm_esprun
   use seq_timemgr_mod, only: seq_timemgr_alarm_barrier

   ! "infodata" gathers various control flags into one datatype
   use seq_infodata_mod, only: seq_infodata_putData, seq_infodata_GetData
   use seq_infodata_mod, only: seq_infodata_init, seq_infodata_exchange
   use seq_infodata_mod, only: seq_infodata_type, seq_infodata_orb_variable_year
   use seq_infodata_mod, only: seq_infodata_print

   ! domain related routines
   use seq_domain_mct, only : seq_domain_check

   ! history file routines
   use seq_hist_mod, only : seq_hist_write, seq_hist_writeavg, seq_hist_writeaux

   ! restart file routines
   use seq_rest_mod, only : seq_rest_read, seq_rest_write

   ! flux calc routines
   use seq_flux_mct, only: seq_flux_init_mct, seq_flux_initexch_mct, seq_flux_ocnalb_mct
   use seq_flux_mct, only: seq_flux_atmocn_mct, seq_flux_atmocnexch_mct

   ! domain fraction routines
   use seq_frac_mct, only : seq_frac_init, seq_frac_set

   ! i/o subroutines
   use seq_io_mod, only : seq_io_cpl_init

   ! rearrange type routines
   use cplcomp_exchange_mod, only: seq_mctext_decomp

   ! diagnostic routines
   use seq_diag_mct, only : seq_diag_zero_mct , seq_diag_avect_mct, seq_diag_lnd_mct
   use seq_diag_mct, only : seq_diag_rof_mct  , seq_diag_ocn_mct  , seq_diag_atm_mct
   use seq_diag_mct, only : seq_diag_ice_mct  , seq_diag_accum_mct, seq_diag_print_mct

   ! list of fields transferred between components
   use seq_flds_mod, only : seq_flds_a2x_fluxes, seq_flds_x2a_fluxes
   use seq_flds_mod, only : seq_flds_i2x_fluxes, seq_flds_x2i_fluxes
   use seq_flds_mod, only : seq_flds_l2x_fluxes, seq_flds_x2l_fluxes
   use seq_flds_mod, only : seq_flds_o2x_fluxes, seq_flds_x2o_fluxes
   use seq_flds_mod, only : seq_flds_g2x_fluxes, seq_flds_x2g_fluxes
   use seq_flds_mod, only : seq_flds_w2x_fluxes, seq_flds_x2w_fluxes
   use seq_flds_mod, only : seq_flds_r2x_fluxes, seq_flds_x2r_fluxes
   use seq_flds_mod, only : seq_flds_set

   ! component type and accessor functions
   use component_type_mod , only: component_get_iamin_compid, component_get_suffix
   use component_type_mod , only: component_get_name, component_get_c2x_cx
   use component_type_mod , only: atm, lnd, ice, ocn, rof, glc, wav, esp
   use component_mod      , only: component_init_pre
   use component_mod      , only: component_init_cc, component_init_cx, component_run, component_final
   use component_mod      , only: component_init_areacor, component_init_aream
   use component_mod      , only: component_exch, component_diag

   ! prep routines (includes mapping routines between components and merging routines)
   use prep_lnd_mod
   use prep_ice_mod
   use prep_wav_mod
   use prep_rof_mod
   use prep_glc_mod
   use prep_ocn_mod
   use prep_atm_mod
   use prep_aoflux_mod

   !--- mapping routines ---
   use seq_map_type_mod
   use seq_map_mod      ! generic mapping

   ! --- timing routines ---
   use t_drv_timers_mod

   implicit none

   private

   public cesm_pre_init1, cesm_pre_init2, cesm_init, cesm_run, cesm_final
   public timing_dir, mpicom_GLOID

#include <mpif.h>

   !----------------------------------------------------------------------------
   ! temporary variables
   !----------------------------------------------------------------------------

   !- from prep routines (arrays of instances)
   type(mct_aVect) , pointer :: a2x_ox(:) => null()
   type(mct_aVect) , pointer :: o2x_ax(:) => null()
   type(mct_aVect) , pointer :: xao_ox(:) => null()
   type(mct_aVect) , pointer :: xao_ax(:) => null()

   !- from component type (single instance inside array of components)
   type(mct_aVect) , pointer :: o2x_ox => null()
   type(mct_aVect) , pointer :: a2x_ax => null()

   character(len=CL) :: suffix
   logical           :: iamin_id
   logical           :: iamroot_id
   integer           :: mpicom
   character(len=seq_comm_namelen) :: compname

   !----------------------------------------------------------------------------
   ! domains & related
   !----------------------------------------------------------------------------

   !--- domain fractions (only defined on cpl pes) ---
   type(mct_aVect) , pointer :: fractions_ax(:)   ! Fractions on atm grid, cpl processes
   type(mct_aVect) , pointer :: fractions_lx(:)   ! Fractions on lnd grid, cpl processes
   type(mct_aVect) , pointer :: fractions_ix(:)   ! Fractions on ice grid, cpl processes
   type(mct_aVect) , pointer :: fractions_ox(:)   ! Fractions on ocn grid, cpl processes
   type(mct_aVect) , pointer :: fractions_gx(:)   ! Fractions on glc grid, cpl processes
   type(mct_aVect) , pointer :: fractions_rx(:)   ! Fractions on rof grid, cpl processes
   type(mct_aVect) , pointer :: fractions_wx(:)   ! Fractions on wav grid, cpl processes

   !--- domain equivalent 2d grid size ---
   integer  :: atm_nx, atm_ny  ! nx, ny of 2d grid, if known
   integer  :: lnd_nx, lnd_ny
   integer  :: ice_nx, ice_ny
   integer  :: ocn_nx, ocn_ny
   integer  :: rof_nx, rof_ny
   integer  :: glc_nx, glc_ny
   integer  :: wav_nx, wav_ny

   !----------------------------------------------------------------------------
   ! Infodata: inter-model control flags, domain info
   !----------------------------------------------------------------------------

   type (seq_infodata_type), target :: infodata ! single instance for cpl and all comps

   !----------------------------------------------------------------------------
   ! time management
   !----------------------------------------------------------------------------

   type (seq_timemgr_type), SAVE :: seq_SyncClock ! array of all clocks & alarm
   type (ESMF_Clock), target :: EClock_d      ! driver clock
   type (ESMF_Clock), target :: EClock_a
   type (ESMF_Clock), target :: EClock_l
   type (ESMF_Clock), target :: EClock_o
   type (ESMF_Clock), target :: EClock_i
   type (ESMF_Clock), target :: EClock_g
   type (ESMF_Clock), target :: EClock_r
   type (ESMF_Clock), target :: EClock_w
   type (ESMF_Clock), target :: EClock_e

   logical  :: restart_alarm          ! restart alarm
   logical  :: history_alarm          ! history alarm
   logical  :: histavg_alarm          ! history alarm
   logical  :: stop_alarm             ! stop alarm
   logical  :: atmrun_alarm           ! atm run alarm
   logical  :: lndrun_alarm           ! lnd run alarm
   logical  :: icerun_alarm           ! ice run alarm
   logical  :: ocnrun_alarm           ! ocn run alarm
   logical  :: ocnnext_alarm          ! ocn run alarm on next timestep
   logical  :: glcrun_alarm           ! glc run alarm
   logical  :: rofrun_alarm           ! rof run alarm
   logical  :: wavrun_alarm           ! wav run alarm
   logical  :: esprun_alarm           ! esp run alarm
   logical  :: tprof_alarm            ! timing profile alarm
   logical  :: barrier_alarm          ! barrier alarm
   logical  :: t1hr_alarm             ! alarm every hour
   logical  :: t2hr_alarm             ! alarm every two hours
   logical  :: t3hr_alarm             ! alarm every three hours
   logical  :: t6hr_alarm             ! alarm every six hours
   logical  :: t12hr_alarm            ! alarm every twelve hours
   logical  :: t24hr_alarm            ! alarm every twentyfour hours
   logical  :: t1yr_alarm             ! alarm every year, at start of year

   real(r8) :: days_per_year = 365.0  ! days per year

   integer  :: dtime                  ! dt of one coupling interval
   integer  :: ncpl                   ! number of coupling intervals per day
   integer  :: ymd                    ! Current date (YYYYMMDD)
   integer  :: year                   ! Current date (YYYY)
   integer  :: month                  ! Current date (MM)
   integer  :: day                    ! Current date (DD)
   integer  :: tod                    ! Current time of day (seconds)
   integer  :: ymdtmp                 ! temporary date (YYYYMMDD)
   integer  :: todtmp                 ! temporary time of day (seconds)
   character(CL) :: orb_mode          ! orbital mode
   character(CS) :: tfreeze_option    ! Freezing point calculation
   integer  :: orb_iyear              ! orbital year
   integer  :: orb_iyear_align        ! associated with model year
   integer  :: orb_cyear              ! orbital year for current orbital computation
   integer  :: orb_nyear              ! orbital year associated with currrent model year
   real(r8) :: orb_eccen              ! orbital eccentricity
   real(r8) :: orb_obliq              ! obliquity in degrees
   real(r8) :: orb_mvelp              ! moving vernal equinox long
   real(r8) :: orb_obliqr             ! Earths obliquity in rad
   real(r8) :: orb_lambm0             ! Mean long of perihelion at vernal equinox (radians)
   real(r8) :: orb_mvelpp             ! moving vernal equinox long
   real(r8) :: wall_time_limit        ! wall time limit in hours
   real(r8) :: wall_time              ! current wall time used
   character(CS) :: force_stop_at     ! force stop at next (month, day, etc)
   logical  :: force_stop             ! force the model to stop
   integer  :: force_stop_ymd         ! force stop ymd
   integer  :: force_stop_tod         ! force stop tod

   !--- for documenting speed of the model ---
   character( 8) :: dstr              ! date string
   character(10) :: tstr              ! time string
   integer       :: begStep, endStep  ! Begining and ending step number
   character(CL) :: calendar          ! calendar name
   real(r8)      :: simDays           ! Number of simulated days
   real(r8)      :: SYPD              ! Simulated years per day
   real(r8)      :: Time_begin        ! Start time
   real(r8)      :: Time_end          ! Ending time
   real(r8)      :: Time_bstep        ! Start time
   real(r8)      :: Time_estep        ! Ending time
   real(r8)      :: time_brun         ! Start time
   real(r8)      :: time_erun         ! Ending time
   real(r8)      :: cktime            ! delta time
   real(r8)      :: cktime_acc(10)    ! cktime accumulator array 1 = all, 2 = atm, etc
   integer       :: cktime_cnt(10)    ! cktime counter array
   real(r8)      :: max_cplstep_time
   character(CL) :: timing_file       ! Local path to tprof filename
   character(CL) :: timing_dir        ! timing directory
   character(CL) :: tchkpt_dir        ! timing checkpoint directory

   !----------------------------------------------------------------------------
   ! control flags
   !----------------------------------------------------------------------------

   logical  :: atm_present            ! .true.  => atm is present
   logical  :: lnd_present            ! .true.  => land is present
   logical  :: ice_present            ! .true.  => ice is present
   logical  :: ocn_present            ! .true.  => ocn is present
   logical  :: glc_present            ! .true.  => glc is present
   logical  :: glclnd_present         ! .true.  => glc is computing land coupling
   logical  :: glcocn_present         ! .true.  => glc is computing ocean runoff
   logical  :: glcice_present         ! .true.  => glc is computing icebergs
   logical  :: rofice_present         ! .true.  => rof is computing icebergs
   logical  :: rof_present            ! .true.  => rof is present
   logical  :: flood_present          ! .true.  => rof is computing flood
   logical  :: wav_present            ! .true.  => wav is present
   logical  :: esp_present            ! .true.  => esp is present

   logical  :: atm_prognostic         ! .true.  => atm comp expects input
   logical  :: lnd_prognostic         ! .true.  => lnd comp expects input
   logical  :: ice_prognostic         ! .true.  => ice comp expects input
   logical  :: iceberg_prognostic     ! .true.  => ice comp can handle iceberg input
   logical  :: ocn_prognostic         ! .true.  => ocn comp expects input
   logical  :: ocnrof_prognostic      ! .true.  => ocn comp expects runoff input
   logical  :: glc_prognostic         ! .true.  => glc comp expects input
   logical  :: rof_prognostic         ! .true.  => rof comp expects input
   logical  :: wav_prognostic         ! .true.  => wav comp expects input
   logical  :: esp_prognostic         ! .true.  => esp comp expects input

   logical  :: atm_c2_lnd             ! .true.  => atm to lnd coupling on
   logical  :: atm_c2_ocn             ! .true.  => atm to ocn coupling on
   logical  :: atm_c2_ice             ! .true.  => atm to ice coupling on
   logical  :: atm_c2_wav             ! .true.  => atm to wav coupling on
   logical  :: lnd_c2_atm             ! .true.  => lnd to atm coupling on
   logical  :: lnd_c2_rof             ! .true.  => lnd to rof coupling on
   logical  :: lnd_c2_glc             ! .true.  => lnd to glc coupling on
   logical  :: ocn_c2_atm             ! .true.  => ocn to atm coupling on
   logical  :: ocn_c2_ice             ! .true.  => ocn to ice coupling on
   logical  :: ocn_c2_wav             ! .true.  => ocn to wav coupling on
   logical  :: ice_c2_atm             ! .true.  => ice to atm coupling on
   logical  :: ice_c2_ocn             ! .true.  => ice to ocn coupling on
   logical  :: ice_c2_wav             ! .true.  => ice to wav coupling on
   logical  :: rof_c2_lnd             ! .true.  => rof to lnd coupling on
   logical  :: rof_c2_ocn             ! .true.  => rof to ocn coupling on
   logical  :: rof_c2_ice             ! .true.  => rof to ice coupling on
   logical  :: glc_c2_lnd             ! .true.  => glc to lnd coupling on
   logical  :: glc_c2_ocn             ! .true.  => glc to ocn coupling on
   logical  :: glc_c2_ice             ! .true.  => glc to ice coupling on
   logical  :: wav_c2_ocn             ! .true.  => wav to ocn coupling on

   logical  :: dead_comps             ! .true.  => dead components
   logical  :: esmf_map_flag          ! .true.  => use esmf for mapping

   logical  :: areafact_samegrid      ! areafact samegrid flag
   logical  :: single_column          ! scm mode logical
   real(r8) :: scmlon                 ! single column lon
   real(r8) :: scmlat                 ! single column lat
   logical  :: aqua_planet            ! aqua planet mode
   real(r8) :: nextsw_cday            ! radiation control
   logical  :: atm_aero               ! atm provides aerosol data

   character(CL) :: cpl_seq_option    ! coupler sequencing option
   logical  :: skip_ocean_run         ! skip the ocean model first pass
   logical  :: cpl2ocn_first          ! use to call initial cpl2ocn timer
   logical  :: run_barriers           ! barrier the component run calls

   character(CS) :: aoflux_grid       ! grid for a/o flux calc: atm xor ocn
   character(CS) :: vect_map          ! vector mapping type

   character(CL) :: atm_gnam          ! atm grid
   character(CL) :: lnd_gnam          ! lnd grid
   character(CL) :: ocn_gnam          ! ocn grid
   character(CL) :: ice_gnam          ! ice grid
   character(CL) :: rof_gnam          ! rof grid
   character(CL) :: glc_gnam          ! glc grid
   character(CL) :: wav_gnam          ! wav grid

   logical  :: samegrid_ao            ! samegrid atm and ocean
   logical  :: samegrid_al            ! samegrid atm and land
   logical  :: samegrid_lr            ! samegrid land and rof
   logical  :: samegrid_oi            ! samegrid ocean and ice
   logical  :: samegrid_ro            ! samegrid runoff and ocean
   logical  :: samegrid_aw            ! samegrid atm and wave
   logical  :: samegrid_ow            ! samegrid ocean and wave
   logical  :: samegrid_lg            ! samegrid glc and land
   logical  :: samegrid_og            ! samegrid glc and ocean
   logical  :: samegrid_ig            ! samegrid glc and ice
   logical  :: samegrid_alo           ! samegrid atm, lnd, ocean

   logical       :: read_restart      ! local read restart flag
   character(CL) :: rest_file         ! restart file path + filename

   logical  :: shr_map_dopole         ! logical for dopole in shr_map_mod
   logical  :: domain_check           ! .true.  => check consistency of domains
   logical  :: reprosum_use_ddpdd     ! setup reprosum, use ddpdd
   real(r8) :: reprosum_diffmax       ! setup reprosum, set rel_diff_max
   logical  :: reprosum_recompute     ! setup reprosum, recompute if tolerance exceeded

   logical  :: output_perf = .false.  ! require timing data output for this pe

   !--- history & budgets ---
   logical :: do_budgets              ! heat/water budgets on
   logical :: do_histinit             ! initial hist file
   logical :: do_histavg              ! histavg on or off
   logical :: do_hist_r2x             ! create aux files: r2x
   logical :: do_hist_l2x             ! create aux files: l2x
   logical :: do_hist_a2x24hr         ! create aux files: a2x
   logical :: do_hist_l2x1yr          ! create aux files: l2x
   logical :: do_hist_a2x             ! create aux files: a2x
   logical :: do_hist_a2x3hrp         ! create aux files: a2x 3hr precip
   logical :: do_hist_a2x3hr          ! create aux files: a2x 3hr states
   logical :: do_hist_a2x1hri         ! create aux files: a2x 1hr instantaneous
   logical :: do_hist_a2x1hr          ! create aux files: a2x 1hr
   integer :: budget_inst             ! instantaneous budget flag
   integer :: budget_daily            ! daily budget flag
   integer :: budget_month            ! monthly budget flag
   integer :: budget_ann              ! annual budget flag
   integer :: budget_ltann            ! long term budget flag for end of year writing
   integer :: budget_ltend            ! long term budget flag for end of run writing

!  character(CL) :: hist_r2x_flds     = 'all'
!  character(CL) :: hist_l2x_flds     = 'all'
!  character(CL) :: hist_a2x24hr_flds = 'all'

   character(CL) :: hist_a2x_flds     = &
        'Faxa_swndr:Faxa_swvdr:Faxa_swndf:Faxa_swvdf'

   character(CL) :: hist_a2x3hrp_flds = &
        'Faxa_rainc:Faxa_rainl:Faxa_snowc:Faxa_snowl'

   character(CL) :: hist_a2x24hr_flds = &
        'Faxa_bcphiwet:Faxa_bcphodry:Faxa_bcphidry:Faxa_ocphiwet:Faxa_ocphidry:&
        &Faxa_ocphodry:Faxa_dstwet1:Faxa_dstdry1:Faxa_dstwet2:Faxa_dstdry2:Faxa_dstwet3:&
        &Faxa_dstdry3:Faxa_dstwet4:Faxa_dstdry4:Sa_co2prog:Sa_co2diag'

   character(CL) :: hist_a2x1hri_flds = &
        'Faxa_swndr:Faxa_swvdr:Faxa_swndf:Faxa_swvdf'

   character(CL) :: hist_a2x1hr_flds  = &
        'Sa_u:Sa_v'

   character(CL) :: hist_a2x3hr_flds  = &
        'Sa_z:Sa_topo:Sa_u:Sa_v:Sa_tbot:Sa_ptem:Sa_shum:Sa_dens:Sa_pbot:Sa_pslv:Faxa_lwdn:&
        &Faxa_rainc:Faxa_rainl:Faxa_snowc:Faxa_snowl:&
        &Faxa_swndr:Faxa_swvdr:Faxa_swndf:Faxa_swvdf:&
        &Sa_co2diag:Sa_co2prog'

   ! --- other ---
   integer  :: ka,km,k1,k2,k3         ! aVect field indices
   integer  :: ocnrun_count           ! number of times ocn run alarm went on
   logical  :: exists                 ! true if file exists
   integer  :: ierr                   ! MPI error return
   integer  :: rc                     ! return code
   logical  :: cdf64                  ! true => use 64 bit addressing in netCDF files

   character(*), parameter :: NLFileName = "drv_in"  ! input namelist filename

   integer  :: info_debug = 0         ! local info_debug level

   !----------------------------------------------------------------------------
   ! memory monitoring
   !----------------------------------------------------------------------------
   real(r8) :: msize,msize0,msize1     ! memory size (high water)
   real(r8) :: mrss ,mrss0 ,mrss1      ! resident size (current memory use)

   !----------------------------------------------------------------------------
   ! threading control
   !----------------------------------------------------------------------------
   integer  :: nthreads_GLOID         ! OMP global number of threads
   integer  :: nthreads_CPLID         ! OMP cpl number of threads
   integer  :: nthreads_ATMID         ! OMP atm number of threads
   integer  :: nthreads_LNDID         ! OMP lnd number of threads
   integer  :: nthreads_ICEID         ! OMP ice number of threads
   integer  :: nthreads_OCNID         ! OMP ocn number of threads
   integer  :: nthreads_GLCID         ! OMP glc number of threads
   integer  :: nthreads_ROFID         ! OMP glc number of threads
   integer  :: nthreads_WAVID         ! OMP wav number of threads
   integer  :: nthreads_ESPID         ! OMP esp number of threads

   integer  :: pethreads_GLOID        ! OMP number of threads per task

   logical  :: drv_threading          ! driver threading control

   !----------------------------------------------------------------------------
   ! communicator groups and related
   !----------------------------------------------------------------------------
   integer  :: Global_Comm

   integer  :: mpicom_GLOID          ! MPI global communicator
   integer  :: mpicom_CPLID          ! MPI cpl communicator
   integer  :: mpicom_OCNID          ! MPI ocn communicator for ensemble member 1

   integer  :: mpicom_CPLALLATMID    ! MPI comm for CPLALLATMID
   integer  :: mpicom_CPLALLLNDID    ! MPI comm for CPLALLLNDID
   integer  :: mpicom_CPLALLICEID    ! MPI comm for CPLALLICEID
   integer  :: mpicom_CPLALLOCNID    ! MPI comm for CPLALLOCNID
   integer  :: mpicom_CPLALLGLCID    ! MPI comm for CPLALLGLCID
   integer  :: mpicom_CPLALLROFID    ! MPI comm for CPLALLROFID
   integer  :: mpicom_CPLALLWAVID    ! MPI comm for CPLALLWAVID

   integer  :: iam_GLOID             ! pe number in global id
   logical  :: iamin_CPLID           ! pe associated with CPLID
   logical  :: iamroot_GLOID         ! GLOID masterproc
   logical  :: iamroot_CPLID         ! CPLID masterproc

   logical  :: iamin_CPLALLATMID     ! pe associated with CPLALLATMID
   logical  :: iamin_CPLALLLNDID     ! pe associated with CPLALLLNDID
   logical  :: iamin_CPLALLICEID     ! pe associated with CPLALLICEID
   logical  :: iamin_CPLALLOCNID     ! pe associated with CPLALLOCNID
   logical  :: iamin_CPLALLGLCID     ! pe associated with CPLALLGLCID
   logical  :: iamin_CPLALLROFID     ! pe associated with CPLALLROFID
   logical  :: iamin_CPLALLWAVID     ! pe associated with CPLALLWAVID

   !----------------------------------------------------------------------------
   ! complist: list of comps on this pe
   !----------------------------------------------------------------------------

   ! allow enough room for names of all physical components + coupler,
   ! where each string can be up to (max_inst_name_len+1) characters
   ! long (+1 allows for a space before each name)
   character(len=(seq_comm_namelen+1)*(num_inst_phys+1)) :: complist

   !----------------------------------------------------------------------------
   ! comp_num_<comp>: unique component number for each component type
   !----------------------------------------------------------------------------
   integer, parameter :: comp_num_atm = 1
   integer, parameter :: comp_num_lnd = 2
   integer, parameter :: comp_num_ice = 3
   integer, parameter :: comp_num_ocn = 4
   integer, parameter :: comp_num_glc = 5
   integer, parameter :: comp_num_rof = 6
   integer, parameter :: comp_num_wav = 7
   integer, parameter :: comp_num_esp = 8

   !----------------------------------------------------------------------------
   ! misc
   !----------------------------------------------------------------------------

   integer, parameter :: ens1=1         ! use first instance of ensemble only
   integer, parameter :: fix1=1         ! temporary hard-coding to first ensemble, needs to be fixed
   integer :: eai, eli, eoi, eii, egi, eri, ewi, eei, exi, efi  ! component instance counters

   !----------------------------------------------------------------------------
   ! formats
   !----------------------------------------------------------------------------
   character(*), parameter :: subname = '(seq_mct_drv)'
   character(*), parameter :: F00 = "('"//subname//" : ', 4A )"
   character(*), parameter :: F0L = "('"//subname//" : ', A, L6 )"
   character(*), parameter :: F0I = "('"//subname//" : ', A, 2i8 )"
   character(*), parameter :: F01 = "('"//subname//" : ', A, 2i8, 3x, A )"
   character(*), parameter :: F0R = "('"//subname//" : ', A, 2g23.15 )"
   character(*), parameter :: FormatA = '(A,": =============== ", A41,          " ===============")'
   character(*), parameter :: FormatD = '(A,": =============== ", A20,2I8,5x,   " ===============")'
   character(*), parameter :: FormatR = '(A,": =============== ", A31,F9.3,1x,  " ===============")'
   character(*), parameter :: FormatQ = '(A,": =============== ", A20,2F10.2,1x," ===============")'
!===============================================================================
contains
!===============================================================================

!===============================================================================
!*******************************************************************************
!===============================================================================

subroutine cesm_pre_init1()
   use shr_pio_mod, only : shr_pio_init1, shr_pio_init2
   implicit none

   !----------------------------------------------------------
   !| Initialize MCT and MPI communicators and IO
   !----------------------------------------------------------

   integer, dimension(num_inst_total) :: comp_id, comp_comm, comp_comm_iam
   logical :: comp_iamin(num_inst_total)
   character(len=seq_comm_namelen) :: comp_name(num_inst_total)
   integer :: i, it

   call mpi_init(ierr)
   call shr_mpi_chkerr(ierr,subname//' mpi_init')

   Global_Comm=MPI_COMM_WORLD
   comp_comm = MPI_COMM_NULL
   time_brun = mpi_wtime()

   call shr_pio_init1(num_inst_total,NLFileName, Global_Comm)
   !
   ! If pio_async_interface is true Global_Comm is MPI_COMM_NULL on the servernodes
   ! and server nodes do not return from shr_pio_init2
   !
   !   if (Global_Comm /= MPI_COMM_NULL) then

   call seq_comm_init(Global_Comm, NLFileName)

   !--- set task based threading counts ---
   call seq_comm_getinfo(GLOID,pethreads=pethreads_GLOID,iam=iam_GLOID)
   call seq_comm_setnthreads(pethreads_GLOID)

   !--- get some general data ---
   it=1
   call seq_comm_getinfo(GLOID,mpicom=mpicom_GLOID,&
        iamroot=iamroot_GLOID,nthreads=nthreads_GLOID)
   if (iamroot_GLOID) output_perf = .true.

   call seq_comm_getinfo(CPLID,mpicom=mpicom_CPLID,&
        iamroot=iamroot_CPLID,nthreads=nthreads_CPLID,&
        iam=comp_comm_iam(it))
   if (iamroot_CPLID) output_perf = .true.

   if (iamin_CPLID) complist = trim(complist)//' cpl'

   comp_id(it)    = CPLID
   comp_comm(it)  = mpicom_CPLID
   iamin_CPLID    = seq_comm_iamin(CPLID)
   comp_iamin(it) = seq_comm_iamin(comp_id(it))
   comp_name(it)  = seq_comm_name(comp_id(it))

   do eai = 1,num_inst_atm
      it=it+1
      comp_id(it)    = ATMID(eai)
      comp_iamin(it) = seq_comm_iamin(comp_id(it))
      comp_name(it)  = seq_comm_name(comp_id(it))
      call seq_comm_getinfo(ATMID(eai), mpicom=comp_comm(it), &
           nthreads=nthreads_ATMID, iam=comp_comm_iam(it))
      if (seq_comm_iamin(ATMID(eai))) then
         complist = trim(complist)//' '//trim(seq_comm_name(ATMID(eai)))
      endif
      if (seq_comm_iamroot(ATMID(eai))) output_perf = .true.
   enddo
   call seq_comm_getinfo(CPLALLATMID, mpicom=mpicom_CPLALLATMID)
   iamin_CPLALLATMID = seq_comm_iamin(CPLALLATMID)

   do eli = 1,num_inst_lnd
      it=it+1
      comp_id(it)    = LNDID(eli)
      comp_iamin(it) = seq_comm_iamin(comp_id(it))
      comp_name(it)  = seq_comm_name(comp_id(it))
      call seq_comm_getinfo(LNDID(eli), mpicom=comp_comm(it), &
           nthreads=nthreads_LNDID, iam=comp_comm_iam(it))
      if (seq_comm_iamin(LNDID(eli))) then
         complist = trim(complist)//' '//trim(seq_comm_name(LNDID(eli)))
      endif
      if (seq_comm_iamroot(LNDID(eli))) output_perf = .true.
   enddo
   call seq_comm_getinfo(CPLALLLNDID, mpicom=mpicom_CPLALLLNDID)
   iamin_CPLALLLNDID = seq_comm_iamin(CPLALLLNDID)

   do eoi = 1,num_inst_ocn
      it=it+1
      comp_id(it)    = OCNID(eoi)
      comp_iamin(it) = seq_comm_iamin(comp_id(it))
      comp_name(it)  = seq_comm_name(comp_id(it))
      call seq_comm_getinfo(OCNID(eoi), mpicom=comp_comm(it), &
           nthreads=nthreads_OCNID, iam=comp_comm_iam(it))
      if (seq_comm_iamin (OCNID(eoi))) then
         complist = trim(complist)//' '//trim(seq_comm_name(OCNID(eoi)))
      endif
      if (seq_comm_iamroot(OCNID(eoi))) output_perf = .true.
   enddo
   call seq_comm_getinfo(CPLALLOCNID, mpicom=mpicom_CPLALLOCNID)
   iamin_CPLALLOCNID = seq_comm_iamin(CPLALLOCNID)

   do eii = 1,num_inst_ice
      it=it+1
      comp_id(it)    = ICEID(eii)
      comp_iamin(it) = seq_comm_iamin(comp_id(it))
      comp_name(it)  = seq_comm_name(comp_id(it))
      call seq_comm_getinfo(ICEID(eii), mpicom=comp_comm(it), &
           nthreads=nthreads_ICEID, iam=comp_comm_iam(it))
      if (seq_comm_iamin (ICEID(eii))) then
         complist = trim(complist)//' '//trim(seq_comm_name(ICEID(eii)))
      endif
      if (seq_comm_iamroot(ICEID(eii))) output_perf = .true.
   enddo
   call seq_comm_getinfo(CPLALLICEID, mpicom=mpicom_CPLALLICEID)
   iamin_CPLALLICEID = seq_comm_iamin(CPLALLICEID)

   do egi = 1,num_inst_glc
      it=it+1
      comp_id(it)    = GLCID(egi)
      comp_iamin(it) = seq_comm_iamin(comp_id(it))
      comp_name(it)  = seq_comm_name(comp_id(it))
      call seq_comm_getinfo(GLCID(egi), mpicom=comp_comm(it), nthreads=nthreads_GLCID, iam=comp_comm_iam(it))
      if (seq_comm_iamin (GLCID(egi))) then
         complist = trim(complist)//' '//trim(seq_comm_name(GLCID(egi)))
      endif
      if (seq_comm_iamroot(GLCID(egi))) output_perf = .true.
   enddo
   call seq_comm_getinfo(CPLALLGLCID, mpicom=mpicom_CPLALLGLCID)
   iamin_CPLALLGLCID = seq_comm_iamin(CPLALLGLCID)

   do eri = 1,num_inst_rof
      it=it+1
      comp_id(it)    = ROFID(eri)
      comp_iamin(it) = seq_comm_iamin(comp_id(it))
      comp_name(it)  = seq_comm_name(comp_id(it))
      call seq_comm_getinfo(ROFID(eri), mpicom=comp_comm(it), &
           nthreads=nthreads_ROFID, iam=comp_comm_iam(it))
      if (seq_comm_iamin(ROFID(eri))) then
         complist = trim(complist)//' '//trim( seq_comm_name(ROFID(eri)))
      endif
      if (seq_comm_iamroot(ROFID(eri))) output_perf = .true.
   enddo
   call seq_comm_getinfo(CPLALLROFID, mpicom=mpicom_CPLALLROFID)
   iamin_CPLALLROFID = seq_comm_iamin(CPLALLROFID)

   do ewi = 1,num_inst_wav
      it=it+1
      comp_id(it)    = WAVID(ewi)
      comp_iamin(it) = seq_comm_iamin(comp_id(it))
      comp_name(it)  = seq_comm_name(comp_id(it))
      call seq_comm_getinfo(WAVID(ewi), mpicom=comp_comm(it), &
           nthreads=nthreads_WAVID, iam=comp_comm_iam(it))
      if (seq_comm_iamin(WAVID(ewi))) then
         complist = trim(complist)//' '//trim(seq_comm_name(WAVID(ewi)))
      endif
      if (seq_comm_iamroot(WAVID(ewi))) output_perf = .true.
   enddo
   call seq_comm_getinfo(CPLALLWAVID, mpicom=mpicom_CPLALLWAVID)
   iamin_CPLALLWAVID = seq_comm_iamin(CPLALLWAVID)

   do eei = 1,num_inst_esp
      it=it+1
      comp_id(it)    = ESPID(eei)
      comp_iamin(it) = seq_comm_iamin(comp_id(it))
      comp_name(it)  = seq_comm_name(comp_id(it))
      call seq_comm_getinfo(ESPID(eei), mpicom=comp_comm(it), &
           nthreads=nthreads_ESPID, iam=comp_comm_iam(it))
      if (seq_comm_iamin (ESPID(eei))) then
         complist = trim(complist)//' '//trim(seq_comm_name(ESPID(eei)))
      endif
   enddo
   ! ESP components do not use the coupler (they are 'external')

   !----------------------------------------------------------
   !| Set logging parameters both for shr code and locally
   !----------------------------------------------------------

   if (iamroot_CPLID) then
      inquire(file='cpl_modelio.nml',exist=exists)
      if (exists) then
         logunit = shr_file_getUnit()
         call shr_file_setIO('cpl_modelio.nml',logunit)
         call shr_file_setLogUnit(logunit)
         loglevel = 1
         call shr_file_setLogLevel(loglevel)
      endif
   else
      loglevel = 0
      call shr_file_setLogLevel(loglevel)
   endif

   !----------------------------------------------------------
   ! Log info about the environment settings
   !----------------------------------------------------------

   if (iamroot_CPLID) then
#ifdef USE_ESMF_LIB
      write(logunit,'(2A)') subname,' USE_ESMF_LIB is set'
#else
      write(logunit,'(2A)') subname,' USE_ESMF_LIB is NOT set, using esmf_wrf_timemgr'
#endif
      write(logunit,'(2A)') subname,' MCT_INTERFACE is set'
   endif

   !
   !  When using io servers (pio_async_interface=.true.) the server tasks do not return from
   !  shr_pio_init2
   !
   call shr_pio_init2(comp_id,comp_name,comp_iamin,comp_comm,comp_comm_iam)

end subroutine cesm_pre_init1

!===============================================================================
!*******************************************************************************
!===============================================================================

subroutine cesm_pre_init2()
   use pio, only : file_desc_t, pio_closefile, pio_file_is_open
   use shr_const_mod, only: shr_const_tkfrz, shr_const_tktrip, &
        shr_const_mwwv, shr_const_mwdair
   use shr_wv_sat_mod, only: shr_wv_sat_set_default, shr_wv_sat_init, &
        ShrWVSatTableSpec, shr_wv_sat_make_tables

   implicit none
   type(file_desc_t) :: pioid
   integer :: maxthreads

   character(CS) :: wv_sat_scheme
   real(r8) :: wv_sat_transition_start
   logical :: wv_sat_use_tables
   real(r8) :: wv_sat_table_spacing
   character(CL) :: errstring

   type(ShrWVSatTableSpec) :: liquid_spec, ice_spec, mixed_spec

   real(r8), parameter :: epsilo = shr_const_mwwv/shr_const_mwdair

   !----------------------------------------------------------
   ! Print Model heading and copyright message
   !----------------------------------------------------------

   if (iamroot_CPLID) call seq_cesm_printlogheader()

   !----------------------------------------------------------
   !| Timer initialization (has to be after mpi init)
   !----------------------------------------------------------
   maxthreads = max(nthreads_GLOID,nthreads_CPLID,nthreads_ATMID, &
        nthreads_LNDID,nthreads_ICEID,nthreads_OCNID,nthreads_GLCID, &
        nthreads_ROFID, nthreads_WAVID, nthreads_ESPID, pethreads_GLOID )

   call t_initf(NLFileName, LogPrint=.true., mpicom=mpicom_GLOID, &
        MasterTask=iamroot_GLOID,MaxThreads=maxthreads)

   if (iamin_CPLID) then
      call seq_io_cpl_init()
   endif

   call t_startf('CPL:INIT')
   call t_adj_detailf(+1)

   call t_startf('CPL:cesm_pre_init2')
   !----------------------------------------------------------
   !| Memory test
   !----------------------------------------------------------

!mt   call shr_mem_init(prt=.true.)
   call shr_mem_init(prt=iamroot_CPLID)

   !----------------------------------------------------------
   !| Initialize infodata
   !----------------------------------------------------------

   call seq_infodata_init(infodata,nlfilename, GLOID, pioid)

   !----------------------------------------------------------
   !| Initialize coupled fields (depends on infodata)
   !----------------------------------------------------------

   call seq_flds_set(nlfilename, GLOID, infodata)

   !----------------------------------------------------------
   !| Obtain infodata info
   !----------------------------------------------------------

   call seq_infodata_GetData(infodata, &
        info_debug=info_debug)

   if (info_debug > 1 .and. iamroot_CPLID) then
      write(logunit,*) ' '
      write(logunit,'(2A)') 'Status of infodata after seq_infodata_init'
      call seq_infodata_print( infodata )
      write(logunit,*) ' '
   endif

   call seq_infodata_GetData(infodata             , &
        read_restart=read_restart                 , &
        restart_file=rest_file                    , &
        timing_dir=timing_dir                     , &
        tchkpt_dir=tchkpt_dir                     , &
        info_debug=info_debug                     , &
        atm_present=atm_present                   , &
        lnd_present=lnd_present                   , &
        ice_present=ice_present                   , &
        ocn_present=ocn_present                   , &
        glc_present=glc_present                   , &
        rof_present=rof_present                   , &
        wav_present=wav_present                   , &
        esp_present=esp_present                   , &
        single_column=single_column               , &
        aqua_planet=aqua_planet                   , &
        cpl_seq_option=cpl_seq_option             , &
        drv_threading=drv_threading               , &
        do_histinit=do_histinit                   , &
        do_budgets=do_budgets                     , &
        budget_inst=budget_inst                   , &
        budget_daily=budget_daily                 , &
        budget_month=budget_month                 , &
        budget_ann=budget_ann                     , &
        budget_ltann=budget_ltann                 , &
        budget_ltend=budget_ltend                 , &
        histaux_a2x=do_hist_a2x                   , &
        histaux_a2x1hri=do_hist_a2x1hri           , &
        histaux_a2x1hr=do_hist_a2x1hr             , &
        histaux_a2x3hr =do_hist_a2x3hr            , &
        histaux_a2x3hrp=do_hist_a2x3hrp           , &
        histaux_a2x24hr=do_hist_a2x24hr           , &
        histaux_l2x=do_hist_l2x                   , &
        histaux_l2x1yr=do_hist_l2x1yr             , &
        histaux_r2x=do_hist_r2x                   , &
        run_barriers=run_barriers                 , &
        mct_usealltoall=mct_usealltoall           , &
        mct_usevector=mct_usevector               , &
        aoflux_grid=aoflux_grid                   , &
        vect_map=vect_map                         , &
        atm_gnam=atm_gnam                         , &
        lnd_gnam=lnd_gnam                         , &
        ocn_gnam=ocn_gnam                         , &
        ice_gnam=ice_gnam                         , &
        rof_gnam=rof_gnam                         , &
        glc_gnam=glc_gnam                         , &
        wav_gnam=wav_gnam                         , &
        tfreeze_option = tfreeze_option           , &
        cpl_decomp=seq_mctext_decomp              , &
        shr_map_dopole=shr_map_dopole             , &
        wall_time_limit=wall_time_limit           , &
        force_stop_at=force_stop_at               , &
        reprosum_use_ddpdd=reprosum_use_ddpdd     , &
        reprosum_diffmax=reprosum_diffmax         , &
        reprosum_recompute=reprosum_recompute, &
        max_cplstep_time=max_cplstep_time)

   ! above - cpl_decomp is set to pass the cpl_decomp value to seq_mctext_decomp
   ! (via a use statement)

   call shr_map_setDopole(shr_map_dopole)

   call shr_reprosum_setopts(&
        repro_sum_use_ddpdd_in    = reprosum_use_ddpdd, &
        repro_sum_rel_diff_max_in = reprosum_diffmax, &
        repro_sum_recompute_in    = reprosum_recompute)

   ! Check cpl_seq_option

   if (trim(cpl_seq_option) /= 'CESM1_ORIG' .and. &
       trim(cpl_seq_option) /= 'CESM1_ORIG_TIGHT' .and. &
       trim(cpl_seq_option) /= 'CESM1_MOD' .and. &
       trim(cpl_seq_option) /= 'CESM1_MOD_TIGHT' .and. &
       trim(cpl_seq_option) /= 'RASM_OPTION1' .and. &
       trim(cpl_seq_option) /= 'RASM_OPTION2' ) then
      call shr_sys_abort(subname//' invalid cpl_seq_option = '//trim(cpl_seq_option))
   endif

   !----------------------------------------------------------
   !| Test Threading Setup in driver
   !  happens to be valid on all pes for all IDs
   !----------------------------------------------------------

   if (drv_threading) then
      if (iamroot_GLOID) write(logunit,*) ' '
      if (iamroot_GLOID) write(logunit,'(2A)    ') subname,' Test Threading in driver'
      call seq_comm_setnthreads(nthreads_GLOID)
      if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_GLOID = ',&
           nthreads_GLOID,seq_comm_getnthreads()
      call seq_comm_setnthreads(nthreads_CPLID)
      if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_CPLID = ',&
           nthreads_CPLID,seq_comm_getnthreads()
      call seq_comm_setnthreads(nthreads_ATMID)
      if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_ATMID = ',&
           nthreads_ATMID,seq_comm_getnthreads()
      call seq_comm_setnthreads(nthreads_LNDID)
      if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_LNDID = ',&
           nthreads_LNDID,seq_comm_getnthreads()
      call seq_comm_setnthreads(nthreads_OCNID)
      if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_OCNID = ',&
           nthreads_OCNID,seq_comm_getnthreads()
      call seq_comm_setnthreads(nthreads_ICEID)
      if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_ICEID = ',&
           nthreads_ICEID,seq_comm_getnthreads()
      call seq_comm_setnthreads(nthreads_GLCID)
      if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_GLCID = ',&
           nthreads_GLCID,seq_comm_getnthreads()
      call seq_comm_setnthreads(nthreads_ROFID)
      if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_ROFID = ',&
           nthreads_ROFID,seq_comm_getnthreads()
      call seq_comm_setnthreads(nthreads_WAVID)
      if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_WAVID = ',&
           nthreads_WAVID,seq_comm_getnthreads()
      call seq_comm_setnthreads(nthreads_ESPID)
      if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_ESPID = ',&
           nthreads_ESPID,seq_comm_getnthreads()
      if (iamroot_GLOID) write(logunit,*) ' '

      call seq_comm_setnthreads(nthreads_GLOID)
   endif

   !----------------------------------------------------------
   !| Initialize time manager
   !----------------------------------------------------------

   call seq_timemgr_clockInit(seq_SyncClock, nlfilename, &
        read_restart, rest_file, pioid, mpicom_gloid,           &
        EClock_d, EClock_a, EClock_l, EClock_o,          &
        EClock_i, Eclock_g, Eclock_r, Eclock_w, Eclock_e)

   if (iamroot_CPLID) then
      call seq_timemgr_clockPrint(seq_SyncClock)
   endif

   call seq_infodata_getData(infodata,   &
        orb_iyear=orb_iyear,             &
        orb_iyear_align=orb_iyear_align, &
        orb_mode=orb_mode)

   !----------------------------------------------------------
   ! Initialize freezing point calculation for all components
   !----------------------------------------------------------

   call shr_frz_freezetemp_init(tfreeze_option)

   if (trim(orb_mode) == trim(seq_infodata_orb_variable_year)) then
      call seq_timemgr_EClockGetData( EClock_d, curr_ymd=ymd)

      call shr_cal_date2ymd(ymd,year,month,day)
      orb_cyear = orb_iyear + (year - orb_iyear_align)

      call shr_orb_params(orb_cyear, orb_eccen, orb_obliq, orb_mvelp, &
           orb_obliqr, orb_lambm0, orb_mvelpp, iamroot_CPLID)

      call seq_infodata_putData(infodata, &
           orb_eccen=orb_eccen,           &
           orb_obliqr=orb_obliqr,         &
           orb_lambm0=orb_lambm0,         &
           orb_mvelpp=orb_mvelpp)
   endif

   call seq_infodata_getData(infodata,                   &
        wv_sat_scheme=wv_sat_scheme,                     &
        wv_sat_transition_start=wv_sat_transition_start, &
        wv_sat_use_tables=wv_sat_use_tables,             &
        wv_sat_table_spacing=wv_sat_table_spacing)

   if (.not. shr_wv_sat_set_default(wv_sat_scheme)) then
      call shr_sys_abort('Invalid wv_sat_scheme.')
   end if

   call shr_wv_sat_init(shr_const_tkfrz, shr_const_tktrip, &
        wv_sat_transition_start, epsilo, errstring)

   if (errstring /= "") then
      call shr_sys_abort('shr_wv_sat_init: '//trim(errstring))
   end if

   ! The below produces internal lookup tables in the range 175-374K for
   ! liquid water, and 125-274K for ice, with a resolution set by the
   ! option wv_sat_table_spacing.
   ! In theory these ranges could be specified in the namelist, but in
   ! practice users will want to change them *very* rarely if ever, which
   ! is why only the spacing is in the namelist.
   if (wv_sat_use_tables) then
      liquid_spec = ShrWVSatTableSpec(ceiling(200._r8/wv_sat_table_spacing), &
           175._r8, wv_sat_table_spacing)
      ice_spec = ShrWVSatTableSpec(ceiling(150._r8/wv_sat_table_spacing), &
           125._r8, wv_sat_table_spacing)
      mixed_spec = ShrWVSatTableSpec(ceiling(250._r8/wv_sat_table_spacing), &
           125._r8, wv_sat_table_spacing)
      call shr_wv_sat_make_tables(liquid_spec, ice_spec, mixed_spec)
   end if

   call seq_infodata_putData(infodata, &
        atm_phase=1,                   &
        lnd_phase=1,                   &
        ocn_phase=1,                   &
        ice_phase=1,                   &
        glc_phase=1,                   &
        wav_phase=1,                   &
        esp_phase=1)

   !----------------------------------------------------------
   !| Set aqua_planet and single_column flags
   !  If in single column mode, overwrite flags according to focndomain file
   !  in ocn_in namelist. SCAM can reset the "present" flags for lnd,
   !  ocn, ice, rof, and flood.
   !----------------------------------------------------------

   if (.not.aqua_planet .and. single_column) then
      call seq_infodata_getData( infodata, &
           scmlon=scmlon, scmlat=scmlat)

      call seq_comm_getinfo(OCNID(ens1), mpicom=mpicom_OCNID)

      call shr_scam_checkSurface(scmlon, scmlat, &
           OCNID(ens1), mpicom_OCNID,            &
           lnd_present=lnd_present,              &
           ocn_present=ocn_present,              &
           ice_present=ice_present,              &
           rof_present=rof_present,              &
           flood_present=flood_present,          &
           rofice_present=rofice_present)

      call seq_infodata_putData(infodata,  &
           lnd_present=lnd_present,        &
           ocn_present=ocn_present,        &
           ice_present=ice_present,        &
           rof_present=rof_present,        &
           flood_present=flood_present,    &
           rofice_present=rofice_present)
   endif
   if(PIO_FILE_IS_OPEN(pioid)) then
      call pio_closefile(pioid)
   endif

   call t_stopf('CPL:cesm_pre_init2')

   call t_adj_detailf(-1)
   call t_stopf('CPL:INIT')

end subroutine cesm_pre_init2

!===============================================================================
!*******************************************************************************
!===============================================================================

subroutine cesm_init()

  implicit none

 101  format( A, 2i8, 12A, A, F8.2, A, F8.2 )
 102  format( A, 2i8, A, 8L3 )
 103  format( 5A )
 104  format( A, 2i8)
 105  format( A, 2i8, A, f10.2, A, f10.2, A, A, i5, A, A)
 106  format( A, f23.12)

   !-----------------------------------------------------------------------------
   !| Component Initialization
   !  Note that within each component initialization, the relevant x_present flag
   !  part of CESMInit can be modified
   !  By default, all these flags are set to true
   !  The atm can reset the lnd_present, ice_present and ocn_present flags based
   !  on aqua_planet, ideal_phys and adiabatic modes
   !  The stub components will reset the present flags to false, all other
   !  components will set them to true for the purposes of symmetry
   !-----------------------------------------------------------------------------

   call t_startf('cesm_init')
   call t_adj_detailf(+1)

   call t_startf('CPL:init_comps')
   if (iamroot_CPLID )then
      write(logunit,*) ' '
      write(logunit,F00) 'Initialize each component: atm, lnd, rof, ocn, ice, glc, wav, esp'
      call shr_sys_flush(logunit)
   endif

   call t_startf('comp_init_pre_all')
   call component_init_pre(atm, ATMID, CPLATMID, CPLALLATMID, infodata, ntype='atm')
   call component_init_pre(lnd, LNDID, CPLLNDID, CPLALLLNDID, infodata, ntype='lnd')
   call component_init_pre(rof, ROFID, CPLROFID, CPLALLROFID, infodata, ntype='rof')
   call component_init_pre(ocn, OCNID, CPLOCNID, CPLALLOCNID, infodata, ntype='ocn')
   call component_init_pre(ice, ICEID, CPLICEID, CPLALLICEID, infodata, ntype='ice')
   call component_init_pre(glc, GLCID, CPLGLCID, CPLALLGLCID, infodata, ntype='glc')
   call component_init_pre(wav, WAVID, CPLWAVID, CPLALLWAVID, infodata, ntype='wav')
   call component_init_pre(esp, ESPID, CPLESPID, CPLALLESPID, infodata, ntype='esp')
   call t_stopf('comp_init_pre_all')

   call t_startf('comp_init_cc_atm')
   call t_adj_detailf(+2)

   call component_init_cc(Eclock_a, atm, atm_init, infodata, NLFilename)
   call t_adj_detailf(-2)
   call t_stopf('comp_init_cc_atm')

   call t_startf('comp_init_cc_lnd')
   call t_adj_detailf(+2)
   call component_init_cc(Eclock_l, lnd, lnd_init, infodata, NLFilename)
   call t_adj_detailf(-2)
   call t_stopf('comp_init_cc_lnd')

   call t_startf('comp_init_cc_rof')
   call t_adj_detailf(+2)
   call component_init_cc(Eclock_r, rof, rof_init, infodata, NLFilename)
   call t_adj_detailf(-2)
   call t_stopf('comp_init_cc_rof')

   call t_startf('comp_init_cc_ocn')
   call t_adj_detailf(+2)
   call component_init_cc(Eclock_o, ocn, ocn_init, infodata, NLFilename)
   call t_adj_detailf(-2)
   call t_stopf('comp_init_cc_ocn')

   call t_startf('comp_init_cc_ice')
   call t_adj_detailf(+2)
   call component_init_cc(Eclock_i, ice, ice_init, infodata, NLFilename)
   call t_adj_detailf(-2)
   call t_stopf('comp_init_cc_ice')

   call t_startf('comp_init_cc_glc')
   call t_adj_detailf(+2)
   call component_init_cc(Eclock_g, glc, glc_init, infodata, NLFilename)
   call t_adj_detailf(-2)
   call t_stopf('comp_init_cc_glc')

   call t_startf('comp_init_cc_wav')
   call t_adj_detailf(+2)
   call component_init_cc(Eclock_w, wav, wav_init, infodata, NLFilename)
   call component_init_cc(Eclock_e, esp, esp_init, infodata, NLFilename)

   call t_adj_detailf(-2)
   call t_stopf('comp_init_cc_wav')

   call t_startf('comp_init_cx_all')
   call t_adj_detailf(+2)
   call component_init_cx(atm, infodata)
   call component_init_cx(lnd, infodata)
   call component_init_cx(rof, infodata)
   call component_init_cx(ocn, infodata)
   call component_init_cx(ice, infodata)
   call component_init_cx(glc, infodata)
   call component_init_cx(wav, infodata)
   call component_init_cx(esp, infodata)
   call t_adj_detailf(-2)
   call t_stopf('comp_init_cx_all')

   ! Determine complist (list of comps for each id)

   call t_startf('comp_list_all')
   call t_adj_detailf(+2)
   complist = " "
   if (iamin_CPLID) complist = trim(complist)//' cpl'

   do eai = 1,num_inst_atm
      iamin_ID = component_get_iamin_compid(atm(eai))
      if (iamin_ID) then
         compname = component_get_name(atm(eai))
         complist = trim(complist)//' '//trim(compname)
      endif
   enddo
   do eli = 1,num_inst_lnd
      iamin_ID = component_get_iamin_compid(lnd(eli))
      if (iamin_ID) then
         compname = component_get_name(lnd(eli))
         complist = trim(complist)//' '//trim(compname)
      endif
   enddo
   do eii = 1,num_inst_ice
      iamin_ID = component_get_iamin_compid(ice(eii))
      if (iamin_ID) then
         compname = component_get_name(ice(eii))
         complist = trim(complist)//' '//trim(compname)
      endif
   enddo
   do eoi = 1,num_inst_ocn
      iamin_ID = component_get_iamin_compid(ocn(eoi))
      if (iamin_ID) then
         compname = component_get_name(ocn(eoi))
         complist = trim(complist)//' '//trim(compname)
      endif
   enddo
   do egi = 1,num_inst_glc
      iamin_ID = component_get_iamin_compid(glc(egi))
      if (iamin_ID) then
         compname = component_get_name(glc(egi))
         complist = trim(complist)//' '//trim(compname)
      endif
   enddo
   do ewi = 1,num_inst_wav
      iamin_ID = component_get_iamin_compid(wav(ewi))
      if (iamin_ID) then
         compname = component_get_name(wav(ewi))
         complist = trim(complist)//' '//trim(compname)
      endif
   enddo

   do eei = 1,num_inst_esp
      iamin_ID = component_get_iamin_compid(esp(eei))
      if (iamin_ID) then
         compname = component_get_name(esp(eei))
         complist = trim(complist)//' '//trim(compname)
      endif
   enddo

   call t_adj_detailf(-2)
   call t_stopf('comp_list_all')

   call t_stopf('CPL:init_comps')
   !----------------------------------------------------------
   !| Determine coupling interactions based on present and prognostic flags
   !----------------------------------------------------------

   if (iamin_CPLALLATMID) call seq_infodata_exchange(infodata,CPLALLATMID,'cpl2atm_init')
   if (iamin_CPLALLLNDID) call seq_infodata_exchange(infodata,CPLALLLNDID,'cpl2lnd_init')
   if (iamin_CPLALLOCNID) call seq_infodata_exchange(infodata,CPLALLOCNID,'cpl2ocn_init')
   if (iamin_CPLALLICEID) call seq_infodata_exchange(infodata,CPLALLICEID,'cpl2ice_init')
   if (iamin_CPLALLGLCID) call seq_infodata_exchange(infodata,CPLALLGLCID,'cpl2glc_init')
   if (iamin_CPLALLROFID) call seq_infodata_exchange(infodata,CPLALLROFID,'cpl2rof_init')
   if (iamin_CPLALLWAVID) call seq_infodata_exchange(infodata,CPLALLWAVID,'cpl2wav_init')

   if (iamroot_CPLID) then
      write(logunit,F00) 'Determine final settings for presence of surface components'
      call shr_sys_flush(logunit)
   endif

   call seq_infodata_getData(infodata,         &
        atm_present=atm_present,               &
        lnd_present=lnd_present,               &
        ice_present=ice_present,               &
        ocn_present=ocn_present,               &
        glc_present=glc_present,               &
        glclnd_present=glclnd_present,         &
        glcocn_present=glcocn_present,         &
        glcice_present=glcice_present,         &
        rof_present=rof_present,               &
        rofice_present=rofice_present,         &
        wav_present=wav_present,               &
        esp_present=esp_present,               &
        flood_present=flood_present,           &
        atm_prognostic=atm_prognostic,         &
        lnd_prognostic=lnd_prognostic,         &
        ice_prognostic=ice_prognostic,         &
        iceberg_prognostic=iceberg_prognostic, &
        ocn_prognostic=ocn_prognostic,         &
        ocnrof_prognostic=ocnrof_prognostic,   &
        glc_prognostic=glc_prognostic,         &
        rof_prognostic=rof_prognostic,         &
        wav_prognostic=wav_prognostic,         &
        esp_prognostic=esp_prognostic,         &
        dead_comps=dead_comps,                 &
        esmf_map_flag=esmf_map_flag,           &
        atm_nx=atm_nx, atm_ny=atm_ny,          &
        lnd_nx=lnd_nx, lnd_ny=lnd_ny,          &
        rof_nx=rof_nx, rof_ny=rof_ny,          &
        ice_nx=ice_nx, ice_ny=ice_ny,          &
        glc_nx=glc_nx, glc_ny=glc_ny,          &
        ocn_nx=ocn_nx, ocn_ny=ocn_ny,          &
        wav_nx=wav_nx, wav_ny=wav_ny,          &
        cpl_cdf64=cdf64,                       &
        atm_aero=atm_aero )

   ! derive samegrid flags

   samegrid_ao  = .true.
   samegrid_al  = .true.
   samegrid_lr  = .true.
   samegrid_oi  = .true.
   samegrid_ro  = .true.
   samegrid_aw  = .true.
   samegrid_ow  = .true.
   samegrid_lg  = .true.
   samegrid_og  = .true.
   samegrid_ig  = .true.
   samegrid_alo = .true.

   ! set samegrid to true for single column
   if (.not. single_column) then
      if (trim(atm_gnam) /= trim(ocn_gnam)) samegrid_ao = .false.
      if (trim(atm_gnam) /= trim(lnd_gnam)) samegrid_al = .false.
      if (trim(lnd_gnam) /= trim(rof_gnam)) samegrid_lr = .false.
      if (trim(rof_gnam) /= trim(ocn_gnam)) samegrid_ro = .false.
      if (trim(ocn_gnam) /= trim(ice_gnam)) samegrid_oi = .false.
      if (trim(atm_gnam) /= trim(wav_gnam)) samegrid_aw = .false.
      if (trim(ocn_gnam) /= trim(wav_gnam)) samegrid_ow = .false.
      if (trim(lnd_gnam) /= trim(glc_gnam)) samegrid_lg = .false.
      if (trim(ocn_gnam) /= trim(glc_gnam)) samegrid_og = .false.
      if (trim(ice_gnam) /= trim(glc_gnam)) samegrid_ig = .false.
      samegrid_alo = (samegrid_al .and. samegrid_ao)
   endif

   ! derive coupling connection flags

   atm_c2_lnd = .false.
   atm_c2_ocn = .false.
   atm_c2_ice = .false.
   atm_c2_wav = .false.
   lnd_c2_atm = .false.
   lnd_c2_rof = .false.
   lnd_c2_glc = .false.
   ocn_c2_atm = .false.
   ocn_c2_ice = .false.
   ocn_c2_wav = .false.
   ice_c2_atm = .false.
   ice_c2_ocn = .false.
   ice_c2_wav = .false.
   rof_c2_lnd = .false.
   rof_c2_ocn = .false.
   rof_c2_ice = .false.
   glc_c2_lnd = .false.
   glc_c2_ocn = .false.
   glc_c2_ice = .false.
   wav_c2_ocn = .false.

   if (atm_present) then
      if (lnd_prognostic) atm_c2_lnd = .true.
      if (ocn_prognostic) atm_c2_ocn = .true.
      if (ocn_present   ) atm_c2_ocn = .true. ! needed for aoflux calc if aoflux=ocn
      if (ice_prognostic) atm_c2_ice = .true.
      if (wav_prognostic) atm_c2_wav = .true.
   endif
   if (lnd_present) then
      if (atm_prognostic) lnd_c2_atm = .true.
      if (rof_prognostic) lnd_c2_rof = .true.
      if (glc_prognostic) lnd_c2_glc = .true.
   endif
   if (ocn_present) then
      if (atm_prognostic) ocn_c2_atm = .true.
      if (atm_present   ) ocn_c2_atm = .true. ! needed for aoflux calc if aoflux=atm
      if (ice_prognostic) ocn_c2_ice = .true.
      if (wav_prognostic) ocn_c2_wav = .true.
   endif
   if (ice_present) then
      if (atm_prognostic) ice_c2_atm = .true.
      if (ocn_prognostic) ice_c2_ocn = .true.
      if (wav_prognostic) ice_c2_wav = .true.
   endif
   if (rof_present) then
      if (lnd_prognostic   ) rof_c2_lnd = .true.
      if (ocnrof_prognostic) rof_c2_ocn = .true.
      if (rofice_present .and. iceberg_prognostic) rof_c2_ice = .true.
   endif
   if (glc_present) then
      if (glclnd_present .and. lnd_prognostic) glc_c2_lnd = .true.
      if (glcocn_present .and. ocn_prognostic) glc_c2_ocn = .true.
      if (glcice_present .and. iceberg_prognostic) glc_c2_ice = .true.
   endif
   if (wav_present) then
      if (ocn_prognostic) wav_c2_ocn = .true.
   endif

   !----------------------------------------------------------
   ! Set domain check and other flag
   !----------------------------------------------------------

   domain_check = .true.
   if (single_column         ) domain_check = .false.
   if (dead_comps            ) domain_check = .false.

   ! set skip_ocean_run flag, used primarily for ocn run on first timestep
   ! use reading a restart as a surrogate from whether this is a startup run

   skip_ocean_run = .true.
   if ( read_restart) skip_ocean_run = .false.
   ocnrun_count = 0
   cpl2ocn_first = .true.

   do_histavg = .true.
   if (seq_timemgr_histavg_type == seq_timemgr_type_never) then
      do_histavg = .false.
   endif

   !----------------------------------------------------------
   !| Write component and coupler setup information
   !----------------------------------------------------------

   if (iamroot_CPLID) then
      write(logunit,*  )' '
      write(logunit,F00)'After component initialization:'
      write(logunit,F0L)'atm model present     = ',atm_present
      write(logunit,F0L)'lnd model present     = ',lnd_present
      write(logunit,F0L)'ocn model present     = ',ocn_present
      write(logunit,F0L)'ice model present     = ',ice_present
      write(logunit,F0L)'glc model present     = ',glc_present
      write(logunit,F0L)'glc/lnd   present     = ',glclnd_present
      write(logunit,F0L)'glc/ocn   present     = ',glcocn_present
      write(logunit,F0L)'glc/ice   present     = ',glcice_present
      write(logunit,F0L)'rof model present     = ',rof_present
      write(logunit,F0L)'rof/ice   present     = ',rofice_present
      write(logunit,F0L)'rof/flood present     = ',flood_present
      write(logunit,F0L)'wav model present     = ',wav_present
      write(logunit,F0L)'esp model present     = ',esp_present

      write(logunit,F0L)'atm model prognostic  = ',atm_prognostic
      write(logunit,F0L)'lnd model prognostic  = ',lnd_prognostic
      write(logunit,F0L)'ocn model prognostic  = ',ocn_prognostic
      write(logunit,F0L)'ice model prognostic  = ',ice_prognostic
      write(logunit,F0L)'iceberg   prognostic  = ',iceberg_prognostic
      write(logunit,F0L)'glc model prognostic  = ',glc_prognostic
      write(logunit,F0L)'rof model prognostic  = ',rof_prognostic
      write(logunit,F0L)'ocn rof   prognostic  = ',ocnrof_prognostic
      write(logunit,F0L)'wav model prognostic  = ',wav_prognostic
      write(logunit,F0L)'esp model prognostic  = ',esp_prognostic

      write(logunit,F0L)'atm_c2_lnd            = ',atm_c2_lnd
      write(logunit,F0L)'atm_c2_ocn            = ',atm_c2_ocn
      write(logunit,F0L)'atm_c2_ice            = ',atm_c2_ice
      write(logunit,F0L)'atm_c2_wav            = ',atm_c2_wav
      write(logunit,F0L)'lnd_c2_atm            = ',lnd_c2_atm
      write(logunit,F0L)'lnd_c2_rof            = ',lnd_c2_rof
      write(logunit,F0L)'lnd_c2_glc            = ',lnd_c2_glc
      write(logunit,F0L)'ocn_c2_atm            = ',ocn_c2_atm
      write(logunit,F0L)'ocn_c2_ice            = ',ocn_c2_ice
      write(logunit,F0L)'ocn_c2_wav            = ',ocn_c2_wav
      write(logunit,F0L)'ice_c2_atm            = ',ice_c2_atm
      write(logunit,F0L)'ice_c2_ocn            = ',ice_c2_ocn
      write(logunit,F0L)'ice_c2_wav            = ',ice_c2_wav
      write(logunit,F0L)'rof_c2_lnd            = ',rof_c2_lnd
      write(logunit,F0L)'rof_c2_ocn            = ',rof_c2_ocn
      write(logunit,F0L)'rof_c2_ice            = ',rof_c2_ice
      write(logunit,F0L)'glc_c2_lnd            = ',glc_c2_lnd
      write(logunit,F0L)'glc_c2_ocn            = ',glc_c2_ocn
      write(logunit,F0L)'glc_c2_ice            = ',glc_c2_ice
      write(logunit,F0L)'wav_c2_ocn            = ',wav_c2_ocn

      write(logunit,F0L)'dead components       = ',dead_comps
      write(logunit,F0L)'domain_check          = ',domain_check
      write(logunit,F01)'atm_nx,atm_ny         = ',atm_nx,atm_ny,trim(atm_gnam)
      write(logunit,F01)'lnd_nx,lnd_ny         = ',lnd_nx,lnd_ny,trim(lnd_gnam)
      write(logunit,F01)'rof_nx,rof_ny         = ',rof_nx,rof_ny,trim(rof_gnam)
      write(logunit,F01)'ice_nx,ice_ny         = ',ice_nx,ice_ny,trim(ice_gnam)
      write(logunit,F01)'ocn_nx,ocn_ny         = ',ocn_nx,ocn_ny,trim(ocn_gnam)
      write(logunit,F01)'glc_nx,glc_ny         = ',glc_nx,glc_ny,trim(glc_gnam)
      write(logunit,F01)'wav_nx,wav_ny         = ',wav_nx,wav_ny,trim(wav_gnam)
      write(logunit,F0L)'samegrid_ao           = ',samegrid_ao
      write(logunit,F0L)'samegrid_al           = ',samegrid_al
      write(logunit,F0L)'samegrid_ro           = ',samegrid_ro
      write(logunit,F0L)'samegrid_aw           = ',samegrid_aw
      write(logunit,F0L)'samegrid_ow           = ',samegrid_ow
      write(logunit,F0L)'skip init ocean run   = ',skip_ocean_run
      write(logunit,F00)'cpl sequence option   = ',trim(cpl_seq_option)
      write(logunit,F0L)'cpl_cdf64             = ',cdf64
      write(logunit,F0L)'do_histavg            = ',do_histavg
      write(logunit,F0L)'atm_aero              = ',atm_aero
      write(logunit,*  )' '
      call shr_sys_flush(logunit)
   endif

   !----------------------------------------------------------
   !| Present and prognostic consistency checks
   !----------------------------------------------------------

   if (atm_prognostic .and. .not.atm_present) then
      call shr_sys_abort(subname//' ERROR: if prognostic atm must also have atm present')
   endif
   if (ocn_prognostic .and. .not.ocn_present) then
      call shr_sys_abort(subname//' ERROR: if prognostic ocn must also have ocn present')
   endif
   if (lnd_prognostic .and. .not.lnd_present) then
      call shr_sys_abort(subname//' ERROR: if prognostic lnd must also have lnd present')
   endif
   if (ice_prognostic .and. .not.ice_present) then
      call shr_sys_abort(subname//' ERROR: if prognostic ice must also have ice present')
   endif
   if (iceberg_prognostic .and. .not.ice_prognostic) then
      call shr_sys_abort(subname//' ERROR: if prognostic iceberg must also have ice prognostic')
   endif
   if (glc_prognostic .and. .not.glc_present) then
      call shr_sys_abort(subname//' ERROR: if prognostic glc must also have glc present')
   endif
   if (rof_prognostic .and. .not.rof_present) then
      call shr_sys_abort(subname//' ERROR: if prognostic rof must also have rof present')
   endif
   if (wav_prognostic .and. .not.wav_present) then
      call shr_sys_abort(subname//' ERROR: if prognostic wav must also have wav present')
   endif
   if (esp_prognostic .and. .not.esp_present) then
      call shr_sys_abort(subname//' ERROR: if prognostic esp must also have esp present')
   endif
#ifndef CPL_BYPASS
   if ((ice_prognostic .or. ocn_prognostic .or. lnd_prognostic) .and. .not. atm_present) then
      call shr_sys_abort(subname//' ERROR: if prognostic surface model must also have atm present')
   endif
#endif
   if ((glclnd_present .or. glcocn_present .or. glcice_present) .and. .not.glc_present) then
      call shr_sys_abort(subname//' ERROR: if glcxxx present must also have glc present')
   endif
   if (rofice_present .and. .not.rof_present) then
      call shr_sys_abort(subname//' ERROR: if rofice present must also have rof present')
   endif
   if (ocnrof_prognostic .and. .not.rof_present) then
      if (iamroot_CPLID) then
         write(logunit,F00) 'WARNING: ocnrof_prognostic is TRUE but rof_present is FALSE'
         call shr_sys_flush(logunit)
      endif
   endif

   !----------------------------------------------------------
   !| Samegrid checks
   !----------------------------------------------------------

   if (.not. samegrid_oi) then
      call shr_sys_abort(subname//' ERROR: samegrid_oi is false')
   endif

   !----------------------------------------------------------
   !| Check instances of prognostic components
   !----------------------------------------------------------

   if (atm_prognostic .and. num_inst_atm /= num_inst_max) &
      call shr_sys_abort(subname//' ERROR: atm_prognostic but num_inst_atm not num_inst_max')
   if (lnd_prognostic .and. num_inst_lnd /= num_inst_max) &
      call shr_sys_abort(subname//' ERROR: lnd_prognostic but num_inst_lnd not num_inst_max')
   if (ocn_prognostic .and. (num_inst_ocn /= num_inst_max .and. num_inst_ocn /= 1)) &
      call shr_sys_abort(subname//' ERROR: ocn_prognostic but num_inst_ocn not 1 or num_inst_max')
   if (ice_prognostic .and. num_inst_ice /= num_inst_max) &
      call shr_sys_abort(subname//' ERROR: ice_prognostic but num_inst_ice not num_inst_max')
   if (glc_prognostic .and. num_inst_glc /= num_inst_max) &
      call shr_sys_abort(subname//' ERROR: glc_prognostic but num_inst_glc not num_inst_max')
   if (rof_prognostic .and. num_inst_rof /= num_inst_max) &
      call shr_sys_abort(subname//' ERROR: rof_prognostic but num_inst_rof not num_inst_max')
   if (wav_prognostic .and. num_inst_wav /= num_inst_max) &
      call shr_sys_abort(subname//' ERROR: wav_prognostic but num_inst_wav not num_inst_max')

   !----------------------------------------------------------
   !| Initialize attribute vectors for prep_c2C_init_avs routines and fractions
   !| Initialize mapping between components
   !----------------------------------------------------------

   if (iamin_CPLID) then

      call t_startf('CPL:init_maps')
      call t_adj_detailf(+2)
      if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)

      call prep_atm_init(infodata, ocn_c2_atm, ice_c2_atm, lnd_c2_atm)

      call prep_lnd_init(infodata, atm_c2_lnd, rof_c2_lnd, glc_c2_lnd)

      call prep_ocn_init(infodata, atm_c2_ocn, atm_c2_ice, ice_c2_ocn, rof_c2_ocn, wav_c2_ocn, glc_c2_ocn)

      call prep_ice_init(infodata, ocn_c2_ice, glc_c2_ice, rof_c2_ice )

      call prep_rof_init(infodata, lnd_c2_rof)

      call prep_glc_init(infodata, lnd_c2_glc)

      call prep_wav_init(infodata, atm_c2_wav, ocn_c2_wav, ice_c2_wav)

      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
      call t_adj_detailf(-2)
      call t_stopf('CPL:init_maps')

   endif

   !----------------------------------------------------------
   !| Update aream in domains where appropriate
   !----------------------------------------------------------

   if (iamin_CPLID) then
      call t_startf ('CPL:init_aream')
      call t_adj_detailf(+2)

      if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)

      call component_init_aream(infodata, rof_c2_ocn, samegrid_ao, samegrid_al, &
           samegrid_ro, samegrid_lg)

      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)

      call t_adj_detailf(-2)
      call t_stopf ('CPL:init_aream')
   endif ! iamin_CPLID

   !----------------------------------------------------------
   !| Check domains
   !  This must be done after the mappers are initialized since
   !  checking is done on each processor and not with a global gather
   !----------------------------------------------------------

   if (iamin_CPLID) then
      call t_startf ('CPL:init_domain_check')
      call t_adj_detailf(+2)

      if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
      if (domain_check) then
         if (iamroot_CPLID) then
            write(logunit,*) ' '
            write(logunit,F00) 'Performing domain checking'
            call shr_sys_flush(logunit)
         endif

         call seq_domain_check( infodata,                                             &
              atm(ens1), ice(ens1), lnd(ens1), ocn(ens1), rof(ens1), glc(ens1),       &
              samegrid_al, samegrid_ao, samegrid_ro, samegrid_lg)

      endif
      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)

      call t_adj_detailf(-2)
      call t_stopf ('CPL:init_domain_check')
   endif ! iamin_CPLID

   !----------------------------------------------------------
   !| Initialize area corrections based on aream (read in map_init) and area
   !| Area correct component initialization output fields
   !| Map initial component AVs from component to coupler pes
   !----------------------------------------------------------

   areafact_samegrid = .false.
#if (defined BFB_CAM_SCAM_IOP )
   if (.not.samegrid_alo) then
      call shr_sys_abort(subname//' ERROR: samegrid_alo is false - Must run with same atm/ocn/lnd grids when configured for scam iop')
   else
      areafact_samegrid = .true.
   endif
#endif
   if (single_column) areafact_samegrid = .true.

   call t_startf ('CPL:init_areacor')
   call t_adj_detailf(+2)

   call mpi_barrier(mpicom_GLOID,ierr)
   if (atm_present) call component_init_areacor(atm, areafact_samegrid, seq_flds_a2x_fluxes)

   call mpi_barrier(mpicom_GLOID,ierr)
   if (lnd_present) call component_init_areacor(lnd, areafact_samegrid, seq_flds_l2x_fluxes)

   call mpi_barrier(mpicom_GLOID,ierr)
   if (rof_present) call component_init_areacor(rof, areafact_samegrid, seq_flds_r2x_fluxes)

   call mpi_barrier(mpicom_GLOID,ierr)
   if (ocn_present) call component_init_areacor(ocn, areafact_samegrid, seq_flds_o2x_fluxes)

   call mpi_barrier(mpicom_GLOID,ierr)
   if (ice_present) call component_init_areacor(ice, areafact_samegrid, seq_flds_i2x_fluxes)

   call mpi_barrier(mpicom_GLOID,ierr)
   if (glc_present) call component_init_areacor(glc, areafact_samegrid, seq_flds_g2x_fluxes)

   call mpi_barrier(mpicom_GLOID,ierr)
   if (wav_present) call component_init_areacor(wav, areafact_samegrid, seq_flds_w2x_fluxes)

   call t_adj_detailf(-2)
   call t_stopf ('CPL:init_areacor')

   !----------------------------------------------------------
   !| global sum diagnostics for IC data
   !----------------------------------------------------------

   if (iamin_CPLID .and. info_debug > 1) then
      call t_startf ('CPL:init_diag')
      call t_adj_detailf(+2)

      if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
      if (atm_present) then
         call component_diag(infodata, atm, flow='c2x', comment='recv IC atm', &
              info_debug=info_debug)
      endif
      if (ice_present) then
         call component_diag(infodata, ice, flow='c2x', comment='recv IC ice', &
              info_debug=info_debug)
      endif
      if (lnd_present) then
         call component_diag(infodata, lnd, flow='c2x', comment='recv IC lnd', &
              info_debug=info_debug)
      endif
      if (rof_present) then
         call component_diag(infodata, rof, flow='c2x', comment='recv IC rof', &
              info_debug=info_debug)
      endif
      if (ocn_present) then
         call component_diag(infodata, ocn, flow='c2x', comment='recv IC ocn', &
              info_debug=info_debug)
      endif
      if (glc_present) then
         call component_diag(infodata, glc, flow='c2x', comment='recv IC glc', &
              info_debug=info_debug)
      endif
      if (wav_present) then
         call component_diag(infodata, wav, flow='c2x', comment='recv IC wav', &
              info_debug=info_debug)
      endif
      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)

      call t_adj_detailf(-2)
      call t_stopf ('CPL:init_diag')
   endif

   !----------------------------------------------------------
   !| Initialize fractions
   !----------------------------------------------------------

   if (iamin_CPLID) then
      call t_startf ('CPL:init_fracs')
      call t_adj_detailf(+2)

      allocate(fractions_ax(num_inst_frc))
      allocate(fractions_lx(num_inst_frc))
      allocate(fractions_ox(num_inst_frc))
      allocate(fractions_ix(num_inst_frc))
      allocate(fractions_gx(num_inst_frc))
      allocate(fractions_rx(num_inst_frc))
      allocate(fractions_wx(num_inst_frc))

      if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
      do efi = 1,num_inst_frc
         eii = mod((efi-1),num_inst_ice) + 1

         if (iamroot_CPLID) then
            write(logunit,*) ' '
            if (efi == 1) write(logunit,F00) 'Initializing fractions'
         endif

         call seq_frac_init(infodata,                                  &
              atm(ens1), ice(ens1), lnd(ens1),                         &
              ocn(ens1), glc(ens1), rof(ens1),                         &
              wav(ens1),                                               &
              fractions_ax(efi), fractions_ix(efi), fractions_lx(efi), &
              fractions_ox(efi), fractions_gx(efi), fractions_rx(efi), &
              fractions_wx(efi))

         if (iamroot_CPLID) then
            write(logunit,*) ' '
            if (efi == 1) write(logunit,F00) 'Setting fractions'
         endif

         call seq_frac_set(infodata, ice(eii), &
              fractions_ax(efi), fractions_ix(efi), fractions_ox(efi))

      enddo
      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)

      call t_adj_detailf(-2)
      call t_stopf ('CPL:init_fracs')
   endif

   !----------------------------------------------------------
   !| Initialize prep_aoflux_mod module variables
   !----------------------------------------------------------

   if (iamin_CPLID) then
      call prep_aoflux_init(infodata, fractions_ox, fractions_ax)
   endif

   !----------------------------------------------------------
   !| Initialize atm/ocn flux component and compute ocean albedos
   !----------------------------------------------------------

   if (iamin_CPLID) then
      if (ocn_present) then
         call t_startf ('CPL:init_aoflux')
         call t_adj_detailf(+2)

         if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
         if (iamroot_CPLID) then
            write(logunit,*) ' '
            write(logunit,F00) 'Initializing atm/ocn flux component'
         endif

         if (trim(aoflux_grid) == 'ocn') then

            call seq_flux_init_mct(ocn(ens1), fractions_ox(ens1))

         elseif (trim(aoflux_grid) == 'atm') then

            call seq_flux_init_mct(atm(ens1), fractions_ax(ens1))

         elseif (trim(aoflux_grid) == 'exch') then

            call shr_sys_abort(subname//' aoflux_grid = exch not validated')
            call seq_flux_initexch_mct(atm(ens1), ocn(ens1), mpicom_cplid, cplid)

         else
            call shr_sys_abort(subname//' aoflux_grid = '//trim(aoflux_grid)//' not available')

         endif

         do exi = 1,num_inst_xao
            !tcx is this correct? relation between xao and frc for ifrad and ofrad
            efi = mod((exi-1),num_inst_frc) + 1
            eai = mod((exi-1),num_inst_atm) + 1
            xao_ox => prep_aoflux_get_xao_ox()        ! array over all instances
            a2x_ox => prep_ocn_get_a2x_ox()
            call seq_flux_ocnalb_mct(infodata, ocn(1), a2x_ox(eai), fractions_ox(efi), xao_ox(exi))
         enddo

         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)

         call t_adj_detailf(-2)
         call t_stopf ('CPL:init_aoflux')
      endif
   endif

   !----------------------------------------------------------
   !| ATM PREP for recalculation of initial solar
   !  Note that ocean albedos are ALWAYS CALCULATED on the ocean grid
   !  If aoflux_grid = 'ocn' , xao_ox is input for atm/ocn fluxes and xao_ax is output
   !  If aoflux_grid = 'atm' , xao_ax is input for atm/ocn fluxes and xao_ox is not used
   !  If aoflux_grid = 'exch', xao_ax is input for atm/ocn /fluxes and xao_ox is not used
   !  Merge atmosphere input state and run atmospheric radiation
   !----------------------------------------------------------

   if (atm_prognostic) then
      if (iamin_CPLID) then

         if (lnd_present) then
            ! Get lnd output on atm grid
            call prep_atm_calc_l2x_ax(fractions_lx, timer='CPL:init_atminit')
         endif

         if (ice_present) then
            ! Get ice output on atm grid
            call prep_atm_calc_i2x_ax(fractions_ix, timer='CPL:init_atminit')
         endif

         if (ocn_present) then
            ! Get ocn output on atm grid
            call prep_atm_calc_o2x_ax(fractions_ox, timer='CPL:init_atminit')
         endif

         if (ocn_present) then
            ! Get albedos on atm grid
            call prep_aoflux_calc_xao_ax(fractions_ox, flds='albedos', timer='CPL:init_atminit')

            ! Get atm/ocn fluxes on atm grid
            if (trim(aoflux_grid) == 'ocn') then
               call prep_aoflux_calc_xao_ax(fractions_ox, flds='states_and_fluxes', &
                    timer='CPL:init_atminit')
            endif
         endif

         if (lnd_present .or. ocn_present) then
            ! Merge input to atmosphere on coupler pes
            xao_ax => prep_aoflux_get_xao_ax()
            if (associated(xao_ax)) then
               call  prep_atm_mrg(infodata, &
                    fractions_ax=fractions_ax, xao_ax=xao_ax, timer_mrg='CPL:init_atminit')
            endif
         endif

         call component_diag(infodata, atm, flow='x2c', comment='send atm', info_debug=info_debug)

      endif

   endif  ! atm_prognostic

   !----------------------------------------------------------
   !| Second phase of atmosphere component initialization
   !  Recalculate solar based on input albedo's from surface components.
   !  Data or dead atmosphere may just return on this phase.
   !----------------------------------------------------------

   if (atm_present) then
      call t_startf('comp_init_cc_atm2')
      call t_adj_detailf(+2)

      if (iamroot_CPLID) then
         write(logunit,F00) 'Calling atm_init_mct phase 2'
      endif

      ! Send atm input data from coupler pes to atm pes
      if (atm_prognostic) then
         call component_exch(atm, flow='x2c', infodata=infodata, &
              infodata_string='cpl2atm_init')
      endif

      ! Set atm init phase to 2 for all atm instances on component instance pes
      do eai = 1,num_inst_atm
         if (component_get_iamin_compid(atm(eai))) then
            call seq_infodata_putData(infodata, atm_phase=2)
         endif
      enddo

      ! Run atm_init_mct with init phase of 2
      call component_init_cc(Eclock_a, atm, atm_init,                   &
           infodata, NLFilename,                                        &
           seq_flds_x2c_fluxes=seq_flds_x2a_fluxes,                     &
           seq_flds_c2x_fluxes=seq_flds_a2x_fluxes)

      ! Map atm output data from atm pes to cpl pes
      call component_exch(atm, flow='c2x', infodata=infodata, &
           infodata_string='atm2cpl_init')

      if (iamin_CPLID) then
         if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
         call component_diag(infodata, atm, flow='c2x', comment= 'recv IC2 atm', &
              info_debug=info_debug)
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
      endif

      call t_adj_detailf(-2)
      call t_stopf('comp_init_cc_atm2')
   endif   ! atm present

   !----------------------------------------------------------
   !| Read driver restart file, overwrite anything previously sent or computed
   !----------------------------------------------------------

   call t_startf('CPL:init_readrestart')
   call t_adj_detailf(+2)

   call seq_diag_zero_mct(mode='all')
   if (read_restart .and. iamin_CPLID) then
      call seq_rest_read(rest_file, infodata, &
           atm, lnd, ice, ocn, rof, glc, wav, esp, &
           fractions_ax, fractions_lx, fractions_ix, fractions_ox, &
           fractions_rx, fractions_gx, fractions_wx)
   endif

   call t_adj_detailf(-2)
   call t_stopf  ('CPL:init_readrestart')

   !----------------------------------------------------------
   !| Map initial r2x_rx and g2x_gx to _ox, _ix and _lx
   !----------------------------------------------------------

   if (iamin_CPLID ) then
      if (rof_c2_ocn) then
         call prep_ocn_calc_r2x_ox(timer='CPL:init_rof2ocn')
      endif
      if (glc_c2_ocn) then
         call prep_ocn_calc_g2x_ox(timer='CPL:init_glc2ocn')
      endif
      if (rof_c2_ice) then
         call prep_ice_calc_r2x_ix(timer='CPL:init_rof2ice')
      endif
      if (glc_c2_ice) then
         call prep_ice_calc_g2x_ix(timer='CPL:init_glc2ice')
      endif
      if (rof_c2_lnd) then
         call prep_lnd_calc_r2x_lx(timer='CPL:init_rof2lnd')
      endif
      if (glc_c2_lnd) then
         call prep_lnd_calc_g2x_lx(timer='CPL:init_gllndnd')
      endif
   endif

   !----------------------------------------------------------
   !| Write histinit output file
   !----------------------------------------------------------

   if (do_histinit) then
      if (iamin_CPLID) then
         call t_startf('CPL:init_histinit')
         call t_adj_detailf(+2)

         if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
         if (iamroot_CPLID) then
            call seq_timemgr_EClockGetData( EClock_d, curr_ymd=ymd, curr_tod=tod )
            write(logunit,104) ' Write history file at ',ymd,tod
            call shr_sys_flush(logunit)
         endif
         call seq_hist_write(infodata, EClock_d, &
              atm, lnd, ice, ocn, rof, glc, wav, &
              fractions_ax, fractions_lx, fractions_ix, fractions_ox, &
              fractions_rx, fractions_gx, fractions_wx)
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)

         call t_adj_detailf(-2)
         call t_stopf('CPL:init_histinit')
      endif
   endif

   if (iamroot_CPLID )then
      write(logunit,*) ' '
      write(logunit,F00) 'Model initialization complete '
      write(logunit,*) ' '
      call shr_sys_flush(logunit)
   endif

   call t_adj_detailf(-1)
   call t_stopf('cesm_init')

end subroutine cesm_init

 !===============================================================================
 !*******************************************************************************
 !===============================================================================

 subroutine cesm_run()
   use seq_comm_mct, only : atm_layout, lnd_layout, ice_layout, glc_layout, rof_layout, &
         ocn_layout, wav_layout, esp_layout

   implicit none
   ! gptl timer lookup variables
   integer, parameter :: hashcnt=7
   integer :: hashint(hashcnt)

101 format( A, 2i8, 12A, A, F8.2, A, F8.2 )
102 format( A, 2i8, A, 8L3 )
103 format( 5A )
104 format( A, 2i8)
105 format( A, 2i8, A, f10.2, A, f10.2, A, A, i5, A, A)
106 format( A, f23.12)
107 format( A, 2i8, A, f12.4, A, f12.4 )
108 format( A, f10.2, A, i8.8)
109 format( A, 2f10.3)
110 format( A, 2i8, A, 9L3 )


   hashint = 0


   call seq_infodata_putData(infodata,atm_phase=1,lnd_phase=1,ocn_phase=1,ice_phase=1)
   call seq_timemgr_EClockGetData( EClock_d, stepno=begstep)
   call seq_timemgr_EClockGetData( EClock_d, dtime=dtime)
   call seq_timemgr_EClockGetData( EClock_d, calendar=calendar)
   ncpl = 86400/dtime
   cktime_acc = 0._r8
   cktime_cnt = 0
   stop_alarm = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_stop)
   if (seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_datestop)) then
      if (iamroot_CPLID) then
         write(logunit,*) ' '
         write(logunit,103) subname,' NOTE: Stopping from alarm STOP DATE'
         write(logunit,*) ' '
      endif
      stop_alarm = .true.
   endif
   force_stop = .false.
   force_stop_ymd = -1
   force_stop_tod = -1

   !|----------------------------------------------------------
   !| Beginning of driver time step loop
   !|----------------------------------------------------------

   call t_startf ('CPL:RUN_LOOP_BSTART')
   call mpi_barrier(mpicom_GLOID,ierr)
   call t_stopf ('CPL:RUN_LOOP_BSTART')
   Time_begin = mpi_wtime()
   Time_bstep = mpi_wtime()
   do while ( .not. stop_alarm)

      call t_startf('CPL:RUN_LOOP', hashint(1))
      call t_startf('CPL:CLOCK_ADVANCE')

      !----------------------------------------------------------
      !| Advance Clock
      !  (this is time that models should have before they return
      !  to the driver).  Write timestamp and run alarm status
      !----------------------------------------------------------

      call seq_timemgr_clockAdvance( seq_SyncClock, force_stop, force_stop_ymd, force_stop_tod)
      call seq_timemgr_EClockGetData( EClock_d, curr_ymd=ymd, curr_tod=tod )
      call shr_cal_date2ymd(ymd,year,month,day)
      stop_alarm    = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_stop)
      atmrun_alarm  = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_atmrun)
      lndrun_alarm  = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_lndrun)
      rofrun_alarm  = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_rofrun)
      icerun_alarm  = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_icerun)
      glcrun_alarm  = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_glcrun)
      wavrun_alarm  = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_wavrun)
      esprun_alarm  = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_esprun)
      ocnrun_alarm  = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_ocnrun)
      ocnnext_alarm = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_ocnnext)
      restart_alarm = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_restart)
      history_alarm = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_history)
      histavg_alarm = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_histavg)
      tprof_alarm   = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_tprof)
      barrier_alarm = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_barrier)

      ! this probably belongs in seq_timemgr somewhere using proper clocks
      t1hr_alarm = .false.
      t2hr_alarm = .false.
      t3hr_alarm = .false.
      t6hr_alarm = .false.
      t12hr_alarm = .false.
      t24hr_alarm = .false.
      t1yr_alarm = .false.
      if (mod(tod, 3600) == 0) t1hr_alarm = .true.
      if (mod(tod, 7200) == 0) t2hr_alarm = .true.
      if (mod(tod,10800) == 0) t3hr_alarm = .true.
      if (mod(tod,21600) == 0) t6hr_alarm = .true.
      if (mod(tod,43200) == 0) t12hr_alarm = .true.
      if (tod            == 0) t24hr_alarm = .true.
      if (month==1 .and. day==1 .and. tod==0) t1yr_alarm = .true.

      call seq_infodata_putData(infodata, glcrun_alarm=glcrun_alarm)

      if (seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_datestop)) then
         if (iamroot_CPLID) then
            write(logunit,*) ' '
            write(logunit,103) subname,' NOTE: Stopping from alarm STOP DATE'
            write(logunit,*) ' '
         endif
         stop_alarm = .true.
      endif

      ! update the orbital data as needed
      if (trim(orb_mode) == trim(seq_infodata_orb_variable_year)) then
         orb_nyear =  orb_iyear + (year - orb_iyear_align)
         if (orb_nyear /= orb_cyear) then
            orb_cyear = orb_nyear
            call shr_orb_params(orb_cyear, orb_eccen, orb_obliq, orb_mvelp, &
                 orb_obliqr, orb_lambm0, orb_mvelpp, iamroot_CPLID)
            call seq_infodata_putData(infodata,orb_eccen=orb_eccen,orb_obliqr=orb_obliqr, &
                 orb_lambm0=orb_lambm0,orb_mvelpp=orb_mvelpp)
         endif
      endif

      ! override ocnrun_alarm and ocnnext_alarm for first ocn run
      ! skip_ocean_run is initialized above to true if it's a startup
      ! if it's not a startup, ignore all of this
      ! stop the overide on the second ocnrun_alarm

      if (ocnrun_alarm) ocnrun_count = ocnrun_count + 1
      if (ocnrun_count > 1) skip_ocean_run = .false.
      if (skip_ocean_run) then
         ocnrun_alarm = .false.
         ocnnext_alarm = .false.
      endif

      if (iamroot_CPLID) then
         if (loglevel > 1) then
            write(logunit,102) ' Alarm_state: model date = ',ymd,tod, &
                 ' aliogrw run alarms = ',  atmrun_alarm, lndrun_alarm, &
                 icerun_alarm, ocnrun_alarm, glcrun_alarm, &
                 rofrun_alarm, wavrun_alarm, esprun_alarm
            write(logunit,102) ' Alarm_state: model date = ',ymd,tod, &
                 ' 1.2.3.6.12.24 run alarms = ',  t1hr_alarm, t2hr_alarm, &
                 t3hr_alarm, t6hr_alarm, t12hr_alarm, t24hr_alarm
            call shr_sys_flush(logunit)
         endif
      endif

      call t_stopf ('CPL:CLOCK_ADVANCE')

      !----------------------------------------------------------
      !| MAP ATM to OCN
      !  Set a2x_ox as a module variable in prep_ocn_mod
      !  This will be used later in the ice prep and in the
      !  atm/ocn flux calculation
      !----------------------------------------------------------

      if (iamin_CPLID .and. (atm_c2_ocn .or. atm_c2_ice)) then
         call cesm_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:OCNPRE1_BARRIER')
         call t_drvstartf ('CPL:OCNPRE1',cplrun=.true.,barrier=mpicom_CPLID,hashint=hashint(3))
         if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)

         call prep_ocn_calc_a2x_ox(timer='CPL:ocnpre1_atm2ocn')

         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
         call t_drvstopf  ('CPL:OCNPRE1',cplrun=.true.,hashint=hashint(3))
      endif

      !----------------------------------------------------------
      !| ATM/OCN SETUP (rasm_option1)
      !----------------------------------------------------------

      if ((trim(cpl_seq_option) == 'RASM_OPTION1') .and. &
           iamin_CPLID .and. ocn_present) then

         call cesm_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:ATMOCN1_BARRIER')
         call t_drvstartf ('CPL:ATMOCN1',cplrun=.true.,barrier=mpicom_CPLID,hashint=hashint(4))
         if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)

         if (ocn_prognostic) then
            ! Map ice to ocn
            if (ice_c2_ocn) call prep_ocn_calc_i2x_ox(timer='CPL:atmocnp_ice2ocn')

            ! Map wav to ocn
            if (wav_c2_ocn) call prep_ocn_calc_w2x_ox(timer='CPL:atmocnp_wav2ocn')
         endif

         !----------------------------------------------------------
         !| atm/ocn flux on atm grid (rasm_option1 and aoflux='atm')
         !----------------------------------------------------------

         if (trim(aoflux_grid) == 'atm') then
            ! compute o2x_ax for flux_atmocn, will be updated before atm merge
            ! do not use fractions because fractions here are NOT consistent with fractions in atm_mrg
            if (ocn_c2_atm) call prep_atm_calc_o2x_ax(timer='CPL:atmoca_ocn2atm')

            call t_drvstartf ('CPL:atmocna_fluxa',barrier=mpicom_CPLID)
            do exi = 1,num_inst_xao
               eai = mod((exi-1),num_inst_atm) + 1
               eoi = mod((exi-1),num_inst_ocn) + 1
               efi = mod((exi-1),num_inst_frc) + 1
               a2x_ax => component_get_c2x_cx(atm(eai))
               o2x_ax => prep_atm_get_o2x_ax()    ! array over all instances
               xao_ax => prep_aoflux_get_xao_ax() ! array over all instances
               call seq_flux_atmocn_mct(infodata, tod, dtime, a2x_ax, o2x_ax(eoi), xao_ax(exi))
            enddo
            call t_drvstopf  ('CPL:atmocna_fluxa')

            if (atm_c2_ocn) call prep_aoflux_calc_xao_ox(timer='CPL:atmocna_atm2ocn')
         endif  ! aoflux_grid

         !----------------------------------------------------------
         !| atm/ocn flux on ocn grid (rasm_option1 and aoflux='ocn')
         !----------------------------------------------------------

         if (trim(aoflux_grid) == 'ocn') then
            call t_drvstartf ('CPL:atmocnp_fluxo',barrier=mpicom_CPLID,hashint=hashint(6))
            do exi = 1,num_inst_xao
               eai = mod((exi-1),num_inst_atm) + 1
               eoi = mod((exi-1),num_inst_ocn) + 1
               efi = mod((exi-1),num_inst_frc) + 1
               a2x_ox => prep_ocn_get_a2x_ox()
               o2x_ox => component_get_c2x_cx(ocn(eoi))
               xao_ox => prep_aoflux_get_xao_ox()
               call seq_flux_atmocn_mct(infodata, tod, dtime, a2x_ox(eai), o2x_ox, xao_ox(exi))
            enddo
            call t_drvstopf  ('CPL:atmocnp_fluxo',hashint=hashint(6))
         endif

         !----------------------------------------------------------
         !| ocn prep-merge (rasm_option1)
         !----------------------------------------------------------

         xao_ox => prep_aoflux_get_xao_ox()
         call prep_ocn_mrg(infodata, fractions_ox, xao_ox=xao_ox, timer_mrg='CPL:atmocnp_mrgx2o')

         ! Accumulate ocn inputs - form partial sum of tavg ocn inputs (virtual "send" to ocn)
         call prep_ocn_accum(timer='CPL:atmocnp_accum')

         !----------------------------------------------------------
         !| ocn albedos (rasm_option1)
         !  (MUST BE AFTER prep_ocn_mrg for swnet to ocn to be computed properly
         !----------------------------------------------------------

         call t_drvstartf ('CPL:atmocnp_ocnalb', barrier=mpicom_CPLID,hashint=hashint(5))
         do exi = 1,num_inst_xao
            efi = mod((exi-1),num_inst_frc) + 1
            eai = mod((exi-1),num_inst_atm) + 1
            xao_ox => prep_aoflux_get_xao_ox()        ! array over all instances
            a2x_ox => prep_ocn_get_a2x_ox()
            call seq_flux_ocnalb_mct(infodata, ocn(1), a2x_ox(eai), fractions_ox(efi), xao_ox(exi))
         enddo
         call t_drvstopf  ('CPL:atmocnp_ocnalb',hashint=hashint(5))

         !----------------------------------------------------------
         !| ocn budget (rasm_option1)
         !----------------------------------------------------------

         if (do_budgets) then
            call cesm_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:BUDGET0_BARRIER')
            call t_drvstartf ('CPL:BUDGET0',cplrun=.true.,budget=.true.,barrier=mpicom_CPLID)
            xao_ox => prep_aoflux_get_xao_ox() ! array over all instances
            call seq_diag_ocn_mct(ocn(ens1), xao_ox(1), fractions_ox(ens1), infodata, &
                 do_o2x=.true., do_x2o=.true., do_xao=.true.)
            call t_drvstopf ('CPL:BUDGET0',cplrun=.true.,budget=.true.)
         endif

         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
         call t_drvstopf  ('CPL:ATMOCN1',cplrun=.true.,hashint=hashint(4))
      endif

      !----------------------------------------------------------
      !| ATM/OCN SETUP-SEND (cesm1_orig, cesm1_orig_tight, cesm1_mod, cesm1_mod_tight, or rasm_option1)
      !----------------------------------------------------------

      if ((trim(cpl_seq_option) == 'CESM1_ORIG' .or. &
           trim(cpl_seq_option) == 'CESM1_ORIG_TIGHT' .or. &
           trim(cpl_seq_option) == 'CESM1_MOD'  .or. &
           trim(cpl_seq_option) == 'CESM1_MOD_TIGHT'  .or. &
           trim(cpl_seq_option) == 'RASM_OPTION1'  ) .and. &
           ocn_present .and. ocnrun_alarm) then

         !----------------------------------------------------
         ! "startup" wait (cesm1_orig, cesm1_mod, or rasm_option1)
         !----------------------------------------------------

         if (iamin_CPLALLOCNID) then
            ! want to know the time the ocean pes waited for the cpl pes
            ! at the first ocnrun_alarm, min ocean wait is wait time
            ! do not use t_barrierf here since it can be "off", use mpi_barrier
            do eoi = 1,num_inst_ocn
               if (ocn(eoi)%iamin_compid) call t_drvstartf ('CPL:C2O_INITWAIT')
            enddo
            call mpi_barrier(mpicom_CPLALLOCNID,ierr)
            do eoi = 1,num_inst_ocn
               if (ocn(eoi)%iamin_compid) call t_drvstopf  ('CPL:C2O_INITWAIT')
            enddo
            cpl2ocn_first = .false.
         endif

         !----------------------------------------------------
         !| ocn average (cesm1_orig, cesm1_orig_tight, cesm1_mod, cesm1_mod_tight, or rasm_option1)
         !----------------------------------------------------

         if (iamin_CPLID .and. ocn_prognostic) then
            call cesm_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:OCNPREP_BARRIER')
            call t_drvstartf ('CPL:OCNPREP',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)

            ! finish accumulating ocean inputs
            ! reset the value of x2o_ox with the value in x2oacc_ox
            ! (module variable in prep_ocn_mod)
            call prep_ocn_accum_avg(timer_accum='CPL:ocnprep_avg')

            call component_diag(infodata, ocn, flow='x2c', comment= 'send ocn', &
                 info_debug=info_debug, timer_diag='CPL:ocnprep_diagav')

            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('CPL:OCNPREP',cplrun=.true.)
         endif

         !----------------------------------------------------
         !| cpl -> ocn (cesm1_orig, cesm1_orig_tight, cesm1_mod, cesm1_mod_tight, or rasm_option1)
         !----------------------------------------------------

         if (iamin_CPLALLOCNID .and. ocn_prognostic) then
            call component_exch(ocn, flow='x2c', &
                 infodata=infodata, infodata_string='cpl2ocn_run', &
                 mpicom_barrier=mpicom_CPLALLOCNID, run_barriers=run_barriers, &
                 timer_barrier='CPL:C2O_BARRIER', timer_comp_exch='CPL:C2O', &
                 timer_map_exch='CPL:c2o_ocnx2ocno', timer_infodata_exch='CPL:c2o_infoexch')
         endif

      endif ! end of OCN SETUP

      !----------------------------------------------------------
      !| LND SETUP-SEND
      !----------------------------------------------------------

      if (lnd_present .and. lndrun_alarm) then

         !----------------------------------------------------
         !| lnd prep-merge
         !----------------------------------------------------

         if (iamin_CPLID) then
            call cesm_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:LNDPREP_BARRIER')
            call t_drvstartf ('CPL:LNDPREP',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)

            if (atm_c2_lnd) then
               call prep_lnd_calc_a2x_lx(timer='CPL:lndprep_atm2lnd')
            endif

            if (lnd_prognostic) then
               call prep_lnd_mrg(infodata, timer_mrg='CPL:lndprep_mrgx2l')

               call component_diag(infodata, lnd, flow='x2c', comment= 'send lnd', &
                    info_debug=info_debug, timer_diag='CPL:lndprep_diagav')
            endif

            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('CPL:LNDPREP',cplrun=.true.)
         endif

         !----------------------------------------------------
         !| cpl -> lnd
         !----------------------------------------------------

         if (iamin_CPLALLLNDID) then
            call component_exch(lnd, flow='x2c', &
                 infodata=infodata, infodata_string='cpl2lnd_run', &
                 mpicom_barrier=mpicom_CPLALLLNDID, run_barriers=run_barriers, &
                 timer_barrier='CPL:C2L_BARRIER', timer_comp_exch='CPL:C2L', &
                 timer_map_exch='CPL:c2l_lndx2lndl', timer_infodata_exch='CPL:c2l_infoexch')
         endif

      endif

      !----------------------------------------------------------
      !| ICE SETUP-SEND
      !  Note that for atm->ice mapping below will leverage the assumption that the
      !  ice and ocn are on the same grid and that mapping of atm to ocean is
      !  done already for use by atmocn flux and ice model prep
      !----------------------------------------------------------

      if (ice_present .and. icerun_alarm) then

         !----------------------------------------------------
         !| ice prep-merge
         !----------------------------------------------------

         if (iamin_CPLID .and. ice_prognostic) then
            call cesm_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:ICEPREP_BARRIER')

            call t_drvstartf ('CPL:ICEPREP',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)


            if (ocn_c2_ice) then
               call prep_ice_calc_o2x_ix(timer='CPL:iceprep_ocn2ice')
            endif

            if (atm_c2_ice) then
               ! This is special to avoid remapping atm to ocn
               ! Note it is constrained that different prep modules cannot
               ! use or call each other
               a2x_ox => prep_ocn_get_a2x_ox() ! array
               call prep_ice_calc_a2x_ix(a2x_ox, timer='CPL:iceprep_atm2ice')
            endif

            call prep_ice_mrg(infodata, timer_mrg='CPL:iceprep_mrgx2i')

            call component_diag(infodata, ice, flow='x2c', comment= 'send ice', &
                 info_debug=info_debug, timer_diag='CPL:iceprep_diagav')

            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('CPL:ICEPREP',cplrun=.true.)
         endif

         !----------------------------------------------------
         !| cpl -> ice
         !----------------------------------------------------

         if (iamin_CPLALLICEID .and. ice_prognostic) then
            call component_exch(ice, flow='x2c', &
                 infodata=infodata, infodata_string='cpl2ice_run', &
                 mpicom_barrier=mpicom_CPLALLICEID, run_barriers=run_barriers, &
                 timer_barrier='CPL:C2I_BARRIER', timer_comp_exch='CPL:C2I', &
                 timer_map_exch='CPL:c2i_icex2icei', timer_infodata_exch='CPL:ice_infoexch')
         endif

      endif

      !----------------------------------------------------------
      !| WAV SETUP-SEND
      !----------------------------------------------------------

      if (wav_present .and. wavrun_alarm) then

         !----------------------------------------------------------
         !| wav prep-merge
         !----------------------------------------------------------

         if (iamin_CPLID .and. wav_prognostic) then
            call cesm_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:WAVPREP_BARRIER')

            call t_drvstartf ('CPL:WAVPREP',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)

            if (atm_c2_wav) then
               call prep_wav_calc_a2x_wx(timer='CPL:wavprep_atm2wav')
            endif

            if (ocn_c2_wav) then
               call prep_wav_calc_o2x_wx(timer='CPL:wavprep_ocn2wav')
            endif

            if (ice_c2_wav) then
               call prep_wav_calc_i2x_wx(timer='CPL:wavprep_ice2wav')
            endif

            call prep_wav_mrg(infodata, fractions_wx, timer_mrg='CPL:wavprep_mrgx2w')

            call component_diag(infodata, wav, flow='x2c', comment= 'send wav', &
                 info_debug=info_debug, timer_diag='CPL:wavprep_diagav')

            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('CPL:WAVPREP',cplrun=.true.)
         endif

         !----------------------------------------------------------
         !| cpl -> wav
         !----------------------------------------------------------

         if (iamin_CPLALLWAVID .and. wav_prognostic) then
            call component_exch(wav, flow='x2c', &
                 infodata=infodata, infodata_string='cpl2wav_run', &
                 mpicom_barrier=mpicom_CPLALLWAVID, run_barriers=run_barriers, &
                 timer_barrier='CPL:C2W_BARRIER', timer_comp_exch='CPL:C2W', &
                 timer_map_exch='CPL:c2w_wavx2wavw', timer_infodata_exch='CPL:c2w_infoexch')
         endif

      endif

      !----------------------------------------------------------
      !| ROF SETUP-SEND
      !----------------------------------------------------------

      if (rof_present .and. rofrun_alarm) then

         !----------------------------------------------------
         !| rof prep-merge
         !----------------------------------------------------

         if (iamin_CPLID .and. rof_prognostic) then
            call cesm_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:ROFPREP_BARRIER')

            call t_drvstartf ('CPL:ROFPREP', cplrun=.true., barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)

            call prep_rof_accum_avg(timer='CPL:rofprep_l2xavg')

            if (lnd_c2_rof) then
               call prep_rof_calc_l2r_rx(fractions_lx, timer='CPL:rofprep_lnd2rof')
            endif

            call prep_rof_mrg(infodata, fractions_rx, timer_mrg='CPL:rofprep_mrgx2r')

            call component_diag(infodata, rof, flow='x2c', comment= 'send rof', &
                 info_debug=info_debug, timer_diag='CPL:rofprep_diagav')

            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('CPL:ROFPREP',cplrun=.true.)
         endif

         !----------------------------------------------------
         !| cpl -> rof
         !----------------------------------------------------

         if (iamin_CPLALLROFID .and. rof_prognostic) then
            call component_exch(rof, flow='x2c', &
                 infodata=infodata, infodata_string='cpl2rof_run', &
                 mpicom_barrier=mpicom_CPLALLLNDID, run_barriers=run_barriers, &
                 timer_barrier='CPL:C2R_BARRIER', timer_comp_exch='CPL:C2R', &
                 timer_map_exch='CPL:c2r_rofx2rofr', timer_infodata_exch='CPL:c2r_infoexch')
         endif

      endif

      !----------------------------------------------------------
      !| RUN ICE MODEL
      !----------------------------------------------------------

      if (ice_present .and. icerun_alarm) then
         call component_run(Eclock_i, ice, ice_run, infodata, &
              seq_flds_x2c_fluxes=seq_flds_x2i_fluxes, &
              seq_flds_c2x_fluxes=seq_flds_i2x_fluxes, &
              comp_prognostic=ice_prognostic, comp_num=comp_num_ice, &
              timer_barrier= 'CPL:ICE_RUN_BARRIER', timer_comp_run='CPL:ICE_RUN', &
              run_barriers=run_barriers, ymd=ymd, tod=tod,comp_layout=ice_layout)
      endif

      !----------------------------------------------------------
      !| RUN LND MODEL
      !----------------------------------------------------------

      if (lnd_present .and. lndrun_alarm) then
         call component_run(Eclock_l, lnd, lnd_run, infodata, &
              seq_flds_x2c_fluxes=seq_flds_x2l_fluxes, &
              seq_flds_c2x_fluxes=seq_flds_l2x_fluxes, &
              comp_prognostic=lnd_prognostic, comp_num=comp_num_lnd, &
              timer_barrier= 'CPL:LND_RUN_BARRIER', timer_comp_run='CPL:LND_RUN', &
              run_barriers=run_barriers, ymd=ymd, tod=tod,comp_layout=lnd_layout)
      endif

      !----------------------------------------------------------
      !| RUN ROF MODEL
      !----------------------------------------------------------

      if (rof_present .and. rofrun_alarm) then
         call component_run(Eclock_r, rof, rof_run, infodata, &
              seq_flds_x2c_fluxes=seq_flds_x2r_fluxes, &
              seq_flds_c2x_fluxes=seq_flds_r2x_fluxes, &
              comp_prognostic=rof_prognostic, comp_num=comp_num_rof, &
              timer_barrier= 'CPL:ROF_RUN_BARRIER', timer_comp_run='CPL:ROF_RUN', &
              run_barriers=run_barriers, ymd=ymd, tod=tod,comp_layout=rof_layout)
      endif

      !----------------------------------------------------------
      !| RUN WAV MODEL
      !----------------------------------------------------------

      if (wav_present .and. wavrun_alarm) then
         call component_run(Eclock_w, wav, wav_run, infodata, &
              seq_flds_x2c_fluxes=seq_flds_x2w_fluxes, &
              seq_flds_c2x_fluxes=seq_flds_w2x_fluxes, &
              comp_prognostic=wav_prognostic, comp_num=comp_num_wav, &
              timer_barrier= 'CPL:WAV_RUN_BARRIER', timer_comp_run='CPL:WAV_RUN', &
              run_barriers=run_barriers, ymd=ymd, tod=tod,comp_layout=wav_layout)
      endif

      !----------------------------------------------------------
      !| RUN OCN MODEL (cesm1_orig_tight or cesm1_mod_tight)
      !----------------------------------------------------------

      if ((trim(cpl_seq_option) == 'CESM1_ORIG_TIGHT' .or. &
           trim(cpl_seq_option) == 'CESM1_MOD_TIGHT'   ) .and. &
          ocn_present .and. ocnrun_alarm) then
         call component_run(Eclock_o, ocn, ocn_run, infodata, &
              seq_flds_x2c_fluxes=seq_flds_x2o_fluxes, &
              seq_flds_c2x_fluxes=seq_flds_o2x_fluxes, &
              comp_prognostic=ocn_prognostic, comp_num=comp_num_ocn, &
              timer_barrier= 'CPL:OCNT_RUN_BARRIER', timer_comp_run='CPL:OCNT_RUN', &
              run_barriers=run_barriers, ymd=ymd, tod=tod,comp_layout=ocn_layout)
      endif

      !----------------------------------------------------------
      !| OCN RECV-POST (cesm1_orig_tight or cesm1_mod_tight)
      !----------------------------------------------------------

      if ((trim(cpl_seq_option) == 'CESM1_ORIG_TIGHT' .or. &
           trim(cpl_seq_option) == 'CESM1_MOD_TIGHT'   ) .and. &
          ocn_present .and. ocnnext_alarm) then

         !----------------------------------------------------------
         !| ocn -> cpl (cesm1_orig_tight or cesm1_mod_tight)
         !----------------------------------------------------------

         if (iamin_CPLALLOCNID) then
            call component_exch(ocn, flow='c2x', &
                 infodata=infodata, infodata_string='ocn2cpl_run', &
                 mpicom_barrier=mpicom_CPLALLOCNID, run_barriers=run_barriers, &
                 timer_barrier='CPL:O2CT_BARRIER', timer_comp_exch='CPL:O2CT', &
                 timer_map_exch='CPL:o2c_ocno2ocnx', timer_infodata_exch='CPL:o2c_infoexch')
         endif

         !----------------------------------------------------------
         !| ocn post (cesm1_orig_tight or cesm1_mod_tight)
         !----------------------------------------------------------

         if (iamin_CPLID) then
            call cesm_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:OCNPOSTT_BARRIER')
            call t_drvstartf  ('CPL:OCNPOSTT',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)

            call component_diag(infodata, ocn, flow='c2x', comment= 'recv ocn', &
                 info_debug=info_debug, timer_diag='CPL:ocnpost_diagav')

            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('CPL:OCNPOSTT',cplrun=.true.)
        endif

      endif

      !----------------------------------------------------------
      !| ATM/OCN SETUP (cesm1_orig, cesm1_orig_tight, cesm1_mod or cesm1_mod_tight)
      !----------------------------------------------------------

      if ((trim(cpl_seq_option) == 'CESM1_ORIG'       .or. &
           trim(cpl_seq_option) == 'CESM1_ORIG_TIGHT' .or. &
           trim(cpl_seq_option) == 'CESM1_MOD'        .or. &
           trim(cpl_seq_option) == 'CESM1_MOD_TIGHT' ) .and. &
           iamin_CPLID .and. ocn_present) then

         call cesm_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:ATMOCNP_BARRIER')
         call t_drvstartf ('CPL:ATMOCNP',cplrun=.true.,barrier=mpicom_CPLID,hashint=hashint(7))
         if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)

         !----------------------------------------------------------
         !| ocn prep-merge (cesm1_orig or cesm1_orig_tight)
         !----------------------------------------------------------

         if (ocn_prognostic) then
            ! Map ice to ocn
            if (ice_c2_ocn) call prep_ocn_calc_i2x_ox(timer='CPL:atmocnp_ice2ocn')

            ! Map wav to ocn
            if (wav_c2_ocn) call prep_ocn_calc_w2x_ox(timer='CPL:atmocnp_wav2ocn')

            if (cpl_seq_option == 'CESM1_ORIG' .or. &
                cpl_seq_option == 'CESM1_ORIG_TIGHT') then
               xao_ox => prep_aoflux_get_xao_ox()
               call prep_ocn_mrg(infodata, fractions_ox, xao_ox=xao_ox, timer_mrg='CPL:atmocnp_mrgx2o')

               ! Accumulate ocn inputs - form partial sum of tavg ocn inputs (virtual "send" to ocn)
               call prep_ocn_accum(timer='CPL:atmocnp_accum')
            endif
         endif

         !----------------------------------------------------------
         !| atm/ocn flux on atm grid ((cesm1_orig, cesm1_orig_tight, cesm1_mod or cesm1_mod_tight) and aoflux='atm')
         !----------------------------------------------------------

         if (trim(aoflux_grid) == 'atm') then
            ! compute o2x_ax for flux_atmocn, will be updated before atm merge
            ! do not use fractions because fractions here are NOT consistent with fractions in atm_mrg
            if (ocn_c2_atm) call prep_atm_calc_o2x_ax(timer='CPL:atmoca_ocn2atm')

            call t_drvstartf ('CPL:atmocna_fluxa',barrier=mpicom_CPLID)
            do exi = 1,num_inst_xao
               eai = mod((exi-1),num_inst_atm) + 1
               eoi = mod((exi-1),num_inst_ocn) + 1
               efi = mod((exi-1),num_inst_frc) + 1
               a2x_ax => component_get_c2x_cx(atm(eai))
               o2x_ax => prep_atm_get_o2x_ax()    ! array over all instances
               xao_ax => prep_aoflux_get_xao_ax() ! array over all instances
               call seq_flux_atmocn_mct(infodata, tod, dtime, a2x_ax, o2x_ax(eoi), xao_ax(exi))
            enddo
            call t_drvstopf  ('CPL:atmocna_fluxa')

            if (atm_c2_ocn) call prep_aoflux_calc_xao_ox(timer='CPL:atmocna_atm2ocn')
         endif  ! aoflux_grid

         !----------------------------------------------------------
         !| atm/ocn flux on ocn grid ((cesm1_orig, cesm1_orig_tight, cesm1_mod or cesm1_mod_tight) and aoflux='ocn')
         !----------------------------------------------------------

         if (trim(aoflux_grid) == 'ocn') then
            call t_drvstartf ('CPL:atmocnp_fluxo',barrier=mpicom_CPLID)
            do exi = 1,num_inst_xao
               eai = mod((exi-1),num_inst_atm) + 1
               eoi = mod((exi-1),num_inst_ocn) + 1
               efi = mod((exi-1),num_inst_frc) + 1
               a2x_ox => prep_ocn_get_a2x_ox()
               o2x_ox => component_get_c2x_cx(ocn(eoi))
               xao_ox => prep_aoflux_get_xao_ox()
               call seq_flux_atmocn_mct(infodata, tod, dtime, a2x_ox(eai), o2x_ox, xao_ox(exi))
            enddo
            call t_drvstopf  ('CPL:atmocnp_fluxo')
!         else if (trim(aoflux_grid) == 'atm') then
!            !--- compute later ---
!
!         else if (trim(aoflux_grid) == 'exch') then
!            xao_ax   => prep_aoflux_get_xao_ax()
!            xao_ox   => prep_aoflux_get_xao_ox()
!
!            call t_drvstartf ('CPL:atmocnp_fluxe',barrier=mpicom_CPLID)
!            call seq_flux_atmocnexch_mct( infodata, atm(eai), ocn(eoi), &
!                 fractions_ax(efi), fractions_ox(efi), xao_ax(exi), xao_ox(exi) )
!            call t_drvstopf  ('CPL:atmocnp_fluxe')
         endif  ! aoflux_grid

         !----------------------------------------------------------
         !| ocn prep-merge (cesm1_mod or cesm1_mod_tight)
         !----------------------------------------------------------

         if (ocn_prognostic) then
            if (cpl_seq_option == 'CESM1_MOD' .or. &
                cpl_seq_option == 'CESM1_MOD_TIGHT') then

               xao_ox => prep_aoflux_get_xao_ox()
               call prep_ocn_mrg(infodata, fractions_ox, xao_ox=xao_ox, timer_mrg='CPL:atmocnp_mrgx2o')

               ! Accumulate ocn inputs - form partial sum of tavg ocn inputs (virtual "send" to ocn)
               call prep_ocn_accum(timer='CPL:atmocnp_accum')
            endif
         endif

         !----------------------------------------------------------
         !| ocn albedos (cesm1_orig, cesm1_orig_tight, cesm1_mod or cesm1_mod_tight)
         !  (MUST BE AFTER prep_ocn_mrg for swnet to ocn to be computed properly
         !----------------------------------------------------------

         call t_drvstartf ('CPL:atmocnp_ocnalb', barrier=mpicom_CPLID)
         do exi = 1,num_inst_xao
            efi = mod((exi-1),num_inst_frc) + 1
            eai = mod((exi-1),num_inst_atm) + 1
            xao_ox => prep_aoflux_get_xao_ox()        ! array over all instances
            a2x_ox => prep_ocn_get_a2x_ox()
            call seq_flux_ocnalb_mct(infodata, ocn(1), a2x_ox(eai), fractions_ox(efi), xao_ox(exi))
         enddo
         call t_drvstopf  ('CPL:atmocnp_ocnalb')

         !----------------------------------------------------------
         !| ocn budget (cesm1_orig, cesm1_orig_tight, cesm1_mod or cesm1_mod_tight)
         !----------------------------------------------------------

         if (do_budgets) then
            call cesm_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:BUDGET0_BARRIER')
            call t_drvstartf ('CPL:BUDGET0',cplrun=.true.,budget=.true.,barrier=mpicom_CPLID)
            xao_ox => prep_aoflux_get_xao_ox() ! array over all instances
            call seq_diag_ocn_mct(ocn(ens1), xao_ox(1), fractions_ox(ens1), infodata, &
                 do_o2x=.true., do_x2o=.true., do_xao=.true.)
            call t_drvstopf ('CPL:BUDGET0',cplrun=.true.,budget=.true.)
         endif

         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
         call t_drvstopf  ('CPL:ATMOCNP',cplrun=.true.,hashint=hashint(7))
      endif

      !----------------------------------------------------------
      !| LND RECV-POST
      !----------------------------------------------------------

      if (lnd_present .and. lndrun_alarm) then

         !----------------------------------------------------------
         !| lnd -> cpl
         !----------------------------------------------------------

         if (iamin_CPLALLLNDID) then
            call component_exch(lnd, flow='c2x', infodata=infodata, infodata_string='lnd2cpl_run', &
                 mpicom_barrier=mpicom_CPLALLLNDID, run_barriers=run_barriers, &
                 timer_barrier='CPL:L2C_BARRIER', timer_comp_exch='CPL:L2C', &
                 timer_map_exch='CPL:l2c_lndl2lndx', timer_infodata_exch='lnd2cpl_run')
         endif

         !----------------------------------------------------------
         !| lnd post
         !----------------------------------------------------------

         if (iamin_CPLID) then
            call cesm_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:LNDPOST_BARRIER')
            call t_drvstartf  ('CPL:LNDPOST',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)

            call component_diag(infodata, lnd, flow='c2x', comment='recv lnd', &
                 info_debug=info_debug, timer_diag='CPL:lndpost_diagav')

            ! Accumulate rof and glc inputs (module variables in prep_rof_mod and prep_glc_mod)
            if (lnd_c2_rof) then
               call prep_rof_accum(timer='CPL:lndpost_accl2r')
            endif
            if (lnd_c2_glc) then
               call prep_glc_accum(timer='CPL:lndpost_accl2g' )
            endif

            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('CPL:LNDPOST',cplrun=.true.)
         endif
      endif

      !----------------------------------------------------------
      !| GLC SETUP-SEND
      !----------------------------------------------------------

      if (glc_present .and. glcrun_alarm) then

         !----------------------------------------------------
         !| glc prep-merge
         !----------------------------------------------------

         if (iamin_CPLID .and. glc_prognostic) then
            call cesm_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:GLCPREP_BARRIER')
            call t_drvstartf ('CPL:GLCPREP',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)

            if (lnd_c2_glc) then
               call prep_glc_accum_avg(timer='CPL:glcprep_avg')

               ! Note that l2x_gx is obtained from mapping the module variable l2gacc_lx
               call prep_glc_calc_l2x_gx(fractions_lx, timer='CPL:glcprep_lnd2glc')

               call prep_glc_mrg(infodata, fractions_gx, timer_mrg='CPL:glcprep_mrgx2g')

               call component_diag(infodata, glc, flow='x2c', comment='send glc', &
                    info_debug=info_debug, timer_diag='CPL:glcprep_diagav')
            endif

            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('CPL:GLCPREP',cplrun=.true.)
         endif

         !----------------------------------------------------
         !| cpl -> glc
         !----------------------------------------------------

         if (iamin_CPLALLGLCID .and. glc_prognostic) then
            call component_exch(glc, flow='x2c', &
                 infodata=infodata, infodata_string='cpl2glc_run', &
                 mpicom_barrier=mpicom_CPLALLGLCID, run_barriers=run_barriers, &
                 timer_barrier='CPL:C2G_BARRIER', timer_comp_exch='CPL:C2G', &
                 timer_map_exch='CPL:c2g_glcx2glcg', timer_infodata_exch='CPL:c2g_infoexch')
         endif

      endif

      !----------------------------------------------------------
      !| ROF RECV-POST
      !----------------------------------------------------------

      if (rof_present .and. rofrun_alarm) then

         !----------------------------------------------------------
         !| rof -> cpl
         !----------------------------------------------------------

         if (iamin_CPLALLROFID) then
            call component_exch(rof, flow='c2x', &
                 infodata=infodata, infodata_string='rof2cpl_run', &
                 mpicom_barrier=mpicom_CPLALLROFID, run_barriers=run_barriers, &
                 timer_barrier='CPL:R2C_BARRIER', timer_comp_exch='CPL:R2C', &
                 timer_map_exch='CPL:r2c_rofr2rofx', timer_infodata_exch='CPL:r2c_infoexch')
         endif

         !----------------------------------------------------------
         !| rof post
         !----------------------------------------------------------

         if (iamin_CPLID) then
            call cesm_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:ROFPOST_BARRIER')
            call t_drvstartf  ('CPL:ROFPOST',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)

            call component_diag(infodata, rof, flow='c2x', comment= 'recv rof', &
                 info_debug=info_debug, timer_diag='CPL:rofpost_diagav')

            if (rof_c2_lnd) then
               call prep_lnd_calc_r2x_lx(timer='CPL:rofpost_rof2lnd')
            endif

            if (rof_c2_ice) then
               call prep_ice_calc_r2x_ix(timer='CPL:rofpost_rof2ice')
            endif

            if (rof_c2_ocn) then
               call prep_ocn_calc_r2x_ox(timer='CPL:rofpost_rof2ocn')
            endif

            call t_drvstopf  ('CPL:ROFPOST', cplrun=.true.)
         endif
      endif

      if (rof_present) then
         if (iamin_CPLID) then
            call cesm_comp_barriers(mpicom=mpicom_CPLID, timer='DRIVER_ROFPOST_BARRIER')
            call t_drvstartf  ('DRIVER_ROFPOST',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
            if (do_hist_r2x) then
               call t_drvstartf ('driver_rofpost_histaux', barrier=mpicom_CPLID)
               do eri = 1,num_inst_rof
                  suffix =  component_get_suffix(rof(eri)) 
                  call seq_hist_writeaux(infodata, EClock_d, rof(eri), flow='c2x', &
                       aname='r2x'//trim(suffix), dname='domrb', &
                       nx=rof_nx, ny=rof_ny, nt=1, write_now=t24hr_alarm)
               enddo
               call t_drvstopf ('driver_rofpost_histaux')
            endif
            call t_drvstopf  ('DRIVER_ROFPOST', cplrun=.true.)
         endif
      endif

      !----------------------------------------------------------
      !| Budget with old fractions
      !----------------------------------------------------------

      ! WJS (2-17-11): I am just using the first instance for the budgets because we
      ! don't expect budgets to be conserved for our case (I case). Also note that we
      ! don't expect budgets to be conserved for the interactive ensemble use case either.
      ! tcraig (aug 2012): put this after rof->cpl so the budget sees the new r2x_rx.
      ! it will also use the current r2x_ox here which is the value from the last timestep
      ! consistent with the ocean coupling

      if (iamin_CPLID .and. do_budgets) then
         call cesm_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:BUDGET1_BARRIER')
         call t_drvstartf ('CPL:BUDGET1',cplrun=.true.,budget=.true.,barrier=mpicom_CPLID)
         if (lnd_present) then
            call seq_diag_lnd_mct(lnd(ens1), fractions_lx(ens1), infodata, &
                 do_l2x=.true., do_x2l=.true.)
         endif
         if (rof_present) then
            call seq_diag_rof_mct(rof(ens1), fractions_rx(ens1), infodata)
         endif
         if (ice_present) then
            call seq_diag_ice_mct(ice(ens1), fractions_ix(ens1), infodata, &
                 do_x2i=.true.)
         endif
         call t_drvstopf  ('CPL:BUDGET1',cplrun=.true.,budget=.true.)
      endif


      !----------------------------------------------------------
      !| ICE RECV-POST
      !----------------------------------------------------------

      if (ice_present .and. icerun_alarm) then

         !----------------------------------------------------------
         !| ice -> cpl
         !----------------------------------------------------------

         if (iamin_CPLALLICEID) then
            call component_exch(ice, flow='c2x', &
                 infodata=infodata, infodata_string='ice2cpl_run', &
                 mpicom_barrier=mpicom_CPLALLICEID, run_barriers=run_barriers, &
                 timer_barrier='CPL:I2C_BARRIER', timer_comp_exch='CPL:I2C', &
                 timer_map_exch='CPL:i2c_icei2icex', timer_infodata_exch='CPL:i2c_infoexch')
         endif

         !----------------------------------------------------------
         !| ice post
         !----------------------------------------------------------

         if (iamin_CPLID) then
            call cesm_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:ICEPOST_BARRIER')
            call t_drvstartf  ('CPL:ICEPOST',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)

            call component_diag(infodata, ice, flow='c2x', comment= 'recv ice', &
                 info_debug=info_debug, timer_diag='CPL:icepost_diagav')

            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('CPL:ICEPOST',cplrun=.true.)
         endif
      endif

      !----------------------------------------------------------
      !| Update fractions based on new ice fractions
      !----------------------------------------------------------

      if (iamin_CPLID) then
         call cesm_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:FRACSET_BARRIER')
         call t_drvstartf ('CPL:FRACSET',cplrun=.true.,barrier=mpicom_CPLID)
         if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
         call t_drvstartf ('CPL:fracset_fracset',barrier=mpicom_CPLID)

         do efi = 1,num_inst_frc
            eii = mod((efi-1),num_inst_ice) + 1

            call seq_frac_set(infodata, ice(eii), &
                 fractions_ax(efi), fractions_ix(efi), fractions_ox(efi))
         enddo
         call t_drvstopf  ('CPL:fracset_fracset')

         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
         call t_drvstopf  ('CPL:FRACSET',cplrun=.true.)
      endif

      !----------------------------------------------------------
      !| ATM/OCN SETUP (rasm_option2)
      !----------------------------------------------------------

      if ((trim(cpl_seq_option) == 'RASM_OPTION2') .and. &
           iamin_CPLID .and. ocn_present) then

         call cesm_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:ATMOCN2_BARRIER')
         call t_drvstartf ('CPL:ATMOCN2',cplrun=.true.,barrier=mpicom_CPLID)
         if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)

         if (ocn_prognostic) then
            ! Map ice to ocn
            if (ice_c2_ocn) call prep_ocn_calc_i2x_ox(timer='CPL:atmocnp_ice2ocn')

            ! Map wav to ocn
            if (wav_c2_ocn) call prep_ocn_calc_w2x_ox(timer='CPL:atmocnp_wav2ocn')
         endif

         !----------------------------------------------------------
         !| atm/ocn flux on atm grid (rasm_option2 and aoflux_grid='atm')
         !----------------------------------------------------------

         if (trim(aoflux_grid) == 'atm') then
            ! compute o2x_ax for flux_atmocn, will be updated before atm merge
            ! can use fractions because fractions here are consistent with fractions in atm_mrg
            if (ocn_c2_atm) call prep_atm_calc_o2x_ax(fractions_ox,timer='CPL:atmoca_ocn2atm')

            call t_drvstartf ('CPL:atmocna_fluxa',barrier=mpicom_CPLID)
            do exi = 1,num_inst_xao
               eai = mod((exi-1),num_inst_atm) + 1
               eoi = mod((exi-1),num_inst_ocn) + 1
               efi = mod((exi-1),num_inst_frc) + 1
               a2x_ax => component_get_c2x_cx(atm(eai))
               o2x_ax => prep_atm_get_o2x_ax()    ! array over all instances
               xao_ax => prep_aoflux_get_xao_ax() ! array over all instances
               call seq_flux_atmocn_mct(infodata, tod, dtime, a2x_ax, o2x_ax(eoi), xao_ax(exi))
            enddo
            call t_drvstopf  ('CPL:atmocna_fluxa')

            if (atm_c2_ocn) call prep_aoflux_calc_xao_ox(timer='CPL:atmocna_atm2ocn')
         endif  ! aoflux_grid

         !----------------------------------------------------------
         !| atm/ocn flux on ocn grid (rasm_option2 and aoflux_grid='ocn')
         !----------------------------------------------------------

         if (trim(aoflux_grid) == 'ocn') then
            call t_drvstartf ('CPL:atmocnp_fluxo',barrier=mpicom_CPLID)
            do exi = 1,num_inst_xao
               eai = mod((exi-1),num_inst_atm) + 1
               eoi = mod((exi-1),num_inst_ocn) + 1
               efi = mod((exi-1),num_inst_frc) + 1
               a2x_ox => prep_ocn_get_a2x_ox()
               o2x_ox => component_get_c2x_cx(ocn(eoi))
               xao_ox => prep_aoflux_get_xao_ox()
               call seq_flux_atmocn_mct(infodata, tod, dtime, a2x_ox(eai), o2x_ox, xao_ox(exi))
            enddo
            call t_drvstopf  ('CPL:atmocnp_fluxo')
         endif  ! aoflux_grid

         !----------------------------------------------------------
         !| ocn prep-merge (rasm_option2)
         !----------------------------------------------------------

         xao_ox => prep_aoflux_get_xao_ox()
         call prep_ocn_mrg(infodata, fractions_ox, xao_ox=xao_ox, timer_mrg='CPL:atmocnp_mrgx2o')

         ! Accumulate ocn inputs - form partial sum of tavg ocn inputs (virtual "send" to ocn)
         call prep_ocn_accum(timer='CPL:atmocnp_accum')

         !----------------------------------------------------------
         !| ocn albedos (rasm_option2)
         !  (MUST BE AFTER prep_ocn_mrg for swnet to ocn to be computed properly
         !----------------------------------------------------------

         call t_drvstartf ('CPL:atmocnp_ocnalb', barrier=mpicom_CPLID)
         do exi = 1,num_inst_xao
            efi = mod((exi-1),num_inst_frc) + 1
            eai = mod((exi-1),num_inst_atm) + 1
            xao_ox => prep_aoflux_get_xao_ox()        ! array over all instances
            a2x_ox => prep_ocn_get_a2x_ox()
            call seq_flux_ocnalb_mct(infodata, ocn(1), a2x_ox(eai), fractions_ox(efi), xao_ox(exi))
         enddo
         call t_drvstopf  ('CPL:atmocnp_ocnalb')

         !----------------------------------------------------------
         !| ocn budget (rasm_option2)
         !----------------------------------------------------------

         if (do_budgets) then
            call cesm_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:BUDGET0_BARRIER')
            call t_drvstartf ('CPL:BUDGET0',cplrun=.true.,budget=.true.,barrier=mpicom_CPLID)
            xao_ox => prep_aoflux_get_xao_ox() ! array over all instances
            call seq_diag_ocn_mct(ocn(ens1), xao_ox(1), fractions_ox(ens1), infodata, &
                 do_o2x=.true., do_x2o=.true., do_xao=.true.)
            call t_drvstopf ('CPL:BUDGET0',cplrun=.true.,budget=.true.)
         endif

         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
         call t_drvstopf  ('CPL:ATMOCN2',cplrun=.true.)
      endif

      !----------------------------------------------------------
      !| OCN SETUP-SEND (rasm_option2)
      !----------------------------------------------------------

      if ((trim(cpl_seq_option) == 'RASM_OPTION2'  ) .and. &
           ocn_present .and. ocnrun_alarm) then

         !----------------------------------------------------
         ! "startup" wait (rasm_option2)
         !----------------------------------------------------

         if (iamin_CPLALLOCNID) then
            ! want to know the time the ocean pes waited for the cpl pes
            ! at the first ocnrun_alarm, min ocean wait is wait time
            ! do not use t_barrierf here since it can be "off", use mpi_barrier
            do eoi = 1,num_inst_ocn
               if (ocn(eoi)%iamin_compid) call t_drvstartf ('CPL:C2O_INITWAIT')
            enddo
            call mpi_barrier(mpicom_CPLALLOCNID,ierr)
            do eoi = 1,num_inst_ocn
               if (ocn(eoi)%iamin_compid) call t_drvstopf  ('CPL:C2O_INITWAIT')
            enddo
            cpl2ocn_first = .false.
         endif

         !----------------------------------------------------
         !| ocn average (rasm_option2)
         !----------------------------------------------------

         if (iamin_CPLID .and. ocn_prognostic) then
            call cesm_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:OCNPRE2_BARRIER')
            call t_drvstartf ('CPL:OCNPRE2',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)

            ! finish accumulating ocean inputs
            ! reset the value of x2o_ox with the value in x2oacc_ox
            ! (module variable in prep_ocn_mod)
            call prep_ocn_accum_avg(timer_accum='CPL:ocnprep_avg')

            call component_diag(infodata, ocn, flow='x2c', comment= 'send ocn', &
                 info_debug=info_debug, timer_diag='CPL:ocnprep_diagav')

            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('CPL:OCNPRE2',cplrun=.true.)
         endif

         !----------------------------------------------------
         !| cpl -> ocn (rasm_option2)
         !----------------------------------------------------

         if (iamin_CPLALLOCNID .and. ocn_prognostic) then
            call component_exch(ocn, flow='x2c', &
                 infodata=infodata, infodata_string='cpl2ocn_run', &
                 mpicom_barrier=mpicom_CPLALLOCNID, run_barriers=run_barriers, &
                 timer_barrier='CPL:C2O2_BARRIER', timer_comp_exch='CPL:C2O2', &
                 timer_map_exch='CPL:c2o2_ocnx2ocno', timer_infodata_exch='CPL:c2o2_infoexch')
         endif

      endif

      !----------------------------------------------------------
      !| ATM SETUP-SEND
      !----------------------------------------------------------

      if (atm_present .and. atmrun_alarm) then

         !----------------------------------------------------------
         !| atm prep-merge
         !----------------------------------------------------------

         if (iamin_CPLID .and. atm_prognostic) then
            call cesm_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:ATMPREP_BARRIER')
            call t_drvstartf ('CPL:ATMPREP',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)

            if (ocn_c2_atm) then
               if (trim(aoflux_grid) == 'ocn') then
                  ! map xao_ox states and fluxes to xao_ax if fluxes were computed on ocn grid
                  call prep_aoflux_calc_xao_ax(fractions_ox, flds='states_and_fluxes', &
                       timer='CPL:atmprep_xao2atm')
               endif

               ! recompute o2x_ax now for the merge with fractions associated with merge
               call prep_atm_calc_o2x_ax(fractions_ox, timer='CPL:atmprep_ocn2atm')

               ! map xao_ox albedos to the atm grid, these are always computed on the ocean grid
               call prep_aoflux_calc_xao_ax(fractions_ox, flds='albedos', timer='CPL:atmprep_alb2atm')
            endif

            if (ice_c2_atm) then
               call prep_atm_calc_i2x_ax(fractions_ix, timer='CPL:atmprep_ice2atm')
            endif

            if (lnd_c2_atm) then
               call prep_atm_calc_l2x_ax(fractions_lx, timer='CPL:atmprep_lnd2atm')
            endif

            if (associated(xao_ax)) then
               call prep_atm_mrg(infodata, fractions_ax, xao_ax=xao_ax, timer_mrg='CPL:atmprep_mrgx2a')
            endif

            call component_diag(infodata, atm, flow='x2c', comment= 'send atm', info_debug=info_debug, &
                 timer_diag='CPL:atmprep_diagav')

            call t_drvstopf  ('CPL:ATMPREP',cplrun=.true.)
            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
         endif

         !----------------------------------------------------------
         !| cpl -> atm
         !----------------------------------------------------------

         if (iamin_CPLALLATMID .and. atm_prognostic) then
            call component_exch(atm, flow='x2c', infodata=infodata, infodata_string='cpl2atm_run', &
                 mpicom_barrier=mpicom_CPLALLATMID, run_barriers=run_barriers, &
                 timer_barrier='CPL:C2A_BARRIER', timer_comp_exch='CPL:C2A', &
                 timer_map_exch='CPL:c2a_atmx2atmg', timer_infodata_exch='CPL:c2a_infoexch')
         endif

      endif

      !----------------------------------------------------------
      !| RUN OCN MODEL (NOT cesm1_orig_tight or cesm1_mod_tight)
      !----------------------------------------------------------

      if ((trim(cpl_seq_option) /= 'CESM1_ORIG_TIGHT' .and. &
           trim(cpl_seq_option) /= 'CESM1_MOD_TIGHT'   ) .and. &
          ocn_present .and. ocnrun_alarm) then
         call component_run(Eclock_o, ocn, ocn_run, infodata, &
              seq_flds_x2c_fluxes=seq_flds_x2o_fluxes, &
              seq_flds_c2x_fluxes=seq_flds_o2x_fluxes, &
              comp_prognostic=ocn_prognostic, comp_num=comp_num_ocn, &
              timer_barrier= 'CPL:OCN_RUN_BARRIER', timer_comp_run='CPL:OCN_RUN', &
              run_barriers=run_barriers, ymd=ymd, tod=tod,comp_layout=ocn_layout)
      endif

      !----------------------------------------------------------
      !| RUN ATM MODEL
      !----------------------------------------------------------

      if (atm_present .and. atmrun_alarm) then
         call component_run(Eclock_a, atm, atm_run, infodata, &
              seq_flds_x2c_fluxes=seq_flds_x2a_fluxes, &
              seq_flds_c2x_fluxes=seq_flds_a2x_fluxes, &
              comp_prognostic=atm_prognostic, comp_num=comp_num_atm, &
              timer_barrier= 'CPL:ATM_RUN_BARRIER', timer_comp_run='CPL:ATM_RUN', &
              run_barriers=run_barriers, ymd=ymd, tod=tod, comp_layout=atm_layout)
      endif

      !----------------------------------------------------------
      !| RUN GLC MODEL
      !----------------------------------------------------------

      if (glc_present .and. glcrun_alarm) then
         call component_run(Eclock_g, glc, glc_run, infodata, &
              seq_flds_x2c_fluxes=seq_flds_x2g_fluxes, &
              seq_flds_c2x_fluxes=seq_flds_g2x_fluxes, &
              comp_prognostic=glc_prognostic, comp_num=comp_num_glc, &
              timer_barrier= 'CPL:GLC_RUN_BARRIER', timer_comp_run='CPL:GLC_RUN', &
              run_barriers=run_barriers, ymd=ymd, tod=tod,comp_layout=glc_layout)
      endif

      !----------------------------------------------------------
      !| RUN ESP MODEL
      !----------------------------------------------------------
      if (esp_present .and. esprun_alarm) then
         call component_run(Eclock_e, esp, esp_run, infodata, &
              comp_prognostic=esp_prognostic, comp_num=comp_num_esp, &
              timer_barrier= 'CPL:ESP_RUN_BARRIER', timer_comp_run='CPL:ESP_RUN', &
              run_barriers=run_barriers, ymd=ymd, tod=tod,comp_layout=esp_layout)
      endif

      !----------------------------------------------------------
      !| WAV RECV-POST
      !----------------------------------------------------------

      if (wav_present .and. wavrun_alarm) then

         !----------------------------------------------------------
         !| wav -> cpl
         !----------------------------------------------------------

         if (iamin_CPLALLWAVID) then
            call component_exch(wav, flow='c2x', infodata=infodata, infodata_string='wav2cpl_run', &
                 mpicom_barrier=mpicom_CPLALLWAVID, run_barriers=run_barriers, &
                 timer_barrier='CPL:W2C_BARRIER', timer_comp_exch='CPL:W2C', &
                 timer_map_exch='CPL:w2c_wavw2wavx', timer_infodata_exch='CPL:w2c_infoexch')
         endif

         !----------------------------------------------------------
         !| wav post
         !----------------------------------------------------------

         if (iamin_CPLID) then
            call cesm_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:WAVPOST_BARRIER')
            call t_drvstartf  ('CPL:WAVPOST',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)

            call component_diag(infodata, wav, flow='c2x', comment= 'recv wav', &
                 info_debug=info_debug, timer_diag='CPL:wavpost_diagav')

            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('CPL:WAVPOST',cplrun=.true.)
         endif
      endif

      !----------------------------------------------------------
      !| GLC RECV-POST
      !----------------------------------------------------------

      if (glc_present .and. glcrun_alarm) then

         !----------------------------------------------------------
         !| glc -> cpl
         !----------------------------------------------------------

         if (iamin_CPLALLGLCID) then
            call component_exch(glc, flow='c2x', infodata=infodata, infodata_string='glc2cpl_run', &
                 mpicom_barrier=mpicom_CPLALLGLCID, run_barriers=run_barriers, &
                 timer_barrier='CPL:G2C_BARRIER', timer_comp_exch='CPL:G2C', &
                 timer_map_exch='CPL:g2c_glcg2glcx', timer_infodata_exch='CPL:g2c_infoexch')
         endif

         !----------------------------------------------------------
         !| glc post
         !----------------------------------------------------------

         if (iamin_CPLID) then
            call cesm_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:GLCPOST_BARRIER')
            call t_drvstartf  ('CPL:GLCPOST',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)

            call component_diag(infodata, glc, flow='c2x', comment= 'recv glc', &
                 info_debug=info_debug, timer_diag='CPL:glcpost_diagav')

            if (glc_c2_lnd) then
               call prep_lnd_calc_g2x_lx(timer='CPL:glcpost_glc2lnd')
            endif

            if (glc_c2_ice) then
               call prep_ice_calc_g2x_ix(timer='CPL:glcpost_glc2ice')
            endif

            if (glc_c2_ocn) then
               call prep_ocn_calc_g2x_ox(timer='CPL:glcpost_glc2ocn')
            endif

            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('CPL:GLCPOST',cplrun=.true.)
         endif
      endif

      !----------------------------------------------------------
      !| ATM RECV-POST
      !----------------------------------------------------------

      if (atm_present .and. atmrun_alarm) then

         !----------------------------------------------------------
         !| atm -> cpl
         !----------------------------------------------------------

         if (iamin_CPLALLATMID) then
            call component_exch(atm, flow='c2x', infodata=infodata, infodata_string='atm2cpl_run', &
                 mpicom_barrier=mpicom_CPLALLATMID, run_barriers=run_barriers, &
                 timer_barrier='CPL:A2C_BARRIER', timer_comp_exch='CPL:A2C', &
                 timer_map_exch='CPL:a2c_atma2atmx', timer_infodata_exch='CPL:a2c_infoexch')
         endif

         !----------------------------------------------------------
         !| atm post
         !----------------------------------------------------------

         if (iamin_CPLID) then
            call cesm_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:ATMPOST_BARRIER')
            call t_drvstartf ('CPL:ATMPOST',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)

            call component_diag(infodata, atm, flow='c2x', comment= 'recv atm', &
                 info_debug=info_debug, timer_diag='CPL:atmpost_diagav')

            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('CPL:ATMPOST',cplrun=.true.)
         endif
      endif

      !----------------------------------------------------------
      !| Budget with new fractions
      !----------------------------------------------------------

      if (iamin_CPLID .and. do_budgets) then
         call cesm_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:BUDGET2_BARRIER')

         call t_drvstartf ('CPL:BUDGET2',cplrun=.true.,budget=.true.,barrier=mpicom_CPLID)
         if (atm_present) then
            call seq_diag_atm_mct(atm(ens1), fractions_ax(ens1), infodata, &
                 do_a2x=.true., do_x2a=.true.)
         endif
         if (ice_present) then
            call seq_diag_ice_mct(ice(ens1), fractions_ix(ens1), infodata, &
                 do_i2x=.true.)
         endif
         call t_drvstopf  ('CPL:BUDGET2',cplrun=.true.,budget=.true.)

         call t_drvstartf ('CPL:BUDGET3',cplrun=.true.,budget=.true.,barrier=mpicom_CPLID)
         call seq_diag_accum_mct()
         call t_drvstopf  ('CPL:BUDGET3',cplrun=.true.,budget=.true.)

         call t_drvstartf ('CPL:BUDGETF',cplrun=.true.,budget=.true.,barrier=mpicom_CPLID)
         if (.not. dead_comps) then
            call seq_diag_print_mct(EClock_d,stop_alarm,budget_inst, &
                 budget_daily, budget_month, budget_ann, budget_ltann, budget_ltend)
         endif
         call seq_diag_zero_mct(EClock=EClock_d)

         call t_drvstopf  ('CPL:BUDGETF',cplrun=.true.,budget=.true.)
      endif

      !----------------------------------------------------------
      !| OCN RECV-POST (NOT cesm1_orig_tight and cesm1_mod_tight)
      !----------------------------------------------------------

      if ((trim(cpl_seq_option) /= 'CESM1_ORIG_TIGHT' .and. &
           trim(cpl_seq_option) /= 'CESM1_MOD_TIGHT'   ) .and. &
          ocn_present .and. ocnnext_alarm) then

         !----------------------------------------------------------
         !| ocn -> cpl (NOT cesm1_orig_tight and cesm1_mod_tight)
         !----------------------------------------------------------

         if (iamin_CPLALLOCNID) then
            call component_exch(ocn, flow='c2x', &
                 infodata=infodata, infodata_string='ocn2cpl_run', &
                 mpicom_barrier=mpicom_CPLALLOCNID, run_barriers=run_barriers, &
                 timer_barrier='CPL:O2C_BARRIER', timer_comp_exch='CPL:O2C', &
                 timer_map_exch='CPL:o2c_ocno2ocnx', timer_infodata_exch='CPL:o2c_infoexch')
         endif

         !----------------------------------------------------------
         !| ocn post (NOT cesm1_orig_tight and cesm1_mod_tight)
         !----------------------------------------------------------

         if (iamin_CPLID) then
            call cesm_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:OCNPOST_BARRIER')
            call t_drvstartf  ('CPL:OCNPOST',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)

            call component_diag(infodata, ocn, flow='c2x', comment= 'recv ocn', &
               info_debug=info_debug, timer_diag='CPL:ocnpost_diagav')

            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('CPL:OCNPOST',cplrun=.true.)
         endif
      endif

      !----------------------------------------------------------
      !| Write driver restart file
      !----------------------------------------------------------

      if ( restart_alarm .and. iamin_CPLID) then
         call cesm_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:RESTART_BARRIER')
         call t_drvstartf ('CPL:RESTART',cplrun=.true.,barrier=mpicom_CPLID)
         if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
         if (iamroot_CPLID) then
            write(logunit,104) ' Write restart file at ',ymd,tod
            call shr_sys_flush(logunit)
         endif

         call seq_rest_write(EClock_d, seq_SyncClock, infodata,       &
              atm, lnd, ice, ocn, rof, glc, wav, esp,                 &
              fractions_ax, fractions_lx, fractions_ix, fractions_ox, &
              fractions_rx, fractions_gx, fractions_wx)

         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
         call t_drvstopf  ('CPL:RESTART',cplrun=.true.)
      endif

      !----------------------------------------------------------
      !| Write history file, only AVs on CPLID
      !----------------------------------------------------------

      if (iamin_CPLID) then

         call cesm_comp_barriers(mpicom=mpicom_CPLID, timer='CPL:HISTORY_BARRIER')
         call t_drvstartf ('CPL:HISTORY',cplrun=.true.,barrier=mpicom_CPLID)
         if ( history_alarm) then
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
            if (iamroot_CPLID) then
               write(logunit,104) ' Write history file at ',ymd,tod
               call shr_sys_flush(logunit)
            endif

            call seq_hist_write(infodata, EClock_d, &
                 atm, lnd, ice, ocn, rof, glc, wav, &
                 fractions_ax, fractions_lx, fractions_ix, fractions_ox,     &
                 fractions_rx, fractions_gx, fractions_wx)

            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
         endif

         if (do_histavg) then
            call seq_hist_writeavg(infodata, EClock_d, &
                 atm, lnd, ice, ocn, rof, glc, wav, histavg_alarm)
         endif

         if (do_hist_a2x) then
            do eai = 1,num_inst_atm
               suffix =  component_get_suffix(atm(eai))
               if (trim(hist_a2x_flds) == 'all') then
                  call seq_hist_writeaux(infodata, EClock_d, atm(eai), flow='c2x', &
                       aname='a2x'//trim(suffix), dname='doma', &
                       nx=atm_nx, ny=atm_ny, nt=ncpl)
               else
                  call seq_hist_writeaux(infodata, EClock_d, atm(eai), flow='c2x', &
                       aname='a2x'//trim(suffix), dname='doma', &
                       nx=atm_nx, ny=atm_ny, nt=ncpl, flds=hist_a2x_flds)
               endif
            enddo
         endif

         if (do_hist_a2x1hri .and. t1hr_alarm) then
            do eai = 1,num_inst_atm
               suffix =  component_get_suffix(atm(eai))
               if (trim(hist_a2x1hri_flds) == 'all') then
                  call seq_hist_writeaux(infodata, EClock_d, atm(eai), flow='c2x', &
                       aname='a2x1hi'//trim(suffix), dname='doma', &
                       nx=atm_nx, ny=atm_ny, nt=24)
               else
                  call seq_hist_writeaux(infodata, EClock_d, atm(eai), flow='c2x', &
                       aname='a2x1hi'//trim(suffix), dname='doma', &
                       nx=atm_nx, ny=atm_ny, nt=24, flds=hist_a2x1hri_flds)
               endif
            enddo
         endif

         if (do_hist_a2x1hr) then
            do eai = 1,num_inst_atm
               suffix =  component_get_suffix(atm(eai))
               if (trim(hist_a2x1hr_flds) == 'all') then
                  call seq_hist_writeaux(infodata, EClock_d, atm(eai), flow='c2x', &
                       aname='a2x1h'//trim(suffix), dname='doma', &
                       nx=atm_nx, ny=atm_ny, nt=24, write_now=t1hr_alarm)
               else
                  call seq_hist_writeaux(infodata, EClock_d, atm(eai), flow='c2x', &
                       aname='a2x1h'//trim(suffix), dname='doma', &
                       nx=atm_nx, ny=atm_ny, nt=24, write_now=t1hr_alarm, flds=hist_a2x1hr_flds)
               endif
            enddo
         endif

         if (do_hist_a2x3hr) then
            do eai = 1,num_inst_atm
               suffix =  component_get_suffix(atm(eai))
               if (trim(hist_a2x3hr_flds) == 'all') then
                  call seq_hist_writeaux(infodata, EClock_d, atm(eai), flow='c2x', &
                       aname='a2x3h'//trim(suffix), dname='doma', &
                       nx=atm_nx, ny=atm_ny, nt=8, write_now=t3hr_alarm)
               else
                  call seq_hist_writeaux(infodata, EClock_d, atm(eai), flow='c2x', &
                       aname='a2x3h'//trim(suffix), dname='doma', &
                       nx=atm_nx, ny=atm_ny, nt=8, write_now=t3hr_alarm, flds=hist_a2x3hr_flds)
               endif
            enddo
         endif

         if (do_hist_a2x3hrp) then
            do eai = 1,num_inst_atm
               suffix = component_get_suffix(atm(eai))
               if (trim(hist_a2x3hrp_flds) == 'all') then
                  call seq_hist_writeaux(infodata, EClock_d, atm(eai), flow='c2x', &
                       aname='a2x3h_prec'//trim(suffix), dname='doma', &
                       nx=atm_nx, ny=atm_ny, nt=8, write_now=t3hr_alarm)
               else
                  call seq_hist_writeaux(infodata, EClock_d, atm(eai), flow='c2x', &
                       aname='a2x3h_prec'//trim(suffix), dname='doma', &
                       nx=atm_nx, ny=atm_ny, nt=8, write_now=t3hr_alarm, flds=hist_a2x3hrp_flds)
               endif
            enddo
         endif

         if (do_hist_a2x24hr) then
            do eai = 1,num_inst_atm
               suffix = component_get_suffix(atm(eai))
               if (trim(hist_a2x24hr_flds) == 'all') then
                  call seq_hist_writeaux(infodata, EClock_d, atm(eai), flow='c2x', &
                       aname='a2x1d'//trim(suffix), dname='doma', &
                       nx=atm_nx, ny=atm_ny, nt=1, write_now=t24hr_alarm)
               else
                  call seq_hist_writeaux(infodata, EClock_d, atm(eai), flow='c2x', &
                       aname='a2x1d'//trim(suffix), dname='doma', &
                       nx=atm_nx, ny=atm_ny, nt=1, write_now=t24hr_alarm, flds=hist_a2x24hr_flds)
               endif
            enddo
         endif

         if (do_hist_l2x1yr .and. glcrun_alarm) then
            ! Use yr_offset=-1 so the file with fields from year 1 has time stamp
            ! 0001-01-01 rather than 0002-01-01, etc.
            do eli = 1,num_inst_lnd
               suffix = component_get_suffix(lnd(eli))
               call seq_hist_writeaux(infodata, EClock_d, lnd(eli), flow='c2x', &
                    aname='l2x'//trim(suffix), dname='doml', &
                    nx=lnd_nx, ny=lnd_ny, nt=1, write_now=t1yr_alarm, yr_offset=-1)
            enddo
         endif

         if (do_hist_l2x) then
            do eli = 1,num_inst_lnd
               suffix =  component_get_suffix(lnd(eli))
               call seq_hist_writeaux(infodata, EClock_d, lnd(eli), flow='c2x', &
                    aname='l2x'//trim(suffix), dname='doml', &
                    nx=lnd_nx, ny=lnd_ny, nt=ncpl)
            enddo
         endif
         call t_drvstopf  ('CPL:HISTORY',cplrun=.true.)

      endif

      !----------------------------------------------------------
      !| Timing and memory diagnostics
      !----------------------------------------------------------

      call t_drvstartf ('CPL:TSTAMP_WRITE',cplrun=.true.)
      if (tod == 0 .or. info_debug > 1) then
         if (iamroot_CPLID) then
            call date_and_time(dstr,tstr)
            Time_estep = mpi_wtime()
            cktime = time_estep-time_bstep
            cktime_acc(1) = cktime_acc(1) + cktime
            cktime_cnt(1) = cktime_cnt(1) + 1
#ifndef CPL_BYPASS
            write(logunit,101) ' tStamp_write: model date = ',ymd,tod, &
                 ' wall clock = ',dstr(1:4),'-',dstr(5:6),'-',dstr(7:8),' ',&
                 tstr(1:2),':',tstr(3:4),':',tstr(5:6), &
                 ' avg dt = ',cktime_acc(1)/cktime_cnt(1),' dt = ',cktime
#endif
            Time_bstep = mpi_wtime()
            call shr_sys_flush(logunit)
            if(cktime > max_cplstep_time .and. max_cplstep_time > 0.0) then
               call shr_sys_abort(subname//'Wall clock time exceeds max_cplstep_time')
            else if(max_cplstep_time < -0.05) then
               ! if max_cplstep_time is < 0 we use abs(max_cplstep_time)
               ! times the initial cktime value as a threshhold
               max_cplstep_time = -(max_cplstep_time)*cktime
            endif
         endif
      end if
      if (tod == 0 .and. wall_time_limit > 0.0_r8 .and. .not. force_stop) then
         time_erun = mpi_wtime()
         ! time_*run is seconds, wall_time_limit is hours
         wall_time = (time_erun - time_brun) / 3600._r8   ! convert secs to hrs
         write(logunit,109) subname//' check wall_time_limit: ',wall_time, wall_time_limit
         if (wall_time > wall_time_limit) then
            force_stop = .true.
            force_stop_tod = 0
            if (trim(force_stop_at) == 'month') then
               call shr_cal_date2ymd(ymd,year,month,day)
               month = month + 1
               do while (month > 12)
                  month = month - 12
                  year = year + 1
               enddo
               call shr_cal_ymd2date(year,month,1,force_stop_ymd)
            elseif (trim(force_stop_at) == 'year') then  ! next year
               call shr_cal_date2ymd(ymd,year,month,day)
               call shr_cal_ymd2date(year+1,1,1,force_stop_ymd)
            elseif (trim(force_stop_at) == 'day') then   ! next day
               ymdtmp = ymd
               call shr_cal_advDateInt(1,'days'  ,ymdtmp,0,force_stop_ymd,todtmp,calendar)
            else    ! day is default
               ymdtmp = ymd
               call shr_cal_advDateInt(1,'days'  ,ymdtmp,0,force_stop_ymd,todtmp,calendar)
            endif
            write(logunit,108) subname//' reached wall_time_limit (hours) =',wall_time_limit, &
                               ' :stop at ',force_stop_ymd
         endif
      endif
#ifndef CPL_BYPASS
      if (tod == 0 .or. info_debug > 1) then
         !! Report on memory usage
         !! For now, just look at the first instance of each component
         if ( iamroot_CPLID .or. &
              ocn(ens1)%iamroot_compid .or. &
              atm(ens1)%iamroot_compid .or. &
              lnd(ens1)%iamroot_compid .or. &
              ice(ens1)%iamroot_compid .or. &
              glc(ens1)%iamroot_compid .or. &
              wav(ens1)%iamroot_compid) then
            call shr_mem_getusage(msize,mrss)

            write(logunit,105) ' memory_write: model date = ',ymd,tod, &
                 ' memory = ',mrss,' MB (highwater)    ',msize,' MB (usage)', &
                 '  (pe=',iam_GLOID,' comps=',trim(complist)//')'
         endif
      endif
#endif
      if (info_debug > 1) then
         if (iamroot_CPLID) then
            call seq_infodata_GetData(infodata,nextsw_cday=nextsw_cday)
            !            write(logunit,106) ' nextsw_cday = ',nextsw_cday
            write(logunit,*) '  nextsw_cday = ',nextsw_cday
         endif
      endif
      call t_drvstopf  ('CPL:TSTAMP_WRITE',cplrun=.true.)

      call t_stopf  ('CPL:RUN_LOOP', hashint(1))

      ! --- Write out performance data
      call t_startf  ('CPL:TPROF_WRITE')
      if (tprof_alarm) then
         call t_adj_detailf(+1)

         call t_startf("sync1_tprof")
         call mpi_barrier(mpicom_GLOID,ierr)
         call t_stopf("sync1_tprof")

         write(timing_file,'(a,i8.8,a1,i5.5)') trim(tchkpt_dir)//"/model_timing_",ymd,"_",tod
         if (output_perf) then
            call t_prf(filename=trim(timing_file), mpicom=mpicom_GLOID, &
                       num_outpe=0, output_thispe=output_perf)
         else
            call t_prf(filename=trim(timing_file), mpicom=mpicom_GLOID, &
                       num_outpe=0)
         endif

         call t_startf("sync2_tprof")
         call mpi_barrier(mpicom_GLOID,ierr)
         call t_stopf("sync2_tprof")

         call t_adj_detailf(-1)
      endif
      call t_stopf  ('CPL:TPROF_WRITE')

      call t_drvstartf  ('CPL:BARRIERALARM',cplrun=.true.)
      if (barrier_alarm) then
         call mpi_barrier(mpicom_GLOID,ierr)
      endif
      call t_drvstopf   ('CPL:BARRIERALARM',cplrun=.true.)

   enddo   ! driver run loop

   !|----------------------------------------------------------
   !| End of driver time step loop
   !|---------------------------------------------------------

   call t_startf ('CPL:RUN_LOOP_BSTOP')
   call mpi_barrier(mpicom_GLOID,ierr)
   call t_stopf ('CPL:RUN_LOOP_BSTOP')

   Time_end = mpi_wtime()

 end subroutine cesm_run

!===============================================================================
!*******************************************************************************
!===============================================================================

 subroutine cesm_final()

   use shr_pio_mod, only : shr_pio_finalize
   use shr_wv_sat_mod, only: shr_wv_sat_final
   implicit none

   !------------------------------------------------------------------------
   ! Finalization of all models
   !------------------------------------------------------------------------

   call t_barrierf ('CPL:FINAL_BARRIER', mpicom_GLOID)
   call t_startf ('CPL:FINAL')
   call t_adj_detailf(+1)

   call t_startf('cesm_final')
   call t_adj_detailf(+1)

   call seq_timemgr_EClockGetData( EClock_d, stepno=endstep)
   call shr_mem_getusage(msize,mrss)

   call component_final(EClock_a, atm, atm_final)
   call component_final(EClock_l, lnd, lnd_final)
   call component_final(EClock_r, rof, rof_final)
   call component_final(EClock_i, ice, ice_final)
   call component_final(EClock_o, ocn, ocn_final)
   call component_final(EClock_g, glc, glc_final)
   call component_final(EClock_w, wav, wav_final)

   !------------------------------------------------------------------------
   ! End the run cleanly
   !------------------------------------------------------------------------

   call shr_wv_sat_final()

   call shr_pio_finalize( )

   call shr_mpi_min(msize ,msize0,mpicom_GLOID,' driver msize0', all=.true.)
   call shr_mpi_max(msize ,msize1,mpicom_GLOID,' driver msize1', all=.true.)
   call shr_mpi_min(mrss  ,mrss0,mpicom_GLOID,'  driver mrss0',  all=.true.)
   call shr_mpi_max(mrss  ,mrss1,mpicom_GLOID,'  driver mrss1',  all=.true.)

   if (iamroot_CPLID )then
      call seq_timemgr_EClockGetData( EClock_d, curr_ymd=ymd, curr_tod=tod, dtime=dtime)
      simDays = (endStep-begStep)*dtime/(24._r8*3600._r8)
      write(logunit,'(//)')
      write(logunit,FormatA) subname, 'SUCCESSFUL TERMINATION OF CPL7-CESM'
      write(logunit,FormatD) subname, '  at YMD,TOD = ',ymd,tod
      write(logunit,FormatR) subname, '# simulated days (this run) = ', simDays
      write(logunit,FormatR) subname, 'compute time (hrs)          = ', (Time_end-Time_begin)/3600._r8
      if ( (Time_end /= Time_begin) .and. (simDays /= 0.0_r8) )then
         SYPD = shr_const_cday*simDays/(days_per_year*(Time_end-Time_begin))
         write(logunit,FormatR) subname, '# simulated years / cmp-day = ', SYPD
      endif
      write(logunit,FormatR) subname,' pes min memory highwater  (MB)  = ',mrss0
      write(logunit,FormatR) subname,' pes max memory highwater  (MB)  = ',mrss1
      write(logunit,FormatR) subname,' pes min memory last usage (MB)  = ',msize0
      write(logunit,FormatR) subname,' pes max memory last usage (MB)  = ',msize1
      write(logunit,'(//)')
      close(logunit)
   endif

   call t_adj_detailf(-1)
   call t_stopf('cesm_final')

   call t_startf("final:sync1_tprof")
   call mpi_barrier(mpicom_GLOID,ierr)
   call t_stopf("final:sync1_tprof")

   call t_adj_detailf(-1)
   call t_stopf  ('CPL:FINAL')

   call t_set_prefixf("final:")
   if (output_perf) then
      call t_prf(trim(timing_dir)//'/model_timing', mpicom=mpicom_GLOID, &
                 output_thispe=output_perf)
   else
      call t_prf(trim(timing_dir)//'/model_timing', mpicom=mpicom_GLOID)
   endif
   call t_unset_prefixf()

   call t_finalizef()

end subroutine cesm_final

!===============================================================================
!*******************************************************************************
!===============================================================================

subroutine seq_cesm_printlogheader()

  !-----------------------------------------------------------------------
  !
  ! Purpose: Print basic information on what this driver program is
  ! to the logfile.
  !
  !-----------------------------------------------------------------------
  !
  ! Local variables
  !
  implicit none

  character(len=8) :: cdate          ! System date
  character(len=8) :: ctime          ! System time
  integer          :: values(8)
  character        :: date*8, time*10, zone*5

!-------------------------------------------------------------------------------

  call date_and_time (date, time, zone, values)
  cdate(1:2) = date(5:6)
  cdate(3:3) = '/'
  cdate(4:5) = date(7:8)
  cdate(6:6) = '/'
  cdate(7:8) = date(3:4)
  ctime(1:2) = time(1:2)
  ctime(3:3) = ':'
  ctime(4:5) = time(3:4)
  ctime(6:6) = ':'
  ctime(7:8) = time(5:6)
  write(logunit,F00) '------------------------------------------------------------'
  write(logunit,F00) '        NCAR CPL7 Community Earth System Model (CESM)  '
  write(logunit,F00) '------------------------------------------------------------'
  write(logunit,F00) '     (Online documentation is available on the CESM         '
  write(logunit,F00) '      Models page: http://www.cesm.ucar.edu/models/         '
  write(logunit,F00) '      License information is available as a link from above '
  write(logunit,F00) '------------------------------------------------------------'
  write(logunit,F00) '                DATE ',cdate, ' TIME ', ctime
  write(logunit,F00) '------------------------------------------------------------'
  write(logunit,*)' '
  write(logunit,*)' '

end subroutine seq_cesm_printlogheader

!===============================================================================

subroutine cesm_comp_barriers(mpicom, timer)
  implicit none
  integer         , intent(in) :: mpicom
  character(len=*), intent(in) :: timer
  integer :: ierr

  if (run_barriers) then
     call t_drvstartf (trim(timer))
     call mpi_barrier(mpicom,ierr)
     call t_drvstopf (trim(timer))
  endif
end subroutine cesm_comp_barriers

end module cesm_comp_mod
