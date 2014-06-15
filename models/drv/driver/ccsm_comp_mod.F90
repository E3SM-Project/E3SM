module ccsm_comp_mod

!-------------------------------------------------------------------------------
!
! Purpose: Main program for NCAR CCSM4/cpl7. Can have different
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
   use shr_cal_mod,       only: shr_cal_date2ymd
   use shr_orb_mod,       only: shr_orb_params
   use shr_reprosum_mod,  only: shr_reprosum_setopts
   use mct_mod            ! mct_ wrappers for mct lib
   use perf_mod
   use ESMF

   !----------------------------------------------------------------------------
   ! component model interfaces (init, run, final methods)
   !----------------------------------------------------------------------------
   use atm_comp_mct, only: atm_init_mct, atm_run_mct, atm_final_mct
   use lnd_comp_mct, only: lnd_init_mct, lnd_run_mct, lnd_final_mct
   use ocn_comp_mct, only: ocn_init_mct, ocn_run_mct, ocn_final_mct
   use ice_comp_mct, only: ice_init_mct, ice_run_mct, ice_final_mct
   use glc_comp_mct, only: glc_init_mct, glc_run_mct, glc_final_mct
   use wav_comp_mct, only: wav_init_mct, wav_run_mct, wav_final_mct
   use rof_comp_mct, only: rof_init_mct, rof_run_mct, rof_final_mct

   !----------------------------------------------------------------------------
   ! cpl7 modules
   !----------------------------------------------------------------------------

   !--- modules with public read/write data ---
   use seq_avdata_mod    ! drv aVects & associated domain, fraction, cdata
   use seq_diag_mct      ! diagnostic routines

   !--- other ---
   use seq_comm_mct      ! mpi comm data & routines, plus logunit and loglevel
   use seq_timemgr_mod   ! clock & alarm routines 
   use seq_infodata_mod  ! "infodata" gathers various control flags into one datatype
   use seq_cdata_mod     ! "cdata" type & methods (domain + decomp + infodata in one datatype)
   use seq_domain_mct    ! domain related routines
   use seq_flux_mct      ! flux calc routines
   use seq_frac_mct      ! domain fraction routines
   use seq_rest_mod      ! restart file routines
   use seq_hist_mod      ! history file routines
   use seq_io_mod        ! i/o subroutines
   use seq_mctext_mod    ! rearrange type routines
   use seq_flds_mod

   !--- merging routines ---
   use mrg_mod          ! merge gridded components

   !--- mapping routines ---
   use seq_map_mod      ! generic mapping

   implicit none

   private

   public ccsm_pre_init, ccsm_init, ccsm_run, ccsm_final
#ifdef USE_ESMF_LIB
   public ccsm_comp_register
#endif
   public timing_dir, mpicom_GLOID

#include <mpif.h>

   !----------------------------------------------------------------------------
   ! domains & related
   !----------------------------------------------------------------------------

   !--- domain decomps (MCT Global Seg Maps) ---
   type(mct_gsMap), target  :: gsMap_aa(num_inst_atm)    ! on component pes
   type(mct_gsMap), target  :: gsMap_ll(num_inst_lnd)
   type(mct_gsMap), target  :: gsMap_oo(num_inst_ocn)
   type(mct_gsMap), target  :: gsMap_ii(num_inst_ice)
   type(mct_gsMap), target  :: gsMap_rr(num_inst_rof)
   type(mct_gsMap), target  :: gsMap_gg(num_inst_glc)
   type(mct_gsMap), target  :: gsMap_ss(num_inst_lnd)
   type(mct_gsMap), target  :: gsMap_ww(num_inst_wav)

   type(mct_gsMap), target  :: gsMap_ax    ! on cpl pes
   type(mct_gsMap), target  :: gsMap_lx
   type(mct_gsMap), target  :: gsMap_ox
   type(mct_gsMap), target  :: gsMap_ix
   type(mct_gsMap), target  :: gsMap_rx
   type(mct_gsMap), target  :: gsMap_gx
   type(mct_gsMap), target  :: gsMap_sx
   type(mct_gsMap), target  :: gsMap_wx

   type(mct_gGrid)  :: dom_tmp   ! temporary

   integer          :: lsize_a   ! local size of atm AV
   integer          :: lsize_l
   integer          :: lsize_r
   integer          :: lsize_s
   integer          :: lsize_o
   integer          :: lsize_i
   integer          :: lsize_g
   integer          :: lsize_w

   !--- domain area correction factors (only defined on cpl pes) ---

   type AreaCorrectFactor
      real(r8), pointer :: drv2mdl(:), mdl2drv(:)
   end type AreaCorrectFactor

   type(AreaCorrectFactor) :: areacor_aa(num_inst_atm)
   type(AreaCorrectFactor) :: areacor_ll(num_inst_lnd)
   type(AreaCorrectFactor) :: areacor_rr(num_inst_rof)
   type(AreaCorrectFactor) :: areacor_ii(num_inst_ice)
   type(AreaCorrectFactor) :: areacor_oo(num_inst_ocn)
   type(AreaCorrectFactor) :: areacor_gg(num_inst_glc)
   type(AreaCorrectFactor) :: areacor_ss(num_inst_lnd)
   type(AreaCorrectFactor) :: areacor_ww(num_inst_wav)

   !--- domain equivalent 2d grid size ---
   integer          :: atm_nx, atm_ny  ! nx, ny of 2d grid, if known
   integer          :: lnd_nx, lnd_ny
   integer          :: ice_nx, ice_ny
   integer          :: ocn_nx, ocn_ny
   integer          :: rof_nx, rof_ny
   integer          :: glc_nx, glc_ny
   integer          :: sno_nx, sno_ny
   integer          :: wav_nx, wav_ny


   !----------------------------------------------------------------------------
   ! mappers 
   !   _C are component/coupler rearrangers
   !   _S are for states
   !   _F are for fluxes
   !   _SF are for states and fluxes
   !----------------------------------------------------------------------------

   type(seq_map), SAVE :: mapper_Ca2x(num_inst_atm)
   type(seq_map), SAVE :: mapper_Cx2a(num_inst_atm)
   type(seq_map), SAVE :: mapper_Cl2x(num_inst_lnd)
   type(seq_map), SAVE :: mapper_Cx2l(num_inst_lnd)
   type(seq_map), SAVE :: mapper_Cs2x(num_inst_lnd)
   type(seq_map), SAVE :: mapper_Cx2s(num_inst_lnd)
   type(seq_map), SAVE :: mapper_Cr2x(num_inst_rof)
   type(seq_map), SAVE :: mapper_Cx2r(num_inst_rof)
   type(seq_map), SAVE :: mapper_Ci2x(num_inst_ice)
   type(seq_map), SAVE :: mapper_Cx2i(num_inst_ice)
   type(seq_map), SAVE :: mapper_Co2x(num_inst_ocn)
   type(seq_map), SAVE :: mapper_Cx2o(num_inst_ocn)
   type(seq_map), SAVE :: mapper_Cg2x(num_inst_glc)
   type(seq_map), SAVE :: mapper_Cx2g(num_inst_glc)
   type(seq_map), SAVE :: mapper_Cw2x(num_inst_wav)
   type(seq_map), SAVE :: mapper_Cx2w(num_inst_wav)

   type(seq_map), SAVE :: mapper_Sa2o
   type(seq_map), SAVE :: mapper_Va2o
   type(seq_map), SAVE :: mapper_Fa2o
   type(seq_map), SAVE :: mapper_So2a
   type(seq_map), SAVE :: mapper_Fo2a
   type(seq_map), SAVE :: mapper_Sa2l
   type(seq_map), SAVE :: mapper_Fa2l
   type(seq_map), SAVE :: mapper_Sl2a
   type(seq_map), SAVE :: mapper_Fl2a
   type(seq_map), SAVE :: mapper_Fl2r
   type(seq_map), SAVE :: mapper_Si2a
   type(seq_map), SAVE :: mapper_Fi2a
   type(seq_map), SAVE :: mapper_Fr2o
   type(seq_map), SAVE :: mapper_Rr2o
   type(seq_map), SAVE :: mapper_Fr2l
   type(seq_map), SAVE :: mapper_Sr2l
   type(seq_map), SAVE :: mapper_SFi2o
   type(seq_map), SAVE :: mapper_SFo2i
   type(seq_map), SAVE :: mapper_SFg2s
   type(seq_map), SAVE :: mapper_SFs2g
   type(seq_map), SAVE :: mapper_Sa2w
   type(seq_map), SAVE :: mapper_So2w
   type(seq_map), SAVE :: mapper_Si2w
   type(seq_map), SAVE :: mapper_Sw2o

   !----------------------------------------------------------------------------
   ! time management
   !----------------------------------------------------------------------------

   type (seq_timemgr_type), SAVE :: seq_SyncClock ! array of all clocks & alarm
   type (ESMF_Clock), target, SAVE       :: EClock_d      ! driver clock
   type (ESMF_Clock), target, SAVE       :: EClock_a
   type (ESMF_Clock), target, SAVE       :: EClock_l
   type (ESMF_Clock), target, SAVE       :: EClock_o
   type (ESMF_Clock), target, SAVE       :: EClock_i
   type (ESMF_Clock), target, SAVE       :: EClock_g
   type (ESMF_Clock), target, SAVE       :: EClock_r
   type (ESMF_Clock), target, SAVE       :: EClock_w

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
   logical  :: tprof_alarm            ! timing profile alarm
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
   character(CL) :: orb_mode          ! orbital mode
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

   !--- for documenting speed of the model ---
   character( 8) :: dstr              ! date string
   character(10) :: tstr              ! time string
   integer       :: begStep, endStep  ! Begining and ending step number
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
   logical  :: rof_present            ! .true.  => rof is present
   logical  :: flood_present          ! .true.  => rof is computing flood
   logical  :: sno_present            ! .true.  => land sno is present
   logical  :: wav_present            ! .true.  => wav is present

   logical  :: atm_prognostic         ! .true.  => atm comp expects input
   logical  :: lnd_prognostic         ! .true.  => lnd comp expects input
   logical  :: ice_prognostic         ! .true.  => ice comp expects input
   logical  :: ocn_prognostic         ! .true.  => ocn comp expects input
   logical  :: ocnrof_prognostic      ! .true.  => ocn comp expects runoff input
   logical  :: glc_prognostic         ! .true.  => glc comp expects input
   logical  :: rof_prognostic         ! .true.  => rof comp expects input
   logical  :: sno_prognostic         ! .true.  => sno comp expects input
   logical  :: wav_prognostic         ! .true.  => wav comp expects input

   logical  :: dead_comps             ! .true.  => dead components 
   logical  :: esmf_map_flag          ! .true.  => use esmf for mapping

   logical  :: single_column          ! scm mode logical
   real(r8) :: scmlon                 ! single column lon
   real(r8) :: scmlat                 ! single column lat
   logical  :: aqua_planet            ! aqua planet mode
   real(r8) :: nextsw_cday            ! radiation control
   logical  :: atm_aero               ! atm provides aerosol data
   real(r8) :: flux_epbalfact         ! precip factor

   logical  :: ocean_tight_coupling   ! couple ocn as frequently as lnd & ice
   logical  :: skip_ocean_run         ! skip the ocean model first pass
   logical  :: cpl2ocn_first          ! use to call initial cpl2ocn timer
   character(CS) :: aoflux_grid       ! grid for a/o flux calc: atm xor ocn 
   character(CS) :: vect_map          ! vector mapping type
   logical  :: run_barriers           ! barrier the component run calls
   logical  :: samegrid_ao            ! samegrid atm and ocean
   logical  :: samegrid_al            ! samegrid atm and land
   logical  :: samegrid_ro            ! samegrid runoff and ocean
   logical  :: samegrid_aw            ! samegrid atm and wave
   logical  :: samegrid_ow            ! samegrid ocean and wave

   logical       :: read_restart      ! local read restart flag
   character(CL) :: rest_file         ! restart file path + filename

   logical  :: domain_check           ! .true.  => check consistency of domains
   logical  :: shr_map_dopole         ! logical for dopole in shr_map_mod
   logical  :: reprosum_use_ddpdd     ! setup reprosum, use ddpdd
   real(r8) :: reprosum_diffmax       ! setup reprosum, set rel_diff_max
   logical  :: reprosum_recompute     ! setup reprosum, recompute if tolerance exceeded

   logical  :: output_perf = .false.  ! require timing data output for this pe

   !--- history & budgets ---
   logical       :: do_budgets        ! heat/water budgets on
   logical       :: do_histinit       ! initial hist file
   logical       :: do_histavg        ! histavg on or off
   logical       :: do_hist_r2x       ! create aux files: r2x
   logical       :: do_hist_l2x       ! create aux files: l2x
   logical       :: do_hist_a2x24hr   ! create aux files: a2x
   logical       :: do_hist_s2x1yr    ! create aux files: s2x
   logical       :: do_hist_a2x       ! create aux files: a2x
   logical       :: do_hist_a2x3hrp   ! create aux files: a2x 3hr precip
   logical       :: do_hist_a2x3hr    ! create aux files: a2x 3hr states
!  character(CL) :: hist_r2x_flds     = 'all'
!  character(CL) :: hist_l2x_flds     = 'all'
   character(CL) :: hist_a2x_flds     = 'Faxa_swndr:Faxa_swvdr:Faxa_swndf:Faxa_swvdf'
!  character(CL) :: hist_a2x24hr_flds = 'all'
   character(CL) :: hist_a2x3hrp_flds = 'Faxa_rainc:Faxa_rainl:Faxa_snowc:Faxa_snowl'
   character(CL) :: hist_a2x3hr_flds  = 'Sa_z:Sa_u:Sa_v:Sa_tbot:Sa_ptem:Sa_shum:Sa_dens:Sa_pbot:Sa_pslv:Faxa_lwdn'
   integer  :: budget_inst            ! instantaneous budget flag
   integer  :: budget_daily           ! daily budget flag
   integer  :: budget_month           ! monthly budget flag
   integer  :: budget_ann             ! annual budget flag
   integer  :: budget_ltann           ! long term budget flag for end of year writing
   integer  :: budget_ltend           ! long term budget flag for end of run writing

   ! --- field indexes used in ccsm_comp ---
   integer  :: index_r2x_Forr_roff
   integer  :: index_r2x_Forr_ioff
   integer  :: index_r2x_Flrr_flood
   integer  :: index_r2x_Slrr_volr
   integer  :: index_l2x_Flrl_rofliq
   integer  :: index_l2x_Flrl_rofice
   integer  :: index_x2r_Flrl_rofliq
   integer  :: index_x2r_Flrl_rofice

   ! --- other ---
   integer  :: ka,km,k1,k2,k3         ! aVect field indices
   integer  :: ocnrun_count           ! number of times ocn run alarm went on
   logical  :: exists                 ! true if file exists
   integer  :: ierr                   ! MPI error return
   integer  :: rc                     ! return code
   logical  :: cdf64                  ! true => use 64 bit addressing in netCDF files
   type(mct_aVect) :: x2r_rx_tmp      ! temporary

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

   integer  :: pethreads_GLOID        ! OMP number of threads per task

   integer  :: nthreads_CPLATMID      ! OMP cpl-atm number of threads
   integer  :: nthreads_CPLLNDID      ! OMP cpl-lnd number of threads
   integer  :: nthreads_CPLICEID      ! OMP cpl-ice number of threads
   integer  :: nthreads_CPLOCNID      ! OMP cpl-ocn number of threads
   integer  :: nthreads_CPLGLCID      ! OMP cpl-glc number of threads
   integer  :: nthreads_CPLROFID      ! OMP cpl-glc number of threads
   integer  :: nthreads_CPLWAVID      ! OMP cpl-wav number of threads

   logical  :: drv_threading          ! driver threading control

   !----------------------------------------------------------------------------
   ! communicator groups and related
   !----------------------------------------------------------------------------
   integer  :: Global_Comm
   integer  :: mpicom_GLOID                         ! MPI global communicator
   integer  :: mpicom_CPLID                         ! MPI cpl communicator

   integer  :: mpicom_CPLALLATMID                   ! MPI comm for CPLALLATMID
   integer  :: mpicom_CPLALLLNDID                   ! MPI comm for CPLALLLNDID
   integer  :: mpicom_CPLALLICEID                   ! MPI comm for CPLALLICEID
   integer  :: mpicom_CPLALLOCNID                   ! MPI comm for CPLALLOCNID
   integer  :: mpicom_CPLALLGLCID                   ! MPI comm for CPLALLGLCID
   integer  :: mpicom_CPLALLROFID                   ! MPI comm for CPLALLROFID
   integer  :: mpicom_CPLALLWAVID                   ! MPI comm for CPLALLWAVID

   integer  :: mpicom_ATMID(num_inst_atm)           ! MPI atm communicator
   integer  :: mpicom_LNDID(num_inst_lnd)           ! MPI lnd communicator
   integer  :: mpicom_ICEID(num_inst_ice)           ! MPI ice communicator
   integer  :: mpicom_OCNID(num_inst_ocn)           ! MPI ocn communicator
   integer  :: mpicom_GLCID(num_inst_glc)           ! MPI glc communicator
   integer  :: mpicom_ROFID(num_inst_rof)           ! MPI rof communicator
   integer  :: mpicom_WAVID(num_inst_wav)           ! MPI wav communicator

   integer  :: mpicom_CPLATMID(num_inst_atm)        ! MPI cpl-atm communicator
   integer  :: mpicom_CPLLNDID(num_inst_lnd)        ! MPI cpl-lnd communicator
   integer  :: mpicom_CPLICEID(num_inst_ice)        ! MPI cpl-ice communicator
   integer  :: mpicom_CPLOCNID(num_inst_ocn)        ! MPI cpl-ocn communicator
   integer  :: mpicom_CPLGLCID(num_inst_glc)        ! MPI cpl-glc communicator
   integer  :: mpicom_CPLROFID(num_inst_rof)        ! MPI cpl-rof communicator
   integer  :: mpicom_CPLWAVID(num_inst_wav)        ! MPI cpl-wav communicator

   logical  :: iamroot_GLOID                        ! GLOID masterproc
   logical  :: iamroot_CPLID                        ! CPLID masterproc
   logical  :: iamroot_ATMID(num_inst_atm)          ! ATMID masterproc
   logical  :: iamroot_LNDID(num_inst_lnd)          ! LNDID masterproc
   logical  :: iamroot_ICEID(num_inst_ice)          ! ICEID masterproc
   logical  :: iamroot_OCNID(num_inst_ocn)          ! OCNID masterproc
   logical  :: iamroot_GLCID(num_inst_glc)          ! GLCID masterproc
   logical  :: iamroot_ROFID(num_inst_rof)          ! ROFID masterproc
   logical  :: iamroot_WAVID(num_inst_wav)          ! WAVID masterproc

   logical  :: iamin_CPLID                          ! pe associated with CPLID
   logical  :: iamin_CPLALLATMID                    ! pe associated with CPLALLATMID
   logical  :: iamin_CPLALLLNDID                    ! pe associated with CPLALLLNDID
   logical  :: iamin_CPLALLICEID                    ! pe associated with CPLALLICEID
   logical  :: iamin_CPLALLOCNID                    ! pe associated with CPLALLOCNID
   logical  :: iamin_CPLALLGLCID                    ! pe associated with CPLALLGLCID
   logical  :: iamin_CPLALLROFID                    ! pe associated with CPLALLROFID
   logical  :: iamin_CPLALLWAVID                    ! pe associated with CPLALLWAVID

   logical  :: iamin_ATMID(num_inst_atm)            ! pe associated with ATMID
   logical  :: iamin_LNDID(num_inst_lnd)            ! pe associated with LNDID
   logical  :: iamin_ICEID(num_inst_ice)            ! pe associated with ICEID
   logical  :: iamin_OCNID(num_inst_ocn)            ! pe associated with OCNID
   logical  :: iamin_GLCID(num_inst_glc)            ! pe associated with GLCID
   logical  :: iamin_ROFID(num_inst_rof)            ! pe associated with ROFID
   logical  :: iamin_WAVID(num_inst_wav)            ! pe associated with WAVID

   logical  :: iamin_CPLATMID(num_inst_atm)         ! pe associated with CPLATMID
   logical  :: iamin_CPLLNDID(num_inst_lnd)         ! pe associated with CPLLNDID
   logical  :: iamin_CPLICEID(num_inst_ice)         ! pe associated with CPLICEID
   logical  :: iamin_CPLOCNID(num_inst_ocn)         ! pe associated with CPLOCNID
   logical  :: iamin_CPLGLCID(num_inst_glc)         ! pe associated with CPLGLCID
   logical  :: iamin_CPLROFID(num_inst_rof)         ! pe associated with CPLROFID
   logical  :: iamin_CPLWAVID(num_inst_wav)         ! pe associated with CPLWAVID

   ! complist: list of comps on this pe
   ! allow enough room for names of all physical components + coupler, 
   ! where each string can be up to (max_inst_name_len+1) characters
   ! long (+1 allows for a space before each name)
   character(len=(seq_comm_namelen+1)*(num_inst_phys+1)) :: complist

   integer  :: iam_GLOID              ! pe number in global id
   integer, pointer :: atm_petlist(:), lnd_petlist(:), ice_petlist(:), ocn_petlist(:), &
                       glc_petlist(:), cpl_petlist(:), rof_petlist(:), wav_petlist(:)

   integer, parameter :: ens1=1        ! use first instance of ensemble only
   integer, parameter :: fix1=1        ! temporary hard-coding to first ensemble, needs to be fixed
   integer :: eai, eli, eoi, eii, egi, eri, ewi, exi, efi  ! component instance counters
   character(len=seq_comm_namelen) :: atm_name(num_inst_atm)  !! For holding component instance names
   character(len=seq_comm_namelen) :: lnd_name(num_inst_lnd)
   character(len=seq_comm_namelen) :: ocn_name(num_inst_ocn)
   character(len=seq_comm_namelen) :: ice_name(num_inst_ice)
   character(len=seq_comm_namelen) :: glc_name(num_inst_glc)
   character(len=seq_comm_namelen) :: rof_name(num_inst_rof)
   character(len=seq_comm_namelen) :: wav_name(num_inst_wav)
   character(CL)                   :: atm_suffix(num_inst_atm) !! for holding per-instance suffix
   character(CL)                   :: lnd_suffix(num_inst_lnd)
   character(CL)                   :: ocn_suffix(num_inst_ocn)
   character(CL)                   :: ice_suffix(num_inst_ice)
   character(CL)                   :: glc_suffix(num_inst_glc)
   character(CL)                   :: rof_suffix(num_inst_rof)
   character(CL)                   :: wav_suffix(num_inst_wav)

   !----------------------------------------------------------------------------
   ! formats
   !----------------------------------------------------------------------------
   character(*), parameter :: subname = '(seq_mct_drv)'
   character(*), parameter :: F00 = "('"//subname//" : ', 4A )"
   character(*), parameter :: F0L = "('"//subname//" : ', A, L6 )"
   character(*), parameter :: F0I = "('"//subname//" : ', A, 2i8 )"
   character(*), parameter :: F0R = "('"//subname//" : ', A, 2g23.15 )"
   character(*), parameter :: FormatA = '(A,": =============== ", A41,          " ===============")'
   character(*), parameter :: FormatD = '(A,": =============== ", A20,2I8,5x,   " ===============")'
   character(*), parameter :: FormatR = '(A,": =============== ", A31,F9.3,1x,  " ===============")'
   character(*), parameter :: FormatQ = '(A,": =============== ", A20,2F10.2,1x," ===============")'

   save

!===============================================================================
contains
!===============================================================================

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

subroutine ccsm_pre_init()
   use shr_pio_mod, only : shr_pio_init1, shr_pio_init2
   implicit none
   !--------------------------------------------------------------------------
   ! Initialize MCT and MPI communicators and IO
   !--------------------------------------------------------------------------

   integer, dimension(num_inst_total) :: comp_id, comp_comm, comp_comm_iam
   logical :: comp_iamin(num_inst_total)
   character(len=seq_comm_namelen) :: comp_name(num_inst_total)
   integer :: i, it

   call mpi_init(ierr)
   call shr_mpi_chkerr(ierr,subname//' mpi_init')

   Global_Comm=MPI_COMM_WORLD
   comp_comm = MPI_COMM_NULL

   call shr_pio_init1(num_inst_total,NLFileName, Global_Comm)
!
! If pio_async_interface is true Global_Comm is MPI_COMM_NULL on the servernodes
! and server nodes do not return from shr_pio_init2
!
!   if (Global_Comm /= MPI_COMM_NULL) then
      call seq_comm_init(Global_Comm, NLFileName)

      !--- set task based threading counts ---
      call seq_comm_setptrs(GLOID,pethreads=pethreads_GLOID,iam=iam_GLOID)
      call seq_comm_setnthreads(pethreads_GLOID)

      !--- get some general data ---
      call seq_comm_setptrs(GLOID,mpicom=mpicom_GLOID,iamroot=iamroot_GLOID,nthreads=nthreads_GLOID)
      if (iamroot_GLOID) output_perf = .true.
      it=1
      call seq_comm_setptrs(CPLID,mpicom=mpicom_CPLID,iamroot=iamroot_CPLID,nthreads=nthreads_CPLID,iam=comp_comm_iam(it))
      if (iamroot_CPLID) output_perf = .true.
      comp_id(it) = CPLID
      comp_comm(it) = mpicom_CPLID
      iamin_CPLID    = seq_comm_iamin(CPLID)
      comp_iamin(it) =  seq_comm_iamin(comp_id(it))
      comp_name(it) = seq_comm_name(comp_id(it))
      complist = " "
      if (iamin_CPLID) complist = trim(complist)//' cpl'
      
      do eai = 1,num_inst_atm
         it=it+1
         comp_id(it) = ATMID(eai)
         comp_iamin(it) =  seq_comm_iamin(comp_id(it))
         comp_name(it) = seq_comm_name(comp_id(it))
         call seq_comm_setptrs(ATMID(eai), &
              mpicom=mpicom_ATMID(eai), &
              iamroot=iamroot_ATMID(eai), &
              nthreads=nthreads_ATMID, iam=comp_comm_iam(it))
         if (iamroot_ATMID(eai)) output_perf = .true.
         iamin_ATMID(eai) = seq_comm_iamin (ATMID(eai))
         atm_name   (eai) = seq_comm_name  (ATMID(eai))
         atm_suffix (eai) = seq_comm_suffix(ATMID(eai))
         comp_comm(it) = mpicom_ATMID(eai)
         if (iamin_ATMID(eai)) then
            complist = trim(complist)//' '//trim(atm_name(eai))
         end if

         call seq_comm_setptrs(CPLATMID(eai), &
              mpicom=mpicom_CPLATMID(eai), &
              nthreads=nthreads_CPLATMID)
         iamin_CPLATMID(eai) = seq_comm_iamin(CPLATMID(eai))
      enddo
      call seq_comm_setptrs(CPLALLATMID, mpicom=mpicom_CPLALLATMID)
      iamin_CPLALLATMID = seq_comm_iamin(CPLALLATMID)

      do eli = 1,num_inst_lnd
         it=it+1
         comp_id(it) = LNDID(eli)
         comp_iamin(it) =  seq_comm_iamin(comp_id(it))
         comp_name(it) = seq_comm_name(comp_id(it))
         call seq_comm_setptrs(LNDID(eli), &
              mpicom=mpicom_LNDID(eli), &
              iamroot=iamroot_LNDID(eli), &
              nthreads=nthreads_LNDID, iam=comp_comm_iam(it))
         if (iamroot_LNDID(eli)) output_perf = .true.
         iamin_LNDID(eli) = seq_comm_iamin (LNDID(eli))
         lnd_name   (eli) = seq_comm_name  (LNDID(eli))
         lnd_suffix (eli) = seq_comm_suffix(LNDID(eli))
         comp_comm(it)=mpicom_lndid(eli)
         if (iamin_LNDID(eli)) then
            complist = trim(complist)//' '//trim(lnd_name(eli))
         end if

         call seq_comm_setptrs(CPLLNDID(eli), &
              mpicom=mpicom_CPLLNDID(eli), &
              nthreads=nthreads_CPLLNDID)
         iamin_CPLLNDID(eli) = seq_comm_iamin(CPLLNDID(eli))
      enddo
      call seq_comm_setptrs(CPLALLLNDID, mpicom=mpicom_CPLALLLNDID)
      iamin_CPLALLLNDID = seq_comm_iamin(CPLALLLNDID)

      do eoi = 1,num_inst_ocn
         it=it+1
         comp_id(it) = OCNID(eoi)
         comp_iamin(it) =  seq_comm_iamin(comp_id(it))
         comp_name(it) = seq_comm_name(comp_id(it))
         call seq_comm_setptrs(OCNID(eoi), &
              mpicom=mpicom_OCNID(eoi), &
              iamroot=iamroot_OCNID(eoi), &
              nthreads=nthreads_OCNID, &
              iam=comp_comm_iam(it))
         if (iamroot_OCNID(eoi)) output_perf = .true.
         iamin_OCNID(eoi) = seq_comm_iamin (OCNID(eoi))
         ocn_name   (eoi) = seq_comm_name  (OCNID(eoi))
         ocn_suffix (eoi) = seq_comm_suffix(OCNID(eoi))
         comp_comm(it) = mpicom_ocnid(eoi)
         if (iamin_OCNID(eoi)) then
            complist = trim(complist)//' '//trim(ocn_name(eoi))
         end if
         
         call seq_comm_setptrs(CPLOCNID(eoi), &
              mpicom=mpicom_CPLOCNID(eoi), &
              nthreads=nthreads_CPLOCNID)
         iamin_CPLOCNID(eoi) = seq_comm_iamin(CPLOCNID(eoi))
      enddo
      call seq_comm_setptrs(CPLALLOCNID, mpicom=mpicom_CPLALLOCNID)
      iamin_CPLALLOCNID = seq_comm_iamin(CPLALLOCNID)

      do eii = 1,num_inst_ice
         it=it+1
         comp_id(it) = ICEID(eii)
         comp_iamin(it) =  seq_comm_iamin(comp_id(it))
         comp_name(it) = seq_comm_name(comp_id(it))
         call seq_comm_setptrs(ICEID(eii), &
              mpicom=mpicom_ICEID(eii), &
              iamroot=iamroot_ICEID(eii), &
              nthreads=nthreads_ICEID, &
              iam=comp_comm_iam(it))
         if (iamroot_ICEID(eii)) output_perf = .true.
         iamin_ICEID(eii) = seq_comm_iamin (ICEID(eii))
         ice_name   (eii) = seq_comm_name  (ICEID(eii))
         ice_suffix (eii) = seq_comm_suffix(ICEID(eii))
         comp_comm(it) = mpicom_iceid(eii)
         if (iamin_ICEID(eii)) then
            complist = trim(complist)//' '//trim(ice_name(eii))
         end if

         call seq_comm_setptrs(CPLICEID(eii), &
              mpicom=mpicom_CPLICEID(eii), &
              nthreads=nthreads_CPLICEID)
         iamin_CPLICEID(eii) = seq_comm_iamin(CPLICEID(eii))
      enddo
      call seq_comm_setptrs(CPLALLICEID, mpicom=mpicom_CPLALLICEID)
      iamin_CPLALLICEID = seq_comm_iamin(CPLALLICEID)
      
      do egi = 1,num_inst_glc
         it=it+1
         comp_id(it) = GLCID(egi)
         comp_iamin(it) =  seq_comm_iamin(comp_id(it))
         comp_name(it) = seq_comm_name(comp_id(it))
         call seq_comm_setptrs(GLCID(egi), &
              mpicom=mpicom_GLCID(egi), &
              iamroot=iamroot_GLCID(egi), &
              nthreads=nthreads_GLCID, &
              iam=comp_comm_iam(it))
         if (iamroot_GLCID(egi)) output_perf = .true.
         comp_comm(it) = mpicom_glcid(egi)
         iamin_GLCID(egi) = seq_comm_iamin (GLCID(egi))
         glc_name   (egi) = seq_comm_name  (GLCID(egi))
         glc_suffix (egi) = seq_comm_suffix(GLCID(egi))
         if (iamin_GLCID(egi)) then
            complist = trim(complist)//' '//trim(glc_name(egi))
         end if

         call seq_comm_setptrs(CPLGLCID(egi), &
              mpicom=mpicom_CPLGLCID(egi), &
              nthreads=nthreads_CPLGLCID)
         iamin_CPLGLCID(egi) = seq_comm_iamin(CPLGLCID(egi))
      enddo
      call seq_comm_setptrs(CPLALLGLCID, mpicom=mpicom_CPLALLGLCID)
      iamin_CPLALLGLCID = seq_comm_iamin(CPLALLGLCID)

      do eri = 1,num_inst_rof
         it=it+1
         comp_id(it) = ROFID(eri)
         comp_iamin(it) = seq_comm_iamin(comp_id(it))
         comp_name(it) = seq_comm_name(comp_id(it))
         call seq_comm_setptrs(ROFID(eri), &
              mpicom=mpicom_ROFID(eri), &
              iamroot=iamroot_ROFID(eri), &
              nthreads=nthreads_ROFID, iam=comp_comm_iam(it))
         if (iamroot_ROFID(eri)) output_perf = .true.
         iamin_ROFID(eri) = seq_comm_iamin (ROFID(eri))
         rof_name   (eri) = seq_comm_name  (ROFID(eri))
         rof_suffix (eri) = seq_comm_suffix(ROFID(eri))
         comp_comm(it)=mpicom_rofid(eri)
         if (iamin_ROFID(eri)) then
            complist = trim(complist)//' '//trim(rof_name(eri))
         end if

         call seq_comm_setptrs(CPLROFID(eri), &
              mpicom=mpicom_CPLROFID(eri), &
              nthreads=nthreads_CPLROFID)
         iamin_CPLROFID(eri) = seq_comm_iamin(CPLROFID(eri))
      enddo
      call seq_comm_setptrs(CPLALLROFID, mpicom=mpicom_CPLALLROFID)
      iamin_CPLALLROFID = seq_comm_iamin(CPLALLROFID)

      do ewi = 1,num_inst_wav
         it=it+1
         comp_id(it) = WAVID(ewi)
         comp_iamin(it) =  seq_comm_iamin(comp_id(it))
         comp_name(it) = seq_comm_name(comp_id(it))
         call seq_comm_setptrs(WAVID(ewi), &
              mpicom=mpicom_WAVID(ewi), &
              iamroot=iamroot_WAVID(ewi), &
              nthreads=nthreads_WAVID, &
              iam=comp_comm_iam(it))
         if (iamroot_WAVID(ewi)) output_perf = .true.
         comp_comm(it) = mpicom_wavid(ewi)
         iamin_WAVID(ewi) = seq_comm_iamin (WAVID(ewi))
         wav_name   (ewi) = seq_comm_name  (WAVID(ewi))
         wav_suffix (ewi) = seq_comm_suffix(WAVID(ewi))
         if (iamin_WAVID(ewi)) then
            complist = trim(complist)//' '//trim(wav_name(ewi))
         end if

         call seq_comm_setptrs(CPLWAVID(ewi), &
              mpicom=mpicom_CPLWAVID(ewi), &
              nthreads=nthreads_CPLWAVID)
         iamin_CPLWAVID(ewi) = seq_comm_iamin(CPLWAVID(ewi))
      enddo
      call seq_comm_setptrs(CPLALLWAVID, mpicom=mpicom_CPLALLWAVID)
      iamin_CPLALLWAVID = seq_comm_iamin(CPLALLWAVID)

   !--------------------------------------------------------------------------
   ! Set logging parameters both for shr code and locally
   !--------------------------------------------------------------------------

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

   !--------------------------------------------------------------------------
   ! Log info about the environment settings
   !--------------------------------------------------------------------------

   if (iamroot_CPLID) then
#ifdef USE_ESMF_LIB
      write(logunit,'(2A)') subname,' USE_ESMF_LIB is set'
#else
      write(logunit,'(2A)') subname,' USE_ESMF_LIB is NOT set, using esmf_wrf_timemgr'
#endif
#ifdef MCT_INTERFACE
      write(logunit,'(2A)') subname,' MCT_INTERFACE is set'
#endif
#ifdef ESMF_INTERFACE
      write(logunit,'(2A)') subname,' ESMF_INTERFACE is set'
#endif
   endif
!
!  When using io servers (pio_async_interface=.true.) the server tasks do not return from 
!  shr_pio_init2 
!

   call shr_pio_init2(comp_id,comp_name,comp_iamin,comp_comm,comp_comm_iam)



end subroutine ccsm_pre_init

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

subroutine ccsm_init()

   implicit none

 101  format( A, 2i8, 12A, A, F8.2, A, F8.2 )
 102  format( A, 2i8, A, 8L3 )
 103  format( 5A )
 104  format( A, 2i8)
 105  format( A, 2i8, A, f10.2, A, f10.2, A, A, i5, A, A)
 106  format( A, f23.12)

   !--------------------------------------------------------------------------
   ! Print Model heading and copyright message
   !--------------------------------------------------------------------------

   if (iamroot_CPLID) call seq_ccsm_printlogheader()

   !-----------------------------------------------------------------------------
   ! Timer initialization (has to be after mpi init)
   !-----------------------------------------------------------------------------

   call t_initf(NLFileName, LogPrint=.false., mpicom=mpicom_GLOID, &
                MasterTask=iamroot_GLOID)

   if (iamin_CPLID) then
      call seq_io_cpl_init()
   endif

   call t_startf('DRIVER_INIT')

   !-----------------------------------------------------------------------------
   ! Memory test
   !-----------------------------------------------------------------------------
!   call shr_mem_init(prt=.true.)
    call shr_mem_init(prt=iamroot_CPLID)

   !-----------------------------------------------------------------------------
   ! Initialize coupled fields
   !-----------------------------------------------------------------------------

   call seq_flds_set(nlfilename,GLOID)

   !-----------------------------------------------------------------------------
   ! Initialize infodata
   !-----------------------------------------------------------------------------

   call seq_infodata_init(infodata,nlfilename,GLOID)
   if (iamroot_CPLID) then
      write(logunit,*) ' '
      write(logunit,'(2A)') 'Status of infodata after seq_infodata_init'
      call seq_infodata_print( infodata )
      write(logunit,*) ' '
   endif

   call seq_infodata_GetData(infodata,read_restart=read_restart, restart_file=rest_file, &
        timing_dir=timing_dir, tchkpt_dir=tchkpt_dir)
   call seq_infodata_GetData(infodata, info_debug=info_debug, atm_present=atm_present, &
        lnd_present=lnd_present, ice_present=ice_present, ocn_present=ocn_present, &
        glc_present=glc_present, sno_present=sno_present, rof_present=rof_present, &
        wav_present=wav_present, &
        single_column=single_column, aqua_planet=aqua_planet, &
        ocean_tight_coupling=ocean_tight_coupling, drv_threading=drv_threading)
   call seq_infodata_GetData(infodata, do_histinit=do_histinit)
   call seq_infodata_GetData(infodata, do_budgets=do_budgets, budget_inst=budget_inst, &
        budget_daily=budget_daily, budget_month=budget_month, budget_ann=budget_ann, &
        budget_ltann=budget_ltann, budget_ltend=budget_ltend)
   call seq_infodata_GetData(infodata, &
        histaux_a2x    =do_hist_a2x    , histaux_a2x3hr =do_hist_a2x3hr , &
        histaux_a2x3hrp=do_hist_a2x3hrp, histaux_a2x24hr=do_hist_a2x24hr, &
        histaux_l2x    =do_hist_l2x    , histaux_r2x    =do_hist_r2x,     &
        histaux_s2x1yr=do_hist_s2x1yr      )
   call seq_infodata_GetData(infodata, run_barriers = run_barriers)
   call seq_infodata_GetData(infodata, mct_usealltoall=mct_usealltoall, &
        mct_usevector=mct_usevector)

   call seq_infodata_GetData(infodata, aoflux_grid=aoflux_grid, vect_map=vect_map)
   call seq_infodata_GetData(infodata, samegrid_ao=samegrid_ao, samegrid_al=samegrid_al, &
        samegrid_ro=samegrid_ro, samegrid_aw=samegrid_aw, samegrid_ow=samegrid_ow)
   ! pass the cpl_decomp value to seq_mctext_decomp
   call seq_infodata_GetData(infodata, cpl_decomp = seq_mctext_decomp)

   call seq_infodata_GetData(infodata, shr_map_dopole=shr_map_dopole)
   call shr_map_setDopole(shr_map_dopole)

   call seq_infodata_GetData(infodata, reprosum_use_ddpdd = reprosum_use_ddpdd, &
        reprosum_diffmax = reprosum_diffmax, reprosum_recompute = reprosum_recompute)
   call shr_reprosum_setopts(repro_sum_use_ddpdd_in = reprosum_use_ddpdd, &
        repro_sum_rel_diff_max_in = reprosum_diffmax, &
        repro_sum_recompute_in = reprosum_recompute)

   !-----------------------------------------------------------------------------
   ! Test Threading Setup in driver, happens to be valid on all pes for all IDs
   !-----------------------------------------------------------------------------

   if (drv_threading) then
      if (iamroot_GLOID) write(logunit,*) ' '
      if (iamroot_GLOID) write(logunit,'(2A)    ') subname,' Test Threading in driver'
      call seq_comm_setnthreads(nthreads_GLOID)
      if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_GLOID = ',nthreads_GLOID,seq_comm_getnthreads()
      call seq_comm_setnthreads(nthreads_CPLID)
      if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_CPLID = ',nthreads_CPLID,seq_comm_getnthreads()
      call seq_comm_setnthreads(nthreads_ATMID)
      if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_ATMID = ',nthreads_ATMID,seq_comm_getnthreads()
      call seq_comm_setnthreads(nthreads_LNDID)
      if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_LNDID = ',nthreads_LNDID,seq_comm_getnthreads()
      call seq_comm_setnthreads(nthreads_OCNID)
      if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_OCNID = ',nthreads_OCNID,seq_comm_getnthreads()
      call seq_comm_setnthreads(nthreads_ICEID)
      if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_ICEID = ',nthreads_ICEID,seq_comm_getnthreads()
      call seq_comm_setnthreads(nthreads_GLCID)
      if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_GLCID = ',nthreads_GLCID,seq_comm_getnthreads()
      call seq_comm_setnthreads(nthreads_ROFID)
      if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_ROFID = ',nthreads_ROFID,seq_comm_getnthreads()
      call seq_comm_setnthreads(nthreads_WAVID)
      if (iamroot_GLOID) write(logunit,'(2A,2I4)') subname,'    nthreads_WAVID = ',nthreads_WAVID,seq_comm_getnthreads()
      if (iamroot_GLOID) write(logunit,*) ' '

      call seq_comm_setnthreads(nthreads_GLOID)
   endif

   !-----------------------------------------------------------------------------
   ! Setup cdata types, call on all pes so the ID is set on all pes even 
   ! though other data may be invalid
   !-----------------------------------------------------------------------------

   call seq_cdata_init(cdata_ax, CPLID, dom_ax, gsMap_ax, infodata, name='cdata_ax' )
   call seq_cdata_init(cdata_lx, CPLID, dom_lx, gsMap_lx, infodata, name='cdata_lx' )
   call seq_cdata_init(cdata_rx, CPLID, dom_rx, gsMap_rx, infodata, name='cdata_rx' )
   call seq_cdata_init(cdata_ix, CPLID, dom_ix, gsMap_ix, infodata, name='cdata_ix' )
   call seq_cdata_init(cdata_ox, CPLID, dom_ox, gsMap_ox, infodata, name='cdata_ox' )
   call seq_cdata_init(cdata_gx, CPLID, dom_gx, gsMap_gx, infodata, name='cdata_gx' )
   call seq_cdata_init(cdata_sx, CPLID, dom_sx, gsMap_sx, infodata, name='cdata_sx' )
   call seq_cdata_init(cdata_wx, CPLID, dom_wx, gsMap_wx, infodata, name='cdata_wx' )

   do eai = 1,num_inst_atm
      call seq_cdata_init(cdata_aa(eai), ATMID(eai), &
                          dom_aa(eai), gsmap_aa(eai), infodata, &
                          name='cdata_aa'//trim(atm_name(eai)))
   enddo

   do eli = 1,num_inst_lnd
      call seq_cdata_init(cdata_ll(eli), LNDID(eli), &
                          dom_ll(eli), gsMap_ll(eli), infodata, &
                          name=('cdata_ll' // trim(lnd_name(eli))))
      call seq_cdata_init(cdata_ss(eli), LNDID(eli), &
                          dom_ss(eli), gsMap_ss(eli), infodata, &
                          name=('cdata_ss' // trim(lnd_name(eli))))
   enddo

   do eri = 1,num_inst_rof
      call seq_cdata_init(cdata_rr(eri), ROFID(eri), &
                          dom_rr(eri), gsMap_rr(eri), infodata, &
                          name=('cdata_rr' // trim(rof_name(eri))))
   end do

   do eii = 1,num_inst_ice
      call seq_cdata_init(cdata_ii(eii), ICEID(eii), &
                          dom_ii(eii), gsmap_ii(eii), infodata, &
                          name='cdata_ii'//trim(ice_name(eii)))
   enddo

   do eoi = 1,num_inst_ocn
      call seq_cdata_init(cdata_oo(eoi), OCNID(eoi), &
                          dom_oo(eoi), gsmap_oo(eoi), infodata, &
                          name='cdata_oo'//trim(ocn_name(eoi)))
   enddo

   do egi = 1,num_inst_glc
      call seq_cdata_init(cdata_gg(egi), GLCID(egi), &
                          dom_gg(egi), gsmap_gg(egi), infodata, &
                          name='cdata_gg'//trim(glc_name(egi)))
   enddo

   do ewi = 1,num_inst_wav
      call seq_cdata_init(cdata_ww(ewi), WAVID(ewi), &
                          dom_ww(ewi), gsmap_ww(ewi), infodata, &
                          name='cdata_ww'//trim(wav_name(ewi)))
   enddo

   !-----------------------------------------------------------------------------
   ! Initialize time manager
   !-----------------------------------------------------------------------------

   call seq_timemgr_clockInit(seq_SyncClock,nlfilename,read_restart,rest_file,mpicom_gloid, &
        EClock_d, EClock_a, EClock_l, EClock_o, EClock_i, Eclock_g, Eclock_r, Eclock_w)
   if (iamroot_CPLID) then
       call seq_timemgr_clockPrint(seq_SyncClock)
   endif

   call seq_infodata_getData(infodata,orb_iyear=orb_iyear,orb_iyear_align=orb_iyear_align, &
      orb_mode=orb_mode)
   if (trim(orb_mode) == trim(seq_infodata_orb_variable_year)) then
      call seq_timemgr_EClockGetData( EClock_d, curr_ymd=ymd)
      call shr_cal_date2ymd(ymd,year,month,day)
      orb_cyear = orb_iyear + (year - orb_iyear_align)
      call shr_orb_params(orb_cyear, orb_eccen, orb_obliq, orb_mvelp, &
                          orb_obliqr, orb_lambm0, orb_mvelpp, iamroot_CPLID)
      call seq_infodata_putData(infodata,orb_eccen=orb_eccen,orb_obliqr=orb_obliqr, &
           orb_lambm0=orb_lambm0,orb_mvelpp=orb_mvelpp)
   endif

   call seq_infodata_putData(infodata,atm_phase=1,lnd_phase=1,ocn_phase=1,ice_phase=1, &
                                      glc_phase=1,wav_phase=1)

   !-----------------------------------------------------------------------------
   ! If in single column mode, overwrite flags according to focndomain file
   ! in ocn_in namelist. SCAM can reset the "present" flags for lnd, sno,
   ! ocn, ice, rof, and flood.
   !-----------------------------------------------------------------------------

   if (.not.aqua_planet .and. single_column) then
      call seq_infodata_getData( infodata, scmlon=scmlon, scmlat=scmlat)
      call shr_scam_checkSurface(scmlon, scmlat, OCNID(ens1), mpicom_OCNID(ens1), &
           lnd_present=lnd_present, sno_present=sno_present, &
           ocn_present=ocn_present, ice_present=ice_present, &
           rof_present=rof_present, flood_present=flood_present)
      call seq_infodata_putData( infodata, &
           lnd_present=lnd_present, sno_present=sno_present, &
           ocn_present=ocn_present, ice_present=ice_present, &
           rof_present=rof_present, flood_present=flood_present)
   endif

   !-----------------------------------------------------------------------------
   ! Component Initialization
   ! Note that within each component initialization, the relevant x_pxresent flag 
   ! part of CCSMInit (contained as a pointer in cdata_xc) can be modified
   ! By default, all these flags are set to true
   ! The atm can reset the lnd_present, ice_present and ocn_present flags based
   ! on aqua_planet, ideal_phys and adiabatic modes
   ! The stub components will reset the present flags to false, all other
   ! components will set them to true for the purposes of symmetry
   !-----------------------------------------------------------------------------

   call t_startf('driver_init_comps')
   if (iamroot_CPLID )then
      write(logunit,*) ' '
      write(logunit,F00) 'Initialize each component: atm, lnd, rof, ocn, ice, glc, wav'
      call shr_sys_flush(logunit)
   endif

   call t_adj_detailf(+2)

   !-----------------------------------------------------------------------------
   ! Initialization atmospheric component
   !-----------------------------------------------------------------------------

   if (iamin_CPLALLATMID) then
      call seq_infodata_exchange(infodata,CPLALLATMID,'cpl2atm_init')
   endif
   do eai = 1,num_inst_atm
      if (iamroot_CPLID) then
         write(logunit,F00) 'Initialize atm component '//trim(atm_name(eai))
         call shr_sys_flush(logunit)
      endif
      if (iamin_ATMID(eai) .and. atm_present) then
         if (drv_threading) call seq_comm_setnthreads(nthreads_ATMID)
         call shr_sys_flush(logunit)
         call atm_init_mct( EClock_a, cdata_aa(eai), x2a_aa(eai), a2x_aa(eai), NLFilename=NLFilename )
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
      endif
   enddo
   if (iamin_CPLATMID(ens1)) then
      call seq_infodata_exchange(infodata,CPLATMID(ens1),'atm2cpl_init')
   endif

   !-----------------------------------------------------------------------------
   ! Initialization land component
   !-----------------------------------------------------------------------------

   if (iamin_CPLALLLNDID) then
      call seq_infodata_exchange(infodata,CPLALLLNDID,'cpl2lnd_init')
   endif
   do eli = 1,num_inst_lnd
      if (iamroot_CPLID) then
         write(logunit,F00) 'Initialize lnd component '//trim(lnd_name(eli))
         call shr_sys_flush(logunit)
      endif
      if (iamin_LNDID(eli) .and. lnd_present) then
         if (drv_threading) call seq_comm_setnthreads(nthreads_LNDID)
         call shr_sys_flush(logunit)
         call lnd_init_mct( EClock_l, &
                            cdata_ll(eli), x2l_ll(eli), l2x_ll(eli), &
                            cdata_ss(eli), x2s_ss(eli), s2x_ss(eli), &
                            NLFilename=NLFilename )
      end if
   end do
   if (iamin_CPLLNDID(ens1)) then
      call seq_infodata_exchange(infodata,CPLLNDID(ens1),'lnd2cpl_init')
   endif

   !----------------------------------------------------
   ! Initialization river runoff component 
   !----------------------------------------------------

   if (iamin_CPLALLROFID) then
      call seq_infodata_exchange(infodata,CPLALLROFID,'cpl2rof_init')
   endif
   do eri = 1,num_inst_rof
      if (iamroot_CPLID) then
         write(logunit,F00) 'Initialize rof component '//trim(rof_name(eri))
         call shr_sys_flush(logunit)
      endif
      if (iamin_ROFID(eri) .and. rof_present) then
         if (drv_threading) call seq_comm_setnthreads(nthreads_ROFID)
         call rof_init_mct( EClock_r, &
                            cdata_rr(eri),  x2r_rr(eri), r2x_rr(eri), &
                            NLFilename=NLFilename)
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
      endif
   enddo
   if (iamin_CPLROFID(ens1)) then
      call seq_infodata_exchange(infodata,CPLROFID(ens1),'rof2cpl_init')
   endif

   !-----------------------------------------------------------------------------
   ! Initialization ocean component
   !-----------------------------------------------------------------------------

   if (iamin_CPLALLOCNID) then
      call seq_infodata_exchange(infodata,CPLALLOCNID,'cpl2ocn_init')
   endif
   do eoi = 1,num_inst_ocn
      if (iamroot_CPLID) then
         write(logunit,F00) 'Initialize ocn component '//trim(ocn_name(eoi))
         call shr_sys_flush(logunit)
      endif
      if (iamin_OCNID(eoi) .and. ocn_present) then
         if (drv_threading) call seq_comm_setnthreads(nthreads_OCNID)
         call shr_sys_flush(logunit)
         call ocn_init_mct( EClock_o, &
                            cdata_oo(eoi), x2o_oo(eoi), o2x_oo(eoi), &
                            NLFilename=NLFilename )
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
      endif
   enddo
   if (iamin_CPLOCNID(ens1)) then
      call seq_infodata_exchange(infodata,CPLOCNID(ens1),'ocn2cpl_init')
   endif

   !-----------------------------------------------------------------------------
   ! Initialization ice component
   !-----------------------------------------------------------------------------

   if (iamin_CPLALLICEID) then
      call seq_infodata_exchange(infodata,CPLALLICEID,'cpl2ice_init')
   endif
   do eii = 1,num_inst_ice
      if (iamroot_CPLID) then
         write(logunit,F00) 'Initialize ice component '//trim(ice_name(eii))
         call shr_sys_flush(logunit)
      endif
      if (iamin_ICEID(eii) .and. ice_present) then
         if (drv_threading) call seq_comm_setnthreads(nthreads_ICEID)
         call shr_sys_flush(logunit)
         call ice_init_mct( EClock_i, &
                            cdata_ii(eii), x2i_ii(eii), i2x_ii(eii), &
                            NLFilename=NLFilename )
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
      endif
   enddo
   if (iamin_CPLICEID(ens1)) then
      call seq_infodata_exchange(infodata,CPLICEID(ens1),'ice2cpl_init')
   endif

   !-----------------------------------------------------------------------------
   ! Initialization glc component
   !-----------------------------------------------------------------------------

   if (iamin_CPLALLGLCID) then
      call seq_infodata_exchange(infodata,CPLALLGLCID,'cpl2glc_init')
   endif
   do egi = 1,num_inst_glc
      if (iamroot_CPLID) then
         write(logunit,F00) 'Initialize glc component '//trim(glc_name(egi))
         call shr_sys_flush(logunit)
      endif
      if (iamin_GLCID(egi) .and. glc_present) then
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLCID)
         call shr_sys_flush(logunit)
         call glc_init_mct( EClock_g, &
                            cdata_gg(egi), x2g_gg(egi), g2x_gg(egi), &
                            NLFilename=NLFilename )
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
      endif
   enddo
   if (iamin_CPLGLCID(ens1)) then
      call seq_infodata_exchange(infodata,CPLGLCID(ens1),'glc2cpl_init')
   endif

   !-----------------------------------------------------------------------------
   ! Initialization wav component
   !-----------------------------------------------------------------------------

   if (iamin_CPLALLWAVID) then
      call seq_infodata_exchange(infodata,CPLALLWAVID,'cpl2wav_init')
   endif
   do ewi = 1,num_inst_wav
      if (iamroot_CPLID) then
         write(logunit,F00) 'Initialize wav component '//trim(wav_name(ewi))
         call shr_sys_flush(logunit)
      endif
      if (iamin_WAVID(ewi) .and. wav_present) then
         if (drv_threading) call seq_comm_setnthreads(nthreads_WAVID)
         call shr_sys_flush(logunit)
         call wav_init_mct( EClock_w, &
                            cdata_ww(ewi), x2w_ww(ewi), w2x_ww(ewi), &
                            NLFilename=NLFilename )
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
      endif
   enddo
   if (iamin_CPLWAVID(ens1)) then
      call seq_infodata_exchange(infodata,CPLWAVID(ens1),'wav2cpl_init')
   endif

   !-----------------------------------------------------------------------------

   call t_adj_detailf(-2)

   call t_stopf  ('driver_init_comps')

   !-----------------------------------------------------------------------------
   ! Determine final settings for presence of land, ice and ocean and the prognostic flags
   !-----------------------------------------------------------------------------

   if (iamin_CPLALLATMID) then
      call seq_infodata_exchange(infodata,CPLALLATMID,'cpl2atm_init')
   endif
   if (iamin_CPLALLLNDID) then
      call seq_infodata_exchange(infodata,CPLALLLNDID,'cpl2lnd_init')
   endif
   if (iamin_CPLALLOCNID) then
      call seq_infodata_exchange(infodata,CPLALLOCNID,'cpl2ocn_init')
   endif
   if (iamin_CPLALLICEID) then
      call seq_infodata_exchange(infodata,CPLALLICEID,'cpl2ice_init')
   endif
   if (iamin_CPLALLGLCID) then
      call seq_infodata_exchange(infodata,CPLALLGLCID,'cpl2glc_init')
   endif
   if (iamin_CPLALLROFID) then
      call seq_infodata_exchange(infodata,CPLALLROFID,'cpl2rof_init')
   endif
   if (iamin_CPLALLWAVID) then
      call seq_infodata_exchange(infodata,CPLALLWAVID,'cpl2wav_init')
   endif

   if (iamroot_CPLID) then
      write(logunit,F00) 'Determine final settings for presence of surface components'
      call shr_sys_flush(logunit)
   endif

   call seq_infodata_getData(infodata, &
        atm_present=atm_present, &
        lnd_present=lnd_present, &
        ice_present=ice_present, &
        ocn_present=ocn_present, & 
        glc_present=glc_present, & 
        rof_present=rof_present, &
        wav_present=wav_present, & 
        flood_present=flood_present, &
        sno_present=sno_present, & 
        atm_prognostic=atm_prognostic, &
        lnd_prognostic=lnd_prognostic, &
        ice_prognostic=ice_prognostic, &
        ocn_prognostic=ocn_prognostic, &
        ocnrof_prognostic=ocnrof_prognostic, &
        glc_prognostic=glc_prognostic, &
        rof_prognostic=rof_prognostic, &
        sno_prognostic=sno_prognostic, &
        wav_prognostic=wav_prognostic, &
        dead_comps=dead_comps, &
        esmf_map_flag=esmf_map_flag, &
        atm_nx=atm_nx, atm_ny=atm_ny, &
        lnd_nx=lnd_nx, lnd_ny=lnd_ny, &
        rof_nx=rof_nx, rof_ny=rof_ny, &
        ice_nx=ice_nx, ice_ny=ice_ny, &
        glc_nx=glc_nx, glc_ny=glc_ny, &
        sno_nx=sno_nx, sno_ny=sno_ny, &
        ocn_nx=ocn_nx, ocn_ny=ocn_ny, &
        wav_nx=wav_nx, wav_ny=wav_ny, &
        cpl_cdf64=cdf64, &
        atm_aero=atm_aero )

   if (ocnrof_prognostic .and. .not.rof_present) then
      if (iamroot_CPLID) then
         write(logunit,F00) 'WARNING: ocnrof_prognostic is TRUE but rof_present is FALSE'
         call shr_sys_flush(logunit)
      endif
   endif
   if (ocn_prognostic .and. .not.ocn_present) then
      call shr_sys_abort('if prognostic ocn must also have ocn present')
   endif
   if (lnd_prognostic .and. .not.lnd_present) then
      call shr_sys_abort('if prognostic lnd must also have lnd present')
   endif
   if (ice_prognostic .and. .not.ice_present) then
      call shr_sys_abort('if prognostic ice must also have ice present')
   endif
   if (glc_prognostic .and. .not.glc_present) then
      call shr_sys_abort('if prognostic glc must also have glc present')
   endif
   if (rof_prognostic .and. .not.rof_present) then
      call shr_sys_abort('if prognostic rof must also have rof present')
   endif
   if (sno_prognostic .and. .not.sno_present) then
      call shr_sys_abort('if prognostic sno must also have sno present')
   endif
   if (wav_prognostic .and. .not.wav_present) then
      call shr_sys_abort('if prognostic wav must also have wav present')
   endif
   if ((ice_prognostic .or. ocn_prognostic .or. lnd_prognostic) .and. .not. atm_present) then
      call shr_sys_abort('if prognostic surface model must also have atm present')
   endif
! tcx remove temporarily for development
!   if (glc_prognostic .and. .not.sno_present) then
!      call shr_sys_abort('if prognostic glc must also have sno present')
!   endif
!   if (sno_prognostic .and. .not.glc_present) then
!      call shr_sys_abort('if prognostic sno must also have glc present')!
!   endif

   ! Prognostic components must be consistent with num_inst_max for coupling 

   if (atm_prognostic .and. num_inst_atm /= num_inst_max) &
      call shr_sys_abort('atm_prognostic but num_inst_atm not num_inst_max')
   if (lnd_prognostic .and. num_inst_lnd /= num_inst_max) &
      call shr_sys_abort('lnd_prognostic but num_inst_lnd not num_inst_max')
   if (ocn_prognostic .and. num_inst_ocn /= num_inst_max) &
      call shr_sys_abort('ocn_prognostic but num_inst_ocn not num_inst_max')
   if (ice_prognostic .and. num_inst_ice /= num_inst_max) &
      call shr_sys_abort('ice_prognostic but num_inst_ice not num_inst_max')
   if (glc_prognostic .and. num_inst_glc /= num_inst_max) &
      call shr_sys_abort('glc_prognostic but num_inst_glc not num_inst_max')
   if (rof_prognostic .and. num_inst_rof /= num_inst_max) &
      call shr_sys_abort('rof_prognostic but num_inst_rof not num_inst_max')
   if (wav_prognostic .and. num_inst_wav /= num_inst_max) &
      call shr_sys_abort('wav_prognostic but num_inst_wav not num_inst_max')

   !-----------------------------------------------------------------------------
   ! Set domain check and other flag
   !-----------------------------------------------------------------------------

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

   !-----------------------------------------------------------------------------
   ! Write output
   ! NOTE- assume that runoff will only be mapped from land to ocean if 
   !       prognostic ocean is true
   !-----------------------------------------------------------------------------

   if (iamroot_CPLID) then
      write(logunit,*  )' '
      write(logunit,F00)'After component initialization:'
      write(logunit,F0L)'atm model present     = ',atm_present
      write(logunit,F0L)'lnd model present     = ',lnd_present
      write(logunit,F0L)'ocn model present     = ',ocn_present
      write(logunit,F0L)'ice model present     = ',ice_present
      write(logunit,F0L)'glc model present     = ',glc_present
      write(logunit,F0L)'sno model present     = ',sno_present
      write(logunit,F0L)'rof model present     = ',rof_present
      write(logunit,F0L)'rof/flood present     = ',flood_present
      write(logunit,F0L)'wav model present     = ',wav_present
      write(logunit,F0L)'atm model prognostic  = ',atm_prognostic
      write(logunit,F0L)'lnd model prognostic  = ',lnd_prognostic
      write(logunit,F0L)'ocn model prognostic  = ',ocn_prognostic
      write(logunit,F0L)'ice model prognostic  = ',ice_prognostic
      write(logunit,F0L)'glc model prognostic  = ',glc_prognostic
      write(logunit,F0L)'sno model prognostic  = ',sno_prognostic
      write(logunit,F0L)'rof model prognostic  = ',rof_prognostic
      write(logunit,F0L)'ocn rof   prognostic  = ',ocnrof_prognostic
      write(logunit,F0L)'wav model prognostic  = ',wav_prognostic
      write(logunit,F0L)'dead components       = ',dead_comps
      write(logunit,F0L)'domain_check          = ',domain_check
      write(logunit,F0I)'atm_nx,atm_ny         = ',atm_nx,atm_ny
      write(logunit,F0I)'lnd_nx,lnd_ny         = ',lnd_nx,lnd_ny
      write(logunit,F0I)'rof_nx,rof_ny         = ',rof_nx,rof_ny
      write(logunit,F0I)'ice_nx,ice_ny         = ',ice_nx,ice_ny
      write(logunit,F0I)'ocn_nx,ocn_ny         = ',ocn_nx,ocn_ny
      write(logunit,F0I)'glc_nx,glc_ny         = ',glc_nx,glc_ny
      write(logunit,F0I)'sno_nx,sno_ny         = ',sno_nx,sno_ny
      write(logunit,F0I)'wav_nx,wav_ny         = ',wav_nx,wav_ny
      write(logunit,F0L)'skip init ocean run   = ',skip_ocean_run
      write(logunit,F0L)'ocean tight coupling  = ',ocean_tight_coupling
      write(logunit,F0L)'cpl_cdf64             = ',cdf64
      write(logunit,F0L)'do_histavg            = ',do_histavg
      write(logunit,F0L)'atm_aero              = ',atm_aero
      write(logunit,*  )' '
      call shr_sys_flush(logunit)
   endif

   !-----------------------------------------------------------------------------
   ! Need to initialize aream, set it to area for now until maps are read
   !   in some cases, maps are not read at all !!
   ! Entire domain must have reasonable values before calling xxx2xxx init
   ! NOTE (tcx) : use cdata%dom instead of dom% due to seg fault on bluevista I, why?
   !-----------------------------------------------------------------------------

   if (atm_present) then
      do eai = 1,num_inst_atm
         if (iamin_ATMID(eai)) then
            if (drv_threading) call seq_comm_setnthreads(nthreads_ATMID)
            k1 = mct_aVect_indexRa(cdata_aa(eai)%dom%data,"area"  ,perrWith='aa area ')
            k2 = mct_aVect_indexRa(cdata_aa(eai)%dom%data,"aream" ,perrWith='aa aream')
            cdata_aa(eai)%dom%data%rAttr(k2,:) = cdata_aa(eai)%dom%data%rAttr(k1,:)
            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
         endif
      enddo
   endif

   if (lnd_present) then
      do eli = 1,num_inst_lnd
         if (iamin_LNDID(eli)) then
            if (drv_threading) call seq_comm_setnthreads(nthreads_LNDID)
            k1 = mct_aVect_indexRa(cdata_ll(eli)%dom%data,"area"  ,perrWith='ll area ')
            k2 = mct_aVect_indexRa(cdata_ll(eli)%dom%data,"aream" ,perrWith='ll aream')
            cdata_ll(eli)%dom%data%rAttr(k2,:) = cdata_ll(eli)%dom%data%rAttr(k1,:)
            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
         endif
      enddo
   end if

   if (sno_present) then
      do eli = 1,num_inst_lnd
         if (iamin_LNDID(eli)) then
            if (drv_threading) call seq_comm_setnthreads(nthreads_LNDID)
            k1 = mct_aVect_indexRa(cdata_ss(eli)%dom%data,"area"  ,perrWith='ss area ')
            k2 = mct_aVect_indexRa(cdata_ss(eli)%dom%data,"aream" ,perrWith='ss aream')
            cdata_ss(eli)%dom%data%rAttr(k2,:) = cdata_ss(eli)%dom%data%rAttr(k1,:)
            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
         endif
      enddo
   end if

   if (rof_present) then
      do eri = 1,num_inst_rof
         if (iamin_ROFID(eri)) then
            if (drv_threading) call seq_comm_setnthreads(nthreads_ROFID)
            k1 = mct_aVect_indexRa(cdata_rr(eri)%dom%data,"area"  ,perrWith='rr area ')
            k2 = mct_aVect_indexRa(cdata_rr(eri)%dom%data,"aream" ,perrWith='rr aream')
            cdata_rr(eri)%dom%data%rAttr(k2,:) = cdata_rr(eri)%dom%data%rAttr(k1,:)
            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
         endif
      enddo
   end if

   if (ocn_present) then
      do eoi = 1,num_inst_ocn
         if (iamin_OCNID(eoi)) then
            if (drv_threading) call seq_comm_setnthreads(nthreads_OCNID)
            k1 = mct_aVect_indexRa(cdata_oo(eoi)%dom%data,"area"  ,perrWith='oo area ')
            k2 = mct_aVect_indexRa(cdata_oo(eoi)%dom%data,"aream" ,perrWith='oo aream')
            cdata_oo(eoi)%dom%data%rAttr(k2,:) = cdata_oo(eoi)%dom%data%rAttr(k1,:)
            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
         endif
      enddo
   endif

   if (ice_present) then
      do eii = 1,num_inst_ice
         if (iamin_ICEID(eii)) then
            if (drv_threading) call seq_comm_setnthreads(nthreads_ICEID)
            k1 = mct_aVect_indexRa(cdata_ii(eii)%dom%data,"area"  ,perrWith='ii area ')
            k2 = mct_aVect_indexRa(cdata_ii(eii)%dom%data,"aream" ,perrWith='ii aream')
            cdata_ii(eii)%dom%data%rAttr(k2,:) = cdata_ii(eii)%dom%data%rAttr(k1,:)
            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
         endif
      enddo
   endif

   if (glc_present) then
      do egi = 1,num_inst_glc
         if (iamin_GLCID(egi)) then
            if (drv_threading) call seq_comm_setnthreads(nthreads_GLCID)
            k1 = mct_aVect_indexRa(cdata_gg(egi)%dom%data,"area"  ,perrWith='gg area ')
            k2 = mct_aVect_indexRa(cdata_gg(egi)%dom%data,"aream" ,perrWith='gg aream')
            cdata_gg(egi)%dom%data%rAttr(k2,:) = cdata_gg(egi)%dom%data%rAttr(k1,:)
            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
         endif
      enddo
   endif

   if (wav_present) then
      do ewi = 1,num_inst_wav
         if (iamin_WAVID(ewi)) then
            if (drv_threading) call seq_comm_setnthreads(nthreads_WAVID)
            k1 = mct_aVect_indexRa(cdata_ww(ewi)%dom%data,"area"  ,perrWith='ww area ')
            k2 = mct_aVect_indexRa(cdata_ww(ewi)%dom%data,"aream" ,perrWith='ww aream')
            cdata_ww(ewi)%dom%data%rAttr(k2,:) = cdata_ww(ewi)%dom%data%rAttr(k1,:)
            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
         endif
      enddo
   endif

   !-----------------------------------------------------------------------------
   ! Initialize driver rearrangers and AVs on driver
   ! Initialize cdata_*x data
   ! Zero out x2*_** in case it never gets used then it'll produce zeros in diags
   ! For ensembles, create only a single dom_*x for the coupler based on the
   !   first ensemble member.  otherwise, just extend the dom_** and dom_*x to
   !   other ensemble members.
   !-----------------------------------------------------------------------------
   call t_startf('driver_init_xxx2xxx')

   call mpi_barrier(mpicom_GLOID,ierr)

   if (atm_present) then
      do eai = 1,num_inst_atm
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) call shr_sys_flush(logunit)
         if (iamin_CPLATMID(eai)) then
            if (eai == 1) then
               if (iamroot_CPLID) write(logunit,F0I) 'creating gsmap_ax'
               if (iamroot_CPLID) call shr_sys_flush(logunit)
               call seq_mctext_gsmapInit(gsmap_aa(eai), ATMID(eai), gsmap_ax, CPLID, CPLATMID(eai))
            endif  
            if (iamroot_CPLID) write(logunit,F0I) 'Initializing mapper_Ca2x',eai
            if (iamroot_CPLID) call shr_sys_flush(logunit)
            call seq_map_init_rearrsplit(mapper_Ca2x(eai), gsmap_aa(eai), ATMID(eai), gsmap_ax, CPLID   , CPLATMID(eai))
            if (iamroot_CPLID) write(logunit,F0I) 'Initializing mapper_Cx2a',eai
            if (iamroot_CPLID) call shr_sys_flush(logunit)
            call seq_map_init_rearrsplit(mapper_Cx2a(eai), gsmap_ax, CPLID   , gsmap_aa(eai), ATMID(eai), CPLATMID(eai))
            call seq_mctext_avInit(x2a_aa(eai), ATMID(eai), x2a_ax(eai), CPLID, gsmap_ax, CPLATMID(eai))
            call seq_mctext_avInit(a2x_aa(eai), ATMID(eai), a2x_ax(eai), CPLID, gsmap_ax, CPLATMID(eai))
            if (eai == 1) then  ! create dom_ax
               if (iamroot_CPLID) write(logunit,F0I) 'creating dom_ax'
               if (iamroot_CPLID) call shr_sys_flush(logunit)
               call seq_mctext_gGridInit(dom_aa(eai), ATMID(eai), dom_ax, CPLID, gsmap_ax, CPLATMID(eai))
               call seq_map_map(mapper_Ca2x(eai), dom_aa(eai)%data, dom_ax%data, msgtag=CPLATMID(eai)*100+eai*10+1)
            else                ! veryify other ensembles have same domain by comparing to dom_ax
               if (iamroot_CPLID) write(logunit,F0I) 'comparing atm domain ensemble number ',eai
               if (iamroot_CPLID) call shr_sys_flush(logunit)
               call seq_mctext_avExtend(dom_ax%data,      CPLID,      CPLATMID(eai))
               call seq_mctext_gGridInit(dom_aa(eai), ATMID(eai), dom_tmp, CPLID, gsmap_ax, CPLATMID(eai))
               call seq_map_map(mapper_Ca2x(eai), dom_aa(eai)%data, dom_tmp%data, msgtag=CPLATMID(eai)*100+eai*10+1)
               if (iamin_CPLID) call seq_domain_compare(dom_ax,dom_tmp,mpicom_CPLID)
               call mct_ggrid_clean(dom_tmp,rc)
            endif
            call mct_avect_zero(x2a_aa(eai))
            call seq_map_map(mapper_Ca2x(eai), x2a_aa(eai), x2a_ax(eai), msgtag=CPLATMID(eai)*100+eai*10+3)
        endif
     enddo
   endif

   call mpi_barrier(mpicom_GLOID,ierr)

   if (lnd_present) then
      do eli = 1,num_inst_lnd
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) call shr_sys_flush(logunit)
         if (iamin_CPLLNDID(eli)) then
            if (eli == 1) then
               if (iamroot_CPLID) write(logunit,F0I) 'creating gsmap_lx'
               if (iamroot_CPLID) call shr_sys_flush(logunit)
               call seq_mctext_gsmapInit(gsmap_ll(eli), LNDID(eli), gsmap_lx, CPLID, CPLLNDID(eli))
            endif
            if (iamroot_CPLID) write(logunit,F0I) 'Initializing mapper_Cl2x',eli
            if (iamroot_CPLID) call shr_sys_flush(logunit)
            call seq_map_init_rearrsplit(mapper_Cl2x(eli), gsmap_ll(eli), LNDID(eli), gsmap_lx, CPLID, CPLLNDID(eli))
            if (iamroot_CPLID) write(logunit,F0I) 'Initializing mapper_Cx2l',eli
            if (iamroot_CPLID) call shr_sys_flush(logunit)
            call seq_map_init_rearrsplit(mapper_Cx2l(eli), gsmap_lx, CPLID, gsmap_ll(eli), LNDID(eli), CPLLNDID(eli))
            call seq_mctext_avInit(x2l_ll(eli), LNDID(eli), x2l_lx(eli), CPLID, gsmap_lx, CPLLNDID(eli))
            call seq_mctext_avInit(l2x_ll(eli), LNDID(eli), l2x_lx(eli), CPLID, gsmap_lx, CPLLNDID(eli))
            if (eli == 1) then
               if (iamroot_CPLID) write(logunit,F0I) 'creating dom_lx'
               if (iamroot_CPLID) call shr_sys_flush(logunit)
               call seq_mctext_gGridInit(dom_ll(eli), LNDID(eli), dom_lx, CPLID, gsmap_lx, CPLLNDID(eli))
               call seq_map_map(mapper_Cl2x(eli), dom_ll(eli)%data, dom_lx%data, msgtag=CPLLNDID(eli)*100+eli*10+1)
            else
               if (iamroot_CPLID) write(logunit,F0I) 'comparing lnd domain ensemble number ',eli
               if (iamroot_CPLID) call shr_sys_flush(logunit)
               call seq_mctext_avExtend(dom_lx%data,      CPLID,      CPLLNDID(eli))
               call seq_mctext_gGridInit(dom_ll(eli), LNDID(eli), dom_tmp, CPLID, gsmap_lx, CPLLNDID(eli))
               call seq_map_map(mapper_Cl2x(eli), dom_ll(eli)%data, dom_tmp%data, msgtag=CPLLNDID(eli)*100+eli*10+1)
               if (iamin_CPLID) call seq_domain_compare(dom_lx,dom_tmp,mpicom_CPLID)
               call mct_ggrid_clean(dom_tmp,rc)
            endif
            call mct_avect_zero(x2l_ll(eli))
            call seq_map_map(mapper_Cl2x(eli), x2l_ll(eli), x2l_lx(eli), msgtag=CPLLNDID(eli)*100+eli*10+3)
         endif
      enddo
   end if

   call mpi_barrier(mpicom_GLOID,ierr)

   if (sno_present) then
      do eli = 1,num_inst_lnd
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) call shr_sys_flush(logunit)
         if (iamin_CPLLNDID(eli)) then
            if (eli == 1) then
               if (iamroot_CPLID) write(logunit,F0I) 'creating gsmap_sx'
               if (iamroot_CPLID) call shr_sys_flush(logunit)
               call seq_mctext_gsmapInit(gsmap_ss(eli), LNDID(eli), gsmap_sx, CPLID, CPLLNDID(eli))
            endif
            if (iamroot_CPLID) write(logunit,F0I) 'Initializing mapper_Cs2x',eli
            if (iamroot_CPLID) call shr_sys_flush(logunit)
            call seq_map_init_rearrsplit(mapper_Cs2x(eli), gsmap_ss(eli), LNDID(eli), gsmap_sx, CPLID    , CPLLNDID(eli))
            if (iamroot_CPLID) write(logunit,F0I) 'Initializing mapper_Cx2s',eli
            if (iamroot_CPLID) call shr_sys_flush(logunit)
            call seq_map_init_rearrsplit(mapper_Cx2s(eli), gsmap_sx, CPLID    , gsmap_ss(eli), LNDID(eli), CPLLNDID(eli))
            call seq_mctext_avInit(x2s_ss(eli), LNDID(eli), x2s_sx(eli), CPLID, gsmap_sx, CPLLNDID(eli))
            call seq_mctext_avInit(s2x_ss(eli), LNDID(eli), s2x_sx(eli), CPLID, gsmap_sx, CPLLNDID(eli))
            if (eli == 1) then
               if (iamroot_CPLID) write(logunit,F0I) 'creating dom_sx'
               if (iamroot_CPLID) call shr_sys_flush(logunit)
               call seq_mctext_gGridInit(dom_ss(eli), LNDID(eli), dom_sx, CPLID, gsmap_sx, CPLLNDID(eli))
               call seq_map_map(mapper_Cs2x(eli), dom_ss(eli)%data, dom_sx%data, msgtag=CPLLNDID(eli)*100+eli*10+1001)
            else
               if (iamroot_CPLID) write(logunit,F0I) 'comparing sno domain ensemble number ',eli
               if (iamroot_CPLID) call shr_sys_flush(logunit)
               call seq_mctext_avExtend(dom_sx%data     , CPLID     , CPLLNDID(eli))
               call seq_mctext_gGridInit(dom_ss(eli), LNDID(eli), dom_tmp, CPLID, gsmap_sx, CPLLNDID(eli))
               call seq_map_map(mapper_Cs2x(eli), dom_ss(eli)%data, dom_tmp%data, msgtag=CPLLNDID(eli)*100+eli*10+1001)
               if (iamin_CPLID) call seq_domain_compare(dom_sx,dom_tmp,mpicom_CPLID)
               call mct_ggrid_clean(dom_tmp,rc)
            endif
            call mct_avect_zero(x2s_ss(eli))
            call seq_map_map(mapper_Cs2x(eli), x2s_ss(eli), x2s_sx(eli), msgtag=CPLLNDID(eli)*100+eli*10+1003)
         endif
      enddo
   end if

   call mpi_barrier(mpicom_GLOID,ierr)

   if (rof_present) then
      do eri = 1,num_inst_rof
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) call shr_sys_flush(logunit)
         if (iamin_CPLROFID(eri)) then
            if (eri == 1) then
               if (iamroot_CPLID) write(logunit,F0I) 'creating gsmap_rx'
               if (iamroot_CPLID) call shr_sys_flush(logunit)
               call seq_mctext_gsmapInit(gsmap_rr(eri), ROFID(eri), gsmap_rx, CPLID, CPLROFID(eri))
            endif
            if (iamroot_CPLID) write(logunit,F0I) 'Initializing mapper_Cr2x',eri
            if (iamroot_CPLID) call shr_sys_flush(logunit)
            call seq_map_init_rearrsplit(mapper_Cr2x(eri), gsmap_rr(eri), ROFID(eri), gsmap_rx, CPLID, CPLROFID(eri))
            if (iamroot_CPLID) write(logunit,F0I) 'Initializing mapper_Cx2r',eri
            if (iamroot_CPLID) call shr_sys_flush(logunit)
            call seq_map_init_rearrsplit(mapper_Cx2r(eri), gsmap_rx, CPLID, gsmap_rr(eri), ROFID(eri), CPLROFID(eri))
            call seq_mctext_avInit(x2r_rr(eri), ROFID(eri), x2r_rx(eri), CPLID, gsmap_rx, CPLROFID(eri))
            call seq_mctext_avInit(r2x_rr(eri), ROFID(eri), r2x_rx(eri), CPLID, gsmap_rx, CPLROFID(eri))
            if (eri == 1) then
               if (iamroot_CPLID) write(logunit,F0I) 'creating dom_rx'
               if (iamroot_CPLID) call shr_sys_flush(logunit)
               call seq_mctext_gGridInit(dom_rr(eri), ROFID(eri), dom_rx, CPLID, gsmap_rx, CPLROFID(eri))
               call seq_map_map(mapper_Cr2x(eri), dom_rr(eri)%data, dom_rx%data, msgtag=CPLROFID(eri)*100+eri*10+1)
            else
               if (iamroot_CPLID) write(logunit,F0I) 'comparing rof domain ensemble number ',eri
               if (iamroot_CPLID) call shr_sys_flush(logunit)
               call seq_mctext_avExtend(dom_rx%data,      CPLID,      CPLROFID(eri))
               call seq_mctext_gGridInit(dom_rr(eri), ROFID(eri), dom_tmp, CPLID, gsmap_rx, CPLROFID(eri))
               call seq_map_map(mapper_Cr2x(eri), dom_rr(eri)%data, dom_tmp%data, msgtag=CPLROFID(eri)*100+eri*10+1)
               if (iamin_CPLID) call seq_domain_compare(dom_rx,dom_tmp,mpicom_CPLID)
               call mct_ggrid_clean(dom_tmp,rc)
            endif
            call mct_avect_zero(x2r_rr(eri))
            call seq_map_map(mapper_Cr2x(eri), x2r_rr(eri), x2r_rx(eri), msgtag=CPLROFID(eri)*100+eri*10+3)
         endif
      enddo
   end if

   call mpi_barrier(mpicom_GLOID,ierr)

   if (ice_present) then
      do eii = 1,num_inst_ice
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) call shr_sys_flush(logunit)
         if (iamin_CPLICEID(eii)) then
            if (eii == 1) then
               if (iamroot_CPLID) write(logunit,F0I) 'creating gsmap_ix'
               if (iamroot_CPLID) call shr_sys_flush(logunit)
               call seq_mctext_gsmapInit(gsmap_ii(eii), ICEID(eii), gsmap_ix, CPLID, CPLICEID(eii))
            endif
            if (iamroot_CPLID) write(logunit,F0I) 'Initializing mapper_Ci2x',eii
            if (iamroot_CPLID) call shr_sys_flush(logunit)
            call seq_map_init_rearrsplit(mapper_Ci2x(eii), gsmap_ii(eii), ICEID(eii), gsmap_ix, CPLID   , CPLICEID(eii))
            if (iamroot_CPLID) write(logunit,F0I) 'Initializing mapper_Cx2i',eii
            if (iamroot_CPLID) call shr_sys_flush(logunit)
            call seq_map_init_rearrsplit(mapper_Cx2i(eii), gsmap_ix, CPLID   , gsmap_ii(eii), ICEID(eii), CPLICEID(eii))
            call seq_mctext_avInit(x2i_ii(eii), ICEID(eii), x2i_ix(eii), CPLID, gsmap_ix, CPLICEID(eii))
            call seq_mctext_avInit(i2x_ii(eii), ICEID(eii), i2x_ix(eii), CPLID, gsmap_ix, CPLICEID(eii))
            if (eii == 1) then
               if (iamroot_CPLID) write(logunit,F0I) 'creating dom_ix'
               if (iamroot_CPLID) call shr_sys_flush(logunit)
               call seq_mctext_gGridInit(dom_ii(eii), ICEID(eii), dom_ix, CPLID, gsmap_ix, CPLICEID(eii))
               call seq_map_map(mapper_Ci2x(eii), dom_ii(eii)%data, dom_ix%data, msgtag=CPLICEID(eii)*100+eii*10+1)
            else
               if (iamroot_CPLID) write(logunit,F0I) 'comparing ice domain ensemble number ',eii
               if (iamroot_CPLID) call shr_sys_flush(logunit)
               call seq_mctext_avExtend(dom_ix%data     , CPLID     , CPLICEID(eii))
               call seq_mctext_gGridInit(dom_ii(eii), ICEID(eii), dom_tmp, CPLID, gsmap_ix, CPLICEID(eii))
               call seq_map_map(mapper_Ci2x(eii), dom_ii(eii)%data, dom_tmp%data, msgtag=CPLICEID(eii)*100+eii*10+1)
               if (iamin_CPLID) call seq_domain_compare(dom_ix,dom_tmp,mpicom_CPLID)
               call mct_ggrid_clean(dom_tmp,rc)
            endif
            call mct_avect_zero(x2i_ii(eii))
            call seq_map_map(mapper_Ci2x(eii), x2i_ii(eii), x2i_ix(eii), msgtag=CPLICEID(eii)*100+eii*10+3)
         endif
      enddo
   endif

   call mpi_barrier(mpicom_GLOID,ierr)

   if (glc_present) then
      do egi = 1,num_inst_glc
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) call shr_sys_flush(logunit)
         if (iamin_CPLGLCID(egi)) then
            if (egi == 1) then
               if (iamroot_CPLID) write(logunit,F0I) 'creating gsmap_gx'
               if (iamroot_CPLID) call shr_sys_flush(logunit)
               call seq_mctext_gsmapInit(gsmap_gg(egi), GLCID(egi), gsmap_gx, CPLID, CPLGLCID(egi))
            endif
            if (iamroot_CPLID) write(logunit,F0I) 'Initializing mapper_Cg2x',egi
            if (iamroot_CPLID) call shr_sys_flush(logunit)
            call seq_map_init_rearrsplit(mapper_Cg2x(egi), gsmap_gg(egi), GLCID(egi), gsmap_gx, CPLID   , CPLGLCID(egi))
            if (iamroot_CPLID) write(logunit,F0I) 'Initializing mapper_Cx2g',egi
            if (iamroot_CPLID) call shr_sys_flush(logunit)
            call seq_map_init_rearrsplit(mapper_Cx2g(egi), gsmap_gx, CPLID   , gsmap_gg(egi), GLCID(egi), CPLGLCID(egi))
            call seq_mctext_avInit(x2g_gg(egi), GLCID(egi), x2g_gx(egi), CPLID, gsmap_gx, CPLGLCID(egi))
            call seq_mctext_avInit(g2x_gg(egi), GLCID(egi), g2x_gx(egi), CPLID, gsmap_gx, CPLGLCID(egi))
            if (egi == 1) then
               if (iamroot_CPLID) write(logunit,F0I) 'creating dom_gx'
               if (iamroot_CPLID) call shr_sys_flush(logunit)
               call seq_mctext_gGridInit(dom_gg(egi), GLCID(egi), dom_gx, CPLID, gsmap_gx, CPLGLCID(egi))
               call seq_map_map(mapper_Cg2x(egi), dom_gg(egi)%data, dom_gx%data, msgtag=CPLGLCID(egi)*100+egi*10+1)
            else
               if (iamroot_CPLID) write(logunit,F0I) 'comparing glc domain ensemble number ',egi
               if (iamroot_CPLID) call shr_sys_flush(logunit)
               call seq_mctext_avExtend(dom_gx%data     , CPLID     , CPLGLCID(egi))
               call seq_mctext_gGridInit(dom_gg(egi), GLCID(egi), dom_tmp, CPLID, gsmap_gx, CPLGLCID(egi))
               call seq_map_map(mapper_Cg2x(egi), dom_gg(egi)%data, dom_tmp%data, msgtag=CPLGLCID(egi)*100+egi*10+1)
               if (iamin_CPLID) call seq_domain_compare(dom_gx,dom_tmp,mpicom_CPLID)
               call mct_ggrid_clean(dom_tmp,rc)
            endif
            call mct_avect_zero(x2g_gg(egi))
            call seq_map_map(mapper_Cg2x(egi), x2g_gg(egi), x2g_gx(egi), msgtag=CPLGLCID(egi)*100+egi*10+3)
         endif
      enddo
   endif

   call mpi_barrier(mpicom_GLOID,ierr)

   if (wav_present) then
      do ewi = 1,num_inst_wav
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) call shr_sys_flush(logunit)
         if (iamin_CPLWAVID(ewi)) then
            if (ewi == 1) then
               if (iamroot_CPLID) write(logunit,F0I) 'creating gsmap_wx'
               if (iamroot_CPLID) call shr_sys_flush(logunit)
               call seq_mctext_gsmapInit(gsmap_ww(ewi), WAVID(ewi), gsmap_wx, CPLID, CPLWAVID(ewi))
            endif
            if (iamroot_CPLID) write(logunit,F0I) 'Initializing mapper_Cw2x',ewi
            if (iamroot_CPLID) call shr_sys_flush(logunit)
            call seq_map_init_rearrsplit(mapper_Cw2x(ewi), gsmap_ww(ewi), WAVID(ewi), gsmap_wx, CPLID   , CPLWAVID(ewi))
            if (iamroot_CPLID) write(logunit,F0I) 'Initializing mapper_Cx2w',ewi
            if (iamroot_CPLID) call shr_sys_flush(logunit)
            call seq_map_init_rearrsplit(mapper_Cx2w(ewi), gsmap_wx, CPLID   , gsmap_ww(ewi), WAVID(ewi), CPLWAVID(ewi))
            call seq_mctext_avInit(x2w_ww(ewi), WAVID(ewi), x2w_wx(ewi), CPLID, gsmap_wx, CPLWAVID(ewi))
            call seq_mctext_avInit(w2x_ww(ewi), WAVID(ewi), w2x_wx(ewi), CPLID, gsmap_wx, CPLWAVID(ewi))
            if (ewi == 1) then
               if (iamroot_CPLID) write(logunit,F0I) 'creating dom_wx'
               if (iamroot_CPLID) call shr_sys_flush(logunit)
               call seq_mctext_gGridInit(dom_ww(ewi), WAVID(ewi), dom_wx, CPLID, gsmap_wx, CPLWAVID(ewi))
               call seq_map_map(mapper_Cw2x(ewi), dom_ww(ewi)%data, dom_wx%data, msgtag=CPLWAVID(ewi)*100+ewi*10+1)
            else
               if (iamroot_CPLID) write(logunit,F0I) 'comparing wav domain ensemble number ',ewi
               if (iamroot_CPLID) call shr_sys_flush(logunit)
               call seq_mctext_avExtend(dom_wx%data     , CPLID     , CPLWAVID(ewi))
               call seq_mctext_gGridInit(dom_ww(ewi), WAVID(ewi), dom_tmp, CPLID, gsmap_wx, CPLWAVID(ewi))
               call seq_map_map(mapper_Cw2x(ewi), dom_ww(ewi)%data, dom_tmp%data, msgtag=CPLWAVID(ewi)*100+ewi*10+1)
               if (iamin_CPLID) call seq_domain_compare(dom_wx,dom_tmp,mpicom_CPLID)
               call mct_ggrid_clean(dom_tmp,rc)
            endif
            call mct_avect_zero(x2w_ww(ewi))
            call seq_map_map(mapper_Cw2x(ewi), x2w_ww(ewi), x2w_wx(ewi), msgtag=CPLWAVID(ewi)*100+ewi*10+3)
         endif
      enddo
   endif

   call mpi_barrier(mpicom_GLOID,ierr)

   if (ocn_present) then
      do eoi = 1,num_inst_ocn
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) call shr_sys_flush(logunit)
         if (iamin_CPLOCNID(eoi)) then
            if (eoi == 1) then
               if (iamroot_CPLID) write(logunit,F0I) 'creating gsmap_ox'
               if (iamroot_CPLID) call shr_sys_flush(logunit)
               call seq_mctext_gsmapInit(gsmap_oo(eoi), OCNID(eoi), gsmap_ox, CPLID, CPLOCNID(eoi))
            endif
            if (iamroot_CPLID) write(logunit,F0I) 'Initializing mapper_Co2x',eoi
            if (iamroot_CPLID) call shr_sys_flush(logunit)
            call seq_map_init_rearrsplit(mapper_Co2x(eoi), gsmap_oo(eoi), OCNID(eoi), gsmap_ox, CPLID   , CPLOCNID(eoi))
            if (iamroot_CPLID) write(logunit,F0I) 'Initializing mapper_Cx2o',eoi
            if (iamroot_CPLID) call shr_sys_flush(logunit)
            call seq_map_init_rearrsplit(mapper_Cx2o(eoi), gsmap_ox, CPLID   , gsmap_oo(eoi), OCNID(eoi), CPLOCNID(eoi))
            call seq_mctext_avInit(x2o_oo(eoi), OCNID(eoi), x2o_ox(eoi), CPLID, gsmap_ox, CPLOCNID(eoi))
            call seq_mctext_avInit(o2x_oo(eoi), OCNID(eoi), o2x_ox(eoi), CPLID, gsmap_ox, CPLOCNID(eoi))
            if (eoi == 1) then
               if (iamroot_CPLID) write(logunit,F0I) 'creating dom_ox'
               if (iamroot_CPLID) call shr_sys_flush(logunit)
               call seq_mctext_gGridInit(dom_oo(eoi), OCNID(eoi), dom_ox, CPLID, gsmap_ox, CPLOCNID(eoi))
               call seq_map_map(mapper_Co2x(eoi), dom_oo(eoi)%data, dom_ox%data, msgtag=CPLOCNID(eoi)*100+eoi*10+1)
            else
               if (iamroot_CPLID) write(logunit,F0I) 'comparing ocn domain ensemble number ',eoi
               if (iamroot_CPLID) call shr_sys_flush(logunit)
               call seq_mctext_avExtend(dom_ox%data     , CPLID     , CPLOCNID(eoi))
               call seq_mctext_gGridInit(dom_oo(eoi), OCNID(eoi), dom_tmp, CPLID, gsmap_ox, CPLOCNID(eoi))
               call seq_map_map(mapper_Co2x(eoi), dom_oo(eoi)%data, dom_tmp%data, msgtag=CPLOCNID(eoi)*100+eoi*10+1)
               if (iamin_CPLID) call seq_domain_compare(dom_ox,dom_tmp,mpicom_CPLID)
               call mct_ggrid_clean(dom_tmp,rc)
            endif
            call mct_avect_zero(x2o_oo(eoi))
            call seq_map_map(mapper_Co2x(eoi), x2o_oo(eoi), x2o_ox(eoi), msgtag=CPLOCNID(eoi)*100+eoi*10+3)

            !--- this needs to be here because it's used on cplocn pes for mapping
            call mct_avect_init(x2oacc_ox(eoi), x2o_ox(eoi), mct_aVect_lsize(x2o_ox(eoi)))
            call mct_aVect_zero(x2oacc_ox(eoi))
            x2oacc_ox_cnt = 0
         endif
      enddo
   endif

   call t_stopf  ('driver_init_xxx2xxx')

   !-----------------------------------------------------------------------------
   ! Initialize Remaining Coupler AVects
   !-----------------------------------------------------------------------------

   if (iamin_CPLID) then
      lsize_a = mct_aVect_lsize(a2x_ax(ens1))
      lsize_l = mct_aVect_lsize(l2x_lx(ens1))
      lsize_s = mct_aVect_lsize(s2x_sx(ens1))
      lsize_r = mct_aVect_lsize(r2x_rx(ens1))
      lsize_o = mct_aVect_lsize(o2x_ox(ens1))
      lsize_i = mct_aVect_lsize(i2x_ix(ens1))
      lsize_g = mct_aVect_lsize(g2x_gx(ens1))
      lsize_w = mct_aVect_lsize(w2x_wx(ens1))

      do eai = 1,num_inst_atm
         call mct_aVect_init(a2x_ix(eai), rList=seq_flds_a2x_fields, lsize=lsize_i)
         call mct_aVect_init(a2x_lx(eai), rList=seq_flds_a2x_fields, lsize=lsize_l)
         call mct_aVect_init(a2x_ox(eai), rList=seq_flds_a2x_fields, lsize=lsize_o)
         call mct_aVect_init(a2x_wx(eai), rList=seq_flds_a2x_fields, lsize=lsize_w)
         call mct_aVect_zero(a2x_ix(eai))
         call mct_aVect_zero(a2x_lx(eai))
         call mct_aVect_zero(a2x_ox(eai))
         call mct_aVect_zero(a2x_wx(eai))
      enddo

      do eli = 1,num_inst_lnd
         call mct_aVect_init(l2x_ax(eli), rList=seq_flds_l2x_fields, lsize=lsize_a)
         call mct_aVect_init(s2x_gx(eli), rList=seq_flds_s2x_fields, lsize=lsize_g)
         call mct_aVect_zero(l2x_ax(eli))
         call mct_aVect_zero(s2x_gx(eli))
         call mct_avect_init(x2racc_lx(eli), rlist=seq_flds_x2r_fields, lsize=lsize_l)
         call mct_aVect_zero(x2racc_lx(eli))
      enddo
      x2racc_lx_cnt = 0
      call mct_avect_init(x2r_rx_tmp, rList=seq_flds_x2r_fields, lsize=lsize_r)
      call mct_avect_zero(x2r_rx_tmp)

      do eri = 1,num_inst_rof
         call mct_aVect_init(r2x_ox(eri), rList=seq_flds_r2x_fields, lsize=lsize_o)
         call mct_aVect_zero(r2x_ox(eri))
         call mct_aVect_init(r2x_lx(eri), rlist=seq_flds_r2x_fields, lsize=lsize_l)
         call mct_aVect_zero(r2x_lx(eri)) 
         call mct_avect_init(r2xacc_rx(eri), rList=seq_flds_r2x_fields, lsize=lsize_r)
         call mct_aVect_zero(r2xacc_rx(eri))
      enddo
      r2xacc_rx_cnt = 0

      do eoi = 1,num_inst_ocn
         call mct_aVect_init(o2x_ax(eoi), rList=seq_flds_o2x_fields, lsize=lsize_a)
         call mct_aVect_init(o2x_ix(eoi), rList=seq_flds_o2x_fields, lsize=lsize_i)
         call mct_aVect_init(o2x_wx(eoi), rList=seq_flds_o2x_fields, lsize=lsize_w)
         call mct_aVect_zero(o2x_ax(eoi))
         call mct_aVect_zero(o2x_ix(eoi))
         call mct_aVect_zero(o2x_wx(eoi))
      enddo

      do eii = 1,num_inst_ice
         call mct_aVect_init(i2x_ax(eii), rList=seq_flds_i2x_fields, lsize=lsize_a)
         call mct_aVect_init(i2x_ox(eii), rList=seq_flds_i2x_fields, lsize=lsize_o)
         call mct_aVect_init(i2x_wx(eii), rList=seq_flds_i2x_fields, lsize=lsize_w)
         call mct_aVect_zero(i2x_ax(eii))
         call mct_aVect_zero(i2x_ox(eii))
         call mct_aVect_zero(i2x_wx(eii))
      enddo

      do egi = 1,num_inst_glc
         call mct_aVect_init(g2x_sx(egi), rList=seq_flds_g2x_fields, lsize=lsize_s)
         call mct_aVect_zero(g2x_sx(egi))
      enddo

      do ewi = 1,num_inst_wav
         call mct_aVect_init(w2x_ox(ewi), rList=seq_flds_w2x_fields, lsize=lsize_o)
         call mct_aVect_zero(w2x_ox(ewi))
      enddo

      allocate(xao_ax(num_inst_xao))
      allocate(xao_ox(num_inst_xao))
      do exi = 1,num_inst_xao
         call mct_aVect_init(xao_ax(exi), rList=seq_flds_xao_fields, lsize=lsize_a)
         call mct_aVect_zero(xao_ax(exi))
         call mct_aVect_init(xao_ox(exi), rList=seq_flds_xao_fields, lsize=lsize_o)
         call mct_aVect_zero(xao_ox(exi))
      enddo

      allocate(fractions_ax(num_inst_frc))
      allocate(fractions_lx(num_inst_frc))
      allocate(fractions_ox(num_inst_frc))
      allocate(fractions_ix(num_inst_frc))
      allocate(fractions_gx(num_inst_frc))
      allocate(fractions_rx(num_inst_frc))
      allocate(fractions_wx(num_inst_frc))

      if (rof_present) then
         index_r2x_Forr_roff  = mct_avect_indexRA(r2x_rx(ens1)   ,"Forr_roff" ,perrWith='indexa Forr_roff')
         index_r2x_Forr_ioff  = mct_avect_indexRA(r2x_rx(ens1)   ,"Forr_ioff" ,perrWith='indexa Forr_ioff')
         index_r2x_Flrr_flood = mct_avect_indexRA(r2x_rx(ens1)   ,"Flrr_flood",perrWith='indexa Flrr_flood')
      endif
      if (lnd_present) then
         index_l2x_Flrl_rofliq = mct_avect_indexRA(l2x_lx(ens1)   ,"Flrl_rofliq",perrWith='indexb Flrl_rofliq')
         index_l2x_Flrl_rofice = mct_avect_indexRA(l2x_lx(ens1)   ,"Flrl_rofice",perrWith='indexb Flrl_rofice')
         index_x2r_Flrl_rofliq = mct_avect_indexRA(x2racc_lx(ens1),"Flrl_rofliq",perrWith='indexa Flrl_rofliq')
         index_x2r_Flrl_rofice = mct_avect_indexRA(x2racc_lx(ens1),"Flrl_rofice",perrWith='indexa Flrl_rofice')
      endif
   endif

   !-----------------------------------------------------------------------------
   ! Remainder of initialization
   !-----------------------------------------------------------------------------

   if (iamin_CPLID) then
      if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)

      !-----------------------------------------------------------------------------
      ! Initialize mapping
      ! Read aream into domains!
      !-----------------------------------------------------------------------------

      call t_startf('driver_init_maps')

      if (ocn_present) then
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) write(logunit,F00) 'Initializing mapper_Sa2o'
         call seq_map_init_rcfile(mapper_Sa2o,gsmap_ax,gsmap_ox,mpicom_CPLID, &
             'seq_maps.rc','atm2ocn_smapname:','atm2ocn_smaptype:',samegrid_ao, &
             'mapper_Sa2o initialization',esmf_map_flag)
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) write(logunit,F00) 'Initializing mapper_Va2o'
         call seq_map_init_rcfile(mapper_Va2o,gsmap_ax,gsmap_ox,mpicom_CPLID, &
             'seq_maps.rc','atm2ocn_vmapname:','atm2ocn_vmaptype:',samegrid_ao, &
             'mapper_Va2o initialization',esmf_map_flag)
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) write(logunit,F00) 'Initializing mapper_Fa2o'
         call seq_map_init_rcfile(mapper_Fa2o,gsmap_ax,gsmap_ox,mpicom_CPLID, &
             'seq_maps.rc','atm2ocn_fmapname:','atm2ocn_fmaptype:',samegrid_ao, &
             'mapper_Fa2o initialization',esmf_map_flag)
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) write(logunit,F00) 'Initializing mapper_So2a'
         call seq_map_init_rcfile(mapper_So2a,gsmap_ox,gsmap_ax,mpicom_CPLID, &
             'seq_maps.rc','ocn2atm_smapname:','ocn2atm_smaptype:',samegrid_ao, &
             'mapper_So2a initialization',esmf_map_flag)
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) write(logunit,F00) 'Initializing mapper_Fo2a'
         call seq_map_init_rcfile(mapper_Fo2a,gsmap_ox,gsmap_ax,mpicom_CPLID, &
             'seq_maps.rc','ocn2atm_fmapname:','ocn2atm_fmaptype:',samegrid_ao, &
             'mapper_Fo2a initialization',esmf_map_flag)

         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) write(logunit,F00) 'Initializing mapper_Va2o vect'
         call seq_map_initvect(mapper_Va2o,vect_map,dom_ax,dom_ox, &
              gsmap_s=gsmap_ax,ni=atm_nx,nj=atm_ny,string='mapper_Va2o initvect')

         if (samegrid_ao) then
            ka = mct_aVect_indexRa(dom_ax%data, "area" )
            km = mct_aVect_indexRa(dom_ax%data, "aream" )
            dom_ax%data%rAttr(km,:) = dom_ax%data%rAttr(ka,:)
            call seq_map_map(mapper_Fa2o, dom_ax%data, dom_ox%data, fldlist='aream')
         else
            call seq_map_readdata('seq_maps.rc','ocn2atm_fmapname:',mpicom_CPLID, CPLID, &
               gsmap_s=gsmap_ox, av_s=dom_ox%data, avfld_s='aream',filefld_s='area_a', &
               gsmap_d=gsmap_ax, av_d=dom_ax%data, avfld_d='aream',filefld_d='area_b', &
               string='ocn2atm aream initialization')
         endif
      endif
      call shr_sys_flush(logunit)
      if (ice_present .and. ocn_present) then
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) write(logunit,F00) 'Initializing mapper_SFo2i'
         call seq_map_init_rearrolap(mapper_SFo2i, gsmap_ox, gsmap_ix, mpicom_CPLID, 'mapper_SFo2i')
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) write(logunit,F00) 'Initializing mapper_SFi2o'
         call seq_map_init_rearrolap(mapper_SFi2o, gsmap_ix, gsmap_ox, mpicom_CPLID, 'mapper_SFi2o')
         call seq_map_map(mapper_SFo2i, dom_ox%data, dom_ix%data, fldlist='aream')
      endif
      call shr_sys_flush(logunit)
      if (ice_present) then
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) write(logunit,F00) 'Initializing mapper_Si2a'
         call seq_map_init_rcfile(mapper_Si2a,gsmap_ix,gsmap_ax,mpicom_CPLID, &
             'seq_maps.rc','ice2atm_smapname:','ice2atm_smaptype:',samegrid_ao, &
             'mapper_Si2a initialization',esmf_map_flag)
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) write(logunit,F00) 'Initializing mapper_Fi2a'
         call seq_map_init_rcfile(mapper_Fi2a,gsmap_ix,gsmap_ax,mpicom_CPLID, &
             'seq_maps.rc','ice2atm_fmapname:','ice2atm_fmaptype:',samegrid_ao, &
             'mapper_Fi2a initialization',esmf_map_flag)
      endif
      call shr_sys_flush(logunit)
      if (rof_present .and. ocnrof_prognostic) then
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) write(logunit,F00) 'Initializing mapper_Rr2o'
         call seq_map_init_rcfile(mapper_Rr2o,gsmap_rx,gsmap_ox,mpicom_CPLID, &
             'seq_maps.rc','rof2ocn_rmapname:','rof2ocn_rmaptype:',samegrid_ro, &
             'mapper_Rr2o initialization',esmf_map_flag)
         if (.not.samegrid_ro) then
            call seq_map_readdata('seq_maps.rc','rof2ocn_rmapname:',mpicom_CPLID, CPLID, &
               gsmap_s=gsmap_rx, av_s=dom_rx%data, avfld_s='aream',filefld_s='area_a', &
               string='rof2ocn aream initialization')
         endif
         if (flood_present) then
            if (iamroot_CPLID) write(logunit,*) ' '
            if (iamroot_CPLID) write(logunit,F00) 'Initializing mapper_Fr2o'
            call seq_map_init_rcfile( mapper_Fr2o,gsmap_rx,gsmap_ox,mpicom_CPLID, &
                 'seq_maps.rc','rof2ocn_fmapname:', 'rof2ocn_fmaptype:',samegrid=.false., &
                 string='mapper_Fr2o initialization', esmf_map=esmf_map_flag)
         endif
      endif
      call shr_sys_flush(logunit)
      if (rof_present .and. lnd_present ) then ! TODO - land and rof might be on the same grid (see below)
         if (rof_prognostic) then 
            if (iamroot_CPLID) write(logunit,*) ' '
            if (iamroot_CPLID) write(logunit,F00) 'Initializing mapper_Fl2r'
            call seq_map_init_rcfile( mapper_Fl2r,gsmap_lx,gsmap_rx,mpicom_CPLID, &
                 'seq_maps.rc','lnd2rof_fmapname:','lnd2rof_fmaptype:',samegrid=.false., &
                 string='mapper_Fl2r initialization', esmf_map=esmf_map_flag)

            if (iamroot_CPLID) write(logunit,*) ' '
            if (iamroot_CPLID) write(logunit,F00) 'Initializing mapper_Sr2l'
            call seq_map_init_rcfile(mapper_Sr2l,gsmap_rx,gsmap_lx,mpicom_CPLID, &
             'seq_maps.rc','rof2lnd_smapname:','rof2lnd_smaptype:',samegrid=.false., &
             string='mapper_Sr2l initialization',esmf_map=esmf_map_flag)


         end if
         if (flood_present .and. lnd_prognostic) then
            if (iamroot_CPLID) write(logunit,*) ' '
            if (iamroot_CPLID) write(logunit,F00) 'Initializing mapper_Fr2l'
            call seq_map_init_rcfile( mapper_Fr2l,gsmap_rx,gsmap_lx,mpicom_CPLID, &
                 'seq_maps.rc','rof2lnd_fmapname:', 'rof2lnd_fmaptype:',samegrid=.false., &
                 string='mapper_Fr2l initialization', esmf_map=esmf_map_flag)
         endif
      end if
      call shr_sys_flush(logunit)
      if (lnd_present) then
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) write(logunit,F00) 'Initializing mapper_Sa2l'
         call seq_map_init_rcfile(mapper_Sa2l,gsmap_ax,gsmap_lx,mpicom_CPLID, &
             'seq_maps.rc','atm2lnd_smapname:','atm2lnd_smaptype:',samegrid_al, &
             'mapper_Sa2l initialization',esmf_map_flag)
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) write(logunit,F00) 'Initializing mapper_Fa2l'
         call seq_map_init_rcfile(mapper_Fa2l,gsmap_ax,gsmap_lx,mpicom_CPLID, &
             'seq_maps.rc','atm2lnd_fmapname:','atm2lnd_fmaptype:',samegrid_al, &
             'mapper_Fa2l initialization',esmf_map_flag)
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) write(logunit,F00) 'Initializing mapper_Sl2a'
         call seq_map_init_rcfile(mapper_Sl2a,gsmap_lx,gsmap_ax,mpicom_CPLID, &
             'seq_maps.rc','lnd2atm_smapname:','lnd2atm_smaptype:',samegrid_al, &
             'mapper_Sl2a initialization',esmf_map_flag)
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) write(logunit,F00) 'Initializing mapper_Fl2a'
         call seq_map_init_rcfile(mapper_Fl2a,gsmap_lx,gsmap_ax,mpicom_CPLID, &
             'seq_maps.rc','lnd2atm_fmapname:','lnd2atm_fmaptype:',samegrid_al, &
             'mapper_Fl2a initialization',esmf_map_flag)
         if (samegrid_al) then
            call seq_map_map(mapper_Sa2l, dom_ax%data, dom_lx%data, fldlist='aream')
         else
            call seq_map_readdata('seq_maps.rc','atm2lnd_fmapname:',mpicom_CPLID, CPLID, &
               gsmap_d=gsmap_lx, av_d=dom_lx%data, avfld_d='aream',filefld_d='area_b', &
               string='atm2lnd aream initialization')
         endif
      endif
      call shr_sys_flush(logunit)
      if (sno_present .and. glc_present) then
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) write(logunit,F00) 'Initializing mapper_SFs2g'
         call seq_map_init_rearrolap(mapper_SFs2g, gsmap_sx, gsmap_gx, mpicom_CPLID, 'mapper_SFs2g')
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) write(logunit,F00) 'Initializing mapper_SFg2s'
         call seq_map_init_rearrolap(mapper_SFg2s, gsmap_gx, gsmap_sx, mpicom_CPLID, 'mapper_SFg2s')
         call seq_map_map(mapper_SFs2g, dom_sx%data, dom_gx%data, fldlist='aream')
      endif
      call shr_sys_flush(logunit)
      if (wav_present) then
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) write(logunit,F00) 'Initializing mapper_Sa2w'
         call seq_map_init_rcfile(mapper_Sa2w,gsmap_ax,gsmap_wx,mpicom_CPLID, &
             'seq_maps.rc','atm2wav_smapname:','atm2wav_smaptype:',samegrid_aw, &
             'mapper_Sa2w initialization')
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) write(logunit,F00) 'Initializing mapper_So2w'
         call seq_map_init_rcfile(mapper_So2w,gsmap_ox,gsmap_wx,mpicom_CPLID, &
             'seq_maps.rc','ocn2wav_smapname:','ocn2wav_smaptype:',samegrid_ow, &
             'mapper_So2w initialization')
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) write(logunit,F00) 'Initializing mapper_Si2w'
         call seq_map_init_rcfile(mapper_Si2w,gsmap_ox,gsmap_wx,mpicom_CPLID, &
             'seq_maps.rc','ice2wav_smapname:','ice2wav_smaptype:',samegrid_ow, &
             'mapper_Si2w initialization')
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) write(logunit,F00) 'Initializing mapper_Sw2o'
         call seq_map_init_rcfile(mapper_Sw2o,gsmap_wx,gsmap_ox,mpicom_CPLID, &
             'seq_maps.rc','wav2ocn_smapname:','wav2ocn_smaptype:',samegrid_ow, &
             'mapper_Sw2o initialization')
      endif
      call shr_sys_flush(logunit)

      call t_stopf  ('driver_init_maps')

      !-----------------------------------------------------------------------------
      ! Check domains if appropriate
      ! This must be done after the mappers are initialized since
      ! checking is done on each processor and not with a global gather
      !-----------------------------------------------------------------------------

      if (domain_check) then
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) write(logunit,F00) 'Performing domain checking'
         call shr_sys_flush(logunit)
         call seq_domain_check( cdata_ax, cdata_ix, cdata_lx, cdata_ox, &
                                cdata_rx, cdata_gx, cdata_sx, &
                                mapper_Fi2a, mapper_SFi2o, mapper_Fo2a, &
                                mapper_SFs2g, mapper_Fa2l, mapper_Fl2a)
      endif

      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
   endif ! iamin_CPLID

   !-----------------------------------------------------------------------------
   ! Map  dom_*x to dom_** in case any domain fields have been updated on cpl pes
   ! Initialize area corrections based on aream (read in map_init) and area
   ! Area correct component initialization output fields
   ! Map initial component AVs from component to coupler pes
   !-----------------------------------------------------------------------------

   call mpi_barrier(mpicom_GLOID,ierr)

   if (atm_present) then
      do eai = 1,num_inst_atm
         if (iamin_CPLATMID(eai)) then
            call seq_map_map(mapper_Cx2a(eai), dom_ax%data, dom_aa(eai)%data, msgtag=CPLATMID(eai)*100+eai*10+5)
            if (iamin_ATMID(eai)) then
               call seq_domain_areafactinit(cdata_aa(eai),areacor_aa(eai)%mdl2drv,areacor_aa(eai)%drv2mdl,&
                    'areafact_a_'//trim(atm_name(eai)))
               call mct_avect_vecmult(a2x_aa(eai),areacor_aa(eai)%mdl2drv,seq_flds_a2x_fluxes)
            endif
            call seq_map_map(mapper_Ca2x(eai), a2x_aa(eai), a2x_ax(eai), msgtag=CPLATMID(eai)*100+eai*10+7)
         endif
      enddo
   endif

   call mpi_barrier(mpicom_GLOID,ierr)

   if (lnd_present) then
      do eli = 1,num_inst_lnd
         if (iamin_CPLLNDID(eli)) then
            call seq_map_map(mapper_Cx2l(eli), dom_lx%data, dom_ll(eli)%data, msgtag=CPLLNDID(eli)*100+eli*10+5)
            if (iamin_LNDID(eli)) then
               call seq_domain_areafactinit(cdata_ll(eli),areacor_ll(eli)%mdl2drv,areacor_ll(eli)%drv2mdl,&
                    'areafact_l_'//trim(lnd_name(eli)))
               call mct_avect_vecmult(l2x_ll(eli),areacor_ll(eli)%mdl2drv,seq_flds_l2x_fluxes)
            endif
            call seq_map_map(mapper_Cl2x(eli), l2x_ll(eli), l2x_lx(eli), msgtag=CPLLNDID(eli)*100+eli*10+7)
         endif
      enddo
   end if

   call mpi_barrier(mpicom_GLOID,ierr)

   if (sno_present) then
      do eli = 1,num_inst_lnd
         if (iamin_CPLLNDID(eli)) then
            call seq_map_map(mapper_Cx2s(eli), dom_sx%data, dom_ss(eli)%data, msgtag=CPLLNDID(eli)*100+eli*10+1005)
            if (iamin_LNDID(eli)) then
               call seq_domain_areafactinit(cdata_ss(eli),areacor_ss(eli)%mdl2drv,areacor_ss(eli)%drv2mdl,&
                    'areafact_s_'//trim(lnd_name(eli)))
               call mct_avect_vecmult(s2x_ss(eli),areacor_ss(eli)%mdl2drv,seq_flds_s2x_fluxes)
            endif
            call seq_map_map(mapper_Cs2x(eli), s2x_ss(eli), s2x_sx(eli), msgtag=CPLLNDID(eli)*100+eli*10+1007)
         endif
      enddo
   end if

   call mpi_barrier(mpicom_GLOID,ierr)

   if (rof_present) then
      do eri = 1,num_inst_rof
         if (iamin_CPLROFID(eri)) then
            call seq_map_map(mapper_Cx2r(eri), dom_rx%data, dom_rr(eri)%data, msgtag=CPLROFID(eri)*100+eri*10+5)
            if (iamin_ROFID(eri)) then
               call seq_domain_areafactinit(cdata_rr(eri),areacor_rr(eri)%mdl2drv,areacor_rr(eri)%drv2mdl,&
                    'areafact_r_'//trim(rof_name(eri)))
               call mct_avect_vecmult(r2x_rr(eri),areacor_rr(eri)%mdl2drv,seq_flds_r2x_fluxes)
            endif
            call seq_map_map(mapper_Cr2x(eri), r2x_rr(eri), r2x_rx(eri), msgtag=CPLROFID(eri)*100+eri*10+7)
         endif
      enddo
   end if

   call mpi_barrier(mpicom_GLOID,ierr)

   if (ocn_present) then
      do eoi = 1,num_inst_ocn
         if (iamin_CPLOCNID(eoi)) then
            call seq_map_map(mapper_Cx2o(eoi), dom_ox%data, dom_oo(eoi)%data, msgtag=CPLOCNID(eoi)*100+eoi*10+5)
            if (iamin_OCNID(eoi)) then
               call seq_domain_areafactinit(cdata_oo(eoi),areacor_oo(eoi)%mdl2drv,areacor_oo(eoi)%drv2mdl,&
                    'areafact_o_'//trim(ocn_name(eoi)))
               call mct_avect_vecmult(o2x_oo(eoi),areacor_oo(eoi)%mdl2drv,seq_flds_o2x_fluxes)
            endif
           call seq_map_map(mapper_Co2x(eoi), o2x_oo(eoi), o2x_ox(eoi), msgtag=CPLOCNID(eoi)*100+eoi*10+7)
        endif
      enddo
   endif

   call mpi_barrier(mpicom_GLOID,ierr)

   if (ice_present) then
      do eii = 1,num_inst_ice
         if (iamin_CPLICEID(eii)) then
            call seq_map_map(mapper_Cx2i(eii), dom_ix%data, dom_ii(eii)%data, msgtag=CPLICEID(eii)*100+eii*10+5)
            if (iamin_ICEID(eii)) then
               call seq_domain_areafactinit(cdata_ii(eii),areacor_ii(eii)%mdl2drv,areacor_ii(eii)%drv2mdl,&
                    'areafact_i_'//trim(ice_name(eii)))
               call mct_avect_vecmult(i2x_ii(eii),areacor_ii(eii)%mdl2drv,seq_flds_i2x_fluxes)
            endif
            call seq_map_map(mapper_Ci2x(eii), i2x_ii(eii), i2x_ix(eii), msgtag=CPLICEID(eii)*100+eii*10+7)
         endif
      enddo
   endif

   call mpi_barrier(mpicom_GLOID,ierr)

   if (glc_present) then
      do egi = 1,num_inst_glc
         if (iamin_CPLGLCID(egi)) then
            call seq_map_map(mapper_Cx2g(egi), dom_gx%data, dom_gg(egi)%data, msgtag=CPLGLCID(egi)*100+egi*10+5)
            if (iamin_GLCID(egi)) then
               call seq_domain_areafactinit(cdata_gg(egi),areacor_gg(egi)%mdl2drv,areacor_gg(egi)%drv2mdl,&
                    'areafact_g_'//trim(glc_name(egi)))
               call mct_avect_vecmult(g2x_gg(egi),areacor_gg(egi)%mdl2drv,seq_flds_g2x_fluxes)
            endif
            call seq_map_map(mapper_Cg2x(egi), g2x_gg(egi), g2x_gx(egi), msgtag=CPLGLCID(egi)*100+egi*10+7)
         endif
      enddo
   endif

   call mpi_barrier(mpicom_GLOID,ierr)

   if (wav_present) then
      do ewi = 1,num_inst_wav
         if (iamin_CPLWAVID(ewi)) then
            call seq_map_map(mapper_Cx2w(ewi), dom_wx%data, dom_ww(ewi)%data, msgtag=CPLWAVID(ewi)*100+ewi*10+5)
            if (iamin_WAVID(ewi)) then
               call seq_domain_areafactinit(cdata_ww(ewi),areacor_ww(ewi)%mdl2drv,areacor_ww(ewi)%drv2mdl,&
                    'areafact_w_'//trim(wav_name(ewi)))
               call mct_avect_vecmult(w2x_ww(ewi),areacor_ww(ewi)%mdl2drv,seq_flds_w2x_fluxes)
            endif
            call seq_map_map(mapper_Cw2x(ewi), w2x_ww(ewi), w2x_wx(ewi), msgtag=CPLWAVID(ewi)*100+ewi*10+7)
         endif
      enddo
   endif

   !-----------------------------------------------------------------------------
   ! global sum diagnostics for IC data
   !-----------------------------------------------------------------------------
   if (iamin_CPLID .and. info_debug > 1) then
      if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
      if (atm_present) then
         do eai = 1,num_inst_atm
            call seq_diag_avect_mct(cdata_ax,a2x_ax(eai),'recv atm'//trim(atm_suffix(eai))//' IC')
         enddo
      endif

      if (ice_present) then
         do eii = 1,num_inst_ice
            call seq_diag_avect_mct(cdata_ix,i2x_ix(eii),'recv ice'//trim(ice_suffix(eii))//' IC')
         enddo
      endif

      if (lnd_present .or. sno_present) then
         do eli = 1,num_inst_lnd
            if (lnd_present) then 
               call seq_diag_avect_mct(cdata_lx, l2x_lx(eli),'recv lnd'//trim(lnd_suffix(eli))//' IC')
            end if
            if (sno_present) then
               call seq_diag_avect_mct(cdata_sx, s2x_sx(eli),'recv sno'//trim(lnd_suffix(eli))//' IC')
            end if
         enddo
      end if

      if (rof_present) then
         do eri = 1,num_inst_rof
            call seq_diag_avect_mct(cdata_rx, r2x_rx(eri),'recv roff'//trim(lnd_suffix(eri))//' IC')
         end do
      end if

      if (ocn_present) then
         do eoi = 1,num_inst_ocn
            call seq_diag_avect_mct(cdata_ox,o2x_ox(eoi),'recv ocn'//trim(ocn_suffix(eoi))//' IC')
         enddo
      endif

      if (glc_present) then
         do egi = 1,num_inst_glc
            call seq_diag_avect_mct(cdata_gx,g2x_gx(egi),'recv glc'//trim(glc_suffix(egi))//' IC')
         enddo
      endif

      if (wav_present) then
         do ewi = 1,num_inst_wav
            call seq_diag_avect_mct(cdata_wx,w2x_wx(ewi),'recv wav'//trim(wav_suffix(ewi))//' IC')
         enddo
      endif

      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
   end if

   !-----------------------------------------------------------------------------
   ! Initialize fractions
   !-----------------------------------------------------------------------------

   if (iamin_CPLID) then
      if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
      do efi = 1,num_inst_frc
         eii = mod((efi-1),num_inst_ice) + 1
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID .and. efi == 1) write(logunit,F00) 'Initializing fractions'
         call seq_frac_init(cdata_ax, cdata_ix, cdata_lx, cdata_ox, cdata_gx, cdata_rx, &
                         cdata_wx, &
                         ice_present, ocn_present, lnd_present, glc_present, rof_present, &
                         wav_present, dead_comps, &
                         fractions_ax(efi), fractions_ix(efi), fractions_lx(efi), &
                         fractions_ox(efi), fractions_gx(efi), fractions_rx(efi), &
                         fractions_wx(efi), &
                         mapper_Fi2a, mapper_SFo2i, mapper_SFi2o, mapper_Fo2a, &
                         mapper_Fa2o, mapper_Fa2l, mapper_Fl2a, mapper_Fl2r, mapper_Fr2l)
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID .and. efi == 1) write(logunit,F00) 'Setting fractions'
         call seq_frac_set(i2x_ix(eii), &
                        cdata_ax, cdata_ix, cdata_lx, cdata_ox, cdata_gx, cdata_rx, cdata_wx, &
                        ice_present, ocn_present, lnd_present, glc_present, rof_present, wav_present, &
                        fractions_ax(efi), fractions_ix(efi), fractions_lx(efi), &
                        fractions_ox(efi), fractions_gx(efi), fractions_rx(efi), &
                        fractions_wx(efi), &
                        mapper_Fi2a, mapper_SFi2o)
      enddo

      !-----------------------------------------------------------------------------
      ! Initialize atm/ocn flux component and compute ocean albedos
      !-----------------------------------------------------------------------------
      if (ocn_present) then
         if (iamroot_CPLID) write(logunit,*) ' '
         if (iamroot_CPLID) write(logunit,F00) 'Initializing atm/ocn flux component'
         ! note: albedo_only mode doesn't use a2x_ox(eai) or o2x_ox(eoi) or a2x_ax(eai) or o2x_ax(eoi)
         ! Initialize attribute vector
         if (trim(aoflux_grid) == 'ocn') then
            call seq_flux_init_mct(cdata_ox,fractions_ox(ens1))
         elseif (trim(aoflux_grid) == 'atm') then
            call seq_flux_init_mct(cdata_ax,fractions_ax(ens1))
         elseif (trim(aoflux_grid) == 'exch') then
            call seq_flux_initexch_mct(cdata_ax,cdata_ox)
         endif
         do exi = 1,num_inst_xao
!tcx is this correct? relation between xao and frc for ifrad and ofrad
            efi = mod((exi-1),num_inst_frc) + 1
            call seq_flux_ocnalb_mct(cdata_ox,xao_ox(exi),fractions_ox(efi))
         enddo
      endif

      if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
   endif

   !-----------------------------------------------------------------------------
   ! Recalculate initial solar. Merge atmosphere input state and run atmospheric radiation
   ! tcx - for initialization only?
   !-----------------------------------------------------------------------------

   call t_startf('driver_init_atminit')

   if (atm_prognostic) then
      if (iamin_CPLID) then
         if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
         if (lnd_present) then
            if (iamroot_CPLID) write(logunit,F00) 'Calling map_lnd2atm_mct'
            do eli = 1,num_inst_lnd
               efi = mod((eli-1),num_inst_frc) + 1
               call seq_map_map(mapper_Sl2a,l2x_lx(eli),l2x_ax(eli), &
                                fldlist=seq_flds_l2x_states, norm=.true., &
                                avwts_s=fractions_lx(efi),avwtsfld_s='lfrin')
               call seq_map_map(mapper_Fl2a,l2x_lx(eli),l2x_ax(eli), &
                                fldlist=seq_flds_l2x_fluxes, norm=.true., &
                                avwts_s=fractions_lx(efi),avwtsfld_s='lfrin')
            enddo
         endif
         if (ocn_present) then
            if (iamroot_CPLID) write(logunit,F00) 'Calling map_ocn2atm_mct for mapping o2x_ox(eoi) to o2x_ax(eoi)'
            do eoi = 1,num_inst_ocn
               efi = mod((eoi-1),num_inst_frc) + 1
               call seq_map_map(mapper_So2a,o2x_ox(eoi),o2x_ax(eoi),fldlist=seq_flds_o2x_states,norm=.true., &
                                avwts_s=fractions_ox(efi),avwtsfld_s='ofrac')
               call seq_map_map(mapper_Fo2a,o2x_ox(eoi),o2x_ax(eoi),fldlist=seq_flds_o2x_fluxes,norm=.true.)
            enddo
            do exi = 1,num_inst_xao
               efi = mod((exi-1),num_inst_frc) + 1
               call seq_map_map(mapper_So2a,xao_ox(exi),xao_ax(exi),fldlist=seq_flds_xao_albedo,norm=.true., &
                                avwts_s=fractions_ox(efi),avwtsfld_s='ofrac')
               if (trim(aoflux_grid) == 'ocn') then
                  if (iamroot_CPLID .and. exi == 1) &
                     write(logunit,F00) 'Calling map_ocn2atm_mct for mapping xao_ox to xao_ax'
                  call seq_map_map(mapper_So2a,xao_ox(exi),xao_ax(exi),fldlist=seq_flds_xao_states,norm=.true., &
                                   avwts_s=fractions_ox(efi),avwtsfld_s='ofrac')
                  call seq_map_map(mapper_Fo2a,xao_ox(exi),xao_ax(exi),fldlist=seq_flds_xao_fluxes,norm=.true., &
                                   avwts_s=fractions_ox(efi),avwtsfld_s='ofrac')
               endif
               if (trim(aoflux_grid) == 'atm') then
                  if (iamroot_CPLID .and. exi == 1) &
                       write(logunit,F00) 'Calling map_atm2ocn_mct for mapping xao_ax to xao_ox'
                  ! tcraig: this mapping has to be done with area overlap mapping for all fields 
                  ! due to the masking of the xao_ax data and the fact that states are mapped with 
                  ! bilinear mapping currently
                  call seq_map_map(mapper_Fa2o, xao_ax(exi), xao_ox(exi), norm=.true.)
               endif
            enddo
         endif
         if (ice_present) then
            if (iamroot_CPLID) write(logunit,F00) 'Calling map_ice2atm_mct for mapping i2x_ix(eii) to i2x_ax(eii)'
            do eii = 1,num_inst_ice
               efi = mod((eii-1),num_inst_frc) + 1
               call seq_map_map(mapper_Si2a, i2x_ix(eii), i2x_ax(eii), fldlist=seq_flds_i2x_states, &
                    avwts_s=fractions_ix(eii), avwtsfld_s='ifrac')
               call seq_map_map(mapper_Fi2a, i2x_ix(eii), i2x_ax(eii), fldlist=seq_flds_i2x_fluxes, &
                    avwts_s=fractions_ix(eii), avwtsfld_s='ifrac')
            enddo
         endif
         if (lnd_present .or. ocn_present) then
            if (iamroot_CPLID) write(logunit,F00) 'Calling mrg_x2a_run_mct'
            ! Use fortran mod to address ensembles in merge
            do eai = 1,num_inst_atm
               eli = mod((eai-1),num_inst_lnd) + 1
               eoi = mod((eai-1),num_inst_ocn) + 1
               eii = mod((eai-1),num_inst_ice) + 1
               exi = mod((eai-1),num_inst_xao) + 1
               efi = mod((eai-1),num_inst_frc) + 1
               call mrg_x2a_run_mct( cdata_ax, l2x_ax(eli), o2x_ax(eoi), xao_ax(exi), i2x_ax(eii), fractions_ax(efi), x2a_ax(eai))
            enddo
         endif
   
         if (info_debug > 1) then
            do eai = 1,num_inst_atm
               call seq_diag_avect_mct(cdata_ax,x2a_ax(eai),'send atm'//trim(atm_suffix(eai))//' IC2')
            enddo
         endif
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
      endif

      do eai = 1,num_inst_atm
         if (iamin_CPLATMID(eai)) then
            call seq_map_map(mapper_Cx2a(eai), x2a_ax(eai), x2a_aa(eai), msgtag=CPLATMID(eai)*100+eai*10+2)
         endif
      enddo
      if (iamin_CPLALLATMID) then
         call seq_infodata_exchange(infodata,CPLALLATMID,'cpl2atm_init')
      endif
   endif  ! atm_prognostic

   !-----------------------------------------------------------------------------
   ! Second phase of atmosphere component initialization, recalculate solar based 
   ! on input albedo's from surface components. Data or dead atmosphere may just
   ! return on this phase.
   !-----------------------------------------------------------------------------

   if (atm_present) then
      if (iamroot_CPLID) write(logunit,F00) 'Calling atm_init_mct phase 2'
      do eai = 1,num_inst_atm
         if (iamin_ATMID(eai)) then
            call t_adj_detailf(+2)
            if (drv_threading) call seq_comm_setnthreads(nthreads_ATMID)
            if (iamroot_ATMID(eai)) write(logunit,F00) 'Initialize atm component phase 2 '//trim(atm_name(eai))
            call seq_infodata_putData(infodata,atm_phase=2)
            call mct_avect_vecmult(x2a_aa(eai),areacor_aa(eai)%drv2mdl,seq_flds_x2a_fluxes)
            call atm_init_mct( EClock_a, cdata_aa(eai), x2a_aa(eai), a2x_aa(eai))
            call mct_avect_vecmult(a2x_aa(eai),areacor_aa(eai)%mdl2drv,seq_flds_a2x_fluxes)
            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_adj_detailf(-2)
         endif
      enddo
      do eai = 1,num_inst_atm
         if (iamin_CPLATMID(eai)) then
            call seq_map_map(mapper_Ca2x(eai), a2x_aa(eai), a2x_ax(eai), msgtag=CPLATMID(eai)*100+eai*10+4)
         endif
      enddo
      if (iamin_CPLATMID(ens1)) then
         call seq_infodata_exchange(infodata,CPLATMID(ens1),'atm2cpl_init')
      endif

      if (iamin_CPLID) then
         if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
         if (info_debug > 1) then
            do eai = 1,num_inst_atm
               call seq_diag_avect_mct(cdata_ax,a2x_ax(eai),'recv atm'//trim(atm_suffix(eai))//' IC2')
            enddo
         endif
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
      endif
   endif   ! atm present

   call t_stopf  ('driver_init_atminit')

   !-----------------------------------------------------------------------------
   ! Read driver restart file, overwrite anything previously sent or computed
   !-----------------------------------------------------------------------------

   call t_startf('driver_init_readrestart')
   call seq_diag_zero_mct(mode='all')
   if (read_restart) call seq_rest_read(rest_file)
   call t_stopf  ('driver_init_readrestart')

   if (do_histinit) then
      if (iamin_CPLID) then
         if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
         if (iamroot_CPLID) then
            call seq_timemgr_EClockGetData( EClock_d, curr_ymd=ymd, curr_tod=tod )
            write(logunit,104) ' Write history file at ',ymd,tod
            call shr_sys_flush(logunit)
         endif
         call seq_hist_write(EClock_d)
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
      endif
   endif

   if (iamroot_CPLID )then
      write(logunit,*) ' '
      write(logunit,F00) 'Model initialization complete '
      write(logunit,*) ' '
      call shr_sys_flush(logunit)
   endif

   call t_stopf  ('DRIVER_INIT')

end subroutine ccsm_init


!===============================================================================

subroutine ccsm_run()

   implicit none

 101  format( A, 2i8, 12A, A, F8.2, A, F8.2 )
 102  format( A, 2i8, A, 8L3 )
 103  format( 5A )
 104  format( A, 2i8)
 105  format( A, 2i8, A, f10.2, A, f10.2, A, A, i5, A, A)
 106  format( A, f23.12)
 107  format( A, 2i8, A, f12.4, A, f12.4 )

   call seq_infodata_putData(infodata,atm_phase=1,lnd_phase=1,ocn_phase=1,ice_phase=1)
   call seq_timemgr_EClockGetData( EClock_d, stepno=begstep)
   call seq_timemgr_EClockGetData( EClock_d, dtime=dtime)
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

   !----------------------------------------------------------
   ! Beginning of basic time step loop
   !----------------------------------------------------------

   call t_startf ('DRIVER_RUN_LOOP_BSTART')
   call mpi_barrier(mpicom_GLOID,ierr)
   call t_stopf ('DRIVER_RUN_LOOP_BSTART')
   Time_begin = mpi_wtime()
   Time_bstep = mpi_wtime()
   do while ( .not. stop_alarm)

      call t_startf('DRIVER_RUN_LOOP')
      call t_drvstartf ('DRIVER_CLOCK_ADVANCE',cplrun=.true.)

      !----------------------------------------------------------
      ! Advance sync clock time (this is time that models should have before 
      ! they return to the driver).  Write timestamp and run alarm status
      !----------------------------------------------------------

      call seq_timemgr_clockAdvance( seq_SyncClock)
      call seq_timemgr_EClockGetData( EClock_d, curr_ymd=ymd, curr_tod=tod )
      call shr_cal_date2ymd(ymd,year,month,day)
      stop_alarm    = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_stop)
      atmrun_alarm  = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_atmrun)
      lndrun_alarm  = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_lndrun)
      rofrun_alarm  = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_rofrun)
      icerun_alarm  = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_icerun)
      glcrun_alarm  = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_glcrun)
      wavrun_alarm  = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_wavrun)
      ocnrun_alarm  = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_ocnrun)
      ocnnext_alarm = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_ocnnext)
      restart_alarm = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_restart)
      history_alarm = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_history)
      histavg_alarm = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_histavg)
      tprof_alarm   = seq_timemgr_alarmIsOn(EClock_d,seq_timemgr_alarm_tprof)

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
                         rofrun_alarm, wavrun_alarm
            write(logunit,102) ' Alarm_state: model date = ',ymd,tod, &
               ' 1.2.3.6.12.24 run alarms = ',  t1hr_alarm, t2hr_alarm, &
                         t3hr_alarm, t6hr_alarm, t12hr_alarm, t24hr_alarm
            call shr_sys_flush(logunit)
         endif
      endif

      call t_drvstopf  ('DRIVER_CLOCK_ADVANCE',cplrun=.true.)

      !----------------------------------------------------------
      ! OCN/ICE PREP
      ! Map for ice prep and atmocn flux
      !----------------------------------------------------------

      if (iamin_CPLID .and. (ice_present.or.ocn_present) .and. atm_present) then
         if (run_barriers) then
            call t_drvstartf ('DRIVER_OCNPREP_BARRIER')
            call mpi_barrier(mpicom_CPLID,ierr)
            call t_drvstopf ('DRIVER_OCNPREP_BARRIER')
         endif
         call t_drvstartf ('DRIVER_OCNPREP',cplrun=.true.,barrier=mpicom_CPLID)
         if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
         call t_drvstartf ('driver_ocnprep_atm2ocn',barrier=mpicom_CPLID)
         do eai = 1,num_inst_atm
            call seq_map_map(mapper_Sa2o, a2x_ax(eai), a2x_ox(eai), fldlist=seq_flds_a2x_states, norm=.true.)
            call seq_map_map(mapper_Fa2o, a2x_ax(eai), a2x_ox(eai), fldlist=seq_flds_a2x_fluxes, norm=.true.)
            !--- tcx this Va2o call will not be necessary when npfix goes away
            call seq_map_map(mapper_Va2o, a2x_ax(eai), a2x_ox(eai), fldlist='Sa_u:Sa_v', norm=.true.)
            !--- tcx the norm should be true below, it's false for bfb backwards compatability
            call seq_map_mapvect(mapper_Va2o,vect_map,a2x_ax(eai),a2x_ox(eai),'Sa_u','Sa_v',norm=.false.)
         enddo
         call t_drvstopf  ('driver_ocnprep_atm2ocn')
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
         call t_drvstopf  ('DRIVER_OCNPREP',cplrun=.true.)
      endif

      !----------------------------------------------------------
      ! OCN SETUP
      !----------------------------------------------------------

      if (ocn_present .and. ocnrun_alarm) then

         !----------------------------------------------------
         ! "startup" wait
         !----------------------------------------------------

         if (iamin_CPLALLOCNID .and. cpl2ocn_first) then
            ! want to know the time the ocean pes waited for the cpl pes
            !   at the first ocnrun_alarm, min ocean wait is wait time
            ! do not use t_barrierf here since it can be "off", use mpi_barrier
            do eoi = 1,num_inst_ocn
               if (iamin_OCNID(eoi)) call t_drvstartf ('DRIVER_C2O_INITWAIT')
            enddo
            call mpi_barrier(mpicom_CPLALLOCNID,ierr)
            do eoi = 1,num_inst_ocn
               if (iamin_OCNID(eoi)) call t_drvstopf  ('DRIVER_C2O_INITWAIT')
            enddo
            cpl2ocn_first = .false.
         endif

         !----------------------------------------------------
         ! ocn prep
         ! note due to x2oacc and r2xacc, need to be careful filling x2o_oo
         ! average x2oacc and r2xacc first, separately, now averages in those AVs
         ! map r2xacc_rx to r2x_ox then copy r2x_ox into x2oacc_ox
         ! finally, rearrange x2oacc_ox to x2o_oo
         !----------------------------------------------------

         if (iamin_CPLID .and. ocn_prognostic) then
            if (run_barriers) then
               call t_drvstartf ('DRIVER_OCNPREP_BARRIER')
               call mpi_barrier(mpicom_CPLID,ierr)
               call t_drvstopf ('DRIVER_OCNPREP_BARRIER')
            endif
            call t_drvstartf ('DRIVER_OCNPREP',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
            call t_drvstartf ('driver_ocnprep_avg',barrier=mpicom_CPLID)
            do eoi = 1,num_inst_ocn
               ! temporary formation of average
!              call mct_aVect_average(x2oacc_ox(eoi))
               if (x2oacc_ox_cnt > 0) then
                  x2oacc_ox(eoi)%rAttr = x2oacc_ox(eoi)%rAttr / (x2oacc_ox_cnt*1.0_r8)
               endif
            enddo
            x2oacc_ox_cnt = 0
            call t_drvstopf  ('driver_ocnprep_avg')
            if (rof_present .and. ocnrof_prognostic) then
               ! Map runoff to ocn, average, put in x2oacc_ox(eoi)
               if (r2xacc_rx_cnt > 0) then
                  call t_drvstartf ('driver_ocnprep_ravg',barrier=mpicom_CPLID)
                  do eri = 1,num_inst_rof
                     r2xacc_rx(eri)%rAttr = r2xacc_rx(eri)%rAttr / (r2xacc_rx_cnt*1.0_r8)
                  enddo
                  r2xacc_rx_cnt = 0
                  call t_drvstopf ('driver_ocnprep_ravg')

                  call t_drvstartf ('driver_ocnprep_rof2ocn',barrier=mpicom_CPLID)
                  do eri = 1,num_inst_rof
                     call seq_map_map(mapper_Rr2o, r2xacc_rx(eri), r2x_ox(eri), &
                                      fldlist='Forr_roff:Forr_ioff', norm=.false.)
                     if (flood_present) then
                        call seq_map_map(mapper_Fr2o, r2xacc_rx(eri), r2x_ox(eri), &
                                         fldlist='Flrr_flood', norm=.true.)
                        ! add flood to roff and zero flood
                        r2x_ox(eri)%rAttr(index_r2x_Forr_roff,:) = r2x_ox(eri)%rAttr(index_r2x_Forr_roff,:) + &
                                                                   r2x_ox(eri)%rAttr(index_r2x_Flrr_flood,:)
                        r2x_ox(eri)%rAttr(index_r2x_Flrr_flood,:) = 0.0_r8
                     endif
                  enddo
                  if (do_hist_r2x) then
                     do eri = 1,num_inst_rof
                        call seq_hist_writeaux(EClock_d,'r2xacc'//trim(lnd_suffix(eri)), &
                                'domr',cdata_rx, r2xacc_rx(eri), rof_nx,rof_ny,1)
                     enddo
                  endif
                  call t_drvstopf  ('driver_ocnprep_rof2ocn')

                  ! Use fortran mod to address ensembles in merge
                  call t_drvstartf ('driver_ocnprep_rofcopy',barrier=mpicom_CPLID)
                  do eoi = 1,num_inst_ocn
                     eri = mod((eoi-1),num_inst_rof) + 1
                     ! The following copies field Forr_roff and Forr_ioff 
                     ! from r2x_ox to x2oacc_ox
                     call mct_aVect_copy(aVin=r2x_ox(eri), aVout=x2oacc_ox(eoi))
                     call mct_aVect_copy(aVin=r2x_ox(eri), aVout=x2o_ox(eoi))   ! for history file
                  enddo
                  call t_drvstopf  ('driver_ocnprep_rofcopy')
               endif
            endif
            if (info_debug > 1) then
               call t_drvstartf ('driver_ocnprep_diagav',barrier=mpicom_CPLID)
               do eoi = 1,num_inst_ocn
                  call seq_diag_avect_mct(cdata_ox,x2oacc_ox(eoi),'send ocn'//trim(ocn_suffix(eoi)))
               enddo
               call t_drvstopf  ('driver_ocnprep_diagav')
            endif
            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('DRIVER_OCNPREP',cplrun=.true.)
         endif

         !----------------------------------------------------
         ! cpl -> ocn
         !----------------------------------------------------
         if (iamin_CPLALLOCNID .and. ocn_prognostic) then
            if (run_barriers) then
               call t_drvstartf ('DRIVER_C2O_BARRIER')
               call mpi_barrier(mpicom_CPLALLOCNID,ierr)
               call t_drvstopf ('DRIVER_C2O_BARRIER')
            endif
            call t_drvstartf ('DRIVER_C2O',cplcom=.true.,barrier=mpicom_CPLALLOCNID)
            do eoi = 1,num_inst_ocn
               if (iamin_CPLOCNID(eoi)) then
                  call t_drvstartf ('driver_c2o_ocnx2ocno',barrier=mpicom_CPLOCNID(eoi))
                  call seq_map_map(mapper_Cx2o(eoi), x2oacc_ox(eoi), x2o_oo(eoi), msgtag=CPLOCNID(eoi)*100+eoi*10+2)
                  call t_drvstopf  ('driver_c2o_ocnx2ocno')
               endif
            enddo
            call t_drvstartf ('driver_c2o_infoexch',barrier=mpicom_CPLALLOCNID)
            call seq_infodata_exchange(infodata,CPLALLOCNID,'cpl2ocn_run')
            call t_drvstopf  ('driver_c2o_infoexch')
            call t_drvstopf  ('DRIVER_C2O',cplcom=.true.)
         endif
      endif
  
      !----------------------------------------------------------
      ! LND SETUP
      !----------------------------------------------------------

      if (lnd_present .and. lndrun_alarm) then

         !----------------------------------------------------
         ! lnd prep
         !----------------------------------------------------

         if (iamin_CPLID) then
            if (run_barriers) then
               call t_drvstartf ('DRIVER_LNDPREP_BARRIER')
               call mpi_barrier(mpicom_CPLID,ierr)
               call t_drvstopf ('DRIVER_LNDPREP_BARRIER')
            endif
            call t_drvstartf ('DRIVER_LNDPREP',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)

            if (lnd_prognostic) then
               call t_drvstartf ('driver_lndprep_atm2lnd',barrier=mpicom_CPLID)
               do eai = 1,num_inst_atm
                  call seq_map_map(mapper_Fa2l, a2x_ax(eai), a2x_lx(eai), norm=.true.)
               enddo
               call t_drvstopf  ('driver_lndprep_atm2lnd')
               if (flood_present) then
                  call t_drvstartf ('driver_lndprep_rof2lnd',barrier=mpicom_CPLID)
                  ! Obtain flooding input from rof to be sent back to land 
                  do eri = 1,num_inst_rof
                     efi = mod((eri-1),num_inst_frc) + 1
                     call seq_map_map(mapper_Fr2l, r2x_rx(eri), r2x_lx(eri), &
                          fldlist=seq_flds_r2x_fluxes, norm=.true.)
                  end do
                  call t_drvstopf  ('driver_lndprep_rof2lnd')
               end if

               if (rof_present) then
                  call t_drvstartf ('driver_lndprep_rof2lnd',barrier=mpicom_CPLID)
                  ! Obtain volr from rof to be sent back to land 
                  do eri = 1,num_inst_rof
                     efi = mod((eri-1),num_inst_frc) + 1
                     call seq_map_map(mapper_Sr2l,r2x_rx(eri),r2x_lx(eri), &
                                fldlist=seq_flds_r2x_states, norm=.true.)
                  end do
                  call t_drvstopf  ('driver_lndprep_rof2lnd')
               end if

               call t_drvstartf ('driver_lndprep_mrgx2l',barrier=mpicom_CPLID)
               ! Use fortran mod to address ensembles in merge
               do eli = 1,num_inst_lnd
                  eai = mod((eli-1),num_inst_atm) + 1
                  eri = mod((eri-1),num_inst_rof) + 1
                  call mrg_x2l_run_mct( cdata_lx, a2x_lx(eai), r2x_lx(eri), x2l_lx(eli))
               enddo
               call t_drvstopf  ('driver_lndprep_mrgx2l')

               if (info_debug > 1) then
                  call t_drvstartf ('driver_lndprep_diagav',barrier=mpicom_CPLID)
                  do eli = 1,num_inst_lnd
                     call seq_diag_avect_mct(cdata_lx,x2l_lx(eli),'send lnd'//trim(lnd_suffix(eli)))
                  enddo
                  call t_drvstopf  ('driver_lndprep_diagav')
               endif
            endif

            if (glc_present .and. sno_prognostic) then
               if (info_debug > 1) then
                  call t_drvstartf ('driver_lndprep_diagav',barrier=mpicom_CPLID)
                  do eli = 1,num_inst_lnd
                     call seq_diag_avect_mct(cdata_sx,x2s_sx(eli),'send sno'//trim(lnd_suffix(eli)))
                  enddo
                  call t_drvstopf  ('driver_lndprep_diagav')
               endif
            endif

            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('DRIVER_LNDPREP',cplrun=.true.)
         endif

         !----------------------------------------------------
         ! cpl -> lnd
         !----------------------------------------------------

         if (iamin_CPLALLLNDID) then
            if (run_barriers) then
               call t_drvstartf ('DRIVER_C2L_BARRIER')
               call mpi_barrier(mpicom_CPLALLLNDID,ierr)
               call t_drvstopf ('DRIVER_C2L_BARRIER')
            endif
            call t_drvstartf ('DRIVER_C2L',cplcom=.true.,barrier=mpicom_CPLALLLNDID)
            do eli = 1,num_inst_lnd
               if (iamin_CPLLNDID(eli)) then
                  if (lnd_prognostic) then
                     call t_drvstartf ('driver_c2l_lndx2lndl', &
                                       barrier=mpicom_CPLLNDID(eli))
                     call seq_map_map(mapper_Cx2l(eli), x2l_lx(eli), x2l_ll(eli), msgtag=CPLLNDID(eli)*100+eli*10+2)
                     call t_drvstopf  ('driver_c2l_lndx2lndl')
                  endif

                  if (glc_present .and. sno_prognostic) then
                     call t_drvstartf ('driver_c2l_snox2snos', &
                                       barrier=mpicom_CPLLNDID(eli))
                     call seq_map_map(mapper_Cx2s(eli), x2s_sx(eli), x2s_ss(eli), msgtag=CPLLNDID(eli)*100+eli*10+4)
                     call t_drvstopf  ('driver_c2l_snox2snos')
                  endif
               endif
            enddo
            if (lnd_prognostic .or. sno_prognostic) then
               call t_drvstartf ('driver_c2l_infoexch',barrier=mpicom_CPLALLLNDID)
               call seq_infodata_exchange(infodata,CPLALLLNDID,'cpl2lnd_run')
               call t_drvstopf  ('driver_c2l_infoexch')
            endif
            call t_drvstopf  ('DRIVER_C2L',cplcom=.true.)
         endif


      endif

      !----------------------------------------------------------
      ! ICE SETUP
      ! Note that for atm->ice mapping below will leverage the assumption that the
      ! ice and ocn are on the same grid and that mapping of atm to ocean is 
      ! done already for use by atmocn flux and ice model prep
      !----------------------------------------------------------

      if (ice_present .and. icerun_alarm) then

         !----------------------------------------------------
         ! ice prep
         !----------------------------------------------------

         if (iamin_CPLID .and. ice_prognostic) then
            if (run_barriers) then
               call t_drvstartf ('DRIVER_ICEPREP_BARRIER')
               call mpi_barrier(mpicom_CPLID,ierr)
               call t_drvstopf ('DRIVER_ICEPREP_BARRIER')
            endif
            call t_drvstartf ('DRIVER_ICEPREP',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
            call t_drvstartf ('driver_iceprep_ocn2ice',barrier=mpicom_CPLID)
            do eoi = 1,num_inst_ocn
               call seq_map_map(mapper_SFo2i, o2x_ox(eoi), o2x_ix(eoi), norm=.true.)
            enddo
            call t_drvstopf  ('driver_iceprep_ocn2ice')
            
            call t_drvstartf ('driver_iceprep_atm2ice',barrier=mpicom_CPLID)
            do eai = 1,num_inst_atm
               call seq_map_map(mapper_SFo2i, a2x_ox(eai), a2x_ix(eai), norm=.true.)
!tcx fails     call map_atm2ice_mct( cdata_ax, a2x_ax(eai), cdata_ix, a2x_ix(eai) )
            enddo
            call t_drvstopf  ('driver_iceprep_atm2ice')
            
            call t_drvstartf ('driver_iceprep_mrgx2i',barrier=mpicom_CPLID)
            ! Use fortran mod to address ensembles in merge
            do eii = 1,num_inst_ice
               eai = mod((eii-1),num_inst_atm) + 1
               eoi = mod((eii-1),num_inst_ocn) + 1
               call mrg_x2i_run_mct( cdata_ix, a2x_ix(eai), o2x_ix(eoi), x2i_ix(eii) )
            enddo
            call t_drvstopf  ('driver_iceprep_mrgx2i')

            if (info_debug > 1) then
               call t_drvstartf ('driver_iceprep_diagav',barrier=mpicom_CPLID)
               do eii = 1,num_inst_ice
                  call seq_diag_avect_mct(cdata_ix,x2i_ix(eii),'send ice'//trim(ice_suffix(eii)))
               enddo
               call t_drvstopf  ('driver_iceprep_diagav')
            endif
            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('DRIVER_ICEPREP',cplrun=.true.)
         endif

         !----------------------------------------------------
         ! cpl -> ice
         !----------------------------------------------------

         if (iamin_CPLALLICEID .and. ice_prognostic) then
            if (run_barriers) then
               call t_drvstartf ('DRIVER_C2I_BARRIER')
               call mpi_barrier(mpicom_CPLALLICEID,ierr)
               call t_drvstopf ('DRIVER_C2I_BARRIER')
            endif
            call t_drvstartf ('DRIVER_C2I',cplcom=.true.,barrier=mpicom_CPLALLICEID)
            do eii = 1,num_inst_ice
               if (iamin_CPLICEID(eii)) then
                  call t_drvstartf ('driver_c2i_icex2icei',barrier=mpicom_CPLICEID(eii))
                  call seq_map_map(mapper_Cx2i(eii), x2i_ix(eii), x2i_ii(eii), msgtag=CPLICEID(eii)*100+eii*10+2)
                  call t_drvstopf  ('driver_c2i_icex2icei')
               endif
            enddo
            call t_drvstartf ('driver_c2i_infoexch',barrier=mpicom_CPLALLICEID)
            call seq_infodata_exchange(infodata,CPLALLICEID,'cpl2ice_run')
            call t_drvstopf  ('driver_c2i_infoexch')
            call t_drvstopf  ('DRIVER_C2I',cplcom=.true.)
         endif

      endif

      !----------------------------------------------------------
      ! WAV Setup
      !----------------------------------------------------------

      if (wav_present .and. wavrun_alarm) then

         !----------------------------------------------------------
         ! wav Prep
         !----------------------------------------------------------

         if (iamin_CPLID .and. wav_prognostic) then
            if (run_barriers) then
               call t_drvstartf ('DRIVER_WAVPREP_BARRIER')
               call mpi_barrier(mpicom_CPLID,ierr)
               call t_drvstopf ('DRIVER_WAVPREP_BARRIER')
            endif
            call t_drvstartf ('DRIVER_WAVPREP',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
            call t_drvstartf ('driver_wavprep_atm2wav',barrier=mpicom_CPLID)
            do eai = 1,num_inst_atm
               call seq_map_map(mapper_Sa2w, a2x_ax(eai), a2x_wx(eai), norm=.true.)
            enddo
            call t_drvstopf  ('driver_wavprep_atm2wav')

            call t_drvstartf ('driver_wavprep_ocn2wav',barrier=mpicom_CPLID)
            do eoi = 1,num_inst_ocn
               call seq_map_map(mapper_So2w, o2x_ox(eoi), o2x_wx(eoi), norm=.true.)
            enddo
            call t_drvstopf  ('driver_wavprep_ocn2wav')

            call t_drvstartf ('driver_wavprep_ice2wav',barrier=mpicom_CPLID)
            do eii = 1,num_inst_ice
               call seq_map_map(mapper_Si2w, i2x_ix(eii), i2x_wx(eii), norm=.true.)
            enddo
            call t_drvstopf  ('driver_wavprep_ice2wav')

            ! Merge wav inputs
            ! Use fortran mod to address ensembles in merge
            call t_drvstartf ('driver_wavprep_mrgx2w',barrier=mpicom_CPLID)
            do ewi = 1,num_inst_wav
               eai = mod((ewi-1),num_inst_atm) + 1
               eoi = mod((ewi-1),num_inst_ocn) + 1
               eii = mod((ewi-1),num_inst_ice) + 1
               efi = mod((ewi-1),num_inst_frc) + 1
               call mrg_x2w_run_mct( cdata_wx, a2x_wx(eai), o2x_wx(eoi), i2x_wx(eii), fractions_wx(efi), x2w_wx(ewi))
            enddo
            call t_drvstopf  ('driver_wavprep_mrgx2w')

            if (info_debug > 1) then
               call t_drvstartf ('driver_wavprep_diagav',barrier=mpicom_CPLID)
               do ewi = 1,num_inst_wav
                  call seq_diag_avect_mct(cdata_wx,x2w_wx(ewi),'send wav'//trim(wav_suffix(ewi)))
               enddo
               call t_drvstopf  ('driver_wavprep_diagav')
            endif
            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('DRIVER_WAVPREP',cplrun=.true.)
         end if

         !----------------------------------------------------------
         ! cpl -> wav
         !----------------------------------------------------------

         if (iamin_CPLALLWAVID .and. wav_prognostic) then
            if (run_barriers) then
               call t_drvstartf ('DRIVER_C2W_BARRIER')
               call mpi_barrier(mpicom_CPLALLWAVID,ierr)
               call t_drvstopf ('DRIVER_C2W_BARRIER')
            endif
            call t_drvstartf ('DRIVER_C2W',cplcom=.true.,barrier=mpicom_CPLALLWAVID)
            do ewi = 1,num_inst_wav
               if (iamin_CPLWAVID(ewi)) then
                  call t_drvstartf ('driver_c2w_wavx2wavw',barrier=mpicom_CPLWAVID(ewi))
                  call seq_map_map(mapper_Cx2w(ewi), x2w_wx(ewi), x2w_ww(ewi))
                  call t_drvstopf  ('driver_c2w_wavx2wavw')
               endif
            enddo
            call t_drvstartf ('driver_c2w_infoexch',barrier=mpicom_CPLALLWAVID)
            call seq_infodata_exchange(infodata,CPLALLWAVID,'cpl2wav_run')
            call t_drvstopf  ('driver_c2w_infoexch')
            call t_drvstopf  ('DRIVER_C2W',cplcom=.true.)
         endif
      end if

      !----------------------------------------------------------
      ! ROF SETUP
      !----------------------------------------------------------

      if (rof_present .and. rofrun_alarm) then

         !----------------------------------------------------
         ! rof prep
         !----------------------------------------------------

         if (iamin_CPLID .and. rof_prognostic) then
            if (run_barriers) then
               call t_drvstartf ('DRIVER_ROFPREP_BARRIER')
               call mpi_barrier(mpicom_CPLID,ierr)
               call t_drvstopf ('DRIVER_ROFPREP_BARRIER')
            endif
            call t_drvstartf ('DRIVER_ROFPREP',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
            do eri = 1,num_inst_rof
               eli = mod((eri-1),num_inst_lnd) + 1
               efi = mod((eri-1),num_inst_frc) + 1
               call t_drvstartf ('driver_rofprep_l2xavg',barrier=mpicom_CPLID)
               if (x2racc_lx_cnt > 0) then
                  x2racc_lx(eli)%rAttr = x2racc_lx(eli)%rAttr / (x2racc_lx_cnt*1.0_r8)
               endif
               call t_drvstopf ('driver_rofprep_l2xavg')

               call t_drvstartf ('driver_rofprep_lnd2rof',barrier=mpicom_CPLID)
               call seq_map_map(mapper_Fl2r, x2racc_lx(eli), x2r_rx_tmp, &
                  fldlist=seq_flds_x2r_fluxes, norm=.true., &
                  avwts_s=fractions_lx(efi), avwtsfld_s='lfrin')
               call t_drvstopf ('driver_rofprep_lnd2rof')

               call t_drvstartf ('driver_rofprep_mrgx2r',barrier=mpicom_CPLID)
               call mrg_x2r_run_mct( cdata_rx, x2r_rx_tmp, fractions_rx(efi), x2r_rx(eri))
               call t_drvstopf ('driver_rofprep_mrgx2r')
            end do
            x2racc_lx_cnt = 0 

            if (info_debug > 1) then
               call t_drvstartf ('driver_rofprep_diagav',barrier=mpicom_CPLID)
               do eri = 1,num_inst_rof
                  call seq_diag_avect_mct(cdata_rx, x2r_rx(eri),'send roff'//trim(rof_suffix(eri)))
               enddo
               call t_drvstopf ('driver_rofprep_diagav')
            end if
            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('DRIVER_ROFPREP',cplrun=.true.)
         end if

         !----------------------------------------------------
         ! cpl -> rof
         !----------------------------------------------------

         if (iamin_CPLALLROFID .and. rof_prognostic) then
            if (run_barriers) then
               call t_drvstartf ('DRIVER_C2R_BARRIER')
               call mpi_barrier(mpicom_CPLALLROFID,ierr)
               call t_drvstopf ('DRIVER_C2R_BARRIER')
            endif
            call t_drvstartf ('DRIVER_C2R',cplcom=.true.,barrier=mpicom_CPLALLLNDID)
            do eri = 1,num_inst_rof
               if (iamin_CPLROFID(eri)) then
                  call t_drvstartf ('driver_c2r_rofx2rofi',barrier=mpicom_CPLROFID(eri))
                  call seq_map_map(mapper_Cx2r(eri), x2r_rx(eri), x2r_rr(eri), msgtag=CPLROFID(eri)*100+eri*10+2)
                  call t_drvstopf  ('driver_c2r_rofx2rofi')
               end if
            end do
            call t_drvstartf ('driver_c2r_infoexch',barrier=mpicom_CPLALLROFID)
            call seq_infodata_exchange(infodata,CPLALLROFID,'cpl2rof_run')
            call t_drvstopf  ('driver_c2r_infoexch')
            call t_drvstopf  ('DRIVER_C2R',cplcom=.true.)
         end if

      end if

      !----------------------------------------------------------
      ! Run Ice Model
      !----------------------------------------------------------

      if (ice_present .and. icerun_alarm) then
         do eii = 1,num_inst_ice
            if (iamin_ICEID(eii)) then
               if (run_barriers) then
                  call t_drvstartf ('DRIVER_ICE_RUN_BARRIER')
                  call mpi_barrier(mpicom_ICEID(eii),ierr)
                  call t_drvstopf ('DRIVER_ICE_RUN_BARRIER')
                  time_brun = mpi_wtime()
               endif
               call t_drvstartf ('DRIVER_ICE_RUN',barrier=mpicom_ICEID(eii))
               if (drv_threading) call seq_comm_setnthreads(nthreads_ICEID)
               if (ice_prognostic) call mct_avect_vecmult(x2i_ii(eii),areacor_ii(eii)%drv2mdl,seq_flds_x2i_fluxes)
               call ice_run_mct( EClock_i, cdata_ii(eii), x2i_ii(eii), i2x_ii(eii))
               call mct_avect_vecmult(i2x_ii(eii),areacor_ii(eii)%mdl2drv,seq_flds_i2x_fluxes)
               if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
               call t_drvstopf  ('DRIVER_ICE_RUN')
               if (run_barriers) then
                  time_erun = mpi_wtime()
                  cktime = time_erun-time_brun
                  cktime_acc(4) = cktime_acc(4) + cktime
                  cktime_cnt(4) = cktime_cnt(4) + 1
                  write(logunit,107) ' rstamp ice_run_time: model date = ', &
                     ymd,tod,' avg dt = ',cktime_acc(4)/cktime_cnt(4),' dt = ',cktime
               endif
            endif
         enddo
      endif

      !----------------------------------------------------------
      ! Run Land Model 
      !----------------------------------------------------------

      if ((lnd_present .or. sno_present) .and. lndrun_alarm) then
         do eli = 1,num_inst_lnd
            if (iamin_LNDID(eli)) then
               if (run_barriers) then
                  call t_drvstartf ('DRIVER_LND_RUN_BARRIER')
                  call mpi_barrier(mpicom_LNDID(eli),ierr)
                  call t_drvstopf ('DRIVER_LND_RUN_BARRIER')
                  time_brun = mpi_wtime()
               endif
               call t_drvstartf ('DRIVER_LND_RUN',barrier=mpicom_LNDID(eli))
               if (drv_threading) call seq_comm_setnthreads(nthreads_LNDID)
               if (lnd_prognostic) call mct_avect_vecmult(x2l_ll(eli),areacor_ll(eli)%drv2mdl,seq_flds_x2l_fluxes)
               if (sno_prognostic) call mct_avect_vecmult(x2s_ss(eli),areacor_ss(eli)%drv2mdl,seq_flds_x2s_fluxes)
               call lnd_run_mct( EClock_l, cdata_ll(eli), x2l_ll(eli), l2x_ll(eli), &
                                           cdata_ss(eli), x2s_ss(eli), s2x_ss(eli))
               if (lnd_present) call mct_avect_vecmult(l2x_ll(eli),areacor_ll(eli)%mdl2drv,seq_flds_l2x_fluxes)
               if (sno_present) call mct_avect_vecmult(s2x_ss(eli),areacor_ss(eli)%mdl2drv,seq_flds_s2x_fluxes)
               if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
               call t_drvstopf  ('DRIVER_LND_RUN')
               if (run_barriers) then
                  time_erun = mpi_wtime()
                  cktime = time_erun-time_brun
                  cktime_acc(3) = cktime_acc(3) + cktime
                  cktime_cnt(3) = cktime_cnt(3) + 1
                  write(logunit,107) ' rstamp lnd_run_time: model date = ', &
                     ymd,tod,' avg dt = ',cktime_acc(3)/cktime_cnt(3),' dt = ',cktime
               endif
            endif   ! iamin_LNDID(eli)
         enddo   ! eli
      endif   ! lnd_present or sno_present and lndrun_alarm

      !----------------------------------------------------------
      ! Run River Runoff model
      !----------------------------------------------------------

      if (rof_present .and. rofrun_alarm) then
         do eri = 1,num_inst_rof
            if (iamin_ROFID(eri)) then
               if (run_barriers) then
                  call t_drvstartf ('DRIVER_ROF_RUN_BARRIER')
                  call mpi_barrier(mpicom_ROFID(eri),ierr)
                  call t_drvstopf ('DRIVER_ROF_RUN_BARRIER')
                  time_brun = mpi_wtime()
               endif
               call t_drvstartf  ('DRIVER_ROF_RUN',barrier=mpicom_ROFID(eri))
               if (drv_threading) call seq_comm_setnthreads(nthreads_ROFID)
               if (rof_prognostic) call mct_avect_vecmult(x2r_rr(eri),areacor_rr(eri)%drv2mdl,seq_flds_x2r_fluxes)
               call rof_run_mct(Eclock_r, cdata_rr(eri), x2r_rr(eri), r2x_rr(eri))
               call mct_avect_vecmult(r2x_rr(eri),areacor_rr(eri)%mdl2drv,seq_flds_r2x_fluxes)
               if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
               call t_drvstopf  ('DRIVER_ROF_RUN')
               if (run_barriers) then
                  time_erun = mpi_wtime()
                  cktime = time_erun-time_brun
                  cktime_acc(7) = cktime_acc(7) + cktime
                  cktime_cnt(7) = cktime_cnt(7) + 1
                  write(logunit,107) ' rstamp rof_run_time: model date = ', &
                     ymd,tod,' avg dt = ',cktime_acc(7)/cktime_cnt(7),' dt = ',cktime
               endif
            end if
         enddo
      end if 

      !----------------------------------------------------------
      ! Run wave model
      !----------------------------------------------------------

      if (wav_present .and. wavrun_alarm) then
         do ewi = 1,num_inst_wav
            if (iamin_WAVID(ewi)) then
               if (run_barriers) then
                  call t_drvstartf ('DRIVER_WAV_RUN_BARRIER')
                  call mpi_barrier(mpicom_WAVID(ewi),ierr)
                  call t_drvstopf ('DRIVER_WAV_RUN_BARRIER')
                  time_brun = mpi_wtime()
               endif
               call t_drvstartf ('DRIVER_WAV_RUN',barrier=mpicom_WAVID(ewi))
               if (drv_threading) call seq_comm_setnthreads(nthreads_WAVID)
               if (wav_prognostic) call mct_avect_vecmult(x2w_ww(ewi),areacor_ww(ewi)%drv2mdl,seq_flds_x2w_fluxes)
               call wav_run_mct( EClock_w, cdata_ww(ewi), x2w_ww(ewi), w2x_ww(ewi))
               call mct_avect_vecmult(w2x_ww(ewi),areacor_ww(ewi)%mdl2drv,seq_flds_w2x_fluxes)
               if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
               call t_drvstopf  ('DRIVER_WAV_RUN')
               if (run_barriers) then
                  time_erun = mpi_wtime()
                  cktime = time_erun-time_brun
                  cktime_acc(8) = cktime_acc(8) + cktime
                  cktime_cnt(8) = cktime_cnt(8) + 1
                  write(logunit,107) ' rstamp wav_run_time: model date = ', &
                     ymd,tod,' avg dt = ',cktime_acc(8)/cktime_cnt(8),' dt = ',cktime
               endif
            endif
         enddo
      end if

      !----------------------------------------------------------
      ! Run Ocn Model HERE, if ocean_tight_coupling
      !----------------------------------------------------------

      if (ocean_tight_coupling) then
      if (ocn_present .and. ocnrun_alarm) then
         do eoi = 1,num_inst_ocn
            if (iamin_OCNID(eoi)) then
               if (run_barriers) then
                  call t_drvstartf ('DRIVER_OCN_RUN_BARRIER')
                  call mpi_barrier(mpicom_OCNID(eoi),ierr)
                  call t_drvstopf ('DRIVER_OCN_RUN_BARRIER')
                  time_brun = mpi_wtime()
               endif
               call t_drvstartf ('DRIVER_OCN_RUN',barrier=mpicom_OCNID(eoi))
               if (drv_threading) call seq_comm_setnthreads(nthreads_OCNID)
               if (ocn_prognostic) call mct_avect_vecmult(x2o_oo(eoi),areacor_oo(eoi)%drv2mdl,seq_flds_x2o_fluxes)
               call ocn_run_mct( EClock_o, cdata_oo(eoi), x2o_oo(eoi), o2x_oo(eoi))
               call mct_avect_vecmult(o2x_oo(eoi),areacor_oo(eoi)%mdl2drv,seq_flds_o2x_fluxes)
               if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
               call t_drvstopf  ('DRIVER_OCN_RUN')
               if (run_barriers) then
                  time_erun = mpi_wtime()
                  cktime = time_erun-time_brun
                  cktime_acc(5) = cktime_acc(5) + cktime
                  cktime_cnt(5) = cktime_cnt(5) + 1
                  write(logunit,107) ' rstamp ocn_run_time: model date = ', &
                     ymd,tod,' avg dt = ',cktime_acc(5)/cktime_cnt(5),' dt = ',cktime
               endif
            endif
         enddo
      endif
      endif
 
      !----------------------------------------------------------
      ! ocn -> cpl, tight coupling (sequential type mode)
      !----------------------------------------------------------

      if (ocean_tight_coupling) then
      if (iamin_CPLALLOCNID) then
      if (ocn_present .and. ocnnext_alarm) then
         if (run_barriers) then
            call t_drvstartf ('DRIVER_O2CT_BARRIER')
            call mpi_barrier(mpicom_CPLALLOCNID,ierr)
            call t_drvstopf ('DRIVER_O2CT_BARRIER')
         endif
         call t_drvstartf ('DRIVER_O2CT',cplcom=.true.,barrier=mpicom_CPLALLOCNID)
         do eoi = 1,num_inst_ocn
            if (iamin_CPLOCNID(eoi)) then
               call t_drvstartf ('driver_o2ct_ocno2ocnx',barrier=mpicom_CPLOCNID(eoi))
               call seq_map_map(mapper_Co2x(eoi), o2x_oo(eoi), o2x_ox(eoi), msgtag=CPLOCNID(eoi)*100+eoi*10+4)
               call t_drvstopf  ('driver_o2ct_ocno2ocnx')
            endif
         enddo
         if (iamin_CPLOCNID(ens1)) then
            call t_drvstartf ('driver_o2ct_infoexch',barrier=mpicom_CPLOCNID(ens1))
            call seq_infodata_exchange(infodata,CPLOCNID(ens1),'ocn2cpl_run')
            call t_drvstopf  ('driver_o2ct_infoexch')
         endif
         call t_drvstopf  ('DRIVER_O2CT',cplcom=.true.)
         if (iamin_CPLID) then
            if (run_barriers) then
               call t_drvstartf ('DRIVER_OCNPOSTT_BARRIER')
               call mpi_barrier(mpicom_CPLID,ierr)
               call t_drvstopf ('DRIVER_OCNPOSTT_BARRIER')
            endif
            call t_drvstartf  ('DRIVER_OCNPOSTT',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
            if (info_debug > 1) then
               call t_drvstartf ('driver_ocnpostt_diagav',barrier=mpicom_CPLID)
               do eoi = 1,num_inst_ocn
                  call seq_diag_avect_mct(cdata_ox,o2x_ox(eoi),'recv ocn'//trim(ocn_suffix(eoi)))
               enddo
               call t_drvstopf  ('driver_ocnpostt_diagav')
            endif
            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('DRIVER_OCNPOSTT',cplrun=.true.)
         endif
      endif
      endif
      endif

      !----------------------------------------------------------
      ! OCN PREP
      !----------------------------------------------------------

      if (ocn_present .and. iamin_CPLID) then
         if (run_barriers) then
            call t_drvstartf ('DRIVER_ATMOCNP_BARRIER')
            call mpi_barrier(mpicom_CPLID,ierr)
            call t_drvstopf ('DRIVER_ATMOCNP_BARRIER')
         endif
         call t_drvstartf ('DRIVER_ATMOCNP',cplrun=.true.,barrier=mpicom_CPLID)
         if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
         if (ocn_prognostic) then

            ! Map ice to ocn
            if (ice_present) then
               call t_drvstartf ('driver_atmocnp_ice2ocn',barrier=mpicom_CPLID)
               do eii = 1,num_inst_ice
                  call seq_map_map(mapper_SFi2o, i2x_ix(eii), i2x_ox(eii), norm=.true.)
               enddo
               call t_drvstopf  ('driver_atmocnp_ice2ocn')
            endif

            ! Map wav to ocn
            if (wav_present) then
               call t_drvstartf ('driver_atmocnp_wav2ocn',barrier=mpicom_CPLID)
               do ewi = 1,num_inst_wav
                  call seq_map_map(mapper_Sw2o, w2x_wx(ewi), w2x_ox(ewi), norm=.true.)
               enddo
               call t_drvstopf  ('driver_atmocnp_wav2ocn')
            endif

            ! Merge ocn inputs
            call t_drvstartf ('driver_atmocnp_mrgx2o',barrier=mpicom_CPLID)
            ! Use fortran mod to address ensembles in merge
            do eoi = 1,num_inst_ocn
               eai = mod((eoi-1),num_inst_atm) + 1
               eii = mod((eoi-1),num_inst_ice) + 1
               ewi = mod((eoi-1),num_inst_wav) + 1
               exi = mod((eoi-1),num_inst_xao) + 1
               efi = mod((eoi-1),num_inst_frc) + 1
               call mrg_x2o_run_mct( cdata_ox, a2x_ox(eai), i2x_ox(eii), w2x_ox(ewi), &
                    xao_ox(exi), fractions_ox(efi), x2o_ox(eoi) )
            enddo
            call t_drvstopf  ('driver_atmocnp_mrgx2o')

            ! Accumulate ocn inputs
            ! Form partial sum of tavg ocn inputs (virtual "send" to ocn) 
            call t_drvstartf ('driver_atmocnp_accum',barrier=mpicom_CPLID)
            do eoi = 1,num_inst_ocn
!     !        call mct_aVect_accumulate(x2o_ox(eoi), x2oacc_ox(eoi))
               if (x2oacc_ox_cnt == 0) then
                  x2oacc_ox(eoi)%rAttr = x2o_ox(eoi)%rAttr
               else
                  x2oacc_ox(eoi)%rAttr = x2oacc_ox(eoi)%rAttr + x2o_ox(eoi)%rAttr
               endif
            enddo
            x2oacc_ox_cnt = x2oacc_ox_cnt + 1
            call t_drvstopf  ('driver_atmocnp_accum')
         endif
 
         ! Compute atm/ocn fluxes (virtual "recv" from ocn)
         do exi = 1,num_inst_xao
            eai = mod((exi-1),num_inst_atm) + 1
            eoi = mod((exi-1),num_inst_ocn) + 1
            efi = mod((exi-1),num_inst_frc) + 1
            if (trim(aoflux_grid) == 'ocn') then
               call t_drvstartf ('driver_atmocnp_fluxo',barrier=mpicom_CPLID)
               call seq_flux_atmocn_mct( cdata_ox, a2x_ox(eai), o2x_ox(eoi), xao_ox(exi))
               call seq_flux_ocnalb_mct(cdata_ox,xao_ox(exi),fractions_ox(efi))
               call t_drvstopf  ('driver_atmocnp_fluxo')
            else if (trim(aoflux_grid) == 'atm') then
               !--- compute later ---
            else if (trim(aoflux_grid) == 'exch') then
               call t_drvstartf ('driver_atmocnp_fluxe',barrier=mpicom_CPLID)
               call seq_flux_atmocnexch_mct( cdata_ax, cdata_ox, a2x_ax(eai), o2x_ox(eoi), xao_ax(exi), xao_ox(exi), &
                                          fractions_ax(efi), fractions_ox(efi))
               call seq_flux_ocnalb_mct(cdata_ox,xao_ox(exi),fractions_ox(efi))
               call t_drvstopf  ('driver_atmocnp_fluxe')
            endif  ! aoflux_grid
         enddo
         
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
         call t_drvstopf  ('DRIVER_ATMOCNP',cplrun=.true.)
      endif

      !----------------------------------------------------------
      ! lnd -> cpl
      !----------------------------------------------------------

      if ((lnd_present.or.sno_present) .and. lndrun_alarm) then
      if (iamin_CPLALLLNDID) then

         do eli = 1,num_inst_lnd
            if (run_barriers) then
               call t_drvstartf ('DRIVER_L2C_BARRIER')
               call mpi_barrier(mpicom_CPLALLLNDID,ierr)
               call t_drvstopf ('DRIVER_L2C_BARRIER')
            endif
            call t_drvstartf ('DRIVER_L2C',cplcom=.true.,barrier=mpicom_CPLALLLNDID)
            if (iamin_CPLLNDID(eli)) then
               if (lnd_present) then
                  call t_drvstartf ('driver_l2c_lndl2lndx',barrier=mpicom_CPLLNDID(eli))
                  call seq_map_map(mapper_Cl2x(eli), l2x_ll(eli), l2x_lx(eli), msgtag=CPLLNDID(eli)*100+eli*10+6)
                  call t_drvstopf  ('driver_l2c_lndl2lndx')
               endif

               if (sno_present .and. glc_prognostic .and. glcrun_alarm) then
                  call t_drvstartf ('driver_l2c_snos2snox',barrier=mpicom_CPLLNDID(eli))
                  call seq_map_map(mapper_Cs2x(eli), s2x_ss(eli), s2x_sx(eli), msgtag=CPLLNDID(eli)*100+eli*10+8)
                  call t_drvstopf  ('driver_l2c_snos2snox')
               endif

            endif
            if (eli == 1 .and. iamin_CPLLNDID(ens1)) then
               call t_drvstartf ('driver_l2c_infoexch', barrier=mpicom_CPLLNDID(ens1))
               call seq_infodata_exchange(infodata,CPLLNDID(ens1),'lnd2cpl_run')
               call t_drvstopf  ('driver_l2c_infoexch')
            endif
            call t_drvstopf  ('DRIVER_L2C',cplcom=.true.)

            if (iamin_CPLID) then
               if (run_barriers) then
                  call t_drvstartf ('DRIVER_LNDPOST_BARRIER')
                  call mpi_barrier(mpicom_CPLID,ierr)
                  call t_drvstopf ('DRIVER_LNDPOST_BARRIER')
               endif
               call t_drvstartf  ('DRIVER_LNDPOST',cplrun=.true.,barrier=mpicom_CPLID)
               if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
               if (info_debug > 1) then
                  call t_drvstartf ('driver_lndpost_diagav',barrier=mpicom_CPLID)

                  if (lnd_present) then
                     call seq_diag_avect_mct(cdata_lx, l2x_lx(eli),'recv lnd'//trim(lnd_suffix(eli)))
                  endif

                  if (sno_present .and. glc_prognostic .and. glcrun_alarm) then
                     call seq_diag_avect_mct(cdata_sx, s2x_sx(eli),'recv sno'//trim(lnd_suffix(eli)))
                  endif

                  call t_drvstopf  ('driver_lndpost_diagav')
               endif

               if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
               call t_drvstopf  ('DRIVER_LNDPOST',cplrun=.true.)
            endif
         enddo   ! eli
      endif   ! CPLALLLNDID

      if (iamin_CPLID .and. lnd_present) then
         do eli = 1,num_inst_lnd
            ! want to do mct_avect_copy but need to accumulate so do it manually
            if (x2racc_lx_cnt == 0) then
               x2racc_lx(eli)%rAttr(index_x2r_Flrl_rofliq,:) = l2x_lx(eli)%rAttr(index_l2x_Flrl_rofliq,:)
               x2racc_lx(eli)%rAttr(index_x2r_Flrl_rofice,:) = l2x_lx(eli)%rAttr(index_l2x_Flrl_rofice,:)
            else
               x2racc_lx(eli)%rAttr(index_x2r_Flrl_rofliq,:) = x2racc_lx(eli)%rAttr(index_x2r_Flrl_rofliq,:) + &
                                                               l2x_lx(eli)%rAttr(index_l2x_Flrl_rofliq,:)
               x2racc_lx(eli)%rAttr(index_x2r_Flrl_rofice,:) = x2racc_lx(eli)%rAttr(index_x2r_Flrl_rofice,:) + &
                                                               l2x_lx(eli)%rAttr(index_l2x_Flrl_rofice,:)
            endif
         end do
         x2racc_lx_cnt = x2racc_lx_cnt + 1
      end if

      endif   ! run alarm, lnd_present

      !----------------------------------------------------------
      ! GLC SETUP
      !----------------------------------------------------------

      if (sno_present .and. glcrun_alarm) then
         if (iamin_CPLID .and. glc_prognostic) then
            if (run_barriers) then
               call t_drvstartf ('DRIVER_GLCPREP_BARRIER')
               call mpi_barrier(mpicom_CPLID,ierr)
               call t_drvstopf ('DRIVER_GLCPREP_BARRIER')
            endif
            call t_drvstartf ('DRIVER_GLCPREP',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)

            ! Map sno to glc
            call t_drvstartf ('driver_glcprep_sno2glc',barrier=mpicom_CPLID)
            do eli = 1,num_inst_lnd
               call seq_map_map(mapper_SFs2g, s2x_sx(eli), s2x_gx(eli), norm=.true.)
            enddo
            call t_drvstopf  ('driver_glcprep_sno2glc')

            ! Merge glc inputs
            call t_drvstartf ('driver_glcprep_mrgx2g',barrier=mpicom_CPLID)
            ! Use fortran mod to address ensembles in merge
            do egi = 1,num_inst_glc
               eli = mod((egi-1),num_inst_lnd) + 1
               call mrg_x2g_run_mct( cdata_gx, s2x_gx(eli), x2g_gx(egi))
            enddo
            call t_drvstopf  ('driver_glcprep_mrgx2g')

            if (info_debug > 1) then
               call t_drvstartf ('driver_glcprep_diagav',barrier=mpicom_CPLID)
               do egi = 1,num_inst_glc
                  call seq_diag_avect_mct(cdata_gx,x2g_gx(egi),'send glc'//trim(glc_suffix(egi)))
               enddo
               call t_drvstopf  ('driver_glcprep_diagav')
            endif

            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('DRIVER_GLCPREP',cplrun=.true.)
         endif

         !----------------------------------------------------
         ! cpl -> glc
         !----------------------------------------------------

         if (iamin_CPLALLGLCID .and. glc_prognostic) then
            if (run_barriers) then
               call t_drvstartf ('DRIVER_C2G_BARRIER')
               call mpi_barrier(mpicom_CPLALLGLCID,ierr)
               call t_drvstopf ('DRIVER_C2G_BARRIER')
            endif
            call t_drvstartf ('DRIVER_C2G',cplcom=.true.,barrier=mpicom_CPLALLGLCID)
            do egi = 1,num_inst_glc
               if (iamin_CPLGLCID(egi)) then
                  call t_drvstartf ('driver_c2g_glcx2glcg',barrier=mpicom_CPLGLCID(egi))
                  call seq_map_map(mapper_Cx2g(egi), x2g_gx(egi), x2g_gg(egi), msgtag=CPLGLCID(egi)*100+egi*10+2)
                  call t_drvstopf  ('driver_c2g_glcx2glcg')
               endif
            enddo
            call t_drvstartf ('driver_c2g_infoexch',barrier=mpicom_CPLALLGLCID)
            call seq_infodata_exchange(infodata,CPLALLGLCID,'cpl2glc_run')
            call t_drvstopf  ('driver_c2g_infoexch')
            call t_drvstopf  ('DRIVER_C2G',cplcom=.true.)
         endif
      endif

      !----------------------------------------------------------
      ! rof -> cpl
      !----------------------------------------------------------

      if (rof_present .and. rofrun_alarm) then
      if (iamin_CPLALLROFID) then
         if (run_barriers) then
            call t_drvstartf ('DRIVER_R2C_BARRIER')
            call mpi_barrier(mpicom_CPLALLROFID,ierr)
            call t_drvstopf ('DRIVER_R2C_BARRIER')
         endif
         call t_drvstartf ('DRIVER_R2C',cplcom=.true., barrier=mpicom_CPLALLROFID)
         do eri = 1,num_inst_rof
            if (iamin_CPLROFID(eri)) then
               call t_drvstartf ('driver_r2c_rofr2rofx',barrier=mpicom_CPLROFID(eri))
               call seq_map_map(mapper_Cr2x(eri), r2x_rr(eri), r2x_rx(eri), msgtag=CPLROFID(eri)*100+eri*10+4)
               call t_drvstopf ('driver_r2c_rofr2rofx')
            end if
         enddo
         if (iamin_CPLROFID(ens1)) then
            call t_drvstartf ('driver_r2c_infoexch',barrier=mpicom_CPLROFID(ens1))
            call seq_infodata_exchange(infodata,CPLROFID(ens1),'rof2cpl_run')
            call t_drvstopf  ('driver_r2c_infoexch')
         endif
         call t_drvstopf  ('DRIVER_R2C',cplcom=.true.)

         if (iamin_CPLID) then
            if (run_barriers) then
               call t_drvstartf ('DRIVER_ROFPOST_BARRIER')
               call mpi_barrier(mpicom_CPLID,ierr)
               call t_drvstopf ('DRIVER_ROFPOST_BARRIER')
            endif
            call t_drvstartf  ('DRIVER_ROFPOST',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
            if (info_debug > 1) then
               call t_drvstartf ('driver_rofpost_diagav',barrier=mpicom_CPLID)
               do eri = 1,num_inst_rof
                  call seq_diag_avect_mct(cdata_rx, r2x_rx(eri),'recv roff'//trim(rof_suffix(eri)))
               enddo
               call t_drvstopf ('driver_rofpost_diagav')
            endif

            if (ocnrof_prognostic) then
               call t_drvstartf ('driver_rofpost_raccum',barrier=mpicom_CPLID)
               ! better to flux correct here if flux_epbalfact varies
               ! over the accumulation period
               call seq_infodata_GetData(infodata, flux_epbalfact = flux_epbalfact)
               ! copy and accumulate only r2o fields 
               do eri = 1,num_inst_rof
                  if (r2xacc_rx_cnt == 0) then
                     r2xacc_rx(eri)%rAttr = r2x_rx(eri)%rAttr * flux_epbalfact
                  else
                     r2xacc_rx(eri)%rAttr = r2xacc_rx(eri)%rAttr + r2x_rx(eri)%rAttr * flux_epbalfact
                  endif
               enddo
               r2xacc_rx_cnt = r2xacc_rx_cnt + 1
               call t_drvstopf ('driver_rofpost_raccum')
            endif
               call t_drvstopf  ('DRIVER_ROFPOST', cplrun=.true.)
         endif ! CPLID
      endif  ! CPLALLROFID
      endif  ! (rof_present .and. rofrun_alarm)
      
      !----------------------------------------------------------
      ! budget with old fractions
      !----------------------------------------------------------

      ! WJS (2-17-11): I am just using the first instance for the budgets because we
      ! don't expect budgets to be conserved for our case (I case). Also note that we
      ! don't expect budgets to be conserved for the interactive ensemble use case either.
      ! tcraig (aug 2012): put this after rof->cpl so the budget sees the new r2x_rx.
      ! it will also use the current r2x_ox here which is the value from the last timestep
      ! consistent with the ocean coupling

      if (iamin_CPLID .and. do_budgets) then
         if (run_barriers) then
            call t_drvstartf ('DRIVER_BUDGET1_BARRIER')
            call mpi_barrier(mpicom_CPLID,ierr)
            call t_drvstopf ('DRIVER_BUDGET1_BARRIER')
         endif
         call t_drvstartf ('DRIVER_BUDGET1',cplrun=.true.,budget=.true.,barrier=mpicom_CPLID)
         if (lnd_present) then
            call seq_diag_lnd_mct(dom_lx, fractions_lx(ens1), l2x_l=l2x_lx(ens1), x2l_l=x2l_lx(ens1))
         endif
         if (rof_present) then
            call seq_diag_rtm_mct(dom_rx, fractions_rx(ens1), r2x_r=r2x_rx(ens1), x2r_r=x2r_rx(ens1))
         endif
         if (ocn_present) then
            call seq_diag_ocn_mct(dom_ox, fractions_ox(ens1), o2x_o=o2x_ox(ens1), x2o_o=x2o_ox(ens1), &
                 xao_o=xao_ox(ens1), r2x_o=r2x_ox(ens1))
         endif
         if (ice_present) then
            call seq_diag_ice_mct(dom_ix, fractions_ix(ens1), x2i_i=x2i_ix(ens1))
         endif
         call t_drvstopf  ('DRIVER_BUDGET1',cplrun=.true.,budget=.true.)
      endif

      !----------------------------------------------------------
      ! ice -> cpl
      !----------------------------------------------------------

      if (ice_present .and. icerun_alarm) then
      if (iamin_CPLALLICEID) then
         if (run_barriers) then
            call t_drvstartf ('DRIVER_I2C_BARRIER')
            call mpi_barrier(mpicom_CPLALLICEID,ierr)
            call t_drvstopf ('DRIVER_I2C_BARRIER')
         endif
         call t_drvstartf ('DRIVER_I2C',cplcom=.true.,barrier=mpicom_CPLALLICEID)
         do eii = 1,num_inst_ice
            if (iamin_CPLICEID(eii)) then
               call t_drvstartf ('driver_i2c_icei2icex',barrier=mpicom_CPLICEID(eii))
               call seq_map_map(mapper_Ci2x(eii), i2x_ii(eii), i2x_ix(eii), msgtag=CPLICEID(eii)*100+eii*10+4)
               call t_drvstopf  ('driver_i2c_icei2icex')
            endif
         enddo
         if (iamin_CPLICEID(ens1)) then
            call t_drvstartf ('driver_i2c_infoexch',barrier=mpicom_CPLICEID(ens1))
            call seq_infodata_exchange(infodata,CPLICEID(ens1),'ice2cpl_run')
            call t_drvstopf  ('driver_i2c_infoexch')
         endif
         call t_drvstopf  ('DRIVER_I2C',cplcom=.true.)

         if (iamin_CPLID) then
            if (run_barriers) then
               call t_drvstartf ('DRIVER_ICEPOST_BARRIER')
               call mpi_barrier(mpicom_CPLID,ierr)
               call t_drvstopf ('DRIVER_ICEPOST_BARRIER')
            endif
            call t_drvstartf  ('DRIVER_ICEPOST',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
            if (info_debug > 1) then
               call t_drvstartf ('driver_icepost_diagav',barrier=mpicom_CPLID)
               do eii = 1,num_inst_ice
                  call seq_diag_avect_mct(cdata_ix,i2x_ix(eii),'recv ice'//trim(ice_suffix(eii)))
               enddo
               call t_drvstopf  ('driver_icepost_diagav')
            endif
            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('DRIVER_ICEPOST',cplrun=.true.)
         endif
      endif   ! CPLALLICEID
      endif   ! run alarm, ice_present

      !----------------------------------------------------------
      ! Update fractions based on new ice fractions
      !----------------------------------------------------------

      if (iamin_CPLID) then
         if (run_barriers) then
            call t_drvstartf ('DRIVER_FRACSET_BARRIER')
            call mpi_barrier(mpicom_CPLID,ierr)
            call t_drvstopf ('DRIVER_FRACSET_BARRIER')
         endif
         call t_drvstartf ('DRIVER_FRACSET',cplrun=.true.,barrier=mpicom_CPLID)
         if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
         call t_drvstartf ('driver_fracset_fracset',barrier=mpicom_CPLID)
         do efi = 1,num_inst_frc
            eii = mod((efi-1),num_inst_ice) + 1
            call seq_frac_set(i2x_ix(eii), &
                        cdata_ax, cdata_ix, cdata_lx, cdata_ox, cdata_gx, cdata_rx, cdata_wx, &
                        ice_present, ocn_present, lnd_present, glc_present, rof_present, wav_present, &
                        fractions_ax(efi), fractions_ix(efi), fractions_lx(efi), &
                        fractions_ox(efi), fractions_gx(efi), fractions_rx(efi), &
                        fractions_wx(efi), &
                        mapper_Fi2a, mapper_SFi2o)
         enddo
         call t_drvstopf  ('driver_fracset_fracset')
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
         call t_drvstopf  ('DRIVER_FRACSET',cplrun=.true.)
      endif

      !----------------------------------------------------------
      ! ATM/OCN FLUX CALC on atm grid with NEW FRACTIONS
      !----------------------------------------------------------

      if (ocn_present .and. iamin_CPLID) then
         ! Compute atm/ocn fluxes (virtual "recv" from ocn)
         if (trim(aoflux_grid) == 'atm') then
            if (run_barriers) then
               call t_drvstartf ('DRIVER_ATMOCNQ_BARRIER')
               call mpi_barrier(mpicom_CPLID,ierr)
               call t_drvstopf ('DRIVER_ATMOCNQ_BARRIER')
            endif
            call t_drvstartf ('DRIVER_ATMOCNQ',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
            do eoi = 1,num_inst_ocn
               efi = mod((eoi-1),num_inst_frc) + 1
               call t_drvstartf ('driver_atmocnq_ocn2atm1',barrier=mpicom_CPLID)
               call seq_map_map(mapper_So2a,o2x_ox(eoi),o2x_ax(eoi),fldlist=seq_flds_o2x_states,norm=.true., &
                                avwts_s=fractions_ox(efi),avwtsfld_s='ofrac')
               call t_drvstopf ('driver_atmocnq_ocn2atm1')
               call t_drvstartf ('driver_atmocnq_ocn2atm2',barrier=mpicom_CPLID)
               call seq_map_map(mapper_Fo2a,o2x_ox(eoi),o2x_ax(eoi),fldlist=seq_flds_o2x_fluxes,norm=.true.)
               call t_drvstopf ('driver_atmocnq_ocn2atm2')
            enddo
            call t_drvstartf ('driver_atmocnq_fluxa',barrier=mpicom_CPLID)
            do exi = 1,num_inst_xao 
               eai = mod((exi-1),num_inst_atm) + 1
               eoi = mod((exi-1),num_inst_ocn) + 1
               efi = mod((exi-1),num_inst_frc) + 1
               call seq_flux_atmocn_mct( cdata_ax, a2x_ax(eai), o2x_ax(eoi), xao_ax(exi))
               call seq_flux_ocnalb_mct(cdata_ox,xao_ox(exi),fractions_ox(efi))
!               call seq_hist_spewav('a2x_ax(eai)',gsmap_ax,a2x_ax(eai),atm_nx,atm_ny,1)
!               call seq_hist_spewav('o2x_ax(eoi)',gsmap_ax,o2x_ax(eoi),atm_nx,atm_ny,1)
!               call seq_hist_spewav('xao_ax',gsmap_ax,xao_ax,atm_nx,atm_ny,1)
            enddo
            call t_drvstopf  ('driver_atmocnq_fluxa')

            call t_drvstartf ('driver_atmocnq_atm2ocnf',barrier=mpicom_CPLID)
! this mapping has to be done with area overlap mapping for all fields 
! due to the masking of the xao_ax data and the fact that a2oS is bilinear
            do exi = 1,num_inst_xao 
               call seq_map_map(mapper_Fa2o,xao_ax(exi),xao_ox(exi),norm=.true.)
            enddo
            call t_drvstopf  ('driver_atmocnq_atm2ocnf')
            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('DRIVER_ATMOCNQ',cplrun=.true.)
         endif
      endif

      !----------------------------------------------------------
      ! ATM SETUP
      !----------------------------------------------------------

      if (atm_present .and. atmrun_alarm) then
 
         !----------------------------------------------------------
         ! atm prep
         !----------------------------------------------------------

         if (iamin_CPLID .and. atm_prognostic) then
            if (run_barriers) then
               call t_drvstartf ('DRIVER_ATMPREP_BARRIER')
               call mpi_barrier(mpicom_CPLID,ierr)
               call t_drvstopf ('DRIVER_ATMPREP_BARRIER')
            endif
            call t_drvstartf ('DRIVER_ATMPREP',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
            if (ocn_present) then
               if (trim(aoflux_grid) == 'ocn') then
                  call t_drvstartf ('driver_atmprep_ocn2atmf',barrier=mpicom_CPLID)
                  do exi = 1,num_inst_xao
                     efi = mod((exi-1),num_inst_frc) + 1
                     call seq_map_map(mapper_So2a,xao_ox(exi),xao_ax(exi),fldlist=seq_flds_xao_states,norm=.true., &
                                      avwts_s=fractions_ox(efi),avwtsfld_s='ofrac')
                     call seq_map_map(mapper_Fo2a,xao_ox(exi),xao_ax(exi),fldlist=seq_flds_xao_fluxes,norm=.true., &
                                      avwts_s=fractions_ox(efi),avwtsfld_s='ofrac')
                  enddo
                  call t_drvstopf  ('driver_atmprep_ocn2atmf')
               endif
            endif  ! ocn_present
            if (ocn_present) then
               if (trim(aoflux_grid) /= 'atm') then
                  do eoi = 1,num_inst_ocn
                     efi = mod((eoi-1),num_inst_frc) + 1
                     call t_drvstartf ('driver_atmprep_ocn2atm1',barrier=mpicom_CPLID)
                     call seq_map_map(mapper_So2a,o2x_ox(eoi),o2x_ax(eoi),fldlist=seq_flds_o2x_states,norm=.true., &
                                      avwts_s=fractions_ox(efi),avwtsfld_s='ofrac')
                     call t_drvstopf  ('driver_atmprep_ocn2atm1')
                     call t_drvstartf ('driver_atmprep_ocn2atm2',barrier=mpicom_CPLID)
                     call seq_map_map(mapper_Fo2a,o2x_ox(eoi),o2x_ax(eoi),fldlist=seq_flds_o2x_fluxes,norm=.true.)
                     call t_drvstopf  ('driver_atmprep_ocn2atm2')
                  enddo
               endif
               call t_drvstartf ('driver_atmprep_ocn2atmb',barrier=mpicom_CPLID)
               do exi = 1,num_inst_xao
                  efi = mod((exi-1),num_inst_frc) + 1
                  call seq_map_map(mapper_So2a,xao_ox(exi),xao_ax(exi),fldlist=seq_flds_xao_albedo,norm=.true., &
                                   avwts_s=fractions_ox(efi),avwtsfld_s='ofrac')
               enddo
               call t_drvstopf  ('driver_atmprep_ocn2atmb')
            endif
            if (ice_present) then
               call t_drvstartf ('driver_atmprep_ice2atm',barrier=mpicom_CPLID)
               do eii = 1,num_inst_ice
                  efi = mod((eii-1),num_inst_frc) + 1
                  call seq_map_map(mapper_Si2a, i2x_ix(eii), i2x_ax(eii), fldlist=seq_flds_i2x_states, &
                       avwts_s=fractions_ix(efi), avwtsfld_s='ifrac')
                  call seq_map_map(mapper_Fi2a, i2x_ix(eii), i2x_ax(eii), fldlist=seq_flds_i2x_fluxes, &
                       avwts_s=fractions_ix(efi), avwtsfld_s='ifrac')
               enddo
               call t_drvstopf  ('driver_atmprep_ice2atm')
            endif
            if (lnd_present) then
               call t_drvstartf ('driver_atmprep_lnd2atm',barrier=mpicom_CPLID)
               do eli = 1,num_inst_lnd
                  efi = mod((eli-1),num_inst_frc) + 1
                  call seq_map_map(mapper_Sl2a,l2x_lx(eli),l2x_ax(eli), &
                                   fldlist=seq_flds_l2x_states, norm=.true., &
                                   avwts_s=fractions_lx(efi),avwtsfld_s='lfrin')
                  call seq_map_map(mapper_Fl2a,l2x_lx(eli),l2x_ax(eli), &
                                   fldlist=seq_flds_l2x_fluxes, norm=.true., &
                                   avwts_s=fractions_lx(efi),avwtsfld_s='lfrin')
               enddo
               call t_drvstopf  ('driver_atmprep_lnd2atm')
            endif
            call t_drvstartf ('driver_atmprep_mrgx2a',barrier=mpicom_CPLID)

            ! Use fortran mod to address ensembles in merge
            do eai = 1,num_inst_atm
               eli = mod((eai-1),num_inst_lnd) + 1
               eoi = mod((eai-1),num_inst_ocn) + 1
               eii = mod((eai-1),num_inst_ice) + 1
               exi = mod((eai-1),num_inst_xao) + 1
               efi = mod((eai-1),num_inst_frc) + 1
               call mrg_x2a_run_mct( cdata_ax, l2x_ax(eli), o2x_ax(eoi), xao_ax(exi), i2x_ax(eii), fractions_ax(efi), x2a_ax(eai))
            enddo
!            call seq_hist_spewav('x2a_ax(eai)',gsmap_ax,x2a_ax(eai),atm_nx,atm_ny,1)
            call t_drvstopf  ('driver_atmprep_mrgx2a')

            if (info_debug > 1) then
               call t_drvstartf ('driver_atmprep_diagav',barrier=mpicom_CPLID)
               do eai = 1,num_inst_atm
                  call seq_diag_avect_mct(cdata_ax,x2a_ax(eai),'send atm'//trim(atm_suffix(eai)))
               enddo
               call t_drvstopf  ('driver_atmprep_diagav')
            endif
            call t_drvstopf  ('DRIVER_ATMPREP',cplrun=.true.)
            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
         endif  ! CPLID

         !----------------------------------------------------------
         ! cpl -> atm
         !----------------------------------------------------------

         if (iamin_CPLALLATMID) then
         if (atm_prognostic) then
            if (run_barriers) then
               call t_drvstartf ('DRIVER_C2A_BARRIER')
               call mpi_barrier(mpicom_CPLALLATMID,ierr)
               call t_drvstopf ('DRIVER_C2A_BARRIER')
            endif
            call t_drvstartf ('DRIVER_C2A',cplcom=.true.,barrier=mpicom_CPLALLATMID)
            do eai = 1,num_inst_atm
               if (iamin_CPLATMID(eai)) then
                  call t_drvstartf ('driver_c2a_atmx2atma',barrier=mpicom_CPLATMID(eai))
                  call seq_map_map(mapper_Cx2a(eai), x2a_ax(eai), x2a_aa(eai), msgtag=CPLATMID(eai)*100+eai*10+2)
                  call t_drvstopf  ('driver_c2a_atmx2atma')
               endif
            enddo
            call t_drvstartf ('driver_c2a_infoexch',barrier=mpicom_CPLALLATMID)
            call seq_infodata_exchange(infodata,CPLALLATMID,'cpl2atm_run')
            call t_drvstopf  ('driver_c2a_infoexch')
            call t_drvstopf  ('DRIVER_C2A',cplcom=.true.)
         endif
         endif

      endif

      !----------------------------------------------------------
      ! Run Ocn Model HERE if NOT ocean_tight_coupling
      !----------------------------------------------------------

      if (.not.ocean_tight_coupling) then
      if (ocn_present .and. ocnrun_alarm) then
         do eoi = 1,num_inst_ocn
            if (iamin_OCNID(eoi)) then
               if (run_barriers) then
                  call t_drvstartf ('DRIVER_OCN_RUN_BARRIER')
                  call mpi_barrier(mpicom_OCNID(eoi),ierr)
                  call t_drvstopf ('DRIVER_OCN_RUN_BARRIER')
                  time_brun = mpi_wtime()
               endif
               call t_drvstartf ('DRIVER_OCN_RUN',barrier=mpicom_OCNID(eoi))
               if (drv_threading) call seq_comm_setnthreads(nthreads_OCNID)
               if (ocn_prognostic) call mct_avect_vecmult(x2o_oo(eoi),areacor_oo(eoi)%drv2mdl,seq_flds_x2o_fluxes)
               call ocn_run_mct( EClock_o, cdata_oo(eoi), x2o_oo(eoi), o2x_oo(eoi))
               call mct_avect_vecmult(o2x_oo(eoi),areacor_oo(eoi)%mdl2drv,seq_flds_o2x_fluxes)
               if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
               call t_drvstopf  ('DRIVER_OCN_RUN')
               if (run_barriers) then
                  time_erun = mpi_wtime()
                  cktime = time_erun-time_brun
                  cktime_acc(5) = cktime_acc(5) + cktime
                  cktime_cnt(5) = cktime_cnt(5) + 1
                  write(logunit,107) ' rstamp ocn_run_time: model date = ', &
                     ymd,tod,' avg dt = ',cktime_acc(5)/cktime_cnt(5),' dt = ',cktime
               endif
            endif
         enddo
      endif
      endif
 
      !----------------------------------------------------------
      ! RUN atm model
      !----------------------------------------------------------

      if (atm_present .and. atmrun_alarm) then
         do eai = 1,num_inst_atm
            if (iamin_ATMID(eai)) then
               if (run_barriers) then
                  call t_drvstartf ('DRIVER_ATM_RUN_BARRIER')
                  call mpi_barrier(mpicom_ATMID(eai),ierr)
                  call t_drvstopf ('DRIVER_ATM_RUN_BARRIER')
                  time_brun = mpi_wtime()
               endif
               call t_drvstartf ('DRIVER_ATM_RUN',barrier=mpicom_ATMID(eai))
               if (drv_threading) call seq_comm_setnthreads(nthreads_ATMID)
               if (atm_prognostic) call mct_avect_vecmult(x2a_aa(eai),areacor_aa(eai)%drv2mdl,seq_flds_x2a_fluxes)
               call atm_run_mct( EClock_a, cdata_aa(eai), x2a_aa(eai), a2x_aa(eai))
               call mct_avect_vecmult(a2x_aa(eai),areacor_aa(eai)%mdl2drv,seq_flds_a2x_fluxes)
               if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
               call t_drvstopf  ('DRIVER_ATM_RUN')
               if (run_barriers) then
                  time_erun = mpi_wtime()
                  cktime = time_erun-time_brun
                  cktime_acc(2) = cktime_acc(2) + cktime
                  cktime_cnt(2) = cktime_cnt(2) + 1
                  write(logunit,107) ' rstamp atm_run_time: model date = ', &
                     ymd,tod,' avg dt = ',cktime_acc(2)/cktime_cnt(2),' dt = ',cktime
               endif
            endif
         enddo
      endif

      !----------------------------------------------------------
      ! Run Glc Model
      !----------------------------------------------------------

      if (glc_present .and. glcrun_alarm) then
         do egi = 1,num_inst_glc
            if (iamin_GLCID(egi)) then
               if (run_barriers) then
                  call t_drvstartf ('DRIVER_GLC_RUN_BARRIER')
                  call mpi_barrier(mpicom_GLCID(egi),ierr)
                  call t_drvstopf ('DRIVER_GLC_RUN_BARRIER')
                  time_brun = mpi_wtime()
               endif
               call t_drvstartf ('DRIVER_GLC_RUN',barrier=mpicom_GLCID(egi))
               if (drv_threading) call seq_comm_setnthreads(nthreads_GLCID)
               if (glc_prognostic) call mct_avect_vecmult(x2g_gg(egi),areacor_gg(egi)%drv2mdl,seq_flds_x2g_fluxes)
               call glc_run_mct( EClock_g, cdata_gg(egi), x2g_gg(egi), g2x_gg(egi))
               call mct_avect_vecmult(g2x_gg(egi),areacor_gg(egi)%mdl2drv,seq_flds_g2x_fluxes)
               if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
               call t_drvstopf  ('DRIVER_GLC_RUN')
               if (run_barriers) then
                  time_erun = mpi_wtime()
                  cktime = time_erun-time_brun
                  cktime_acc(6) = cktime_acc(6) + cktime
                  cktime_cnt(6) = cktime_cnt(6) + 1
                  write(logunit,107) ' rstamp glc_run_time: model date = ', &
                     ymd,tod,' avg dt = ',cktime_acc(6)/cktime_cnt(6),' dt = ',cktime
               endif
            endif
         enddo
      endif
 
      !----------------------------------------------------------
      ! wav -> cpl
      !----------------------------------------------------------

      if (wav_present .and. wavrun_alarm) then
      if (iamin_CPLALLWAVID) then
         if (run_barriers) then
            call t_drvstartf ('DRIVER_W2C_BARRIER')
            call mpi_barrier(mpicom_CPLALLWAVID,ierr)
            call t_drvstopf ('DRIVER_W2C_BARRIER')
         endif
         call t_drvstartf ('DRIVER_W2C',cplcom=.true.,barrier=mpicom_CPLALLWAVID)
         do ewi = 1,num_inst_wav
            if (iamin_CPLWAVID(ewi)) then
               call t_drvstartf ('driver_w2c_wavw2wavx',barrier=mpicom_CPLWAVID(ewi))
               call seq_map_map(mapper_Cw2x(ewi), w2x_ww(ewi), w2x_wx(ewi))
               call t_drvstopf  ('driver_w2c_wavw2wavx')
            endif
         enddo
         if (iamin_CPLWAVID(ens1)) then
            call t_drvstartf ('driver_w2c_infoexch',barrier=mpicom_CPLWAVID(ens1))
            call seq_infodata_exchange(infodata,CPLWAVID(ens1),'wav2cpl_run')
            call t_drvstopf  ('driver_w2c_infoexch')
         endif
         call t_drvstopf  ('DRIVER_W2C',cplcom=.true.)

         if (iamin_CPLID) then
            if (run_barriers) then
               call t_drvstartf ('DRIVER_WAVPOST_BARRIER')
               call mpi_barrier(mpicom_CPLID,ierr)
               call t_drvstopf ('DRIVER_WAVPOST_BARRIER')
            endif
            call t_drvstartf  ('DRIVER_WAVPOST',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
            if (info_debug > 1) then
               call t_drvstartf ('driver_wavpost_diagav',barrier=mpicom_CPLID)
               do ewi = 1,num_inst_wav
                  call seq_diag_avect_mct(cdata_wx,w2x_wx(ewi),'recv wav'//trim(wav_suffix(ewi)))
               enddo
               call t_drvstopf  ('driver_wavpost_diagav')
            endif
            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('DRIVER_WAVPOST',cplrun=.true.)
         endif
      endif   ! CPLALLWAVID
      endif   ! run alarm, wav_present

      !----------------------------------------------------------
      ! glc -> cpl
      !----------------------------------------------------------

      if (glc_present .and. glcrun_alarm) then
      if (iamin_CPLALLGLCID) then
         if (run_barriers) then
            call t_drvstartf ('DRIVER_G2C_BARRIER')
            call mpi_barrier(mpicom_CPLALLGLCID,ierr)
            call t_drvstopf ('DRIVER_G2C_BARRIER')
         endif
         call t_drvstartf ('DRIVER_G2C',cplcom=.true.,barrier=mpicom_CPLALLGLCID)
         do egi = 1,num_inst_glc
            if (iamin_CPLGLCID(egi)) then
               call t_drvstartf ('driver_g2c_glcg2glcx',barrier=mpicom_CPLGLCID(egi))
               call seq_map_map(mapper_Cg2x(egi), g2x_gg(egi), g2x_gx(egi), msgtag=CPLGLCID(egi)*100+egi*10+4)
               call t_drvstopf  ('driver_g2c_glcg2glcx')
            endif
         enddo
         if (iamin_CPLGLCID(ens1)) then
            call t_drvstartf ('driver_g2c_infoexch',barrier=mpicom_CPLGLCID(ens1))
            call seq_infodata_exchange(infodata,CPLGLCID(ens1),'glc2cpl_run')
            call t_drvstopf  ('driver_g2c_infoexch')
         endif
         call t_drvstopf  ('DRIVER_G2C',cplcom=.true.)

         if (iamin_CPLID) then
            if (run_barriers) then
               call t_drvstartf ('DRIVER_GLCPOST_BARRIER')
               call mpi_barrier(mpicom_CPLID,ierr)
               call t_drvstopf ('DRIVER_GLCPOST_BARRIER')
            endif
            call t_drvstartf  ('DRIVER_GLCPOST',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
            if (info_debug > 1) then
               call t_drvstartf ('driver_glcpost_diagav',barrier=mpicom_CPLID)
               do egi = 1,num_inst_glc
                  call seq_diag_avect_mct(cdata_gx,g2x_gx(egi),'recv glc'//trim(glc_suffix(egi)))
               enddo
               call t_drvstopf  ('driver_glcpost_diagav')
            endif
            
            if (sno_prognostic) then
               call t_drvstartf ('driver_glcpost_glc2sno',barrier=mpicom_CPLID)
               do egi = 1,num_inst_glc
                  call seq_map_map(mapper_SFg2s, g2x_gx(egi), g2x_sx(egi), norm=.true.)
               enddo
               call t_drvstopf  ('driver_glcpost_glc2sno')

               call t_drvstartf ('driver_glcpost_mrgx2s',barrier=mpicom_CPLID)
               ! Use fortran mod to address ensembles in merge
               do eli = 1,num_inst_lnd
                  egi = mod((eli-1),num_inst_glc) + 1
                  call mrg_x2s_run_mct( cdata_sx, g2x_sx(egi), x2s_sx(eli) )
               enddo
               call t_drvstopf  ('driver_glcpost_mrgx2s')
            end if

            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('DRIVER_GLCPOST',cplrun=.true.)
         endif
      endif   ! CPLALLGLCID
      endif   ! run alarm, glc_present

      !----------------------------------------------------------
      ! atm -> cpl
      !----------------------------------------------------------

      if (atm_present .and. atmrun_alarm) then
      if (iamin_CPLALLATMID) then
         if (run_barriers) then
            call t_drvstartf ('DRIVER_A2C_BARRIER')
            call mpi_barrier(mpicom_CPLALLATMID,ierr)
            call t_drvstopf ('DRIVER_A2C_BARRIER')
         endif
         call t_drvstartf ('DRIVER_A2C',cplcom=.true.,barrier=mpicom_CPLALLATMID)
         do eai = 1,num_inst_atm
            if (iamin_CPLATMID(eai)) then
               call t_drvstartf ('driver_a2c_atma2atmx',barrier=mpicom_CPLATMID(eai))
               call seq_map_map(mapper_Ca2x(eai), a2x_aa(eai), a2x_ax(eai), msgtag=CPLATMID(eai)*100+eai*10+4)
               call t_drvstopf  ('driver_a2c_atma2atmx')
            endif
         enddo
         if (iamin_CPLATMID(ens1)) then
            call t_drvstartf ('driver_a2c_infoexch',barrier=mpicom_CPLATMID(ens1))
            call seq_infodata_exchange(infodata,CPLATMID(ens1),'atm2cpl_run')
            call t_drvstopf  ('driver_a2c_infoexch')
         endif
         call t_drvstopf  ('DRIVER_A2C',cplcom=.true.)

         if (iamin_CPLID) then
            if (run_barriers) then
               call t_drvstartf ('DRIVER_ATMPOST_BARRIER')
               call mpi_barrier(mpicom_CPLID,ierr)
               call t_drvstopf ('DRIVER_ATMPOST_BARRIER')
            endif
            call t_drvstartf ('DRIVER_ATMPOST',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
            if (info_debug > 1) then
               call t_drvstartf ('driver_atmpost_diagav',barrier=mpicom_CPLID)
               do eai = 1,num_inst_atm
                  call seq_diag_avect_mct(cdata_ax,a2x_ax(eai),'recv atm'//trim(atm_suffix(eai)))
               enddo
               call t_drvstopf  ('driver_atmpost_diagav')
            endif
            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('DRIVER_ATMPOST',cplrun=.true.)
         endif
      endif ! CPLALLATMID
      endif ! run alarm

      !----------------------------------------------------------
      ! budget with new fractions
      !----------------------------------------------------------

      if (iamin_CPLID .and. do_budgets) then
         if (run_barriers) then
            call t_drvstartf ('DRIVER_BUDGET2_BARRIER')
            call mpi_barrier(mpicom_CPLID,ierr)
            call t_drvstopf ('DRIVER_BUDGET2_BARRIER')
         endif
         call t_drvstartf ('DRIVER_BUDGET2',cplrun=.true.,budget=.true.,barrier=mpicom_CPLID)
         if (atm_present) then
            call seq_diag_atm_mct(dom_ax, fractions_ax(ens1), a2x_a=a2x_ax(ens1), x2a_a=x2a_ax(ens1))
         endif
         if (ice_present) then
            call seq_diag_ice_mct(dom_ix, fractions_ix(ens1), i2x_i=i2x_ix(ens1))
         endif
         call t_drvstopf  ('DRIVER_BUDGET2',cplrun=.true.,budget=.true.)

         call t_drvstartf ('DRIVER_BUDGET3',cplrun=.true.,budget=.true.,barrier=mpicom_CPLID)
         call seq_diag_accum_mct()
         call t_drvstopf  ('DRIVER_BUDGET3',cplrun=.true.,budget=.true.)

         call t_drvstartf ('DRIVER_BUDGETF',cplrun=.true.,budget=.true.,barrier=mpicom_CPLID)
         if (.not. dead_comps) call seq_diag_print_mct(EClock_d,stop_alarm,budget_inst, &
            budget_daily, budget_month, budget_ann, budget_ltann, budget_ltend)
         call seq_diag_zero_mct(EClock=EClock_d)
         call t_drvstopf  ('DRIVER_BUDGETF',cplrun=.true.,budget=.true.)
      endif

      !----------------------------------------------------------
      ! ocn -> cpl, loose coupling (concurrent type mode)
      !----------------------------------------------------------

      if (.not.ocean_tight_coupling) then
      if (iamin_CPLALLOCNID) then
      if (ocn_present .and. ocnnext_alarm) then
         if (run_barriers) then
            call t_drvstartf ('DRIVER_O2C_BARRIER')
            call mpi_barrier(mpicom_CPLALLOCNID,ierr)
            call t_drvstopf ('DRIVER_O2C_BARRIER')
         endif
         call t_drvstartf ('DRIVER_O2C',cplcom=.true.,barrier=mpicom_CPLALLOCNID)
         do eoi = 1,num_inst_ocn
            if (iamin_CPLOCNID(eoi)) then
               call t_drvstartf ('driver_o2c_ocno2ocnx',barrier=mpicom_CPLOCNID(eoi))
               call seq_map_map(mapper_Co2x(eoi), o2x_oo(eoi), o2x_ox(eoi), msgtag=CPLOCNID(eoi)*100+eoi*10+6)
               call t_drvstopf  ('driver_o2c_ocno2ocnx')
            endif
         enddo
         if (iamin_CPLOCNID(ens1)) then
            call t_drvstartf ('driver_o2c_infoexch',barrier=mpicom_CPLOCNID(ens1))
            call seq_infodata_exchange(infodata,CPLOCNID(ens1),'ocn2cpl_run')
            call t_drvstopf  ('driver_o2c_infoexch')
         endif
         call t_drvstopf  ('DRIVER_O2C',cplcom=.true.)

         if (iamin_CPLID) then
            if (run_barriers) then
               call t_drvstartf ('DRIVER_OCNPOST_BARRIER')
               call mpi_barrier(mpicom_CPLID,ierr)
               call t_drvstopf ('DRIVER_OCNPOST_BARRIER')
            endif
            call t_drvstartf  ('DRIVER_OCNPOST',cplrun=.true.,barrier=mpicom_CPLID)
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
            if (info_debug > 1) then
               call t_drvstartf ('driver_ocnpost_diagav',barrier=mpicom_CPLID)
               do eoi = 1,num_inst_ocn
                  call seq_diag_avect_mct(cdata_ox,o2x_ox(eoi),'recv ocn'//trim(ocn_suffix(eoi)))
               enddo
               call t_drvstopf  ('driver_ocnpost_diagav')
            endif
            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
            call t_drvstopf  ('DRIVER_OCNPOST',cplrun=.true.)
         endif
      endif
      endif
      endif

      !----------------------------------------------------------
      ! Save driver level restart information
      !----------------------------------------------------------

      if ( restart_alarm .and. iamin_CPLID) then
         if (run_barriers) then
            call t_drvstartf ('DRIVER_RESTART_BARRIER')
            call mpi_barrier(mpicom_CPLID,ierr)
            call t_drvstopf ('DRIVER_RESTART_BARRIER')
         endif
         call t_drvstartf ('DRIVER_RESTART',cplrun=.true.,barrier=mpicom_CPLID)
         if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
         if (iamroot_CPLID) then
            write(logunit,104) ' Write restart file at ',ymd,tod
            call shr_sys_flush(logunit)
         endif
         call seq_rest_write(EClock_d,seq_SyncClock)
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
         call t_drvstopf  ('DRIVER_RESTART',cplrun=.true.)
      endif

      !----------------------------------------------------------
      ! Write history file, only AVs on CPLID
      !----------------------------------------------------------

      if (iamin_CPLID) then

         if (run_barriers) then
            call t_drvstartf ('DRIVER_HISTORY_BARRIER')
            call mpi_barrier(mpicom_CPLID,ierr)
            call t_drvstopf ('DRIVER_HISTORY_BARRIER')
         endif
         call t_drvstartf ('DRIVER_HISTORY',cplrun=.true.,barrier=mpicom_CPLID)
         if ( history_alarm) then
            if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
            if (iamroot_CPLID) then
               write(logunit,104) ' Write history file at ',ymd,tod
               call shr_sys_flush(logunit)
            endif
            call seq_hist_write(EClock_d)
            if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
         endif

         if (do_histavg) then
            call seq_hist_writeavg(EClock_d,histavg_alarm)
         endif

         if (do_hist_a2x) then
            do eai = 1,num_inst_atm
               if (trim(hist_a2x_flds) == 'all') then
                  call seq_hist_writeaux(EClock_d,'a2x'//trim(atm_suffix(eai)),'doma', &
                       cdata_ax,a2x_ax(eai),atm_nx,atm_ny,ncpl)
               else
                  call seq_hist_writeaux(EClock_d,'a2x'//trim(atm_suffix(eai)),'doma', &
                       cdata_ax,a2x_ax(eai),atm_nx,atm_ny,ncpl,flds=hist_a2x_flds)
               endif
            enddo
         endif

         if (do_hist_a2x3hr) then
            do eai = 1,num_inst_atm
               if (trim(hist_a2x3hr_flds) == 'all') then
                  call seq_hist_writeaux(EClock_d,'a2x3h'//trim(atm_suffix(eai)),'doma', &
                       cdata_ax,a2x_ax(eai),atm_nx,atm_ny,8,t3hr_alarm)
               else
                  call seq_hist_writeaux(EClock_d,'a2x3h'//trim(atm_suffix(eai)),'doma', &
                       cdata_ax,a2x_ax(eai),atm_nx,atm_ny,8,t3hr_alarm,flds=hist_a2x3hr_flds)
               end if
            enddo
         endif

         if (do_hist_a2x3hrp) then
            do eai = 1,num_inst_atm
               if (trim(hist_a2x3hrp_flds) == 'all') then
                  call seq_hist_writeaux(EClock_d,'a2x3h_prec'//trim(atm_suffix(eai)),'doma', &
                       cdata_ax,a2x_ax(eai),atm_nx,atm_ny,8,t3hr_alarm)
               else
                  call seq_hist_writeaux(EClock_d,'a2x3h_prec'//trim(atm_suffix(eai)),'doma', &
                       cdata_ax,a2x_ax(eai),atm_nx,atm_ny,8,t3hr_alarm,flds=hist_a2x3hrp_flds)
               end if
            enddo
         endif

         if (do_hist_a2x24hr) then
            do eai = 1,num_inst_atm
               call seq_hist_writeaux(EClock_d,'a2x1d'//trim(atm_suffix(eai)),'doma', &
                    cdata_ax,a2x_ax(eai),atm_nx,atm_ny,1,t24hr_alarm)
            enddo
         endif

         if (do_hist_s2x1yr .and. glcrun_alarm) then
            do eli = 1,num_inst_lnd
               ! Use yr_offset=-1 so the file with fields from year 1 has time stamp 
               ! 0001-01-01 rather than 0002-01-01, etc.
               call seq_hist_writeaux(EClock_d,'s2x'//trim(lnd_suffix(eli)),'doml', &
                    cdata_sx,s2x_sx(eli),lnd_nx,lnd_ny,1,t1yr_alarm,yr_offset=-1)
            enddo
         endif

         if (do_hist_l2x) then
            do eli = 1,num_inst_lnd
               call seq_hist_writeaux(EClock_d,'l2x'//trim(lnd_suffix(eli)),'doml', &
                    cdata_lx,l2x_lx(eli),lnd_nx,lnd_ny,ncpl)
            enddo
         endif
         call t_drvstopf  ('DRIVER_HISTORY',cplrun=.true.)

      end if

      ! --- End timestep clock/timing diagnostics
      call t_drvstartf ('DRIVER_TSTAMP_WRITE',cplrun=.true.)
      if (tod == 0 .or. info_debug > 1) then
         if (iamroot_CPLID) then
            call date_and_time(dstr,tstr)
            Time_estep = mpi_wtime()
            cktime = time_estep-time_bstep
            cktime_acc(1) = cktime_acc(1) + cktime
            cktime_cnt(1) = cktime_cnt(1) + 1
            write(logunit,101) ' tStamp_write: model date = ',ymd,tod, &
               ' wall clock = ',dstr(1:4),'-',dstr(5:6),'-',dstr(7:8),' ',&
                                tstr(1:2),':',tstr(3:4),':',tstr(5:6), &
               ' avg dt = ',cktime_acc(1)/cktime_cnt(1),' dt = ',cktime 
            Time_bstep = mpi_wtime()
            call shr_sys_flush(logunit)
         endif
      endif
      if (tod == 0 .or. info_debug > 1) then
         !! Report on memory usage
         !! For now, just look at the first instance of each component
         if (iamroot_CPLID .or. iamroot_OCNID(ens1) .or. iamroot_ATMID(ens1) .or. &
             iamroot_LNDID(ens1) .or. iamroot_ICEID(ens1) .or. iamroot_GLCID(ens1) .or. &
             iamroot_WAVID(ens1)) then
            call shr_mem_getusage(msize,mrss)
            write(logunit,105) ' memory_write: model date = ',ymd,tod, &
               ' memory = ',mrss,' MB (highwater)    ',msize,' MB (usage)', &
               '  (pe=',iam_GLOID,' comps=',trim(complist)//')'
         endif
      endif
      if (info_debug > 1) then
         if (iamroot_CPLID) then
            call seq_infodata_GetData(infodata,nextsw_cday=nextsw_cday)
!            write(logunit,106) ' nextsw_cday = ',nextsw_cday
            write(logunit,*) '  nextsw_cday = ',nextsw_cday
         endif
      endif
      call t_drvstopf  ('DRIVER_TSTAMP_WRITE',cplrun=.true.)

      call t_stopf  ('DRIVER_RUN_LOOP')
      ! --- Write out performance data 
      call t_drvstartf  ('DRIVER_TPROF_WRITE',cplrun=.true.)
      if (tprof_alarm) then
         call t_startf("sync1_tprof")
         call mpi_barrier(mpicom_GLOID,ierr)
         call t_stopf("sync1_tprof")

         write(timing_file,'(a,i8.8,a1,i5.5)') trim(tchkpt_dir)//"/ccsm_timing_",ymd,"_",tod
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
      endif
      call t_drvstopf  ('DRIVER_TPROF_WRITE',cplrun=.true.)

   enddo   ! driver run loop

   call t_startf ('DRIVER_RUN_LOOP_BSTOP')
   call mpi_barrier(mpicom_GLOID,ierr)
   call t_stopf ('DRIVER_RUN_LOOP_BSTOP')

   Time_end = mpi_wtime()

   !----------------------------------------------------------
   ! Ending of basic time step loop
   !----------------------------------------------------------

end subroutine ccsm_run


!===============================================================================

subroutine ccsm_final()
  use shr_pio_mod, only : shr_pio_finalize
   implicit none
   !------------------------------------------------------------------------
   ! Finalization of all models
   !------------------------------------------------------------------------

 103  format( 5A )
   ! TODO finalize routines need to be cleaned up 

   call t_barrierf ('DRIVER_FINAL_BARRIER', mpicom_GLOID)
   call t_startf ('DRIVER_FINAL')

   call seq_timemgr_EClockGetData( EClock_d, stepno=endstep)
   call shr_mem_getusage(msize,mrss)

   do ewi = 1,num_inst_wav
      if (iamin_WAVID(ewi)) then
         if (drv_threading) call seq_comm_setnthreads(nthreads_WAVID)
         call wav_final_mct( EClock_w, cdata_ww(ewi), x2w_ww(ewi), w2x_ww(ewi))
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
      end if
   enddo

   do eai = 1,num_inst_atm
      if (iamin_ATMID(eai)) then
         if (drv_threading) call seq_comm_setnthreads(nthreads_ATMID)
         call atm_final_mct( EClock_a, cdata_aa(eai), x2a_aa(eai), a2x_aa(eai))
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
      end if
   enddo

   do eli = 1,num_inst_lnd
      if (iamin_LNDID(eli)) then
         if (drv_threading) call seq_comm_setnthreads(nthreads_LNDID)
         call lnd_final_mct( EClock_l, cdata_ll(eli), x2l_ll(eli), l2x_ll(eli), &
                                       cdata_ss(eli), x2s_ss(eli), s2x_ss(eli))
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
      end if
   enddo

   do eri = 1,num_inst_rof
      if (iamin_ROFID(eri)) then
         if (drv_threading) call seq_comm_setnthreads(nthreads_ROFID)
         call rof_final_mct( EClock_l, cdata_rr(eri), x2r_rr(eri), r2x_rr(eri))
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
      end if
   enddo

   do eii = 1,num_inst_ice
      if (iamin_ICEID(eii)) then
         if (drv_threading) call seq_comm_setnthreads(nthreads_ICEID)
         call ice_final_mct( EClock_i, cdata_ii(eii), x2i_ii(eii), i2x_ii(eii))
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
      end if
   enddo

   do eoi = 1,num_inst_ocn
      if (iamin_OCNID(eoi)) then
         if (drv_threading) call seq_comm_setnthreads(nthreads_OCNID)
         call ocn_final_mct( EClock_o, cdata_oo(eoi), x2o_oo(eoi), o2x_oo(eoi))
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
      end if
   enddo

   do egi = 1,num_inst_glc
      if (iamin_GLCID(egi)) then
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLCID)
         call glc_final_mct( EClock_g, cdata_gg(egi), x2g_gg(egi), g2x_gg(egi))
         if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
      end if
   enddo

   !------------------------------------------------------------------------
   ! End the run cleanly
   !------------------------------------------------------------------------

   call shr_pio_finalize( )
   
   call shr_mpi_min(msize,msize0,mpicom_GLOID,'driver msize0',all=.true.)
   call shr_mpi_max(msize,msize1,mpicom_GLOID,'driver msize1',all=.true.)
   call shr_mpi_min(mrss,mrss0,mpicom_GLOID,'driver mrss0',all=.true.)
   call shr_mpi_max(mrss,mrss1,mpicom_GLOID,'driver mrss1',all=.true.)
   if (iamroot_CPLID )then
      call seq_timemgr_EClockGetData( EClock_d, curr_ymd=ymd, curr_tod=tod, dtime=dtime)
      write(logunit,'(//)')
      write(logunit,FormatA) subname, 'SUCCESSFUL TERMINATION OF CPL7-CCSM'
      write(logunit,FormatD) subname, '  at YMD,TOD = ',ymd,tod
      simDays = (endStep-begStep)*dtime/(24._r8*3600._r8)
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

   call t_stopf  ('DRIVER_FINAL')
   if (output_perf) then
      call t_prf(trim(timing_dir)//'/ccsm_timing', mpicom=mpicom_GLOID, &
                 output_thispe=output_perf)
   else
      call t_prf(trim(timing_dir)//'/ccsm_timing', mpicom=mpicom_GLOID)
   endif
   call t_finalizef()

end subroutine ccsm_final


!===============================================================================

subroutine t_drvstartf(string,cplrun,cplcom,budget,barrier)

   implicit none

   character(len=*),intent(in) :: string
   logical,intent(in),optional :: cplrun
   logical,intent(in),optional :: cplcom
   logical,intent(in),optional :: budget
   integer,intent(in),optional :: barrier

   character(len=128) :: strbar
   character(len=*),parameter :: strcpl = 'DRIVER_CPL_RUN'
   character(len=*),parameter :: strcom = 'DRIVER_CPL_COMM'
   character(len=*),parameter :: strbud = 'DRIVER_BUDGET'
   logical :: lcplrun,lcplcom,lbudget

!-------------------------------------------------------------------------------

   lcplrun  = .false.
   lcplcom  = .false.
   lbudget  = .false.
   if (present(cplrun)) then
      lcplrun = cplrun
   endif
   if (present(cplcom)) then
      lcplcom = cplcom
   endif
   if (present(budget)) then
      lbudget = budget
   endif

   if (present(barrier)) then
      strbar = trim(string)//'_BARRIER'
      call t_barrierf (trim(strbar), barrier)
   endif

   if (lcplrun) then
      call t_startf   (trim(strcpl))
   endif

   if (lcplcom) then
      call t_startf   (trim(strcom))
   endif

   if (lbudget) then
      call t_startf   (trim(strbud))
   endif

   call t_startf   (trim(string))

end subroutine t_drvstartf

!===============================================================================

subroutine t_drvstopf(string,cplrun,cplcom,budget)

   implicit none

   character(len=*),intent(in) :: string
   logical,intent(in),optional :: cplrun
   logical,intent(in),optional :: cplcom
   logical,intent(in),optional :: budget

   character(len=128) :: strbar
   character(len=*),parameter :: strcpl = 'DRIVER_CPL_RUN'
   character(len=*),parameter :: strcom = 'DRIVER_CPL_COMM'
   character(len=*),parameter :: strbud = 'DRIVER_BUDGET'
   logical :: lcplrun,lcplcom,lbudget

!-------------------------------------------------------------------------------

   lcplrun = .false.
   lcplcom = .false.
   lbudget = .false.
   if (present(cplrun)) then
      lcplrun = cplrun
   endif
   if (present(cplcom)) then
      lcplcom = cplcom
   endif
   if (present(budget)) then
      lbudget = budget
   endif

!  strbar = trim(string)//'_BARRIER'

   call t_stopf   (trim(string))

   if (lbudget) then
      call t_stopf   (trim(strbud))
   endif

   if (lcplrun) then
      call t_stopf   (trim(strcpl))
   endif

   if (lcplcom) then
      call t_stopf   (trim(strcom))
   endif

end subroutine t_drvstopf

!===============================================================================

subroutine seq_ccsm_printlogheader()

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
  write(logunit,F00) '     (Online documentation is available on the CCSM         '
  write(logunit,F00) '      Models page: http://www.cesm.ucar.edu/models/         '
  write(logunit,F00) '      License information is available as a link from above '
  write(logunit,F00) '------------------------------------------------------------'
  write(logunit,F00) '                DATE ',cdate, ' TIME ', ctime
  write(logunit,F00) '------------------------------------------------------------'
  write(logunit,*)' '
  write(logunit,*)' '

end subroutine seq_ccsm_printlogheader


#ifdef USE_ESMF_LIB

!===============================================================================

subroutine ccsm_comp_init(comp, importState, exportState, clock, rc)
   use esmfshr_attribute_mod
#ifdef ESMF_INTERFACE
   use atm_comp_mct, only: atm_register
   use lnd_comp_mct, only: lnd_register
   use ocn_comp_mct, only: ocn_register
   use ice_comp_mct, only: ice_register
   use glc_comp_mct, only: glc_register
   use rof_comp_mct, only: rof_register
   use wav_comp_mct, only: wav_register
#endif
   implicit none

   type(ESMF_CplComp)   :: comp
   type(ESMF_State)     :: importState, exportState
   type(ESMF_Clock)     :: clock
   integer, intent(out) :: rc

   ! Local variables
   type(ESMF_State)    :: attState, imp_state, exp_state
   integer             :: localrc
   type(ESMF_GridComp) :: atmComp, lndComp, iceComp, ocnComp, glcComp, &
                          rofComp, wavComp, cplComp
   type(ESMF_VM)       :: vm
   character(len=80)   :: str
   integer             :: rc2
   integer, dimension(1) :: rootList

   rc = ESMF_SUCCESS

   !------
   ! Create a state object to which the field level attributes will be 
   ! attached, and link the state to the specified component
   !------
   attState = ESMF_StateCreate(name="ccsm_atts", rc=localrc)
   if (localrc /= ESMF_SUCCESS) call shr_sys_abort('failed to create state for attributes')

   call ESMF_AttributeLink(comp, attState, rc=localrc)
   if (localrc /= ESMF_SUCCESS) call shr_sys_abort('failed to link attributes')

   !------
   ! Create and setup the model components
   ! import and export states are inout variables to register subroutines and their
   ! values are changed in each iteration and saved in the seq_comm_type array.
   !------

   call seq_comm_petlist(CPLID, cpl_petlist)
   call seq_map_register(cpl_petlist, comp, cplComp, imp_state, exp_state)
   call seq_comm_setcompstates(CPLID, cplComp, imp_state, exp_state)

#ifdef ESMF_INTERFACE
   do eai = 1, num_inst_atm
      call seq_comm_petlist(ATMID(eai),atm_petlist)
      call atm_register(atm_petlist, comp, atmComp, imp_state, exp_state)
      call seq_comm_setcompstates(ATMID(eai), atmComp, imp_state, exp_state)
   enddo
   do eli = 1, num_inst_lnd
      call seq_comm_petlist(LNDID(eli),lnd_petlist)
      call lnd_register(lnd_petlist, comp, lndComp, imp_state, exp_state)
      call seq_comm_setcompstates(LNDID(eli), lndComp, imp_state, exp_state)
   enddo
   do eii = 1, num_inst_ice
      call seq_comm_petlist(ICEID(eii),ice_petlist)
      call ice_register(ice_petlist, comp, iceComp, imp_state, exp_state)
      call seq_comm_setcompstates(ICEID(eii), iceComp, imp_state, exp_state)
   enddo
   do eoi = 1, num_inst_ocn
      call seq_comm_petlist(OCNID(eoi),ocn_petlist)
      call ocn_register(ocn_petlist, comp, ocnComp, imp_state, exp_state)
      call seq_comm_setcompstates(OCNID(eoi), ocnComp, imp_state, exp_state)
   enddo
   do egi = 1, num_inst_glc
      call seq_comm_petlist(GLCID(egi),glc_petlist)
      call glc_register(glc_petlist, comp, glcComp, imp_state, exp_state)
      call seq_comm_setcompstates(GLCID(egi), glcComp, imp_state, exp_state)
   enddo
   do eri = 1, num_inst_rof
      call seq_comm_petlist(ROFID(eri),rof_petlist)
      call rof_register(rof_petlist, comp, rofComp, imp_state, exp_state)
      call seq_comm_setcompstates(ROFID(eri), rofComp, imp_state, exp_state)
   enddo
   do ewi = 1, num_inst_wav
      call seq_comm_petlist(WAVID(ewi),wav_petlist)
      call wav_register(wav_petlist, comp, wavComp, imp_state, exp_state)
      call seq_comm_setcompstates(WAVID(ewi), wavComp, imp_state, exp_state)
   enddo
#endif

   !------
   ! Process the CESM initialization
   !------
   call ccsm_init()

   !------
   ! Set the application and field level attributes
   !------
   call esmfshr_attribute_appl_init(comp, rc=localrc)
   !call esmfshr_attribute_fields_init(attState, rc=localrc)

   !------
   ! Get the VM and root pet list to be used for the AttributeUpdate call
   !------
! get current
   !call ESMF_VMGetGlobal(vm, rc=localrc)
   call ESMF_VMGetCurrent(vm, rc=localrc)
   if (localrc /= 0) call shr_sys_abort('failed to get VM')

#ifdef ESMF_INTERFACE
!-------------------------------------------------------------------------
! The attribute handling part of the code is updated to loop
! through ensemble instances.
!-------------------------------------------------------------------------
   do eai = 1, num_inst_atm
      call seq_comm_petlist(ATMID(eai),atm_petlist)
      call seq_comm_getcompstates(ATMID(eai), atmComp)
      call ESMF_AttributeUpdate(atmComp, vm, rootList=atm_petlist, rc=localrc)
      if (localrc /= ESMF_SUCCESS) call shr_sys_abort('failed to update atm attributes')
   enddo 
   do eli = 1, num_inst_lnd
      call seq_comm_petlist(LNDID(eli),lnd_petlist)
      call seq_comm_getcompstates(LNDID(eli), lndComp)
      call ESMF_AttributeUpdate(lndComp, vm, rootList=lnd_petlist, rc=localrc)
      if (localrc /= ESMF_SUCCESS) call shr_sys_abort('failed to update lnd attributes')
   enddo 
   do eii = 1, num_inst_ice
      call seq_comm_petlist(ICEID(eii),ice_petlist)
      call seq_comm_getcompstates(ICEID(eii), iceComp)
      call ESMF_AttributeUpdate(iceComp, vm, rootList=ice_petlist, rc=localrc)
      if (localrc /= ESMF_SUCCESS) call shr_sys_abort('failed to update ice attributes')
   enddo
   do eoi = 1, num_inst_ocn
      call seq_comm_petlist(OCNID(eoi),ocn_petlist)
      call seq_comm_getcompstates(OCNID(eoi), ocnComp)
      call ESMF_AttributeUpdate(ocnComp, vm, rootList=ocn_petlist, rc=localrc)
      if (localrc /= ESMF_SUCCESS) call shr_sys_abort('failed to update ocn attributes')
   enddo
   do egi = 1, num_inst_glc
      call seq_comm_petlist(GLCID(egi),glc_petlist)
      call seq_comm_getcompstates(GLCID(egi), glcComp)
      call ESMF_AttributeUpdate(glcComp, vm, rootList=glc_petlist, rc=localrc)
      if (localrc /= ESMF_SUCCESS) call shr_sys_abort('failed to update glc attributes')
   enddo
   do eri = 1, num_inst_rof
      call seq_comm_petlist(ROFID(eri),rof_petlist)
      call seq_comm_getcompstates(ROFID(eri), rofComp)
      call ESMF_AttributeUpdate(rofComp, vm, rootList=rof_petlist, rc=localrc)
      if (localrc /= ESMF_SUCCESS) call shr_sys_abort('failed to update rof attributes')
   enddo 
   do ewi = 1, num_inst_wav
      call seq_comm_petlist(WAVID(ewi),wav_petlist)
      call seq_comm_getcompstates(WAVID(ewi), wavComp)
      call ESMF_AttributeUpdate(wavComp, vm, rootList=wav_petlist, rc=localrc)
      if (localrc /= ESMF_SUCCESS) call shr_sys_abort('failed to update wav attributes')
   enddo

   !------
   ! Write out all of the attributes to the CIM compliant XML file
   !------
   if (iamroot_GLOID) then

      call ESMF_AttributeWrite( &
              comp, &
              convention='CIM', &
              purpose='Model Component Simulation Description', &
              attwriteflag=ESMF_ATTWRITE_XML, rc=localrc)

   endif
#endif

   rc = localrc

end subroutine ccsm_comp_init

!===============================================================================

subroutine ccsm_comp_run(comp, importState, exportState, clock, rc)
   implicit none
   type(ESMF_CplComp)   :: comp
   type(ESMF_State)     :: importState, exportState
   type(ESMF_Clock)     :: clock
   integer, intent(out) :: rc

   rc = ESMF_SUCCESS

   call ccsm_run()

end subroutine ccsm_comp_run

!===============================================================================

subroutine ccsm_comp_final(comp, importState, exportState, clock, rc)
   implicit none
   type(ESMF_CplComp)   :: comp
   type(ESMF_State)     :: importState, exportState
   type(ESMF_Clock)     :: clock
   integer, intent(out) :: rc

   rc = ESMF_SUCCESS

   call ccsm_final()

end subroutine ccsm_comp_final


!===============================================================================
!
! This subroutine registers the initialization, run and finalization routines
! for the specified coupler component.  
!
subroutine ccsm_comp_register(comp, rc)
   implicit none
   type(ESMF_CplComp) :: comp
   integer, intent(out) :: rc

   rc = ESMF_SUCCESS

   call ESMF_CplCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
                                  userRoutine=ccsm_comp_init, rc=rc)
   call ESMF_CplCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
                                  userRoutine=ccsm_comp_run, rc=rc)
   call ESMF_CplCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, &
                                  userRoutine=ccsm_comp_final, rc=rc)

end subroutine ccsm_comp_register

!===============================================================================

#endif

end module ccsm_comp_mod

