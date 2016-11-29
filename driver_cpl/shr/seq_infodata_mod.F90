!===============================================================================
! SVN $Id: seq_infodata_mod.F90 68253 2015-02-18 22:24:57Z mvertens $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/drv/seq_mct/trunk_tags/drvseq5_1_15/shr/seq_infodata_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: seq_infodata_mod --- Module for input data shared between CCSM components
!
! !DESCRIPTION:
!
!     A module to get, put, and store some standard scalar data
!
! Typical usage:
!
!
! !REMARKS:
!
! !REVISION HISTORY:
!     2005-Nov-11 - E. Kluzek - creation of shr_inputinfo_mod
!     2007-Nov-15 - T. Craig - refactor for ccsm4 system and move to seq_infodata_mod
!
! !INTERFACE: ------------------------------------------------------------------

MODULE seq_infodata_mod

! !USES:

   use shr_kind_mod,      only : SHR_KIND_CS, SHR_KIND_CL, SHR_KIND_IN,      &
                                 SHR_KIND_R8, SHR_KIND_I8
   use shr_sys_mod,       only : shr_sys_flush, shr_sys_abort, shr_sys_getenv
   use seq_comm_mct,      only : logunit, loglevel, CPLID, seq_comm_gloroot
   use seq_comm_mct,      only : seq_comm_setptrs, seq_comm_iamroot, seq_comm_iamin
   use shr_orb_mod

   implicit none

   private  ! default private

! !PUBLIC TYPES:

   public :: seq_infodata_type

! !PUBLIC MEMBER FUNCTIONS

   public :: seq_infodata_Init            ! Initialize
   public :: seq_infodata_GetData         ! Get values from object
   public :: seq_infodata_PutData         ! Change values
   public :: seq_infodata_Print           ! print current info
   public :: seq_infodata_Exchange        ! exchange data across pes

! !PUBLIC DATA MEMBERS:

!EOP

   ! Strings of valid start_type options
   character(len=*), public, parameter :: seq_infodata_start_type_start = "startup"
   character(len=*), public, parameter :: seq_infodata_start_type_cont  = "continue"
   character(len=*), public, parameter :: seq_infodata_start_type_brnch = "branch"
   character(len=*), public, parameter :: seq_infodata_orb_fixed_year       = 'fixed_year'
   character(len=*), public, parameter :: seq_infodata_orb_variable_year    = 'variable_year'
   character(len=*), public, parameter :: seq_infodata_orb_fixed_parameters = 'fixed_parameters'

   ! InputInfo derived type

   type seq_infodata_type
      private     ! This type is opaque

      !--- set via namelist and held fixed ----
      character(SHR_KIND_CS)  :: cime_model      ! acme or cesm
      character(SHR_KIND_CL)  :: start_type      ! Type of startup
      character(SHR_KIND_CL)  :: case_name       ! Short case identification
      character(SHR_KIND_CL)  :: case_desc       ! Long description of this case
      character(SHR_KIND_CL)  :: model_version   ! Model version
      character(SHR_KIND_CS)  :: username        ! Current user
      character(SHR_KIND_CS)  :: hostname        ! Current machine
      character(SHR_KIND_CL)  :: timing_dir      ! Dir for timing files
      character(SHR_KIND_CL)  :: tchkpt_dir      ! Dir for timing checkpoint files
      logical                 :: atm_adiabatic   ! No surface models and atm adiabatic mode
      logical                 :: atm_ideal_phys  ! No surface models and atm ideal-physics
      logical                 :: aqua_planet     ! No ice/lnd, analytic ocn, perpetual time
      integer(SHR_KIND_IN)    :: aqua_planet_sst ! aqua planet analytic sst type
      logical                 :: run_barriers    ! barrier component run calls
      logical                 :: brnch_retain_casename   ! If branch and can use same casename
      logical                 :: read_restart    ! read the restart file, based on start_type
      character(SHR_KIND_CL)  :: restart_pfile   ! Restart pointer file
      character(SHR_KIND_CL)  :: restart_file    ! Full archive path to restart file
      logical                 :: single_column   ! single column mode
      real (SHR_KIND_R8)      :: scmlat          ! single column lat
      real (SHR_KIND_R8)      :: scmlon          ! single column lon
      character(SHR_KIND_CS)  :: logFilePostFix  ! postfix for output log files
      character(SHR_KIND_CL)  :: outPathRoot     ! root for output log files
      logical                 :: perpetual       ! perpetual flag
      integer(SHR_KIND_IN)    :: perpetual_ymd   ! perpetual date
      integer(SHR_KIND_IN)    :: orb_iyear       ! orbital year
      integer(SHR_KIND_IN)    :: orb_iyear_align ! model year associated with orb year
      character(SHR_KIND_CL)  :: orb_mode        ! orbital mode
      real(SHR_KIND_R8)       :: orb_eccen       ! See shr_orb_mod
      real(SHR_KIND_R8)       :: orb_obliq       ! See shr_orb_mod
      real(SHR_KIND_R8)       :: orb_mvelp       ! See shr_orb_mod
      real(SHR_KIND_R8)       :: orb_obliqr      ! See shr_orb_mod
      real(SHR_KIND_R8)       :: orb_lambm0      ! See shr_orb_mod
      real(SHR_KIND_R8)       :: orb_mvelpp      ! See shr_orb_mod
      character(SHR_KIND_CS)  :: wv_sat_scheme   ! Water vapor saturation pressure scheme
      real(SHR_KIND_R8)       :: wv_sat_transition_start ! Saturation transition range
      logical                 :: wv_sat_use_tables   ! Saturation pressure lookup tables
      real(SHR_KIND_R8)       :: wv_sat_table_spacing! Saturation pressure table resolution
      character(SHR_KIND_CS)  :: tfreeze_option  ! Freezing point calculation
      character(SHR_KIND_CL)  :: flux_epbal      ! selects E,P,R adjustment technique
      logical                 :: flux_albav      ! T => no diurnal cycle in ocn albedos
      logical                 :: flux_diurnal    ! T => diurnal cycle in atm/ocn fluxes
      real(SHR_KIND_R8)       :: wall_time_limit ! force stop time limit (hours)
      character(SHR_KIND_CS)  :: force_stop_at   ! when to force a stop (month, day, etc)
      character(SHR_KIND_CL)  :: atm_gnam        ! atm grid
      character(SHR_KIND_CL)  :: lnd_gnam        ! lnd grid
      character(SHR_KIND_CL)  :: ocn_gnam        ! ocn grid
      character(SHR_KIND_CL)  :: ice_gnam        ! ice grid
      character(SHR_KIND_CL)  :: rof_gnam        ! rof grid
      character(SHR_KIND_CL)  :: glc_gnam        ! glc grid
      character(SHR_KIND_CL)  :: wav_gnam        ! wav grid
      logical                 :: shr_map_dopole  ! pole corrections in shr_map_mod
      character(SHR_KIND_CL)  :: vect_map        ! vector mapping option, none, cart3d, cart3d_diag, cart3d_uvw, cart3d_uvw_diag
      character(SHR_KIND_CS)  :: aoflux_grid     ! grid for atm ocn flux calc
      integer                 :: cpl_decomp      ! coupler decomp
      character(SHR_KIND_CL)  :: cpl_seq_option  ! coupler sequencing option
      logical                 :: cpl_cdf64       ! use netcdf 64 bit offset, large file support
      logical                 :: do_budgets      ! do heat/water budgets diagnostics
      logical                 :: do_histinit     ! write out initial history file
      integer                 :: budget_inst     ! instantaneous budget level
      integer                 :: budget_daily    ! daily budget level
      integer                 :: budget_month    ! monthly budget level
      integer                 :: budget_ann      ! annual budget level
      integer                 :: budget_ltann    ! long term budget level written at end of year
      integer                 :: budget_ltend    ! long term budget level written at end of run
      logical                 :: drv_threading   ! is threading control in driver turned on
      logical                 :: histaux_a2x     ! cpl writes aux hist files: a2x every c2a comm
      logical                 :: histaux_a2x1hri ! cpl writes aux hist files: a2x 1hr instaneous values
      logical                 :: histaux_a2x1hr  ! cpl writes aux hist files: a2x 1hr
      logical                 :: histaux_a2x3hr  ! cpl writes aux hist files: a2x 3hr states
      logical                 :: histaux_a2x3hrp ! cpl writes aux hist files: a2x 3hr precip
      logical                 :: histaux_a2x24hr ! cpl writes aux hist files: a2x daily all
      logical                 :: histaux_l2x1yr  ! cpl writes aux hist files: l2x annual all
      logical                 :: histaux_l2x     ! cpl writes aux hist files: l2x every c2l comm
      logical                 :: histaux_r2x     ! cpl writes aux hist files: r2x every c2o comm
      logical                 :: histavg_atm     ! cpl writes atm fields in average history file
      logical                 :: histavg_lnd     ! cpl writes lnd fields in average history file
      logical                 :: histavg_ocn     ! cpl writes ocn fields in average history file
      logical                 :: histavg_ice     ! cpl writes ice fields in average history file
      logical                 :: histavg_rof     ! cpl writes rof fields in average history file
      logical                 :: histavg_glc     ! cpl writes glc fields in average history file
      logical                 :: histavg_wav     ! cpl writes wav fields in average history file
      logical                 :: histavg_xao     ! cpl writes flux xao fields in average history file
      real(SHR_KIND_R8)       :: eps_frac        ! fraction error tolerance
      real(SHR_KIND_R8)       :: eps_amask       ! atm mask error tolerance
      real(SHR_KIND_R8)       :: eps_agrid       ! atm grid error tolerance
      real(SHR_KIND_R8)       :: eps_aarea       ! atm area error tolerance
      real(SHR_KIND_R8)       :: eps_omask       ! ocn mask error tolerance
      real(SHR_KIND_R8)       :: eps_ogrid       ! ocn grid error tolerance
      real(SHR_KIND_R8)       :: eps_oarea       ! ocn area error tolerance
      logical                 :: mct_usealltoall ! flag for mct alltoall
      logical                 :: mct_usevector   ! flag for mct vector

      logical                 :: reprosum_use_ddpdd  ! use ddpdd algorithm
      real(SHR_KIND_R8)       :: reprosum_diffmax    ! maximum difference tolerance
      logical                 :: reprosum_recompute  ! recompute reprosum with nonscalable algorithm
                                                     ! if reprosum_diffmax is exceeded

      !--- set via namelist and may be time varying ---
      integer(SHR_KIND_IN)    :: info_debug      ! debug level
      logical                 :: bfbflag         ! turn on bfb option
      logical                 :: esmf_map_flag   ! do we use esmf mapping

      !--- set via components and held fixed ---
      logical                 :: atm_present     ! does component model exist
      logical                 :: atm_prognostic  ! does component model need input data from driver
      logical                 :: lnd_present     ! does component model exist
      logical                 :: lnd_prognostic  ! does component model need input data from driver
      logical                 :: rof_present     ! does rof component exist
      logical                 :: rofice_present  ! does rof have iceberg coupling on
      logical                 :: rof_prognostic  ! does rof component need input data
      logical                 :: flood_present   ! does rof have flooding on
      logical                 :: ocn_present     ! does component model exist
      logical                 :: ocn_prognostic  ! does component model need input data from driver
      logical                 :: ocnrof_prognostic ! does component need rof data
      logical                 :: ice_present     ! does component model exist
      logical                 :: ice_prognostic  ! does component model need input data from driver
      logical                 :: iceberg_prognostic ! does the ice model support icebergs
      logical                 :: glc_present     ! does component model exist
      logical                 :: glclnd_present  ! does glc have land coupling fields on
      logical                 :: glcocn_present  ! does glc have ocean runoff on
      logical                 :: glcice_present  ! does glc have iceberg coupling on
      logical                 :: glc_prognostic  ! does component model need input data from driver
      logical                 :: wav_present     ! does component model exist
      logical                 :: wav_prognostic  ! does component model need input data from driver
      logical                 :: esp_present     ! does component model exist
      logical                 :: esp_prognostic  ! does component model need input data from driver
      logical                 :: dead_comps      ! do we have dead models
      integer(SHR_KIND_IN)    :: atm_nx          ! nx, ny of "2d" grid
      integer(SHR_KIND_IN)    :: atm_ny          ! nx, ny of "2d" grid
      integer(SHR_KIND_IN)    :: lnd_nx          ! nx, ny of "2d" grid
      integer(SHR_KIND_IN)    :: lnd_ny          ! nx, ny of "2d" grid
      integer(SHR_KIND_IN)    :: ice_nx          ! nx, ny of "2d" grid
      integer(SHR_KIND_IN)    :: ice_ny          ! nx, ny of "2d" grid
      integer(SHR_KIND_IN)    :: ocn_nx          ! nx, ny of "2d" grid
      integer(SHR_KIND_IN)    :: ocn_ny          ! nx, ny of "2d" grid
      integer(SHR_KIND_IN)    :: rof_nx          ! nx, ny of "2d" grid
      integer(SHR_KIND_IN)    :: rof_ny          ! nx, ny of "2d" grid
      integer(SHR_KIND_IN)    :: glc_nx          ! nx, ny of "2d" grid
      integer(SHR_KIND_IN)    :: glc_ny          ! nx, ny of "2d" grid
      integer(SHR_KIND_IN)    :: wav_nx          ! nx, ny of "2d" grid
      integer(SHR_KIND_IN)    :: wav_ny          ! nx, ny of "2d" grid

      !--- set via components and may be time varying ---
      real(SHR_KIND_R8)       :: nextsw_cday     ! calendar of next atm shortwave
      real(SHR_KIND_R8)       :: precip_fact     ! precip factor
      integer(SHR_KIND_IN)    :: atm_phase       ! atm phase
      integer(SHR_KIND_IN)    :: lnd_phase       ! lnd phase
      integer(SHR_KIND_IN)    :: ice_phase       ! ice phase
      integer(SHR_KIND_IN)    :: ocn_phase       ! ocn phase
      integer(SHR_KIND_IN)    :: glc_phase       ! glc phase
      integer(SHR_KIND_IN)    :: rof_phase       ! rof phase
      integer(SHR_KIND_IN)    :: wav_phase       ! wav phase
      integer(SHR_KIND_IN)    :: esp_phase       ! esp phase
      logical                 :: atm_aero        ! atmosphere aerosols
      logical                 :: glcrun_alarm    ! glc run alarm
      logical                 :: glc_g2lupdate   ! update glc2lnd fields in lnd model
      logical                 :: pause_alarm     ! active components should write restart files
      logical                 :: resume_alarm    ! active components should 'resume' from provided files
      real(shr_kind_r8) :: max_cplstep_time  ! abort if cplstep time exceeds this value
      !--- set from restart file ---
      character(SHR_KIND_CL)  :: rest_case_name  ! Short case identification
   end type seq_infodata_type

   ! --- public interfaces --------------------------------------------------------
   interface seq_infodata_GetData
     module procedure seq_infodata_GetData_explicit
#ifndef CPRPGI
     module procedure seq_infodata_GetData_bytype
#endif
! ^ ifndef CPRPGI
   end interface

   interface seq_infodata_PutData
     module procedure seq_infodata_PutData_explicit
#ifndef CPRPGI
     module procedure seq_infodata_PutData_bytype
#endif
! ^ ifndef CPRPGI
   end interface

   ! --- Private local data -------------------------------------------------------

   character(len=*),parameter :: sp_str = 'str_undefined'

!===============================================================================
CONTAINS
!===============================================================================

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_infodata_Init -- read in CCSM shared namelist
!
! !DESCRIPTION:
!
!     Read in input from seq_infodata_inparm namelist, output ccsm derived type for
!     miscillaneous info.
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE seq_infodata_Init( infodata, nmlfile, ID, pioid)

! !USES:

   use shr_file_mod,   only : shr_file_getUnit, shr_file_freeUnit
   use shr_string_mod, only : shr_string_toUpper, shr_string_listAppend
   use shr_mpi_mod,    only : shr_mpi_bcast
   use seq_io_read_mod
   use pio, only : file_desc_t
   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(seq_infodata_type), intent(INOUT) :: infodata  ! infodata object
   character(len=*),        intent(IN)    :: nmlfile   ! Name-list filename
   integer(SHR_KIND_IN),    intent(IN)    :: ID        ! seq_comm ID
   type(file_desc_T) :: pioid
!EOP

    !----- local -----
    character(len=*),    parameter :: subname = '(seq_infodata_Init) '
    integer(SHR_KIND_IN),parameter :: aqua_perpetual_ymd = 321

    integer :: mpicom             ! MPI communicator
    integer :: ierr               ! I/O error code
    integer :: unitn              ! Namelist unit number to read

    !------ namelist -----
    character(SHR_KIND_CS) :: cime_model         ! acme or cesm
    character(SHR_KIND_CL) :: case_desc          ! Case long description
    character(SHR_KIND_CL) :: case_name          ! Case short name
    character(SHR_KIND_CL) :: model_version      ! Model version
    character(SHR_KIND_CS) :: username           ! Current user
    character(SHR_KIND_CS) :: hostname           ! Current machine
    character(SHR_KIND_CL) :: start_type         ! Startup-type: startup, continue, branch
    character(SHR_KIND_CL) :: timing_dir         ! Dir for timing files
    character(SHR_KIND_CL) :: tchkpt_dir         ! Dir for timing checkpoint files
    logical                :: atm_adiabatic      ! Atmosphere adiabatic physics mode
    logical                :: atm_ideal_phys     ! Atmosphere idealized physics mode
    logical                :: aqua_planet        ! Aqua-planet mode (surface is all ocean)
    integer(SHR_KIND_IN)   :: aqua_planet_sst    ! analytic sst field
    logical                :: run_barriers       ! barrier component run calls
    logical                :: brnch_retain_casename ! If retain casename for branch
    integer(SHR_KIND_IN)   :: info_debug         ! debug flag
    logical                :: bfbflag            ! bit for bit flag
    logical                :: esmf_map_flag      ! esmf mapping flag
    character(SHR_KIND_CL) :: restart_pfile      ! Restart pointer filename
    character(SHR_KIND_CL) :: restart_file       ! Restart filename

    logical                :: single_column      ! single column mode
    real (SHR_KIND_R8)     :: scmlat             ! single column mode latitude
    real (SHR_KIND_R8)     :: scmlon             ! single column mode longitude
    character(SHR_KIND_CS) :: logFilePostFix     ! postfix for output log files
    character(SHR_KIND_CL) :: outPathRoot        ! root output files
    logical                :: perpetual          ! perpetual mode
    integer(SHR_KIND_IN)   :: perpetual_ymd      ! perpetual ymd
    integer(SHR_KIND_IN)   :: orb_iyear          ! orbital year
    integer(SHR_KIND_IN)   :: orb_iyear_align    ! model year associated with orb year
    character(SHR_KIND_CL) :: orb_mode           ! orbital mode
    real(SHR_KIND_R8)      :: orb_obliq          ! Obliquity of orbit
    real(SHR_KIND_R8)      :: orb_eccen          ! Eccentricity of orbit
    real(SHR_KIND_R8)      :: orb_mvelp          ! Location of vernal equinox
    real(SHR_KIND_R8)      :: orb_obliqr         ! Obliquity in radians
    real(SHR_KIND_R8)      :: orb_lambm0         ! lon of per at vernal equ
    real(SHR_KIND_R8)      :: orb_mvelpp         ! mvelp plus pi
    character(SHR_KIND_CS) :: wv_sat_scheme      ! Water vapor saturation pressure scheme
    real(SHR_KIND_R8)      :: wv_sat_transition_start! Saturation transition range
    logical                :: wv_sat_use_tables  ! Saturation pressure lookup tables
    real(SHR_KIND_R8)      :: wv_sat_table_spacing   ! Saturation pressure table resolution
    character(SHR_KIND_CS) :: tfreeze_option     ! Freezing point calculation
    character(SHR_KIND_CL) :: flux_epbal         ! selects E,P,R adjustment technique
    logical                :: flux_albav         ! T => no diurnal cycle in ocn albedos
    logical                :: flux_diurnal       ! T => diurnal cycle in atm/ocn fluxes
    real(SHR_KIND_R8)      :: wall_time_limit    ! force stop time limit (hours)
    character(SHR_KIND_CS) :: force_stop_at      ! when to force a stop (month, day, etc)
    character(SHR_KIND_CL) :: atm_gnam           ! atm grid
    character(SHR_KIND_CL) :: lnd_gnam           ! lnd grid
    character(SHR_KIND_CL) :: ocn_gnam           ! ocn grid
    character(SHR_KIND_CL) :: ice_gnam           ! ice grid
    character(SHR_KIND_CL) :: rof_gnam           ! rof grid
    character(SHR_KIND_CL) :: glc_gnam           ! glc grid
    character(SHR_KIND_CL) :: wav_gnam           ! wav grid
    logical                :: shr_map_dopole     ! pole corrections in shr_map_mod
    character(SHR_KIND_CL) :: vect_map           ! vector mapping option
    character(SHR_KIND_CS) :: aoflux_grid        ! grid for atm ocn flux calc
    integer                :: cpl_decomp         ! coupler decomp
    character(SHR_KIND_CL) :: cpl_seq_option     ! coupler sequencing option
    logical                :: cpl_cdf64          ! use netcdf 64 bit offset, large file support
    logical                :: do_budgets         ! do heat/water budgets diagnostics
    logical                :: do_histinit        ! write out initial history file
    integer                :: budget_inst        ! instantaneous budget level
    integer                :: budget_daily       ! daily budget level
    integer                :: budget_month       ! monthly budget level
    integer                :: budget_ann         ! annual budget level
    integer                :: budget_ltann       ! long term budget level written at end of year
    integer                :: budget_ltend       ! long term budget level written at end of run
    logical                :: histaux_a2x        ! cpl writes aux hist files: a2x every c2a comm
    logical                :: histaux_a2x1hri    ! cpl writes aux hist files: a2x 1hr instaneous values
    logical                :: histaux_a2x1hr     ! cpl writes aux hist files: a2x 1hr
    logical                :: histaux_a2x3hr     ! cpl writes aux hist files: a2x 3hr states
    logical                :: histaux_a2x3hrp    ! cpl writes aux hist files: a2x 2hr precip
    logical                :: histaux_a2x24hr    ! cpl writes aux hist files: a2x daily all
    logical                :: histaux_l2x1yr     ! cpl writes aux hist files: l2x annual all
    logical                :: histaux_l2x        ! cpl writes aux hist files: l2x every c2l comm
    logical                :: histaux_r2x        ! cpl writes aux hist files: r2x every c2o comm
    logical                :: histavg_atm        ! cpl writes atm fields in average history file
    logical                :: histavg_lnd        ! cpl writes lnd fields in average history file
    logical                :: histavg_ocn        ! cpl writes ocn fields in average history file
    logical                :: histavg_ice        ! cpl writes ice fields in average history file
    logical                :: histavg_rof        ! cpl writes rof fields in average history file
    logical                :: histavg_glc        ! cpl writes glc fields in average history file
    logical                :: histavg_wav        ! cpl writes wav fields in average history file
    logical                :: histavg_xao        ! cpl writes flux xao fields in average history file
    logical                :: drv_threading      ! is threading control in driver turned on
    real(SHR_KIND_R8)      :: eps_frac           ! fraction error tolerance
    real(SHR_KIND_R8)      :: eps_amask          ! atm mask error tolerance
    real(SHR_KIND_R8)      :: eps_agrid          ! atm grid error tolerance
    real(SHR_KIND_R8)      :: eps_aarea          ! atm area error tolerance
    real(SHR_KIND_R8)      :: eps_omask          ! ocn mask error tolerance
    real(SHR_KIND_R8)      :: eps_ogrid          ! ocn grid error tolerance
    real(SHR_KIND_R8)      :: eps_oarea          ! ocn area error tolerance
    logical                :: reprosum_use_ddpdd ! use ddpdd algorithm
    real(SHR_KIND_R8)      :: reprosum_diffmax   ! maximum difference tolerance
    logical                :: reprosum_recompute ! recompute reprosum with nonscalable algorithm
                                                 ! if reprosum_diffmax is exceeded
    logical                :: mct_usealltoall    ! flag for mct alltoall
    logical                :: mct_usevector      ! flag for mct vector
    real(shr_kind_r8) :: max_cplstep_time  ! abort if cplstep time exceeds this value

    namelist /seq_infodata_inparm/  &
         cime_model, case_desc, case_name, start_type, tchkpt_dir,     &
         model_version, username, hostname, timing_dir,    &
         atm_adiabatic, atm_ideal_phys, aqua_planet,aqua_planet_sst, &
         brnch_retain_casename, info_debug, bfbflag,       &
         restart_pfile, restart_file, run_barriers,        &
         single_column, scmlat, force_stop_at,             &
         scmlon, logFilePostFix, outPathRoot, flux_diurnal,&
         perpetual, perpetual_ymd, flux_epbal, flux_albav, &
         orb_iyear_align, orb_mode, wall_time_limit,       &
         orb_iyear, orb_obliq, orb_eccen, orb_mvelp,       &
         wv_sat_scheme, wv_sat_transition_start,           &
         wv_sat_use_tables, wv_sat_table_spacing,          &
         tfreeze_option,                                      &
         ice_gnam, rof_gnam, glc_gnam, wav_gnam,           &
         atm_gnam, lnd_gnam, ocn_gnam, cpl_decomp,         &
         shr_map_dopole, vect_map, aoflux_grid, do_histinit,  &
         do_budgets, drv_threading,                        &
         budget_inst, budget_daily, budget_month,          &
         budget_ann, budget_ltann, budget_ltend,           &
         histaux_a2x,histaux_a2x1hri,histaux_a2x1hr,       &
         histaux_a2x3hr,histaux_a2x3hrp,                   &
         histaux_a2x24hr,histaux_l2x   ,histaux_r2x,       &
         histavg_atm, histavg_lnd, histavg_ocn, histavg_ice, &
         histavg_rof, histavg_glc, histavg_wav, histavg_xao, &
         histaux_l2x1yr, cpl_seq_option,                   &
         cpl_cdf64, eps_frac, eps_amask,                   &
         eps_agrid, eps_aarea, eps_omask, eps_ogrid,       &
         eps_oarea, esmf_map_flag,                         &
         reprosum_use_ddpdd, reprosum_diffmax, reprosum_recompute, &
         mct_usealltoall, mct_usevector, max_cplstep_time

!-------------------------------------------------------------------------------

    call seq_comm_setptrs(ID,mpicom=mpicom)

    !---------------------------------------------------------------------------
    ! Set infodata on root pe
    !---------------------------------------------------------------------------
    if (seq_comm_iamroot(ID)) then

       !---------------------------------------------------------------------------
       ! Set namelist defaults
       !---------------------------------------------------------------------------
       cime_model            = 'unknown'
       case_desc             = ' '
       case_name             = ' '
       model_version         = 'unknown'
       username              = 'unknown'
       hostname              = 'unknown'
       timing_dir            = '.'
       tchkpt_dir            = '.'
       start_type            = ' '
       atm_ideal_phys        = .false.
       atm_adiabatic         = .false.
       aqua_planet           = .false.
       aqua_planet_sst       = 1
       run_barriers          = .false.
       brnch_retain_casename = .false.
       info_debug            = 1
       bfbflag               = .false.
       esmf_map_flag         = .false.
       restart_pfile         = 'rpointer.drv'
       restart_file          = trim(sp_str)
       single_column         = .false.
       scmlat                = -999.
       scmlon                = -999.
       logFilePostFix        = '.log'
       outPathRoot           = './'
       perpetual             = .false.
       perpetual_ymd         = -999
       orb_mode              = seq_infodata_orb_fixed_year
       orb_iyear             = SHR_ORB_UNDEF_INT
       orb_iyear_align       = SHR_ORB_UNDEF_INT
       orb_obliq             = SHR_ORB_UNDEF_REAL
       orb_eccen             = SHR_ORB_UNDEF_REAL
       orb_mvelp             = SHR_ORB_UNDEF_REAL
       wv_sat_scheme         = "GoffGratch"
       wv_sat_transition_start = 20.0
       wv_sat_use_tables     = .false.
       wv_sat_table_spacing  = 1.0
       tfreeze_option        = 'minus1p8'
       flux_epbal            = 'off'
       flux_albav            = .false.
       flux_diurnal          = .false.
       wall_time_limit       = -1.0
       force_stop_at         = 'month'
       atm_gnam              = 'undefined'
       lnd_gnam              = 'undefined'
       ocn_gnam              = 'undefined'
       ice_gnam              = 'undefined'
       rof_gnam              = 'undefined'
       glc_gnam              = 'undefined'
       wav_gnam              = 'undefined'
       shr_map_dopole        = .true.
       vect_map              = 'cart3d'
       aoflux_grid           = 'ocn'
       cpl_decomp            = 0
       cpl_seq_option        = 'CESM1_MOD'
       cpl_cdf64             = .true.
       do_budgets            = .false.
       do_histinit           = .false.
       budget_inst           = 0
       budget_daily          = 0
       budget_month          = 1
       budget_ann            = 1
       budget_ltann          = 1
       budget_ltend          = 0
       histaux_a2x           = .false.
       histaux_a2x1hri       = .false.
       histaux_a2x1hr        = .false.
       histaux_a2x3hr        = .false.
       histaux_a2x3hrp       = .false.
       histaux_a2x24hr       = .false.
       histaux_l2x1yr        = .false.
       histaux_l2x           = .false.
       histaux_r2x           = .false.
       histavg_atm           = .true.
       histavg_lnd           = .true.
       histavg_ocn           = .true.
       histavg_ice           = .true.
       histavg_rof           = .true.
       histavg_glc           = .true.
       histavg_wav           = .true.
       histavg_xao           = .true.
       drv_threading         = .false.
       eps_frac              = 1.0e-02_SHR_KIND_R8
       eps_amask             = 1.0e-13_SHR_KIND_R8
       eps_agrid             = 1.0e-12_SHR_KIND_R8
       eps_aarea             = 9.0e-07_SHR_KIND_R8
       eps_omask             = 1.0e-06_SHR_KIND_R8
       eps_ogrid             = 1.0e-02_SHR_KIND_R8
       eps_oarea             = 1.0e-01_SHR_KIND_R8
       reprosum_use_ddpdd    = .false.
       reprosum_diffmax      = -1.0e-8
       reprosum_recompute    = .false.
       mct_usealltoall       = .false.
       mct_usevector         = .false.
       max_cplstep_time      = 0.0
       !---------------------------------------------------------------------------
       ! Read in namelist
       !---------------------------------------------------------------------------
       unitn = shr_file_getUnit()
       write(logunit,"(A)") subname,' read seq_infodata_inparm namelist from: '//trim(nmlfile)
       open( unitn, file=trim(nmlfile), status='old' )
       ierr = 1
       do while( ierr /= 0 )
          read(unitn,nml=seq_infodata_inparm,iostat=ierr)
          if (ierr < 0) then
             call shr_sys_abort( subname//':: namelist read returns an'// &
                                 ' end of file or end of record condition' )
          end if
       end do
       close(unitn)
       call shr_file_freeUnit( unitn )

       !---------------------------------------------------------------------------
       ! Set infodata on root pe
       !---------------------------------------------------------------------------
       infodata%cime_model            = cime_model
       infodata%case_desc             = case_desc
       infodata%case_name             = case_name
       infodata%model_version         = model_version
       infodata%username              = username
       infodata%hostname              = hostname
       infodata%start_type            = start_type
       infodata%timing_dir            = timing_dir
       infodata%tchkpt_dir            = tchkpt_dir
       infodata%atm_ideal_phys        = atm_ideal_phys
       infodata%atm_adiabatic         = atm_adiabatic
       infodata%aqua_planet           = aqua_planet
       infodata%aqua_planet_sst       = aqua_planet_sst
       infodata%run_barriers          = run_barriers
       infodata%brnch_retain_casename = brnch_retain_casename
       infodata%restart_pfile         = restart_pfile
       infodata%restart_file          = restart_file
       infodata%single_column         = single_column
       infodata%scmlat                = scmlat
       infodata%scmlon                = scmlon
       infodata%logFilePostFix        = logFilePostFix
       infodata%outPathRoot           = outPathRoot
       infodata%perpetual             = perpetual
       infodata%perpetual_ymd         = perpetual_ymd
       infodata%wv_sat_scheme         = wv_sat_scheme
       infodata%wv_sat_transition_start = wv_sat_transition_start
       infodata%wv_sat_use_tables     = wv_sat_use_tables
       infodata%wv_sat_table_spacing  = wv_sat_table_spacing
       infodata%tfreeze_option        = tfreeze_option
       infodata%flux_epbal            = flux_epbal
       infodata%flux_albav            = flux_albav
       infodata%flux_diurnal          = flux_diurnal
       infodata%wall_time_limit       = wall_time_limit
       infodata%force_stop_at         = force_stop_at
       infodata%atm_gnam              = atm_gnam
       infodata%lnd_gnam              = lnd_gnam
       infodata%ocn_gnam              = ocn_gnam
       infodata%ice_gnam              = ice_gnam
       infodata%rof_gnam              = rof_gnam
       infodata%glc_gnam              = glc_gnam
       infodata%wav_gnam              = wav_gnam
       infodata%shr_map_dopole        = shr_map_dopole
       infodata%vect_map              = vect_map
       infodata%aoflux_grid           = aoflux_grid
       infodata%cpl_decomp            = cpl_decomp
       infodata%cpl_seq_option        = cpl_seq_option
       infodata%cpl_cdf64             = cpl_cdf64
       infodata%do_budgets            = do_budgets
       infodata%do_histinit           = do_histinit
       infodata%budget_inst           = budget_inst
       infodata%budget_daily          = budget_daily
       infodata%budget_month          = budget_month
       infodata%budget_ann            = budget_ann
       infodata%budget_ltann          = budget_ltann
       infodata%budget_ltend          = budget_ltend
       infodata%histaux_a2x           = histaux_a2x
       infodata%histaux_a2x1hri       = histaux_a2x1hri
       infodata%histaux_a2x1hr        = histaux_a2x1hr
       infodata%histaux_a2x3hr        = histaux_a2x3hr
       infodata%histaux_a2x3hrp       = histaux_a2x3hrp
       infodata%histaux_a2x24hr       = histaux_a2x24hr
       infodata%histaux_l2x1yr        = histaux_l2x1yr
       infodata%histaux_l2x           = histaux_l2x
       infodata%histaux_r2x           = histaux_r2x
       infodata%histavg_atm           = histavg_atm
       infodata%histavg_lnd           = histavg_lnd
       infodata%histavg_ocn           = histavg_ocn
       infodata%histavg_ice           = histavg_ice
       infodata%histavg_rof           = histavg_rof
       infodata%histavg_glc           = histavg_glc
       infodata%histavg_wav           = histavg_wav
       infodata%histavg_xao           = histavg_xao
       infodata%drv_threading         = drv_threading
       infodata%eps_frac              = eps_frac
       infodata%eps_amask             = eps_amask
       infodata%eps_agrid             = eps_agrid
       infodata%eps_aarea             = eps_aarea
       infodata%eps_omask             = eps_omask
       infodata%eps_ogrid             = eps_ogrid
       infodata%eps_oarea             = eps_oarea
       infodata%reprosum_use_ddpdd    = reprosum_use_ddpdd
       infodata%reprosum_diffmax      = reprosum_diffmax
       infodata%reprosum_recompute    = reprosum_recompute
       infodata%mct_usealltoall       = mct_usealltoall
       infodata%mct_usevector         = mct_usevector

       infodata%info_debug            = info_debug
       infodata%bfbflag               = bfbflag
       infodata%esmf_map_flag         = esmf_map_flag

       infodata%atm_present = .true.
       infodata%lnd_present = .true.
       infodata%rof_present = .true.
       infodata%rofice_present = .true.
       infodata%flood_present = .true.
       infodata%ocn_present = .true.
       infodata%ice_present = .true.
       infodata%glc_present = .true.
       infodata%wav_present = .true.
       infodata%glclnd_present = .true.
       infodata%glcocn_present = .true.
       infodata%glcice_present = .true.
       infodata%esp_present = .true.

       infodata%atm_prognostic = .false.
       infodata%lnd_prognostic = .false.
       infodata%rof_prognostic = .false.
       infodata%ocn_prognostic = .false.
       infodata%ocnrof_prognostic = .false.
       infodata%ice_prognostic = .false.
       infodata%glc_prognostic = .false.
       infodata%wav_prognostic = .false.
       infodata%iceberg_prognostic = .false.
       infodata%esp_prognostic = .false.
       infodata%dead_comps = .false.

       infodata%atm_nx = 0
       infodata%atm_ny = 0
       infodata%lnd_nx = 0
       infodata%lnd_ny = 0
       infodata%rof_nx = 0
       infodata%rof_ny = 0
       infodata%ice_nx = 0
       infodata%ice_ny = 0
       infodata%ocn_nx = 0
       infodata%ocn_ny = 0
       infodata%glc_nx = 0
       infodata%glc_ny = 0
       infodata%wav_nx = 0
       infodata%wav_ny = 0

       infodata%nextsw_cday = -1.0_SHR_KIND_R8
       infodata%precip_fact =  1.0_SHR_KIND_R8
       infodata%atm_phase = 1
       infodata%lnd_phase = 1
       infodata%ocn_phase = 1
       infodata%ice_phase = 1
       infodata%glc_phase = 1
       infodata%rof_phase = 1
       infodata%wav_phase = 1
       infodata%atm_aero     = .false.
       infodata%glcrun_alarm = .false.
       infodata%glc_g2lupdate= .false.
       infodata%max_cplstep_time = max_cplstep_time
       infodata%pause_alarm  = .false.
       infodata%resume_alarm = .false.


       !---------------------------------------------------------------
       ! check orbital mode, reset unused parameters, validate settings
       !---------------------------------------------------------------
       if (trim(orb_mode) == trim(seq_infodata_orb_fixed_year)) then
          orb_obliq = SHR_ORB_UNDEF_REAL
          orb_eccen = SHR_ORB_UNDEF_REAL
          orb_mvelp = SHR_ORB_UNDEF_REAL
          if (orb_iyear == SHR_ORB_UNDEF_INT) then
             write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
             write(logunit,*) trim(subname),' ERROR: fixed_year settings = ',orb_iyear
             call shr_sys_abort(subname//' ERROR: invalid settings for orb_mode '//trim(orb_mode))
          endif
       elseif (trim(orb_mode) == trim(seq_infodata_orb_variable_year)) then
          orb_obliq = SHR_ORB_UNDEF_REAL
          orb_eccen = SHR_ORB_UNDEF_REAL
          orb_mvelp = SHR_ORB_UNDEF_REAL
          if (orb_iyear       == SHR_ORB_UNDEF_INT .or. &
              orb_iyear_align == SHR_ORB_UNDEF_INT) then
             write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
             write(logunit,*) trim(subname),' ERROR: variable_year settings = ',orb_iyear,orb_iyear_align
             call shr_sys_abort(subname//' ERROR: invalid settings for orb_mode '//trim(orb_mode))
          endif
       elseif (trim(orb_mode) == trim(seq_infodata_orb_fixed_parameters)) then
          !-- force orb_iyear to undef to make sure shr_orb_params works properly
          orb_iyear = SHR_ORB_UNDEF_INT
          orb_iyear_align = SHR_ORB_UNDEF_INT
          if (orb_eccen == SHR_ORB_UNDEF_REAL .or. &
              orb_obliq == SHR_ORB_UNDEF_REAL .or. &
              orb_mvelp == SHR_ORB_UNDEF_REAL) then
             write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
             write(logunit,*) trim(subname),' ERROR: orb_eccen = ',orb_eccen
             write(logunit,*) trim(subname),' ERROR: orb_obliq = ',orb_obliq
             write(logunit,*) trim(subname),' ERROR: orb_mvelp = ',orb_mvelp
             call shr_sys_abort(subname//' ERROR: invalid settings for orb_mode '//trim(orb_mode))
          endif
       else
          call shr_sys_abort(subname//' ERROR: invalid orb_mode '//trim(orb_mode))
       endif

       call shr_orb_params(orb_iyear, orb_eccen, orb_obliq, orb_mvelp, &
                           orb_obliqr, orb_lambm0, orb_mvelpp, .true.)

       infodata%orb_mode   = orb_mode
       infodata%orb_iyear  = orb_iyear
       infodata%orb_iyear_align = orb_iyear_align
       infodata%orb_eccen  = orb_eccen
       infodata%orb_obliq  = orb_obliq
       infodata%orb_mvelp  = orb_mvelp
       infodata%orb_obliqr = orb_obliqr
       infodata%orb_lambm0 = orb_lambm0
       infodata%orb_mvelpp = orb_mvelpp

       !--- Derive a few things ---
       infodata%rest_case_name = ' '
       infodata%read_restart = .false.
       if (trim(start_type) == trim(seq_infodata_start_type_cont) .or. &
           trim(start_type) == trim(seq_infodata_start_type_brnch)) then
          infodata%read_restart = .true.
       endif

    end if

    !-----------------------------------------------------
    ! Read Restart (seq_io_read must be called on all pes)
    !-----------------------------------------------------
    call shr_mpi_bcast(infodata%read_restart,mpicom)
    if (infodata%read_restart) then
       !--- read rpointer if restart_file is set to sp_str ---
       if (seq_comm_iamroot(ID)) then
          if (trim(infodata%restart_file) == trim(sp_str)) then
             unitn = shr_file_getUnit()
             if (loglevel > 0) write(logunit,"(3A)") subname," read rpointer file ", &
                trim(infodata%restart_pfile)
             open(unitn, file=infodata%restart_pfile, form='FORMATTED', status='old',iostat=ierr)
             if (ierr < 0) then
                call shr_sys_abort( subname//':: rpointer file open returns an'// &
                     ' error condition' )
             end if
             read(unitn,'(a)', iostat=ierr) infodata%restart_file
             if (ierr < 0) then
                call shr_sys_abort( subname//':: rpointer file read returns an'// &
                     ' error condition' )
             end if
             close(unitn)
             call shr_file_freeUnit( unitn )
             write(logunit,"(3A)") subname,' restart file from rpointer= ', &
                trim(infodata%restart_file)
          endif
       endif
       call shr_mpi_bcast(infodata%restart_file,mpicom)
       !--- NOTE: use CPLID here because seq_io is only value on CPLID
       if (seq_comm_iamin(CPLID)) then
          call seq_io_read(infodata%restart_file,pioid,infodata%nextsw_cday   ,'seq_infodata_nextsw_cday')
          call seq_io_read(infodata%restart_file,pioid,infodata%precip_fact   ,'seq_infodata_precip_fact')
          call seq_io_read(infodata%restart_file,pioid,infodata%rest_case_name,'seq_infodata_case_name')
       endif
       !--- Send from CPLID ROOT to GLOBALID ROOT, use bcast as surrogate
       call shr_mpi_bcast(infodata%nextsw_cday,mpicom,pebcast=seq_comm_gloroot(CPLID))
       call shr_mpi_bcast(infodata%precip_fact,mpicom,pebcast=seq_comm_gloroot(CPLID))
       call shr_mpi_bcast(infodata%rest_case_name,mpicom,pebcast=seq_comm_gloroot(CPLID))
    endif

    if (seq_comm_iamroot(ID)) then
       if (infodata%aqua_planet) then
          infodata%atm_present = .true.
          infodata%lnd_present = .false.
          infodata%rof_present = .false.
          infodata%rofice_present = .false.
          infodata%flood_present = .false.
          infodata%ice_present = .false.
          infodata%ocn_present = .true.
          infodata%glc_present = .false.
          infodata%wav_present = .false.
          infodata%glclnd_present = .false.
          infodata%glcocn_present = .false.
          infodata%glcice_present = .false.
          infodata%esp_present = .false.
       end if
       if (infodata%atm_adiabatic .or. infodata%atm_ideal_phys) then
          infodata%atm_present = .true.
          infodata%lnd_present = .false.
          infodata%rof_present = .false.
          infodata%rofice_present = .false.
          infodata%flood_present = .false.
          infodata%ice_present = .false.
          infodata%ocn_present = .false.
          infodata%glc_present = .false.
          infodata%wav_present = .false.
          infodata%glclnd_present = .false.
          infodata%glcocn_present = .false.
          infodata%glcice_present = .false.
          infodata%esp_present = .false.
       end if

       if ( infodata%aqua_planet ) then
          infodata%aqua_planet_sst = 1
          infodata%perpetual      = .true.
          infodata%perpetual_ymd  = aqua_perpetual_ymd
       endif

       ! --- Error check the input values ------
       call seq_infodata_Check( infodata )

    end if

    call seq_infodata_bcast(infodata,mpicom)

END SUBROUTINE seq_infodata_Init

!===============================================================================
!===============================================================================
! !IROUTINE: seq_infodata_GetData_explicit -- Get values from infodata object
!
! !DESCRIPTION:
!
!     Get values out of the infodata object.
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE seq_infodata_GetData_explicit( infodata, cime_model, case_name, case_desc, timing_dir,  &
           model_version, username, hostname, rest_case_name, tchkpt_dir,     &
           start_type, restart_pfile, restart_file, perpetual, perpetual_ymd, &
           aqua_planet,aqua_planet_sst, atm_ideal_phys, atm_adiabatic, brnch_retain_casename, &
           single_column, scmlat,scmlon,logFilePostFix, outPathRoot,          &
           atm_present, atm_prognostic, lnd_present, lnd_prognostic, rof_prognostic, &
           rof_present, ocn_present, ocn_prognostic, ocnrof_prognostic,       &
           ice_present, ice_prognostic, glc_present, glc_prognostic,          &
           flood_present, wav_present, wav_prognostic, rofice_present,        &
           glclnd_present, glcocn_present, glcice_present, iceberg_prognostic,&
           esp_present, esp_prognostic,                                       &
           bfbflag, lnd_gnam, cpl_decomp, cpl_seq_option,                     &
           ice_gnam, rof_gnam, glc_gnam, wav_gnam,                            &
           atm_gnam, ocn_gnam, info_debug, dead_comps, read_restart,          &
           shr_map_dopole, vect_map, aoflux_grid, flux_epbalfact,             &
           nextsw_cday, precip_fact, flux_epbal, flux_albav, glcrun_alarm,    &
           glc_g2lupdate, atm_aero, run_barriers, esmf_map_flag,              &
           do_budgets, do_histinit, drv_threading, flux_diurnal,              &
           budget_inst, budget_daily, budget_month, wall_time_limit,          &
           budget_ann, budget_ltann, budget_ltend , force_stop_at,            &
           histaux_a2x    , histaux_a2x1hri, histaux_a2x1hr,                  &
           histaux_a2x3hr, histaux_a2x3hrp , histaux_l2x1yr,                  &
           histaux_a2x24hr, histaux_l2x   , histaux_r2x     , orb_obliq,      &
           histavg_atm, histavg_lnd, histavg_ocn, histavg_ice,                &
           histavg_rof, histavg_glc, histavg_wav, histavg_xao,                &
           cpl_cdf64, orb_iyear, orb_iyear_align, orb_mode, orb_mvelp,        &
           orb_eccen, orb_obliqr, orb_lambm0, orb_mvelpp, wv_sat_scheme,      &
           wv_sat_transition_start, wv_sat_use_tables, wv_sat_table_spacing,  &
           tfreeze_option,                                                       &
           glc_phase, rof_phase, atm_phase, lnd_phase, ocn_phase, ice_phase,  &
           wav_phase, esp_phase, wav_nx, wav_ny, atm_nx, atm_ny,              &
           lnd_nx, lnd_ny, rof_nx, rof_ny, ice_nx, ice_ny, ocn_nx, ocn_ny,    &
           glc_nx, glc_ny, eps_frac, eps_amask,                               &
           eps_agrid, eps_aarea, eps_omask, eps_ogrid, eps_oarea,             &
           reprosum_use_ddpdd, reprosum_diffmax, reprosum_recompute,          &
           mct_usealltoall, mct_usevector, max_cplstep_time,                  &
           pause_alarm, resume_alarm)


   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(seq_infodata_type),       intent(IN)  :: infodata      ! Input CCSM structure
   character(len=*),    optional, intent(OUT) :: cime_model    ! CIME model (acme or cesm)
   character(len=*),    optional, intent(OUT) :: start_type    ! Start type
   character(len=*),    optional, intent(OUT) :: case_name     ! Short case identification
   character(len=*),    optional, intent(OUT) :: case_desc     ! Long case description
   character(len=*),    optional, intent(OUT) :: model_version ! Model version
   character(len=*),    optional, intent(OUT) :: username      ! Username
   character(len=*),    optional, intent(OUT) :: hostname      ! Hostname
   character(len=*),    optional, intent(OUT) :: rest_case_name ! restart casename
   character(len=*),    optional, intent(OUT) :: timing_dir    ! timing dir name
   character(len=*),    optional, intent(OUT) :: tchkpt_dir    ! timing checkpoint dir name
   logical,             optional, intent(OUT) :: atm_adiabatic ! atm adiabatic mode
   logical,             optional, intent(OUT) :: atm_ideal_phys! atm idealized-physics mode
   logical,             optional, intent(OUT) :: aqua_planet   ! aqua_planet mode
   integer(SHR_KIND_IN),optional, intent(OUT) :: aqua_planet_sst! aqua_planet sst_type
   logical,             optional, intent(OUT) :: run_barriers  ! barrier component run calls
   logical,             optional, intent(OUT) :: brnch_retain_casename
   logical,             optional, intent(OUT) :: read_restart  ! read restart flag
   character(len=*),    optional, intent(OUT) :: restart_pfile ! Restart pointer file
   character(len=*),    optional, intent(OUT) :: restart_file  ! Restart file pathname
   logical,             optional, intent(OUT) :: single_column
   real (SHR_KIND_R8),  optional, intent(OUT) :: scmlat
   real (SHR_KIND_R8),  optional, intent(OUT) :: scmlon
   character(len=*),    optional, intent(OUT) :: logFilePostFix! output log file postfix
   character(len=*),    optional, intent(OUT) :: outPathRoot   ! output file root
   logical,             optional, intent(OUT) :: perpetual     ! If this is perpetual
   integer,             optional, intent(OUT) :: perpetual_ymd ! If perpetual, date
   character(len=*),    optional, intent(OUT) :: orb_mode      ! orbital mode
   integer,             optional, intent(OUT) :: orb_iyear     ! orbital year
   integer,             optional, intent(OUT) :: orb_iyear_align  ! orbital year model year align
   real(SHR_KIND_R8)   ,optional, intent(OUT) :: orb_eccen     ! See shr_orb_mod
   real(SHR_KIND_R8)   ,optional, intent(OUT) :: orb_obliqr    ! See shr_orb_mod
   real(SHR_KIND_R8)   ,optional, intent(OUT) :: orb_obliq     ! See shr_orb_mod
   real(SHR_KIND_R8)   ,optional, intent(OUT) :: orb_lambm0    ! See shr_orb_mod
   real(SHR_KIND_R8)   ,optional, intent(OUT) :: orb_mvelpp    ! See shr_orb_mod
   real(SHR_KIND_R8)   ,optional, intent(OUT) :: orb_mvelp     ! See shr_orb_mod
   character(len=*)    ,optional, intent(OUT) :: wv_sat_scheme ! Water vapor saturation pressure scheme
   real(SHR_KIND_R8)   ,optional, intent(OUT) :: wv_sat_transition_start   ! Saturation transition range
   logical             ,optional, intent(OUT) :: wv_sat_use_tables ! Saturation pressure lookup tables
   real(SHR_KIND_R8)   ,optional, intent(OUT) :: wv_sat_table_spacing  ! Saturation pressure table resolution
   character(len=*)    ,optional, intent(OUT) :: tfreeze_option   ! Freezing point of salt water
   character(len=*)    ,optional, intent(OUT) :: flux_epbal    ! selects E,P,R adjustment technique
   logical             ,optional, intent(OUT) :: flux_albav    ! T => no diurnal cycle in ocn albedos
   logical             ,optional, intent(OUT) :: flux_diurnal  ! T => diurnal cycle in atm/ocn flux
   real(SHR_KIND_R8)   ,optional, intent(OUT) :: wall_time_limit ! force stop wall time (hours)
   character(len=*)    ,optional, intent(OUT) :: force_stop_at ! force stop at next (month, day, etc)
   character(len=*)    ,optional, intent(OUT) :: atm_gnam      ! atm grid
   character(len=*)    ,optional, intent(OUT) :: lnd_gnam      ! lnd grid
   character(len=*)    ,optional, intent(OUT) :: ocn_gnam      ! ocn grid
   character(len=*)    ,optional, intent(OUT) :: ice_gnam      ! ice grid
   character(len=*)    ,optional, intent(OUT) :: rof_gnam      ! rof grid
   character(len=*)    ,optional, intent(OUT) :: glc_gnam      ! glc grid
   character(len=*)    ,optional, intent(OUT) :: wav_gnam      ! wav grid
   logical             ,optional, intent(OUT) :: shr_map_dopole  ! pole corrections in shr_map_mod
   character(len=*)    ,optional, intent(OUT) :: vect_map      ! vector mapping option
   character(len=*)    ,optional, intent(OUT) :: aoflux_grid   ! grid for atm ocn flux calc
   integer             ,optional, intent(OUT) :: cpl_decomp    ! coupler decomp
   character(len=*)    ,optional, intent(OUT) :: cpl_seq_option! coupler sequencing option
   logical             ,optional, intent(OUT) :: cpl_cdf64     ! netcdf large file setting
   logical             ,optional, intent(OUT) :: do_budgets    ! heat/water budgets
   logical             ,optional, intent(OUT) :: do_histinit   ! initial history file
   integer             ,optional, intent(OUT) :: budget_inst   ! inst budget
   integer             ,optional, intent(OUT) :: budget_daily  ! daily budget
   integer             ,optional, intent(OUT) :: budget_month  ! month budget
   integer             ,optional, intent(OUT) :: budget_ann    ! ann budget
   integer             ,optional, intent(OUT) :: budget_ltann  ! ltann budget
   integer             ,optional, intent(OUT) :: budget_ltend  ! ltend budget
   logical             ,optional, intent(OUT) :: histaux_a2x
   logical             ,optional, intent(OUT) :: histaux_a2x1hri
   logical             ,optional, intent(OUT) :: histaux_a2x1hr
   logical             ,optional, intent(OUT) :: histaux_a2x3hr
   logical             ,optional, intent(OUT) :: histaux_a2x3hrp
   logical             ,optional, intent(OUT) :: histaux_a2x24hr
   logical             ,optional, intent(OUT) :: histaux_l2x1yr
   logical             ,optional, intent(OUT) :: histaux_l2x
   logical             ,optional, intent(OUT) :: histaux_r2x
   logical             ,optional, intent(OUT) :: histavg_atm
   logical             ,optional, intent(OUT) :: histavg_lnd
   logical             ,optional, intent(OUT) :: histavg_ocn
   logical             ,optional, intent(OUT) :: histavg_ice
   logical             ,optional, intent(OUT) :: histavg_rof
   logical             ,optional, intent(OUT) :: histavg_glc
   logical             ,optional, intent(OUT) :: histavg_wav
   logical             ,optional, intent(OUT) :: histavg_xao
   logical             ,optional, intent(OUT) :: drv_threading ! driver threading control flag
   real(SHR_KIND_R8)   ,optional, intent(OUT) :: eps_frac      ! fraction error tolerance
   real(SHR_KIND_R8)   ,optional, intent(OUT) :: eps_amask     ! atm mask error tolerance
   real(SHR_KIND_R8)   ,optional, intent(OUT) :: eps_agrid     ! atm grid error tolerance
   real(SHR_KIND_R8)   ,optional, intent(OUT) :: eps_aarea     ! atm area error tolerance
   real(SHR_KIND_R8)   ,optional, intent(OUT) :: eps_omask     ! ocn mask error tolerance
   real(SHR_KIND_R8)   ,optional, intent(OUT) :: eps_ogrid     ! ocn grid error tolerance
   real(SHR_KIND_R8)   ,optional, intent(OUT) :: eps_oarea     ! ocn area error tolerance
   logical             ,optional, intent(OUT) :: reprosum_use_ddpdd ! use ddpdd algorithm
   real(SHR_KIND_R8)   ,optional, intent(OUT) :: reprosum_diffmax   ! maximum difference tolerance
   logical             ,optional, intent(OUT) :: reprosum_recompute ! recompute if tolerance exceeded
   logical             ,optional, intent(OUT) :: mct_usealltoall ! flag for mct alltoall
   logical             ,optional, intent(OUT) :: mct_usevector   ! flag for mct vector

   integer(SHR_KIND_IN),optional, intent(OUT) :: info_debug
   logical             ,optional, intent(OUT) :: bfbflag
   logical             ,optional, intent(OUT) :: esmf_map_flag
   logical             ,optional, intent(OUT) :: dead_comps    ! do we have dead models

   logical             ,optional, intent(OUT) :: atm_present    ! provide data
   logical             ,optional, intent(OUT) :: atm_prognostic ! need data
   logical             ,optional, intent(OUT) :: lnd_present
   logical             ,optional, intent(OUT) :: lnd_prognostic
   logical             ,optional, intent(OUT) :: rof_present
   logical             ,optional, intent(OUT) :: rofice_present
   logical             ,optional, intent(OUT) :: rof_prognostic
   logical             ,optional, intent(OUT) :: flood_present
   logical             ,optional, intent(OUT) :: ocn_present
   logical             ,optional, intent(OUT) :: ocn_prognostic
   logical             ,optional, intent(OUT) :: ocnrof_prognostic
   logical             ,optional, intent(OUT) :: ice_present
   logical             ,optional, intent(OUT) :: ice_prognostic
   logical             ,optional, intent(OUT) :: iceberg_prognostic
   logical             ,optional, intent(OUT) :: glc_present
   logical             ,optional, intent(OUT) :: glclnd_present
   logical             ,optional, intent(OUT) :: glcocn_present
   logical             ,optional, intent(OUT) :: glcice_present
   logical             ,optional, intent(OUT) :: glc_prognostic
   logical             ,optional, intent(OUT) :: wav_present
   logical             ,optional, intent(OUT) :: wav_prognostic
   logical             ,optional, intent(OUT) :: esp_present
   logical             ,optional, intent(OUT) :: esp_prognostic
   integer(SHR_KIND_IN),optional, intent(OUT) :: atm_nx        ! nx,ny 2d grid size global
   integer(SHR_KIND_IN),optional, intent(OUT) :: atm_ny        ! nx,ny 2d grid size global
   integer(SHR_KIND_IN),optional, intent(OUT) :: lnd_nx
   integer(SHR_KIND_IN),optional, intent(OUT) :: lnd_ny
   integer(SHR_KIND_IN),optional, intent(OUT) :: rof_nx
   integer(SHR_KIND_IN),optional, intent(OUT) :: rof_ny
   integer(SHR_KIND_IN),optional, intent(OUT) :: ice_nx
   integer(SHR_KIND_IN),optional, intent(OUT) :: ice_ny
   integer(SHR_KIND_IN),optional, intent(OUT) :: ocn_nx
   integer(SHR_KIND_IN),optional, intent(OUT) :: ocn_ny
   integer(SHR_KIND_IN),optional, intent(OUT) :: glc_nx
   integer(SHR_KIND_IN),optional, intent(OUT) :: glc_ny
   integer(SHR_KIND_IN),optional, intent(OUT) :: wav_nx
   integer(SHR_KIND_IN),optional, intent(OUT) :: wav_ny

   real(SHR_KIND_R8)   ,optional, intent(OUT) :: nextsw_cday   ! calendar of next atm shortwave
   real(SHR_KIND_R8)   ,optional, intent(OUT) :: precip_fact   ! precip factor
   real(SHR_KIND_R8)   ,optional, intent(OUT) :: flux_epbalfact ! adjusted precip factor
   integer(SHR_KIND_IN),optional, intent(OUT) :: atm_phase     ! atm phase
   integer(SHR_KIND_IN),optional, intent(OUT) :: lnd_phase     ! lnd phase
   integer(SHR_KIND_IN),optional, intent(OUT) :: ice_phase     ! ice phase
   integer(SHR_KIND_IN),optional, intent(OUT) :: ocn_phase     ! ocn phase
   integer(SHR_KIND_IN),optional, intent(OUT) :: glc_phase     ! glc phase
   integer(SHR_KIND_IN),optional, intent(OUT) :: rof_phase     ! rof phase
   integer(SHR_KIND_IN),optional, intent(OUT) :: wav_phase     ! wav phase
   integer(SHR_KIND_IN),optional, intent(OUT) :: esp_phase     ! wav phase
   logical             ,optional, intent(OUT) :: atm_aero      ! atmosphere aerosols
   logical             ,optional, intent(OUT) :: glcrun_alarm  ! glc run alarm
   logical             ,optional, intent(OUT) :: glc_g2lupdate ! update glc2lnd fields in lnd model
   real(shr_kind_r8), optional, intent(out) :: max_cplstep_time
   logical             ,optional, intent(OUT) :: pause_alarm
   logical             ,optional, intent(OUT) :: resume_alarm

    !----- local -----
    character(len=*), parameter :: subname = '(seq_infodata_GetData_explicit) '

!-------------------------------------------------------------------------------

    if ( present(cime_model)     ) cime_model     = infodata%cime_model
    if ( present(start_type)     ) start_type     = infodata%start_type
    if ( present(case_name)      ) case_name      = infodata%case_name
    if ( present(case_desc)      ) case_desc      = infodata%case_desc
    if ( present(model_version)  ) model_version  = infodata%model_version
    if ( present(username)       ) username       = infodata%username
    if ( present(hostname)       ) hostname       = infodata%hostname
    if ( present(rest_case_name) ) rest_case_name = infodata%rest_case_name
    if ( present(timing_dir)     ) timing_dir     = infodata%timing_dir
    if ( present(tchkpt_dir)     ) tchkpt_dir     = infodata%tchkpt_dir
    if ( present(atm_adiabatic)  ) atm_adiabatic  = infodata%atm_adiabatic
    if ( present(atm_ideal_phys) ) atm_ideal_phys = infodata%atm_ideal_phys
    if ( present(aqua_planet)    ) aqua_planet    = infodata%aqua_planet
    if ( present(aqua_planet_sst)) aqua_planet_sst= infodata%aqua_planet_sst
    if ( present(run_barriers)   ) run_barriers   = infodata%run_barriers
    if ( present(brnch_retain_casename) ) &
         brnch_retain_casename =  infodata%brnch_retain_casename
    if ( present(read_restart)   ) read_restart   = infodata%read_restart
    if ( present(restart_pfile)  ) restart_pfile  = infodata%restart_pfile
    if ( present(restart_file)   ) restart_file   = infodata%restart_file
    if ( present(single_column)  ) single_column  = infodata%single_column
    if ( present(scmlat)         ) scmlat         = infodata%scmlat
    if ( present(scmlon)         ) scmlon         = infodata%scmlon
    if ( present(logFilePostFix) ) logFilePostFix = infodata%logFilePostFix
    if ( present(outPathRoot)    ) outPathRoot    = infodata%outPathRoot
    if ( present(perpetual)      ) perpetual      = infodata%perpetual
    if ( present(perpetual_ymd)  ) perpetual_ymd  = infodata%perpetual_ymd
    if ( present(orb_iyear)      ) orb_iyear      = infodata%orb_iyear
    if ( present(orb_iyear_align)) orb_iyear_align= infodata%orb_iyear_align
    if ( present(orb_mode)       ) orb_mode       = infodata%orb_mode
    if ( present(orb_eccen)      ) orb_eccen      = infodata%orb_eccen
    if ( present(orb_obliqr)     ) orb_obliqr     = infodata%orb_obliqr
    if ( present(orb_obliq)      ) orb_obliq      = infodata%orb_obliq
    if ( present(orb_lambm0)     ) orb_lambm0     = infodata%orb_lambm0
    if ( present(orb_mvelpp)     ) orb_mvelpp     = infodata%orb_mvelpp
    if ( present(orb_mvelp)      ) orb_mvelp      = infodata%orb_mvelp
    if ( present(wv_sat_scheme)  ) wv_sat_scheme  = infodata%wv_sat_scheme
    if ( present(wv_sat_transition_start)) &
         wv_sat_transition_start = infodata%wv_sat_transition_start
    if ( present(wv_sat_use_tables)) wv_sat_use_tables = infodata%wv_sat_use_tables
    if ( present(wv_sat_table_spacing)) wv_sat_table_spacing = infodata%wv_sat_table_spacing
    if ( present(tfreeze_option) ) tfreeze_option = infodata%tfreeze_option
    if ( present(flux_epbal)     ) flux_epbal     = infodata%flux_epbal
    if ( present(flux_albav)     ) flux_albav     = infodata%flux_albav
    if ( present(flux_diurnal)   ) flux_diurnal   = infodata%flux_diurnal
    if ( present(wall_time_limit)) wall_time_limit= infodata%wall_time_limit
    if ( present(force_stop_at)  ) force_stop_at  = infodata%force_stop_at
    if ( present(atm_gnam)       ) atm_gnam       = infodata%atm_gnam
    if ( present(lnd_gnam)       ) lnd_gnam       = infodata%lnd_gnam
    if ( present(ocn_gnam)       ) ocn_gnam       = infodata%ocn_gnam
    if ( present(ice_gnam)       ) ice_gnam       = infodata%ice_gnam
    if ( present(rof_gnam)       ) rof_gnam       = infodata%rof_gnam
    if ( present(glc_gnam)       ) glc_gnam       = infodata%glc_gnam
    if ( present(wav_gnam)       ) wav_gnam       = infodata%wav_gnam
    if ( present(shr_map_dopole) ) shr_map_dopole = infodata%shr_map_dopole
    if ( present(vect_map)       ) vect_map       = infodata%vect_map
    if ( present(aoflux_grid)    ) aoflux_grid    = infodata%aoflux_grid
    if ( present(cpl_decomp)     ) cpl_decomp     = infodata%cpl_decomp
    if ( present(cpl_seq_option) ) cpl_seq_option = infodata%cpl_seq_option
    if ( present(cpl_cdf64)      ) cpl_cdf64      = infodata%cpl_cdf64
    if ( present(do_budgets)     ) do_budgets     = infodata%do_budgets
    if ( present(do_histinit)    ) do_histinit    = infodata%do_histinit
    if ( present(budget_inst)    ) budget_inst    = infodata%budget_inst
    if ( present(budget_daily)   ) budget_daily   = infodata%budget_daily
    if ( present(budget_month)   ) budget_month   = infodata%budget_month
    if ( present(budget_ann)     ) budget_ann     = infodata%budget_ann
    if ( present(budget_ltann)   ) budget_ltann   = infodata%budget_ltann
    if ( present(budget_ltend)   ) budget_ltend   = infodata%budget_ltend
    if ( present(histaux_a2x)    ) histaux_a2x    = infodata%histaux_a2x
    if ( present(histaux_a2x1hri)) histaux_a2x1hri= infodata%histaux_a2x1hri
    if ( present(histaux_a2x1hr) ) histaux_a2x1hr = infodata%histaux_a2x1hr
    if ( present(histaux_a2x3hr) ) histaux_a2x3hr = infodata%histaux_a2x3hr
    if ( present(histaux_a2x3hrp)) histaux_a2x3hrp= infodata%histaux_a2x3hrp
    if ( present(histaux_a2x24hr)) histaux_a2x24hr= infodata%histaux_a2x24hr
    if ( present(histaux_l2x1yr) ) histaux_l2x1yr = infodata%histaux_l2x1yr
    if ( present(histaux_l2x)    ) histaux_l2x    = infodata%histaux_l2x
    if ( present(histaux_r2x)    ) histaux_r2x    = infodata%histaux_r2x
    if ( present(histavg_atm)    ) histavg_atm    = infodata%histavg_atm
    if ( present(histavg_lnd)    ) histavg_lnd    = infodata%histavg_lnd
    if ( present(histavg_ocn)    ) histavg_ocn    = infodata%histavg_ocn
    if ( present(histavg_ice)    ) histavg_ice    = infodata%histavg_ice
    if ( present(histavg_rof)    ) histavg_rof    = infodata%histavg_rof
    if ( present(histavg_glc)    ) histavg_glc    = infodata%histavg_glc
    if ( present(histavg_wav)    ) histavg_wav    = infodata%histavg_wav
    if ( present(histavg_xao)    ) histavg_xao    = infodata%histavg_xao
    if ( present(drv_threading)  ) drv_threading  = infodata%drv_threading
    if ( present(eps_frac)       ) eps_frac       = infodata%eps_frac
    if ( present(eps_amask)      ) eps_amask      = infodata%eps_amask
    if ( present(eps_agrid)      ) eps_agrid      = infodata%eps_agrid
    if ( present(eps_aarea)      ) eps_aarea      = infodata%eps_aarea
    if ( present(eps_omask)      ) eps_omask      = infodata%eps_omask
    if ( present(eps_ogrid)      ) eps_ogrid      = infodata%eps_ogrid
    if ( present(eps_oarea)      ) eps_oarea      = infodata%eps_oarea
    if ( present(reprosum_use_ddpdd)) reprosum_use_ddpdd = infodata%reprosum_use_ddpdd
    if ( present(reprosum_diffmax)  ) reprosum_diffmax   = infodata%reprosum_diffmax
    if ( present(reprosum_recompute)) reprosum_recompute = infodata%reprosum_recompute
    if ( present(mct_usealltoall)) mct_usealltoall = infodata%mct_usealltoall
    if ( present(mct_usevector)  ) mct_usevector  = infodata%mct_usevector

    if ( present(info_debug)     ) info_debug     = infodata%info_debug
    if ( present(bfbflag)        ) bfbflag        = infodata%bfbflag
    if ( present(esmf_map_flag)  ) esmf_map_flag  = infodata%esmf_map_flag
    if ( present(dead_comps)     ) dead_comps     = infodata%dead_comps

    if ( present(atm_present)    ) atm_present    = infodata%atm_present
    if ( present(atm_prognostic) ) atm_prognostic = infodata%atm_prognostic
    if ( present(lnd_present)    ) lnd_present    = infodata%lnd_present
    if ( present(lnd_prognostic) ) lnd_prognostic = infodata%lnd_prognostic
    if ( present(rof_present)    ) rof_present    = infodata%rof_present
    if ( present(rofice_present) ) rofice_present = infodata%rofice_present
    if ( present(rof_prognostic) ) rof_prognostic = infodata%rof_prognostic
    if ( present(flood_present)  ) flood_present  = infodata%flood_present
    if ( present(ocn_present)    ) ocn_present    = infodata%ocn_present
    if ( present(ocn_prognostic) ) ocn_prognostic = infodata%ocn_prognostic
    if ( present(ocnrof_prognostic) ) ocnrof_prognostic = infodata%ocnrof_prognostic
    if ( present(ice_present)    ) ice_present    = infodata%ice_present
    if ( present(ice_prognostic) ) ice_prognostic = infodata%ice_prognostic
    if ( present(iceberg_prognostic)) iceberg_prognostic = infodata%iceberg_prognostic
    if ( present(glc_present)    ) glc_present    = infodata%glc_present
    if ( present(glclnd_present) ) glclnd_present = infodata%glclnd_present
    if ( present(glcocn_present) ) glcocn_present = infodata%glcocn_present
    if ( present(glcice_present) ) glcice_present = infodata%glcice_present
    if ( present(glc_prognostic) ) glc_prognostic = infodata%glc_prognostic
    if ( present(wav_present)    ) wav_present    = infodata%wav_present
    if ( present(wav_prognostic) ) wav_prognostic = infodata%wav_prognostic
    if ( present(esp_present)    ) esp_present    = infodata%esp_present
    if ( present(esp_prognostic) ) esp_prognostic = infodata%esp_prognostic
    if ( present(atm_nx)         ) atm_nx         = infodata%atm_nx
    if ( present(atm_ny)         ) atm_ny         = infodata%atm_ny
    if ( present(lnd_nx)         ) lnd_nx         = infodata%lnd_nx
    if ( present(lnd_ny)         ) lnd_ny         = infodata%lnd_ny
    if ( present(rof_nx)         ) rof_nx         = infodata%rof_nx
    if ( present(rof_ny)         ) rof_ny         = infodata%rof_ny
    if ( present(ice_nx)         ) ice_nx         = infodata%ice_nx
    if ( present(ice_ny)         ) ice_ny         = infodata%ice_ny
    if ( present(ocn_nx)         ) ocn_nx         = infodata%ocn_nx
    if ( present(ocn_ny)         ) ocn_ny         = infodata%ocn_ny
    if ( present(glc_nx)         ) glc_nx         = infodata%glc_nx
    if ( present(glc_ny)         ) glc_ny         = infodata%glc_ny
    if ( present(wav_nx)         ) wav_nx         = infodata%wav_nx
    if ( present(wav_ny)         ) wav_ny         = infodata%wav_ny

    if ( present(nextsw_cday)    ) nextsw_cday    = infodata%nextsw_cday
    if ( present(precip_fact)    ) precip_fact    = infodata%precip_fact
    if ( present(flux_epbalfact) ) then
       flux_epbalfact = 1.0_SHR_KIND_R8
       if (trim(infodata%flux_epbal) == 'ocn') then
          flux_epbalfact = infodata%precip_fact
       end if
       if (flux_epbalfact <= 0.0_SHR_KIND_R8) then
          if (loglevel > 0) write(logunit,'(2a,e16.6)') &
             trim(subname),' WARNING: factor from ocn = ',flux_epbalfact
          if (loglevel > 0) write(logunit,'(2a)') &
             trim(subname),' WARNING: resetting flux_epbalfact to 1.0'
          flux_epbalfact = 1.0_SHR_KIND_R8
       end if
    endif
    if ( present(atm_phase)      ) atm_phase      = infodata%atm_phase
    if ( present(lnd_phase)      ) lnd_phase      = infodata%lnd_phase
    if ( present(ice_phase)      ) ice_phase      = infodata%ice_phase
    if ( present(ocn_phase)      ) ocn_phase      = infodata%ocn_phase
    if ( present(glc_phase)      ) glc_phase      = infodata%glc_phase
    if ( present(rof_phase)      ) rof_phase      = infodata%rof_phase
    if ( present(wav_phase)      ) wav_phase      = infodata%wav_phase
    if ( present(esp_phase)      ) esp_phase      = infodata%esp_phase
    if ( present(atm_aero)       ) atm_aero       = infodata%atm_aero
    if ( present(glcrun_alarm)   ) glcrun_alarm   = infodata%glcrun_alarm
    if ( present(glc_g2lupdate)  ) glc_g2lupdate  = infodata%glc_g2lupdate
    if ( present(max_cplstep_time) ) max_cplstep_time = infodata%max_cplstep_time
    if ( present(pause_alarm)  ) pause_alarm  = infodata%pause_alarm
    if ( present(resume_alarm)  ) resume_alarm  = infodata%resume_alarm

END SUBROUTINE seq_infodata_GetData_explicit

#ifndef CPRPGI
!===============================================================================
! !IROUTINE: seq_infodata_GetData_bytype -- Get values from infodata object
!
! !DESCRIPTION:
!
!     Get values out of the infodata object.
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE seq_infodata_GetData_bytype( component_firstletter, infodata,      &
           comp_present, comp_prognostic, comp_gnam,                          &
           histavg_comp, comp_phase, comp_nx, comp_ny)


   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(len=1),              intent(IN)  :: component_firstletter
   type(seq_infodata_type),       intent(IN)  :: infodata      ! Input CCSM structure
   logical             ,optional, intent(OUT) :: comp_present    ! provide data
   logical             ,optional, intent(OUT) :: comp_prognostic ! need data
   character(len=*)    ,optional, intent(OUT) :: comp_gnam      ! comp grid
   integer(SHR_KIND_IN),optional, intent(OUT) :: comp_nx        ! nx,ny 2d grid size global
   integer(SHR_KIND_IN),optional, intent(OUT) :: comp_ny        ! nx,ny 2d grid size global
   integer(SHR_KIND_IN),optional, intent(OUT) :: comp_phase
   logical             ,optional, intent(OUT) :: histavg_comp

    !----- local -----
    character(len=*), parameter :: subname = '(seq_infodata_GetData_bytype) '

!-------------------------------------------------------------------------------

    if (component_firstletter == 'a') then
      call seq_infodata_GetData(infodata, atm_present=comp_present,           &
           atm_prognostic=comp_prognostic, atm_gnam=comp_gnam,                &
           atm_phase=comp_phase, atm_nx=comp_nx, atm_ny=comp_ny, histavg_atm=histavg_comp)
    else if (component_firstletter == 'l') then
      call seq_infodata_GetData(infodata, lnd_present=comp_present,           &
           lnd_prognostic=comp_prognostic, lnd_gnam=comp_gnam,                &
           lnd_phase=comp_phase, lnd_nx=comp_nx, lnd_ny=comp_ny, histavg_lnd=histavg_comp)
    else if (component_firstletter == 'i') then
      call seq_infodata_GetData(infodata, ice_present=comp_present,           &
           ice_prognostic=comp_prognostic, ice_gnam=comp_gnam,                &
           ice_phase=comp_phase, ice_nx=comp_nx, ice_ny=comp_ny, histavg_ice=histavg_comp)
    else if (component_firstletter == 'o') then
      call seq_infodata_GetData(infodata, ocn_present=comp_present,           &
           ocn_prognostic=comp_prognostic, ocn_gnam=comp_gnam,                &
           ocn_phase=comp_phase, ocn_nx=comp_nx, ocn_ny=comp_ny, histavg_ocn=histavg_comp)
    else if (component_firstletter == 'r') then
      call seq_infodata_GetData(infodata, rof_present=comp_present,           &
           rof_prognostic=comp_prognostic, rof_gnam=comp_gnam,                &
           rof_phase=comp_phase, rof_nx=comp_nx, rof_ny=comp_ny, histavg_rof=histavg_comp)
    else if (component_firstletter == 'g') then
      call seq_infodata_GetData(infodata, glc_present=comp_present,           &
           glc_prognostic=comp_prognostic, glc_gnam=comp_gnam,                &
           glc_phase=comp_phase, glc_nx=comp_nx, glc_ny=comp_ny, histavg_glc=histavg_comp)
    else if (component_firstletter == 'w') then
      call seq_infodata_GetData(infodata, wav_present=comp_present,           &
           wav_prognostic=comp_prognostic, wav_gnam=comp_gnam,                &
           wav_phase=comp_phase, wav_nx=comp_nx, wav_ny=comp_ny, histavg_wav=histavg_comp)
    else if (component_firstletter == 'e') then
      if (present(comp_gnam)) then
        comp_gnam = ''
        if ((loglevel > 1) .and. seq_comm_iamroot(1)) then
          write(logunit,*) trim(subname),' Note: ESP type has no gnam property'
        end if
      end if
      if (present(comp_nx)) then
        comp_nx = 1
        if ((loglevel > 1) .and. seq_comm_iamroot(1)) then
          write(logunit,*) trim(subname),' Note: ESP type has no nx property'
        end if
      end if
      if (present(comp_ny)) then
        comp_ny = 1
        if ((loglevel > 1) .and. seq_comm_iamroot(1)) then
          write(logunit,*) trim(subname),' Note: ESP type has no ny property'
        end if
      end if
      if (present(histavg_comp)) then
        histavg_comp = .false.
        if ((loglevel > 1) .and. seq_comm_iamroot(1)) then
          write(logunit,*) trim(subname),' Note: ESP type has no histavg property'
        end if
      end if
     
      call seq_infodata_GetData(infodata, esp_present=comp_present,           &
           esp_prognostic=comp_prognostic, esp_phase=comp_phase)
    else
       call shr_sys_abort( subname//": unknown component-type first letter,'"//component_firstletter//"', aborting")
    end if
  END SUBROUTINE seq_infodata_GetData_bytype
#endif
! ^ ifndef CPRPGI

!===============================================================================
! !IROUTINE: seq_infodata_PutData_explicit -- Put values into infodata object
!
! !DESCRIPTION:
!
!     Put values into the infodata object.
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE seq_infodata_PutData_explicit( infodata, cime_model, case_name, case_desc, timing_dir,  &
           model_version, username, hostname, rest_case_name, tchkpt_dir,     &
           start_type, restart_pfile, restart_file, perpetual, perpetual_ymd, &
           aqua_planet,aqua_planet_sst, atm_ideal_phys, atm_adiabatic, brnch_retain_casename, &
           single_column, scmlat,scmlon,logFilePostFix, outPathRoot,          &
           atm_present, atm_prognostic, lnd_present, lnd_prognostic, rof_prognostic, &
           rof_present, ocn_present, ocn_prognostic, ocnrof_prognostic,       &
           ice_present, ice_prognostic, glc_present, glc_prognostic,          &
           flood_present, wav_present, wav_prognostic, rofice_present,        &
           glclnd_present, glcocn_present, glcice_present, iceberg_prognostic,&
           esp_present, esp_prognostic,                                       &
           bfbflag, lnd_gnam, cpl_decomp, cpl_seq_option,                     &
           ice_gnam, rof_gnam, glc_gnam, wav_gnam,                            &
           atm_gnam, ocn_gnam, info_debug, dead_comps, read_restart,          &
           shr_map_dopole, vect_map, aoflux_grid, run_barriers,               &
           nextsw_cday, precip_fact, flux_epbal, flux_albav, glcrun_alarm,    &
           glc_g2lupdate, atm_aero, esmf_map_flag, wall_time_limit,           &
           do_budgets, do_histinit, drv_threading, flux_diurnal,              &
           budget_inst, budget_daily, budget_month, force_stop_at,            &
           budget_ann, budget_ltann, budget_ltend ,                           &
           histaux_a2x    , histaux_a2x1hri, histaux_a2x1hr,                  &
           histaux_a2x3hr, histaux_a2x3hrp , histaux_l2x1yr,                  &
           histaux_a2x24hr, histaux_l2x   , histaux_r2x     , orb_obliq,      &
           histavg_atm, histavg_lnd, histavg_ocn, histavg_ice,                &
           histavg_rof, histavg_glc, histavg_wav, histavg_xao,                &
           cpl_cdf64, orb_iyear, orb_iyear_align, orb_mode, orb_mvelp,        &
           orb_eccen, orb_obliqr, orb_lambm0, orb_mvelpp, wv_sat_scheme,      &
           wv_sat_transition_start, wv_sat_use_tables, wv_sat_table_spacing,  &
           tfreeze_option, &
           glc_phase, rof_phase, atm_phase, lnd_phase, ocn_phase, ice_phase,  &
           wav_phase, esp_phase, wav_nx, wav_ny, atm_nx, atm_ny,              &
           lnd_nx, lnd_ny, rof_nx, rof_ny, ice_nx, ice_ny, ocn_nx, ocn_ny,    &
           glc_nx, glc_ny, eps_frac, eps_amask,                               &
           eps_agrid, eps_aarea, eps_omask, eps_ogrid, eps_oarea,             &
           reprosum_use_ddpdd, reprosum_diffmax, reprosum_recompute,          &
           mct_usealltoall, mct_usevector, pause_alarm, resume_alarm )


   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(seq_infodata_type),    intent(INOUT) :: infodata      ! Input CCSM structure
   character(len=*),    optional, intent(IN) :: cime_model    ! CIME model (acme or cesm)
   character(len=*),    optional, intent(IN) :: start_type    ! Start type
   character(len=*),    optional, intent(IN) :: case_name     ! Short case identification
   character(len=*),    optional, intent(IN) :: case_desc     ! Long case description
   character(len=*),    optional, intent(IN) :: model_version ! Model version
   character(len=*),    optional, intent(IN) :: username      ! Username
   character(len=*),    optional, intent(IN) :: hostname      ! Hostname
   character(len=*),    optional, intent(IN) :: rest_case_name ! restart casename
   character(len=*),    optional, intent(IN) :: timing_dir    ! timing dir name
   character(len=*),    optional, intent(IN) :: tchkpt_dir    ! timing checkpoint dir name
   logical,             optional, intent(IN) :: atm_adiabatic ! atm adiabatic mode
   logical,             optional, intent(IN) :: atm_ideal_phys! atm idealized-physics mode
   logical,             optional, intent(IN) :: aqua_planet   ! aqua_planet mode
   integer(SHR_KIND_IN),optional, intent(IN) :: aqua_planet_sst ! aqua_planet sst type
   logical,             optional, intent(IN) :: run_barriers  ! barrier component run calls
   logical,             optional, intent(IN) :: brnch_retain_casename
   logical,             optional, intent(IN) :: read_restart  ! read restart flag
   character(len=*),    optional, intent(IN) :: restart_pfile ! Restart pointer file
   character(len=*),    optional, intent(IN) :: restart_file  ! Restart file pathname
   logical,             optional, intent(IN) :: single_column
   real (SHR_KIND_R8),  optional, intent(IN) :: scmlat
   real (SHR_KIND_R8),  optional, intent(IN) :: scmlon
   character(len=*),    optional, intent(IN) :: logFilePostFix! output log file postfix
   character(len=*),    optional, intent(IN) :: outPathRoot   ! output file root
   logical,             optional, intent(IN) :: perpetual     ! If this is perpetual
   integer,             optional, intent(IN) :: perpetual_ymd ! If perpetual, date
   character(len=*),    optional, intent(IN) :: orb_mode      ! orbital mode
   integer,             optional, intent(IN) :: orb_iyear     ! orbital year
   integer,             optional, intent(IN) :: orb_iyear_align  ! orbital year model year align
   real(SHR_KIND_R8)   ,optional, intent(IN) :: orb_eccen     ! See shr_orb_mod
   real(SHR_KIND_R8)   ,optional, intent(IN) :: orb_obliqr    ! See shr_orb_mod
   real(SHR_KIND_R8)   ,optional, intent(IN) :: orb_obliq     ! See shr_orb_mod
   real(SHR_KIND_R8)   ,optional, intent(IN) :: orb_lambm0    ! See shr_orb_mod
   real(SHR_KIND_R8)   ,optional, intent(IN) :: orb_mvelpp    ! See shr_orb_mod
   real(SHR_KIND_R8)   ,optional, intent(IN) :: orb_mvelp     ! See shr_orb_mod
   character(len=*)    ,optional, intent(IN) :: wv_sat_scheme ! Water vapor saturation pressure scheme
   real(SHR_KIND_R8)   ,optional, intent(IN) :: wv_sat_transition_start  ! Saturation transition range
   logical             ,optional, intent(IN) :: wv_sat_use_tables ! Saturation pressure lookup tables
   real(SHR_KIND_R8)   ,optional, intent(IN) :: wv_sat_table_spacing  ! Saturation pressure table resolution
   character(len=*)    ,optional, intent(IN) :: tfreeze_option   ! Freezing point of salt water
   character(len=*)    ,optional, intent(IN) :: flux_epbal    ! selects E,P,R adjustment technique
   logical             ,optional, intent(IN) :: flux_albav    ! T => no diurnal cycle in ocn albedos
   logical             ,optional, intent(IN) :: flux_diurnal  ! T => diurnal cycle in atm/ocn flux
   real(SHR_KIND_R8)   ,optional, intent(IN) :: wall_time_limit ! force stop wall time (hours)
   character(len=*)    ,optional, intent(IN) :: force_stop_at ! force a stop at next (month, day, etc)
   character(len=*)    ,optional, intent(IN) :: atm_gnam   ! atm grid
   character(len=*)    ,optional, intent(IN) :: lnd_gnam   ! lnd grid
   character(len=*)    ,optional, intent(IN) :: ocn_gnam   ! ocn grid
   character(len=*)    ,optional, intent(IN) :: ice_gnam   ! ice grid
   character(len=*)    ,optional, intent(IN) :: rof_gnam   ! rof grid
   character(len=*)    ,optional, intent(IN) :: glc_gnam   ! glc grid
   character(len=*)    ,optional, intent(IN) :: wav_gnam   ! wav grid
   logical             ,optional, intent(IN) :: shr_map_dopole  ! pole corrections in shr_map_mod
   character(len=*)    ,optional, intent(IN) :: vect_map      ! vector mapping option
   character(len=*)    ,optional, intent(IN) :: aoflux_grid   ! grid for atm ocn flux calc
   integer             ,optional, intent(IN) :: cpl_decomp    ! coupler decomp
   character(len=*)    ,optional, intent(IN) :: cpl_seq_option! coupler sequencing option
   logical             ,optional, intent(IN) :: cpl_cdf64     ! netcdf large file setting
   logical             ,optional, intent(IN) :: do_budgets    ! heat/water budgets
   logical             ,optional, intent(IN) :: do_histinit   ! initial history file
   integer             ,optional, intent(IN) :: budget_inst   ! inst budget
   integer             ,optional, intent(IN) :: budget_daily  ! daily budget
   integer             ,optional, intent(IN) :: budget_month  ! month budget
   integer             ,optional, intent(IN) :: budget_ann    ! ann budget
   integer             ,optional, intent(IN) :: budget_ltann  ! ltann budget
   integer             ,optional, intent(IN) :: budget_ltend  ! ltend budget
   logical             ,optional, intent(IN) :: histaux_a2x
   logical             ,optional, intent(IN) :: histaux_a2x1hri
   logical             ,optional, intent(IN) :: histaux_a2x1hr
   logical             ,optional, intent(IN) :: histaux_a2x3hr
   logical             ,optional, intent(IN) :: histaux_a2x3hrp
   logical             ,optional, intent(IN) :: histaux_a2x24hr
   logical             ,optional, intent(IN) :: histaux_l2x1yr
   logical             ,optional, intent(IN) :: histaux_l2x
   logical             ,optional, intent(IN) :: histaux_r2x
   logical             ,optional, intent(IN) :: histavg_atm
   logical             ,optional, intent(IN) :: histavg_lnd
   logical             ,optional, intent(IN) :: histavg_ocn
   logical             ,optional, intent(IN) :: histavg_ice
   logical             ,optional, intent(IN) :: histavg_rof
   logical             ,optional, intent(IN) :: histavg_glc
   logical             ,optional, intent(IN) :: histavg_wav
   logical             ,optional, intent(IN) :: histavg_xao
   logical             ,optional, intent(IN) :: drv_threading ! driver threading control flag
   real(SHR_KIND_R8)   ,optional, intent(IN) :: eps_frac      ! fraction error tolerance
   real(SHR_KIND_R8)   ,optional, intent(IN) :: eps_amask     ! atm mask error tolerance
   real(SHR_KIND_R8)   ,optional, intent(IN) :: eps_agrid     ! atm grid error tolerance
   real(SHR_KIND_R8)   ,optional, intent(IN) :: eps_aarea     ! atm area error tolerance
   real(SHR_KIND_R8)   ,optional, intent(IN) :: eps_omask     ! ocn mask error tolerance
   real(SHR_KIND_R8)   ,optional, intent(IN) :: eps_ogrid     ! ocn grid error tolerance
   real(SHR_KIND_R8)   ,optional, intent(IN) :: eps_oarea     ! ocn area error tolerance
   logical             ,optional, intent(IN) :: reprosum_use_ddpdd ! use ddpdd algorithm
   real(SHR_KIND_R8)   ,optional, intent(IN) :: reprosum_diffmax   ! maximum difference tolerance
   logical             ,optional, intent(IN) :: reprosum_recompute ! recompute if tolerance exceeded
   logical             ,optional, intent(IN) :: mct_usealltoall ! flag for mct alltoall
   logical             ,optional, intent(IN) :: mct_usevector   ! flag for mct vector

   integer(SHR_KIND_IN),optional, intent(IN) :: info_debug
   logical             ,optional, intent(IN) :: bfbflag
   logical             ,optional, intent(IN) :: esmf_map_flag
   logical             ,optional, intent(IN) :: dead_comps    ! do we have dead models

   logical             ,optional, intent(IN) :: atm_present    ! provide data
   logical             ,optional, intent(IN) :: atm_prognostic ! need data
   logical             ,optional, intent(IN) :: lnd_present
   logical             ,optional, intent(IN) :: lnd_prognostic
   logical             ,optional, intent(IN) :: rof_present
   logical             ,optional, intent(IN) :: rofice_present
   logical             ,optional, intent(IN) :: rof_prognostic
   logical             ,optional, intent(IN) :: flood_present
   logical             ,optional, intent(IN) :: ocn_present
   logical             ,optional, intent(IN) :: ocn_prognostic
   logical             ,optional, intent(IN) :: ocnrof_prognostic
   logical             ,optional, intent(IN) :: ice_present
   logical             ,optional, intent(IN) :: ice_prognostic
   logical             ,optional, intent(IN) :: iceberg_prognostic
   logical             ,optional, intent(IN) :: glc_present
   logical             ,optional, intent(IN) :: glclnd_present
   logical             ,optional, intent(IN) :: glcocn_present
   logical             ,optional, intent(IN) :: glcice_present
   logical             ,optional, intent(IN) :: glc_prognostic
   logical             ,optional, intent(IN) :: wav_present
   logical             ,optional, intent(IN) :: wav_prognostic
   logical             ,optional, intent(IN) :: esp_present
   logical             ,optional, intent(IN) :: esp_prognostic
   integer(SHR_KIND_IN),optional, intent(IN) :: atm_nx        ! nx,ny 2d grid size global
   integer(SHR_KIND_IN),optional, intent(IN) :: atm_ny        ! nx,ny 2d grid size global
   integer(SHR_KIND_IN),optional, intent(IN) :: lnd_nx
   integer(SHR_KIND_IN),optional, intent(IN) :: lnd_ny
   integer(SHR_KIND_IN),optional, intent(IN) :: rof_nx
   integer(SHR_KIND_IN),optional, intent(IN) :: rof_ny
   integer(SHR_KIND_IN),optional, intent(IN) :: ice_nx
   integer(SHR_KIND_IN),optional, intent(IN) :: ice_ny
   integer(SHR_KIND_IN),optional, intent(IN) :: ocn_nx
   integer(SHR_KIND_IN),optional, intent(IN) :: ocn_ny
   integer(SHR_KIND_IN),optional, intent(IN) :: glc_nx
   integer(SHR_KIND_IN),optional, intent(IN) :: glc_ny
   integer(SHR_KIND_IN),optional, intent(IN) :: wav_nx
   integer(SHR_KIND_IN),optional, intent(IN) :: wav_ny

   real(SHR_KIND_R8)   ,optional, intent(IN) :: nextsw_cday   ! calendar of next atm shortwave
   real(SHR_KIND_R8)   ,optional, intent(IN) :: precip_fact   ! precip factor
   integer(SHR_KIND_IN),optional, intent(IN) :: atm_phase     ! atm phase
   integer(SHR_KIND_IN),optional, intent(IN) :: lnd_phase     ! lnd phase
   integer(SHR_KIND_IN),optional, intent(IN) :: ice_phase     ! ice phase
   integer(SHR_KIND_IN),optional, intent(IN) :: ocn_phase     ! ocn phase
   integer(SHR_KIND_IN),optional, intent(IN) :: glc_phase     ! glc phase
   integer(SHR_KIND_IN),optional, intent(IN) :: rof_phase     ! glc phase
   integer(SHR_KIND_IN),optional, intent(IN) :: wav_phase     ! wav phase
   integer(SHR_KIND_IN),optional, intent(IN) :: esp_phase     ! esp phase
   logical             ,optional, intent(IN) :: atm_aero      ! atm aerosols
   logical             ,optional, intent(IN) :: glcrun_alarm  ! glc run alarm
   logical             ,optional, intent(IN) :: glc_g2lupdate ! update glc2lnd fields in lnd model
   logical             ,optional, intent(IN) :: pause_alarm
   logical             ,optional, intent(IN) :: resume_alarm

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(seq_infodata_PutData_explicit) '

!-------------------------------------------------------------------------------

    if ( present(cime_model)     ) infodata%cime_model     = cime_model
    if ( present(start_type)     ) infodata%start_type     = start_type
    if ( present(case_name)      ) infodata%case_name      = case_name
    if ( present(case_desc)      ) infodata%case_desc      = case_desc
    if ( present(model_version)  ) infodata%model_version  = model_version
    if ( present(username)       ) infodata%username       = username
    if ( present(hostname)       ) infodata%hostname       = hostname
    if ( present(rest_case_name) ) infodata%rest_case_name = rest_case_name
    if ( present(timing_dir)     ) infodata%timing_dir     = timing_dir
    if ( present(tchkpt_dir)     ) infodata%tchkpt_dir     = tchkpt_dir
    if ( present(atm_adiabatic)  ) infodata%atm_adiabatic  = atm_adiabatic
    if ( present(atm_ideal_phys) ) infodata%atm_ideal_phys = atm_ideal_phys
    if ( present(aqua_planet)    ) infodata%aqua_planet    = aqua_planet
    if ( present(aqua_planet_sst)) infodata%aqua_planet_sst= aqua_planet_sst
    if ( present(run_barriers)   ) infodata%run_barriers   = run_barriers
    if ( present(brnch_retain_casename)) infodata%brnch_retain_casename = brnch_retain_casename
    if ( present(read_restart)   ) infodata%read_restart   = read_restart
    if ( present(restart_pfile)  ) infodata%restart_pfile  = restart_pfile
    if ( present(restart_file)   ) infodata%restart_file   = restart_file
    if ( present(single_column)  ) infodata%single_column  = single_column
    if ( present(scmlat)         ) infodata%scmlat         = scmlat
    if ( present(scmlon)         ) infodata%scmlon         = scmlon
    if ( present(logFilePostFix) ) infodata%logFilePostFix = logFilePostFix
    if ( present(outPathRoot)    ) infodata%outPathRoot    = outPathRoot
    if ( present(perpetual)      ) infodata%perpetual      = perpetual
    if ( present(perpetual_ymd)  ) infodata%perpetual_ymd  = perpetual_ymd
    if ( present(orb_iyear)      ) infodata%orb_iyear      = orb_iyear
    if ( present(orb_iyear_align)) infodata%orb_iyear_align= orb_iyear_align
    if ( present(orb_mode)       ) infodata%orb_mode       = orb_mode
    if ( present(orb_eccen)      ) infodata%orb_eccen      = orb_eccen
    if ( present(orb_obliqr)     ) infodata%orb_obliqr     = orb_obliqr
    if ( present(orb_obliq)      ) infodata%orb_obliq      = orb_obliq
    if ( present(orb_lambm0)     ) infodata%orb_lambm0     = orb_lambm0
    if ( present(orb_mvelpp)     ) infodata%orb_mvelpp     = orb_mvelpp
    if ( present(orb_mvelp)      ) infodata%orb_mvelp      = orb_mvelp
    if ( present(wv_sat_scheme)  ) infodata%wv_sat_scheme  = wv_sat_scheme
    if ( present(wv_sat_transition_start)) &
         infodata%wv_sat_transition_start = wv_sat_transition_start
    if ( present(wv_sat_use_tables)) infodata%wv_sat_use_tables = wv_sat_use_tables
    if ( present(wv_sat_table_spacing)) infodata%wv_sat_table_spacing = wv_sat_table_spacing
    if ( present(tfreeze_option)    ) infodata%tfreeze_option    = tfreeze_option
    if ( present(flux_epbal)     ) infodata%flux_epbal     = flux_epbal
    if ( present(flux_albav)     ) infodata%flux_albav     = flux_albav
    if ( present(flux_diurnal)   ) infodata%flux_diurnal   = flux_diurnal
    if ( present(wall_time_limit)) infodata%wall_time_limit= wall_time_limit
    if ( present(force_stop_at)  ) infodata%force_stop_at  = force_stop_at
    if ( present(atm_gnam)       ) infodata%atm_gnam       = atm_gnam
    if ( present(lnd_gnam)       ) infodata%lnd_gnam       = lnd_gnam
    if ( present(ocn_gnam)       ) infodata%ocn_gnam       = ocn_gnam
    if ( present(ice_gnam)       ) infodata%ice_gnam       = ice_gnam
    if ( present(rof_gnam)       ) infodata%rof_gnam       = rof_gnam
    if ( present(glc_gnam)       ) infodata%glc_gnam       = glc_gnam
    if ( present(wav_gnam)       ) infodata%wav_gnam       = wav_gnam
    if ( present(shr_map_dopole) ) infodata%shr_map_dopole = shr_map_dopole
    if ( present(vect_map)       ) infodata%vect_map       = vect_map
    if ( present(aoflux_grid)    ) infodata%aoflux_grid    = aoflux_grid
    if ( present(cpl_decomp)     ) infodata%cpl_decomp     = cpl_decomp
    if ( present(cpl_seq_option) ) infodata%cpl_seq_option = cpl_seq_option
    if ( present(cpl_cdf64)      ) infodata%cpl_cdf64      = cpl_cdf64
    if ( present(do_budgets)     ) infodata%do_budgets     = do_budgets
    if ( present(do_histinit)    ) infodata%do_histinit    = do_histinit
    if ( present(budget_inst)    ) infodata%budget_inst    = budget_inst
    if ( present(budget_daily)   ) infodata%budget_daily   = budget_daily
    if ( present(budget_month)   ) infodata%budget_month   = budget_month
    if ( present(budget_ann)     ) infodata%budget_ann     = budget_ann
    if ( present(budget_ltann)   ) infodata%budget_ltann   = budget_ltann
    if ( present(budget_ltend)   ) infodata%budget_ltend   = budget_ltend
    if ( present(histaux_a2x)    ) infodata%histaux_a2x    = histaux_a2x
    if ( present(histaux_a2x1hri)) infodata%histaux_a2x1hri= histaux_a2x1hri
    if ( present(histaux_a2x1hr) ) infodata%histaux_a2x1hr = histaux_a2x1hr
    if ( present(histaux_a2x3hr) ) infodata%histaux_a2x3hr = histaux_a2x3hr
    if ( present(histaux_a2x3hrp)) infodata%histaux_a2x3hrp= histaux_a2x3hrp
    if ( present(histaux_a2x24hr)) infodata%histaux_a2x24hr= histaux_a2x24hr
    if ( present(histaux_l2x1yr) ) infodata%histaux_l2x1yr = histaux_l2x1yr
    if ( present(histaux_l2x)    ) infodata%histaux_l2x    = histaux_l2x
    if ( present(histaux_r2x)    ) infodata%histaux_r2x    = histaux_r2x
    if ( present(histavg_atm)    ) infodata%histavg_atm    = histavg_atm
    if ( present(histavg_lnd)    ) infodata%histavg_lnd    = histavg_lnd
    if ( present(histavg_ocn)    ) infodata%histavg_ocn    = histavg_ocn
    if ( present(histavg_ice)    ) infodata%histavg_ice    = histavg_ice
    if ( present(histavg_rof)    ) infodata%histavg_rof    = histavg_rof
    if ( present(histavg_glc)    ) infodata%histavg_glc    = histavg_glc
    if ( present(histavg_wav)    ) infodata%histavg_wav    = histavg_wav
    if ( present(histavg_xao)    ) infodata%histavg_xao    = histavg_xao
    if ( present(drv_threading)  ) infodata%drv_threading  = drv_threading
    if ( present(eps_frac)       ) infodata%eps_frac       = eps_frac
    if ( present(eps_amask)      ) infodata%eps_amask      = eps_amask
    if ( present(eps_agrid)      ) infodata%eps_agrid      = eps_agrid
    if ( present(eps_aarea)      ) infodata%eps_aarea      = eps_aarea
    if ( present(eps_omask)      ) infodata%eps_omask      = eps_omask
    if ( present(eps_ogrid)      ) infodata%eps_ogrid      = eps_ogrid
    if ( present(eps_oarea)      ) infodata%eps_oarea      = eps_oarea
    if ( present(reprosum_use_ddpdd)) infodata%reprosum_use_ddpdd = reprosum_use_ddpdd
    if ( present(reprosum_diffmax)  ) infodata%reprosum_diffmax   = reprosum_diffmax
    if ( present(reprosum_recompute)) infodata%reprosum_recompute = reprosum_recompute
    if ( present(mct_usealltoall)) infodata%mct_usealltoall = mct_usealltoall
    if ( present(mct_usevector)  ) infodata%mct_usevector  = mct_usevector

    if ( present(info_debug)     ) infodata%info_debug     = info_debug
    if ( present(bfbflag)        ) infodata%bfbflag        = bfbflag
    if ( present(esmf_map_flag)  ) infodata%esmf_map_flag  = esmf_map_flag
    if ( present(dead_comps)     ) infodata%dead_comps     = dead_comps

    if ( present(atm_present)    ) infodata%atm_present    = atm_present
    if ( present(atm_prognostic) ) infodata%atm_prognostic = atm_prognostic
    if ( present(lnd_present)    ) infodata%lnd_present    = lnd_present
    if ( present(lnd_prognostic) ) infodata%lnd_prognostic = lnd_prognostic
    if ( present(rof_present)    ) infodata%rof_present    = rof_present
    if ( present(rofice_present) ) infodata%rofice_present = rofice_present
    if ( present(rof_prognostic) ) infodata%rof_prognostic = rof_prognostic
    if ( present(flood_present)  ) infodata%flood_present  = flood_present
    if ( present(ocn_present)    ) infodata%ocn_present    = ocn_present
    if ( present(ocn_prognostic) ) infodata%ocn_prognostic = ocn_prognostic
    if ( present(ocnrof_prognostic)) infodata%ocnrof_prognostic = ocnrof_prognostic
    if ( present(ice_present)    ) infodata%ice_present    = ice_present
    if ( present(ice_prognostic) ) infodata%ice_prognostic = ice_prognostic
    if ( present(iceberg_prognostic)) infodata%iceberg_prognostic = iceberg_prognostic
    if ( present(glc_present)    ) infodata%glc_present    = glc_present
    if ( present(glclnd_present) ) infodata%glclnd_present = glclnd_present
    if ( present(glcocn_present) ) infodata%glcocn_present = glcocn_present
    if ( present(glcice_present) ) infodata%glcice_present = glcice_present
    if ( present(glc_prognostic) ) infodata%glc_prognostic = glc_prognostic
    if ( present(wav_present)    ) infodata%wav_present    = wav_present
    if ( present(wav_prognostic) ) infodata%wav_prognostic = wav_prognostic
    if ( present(esp_present)    ) infodata%esp_present    = esp_present
    if ( present(esp_prognostic) ) infodata%esp_prognostic = esp_prognostic
    if ( present(atm_nx)         ) infodata%atm_nx         = atm_nx
    if ( present(atm_ny)         ) infodata%atm_ny         = atm_ny
    if ( present(lnd_nx)         ) infodata%lnd_nx         = lnd_nx
    if ( present(lnd_ny)         ) infodata%lnd_ny         = lnd_ny
    if ( present(rof_nx)         ) infodata%rof_nx         = rof_nx
    if ( present(rof_ny)         ) infodata%rof_ny         = rof_ny
    if ( present(ice_nx)         ) infodata%ice_nx         = ice_nx
    if ( present(ice_ny)         ) infodata%ice_ny         = ice_ny
    if ( present(ocn_nx)         ) infodata%ocn_nx         = ocn_nx
    if ( present(ocn_ny)         ) infodata%ocn_ny         = ocn_ny
    if ( present(glc_nx)         ) infodata%glc_nx         = glc_nx
    if ( present(glc_ny)         ) infodata%glc_ny         = glc_ny
    if ( present(wav_nx)         ) infodata%wav_nx         = wav_nx
    if ( present(wav_ny)         ) infodata%wav_ny         = wav_ny

    if ( present(nextsw_cday)    ) infodata%nextsw_cday    = nextsw_cday
    if ( present(precip_fact)    ) infodata%precip_fact    = precip_fact
    if ( present(atm_phase)      ) infodata%atm_phase      = atm_phase
    if ( present(lnd_phase)      ) infodata%lnd_phase      = lnd_phase
    if ( present(ice_phase)      ) infodata%ice_phase      = ice_phase
    if ( present(ocn_phase)      ) infodata%ocn_phase      = ocn_phase
    if ( present(glc_phase)      ) infodata%glc_phase      = glc_phase
    if ( present(rof_phase)      ) infodata%rof_phase      = rof_phase
    if ( present(wav_phase)      ) infodata%wav_phase      = wav_phase
    if ( present(esp_phase)      ) infodata%esp_phase      = esp_phase
    if ( present(atm_aero)       ) infodata%atm_aero       = atm_aero
    if ( present(glcrun_alarm)   ) infodata%glcrun_alarm   = glcrun_alarm
    if ( present(glc_g2lupdate)  ) infodata%glc_g2lupdate  = glc_g2lupdate
    if ( present(pause_alarm)    ) infodata%pause_alarm    = pause_alarm
    if ( present(resume_alarm)   ) infodata%resume_alarm   = resume_alarm

END SUBROUTINE seq_infodata_PutData_explicit

#ifndef CPRPGI
!===============================================================================
! !IROUTINE: seq_infodata_PutData_bytype -- Put values into infodata object
!
! !DESCRIPTION:
!
!     Put values into the infodata object.
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE seq_infodata_PutData_bytype( component_firstletter, infodata,      &
           comp_present, comp_prognostic, comp_gnam,                          &
           histavg_comp, comp_phase, comp_nx, comp_ny)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(len=1),              intent(IN)  :: component_firstletter
   type(seq_infodata_type),       intent(INOUT) :: infodata      ! Input CCSM structure
   logical             ,optional, intent(IN) :: comp_present    ! provide data
   logical             ,optional, intent(IN) :: comp_prognostic ! need data
   character(len=*)    ,optional, intent(IN) :: comp_gnam      ! comp grid
   integer(SHR_KIND_IN),optional, intent(IN) :: comp_nx        ! nx,ny 2d grid size global
   integer(SHR_KIND_IN),optional, intent(IN) :: comp_ny        ! nx,ny 2d grid size global
   integer(SHR_KIND_IN),optional, intent(IN) :: comp_phase
   logical             ,optional, intent(IN) :: histavg_comp

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(seq_infodata_PutData_bytype) '

!-------------------------------------------------------------------------------

    if (component_firstletter == 'a') then
      call seq_infodata_PutData(infodata, atm_present=comp_present,           &
           atm_prognostic=comp_prognostic, atm_gnam=comp_gnam,                &
           atm_phase=comp_phase, atm_nx=comp_nx, atm_ny=comp_ny, histavg_atm=histavg_comp)
    else if (component_firstletter == 'l') then
      call seq_infodata_PutData(infodata, lnd_present=comp_present,           &
           lnd_prognostic=comp_prognostic, lnd_gnam=comp_gnam,                &
           lnd_phase=comp_phase, lnd_nx=comp_nx, lnd_ny=comp_ny, histavg_lnd=histavg_comp)
    else if (component_firstletter == 'i') then
      call seq_infodata_PutData(infodata, ice_present=comp_present,           &
           ice_prognostic=comp_prognostic, ice_gnam=comp_gnam,                &
           ice_phase=comp_phase, ice_nx=comp_nx, ice_ny=comp_ny, histavg_ice=histavg_comp)
    else if (component_firstletter == 'o') then
      call seq_infodata_PutData(infodata, ocn_present=comp_present,           &
           ocn_prognostic=comp_prognostic, ocn_gnam=comp_gnam,                &
           ocn_phase=comp_phase, ocn_nx=comp_nx, ocn_ny=comp_ny, histavg_ocn=histavg_comp)
    else if (component_firstletter == 'r') then
      call seq_infodata_PutData(infodata, rof_present=comp_present,           &
           rof_prognostic=comp_prognostic, rof_gnam=comp_gnam,                &
           rof_phase=comp_phase, rof_nx=comp_nx, rof_ny=comp_ny, histavg_rof=histavg_comp)
    else if (component_firstletter == 'g') then
      call seq_infodata_PutData(infodata, glc_present=comp_present,           &
           glc_prognostic=comp_prognostic, glc_gnam=comp_gnam,                &
           glc_phase=comp_phase, glc_nx=comp_nx, glc_ny=comp_ny, histavg_glc=histavg_comp)
    else if (component_firstletter == 'w') then
      call seq_infodata_PutData(infodata, wav_present=comp_present,           &
           wav_prognostic=comp_prognostic, wav_gnam=comp_gnam,                &
           wav_phase=comp_phase, wav_nx=comp_nx, wav_ny=comp_ny, histavg_wav=histavg_comp)
    else if (component_firstletter == 'e') then
      if ((loglevel > 1) .and. seq_comm_iamroot(1)) then
        if (present(comp_gnam)) then
          write(logunit,*) trim(subname),' Note: ESP type has no gnam property'
        end if
        if (present(comp_nx)) then
          write(logunit,*) trim(subname),' Note: ESP type has no nx property'
        end if
        if (present(comp_ny)) then
          write(logunit,*) trim(subname),' Note: ESP type has no ny property'
        end if
        if (present(histavg_comp)) then
          write(logunit,*) trim(subname),' Note: ESP type has no histavg property'
        end if
      end if
     
      call seq_infodata_PutData(infodata, esp_present=comp_present,           &
           esp_prognostic=comp_prognostic, esp_phase=comp_phase)
    else
       call shr_sys_abort( subname//": unknown component-type first letter,'"//component_firstletter//"', aborting")
    end if

END SUBROUTINE seq_infodata_PutData_bytype
#endif
! ^ ifndef CPRPGI

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_infodata_bcast -- Broadcast an infodata from root pe
!
! !DESCRIPTION:
!
! Broadcast an infodata across pes
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_infodata_bcast(infodata,mpicom)

   use shr_mpi_mod, only : shr_mpi_bcast

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(seq_infodata_type), intent(INOUT) :: infodata    ! assume valid on root pe
  integer(SHR_KIND_IN),    intent(IN)    :: mpicom      ! mpi comm

!EOP

    !----- local -----

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    call shr_mpi_bcast(infodata%cime_model,            mpicom)
    call shr_mpi_bcast(infodata%start_type,            mpicom)
    call shr_mpi_bcast(infodata%case_desc,             mpicom)
    call shr_mpi_bcast(infodata%model_version,         mpicom)
    call shr_mpi_bcast(infodata%username,              mpicom)
    call shr_mpi_bcast(infodata%hostname,              mpicom)
    call shr_mpi_bcast(infodata%case_name,             mpicom)
    call shr_mpi_bcast(infodata%timing_dir,            mpicom)
    call shr_mpi_bcast(infodata%tchkpt_dir,            mpicom)
    call shr_mpi_bcast(infodata%atm_ideal_phys,        mpicom)
    call shr_mpi_bcast(infodata%atm_adiabatic,         mpicom)
    call shr_mpi_bcast(infodata%aqua_planet,           mpicom)
    call shr_mpi_bcast(infodata%aqua_planet_sst,       mpicom)
    call shr_mpi_bcast(infodata%run_barriers,          mpicom)
    call shr_mpi_bcast(infodata%brnch_retain_casename, mpicom)
    call shr_mpi_bcast(infodata%read_restart,          mpicom)
    call shr_mpi_bcast(infodata%restart_pfile,         mpicom)
    call shr_mpi_bcast(infodata%restart_file,          mpicom)
    call shr_mpi_bcast(infodata%single_column,         mpicom)
    call shr_mpi_bcast(infodata%scmlat,                mpicom)
    call shr_mpi_bcast(infodata%scmlon,                mpicom)
    call shr_mpi_bcast(infodata%logFilePostFix,        mpicom)
    call shr_mpi_bcast(infodata%outPathRoot,           mpicom)
    call shr_mpi_bcast(infodata%perpetual,             mpicom)
    call shr_mpi_bcast(infodata%perpetual_ymd,         mpicom)
    call shr_mpi_bcast(infodata%orb_iyear,             mpicom)
    call shr_mpi_bcast(infodata%orb_iyear_align,       mpicom)
    call shr_mpi_bcast(infodata%orb_mode,              mpicom)
    call shr_mpi_bcast(infodata%orb_eccen,             mpicom)
    call shr_mpi_bcast(infodata%orb_obliq,             mpicom)
    call shr_mpi_bcast(infodata%orb_mvelp,             mpicom)
    call shr_mpi_bcast(infodata%orb_obliqr,            mpicom)
    call shr_mpi_bcast(infodata%orb_lambm0,            mpicom)
    call shr_mpi_bcast(infodata%orb_mvelpp,            mpicom)
    call shr_mpi_bcast(infodata%wv_sat_scheme,         mpicom)
    call shr_mpi_bcast(infodata%wv_sat_transition_start, mpicom)
    call shr_mpi_bcast(infodata%wv_sat_use_tables,     mpicom)
    call shr_mpi_bcast(infodata%wv_sat_table_spacing,  mpicom)
    call shr_mpi_bcast(infodata%tfreeze_option,        mpicom)
    call shr_mpi_bcast(infodata%flux_epbal,            mpicom)
    call shr_mpi_bcast(infodata%flux_albav,            mpicom)
    call shr_mpi_bcast(infodata%flux_diurnal,          mpicom)
    call shr_mpi_bcast(infodata%wall_time_limit,       mpicom)
    call shr_mpi_bcast(infodata%force_stop_at,         mpicom)
    call shr_mpi_bcast(infodata%atm_gnam,              mpicom)
    call shr_mpi_bcast(infodata%lnd_gnam,              mpicom)
    call shr_mpi_bcast(infodata%ocn_gnam,              mpicom)
    call shr_mpi_bcast(infodata%ice_gnam,              mpicom)
    call shr_mpi_bcast(infodata%rof_gnam,              mpicom)
    call shr_mpi_bcast(infodata%glc_gnam,              mpicom)
    call shr_mpi_bcast(infodata%wav_gnam,              mpicom)
    call shr_mpi_bcast(infodata%shr_map_dopole,        mpicom)
    call shr_mpi_bcast(infodata%vect_map,              mpicom)
    call shr_mpi_bcast(infodata%aoflux_grid,           mpicom)
    call shr_mpi_bcast(infodata%cpl_decomp,            mpicom)
    call shr_mpi_bcast(infodata%cpl_seq_option,        mpicom)
    call shr_mpi_bcast(infodata%cpl_cdf64,             mpicom)
    call shr_mpi_bcast(infodata%do_budgets,            mpicom)
    call shr_mpi_bcast(infodata%do_histinit,           mpicom)
    call shr_mpi_bcast(infodata%budget_inst,           mpicom)
    call shr_mpi_bcast(infodata%budget_daily,          mpicom)
    call shr_mpi_bcast(infodata%budget_month,          mpicom)
    call shr_mpi_bcast(infodata%budget_ann,            mpicom)
    call shr_mpi_bcast(infodata%budget_ltann,          mpicom)
    call shr_mpi_bcast(infodata%budget_ltend,          mpicom)
    call shr_mpi_bcast(infodata%histaux_a2x           ,mpicom)
    call shr_mpi_bcast(infodata%histaux_a2x1hri       ,mpicom)
    call shr_mpi_bcast(infodata%histaux_a2x1hr        ,mpicom)
    call shr_mpi_bcast(infodata%histaux_a2x3hr        ,mpicom)
    call shr_mpi_bcast(infodata%histaux_a2x3hrp       ,mpicom)
    call shr_mpi_bcast(infodata%histaux_a2x24hr       ,mpicom)
    call shr_mpi_bcast(infodata%histaux_l2x1yr        ,mpicom)
    call shr_mpi_bcast(infodata%histaux_l2x           ,mpicom)
    call shr_mpi_bcast(infodata%histaux_r2x           ,mpicom)
    call shr_mpi_bcast(infodata%histavg_atm           ,mpicom)
    call shr_mpi_bcast(infodata%histavg_lnd           ,mpicom)
    call shr_mpi_bcast(infodata%histavg_ocn           ,mpicom)
    call shr_mpi_bcast(infodata%histavg_ice           ,mpicom)
    call shr_mpi_bcast(infodata%histavg_rof           ,mpicom)
    call shr_mpi_bcast(infodata%histavg_glc           ,mpicom)
    call shr_mpi_bcast(infodata%histavg_wav           ,mpicom)
    call shr_mpi_bcast(infodata%histavg_xao           ,mpicom)
    call shr_mpi_bcast(infodata%drv_threading,         mpicom)
    call shr_mpi_bcast(infodata%eps_frac,              mpicom)
    call shr_mpi_bcast(infodata%eps_amask,             mpicom)
    call shr_mpi_bcast(infodata%eps_agrid,             mpicom)
    call shr_mpi_bcast(infodata%eps_aarea,             mpicom)
    call shr_mpi_bcast(infodata%eps_omask,             mpicom)
    call shr_mpi_bcast(infodata%eps_ogrid,             mpicom)
    call shr_mpi_bcast(infodata%eps_oarea,             mpicom)
    call shr_mpi_bcast(infodata%reprosum_use_ddpdd,    mpicom)
    call shr_mpi_bcast(infodata%reprosum_diffmax,      mpicom)
    call shr_mpi_bcast(infodata%reprosum_recompute,    mpicom)
    call shr_mpi_bcast(infodata%mct_usealltoall,       mpicom)
    call shr_mpi_bcast(infodata%mct_usevector,         mpicom)

    call shr_mpi_bcast(infodata%info_debug,            mpicom)
    call shr_mpi_bcast(infodata%bfbflag,               mpicom)
    call shr_mpi_bcast(infodata%esmf_map_flag,         mpicom)
    call shr_mpi_bcast(infodata%dead_comps,            mpicom)

    call shr_mpi_bcast(infodata%atm_present,           mpicom)
    call shr_mpi_bcast(infodata%atm_prognostic,        mpicom)
    call shr_mpi_bcast(infodata%lnd_present,           mpicom)
    call shr_mpi_bcast(infodata%lnd_prognostic,        mpicom)
    call shr_mpi_bcast(infodata%rof_present,           mpicom)
    call shr_mpi_bcast(infodata%rofice_present,        mpicom)
    call shr_mpi_bcast(infodata%rof_prognostic,        mpicom)
    call shr_mpi_bcast(infodata%flood_present,         mpicom)
    call shr_mpi_bcast(infodata%ocn_present,           mpicom)
    call shr_mpi_bcast(infodata%ocn_prognostic,        mpicom)
    call shr_mpi_bcast(infodata%ocnrof_prognostic,     mpicom)
    call shr_mpi_bcast(infodata%ice_present,           mpicom)
    call shr_mpi_bcast(infodata%ice_prognostic,        mpicom)
    call shr_mpi_bcast(infodata%iceberg_prognostic,    mpicom)
    call shr_mpi_bcast(infodata%glc_present,           mpicom)
    call shr_mpi_bcast(infodata%glclnd_present,        mpicom)
    call shr_mpi_bcast(infodata%glcocn_present,        mpicom)
    call shr_mpi_bcast(infodata%glcice_present,        mpicom)
    call shr_mpi_bcast(infodata%glc_prognostic,        mpicom)
    call shr_mpi_bcast(infodata%wav_present,           mpicom)
    call shr_mpi_bcast(infodata%wav_prognostic,        mpicom)
    call shr_mpi_bcast(infodata%esp_present,           mpicom)
    call shr_mpi_bcast(infodata%esp_prognostic,        mpicom)

    call shr_mpi_bcast(infodata%atm_nx,                mpicom)
    call shr_mpi_bcast(infodata%atm_ny,                mpicom)
    call shr_mpi_bcast(infodata%lnd_nx,                mpicom)
    call shr_mpi_bcast(infodata%lnd_ny,                mpicom)
    call shr_mpi_bcast(infodata%rof_nx,                mpicom)
    call shr_mpi_bcast(infodata%rof_ny,                mpicom)
    call shr_mpi_bcast(infodata%ice_nx,                mpicom)
    call shr_mpi_bcast(infodata%ice_ny,                mpicom)
    call shr_mpi_bcast(infodata%ocn_nx,                mpicom)
    call shr_mpi_bcast(infodata%ocn_ny,                mpicom)
    call shr_mpi_bcast(infodata%glc_nx,                mpicom)
    call shr_mpi_bcast(infodata%glc_ny,                mpicom)
    call shr_mpi_bcast(infodata%wav_nx,                mpicom)
    call shr_mpi_bcast(infodata%wav_ny,                mpicom)

    call shr_mpi_bcast(infodata%nextsw_cday,           mpicom)
    call shr_mpi_bcast(infodata%precip_fact,           mpicom)
    call shr_mpi_bcast(infodata%atm_phase,             mpicom)
    call shr_mpi_bcast(infodata%lnd_phase,             mpicom)
    call shr_mpi_bcast(infodata%ice_phase,             mpicom)
    call shr_mpi_bcast(infodata%ocn_phase,             mpicom)
    call shr_mpi_bcast(infodata%glc_phase,             mpicom)
    call shr_mpi_bcast(infodata%rof_phase,             mpicom)
    call shr_mpi_bcast(infodata%wav_phase,             mpicom)
    call shr_mpi_bcast(infodata%atm_aero,              mpicom)
    call shr_mpi_bcast(infodata%glcrun_alarm,          mpicom)
    call shr_mpi_bcast(infodata%glc_g2lupdate,         mpicom)
    call shr_mpi_bcast(infodata%pause_alarm,           mpicom)
    call shr_mpi_bcast(infodata%resume_alarm,          mpicom)

end subroutine seq_infodata_bcast

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_infodata_Exchange -- Broadcast a subset of infodata between pes
!
! !DESCRIPTION:
!
! Broadcast a subset of infodata data between pes to support "exchange" of information
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_infodata_Exchange(infodata,ID,type)

   use shr_mpi_mod, only : shr_mpi_bcast

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(seq_infodata_type), intent(INOUT) :: infodata    ! assume valid on root pe
  integer(SHR_KIND_IN),    intent(IN)    :: ID          ! mpi comm
  character(len=*),        intent(IN)    :: type        ! type

!EOP

  !----- local -----
  integer(SHR_KIND_IN) :: mpicom     ! mpicom
  integer(SHR_KIND_IN) :: pebcast    ! pe sending
  logical :: atm2cpli,atm2cplr
  logical :: lnd2cpli,lnd2cplr
  logical :: rof2cpli,rof2cplr
  logical :: ocn2cpli,ocn2cplr
  logical :: ice2cpli,ice2cplr
  logical :: glc2cpli,glc2cplr
  logical :: wav2cpli,wav2cplr
  logical :: cpl2i,cpl2r
  logical :: logset
  logical :: deads   ! local variable to hold info temporarily
  character(len=*), parameter :: subname = '(seq_infodata_Exchange) '

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

  ! assume the comp pe is going to broadcast, change to cplpe below if appropriate
  call seq_comm_setptrs(ID,mpicom=mpicom,cmppe=pebcast)

  logset = .false.

  atm2cpli = .false.
  atm2cplr = .false.
  lnd2cpli = .false.
  lnd2cplr = .false.
  rof2cpli = .false.
  rof2cplr = .false.
  ocn2cpli = .false.
  ocn2cplr = .false.
  ice2cpli = .false.
  ice2cplr = .false.
  glc2cpli = .false.
  glc2cplr = .false.
  wav2cpli = .false.
  wav2cplr = .false.
  cpl2i = .false.
  cpl2r = .false.

  ! --- translate type into logicals ---

  if (trim(type) == 'atm2cpl_init') then
     atm2cpli = .true.
     atm2cplr = .true.
     logset = .true.
  endif
  if (trim(type) == 'atm2cpl_run') then
     atm2cplr = .true.
     logset = .true.
  endif

  if (trim(type) == 'lnd2cpl_init') then
     lnd2cpli = .true.
     lnd2cplr = .true.
     logset = .true.
  endif
  if (trim(type) == 'lnd2cpl_run') then
     lnd2cplr = .true.
     logset = .true.
  endif

  if (trim(type) == 'rof2cpl_init') then
     rof2cpli = .true.
     rof2cplr = .true.
     logset = .true.
  endif
  if (trim(type) == 'rof2cpl_run') then
     rof2cplr = .true.
     logset = .true.
  endif

  if (trim(type) == 'ocn2cpl_init') then
     ocn2cpli = .true.
     ocn2cplr = .true.
     logset = .true.
  endif
  if (trim(type) == 'ocn2cpl_run') then
     ocn2cplr = .true.
     logset = .true.
  endif

  if (trim(type) == 'ice2cpl_init') then
     ice2cpli = .true.
     ice2cplr = .true.
     logset = .true.
  endif
  if (trim(type) == 'ice2cpl_run') then
     ice2cplr = .true.
     logset = .true.
  endif

  if (trim(type) == 'glc2cpl_init') then
     glc2cpli = .true.
     glc2cplr = .true.
     logset = .true.
  endif
  if (trim(type) == 'glc2cpl_run') then
     glc2cplr = .true.
     logset = .true.
  endif

  if (trim(type) == 'wav2cpl_init') then
     wav2cpli = .true.
     wav2cplr = .true.
     logset = .true.
  endif
  if (trim(type) == 'wav2cpl_run') then
     wav2cplr = .true.
     logset = .true.
  endif

  if (trim(type) == 'cpl2atm_init' .or. &
      trim(type) == 'cpl2lnd_init' .or. &
      trim(type) == 'cpl2rof_init' .or. &
      trim(type) == 'cpl2ocn_init' .or. &
      trim(type) == 'cpl2glc_init' .or. &
      trim(type) == 'cpl2wav_init' .or. &
      trim(type) == 'cpl2ice_init') then
     cpl2i = .true.
     cpl2r = .true.
     call seq_comm_setptrs(ID,cplpe=pebcast)
     logset = .true.
  endif

  if (trim(type) == 'cpl2atm_run' .or. &
      trim(type) == 'cpl2lnd_run' .or. &
      trim(type) == 'cpl2rof_run' .or. &
      trim(type) == 'cpl2ocn_run' .or. &
      trim(type) == 'cpl2glc_run' .or. &
      trim(type) == 'cpl2wav_run' .or. &
      trim(type) == 'cpl2ice_run') then
     cpl2r = .true.
     call seq_comm_setptrs(ID,cplpe=pebcast)
     logset = .true.
  endif

  ! --- make sure the type was valid ---

  if (.not. logset) then
     write(logunit,*) trim(subname),' ERROR: type invalid ',trim(type)
     call shr_sys_abort()
  endif

  ! --- now execute exchange ---

  if (atm2cpli) then
    call shr_mpi_bcast(infodata%atm_present,      mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%atm_prognostic,   mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%atm_nx,           mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%atm_ny,           mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%atm_aero,         mpicom,pebcast=pebcast)
    ! dead_comps is true if it's ever set to true
    deads = infodata%dead_comps
    call shr_mpi_bcast(deads,                     mpicom,pebcast=pebcast)
    if (deads .or. infodata%dead_comps) infodata%dead_comps = .true.
  endif

  if (lnd2cpli) then
    call shr_mpi_bcast(infodata%lnd_present,      mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%lnd_prognostic,   mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%lnd_nx,           mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%lnd_ny,           mpicom,pebcast=pebcast)
    ! dead_comps is true if it's ever set to true
    deads = infodata%dead_comps
    call shr_mpi_bcast(deads,                     mpicom,pebcast=pebcast)
    if (deads .or. infodata%dead_comps) infodata%dead_comps = .true.
  endif

  if (rof2cpli) then
    call shr_mpi_bcast(infodata%rof_present,      mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%rofice_present,   mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%rof_prognostic,   mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%rof_nx,           mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%rof_ny,           mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%flood_present,    mpicom,pebcast=pebcast)
    ! dead_comps is true if it's ever set to true
    deads = infodata%dead_comps
    call shr_mpi_bcast(deads,                     mpicom,pebcast=pebcast)
    if (deads .or. infodata%dead_comps) infodata%dead_comps = .true.
  endif

  if (ocn2cpli) then
    call shr_mpi_bcast(infodata%ocn_present,      mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%ocn_prognostic,   mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%ocnrof_prognostic,mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%ocn_nx,           mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%ocn_ny,           mpicom,pebcast=pebcast)
    ! dead_comps is true if it's ever set to true
    deads = infodata%dead_comps
    call shr_mpi_bcast(deads,                     mpicom,pebcast=pebcast)
    if (deads .or. infodata%dead_comps) infodata%dead_comps = .true.
  endif

  if (ice2cpli) then
    call shr_mpi_bcast(infodata%ice_present,      mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%ice_prognostic,   mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%iceberg_prognostic,mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%ice_nx,           mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%ice_ny,           mpicom,pebcast=pebcast)
    ! dead_comps is true if it's ever set to true
    deads = infodata%dead_comps
    call shr_mpi_bcast(deads,                     mpicom,pebcast=pebcast)
    if (deads .or. infodata%dead_comps) infodata%dead_comps = .true.
  endif

  if (glc2cpli) then
    call shr_mpi_bcast(infodata%glc_present,      mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%glclnd_present,   mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%glcocn_present,   mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%glcice_present,   mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%glc_prognostic,   mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%glc_nx,           mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%glc_ny,           mpicom,pebcast=pebcast)
    ! dead_comps is true if it's ever set to true
    deads = infodata%dead_comps
    call shr_mpi_bcast(deads,                     mpicom,pebcast=pebcast)
    if (deads .or. infodata%dead_comps) infodata%dead_comps = .true.
  endif

  if (wav2cpli) then
    call shr_mpi_bcast(infodata%wav_present,      mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%wav_prognostic,   mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%wav_nx,           mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%wav_ny,           mpicom,pebcast=pebcast)
    ! dead_comps is true if it's ever set to true
    deads = infodata%dead_comps
    call shr_mpi_bcast(deads,                     mpicom,pebcast=pebcast)
    if (deads .or. infodata%dead_comps) infodata%dead_comps = .true.
  endif

  if (cpl2i) then
    call shr_mpi_bcast(infodata%atm_present,      mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%atm_prognostic,   mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%lnd_present,      mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%lnd_prognostic,   mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%rof_present,      mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%rofice_present,   mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%rof_prognostic,   mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%flood_present,    mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%ocn_present,      mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%ocn_prognostic,   mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%ocnrof_prognostic,mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%ice_present,      mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%ice_prognostic,   mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%iceberg_prognostic,mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%glc_present,      mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%glclnd_present,   mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%glcocn_present,   mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%glcice_present,   mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%glc_prognostic,   mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%wav_present,      mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%wav_prognostic,   mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%dead_comps,       mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%atm_aero,         mpicom,pebcast=pebcast)
  endif

  if (atm2cplr) then
    call shr_mpi_bcast(infodata%nextsw_cday,      mpicom,pebcast=pebcast)
  endif

  if (ocn2cplr) then
    call shr_mpi_bcast(infodata%precip_fact,      mpicom,pebcast=pebcast)
  endif

  if (cpl2r) then
    call shr_mpi_bcast(infodata%nextsw_cday,      mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%precip_fact,      mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%glcrun_alarm,     mpicom,pebcast=pebcast)
    call shr_mpi_bcast(infodata%glc_g2lupdate,    mpicom,pebcast=pebcast)
  endif

end subroutine seq_infodata_Exchange

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: seq_infodata_Check  -- Check that input InputInfo derived type is valid
!
! !DESCRIPTION:
!
! Check that input infodata object has reasonable values
!
! !INTERFACE: ------------------------------------------------------------------

subroutine seq_infodata_Check( infodata )

! !USES:

  use shr_assert_mod,   only: shr_assert_in_domain
  use shr_string_mod,   only: shr_string_listIntersect
  use shr_wv_sat_mod,   only: shr_wv_sat_get_scheme_idx, shr_wv_sat_valid_idx

  implicit none

! !INPUT/OUTPUT PARAMETERS:

    type(seq_infodata_type), intent(INOUT) :: infodata    ! Output CCSM structure

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(seq_infodata_Check) '
    integer :: lastchar                        ! Last character index
    integer :: rc                              ! Return code

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

    ! --- CIME model ------
    if ( trim(infodata%cime_model) /= 'acme' .and. trim(infodata%cime_model) /= 'cesm') then
       call shr_sys_abort( subname//': cime_model must be set to acme or cesm, aborting')
    end if

    ! --- Case name ------
    lastchar = len(infodata%case_name)
    if ( len_trim(infodata%case_name) == 0) then
       call shr_sys_abort( subname//': variable case_name must be set, aborting')
    end if
    if (infodata%case_name(lastchar:lastchar) /= ' ') then
       write(logunit,"(A,I4,A)")'ERROR: case_name must not exceed ', len(infodata%case_name)-1, &
                 ' characters'
       call shr_sys_abort( subname//': variable case_name must be set, aborting')
    end if

    ! --- Special configurations -----
    if ( infodata%atm_adiabatic .and. (infodata%atm_ideal_phys .or. &
         infodata%aqua_planet) )then
       call shr_sys_abort( subname//': only one of atm_adiabatic, ' // &
                           'atm_ideal_phys or aqua_planet can be set' )
    end if

    ! --- Restart pointer file -----
    if ( len_trim(infodata%restart_pfile) == 0 ) then
       call shr_sys_abort( subname//': restart_pfile must be set' )
    end if

    ! --- LogFile ending name -----
    if ( len_trim(infodata%logFilePostFix) == 0 ) then
       call shr_sys_abort( subname//': logFilePostFix  must be set to something not blank' )
    end if

    ! --- Output path root directory -----
    if ( len_trim(infodata%outPathRoot) == 0 ) then
       call shr_sys_abort( subname//': outPathRoot  must be set' )
    end if
    if ( index(infodata%outPathRoot,"/",back=.true.) /= &
         len_trim(infodata%outPathRoot) ) then
       call shr_sys_abort( subname//': outPathRoot must end with a slash' )
    end if

    ! --- Start-type ------
    if ((trim(infodata%start_type) /= seq_infodata_start_type_start) .and.  &
        (trim(infodata%start_type) /= seq_infodata_start_type_cont ) .and.  &
        (trim(infodata%start_type) /= seq_infodata_start_type_brnch)) then
       call shr_sys_abort(subname//': start_type invalid = '//trim(infodata%start_type))
    end if

    if ((trim(infodata%start_type) == seq_infodata_start_type_cont ) .and.  &
        (trim(infodata%case_name)  /= trim(infodata%rest_case_name))) then
       write(logunit,'(10a)') subname,' case_name =',trim(infodata%case_name),':', &
                                      ' rest_case_name =',trim(infodata%rest_case_name),':'
       call shr_sys_abort(subname//': invalid continue restart case name = '//trim(infodata%rest_case_name))
    endif

    if (infodata%orb_eccen  == SHR_ORB_UNDEF_REAL .or. &
        infodata%orb_obliqr == SHR_ORB_UNDEF_REAL .or. &
        infodata%orb_mvelpp == SHR_ORB_UNDEF_REAL .or. &
        infodata%orb_lambm0 == SHR_ORB_UNDEF_REAL) then
       call shr_sys_abort(subname//': orb params incorrect')
    endif

    if (.not. shr_wv_sat_valid_idx(shr_wv_sat_get_scheme_idx(trim(infodata%wv_sat_scheme)))) then
       call shr_sys_abort(subname//': "'//trim(infodata%wv_sat_scheme)//'" &
            &is not a recognized saturation vapor pressure scheme name')
    end if

    ! A transition range averaging method in CAM is only valid for:
    !
    ! -40 deg C <= T <= 0 deg C
    !
    ! shr_wv_sat_mod itself checks for values with the wrong sign, but we
    ! have to check that the range is no more than 40 deg C here. Even
    ! though this is a CAM-specific restriction, it's not really likely
    ! that any other parameterization will be dealing with mixed-phase
    ! water below 40 deg C anyway.
    call shr_assert_in_domain(infodata%wv_sat_transition_start, &
         ge=0._SHR_KIND_R8, le=40._SHR_KIND_R8, &
         varname="wv_sat_transition_start",&
         msg="Invalid transition temperature range.")

    if ((trim(infodata%aoflux_grid) /= 'ocn') .and. &
        (trim(infodata%aoflux_grid) /= 'atm') .and. &
        (trim(infodata%aoflux_grid) /= 'exch')) then
       write(logunit,'(2a)') 'ERROR aoflux_grid not supported = ',trim(infodata%aoflux_grid)
       call shr_sys_abort(subname//': aoflux_grid invalid = '//trim(infodata%aoflux_grid))
    endif

    if ((trim(infodata%vect_map) /= 'none') .and. &
        (trim(infodata%vect_map) /= 'cart3d') .and. &
        (trim(infodata%vect_map) /= 'cart3d_diag') .and. &
        (trim(infodata%vect_map) /= 'cart3d_uvw') .and. &
        (trim(infodata%vect_map) /= 'cart3d_uvw_diag')) then
       write(logunit,'(2a)') 'ERROR vect_map not supported = ',trim(infodata%vect_map)
       call shr_sys_abort(subname//': vect_map invalid = '//trim(infodata%vect_map))
    endif

    if (infodata%pause_alarm .and. infodata%resume_alarm) then
       call shr_sys_abort(subname//': pause_alarm and resume_alarm should not both be .true.')
    endif

END SUBROUTINE seq_infodata_Check

!===============================================================================
!===============================================================================
! !IROUTINE: seq_infodata_print -- Print out values to log file
!
! !DESCRIPTION:
!
!     Print derivied type out to screen.
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE seq_infodata_print( infodata )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(seq_infodata_type), intent(IN) :: infodata  ! Output CCSM structure

!EOP

    !----- local -----
    character(len=*), parameter :: subname = '(seq_infodata_print) '
    character(len=*), parameter ::  F0A = "(2A,A)"
    character(len=*), parameter ::  F0L = "(2A,L3)"
    character(len=*), parameter ::  F0I = "(2A,I10)"
    character(len=*), parameter ::  F0S = "(2A,I4)"
    character(len=*), parameter ::  F0R = "(2A,g22.14)"

!-------------------------------------------------------------------------------
! Notes:
!-------------------------------------------------------------------------------

!    if (loglevel > 0) then
       write(logunit,F0A) subname,'CIME model               = ', trim(infodata%cime_model)
       write(logunit,F0A) subname,'Start type               = ', trim(infodata%start_type)
       write(logunit,F0A) subname,'Case name                = ', trim(infodata%case_name)
       write(logunit,F0A) subname,'Case description         = ', trim(infodata%case_desc)
       write(logunit,F0A) subname,'Model version            = ', trim(infodata%model_version)
       write(logunit,F0A) subname,'Username                 = ', trim(infodata%username)
       write(logunit,F0A) subname,'Hostname                 = ', trim(infodata%hostname)
       write(logunit,F0A) subname,'Timing Dir               = ', trim(infodata%timing_dir)
       write(logunit,F0A) subname,'Timing Checkpoint Dir    = ', trim(infodata%tchkpt_dir)
       write(logunit,F0A) subname,'Restart case name        = ', trim(infodata%rest_case_name)

       write(logunit,F0L) subname,'atm_ideal_phys           = ', infodata%atm_ideal_phys
       write(logunit,F0L) subname,'atm adiabatic mode       = ', infodata%atm_adiabatic
       write(logunit,F0L) subname,'aqua_planet mode         = ', infodata%aqua_planet
       write(logunit,F0I) subname,'aqua_planet analytic sst = ', infodata%aqua_planet_sst
       write(logunit,F0L) subname,'brnch_retain_casename    = ', infodata%brnch_retain_casename

       write(logunit,F0L) subname,'read_restart flag        = ', infodata%read_restart
       write(logunit,F0A) subname,'Restart pointer file     = ', trim(infodata%restart_pfile)
       write(logunit,F0A) subname,'Restart file (full path) = ', trim(infodata%restart_file)

       write(logunit,F0L) subname,'single_column            = ', infodata%single_column
       write(logunit,F0R) subname,'scmlat                   = ', infodata%scmlat
       write(logunit,F0R) subname,'scmlon                   = ', infodata%scmlon

       write(logunit,F0A) subname,'Log output end name      = ', trim(infodata%logFilePostFix)
       write(logunit,F0A) subname,'Output path dir          = ', trim(infodata%outPathRoot)

       write(logunit,F0L) subname,'perpetual                = ', infodata%perpetual
       write(logunit,F0I) subname,'perpetual_ymd            = ', infodata%perpetual_ymd

       write(logunit,F0A) subname,'orb_mode                 = ', trim(infodata%orb_mode)
     if (trim(infodata%orb_mode) == trim(seq_infodata_orb_fixed_parameters)) then
       write(logunit,F0R) subname,'orb_eccen                = ', infodata%orb_eccen
       write(logunit,F0R) subname,'orb_obliq                = ', infodata%orb_obliq
       write(logunit,F0R) subname,'orb_mvelp                = ', infodata%orb_mvelp
       write(logunit,F0R) subname,'orb_obliqr               = ', infodata%orb_obliqr
       write(logunit,F0R) subname,'orb_mvelpp               = ', infodata%orb_mvelpp
       write(logunit,F0R) subname,'orb_lambm0               = ', infodata%orb_lambm0
     elseif (trim(infodata%orb_mode) == trim(seq_infodata_orb_fixed_year)) then
       write(logunit,F0I) subname,'orb_iyear                = ', infodata%orb_iyear
       write(logunit,F0R) subname,'orb_eccen                = ', infodata%orb_eccen
       write(logunit,F0R) subname,'orb_obliq                = ', infodata%orb_obliq
       write(logunit,F0R) subname,'orb_mvelp                = ', infodata%orb_mvelp
       write(logunit,F0R) subname,'orb_obliqr               = ', infodata%orb_obliqr
       write(logunit,F0R) subname,'orb_mvelpp               = ', infodata%orb_mvelpp
       write(logunit,F0R) subname,'orb_lambm0               = ', infodata%orb_lambm0
     elseif (trim(infodata%orb_mode) == trim(seq_infodata_orb_variable_year)) then
       write(logunit,F0I) subname,'orb_iyear                = ', infodata%orb_iyear
       write(logunit,F0I) subname,'orb_iyear_align          = ', infodata%orb_iyear_align
     endif

       write(logunit,F0A) subname,'wv_sat_scheme            = ', trim(infodata%wv_sat_scheme)
       write(logunit,F0R) subname,'wv_sat_transition_start  = ', infodata%wv_sat_transition_start
       write(logunit,F0L) subname,'wv_sat_use_tables        = ', infodata%wv_sat_use_tables
       write(logunit,F0R) subname,'wv_sat_table_spacing     = ', infodata%wv_sat_table_spacing

       write(logunit,F0A) subname,'tfreeze_option           = ', trim(infodata%tfreeze_option)
       write(logunit,F0A) subname,'flux_epbal               = ', trim(infodata%flux_epbal)
       write(logunit,F0L) subname,'flux_albav               = ', infodata%flux_albav
       write(logunit,F0L) subname,'flux_diurnal             = ', infodata%flux_diurnal
       write(logunit,F0R) subname,'wall_time_limit          = ', infodata%wall_time_limit
       write(logunit,F0A) subname,'force_stop_at            = ', trim(infodata%force_stop_at)
       write(logunit,F0A) subname,'atm_gridname             = ', trim(infodata%atm_gnam)
       write(logunit,F0A) subname,'lnd_gridname             = ', trim(infodata%lnd_gnam)
       write(logunit,F0A) subname,'ocn_gridname             = ', trim(infodata%ocn_gnam)
       write(logunit,F0A) subname,'ice_gridname             = ', trim(infodata%ice_gnam)
       write(logunit,F0A) subname,'rof_gridname             = ', trim(infodata%rof_gnam)
       write(logunit,F0A) subname,'glc_gridname             = ', trim(infodata%glc_gnam)
       write(logunit,F0A) subname,'wav_gridname             = ', trim(infodata%wav_gnam)
       write(logunit,F0L) subname,'shr_map_dopole           = ', infodata%shr_map_dopole
       write(logunit,F0A) subname,'vect_map                 = ', trim(infodata%vect_map)
       write(logunit,F0A) subname,'aoflux_grid              = ', trim(infodata%aoflux_grid)
       write(logunit,F0A) subname,'cpl_seq_option           = ', trim(infodata%cpl_seq_option)
       write(logunit,F0S) subname,'cpl_decomp               = ', infodata%cpl_decomp
       write(logunit,F0L) subname,'cpl_cdf64                = ', infodata%cpl_cdf64
       write(logunit,F0L) subname,'do_budgets               = ', infodata%do_budgets
       write(logunit,F0L) subname,'do_histinit              = ', infodata%do_histinit
       write(logunit,F0S) subname,'budget_inst              = ', infodata%budget_inst
       write(logunit,F0S) subname,'budget_daily             = ', infodata%budget_daily
       write(logunit,F0S) subname,'budget_month             = ', infodata%budget_month
       write(logunit,F0S) subname,'budget_ann               = ', infodata%budget_ann
       write(logunit,F0S) subname,'budget_ltann             = ', infodata%budget_ltann
       write(logunit,F0S) subname,'budget_ltend             = ', infodata%budget_ltend
       write(logunit,F0L) subname,'histaux_a2x              = ', infodata%histaux_a2x
       write(logunit,F0L) subname,'histaux_a2x1hri          = ', infodata%histaux_a2x1hri
       write(logunit,F0L) subname,'histaux_a2x1hr           = ', infodata%histaux_a2x1hr
       write(logunit,F0L) subname,'histaux_a2x3hr           = ', infodata%histaux_a2x3hr
       write(logunit,F0L) subname,'histaux_a2x3hrp          = ', infodata%histaux_a2x3hrp
       write(logunit,F0L) subname,'histaux_a2x24hr          = ', infodata%histaux_a2x24hr
       write(logunit,F0L) subname,'histaux_l2x1yr           = ', infodata%histaux_l2x1yr
       write(logunit,F0L) subname,'histaux_l2x              = ', infodata%histaux_l2x
       write(logunit,F0L) subname,'histaux_r2x              = ', infodata%histaux_r2x
       write(logunit,F0L) subname,'histavg_atm              = ', infodata%histavg_atm
       write(logunit,F0L) subname,'histavg_lnd              = ', infodata%histavg_lnd
       write(logunit,F0L) subname,'histavg_ocn              = ', infodata%histavg_ocn
       write(logunit,F0L) subname,'histavg_ice              = ', infodata%histavg_ice
       write(logunit,F0L) subname,'histavg_rof              = ', infodata%histavg_rof
       write(logunit,F0L) subname,'histavg_glc              = ', infodata%histavg_glc
       write(logunit,F0L) subname,'histavg_wav              = ', infodata%histavg_wav
       write(logunit,F0L) subname,'histavg_xao              = ', infodata%histavg_xao
       write(logunit,F0L) subname,'drv_threading            = ', infodata%drv_threading

       write(logunit,F0R) subname,'eps_frac                 = ', infodata%eps_frac
       write(logunit,F0R) subname,'eps_amask                = ', infodata%eps_amask
       write(logunit,F0R) subname,'eps_agrid                = ', infodata%eps_agrid
       write(logunit,F0R) subname,'eps_aarea                = ', infodata%eps_aarea
       write(logunit,F0R) subname,'eps_omask                = ', infodata%eps_omask
       write(logunit,F0R) subname,'eps_ogrid                = ', infodata%eps_ogrid
       write(logunit,F0R) subname,'eps_oarea                = ', infodata%eps_oarea

       write(logunit,F0L) subname,'reprosum_use_ddpdd       = ', infodata%reprosum_use_ddpdd
       write(logunit,F0R) subname,'reprosum_diffmax         = ', infodata%reprosum_diffmax
       write(logunit,F0L) subname,'reprosum_recompute       = ', infodata%reprosum_recompute

       write(logunit,F0L) subname,'mct_usealltoall          = ', infodata%mct_usealltoall
       write(logunit,F0L) subname,'mct_usevector            = ', infodata%mct_usevector

       write(logunit,F0S) subname,'info_debug               = ', infodata%info_debug
       write(logunit,F0L) subname,'bfbflag                  = ', infodata%bfbflag
       write(logunit,F0L) subname,'esmf_map_flag            = ', infodata%esmf_map_flag
       write(logunit,F0L) subname,'dead_comps               = ', infodata%dead_comps
       write(logunit,F0L) subname,'run_barriers             = ', infodata%run_barriers

       write(logunit,F0L) subname,'atm_present              = ', infodata%atm_present
       write(logunit,F0L) subname,'atm_prognostic           = ', infodata%atm_prognostic
       write(logunit,F0L) subname,'lnd_present              = ', infodata%lnd_present
       write(logunit,F0L) subname,'lnd_prognostic           = ', infodata%lnd_prognostic
       write(logunit,F0L) subname,'rof_present              = ', infodata%rof_present
       write(logunit,F0L) subname,'rofice_present           = ', infodata%rofice_present
       write(logunit,F0L) subname,'rof_prognostic           = ', infodata%rof_prognostic
       write(logunit,F0L) subname,'flood_present            = ', infodata%flood_present
       write(logunit,F0L) subname,'ocn_present              = ', infodata%ocn_present
       write(logunit,F0L) subname,'ocn_prognostic           = ', infodata%ocn_prognostic
       write(logunit,F0L) subname,'ocnrof_prognostic        = ', infodata%ocnrof_prognostic
       write(logunit,F0L) subname,'ice_present              = ', infodata%ice_present
       write(logunit,F0L) subname,'ice_prognostic           = ', infodata%ice_prognostic
       write(logunit,F0L) subname,'iceberg_prognostic       = ', infodata%iceberg_prognostic
       write(logunit,F0L) subname,'glc_present              = ', infodata%glc_present
       write(logunit,F0L) subname,'glclnd_present           = ', infodata%glclnd_present
       write(logunit,F0L) subname,'glcocn_present           = ', infodata%glcocn_present
       write(logunit,F0L) subname,'glcice_present           = ', infodata%glcice_present
       write(logunit,F0L) subname,'glc_prognostic           = ', infodata%glc_prognostic
       write(logunit,F0L) subname,'wav_present              = ', infodata%wav_present
       write(logunit,F0L) subname,'wav_prognostic           = ', infodata%wav_prognostic
       write(logunit,F0L) subname,'esp_present              = ', infodata%esp_present
       write(logunit,F0L) subname,'esp_prognostic           = ', infodata%esp_prognostic

       write(logunit,F0I) subname,'atm_nx                   = ', infodata%atm_nx
       write(logunit,F0I) subname,'atm_ny                   = ', infodata%atm_ny
       write(logunit,F0I) subname,'lnd_nx                   = ', infodata%lnd_nx
       write(logunit,F0I) subname,'lnd_ny                   = ', infodata%lnd_ny
       write(logunit,F0I) subname,'rof_nx                   = ', infodata%rof_nx
       write(logunit,F0I) subname,'rof_ny                   = ', infodata%rof_ny
       write(logunit,F0I) subname,'ice_nx                   = ', infodata%ice_nx
       write(logunit,F0I) subname,'ice_ny                   = ', infodata%ice_ny
       write(logunit,F0I) subname,'ocn_nx                   = ', infodata%ocn_nx
       write(logunit,F0I) subname,'ocn_ny                   = ', infodata%ocn_ny
       write(logunit,F0I) subname,'glc_nx                   = ', infodata%glc_nx
       write(logunit,F0I) subname,'glc_ny                   = ', infodata%glc_ny
       write(logunit,F0I) subname,'wav_nx                   = ', infodata%wav_nx
       write(logunit,F0I) subname,'wav_ny                   = ', infodata%wav_ny

       write(logunit,F0R) subname,'nextsw_cday              = ', infodata%nextsw_cday
       write(logunit,F0R) subname,'precip_fact              = ', infodata%precip_fact
       write(logunit,F0L) subname,'atm_aero                 = ', infodata%atm_aero

       write(logunit,F0S) subname,'atm_phase                = ', infodata%atm_phase
       write(logunit,F0S) subname,'lnd_phase                = ', infodata%lnd_phase
       write(logunit,F0S) subname,'ocn_phase                = ', infodata%ocn_phase
       write(logunit,F0S) subname,'ice_phase                = ', infodata%ice_phase
       write(logunit,F0S) subname,'glc_phase                = ', infodata%glc_phase
       write(logunit,F0S) subname,'rof_phase                = ', infodata%rof_phase
       write(logunit,F0S) subname,'wav_phase                = ', infodata%wav_phase

       write(logunit,F0L) subname,'glcrun_alarm             = ', infodata%glcrun_alarm
       write(logunit,F0L) subname,'glc_g2lupdate            = ', infodata%glc_g2lupdate
       write(logunit,F0L) subname,'pause_alarm              = ', infodata%pause_alarm
       write(logunit,F0L) subname,'resume_alarm             = ', infodata%resume_alarm
!     endif

END SUBROUTINE seq_infodata_print

!===============================================================================
!===============================================================================

END MODULE seq_infodata_mod
