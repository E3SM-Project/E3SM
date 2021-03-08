#if defined( MMF_ORIENT_RAND ) && defined( MMF_DIR_NS )
#undef MMF_DIR_NS
#endif

module crm_physics

!---------------------------------------------------------------------------------------------------
! 
! Purpose: Provides the interface to the crm code for the MMF configuration
! 
!---------------------------------------------------------------------------------------------------
   use shr_kind_mod,    only: r8 => shr_kind_r8
   use shr_sys_mod,     only: shr_sys_flush
   use cam_abortutils,  only: endrun
   use cam_logfile,     only: iulog
   use physics_types,   only: physics_state, physics_tend
   use ppgrid,          only: pcols, pver, pverp
   use constituents,    only: pcnst
#ifdef MODAL_AERO
   use modal_aero_data, only: ntot_amode
#endif

   implicit none 
   private
   save

   public :: crm_physics_register
   public :: crm_physics_init
   public :: crm_physics_tend
   public :: crm_surface_flux_bypass_tend
   public :: crm_recall_state_tend
   public :: crm_save_state_tend
   public :: m2005_effradius

   integer :: crm_qaerwat_idx, crm_dgnumwet_idx
   ! integer :: prec_dp_idx, snow_dp_idx, prec_sh_idx, snow_sh_idx
   integer :: prec_sed_idx, snow_sed_idx, snow_str_idx, prec_pcw_idx, snow_pcw_idx
   integer :: cldo_idx, cld_idx, concld_idx, iciwp_idx, iclwp_idx

   ! Physics buffer fields normally registered by offline parameterizations
   integer :: qme_idx
   integer :: prain_idx
   integer :: nevapr_idx
   integer :: wsedl_idx
   integer :: rei_idx
   integer :: rel_idx
   integer :: dei_idx
   integer :: mu_idx
   integer :: prer_evap_idx
   integer :: lambdac_idx
   integer :: iciwpst_idx
   integer :: iclwpst_idx
   integer :: des_idx
   integer :: icswp_idx
   integer :: cldfsnow_idx

   ! Constituent names
   character(len=8), parameter :: cnst_names(8) = (/'CLDLIQ', 'CLDICE','NUMLIQ','NUMICE', &
                                                    'RAINQM', 'SNOWQM','NUMRAI','NUMSNO'/)

   integer :: ixcldliq  = -1   ! cloud liquid amount index
   integer :: ixcldice  = -1   ! cloud ice amount index
   integer :: ixnumliq  = -1   ! cloud liquid number index
   integer :: ixnumice  = -1   ! cloud ice water index
   integer :: ixrain    = -1   ! rain index
   integer :: ixsnow    = -1   ! snow index
   integer :: ixnumrain = -1   ! rain number index
   integer :: ixnumsnow = -1   ! snow number index

   ! Physics buffer indices 
   integer :: icwmrdp_idx     = -1 
   integer :: rprddp_idx      = -1 
   integer :: nevapr_dpcu_idx = -1 
   integer :: prec_dp_idx     = -1
   integer :: snow_dp_idx     = -1
   integer :: ttend_dp_idx    = -1

   integer :: icwmrsh_idx     = -1
   integer :: rprdsh_idx      = -1
   integer :: rprdtot_idx     = -1
   integer :: cldtop_idx      = -1
   integer :: cldbot_idx      = -1
   integer :: nevapr_shcu_idx = -1
   integer :: prec_sh_idx     = -1
   integer :: snow_sh_idx     = -1

   integer :: ast_idx         = -1 ! stratiform cloud fraction index in physics buffer
   integer :: aist_idx        = -1 ! ice stratiform cloud fraction index in physics buffer
   integer :: alst_idx        = -1 ! liquid stratiform cloud fraction index in physics buffer
   integer :: qist_idx        = -1 ! ice stratiform in-cloud IWC 
   integer :: qlst_idx        = -1 ! liquid stratiform in-cloud LWC  
   integer :: fice_idx        = -1

   integer :: qcwat_idx       = -1 ! qcwat index in physics buffer
   integer :: lcwat_idx       = -1 ! lcwat index in physics buffer
   integer :: iccwat_idx      = -1 ! iccwat index in physics buffer
   integer :: nlwat_idx       = -1 ! nlwat index in physics buffer
   integer :: niwat_idx       = -1 ! niwat index in physics buffer
   integer :: tcwat_idx       = -1 ! tcwat index in physics buffer

   real(r8),pointer                        :: acldy_cen_tbeg(:,:)        ! cloud fraction
   real(r8), pointer, dimension(:,:)       :: cldo

   type(physics_state)                     :: state_save
   type(physics_tend)                      :: tend_save
   real(r8), dimension(:,:), allocatable   :: cldo_save     ! saved cloud fraction
   real(r8), dimension(:,:,:), allocatable :: q_aero        ! used to keep aerosol changes from offline physics
#ifdef MODAL_AERO
   real(r8), dimension(:,:,:), allocatable :: qqcw_save
   real(r8), dimension(:,:,:), allocatable :: qqcw_all
   real(r8), dimension(:,:,:), allocatable :: dgnumwet_save
   real(r8),pointer  :: dgnumwet(:,:,:)
#endif


contains
!===================================================================================================
!===================================================================================================

subroutine crm_physics_register()
!---------------------------------------------------------------------------------------------------
! 
! Purpose:  add necessary fields into physics buffer
!
!---------------------------------------------------------------------------------------------------
   use spmd_utils,          only: masterproc
   use physconst,           only: cpair, mwh2o, mwdry
   use ppgrid,              only: pcols, pver, pverp
   use physics_buffer,      only: dyn_time_lvls, pbuf_add_field, dtype_r8, pbuf_get_index
   use phys_control,        only: phys_getopts, use_gw_convect
   use constituents,        only: cnst_add
   use crmdims,             only: crm_nx, crm_ny, crm_nz, &
                                  crm_dx, crm_dy, crm_dt, &
                                  crm_nx_rad, crm_ny_rad
   use constituents,        only: cnst_add
   use crm_history,         only: crm_history_register
#if defined(MMF_SAMXX)
   use cpp_interface_mod,   only: setparm
   use gator_mod, only: gator_init
#elif defined(MMF_SAM)
   use setparm_mod      ,   only: setparm
#endif
#ifdef MODAL_AERO
   use modal_aero_data, only: ntot_amode
#endif
   !----------------------------------------------------------------------------
   ! local variables
   integer idx
   logical           :: use_ECPP
   logical           :: use_MMF_VT
   character(len=16) :: MMF_microphysics_scheme
   logical           :: prog_modal_aero
   integer, dimension(1) :: dims_gcm_1D
   integer, dimension(2) :: dims_gcm_2D
   integer, dimension(3) :: dims_crm_2D
   integer, dimension(4) :: dims_crm_3D
   integer, dimension(4) :: dims_crm_rad
#ifdef MODAL_AERO
   integer, dimension(5) :: dims_crm_aer
#endif
   integer :: cnst_ind ! dummy for adding new constituents for variance transport
   !----------------------------------------------------------------------------
#ifdef MMF_SAMXX
   call gator_init()
#endif
   dims_gcm_1D  = (/pcols/)
   dims_gcm_2D  = (/pcols, pver/)
   dims_crm_2D  = (/pcols, crm_nx, crm_ny/)
   dims_crm_3D  = (/pcols, crm_nx, crm_ny, crm_nz/)
   dims_crm_rad = (/pcols, crm_nx_rad, crm_ny_rad, crm_nz/)
#ifdef MODAL_AERO
   dims_crm_aer = (/pcols, crm_nx_rad, crm_ny_rad, crm_nz, ntot_amode/)
#endif

   call phys_getopts( use_ECPP_out = use_ECPP )
   call phys_getopts( use_MMF_VT_out = use_MMF_VT )
   call phys_getopts( MMF_microphysics_scheme_out = MMF_microphysics_scheme )
   call phys_getopts( prog_modal_aero_out = prog_modal_aero )

   if(masterproc) then
      print*,'_____________________________________________________________'
      print*,'____ Multi-Scale Modelling Framework (MMF) Configuration ____'
      print*,'crm_nx     = ',crm_nx
      print*,'crm_ny     = ',crm_ny
      print*,'crm_nz     = ',crm_nz
      print*,'crm_dx     = ',crm_dx
      print*,'crm_dy     = ',crm_dy
      print*,'crm_dt     = ',crm_dt
      print*,'crm_nx_rad = ',crm_nx_rad
      print*,'crm_ny_rad = ',crm_ny_rad
      print*,'use_ECPP   = ',use_ECPP
      print*,'CRM Microphysics = ',MMF_microphysics_scheme
      print*,'_____________________________________________________________'
   end if

   !----------------------------------------------------------------------------
   ! Setup CRM internal parameters
   !----------------------------------------------------------------------------
   call setparm()

   if (use_MMF_VT) then
      ! add variance tracers
      call cnst_add('VT_T', real(0,r8), real(0,r8), real(0,r8), cnst_ind, &
                    longname='VT_T', readiv=.false., mixtype='dry',cam_outfld=.false.)
      call cnst_add('VT_Q', real(0,r8), real(0,r8), real(0,r8), cnst_ind, &
                    longname='VT_Q', readiv=.false., mixtype='dry',cam_outfld=.false.)
   end if

   !----------------------------------------------------------------------------
   ! constituents
   !----------------------------------------------------------------------------
   call cnst_add(cnst_names(1), mwdry, cpair, 0._r8, ixcldliq, longname='Grid box averaged cld liquid amount',is_convtran1=.true.)
   call cnst_add(cnst_names(2), mwdry, cpair, 0._r8, ixcldice, longname='Grid box averaged cld ice amount',   is_convtran1=.true.)
   call cnst_add(cnst_names(3), mwh2o, cpair, 0._r8, ixnumliq, longname='Grid box averaged cld liquid number',is_convtran1=.true.)
   call cnst_add(cnst_names(4), mwh2o, cpair, 0._r8, ixnumice, longname='Grid box averaged cld ice number',   is_convtran1=.true.)
   call cnst_add(cnst_names(5), mwh2o, cpair, 0._r8, ixrain,   longname='Grid box averaged rain amount',      is_convtran1=.true.)
   call cnst_add(cnst_names(6), mwh2o, cpair, 0._r8, ixsnow,   longname='Grid box averaged snow amount',      is_convtran1=.true.)
   call cnst_add(cnst_names(7), mwh2o, cpair, 0._r8, ixnumrain,longname='Grid box averaged rain number',      is_convtran1=.true.)
   call cnst_add(cnst_names(8), mwh2o, cpair, 0._r8, ixnumsnow,longname='Grid box averaged snow number',      is_convtran1=.true.)

   !----------------------------------------------------------------------------
   ! Register MMF history variables
   !----------------------------------------------------------------------------
   call crm_history_register()

   !----------------------------------------------------------------------------
   ! pbuf fields
   !----------------------------------------------------------------------------

   ! CRM state 
   call pbuf_add_field('CRM_U', 'global', dtype_r8, dims_crm_3D, idx)
   call pbuf_add_field('CRM_V', 'global', dtype_r8, dims_crm_3D, idx)
   call pbuf_add_field('CRM_W', 'global', dtype_r8, dims_crm_3D, idx)
   call pbuf_add_field('CRM_T', 'global', dtype_r8, dims_crm_3D, idx)

   ! Radiation
   call pbuf_add_field('CRM_T_RAD',   'physpkg', dtype_r8, dims_crm_rad, idx)
   call pbuf_add_field('CRM_QV_RAD',  'physpkg', dtype_r8, dims_crm_rad, idx)
   call pbuf_add_field('CRM_QC_RAD',  'physpkg', dtype_r8, dims_crm_rad, idx)
   call pbuf_add_field('CRM_QI_RAD',  'physpkg', dtype_r8, dims_crm_rad, idx)
   call pbuf_add_field('CRM_CLD_RAD', 'physpkg', dtype_r8, dims_crm_rad, idx)
   call pbuf_add_field('CRM_QRAD',    'global',  dtype_r8, dims_crm_rad, idx)
   
   call pbuf_add_field('CLDO',  'global', dtype_r8, (/pcols, pver, dyn_time_lvls/), cldo_idx  )
   call pbuf_add_field('CLD',   'global', dtype_r8, (/pcols, pver, dyn_time_lvls/), cld_idx)
   call pbuf_add_field('CONCLD','global', dtype_r8, (/pcols, pver, dyn_time_lvls/), concld_idx)
   
   call pbuf_add_field('QME',   'physpkg',dtype_r8, (/pcols, pver/), qme_idx)
   call pbuf_add_field('PRAIN', 'physpkg',dtype_r8, (/pcols, pver/), prain_idx)
   call pbuf_add_field('NEVAPR',     'physpkg',dtype_r8,(/pcols,pver/), nevapr_idx)
   call pbuf_add_field('PRER_EVAP',  'global', dtype_r8,(/pcols,pver/), prer_evap_idx)

   call pbuf_add_field('WSEDL',      'physpkg',dtype_r8,(/pcols,pver/), wsedl_idx)

   call pbuf_add_field('REI',     'physpkg',dtype_r8,(/pcols,pver/), rei_idx)
   call pbuf_add_field('REL',     'physpkg',dtype_r8,(/pcols,pver/), rel_idx)
   call pbuf_add_field('DEI',     'physpkg',dtype_r8,(/pcols,pver/), dei_idx)        ! Mitchell ice effective diameter for radiation
   call pbuf_add_field('MU',      'physpkg',dtype_r8,(/pcols,pver/), mu_idx)         ! Size distribution shape parameter for radiation
   call pbuf_add_field('LAMBDAC', 'physpkg',dtype_r8,(/pcols,pver/), lambdac_idx)    ! Size distribution shape parameter for radiation
   call pbuf_add_field('ICIWPST', 'physpkg',dtype_r8,(/pcols,pver/), iciwpst_idx)    ! Stratiform only in cloud ice water path for radiation
   call pbuf_add_field('ICLWPST', 'physpkg',dtype_r8,(/pcols,pver/), iclwpst_idx)    ! Stratiform in cloud liquid water path for radiation
   call pbuf_add_field('DES',     'physpkg',dtype_r8,(/pcols,pver/), des_idx)        ! Snow effective diameter for radiation
   call pbuf_add_field('ICIWP',   'physpkg',dtype_r8,(/pcols,pver/), iciwp_idx)      ! In cloud ice water path for radiation
   call pbuf_add_field('ICLWP',   'physpkg',dtype_r8,(/pcols,pver/), iclwp_idx)      ! In cloud liquid water path for radiation
   call pbuf_add_field('ICSWP',   'physpkg',dtype_r8,(/pcols,pver/), icswp_idx)      ! In cloud snow water path for radiation
   call pbuf_add_field('CLDFSNOW','physpkg',dtype_r8,(/pcols,pver,dyn_time_lvls/), cldfsnow_idx) ! Cloud fraction for liquid drops + snow

   if (MMF_microphysics_scheme .eq. 'm2005') then
      call pbuf_add_field('CRM_NC_RAD','physpkg', dtype_r8, dims_crm_rad, idx)
      call pbuf_add_field('CRM_NI_RAD','physpkg', dtype_r8, dims_crm_rad, idx)
      call pbuf_add_field('CRM_QS_RAD','physpkg', dtype_r8, dims_crm_rad, idx)
      call pbuf_add_field('CRM_NS_RAD','physpkg', dtype_r8, dims_crm_rad, idx)

      call pbuf_add_field('CRM_QT', 'global', dtype_r8, dims_crm_3D, idx)
      call pbuf_add_field('CRM_NC', 'global', dtype_r8, dims_crm_3D, idx)
      call pbuf_add_field('CRM_QR', 'global', dtype_r8, dims_crm_3D, idx)
      call pbuf_add_field('CRM_NR', 'global', dtype_r8, dims_crm_3D, idx)
      call pbuf_add_field('CRM_QI', 'global', dtype_r8, dims_crm_3D, idx)
      call pbuf_add_field('CRM_NI', 'global', dtype_r8, dims_crm_3D, idx)
      call pbuf_add_field('CRM_QS', 'global', dtype_r8, dims_crm_3D, idx)
      call pbuf_add_field('CRM_NS', 'global', dtype_r8, dims_crm_3D, idx)
      call pbuf_add_field('CRM_QG', 'global', dtype_r8, dims_crm_3D, idx)
      call pbuf_add_field('CRM_NG', 'global', dtype_r8, dims_crm_3D, idx)
      call pbuf_add_field('CRM_QC', 'global', dtype_r8, dims_crm_3D, idx)

      if (prog_modal_aero) then
         call pbuf_add_field('RATE1_CW2PR_ST','physpkg',dtype_r8,(/pcols,pver/), idx)
      end if

   else
      call pbuf_add_field('CRM_QT', 'global', dtype_r8, dims_crm_3D, idx)
      call pbuf_add_field('CRM_QP', 'global', dtype_r8, dims_crm_3D, idx)
      call pbuf_add_field('CRM_QN', 'global', dtype_r8, dims_crm_3D, idx)
   end if

   ! CRM mass flux
   call pbuf_add_field('MU_CRM', 'physpkg', dtype_r8, dims_gcm_2D, idx) ! mass flux up
   call pbuf_add_field('MD_CRM', 'physpkg', dtype_r8, dims_gcm_2D, idx) ! mass flux down
   call pbuf_add_field('DU_CRM', 'physpkg', dtype_r8, dims_gcm_2D, idx) ! mass detrainment from updraft
   call pbuf_add_field('EU_CRM', 'physpkg', dtype_r8, dims_gcm_2D, idx) ! mass detrainment from updraft
   call pbuf_add_field('ED_CRM', 'physpkg', dtype_r8, dims_gcm_2D, idx) ! mass detrainment from downdraft

   call pbuf_add_field('JT_CRM',    'physpkg', dtype_r8, dims_gcm_1D, idx) ! index of cloud (convection) top for each column
   call pbuf_add_field('MX_CRM',    'physpkg', dtype_r8, dims_gcm_1D, idx) ! index of cloud (convection) bottom for each column
   call pbuf_add_field('IDEEP_CRM', 'physpkg', dtype_r8, dims_gcm_1D, idx) ! Gathering array for convective columns

   ! CRM turbulence
   call pbuf_add_field('TKE_CRM', 'physpkg', dtype_r8, dims_gcm_2D, idx) ! TKE from CRM  (m2/s2)
   call pbuf_add_field('TK_CRM',  'physpkg', dtype_r8, dims_gcm_2D, idx) ! TK from CRM (m2/s)

   ! ACLDY_CEN has to be global in the physcal buffer to be saved in the restart file
   ! total (all sub-classes) cloudy fractional area in previous time step 
   call pbuf_add_field('ACLDY_CEN','global', dtype_r8, dims_gcm_2D, idx) 
   
   ! CRM orientation angle needs to persist across time steps
   call pbuf_add_field('CRM_ANGLE', 'global', dtype_r8, dims_gcm_1D, idx)

   !----------------------------------------------------------------------------
   ! pbuf fields previously added by offline parameterizations
   !----------------------------------------------------------------------------

   call pbuf_add_field('ICWMRDP',    'physpkg',dtype_r8,(/pcols,pver/),icwmrdp_idx)
   call pbuf_add_field('RPRDDP',     'physpkg',dtype_r8,(/pcols,pver/),rprddp_idx)
   call pbuf_add_field('NEVAPR_DPCU','physpkg',dtype_r8,(/pcols,pver/),nevapr_dpcu_idx)
   call pbuf_add_field('PREC_DP',    'physpkg',dtype_r8,(/pcols/),     prec_dp_idx)
   call pbuf_add_field('SNOW_DP',    'physpkg',dtype_r8,(/pcols/),     snow_dp_idx)

   call pbuf_add_field('SH_FRAC', 'physpkg', dtype_r8, (/pcols,pver/), idx) 
   call pbuf_add_field('DP_FRAC', 'physpkg', dtype_r8, (/pcols,pver/), idx) 

   if (use_gw_convect) call pbuf_add_field('TTEND_DP','physpkg',dtype_r8,(/pcols,pver/),ttend_dp_idx)

   call pbuf_add_field('ICWMRSH',    'physpkg',dtype_r8,(/pcols,pver/),icwmrsh_idx )
   call pbuf_add_field('RPRDSH',     'physpkg',dtype_r8,(/pcols,pver/),rprdsh_idx )
   call pbuf_add_field('RPRDTOT',    'physpkg',dtype_r8,(/pcols,pver/),rprdtot_idx )
   call pbuf_add_field('CLDTOP',     'physpkg' ,dtype_r8,(/pcols,1/),  cldtop_idx )
   call pbuf_add_field('CLDBOT',     'physpkg' ,dtype_r8,(/pcols,1/),  cldbot_idx )
   call pbuf_add_field('NEVAPR_SHCU','physpkg',dtype_r8,(/pcols,pver/),nevapr_shcu_idx )
   call pbuf_add_field('PREC_SH',    'physpkg',dtype_r8,(/pcols/),     prec_sh_idx )
   call pbuf_add_field('SNOW_SH',    'physpkg',dtype_r8,(/pcols/),     snow_sh_idx )

   call pbuf_add_field('cush',       'global'  ,dtype_r8,(/pcols,dyn_time_lvls/), idx )  

   call pbuf_add_field('AST',      'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), ast_idx)
   call pbuf_add_field('FICE',     'physpkg', dtype_r8, (/pcols,pver/), fice_idx)

   !----------------------------------------------------------------------------
   !----------------------------------------------------------------------------

end subroutine crm_physics_register

!===================================================================================================
!===================================================================================================

subroutine crm_physics_init(state, pbuf2d, species_class)
!---------------------------------------------------------------------------------------------------
! 
! Purpose: initialize some variables, and add necessary fields into output fields 
!
!---------------------------------------------------------------------------------------------------
   use time_manager,          only: is_first_step
   use physics_buffer,        only: physics_buffer_desc, pbuf_get_index, pbuf_set_field
   use phys_control,          only: phys_getopts
   use crm_history,           only: crm_history_init
   use cam_history,    only: addfld
   use constituents,   only: apcnst, bpcnst, cnst_name, cnst_longname, cnst_get_ind
#ifdef ECPP
   use module_ecpp_ppdriver2, only: papampollu_init
#endif
   use ppgrid,                only: begchunk, endchunk
   use constituents,          only: pcnst, cnst_get_ind
   !----------------------------------------------------------------------------
   ! interface variables
   ! NOTE - species_class is an input so it needs to be outside of ifdef MODAL_AERO for 1-mom micro
   type(physics_state), intent(inout), dimension(begchunk:endchunk) :: state
   type(physics_buffer_desc), pointer :: pbuf2d(:,:)
   integer, intent(inout) :: species_class(:) 
   !----------------------------------------------------------------------------
   ! local variables
   integer :: m, mm
   integer :: ncnst = 8
   integer :: ierror   ! Error code
   logical :: use_ECPP
   logical :: use_MMF_VT
   character(len=16) :: MMF_microphysics_scheme
   integer :: idx_vt_t, idx_vt_q
   integer :: lchnk, ncol
   !----------------------------------------------------------------------------
   call phys_getopts(use_ECPP_out = use_ECPP)
   call phys_getopts(use_MMF_VT_out = use_MMF_VT)
   call phys_getopts(MMF_microphysics_scheme_out = MMF_microphysics_scheme)

#ifdef ECPP
   if (use_ECPP) then
      ! Initialize ECPP driver
      call papampollu_init()
   end if
#endif

   allocate (cldo_save(pcols,pver), stat=ierror)
   if ( ierror /= 0 ) call endrun('CRM_PHYSICS_INIT error: allocation error cldo_save')

   allocate (q_aero(pcols,pver,pcnst), stat=ierror)
   if ( ierror /= 0 ) call endrun('CRM_PHYSICS_INIT error: allocation error q_aero')

#ifdef MODAL_AERO
   allocate (qqcw_save(pcols,pver,pcnst), stat=ierror)
   if ( ierror /= 0 ) call endrun('CRM_PHYSICS_INIT error: allocation error qqcw_save')

   allocate (qqcw_all(pcols,pver,pcnst), stat=ierror)
   if ( ierror /= 0 ) call endrun('CRM_PHYSICS_INIT error: allocation error qqcw_all')

   allocate (dgnumwet_save(pcols,pver,ntot_amode), stat=ierror)
   if ( ierror /= 0 ) call endrun('CRM_PHYSICS_INIT error: allocation error dgnumwet')
#endif

   call crm_history_init(species_class)

   ! Register history variables
   do m = 1, ncnst
      call cnst_get_ind(cnst_names(m), mm)
      if ( any(mm == (/ ixcldliq, ixcldice, ixrain, ixsnow /)) ) then
         ! mass mixing ratios
         call addfld(cnst_name(mm), (/ 'lev' /), 'A', 'kg/kg', cnst_longname(mm)                   )
         ! call addfld(sflxnam(mm),    horiz_only, 'A',   'kg/m2/s', trim(cnst_name(mm))//' surface flux')
      else if ( any(mm == (/ ixnumliq, ixnumice, ixnumrain, ixnumsnow /)) ) then
         ! number concentrations
         call addfld(cnst_name(mm), (/ 'lev' /), 'A', '1/kg', cnst_longname(mm)                   )
         ! call addfld(sflxnam(mm),    horiz_only, 'A',   '1/m2/s', trim(cnst_name(mm))//' surface flux')
      else
         call endrun( "crm_physics_init: Could not call addfld for constituent with unknown units.")
      endif
   end do

   call addfld(apcnst(ixcldliq), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldliq))//' after physics'  )
   call addfld(apcnst(ixcldice), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldice))//' after physics'  )
   call addfld(bpcnst(ixcldliq), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldliq))//' before physics' )
   call addfld(bpcnst(ixcldice), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldice))//' before physics' )
   
   call addfld(apcnst(ixrain), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixrain))//' after physics'  )
   call addfld(apcnst(ixsnow), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixsnow))//' after physics'  )
   call addfld(bpcnst(ixrain), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixrain))//' before physics' )
   call addfld(bpcnst(ixsnow), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixsnow))//' before physics' )

   prec_dp_idx  = pbuf_get_index('PREC_DP')
   snow_dp_idx  = pbuf_get_index('SNOW_DP')
   prec_sh_idx  = pbuf_get_index('PREC_SH')
   snow_sh_idx  = pbuf_get_index('SNOW_SH')
   prec_sed_idx = pbuf_get_index('PREC_SED')
   snow_sed_idx = pbuf_get_index('SNOW_SED')
   snow_str_idx = pbuf_get_index('SNOW_STR')
   prec_pcw_idx = pbuf_get_index('PREC_PCW')
   snow_pcw_idx = pbuf_get_index('SNOW_PCW')

   if (use_MMF_VT) then
      ! initialize variance transport tracers
      call cnst_get_ind( 'VT_T', idx_vt_t )
      call cnst_get_ind( 'VT_Q', idx_vt_q )
      do lchnk = begchunk, endchunk
         ncol  = state(lchnk)%ncol
         state(lchnk)%q(:ncol,:pver,idx_vt_t) = 0
         state(lchnk)%q(:ncol,:pver,idx_vt_q) = 0
      end do
   end if

   if (is_first_step()) then
      call pbuf_set_field(pbuf2d, cldo_idx,   0._r8)
      call pbuf_set_field(pbuf2d, prer_evap_idx,  0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('ICWMRDP')    , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('RPRDDP')     , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('NEVAPR_DPCU'), 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('PREC_DP')    , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('SNOW_DP')    , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('TTEND_DP')   , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('ICWMRSH')    , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('RPRDSH')     , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('RPRDTOT')    , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('NEVAPR_SHCU'), 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('PREC_SH')    , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('SNOW_SH')    , 0._r8)

      call pbuf_set_field(pbuf2d, pbuf_get_index('ICIWPST')    , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('ICLWPST')    , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('ICSWP')      , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('CLDFSNOW')      , 0._r8)
      
   end if

end subroutine crm_physics_init

!===================================================================================================
!===================================================================================================

subroutine crm_physics_tend(ztodt, state, tend, ptend, pbuf, cam_in, cam_out, &
                            species_class, crm_ecpp_output, mmf_clear_rh, &
                            mmf_qchk_prec_dp, mmf_qchk_snow_dp, mmf_rad_flux)

   !------------------------------------------------------------------------------------------------
   !
   ! Purpose: CRM interface for the MMF configuration to update GCM state
   ! 
   ! Original Author: Marat Khairoutdinov
   ! 
   !------------------------------------------------------------------------------------------------
   use ppgrid
   use perf_mod
   use physics_buffer,  only: physics_buffer_desc, pbuf_old_tim_idx, pbuf_get_index, &
                              dyn_time_lvls, pbuf_get_field, pbuf_set_field
   use physics_types,   only: physics_state, physics_tend, physics_ptend, physics_ptend_init
   use camsrfexch,      only: cam_in_t, cam_out_t
   use time_manager,    only: is_first_step, get_nstep
   use crmdims,         only: crm_nx, crm_ny, crm_nz, crm_nx_rad, crm_ny_rad
   use physconst,       only: cpair, latvap, latice, gravit, cappa
   use constituents,    only: pcnst, cnst_get_ind
#if defined(MMF_SAMXX)
   use cpp_interface_mod, only: crm
#elif defined(MMF_SAM)
   use crm_module       , only: crm
#endif
   use params_kind,          only: crm_rknd
   use phys_control,    only: phys_getopts, phys_do_flux_avg
   use crm_history,     only: crm_history_out
   use wv_saturation,   only: qsat_water
#if (defined  m2005 && defined MODAL_AERO)  
   ! modal_aero_data only exists if MODAL_AERO
   use modal_aero_data, only: ntot_amode, ntot_amode
   use ndrop,           only: loadaer
#endif

#if defined( MMF_ORIENT_RAND )
   use RNG_MT ! random number generator for randomly rotating CRM orientation (MMF_ORIENT_RAND)
#endif

   use crm_state_module,       only: crm_state_type
   use crm_rad_module,         only: crm_rad_type, crm_rad_initialize, crm_rad_finalize
   use crm_input_module,       only: crm_input_type
   use crm_output_module,      only: crm_output_type, crm_output_initialize, crm_output_finalize
   use crm_ecpp_output_module, only: crm_ecpp_output_type

   use iso_c_binding, only: c_bool
   use phys_grid    , only: get_rlon_p, get_rlat_p, get_gcol_p  
   use spmd_utils,          only: masterproc

   real(r8),                        intent(in   ) :: ztodt            ! global model time increment
   type(physics_state),             intent(in   ) :: state            ! Global model state 
   type(physics_tend),              intent(in   ) :: tend             ! 
   type(physics_ptend),             intent(  out) :: ptend            ! output tendencies
   type(physics_buffer_desc),       pointer       :: pbuf(:)          ! physics buffer
   type(cam_in_t),                  intent(in   ) :: cam_in           ! atm input from coupler
   type(cam_out_t),                 intent(inout) :: cam_out          ! atm output to coupler
   integer,                         intent(in   ) :: species_class(:) ! aerosol species type
   type(crm_ecpp_output_type),      intent(inout) :: crm_ecpp_output  ! output data for ECPP calculations
   real(r8), dimension(pcols,pver), intent(  out) :: mmf_clear_rh     ! clear air relative humidity used for aerosol water uptake
   real(r8), dimension(pcols),      intent(  out) :: mmf_qchk_prec_dp ! precipitation diagostic (liq+ice)  used for check_energy_chng
   real(r8), dimension(pcols),      intent(  out) :: mmf_qchk_snow_dp ! precipitation diagostic (ice only) used for check_energy_chng
   real(r8), dimension(pcols),      intent(  out) :: mmf_rad_flux     ! radiative flux diagnostic used for check_energy_chng

   !------------------------------------------------------------------------------------------------
   ! Local variables 
   !------------------------------------------------------------------------------------------------

   ! convective precipitation variables
   real(r8), pointer :: prec_dp(:)          ! total precip from deep convection (ZM)    [m/s]
   real(r8), pointer :: snow_dp(:)          ! snow from deep convection (ZM)            [m/s]
   real(r8), pointer :: prec_sh(:)          ! total precip from shallow convection      [m/s]
   real(r8), pointer :: snow_sh(:)          ! snow from shallow convection              [m/s]

   ! stratiform precipitation variables
   real(r8), pointer :: prec_pcw(:)         ! total precip from prognostic cloud scheme [m/s]
   real(r8), pointer :: snow_pcw(:)         ! snow from prognostic cloud scheme         [m/s]
   real(r8), pointer :: prec_sed(:)         ! total precip from cloud sedimentation     [m/s]
   real(r8), pointer :: snow_sed(:)         ! snow from cloud ice sedimentation         [m/s]
   real(r8), pointer :: snow_str(:)         ! snow from stratiform cloud                [m/s]

   integer lchnk                                   ! chunk identifier
   integer ncol                                    ! number of atmospheric columns
   integer nstep                                   ! time steps
   real(r8) crm_run_time                           ! length of CRM integration

   ! physics buffer fields to compute tendencies for stratiform package
   real(r8), pointer, dimension(:,:) :: cld        ! cloud fraction

#if (defined m2005 && defined MODAL_AERO)
   real(r8), dimension(pcols) :: aerosol_num       ! aerosol number concentration      [/m3]
   real(r8), dimension(pcols) :: aerosol_vol       ! aerosol voume concentration       [m3/m3]
   real(r8), dimension(pcols) :: aerosol_hygro     ! aerosol bulk hygroscopicity
   integer  :: phase ! phase of aerosol - interstitial, cloud-borne, or the sum
#endif

   real(r8) :: dp_g                                ! = state%pdel / gravit
   real(r8), dimension(pcols,pver) :: air_density  ! air density                       [kg/m3]
   real(r8), dimension(pcols,pver) :: TKE_tmp      ! temporary TKE value used for ECPP

   ! CRM column radiation stuff:
   real(r8), pointer, dimension(:,:) :: qrs        ! shortwave radiative heating rate
   real(r8), pointer, dimension(:,:) :: qrl        ! shortwave radiative heating rate

   real(r8), dimension(pcols) :: qli_hydro_before  ! column-integraetd rain + snow + graupel 
   real(r8), dimension(pcols) ::  qi_hydro_before  ! column-integrated snow water + graupel water
   real(r8), dimension(pcols) :: qli_hydro_after   ! column-integraetd rain + snow + graupel 
   real(r8), dimension(pcols) ::  qi_hydro_after   ! column-integrated snow water + graupel water
   real(r8) :: sfactor                             ! used to determine precip type for sam1mom

   integer  :: i, icrm, icol, k, m, ii, jj         ! loop iterators
   integer  :: ixcldliq, ixcldice                  ! constituent indices
   integer  :: ixnumliq, ixnumice                  ! constituent indices
   integer  :: ixrain, ixsnow                      ! constituent indices
   integer  :: ixnumrain, ixnumsnow                ! constituent indices
   integer  :: ifld, itim                          ! pbuf field and "old time" indices
   real(r8) :: ideep_crm(pcols)                    ! gathering array for convective columns
   logical  :: ls, lu, lv, lq(pcnst)               ! flags for updating ptend
   logical  :: use_ECPP                            ! flag for ECPP mode
   character(len=16) :: microp_scheme              ! GCM microphysics scheme
   character(len=16) :: MMF_microphysics_scheme    ! CRM microphysics scheme

   logical(c_bool):: use_MMF_VT                    ! flag for MMF variance transport (for C++ CRM)
   logical        :: use_MMF_VT_tmp                ! flag for MMF variance transport (for Fortran CRM)
   integer        :: MMF_VT_wn_max                 ! wavenumber cutoff for filtered variance transport

   real(r8) :: tmp_e_sat                           ! temporary saturation vapor pressure
   real(r8) :: tmp_q_sat                           ! temporary saturation specific humidity
   real(r8) :: tmp_rh_sum                          ! temporary relative humidity sum
   real(r8) :: tmp_rh_cnt                          ! temporary relative humidity count

   ! variables for changing CRM orientation
   real(crm_rknd), parameter        :: pi   = 3.14159265359
   real(crm_rknd), parameter        :: pix2 = 6.28318530718
   real(crm_rknd), dimension(pcols) :: crm_angle

   ! surface flux variables for using adjusted fluxes from flux_avg_run
   real(crm_rknd), pointer, dimension(:) :: shf_tmp
   real(crm_rknd), pointer, dimension(:) :: lhf_tmp
   real(crm_rknd), pointer, dimension(:) :: wsx_tmp
   real(crm_rknd), pointer, dimension(:) :: wsy_tmp

   ! CRM types
   type(crm_state_type)  :: crm_state
   type(crm_rad_type)    :: crm_rad
   type(crm_input_type)  :: crm_input
   type(crm_output_type) :: crm_output

#if defined( MMF_ORIENT_RAND )
   real(crm_rknd) :: unif_rand1           ! uniform random number 
   real(crm_rknd) :: unif_rand2           ! uniform random number 
   real(crm_rknd) :: norm_rand            ! normally distributed random number - Box-Muller (1958)
   real(crm_rknd) :: crm_rotation_std     ! scaling factor for rotation (std dev of rotation angle)
   real(crm_rknd) :: crm_rotation_offset  ! offset to specify preferred rotation direction 
#endif

   real(crm_rknd), allocatable :: crm_clear_rh(:,:) ! clear air relative humidity for aerosol wateruptake

   real(crm_rknd), allocatable :: longitude0(:)
   real(crm_rknd), allocatable :: latitude0 (:)
   integer       , allocatable :: gcolp     (:)
   real(crm_rknd)              :: crm_accel_factor
   logical                     :: use_crm_accel_tmp
   logical                     :: crm_accel_uv_tmp
   logical(c_bool)             :: use_crm_accel
   logical(c_bool)             :: crm_accel_uv
   integer                     :: igstep
   integer                     :: idx_vt_t, idx_vt_q  ! variance transport constituent indices
   !------------------------------------------------------------------------------------------------
   !------------------------------------------------------------------------------------------------
#if defined( MMF_ORIENT_RAND )
   crm_rotation_std    = 20. * pi/180.                 ! std deviation of normal distribution for CRM rotation [radians]
   crm_rotation_offset = 90. * pi/180. * ztodt/86400.  ! This means that a CRM should rotate 90 deg / day on average
#endif

   crm_run_time = ztodt

   call phys_getopts(use_ECPP_out = use_ECPP)
   call phys_getopts(MMF_microphysics_scheme_out = MMF_microphysics_scheme)

   use_MMF_VT = .false.
   call phys_getopts(use_MMF_VT_out = use_MMF_VT_tmp)
   call phys_getopts(MMF_VT_wn_max_out = MMF_VT_wn_max)
   use_MMF_VT = use_MMF_VT_tmp

   nstep = get_nstep()
   lchnk = state%lchnk
   ncol  = state%ncol

   call t_startf ('crm')

   !------------------------------------------------------------------------------------------------
   ! Initialize ptend
   !------------------------------------------------------------------------------------------------
   lu = .true. 
   lv = .true.
   ls = .true.
   lq(:) = .true.
   call physics_ptend_init(ptend, state%psetcols, 'crm', lu=lu, lv=lv, ls=ls, lq=lq)
   
   !------------------------------------------------------------------------------------------------
   ! Initialize CRM state (nullify pointers, allocate memory, etc)
   !------------------------------------------------------------------------------------------------
   call crm_state%initialize()
   call crm_rad_initialize(crm_rad)
   call crm_input%initialize(pcols,pver)
   call crm_output_initialize(crm_output,pcols,pver)

   !------------------------------------------------------------------------------------------------
   ! Set CRM orientation angle
   !------------------------------------------------------------------------------------------------
   crm_angle(:) = 0

#if defined( MMF_ORIENT_RAND )
   !------------------------------------------------------------------------------------------------
   ! Rotate the CRM using a random walk
   !------------------------------------------------------------------------------------------------
   if ( (crm_ny.eq.1) .or. (crm_nx.eq.1) ) then

      if (.not. is_first_step()) then
         ! get current crm angle from pbuf, except on first step
         call pbuf_get_field(pbuf, pbuf_get_index('CRM_ANGLE'), crm_angle)
      end if

      do i = 1,ncol

         ! set the seed based on the chunk and column index (duplicate seeds are ok)
         call RNG_MT_set_seed( lchnk + i + nstep )

         ! Generate a pair of uniform random numbers
         call RNG_MT_gen_rand(unif_rand1)
         call RNG_MT_gen_rand(unif_rand2)

         ! Box-Muller (1958) method of obtaining a Gaussian distributed random number
         norm_rand = sqrt(-2.*log(unif_rand1))*cos(pix2*unif_rand2)
         crm_angle(i) = crm_angle(i) + norm_rand * crm_rotation_std + crm_rotation_offset

         ! Adjust CRM orientation angle to be between 0 and 2*pi
         if ( crm_angle(i) .lt. 0. )   crm_angle(i) = crm_angle(i) + pix2
         if ( crm_angle(i) .gt. pix2 ) crm_angle(i) = crm_angle(i) - pix2

      end do ! i

      ! write current crm_angle to pbuf
      call pbuf_set_field(pbuf, pbuf_get_index('CRM_ANGLE'), crm_angle)
   end if
#else /* MMF_ORIENT_RAND */
   !------------------------------------------------------------------------------------------------
   ! use static CRM orientation (no rotation)
   !------------------------------------------------------------------------------------------------
#if defined( MMF_DIR_NS )
    if (crm_ny.eq.1) crm_angle(:ncol) = pi/2.
#endif /* MMF_DIR_NS */

#endif /* MMF_ORIENT_RAND */

   !------------------------------------------------------------------------------------------------
   ! Retrieve pbuf fields
   !------------------------------------------------------------------------------------------------
   if (MMF_microphysics_scheme .eq. 'm2005') then
      call pbuf_get_field(pbuf, pbuf_get_index('CRM_NC_RAD'), crm_rad%nc, start=(/1,1,1,1/), kount=(/pcols,crm_nx_rad, crm_ny_rad, crm_nz/))
      call pbuf_get_field(pbuf, pbuf_get_index('CRM_NI_RAD'), crm_rad%ni, start=(/1,1,1,1/), kount=(/pcols,crm_nx_rad, crm_ny_rad, crm_nz/))
      call pbuf_get_field(pbuf, pbuf_get_index('CRM_QS_RAD'), crm_rad%qs, start=(/1,1,1,1/), kount=(/pcols,crm_nx_rad, crm_ny_rad, crm_nz/))
      call pbuf_get_field(pbuf, pbuf_get_index('CRM_NS_RAD'), crm_rad%ns, start=(/1,1,1,1/), kount=(/pcols,crm_nx_rad, crm_ny_rad, crm_nz/))
   end if

   call pbuf_get_field (pbuf, pbuf_get_index('CRM_QRAD'),    crm_rad%qrad)
   call pbuf_get_field (pbuf, pbuf_get_index('CRM_T_RAD'),   crm_rad%temperature)
   call pbuf_get_field (pbuf, pbuf_get_index('CRM_QV_RAD'),  crm_rad%qv)
   call pbuf_get_field (pbuf, pbuf_get_index('CRM_QC_RAD'),  crm_rad%qc)
   call pbuf_get_field (pbuf, pbuf_get_index('CRM_QI_RAD'),  crm_rad%qi)
   call pbuf_get_field (pbuf, pbuf_get_index('CRM_CLD_RAD'), crm_rad%cld)

   call pbuf_get_field(pbuf, prec_dp_idx,  prec_dp  )
   call pbuf_get_field(pbuf, prec_sh_idx,  prec_sh  )
   call pbuf_get_field(pbuf, prec_sed_idx, prec_sed )
   call pbuf_get_field(pbuf, prec_pcw_idx, prec_pcw )
   call pbuf_get_field(pbuf, snow_dp_idx,  snow_dp  )
   call pbuf_get_field(pbuf, snow_sh_idx,  snow_sh  )
   call pbuf_get_field(pbuf, snow_sed_idx, snow_sed )
   call pbuf_get_field(pbuf, snow_str_idx, snow_str )
   call pbuf_get_field(pbuf, snow_pcw_idx, snow_pcw )

   !!! total clouds and precipiation - initialize here to be safe 
   !!! WARNING - this disables aerosol scavenging!
   ! call pbuf_set_field(pbuf, pbuf_get_index('AST'   ), 0.0_r8 )
   ! call pbuf_set_field(pbuf, pbuf_get_index('QME'   ), 0.0_r8 )
   ! call pbuf_set_field(pbuf, pbuf_get_index('PRAIN' ), 0.0_r8 )
   ! call pbuf_set_field(pbuf, pbuf_get_index('NEVAPR'), 0.0_r8 )

   ! set convective rain to be zero for PRAIN already includes precipitation production from convection. 
   call pbuf_set_field(pbuf, pbuf_get_index('RPRDTOT'), 0.0_r8 )
   call pbuf_set_field(pbuf, pbuf_get_index('RPRDDP' ), 0.0_r8 )
   call pbuf_set_field(pbuf, pbuf_get_index('RPRDSH' ), 0.0_r8 )
   call pbuf_set_field(pbuf, pbuf_get_index('ICWMRDP'), 0.0_r8 )
   call pbuf_set_field(pbuf, pbuf_get_index('ICWMRSH'), 0.0_r8 )
   
   prec_dp  = 0.
   snow_dp  = 0.
   prec_sh  = 0.
   snow_sh  = 0. 
   prec_sed = 0.
   snow_sed = 0.
   snow_str = 0.
   prec_pcw = 0
   snow_pcw = 0.
   
   ! Initialize stuff:
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)

   !------------------------------------------------------------------------------------------------
   ! Retreive CRM state data from pbuf
   !------------------------------------------------------------------------------------------------

   ! Set pointers from crm_state to fields that persist on physics buffer
   call pbuf_get_field (pbuf, pbuf_get_index('CRM_U'), crm_state%u_wind)
   call pbuf_get_field (pbuf, pbuf_get_index('CRM_V'), crm_state%v_wind)
   call pbuf_get_field (pbuf, pbuf_get_index('CRM_W'), crm_state%w_wind)
   call pbuf_get_field (pbuf, pbuf_get_index('CRM_T'), crm_state%temperature)

   ! Set pointers to microphysics fields in crm_state
   call pbuf_get_field(pbuf, pbuf_get_index('CRM_QT'), crm_state%qt)
   if (MMF_microphysics_scheme .eq. 'sam1mom') then
      call pbuf_get_field(pbuf, pbuf_get_index('CRM_QP'), crm_state%qp)
      call pbuf_get_field(pbuf, pbuf_get_index('CRM_QN'), crm_state%qn)
   else
      call pbuf_get_field(pbuf, pbuf_get_index('CRM_NC'), crm_state%nc)
      call pbuf_get_field(pbuf, pbuf_get_index('CRM_QR'), crm_state%qr)
      call pbuf_get_field(pbuf, pbuf_get_index('CRM_NR'), crm_state%nr)
      call pbuf_get_field(pbuf, pbuf_get_index('CRM_QI'), crm_state%qi)
      call pbuf_get_field(pbuf, pbuf_get_index('CRM_NI'), crm_state%ni)
      call pbuf_get_field(pbuf, pbuf_get_index('CRM_QS'), crm_state%qs)
      call pbuf_get_field(pbuf, pbuf_get_index('CRM_NS'), crm_state%ns)
      call pbuf_get_field(pbuf, pbuf_get_index('CRM_QG'), crm_state%qg)
      call pbuf_get_field(pbuf, pbuf_get_index('CRM_NG'), crm_state%ng)
      call pbuf_get_field(pbuf, pbuf_get_index('CRM_QC'), crm_state%qc)
   end if

   ! "Old" cloud fraction (what does all this mean?)
   itim = pbuf_old_tim_idx()
   ifld = pbuf_get_index('CLD')
   call pbuf_get_field(pbuf, ifld, cld, start=(/1,1,itim/), kount=(/pcols,pver,1/) )

   !------------------------------------------------------------------------------------------------
   !------------------------------------------------------------------------------------------------

   if (is_first_step()) then
      ! call check_energy_timestep_init(state, tend, pbuf)
      ! initialize crm_state%qt to zero (needed for ncol < i <= pcols)
      crm_state%qt(:,:,:,:) = 0.0_r8
      do i = 1,ncol
         do k = 1,crm_nz
            m = pver-k+1

            ! Initialize CRM state
            crm_state%u_wind(i,:,:,k) = state%u(i,m) * cos( crm_angle(i) ) + state%v(i,m) * sin( crm_angle(i) )
            crm_state%v_wind(i,:,:,k) = state%v(i,m) * cos( crm_angle(i) ) - state%u(i,m) * sin( crm_angle(i) )
            crm_state%w_wind(i,:,:,k) = 0.
            crm_state%temperature(i,:,:,k) = state%t(i,m)

            ! Initialize microphysics arrays
            if (MMF_microphysics_scheme .eq. 'sam1mom') then
               crm_state%qt(i,:,:,k) = state%q(i,m,1)+state%q(i,m,ixcldliq)+state%q(i,m,ixcldice)
               crm_state%qp(i,:,:,k) = 0.0_r8
               crm_state%qn(i,:,:,k) = state%q(i,m,ixcldliq)+state%q(i,m,ixcldice)
            else if (MMF_microphysics_scheme .eq. 'm2005') then
               crm_state%qt(i,:,:,k) = state%q(i,m,1)+state%q(i,m,ixcldliq)
               crm_state%qc(i,:,:,k) = state%q(i,m,ixcldliq)
               crm_state%qi(i,:,:,k) = state%q(i,m,ixcldice)
               crm_state%nc(i,:,:,k) = 0.0_r8
               crm_state%qr(i,:,:,k) = 0.0_r8
               crm_state%nr(i,:,:,k) = 0.0_r8
               crm_state%ni(i,:,:,k) = 0.0_r8
               crm_state%qs(i,:,:,k) = 0.0_r8
               crm_state%ns(i,:,:,k) = 0.0_r8
               crm_state%qg(i,:,:,k) = 0.0_r8
               crm_state%ng(i,:,:,k) = 0.0_r8
            end if

         end do
      end do

      do k = 1,crm_nz
         m = pver-k+1
         do i = 1,ncol
            crm_rad%qrad         (i,:,:,k) = 0.
            crm_rad%temperature  (i,:,:,k) = state%t(i,m)
            crm_rad%qv           (i,:,:,k) = state%q(i,m,1)
            crm_rad%qc           (i,:,:,k) = 0.
            crm_rad%qi           (i,:,:,k) = 0.
            crm_rad%cld          (i,:,:,k) = 0.
#ifdef m2005
            if (MMF_microphysics_scheme .eq. 'm2005') then
               crm_rad%nc(i,:,:,k) = 0.0
               crm_rad%ni(i,:,:,k) = 0.0       
               crm_rad%qs(i,:,:,k) = 0.0
               crm_rad%ns(i,:,:,k) = 0.0
            end if
#endif
         end do
      end do

      ! use radiation from grid-cell mean radctl on first time step
      cld(:,:) = 0.
      ptend%s(:,:) = 0.
      ptend%q(:,:,1) = 0.
      ptend%q(:,:,ixcldliq) = 0.
      ptend%q(:,:,ixcldice) = 0.

#ifdef ECPP
      if (use_ECPP) then
         
         air_density(1:ncol,1:pver) = state%pmid(1:ncol,1:pver) / (287.15*state%t(1:ncol,1:pver))

         ! initialize turbulence for ECPP calculations
         call pbuf_set_field(pbuf, pbuf_get_index('TKE_CRM'), 0.0_r8 )
         call pbuf_set_field(pbuf, pbuf_get_index('TK_CRM'), 0.0_r8 )
      end if 
#endif

      ! only need to do this once when crm_angle is static
      call pbuf_set_field(pbuf, pbuf_get_index('CRM_ANGLE'), crm_angle)

      ! Set this output to zero on first step
      mmf_clear_rh(1:ncol,1:pver) = 0

   else  ! not is_first_step

      ptend%s(:,:) = 0.    ! necessary?
      ptend%q(:,:,1) = 0.  ! necessary?
      ptend%q(:,:,ixcldliq) = 0.
      ptend%q(:,:,ixcldice) = 0.

      do m = 1,crm_nz
         k = pver-m+1
         do i = 1,ncol
            crm_rad%qrad(i,:,:,m) = crm_rad%qrad(i,:,:,m) / state%pdel(i,k) ! for energy conservation
         end do
      end do

#if (defined m2005 && defined MODAL_AERO)
      air_density(1:ncol,1:pver) = state%pmid(1:ncol,1:pver) / (287.15*state%t(1:ncol,1:pver))
#endif

      !---------------------------------------------------------------------------------------------
      ! calculate total water before calling crm - used for check_energy_chng() after CRM
      !---------------------------------------------------------------------------------------------
      do i = 1,ncol
         qli_hydro_before(i) = 0.0_r8
         qi_hydro_before(i) = 0.0_r8

         do m = 1,crm_nz
            k = pver-m+1
            dp_g = state%pdel(i,k)/gravit
            do jj = 1,crm_ny
               do ii = 1,crm_nx
                  if (MMF_microphysics_scheme .eq. 'm2005') then
                     qli_hydro_before(i) = qli_hydro_before(i)+(crm_state%qr(i,ii,jj,m)+ &
                                                                crm_state%qs(i,ii,jj,m)+ &
                                                                crm_state%qg(i,ii,jj,m)) * dp_g
                     qi_hydro_before(i)  =  qi_hydro_before(i)+(crm_state%qs(i,ii,jj,m)+ &
                                                                crm_state%qg(i,ii,jj,m)) * dp_g
                  else if (MMF_microphysics_scheme .eq. 'sam1mom') then
                     sfactor = max(0._r8,min(1._r8,(crm_state%temperature(i,ii,jj,m)-268.16)*1./(283.16-268.16)))
                     qli_hydro_before(i) = qli_hydro_before(i)+crm_state%qp(i,ii,jj,m) * dp_g
                     qi_hydro_before(i)  =  qi_hydro_before(i)+crm_state%qp(i,ii,jj,m) * (1-sfactor) * dp_g
                  end if ! MMF_microphysics_scheme
               end do ! ii
            end do ! jj
         end do ! m

         qli_hydro_before(i) = qli_hydro_before(i)/(crm_nx*crm_ny)
         qi_hydro_before(i)  =  qi_hydro_before(i)/(crm_nx*crm_ny)
      end do ! i = 1,ncol

      ! Set CRM inputs
      ! TODO: move this to a routine and call like:
      !    call set_crm_input(state, cam_in, pbuf, crm_input)
      crm_input%zmid(1:ncol,1:pver) = state%zm(1:ncol,1:pver)
      crm_input%zint(1:ncol,1:pver+1) = state%zi(1:ncol,1:pver+1)
      crm_input%tl(1:ncol,1:pver) = state%t(1:ncol,1:pver)
      crm_input%ql(1:ncol,1:pver) = state%q(1:ncol,1:pver,1)
      crm_input%qccl(1:ncol,1:pver) = state%q(1:ncol,1:pver,ixcldliq)
      crm_input%qiil(1:ncol,1:pver) = state%q(1:ncol,1:pver,ixcldice)
      crm_input%ps(1:ncol) = state%ps(1:ncol)
      crm_input%pmid(1:ncol,1:pver) = state%pmid(1:ncol,1:pver)
      crm_input%pint(1:ncol,1:pver+1) = state%pint(1:ncol,1:pver+1)
      crm_input%pdel(1:ncol,1:pver) = state%pdel(1:ncol,1:pver)
      crm_input%phis(1:ncol) = state%phis(1:ncol)
      crm_input%ul(1:ncol,1:pver) = state%u(1:ncol,1:pver)
      crm_input%vl(1:ncol,1:pver) = state%v(1:ncol,1:pver)
      crm_input%ocnfrac(1:ncol) = cam_in%ocnfrac(1:ncol)
#if defined( MMF_ESMT )
      ! Set the input wind for ESMT
      crm_input%ul_esmt(1:ncol,1:pver) = state%u(1:ncol,1:pver)
      crm_input%vl_esmt(1:ncol,1:pver) = state%v(1:ncol,1:pver)
#endif /* MMF_ESMT */

      ! Set surface flux variables
      if (phys_do_flux_avg()) then
         call pbuf_get_field(pbuf, pbuf_get_index('LHFLX'), shf_tmp)
         call pbuf_get_field(pbuf, pbuf_get_index('SHFLX'), lhf_tmp)
         call pbuf_get_field(pbuf, pbuf_get_index('TAUX'),  wsx_tmp)
         call pbuf_get_field(pbuf, pbuf_get_index('TAUY'),  wsy_tmp)
      else
         shf_tmp = cam_in%shf
         lhf_tmp = cam_in%lhf
         wsx_tmp = cam_in%wsx
         wsy_tmp = cam_in%wsy
      end if
      do i = 1,ncol
         crm_input%tau00(i)   = sqrt(wsx_tmp(i)**2 + wsy_tmp(i)**2)
         crm_input%bflxls(i)  = shf_tmp(i)/cpair + 0.61*state%t(i,pver)*lhf_tmp(i)/latvap
         crm_input%fluxu00(i) = wsx_tmp(i)         ! N/m2
         crm_input%fluxv00(i) = wsy_tmp(i)         ! N/m2
         crm_input%fluxt00(i) = shf_tmp(i)/cpair   ! K Kg/ (m2 s)
         crm_input%fluxq00(i) = lhf_tmp(i)/latvap  ! Kg/(m2 s)
         crm_input%wndls(i)   = sqrt(state%u(i,pver)**2 + state%v(i,pver)**2)
      end do

#if (defined m2005 && defined MODAL_AERO)
      ! Set aerosol
      phase = 1  ! interstital aerosols only
      do i = 1,ncol
         do k = 1, pver
            do m = 1, ntot_amode
               call loadaer( state, pbuf, i, i, k, m, air_density, phase, &
                             aerosol_num, aerosol_vol, aerosol_hygro)
               crm_input%naermod (i,k,m) = aerosol_num(i)
               crm_input%vaerosol(i,k,m) = aerosol_vol(i)
               crm_input%hygro   (i,k,m) = aerosol_hygro(i)
            end do    
         end do
      end do
#endif
      !---------------------------------------------------------------------------------------------
      ! Variance transport
      !---------------------------------------------------------------------------------------------
      if (use_MMF_VT) then
         call cnst_get_ind( 'VT_T', idx_vt_t )
         call cnst_get_ind( 'VT_Q', idx_vt_q )
         crm_input%t_vt(:ncol,:pver) = state%q(:ncol,:pver,idx_vt_t)
         crm_input%q_vt(:ncol,:pver) = state%q(:ncol,:pver,idx_vt_q)
      end if
      !---------------------------------------------------------------------------------------------
      ! Set the input wind (also sets CRM orientation)
      !---------------------------------------------------------------------------------------------
      do i = 1,ncol
         do k = 1,pver
            crm_input%ul(i,k) = state%u(i,k) * cos( crm_angle(i) ) + state%v(i,k) * sin( crm_angle(i) )
            crm_input%vl(i,k) = state%v(i,k) * cos( crm_angle(i) ) - state%u(i,k) * sin( crm_angle(i) )
         end do ! k=1,pver
      end do ! i=1,ncol

      !---------------------------------------------------------------------------------------------
      ! Run the CRM
      !---------------------------------------------------------------------------------------------
      if (.not.allocated(ptend%q)) write(*,*) '=== ptend%q not allocated ==='
      if (.not.allocated(ptend%s)) write(*,*) '=== ptend%s not allocated ==='

      if (.not.allocated(crm_clear_rh)) allocate(crm_clear_rh(ncol,crm_nz))

      ! Load latitude, longitude, and unique column ID for all CRMs
      allocate(longitude0(ncol))
      allocate(latitude0 (ncol))
      allocate(gcolp     (ncol))
      do icrm = 1 , ncol
        latitude0 (icrm) = get_rlat_p(lchnk,icrm) * 57.296_r8
        longitude0(icrm) = get_rlon_p(lchnk,icrm) * 57.296_r8
        gcolp     (icrm) = get_gcol_p(lchnk,icrm)
      enddo
      ! set CRM mean state acceleration (MSA) parameters from namelist
      use_crm_accel = .false.
      crm_accel_factor = 0.
      crm_accel_uv = .false.
      call phys_getopts(use_crm_accel_out    = use_crm_accel_tmp, &
                        crm_accel_factor_out = crm_accel_factor, &
                        crm_accel_uv_out     = crm_accel_uv_tmp)
      use_crm_accel = use_crm_accel_tmp
      crm_accel_uv = crm_accel_uv_tmp

      if (masterproc) then
        if (use_crm_accel) then
#if !defined(sam1mom)
          write(0,*) "CRM time step relaxation is only compatible with sam1mom microphysics"
          call endrun('crm main')
#endif
        endif
      endif

      ! Load the nstep
      igstep = get_nstep()

#if defined(MMF_SAM)

      call t_startf ('crm_call')

      call crm(lchnk, ncol, ztodt, pver, &
               crm_input, crm_state, crm_rad, &
               crm_ecpp_output, crm_output, crm_clear_rh, &
               latitude0, longitude0, gcolp, igstep, &
               use_MMF_VT_tmp, MMF_VT_wn_max, &
               use_crm_accel_tmp, crm_accel_factor, crm_accel_uv_tmp)

      call t_stopf('crm_call')

#elif defined(MMF_SAMXX)

      call t_startf ('crm_call')

      ! Fortran classes don't translate to C++ classes, we we have to separate
      ! this stuff out when calling the C++ routinte crm(...)
      call crm(ncol, pcols, ztodt, pver, crm_input%bflxls, crm_input%wndls, crm_input%zmid, crm_input%zint, &
               crm_input%pmid, crm_input%pint, crm_input%pdel, crm_input%ul, crm_input%vl, &
               crm_input%tl, crm_input%qccl, crm_input%qiil, crm_input%ql, crm_input%tau00, &
#ifdef MMF_ESMT
               crm_input%ul_esmt, crm_input%vl_esmt,                                        &
#endif
               crm_input%t_vt, crm_input%q_vt, &
               crm_state%u_wind, crm_state%v_wind, crm_state%w_wind, crm_state%temperature, &
               crm_state%qt, crm_state%qp, crm_state%qn, crm_rad%qrad, crm_rad%temperature, &
               crm_rad%qv, crm_rad%qc, crm_rad%qi, crm_rad%cld, crm_output%subcycle_factor, &
               crm_output%prectend, crm_output%precstend, crm_output%cld, crm_output%cldtop, &
               crm_output%gicewp, crm_output%gliqwp, crm_output%mctot, crm_output%mcup, crm_output%mcdn, &
               crm_output%mcuup, crm_output%mcudn, crm_output%qc_mean, crm_output%qi_mean, crm_output%qs_mean, &
               crm_output%qg_mean, crm_output%qr_mean, crm_output%mu_crm, crm_output%md_crm, crm_output%eu_crm, &
               crm_output%du_crm, crm_output%ed_crm, crm_output%flux_qt, crm_output%flux_u, crm_output%flux_v, &
               crm_output%fluxsgs_qt, crm_output%tkez, crm_output%tkew, crm_output%tkesgsz, crm_output%tkz, crm_output%flux_qp, &
               crm_output%precflux, crm_output%qt_trans, crm_output%qp_trans, crm_output%qp_fall, crm_output%qp_evp, &
               crm_output%qp_src, crm_output%qt_ls, crm_output%t_ls, crm_output%jt_crm, crm_output%mx_crm, crm_output%cltot, &
               crm_output%clhgh, crm_output%clmed, crm_output%cllow, &
               crm_output%sltend, crm_output%qltend, crm_output%qcltend, crm_output%qiltend, &
               crm_output%t_vt_tend, crm_output%q_vt_tend, crm_output%t_vt_ls, crm_output%q_vt_ls, &
#if defined(MMF_MOMENTUM_FEEDBACK)
               crm_output%ultend, crm_output%vltend, &
#endif /* MMF_MOMENTUM_FEEDBACK */
               crm_output%tk, crm_output%tkh, crm_output%qcl, crm_output%qci, crm_output%qpl, crm_output%qpi, &
               crm_output%z0m, crm_output%taux, crm_output%tauy, crm_output%precc, crm_output%precl, crm_output%precsc, &
               crm_output%precsl, crm_output%prec_crm,                         &
#ifdef MMF_ESMT
               crm_output%u_tend_esmt, crm_output%v_tend_esmt,                 &
#endif
               crm_clear_rh, &
               latitude0, longitude0, gcolp, igstep, &
               use_MMF_VT, MMF_VT_wn_max, &
               use_crm_accel, crm_accel_factor, crm_accel_uv)

      call t_stopf('crm_call')
      
#endif

      deallocate(longitude0)
      deallocate(latitude0 )
      deallocate(gcolp     )

      !---------------------------------------------------------------------------------------------
      ! Copy tendencies from CRM output to ptend
      !---------------------------------------------------------------------------------------------
      ptend%s(:ncol,:pver)          = crm_output%sltend(1:ncol,1:pver)
      ptend%q(:ncol,:pver,1)        = crm_output%qltend(1:ncol,1:pver)
      ptend%q(:ncol,:pver,ixcldliq) = crm_output%qcltend(1:ncol,1:pver)
      ptend%q(:ncol,:pver,ixcldice) = crm_output%qiltend(1:ncol,1:pver)

      if (use_MMF_VT) then
         call cnst_get_ind( 'VT_T', idx_vt_t )
         call cnst_get_ind( 'VT_Q', idx_vt_q )
         ptend%q(1:ncol,1:pver,idx_vt_t) = crm_output%t_vt_tend(1:ncol,1:pver)
         ptend%q(1:ncol,1:pver,idx_vt_q) = crm_output%q_vt_tend(1:ncol,1:pver)
      end if

      !---------------------------------------------------------------------------------------------
      ! Add radiative heating tendency above CRM
      !---------------------------------------------------------------------------------------------

      call pbuf_get_field(pbuf, pbuf_get_index('QRL'), qrl)
      call pbuf_get_field(pbuf, pbuf_get_index('QRS'), qrs)

      do k = 1,pver
         do i = 1,ncol
            qrs(i,k) = qrs(i,k)/state%pdel(i,k)
            qrl(i,k) = qrl(i,k)/state%pdel(i,k)
         end do
      end do

      ! The radiation tendencies in the GCM levels above the CRM and the top 2 CRM levels are set to
      ! be zero in the CRM, So add radiation tendencies to these levels 
      ptend%s(:ncol, :pver-crm_nz+2) = qrs(:ncol,:pver-crm_nz+2) + qrl(:ncol,:pver-crm_nz+2)

      ! This will be used to check energy conservation
      mmf_rad_flux(:ncol) = 0.0_r8
      do k = 1,pver
         do i = 1,ncol
            mmf_rad_flux(i) = mmf_rad_flux(i) + ( qrs(i,k) + qrl(i,k) ) * state%pdel(i,k)/gravit
         end do
      end do

      !---------------------------------------------------------------------------------------------
      ! CRM cloud/precip output
      !---------------------------------------------------------------------------------------------
      ! We need to do this here because we did not set crm_output%cld as a 
      ! pointer to pbuf, so we need to copy the data over. NOTE: I think this 
      ! can be done using pbuf_set_field without making an extra pointer for 
      ! cld, but I do not think we would be able to zero-out the rest of cld 
      ! beyond pcols that way.
      !call pbuf_get_field(pbuf, pbuf_get_index('CLD'), cld)
      !cld(:,:) = 0
      cld(1:ncol,1:pver) = crm_output%cld(1:ncol,1:pver)

      ! There is no separate convective and stratiform precip for CRM,
      ! so make it all convective and zero out the stratform
      crm_output%precc(:ncol)  = crm_output%precc(:ncol) + crm_output%precl(:ncol)
      crm_output%precsc(:ncol) = crm_output%precsc(:ncol) + crm_output%precsl(:ncol)
      crm_output%precl(:ncol)  = 0
      crm_output%precsl(:ncol) = 0

      !---------------------------------------------------------------------------------------------
      ! Set pbuf variables needed by the coupler
      !---------------------------------------------------------------------------------------------

      ! TODO: do we need to zero out the elements beyond ncol here?
      prec_dp(1:ncol) = crm_output%precc(1:ncol)
      snow_dp(1:ncol) = crm_output%precsc(1:ncol)

      !---------------------------------------------------------------------------------------------
      ! Output for ECPP
      !---------------------------------------------------------------------------------------------
      if (use_ECPP) then

         ! turbulence
         air_density(1:ncol,1:pver) = state%pmid(1:ncol,1:pver) / (287.15*state%t(1:ncol,1:pver))
         TKE_tmp(1:ncol,1:pver) = crm_output%tkez(1:ncol,1:pver) / air_density(1:ncol,1:pver)
         call pbuf_set_field(pbuf, pbuf_get_index('TKE_CRM'), TKE_tmp )
         call pbuf_set_field(pbuf, pbuf_get_index('TK_CRM'), crm_output%tkz )

         ! For convective transport
         do i = 1,ncol
            ideep_crm(i) = i*1.0 
         end do
         call pbuf_set_field(pbuf, pbuf_get_index('MU_CRM'), crm_output%mu_crm )
         call pbuf_set_field(pbuf, pbuf_get_index('MD_CRM'), crm_output%md_crm )
         call pbuf_set_field(pbuf, pbuf_get_index('EU_CRM'), crm_output%eu_crm )
         call pbuf_set_field(pbuf, pbuf_get_index('DU_CRM'), crm_output%du_crm )
         call pbuf_set_field(pbuf, pbuf_get_index('ED_CRM'), crm_output%eu_crm )
         call pbuf_set_field(pbuf, pbuf_get_index('JT_CRM'), crm_output%jt_crm )
         call pbuf_set_field(pbuf, pbuf_get_index('MX_CRM'), crm_output%mx_crm )
         call pbuf_set_field(pbuf, pbuf_get_index('IDEEP_CRM'), ideep_crm )

         crm_ecpp_output%qlsinkcen = crm_ecpp_output%qlsink_avgcen
      
      end if ! use_ECPP

      !---------------------------------------------------------------------------------------------
      ! Set ptend logicals with physics tendencies from CRM
      !---------------------------------------------------------------------------------------------
      ptend%name = 'crm'
      ptend%ls           = .TRUE.
      ptend%lq(1)        = .TRUE.
      ptend%lq(ixcldliq) = .TRUE.
      ptend%lq(ixcldice) = .TRUE.
      ptend%lu           = .FALSE.
      ptend%lv           = .FALSE.

      if (use_MMF_VT) then
         call cnst_get_ind( 'VT_T', idx_vt_t )
         call cnst_get_ind( 'VT_Q', idx_vt_q )
         ptend%lq(idx_vt_t) = .TRUE.
         ptend%lq(idx_vt_q) = .TRUE.
      end if

      !---------------------------------------------------------------------------------------------
      ! CRM momentum tendencies
      !---------------------------------------------------------------------------------------------

#if defined( MMF_USE_ESMT )
      ptend%lu = .TRUE.
      ptend%lv = .TRUE.
      ptend%u  = crm_output%u_tend_esmt
      ptend%v  = crm_output%v_tend_esmt
#else /* MMF_USE_ESMT not defined */  

#if defined(MMF_MOMENTUM_FEEDBACK)
      ptend%lu = .TRUE.
      ptend%lv = .TRUE.
      ! rotate resolved CRM momentum tendencies back
      do i = 1, ncol 
         ptend%u(i,1:pver) = crm_output%ultend(i,1:pver) * cos( -1.*crm_angle(i) ) + crm_output%vltend(i,1:pver) * sin( -1.*crm_angle(i) )
         ptend%v(i,1:pver) = crm_output%vltend(i,1:pver) * cos( -1.*crm_angle(i) ) - crm_output%ultend(i,1:pver) * sin( -1.*crm_angle(i) )
      end do
#endif /* MMF_MOMENTUM_FEEDBACK */

#endif /* MMF_USE_ESMT */

      !---------------------------------------------------------------------------------------------
      ! Write out data for history files
      !---------------------------------------------------------------------------------------------

      call crm_history_out(state, ptend, crm_state, crm_rad, crm_output, crm_ecpp_output, qrs, qrl)

      ! Convert heating rate to Q*dp to conserve energy across timesteps
      do m = 1,crm_nz
         k = pver-m+1
         do i = 1,ncol
            crm_rad%qrad(i,:,:,m) = crm_rad%qrad(i,:,:,m) * state%pdel(i,k) ! for energy conservation
         end do
      end do

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

      call phys_getopts(microp_scheme_out=microp_scheme)
      if (microp_scheme .eq. 'MG' ) then
#ifdef m2005
         if (MMF_microphysics_scheme .eq. 'm2005') then
            call cnst_get_ind('NUMLIQ', ixnumliq)
            call cnst_get_ind('NUMICE', ixnumice)
            call cnst_get_ind('RAINQM', ixrain)   
            call cnst_get_ind('SNOWQM', ixsnow)   
            call cnst_get_ind('NUMRAI', ixnumrain)
            call cnst_get_ind('NUMSNO', ixnumsnow)
            ptend%lq(ixnumliq)  = .TRUE.
            ptend%lq(ixnumice)  = .TRUE.
            ptend%q(:,:,ixnumliq)  = 0._r8
            ptend%q(:,:,ixnumice)  = 0._r8
#ifdef ECPP
            ptend%lq(ixrain)    = .TRUE. 
            ptend%lq(ixsnow)    = .TRUE. 
            ptend%lq(ixnumrain) = .TRUE. 
            ptend%lq(ixnumsnow) = .TRUE. 
            ptend%q(:,:,ixrain)    = 0._r8 
            ptend%q(:,:,ixsnow)    = 0._r8 
            ptend%q(:,:,ixnumrain) = 0._r8 
            ptend%q(:,:,ixnumsnow) = 0._r8 
#endif /* ECPP */

            do i = 1, ncol
            do k = 1, crm_nz 
               m = pver-k+1
               do ii = 1, crm_nx
               do jj = 1, crm_ny
                 ptend%q(i,m,ixnumliq)  = ptend%q(i,m,ixnumliq)  + crm_state%nc(i,ii,jj,k) 
                 ptend%q(i,m,ixnumice)  = ptend%q(i,m,ixnumice)  + crm_state%ni(i,ii,jj,k)
#ifdef ECPP
                 ptend%q(i,m,ixrain)    = ptend%q(i,m,ixrain)    + crm_state%qr(i,ii,jj,k)
                 ptend%q(i,m,ixsnow)    = ptend%q(i,m,ixsnow)    + crm_state%qs(i,ii,jj,k)
                 ptend%q(i,m,ixnumrain) = ptend%q(i,m,ixnumrain) + crm_state%nr(i,ii,jj,k)
                 ptend%q(i,m,ixnumsnow) = ptend%q(i,m,ixnumsnow) + crm_state%ns(i,ii,jj,k)
#endif /* ECPP */
               end do
               end do
               ptend%q(i,m,ixnumliq)  = (ptend%q(i,m,ixnumliq) /(crm_nx*crm_ny) - state%q(i,m,ixnumliq)) /crm_run_time
               ptend%q(i,m,ixnumice)  = (ptend%q(i,m,ixnumice) /(crm_nx*crm_ny) - state%q(i,m,ixnumice)) /crm_run_time
#ifdef ECPP
               ptend%q(i,m,ixrain)    = (ptend%q(i,m,ixrain)   /(crm_nx*crm_ny) - state%q(i,m,ixrain))   /crm_run_time
               ptend%q(i,m,ixsnow)    = (ptend%q(i,m,ixsnow)   /(crm_nx*crm_ny) - state%q(i,m,ixsnow))   /crm_run_time
               ptend%q(i,m,ixnumrain) = (ptend%q(i,m,ixnumrain)/(crm_nx*crm_ny) - state%q(i,m,ixnumrain))/crm_run_time
               ptend%q(i,m,ixnumsnow) = (ptend%q(i,m,ixnumsnow)/(crm_nx*crm_ny) - state%q(i,m,ixnumsnow))/crm_run_time
#endif /* ECPP */

            end do
            end do
         end if
#endif /* m2005 */
      end if ! microp_scheme .eq. 'MG'

      !---------------------------------------------------------------------------------------------
      ! calculate column integrated water for energy check
      !---------------------------------------------------------------------------------------------
      do i = 1,ncol
         qli_hydro_after(i) = 0.0_r8
         qi_hydro_after(i) = 0.0_r8
         do m = 1,crm_nz
            k = pver-m+1
            dp_g = state%pdel(i,k)/gravit
            do jj = 1,crm_ny
               do ii = 1,crm_nx
                  if(MMF_microphysics_scheme .eq. 'm2005') then
                     qli_hydro_after(i) = qli_hydro_after(i)+(crm_state%qr(i,ii,jj,m)+ &
                                                              crm_state%qs(i,ii,jj,m)+ &
                                                              crm_state%qg(i,ii,jj,m)) * dp_g
                     qi_hydro_after(i)  =  qi_hydro_after(i)+(crm_state%qs(i,ii,jj,m)+ &
                                                              crm_state%qg(i,ii,jj,m)) * dp_g
                  else if(MMF_microphysics_scheme .eq. 'sam1mom') then 
                     sfactor = max(0._r8,min(1._r8,(crm_state%temperature(i,ii,jj,m)-268.16)*1./(283.16-268.16)))
                     qli_hydro_after(i) = qli_hydro_after(i)+crm_state%qp(i,ii,jj,m) * dp_g
                     qi_hydro_after(i)  =  qi_hydro_after(i)+crm_state%qp(i,ii,jj,m) * (1-sfactor) * dp_g
                  end if ! MMF_microphysics_scheme
               end do ! ii
            end do ! jj
         end do ! m = 1,crm_nz
         qli_hydro_after(i) = qli_hydro_after(i)/(crm_nx*crm_ny)
         qi_hydro_after(i)  =  qi_hydro_after(i)/(crm_nx*crm_ny)
      end do ! i = 1,ncold

      mmf_qchk_prec_dp(:ncol) = prec_dp(:ncol) + (qli_hydro_after (:ncol) - &
                                                  qli_hydro_before(:ncol))/crm_run_time/1000._r8
      mmf_qchk_snow_dp(:ncol) = snow_dp(:ncol) + ( qi_hydro_after (:ncol) - &
                                                   qi_hydro_before(:ncol))/crm_run_time/1000._r8

      !---------------------------------------------------------------------------------------------
      ! copy clear air relative humdity for aerosol water uptake
      !---------------------------------------------------------------------------------------------
      ! initialize to zero, so no aerosol water uptake occurs by default
      mmf_clear_rh(1:ncol,1:pver) = 0
      do icol = 1,ncol
         do m = 1,crm_nz
            k = pver-m+1
            mmf_clear_rh(icol,k) = crm_clear_rh(icol,m)
         end do ! m = 1,crm_nz
      end do ! i = 1,ncol
      deallocate(crm_clear_rh)
      !---------------------------------------------------------------------------------------------
      !---------------------------------------------------------------------------------------------

   end if ! (is_first_step())

   call t_stopf('crm')
   
   !------------------------------------------------------------------------------------------------
   ! Free memory in derived types
   !------------------------------------------------------------------------------------------------

   call crm_state%finalize()
   call crm_rad_finalize(crm_rad)
   call crm_input%finalize()
   call crm_output_finalize(crm_output)

end subroutine crm_physics_tend

!==================================================================================================
!==================================================================================================

subroutine crm_surface_flux_bypass_tend(state, cam_in, ptend)
   !------------------------------------------------------------------------------------------------
   ! This subroutine is used to apply the fluxes when MMF_FLUX_BYPASS is used.
   ! The surface flux bypass option was originally used by Mike Pritchard (UCI)
   ! Without this bypass the surface flux tendencies are applied to the lowest 
   ! layer of the GCM state without being diffused vertically by PBL turbulence. 
   ! This was intentional (confirmed by Marat). This bypass applies the surface 
   ! fluxes after the dycor and prior to running the CRM (the tendency addition 
   ! in diffusion_solver.F90 is disabled). This is a more natural progression 
   ! and does not expose the GCM dynamical core to unrealistic gradients.
   ! (only sensible and latent heat fluxes are affected)
   !------------------------------------------------------------------------------------------------
   use physics_types,   only: physics_state, physics_ptend, physics_ptend_init
   use physics_buffer,  only: physics_buffer_desc
   use camsrfexch,      only: cam_in_t
   use ppgrid,          only: begchunk, endchunk, pcols, pver
   use constituents,    only: pcnst
   use physconst,       only: gravit

   implicit none

   type(physics_state), intent(in   ) :: state
   type(cam_in_t),      intent(in   ) :: cam_in
   type(physics_ptend), intent(  out) :: ptend 

   integer  :: icol     ! loop iterator
   integer  :: ncol     ! number of columns
   real(r8) :: g_dp     ! temporary variable for unit conversion
   logical, dimension(pcnst) :: lq

   ncol  = state%ncol

   ! initialize ptend
   lq(:) = .false.
   lq(1) = .true.
   call physics_ptend_init(ptend, state%psetcols, 'MMF_FLUX_BYPASS', &
                           lu=.false., lv=.false., ls=.true., lq=lq)

   ! apply fluxes to bottom layer
   do icol = 1,ncol
      g_dp = gravit * state%rpdel(icol,pver)             ! note : rpdel = 1./pdel
      ! ptend%s(icol,:)   = 0.
      ! ptend%q(icol,:,1) = 0.
      ptend%s(icol,pver)   = g_dp * cam_in%shf(icol)
      ptend%q(icol,pver,1) = g_dp * cam_in%cflx(icol,1)
   end do

end subroutine crm_surface_flux_bypass_tend

!==================================================================================================
!==================================================================================================

subroutine crm_save_state_tend(state,tend,pbuf)
   !------------------------------------------------------------------------------------------------
   ! This subroutine is used to save state variables at the beginning of tphysbc so they can be 
   ! recalled after they have been changed by conventional physics.
   !------------------------------------------------------------------------------------------------
   use physics_types,   only: physics_state, physics_tend, physics_tend_dealloc, &
                              physics_state_copy, physics_tend_copy, physics_state_dealloc
   use time_manager,    only: is_first_step
   use physics_buffer,  only: pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field, physics_buffer_desc
   use phys_control,    only: phys_getopts
#ifdef MODAL_AERO
   use modal_aero_data, only: ntot_amode, qqcw_get_field
#endif

   implicit none

   type(physics_state),       intent(in   ) :: state
   type(physics_tend),        intent(in   ) :: tend 
   type(physics_buffer_desc), pointer       :: pbuf(:)

   integer itim, ifld, ncol, i, lchnk
   real(r8), pointer, dimension(:,:) :: cld        ! cloud fraction
   real(r8), pointer, dimension(:,:) :: qqcw
   logical :: use_ECPP

   call phys_getopts( use_ECPP_out  = use_ECPP  )

   lchnk = state%lchnk
   ncol  = state%ncol
   itim  = pbuf_old_tim_idx()

   ! Save the state and tend variables 
   ! Overwrite conventional physics effects before calling the crm. 
   ! Non-CRM physics routines are allowed to compute diagnostic tendencies.
   ! Note that state_save and tend_save get allocated in the copy routines
   if ( allocated(state_save%t) ) call physics_state_dealloc(state_save)
   if ( allocated(tend_save%dtdt) ) call physics_tend_dealloc(tend_save)

   call physics_state_copy(state,state_save)
   call physics_tend_copy(tend,tend_save)
   
   ! save the old cloud fraction
   call pbuf_get_field(pbuf, cldo_idx, cldo, start=(/1,1,itim/), kount=(/pcols,pver,1/) )
   cldo_save(:ncol, :) = cldo(:ncol, :)

   ! other things relevant for aerosols
   if (use_ECPP) then
      call pbuf_get_field(pbuf, pbuf_get_index('CLD'), cld )
      call pbuf_get_field(pbuf, pbuf_get_index('ACLDY_CEN'), acldy_cen_tbeg)
      if(is_first_step())then
        acldy_cen_tbeg(:ncol,:) = cld(:ncol,:)
      end if
   end if

#if ( defined MODAL_AERO )
   qqcw_all=0_r8
   do i = 1,pcnst
      qqcw   =>  qqcw_get_field(pbuf, i,lchnk, .true.)
      if (associated(qqcw)) qqcw_all(:,:,i) = qqcw(:,:)
   end do

   ifld = pbuf_get_index('DGNUMWET')
   call pbuf_get_field(pbuf, ifld, dgnumwet, start=(/1,1,1/), kount=(/pcols,pver,ntot_amode/))

   if(is_first_step())then
      do i = 1,pcnst
         qqcw_all(:,:,i) = 1.e-38_r8
      end do
      dgnumwet(1:pcols,1:pver,1:ntot_amode) = 0.0_r8
   end if

   qqcw_save = qqcw_all
   dgnumwet_save = dgnumwet
#endif

end subroutine crm_save_state_tend

!==================================================================================================
!==================================================================================================

subroutine crm_recall_state_tend(state,tend,pbuf)
   !------------------------------------------------------------------------------------------------
   ! This subroutine is used to recall the state that was saved prior to running the conventional 
   ! GCM physics routines.
   !------------------------------------------------------------------------------------------------
   use physics_types,   only: physics_state, physics_tend, physics_tend_dealloc,&
                              physics_state_copy, physics_tend_copy, physics_state_dealloc
   use time_manager,    only: is_first_step
   use physics_buffer,  only: pbuf_get_field, physics_buffer_desc, dyn_time_lvls
   use phys_control,    only: phys_getopts
#ifdef MODAL_AERO
   use modal_aero_data, only: qqcw_get_field
#endif

   implicit none

   type(physics_state),       intent(inout) :: state
   type(physics_tend),        intent(inout) :: tend 
   type(physics_buffer_desc), pointer       :: pbuf(:)

   integer ncol, lchnk, i, m
   logical :: use_ECPP
   character(len=16) :: microp_scheme           ! host model microphysics scheme
   real(r8), pointer, dimension(:,:) :: cld     ! cloud fraction
   real(r8), pointer, dimension(:,:) :: qqcw    ! 
   
   call phys_getopts( use_ECPP_out      = use_ECPP  )
   call phys_getopts( microp_scheme_out = microp_scheme )

   lchnk = state%lchnk
   ncol  = state%ncol

   ! gas and aerosol species are updated by shallow convective transport.
   ! Aerosol changes from dropmixnuc in cldwat2m.F90 are discarded. 
   ! (i.e. when dropmixnuc is called in cldwat2m.F90, the tendency is set to zero)
   q_aero = state%q

   ! Restore state and tend (from beginning of tphysbc)
   if ( allocated(state%t) ) call physics_state_dealloc(state)
   if ( allocated(tend%dtdt) ) call physics_tend_dealloc(tend)

   call physics_state_copy(state_save,state)
   call physics_tend_copy(tend_save,tend)

   call physics_state_dealloc(state_save)
   call physics_tend_dealloc(tend_save)

   ! tracer species other than water vapor and cloud water are updated in convetional CAM.
   ! When ECPP is used, dropmixnuc and all transport(deep and shallow) are done in ECPP.
   ! So all change in aerosol and gas species in the conventinal CAM are discarded.
   ! Minghuai Wang, 2010-01 (Minghuai.Wang@pnl.gov)
   if (.not. use_ECPP) then
      if ( microp_scheme .eq. 'MG' ) then
         state%q(:ncol,:pver,6:pcnst) = q_aero(:ncol,:pver,6:pcnst)
      else if ( microp_scheme .eq. 'RK' ) then
         state%q(:ncol,:pver,4:pcnst) = q_aero(:ncol,:pver,4:pcnst)
      end if
   end if

   ! whannah - not sure why we do this...
   if(is_first_step())then
      do m = 1,dyn_time_lvls
         call pbuf_get_field(pbuf, cldo_idx, cldo, start=(/1,1,m/), kount=(/pcols,pver,1/) )
         cldo(:ncol,:) = 0
      end do
   end if

   ! Restore cloud fraction
   cldo(:ncol, :) = cldo_save(:ncol, :)

#if ( defined MODAL_AERO )
   do i = 1,pcnst
      qqcw => qqcw_get_field(pbuf,i,lchnk,.true.)
      if (associated(qqcw)) qqcw(:,:) = qqcw_save(:,:,i)
   end do
   dgnumwet = dgnumwet_save
#endif

end subroutine crm_recall_state_tend

!==================================================================================================
!==================================================================================================

subroutine m2005_effradius(ql, nl,qi,ni,qs, ns, cld, pres, tk, &
                           effl, effi, effl_fn, deffi,         &
                           lamcrad, pgamrad, des)
   !------------------------------------------------------------------------------------------------
   ! This subroutine is used to calculate droplet and ice crystal effective radius, which will be
   ! used in the CAM radiation code. The method to calculate effective radius is taken out of the
   ! Morrision two moment scheme from M2005MICRO_GRAUPEL. It is also very similar to the subroutine
   ! effradius in the module of cldwat2m in the CAM source codes. 
   ! Adopted by Minghuai Wang (Minghuai.Wang@pnl.gov). 
   !------------------------------------------------------------------------------------------------
   ! Calculate effective radius for radiation code
   ! If no cloud water, default value is:
   !   10 micron for droplets,
   !   25 micron for cloud ice.
   ! Be careful of the unit of effective radius : [micro meter]
   !------------------------------------------------------------------------------------------------
   use shr_spfn_mod,    only: gamma => shr_spfn_gamma
   implicit none

   ! input arguments
   real(r8), intent(in)    :: ql          ! Mean LWC of pixels [ kg/kg ]
   real(r8), intent(in)    :: nl          ! Grid-mean number concentration of cloud liquid droplet [#/kg]
   real(r8), intent(in)    :: qi          ! Mean IWC of pixels [ kg/kg ]
   real(r8), intent(in)    :: ni          ! Grid-mean number concentration of cloud ice    droplet [#/kg]
   real(r8), intent(in)    :: qs          ! mean snow water content [kg/kg]
   real(r8), intent(in)    :: ns          ! Mean snow crystal number concnetration [#/kg]
   real(r8), intent(in)    :: cld         ! Physical stratus fraction
   real(r8), intent(in)    :: pres        ! Air pressure [Pa] 
   real(r8), intent(in)    :: tk          ! air temperature [K]

   ! output arguments
   real(r8), intent(out)   :: effl        ! Effective radius of cloud liquid droplet [micro-meter]
   real(r8), intent(out)   :: effi        ! Effective radius of cloud ice    droplet [micro-meter]
   real(r8), intent(out)   :: effl_fn     ! effl for fixed number concentration of nlic = 1.e8
   real(r8), intent(out)   :: deffi       ! ice effective diameter for optics (radiation)
   real(r8), intent(out)   :: pgamrad     ! gamma parameter for optics (radiation)
   real(r8), intent(out)   :: lamcrad     ! slope of droplet distribution for optics (radiation)
   real(r8), intent(out)   :: des         ! snow effective diameter for optics (radiation) [micro-meter]

   ! local variables
   real(r8)  qlic        ! In-cloud LWC [kg/m3]
   real(r8)  qiic        ! In-cloud IWC [kg/m3]
   real(r8)  nlic        ! In-cloud liquid number concentration [#/kg]
   real(r8)  niic        ! In-cloud ice    number concentration [#/kg]
   real(r8)  mtime       ! Factor to account for droplet activation timescale [no]
   real(r8)  cldm        ! Constrained stratus fraction [no]
   real(r8)  mincld      ! Minimum stratus fraction [no]

   real(r8)  lami, laml, lammax, lammin, pgam, lams, lammaxs, lammins

   real(r8)  dcs         ! autoconversion size threshold   [meter]
   real(r8)  di, ci      ! cloud ice mass-diameter relationship
   real(r8)  ds, cs      ! snow crystal mass-diameter relationship 
   real(r8)  qsmall      !
   real(r8)  rho         ! air density [kg/m3]
   real(r8)  rhow        ! liquid water density [kg/m3]
   real(r8)  rhoi        ! ice density [kg/m3]
   real(r8)  rhos        ! snow density [kg/m3]
   real(r8)  res         ! effective snow diameters
   real(r8)  pi          !
   real(r8)  tempnc      !

   !------------------------------------------------------------------------------------------------
   ! Main computation 
   !------------------------------------------------------------------------------------------------

   pi = 3.1415926535897932384626434
   ! qsmall = 1.0e-18  ! in the CAM source code (cldwat2m)
   qsmall = 1.0e-14  ! in the SAM source code (module_mp_graupel)
   ! rhow = 1000.      ! in cldwat2m, CAM 
   rhow = 997.       ! in module_mp_graupel, SAM
   rhoi = 500.       ! in both CAM and SAM

   ! dcs = 70.e-6_r8    ! in cldwat2m, CAM 
   dcs = 125.e-6_r8   ! in module_mp_graupel, SAM 
   ci = rhoi * pi/6.
   di = 3.

   ! for snow water
   rhos = 100.      ! in both SAM and CAM5 
   cs = rhos*pi/6.
   ds = 3.


   rho = pres / (287.15*tk)    ! air density [kg/m3]

   mincld  = 0.0001_r8
   cldm    = max(cld,mincld)
   qlic    = min(5.e-3_r8,max(0._r8,ql/cldm))
   qiic    = min(5.e-3_r8,max(0._r8,qi/cldm))
   nlic    = max(nl,0._r8)/cldm
   niic    = max(ni,0._r8)/cldm

   !------------------------------------------------------------------------------------------------
   ! Effective diameters of snow crystals
   !------------------------------------------------------------------------------------------------
   if(qs.gt.1.0e-7) then 
      lammaxs=1._r8/10.e-6_r8
      lammins=1._r8/2000.e-6_r8
      lams = (gamma(1._r8+ds)*cs * ns/qs)**(1._r8/ds)
      lams = min(lammaxs,max(lams,lammins))
      res = 1.5/lams*1.0e6_r8
   else
      res = 500._r8 
   end if 

   !
   ! from Hugh Morrision: rhos/917 accouts for assumptions about 
   ! ice density in the Mitchell optics. 
   !

   des = res * rhos/917._r8 *2._r8

   !------------------------------------------------------------------------------------------------
   ! Effective radius of cloud ice droplet 
   !------------------------------------------------------------------------------------------------

   if( qiic.ge.qsmall ) then
      niic   = min(niic,qiic*1.e20_r8)
      ! lammax = 1._r8/10.e-6_r8      ! in cldwat2m, CAM
      ! lammin = 1._r8/(2._r8*dcs)    ! in cldwat2m, CAM
      lammax = 1._r8/1.e-6_r8      ! in module_mp_graupel, SAM 
      lammin = 1._r8/(2._r8*dcs+100.e-6_r8)    ! in module_mp_graupel, SAM 
      lami   = (gamma(1._r8+di)*ci*niic/qiic)**(1._r8/di)
      lami   = min(lammax,max(lami,lammin))
      effi   = 1.5_r8/lami*1.e6_r8
   else
      effi   = 25._r8
   end if

   !--hm ice effective radius for david mitchell's optics
   !--ac morrison indicates that this is effective diameter
   !--ac morrison indicates 917 (for the density of pure ice..)
   deffi  = effi *rhoi/917._r8*2._r8

   !------------------------------------------------------------------------------------------------
   ! Effective radius of cloud liquid droplet 
   !------------------------------------------------------------------------------------------------

   if( qlic.ge.qsmall ) then
      ! Matin et al., 1994 (JAS) formula for pgam (the same is used in both CAM and SAM).
      ! See also Morrison and Grabowski (2007, JAS, Eq. (2))
      nlic   = min(nlic,qlic*1.e20_r8)

      ! set the minimum droplet number as 20/cm3.
      ! nlic   = max(nlic,20.e6_r8/rho) ! sghan minimum in #/cm3
      tempnc = nlic/rho/1.0e6    ! #/kg --> #/cm3
      ! if (tempnc.gt.100._r8) then 
      !   write(0, *) 'nc larger than 100  ', tempnc, rho
      ! end if

      !!!!!! ????? Should be the in-cloud dropelt number calculated as nlic*rho/1.0e6_r8 ????!!!! +++mhwang
      ! pgam   = 0.0005714_r8*(nlic/1.e6_r8/rho) + 0.2714_r8  !wrong, confirmed with Hugh Morrison. fixed in the latest SAM. 
      pgam   = 0.0005714_r8*(nlic*rho/1.e6_r8) + 0.2714_r8
      pgam   = 1._r8/(pgam**2)-1._r8
      ! pgam   = min(15._r8,max(pgam,2._r8))   ! in cldwat2m, CAM
      pgam   = min(10._r8,max(pgam,2._r8))   ! in module_mp_graupel, SAM
      ! if(pgam.gt.2.01_r8 .and.pgam.lt.9.99_r8) then
      !   write(0, *) 'pgam', pgam
      ! end if
      laml   = (pi/6._r8*rhow*nlic*gamma(pgam+4._r8)/(qlic*gamma(pgam+1._r8)))**(1._r8/3._r8)
      lammin = (pgam+1._r8)/50.e-6_r8    ! in cldwat2m, CAM
      lammax = (pgam+1._r8)/2.e-6_r8     ! in cldwat2m, CAM   ! cldwat2m should be used, 
                                                              ! if lammax is too large, this will lead to crash in 
                                                              ! src/physics/rrtmg/cloud_rad_props.F90 because 
                                                              ! klambda-1 can be zero in gam_liquid_lw and gam_liquid_sw
                                                              !  and g_lambda(kmu,klambda-1) will not be defined. 
      ! lammin = (pgam+1._r8)/60.e-6_r8    ! in module_mp_graupel, SAM
      ! lammax = (pgam+1._r8)/1.e-6_r8     ! in module_mp_graupel, SAM

      laml   = min(max(laml,lammin),lammax)
      ! effl   = gamma(qcvar+1._r8/3._r8)/(gamma(qcvar)*qcvar**(1._r8/3._r8))* &
      !          gamma(pgam+4._r8)/gamma(pgam+3._r8)/laml/2._r8*1.e6_r8      ! in cldwat2m, CAM
      effl   =  gamma(pgam+4._r8)/gamma(pgam+3._r8)/laml/2._r8*1.e6_r8  ! in module_mp_graupel, SAM
      lamcrad  = laml 
      pgamrad  = pgam
   else
      ! we chose 10. over 25, since 10 is a more reasonable value for liquid droplet. +++mhwang
      effl   = 10._r8     ! in cldwat2m, CAM
      ! effl   = 25._r8     ! in module_mp_graupel, SAM
      lamcrad  = 0.0_r8
      pgamrad  = 0.0_r8
   end if

   !------------------------------------------------------------------------------------------------
   ! Recalculate effective radius for constant number, in order to separate first and second 
   ! indirect effects. Assume constant number of 10^8 kg-1 
   !------------------------------------------------------------------------------------------------

   nlic = 1.e8
   if( qlic.ge.qsmall ) then
      ! Matin et al., 1994 (JAS) formula for pgam (the same is used in both CAM and SAM). 
      ! See also Morrison and Grabowski (2007, JAS, Eq. (2))  
      nlic   = min(nlic,qlic*1.e20_r8)
      pgam   = 0.0005714_r8*(nlic/1.e6_r8/rho) + 0.2714_r8
      pgam   = 1._r8/(pgam**2)-1._r8
      ! pgam   = min(15._r8,max(pgam,2._r8))   ! in cldwat2m, CAM
      pgam   = min(10._r8,max(pgam,2._r8))   ! in module_mp_graupel, SAM
      laml   = (pi/6._r8*rhow*nlic*gamma(pgam+4._r8)/(qlic*gamma(pgam+1._r8)))**(1._r8/3._r8)
      ! lammin = (pgam+1._r8)/50.e-6_r8    ! in cldwat2m, CAM
      ! lammax = (pgam+1._r8)/2.e-6_r8     ! in cldwat2m, CAM
      lammin = (pgam+1._r8)/60.e-6_r8    ! in module_mp_graupel, SAM
      lammax = (pgam+1._r8)/1.e-6_r8     ! in module_mp_graupel, SAM

      laml   = min(max(laml,lammin),lammax)
      ! effl_fn   = gamma(qcvar+1._r8/3._r8)/(gamma(qcvar)*qcvar**(1._r8/3._r8))* &
      !          gamma(pgam+4._r8)/gamma(pgam+3._r8)/laml/2._r8*1.e6_r8      ! in cldwat2m, CAM
      effl_fn   =  gamma(pgam+4._r8)/gamma(pgam+3._r8)/laml/2._r8*1.e6_r8  ! in module_mp_graupel, SAM
   else
      ! we chose 10. over 25, since 10 is a more reasonable value for liquid droplet. +++mhwang
      effl_fn   = 10._r8     ! in cldwat2m, CAM
      ! effl_fn   = 25._r8     ! in module_mp_graupel, SAM
   end if
   !------------------------------------------------------------------------------------------------
   !------------------------------------------------------------------------------------------------
   return
end subroutine m2005_effradius

!==================================================================================================
!==================================================================================================

end module crm_physics
