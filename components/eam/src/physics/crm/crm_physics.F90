module crm_physics
!---------------------------------------------------------------------------------------------------
! Purpose: Provides the interface to the crm code for the MMF configuration
!---------------------------------------------------------------------------------------------------
   use shr_kind_mod,    only: r8 => shr_kind_r8
   use shr_sys_mod,     only: shr_sys_flush   
   use shr_const_mod,   only: SHR_CONST_RDAIR,SHR_CONST_RWV
   use spmd_utils,      only: masterproc
   use cam_abortutils,  only: endrun
   use cam_control_mod, only: nsrest  ! restart flag
   use cam_logfile,     only: iulog
   use physics_types,   only: physics_state, physics_tend
   use ppgrid,          only: begchunk, endchunk, pcols, pver, pverp
   use constituents,    only: pcnst
#if defined(MODAL_AERO)
   use modal_aero_data, only: ntot_amode
#endif

   implicit none 
   private
   save

   public :: crm_physics_register
   public :: crm_physics_init
   public :: crm_physics_final
   public :: crm_physics_tend
   public :: crm_surface_flux_bypass_tend

   integer, public :: ncrms = -1 ! total number of CRMs summed over all chunks in task

   ! Constituent names - assigned according to MMF_microphysics_scheme
   character(len=8) :: cnst_names(8)

   integer :: ixcldliq  = -1   ! cloud liquid amount index
   integer :: ixcldice  = -1   ! cloud ice amount index
   integer :: ixnumliq  = -1   ! cloud liquid number index
   integer :: ixnumice  = -1   ! cloud ice water index
   integer :: ixrain    = -1   ! rain index
   integer :: ixsnow    = -1   ! snow index
   integer :: ixnumrain = -1   ! rain number index
   integer :: ixnumsnow = -1   ! snow number index
   integer :: ixcldrim  = -1   ! ice rime mass mixing ratio index
   integer :: ixrimvol  = -1   ! ice rime volume mixing ratio index
   integer :: idx_vt_t  = -1   ! CRM variance transport - liquid static energy
   integer :: idx_vt_q  = -1   ! CRM variance transport - total water
   integer :: idx_vt_u  = -1   ! CRM variance transport - horizontal momentum

   ! Physics buffer indices  
   integer :: ttend_dp_idx     = -1
   integer :: mmf_clear_rh_idx = -1
   integer :: crm_angle_idx    = -1
   integer :: cld_idx          = -1
   integer :: prec_dp_idx      = -1
   integer :: snow_dp_idx      = -1

   integer :: crm_t_rad_idx    = -1
   integer :: crm_qv_rad_idx   = -1
   integer :: crm_qc_rad_idx   = -1
   integer :: crm_qi_rad_idx   = -1
   integer :: crm_cld_rad_idx  = -1
   integer :: crm_qrad_idx     = -1

   integer :: crm_nc_rad_idx   = -1
   integer :: crm_ni_rad_idx   = -1
   integer :: crm_qs_rad_idx   = -1
   integer :: crm_ns_rad_idx   = -1

   integer :: crm_t_prev_idx   = -1
   integer :: crm_q_prev_idx   = -1

   integer :: crm_shoc_tk_idx       = -1
   integer :: crm_shoc_tkh_idx      = -1
   integer :: crm_shoc_wthv_idx     = -1
   integer :: crm_shoc_relvar_idx   = -1
   integer :: crm_shoc_cldfrac_idx  = -1

contains
!===================================================================================================
!===================================================================================================

subroutine crm_physics_register()
!---------------------------------------------------------------------------------------------------
! Purpose:  add necessary fields into physics buffer
!---------------------------------------------------------------------------------------------------
   use physconst,           only: cpair, mwh2o, mwdry
   use physics_buffer,      only: dyn_time_lvls, pbuf_add_field, dtype_r8, pbuf_get_index
   use phys_control,        only: phys_getopts, use_gw_convect
   use constituents,        only: cnst_add
   use crmdims,             only: crm_nx, crm_ny, crm_nz, &
                                  crm_dx, crm_dy, crm_dt, &
                                  crm_nx_rad, crm_ny_rad
   use constituents,        only: cnst_add
   use crm_history,         only: crm_history_register
#if defined(MMF_SAMXX) || defined(MMF_PAM)
   use gator_mod,           only: gator_init
#endif
#if defined(MMF_SAMXX)
   use cpp_interface_mod,   only: setparm
#elif defined(MMF_SAM) || defined(MMF_SAMOMP)
   use setparm_mod      ,   only: setparm
#endif
#if defined(MODAL_AERO)
   use modal_aero_data, only: ntot_amode
#endif
   !----------------------------------------------------------------------------
   ! local variables
   integer :: idx, c
   logical           :: use_ECPP
   logical           :: use_MMF_VT
   character(len=16) :: MMF_microphysics_scheme
   logical           :: prog_modal_aero
   integer, dimension(1) :: dims_gcm_1D
   integer, dimension(2) :: dims_gcm_2D
   integer, dimension(3) :: dims_gcm_3D
   integer, dimension(4) :: dims_crm_3D
   integer, dimension(4) :: dims_crm_rad
   integer :: cnst_ind ! dummy for adding new constituents for variance transport
   !----------------------------------------------------------------------------
#if defined(MMF_SAMXX) || defined(MMF_PAM)
   call gator_init()
#endif

   dims_gcm_1D  = (/pcols/)
   dims_gcm_2D  = (/pcols, pver/)
   dims_gcm_3D  = (/pcols, pver, dyn_time_lvls/)
   dims_crm_3D  = (/pcols, crm_nx, crm_ny, crm_nz/)
   dims_crm_rad = (/pcols, crm_nx_rad, crm_ny_rad, crm_nz/)

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
#if defined(MMF_SAMXX) || defined(MMF_SAM) || defined(MMF_SAMOMP)
   call setparm()
#endif

   if (use_MMF_VT) then
      ! add variance tracers
      call cnst_add('VT_T', real(0,r8), real(0,r8), real(0,r8), idx_vt_t, &
                    longname='CRM variance transport tracer for T', &
                    readiv=.false., mixtype='dry',cam_outfld=.false.)
      call cnst_add('VT_Q', real(0,r8), real(0,r8), real(0,r8), idx_vt_q, &
                    longname='CRM variance transport tracer for Q', &
                    readiv=.false., mixtype='dry',cam_outfld=.false.)
      call cnst_add('VT_U', real(0,r8), real(0,r8), real(0,r8), idx_vt_u, &
                    longname='CRM variance transport tracer for U', &
                    readiv=.false., mixtype='dry',cam_outfld=.false.)
   end if

   !----------------------------------------------------------------------------
   ! constituents
   !----------------------------------------------------------------------------
   if ( MMF_microphysics_scheme .eq. 'p3' ) then
      cnst_names(:) = (/'CLDLIQ','CLDICE','NUMLIQ','NUMICE','RAINQM','NUMRAI', 'CLDRIM','BVRIM '/)
      call cnst_add(cnst_names(1),  mwdry, cpair, 0._r8, ixcldliq, longname='Grid box averaged cld liquid amount',is_convtran1=.true.)
      call cnst_add(cnst_names(2),  mwdry, cpair, 0._r8, ixcldice, longname='Grid box averaged cld ice amount',   is_convtran1=.true.)
      call cnst_add(cnst_names(3),  mwh2o, cpair, 0._r8, ixnumliq, longname='Grid box averaged cld liquid number',is_convtran1=.true.)
      call cnst_add(cnst_names(4),  mwh2o, cpair, 0._r8, ixnumice, longname='Grid box averaged cld ice number',   is_convtran1=.true.)
      call cnst_add(cnst_names(5),  mwh2o, cpair, 0._r8, ixrain,   longname='Grid box averaged rain amount',      is_convtran1=.true.)
      call cnst_add(cnst_names(6),  mwh2o, cpair, 0._r8, ixnumrain,longname='Grid box averaged rain number',      is_convtran1=.true.)
      call cnst_add(cnst_names(7),  mwh2o, cpair, 0._r8, ixcldrim, longname='Grid box averaged rime amount',      is_convtran1=.true.)
      call cnst_add(cnst_names(8),  mwh2o, cpair, 0._r8, ixrimvol, longname='Grid box averaged rime volume',      is_convtran1=.true.)
   else
      cnst_names(:) = (/'CLDLIQ','CLDICE','NUMLIQ','NUMICE','RAINQM','SNOWQM','NUMRAI','NUMSNO'/)
      call cnst_add(cnst_names(1),  mwdry, cpair, 0._r8, ixcldliq, longname='Grid box averaged cld liquid amount',is_convtran1=.true.)
      call cnst_add(cnst_names(2),  mwdry, cpair, 0._r8, ixcldice, longname='Grid box averaged cld ice amount',   is_convtran1=.true.)
      call cnst_add(cnst_names(3),  mwh2o, cpair, 0._r8, ixnumliq, longname='Grid box averaged cld liquid number',is_convtran1=.true.)
      call cnst_add(cnst_names(4),  mwh2o, cpair, 0._r8, ixnumice, longname='Grid box averaged cld ice number',   is_convtran1=.true.)
      call cnst_add(cnst_names(5),  mwh2o, cpair, 0._r8, ixrain,   longname='Grid box averaged rain amount',      is_convtran1=.true.)
      call cnst_add(cnst_names(6),  mwh2o, cpair, 0._r8, ixsnow,   longname='Grid box averaged snow amount',      is_convtran1=.true.)
      call cnst_add(cnst_names(7),  mwh2o, cpair, 0._r8, ixnumrain,longname='Grid box averaged rain number',      is_convtran1=.true.)
      call cnst_add(cnst_names(8),  mwh2o, cpair, 0._r8, ixnumsnow,longname='Grid box averaged snow number',      is_convtran1=.true.)
   end if

   !----------------------------------------------------------------------------
   ! Register MMF history variables
   !----------------------------------------------------------------------------
   call crm_history_register()

   !----------------------------------------------------------------------------
   ! pbuf fields
   !----------------------------------------------------------------------------

   ! CRM state 
   call pbuf_add_field('CRM_U',        'global',dtype_r8,dims_crm_3D,idx)
   call pbuf_add_field('CRM_V',        'global',dtype_r8,dims_crm_3D,idx)
   call pbuf_add_field('CRM_W',        'global',dtype_r8,dims_crm_3D,idx)
   call pbuf_add_field('CRM_T',        'global',dtype_r8,dims_crm_3D,idx)
   call pbuf_add_field('CRM_RHO',      'global',dtype_r8,dims_crm_3D,idx)

   ! Radiation
   call pbuf_add_field('CRM_T_RAD',    'physpkg',dtype_r8,dims_crm_rad,crm_t_rad_idx)
   call pbuf_add_field('CRM_QV_RAD',   'physpkg',dtype_r8,dims_crm_rad,crm_qv_rad_idx)
   call pbuf_add_field('CRM_QC_RAD',   'physpkg',dtype_r8,dims_crm_rad,crm_qc_rad_idx)
   call pbuf_add_field('CRM_QI_RAD',   'physpkg',dtype_r8,dims_crm_rad,crm_qi_rad_idx)
   call pbuf_add_field('CRM_CLD_RAD',  'physpkg',dtype_r8,dims_crm_rad,crm_cld_rad_idx)
   call pbuf_add_field('CRM_QRAD',     'global', dtype_r8,dims_crm_rad,crm_qrad_idx)

   call pbuf_add_field('REI',          'physpkg',dtype_r8,dims_gcm_2D,idx)  ! Effective radius (ice)
   call pbuf_add_field('REL',          'physpkg',dtype_r8,dims_gcm_2D,idx)  ! Effective radius (liquid)
   call pbuf_add_field('DEI',          'physpkg',dtype_r8,dims_gcm_2D,idx)  ! Mitchell ice effective diameter for radiation
   call pbuf_add_field('MU',           'physpkg',dtype_r8,dims_gcm_2D,idx)  ! Size distribution shape parameter for radiation
   call pbuf_add_field('LAMBDAC',      'physpkg',dtype_r8,dims_gcm_2D,idx)  ! Size distribution shape parameter for radiation
   call pbuf_add_field('ICIWPST',      'physpkg',dtype_r8,dims_gcm_2D,idx)  ! Stratiform only in cloud ice water path for radiation
   call pbuf_add_field('ICLWPST',      'physpkg',dtype_r8,dims_gcm_2D,idx)  ! Stratiform in cloud liquid water path for radiation
   call pbuf_add_field('DES',          'physpkg',dtype_r8,dims_gcm_2D,idx)  ! Snow effective diameter for radiation
   call pbuf_add_field('ICSWP',        'physpkg',dtype_r8,dims_gcm_2D,idx)  ! In cloud snow water path for radiation
   call pbuf_add_field('CLDFSNOW',     'physpkg',dtype_r8,dims_gcm_3D,idx)  ! Cloud fraction for liquid drops + snow
   call pbuf_add_field('CLDO',         'global', dtype_r8,dims_gcm_3D,idx)  ! "old" cloud fraction
   call pbuf_add_field('CLD',          'global', dtype_r8,dims_gcm_3D,idx)  ! cloud fraction
   call pbuf_add_field('CONCLD',       'global', dtype_r8,dims_gcm_3D,idx)  ! convective cloud fraction

   call pbuf_add_field('CRM_QV',       'global', dtype_r8,dims_crm_3D,idx)

   if (MMF_microphysics_scheme .eq. 'sam1mom') then
      call pbuf_add_field('CRM_QP',    'global', dtype_r8,dims_crm_3D,idx)
      call pbuf_add_field('CRM_QN',    'global', dtype_r8,dims_crm_3D,idx)
   end if

   if (MMF_microphysics_scheme .eq. 'p3') then
      call pbuf_add_field('CRM_NC_RAD','physpkg',dtype_r8,dims_crm_rad,crm_nc_rad_idx)
      call pbuf_add_field('CRM_NI_RAD','physpkg',dtype_r8,dims_crm_rad,crm_ni_rad_idx)
      call pbuf_add_field('CRM_QC',    'global', dtype_r8,dims_crm_3D,idx)
      call pbuf_add_field('CRM_NC',    'global', dtype_r8,dims_crm_3D,idx)
      call pbuf_add_field('CRM_QR',    'global', dtype_r8,dims_crm_3D,idx)
      call pbuf_add_field('CRM_NR',    'global', dtype_r8,dims_crm_3D,idx)
      call pbuf_add_field('CRM_QI',    'global', dtype_r8,dims_crm_3D,idx)
      call pbuf_add_field('CRM_NI',    'global', dtype_r8,dims_crm_3D,idx)
      call pbuf_add_field('CRM_QM',    'global', dtype_r8,dims_crm_3D,idx)
      call pbuf_add_field('CRM_BM',    'global', dtype_r8,dims_crm_3D,idx)
      call pbuf_add_field('CRM_T_PREV','global', dtype_r8,dims_crm_3D,crm_t_prev_idx)
      call pbuf_add_field('CRM_Q_PREV','global', dtype_r8,dims_crm_3D,crm_q_prev_idx)
      if (prog_modal_aero) call pbuf_add_field('RATE1_CW2PR_ST','physpkg',dtype_r8,dims_gcm_2D,idx)
      call pbuf_add_field('CRM_SHOC_TK'     ,'global', dtype_r8,dims_crm_3D,crm_shoc_tk_idx)
      call pbuf_add_field('CRM_SHOC_THH'    ,'global', dtype_r8,dims_crm_3D,crm_shoc_tkh_idx)
      call pbuf_add_field('CRM_SHOC_WTHV'   ,'global', dtype_r8,dims_crm_3D,crm_shoc_wthv_idx)
      call pbuf_add_field('CRM_SHOC_RELVAR' ,'global', dtype_r8,dims_crm_3D,crm_shoc_relvar_idx)
      call pbuf_add_field('CRM_SHOC_CLDFRAC','global', dtype_r8,dims_crm_3D,crm_shoc_cldfrac_idx)
   end if

   ! CRM rad stuff specific to COSP; this does not strictly need to be in
   ! pbuf, we could compute it in rad and pass as optional arguments, but
   ! this seemed to be the cleanest solution for the time being (in other
   ! words, this is probably a lazy hack because I could not think of a
   ! cleaner way of passing these)
   call pbuf_add_field('CRM_REL' , 'physpkg', dtype_r8,dims_crm_3D,idx)
   call pbuf_add_field('CRM_REI' , 'physpkg', dtype_r8,dims_crm_3D,idx)

   ! CRM mass flux
   call pbuf_add_field('MU_CRM',       'physpkg',dtype_r8,dims_gcm_2D,idx) ! mass flux up
   call pbuf_add_field('MD_CRM',       'physpkg',dtype_r8,dims_gcm_2D,idx) ! mass flux down
   call pbuf_add_field('DU_CRM',       'physpkg',dtype_r8,dims_gcm_2D,idx) ! mass detrainment from updraft
   call pbuf_add_field('EU_CRM',       'physpkg',dtype_r8,dims_gcm_2D,idx) ! mass detrainment from updraft
   call pbuf_add_field('ED_CRM',       'physpkg',dtype_r8,dims_gcm_2D,idx) ! mass detrainment from downdraft

   call pbuf_add_field('JT_CRM',       'physpkg',dtype_r8,dims_gcm_1D,idx) ! index of cloud (convection) top for each column
   call pbuf_add_field('MX_CRM',       'physpkg',dtype_r8,dims_gcm_1D,idx) ! index of cloud (convection) bottom for each column
   call pbuf_add_field('IDEEP_CRM',    'physpkg',dtype_r8,dims_gcm_1D,idx) ! Gathering array for convective columns

   ! CRM turbulence
   call pbuf_add_field('TKE_CRM',      'physpkg',dtype_r8,dims_gcm_2D,idx) ! TKE from CRM  (m2/s2)
   call pbuf_add_field('TK_CRM',       'physpkg',dtype_r8,dims_gcm_2D,idx) ! TK from CRM (m2/s)

   ! ACLDY_CEN has to be global in the physcal buffer to be saved in the restart file
   ! total (all sub-classes) cloudy fractional area in previous time step 
   call pbuf_add_field('ACLDY_CEN',    'global', dtype_r8,dims_gcm_2D,idx) 
   
   ! CRM orientation angle needs to persist across time steps
   call pbuf_add_field('CRM_ANGLE',    'global', dtype_r8,dims_gcm_1D,crm_angle_idx)

   ! top and bottom levels of convective activity for chemistry
   call pbuf_add_field('CLDTOP',       'physpkg',dtype_r8,(/pcols,1/),idx)
   call pbuf_add_field('CLDBOT',       'physpkg',dtype_r8,(/pcols,1/),idx)

   ! Deep convective heating for convective gravity wave source
   if (use_gw_convect) call pbuf_add_field('TTEND_DP','physpkg',dtype_r8,dims_gcm_2D,ttend_dp_idx)

   !----------------------------------------------------------------------------
   ! miscellaneous fields previously added by offline parameterizations
   !----------------------------------------------------------------------------
   call pbuf_add_field('QME',          'physpkg',dtype_r8,dims_gcm_2D,idx) ! net condensation/evaporation of cloud water
   call pbuf_add_field('PRAIN',        'physpkg',dtype_r8,dims_gcm_2D,idx) ! total precip rate?
   call pbuf_add_field('NEVAPR',       'physpkg',dtype_r8,dims_gcm_2D,idx) ! total precip evaporation rate (rain + snow)
   call pbuf_add_field('PRER_EVAP',    'global' ,dtype_r8,dims_gcm_2D,idx) ! rain evaporation rate
   call pbuf_add_field('WSEDL',        'physpkg',dtype_r8,dims_gcm_2D,idx) ! Sed. velocity of liq stratus cloud droplet [m/s]
   call pbuf_add_field('ICWMRDP',      'physpkg',dtype_r8,dims_gcm_2D,idx) ! in-cloud deep conv water+ice mixing ratio
   call pbuf_add_field('ICWMRSH',      'physpkg',dtype_r8,dims_gcm_2D,idx) ! in-cloud shallow conv water+ice mixing ratio
   call pbuf_add_field('RPRDDP',       'physpkg',dtype_r8,dims_gcm_2D,idx) ! dq/dt due to deep convective rainout
   call pbuf_add_field('RPRDSH',       'physpkg',dtype_r8,dims_gcm_2D,idx) ! dq/dt due to shallow convective rainout
   call pbuf_add_field('RPRDTOT',      'physpkg',dtype_r8,dims_gcm_2D,idx) ! dq/dt due to total convective rainout
   call pbuf_add_field('NEVAPR_DPCU',  'physpkg',dtype_r8,dims_gcm_2D,idx) ! evaporation of deep convective precipitation
   call pbuf_add_field('NEVAPR_SHCU',  'physpkg',dtype_r8,dims_gcm_2D,idx) ! evaporation of shallow convective precipitation
   call pbuf_add_field('PREC_DP',      'physpkg',dtype_r8,dims_gcm_1D,idx) ! total precip from deep convection
   call pbuf_add_field('SNOW_DP',      'physpkg',dtype_r8,dims_gcm_1D,idx) ! snow from deep convection
   call pbuf_add_field('PREC_SH',      'physpkg',dtype_r8,dims_gcm_1D,idx) ! total precip from shallow convection
   call pbuf_add_field('SNOW_SH',      'physpkg',dtype_r8,dims_gcm_1D,idx) ! snow from shallow convection
   call pbuf_add_field('PREC_PCW',     'physpkg', dtype_r8,dims_gcm_1D,idx) ! total precip from "prognostic cloud scheme" (stratiform)
   call pbuf_add_field('SNOW_PCW',     'physpkg', dtype_r8,dims_gcm_1D,idx) ! snow from "prognostic cloud scheme" (stratiform)
   call pbuf_add_field('PREC_SED',     'physpkg', dtype_r8,dims_gcm_1D,idx) ! total precip from cloud sedimentation
   call pbuf_add_field('SNOW_SED',     'physpkg', dtype_r8,dims_gcm_1D,idx) ! snow from cloud sedimentation
   call pbuf_add_field('SH_FRAC',      'physpkg',dtype_r8,dims_gcm_2D,idx) ! shallow cloud fraction
   call pbuf_add_field('DP_FRAC',      'physpkg',dtype_r8,dims_gcm_2D,idx) ! deep cloud fraction
   call pbuf_add_field('FICE',         'physpkg',dtype_r8,dims_gcm_2D,idx) ! fraction of cloud water that is ice
   call pbuf_add_field('cush',         'global' ,dtype_r8,(/pcols,dyn_time_lvls/),idx ) !  Convective scale height
   call pbuf_add_field('AST',          'global' ,dtype_r8,(/pcols,pver,dyn_time_lvls/),idx) ! Stratiform cloud fraction

   !----------------------------------------------------------------------------
   !----------------------------------------------------------------------------

end subroutine crm_physics_register

!===================================================================================================
!===================================================================================================

subroutine crm_physics_init(state, pbuf2d, species_class)
!---------------------------------------------------------------------------------------------------
! Purpose: initialize some variables, and add necessary fields into output fields 
!---------------------------------------------------------------------------------------------------
   use time_manager,          only: is_first_step
   use physics_buffer,        only: physics_buffer_desc, pbuf_get_index, pbuf_set_field, pbuf_get_chunk
   use phys_control,          only: phys_getopts, use_gw_convect
   use phys_grid,             only: get_ncols_p
   use crm_history,           only: crm_history_init
   use cam_history,           only: addfld, hist_fld_active
   use constituents,          only: apcnst, bpcnst, cnst_name, cnst_longname, cnst_get_ind
   use constituents,          only: pcnst, cnst_get_ind
#ifdef ECPP
   use module_ecpp_ppdriver2, only: papampollu_init
#endif
   !----------------------------------------------------------------------------
   ! interface variables
   ! NOTE - species_class is an input so it needs to be outside of ifdef MODAL_AERO for 1-mom micro
   type(physics_state), intent(inout), dimension(begchunk:endchunk) :: state
   type(physics_buffer_desc), pointer :: pbuf2d(:,:)
   integer, intent(inout) :: species_class(:) 
   !----------------------------------------------------------------------------
   ! local variables
   integer :: m, mm, c
   integer :: ncnst
   logical :: use_ECPP
   logical :: use_MMF_VT
   character(len=16) :: MMF_microphysics_scheme
   integer :: ncol
   logical :: pam_stat_fields_active
   !----------------------------------------------------------------------------
   call phys_getopts(use_ECPP_out = use_ECPP)
   call phys_getopts(use_MMF_VT_out = use_MMF_VT)
   call phys_getopts(MMF_microphysics_scheme_out = MMF_microphysics_scheme)

   ! Determine total number of CRMs per task
   ncrms = 0
   do c=begchunk, endchunk
      ncrms = ncrms + state(c)%ncol
   end do
   
#ifdef ECPP
   ! Initialize ECPP driver
   if (use_ECPP) call papampollu_init()
#endif

   call crm_history_init(species_class)

   ! Register contituent history variables (previously added by micro_mg_cam.F90)
   ncnst = size(cnst_names)
   do m = 1, ncnst
      call cnst_get_ind(cnst_names(m), mm)
      if ( any(mm == (/ ixcldliq, ixcldice, ixrain, ixsnow, ixcldrim /)) ) then
         ! mass mixing ratios
         call addfld(cnst_name(mm), (/ 'lev' /), 'A', 'kg/kg', cnst_longname(mm)                   )
      else if ( any(mm == (/ ixnumliq, ixnumice, ixnumrain, ixnumsnow /)) ) then
         ! number concentrations
         call addfld(cnst_name(mm), (/ 'lev' /), 'A', '1/kg', cnst_longname(mm)                   )
      else if ( any(mm == (/ ixrimvol /)) ) then
         ! rime volume
         call addfld(cnst_name(mm), (/ 'lev' /), 'A', 'm3/kg', cnst_longname(mm)                   )
      else
         call endrun( "crm_physics_init: Could not call addfld for constituent with unknown units.")
      endif
   end do

   if ( MMF_microphysics_scheme .eq. 'sam1mom' ) then
      call addfld(apcnst(ixcldliq), (/'lev'/), 'A', 'kg/kg', trim(cnst_name(ixcldliq))//' after physics'  )
      call addfld(apcnst(ixcldice), (/'lev'/), 'A', 'kg/kg', trim(cnst_name(ixcldice))//' after physics'  )
      call addfld(bpcnst(ixcldliq), (/'lev'/), 'A', 'kg/kg', trim(cnst_name(ixcldliq))//' before physics' )
      call addfld(bpcnst(ixcldice), (/'lev'/), 'A', 'kg/kg', trim(cnst_name(ixcldice))//' before physics' )
      call addfld(apcnst(ixrain),   (/'lev'/), 'A', 'kg/kg', trim(cnst_name(ixrain))//' after physics'  )
      call addfld(apcnst(ixsnow),   (/'lev'/), 'A', 'kg/kg', trim(cnst_name(ixsnow))//' after physics'  )
      call addfld(bpcnst(ixrain),   (/'lev'/), 'A', 'kg/kg', trim(cnst_name(ixrain))//' before physics' )
      call addfld(bpcnst(ixsnow),   (/'lev'/), 'A', 'kg/kg', trim(cnst_name(ixsnow))//' before physics' )
   end if
   if ( MMF_microphysics_scheme .eq. 'p3' ) then
      call addfld(apcnst(ixcldliq ),(/'lev'/),'A','kg/kg',trim(cnst_name(ixcldliq ))//' after physics'  )
      call addfld(bpcnst(ixcldliq ),(/'lev'/),'A','kg/kg',trim(cnst_name(ixcldliq ))//' before physics' )
      call addfld(apcnst(ixcldice ),(/'lev'/),'A','kg/kg',trim(cnst_name(ixcldice ))//' after physics'  )
      call addfld(bpcnst(ixcldice ),(/'lev'/),'A','kg/kg',trim(cnst_name(ixcldice ))//' before physics' )
      call addfld(apcnst(ixnumliq ),(/'lev'/),'A',  '/kg',trim(cnst_name(ixnumliq ))//' after physics'  )
      call addfld(bpcnst(ixnumliq ),(/'lev'/),'A',  '/kg',trim(cnst_name(ixnumliq ))//' before physics' )
      call addfld(apcnst(ixnumice ),(/'lev'/),'A',  '/kg',trim(cnst_name(ixnumice ))//' after physics'  )
      call addfld(bpcnst(ixnumice ),(/'lev'/),'A',  '/kg',trim(cnst_name(ixnumice ))//' before physics' )
      call addfld(apcnst(ixrain   ),(/'lev'/),'A','kg/kg',trim(cnst_name(ixrain   ))//' after physics'  )
      call addfld(bpcnst(ixrain   ),(/'lev'/),'A','kg/kg',trim(cnst_name(ixrain   ))//' before physics' )
      call addfld(apcnst(ixnumrain),(/'lev'/),'A',  '/kg',trim(cnst_name(ixnumrain))//' after physics'  )
      call addfld(bpcnst(ixnumrain),(/'lev'/),'A',  '/kg',trim(cnst_name(ixnumrain))//' before physics' )
      call addfld(apcnst(ixcldrim ),(/'lev'/),'A','kg/kg',trim(cnst_name(ixcldrim ))//' after physics'  )
      call addfld(bpcnst(ixcldrim ),(/'lev'/),'A','kg/kg',trim(cnst_name(ixcldrim ))//' before physics' )
      call addfld(apcnst(ixrimvol ),(/'lev'/),'A','m3/kg',trim(cnst_name(ixrimvol ))//' after physics'  )
      call addfld(bpcnst(ixrimvol ),(/'lev'/),'A','m3/kg',trim(cnst_name(ixrimvol ))//' before physics' )
   end if

   if (use_MMF_VT) then
      ! initialize variance transport tracers
      do c=begchunk, endchunk
         ncol = state(c)%ncol
         state(c)%q(:ncol,:pver,idx_vt_t) = 0
         state(c)%q(:ncol,:pver,idx_vt_q) = 0
         state(c)%q(:ncol,:pver,idx_vt_u) = 0
      end do
   end if

   ! set pbuf indices
   mmf_clear_rh_idx = pbuf_get_index('MMF_CLEAR_RH')
   cld_idx          = pbuf_get_index('CLD')
   prec_dp_idx      = pbuf_get_index('PREC_DP')
   snow_dp_idx      = pbuf_get_index('SNOW_DP')

   ! Initialize pbuf variables
   if (is_first_step()) then
      call pbuf_set_field(pbuf2d, crm_t_rad_idx,  0._r8)
      call pbuf_set_field(pbuf2d, crm_qv_rad_idx, 0._r8)
      call pbuf_set_field(pbuf2d, crm_qc_rad_idx, 0._r8)
      call pbuf_set_field(pbuf2d, crm_qi_rad_idx, 0._r8)
      call pbuf_set_field(pbuf2d, crm_cld_rad_idx,0._r8)
      call pbuf_set_field(pbuf2d, crm_qrad_idx,   0._r8)

      call pbuf_set_field(pbuf2d, pbuf_get_index('CLDO')       , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('PRER_EVAP')  , 0._r8)

      call pbuf_set_field(pbuf2d, pbuf_get_index('ICWMRDP')    , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('RPRDDP')     , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('NEVAPR_DPCU'), 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('PREC_DP')    , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('SNOW_DP')    , 0._r8)

      call pbuf_set_field(pbuf2d, pbuf_get_index('ICWMRSH')    , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('RPRDSH')     , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('RPRDTOT')    , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('NEVAPR_SHCU'), 0._r8)
      call pbuf_set_field(pbuf2d, prec_dp_idx                  , 0._r8)
      call pbuf_set_field(pbuf2d, snow_dp_idx                  , 0._r8)

      call pbuf_set_field(pbuf2d, pbuf_get_index('ICIWPST')    , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('ICLWPST')    , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('ICSWP')      , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('CLDFSNOW')   , 0._r8)

      call pbuf_set_field(pbuf2d, pbuf_get_index('SH_FRAC')    , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('DP_FRAC')    , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('AST')        , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('FICE')       , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('REL')        , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('REI')        , 0._r8)

      call pbuf_set_field(pbuf2d, pbuf_get_index('PREC_PCW')   , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('SNOW_PCW')   , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('PREC_SED')   , 0._r8)
      call pbuf_set_field(pbuf2d, pbuf_get_index('SNOW_SED')   , 0._r8)

      if (use_gw_convect) call pbuf_set_field(pbuf2d, ttend_dp_idx, 0._r8)
   end if

end subroutine crm_physics_init

!===================================================================================================
!===================================================================================================

subroutine crm_physics_final()
#if defined(MMF_PAM)
   use gator_mod,       only: gator_finalize
   use pam_driver_mod,  only: pam_finalize
   call gator_finalize()
   call pam_finalize()
#endif
#if defined(MMF_SAMXX)
   use gator_mod, only: gator_finalize
   call gator_finalize()
#endif
end subroutine crm_physics_final

!===================================================================================================
!===================================================================================================

subroutine crm_physics_tend(ztodt, state, tend, ptend, pbuf2d, cam_in, cam_out, &
                            species_class, crm_ecpp_output, &
                            mmf_qchk_prec_dp, mmf_qchk_snow_dp, mmf_rad_flux)
   !------------------------------------------------------------------------------------------------
   ! Purpose: CRM interface for the MMF configuration to update GCM state
   ! Original Author: Marat Khairoutdinov
   !------------------------------------------------------------------------------------------------
   use perf_mod
   use physics_buffer,        only: physics_buffer_desc, pbuf_old_tim_idx, pbuf_get_index, &
                                    dyn_time_lvls, pbuf_get_field, pbuf_set_field, pbuf_get_chunk
   use physics_types,         only: physics_state, physics_tend, physics_ptend, physics_ptend_init
   use camsrfexch,            only: cam_in_t, cam_out_t
   use time_manager,          only: is_first_step, get_nstep
   use crmdims,               only: crm_nx, crm_ny, crm_nz, crm_nx_rad, crm_ny_rad, crm_dx, crm_dy, crm_dt
   use physconst,             only: cpair, latvap, latice, gravit, cappa, pi
   use constituents,          only: pcnst, cnst_get_ind
#if defined(MMF_SAMXX)
   use cpp_interface_mod,     only: crm
#elif defined(MMF_PAM)
   use pam_fortran_interface
   use pam_driver_mod,        only: pam_driver
#elif defined(MMF_SAM) || defined(MMF_SAMOMP)
   use crm_module,            only: crm
#endif
   use params_kind,           only: crm_rknd
   use phys_control,          only: phys_getopts, phys_do_flux_avg
   use crm_history,           only: crm_history_out
   use wv_saturation,         only: qsat_water
#if defined(MODAL_AERO)
   ! modal_aero_data only exists if MODAL_AERO
   use modal_aero_data,       only: ntot_amode, ntot_amode
   use ndrop,                 only: loadaer
#endif

   use RNG_MT ! random number generator for randomly rotating CRM orientation

   use crm_state_module,      only: crm_state_type, crm_state_initialize, crm_state_finalize
   use crm_rad_module,        only: crm_rad_type, crm_rad_initialize, crm_rad_finalize
   use crm_input_module,      only: crm_input_type, crm_input_initialize, crm_input_finalize
   use crm_output_module,     only: crm_output_type, crm_output_initialize, crm_output_finalize
   use crm_ecpp_output_module,only: crm_ecpp_output_type

   use iso_c_binding,         only: c_bool
   use phys_grid,             only: get_rlon_p, get_rlat_p, get_gcol_p  
   use openacc_utils,         only: prefetch

   real(r8),                                        intent(in   ) :: ztodt            ! global model time increment and CRM run length
   type(physics_state),dimension(begchunk:endchunk),intent(in   ) :: state            ! Global model state 
   type(physics_tend), dimension(begchunk:endchunk),intent(in   ) :: tend             ! 
   type(physics_ptend),dimension(begchunk:endchunk),intent(  out) :: ptend            ! output tendencies
   type(physics_buffer_desc), pointer                             :: pbuf2d(:,:)      ! physics buffer
   type(cam_in_t),     dimension(begchunk:endchunk),intent(in   ) :: cam_in           ! atm input from coupler
   type(cam_out_t),    dimension(begchunk:endchunk),intent(inout) :: cam_out          ! atm output to coupler
   integer,                                         intent(in   ) :: species_class(:) ! aerosol species type
   type(crm_ecpp_output_type),                      intent(inout) :: crm_ecpp_output  ! output data for ECPP calculations
   real(r8), dimension(begchunk:endchunk,pcols),    intent(  out) :: mmf_qchk_prec_dp ! precipitation diagostic (liq+ice)  used for check_energy_chng
   real(r8), dimension(begchunk:endchunk,pcols),    intent(  out) :: mmf_qchk_snow_dp ! precipitation diagostic (ice only) used for check_energy_chng
   real(r8), dimension(begchunk:endchunk,pcols),    intent(  out) :: mmf_rad_flux     ! radiative flux diagnostic used for check_energy_chng

   !------------------------------------------------------------------------------------------------
   ! Local variables 
   !------------------------------------------------------------------------------------------------
   integer lchnk                                   ! chunk identifier
   integer ncol                                    ! number of atmospheric columns
   integer nstep                                   ! time steps

   type(physics_buffer_desc), pointer :: pbuf_chunk(:) ! temporary pbuf pointer for single chunk

   ! convective precipitation variables on pbuf
   real(r8), pointer :: prec_dp(:)                 ! total precip from deep convection        [m/s]
   real(r8), pointer :: snow_dp(:)                 ! snow from deep convection                [m/s]
   real(r8), pointer :: ttend_dp(:,:)              ! Convective heating for gravity wave drag
   real(r8), pointer :: mmf_clear_rh(:,:)          ! clear air RH for aerosol water uptake
   real(r8), pointer :: cld(:,:)                   ! cloud fraction

#if defined(MODAL_AERO)
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

   real(r8), dimension(begchunk:endchunk,pcols) :: qli_hydro_before ! column-integraetd initial precipitating water+ice
   real(r8), dimension(begchunk:endchunk,pcols) ::  qi_hydro_before ! column-integrated initial precipitating ice
   real(r8), dimension(begchunk:endchunk,pcols) :: qli_hydro_after  ! column-integraetd final precipitating water+ice
   real(r8), dimension(begchunk:endchunk,pcols) ::  qi_hydro_after  ! column-integrated final precipitating ice

   real(r8) :: sfactor                             ! used to determine precip type for sam1mom

   integer  :: i, icrm, icol, k, m, ii, jj, c      ! loop iterators
   integer  :: ncol_sum                            ! ncol sum for chunk loops
   integer  :: icrm_beg, icrm_end                  ! CRM column index range for crm_history_out
   integer  :: itim                                ! pbuf field and "old time" indices
   real(r8) :: ideep_crm(pcols)                    ! gathering array for convective columns
   logical  :: lq(pcnst)                           ! flags for initializing ptend
   logical  :: use_ECPP                            ! flag for ECPP mode
   character(len=16) :: MMF_microphysics_scheme    ! CRM microphysics scheme

   logical(c_bool):: use_MMF_VT                    ! flag for MMF variance transport (for C++ CRM)
   logical(c_bool):: use_MMF_ESMT                  ! flag for MMF scalar momentum transport (for C++ CRM)
   logical        :: use_MMF_VT_tmp                ! flag for MMF variance transport (for Fortran CRM)
   logical        :: use_MMF_ESMT_tmp              ! flag for MMF scalar momentum transport (for Fortran CRM)
   integer        :: MMF_VT_wn_max                 ! wavenumber cutoff for filtered variance transport

   real(r8) :: tmp_e_sat                           ! temporary saturation vapor pressure
   real(r8) :: tmp_q_sat                           ! temporary saturation specific humidity
   real(r8) :: tmp_rh_sum                          ! temporary relative humidity sum
   real(r8) :: tmp_rh_cnt                          ! temporary relative humidity count

   ! variables for changing CRM orientation
   real(crm_rknd) :: MMF_orientation_angle         ! CRM orientation [deg] (convert to radians)
   real(crm_rknd) :: unif_rand1, unif_rand2        ! uniform random numbers 
   real(crm_rknd) :: norm_rand                     ! normally distributed random number - Box-Muller (1958)
   real(crm_rknd) :: crm_rotation_std              ! scaling factor for rotation (std dev of rotation angle)
   real(crm_rknd) :: crm_rotation_offset           ! offset to specify preferred rotation direction 
   real(crm_rknd), pointer :: crm_angle(:)         ! CRM orientation angle (pbuf)

   ! surface flux variables for using adjusted fluxes from flux_avg_run
   real(crm_rknd), pointer, dimension(:) :: shf_ptr
   real(crm_rknd), pointer, dimension(:) :: lhf_ptr
   real(crm_rknd), pointer, dimension(:) :: wsx_ptr
   real(crm_rknd), pointer, dimension(:) :: wsy_ptr

   real(crm_rknd), dimension(pcols) :: shf_tmp
   real(crm_rknd), dimension(pcols) :: lhf_tmp
   real(crm_rknd), dimension(pcols) :: wsx_tmp
   real(crm_rknd), dimension(pcols) :: wsy_tmp

   real(crm_rknd) :: tmp_dz
   real(crm_rknd) :: tmp_dp

   ! CRM types
   type(crm_state_type)  :: crm_state
   type(crm_rad_type)    :: crm_rad
   type(crm_input_type)  :: crm_input
   type(crm_output_type) :: crm_output

   real(crm_rknd), allocatable :: crm_clear_rh(:,:) ! clear air relative humidity for aerosol wateruptake

   real(crm_rknd), allocatable :: longitude0(:)
   real(crm_rknd), allocatable :: latitude0 (:)
   integer       , allocatable :: gcolp     (:)
   real(crm_rknd)              :: crm_accel_factor
   logical                     :: use_crm_accel_tmp
   logical                     :: crm_accel_uv_tmp
   logical(c_bool)             :: use_crm_accel
   logical(c_bool)             :: crm_accel_uv

   ! pointers for crm_rad data on pbuf
   real(crm_rknd), pointer :: crm_qrad   (:,:,:,:) ! rad heating
   real(crm_rknd), pointer :: crm_t_rad  (:,:,:,:) ! rad temperature
   real(crm_rknd), pointer :: crm_qv_rad (:,:,:,:) ! rad vapor
   real(crm_rknd), pointer :: crm_qc_rad (:,:,:,:) ! rad cloud water
   real(crm_rknd), pointer :: crm_qi_rad (:,:,:,:) ! rad cloud ice
   real(crm_rknd), pointer :: crm_cld_rad(:,:,:,:) ! rad cloud fraction
   real(crm_rknd), pointer :: crm_nc_rad (:,:,:,:) ! rad cloud droplet number (#/kg)
   real(crm_rknd), pointer :: crm_ni_rad (:,:,:,:) ! rad cloud ice crystal number (#/kg)
   real(crm_rknd), pointer :: crm_qs_rad (:,:,:,:) ! rad cloud snow (kg/kg)
   real(crm_rknd), pointer :: crm_ns_rad (:,:,:,:) ! rad cloud snow crystal number (#/kg)

   ! pointers for crm state data on pbuf
   real(crm_rknd), pointer :: crm_u (:,:,:,:) ! CRM u-wind component
   real(crm_rknd), pointer :: crm_v (:,:,:,:) ! CRM v-wind component
   real(crm_rknd), pointer :: crm_w (:,:,:,:) ! CRM w-wind component
   real(crm_rknd), pointer :: crm_t (:,:,:,:) ! CRM temperature
   real(crm_rknd), pointer :: crm_rho_d(:,:,:,:) ! CRM dry density
   real(crm_rknd), pointer :: crm_qv(:,:,:,:) ! CRM water vapor mixing ratio (kg/kg of wet air)

   real(crm_rknd), pointer :: crm_qp(:,:,:,:) ! 1-mom mass mixing ratio of precipitating condensate
   real(crm_rknd), pointer :: crm_qn(:,:,:,:) ! 1-mom mass mixing ratio of cloud condensate

   real(crm_rknd), pointer :: crm_qc(:,:,:,:)      ! p3 mass mixing ratio of cloud liquid
   real(crm_rknd), pointer :: crm_qr(:,:,:,:)      ! p3 mass mixing ratio of rain
   real(crm_rknd), pointer :: crm_qi(:,:,:,:)      ! p3 mass mixing ratio of cloud ice
   real(crm_rknd), pointer :: crm_nc(:,:,:,:)      ! p3 number concentration of cloud water
   real(crm_rknd), pointer :: crm_nr(:,:,:,:)      ! p3 number concentration of rain
   real(crm_rknd), pointer :: crm_ni(:,:,:,:)      ! p3 number concentration of cloud ice
   real(crm_rknd), pointer :: crm_qm(:,:,:,:)      ! p3 rime density
   real(crm_rknd), pointer :: crm_bm(:,:,:,:)      ! p3 rime volume
   real(crm_rknd), pointer :: crm_q_prev(:,:,:,:)  ! p3 previous qv
   real(crm_rknd), pointer :: crm_t_prev(:,:,:,:)  ! p3 previous t 
   real(crm_rknd), pointer :: crm_shoc_tk(:,:,:,:)      ! SHOC Eddy coefficient for momentum [m2/s]
   real(crm_rknd), pointer :: crm_shoc_tkh(:,:,:,:)     ! SHOC Eddy coefficent for heat [m2/s]
   real(crm_rknd), pointer :: crm_shoc_wthv(:,:,:,:)    ! SHOC Buoyancy flux [K m/s]
   real(crm_rknd), pointer :: crm_shoc_relvar(:,:,:,:)  ! SHOC Relative cloud water variance
   real(crm_rknd), pointer :: crm_shoc_cldfrac(:,:,:,:) ! SHOC Cloud fraction [-]

   !------------------------------------------------------------------------------------------------
   !------------------------------------------------------------------------------------------------

   call phys_getopts(use_ECPP_out = use_ECPP)
   call phys_getopts(MMF_microphysics_scheme_out = MMF_microphysics_scheme)
   call phys_getopts(MMF_orientation_angle_out = MMF_orientation_angle)

   ! CRM variance transport
   use_MMF_VT = .false.
   call phys_getopts(use_MMF_VT_out = use_MMF_VT_tmp)
   call phys_getopts(MMF_VT_wn_max_out = MMF_VT_wn_max)
   use_MMF_VT = use_MMF_VT_tmp

   ! CRM explicit scalar momentum transport (ESMT)
   use_MMF_ESMT = .false.
   call phys_getopts(use_MMF_ESMT_out = use_MMF_ESMT_tmp)
   use_MMF_ESMT = use_MMF_ESMT_tmp

   ! CRM mean state acceleration (MSA) parameters
   use_crm_accel = .false.
   crm_accel_factor = 0.
   crm_accel_uv = .false.
   call phys_getopts(use_crm_accel_out    = use_crm_accel_tmp)
   call phys_getopts(crm_accel_factor_out = crm_accel_factor)
   call phys_getopts(crm_accel_uv_out     = crm_accel_uv_tmp)
   use_crm_accel = use_crm_accel_tmp
   crm_accel_uv = crm_accel_uv_tmp

   nstep = get_nstep()
   itim = pbuf_old_tim_idx() ! "Old" pbuf time index (what does all this mean?)

   call t_startf ('crm')

   !------------------------------------------------------------------------------------------------
   ! Initialize ptend
   !------------------------------------------------------------------------------------------------
   lq(:) = .true.
   do c=begchunk, endchunk
      call physics_ptend_init(ptend(c), state(c)%psetcols, 'crm', lu=.true., lv=.true., ls=.true., lq=lq)
   end do
   
   !------------------------------------------------------------------------------------------------
   ! Initialize CRM state (nullify pointers, allocate memory, etc)
   !------------------------------------------------------------------------------------------------
   call crm_state_initialize(crm_state, ncrms, crm_nx, crm_ny, crm_nz, MMF_microphysics_scheme)
   call crm_rad_initialize(crm_rad, ncrms, crm_nx_rad, crm_ny_rad, crm_nz, MMF_microphysics_scheme)
   call crm_input_initialize(crm_input, ncrms, pver, MMF_microphysics_scheme)
   call crm_output_initialize(crm_output, ncrms, pver, crm_nx, crm_ny, crm_nz, MMF_microphysics_scheme)

   !------------------------------------------------------------------------------------------------
   ! Set CRM orientation angle
   !------------------------------------------------------------------------------------------------
   do c=begchunk, endchunk
      lchnk = state(c)%lchnk
      ncol = state(c)%ncol
      pbuf_chunk => pbuf_get_chunk(pbuf2d, c)

      call pbuf_get_field(pbuf_chunk, crm_angle_idx, crm_angle)

      if (MMF_orientation_angle==-1) then

         crm_rotation_std    = 20. * pi/180.                 ! std deviation of normal distribution for CRM rotation [radians]
         crm_rotation_offset = 90. * pi/180. * ztodt/86400.  ! This means that a CRM should rotate 90 deg / day on average

         if (is_first_step()) crm_angle(1:ncol) = 0

         ! Rotate the CRM using a random walk
         if ( (crm_ny.eq.1) .or. (crm_nx.eq.1) ) then
            do i = 1,ncol
               ! set the seed based on the chunk and column index (duplicate seeds are ok)
               call RNG_MT_set_seed( lchnk + i + nstep )
               ! Generate a pair of uniform random numbers
               call RNG_MT_gen_rand(unif_rand1)
               call RNG_MT_gen_rand(unif_rand2)
               ! Box-Muller (1958) method of obtaining a Gaussian distributed random number
               norm_rand = sqrt(-2.*log(unif_rand1))*cos(pi*2*unif_rand2)
               crm_angle(i) = crm_angle(i) + norm_rand * crm_rotation_std + crm_rotation_offset
               ! Adjust CRM orientation angle to be between 0 and 2*pi
               if ( crm_angle(i).lt. 0.   ) crm_angle(i) = crm_angle(i) + pi*2
               if ( crm_angle(i).gt.(pi*2)) crm_angle(i) = crm_angle(i) - pi*2
            end do ! i
         end if

      else

         ! use static CRM orientation (no rotation) - only set pbuf values once
         if (is_first_step()) crm_angle(1:ncol) = MMF_orientation_angle * pi/180.

      end if

   end do ! c=begchunk, endchunk

   !------------------------------------------------------------------------------------------------
   ! pbuf initialization
   !------------------------------------------------------------------------------------------------
   do c=begchunk, endchunk
      pbuf_chunk => pbuf_get_chunk(pbuf2d, c)
      
      ! Zero these fields to ensure balanced water in land input
      call pbuf_set_field(pbuf_chunk, pbuf_get_index('PREC_SED'), 0._r8 )
      call pbuf_set_field(pbuf_chunk, pbuf_get_index('SNOW_SED'), 0._r8 )
      call pbuf_set_field(pbuf_chunk, pbuf_get_index('PREC_PCW'), 0._r8 )
      call pbuf_set_field(pbuf_chunk, pbuf_get_index('SNOW_PCW'), 0._r8 )
      call pbuf_set_field(pbuf_chunk, pbuf_get_index('PREC_SH'),  0._r8)
      call pbuf_set_field(pbuf_chunk, pbuf_get_index('SNOW_SH'),  0._r8)

      ! set convective rain to be zero for PRAIN already includes precipitation production from convection. 
      call pbuf_set_field(pbuf_chunk, pbuf_get_index('RPRDTOT'), 0.0_r8 )
      call pbuf_set_field(pbuf_chunk, pbuf_get_index('RPRDDP' ), 0.0_r8 )
      call pbuf_set_field(pbuf_chunk, pbuf_get_index('RPRDSH' ), 0.0_r8 )
      call pbuf_set_field(pbuf_chunk, pbuf_get_index('ICWMRDP'), 0.0_r8 )
      call pbuf_set_field(pbuf_chunk, pbuf_get_index('ICWMRSH'), 0.0_r8 )

   end do ! c=begchunk, endchunk
   !------------------------------------------------------------------------------------------------
   !------------------------------------------------------------------------------------------------

   if (is_first_step()) then

      do c=begchunk, endchunk
         pbuf_chunk => pbuf_get_chunk(pbuf2d, c)
         ncol = state(c)%ncol

         call pbuf_get_field(pbuf_chunk, crm_angle_idx, crm_angle)

         !------------------------------------------------------------------------------------------
         ! initialize CRM state stored in pbuf
         !------------------------------------------------------------------------------------------
         
         ! Set pointers to crm_state fields that persist on physics buffer
         call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_U'),  crm_u)
         call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_V'),  crm_v)
         call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_W'),  crm_w)
         call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_T'),  crm_t)
         call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_RHO'),crm_rho_d)
         call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QV'), crm_qv)
         if (MMF_microphysics_scheme .eq. 'sam1mom') then
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QP'), crm_qp)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QN'), crm_qn)
         else if (MMF_microphysics_scheme .eq. 'p3') then
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QC'), crm_qc)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_NC'), crm_nc)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QR'), crm_qr)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_NR'), crm_nr)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QI'), crm_qi)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_NI'), crm_ni)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QM'), crm_qm)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_BM'), crm_bm)
            call pbuf_get_field(pbuf_chunk, crm_t_prev_idx, crm_t_prev)
            call pbuf_get_field(pbuf_chunk, crm_q_prev_idx, crm_q_prev)
            call pbuf_get_field(pbuf_chunk, crm_shoc_tk_idx     , crm_shoc_tk)
            call pbuf_get_field(pbuf_chunk, crm_shoc_tkh_idx    , crm_shoc_tkh)
            call pbuf_get_field(pbuf_chunk, crm_shoc_wthv_idx   , crm_shoc_wthv)
            call pbuf_get_field(pbuf_chunk, crm_shoc_relvar_idx , crm_shoc_relvar)
            call pbuf_get_field(pbuf_chunk, crm_shoc_cldfrac_idx, crm_shoc_cldfrac)
         end if

         ! initialize all water to zero (needed for ncol < i <= pcols)
         crm_qv(1:pcols,:,:,:) = 0.0_r8
         if (MMF_microphysics_scheme .eq. 'sam1mom') then
            crm_qp(1:pcols,:,:,:) = 0.0_r8
            crm_qn(1:pcols,:,:,:) = 0.0_r8
         else if (MMF_microphysics_scheme .eq. 'p3') then
            crm_qc(1:pcols,:,:,:) = 0.0_r8
            crm_nc(1:pcols,:,:,:) = 0.0_r8
            crm_qr(1:pcols,:,:,:) = 0.0_r8
            crm_nr(1:pcols,:,:,:) = 0.0_r8
            crm_qi(1:pcols,:,:,:) = 0.0_r8
            crm_ni(1:pcols,:,:,:) = 0.0_r8
            crm_qm(1:pcols,:,:,:) = 0.0_r8
            crm_bm(1:pcols,:,:,:) = 0.0_r8
         end if

         do i = 1,ncol
            do k = 1,crm_nz
               m = pver-k+1

               crm_u(i,:,:,k) = state(c)%u(i,m) * cos( crm_angle(i) ) + state(c)%v(i,m) * sin( crm_angle(i) )
               crm_v(i,:,:,k) = state(c)%v(i,m) * cos( crm_angle(i) ) - state(c)%u(i,m) * sin( crm_angle(i) )
               crm_w(i,:,:,k) = 0.
               crm_t(i,:,:,k) = state(c)%t(i,m)

               ! diagnose dry density for PAM using:
               !   pdel = pdel_dry + pdel*qv
               !   pdel_dry = rho_dry*dz*g
               tmp_dp = state(c)%pint(i,m+1) - state(c)%pint(i,m)
               tmp_dz = state(c)%zi(i,m+1) - state(c)%zi(i,m)
               crm_rho_d(i,:,:,k) = -1 * tmp_dp * ( 1 - state(c)%q(i,m,1) ) / ( tmp_dz * gravit )

               ! Initialize microphysics arrays
               if (MMF_microphysics_scheme .eq. 'sam1mom') then
                  crm_qv(i,:,:,k) = state(c)%q(i,m,1)!+state(c)%q(i,m,ixcldliq)+state(c)%q(i,m,ixcldice)
                  crm_qp(i,:,:,k) = 0.0_r8
                  crm_qn(i,:,:,k) = state(c)%q(i,m,ixcldliq)+state(c)%q(i,m,ixcldice)
               else if (MMF_microphysics_scheme .eq. 'p3') then
                  crm_qv(i,:,:,k) = state(c)%q(i,m,1)!+state(c)%q(i,m,ixcldliq)
                  crm_qc(i,:,:,k) = state(c)%q(i,m,ixcldliq)
                  crm_qi(i,:,:,k) = state(c)%q(i,m,ixcldice)
                  crm_qr(i,:,:,k) = 0.0_r8
                  crm_nc(i,:,:,k) = 0.0_r8
                  crm_ni(i,:,:,k) = 0.0_r8
                  crm_nr(i,:,:,k) = 0.0_r8
                  crm_qm(i,:,:,k) = 0.0_r8
                  crm_bm(i,:,:,k) = 0.0_r8
                  crm_t_prev(i,:,:,k) = state(c)%t(i,m)
                  crm_q_prev(i,:,:,k) = state(c)%q(i,m,1)
                  crm_shoc_tk     (i,:,:,k) = 0.0_r8
                  crm_shoc_tkh    (i,:,:,k) = 0.0_r8
                  crm_shoc_wthv   (i,:,:,k) = 0.0_r8
                  crm_shoc_relvar (i,:,:,k) = 0.0_r8
                  crm_shoc_cldfrac(i,:,:,k) = 0.0_r8
               end if

            end do
         end do

         !------------------------------------------------------------------------------------------
         ! Initialize radiation variables
         !------------------------------------------------------------------------------------------
         call pbuf_get_field(pbuf_chunk, crm_t_rad_idx,  crm_t_rad)
         call pbuf_get_field(pbuf_chunk, crm_qv_rad_idx, crm_qv_rad)
         call pbuf_get_field(pbuf_chunk, crm_qc_rad_idx, crm_qc_rad)
         call pbuf_get_field(pbuf_chunk, crm_qi_rad_idx, crm_qi_rad)
         call pbuf_get_field(pbuf_chunk, crm_cld_rad_idx,crm_cld_rad)
         call pbuf_get_field(pbuf_chunk, crm_qrad_idx,   crm_qrad)
         if (MMF_microphysics_scheme .eq. 'p3') then
            call pbuf_get_field(pbuf_chunk, crm_nc_rad_idx,crm_nc_rad)
            call pbuf_get_field(pbuf_chunk, crm_ni_rad_idx,crm_ni_rad)
         end if

         do k = 1,crm_nz
            m = pver-k+1
            do i = 1,ncol
               crm_qrad   (i,:,:,k) = 0.
               crm_t_rad  (i,:,:,k) = state(c)%t(i,m)
               crm_qv_rad (i,:,:,k) = state(c)%q(i,m,1)
               crm_qc_rad (i,:,:,k) = 0.
               crm_qi_rad (i,:,:,k) = 0.
               crm_cld_rad(i,:,:,k) = 0.
               if (MMF_microphysics_scheme .eq. 'p3') then
                  crm_nc_rad(i,:,:,k) = 0.0
                  crm_ni_rad(i,:,:,k) = 0.0
               end if
            end do
         end do

         ! use radiation from grid-cell mean radctl on first time step
         call pbuf_get_field(pbuf_chunk, cld_idx, cld, start=(/1,1,itim/), kount=(/pcols,pver,1/) )
         cld(1:ncol,:) = 0.

         ! Set clear air RH to zero on first step
         call pbuf_get_field(pbuf_chunk, mmf_clear_rh_idx, mmf_clear_rh )
         mmf_clear_rh(1:ncol,1:pver) = 0

         !------------------------------------------------------------------------------------------
         ! initialize ECPP variables
         !------------------------------------------------------------------------------------------
#ifdef ECPP
         if (use_ECPP) then
            ! initialize turbulence for ECPP calculations
            call pbuf_set_field(pbuf_chunk, pbuf_get_index('TKE_CRM'), 0.0_r8 )
            call pbuf_set_field(pbuf_chunk, pbuf_get_index('TK_CRM'), 0.0_r8 )
         end if 
#endif
         !------------------------------------------------------------------------------------------
         ! zero output fluxes to avoid error in check_energy_chng()
         !------------------------------------------------------------------------------------------
         do i = 1,ncol
            mmf_qchk_prec_dp(c,i) = 0.0
            mmf_qchk_snow_dp(c,i) = 0.0
            mmf_rad_flux(c,i) = 0.0
         end do ! i = 1,ncol

         !------------------------------------------------------------------------------------------
         !------------------------------------------------------------------------------------------
      end do ! c=begchunk, endchunk

   else  ! not is_first_step

      ncol_sum = 0
      do c=begchunk, endchunk
         pbuf_chunk => pbuf_get_chunk(pbuf2d, c)
         ncol = state(c)%ncol

         call pbuf_get_field(pbuf_chunk, crm_angle_idx, crm_angle)

         !------------------------------------------------------------------------------------------
         ! Retreive CRM state data from pbuf
         !------------------------------------------------------------------------------------------

         ! Set pointers to crm_state fields that persist on physics buffer
         call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_U'),  crm_u)
         call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_V'),  crm_v)
         call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_W'),  crm_w)
         call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_T'),  crm_t)
         call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_RHO'),crm_rho_d)
         call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QV'), crm_qv)
         if (MMF_microphysics_scheme .eq. 'sam1mom') then
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QP'), crm_qp)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QN'), crm_qn)
         end if
         if (MMF_microphysics_scheme .eq. 'p3') then
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_NC'), crm_nc)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QR'), crm_qr)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_NR'), crm_nr)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QI'), crm_qi)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_NI'), crm_ni)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QC'), crm_qc)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QM'), crm_qm)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_BM'), crm_bm)
            call pbuf_get_field(pbuf_chunk, crm_t_prev_idx, crm_t_prev)
            call pbuf_get_field(pbuf_chunk, crm_q_prev_idx, crm_q_prev)
            call pbuf_get_field(pbuf_chunk, crm_shoc_tk_idx     , crm_shoc_tk)
            call pbuf_get_field(pbuf_chunk, crm_shoc_tkh_idx    , crm_shoc_tkh)
            call pbuf_get_field(pbuf_chunk, crm_shoc_wthv_idx   , crm_shoc_wthv)
            call pbuf_get_field(pbuf_chunk, crm_shoc_relvar_idx , crm_shoc_relvar)
            call pbuf_get_field(pbuf_chunk, crm_shoc_cldfrac_idx, crm_shoc_cldfrac)
         end if

         ! copy pbuf data into crm_state
         do i = 1,ncol
            icrm = ncol_sum + i
            crm_state%u_wind     (icrm,:,:,:) = crm_u (i,:,:,:)
            crm_state%v_wind     (icrm,:,:,:) = crm_v (i,:,:,:)
            crm_state%w_wind     (icrm,:,:,:) = crm_w (i,:,:,:)
            crm_state%temperature(icrm,:,:,:) = crm_t (i,:,:,:)
            crm_state%rho_dry    (icrm,:,:,:) = crm_rho_d(i,:,:,:)
            crm_state%qv         (icrm,:,:,:) = crm_qv(i,:,:,:)
            if (MMF_microphysics_scheme .eq. 'sam1mom') then
               crm_state%qp      (icrm,:,:,:) = crm_qp(i,:,:,:)
               crm_state%qn      (icrm,:,:,:) = crm_qn(i,:,:,:)
            end if
            if (MMF_microphysics_scheme .eq. 'p3') then
               crm_state%qc      (icrm,:,:,:) = crm_qc(i,:,:,:)
               crm_state%qi      (icrm,:,:,:) = crm_qi(i,:,:,:)
               crm_state%qr      (icrm,:,:,:) = crm_qr(i,:,:,:)
               crm_state%nc      (icrm,:,:,:) = crm_nc(i,:,:,:)
               crm_state%ni      (icrm,:,:,:) = crm_ni(i,:,:,:)
               crm_state%nr      (icrm,:,:,:) = crm_nr(i,:,:,:)
               crm_state%qm      (icrm,:,:,:) = crm_qm(i,:,:,:)
               crm_state%bm      (icrm,:,:,:) = crm_bm(i,:,:,:)
               crm_state%t_prev  (icrm,:,:,:) = crm_t_prev(i,:,:,:)
               crm_state%q_prev  (icrm,:,:,:) = crm_q_prev(i,:,:,:)
               crm_state%shoc_tk     (icrm,:,:,:) = crm_shoc_tk     (i,:,:,:)
               crm_state%shoc_tkh    (icrm,:,:,:) = crm_shoc_tkh    (i,:,:,:)
               crm_state%shoc_wthv   (icrm,:,:,:) = crm_shoc_wthv   (i,:,:,:)
               crm_state%shoc_relvar (icrm,:,:,:) = crm_shoc_relvar (i,:,:,:)
               crm_state%shoc_cldfrac(icrm,:,:,:) = crm_shoc_cldfrac(i,:,:,:)
            end if
         end do ! i=1,ncol

         ! Retrieve radiative heating tendency
         ! Rad tendency was normalized for energy conservation 
         ! un-normalize before sending to CRM - this yields units of [K/s]
         call pbuf_get_field(pbuf_chunk, crm_qrad_idx, crm_qrad)
         do m = 1,crm_nz
            k = pver-m+1
            do i = 1,ncol
               icrm = ncol_sum + i
               crm_rad%qrad(icrm,:,:,m) = crm_qrad(i,:,:,m) / state(c)%pdel(i,k)
            end do
         end do

         !------------------------------------------------------------------------------------------
         ! calculate total water before calling crm - used for check_energy_chng() after CRM
         !------------------------------------------------------------------------------------------
         do i = 1,ncol
            icrm = ncol_sum + i

            qli_hydro_before(c,i) = 0.0_r8
            qi_hydro_before(c,i) = 0.0_r8

            do m = 1,crm_nz
               k = pver-m+1
               dp_g = state(c)%pdel(i,k)/gravit
               do jj = 1,crm_ny
                  do ii = 1,crm_nx
                     if (MMF_microphysics_scheme .eq. 'p3') then
                        qli_hydro_before(c,i) = qli_hydro_before(c,i)+(crm_state%qr(icrm,ii,jj,m)) * dp_g
                        ! P3 does not have a separate category for snow, so qi_hydro_before should be zero
                     end if
                     if (MMF_microphysics_scheme .eq. 'sam1mom') then
                        sfactor = max(0._r8,min(1._r8,(crm_state%temperature(icrm,ii,jj,m)-268.16)*1./(283.16-268.16)))
                        qli_hydro_before(c,i) = qli_hydro_before(c,i)+crm_state%qp(icrm,ii,jj,m) * dp_g
                        qi_hydro_before(c,i)  =  qi_hydro_before(c,i)+crm_state%qp(icrm,ii,jj,m) * (1-sfactor) * dp_g
                     end if ! MMF_microphysics_scheme
                  end do ! ii
               end do ! jj
            end do ! m

            qli_hydro_before(c,i) = qli_hydro_before(c,i)/(crm_nx*crm_ny)
            qi_hydro_before(c,i)  =  qi_hydro_before(c,i)/(crm_nx*crm_ny)
         end do ! i = 1,ncol

         !------------------------------------------------------------------------------------------
         ! Set CRM inputs
         !------------------------------------------------------------------------------------------
         ! TODO: move this to a routine and call like: call set_crm_input(...)
         do i = 1,ncol
            icrm = ncol_sum + i
            crm_input%zmid(icrm,1:pver)   = state(c)%zm(i,1:pver)
            crm_input%zint(icrm,1:pver+1) = state(c)%zi(i,1:pver+1)
            crm_input%tl(icrm,1:pver)     = state(c)%t(i,1:pver)
            crm_input%ql(icrm,1:pver)     = state(c)%q(i,1:pver,1)
            crm_input%qccl(icrm,1:pver)   = state(c)%q(i,1:pver,ixcldliq)
            crm_input%qiil(icrm,1:pver)   = state(c)%q(i,1:pver,ixcldice)
            crm_input%ps(icrm)            = state(c)%ps(i)
            crm_input%pmid(icrm,1:pver)   = state(c)%pmid(i,1:pver)
            crm_input%pint(icrm,1:pver+1) = state(c)%pint(i,1:pver+1)
            crm_input%pdel(icrm,1:pver)   = state(c)%pdel(i,1:pver)
            crm_input%phis(icrm)          = state(c)%phis(i)
            crm_input%ul(icrm,1:pver)     = state(c)%u(i,1:pver)
            crm_input%vl(icrm,1:pver)     = state(c)%v(i,1:pver)
            crm_input%ocnfrac(icrm)       = cam_in(c)%ocnfrac(i)
            ! input wind for ESMT
            crm_input%ul_esmt(icrm,1:pver) = state(c)%u(i,1:pver)
            crm_input%vl_esmt(icrm,1:pver) = state(c)%v(i,1:pver)
            ! Variance transport
            if (use_MMF_VT) then
               crm_input%t_vt(icrm,:pver) = state(c)%q(i,:pver,idx_vt_t)
               crm_input%q_vt(icrm,:pver) = state(c)%q(i,:pver,idx_vt_q)
               crm_input%u_vt(icrm,:pver) = state(c)%q(i,:pver,idx_vt_u)
            end if
            ! Set the input wind (also sets CRM orientation)
            do k = 1,pver
               crm_input%ul(icrm,k) = state(c)%u(i,k) * cos( crm_angle(i) ) + state(c)%v(i,k) * sin( crm_angle(i) )
               crm_input%vl(icrm,k) = state(c)%v(i,k) * cos( crm_angle(i) ) - state(c)%u(i,k) * sin( crm_angle(i) )
            end do ! k=1,pver
         end do ! i=1,ncol

         !------------------------------------------------------------------------------------------
         ! P3 input data
         !------------------------------------------------------------------------------------------
         if (MMF_microphysics_scheme .eq. 'p3') then
            ! NOTE - nccn_prescribed is only used when do_prescribed_CCN=true
            ! if do_prescribed_CCN=false then nc_nuceat_tend and ni_activated
            ! must be provided by an aerosol activation scheme
            ! TODO: build interface for SPA to get nccn_prescribed consistent with SCREAM
            do i = 1,ncol
               icrm = ncol_sum + i
               do k = 1, pver
                  crm_input%nccn_prescribed(icrm,k) = 50e6 ! 50 [#/cm3] 1e6 [cm3/m3] / (kg/m3) => 50e6 [#/kg]
                  crm_input%nc_nuceat_tend (icrm,k) = 0.0 ! [#/kg/s]
                  crm_input%ni_activated   (icrm,k) = 0.0 ! [#/kg]
               end do
            end do
         end if

         !------------------------------------------------------------------------------------------
         ! Set surface flux variables
         !------------------------------------------------------------------------------------------
         if (phys_do_flux_avg()) then
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('SHFLX'), shf_ptr)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('LHFLX'), lhf_ptr)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('TAUX'),  wsx_ptr)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('TAUY'),  wsy_ptr)
            shf_tmp = shf_ptr
            lhf_tmp = lhf_ptr
            wsx_tmp = wsx_ptr
            wsy_tmp = wsy_ptr
         else
            shf_tmp = cam_in(c)%shf
            lhf_tmp = cam_in(c)%lhf
            wsx_tmp = cam_in(c)%wsx
            wsy_tmp = cam_in(c)%wsy
         end if

         do i = 1,ncol
            icrm = ncol_sum + i
            crm_input%tau00(icrm)   = sqrt(wsx_tmp(i)**2 + wsy_tmp(i)**2)
            crm_input%bflxls(icrm)  = shf_tmp(i)/cpair + 0.61*state(c)%t(i,pver)*lhf_tmp(i)/latvap
            crm_input%fluxu00(icrm) = wsx_tmp(i) ! N/m2
            crm_input%fluxv00(icrm) = wsy_tmp(i) ! N/m2
            crm_input%fluxt00(icrm) = shf_tmp(i) ! W/m2  ( divide by cpair  to get K kg/(m2 s) )
            crm_input%fluxq00(icrm) = lhf_tmp(i) ! W/m2  ( divide by latvap to get   kg/(m2 s) )
            crm_input%wndls(icrm)   = sqrt(state(c)%u(i,pver)**2 + state(c)%v(i,pver)**2)
         end do

         !------------------------------------------------------------------------------------------
         ! Set aerosol
         !------------------------------------------------------------------------------------------
#if defined(MODAL_AERO)
         phase = 1  ! interstital aerosols only
         do i = 1,ncol
            icrm = ncol_sum + i
            air_density(i,1:pver) = state(c)%pmid(i,1:pver) / (287.15*state(c)%t(i,1:pver))
            do k = 1, pver
               do m = 1, ntot_amode
                 call loadaer( state(c), pbuf_chunk, i, i, k, m, air_density, phase, &
                               aerosol_num, aerosol_vol, aerosol_hygro)
                 crm_input%naermod (icrm,k,m) = aerosol_num(i)
                 crm_input%vaerosol(icrm,k,m) = aerosol_vol(i)
                 crm_input%hygro   (icrm,k,m) = aerosol_hygro(i)
               end do    
            end do
         end do
#endif
         !------------------------------------------------------------------------------------------
         !------------------------------------------------------------------------------------------

         ncol_sum = ncol_sum + ncol

      end do ! c=begchunk, endchunk

      !---------------------------------------------------------------------------------------------
      ! Run the CRM
      !---------------------------------------------------------------------------------------------
      if (.not.allocated(crm_clear_rh)) allocate(crm_clear_rh(ncrms,crm_nz))

      ! Load latitude, longitude, and unique column ID for all CRMs
      allocate(longitude0(ncrms))
      allocate(latitude0 (ncrms))
      allocate(gcolp     (ncrms))
      ncol_sum = 0
      do c=begchunk, endchunk
         ncol = state(c)%ncol
         do i = 1 , ncol
           icrm = ncol_sum + i
           latitude0 (icrm) = get_rlat_p(c,i) * 57.296_r8
           longitude0(icrm) = get_rlon_p(c,i) * 57.296_r8
           gcolp     (icrm) = get_gcol_p(c,i)
         enddo
         ncol_sum = ncol_sum + ncol
      end do ! c=begchunk, endchunk

#if defined(MMF_SAM) || defined(MMF_SAMOMP)
      
      call t_startf ('crm_call')
      call crm(ncrms, ztodt, pver, &
               crm_input, crm_state, crm_rad, &
               crm_ecpp_output, crm_output, crm_clear_rh, &
               latitude0, longitude0, gcolp, nstep, &
               use_MMF_VT_tmp, MMF_VT_wn_max, &
               use_crm_accel_tmp, crm_accel_factor, crm_accel_uv_tmp)
      call t_stopf('crm_call')
      
#elif defined(MMF_SAMXX)

      ! Fortran classes don't translate to C++ classes, we we have to separate
      ! this stuff out when calling the C++ routinte crm(...)
      call t_startf ('crm_call')
      call crm(ncrms, ncrms, ztodt, pver, crm_input%bflxls, crm_input%wndls, crm_input%zmid, crm_input%zint, &
               crm_input%pmid, crm_input%pint, crm_input%pdel, crm_input%ul, crm_input%vl, &
               crm_input%tl, crm_input%qccl, crm_input%qiil, crm_input%ql, crm_input%tau00, &
               crm_input%ul_esmt, crm_input%vl_esmt,                                        &
               crm_input%t_vt, crm_input%q_vt, crm_input%u_vt, &
               crm_state%u_wind, crm_state%v_wind, crm_state%w_wind, crm_state%temperature, &
               crm_state%qv, crm_state%qp, crm_state%qn, crm_rad%qrad, crm_rad%temperature, &
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
               crm_output%t_vt_tend, crm_output%q_vt_tend, crm_output%u_vt_tend, &
               crm_output%t_vt_ls, crm_output%q_vt_ls, crm_output%u_vt_ls, &
               crm_output%ultend, crm_output%vltend, &
               crm_output%tk, crm_output%tkh, crm_output%qcl, crm_output%qci, crm_output%qpl, crm_output%qpi, &
               crm_output%z0m, crm_output%taux, crm_output%tauy, crm_output%precc, crm_output%precl, crm_output%precsc, &
               crm_output%precsl, crm_output%prec_crm,                         &
               crm_clear_rh, &
               latitude0, longitude0, gcolp, nstep, &
               use_MMF_VT, MMF_VT_wn_max, use_MMF_ESMT, &
               use_crm_accel, crm_accel_factor, crm_accel_uv)
      call t_stopf('crm_call')

#elif defined(MMF_PAM)

      call pam_mirror_array_readonly( 'latitude',      latitude0   )
      call pam_mirror_array_readonly( 'longitude',     longitude0  )

      call pam_mirror_array_readonly( 'input_bflxls',  crm_input%bflxls  )
      call pam_mirror_array_readonly( 'input_wndls',   crm_input%wndls   )

      call pam_mirror_array_readonly( 'input_shf',     crm_input%fluxt00 )
      call pam_mirror_array_readonly( 'input_lhf',     crm_input%fluxq00 )

      call pam_mirror_array_readonly( 'input_zmid',    crm_input%zmid    )
      call pam_mirror_array_readonly( 'input_zint',    crm_input%zint    )
      call pam_mirror_array_readonly( 'input_pmid',    crm_input%pmid    )
      call pam_mirror_array_readonly( 'input_pint',    crm_input%pint    )
      call pam_mirror_array_readonly( 'input_pdel',    crm_input%pdel    )
      call pam_mirror_array_readonly( 'input_phis',    crm_input%phis    )

      call pam_mirror_array_readonly( 'input_ul',      crm_input%ul      )
      call pam_mirror_array_readonly( 'input_vl',      crm_input%vl      )
      call pam_mirror_array_readonly( 'input_tl',      crm_input%tl      )
      call pam_mirror_array_readonly( 'input_qccl',    crm_input%qccl    )
      call pam_mirror_array_readonly( 'input_qiil',    crm_input%qiil    )
      call pam_mirror_array_readonly( 'input_ql',      crm_input%ql      )
      call pam_mirror_array_readonly( 'input_tau00',   crm_input%tau00   )
      ! call pam_mirror_array_readonly( 'input_ul_esmt', crm_input%ul_esmt )
      ! call pam_mirror_array_readonly( 'input_vl_esmt', crm_input%vl_esmt )
      call pam_mirror_array_readonly( 'input_nccn_prescribed',crm_input%nccn_prescribed )
      call pam_mirror_array_readonly( 'input_nc_nuceat_tend', crm_input%nc_nuceat_tend  )
      call pam_mirror_array_readonly( 'input_ni_activated',   crm_input%ni_activated    )

      ! Variance transport inputs and outputs
      if (use_MMF_VT) then
         call pam_mirror_array_readonly(  'input_vt_t',       crm_input%t_vt )
         call pam_mirror_array_readonly(  'input_vt_q',       crm_input%q_vt )
         call pam_mirror_array_readonly(  'input_vt_u',       crm_input%u_vt )
         call pam_mirror_array_readwrite( 'output_t_vt_tend', crm_output%t_vt_tend )
         call pam_mirror_array_readwrite( 'output_q_vt_tend', crm_output%q_vt_tend )
         call pam_mirror_array_readwrite( 'output_u_vt_tend', crm_output%u_vt_tend )
         call pam_mirror_array_readwrite( 'output_t_vt_ls',   crm_output%t_vt_ls )
         call pam_mirror_array_readwrite( 'output_q_vt_ls',   crm_output%q_vt_ls )
         call pam_mirror_array_readwrite( 'output_u_vt_ls',   crm_output%u_vt_ls )
      end if

      call pam_mirror_array_readwrite( 'state_u_wind',      crm_state%u_wind      )
      call pam_mirror_array_readwrite( 'state_v_wind',      crm_state%v_wind      )
      call pam_mirror_array_readwrite( 'state_w_wind',      crm_state%w_wind      )
      call pam_mirror_array_readwrite( 'state_temperature', crm_state%temperature )
      call pam_mirror_array_readwrite( 'state_rho_dry',     crm_state%rho_dry     )
      call pam_mirror_array_readwrite( 'state_qv',          crm_state%qv          )
      call pam_mirror_array_readwrite( 'state_qc',          crm_state%qc          )
      call pam_mirror_array_readwrite( 'state_nc',          crm_state%nc          )
      call pam_mirror_array_readwrite( 'state_qr',          crm_state%qr          )
      call pam_mirror_array_readwrite( 'state_nr',          crm_state%nr          )
      call pam_mirror_array_readwrite( 'state_qi',          crm_state%qi          )
      call pam_mirror_array_readwrite( 'state_ni',          crm_state%ni          )
      call pam_mirror_array_readwrite( 'state_qm',          crm_state%qm          )
      call pam_mirror_array_readwrite( 'state_bm',          crm_state%bm          )
      call pam_mirror_array_readwrite( 'state_t_prev',      crm_state%t_prev      )
      call pam_mirror_array_readwrite( 'state_q_prev',      crm_state%q_prev      )

      ! SHOC variables
      call pam_mirror_array_readwrite( 'state_shoc_tk',      crm_state%shoc_tk       )
      call pam_mirror_array_readwrite( 'state_shoc_tkh',     crm_state%shoc_tkh      )
      call pam_mirror_array_readwrite( 'state_shoc_wthv',    crm_state%shoc_wthv )
      call pam_mirror_array_readwrite( 'state_shoc_relvar',  crm_state%shoc_relvar   )
      call pam_mirror_array_readwrite( 'state_shoc_cldfrac', crm_state%shoc_cldfrac  )

      ! Radiation tendency and output conditions
      call pam_mirror_array_readwrite( 'rad_qrad',        crm_rad%qrad        )
      call pam_mirror_array_readwrite( 'rad_temperature', crm_rad%temperature )
      call pam_mirror_array_readwrite( 'rad_qv',          crm_rad%qv          )
      call pam_mirror_array_readwrite( 'rad_qc',          crm_rad%qc          )
      call pam_mirror_array_readwrite( 'rad_qi',          crm_rad%qi          )
      call pam_mirror_array_readwrite( 'rad_nc',          crm_rad%nc          )
      call pam_mirror_array_readwrite( 'rad_ni',          crm_rad%ni          )
      call pam_mirror_array_readwrite( 'rad_cld',         crm_rad%cld         )

      call pam_mirror_array_readwrite( 'output_clear_rh', crm_clear_rh        )

      ! call pam_mirror_array_readwrite( 'output_prectend',        crm_output%prectend,        '' )
      ! call pam_mirror_array_readwrite( 'output_precstend',       crm_output%precstend,       '' )
      call pam_mirror_array_readwrite( 'output_cld',             crm_output%cld,             '' )
      ! call pam_mirror_array_readwrite( 'output_cldtop',          crm_output%cldtop,          '' )
      call pam_mirror_array_readwrite( 'output_gicewp',          crm_output%gicewp,          '' )
      call pam_mirror_array_readwrite( 'output_gliqwp',          crm_output%gliqwp,          '' )
      call pam_mirror_array_readwrite( 'output_liq_ice_exchange',crm_output%liq_ice_exchange,'' )
      call pam_mirror_array_readwrite( 'output_vap_liq_exchange',crm_output%vap_liq_exchange,'' )
      call pam_mirror_array_readwrite( 'output_vap_ice_exchange',crm_output%vap_ice_exchange,'' )
      ! call pam_mirror_array_readwrite( 'output_mctot',           crm_output%mctot,           '' )
      ! call pam_mirror_array_readwrite( 'output_mcup',            crm_output%mcup,            '' )
      ! call pam_mirror_array_readwrite( 'output_mcdn',            crm_output%mcdn,            '' )
      ! call pam_mirror_array_readwrite( 'output_mcuup',           crm_output%mcuup,           '' )
      ! call pam_mirror_array_readwrite( 'output_mcudn',           crm_output%mcudn,           '' )

      call pam_mirror_array_readwrite( 'output_qv_mean',         crm_output%qv_mean,         '' )
      call pam_mirror_array_readwrite( 'output_qc_mean',         crm_output%qc_mean,         '' )
      call pam_mirror_array_readwrite( 'output_qr_mean',         crm_output%qr_mean,         '' )
      call pam_mirror_array_readwrite( 'output_qi_mean',         crm_output%qi_mean,         '' )
      call pam_mirror_array_readwrite( 'output_nc_mean',         crm_output%nc_mean,         '' )
      call pam_mirror_array_readwrite( 'output_nr_mean',         crm_output%nr_mean,         '' )
      call pam_mirror_array_readwrite( 'output_ni_mean',         crm_output%ni_mean,         '' )
      call pam_mirror_array_readwrite( 'output_qm_mean',         crm_output%qm_mean,         '' )
      call pam_mirror_array_readwrite( 'output_bm_mean',         crm_output%bm_mean,         '' )
      call pam_mirror_array_readwrite( 'output_rho_d_mean',      crm_output%rho_d_mean,         '' )
      call pam_mirror_array_readwrite( 'output_rho_v_mean',      crm_output%rho_v_mean,         '' )

      call pam_mirror_array_readwrite( 'output_qt_ls',       crm_output%qt_ls,       '' )
      call pam_mirror_array_readwrite( 'output_t_ls',        crm_output%t_ls,        '' )
      call pam_mirror_array_readwrite( 'output_rho_v_ls',    crm_output%rho_v_ls,    '' )
      call pam_mirror_array_readwrite( 'output_rho_d_ls',    crm_output%rho_d_ls,    '' )
      call pam_mirror_array_readwrite( 'output_rho_l_ls',    crm_output%rho_l_ls,    '' )
      call pam_mirror_array_readwrite( 'output_rho_i_ls',    crm_output%rho_i_ls,    '' )
      ! call pam_mirror_array_readwrite( 'output_jt_crm',      crm_output%jt_crm,      '' )
      ! call pam_mirror_array_readwrite( 'output_mx_crm',      crm_output%mx_crm,      '' )
      ! call pam_mirror_array_readwrite( 'output_cltot',       crm_output%cltot,       '' )
      ! call pam_mirror_array_readwrite( 'output_clhgh',       crm_output%clhgh,       '' )
      ! call pam_mirror_array_readwrite( 'output_clmed',       crm_output%clmed,       '' )
      ! call pam_mirror_array_readwrite( 'output_cllow',       crm_output%cllow,       '' )
      call pam_mirror_array_readwrite( 'output_sltend',      crm_output%sltend,      '' )
      call pam_mirror_array_readwrite( 'output_qltend',      crm_output%qltend,      '' )
      call pam_mirror_array_readwrite( 'output_qcltend',     crm_output%qcltend,     '' )
      call pam_mirror_array_readwrite( 'output_qiltend',     crm_output%qiltend,     '' )
      call pam_mirror_array_readwrite( 'output_ultend',      crm_output%ultend,      '' )
      call pam_mirror_array_readwrite( 'output_vltend',      crm_output%vltend,      '' )
      ! call pam_mirror_array_readwrite( 'output_tk',          crm_output%tk,          '' )
      ! call pam_mirror_array_readwrite( 'output_tkh',         crm_output%tkh,         '' )
      ! call pam_mirror_array_readwrite( 'output_qcl',         crm_output%qcl,         '' )
      ! call pam_mirror_array_readwrite( 'output_qci',         crm_output%qci,         '' )
      ! call pam_mirror_array_readwrite( 'output_qpl',         crm_output%qpl,         '' )
      ! call pam_mirror_array_readwrite( 'output_qpi',         crm_output%qpi,         '' )
      call pam_mirror_array_readwrite( 'output_precc',       crm_output%precc,       '' )
      call pam_mirror_array_readwrite( 'output_precl',       crm_output%precl,       '' )
      call pam_mirror_array_readwrite( 'output_precsc',      crm_output%precsc,      '' )
      call pam_mirror_array_readwrite( 'output_precsl',      crm_output%precsl,      '' )
      call pam_mirror_array_readwrite( 'output_prec_crm',    crm_output%prec_crm,    '' )

      call pam_mirror_array_readwrite( 'output_dt_sgs',      crm_output%dt_sgs,      '' )
      call pam_mirror_array_readwrite( 'output_dqv_sgs',     crm_output%dqv_sgs,     '' )
      call pam_mirror_array_readwrite( 'output_dqc_sgs',     crm_output%dqc_sgs,     '' )
      call pam_mirror_array_readwrite( 'output_dqi_sgs',     crm_output%dqi_sgs,     '' )
      call pam_mirror_array_readwrite( 'output_dqr_sgs',     crm_output%dqr_sgs,     '' )

      call pam_mirror_array_readwrite( 'output_dt_micro',    crm_output%dt_micro,    '' )
      call pam_mirror_array_readwrite( 'output_dqv_micro',   crm_output%dqv_micro,   '' )
      call pam_mirror_array_readwrite( 'output_dqc_micro',   crm_output%dqc_micro,   '' )
      call pam_mirror_array_readwrite( 'output_dqi_micro',   crm_output%dqi_micro,   '' )
      call pam_mirror_array_readwrite( 'output_dqr_micro',   crm_output%dqr_micro,   '' )

      call pam_mirror_array_readwrite( 'output_dt_dycor',   crm_output%dt_dycor,     '' )
      call pam_mirror_array_readwrite( 'output_dqv_dycor',  crm_output%dqv_dycor,    '' )
      call pam_mirror_array_readwrite( 'output_dqc_dycor',  crm_output%dqc_dycor,    '' )
      call pam_mirror_array_readwrite( 'output_dqi_dycor',  crm_output%dqi_dycor,    '' )
      call pam_mirror_array_readwrite( 'output_dqr_dycor',  crm_output%dqr_dycor,    '' )

      call pam_mirror_array_readwrite( 'output_dt_sponge',  crm_output%dt_sponge,    '' )
      call pam_mirror_array_readwrite( 'output_dqv_sponge', crm_output%dqv_sponge,   '' )
      call pam_mirror_array_readwrite( 'output_dqc_sponge', crm_output%dqc_sponge,   '' )
      call pam_mirror_array_readwrite( 'output_dqi_sponge', crm_output%dqi_sponge,   '' )
      call pam_mirror_array_readwrite( 'output_dqr_sponge', crm_output%dqr_sponge,   '' )

      call pam_mirror_array_readonly( 'global_column_id', gcolp )

      call pam_set_option('ncrms', ncrms )
      call pam_set_option('gcm_nlev', pver )
      call pam_set_option('crm_nz',crm_nz )
      call pam_set_option('crm_nx',crm_nx )
      call pam_set_option('crm_ny',crm_ny )
      call pam_set_option('rad_nx',crm_nx_rad)
      call pam_set_option('rad_ny',crm_ny_rad)
      call pam_set_option('crm_dx',crm_dx )
      call pam_set_option('crm_dy',crm_dy )
      call pam_set_option('gcm_dt',ztodt )
      call pam_set_option('crm_dt',crm_dt )

      call pam_register_dimension('gcm_lev',pver)

      call pam_set_option('use_MMF_VT', use_MMF_VT_tmp )
      call pam_set_option('use_MMF_ESMT', use_MMF_ESMT_tmp )

      call pam_set_option('use_crm_accel', use_crm_accel_tmp )
      call pam_set_option('crm_accel_uv', crm_accel_uv_tmp)
      call pam_set_option('crm_accel_factor', crm_accel_factor )

      call pam_set_option('enable_physics_tend_stats', .false. )

      call pam_set_option('is_first_step', (nstep<=1) )
      call pam_set_option('is_restart', (nsrest>0) )
      call pam_set_option('am_i_root', masterproc )

      call t_startf ('crm_call')
      call pam_driver()
      call t_stopf('crm_call')

#endif

      deallocate(longitude0)
      deallocate(latitude0 )
      deallocate(gcolp     )

      !---------------------------------------------------------------------------------------------
      ! Deal with CRM outputs
      !---------------------------------------------------------------------------------------------

      ! There is no separate convective and stratiform precip for CRM,
      ! so make it all convective and zero out the stratiform
      crm_output%precc (1:ncrms) = crm_output%precc (1:ncrms) + crm_output%precl (1:ncrms)
      crm_output%precsc(1:ncrms) = crm_output%precsc(1:ncrms) + crm_output%precsl(1:ncrms)
      crm_output%precl (1:ncrms) = 0
      crm_output%precsl(1:ncrms) = 0

      ncol_sum = 0
      do c=begchunk, endchunk
         pbuf_chunk => pbuf_get_chunk(pbuf2d, c)
         ncol = state(c)%ncol

         call pbuf_get_field(pbuf_chunk, crm_angle_idx, crm_angle)

         !------------------------------------------------------------------------------------------
         ! Set ptend logicals with physics tendencies from CRM
         !------------------------------------------------------------------------------------------
         ptend(c)%name = 'crm'
         ptend(c)%ls           = .TRUE.
         ptend(c)%lq(1)        = .TRUE.
         ptend(c)%lq(ixcldliq) = .TRUE.
         ptend(c)%lq(ixcldice) = .TRUE.
         ptend(c)%lu           = .FALSE.
         ptend(c)%lv           = .FALSE.

         if (use_MMF_VT) then
            ptend(c)%lq(idx_vt_t) = .TRUE.
            ptend(c)%lq(idx_vt_q) = .TRUE.
            ptend(c)%lq(idx_vt_u) = .TRUE.
         end if
         !------------------------------------------------------------------------------------------
         ! Populate output tendencies from CRM
         !------------------------------------------------------------------------------------------
         do i = 1,ncol
            icrm = ncol_sum + i
            ptend(c)%s(i,1:pver)          = crm_output%sltend (icrm,1:pver)
            ptend(c)%q(i,1:pver,1)        = crm_output%qltend (icrm,1:pver)
            ptend(c)%q(i,1:pver,ixcldliq) = crm_output%qcltend(icrm,1:pver)
            ptend(c)%q(i,1:pver,ixcldice) = crm_output%qiltend(icrm,1:pver)

            if (use_MMF_VT) then
               ptend(c)%q(i,1:pver,idx_vt_t) = crm_output%t_vt_tend(icrm,1:pver)
               ptend(c)%q(i,1:pver,idx_vt_q) = crm_output%q_vt_tend(icrm,1:pver)
               ptend(c)%q(i,1:pver,idx_vt_u) = crm_output%u_vt_tend(icrm,1:pver)
            end if
         end do ! i = 1,ncol

         !------------------------------------------------------------------------------------------
         ! Populate output tendencies for P3 microphysics
         !------------------------------------------------------------------------------------------
         if ( MMF_microphysics_scheme .eq. 'p3' ) then

            ptend(c)%lq(ixnumliq)  = .TRUE.
            ptend(c)%lq(ixnumice)  = .TRUE.
            ptend(c)%lq(ixrain)    = .TRUE. 
            ptend(c)%lq(ixnumrain) = .TRUE. 
            ptend(c)%lq(ixcldrim)  = .TRUE. 
            ptend(c)%lq(ixrimvol)  = .TRUE. 

            do i = 1, ncol
               do k = 1, crm_nz 
                  m = pver-k+1
                  ptend(c)%q(i,m,ixnumliq)  = 0
                  ptend(c)%q(i,m,ixnumice)  = 0
                  ptend(c)%q(i,m,ixrain)    = 0
                  ptend(c)%q(i,m,ixnumrain) = 0
                  ptend(c)%q(i,m,ixcldrim)  = 0
                  ptend(c)%q(i,m,ixrimvol)  = 0
               end do
            end do

            do i = 1, ncol
               icrm = ncol_sum + i
               do k = 1, crm_nz 
                  m = pver-k+1
                  do ii = 1, crm_nx
                  do jj = 1, crm_ny
                     ptend(c)%q(i,m,ixnumliq)  = ptend(c)%q(i,m,ixnumliq ) + crm_state%nc(icrm,ii,jj,k) 
                     ptend(c)%q(i,m,ixnumice)  = ptend(c)%q(i,m,ixnumice ) + crm_state%ni(icrm,ii,jj,k)
                     ptend(c)%q(i,m,ixrain)    = ptend(c)%q(i,m,ixrain   ) + crm_state%qr(icrm,ii,jj,k)
                     ptend(c)%q(i,m,ixnumrain) = ptend(c)%q(i,m,ixnumrain) + crm_state%nr(icrm,ii,jj,k)
                     ptend(c)%q(i,m,ixcldrim)  = ptend(c)%q(i,m,ixcldrim ) + crm_state%qm(icrm,ii,jj,k)
                     ptend(c)%q(i,m,ixrimvol)  = ptend(c)%q(i,m,ixrimvol ) + crm_state%bm(icrm,ii,jj,k)
                  end do
                  end do
                  ptend(c)%q(i,m,ixnumliq ) = ( ptend(c)%q(i,m,ixnumliq )/(crm_nx*crm_ny) - state(c)%q(i,m,ixnumliq ) )/ztodt
                  ptend(c)%q(i,m,ixnumice ) = ( ptend(c)%q(i,m,ixnumice )/(crm_nx*crm_ny) - state(c)%q(i,m,ixnumice ) )/ztodt
                  ptend(c)%q(i,m,ixrain   ) = ( ptend(c)%q(i,m,ixrain   )/(crm_nx*crm_ny) - state(c)%q(i,m,ixrain   ) )/ztodt
                  ptend(c)%q(i,m,ixnumrain) = ( ptend(c)%q(i,m,ixnumrain)/(crm_nx*crm_ny) - state(c)%q(i,m,ixnumrain) )/ztodt
                  ptend(c)%q(i,m,ixcldrim ) = ( ptend(c)%q(i,m,ixcldrim )/(crm_nx*crm_ny) - state(c)%q(i,m,ixcldrim ) )/ztodt
                  ptend(c)%q(i,m,ixrimvol ) = ( ptend(c)%q(i,m,ixrimvol )/(crm_nx*crm_ny) - state(c)%q(i,m,ixrimvol ) )/ztodt
               end do
            end do

         end if

         !------------------------------------------------------------------------------------------
         ! Add radiative heating tendency above CRM
         !------------------------------------------------------------------------------------------
         call pbuf_get_field(pbuf_chunk, pbuf_get_index('QRL'), qrl)
         call pbuf_get_field(pbuf_chunk, pbuf_get_index('QRS'), qrs)

         do k = 1,pver
            do i = 1,ncol
               qrs(i,k) = qrs(i,k)/state(c)%pdel(i,k)
               qrl(i,k) = qrl(i,k)/state(c)%pdel(i,k)
            end do
         end do

#if defined(MMF_PAM)
         ! Currently PAM is coupled at all levels, so only add rad to levels above CRM top
         ptend(c)%s(1:ncol, 1:pver-crm_nz) = qrs(1:ncol,1:pver-crm_nz) + qrl(1:ncol,1:pver-crm_nz)
#else
         ! The radiation tendencies in the GCM levels above the CRM and the top 2 CRM levels are set to
         ! be zero in the CRM, So add radiation tendencies to these levels 
         ptend(c)%s(1:ncol, 1:pver-crm_nz+2) = qrs(1:ncol,1:pver-crm_nz+2) + qrl(1:ncol,1:pver-crm_nz+2)
#endif

         ! This will be used to check energy conservation
         mmf_rad_flux(c,:ncol) = 0.0_r8
         do k = 1,pver
            do i = 1,ncol
               mmf_rad_flux(c,i) = mmf_rad_flux(c,i) + ( qrs(i,k) + qrl(i,k) ) * state(c)%pdel(i,k)/gravit
            end do
         end do

         !------------------------------------------------------------------------------------------
         ! set convective heating tendency for gravity wave drag
         !------------------------------------------------------------------------------------------
         ! Subtract radiation from all levels so only convective heating is used for the
         ! convective gravityw ave scheme, consistent with standard EAM w/ ZM
         if (ttend_dp_idx > 0) then
            call pbuf_get_field(pbuf_chunk, ttend_dp_idx, ttend_dp)
            ttend_dp(1:ncol,1:pver) = ( ptend(c)%s(1:ncol,:pver) - qrs(1:ncol,:pver) - qrl(1:ncol,:pver) )/cpair
         end if

         !------------------------------------------------------------------------------------------
         ! CRM cloud/precip output
         !------------------------------------------------------------------------------------------
         ! We need to do this here because we did not set crm_output%cld as a 
         ! pointer to pbuf, so we need to copy the data over. NOTE: I think this 
         ! can be done using pbuf_set_field without making an extra pointer for 
         ! cld, but I do not think we would be able to zero-out the rest of cld 
         ! beyond pcols that way.
         call pbuf_get_field(pbuf_chunk, cld_idx, cld, start=(/1,1,itim/), kount=(/pcols,pver,1/) )
         call pbuf_get_field(pbuf_chunk, prec_dp_idx,  prec_dp  )
         call pbuf_get_field(pbuf_chunk, snow_dp_idx,  snow_dp  )
         do i = 1,ncol
            icrm = ncol_sum + i
            cld(i,1:pver) = crm_output%cld(icrm,1:pver)
            prec_dp(i) = crm_output%precc(icrm)
            snow_dp(i) = crm_output%precsc(icrm)
         end do

         !------------------------------------------------------------------------------------------
         ! Output for ECPP
         !------------------------------------------------------------------------------------------
         if (use_ECPP) then

            do i = 1,ncol
               icrm = ncol_sum + i
               air_density(i,1:pver) = state(c)%pmid(i,1:pver) / (287.15*state(c)%t(i,1:pver))
               TKE_tmp(i,1:pver) = crm_output%tkez(icrm,1:pver) / air_density(i,1:pver)
               ideep_crm(i) = i*1.0  ! For convective transport
            end do

            call pbuf_set_field(pbuf_chunk, pbuf_get_index('TKE_CRM'), TKE_tmp )
            call pbuf_set_field(pbuf_chunk, pbuf_get_index('TK_CRM'), crm_output%tkz   (ncol_sum+1:ncol_sum+ncol,1:pver) )
            call pbuf_set_field(pbuf_chunk, pbuf_get_index('MU_CRM'), crm_output%mu_crm(ncol_sum+1:ncol_sum+ncol,1:pver) )
            call pbuf_set_field(pbuf_chunk, pbuf_get_index('MD_CRM'), crm_output%md_crm(ncol_sum+1:ncol_sum+ncol,1:pver) )
            call pbuf_set_field(pbuf_chunk, pbuf_get_index('EU_CRM'), crm_output%eu_crm(ncol_sum+1:ncol_sum+ncol,1:pver) )
            call pbuf_set_field(pbuf_chunk, pbuf_get_index('DU_CRM'), crm_output%du_crm(ncol_sum+1:ncol_sum+ncol,1:pver) )
            call pbuf_set_field(pbuf_chunk, pbuf_get_index('ED_CRM'), crm_output%eu_crm(ncol_sum+1:ncol_sum+ncol,1:pver) )
            call pbuf_set_field(pbuf_chunk, pbuf_get_index('JT_CRM'), crm_output%jt_crm(ncol_sum+1:ncol_sum+ncol) )
            call pbuf_set_field(pbuf_chunk, pbuf_get_index('MX_CRM'), crm_output%mx_crm(ncol_sum+1:ncol_sum+ncol) )
            call pbuf_set_field(pbuf_chunk, pbuf_get_index('IDEEP_CRM'), ideep_crm )

            crm_ecpp_output%qlsinkcen = crm_ecpp_output%qlsink_avgcen
         
         end if ! use_ECPP

         !------------------------------------------------------------------------------------------
         ! CRM momentum tendencies
         !------------------------------------------------------------------------------------------

         if (use_MMF_ESMT) then
            ptend(c)%lu = .TRUE.
            ptend(c)%lv = .TRUE.
            do i = 1,ncol
               icrm = ncol_sum + i
               ptend(c)%u(i,1:pver)  = crm_output%ultend(icrm,1:pver)
               ptend(c)%v(i,1:pver)  = crm_output%vltend(icrm,1:pver)
            end do
         end if

#if defined(MMF_MOMENTUM_FEEDBACK)
         ptend(c)%lu = .TRUE.
         ptend(c)%lv = .TRUE.
         ! rotate resolved CRM momentum tendencies back
         do i = 1, ncol 
            icrm = ncol_sum + i
            ptend(c)%u(i,1:pver) = crm_output%ultend(icrm,1:pver) * cos( -1.*crm_angle(i) ) + crm_output%vltend(icrm,1:pver) * sin( -1.*crm_angle(i) )
            ptend(c)%v(i,1:pver) = crm_output%vltend(icrm,1:pver) * cos( -1.*crm_angle(i) ) - crm_output%ultend(icrm,1:pver) * sin( -1.*crm_angle(i) )
         end do
#endif /* MMF_MOMENTUM_FEEDBACK */

         !------------------------------------------------------------------------------------------
         ! Write out data for history files
         !------------------------------------------------------------------------------------------
         icrm_beg = ncol_sum + 1
         icrm_end = ncol_sum + ncol
         call crm_history_out(state(c), ptend(c), crm_state, crm_rad, crm_output, &
                              crm_ecpp_output, qrs, qrl, icrm_beg, icrm_end)

         !------------------------------------------------------------------------------------------
         ! Convert heating rate to Q*dp to conserve energy across timesteps
         ! and put rad data back in pbuf
         !------------------------------------------------------------------------------------------

         call pbuf_get_field(pbuf_chunk, crm_t_rad_idx,    crm_t_rad)
         call pbuf_get_field(pbuf_chunk, crm_qv_rad_idx,   crm_qv_rad)
         call pbuf_get_field(pbuf_chunk, crm_qc_rad_idx,   crm_qc_rad)
         call pbuf_get_field(pbuf_chunk, crm_qi_rad_idx,   crm_qi_rad)
         call pbuf_get_field(pbuf_chunk, crm_cld_rad_idx,  crm_cld_rad)
         call pbuf_get_field(pbuf_chunk, crm_qrad_idx, crm_qrad)

         if (MMF_microphysics_scheme .eq. 'p3') then
            call pbuf_get_field(pbuf_chunk, crm_nc_rad_idx, crm_nc_rad)
            call pbuf_get_field(pbuf_chunk, crm_ni_rad_idx, crm_ni_rad)
         end if

         do i = 1,ncol
            icrm = ncol_sum + i
            do m = 1,crm_nz
               k = pver-m+1
               crm_rad%qrad(icrm,:,:,m) = crm_rad%qrad(icrm,:,:,m) * state(c)%pdel(i,k) ! for energy conservation
            end do
            crm_qrad     (i,:,:,:) = crm_rad%qrad       (icrm,:,:,:)
            crm_t_rad    (i,:,:,:) = crm_rad%temperature(icrm,:,:,:)
            crm_qv_rad   (i,:,:,:) = crm_rad%qv         (icrm,:,:,:)
            crm_qc_rad   (i,:,:,:) = crm_rad%qc         (icrm,:,:,:)
            crm_qi_rad   (i,:,:,:) = crm_rad%qi         (icrm,:,:,:)
            crm_cld_rad  (i,:,:,:) = crm_rad%cld        (icrm,:,:,:)
            if (MMF_microphysics_scheme .eq. 'p3') then
               crm_nc_rad(i,:,:,:) = crm_rad%nc         (icrm,:,:,:)
               crm_ni_rad(i,:,:,:) = crm_rad%ni         (icrm,:,:,:)
            end if
         end do

         !------------------------------------------------------------------------------------------
         ! put CRM state data back in pbuf
         !------------------------------------------------------------------------------------------
         call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_U'),  crm_u)
         call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_V'),  crm_v)
         call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_W'),  crm_w)
         call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_T'),  crm_t)
         call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_RHO'),  crm_rho_d)
         call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QV'), crm_qv)
         if (MMF_microphysics_scheme .eq. 'sam1mom') then
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QP'), crm_qp)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QN'), crm_qn)
         end if
         if (MMF_microphysics_scheme .eq. 'p3') then
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_NC'), crm_nc)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QR'), crm_qr)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_NR'), crm_nr)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QI'), crm_qi)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_NI'), crm_ni)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QC'), crm_qc)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QM'), crm_qm)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_BM'), crm_bm)
            call pbuf_get_field(pbuf_chunk, crm_t_prev_idx, crm_t_prev)
            call pbuf_get_field(pbuf_chunk, crm_q_prev_idx, crm_q_prev)
            call pbuf_get_field(pbuf_chunk, crm_shoc_tk_idx     , crm_shoc_tk)
            call pbuf_get_field(pbuf_chunk, crm_shoc_tkh_idx    , crm_shoc_tkh)
            call pbuf_get_field(pbuf_chunk, crm_shoc_wthv_idx   , crm_shoc_wthv)
            call pbuf_get_field(pbuf_chunk, crm_shoc_relvar_idx , crm_shoc_relvar)
            call pbuf_get_field(pbuf_chunk, crm_shoc_cldfrac_idx, crm_shoc_cldfrac)
         end if

         do i = 1,ncol
            icrm = ncol_sum + i
            crm_u (i,:,:,:) = crm_state%u_wind     (icrm,:,:,:)
            crm_v (i,:,:,:) = crm_state%v_wind     (icrm,:,:,:)
            crm_w (i,:,:,:) = crm_state%w_wind     (icrm,:,:,:)
            crm_t (i,:,:,:) = crm_state%temperature(icrm,:,:,:)
            crm_rho_d(i,:,:,:) = crm_state%rho_dry (icrm,:,:,:)
            crm_qv(i,:,:,:) = crm_state%qv         (icrm,:,:,:)
            if (MMF_microphysics_scheme .eq. 'sam1mom') then
               crm_qp(i,:,:,:) = crm_state%qp(icrm,:,:,:)
               crm_qn(i,:,:,:) = crm_state%qn(icrm,:,:,:)
            end if
            if (MMF_microphysics_scheme .eq. 'p3') then
               crm_qc(i,:,:,:)     = crm_state%qc    (icrm,:,:,:)
               crm_qi(i,:,:,:)     = crm_state%qi    (icrm,:,:,:)
               crm_qr(i,:,:,:)     = crm_state%qr    (icrm,:,:,:)
               crm_nc(i,:,:,:)     = crm_state%nc    (icrm,:,:,:)
               crm_ni(i,:,:,:)     = crm_state%ni    (icrm,:,:,:)
               crm_nr(i,:,:,:)     = crm_state%nr    (icrm,:,:,:)
               crm_qm(i,:,:,:)     = crm_state%qm    (icrm,:,:,:)
               crm_bm(i,:,:,:)     = crm_state%bm    (icrm,:,:,:)
               crm_t_prev(i,:,:,:) = crm_state%t_prev(icrm,:,:,:)
               crm_q_prev(i,:,:,:) = crm_state%q_prev(icrm,:,:,:)
               crm_shoc_tk     (i,:,:,:) = crm_state%shoc_tk     (icrm,:,:,:)
               crm_shoc_tkh    (i,:,:,:) = crm_state%shoc_tkh    (icrm,:,:,:)
               crm_shoc_wthv   (i,:,:,:) = crm_state%shoc_wthv   (icrm,:,:,:)
               crm_shoc_relvar (i,:,:,:) = crm_state%shoc_relvar (icrm,:,:,:)
               crm_shoc_cldfrac(i,:,:,:) = crm_state%shoc_cldfrac(icrm,:,:,:)
            end if
         end do

         !------------------------------------------------------------------------------------------
         ! calculate column integrated water for energy check
         !------------------------------------------------------------------------------------------
         call pbuf_get_field(pbuf_chunk, prec_dp_idx,  prec_dp  )
         call pbuf_get_field(pbuf_chunk, snow_dp_idx,  snow_dp  )
         do i = 1,ncol
            icrm = ncol_sum + i
            qli_hydro_after(c,i) = 0.0_r8
            qi_hydro_after(c,i) = 0.0_r8
            do m = 1,crm_nz
               k = pver-m+1
               dp_g = state(c)%pdel(i,k)/gravit
               do jj = 1,crm_ny
                  do ii = 1,crm_nx
                     if(MMF_microphysics_scheme .eq. 'p3') then 
                        qli_hydro_after(c,i) = qli_hydro_after(c,i)+(crm_state%qr(icrm,ii,jj,m)) * dp_g
                        ! P3 does not have a separate category for snow, so qi_hydro_after should be zero
                     end if
                     if(MMF_microphysics_scheme .eq. 'sam1mom') then 
                        sfactor = max(0._r8,min(1._r8,(crm_state%temperature(icrm,ii,jj,m)-268.16)*1./(283.16-268.16)))
                        qli_hydro_after(c,i) = qli_hydro_after(c,i)+crm_state%qp(icrm,ii,jj,m) * dp_g
                        qi_hydro_after(c,i)  =  qi_hydro_after(c,i)+crm_state%qp(icrm,ii,jj,m) * (1-sfactor) * dp_g
                     end if ! MMF_microphysics_scheme
                  end do ! ii
               end do ! jj
            end do ! m = 1,crm_nz
            qli_hydro_after(c,i) = qli_hydro_after(c,i)/(crm_nx*crm_ny)
            qi_hydro_after(c,i)  =  qi_hydro_after(c,i)/(crm_nx*crm_ny)

            mmf_qchk_prec_dp(c,i) = prec_dp(i) + (qli_hydro_after(c,i) - qli_hydro_before(c,i))/ztodt/1000._r8
            mmf_qchk_snow_dp(c,i) = snow_dp(i) + ( qi_hydro_after(c,i) -  qi_hydro_before(c,i))/ztodt/1000._r8

         end do ! i = 1,ncol

         !------------------------------------------------------------------------------------------
         ! copy clear air relative humdity for aerosol water uptake
         !------------------------------------------------------------------------------------------
         call pbuf_get_field(pbuf_chunk, mmf_clear_rh_idx, mmf_clear_rh )
         ! initialize to zero, so no aerosol water uptake occurs by default
         mmf_clear_rh(1:ncol,1:pver) = 0
         do i = 1,ncol
            icrm = ncol_sum+i
            do m = 1,crm_nz
               k = pver-m+1
               mmf_clear_rh(i,k) = crm_clear_rh(icrm,m)
            end do ! m = 1,crm_nz
         end do ! i = 1,ncol

         !------------------------------------------------------------------------------------------
         !------------------------------------------------------------------------------------------
         ncol_sum = ncol_sum + ncol
      end do ! c=begchunk, endchunk

      deallocate(crm_clear_rh)

   end if ! (is_first_step())

   call t_stopf('crm')
   
   !------------------------------------------------------------------------------------------------
   ! Free memory in derived types
   !------------------------------------------------------------------------------------------------

   call crm_state_finalize(crm_state, MMF_microphysics_scheme)
   call crm_rad_finalize(crm_rad, MMF_microphysics_scheme)
   call crm_input_finalize(crm_input, MMF_microphysics_scheme)
   call crm_output_finalize(crm_output, MMF_microphysics_scheme)

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
      ptend%s(icol,pver)   = g_dp * cam_in%shf(icol)
      ptend%q(icol,pver,1) = g_dp * cam_in%cflx(icol,1)
   end do

end subroutine crm_surface_flux_bypass_tend

!==================================================================================================
!==================================================================================================

end module crm_physics

