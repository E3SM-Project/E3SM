module crm_physics
!---------------------------------------------------------------------------------------------------
! Purpose: Provides the interface to the crm code for the MMF configuration
!---------------------------------------------------------------------------------------------------
   use shr_kind_mod,    only: r8 => shr_kind_r8
   use shr_sys_mod,     only: shr_sys_flush   
   use cam_abortutils,  only: endrun
   use cam_logfile,     only: iulog
   use physics_types,   only: physics_state, physics_tend
   use ppgrid,          only: begchunk, endchunk, pcols, pver, pverp
   use constituents,    only: pcnst
#ifdef MODAL_AERO
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
   public :: m2005_effradius

   integer, public :: ncrms = -1 ! total number of CRMs summed over all chunks in task

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

contains
!===================================================================================================
!===================================================================================================

subroutine crm_physics_register()
!---------------------------------------------------------------------------------------------------
! Purpose:  add necessary fields into physics buffer
!---------------------------------------------------------------------------------------------------
   use spmd_utils,          only: masterproc
   use physconst,           only: cpair, mwh2o, mwdry
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
#elif defined(MMF_SAM) || defined(MMF_SAMOMP)
   use setparm_mod      ,   only: setparm
#endif
#ifdef MODAL_AERO
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
#ifdef MMF_SAMXX
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
   call setparm()

   if (use_MMF_VT) then
      ! add variance tracers
      call cnst_add('VT_T', real(0,r8), real(0,r8), real(0,r8), idx_vt_t, &
                    longname='VT_T', readiv=.false., mixtype='dry',cam_outfld=.false.)
      call cnst_add('VT_Q', real(0,r8), real(0,r8), real(0,r8), idx_vt_q, &
                    longname='VT_Q', readiv=.false., mixtype='dry',cam_outfld=.false.)
      call cnst_add('VT_U', real(0,r8), real(0,r8), real(0,r8), idx_vt_u, &
                    longname='VT_U', readiv=.false., mixtype='dry',cam_outfld=.false.)
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
   call pbuf_add_field('CRM_U',        'global',dtype_r8,dims_crm_3D,idx)
   call pbuf_add_field('CRM_V',        'global',dtype_r8,dims_crm_3D,idx)
   call pbuf_add_field('CRM_W',        'global',dtype_r8,dims_crm_3D,idx)
   call pbuf_add_field('CRM_T',        'global',dtype_r8,dims_crm_3D,idx)

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

   if (MMF_microphysics_scheme .eq. 'm2005') then
      call pbuf_add_field('CRM_NC_RAD','physpkg',dtype_r8,dims_crm_rad,crm_nc_rad_idx)
      call pbuf_add_field('CRM_NI_RAD','physpkg',dtype_r8,dims_crm_rad,crm_ni_rad_idx)
      call pbuf_add_field('CRM_QS_RAD','physpkg',dtype_r8,dims_crm_rad,crm_qs_rad_idx)
      call pbuf_add_field('CRM_NS_RAD','physpkg',dtype_r8,dims_crm_rad,crm_ns_rad_idx)

      call pbuf_add_field('CRM_QT',    'global', dtype_r8,dims_crm_3D,idx)
      call pbuf_add_field('CRM_NC',    'global', dtype_r8,dims_crm_3D,idx)
      call pbuf_add_field('CRM_QR',    'global', dtype_r8,dims_crm_3D,idx)
      call pbuf_add_field('CRM_NR',    'global', dtype_r8,dims_crm_3D,idx)
      call pbuf_add_field('CRM_QI',    'global', dtype_r8,dims_crm_3D,idx)
      call pbuf_add_field('CRM_NI',    'global', dtype_r8,dims_crm_3D,idx)
      call pbuf_add_field('CRM_QS',    'global', dtype_r8,dims_crm_3D,idx)
      call pbuf_add_field('CRM_NS',    'global', dtype_r8,dims_crm_3D,idx)
      call pbuf_add_field('CRM_QG',    'global', dtype_r8,dims_crm_3D,idx)
      call pbuf_add_field('CRM_NG',    'global', dtype_r8,dims_crm_3D,idx)
      call pbuf_add_field('CRM_QC',    'global', dtype_r8,dims_crm_3D,idx)

      if (prog_modal_aero) then
         call pbuf_add_field('RATE1_CW2PR_ST','physpkg',dtype_r8,dims_gcm_2D,idx)
      end if

   else
      call pbuf_add_field('CRM_QT',    'global', dtype_r8,dims_crm_3D,idx)
      call pbuf_add_field('CRM_QP',    'global', dtype_r8,dims_crm_3D,idx)
      call pbuf_add_field('CRM_QN',    'global', dtype_r8,dims_crm_3D,idx)
   end if

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
   use cam_history,    only: addfld
   use constituents,   only: apcnst, bpcnst, cnst_name, cnst_longname, cnst_get_ind
#ifdef ECPP
   use module_ecpp_ppdriver2, only: papampollu_init
#endif
   use constituents,          only: pcnst, cnst_get_ind
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
   integer :: ierror   ! Error code
   logical :: use_ECPP
   logical :: use_MMF_VT
   character(len=16) :: MMF_microphysics_scheme
   integer :: lchnk, ncol
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
      if ( any(mm == (/ ixcldliq, ixcldice, ixrain, ixsnow /)) ) then
         ! mass mixing ratios
         call addfld(cnst_name(mm), (/ 'lev' /), 'A', 'kg/kg', cnst_longname(mm)                   )
      else if ( any(mm == (/ ixnumliq, ixnumice, ixnumrain, ixnumsnow /)) ) then
         ! number concentrations
         call addfld(cnst_name(mm), (/ 'lev' /), 'A', '1/kg', cnst_longname(mm)                   )
      else
         call endrun( "crm_physics_init: Could not call addfld for constituent with unknown units.")
      endif
   end do

   call addfld(apcnst(ixcldliq), (/'lev'/), 'A', 'kg/kg', trim(cnst_name(ixcldliq))//' after physics'  )
   call addfld(apcnst(ixcldice), (/'lev'/), 'A', 'kg/kg', trim(cnst_name(ixcldice))//' after physics'  )
   call addfld(bpcnst(ixcldliq), (/'lev'/), 'A', 'kg/kg', trim(cnst_name(ixcldliq))//' before physics' )
   call addfld(bpcnst(ixcldice), (/'lev'/), 'A', 'kg/kg', trim(cnst_name(ixcldice))//' before physics' )
   call addfld(apcnst(ixrain),   (/'lev'/), 'A', 'kg/kg', trim(cnst_name(ixrain))//' after physics'  )
   call addfld(apcnst(ixsnow),   (/'lev'/), 'A', 'kg/kg', trim(cnst_name(ixsnow))//' after physics'  )
   call addfld(bpcnst(ixrain),   (/'lev'/), 'A', 'kg/kg', trim(cnst_name(ixrain))//' before physics' )
   call addfld(bpcnst(ixsnow),   (/'lev'/), 'A', 'kg/kg', trim(cnst_name(ixsnow))//' before physics' )

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
      if (MMF_microphysics_scheme .eq. 'm2005') then
         call pbuf_set_field(pbuf2d, crm_nc_rad_idx,0._r8)
         call pbuf_set_field(pbuf2d, crm_ni_rad_idx,0._r8)
         call pbuf_set_field(pbuf2d, crm_qs_rad_idx,0._r8)
         call pbuf_set_field(pbuf2d, crm_ns_rad_idx,0._r8)
      end if

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
   use crmdims,               only: crm_nx, crm_ny, crm_nz, crm_nx_rad, crm_ny_rad
   use physconst,             only: cpair, latvap, latice, gravit, cappa, pi
   use constituents,          only: pcnst, cnst_get_ind
#if defined(MMF_SAMXX)
   use cpp_interface_mod,     only: crm
#elif defined(MMF_SAM) || defined(MMF_SAMOMP)
   use crm_module,            only: crm
#endif
   use params_kind,           only: crm_rknd
   use phys_control,          only: phys_getopts, phys_do_flux_avg
   use crm_history,           only: crm_history_out
   use wv_saturation,         only: qsat_water
#if (defined  m2005 && defined MODAL_AERO)  
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
   use spmd_utils,            only: masterproc
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
   real(r8), pointer :: prec_dp(:)                 ! total precip from deep convection (ZM)    [m/s]
   real(r8), pointer :: snow_dp(:)                 ! snow from deep convection (ZM)            [m/s]
   real(r8), pointer :: ttend_dp(:,:)              ! Convective heating for gravity wave drag
   real(r8), pointer :: mmf_clear_rh(:,:)          ! clear air RH for aerosol water uptake
   real(r8), pointer :: cld(:,:)                   ! cloud fraction

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

   real(r8), dimension(begchunk:endchunk,pcols) :: qli_hydro_before ! column-integraetd rain + snow + graupel 
   real(r8), dimension(begchunk:endchunk,pcols) ::  qi_hydro_before ! column-integrated snow water + graupel water
   real(r8), dimension(begchunk:endchunk,pcols) :: qli_hydro_after  ! column-integraetd rain + snow + graupel 
   real(r8), dimension(begchunk:endchunk,pcols) ::  qi_hydro_after  ! column-integrated snow water + graupel water

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
   logical        :: use_MMF_VT_tmp                ! flag for MMF variance transport (for Fortran CRM)
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
   real(crm_rknd), pointer :: crm_qt(:,:,:,:) ! CRM total water

   real(crm_rknd), pointer :: crm_qp(:,:,:,:) ! 1-mom mass mixing ratio of precipitating condensate
   real(crm_rknd), pointer :: crm_qn(:,:,:,:) ! 1-mom mass mixing ratio of cloud condensate

   real(crm_rknd), pointer :: crm_nc(:,:,:,:) ! 2-mom mass mixing ratio of cloud water
   real(crm_rknd), pointer :: crm_qr(:,:,:,:) ! 2-mom number concentration of cloud water
   real(crm_rknd), pointer :: crm_nr(:,:,:,:) ! 2-mom mass mixing ratio of rain
   real(crm_rknd), pointer :: crm_qi(:,:,:,:) ! 2-mom number concentration of rain
   real(crm_rknd), pointer :: crm_ni(:,:,:,:) ! 2-mom mass mixing ratio of cloud ice
   real(crm_rknd), pointer :: crm_qs(:,:,:,:) ! 2-mom number concentration of cloud ice
   real(crm_rknd), pointer :: crm_ns(:,:,:,:) ! 2-mom mass mixing ratio of snow
   real(crm_rknd), pointer :: crm_qg(:,:,:,:) ! 2-mom number concentration of snow
   real(crm_rknd), pointer :: crm_ng(:,:,:,:) ! 2-mom mass mixing ratio of graupel
   real(crm_rknd), pointer :: crm_qc(:,:,:,:) ! 2-mom number concentration of graupel

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

   ! CRM mean state acceleration (MSA) parameters
   use_crm_accel = .false.
   crm_accel_factor = 0.
   crm_accel_uv = .false.
   call phys_getopts(use_crm_accel_out    = use_crm_accel_tmp)
   call phys_getopts(crm_accel_factor_out = crm_accel_factor)
   call phys_getopts(crm_accel_uv_out     = crm_accel_uv_tmp)
   use_crm_accel = use_crm_accel_tmp
   crm_accel_uv = crm_accel_uv_tmp

   if (masterproc) then
     if (use_crm_accel .and. trim(MMF_microphysics_scheme)/='sam1mom') then
       write(0,*) "CRM time step relaxation is only compatible with sam1mom microphysics"
       call endrun('crm main')
     endif
   endif

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
         call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QT'), crm_qt)
         if (MMF_microphysics_scheme .eq. 'sam1mom') then
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QP'), crm_qp)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QN'), crm_qn)
         else
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_NC'), crm_nc)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QR'), crm_qr)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_NR'), crm_nr)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QI'), crm_qi)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_NI'), crm_ni)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QS'), crm_qs)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_NS'), crm_ns)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QG'), crm_qg)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_NG'), crm_ng)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QC'), crm_qc)
         end if

         ! initialize all of total water to zero (needed for ncol < i <= pcols)
         crm_qt(1:pcols,:,:,:) = 0.0_r8

         do i = 1,ncol
            do k = 1,crm_nz
               m = pver-k+1

               crm_u(i,:,:,k) = state(c)%u(i,m) * cos( crm_angle(i) ) + state(c)%v(i,m) * sin( crm_angle(i) )
               crm_v(i,:,:,k) = state(c)%v(i,m) * cos( crm_angle(i) ) - state(c)%u(i,m) * sin( crm_angle(i) )
               crm_w(i,:,:,k) = 0.
               crm_t(i,:,:,k) = state(c)%t(i,m)

               ! Initialize microphysics arrays
               if (MMF_microphysics_scheme .eq. 'sam1mom') then
                  crm_qt(i,:,:,k) = state(c)%q(i,m,1)+state(c)%q(i,m,ixcldliq)+state(c)%q(i,m,ixcldice)
                  crm_qp(i,:,:,k) = 0.0_r8
                  crm_qn(i,:,:,k) = state(c)%q(i,m,ixcldliq)+state(c)%q(i,m,ixcldice)
               else if (MMF_microphysics_scheme .eq. 'm2005') then
                  crm_qt(i,:,:,k) = state(c)%q(i,m,1)+state(c)%q(i,m,ixcldliq)
                  crm_qc(i,:,:,k) = state(c)%q(i,m,ixcldliq)
                  crm_qi(i,:,:,k) = state(c)%q(i,m,ixcldice)
                  crm_nc(i,:,:,k) = 0.0_r8
                  crm_qr(i,:,:,k) = 0.0_r8
                  crm_nr(i,:,:,k) = 0.0_r8
                  crm_ni(i,:,:,k) = 0.0_r8
                  crm_qs(i,:,:,k) = 0.0_r8
                  crm_ns(i,:,:,k) = 0.0_r8
                  crm_qg(i,:,:,k) = 0.0_r8
                  crm_ng(i,:,:,k) = 0.0_r8
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
         if (MMF_microphysics_scheme .eq. 'm2005') then
            call pbuf_get_field(pbuf_chunk, crm_nc_rad_idx,crm_nc_rad)
            call pbuf_get_field(pbuf_chunk, crm_ni_rad_idx,crm_ni_rad)
            call pbuf_get_field(pbuf_chunk, crm_qs_rad_idx,crm_qs_rad)
            call pbuf_get_field(pbuf_chunk, crm_ns_rad_idx,crm_ns_rad)
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
               if (MMF_microphysics_scheme .eq. 'm2005') then
                  crm_nc_rad(i,:,:,k) = 0.0
                  crm_ni_rad(i,:,:,k) = 0.0
                  crm_qs_rad(i,:,:,k) = 0.0
                  crm_ns_rad(i,:,:,k) = 0.0
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
         call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QT'), crm_qt)
         if (MMF_microphysics_scheme .eq. 'sam1mom') then
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QP'), crm_qp)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QN'), crm_qn)
         else
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_NC'), crm_nc)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QR'), crm_qr)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_NR'), crm_nr)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QI'), crm_qi)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_NI'), crm_ni)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QS'), crm_qs)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_NS'), crm_ns)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QG'), crm_qg)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_NG'), crm_ng)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QC'), crm_qc)
         end if

         ! copy pbuf data into crm_state
         do i = 1,ncol
            icrm = ncol_sum + i
            crm_state%u_wind     (icrm,:,:,:) = crm_u (i,:,:,:)
            crm_state%v_wind     (icrm,:,:,:) = crm_v (i,:,:,:)
            crm_state%w_wind     (icrm,:,:,:) = crm_w (i,:,:,:)
            crm_state%temperature(icrm,:,:,:) = crm_t (i,:,:,:)
            crm_state%qt         (icrm,:,:,:) = crm_qt(i,:,:,:)
            if (MMF_microphysics_scheme .eq. 'sam1mom') then
               crm_state%qp      (icrm,:,:,:) = crm_qp(i,:,:,:)
               crm_state%qn      (icrm,:,:,:) = crm_qn(i,:,:,:)
            else
               crm_state%nc      (icrm,:,:,:) = crm_nc(i,:,:,:)
               crm_state%qr      (icrm,:,:,:) = crm_qr(i,:,:,:)
               crm_state%nr      (icrm,:,:,:) = crm_nr(i,:,:,:)
               crm_state%qi      (icrm,:,:,:) = crm_qi(i,:,:,:)
               crm_state%ni      (icrm,:,:,:) = crm_ni(i,:,:,:)
               crm_state%qs      (icrm,:,:,:) = crm_qs(i,:,:,:)
               crm_state%ns      (icrm,:,:,:) = crm_ns(i,:,:,:)
               crm_state%qg      (icrm,:,:,:) = crm_qg(i,:,:,:)
               crm_state%ng      (icrm,:,:,:) = crm_ng(i,:,:,:)
               crm_state%qc      (icrm,:,:,:) = crm_qc(i,:,:,:)
            end if
         end do ! i=1,ncol

         ! Retrieve radiative heating tendency
         call pbuf_get_field(pbuf_chunk, crm_qrad_idx, crm_qrad)
         do m = 1,crm_nz
            k = pver-m+1
            do i = 1,ncol
               icrm = ncol_sum + i
               crm_rad%qrad(icrm,:,:,m) = crm_qrad(i,:,:,m) / state(c)%pdel(i,k) ! normalize for energy conservation
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
                     if (MMF_microphysics_scheme .eq. 'm2005') then
                        qli_hydro_before(c,i) = qli_hydro_before(c,i)+(crm_state%qr(icrm,ii,jj,m)+ &
                                                                       crm_state%qs(icrm,ii,jj,m)+ &
                                                                       crm_state%qg(icrm,ii,jj,m)) * dp_g
                        qi_hydro_before(c,i)  =  qi_hydro_before(c,i)+(crm_state%qs(icrm,ii,jj,m)+ &
                                                                       crm_state%qg(icrm,ii,jj,m)) * dp_g
                     else if (MMF_microphysics_scheme .eq. 'sam1mom') then
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
#if defined( MMF_ESMT )
            ! Set the input wind for ESMT
            crm_input%ul_esmt(icrm,1:pver) = state(c)%u(i,1:pver)
            crm_input%vl_esmt(icrm,1:pver) = state(c)%v(i,1:pver)
#endif /* MMF_ESMT */
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
         ! Set surface flux variables
         !------------------------------------------------------------------------------------------
         if (phys_do_flux_avg()) then
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('LHFLX'), shf_ptr)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('SHFLX'), lhf_ptr)
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
            crm_input%fluxu00(icrm) = wsx_tmp(i)         ! N/m2
            crm_input%fluxv00(icrm) = wsy_tmp(i)         ! N/m2
            crm_input%fluxt00(icrm) = shf_tmp(i)/cpair   ! K Kg/ (m2 s)
            crm_input%fluxq00(icrm) = lhf_tmp(i)/latvap  ! Kg/(m2 s)
            crm_input%wndls(icrm)   = sqrt(state(c)%u(i,pver)**2 + state(c)%v(i,pver)**2)
         end do

         !------------------------------------------------------------------------------------------
         ! Set aerosol
         !------------------------------------------------------------------------------------------
#if (defined m2005 && defined MODAL_AERO)
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

      call t_startf ('crm_call')

      ! Fortran classes don't translate to C++ classes, we we have to separate
      ! this stuff out when calling the C++ routinte crm(...)
      call crm(ncrms, ncrms, ztodt, pver, crm_input%bflxls, crm_input%wndls, crm_input%zmid, crm_input%zint, &
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
               crm_input%u_vt, crm_output%u_vt_tend, crm_output%u_vt_ls, &
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
               latitude0, longitude0, gcolp, nstep, &
               use_MMF_VT, MMF_VT_wn_max, &
               use_crm_accel, crm_accel_factor, crm_accel_uv)

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
         ! Populate output tendencies for 2-mom microphysics
         !------------------------------------------------------------------------------------------

         if (MMF_microphysics_scheme .eq. 'm2005') then
            ptend(c)%lq(ixnumliq)  = .TRUE.
            ptend(c)%lq(ixnumice)  = .TRUE.
            if (use_ECPP) then
               ptend(c)%lq(ixrain)    = .TRUE. 
               ptend(c)%lq(ixsnow)    = .TRUE. 
               ptend(c)%lq(ixnumrain) = .TRUE. 
               ptend(c)%lq(ixnumsnow) = .TRUE. 
            end if

            do i = 1, ncol
               icrm = ncol_sum + i
               do k = 1, crm_nz 
                  m = pver-k+1
                  do ii = 1, crm_nx
                  do jj = 1, crm_ny
                     ptend(c)%q(i,m,ixnumliq)  = ptend(c)%q(i,m,ixnumliq)  + crm_state%nc(icrm,ii,jj,k) 
                     ptend(c)%q(i,m,ixnumice)  = ptend(c)%q(i,m,ixnumice)  + crm_state%ni(icrm,ii,jj,k)
                     if (use_ECPP) then
                        ptend(c)%q(i,m,ixrain)    = ptend(c)%q(i,m,ixrain)    + crm_state%qr(icrm,ii,jj,k)
                        ptend(c)%q(i,m,ixsnow)    = ptend(c)%q(i,m,ixsnow)    + crm_state%qs(icrm,ii,jj,k)
                        ptend(c)%q(i,m,ixnumrain) = ptend(c)%q(i,m,ixnumrain) + crm_state%nr(icrm,ii,jj,k)
                        ptend(c)%q(i,m,ixnumsnow) = ptend(c)%q(i,m,ixnumsnow) + crm_state%ns(icrm,ii,jj,k)
                     end if
                  end do
                  end do
                  ptend(c)%q(i,m,ixnumliq)  = (ptend(c)%q(i,m,ixnumliq) /(crm_nx*crm_ny) - state(c)%q(i,m,ixnumliq)) /ztodt
                  ptend(c)%q(i,m,ixnumice)  = (ptend(c)%q(i,m,ixnumice) /(crm_nx*crm_ny) - state(c)%q(i,m,ixnumice)) /ztodt
                  if (use_ECPP) then
                     ptend(c)%q(i,m,ixrain)    = (ptend(c)%q(i,m,ixrain)   /(crm_nx*crm_ny) - state(c)%q(i,m,ixrain))   /ztodt
                     ptend(c)%q(i,m,ixsnow)    = (ptend(c)%q(i,m,ixsnow)   /(crm_nx*crm_ny) - state(c)%q(i,m,ixsnow))   /ztodt
                     ptend(c)%q(i,m,ixnumrain) = (ptend(c)%q(i,m,ixnumrain)/(crm_nx*crm_ny) - state(c)%q(i,m,ixnumrain))/ztodt
                     ptend(c)%q(i,m,ixnumsnow) = (ptend(c)%q(i,m,ixnumsnow)/(crm_nx*crm_ny) - state(c)%q(i,m,ixnumsnow))/ztodt
                  end if
               end do
            end do
         end if

         !------------------------------------------------------------------------------------------
         ! set convective heating tendency for gravity wave drag
         !------------------------------------------------------------------------------------------
         if (ttend_dp_idx > 0) then
            call pbuf_get_field(pbuf_chunk, ttend_dp_idx, ttend_dp)
            ttend_dp(1:ncol,1:pver) = ptend(c)%s(1:ncol,1:pver)/cpair
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

         ! The radiation tendencies in the GCM levels above the CRM and the top 2 CRM levels are set to
         ! be zero in the CRM, So add radiation tendencies to these levels 
         ptend(c)%s(1:ncol, 1:pver-crm_nz+2) = qrs(1:ncol,1:pver-crm_nz+2) + qrl(1:ncol,1:pver-crm_nz+2)

         ! This will be used to check energy conservation
         mmf_rad_flux(c,:ncol) = 0.0_r8
         do k = 1,pver
            do i = 1,ncol
               mmf_rad_flux(c,i) = mmf_rad_flux(c,i) + ( qrs(i,k) + qrl(i,k) ) * state(c)%pdel(i,k)/gravit
            end do
         end do

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

#if defined( MMF_USE_ESMT )
         ptend(c)%lu = .TRUE.
         ptend(c)%lv = .TRUE.
         do i = 1,ncol
            icrm = ncol_sum + i
            ptend(c)%u(i,1:pver)  = crm_output%u_tend_esmt(icrm,1:pver)
            ptend(c)%v(i,1:pver)  = crm_output%v_tend_esmt(icrm,1:pver)
         end do
#else /* MMF_USE_ESMT not defined */  

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

#endif /* MMF_USE_ESMT */

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

         if (MMF_microphysics_scheme .eq. 'm2005') then
            call pbuf_get_field(pbuf_chunk, crm_nc_rad_idx, crm_nc_rad)
            call pbuf_get_field(pbuf_chunk, crm_ni_rad_idx, crm_ni_rad)
            call pbuf_get_field(pbuf_chunk, crm_qs_rad_idx, crm_qs_rad)
            call pbuf_get_field(pbuf_chunk, crm_ns_rad_idx, crm_ns_rad)
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
            if (MMF_microphysics_scheme .eq. 'm2005') then
               crm_nc_rad(i,:,:,:) = crm_rad%nc         (icrm,:,:,:)
               crm_ni_rad(i,:,:,:) = crm_rad%ni         (icrm,:,:,:)
               crm_qs_rad(i,:,:,:) = crm_rad%qs         (icrm,:,:,:)
               crm_ns_rad(i,:,:,:) = crm_rad%ns         (icrm,:,:,:)
            end if
         end do

         !------------------------------------------------------------------------------------------
         ! put CRM state data back in pbuf
         !------------------------------------------------------------------------------------------
         call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_U'),  crm_u)
         call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_V'),  crm_v)
         call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_W'),  crm_w)
         call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_T'),  crm_t)
         call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QT'), crm_qt)
         if (MMF_microphysics_scheme .eq. 'sam1mom') then
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QP'), crm_qp)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QN'), crm_qn)
         else
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_NC'), crm_nc)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QR'), crm_qr)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_NR'), crm_nr)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QI'), crm_qi)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_NI'), crm_ni)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QS'), crm_qs)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_NS'), crm_ns)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QG'), crm_qg)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_NG'), crm_ng)
            call pbuf_get_field(pbuf_chunk, pbuf_get_index('CRM_QC'), crm_qc)
         end if

         do i = 1,ncol
            icrm = ncol_sum + i
            crm_u (i,:,:,:) = crm_state%u_wind     (icrm,:,:,:)
            crm_v (i,:,:,:) = crm_state%v_wind     (icrm,:,:,:)
            crm_w (i,:,:,:) = crm_state%w_wind     (icrm,:,:,:)
            crm_t (i,:,:,:) = crm_state%temperature(icrm,:,:,:)
            crm_qt(i,:,:,:) = crm_state%qt         (icrm,:,:,:)
            if (MMF_microphysics_scheme .eq. 'sam1mom') then
               crm_qp(i,:,:,:) = crm_state%qp(icrm,:,:,:)
               crm_qn(i,:,:,:) = crm_state%qn(icrm,:,:,:)
            else if (MMF_microphysics_scheme .eq. 'm2005') then
               crm_qc(i,:,:,:) = crm_state%qc(icrm,:,:,:)
               crm_qi(i,:,:,:) = crm_state%qi(icrm,:,:,:)
               crm_nc(i,:,:,:) = crm_state%nc(icrm,:,:,:)
               crm_qr(i,:,:,:) = crm_state%qr(icrm,:,:,:)
               crm_nr(i,:,:,:) = crm_state%nr(icrm,:,:,:)
               crm_ni(i,:,:,:) = crm_state%ni(icrm,:,:,:)
               crm_qs(i,:,:,:) = crm_state%qs(icrm,:,:,:)
               crm_ns(i,:,:,:) = crm_state%ns(icrm,:,:,:)
               crm_qg(i,:,:,:) = crm_state%qg(icrm,:,:,:)
               crm_ng(i,:,:,:) = crm_state%ng(icrm,:,:,:)
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
                     if(MMF_microphysics_scheme .eq. 'm2005') then
                        qli_hydro_after(c,i) = qli_hydro_after(c,i)+(crm_state%qr(icrm,ii,jj,m)+ &
                                                                     crm_state%qs(icrm,ii,jj,m)+ &
                                                                     crm_state%qg(icrm,ii,jj,m)) * dp_g
                        qi_hydro_after(c,i)  =  qi_hydro_after(c,i)+(crm_state%qs(icrm,ii,jj,m)+ &
                                                                     crm_state%qg(icrm,ii,jj,m)) * dp_g
                     else if(MMF_microphysics_scheme .eq. 'sam1mom') then 
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
