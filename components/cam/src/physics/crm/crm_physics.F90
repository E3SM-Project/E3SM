#if defined( SP_ORIENT_RAND ) && defined( SP_DIR_NS )
#undef SP_DIR_NS
#endif

module crm_physics

!---------------------------------------------------------------------------------------------------
! Purpose: 
! 
!    Provides the CAM interface to the crm code.  
!
! Revision history: 
! June, 2009, Minghuai Wang:  
!          crm_physics_tend 
! July, 2009, Minghuai Wang: m2005_effradius
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
   integer :: prec_dp_idx, snow_dp_idx, prec_sh_idx, snow_sh_idx
   integer :: prec_sed_idx, snow_sed_idx, snow_str_idx, prec_pcw_idx, snow_pcw_idx
   integer :: cldo_idx
#ifdef MAML
   integer :: crm_pcp_idx,crm_snw_idx
#endif
   real(r8),pointer                        :: acldy_cen_tbeg(:,:)        ! cloud fraction
   real(r8), pointer, dimension(:,:)       :: cldo

   type(physics_state)                     :: state_save
   type(physics_tend)                      :: tend_save
   real(r8), dimension(pcols, pver)        :: cldo_save     ! saved cloud fraction
   real(r8), dimension(pcols, pver, pcnst) :: q_aero        ! used to keep aerosol changes from offline physics
#ifdef MODAL_AERO
   real(r8)          :: qqcw_save(pcols,pver,pcnst)
   real(r8)          :: qqcw_all(pcols,pver,pcnst)
   real(r8)          :: dgnumwet_save(pcols, pver, ntot_amode)
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
  use spmd_utils,      only: masterproc
  use physconst,       only: mwdry, cpair
  use ppgrid,          only: pcols, pver, pverp
  use physics_buffer,  only: dyn_time_lvls, pbuf_add_field, dtype_r8, pbuf_get_index
  use phys_control,    only: phys_getopts
  use crmdims,         only: crm_nx, crm_ny, crm_nz, crm_dx, crm_dy, crm_dt, crm_nx_rad, crm_ny_rad
  use setparm_mod,         only: setparm
  use cam_history_support, only: add_hist_coord

  ! local variables
  integer idx
  logical           :: use_ECPP, use_SPCAM
  character(len=16) :: SPCAM_microp_scheme

  call phys_getopts( use_ECPP_out            = use_ECPP)
  call phys_getopts( SPCAM_microp_scheme_out = SPCAM_microp_scheme)
  call phys_getopts( use_SPCAM_out           = use_SPCAM)

  if(masterproc) then
      print*,'_________________________________________'
      print*,'_ Super-parameterization run ____________'
      print*,'crm_nx=',crm_nx,'   crm_ny=',crm_ny,'   crm_nz=',crm_nz
      print*,'crm_dx=',crm_dx,'   crm_dy=',crm_dy,'   crm_dt=',crm_dt
      if (SPCAM_microp_scheme .eq. 'sam1mom') print*,'Microphysics: SAM1MOM'
      if (SPCAM_microp_scheme .eq. 'm2005') print*,'Microphysics: M2005'
      print*,'_________________________________________'
  end if

  if (use_SPCAM) then
     call setparm()
  end if

  ! Add crm dimensions to cam history 
  call add_hist_coord('crm_nx'       ,crm_nx,  'CRM NX')
  call add_hist_coord('crm_ny'       ,crm_ny,  'CRM NY')
  call add_hist_coord('crm_nz'       ,crm_nz,  'CRM NZ')
  call add_hist_coord('crm_nx_rad', crm_nx_rad, 'Number of x columns for radiation')
  call add_hist_coord('crm_ny_rad', crm_ny_rad, 'Number of y columns for radiation')

  call pbuf_add_field('CRM_U',     'global', dtype_r8, (/pcols,crm_nx,crm_ny,crm_nz/), idx)
  call pbuf_add_field('CRM_V',     'global', dtype_r8, (/pcols,crm_nx,crm_ny,crm_nz/), idx)
  call pbuf_add_field('CRM_W',     'global', dtype_r8, (/pcols,crm_nx,crm_ny,crm_nz/), idx)
  call pbuf_add_field('CRM_T',     'global', dtype_r8, (/pcols,crm_nx,crm_ny,crm_nz/), idx)

  call pbuf_add_field('CRM_T_RAD',   'physpkg', dtype_r8, (/pcols,crm_nx_rad,crm_ny_rad,crm_nz/), idx)
  call pbuf_add_field('CRM_QV_RAD',  'physpkg', dtype_r8, (/pcols,crm_nx_rad,crm_ny_rad,crm_nz/), idx)
  call pbuf_add_field('CRM_QC_RAD',  'physpkg', dtype_r8, (/pcols,crm_nx_rad,crm_ny_rad,crm_nz/), idx)
  call pbuf_add_field('CRM_QI_RAD',  'physpkg', dtype_r8, (/pcols,crm_nx_rad,crm_ny_rad,crm_nz/), idx)
  call pbuf_add_field('CRM_CLD_RAD', 'physpkg', dtype_r8, (/pcols,crm_nx_rad,crm_ny_rad,crm_nz/), idx)
  call pbuf_add_field('CRM_QRAD',    'global',  dtype_r8, (/pcols,crm_nx_rad,crm_ny_rad,crm_nz/), idx)

#ifdef MODAL_AERO
  call pbuf_add_field('CRM_QAERWAT', 'physpkg', dtype_r8, (/pcols,crm_nx_rad,crm_ny_rad,crm_nz,ntot_amode/),  crm_qaerwat_idx)
  call pbuf_add_field('CRM_DGNUMWET','physpkg', dtype_r8, (/pcols,crm_nx_rad,crm_ny_rad,crm_nz,ntot_amode/),  crm_dgnumwet_idx)
#endif
   
   cldo_idx = pbuf_get_index('CLDO')
   call pbuf_add_field('CLDO',      'global',  dtype_r8, (/pcols ,pver, dyn_time_lvls/),                  cldo_idx  )

  if (SPCAM_microp_scheme .eq. 'm2005') then
    call pbuf_add_field('CRM_NC_RAD','physpkg', dtype_r8, (/pcols, crm_nx_rad, crm_ny_rad, crm_nz/), idx)
    call pbuf_add_field('CRM_NI_RAD','physpkg', dtype_r8, (/pcols, crm_nx_rad, crm_ny_rad, crm_nz/), idx)
    call pbuf_add_field('CRM_QS_RAD','physpkg', dtype_r8, (/pcols, crm_nx_rad, crm_ny_rad, crm_nz/), idx)
    call pbuf_add_field('CRM_NS_RAD','physpkg', dtype_r8, (/pcols, crm_nx_rad, crm_ny_rad, crm_nz/), idx)

    call pbuf_add_field('CRM_QT',    'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/), idx)
    call pbuf_add_field('CRM_NC',    'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/), idx)
    call pbuf_add_field('CRM_QR',    'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/), idx)
    call pbuf_add_field('CRM_NR',    'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/), idx)
    call pbuf_add_field('CRM_QI',    'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/), idx)
    call pbuf_add_field('CRM_NI',    'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/), idx)
    call pbuf_add_field('CRM_QS',    'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/), idx)
    call pbuf_add_field('CRM_NS',    'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/), idx)
    call pbuf_add_field('CRM_QG',    'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/), idx)
    call pbuf_add_field('CRM_NG',    'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/), idx)
    call pbuf_add_field('CRM_QC',    'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/), idx)
  else
    call pbuf_add_field('CRM_QT',    'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/), idx)
    call pbuf_add_field('CRM_QP',    'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/), idx)
    call pbuf_add_field('CRM_QN',    'global',  dtype_r8, (/pcols, crm_nx, crm_ny, crm_nz/), idx)
  endif
   
  call pbuf_add_field('MU_CRM',    'physpkg', dtype_r8, (/pcols,pver/), idx)  ! mass flux up
  call pbuf_add_field('MD_CRM',    'physpkg', dtype_r8, (/pcols,pver/), idx)  ! mass flux down
  call pbuf_add_field('DU_CRM',    'physpkg', dtype_r8, (/pcols,pver/), idx)  ! mass detrainment from updraft
  call pbuf_add_field('EU_CRM',    'physpkg', dtype_r8, (/pcols,pver/), idx)  ! mass detrainment from updraft
  call pbuf_add_field('ED_CRM',    'physpkg', dtype_r8, (/pcols,pver/), idx)  ! mass detrainment from downdraft
  call pbuf_add_field('JT_CRM',    'physpkg', dtype_r8, (/pcols/),      idx)  ! index of cloud (convection) top for each column
  call pbuf_add_field('MX_CRM',    'physpkg', dtype_r8, (/pcols/),      idx)  ! index of cloud (convection) bottom for each column
  call pbuf_add_field('IDEEP_CRM', 'physpkg', dtype_r8, (/pcols/),      idx)  ! Gathering array for convective columns

  call pbuf_add_field('TKE_CRM',   'physpkg', dtype_r8, (/pcols,pver/), idx) ! TKE from CRM  (m2/s2)
  call pbuf_add_field('TK_CRM',    'physpkg', dtype_r8, (/pcols,pver/), idx) ! TK from CRM (m2/s)

  ! ACLDY_CEN has to be global in the physcal buffer to be saved in the restart file??
  call pbuf_add_field('ACLDY_CEN','global', dtype_r8, (/pcols,pver/), idx) ! total (all sub-classes) cloudy fractional area in previous time step 
#ifdef MAML
  !adding new variables for passing CRM-scale precipition/snow into CLM. Added new variables
  call pbuf_add_field('CRM_PCP',     'physpkg', dtype_r8, (/pcols,crm_nx, crm_ny/),                crm_pcp_idx)
  call pbuf_add_field('CRM_SNW',     'physpkg', dtype_r8, (/pcols,crm_nx, crm_ny/),                crm_snw_idx)
#endif

end subroutine crm_physics_register

!===================================================================================================
!===================================================================================================

subroutine crm_physics_init(species_class)
!---------------------------------------------------------------------------------------------------
! 
! Purpose: initialize some variables, and add necessary fields into output fields 
!
!---------------------------------------------------------------------------------------------------
   use physics_buffer,  only: pbuf_get_index
   use physconst,       only: mwdry, cpair, spec_class_gas
   use ppgrid,          only: pcols, pver, pverp
   use constituents,    only: pcnst, cnst_name
   use cam_history,     only: addfld, add_default, horiz_only
   use cam_history_support, only: add_hist_coord
   use phys_control,    only: phys_getopts
   use crmdims,         only: crm_nx, crm_ny, crm_nz, crm_nx_rad, crm_ny_rad
#ifdef ECPP
   use module_ecpp_ppdriver2,  only: papampollu_init
   use ecppvars,               only: NCLASS_CL,ncls_ecpp_in,NCLASS_PR
#endif
   use accelerate_crm_mod, only: crm_accel_init
#ifdef MODAL_AERO
   use cam_history,     only: fieldname_len
   use modal_aero_data, only: cnst_name_cw, &
                              lmassptr_amode, lmassptrcw_amode, &
                              nspec_amode, ntot_amode, numptr_amode, numptrcw_amode
                              
   integer :: l, lphase, lspec
   character(len=fieldname_len)   :: tmpname
   character(len=fieldname_len+3) :: fieldname
   character(128)                 :: long_name
   character(8)                   :: unit
#endif
    
   ! species_class is defined as input so it needs to be outside of MODAL_AERO condition for 1-moment micro to work
   integer, intent(in) :: species_class(:)

   ! local variables
   integer :: m
   logical :: use_ECPP
   character(len=16) :: SPCAM_microp_scheme

   call phys_getopts(use_ECPP_out = use_ECPP)
   call phys_getopts(SPCAM_microp_scheme_out = SPCAM_microp_scheme)

#ifdef ECPP
   if (use_ECPP) then
      ! Initialize ECPP driver
      call papampollu_init()
   end if
#endif

   if (SPCAM_microp_scheme .eq. 'm2005') then
      call addfld ('SPNC    ',(/ 'lev' /), 'A', '/kg   ','Cloud water dropet number from CRM')
      call addfld ('SPNI    ',(/ 'lev' /), 'A', '/kg   ','Cloud ice crystal number from CRM')
      call addfld ('SPNS    ',(/ 'lev' /), 'A', '/kg   ','Snow particle number from CRM')
      call addfld ('SPNG    ',(/ 'lev' /), 'A', '/kg   ','Graupel particle number from CRM')
      call addfld ('SPNR    ',(/ 'lev' /), 'A', '/kg   ','Rain particle number from CRM')
      call add_default ('SPNC    ', 1, ' ')
      call add_default ('SPNI    ', 1, ' ')
      call add_default ('SPNS    ', 1, ' ')
      call add_default ('SPNG    ', 1, ' ')
      call add_default ('SPNR    ', 1, ' ')
  
      call addfld ('CRM_FLIQ ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '1      ','Frequency of Occurrence of Liquid'      )
      call addfld ('CRM_FICE ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '1      ','Frequency of Occurrence of Ice'         )
      call addfld ('CRM_FRAIN',(/'crm_nx','crm_ny','crm_nz'/), 'A', '1      ','Frequency of Occurrence of Rain'        )
      call addfld ('CRM_FSNOW',(/'crm_nx','crm_ny','crm_nz'/), 'A', '1      ','Frequency of Occurrence of Snow'        )
      call addfld ('CRM_FGRAP',(/'crm_nx','crm_ny','crm_nz'/), 'A', '1      ','Frequency of Occurrence of Graupel'     )
  
      call addfld ('CRM_QS  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', 'kg/kg   ','Snow mixing ratio from CRM'             )
      call addfld ('CRM_QG  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', 'kg/kg   ','Graupel mixing ratio from CRM'          )
      call addfld ('CRM_QR  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', 'kg/kg   ','Rain mixing ratio from CRM'             )
  
      call addfld ('CRM_NC  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/kg     ','Cloud water dropet number from CRM'     )
      call addfld ('CRM_NI  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/kg     ','Cloud ice crystal number from CRM'      )
      call addfld ('CRM_NS  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/kg     ','Snow particle number from CRM'          )
      call addfld ('CRM_NG  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/kg     ','Graupel particle number from CRM'       )
      call addfld ('CRM_NR  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/kg     ','Rain particle number from CRM'          )
  
      ! hm 7/26/11, add new output
      ! below is for *instantaneous* crm output
      call addfld ('CRM_AUT  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/s     ','Autoconversion cloud waterfrom CRM'     )
      call addfld ('CRM_ACC  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/s     ','Accretion cloud water from CRM'         )
      call addfld ('CRM_EVPC ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/s     ','Evaporation cloud water from CRM'       )
      call addfld ('CRM_EVPR ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/s     ','Evaporation rain from CRM'              )
      call addfld ('CRM_MLT  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/s     ','Melting ice snow graupel from CRM'      )
      call addfld ('CRM_SUB  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/s     ','Sublimation ice snow graupel from CRM'  )
      call addfld ('CRM_DEP  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/s     ','Deposition ice snow graupel from CRM'   )
      call addfld ('CRM_CON  ',(/'crm_nx','crm_ny','crm_nz'/), 'A', '/s     ','Condensation cloud water from CRM'      )
  
      ! hm 8/31/11 -  *gcm-grid and time-step-avg* process output
      call addfld ('A_AUT  ',(/ 'lev' /), 'A', '/s   ','Avg autoconversion cloud water from CRM'            )
      call addfld ('A_ACC  ',(/ 'lev' /), 'A', '/s   ','Avg accretion cloud water from CRM'                 )
      call addfld ('A_EVPC ',(/ 'lev' /), 'A', '/s   ','Avg evaporation cloud water from CRM'               )
      call addfld ('A_EVPR ',(/ 'lev' /), 'A', '/s   ','Avg evaporation rain from CRM'                      )
      call addfld ('A_MLT  ',(/ 'lev' /), 'A', '/s   ','Avg melting ice snow graupel from CRM'              )
      call addfld ('A_SUB  ',(/ 'lev' /), 'A', '/s   ','Avg sublimation ice snow graupel from CRM'          )
      call addfld ('A_DEP  ',(/ 'lev' /), 'A', '/s   ','Avg deposition ice snow graupel from CRM'           )
      call addfld ('A_CON  ',(/ 'lev' /), 'A', '/s   ','Avg condensation cloud water from CRM'              )
  
      call addfld ('CRM_DES  ', (/'crm_nx','crm_ny','crm_nz'/), 'A', 'micrometers', 'cloud scale snow effective diameter')
      call addfld ('CRM_MU   ', (/'crm_nx','crm_ny','crm_nz'/), 'A', 'micrometers', &
                                                 'cloud scale droplet size distribution shape parameter for radiation')
      call addfld ('CRM_LAMBDA',(/'crm_nx','crm_ny','crm_nz'/), 'A', 'micrometers',  &
                                                 'cloud scale slope of droplet distribution for radiation')
      call addfld ('CRM_TAU  ', (/'crm_nx','crm_ny','crm_nz'/), 'A', '1',           'cloud scale cloud optical depth'  )
      call addfld ('CRM_WVAR' , (/'crm_nx','crm_ny','crm_nz'/), 'A', 'm/s',         'vertical velocity variance from CRM')
   end if  ! SPCAM_microp_scheme .eq. 'm2005'

   call addfld ('CRM_DEI  ', (/'crm_nx','crm_ny','crm_nz'/), 'A', 'micrometers', 'cloud scale Mitchell ice effective diameter')
   call addfld ('CRM_REL  ', (/'crm_nx','crm_ny','crm_nz'/), 'A', 'micrometers', 'cloud scale droplet effective radius')
   call addfld ('CRM_REI  ', (/'crm_nx','crm_ny','crm_nz'/), 'A', 'micrometers', 'cloud scale ice crystal effective radius')
 
   call addfld ('CRM_FSNT',  (/'crm_nx','crm_ny'/), 'A',  'unitless', 'net TOA shortwave fluxes at CRM grids')
   call addfld ('CRM_FSNTC', (/'crm_nx','crm_ny'/), 'A',  'unitless', 'net TOA clear-sky shortwave fluxes at CRM grids')
   call addfld ('CRM_FSNS',  (/'crm_nx','crm_ny'/), 'A',  'unitless', 'net surface shortwave fluxes at CRM grids')
   call addfld ('CRM_FSNSC', (/'crm_nx','crm_ny'/), 'A',  'unitless',  'net surface clear-sky shortwave fluxes at CRM grids')
   call addfld ('CRM_FLNT',  (/'crm_nx','crm_ny'/), 'A',  'unitless', 'net TOA longwave fluxes at CRM grids')
   call addfld ('CRM_FLNTC', (/'crm_nx','crm_ny'/), 'A',  'unitless', 'net TOA clear-sky longwave fluxes at CRM grids')
   call addfld ('CRM_FLNS',  (/'crm_nx','crm_ny'/), 'A',  'unitless', 'net surface longwave fluxes at CRM grids')
   call addfld ('CRM_FLNSC', (/'crm_nx','crm_ny'/), 'A',  'unitless', 'net surface clear-sky longwave fluxes at CRM grids')
 
   call addfld ('CRM_AODVIS', (/'crm_nx','crm_ny'/),          'A', 'unitless', 'Aerosol optical depth at 550nm in CRM grids',&
                                                                                                            flag_xyfill=.true.)
   call addfld ('CRM_AOD400', (/'crm_nx','crm_ny'/),          'A', 'unitless', 'Aerosol optical depth at 400nm in CRM grids',&
                                                                                                            flag_xyfill=.true.)
   call addfld ('CRM_AOD700', (/'crm_nx','crm_ny'/),          'A', 'unitless', 'Aerosol optical depth at 700nm in CRM grids', &
                                                                                                            flag_xyfill=.true.)
   call addfld ('CRM_AODVISZ',(/'crm_nx','crm_ny','crm_nz'/), 'A', 'unitless',  &
                                               'Aerosol optical depth at each layer at 500nm in CRM grids', flag_xyfill=.true.)
   call addfld ('AOD400',      horiz_only,                    'A', 'unitless', 'Aerosol optical depth at 400nm', &
                                                                                                            flag_xyfill=.true.)
   call addfld ('AOD700',      horiz_only,                    'A', 'unitless', 'Aerosol optical depth at 700nm', &
                                                                                                            flag_xyfill=.true.)
   call add_default ('AOD400',  1, ' ')
   call add_default ('AOD700',  1, ' ')
    
   call addfld ('CRM_TK',  (/'crm_nx','crm_ny','crm_nz'/), 'A','m^2/s',   'Eddy viscosity from CRM')
   call addfld ('CRM_TKH', (/'crm_nx','crm_ny','crm_nz'/), 'A','m^2/s',   'Eddy viscosity from CRM')

#ifdef ECPP
   if (use_ECPP) then
      ! Add history coordinate variables for ECPP class dimensions
      call add_hist_coord('NCLASS_CL'    ,NCLASS_CL,'NCLASS_CL')
      call add_hist_coord('ncls_ecpp_in' ,ncls_ecpp_in,'ncls_ecpp_in')
      call add_hist_coord('NCLASS_PR'    ,NCLASS_PR,'NCLASS_PR')

      ! Add ECPP output variables
      call addfld ('ABND    ', (/'ilev        ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'fraction', &
                   'cloud fraction for each sub-sub class for full time period at layer boundary')
      call addfld ('ABND_TF ', (/'ilev        ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'fraction', &
                   'cloud fraction for each sub-sub class for end-portion of time period at layer boundary')
      call addfld ('MASFBND ', (/'ilev        ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'kg/m2/s',  &
                   'sub-class vertical mass flux (kg/m2/s) at layer boundary')
      call addfld ('ACEN    ', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'fraction', &
                   'cloud fraction for each sub-sub class for full time period at layer center')
      call addfld ('ACEN_TF ', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'fraction', &
                   'cloud fraction for each sub-sub class for end-portion of time period at layer center')
      call addfld ('RHCEN   ', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'fraction', &
                   'relative humidity for each sub-sub calss at layer center')
      call addfld ('QCCEN   ', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'kg/kg',    &
                   'cloud water for each sub-sub class at layer center')
      call addfld ('QICEN   ', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'kg/kg',    &
                   'cloud ice for each sub-sub class at layer center')
      call addfld ('QSINK_AFCEN', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', '/s',    &
                   'cloud water loss rate from precip. using cloud water after precip. for each sub-sub class at layer center')
      call addfld ('QSINK_BFCEN', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', '/s',    &
                  'cloud water loss rate from precip. using cloud water before precip. for each sub-sub class at layer center')
      call addfld ('QSINK_AVGCEN', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', '/s',   &
       'cloud water loss rate from precip. using averaged cloud water and precip. rate for each sub-sub class at layer center')
      call addfld ('PRAINCEN', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'kg/kg/s',  &
                   ' cloud water loss rate from precipitation (kg/kg/s) for each sub-sub class at layer center')
      call addfld ('PRECRCEN', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'kg/m2/s',  &
                   'liquid (rain) precipitation rate for each sub-sub class at layer center')
      call addfld ('PRECSCEN', (/'lev         ','NCLASS_CL   ','ncls_ecpp_in','NCLASS_PR   '/), 'A', 'kg/m2/s',  &
                   'solid (snow, graupel,...) precipitation rate for each sub-sub class at layer center')
      call addfld ('WUPTHRES',      (/ 'ilev' /), 'A', 'm/s',    'vertical velocity threshold for updraft')
      call addfld ('WDNTHRES',      (/ 'ilev' /), 'A', 'm/s',    'vertical velocity threshold for dndraft')
      call addfld ('WWQUI_CEN',     (/ 'lev' /),  'A', 'm2/s2',  'vertical velocity variance in the quiescent class, layer center')
      call addfld ('WWQUI_CLD_CEN', (/ 'lev' /),  'A', 'm2/s2',  &
                                                       'vertical velocity variance in the cloudy quiescent class, layer center')
      call addfld ('WWQUI_BND',     (/ 'ilev' /), 'A', 'm2/s2',  &
                                                       'vertical velocity variance in the quiescent class, layer boundary')
      call addfld ('WWQUI_CLD_BND', (/ 'ilev' /), 'A', 'm2/s2',  &
                                                       'vertical velocity variance in the cloudy quiescent class, layer boundary')
   endif
#endif /* ECPP */

   call addfld ('MU_CRM   ', (/ 'lev' /), 'A', 'Pa/s',     'mass flux up from CRM')
   call addfld ('MD_CRM   ', (/ 'lev' /), 'A', 'Pa/s',     'mass flux down from CRM')
   call addfld ('DU_CRM   ', (/ 'lev' /), 'A', '/s',       'detrainment from updraft from CRM')
   call addfld ('EU_CRM   ', (/ 'lev' /), 'A', '/s',       'entraiment rate from updraft')
   call addfld ('ED_CRM   ', (/ 'lev' /), 'A', '/s',       'entraiment rate from downdraft')
 
   do m=1, pcnst 
     if(cnst_name(m) == 'DMS') then 
        call addfld('DMSCONV',   (/ 'lev' /), 'A', 'kg/kg/s',  'DMS tendency from ZM convection')
     end if
     if(cnst_name(m) == 'SO2') then 
        call addfld('SO2CONV',   (/ 'lev' /), 'A', 'kg/kg/s',  'SO2 tendency from ZM convection')
      end if
   end do
 
   call addfld ('SPQRL    ', (/ 'lev' /), 'A', 'K/s',      'long-wave heating rate')
   call addfld ('SPQRS    ', (/ 'lev' /), 'A', 'K/s',      'short-wave heating rate')
   call addfld ('LENGC    ', (/ 'ilev' /), 'A', 'm  ',      'Mixing length scale for the calcuation of vertical difusivity')
 
   call addfld ('SPKVH     ',(/ 'ilev' /), 'A', 'm2/s    ', 'Vertical diffusivity used in dropmixnuc in the MMF call')
   call addfld ('SPWTKE   ', (/ 'lev' /), 'A', 'm/s',      'Standard deviation of updraft velocity')
   call addfld ('SPLCLOUD  ',(/ 'lev' /), 'A', '        ', 'Liquid cloud fraction')
#ifdef MAML
   call addfld ('CRM_SHF ',(/'crm_nx','crm_ny'/),           'I', 'W/m2    ', 'CRM Sfc sensible heat flux'          )
   call addfld ('CRM_LHF ',(/'crm_nx','crm_ny'/),           'I', 'W/m2    ', 'CRM Sfc latent heat flux'            )
   call addfld ('CRM_SNOW',(/'crm_nx','crm_ny'/),           'I', 'm/s     ', 'CRM Snow Rate'                       )
   call addfld ('CRM_PCP ',(/'crm_nx','crm_ny'/),           'I', 'm/s     ', 'CRM Precipitation Rate'              )
#endif

   call addfld ('SPNDROPMIX',(/ 'lev' /),'A','#/kg/s  ','Droplet number mixing')
   call addfld ('SPNDROPSRC',(/ 'lev' /),'A','#/kg/s  ','Droplet number source')
   call addfld ('SPNDROPCOL',horiz_only,'A', '#/m2    ','Column droplet number')

   call add_default ('SPKVH     ', 1, ' ')
   call add_default ('SPWTKE    ', 1, ' ')
   call add_default ('SPLCLOUD  ', 1, ' ')
   call add_default ('SPNDROPSRC', 1, ' ')
   call add_default ('SPNDROPMIX', 1, ' ')
   call add_default ('SPNDROPCOL', 1, ' ')

#ifdef MODAL_AERO
   ! add dropmixnuc tendencies for all modal aerosol species
   do m = 1, ntot_amode
      do lphase = 1, 2
         do lspec = 0, nspec_amode(m)+1   ! loop over number + chem constituents + water
            unit = 'kg/m2/s'
            if (lspec == 0) then   ! number
               unit = '#/m2/s'
               if (lphase == 1) then
                  l = numptr_amode(m)
               else
                  l = numptrcw_amode(m)
               endif
            else if (lspec <= nspec_amode(m)) then   ! non-water mass
               if (lphase == 1) then
                  l = lmassptr_amode(lspec,m)
               else
                  l = lmassptrcw_amode(lspec,m)
               endif
            else   ! water mass
               cycle
            end if
            if (lphase == 1) then
               tmpname = cnst_name(l)
            else
               tmpname = cnst_name_cw(l)
            end if

            fieldname = trim(tmpname) // '_mixnuc1sp'
            long_name = trim(tmpname) // ' dropmixnuc mixnuc column tendency in the mmf one '
            call addfld( fieldname,  horiz_only, 'A', unit, long_name)
            call add_default( fieldname, 1, ' ' )

         end do   ! lspec
      end do   ! lphase
   end do   ! m  

   do m=1, pcnst
      if(species_class(m).eq.spec_class_gas) then
         fieldname = trim(cnst_name(m)) // '_mixnuc1sp'
         long_name = trim(cnst_name(m)) // ' dropmixnuc mixnuc column tendency in the mmf one '
         call addfld( fieldname,  horiz_only, 'A', unit, long_name)
         call add_default( fieldname, 1, ' ' )
      end if
   end do
#endif

   call addfld ('PRES    ',(/ 'lev' /), 'A', 'Pa      ','Pressure'                                )
   call addfld ('DPRES   ',(/ 'lev' /), 'A', 'Pa      ','Pressure thickness of layer'             )
   call addfld ('SPDT    ',(/ 'lev' /), 'A', 'K/s     ','T tendency due to CRM'                   )
   call addfld ('SPDQ    ',(/ 'lev' /), 'A', 'kg/kg/s ','Q tendency due to CRM'                   )
   call addfld ('SPDQC   ',(/ 'lev' /), 'A', 'kg/kg/s ','QC tendency due to CRM'                  )
   call addfld ('SPDQI   ',(/ 'lev' /), 'A', 'kg/kg/s ','QI tendency due to CRM'                  )
   call addfld ('SPMC    ',(/ 'lev' /), 'A', 'kg/m2/s ','Total mass flux from CRM'                )
   call addfld ('SPMCUP  ',(/ 'lev' /), 'A', 'kg/m2/s ','Updraft mass flux from CRM'              )
   call addfld ('SPMCDN  ',(/ 'lev' /), 'A', 'kg/m2/s ','Downdraft mass flux from CRM'            )
   call addfld ('SPMCUUP ',(/ 'lev' /), 'A', 'kg/m2/s ','Unsaturated updraft mass flux from CRM'  )
   call addfld ('SPMCUDN ',(/ 'lev' /), 'A', 'kg/m2/s ','Unsaturated downdraft mass flux from CRM')
   call addfld ('SPQC    ',(/ 'lev' /), 'A', 'kg/kg   ','Cloud water from CRM'                    )
   call addfld ('SPQI    ',(/ 'lev' /), 'A', 'kg/kg   ','Cloud ice from CRM'                      )
   call addfld ('SPQS    ',(/ 'lev' /), 'A', 'kg/kg   ','Snow from CRM'                           )
   call addfld ('SPQG    ',(/ 'lev' /), 'A', 'kg/kg   ','Graupel from CRM'                        )
   call addfld ('SPQR    ',(/ 'lev' /), 'A', 'kg/kg   ','Rain from CRM'                           )
   call addfld ('SPQTFLX ',(/ 'lev' /), 'A', 'kg/m2/s ','Nonprecip. water flux from CRM'          )
   call addfld ('SPUFLX  ',(/ 'lev' /), 'A', 'm2/s2   ','x-momentum flux from CRM'                )
   call addfld ('SPVFLX  ',(/ 'lev' /), 'A', 'm2/s2   ','y-momentum flux from CRM'                )
   call addfld ('SPQTFLXS',(/ 'lev' /), 'A', 'kg/m2/s ','SGS Nonprecip. water flux from CRM'      )
   call addfld ('SPTKE   ',(/ 'lev' /), 'A', 'kg/m/s2 ','Total TKE in CRM'                        )
   call addfld ('SPTKES  ',(/ 'lev' /), 'A', 'kg/m/s2 ','SGS TKE in CRM'                          )
   call addfld ('SPTK    ',(/ 'lev' /), 'A', 'm2/s    ','SGS TK in CRM'                           )
   call addfld ('SPQPFLX ',(/ 'lev' /), 'A', 'kg/m2/s ','Precip. water flux from CRM'             )
   call addfld ('SPPFLX  ',(/ 'lev' /), 'A', 'm/s     ','Precipitation flux from CRM'             )
   call addfld ('SPQTLS  ',(/ 'lev' /), 'A', 'kg/kg/s ','L.S. Vapor Tendency from CRM'            )
   call addfld ('SPQTTR  ',(/ 'lev' /), 'A', 'kg/kg/s ','Nonprec. water transport from CRM'       )
   call addfld ('SPQPTR  ',(/ 'lev' /), 'A', 'kg/kg/s ','Prec. water transport from CRM'          )
   call addfld ('SPQPEVP ',(/ 'lev' /), 'A', 'kg/kg/s ','Prec. water evaporation from CRM'        )
   call addfld ('SPQPFALL',(/ 'lev' /), 'A', 'kg/kg/s ','Prec. water fall-out from CRM'           )
   call addfld ('SPQPSRC ',(/ 'lev' /), 'A', 'kg/kg/s ','Prec. water source from CRM'             )
   call addfld ('SPTLS   ',(/ 'lev' /), 'A', 'kg/kg/s ','L.S. LIWSE Tendency from CRM'            )
   call addfld ('TIMINGF ', horiz_only, 'A', '        ','CRM CPU usage efficiency: 1 - ideal'     )
   call addfld ('CLOUDTOP',(/ 'lev' /), 'A', '        ','Cloud Top PDF'                           )
#if defined(SPMOMTRANS) || defined(SP_ESMT)
   call addfld ('UCONVMOM',(/ 'lev' /), 'A', 'm/s2 ','U tendency due to CRM'            )
   call addfld ('VCONVMOM',(/ 'lev' /), 'A', 'm/s2 ','V tendency due to CRM'            )
#endif
#if defined(SP_ESMT)
   call addfld ('U_ESMT',(/ 'lev' /), 'A', 'm/s2 ','U tendency due to CRM (ESMT)'            )
   call addfld ('V_ESMT',(/ 'lev' /), 'A', 'm/s2 ','V tendency due to CRM (ESMT)'            )
#endif

   call addfld ('CRM_U   ',(/'crm_nx','crm_ny', 'crm_nz'/), 'I', 'm/s     ', 'CRM x-wind'                          )
   call addfld ('CRM_V   ',(/'crm_nx','crm_ny', 'crm_nz'/), 'I', 'm/s     ', 'CRM y-wind'                          )
   call addfld ('CRM_W   ',(/'crm_nx','crm_ny', 'crm_nz'/), 'I', 'm/s     ', 'CRM z-wind'                          )
   call addfld ('CRM_T   ',(/'crm_nx','crm_ny', 'crm_nz'/), 'I', 'K       ', 'CRM Temperature'                     )
   call addfld ('CRM_QV  ',(/'crm_nx','crm_ny', 'crm_nz'/), 'I', 'kg/kg   ', 'CRM Water Vapor'                     )
   call addfld ('CRM_QC  ',(/'crm_nx','crm_ny', 'crm_nz'/), 'I', 'kg/kg   ', 'CRM Cloud Water'                     )
   call addfld ('CRM_QI  ',(/'crm_nx','crm_ny', 'crm_nz'/), 'I', 'kg/kg   ', 'CRM Cloud Ice'                       )
   call addfld ('CRM_QPC ',(/'crm_nx','crm_ny', 'crm_nz'/), 'I', 'kg/kg   ', 'CRM Precipitating Water'             )
   call addfld ('CRM_QPI ',(/'crm_nx','crm_ny', 'crm_nz'/), 'I', 'kg/kg   ', 'CRM Precipitating Ice'               )
   call addfld ('CRM_PREC',(/'crm_nx','crm_ny'/),           'I', 'm/s     ', 'CRM Precipitation Rate'              )

   !-- MDB 8/2013
   call addfld ('SPTVFLUX ',(/ 'lev' /), 'A', 'W/m2  ','Buoyancy Flux from CRM'             )
   call addfld ('SPBUOY   ',(/ 'lev' /), 'A', 'W/m3  ','Buoyancy Term from CRM'             )
   call addfld ('SPBUOYSD ',(/ 'lev' /), 'A', 'W/m3  ','Std Dev of Buoyancy Term from CRM'  )
   call addfld ('SPMSEF   ',(/ 'lev' /), 'A', 'W/m2  ','Moist Static Energy Flux from CRM'  )
   call addfld ('SPQVFLUX ',(/ 'lev' /), 'A', 'W/m2  ','Water Wapor Flux from CRM'          )

   call add_default ('SPDT    ', 1, ' ')
   call add_default ('SPDQ    ', 1, ' ')
   call add_default ('SPDQC   ', 1, ' ')
   call add_default ('SPDQI   ', 1, ' ')
   call add_default ('SPMC    ', 1, ' ')
   call add_default ('SPMCUP  ', 1, ' ')
   call add_default ('SPMCDN  ', 1, ' ')
   call add_default ('SPMCUUP ', 1, ' ')
   call add_default ('SPMCUDN ', 1, ' ')
   call add_default ('SPQC    ', 1, ' ')
   call add_default ('SPQI    ', 1, ' ')
   call add_default ('SPQS    ', 1, ' ')
   call add_default ('SPQG    ', 1, ' ')
   call add_default ('SPQR    ', 1, ' ')
   call add_default ('SPQTFLX ', 1, ' ')
   call add_default ('SPQTFLXS', 1, ' ')
   call add_default ('SPTKE   ', 1, ' ')
   call add_default ('SPTKES  ', 1, ' ')
   call add_default ('SPTK    ', 1, ' ')
   call add_default ('SPQPFLX ', 1, ' ')
   call add_default ('SPPFLX  ', 1, ' ')
   call add_default ('SPQTLS  ', 1, ' ')
   call add_default ('SPQTTR  ', 1, ' ')
   call add_default ('SPQPTR  ', 1, ' ')
   call add_default ('SPQPEVP ', 1, ' ')
   call add_default ('SPQPFALL', 1, ' ')
   call add_default ('SPQPSRC ', 1, ' ')
   call add_default ('SPTLS   ', 1, ' ')
   call add_default ('CLOUDTOP', 1, ' ')
   call add_default ('TIMINGF ', 1, ' ')

   call add_default ('SPTVFLUX  ', 1, ' ')
   call add_default ('SPBUOY    ', 1, ' ')
   call add_default ('SPBUOYSD  ', 1, ' ')
   call add_default ('SPMSEF    ', 1, ' ')
   call add_default ('SPQVFLUX  ', 1, ' ')

#if defined(SPMOMTRANS) || defined(SP_ESMT)
   call add_default ('UCONVMOM', 1, ' ')
   call add_default ('VCONVMOM', 1, ' ')
#endif

#if defined(SP_ESMT)
   call add_default ('U_ESMT', 1, ' ')
   call add_default ('V_ESMT', 1, ' ')
#endif

   prec_dp_idx  =  pbuf_get_index('PREC_DP')
   snow_dp_idx  =  pbuf_get_index('SNOW_DP')
   prec_sh_idx  =  pbuf_get_index('PREC_SH')
   snow_sh_idx  =  pbuf_get_index('SNOW_SH')
   prec_sed_idx =  pbuf_get_index('PREC_SED')
   snow_sed_idx =  pbuf_get_index('SNOW_SED')
   snow_str_idx =  pbuf_get_index('SNOW_STR')
   prec_pcw_idx =  pbuf_get_index('PREC_PCW')
   snow_pcw_idx =  pbuf_get_index('SNOW_PCW')

   call crm_accel_init()

end subroutine crm_physics_init

!===================================================================================================
!===================================================================================================

subroutine crm_physics_tend(ztodt, state, tend, ptend, pbuf, cam_in, cam_out,    &
                            species_class, crm_ecpp_output,          &
                            sp_qchk_prec_dp, sp_qchk_snow_dp, sp_rad_flux)

   !------------------------------------------------------------------------------------------------
   !
   ! CRM interface for Super-Parameterization
   !
   ! Author: Marat Khairoutdinov (mkhairoutdin@ms.cc.sunysb.edu)
   !
   ! Purpose: to update state from CRM physics. 
   ! 
   ! Revision history: 
   !
   ! June, 2009, Minghuai Wang: 
   !    These codes are taken out from tphysbc.F90 in the spcam3.5, developed by Marat Khairoutdinov 
   !    (mkhairoutdin@ms.cc.sunysb.edu). Here we try to follow the procedure in 'Interface to Column
   !    Physics and Chemistry packages' to implement the CRM physics.
   !    
   ! July, 13, 2009, Minghuai Wang: 
   !    Hydrometer numbers are output from SAM when Morrison's microphysics is used, and will be
   !    used in the radiative transfer code to calculate effective radius. 
   !      
   ! July, 15, 2009, Minghuai Wang: 
   !    Get modal aerosol, and use it in the SAM. 
   ! 
   !------------------------------------------------------------------------------------------------
   use ppgrid
   use physics_buffer,  only: physics_buffer_desc, pbuf_old_tim_idx, pbuf_get_index, &
                              dyn_time_lvls, pbuf_get_field, pbuf_set_field
   use physics_types,   only: physics_state, physics_tend, physics_ptend, physics_ptend_init
   use camsrfexch,      only: cam_in_t, cam_out_t
   use time_manager,    only: is_first_step, get_nstep
   use cam_history,     only: outfld
   use perf_mod
   use crmdims,         only: crm_nx, crm_ny, crm_nz, crm_nx_rad, crm_ny_rad
   use physconst,       only: cpair, latvap, latice, gravit, cappa
   use constituents,    only: pcnst, qmin, cnst_get_ind, cnst_cam_outfld, bpcnst, cnst_name
   use crm_module,      only: crm
   use params,          only: crm_rknd
   use physconst,       only: latvap
   use phys_control,    only: phys_getopts
   use check_energy,    only: check_energy_chng

! modal_aero_data only exists if MODAL_AERO
#if (defined  m2005 && defined MODAL_AERO)  
   use modal_aero_data, only: ntot_amode, ntot_amode
   use ndrop,  only: loadaer
   use microphysics,  only: iqv, iqci, iqr, iqs, iqg, incl, inci, inr, ing, ins   !!!!!! BE CAUIOUS, these indices can only defined before call to crm.
#endif
   use phys_grid,       only: get_rlat_all_p, get_rlon_all_p, get_gcol_p
   !!!use aerosol_intr,    only: aerosol_wet_intr

#if defined( SP_ORIENT_RAND )
   use RNG_MT            ! random number generator for randomly rotating CRM orientation (SP_ORIENT_RAND)
#endif

   use crm_state_module,       only: crm_state_type
   use crm_rad_module,         only: crm_rad_type, crm_rad_initialize, crm_rad_finalize
   use crm_input_module,       only: crm_input_type
   use crm_output_module,      only: crm_output_type, crm_output_initialize, crm_output_finalize
   use crm_ecpp_output_module, only: crm_ecpp_output_type

   real(r8),                   intent(in   ) :: ztodt            ! global model time increment
   type(physics_state),        intent(in   ) :: state            ! Global model state 
   type(physics_tend),         intent(in   ) :: tend             ! 
   type(physics_ptend),        intent(  out) :: ptend            ! output tendencies
   type(physics_buffer_desc),  pointer       :: pbuf(:)          ! physics buffer
   type(cam_in_t),             intent(in   ) :: cam_in           ! atm input from coupler
   type(cam_out_t),            intent(inout) :: cam_out          ! atm output to coupler
   integer,                    intent(in   ) :: species_class(:) ! aerosol species type
   type(crm_ecpp_output_type), intent(inout) :: crm_ecpp_output  ! output data for ECPP calculations
   real(r8), dimension(pcols), intent(out  ) :: sp_qchk_prec_dp  ! precipitation diagostic (liq+ice)  used for check_energy_chng
   real(r8), dimension(pcols), intent(out  ) :: sp_qchk_snow_dp  ! precipitation diagostic (ice only) used for check_energy_chng
   real(r8), dimension(pcols), intent(out  ) :: sp_rad_flux      ! radiative flux diagnostic used for check_energy_chng

   
   ! real(r8), intent(in) :: dlf(pcols,pver)  ! shallow+deep convective detrainment [kg/kg/s] - used for aerosol_wet_intr - no longer needed

   !------------------------------------------------------------------------------------------------
   ! Local variables 
   !------------------------------------------------------------------------------------------------

   ! convective precipitation variables
   real(r8), pointer :: prec_dp(:)          ! total precip from deep convection (ZM)    [m/s]
   real(r8), pointer :: snow_dp(:)          ! snow from deep convection (ZM)            [m/s]
   real(r8), pointer :: prec_sh(:)          ! total precip from shallow convection      [m/s]
   real(r8), pointer :: snow_sh(:)          ! snow from shallow convection              [m/s]

   ! stratiform precipitation variables
   real(r8), pointer :: prec_pcw(:)         ! total precip from prognostic cloud scheme   [m/s]
   real(r8), pointer :: snow_pcw(:)         ! snow from prognostic cloud scheme           [m/s]
   real(r8), pointer :: prec_sed(:)         ! total precip from cloud sedimentation       [m/s]
   real(r8), pointer :: snow_sed(:)         ! snow from cloud ice sedimentation           [m/s]
   real(r8), pointer :: snow_str(:)         ! snow from stratiform cloud                  [m/s]

   integer lchnk                    ! chunk identifier
   integer ncol                     ! number of atmospheric columns
   integer  nstep                   ! time steps
   real(r8) crm_run_time            ! length of CRM integration
   
   character(len=16) :: microp_scheme  ! microphysics scheme
   real(r8) :: cwp      (pcols,pver)   ! in-cloud cloud (total) water path (kg/m2)
   real(r8) :: gwp      (pcols,pver)   ! grid-box cloud (total) water path (kg/m2)
   real(r8) :: tgicewp  (pcols)        ! Vertically integrated ice water path (kg/m2
   real(r8) :: tgliqwp  (pcols)        ! Vertically integrated liquid water path (kg/m2)
   real(r8) :: cicewp   (pcols,pver)   ! in-cloud cloud ice water path (kg/m2)
   real(r8) :: cliqwp   (pcols,pver)   ! in-cloud cloud liquid water path (kg/m2)
   real(r8) :: tgwp     (pcols)        ! Vertically integrated (total) cloud water path  (kg/m2)
   real(r8) :: ftem     (pcols,pver)   ! Temporary workspace for outfld variables
   real(r8) :: ideep_crm(pcols)

   ! physics buffer fields to compute tendencies for stratiform package
   real(r8), pointer, dimension(:,:) :: cld        ! cloud fraction

#if (defined m2005 && defined MODAL_AERO)
   real(r8) :: na(pcols)  ! aerosol number concentration [/m3]
   real(r8) :: va(pcols)  ! aerosol voume concentration [m3/m3]
   real(r8) :: hy(pcols)  ! aerosol bulk hygroscopicity
   integer  :: phase      ! phase to determine whether it is interstitial, cloud-borne, or the sum. 
#endif

   real(r8) :: cs(pcols, pver)  ! air density  [kg/m3]

   ! CRM column radiation stuff:
   real(r8), pointer, dimension(:,:) :: qrs        ! shortwave radiative heating rate
   real(r8), pointer, dimension(:,:) :: qrl        ! shortwave radiative heating rate

   ! TODO: this should not be a pointer, it is being used as an alloctable (and
   ! is only used in one place, so should get a more descriptive variable name)
   real(r8), pointer, dimension(:,:)     :: tempPtr

   real(r8), dimension(pcols) :: qli_hydro_before    ! column-integraetd rain + snow + graupel 
   real(r8), dimension(pcols) ::  qi_hydro_before    ! column-integrated snow water + graupel water
   real(r8), dimension(pcols) :: qli_hydro_after     ! column-integraetd rain + snow + graupel 
   real(r8), dimension(pcols) ::  qi_hydro_after     ! column-integrated snow water + graupel water
   real(r8) :: sfactor                               ! used to determine precip type for sam1mom

   integer :: ii, jj
   integer :: ixcldliq, ixcldice, ixnumliq, ixnumice
   integer :: ixrain, ixsnow, ixnumrain, ixnumsnow
   integer :: i, k, m
   integer :: ifld, itim
   logical :: use_ECPP, use_SPCAM
   character(len=16) :: SPCAM_microp_scheme

   logical :: ls, lu, lv, lq(pcnst)

   real(r8) dp_g   ! = state%pdel / gravit

   integer :: icol(pcols)

   !!! variables for changing CRM orientation
   real(crm_rknd), parameter        :: pi   = 3.14159265359
   real(crm_rknd), parameter        :: pix2 = 6.28318530718
   real(crm_rknd), dimension(pcols) :: crm_angle

   ! CRM types
   type(crm_state_type)  :: crm_state
   type(crm_rad_type)    :: crm_rad
   type(crm_input_type)  :: crm_input
   type(crm_output_type) :: crm_output
#ifdef MAML
   real(r8), pointer, dimension(:,:,:)   :: crm_pcp
   real(r8), pointer, dimension(:,:,:)   :: crm_snw
   real(r8) :: factor_xy
   real(r8) :: tau00_avg, bflxls_avg, fluxu00_avg, fluxv00_avg 
   real(r8) :: fluxt00_avg, fluxq00_avg 
   factor_xy = 1._r8/dble(crm_nx*crm_ny)
#endif

#if defined( SP_ORIENT_RAND )
   real(crm_rknd) :: unif_rand1           ! uniform random number 
   real(crm_rknd) :: unif_rand2           ! uniform random number 
   real(crm_rknd) :: norm_rand            ! normally distributed random number using the Box-Muller (1958) method
   real(crm_rknd) :: crm_rotation_std     ! scaling factor for CRM rotation (std deviation of rotation angle)
   real(crm_rknd) :: crm_rotation_offset  ! offset to specify preferred rotation direction 
   integer :: seed

   crm_rotation_std    = 20. * pi/180.                 ! std deviation of normal distribution for CRM rotation [radians]
   crm_rotation_offset = 90. * pi/180. * ztodt/86400.  ! This means that a CRM should rotate 90 deg / day on average
#endif

   crm_run_time = ztodt

   
   call phys_getopts(use_ECPP_out            = use_ECPP)
   call phys_getopts(use_SPCAM_out           = use_SPCAM)
   call phys_getopts(SPCAM_microp_scheme_out = SPCAM_microp_scheme)

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

#if defined( SP_ORIENT_RAND )
   !------------------------------------------------------------------------------------------------
   ! Rotate the CRM using a random walk
   !------------------------------------------------------------------------------------------------
   if ( (crm_ny.eq.1) .or. (crm_nx.eq.1) ) then

      do i=1,ncol

         !!! set the seed based on the chunk and column index (duplicate seeds are ok)
         seed = lchnk + i + nstep

         call RNG_MT_set_seed(seed)

         !!! Generate a pair of uniform random numbers
         call RNG_MT_gen_rand(unif_rand1)
         call RNG_MT_gen_rand(unif_rand2)

         !!! Box-Muller (1958) method of obtaining a Gaussian distributed random number
         norm_rand = sqrt(-2.*log(unif_rand1))*cos(pix2*unif_rand2)
         crm_angle(i) = crm_angle(i) + norm_rand * crm_rotation_std + crm_rotation_offset

         !!! Adjust CRM orientation angle to be between 0 and 2*pi
         if ( crm_angle(i) .lt. 0. ) then
            crm_angle(i) = crm_angle(i) + pix2
         endif
         if ( crm_angle(i) .gt. pix2 ) then
            crm_angle(i) = crm_angle(i) - pix2
         endif

      enddo

   endif
#else /* SP_ORIENT_RAND */
   !------------------------------------------------------------------------------------------------
   ! use static CRM orientation (no rotation)
   !------------------------------------------------------------------------------------------------
#if defined( SP_DIR_NS )
    if (crm_ny.eq.1) then
       crm_angle(:ncol) = pi/2.
    else 
       crm_angle(:ncol) = 0.
    endif
#else
      crm_angle(:ncol) = 0.
#endif /* SP_DIR_NS */

#endif /* SP_ORIENT_RAND */

   !------------------------------------------------------------------------------------------------
   ! Retreive pbuf fields
   !------------------------------------------------------------------------------------------------
   if (SPCAM_microp_scheme .eq. 'm2005') then
     call pbuf_get_field(pbuf, pbuf_get_index('CRM_NC_RAD'), crm_rad%nc, start=(/1,1,1,1/), kount=(/pcols,crm_nx_rad, crm_ny_rad, crm_nz/))
     call pbuf_get_field(pbuf, pbuf_get_index('CRM_NI_RAD'), crm_rad%ni, start=(/1,1,1,1/), kount=(/pcols,crm_nx_rad, crm_ny_rad, crm_nz/))
     call pbuf_get_field(pbuf, pbuf_get_index('CRM_QS_RAD'), crm_rad%qs, start=(/1,1,1,1/), kount=(/pcols,crm_nx_rad, crm_ny_rad, crm_nz/))
     call pbuf_get_field(pbuf, pbuf_get_index('CRM_NS_RAD'), crm_rad%ns, start=(/1,1,1,1/), kount=(/pcols,crm_nx_rad, crm_ny_rad, crm_nz/))
   endif

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

   !!! set convective rain to be zero for PRAIN already includes precipitation production from convection. 
   ! call pbuf_set_field(pbuf, pbuf_get_index('RPRDTOT'), 0.0_r8 )
   ! call pbuf_set_field(pbuf, pbuf_get_index('RPRDDP' ), 0.0_r8 )
   ! call pbuf_set_field(pbuf, pbuf_get_index('RPRDSH' ), 0.0_r8 )
   ! call pbuf_set_field(pbuf, pbuf_get_index('ICWMRDP'), 0.0_r8 )
   ! call pbuf_set_field(pbuf, pbuf_get_index('ICWMRSH'), 0.0_r8 )
   
   prec_dp  = 0.
   snow_dp  = 0.
   prec_sh  = 0.
   snow_sh  = 0. 
   prec_sed = 0.
   snow_sed = 0.
   snow_str = 0.
   prec_pcw = 0
   snow_pcw = 0.
   
   !!! Initialize stuff:
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
   if (SPCAM_microp_scheme .eq. 'sam1mom') then
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
#ifdef MAML
   call pbuf_get_field (pbuf, crm_pcp_idx, crm_pcp)
   call pbuf_get_field (pbuf, crm_snw_idx, crm_snw)
#endif

   !------------------------------------------------------------------------------------------------
   !------------------------------------------------------------------------------------------------

   if(is_first_step()) then
      ! call check_energy_timestep_init(state, tend, pbuf)
      ! initialize crm_state%qt to zero (needed for ncol < i <= pcols)
      crm_state%qt(:,:,:,:) = 0.0_r8
      do i=1,ncol
         do k=1,crm_nz
            m = pver-k+1

            ! Initialize CRM state
            crm_state%u_wind(i,:,:,k) = state%u(i,m) * cos( crm_angle(i) ) + state%v(i,m) * sin( crm_angle(i) )
            crm_state%v_wind(i,:,:,k) = state%v(i,m) * cos( crm_angle(i) ) - state%u(i,m) * sin( crm_angle(i) )
            crm_state%w_wind(i,:,:,k) = 0.
            crm_state%temperature(i,:,:,k) = state%t(i,m)

            ! Initialize microphysics arrays
            if (SPCAM_microp_scheme .eq. 'sam1mom') then
               crm_state%qt(i,:,:,k) = state%q(i,m,1)+state%q(i,m,ixcldliq)+state%q(i,m,ixcldice)
               crm_state%qp(i,:,:,k) = 0.0_r8
               crm_state%qn(i,:,:,k) = state%q(i,m,ixcldliq)+state%q(i,m,ixcldice)
            else if (SPCAM_microp_scheme .eq. 'm2005') then
               crm_state%qt(i,:,:,k) = state%q(i,m,1)+state%q(i,m,ixcldliq)
               crm_state%nc(i,:,:,k) = 0.0_r8
               crm_state%qr(i,:,:,k) = 0.0_r8
               crm_state%nr(i,:,:,k) = 0.0_r8
               crm_state%qi(i,:,:,k) = state%q(i,m,ixcldice)
               crm_state%ni(i,:,:,k) = 0.0_r8
               crm_state%qs(i,:,:,k) = 0.0_r8
               crm_state%ns(i,:,:,k) = 0.0_r8
               crm_state%qg(i,:,:,k) = 0.0_r8
               crm_state%ng(i,:,:,k) = 0.0_r8
               crm_state%qc(i,:,:,k) = state%q(i,m,ixcldliq)
            endif

         end do
      end do

      do k=1,crm_nz
         m = pver-k+1
         do i=1,ncol
            crm_rad%qrad         (i,:,:,k) = 0.
            crm_rad%temperature  (i,:,:,k) = state%t(i,m)
            crm_rad%qv           (i,:,:,k) = state%q(i,m,1)
            crm_rad%qc           (i,:,:,k) = 0.
            crm_rad%qi           (i,:,:,k) = 0.
            crm_rad%cld          (i,:,:,k) = 0.
#ifdef m2005
            if (SPCAM_microp_scheme .eq. 'm2005') then
               crm_rad%nc(i,:,:,k) = 0.0
               crm_rad%ni(i,:,:,k) = 0.0       
               crm_rad%qs(i,:,:,k) = 0.0
               crm_rad%ns(i,:,:,k) = 0.0
            endif
#endif
         end do
      end do
#ifdef MAML
      crm_pcp(:,:,:) = 0.
      crm_snw(:,:,:) = 0.
#endif

! use radiation from grid-cell mean radctl on first time step
      cld(:,:) = 0.
      ptend%q(:,:,1) = 0.
      ptend%q(:,:,ixcldliq) = 0.
      ptend%q(:,:,ixcldice) = 0.
      ptend%s(:,:) = 0.

#ifdef ECPP
      if (use_ECPP) then
         
         cs(1:ncol,1:pver) = state%pmid(1:ncol,1:pver) / (287.15*state%t(1:ncol,1:pver))

         !!! initialize turbulence for ECPP calculations
         ifld = pbuf_get_index('TKE_CRM')
         call pbuf_set_field(pbuf, ifld, 0.0_r8, start=(/1,1/), kount=(/pcols, pver/) )

         ifld = pbuf_get_index('TK_CRM')
         call pbuf_set_field(pbuf, ifld, 0.0_r8, start=(/1,1/), kount=(/pcols, pver/) )
      endif 
#endif


   else  ! not is_first_step

      ptend%q(:,:,1) = 0.  ! necessary?
      ptend%q(:,:,ixcldliq) = 0.
      ptend%q(:,:,ixcldice) = 0.
      ptend%s(:,:) = 0. ! necessary?
      cwp    = 0.
#ifdef MAML
      crm_pcp = 0.
      crm_snw = 0.
#endif

      do m=1,crm_nz
         k = pver-m+1
         do i = 1,ncol
            crm_rad%qrad(i,:,:,m) = crm_rad%qrad(i,:,:,m) / state%pdel(i,k) ! for energy conservation
         end do
      end do

#if (defined m2005 && defined MODAL_AERO)
      cs(1:ncol,1:pver) = state%pmid(1:ncol,1:pver) / (287.15*state%t(1:ncol,1:pver))
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
                  if (SPCAM_microp_scheme .eq. 'm2005') then
                     qli_hydro_before(i) = qli_hydro_before(i)+(crm_state%qr(i,ii,jj,m)+ &
                                                                crm_state%qs(i,ii,jj,m)+ &
                                                                crm_state%qg(i,ii,jj,m)) * dp_g
                     qi_hydro_before(i)  =  qi_hydro_before(i)+(crm_state%qs(i,ii,jj,m)+ &
                                                                crm_state%qg(i,ii,jj,m)) * dp_g
                  else if (SPCAM_microp_scheme .eq. 'sam1mom') then
                     sfactor = max(0._r8,min(1._r8,(crm_state%temperature(i,ii,jj,m)-268.16)*1./(283.16-268.16)))
                     qli_hydro_before(i) = qli_hydro_before(i)+crm_state%qp(i,ii,jj,m) * dp_g
                     qi_hydro_before(i)  =  qi_hydro_before(i)+crm_state%qp(i,ii,jj,m) * (1-sfactor) * dp_g
                  end if ! SPCAM_microp_scheme
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
      do i = 1,ncol
#ifdef MAML
         tau00_avg =0._r8
         bflxls_avg =0._r8
         fluxu00_avg =0._r8
         fluxv00_avg =0._r8
         fluxt00_avg =0._r8
         fluxq00_avg =0._r8
         do ii = 1, crm_nx*crm_ny
            !seems none of variables below is used, so I don't use the CRM-level
            !input. Later on, they can be changed if needed   
            tau00_avg = tau00_avg + sqrt(cam_in%wsx(i,ii)**2 + cam_in%wsy(i,ii)**2)
            bflxls_avg = bflxls_avg + cam_in%shf(i,ii)/cpair + 0.61*state%t(i,pver)*cam_in%lhf(i,ii)/latvap
            fluxu00_avg = fluxu00_avg + cam_in%wsx(i,ii)     !N/m2
            fluxv00_avg = fluxv00_avg + cam_in%wsy(i,ii)     !N/m2
            fluxt00_avg = fluxt00_avg + cam_in%shf(i,ii)/cpair  ! K Kg/ (m2 s)
            fluxq00_avg = fluxq00_avg + cam_in%lhf(i,ii)/latvap ! Kg/(m2 s)
         end do
         crm_input%tau00(i) = tau00_avg*factor_xy 
         crm_input%bflxls(i) = bflxls_avg*factor_xy 
         crm_input%fluxu00(i) = fluxu00_avg*factor_xy     !N/m2
         crm_input%fluxv00(i) = fluxv00_avg*factor_xy     !N/m2
         crm_input%fluxt00(i) = fluxt00_avg*factor_xy  ! K Kg/ (m2 s)
         crm_input%fluxq00(i) = fluxq00_avg*factor_xy ! Kg/(m2 s)

#else
         crm_input%tau00(i) = sqrt(cam_in%wsx(i)**2 + cam_in%wsy(i)**2)
         crm_input%bflxls(i) = cam_in%shf(i)/cpair + 0.61*state%t(i,pver)*cam_in%lhf(i)/latvap
         crm_input%fluxu00(i) = cam_in%wsx(i)     !N/m2
         crm_input%fluxv00(i) = cam_in%wsy(i)     !N/m2
         crm_input%fluxt00(i) = cam_in%shf(i)/cpair  ! K Kg/ (m2 s)
         crm_input%fluxq00(i) = cam_in%lhf(i)/latvap ! Kg/(m2 s)
#endif  
         crm_input%wndls(i) = sqrt(state%u(i,pver)**2 + state%v(i,pver)**2)
      end do
#if (defined m2005 && defined MODAL_AERO)
      ! Set aerosol
      do i = 1,ncol
         do k=1, pver
            phase = 1  ! interstital aerosols only
            do m=1, ntot_amode
               call loadaer( &
                  state, pbuf, i, i, k, &
                  m, cs, phase, na, va, &
                  hy)
               crm_input%naermod (i,k,m) = na(i)
               crm_input%vaerosol(i,k,m) = va(i)
               crm_input%hygro   (i,k,m) = hy(i)
            end do    
         end do
      end do
#endif
      !---------------------------------------------------------------------------------------------
      ! Set the input wind (also sets CRM orientation)
      !---------------------------------------------------------------------------------------------
      do i = 1,ncol
         do k=1,pver
            crm_input%ul(i,k) = state%u(i,k) * cos( crm_angle(i) ) + state%v(i,k) * sin( crm_angle(i) )
            crm_input%vl(i,k) = state%v(i,k) * cos( crm_angle(i) ) - state%u(i,k) * sin( crm_angle(i) )
#if defined( SP_ESMT )
            ! Set the input wind for ESMT
            crm_input%ul_esmt(i,k) = state%u(i,k)
            crm_input%vl_esmt(i,k) = state%v(i,k)
#endif /* SP_ESMT */
         enddo ! k=1,pver
      enddo ! i=1,ncol

      ! Set the global model column mapping
      ! TODO: get rid of this; we do not need this anymore if we are not
      ! supporting the stand-alone CRM dump routines in crm_module.
      do i = 1,ncol
         icol(i) = i
      end do

      !---------------------------------------------------------------------------------------------
      ! Run the CRM
      !---------------------------------------------------------------------------------------------
      if (.not.allocated(ptend%q)) write(*,*) '=== ptend%q not allocated ==='
      if (.not.allocated(ptend%s)) write(*,*) '=== ptend%s not allocated ==='
      call t_startf ('crm_call')
      call crm( lchnk, icol(:ncol), ncol, ztodt, pver,                    &
                crm_input, crm_state, crm_rad,                            &
#ifdef MAML
                crm_pcp(:ncol,:,:),     crm_snw(:ncol,:,:),               &
#endif
                crm_ecpp_output, crm_output )
      call t_stopf('crm_call')

      ! Copy tendencies from CRM output
      ptend%q(:ncol,:pver,1) = crm_output%qltend(1:ncol,1:pver)
      ptend%q(:ncol,:pver,ixcldliq) = crm_output%qcltend(1:ncol,1:pver)
      ptend%q(:ncol,:pver,ixcldice) = crm_output%qiltend(1:ncol,1:pver)
      ptend%s(:ncol,:pver) = crm_output%sltend(1:ncol,1:pver)

      ! There is no separate convective and stratiform precip for CRM:
      crm_output%precc(:ncol)  = crm_output%precc(:ncol) + crm_output%precl(:ncol)
      crm_output%precsc(:ncol) = crm_output%precsc(:ncol) + crm_output%precsl(:ncol)
      crm_output%precl(:ncol)  = 0
      crm_output%precsl(:ncol) = 0

      !!! these precip pointer variables are used by coupler
      !!! TODO: do we need to zero out the elements beyond ncol here?
      prec_dp(1:ncol) = crm_output%precc(1:ncol)
      snow_dp(1:ncol) = crm_output%precsc(1:ncol)

      !!! These are needed elsewhere in the model when SP_PHYS_BYPASS is used
      ! ifld = pbuf_get_index('AST'   )
      ! call pbuf_set_field(pbuf,ifld, cld   (:ncol,:pver),start=(/1,1/), kount=(/pcols,pver/) )
      ! ifld = pbuf_get_index('PRAIN' )
      ! call pbuf_set_field(pbuf,ifld, crm_output%qp_src(:ncol,:pver),start=(/1,1/), kount=(/pcols,pver/) )
      ! ifld = pbuf_get_index('NEVAPR')
      ! call pbuf_set_field(pbuf,ifld, crm_output%qp_evp(:ncol,:pver),start=(/1,1/), kount=(/pcols,pver/) )

      ! Set cloud field to output from CRM. We need to do this here because we
      ! did not set crm_output%cld as a pointer to the pbuf field, so we need to
      ! copy the data over. NOTE: I think this can be done using pbuf_set_field
      ! without making an extra pointer for cld, but I do not think we would be
      ! able to zero-out the rest of cld beyond pcols that way.
      !call pbuf_get_field(pbuf, pbuf_get_index('CLD'), cld)
      !cld(:,:) = 0
      cld(1:ncol,1:pver) = crm_output%cld(1:ncol,1:pver)

      ! Output CRM heating rate
      ! NOTE: We output the heating rate here because this is the heating that
      ! is applied to the CRM at this GCM timestep, but note that this heating
      ! rate is offset in time from when the radiative heating rate is
      ! calculated in the radiative transfer interface, which comes AFTER this
      ! routine call (because we want it to use the clouds simulated from this
      ! CRM call). Thus, comparing this heating rate with CRM_QRS + CRM_QRL
      ! output in radiation_tend will show a time lag.
      call outfld('CRM_QRAD', crm_rad%qrad(1:ncol,:,:,:), ncol, lchnk)

      ! Convert heating rate to Q*dp to conserve energy across timesteps
      do m=1,crm_nz
         k = pver-m+1
         do i = 1,ncol
            crm_rad%qrad(i,:,:,m) = crm_rad%qrad(i,:,:,m) * state%pdel(i,k) ! for energy conservation
         end do
      end do

      call outfld('PRES    ',state%pmid ,pcols   ,lchnk   )
      call outfld('DPRES   ',state%pdel ,pcols   ,lchnk   )

      call outfld('CRM_U   ',crm_state%u_wind,      pcols, lchnk)
      call outfld('CRM_V   ',crm_state%v_wind,      pcols, lchnk)
      call outfld('CRM_W   ',crm_state%w_wind,      pcols, lchnk)
      call outfld('CRM_T   ',crm_state%temperature, pcols, lchnk)

      if (SPCAM_microp_scheme .eq. 'sam1mom') then
         call outfld('CRM_QV  ',(crm_state%qt-crm_output%qcl-crm_output%qci), pcols, lchnk)
      else if (SPCAM_microp_scheme .eq. 'm2005') then 
         call outfld('CRM_QV  ',crm_state%qt-crm_output%qcl, pcols, lchnk)
      endif
      call outfld('CRM_QC  ',crm_output%qcl   ,pcols   ,lchnk   )
      call outfld('CRM_QI  ',crm_output%qci   ,pcols   ,lchnk   )
      call outfld('CRM_QPC ',crm_output%qpl  ,pcols   ,lchnk   )
      call outfld('CRM_QPI ',crm_output%qpi  ,pcols   ,lchnk   )
      call outfld('CRM_PREC',crm_output%prec_crm       ,pcols   ,lchnk   )
      call outfld('CRM_TK ', crm_output%tk(:, :, :, :)  ,pcols   ,lchnk   )  
      call outfld('CRM_TKH', crm_output%tkh(:, :, :, :)  ,pcols   ,lchnk   ) 
#ifdef MAML
       call outfld('CRM_SHF ', cam_in%shf ,pcols ,lchnk)
       call outfld('CRM_LHF ', cam_in%lhf ,pcols ,lchnk)
       call outfld('CRM_SNOW', crm_snw    ,pcols ,lchnk)
       call outfld('CRM_PCP',  crm_pcp    ,pcols ,lchnk)
#endif

#ifdef m2005
      if (SPCAM_microp_scheme .eq. 'm2005') then
         ! index is defined in ./crm/MICRO_M2005/microphysics.F90
         ! Be cautious to use them here. They are defined in crm codes, and these codes are called only 
         ! after the subroutine of crm is called. So they can only be used after the 'crm' subroutine. 
         ! some initializaiton part of crm codes should be called in the initializaation part of cam 
         ! in the future.
         ! incl, inci, ... can not be used here, for they are defined before we call them???
         ! +++mhwang
         call outfld('CRM_NC ',crm_state%nc(:, :, :, :)   ,pcols   ,lchnk   )
         call outfld('CRM_NI ',crm_state%ni(:, :, :, :)   ,pcols   ,lchnk   )
         call outfld('CRM_NR ',crm_state%nr(:, :, :, :)   ,pcols   ,lchnk   )
         call outfld('CRM_NS ',crm_state%ns(:, :, :, :)   ,pcols   ,lchnk   )
         call outfld('CRM_NG ',crm_state%ng(:, :, :, :)   ,pcols   ,lchnk   )

         call outfld('CRM_WVAR', crm_output%wvar, pcols, lchnk)

         call outfld('CRM_QR ',crm_state%qr(:, :, :, :)   ,pcols   ,lchnk   )
         call outfld('CRM_QS ',crm_state%qs(:, :, :, :)   ,pcols   ,lchnk   )
         call outfld('CRM_QG ',crm_state%qg(:, :, :, :)   ,pcols   ,lchnk   )

         ! hm 7/26/11, add new output
         call outfld('CRM_AUT', crm_output%aut, pcols, lchnk)
         call outfld('CRM_ACC', crm_output%acc, pcols, lchnk)
         call outfld('CRM_EVPC', crm_output%evpc, pcols, lchnk)
         call outfld('CRM_EVPR', crm_output%evpr, pcols, lchnk)
         call outfld('CRM_MLT', crm_output%mlt, pcols, lchnk)
         call outfld('CRM_SUB', crm_output%sub, pcols, lchnk)
         call outfld('CRM_DEP', crm_output%dep, pcols, lchnk)
         call outfld('CRM_CON', crm_output%con, pcols, lchnk)
         ! hm 8/31/11, add new output for time-mean-avg
         call outfld('A_AUT', crm_output%aut_a, pcols, lchnk)
         call outfld('A_ACC', crm_output%acc_a, pcols, lchnk)
         call outfld('A_EVPC', crm_output%evpc_a, pcols, lchnk)
         call outfld('A_EVPR', crm_output%evpr_a, pcols, lchnk)
         call outfld('A_MLT', crm_output%mlt_a, pcols, lchnk)
         call outfld('A_SUB', crm_output%sub_a, pcols, lchnk)
         call outfld('A_DEP', crm_output%dep_a, pcols, lchnk)
         call outfld('A_CON', crm_output%con_a, pcols, lchnk)
      endif ! m2005
#endif /* m2005 */

      !---------------------------------------------------------------------------------------------
      ! Add radiative heating tendency above CRM
      !---------------------------------------------------------------------------------------------

      ifld = pbuf_get_index('QRL')
      call pbuf_get_field(pbuf, ifld, qrl)

      ifld = pbuf_get_index('QRS')
      call pbuf_get_field(pbuf, ifld, qrs)

      do k =1 , pver
         do i = 1, ncol
            qrs(i,k) = qrs(i,k)/state%pdel(i,k)
            qrl(i,k) = qrl(i,k)/state%pdel(i,k)
         end do
      end do

      call outfld('SPQRL   ',qrl/cpair      ,pcols   ,lchnk   )
      call outfld('SPQRS   ',qrs/cpair      ,pcols   ,lchnk   )

      ! The radiation tendencies in the GCM levels above the CRM and the top 2 CRM levels are set to
      ! be zero in the CRM, So add radiation tendencies to these levels 
      ptend%s(:ncol, :pver-crm_nz+2) = qrs(:ncol,:pver-crm_nz+2) + qrl(:ncol,:pver-crm_nz+2)

      !!! This will be used to check energy conservation
      sp_rad_flux(:ncol) = 0.0_r8
      do k=1, pver
         do i=1, ncol
            sp_rad_flux(i) = sp_rad_flux(i) + ( qrs(i,k) + qrl(i,k) ) * state%pdel(i,k)/gravit
         end do
      end do

      !!! Subtract radiative heating for SPDT output
      ftem(:ncol,:pver) = ( ptend%s(:ncol,:pver) - qrs(:ncol,:pver) - qrl(:ncol,:pver) )/cpair

      !---------------------------------------------------------------------------------------------
      ! Output SP-specific fields
      !---------------------------------------------------------------------------------------------
      call outfld('SPDQ    ',ptend%q(1,1,1)        ,pcols ,lchnk )
      call outfld('SPDQC   ',ptend%q(1,1,ixcldliq) ,pcols ,lchnk )
      call outfld('SPDQI   ',ptend%q(1,1,ixcldice) ,pcols ,lchnk )
      call outfld('SPDT    ',ftem   ,pcols ,lchnk )
      call outfld('SPMC    ',crm_output%mctot  ,pcols ,lchnk )
      call outfld('SPMCUP  ',crm_output%mcup   ,pcols ,lchnk )
      call outfld('SPMCDN  ',crm_output%mcdn   ,pcols ,lchnk )
      call outfld('SPMCUUP ',crm_output%mcuup  ,pcols ,lchnk )
      call outfld('SPMCUDN ',crm_output%mcudn  ,pcols ,lchnk )
      call outfld('SPQC    ',crm_output%qc_mean   ,pcols ,lchnk )
      call outfld('SPQI    ',crm_output%qi_mean   ,pcols ,lchnk )
      call outfld('SPQS    ',crm_output%qs_mean   ,pcols ,lchnk )
      call outfld('SPQG    ',crm_output%qg_mean   ,pcols ,lchnk )
      call outfld('SPQR    ',crm_output%qr_mean   ,pcols ,lchnk )

#ifdef m2005
      if (SPCAM_microp_scheme .eq. 'm2005') then
         call outfld('SPNC    ',crm_output%nc_mean         ,pcols ,lchnk )
         call outfld('SPNI    ',crm_output%ni_mean         ,pcols ,lchnk )
         call outfld('SPNS    ',crm_output%ns_mean         ,pcols ,lchnk )
         call outfld('SPNG    ',crm_output%ng_mean         ,pcols ,lchnk )
         call outfld('SPNR    ',crm_output%nr_mean         ,pcols ,lchnk )
      endif
#endif

      call outfld('SPQTFLX ',crm_output%flux_qt        ,pcols ,lchnk )
      call outfld('SPUFLX  ',crm_output%flux_u         ,pcols ,lchnk )
      call outfld('SPVFLX  ',crm_output%flux_v         ,pcols ,lchnk )
      call outfld('SPTKE   ',crm_output%tkez           ,pcols ,lchnk )
      call outfld('SPTKES  ',crm_output%tkesgsz        ,pcols ,lchnk )
      call outfld('SPTK    ',crm_output%tkz            ,pcols ,lchnk )
      call outfld('SPQTFLXS',crm_output%fluxsgs_qt     ,pcols ,lchnk )
      call outfld('SPQPFLX ',crm_output%flux_qp        ,pcols ,lchnk )
      call outfld('SPPFLX  ',crm_output%precflux       ,pcols ,lchnk )
      call outfld('SPQTLS  ',crm_output%qt_ls          ,pcols ,lchnk )
      call outfld('SPQTTR  ',crm_output%qt_trans       ,pcols ,lchnk )
      call outfld('SPQPTR  ',crm_output%qp_trans       ,pcols ,lchnk )
      call outfld('SPQPEVP ',crm_output%qp_evp         ,pcols ,lchnk )
      call outfld('SPQPFALL',crm_output%qp_fall        ,pcols ,lchnk )
      call outfld('SPQPSRC ',crm_output%qp_src         ,pcols ,lchnk )
      call outfld('SPTLS   ',crm_output%t_ls           ,pcols ,lchnk )

      ! NOTE: these should overwrite cloud outputs from CAM routines
      call outfld('CLOUD   ',crm_output%cld,  pcols,lchnk)
      call outfld('CLDTOT  ',crm_output%cltot  ,pcols,lchnk)
      call outfld('CLDHGH  ',crm_output%clhgh  ,pcols,lchnk)
      call outfld('CLDMED  ',crm_output%clmed  ,pcols,lchnk)
      call outfld('CLDLOW  ',crm_output%cllow  ,pcols,lchnk)
      call outfld('CLOUDTOP',crm_output%cldtop, pcols,lchnk)

      call outfld('TIMINGF ',crm_output%timing_factor  ,pcols,lchnk)

      !---------------------------------------------------------------------------------------------
      ! Compute liquid water paths (for diagnostics only)
      !---------------------------------------------------------------------------------------------
       tgicewp(:ncol) = 0.
       tgliqwp(:ncol) = 0.
       do k=1,pver
          do i = 1,ncol
             cicewp(i,k) = crm_output%gicewp(i,k) * 1.0e-3 / max(0.01_r8,crm_output%cld(i,k)) ! In-cloud ice water path.  g/m2 --> kg/m2
             cliqwp(i,k) = crm_output%gliqwp(i,k) * 1.0e-3 / max(0.01_r8,crm_output%cld(i,k)) ! In-cloud liquid water path. g/m2 --> kg/m2
             tgicewp(i)  = tgicewp(i) + crm_output%gicewp(i,k) *1.0e-3 ! grid cell mean ice water path.  g/m2 --> kg/m2
             tgliqwp(i)  = tgliqwp(i) + crm_output%gliqwp(i,k) *1.0e-3 ! grid cell mean ice water path.  g/m2 --> kg/m2
          end do
       end do
       tgwp(:ncol) = tgicewp(:ncol) + tgliqwp(:ncol)
       gwp(:ncol,:pver) = crm_output%gicewp(:ncol,:pver) + crm_output%gliqwp(:ncol,:pver)
       cwp(:ncol,:pver) = cicewp(:ncol,:pver) + cliqwp(:ncol,:pver)

       call outfld('GCLDLWP' ,gwp    , pcols,lchnk)
       call outfld('TGCLDCWP',tgwp   , pcols,lchnk)
       call outfld('TGCLDLWP',tgliqwp, pcols,lchnk)
       call outfld('TGCLDIWP',tgicewp, pcols,lchnk)
       call outfld('ICLDTWP' ,cwp    , pcols,lchnk)
       call outfld('ICLDIWP' ,cicewp , pcols,lchnk)

       if (use_ECPP) then

          ! turbulence
          allocate(tempPtr(pcols,pver))
          ifld = pbuf_get_index('TKE_CRM')
          cs(1:ncol,1:pver) = state%pmid(1:ncol,1:pver) / (287.15*state%t(1:ncol,1:pver))
          tempPtr(:ncol, 1:pver) = crm_output%tkez(:ncol, 1:pver)/cs(:ncol, 1:pver)
          call pbuf_set_field(pbuf, ifld, tempPtr, start=(/1,1/), kount=(/pcols, pver/) )
          deallocate(tempPtr)

          ifld = pbuf_get_index('TK_CRM')
          call pbuf_set_field(pbuf, ifld, crm_output%tkz, start=(/1,1/), kount=(/pcols, pver/) )

          ! For convective transport
          do i=1, ncol
           ideep_crm(i) = i*1.0 
          end do
          ifld = pbuf_get_index( 'MU_CRM' )
          call pbuf_set_field(pbuf, ifld, crm_output%mu_crm, start=(/1,1/), kount=(/pcols, pver/) )
          ifld = pbuf_get_index( 'MD_CRM' )
          call pbuf_set_field(pbuf, ifld, crm_output%md_crm, start=(/1,1/), kount=(/pcols, pver/) )
          ifld = pbuf_get_index( 'EU_CRM' )
          call pbuf_set_field(pbuf, ifld, crm_output%eu_crm, start=(/1,1/), kount=(/pcols, pver/) )
          ifld = pbuf_get_index( 'DU_CRM' )
          call pbuf_set_field(pbuf, ifld, crm_output%du_crm, start=(/1,1/), kount=(/pcols, pver/) )
          ifld = pbuf_get_index( 'ED_CRM' )
          call pbuf_set_field(pbuf, ifld, crm_output%eu_crm, start=(/1,1/), kount=(/pcols, pver/) )
          ifld = pbuf_get_index( 'JT_CRM' )
          call pbuf_set_field(pbuf, ifld, crm_output%jt_crm, start=(/1/), kount=(/pcols/) )
          ifld = pbuf_get_index( 'MX_CRM' )
          call pbuf_set_field(pbuf, ifld, crm_output%mx_crm, start=(/1/), kount=(/pcols/) )
          ifld = pbuf_get_index( 'IDEEP_CRM' )
          call pbuf_set_field(pbuf, ifld, ideep_crm, start=(/1/), kount=(/pcols/) )
       endif
       call outfld('MU_CRM  ', crm_output%mu_crm, pcols, lchnk)
       call outfld('MD_CRM  ', crm_output%md_crm, pcols, lchnk)
       call outfld('EU_CRM  ', crm_output%eu_crm, pcols, lchnk)
       call outfld('DU_CRM  ', crm_output%du_crm, pcols, lchnk)
       call outfld('ED_CRM  ', crm_output%ed_crm, pcols, lchnk)

#ifdef ECPP
       if (use_ECPP) then

         crm_ecpp_output%qlsinkcen = crm_ecpp_output%qlsink_avgcen

         call outfld('ACEN    ',      crm_ecpp_output%acen,             pcols, lchnk)
         call outfld('ABND    ',      crm_ecpp_output%abnd,             pcols, lchnk)
         call outfld('ACEN_TF ',      crm_ecpp_output%acen_tf,          pcols, lchnk)
         call outfld('ABND_TF ',      crm_ecpp_output%abnd_tf,          pcols, lchnk)
         call outfld('MASFBND ',      crm_ecpp_output%massflxbnd,       pcols, lchnk)
         call outfld('RHCEN   ',      crm_ecpp_output%rhcen,            pcols, lchnk)
         call outfld('QCCEN   ',      crm_ecpp_output%qcloudcen,        pcols, lchnk)
         call outfld('QICEN   ',      crm_ecpp_output%qicecen,          pcols, lchnk)
         call outfld('QSINK_AFCEN',   crm_ecpp_output%qlsink_afcen,     pcols, lchnk)
         call outfld('PRECRCEN',      crm_ecpp_output%precrcen,         pcols, lchnk)
         call outfld('PRECSCEN',      crm_ecpp_output%precsolidcen,     pcols, lchnk)
         call outfld('WUPTHRES',      crm_ecpp_output%wupthresh_bnd,    pcols, lchnk)
         call outfld('WDNTHRES',      crm_ecpp_output%wdownthresh_bnd,  pcols, lchnk)
         call outfld('WWQUI_CEN',     crm_ecpp_output%wwqui_cen,        pcols, lchnk)
         call outfld('WWQUI_CLD_CEN', crm_ecpp_output%wwqui_cloudy_cen, pcols, lchnk)
         call outfld('WWQUI_BND',     crm_ecpp_output%wwqui_cen,        pcols, lchnk)
         call outfld('WWQUI_CLD_BND', crm_ecpp_output%wwqui_cloudy_cen, pcols, lchnk)
         call outfld('QSINK_BFCEN',   crm_ecpp_output%qlsink_bfcen,     pcols, lchnk)
         call outfld('QSINK_AVGCEN',  crm_ecpp_output%qlsink_avgcen,    pcols, lchnk)
         call outfld('PRAINCEN',      crm_ecpp_output%praincen,         pcols, lchnk)

       end if ! use_ECPP
#endif /* ECPP */

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

      !---------------------------------------------------------------------------------------------
      ! Output CRM momentum tendencies
      !---------------------------------------------------------------------------------------------

#if defined( SP_ESMT )
      call outfld('U_ESMT',crm_output%u_tend_esmt,pcols   ,lchnk   )
      call outfld('V_ESMT',crm_output%v_tend_esmt,pcols   ,lchnk   )
#endif /* SP_ESMT */

#if defined( SP_USE_ESMT )
      ptend%lu = .TRUE.
      ptend%lv = .TRUE.
      ptend%u  = crm_output%u_tend_esmt
      ptend%v  = crm_output%v_tend_esmt
#else /* SP_USE_ESMT not defined */  

#if defined(SPMOMTRANS)
      ptend%lu = .TRUE.
      ptend%lv = .TRUE.
      
      !!! rotate resolved CRM momentum tendencies back
      do i = 1, ncol 
         ptend%u(i) = crm_output%ultend(i) * cos( -1.*crm_angle(i) ) + crm_output%vltend(i) * sin( -1.*crm_angle(i) )
         ptend%v(i) = crm_output%vltend(i) * cos( -1.*crm_angle(i) ) - crm_output%ultend(i) * sin( -1.*crm_angle(i) )
      enddo

      call outfld('UCONVMOM',ptend%u,pcols   ,lchnk   )
      call outfld('VCONVMOM',ptend%v,pcols   ,lchnk   )
#endif /* SPMOMTRANS */

#endif /* SP_USE_ESMT */

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

       call phys_getopts(microp_scheme_out=microp_scheme)
       if(microp_scheme .eq. 'MG' ) then
#ifdef m2005
         if (SPCAM_microp_scheme .eq. 'm2005') then
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
             do k=1, crm_nz 
               m= pver-k+1
               do ii=1, crm_nx
               do jj=1, crm_ny
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
         endif
#endif /* m2005 */
       end if

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
                  if(SPCAM_microp_scheme .eq. 'm2005') then
                     qli_hydro_after(i) = qli_hydro_after(i)+(crm_state%qr(i,ii,jj,m)+ &
                                                              crm_state%qs(i,ii,jj,m)+ &
                                                              crm_state%qg(i,ii,jj,m)) * dp_g
                     qi_hydro_after(i)  =  qi_hydro_after(i)+(crm_state%qs(i,ii,jj,m)+ &
                                                              crm_state%qg(i,ii,jj,m)) * dp_g
                  else if(SPCAM_microp_scheme .eq. 'sam1mom') then 
                     sfactor = max(0._r8,min(1._r8,(crm_state%temperature(i,ii,jj,m)-268.16)*1./(283.16-268.16)))
                     qli_hydro_after(i) = qli_hydro_after(i)+crm_state%qp(i,ii,jj,m) * dp_g
                     qi_hydro_after(i)  =  qi_hydro_after(i)+crm_state%qp(i,ii,jj,m) * (1-sfactor) * dp_g
                  end if ! SPCAM_microp_scheme
               end do ! ii
            end do ! jj
         end do ! m = 1,crm_nz
         qli_hydro_after(i) = qli_hydro_after(i)/(crm_nx*crm_ny)
         qi_hydro_after(i)  =  qi_hydro_after(i)/(crm_nx*crm_ny)
      end do ! i = 1,ncold

      sp_qchk_prec_dp(:ncol) = prec_dp(:ncol) + (qli_hydro_after (:ncol) - &
                                                 qli_hydro_before(:ncol))/crm_run_time/1000._r8
      sp_qchk_snow_dp(:ncol) = snow_dp(:ncol) + ( qi_hydro_after (:ncol) - &
                                                  qi_hydro_before(:ncol))/crm_run_time/1000._r8

   end if ! (is_first_step())


   call t_stopf('crm')
   
   !------------------------------------------------------------------------------------------------
   ! Aerosol stuff
   !------------------------------------------------------------------------------------------------
    
   !!! Calculate aerosol water at CRM domain using water vapor at CRM domain
   !!!call aerosol_wet_intr (state, ptend, ztodt, pbuf, cam_out, dlf)

   !------------------------------------------------------------------------------------------------

   !!! Free memory in derived types
   call crm_state%finalize()
   call crm_rad_finalize(crm_rad)
   call crm_input%finalize()
   call crm_output_finalize(crm_output)

end subroutine crm_physics_tend

!==================================================================================================
!==================================================================================================

subroutine crm_surface_flux_bypass_tend(state, cam_in, ptend)
   !------------------------------------------------------------------------------------------------
   ! This subroutine is used to apply the fluxes when SP_FLUX_BYPASS is used.
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

   integer  :: ii       ! loop iterator
   integer  :: ncol     ! number of columns
   real(r8) :: g_dp     ! temporary variable for unit conversion
   logical, dimension(pcnst) :: lq

   ncol  = state%ncol

   !!! initialize ptend
   lq(:) = .false.
   lq(1) = .true.
   call physics_ptend_init(ptend, state%psetcols, 'SP_FLUX_BYPASS', &
                           lu=.false., lv=.false., ls=.true., lq=lq)

   !!! apply fluxes to bottom layer
   do ii = 1,ncol
      g_dp = gravit * state%rpdel(ii,pver)             ! note : rpdel = 1./pdel
      ptend%s(ii,:)   = 0.
      ptend%q(ii,:,1) = 0.
#ifdef MAML
      ptend%s(ii,pver)   = g_dp * cam_in%shf(ii,1)
#else
      ptend%s(ii,pver)   = g_dp * cam_in%shf(ii)
#endif
      ptend%q(ii,pver,1) = g_dp * cam_in%cflx(ii,1)
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
   logical                           :: use_SPCAM, use_ECPP

   call phys_getopts( use_SPCAM_out = use_SPCAM )
   call phys_getopts( use_ECPP_out  = use_ECPP  )

   lchnk = state%lchnk
   ncol  = state%ncol
   itim  = pbuf_old_tim_idx()

   !!! Save the state and tend variables 
   !!! Overwrite conventional physics effects before calling the crm. 
   !!! Non-CRM physics routines are allowed to compute diagnostic tendencies.
   !!! Note that state_save and tend_save get allocated in the copy routines
   if ( allocated(state_save%t) ) call physics_state_dealloc(state_save)
   if ( allocated(tend_save%dtdt) ) call physics_tend_dealloc(tend_save)

   call physics_state_copy(state,state_save)
   call physics_tend_copy(tend,tend_save)
   
   !!! save the old cloud fraction
   call pbuf_get_field(pbuf, cldo_idx, cldo, start=(/1,1,itim/), kount=(/pcols,pver,1/) )
   cldo_save(:ncol, :) = cldo(:ncol, :)

   !!! other things relevant for aerosols
   if (use_ECPP) then
      ifld = pbuf_get_index('CLD')
      call pbuf_get_field(pbuf, ifld, cld, (/1,1,itim/),(/pcols,pver,1/))
      ifld = pbuf_get_index('ACLDY_CEN')
      call pbuf_get_field(pbuf, ifld, acldy_cen_tbeg)
      if(is_first_step())then
        acldy_cen_tbeg(:ncol,:) = cld(:ncol, :)
      end if
   endif

#if ( defined MODAL_AERO )
   qqcw_all=0_r8
   do i = 1,pcnst
      qqcw   =>  qqcw_get_field(pbuf, i,lchnk, .true.)
      if (associated(qqcw)) qqcw_all(:,:,i) = qqcw(:,:)
   end do

   ifld = pbuf_get_index('DGNUMWET')
   call pbuf_get_field(pbuf, ifld, dgnumwet, start=(/1,1,1/), kount=(/pcols,pver,ntot_amode/))

   if(is_first_step())then
      do i=1,pcnst
         qqcw_all(:,:,i) = 1.e-38_r8
      end do
      dgnumwet(1:pcols,1:pver,1:ntot_amode) = 0.0_r8
   endif

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
   real(r8), pointer, dimension(:,:) :: cld                 ! cloud fraction
   real(r8), pointer, dimension(:,:) :: qqcw                ! 
   character(len=16)                 :: microp_scheme       ! host model microphysics scheme
   logical                           :: use_SPCAM, use_ECPP
   
   call phys_getopts( use_SPCAM_out     = use_SPCAM )
   call phys_getopts( use_ECPP_out      = use_ECPP  )
   call phys_getopts( microp_scheme_out = microp_scheme )

   lchnk = state%lchnk
   ncol  = state%ncol

   ! gas and aerosol species are updated by shallow convective transport.
   ! Aerosol changes from dropmixnuc in cldwat2m.F90 are discarded. 
   ! (i.e. when dropmixnuc is called in cldwat2m.F90, the tendency is set to zero)
   q_aero = state%q

   !!! Restore state and tend (from beginning of tphysbc)
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
   endif

   !!! whannah - not sure why we do this...
   if(is_first_step())then
      do m=1,dyn_time_lvls
         call pbuf_get_field(pbuf, cldo_idx, cldo, start=(/1,1,m/), kount=(/pcols,pver,1/) )
         cldo(:ncol,:) = 0
      enddo
   endif

   !!! Restore cloud fraction
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

   !!! input arguments
   real(r8), intent(in)    :: ql          ! Mean LWC of pixels [ kg/kg ]
   real(r8), intent(in)    :: nl          ! Grid-mean number concentration of cloud liquid droplet [#/kg]
   real(r8), intent(in)    :: qi          ! Mean IWC of pixels [ kg/kg ]
   real(r8), intent(in)    :: ni          ! Grid-mean number concentration of cloud ice    droplet [#/kg]
   real(r8), intent(in)    :: qs          ! mean snow water content [kg/kg]
   real(r8), intent(in)    :: ns          ! Mean snow crystal number concnetration [#/kg]
   real(r8), intent(in)    :: cld         ! Physical stratus fraction
   real(r8), intent(in)    :: pres        ! Air pressure [Pa] 
   real(r8), intent(in)    :: tk          ! air temperature [K]

   !!! output arguments
   real(r8), intent(out)   :: effl        ! Effective radius of cloud liquid droplet [micro-meter]
   real(r8), intent(out)   :: effi        ! Effective radius of cloud ice    droplet [micro-meter]
   real(r8), intent(out)   :: effl_fn     ! effl for fixed number concentration of nlic = 1.e8
   real(r8), intent(out)   :: deffi       ! ice effective diameter for optics (radiation)
   real(r8), intent(out)   :: pgamrad     ! gamma parameter for optics (radiation)
   real(r8), intent(out)   :: lamcrad     ! slope of droplet distribution for optics (radiation)
   real(r8), intent(out)   :: des         ! snow effective diameter for optics (radiation) [micro-meter]

   !!! local variables
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

   !!! for snow water
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
   endif

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
   endif

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
   endif
   !------------------------------------------------------------------------------------------------
   !------------------------------------------------------------------------------------------------
   return
end subroutine m2005_effradius

!==================================================================================================
!==================================================================================================

end module crm_physics
