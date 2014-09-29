   module convect_shallow

   !----------------------------------------------- !
   ! Purpose:                                       !
   !                                                !
   ! CAM interface to the shallow convection scheme !
   !                                                !
   ! Author: D.B. Coleman                           !
   !         Sungsu Park. Jan. 2010.                !
   !                                                !
   !----------------------------------------------- !

   use shr_kind_mod,      only : r8=>shr_kind_r8
   use physconst,         only : cpair, zvir
   use ppgrid,            only : pver, pcols, pverp
   use zm_conv,           only : zm_conv_evap
   use cam_history,       only : outfld, addfld, add_default, phys_decomp
   use cam_logfile,       only : iulog
   use phys_control,      only : phys_getopts

   implicit none
   private                 
   save

   public :: &
             convect_shallow_register,      & ! Register fields in physics buffer
             convect_shallow_init,          & ! Initialize shallow module
             convect_shallow_tend,          & ! Return tendencies
             convect_shallow_use_shfrc	      ! 

   ! The following namelist variable controls which shallow convection package is used.
   !        'Hack'   = Hack shallow convection (default)
   !        'UW'     = UW shallow convection by Sungsu Park and Christopher S. Bretherton
   !        'off'    = No shallow convection

   character(len=16) :: shallow_scheme      ! Default set in phys_control.F90, use namelist to change
   character(len=16) :: microp_scheme       ! Microphysics scheme
   logical           :: history_amwg        ! output the variables used by the AMWG diag package
   logical           :: history_budget      ! Output tendencies and state variables for CAM4 T, qv, ql, qi
   integer           :: history_budget_histfile_num ! output history file number for budget fields

   ! Physics buffer indices 
   integer    ::     icwmrsh_idx    = 0  
   integer    ::      rprdsh_idx    = 0 
   integer    ::     rprdtot_idx    = 0 
   integer    ::      cldtop_idx    = 0 
   integer    ::      cldbot_idx    = 0 
   integer    ::        cush_idx    = 0 
   integer    :: nevapr_shcu_idx    = 0
   integer    ::       shfrc_idx    = 0 
   integer    ::         cld_idx    = 0 
   integer    ::      concld_idx    = 0
   integer    ::      rprddp_idx    = 0
   integer    ::         tke_idx    = 0

   integer    ::         qpert_idx    = 0
   integer    ::         pblh_idx     = 0
   integer    ::      prec_sh_idx     = 0
   integer    ::      snow_sh_idx     = 0

   integer :: & ! field index in physics buffer
      sh_flxprc_idx, &
      sh_flxsnw_idx, &
      sh_cldliq_idx, &
      sh_cldice_idx
   logical :: convproc_do_aer, convproc_do_gas !BSINGH(09/22/2014): Added for unified convective transport
   contains

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine convect_shallow_register

  !-------------------------------------------------- !
  ! Purpose : Register fields with the physics buffer !
  !-------------------------------------------------- !

  use physics_buffer, only : pbuf_add_field, dtype_r8, dyn_time_lvls
  
  implicit none

  call phys_getopts( shallow_scheme_out = shallow_scheme, microp_scheme_out = microp_scheme, &
       convproc_do_aer_out = convproc_do_aer, convproc_do_gas_out=convproc_do_gas) !BSINGH(09/22/2014): Added for unified convective transport
                     

  call pbuf_add_field('ICWMRSH',    'physpkg' ,dtype_r8,(/pcols,pver/),       icwmrsh_idx )
  call pbuf_add_field('RPRDSH',     'physpkg' ,dtype_r8,(/pcols,pver/),       rprdsh_idx )
  call pbuf_add_field('RPRDTOT',    'physpkg' ,dtype_r8,(/pcols,pver/),       rprdtot_idx )
  call pbuf_add_field('CLDTOP',     'physpkg' ,dtype_r8,(/pcols,1/),          cldtop_idx )
  call pbuf_add_field('CLDBOT',     'physpkg' ,dtype_r8,(/pcols,1/),          cldbot_idx )
  call pbuf_add_field('cush',       'global'  ,dtype_r8,(/pcols,dyn_time_lvls/), cush_idx ) 	
  call pbuf_add_field('NEVAPR_SHCU','physpkg' ,dtype_r8,(/pcols,pver/),       nevapr_shcu_idx )
  call pbuf_add_field('PREC_SH',    'physpkg' ,dtype_r8,(/pcols/),            prec_sh_idx )
  call pbuf_add_field('SNOW_SH',    'physpkg' ,dtype_r8,(/pcols/),            snow_sh_idx )

  if( shallow_scheme .eq. 'UW' ) then
      call pbuf_add_field('shfrc','physpkg' ,dtype_r8,(/pcols,pver/),shfrc_idx )
  endif

! shallow interface gbm flux_convective_cloud_rain+snow (kg/m2/s)
  call pbuf_add_field('SH_FLXPRC','physpkg',dtype_r8,(/pcols,pverp/),sh_flxprc_idx)  

! shallow interface gbm flux_convective_cloud_snow (kg/m2/s)
  call pbuf_add_field('SH_FLXSNW','physpkg',dtype_r8,(/pcols,pverp/),sh_flxsnw_idx)  

! shallow gbm cloud liquid water (kg/kg)
  call pbuf_add_field('SH_CLDLIQ','physpkg',dtype_r8,(/pcols,pver/),sh_cldliq_idx)  

! shallow gbm cloud ice water (kg/kg)
  call pbuf_add_field('SH_CLDICE','physpkg',dtype_r8,(/pcols,pver/),sh_cldice_idx)  

  end subroutine convect_shallow_register

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !


  subroutine convect_shallow_init(pref_edge)

  !------------------------------------------------------------------------------- !
  ! Purpose : Declare output fields, and initialize variables needed by convection !
  !------------------------------------------------------------------------------- !

  use cam_history,       only : addfld, add_default, phys_decomp
  use ppgrid,            only : pcols, pver
  use hk_conv,           only : mfinti
  use uwshcu,            only : init_uwshcu
  use physconst,         only : rair, gravit, latvap, rhoh2o, zvir, &
                                cappa, latice, mwdry, mwh2o
  use pmgrid,            only : plev, plevp
  use spmd_utils,        only : masterproc
  use abortutils,        only : endrun
  use phys_control,      only : cam_physpkg_is
  
  use physics_buffer,            only : pbuf_get_index, physics_buffer_desc, pbuf_set_field
   use time_manager,    only :  is_first_step

  implicit none

  real(r8), intent(in)       :: pref_edge(plevp)        ! Reference pressures at interfaces

  integer limcnv                                   ! Top interface level limit for convection
  integer k
  character(len=16)          :: eddy_scheme
    
  ! ------------------------------------------------- !
  ! Variables for detailed abalysis of UW-ShCu scheme !
  ! ------------------------------------------------- !

  call addfld( 'qt_pre_Cu    ', 'kg/kg'   ,  pver ,  'I' , 'qt_preCU'                                         ,  phys_decomp )
  call addfld( 'sl_pre_Cu    ', 'J/kg'    ,  pver ,  'I' , 'sl_preCU'                                         ,  phys_decomp )
  call addfld( 'slv_pre_Cu   ', 'J/kg'    ,  pver ,  'I' , 'slv_preCU'                                        ,  phys_decomp )
  call addfld( 'u_pre_Cu     ', 'm/s'     ,  pver ,  'I' , 'u_preCU'                                          ,  phys_decomp )
  call addfld( 'v_pre_Cu     ', 'm/s'     ,  pver ,  'I' , 'v_preCU'                                          ,  phys_decomp )
  call addfld( 'qv_pre_Cu    ', 'kg/kg'   ,  pver ,  'I' , 'qv_preCU'                                         ,  phys_decomp )
  call addfld( 'ql_pre_Cu    ', 'kg/kg'   ,  pver ,  'I' , 'ql_preCU'                                         ,  phys_decomp )
  call addfld( 'qi_pre_Cu    ', 'kg/kg'   ,  pver ,  'I' , 'qi_preCU'                                         ,  phys_decomp )
  call addfld( 't_pre_Cu     ', 'K'       ,  pver ,  'I' , 't_preCU'                                          ,  phys_decomp )
  call addfld( 'rh_pre_Cu    ', '%'       ,  pver ,  'I' , 'rh_preCU'                                         ,  phys_decomp )

  call addfld( 'qt_aft_Cu    ', 'kg/kg'   ,  pver ,  'I' , 'qt_afterCU'                                       ,  phys_decomp )
  call addfld( 'sl_aft_Cu    ', 'J/kg'    ,  pver ,  'I' , 'sl_afterCU'                                       ,  phys_decomp )
  call addfld( 'slv_aft_Cu   ', 'J/kg'    ,  pver ,  'I' , 'slv_afterCU'                                      ,  phys_decomp )
  call addfld( 'u_aft_Cu     ', 'm/s'     ,  pver ,  'I' , 'u_afterCU'                                        ,  phys_decomp )
  call addfld( 'v_aft_Cu     ', 'm/s'     ,  pver ,  'I' , 'v_afterCU'                                        ,  phys_decomp )
  call addfld( 'qv_aft_Cu    ', 'kg/kg'   ,  pver ,  'I' , 'qv_afterCU'                                       ,  phys_decomp )
  call addfld( 'ql_aft_Cu    ', 'kg/kg'   ,  pver ,  'I' , 'ql_afterCU'                                       ,  phys_decomp )
  call addfld( 'qi_aft_Cu    ', 'kg/kg'   ,  pver ,  'I' , 'qi_afterCU'                                       ,  phys_decomp )
  call addfld( 't_aft_Cu     ', 'K'       ,  pver ,  'I' , 't_afterCU'                                        ,  phys_decomp )
  call addfld( 'rh_aft_Cu    ', '%'       ,  pver ,  'I' , 'rh_afterCU'                                       ,  phys_decomp )

  call addfld( 'tten_Cu      ', 'K/s'     ,  pver ,  'I' , 'Temprtaure tendency by cumulus convection'        ,  phys_decomp )
  call addfld( 'rhten_Cu     ', '%/s'     ,  pver ,  'I' , 'RH tendency by cumumus convection'                ,  phys_decomp )

  ! ------------------------------------------- !
  ! Common Output for Shallow Convection Scheme !
  ! ------------------------------------------- !

  call addfld( 'CMFDT   '     , 'K/s     ',  pver ,  'A' , &
       'T tendency - shallow convection'                           ,  phys_decomp )
  call addfld( 'CMFDQ   '     , 'kg/kg/s ',  pver ,  'A' , &
       'QV tendency - shallow convection'                          ,  phys_decomp )
  call addfld( 'CMFDLIQ '     , 'kg/kg/s ',  pver ,  'A' , &
       'Cloud liq tendency - shallow convection'                   ,  phys_decomp )
  call addfld( 'CMFDICE '     , 'kg/kg/s ',  pver ,  'A' , &
       'Cloud ice tendency - shallow convection'                   ,  phys_decomp )
  call addfld( 'CMFDQR  '     , 'kg/kg/s ',  pver ,  'A' , &
       'Q tendency - shallow convection rainout'                   ,  phys_decomp )
  call addfld( 'EVAPTCM '     , 'K/s     ',  pver ,  'A' , &
       'T tendency - Evaporation/snow prod from Hack convection'   ,  phys_decomp )
  call addfld( 'FZSNTCM '     , 'K/s     ',  pver ,  'A' , &
       'T tendency - Rain to snow conversion from Hack convection' ,  phys_decomp )
  call addfld( 'EVSNTCM '     , 'K/s     ',  pver ,  'A' , &
       'T tendency - Snow to rain prod from Hack convection'       ,  phys_decomp )
  call addfld( 'EVAPQCM '     , 'kg/kg/s ',  pver ,  'A' , &
       'Q tendency - Evaporation from Hack convection'             ,  phys_decomp )
  call addfld( 'QC      '     , 'kg/kg/s ',  pver ,  'A' , &
       'Q tendency - shallow convection LW export'                 ,  phys_decomp )
  call addfld( 'PRECSH  '     , 'm/s     ',  1,      'A' , &
       'Shallow Convection precipitation rate'                     ,  phys_decomp )
  call addfld( 'CMFMC   '     , 'kg/m2/s ',  pverp,  'A' , &
       'Moist shallow convection mass flux'                        ,  phys_decomp )
  call addfld( 'CMFSL   '     , 'W/m2    ',  pverp,  'A' , &
       'Moist shallow convection liquid water static energy flux'  ,  phys_decomp )
  call addfld( 'CMFLQ   '     , 'W/m2    ',  pverp,  'A' , &
       'Moist shallow convection total water flux'                 ,  phys_decomp )
  call addfld( 'CIN     '     , 'J/kg    ',  1    ,  'A' , &
       'Convective inhibition'                                     ,  phys_decomp )
  call addfld( 'CBMF    '     , 'kg/m2/s ',  1    ,  'A' , &
       'Cloud base mass flux'                                      ,  phys_decomp )
  call addfld( 'CLDTOP  '     , '1       ',  1    ,  'I' , &
       'Vertical index of cloud top'                               ,  phys_decomp )
  call addfld( 'CLDBOT  '     , '1       ',  1    ,  'I' , &
       'Vertical index of cloud base'                              ,  phys_decomp )
  call addfld( 'PCLDTOP '     , '1       ',  1    ,  'A' , &
       'Pressure of cloud top'                                     ,  phys_decomp )
  call addfld( 'PCLDBOT '     , '1       ',  1    ,  'A' , &
       'Pressure of cloud base'                                    ,  phys_decomp )

  call addfld( 'FREQSH '      , 'fraction',  1    ,  'A' , &
       'Fractional occurance of shallow convection'                ,  phys_decomp )
                                                                                                                    
  call addfld( 'HKFLXPRC'     , 'kg/m2/s ',  pverp,  'A' , &
       'Flux of precipitation from HK convection'                  ,  phys_decomp )
  call addfld( 'HKFLXSNW'     , 'kg/m2/s ',  pverp,  'A' , &
       'Flux of snow from HK convection'                           ,  phys_decomp )
  call addfld( 'HKNTPRPD'     , 'kg/kg/s ',  pver ,  'A' , &
       'Net precipitation production from HK convection'           ,  phys_decomp )
  call addfld( 'HKNTSNPD'     , 'kg/kg/s ',  pver ,  'A' , &
       'Net snow production from HK convection'                    ,  phys_decomp )
  call addfld( 'HKEIHEAT'     , 'W/kg'    ,  pver ,  'A' , &
       'Heating by ice and evaporation in HK convection'           ,  phys_decomp )

  call addfld ('ICWMRSH  '    , 'kg/kg   ',  pver,   'A' , &
       'Shallow Convection in-cloud water mixing ratio '           ,  phys_decomp )

  if( shallow_scheme .eq. 'UW' ) then
     call addfld( 'UWFLXPRC'     , 'kg/m2/s ',  pverp,  'A' , &
          'Flux of precipitation from UW shallow convection'          ,  phys_decomp )
     call addfld( 'UWFLXSNW'     , 'kg/m2/s ',  pverp,  'A' , &
          'Flux of snow from UW shallow convection'                   ,  phys_decomp )
  end if



  call phys_getopts( eddy_scheme_out = eddy_scheme      , &
                     history_amwg_out = history_amwg    , &
                     history_budget_out = history_budget, &
                     history_budget_histfile_num_out = history_budget_histfile_num)


  if( history_budget ) then
      call add_default( 'CMFDLIQ  ', history_budget_histfile_num, ' ' )
      call add_default( 'CMFDICE  ', history_budget_histfile_num, ' ' )
      call add_default( 'CMFDT   ', history_budget_histfile_num, ' ' )
      call add_default( 'CMFDQ   ', history_budget_histfile_num, ' ' )
      if( cam_physpkg_is('cam3') .or. cam_physpkg_is('cam4') ) then
         call add_default( 'EVAPQCM  ', history_budget_histfile_num, ' ' )
         call add_default( 'EVAPTCM  ', history_budget_histfile_num, ' ' )
      end if
  end if
  pblh_idx  = pbuf_get_index('pblh')


  select case (shallow_scheme)

  case('off')  ! None

     if( masterproc ) write(iulog,*) 'convect_shallow_init: shallow convection OFF'
     continue

  case('Hack') ! Hack scheme

     qpert_idx = pbuf_get_index('qpert')

     if( masterproc ) write(iulog,*) 'convect_shallow_init: Hack shallow convection'
   ! Limit shallow convection to regions below 40 mb
   ! Note this calculation is repeated in the deep convection interface
     if( pref_edge(1) >= 4.e3_r8 ) then
         limcnv = 1
     else
         do k = 1, plev
            if( pref_edge(k) < 4.e3_r8 .and. pref_edge(k+1) >= 4.e3_r8 ) then
                limcnv = k
                goto 10
            end if
         end do
         limcnv = plevp
     end if
10   continue

     if( masterproc ) then
         write(iulog,*) 'MFINTI: Convection will be capped at intfc ', limcnv, ' which is ', pref_edge(limcnv), ' pascals'
     end if
     
     call mfinti( rair, cpair, gravit, latvap, rhoh2o, limcnv) ! Get args from inti.F90

  case('UW') ! Park and Bretherton shallow convection scheme

     if( masterproc ) write(iulog,*) 'convect_shallow_init: UW shallow convection scheme (McCaa)'
     if( eddy_scheme .ne. 'diag_TKE' ) then
         write(iulog,*) 'ERROR: shallow convection scheme ', shallow_scheme, ' is incompatible with eddy scheme ', eddy_scheme
         call endrun( 'convect_shallow_init: shallow_scheme and eddy_scheme are incompatible' )
     endif
     call init_uwshcu( r8, latvap, cpair, latice, zvir, rair, gravit, mwh2o/mwdry )

     tke_idx = pbuf_get_index('tke')


  end select

  cld_idx      = pbuf_get_index('CLD')
  concld_idx   = pbuf_get_index('CONCLD')
  rprddp_idx   = pbuf_get_index('RPRDDP')
  rprdsh_idx   = pbuf_get_index('RPRDSH')
  rprdtot_idx  = pbuf_get_index('RPRDTOT')

  end subroutine convect_shallow_init

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  function convect_shallow_use_shfrc()
  !----------------------------------------------------------------------- !
  ! Purpose: Return true cloud fraction should use shallow convection      !
  !          calculated convective clouds.                                 !
  !   convect_shallow_use_shfrc() = .true.   for     shallow_scheme = 'UW' !
  !   convect_shallow_use_shfrc() = .false.  for     all other schemes     !
  !                                                                        !
  ! Author: D. B. Coleman                                                  !
  !----------------------------------------------------------------------- !
     implicit none
     logical :: convect_shallow_use_shfrc     ! Return value

     if ( shallow_scheme .eq. 'UW' ) then
          convect_shallow_use_shfrc = .true.
     else
	  convect_shallow_use_shfrc = .false.
     endif

     return

  end function convect_shallow_use_shfrc

  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine convect_shallow_tend( ztodt  , cmfmc   , cmfmc2   , &
                                   qc     , qc2     , rliq     , rliq2    , & 
                                   state  , ptend_all, pbuf    , sh_e_ed_ratio) !BSINGH(09/22/2014): Added sh_e_ed_ratio for unified convective transport
   use physics_buffer,  only : physics_buffer_desc, pbuf_get_field, pbuf_set_field, pbuf_old_tim_idx
   use cam_history,     only : outfld
   use physics_types,   only : physics_state, physics_ptend
   use physics_types,   only : physics_ptend_init, physics_update
   use physics_types,   only : physics_state_copy, physics_state_dealloc
   use physics_types,   only : physics_ptend_dealloc
   use physics_types,   only : physics_ptend_sum
   
   use constituents,    only : pcnst, cnst_get_ind, cnst_get_type_byind
   use hk_conv,         only : cmfmca
   use uwshcu,          only : compute_uwshcu_inv
   use time_manager,    only : get_nstep, is_first_step
   use wv_saturation,   only : qsat
   use physconst,       only : latice, latvap, rhoh2o

   use spmd_utils, only : iam
   implicit none

   ! ---------------------- !
   ! Input-Output Arguments !
   ! ---------------------- !
   type(physics_buffer_desc), pointer :: pbuf(:)
   type(physics_state), intent(in)    :: state                           ! Physics state variables
   real(r8),            intent(in)    :: ztodt                           ! 2 delta-t  [ s ]

   type(physics_ptend), intent(out)   :: ptend_all                       ! Indivdual parameterization tendencies
   real(r8),            intent(out)   :: cmfmc2(pcols,pverp)             ! Updraft mass flux by shallow convection [ kg/s/m2 ]
   real(r8),            intent(out)   :: rliq2(pcols)                    ! Vertically-integrated reserved cloud condensate [ m/s ]
   real(r8),            intent(out)   :: qc2(pcols,pver)                 ! Same as qc but only from shallow convection scheme
   real(r8),            intent(out)   :: sh_e_ed_ratio(pcols,pver)       ! fer/(fer+fdr) from uwschu  !RCE !BSINGH(09/22/2014): Added for unified convective transport
   

   real(r8),            intent(inout) :: cmfmc(pcols,pverp)              ! Moist deep + shallow convection cloud mass flux [ kg/s/m2 ]
   real(r8),            intent(inout) :: qc(pcols,pver)                  ! dq/dt due to export of cloud water into environment by shallow and deep convection [ kg/kg/s ]
   real(r8),            intent(inout) :: rliq(pcols)                     ! Vertical integral of qc [ m/s ]


   ! --------------- !
   ! Local Variables ! 
   ! --------------- !
   integer  :: i, k, m
   integer  :: n, x
   integer  :: ilon                                                      ! Global longitude index of a column
   integer  :: ilat                                                      ! Global latitude  index of a column
   integer  :: lchnk                                                     ! Chunk identifier
   integer  :: ncol                                                      ! Number of atmospheric columns
   integer  :: nstep                                                     ! Current time step index
   integer  :: ixcldice, ixcldliq                                        ! Constituent indices for cloud liquid and ice water.
   integer  :: ixnumice, ixnumliq                                        ! Constituent indices for cloud liquid and ice number concentration

   real(r8),  pointer   :: precc(:)                                      ! Shallow convective precipitation (rain+snow) rate at surface [ m/s ]
   real(r8),  pointer   :: snow(:)                                       ! Shallow convective snow rate at surface [ m/s ]

   real(r8) :: ftem(pcols,pver)                                          ! Temporary workspace for outfld variables
   real(r8) :: cnt2(pcols)                                               ! Top level of shallow convective activity
   real(r8) :: cnb2(pcols)                                               ! Bottom level of convective activity
   real(r8) :: tpert(pcols)                                              ! PBL perturbation theta

   real(r8), pointer   :: pblh(:)                                        ! PBL height [ m ]
   real(r8), pointer   :: qpert(:,:)                                     ! PBL perturbation specific humidity

   real(r8) :: ntprprd(pcols,pver)                                       ! Net precip production in layer
   real(r8) :: ntsnprd(pcols,pver)                                       ! Net snow   production in layer
   real(r8) :: tend_s_snwprd(pcols,pver)                                 ! Heating rate of snow production
   real(r8) :: tend_s_snwevmlt(pcols,pver)                               ! Heating rate of evap/melting of snow
   real(r8) :: slflx(pcols,pverp)                                        ! Shallow convective liquid water static energy flux
   real(r8) :: qtflx(pcols,pverp)                                        ! Shallow convective total water flux
   real(r8) :: cmfdqs(pcols, pver)                                       ! Shallow convective snow production
   real(r8) :: zero(pcols)                                               ! Array of zeros
   real(r8) :: cbmf(pcols)                                               ! Shallow cloud base mass flux [ kg/s/m2 ]
   real(r8) :: freqsh(pcols)                                             ! Frequency of shallow convection occurence
   real(r8) :: pcnt(pcols)                                               ! Top    pressure level of shallow + deep convective activity
   real(r8) :: pcnb(pcols)                                               ! Bottom pressure level of shallow + deep convective activity
   real(r8) :: cmfsl(pcols,pverp )                                       ! Convective flux of liquid water static energy
   real(r8) :: cmflq(pcols,pverp )                                       ! Convective flux of total water in energy unit
   
   real(r8) :: ftem_preCu(pcols,pver)                                    ! Saturation vapor pressure after shallow Cu convection
   real(r8) :: tem2(pcols,pver)                                          ! Saturation specific humidity and RH
   real(r8) :: t_preCu(pcols,pver)                                       ! Temperature after shallow Cu convection
   real(r8) :: tten(pcols,pver)                                          ! Temperature tendency after shallow Cu convection
   real(r8) :: rhten(pcols,pver)                                         ! RH tendency after shallow Cu convection
   real(r8) :: iccmr_UW(pcols,pver)                                      ! In-cloud Cumulus LWC+IWC [ kg/m2 ]
   real(r8) :: icwmr_UW(pcols,pver)                                      ! In-cloud Cumulus LWC     [ kg/m2 ]
   real(r8) :: icimr_UW(pcols,pver)                                      ! In-cloud Cumulus IWC     [ kg/m2 ]
   real(r8) :: ptend_tracer(pcols,pver,pcnst)                            ! Tendencies of tracers
   real(r8) :: sum1, sum2, sum3, pdelx 

   real(r8) :: fer_out(pcols,pver), fdr_out(pcols,pver)                  ! fractional entrainment/detrainment rates rom uwschu  !RCE !BSINGH(09/22/2014): Added for unified convective transport
   real(r8), dimension(pcols,pver) :: sl, qt, slv
   real(r8), dimension(pcols,pver) :: sl_preCu, qt_preCu, slv_preCu

   type(physics_state) :: state1                                         ! Locally modify for evaporation to use, not returned
   type(physics_ptend) :: ptend_loc                                      ! Local tendency from processes, added up to return as ptend_all

   integer itim_old, ifld
   real(r8), pointer, dimension(:,:) :: cld
   real(r8), pointer, dimension(:,:) :: concld
   real(r8), pointer, dimension(:,:) :: icwmr                            ! In cloud water + ice mixing ratio
   real(r8), pointer, dimension(:,:) :: rprddp                           ! dq/dt due to deep convective rainout
   real(r8), pointer, dimension(:,:) :: rprdsh                           ! dq/dt due to deep and shallow convective rainout
   real(r8), pointer, dimension(:,:) :: evapcsh                          ! Evaporation of shallow convective precipitation >= 0.
   real(r8), pointer, dimension(:)   :: cnt
   real(r8), pointer, dimension(:)   :: cnb
   real(r8), pointer, dimension(:)   :: cush
   real(r8), pointer, dimension(:,:) :: tke
   real(r8), pointer, dimension(:,:) :: shfrc
   real(r8), pointer, dimension(:,:) :: flxprec                          ! Shallow convective-scale flux of precip (rain+snow) at interfaces [ kg/m2/s ]
   real(r8), pointer, dimension(:,:) :: flxsnow                          ! Shallow convective-scale flux of snow at interfaces [ kg/m2/s ]
   real(r8), pointer, dimension(:,:) :: sh_cldliq
   real(r8), pointer, dimension(:,:) :: sh_cldice

   logical                           :: lq(pcnst)

   ! ----------------------- !
   ! Main Computation Begins ! 
   ! ----------------------- !

   zero  = 0._r8
   nstep = get_nstep()
   lchnk = state%lchnk
   ncol  = state%ncol
  
   call physics_state_copy( state, state1 )          ! Copy state to local state1.

   ! Associate pointers with physics buffer fields


   itim_old   =  pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, cld_idx,         cld,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, concld_idx,      concld, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

   call pbuf_get_field(pbuf, icwmrsh_idx,     icwmr)

   call pbuf_get_field(pbuf, rprddp_idx,      rprddp )

   call pbuf_get_field(pbuf, rprdsh_idx,      rprdsh )

   call pbuf_get_field(pbuf, nevapr_shcu_idx, evapcsh  )

   call pbuf_get_field(pbuf, cldtop_idx,      cnt )

   call pbuf_get_field(pbuf, cldbot_idx,      cnb )

   call pbuf_get_field(pbuf, prec_sh_idx,   precc )

   call pbuf_get_field(pbuf, snow_sh_idx,    snow )

   if( convect_shallow_use_shfrc() ) then ! Park-Bretherton UW Shallow Convection Schemes
       call pbuf_get_field(pbuf, shfrc_idx,  shfrc  )
   endif

   ! Initialization


   call cnst_get_ind( 'CLDLIQ', ixcldliq )
   call cnst_get_ind( 'CLDICE', ixcldice )

   call pbuf_get_field(pbuf, pblh_idx, pblh)

   !  This field probably should reference the pbuf tpert field but it doesnt
   tpert(:ncol)         = 0._r8

   sh_e_ed_ratio(:ncol,:pver) = -1.0_r8  !RCE !BSINGH(09/22/2014): Added for unified convective transport

   select case (shallow_scheme)

   case('off', 'CLUBB_SGS') ! None

      lq(:) = .TRUE.
      call physics_ptend_init( ptend_loc, state%psetcols, 'convect_shallow (off)', ls=.true., lq=lq ) ! Initialize local ptend type

      cmfmc2      = 0._r8
      ptend_loc%q = 0._r8
      ptend_loc%s = 0._r8
      rprdsh      = 0._r8
      cmfdqs      = 0._r8
      precc       = 0._r8
      slflx       = 0._r8
      qtflx       = 0._r8
      icwmr       = 0._r8
      rliq2       = 0._r8
      qc2         = 0._r8
      cmfsl       = 0._r8
      cmflq       = 0._r8
      cnt2        = pver
      cnb2        = 1._r8
      evapcsh     = 0._r8
      snow        = 0._r8

   case('Hack') ! Hack scheme
                                   
      lq(:) = .TRUE.
      call physics_ptend_init( ptend_loc, state%psetcols, 'cmfmca', ls=.true., lq=lq  ) ! Initialize local ptend type

      call pbuf_get_field(pbuf, qpert_idx, qpert)
      qpert(:ncol,2:pcnst) = 0._r8

      call cmfmca( lchnk        ,  ncol         ,                                               &
                   nstep        ,  ztodt        ,  state%pmid ,  state%pdel  ,                  &
                   state%rpdel  ,  state%zm     ,  tpert      ,  qpert       ,  state%phis  ,   &
                   pblh         ,  state%t      ,  state%q    ,  ptend_loc%s ,  ptend_loc%q ,   &
                   cmfmc2       ,  rprdsh       ,  cmfsl      ,  cmflq       ,  precc       ,   &
                   qc2          ,  cnt2         ,  cnb2       ,  icwmr       ,  rliq2       ,   & 
                   state%pmiddry,  state%pdeldry,  state%rpdeldry )

   case('UW')   ! UW shallow convection scheme

      ! -------------------------------------- !
      ! uwshcu does momentum transport as well !
      ! -------------------------------------- !

      ! Initialize local ptend type
      lq(:) = .TRUE.
      call physics_ptend_init( ptend_loc, state%psetcols, 'UWSHCU', ls=.true., lu=.true., lv=.true., lq=lq  ) 

      call pbuf_get_field(pbuf, cush_idx, cush  ,(/1,itim_old/),  (/pcols,1/))
      call pbuf_get_field(pbuf, tke_idx,  tke)


      call pbuf_get_field(pbuf, sh_flxprc_idx, flxprec)
      call pbuf_get_field(pbuf, sh_flxsnw_idx, flxsnow)

      call compute_uwshcu_inv( pcols     , pver    , ncol           , pcnst         , ztodt         ,                   &
                               state%pint, state%zi, state%pmid     , state%zm      , state%pdel    ,                   & 
                               state%u   , state%v , state%q(:,:,1) , state%q(:,:,ixcldliq), state%q(:,:,ixcldice),     &
                               state%t   , state%s , state%q(:,:,:) ,                                                   &
                               tke       , cld     , concld         , pblh          , cush          ,                   &
                               cmfmc2    , slflx   , qtflx          , 							&
			       flxprec, flxsnow, 			         					&
                               ptend_loc%q(:,:,1)  , ptend_loc%q(:,:,ixcldliq), ptend_loc%q(:,:,ixcldice),              &
                               ptend_loc%s         , ptend_loc%u    , ptend_loc%v   , ptend_tracer  ,                   &
                               rprdsh              , cmfdqs         , precc         , snow          ,                   &
                               evapcsh             , shfrc          , iccmr_UW      , icwmr_UW      ,                   &
                               icimr_UW            , cbmf           , qc2           , rliq2         ,                   &
                               cnt2                , cnb2           , lchnk         , state%pdeldry ,                   &
                               fer_out             , fdr_out                                                            )!BSINGH(09/22/2014): Added fer_out and fdr_out for unified convective transport

      if(convproc_do_aer .or. convproc_do_gas) then
         !RCE mods for modal_aero_convproc
         do k = 1, pver
            do i = 1, ncol
               if ( max(fer_out(i,k),fdr_out(i,k)) > 1.0e-10_r8) then
                  sh_e_ed_ratio(i,k) = max(fer_out(i,k),0.0_r8) &
                       / (max(fer_out(i,k),0.0_r8) + max(fdr_out(i,k),0.0_r8))
               end if
            end do
         end do
      endif

      ! --------------------------------------------------------------------- !
      ! Here, 'rprdsh = qrten', 'cmfdqs = qsten' both in unit of [ kg/kg/s ]  !
      ! In addition, define 'icwmr' which includes both liquid and ice.       !
      ! --------------------------------------------------------------------- !

      icwmr(:ncol,:)  = iccmr_UW(:ncol,:) 
      rprdsh(:ncol,:) = rprdsh(:ncol,:) + cmfdqs(:ncol,:)
      do m = 4, pcnst
         ptend_loc%q(:ncol,:pver,m) = ptend_tracer(:ncol,:pver,m)
      enddo

      ! Conservation check
      
      !  do i = 1, ncol
      !  do m = 1, pcnst
      !     sum1 = 0._r8
      !     sum2 = 0._r8
      !     sum3 = 0._r8
      !  do k = 1, pver
      !       if(cnst_get_type_byind(m).eq.'wet') then
      !          pdelx = state%pdel(i,k)
      !       else
      !          pdelx = state%pdeldry(i,k)
      !       endif
      !       sum1 = sum1 + state%q(i,k,m)*pdelx
      !       sum2 = sum2 +(state%q(i,k,m)+ptend_loc%q(i,k,m)*ztodt)*pdelx  
      !       sum3 = sum3 + ptend_loc%q(i,k,m)*pdelx 
      !  enddo
      !  if( m .gt. 3 .and. abs(sum1) .gt. 1.e-13_r8 .and. abs(sum2-sum1)/sum1 .gt. 1.e-12_r8 ) then
      !! if( m .gt. 3 .and. abs(sum3) .gt. 1.e-13_r8 ) then
      !      write(iulog,*) 'Sungsu : convect_shallow.F90 does not conserve tracers : ', m, sum1, sum2, abs(sum2-sum1)/sum1
      !!     write(iulog,*) 'Sungsu : convect_shallow.F90 does not conserve tracers : ', m, sum3
      !  endif
      !  enddo
      !  enddo

      ! ------------------------------------------------- !
      ! Convective fluxes of 'sl' and 'qt' in energy unit !
      ! ------------------------------------------------- !

      cmfsl(:ncol,:pverp) = slflx(:ncol,:pverp)
      cmflq(:ncol,:pverp) = qtflx(:ncol,:pverp) * latvap

      call outfld( 'PRECSH' , precc  , pcols, lchnk )

   end select

   ! --------------------------------------------------------!     
   ! Calculate fractional occurance of shallow convection    !
   ! --------------------------------------------------------!

 ! Modification : I should check whether below computation of freqsh is correct.

   freqsh(:) = 0._r8
   do i = 1, ncol
      if( maxval(cmfmc2(i,:pver)) <= 0._r8 ) then
          freqsh(i) = 1._r8
      end if
   end do

   ! ------------------------------------------------------------------------------ !
   ! Merge shallow convection output with prior results from deep convection scheme !
   ! ------------------------------------------------------------------------------ !

   ! ----------------------------------------------------------------------- !
   ! Combine cumulus updraft mass flux : 'cmfmc2'(shallow) + 'cmfmc'(deep)   !
   ! ----------------------------------------------------------------------- !

   cmfmc(:ncol,:pver) = cmfmc(:ncol,:pver) + cmfmc2(:ncol,:pver)

   ! -------------------------------------------------------------- !
   ! 'cnt2' & 'cnb2' are from shallow, 'cnt' & 'cnb' are from deep  !
   ! 'cnt2' & 'cnb2' are the interface indices of cloud top & base: ! 
   !        cnt2 = float(kpen)                                      !
   !        cnb2 = float(krel - 1)                                  !
   ! Note that indices decreases with height.                       !
   ! -------------------------------------------------------------- !

   do i = 1, ncol
      if( cnt2(i) < cnt(i)) cnt(i) = cnt2(i)
      if( cnb2(i) > cnb(i)) cnb(i) = cnb2(i)
      pcnt(i) = state%pmid(i,int(cnt(i)))
      pcnb(i) = state%pmid(i,int(cnb(i)))     
   end do
   
   ! ----------------------------------------------- !
   ! This quantity was previously known as CMFDQR.   !
   ! Now CMFDQR is the shallow rain production only. !
   ! ----------------------------------------------- !

   
   call pbuf_set_field(pbuf, rprdtot_idx, rprdsh(:ncol,:pver) + rprddp(:ncol,:pver), start=(/1,1/), kount=(/ncol,pver/))
 
   ! ----------------------------------------------------------------------- ! 
   ! Add shallow reserved cloud condensate to deep reserved cloud condensate !
   !     qc [ kg/kg/s] , rliq [ m/s ]                                        !
   ! ----------------------------------------------------------------------- !

   qc(:ncol,:pver) = qc(:ncol,:pver) + qc2(:ncol,:pver)
   rliq(:ncol)     = rliq(:ncol) + rliq2(:ncol)    

   ! ---------------------------------------------------------------------------- !
   ! Output new partition of cloud condensate variables, as well as precipitation !
   ! ---------------------------------------------------------------------------- ! 

   if( microp_scheme == 'MG' ) then
       call cnst_get_ind( 'NUMLIQ', ixnumliq )
       call cnst_get_ind( 'NUMICE', ixnumice )
   endif

   ftem(:ncol,:pver) = ptend_loc%s(:ncol,:pver)/cpair

   call outfld( 'ICWMRSH ', icwmr                    , pcols   , lchnk )

   call outfld( 'CMFDT  ', ftem                      , pcols   , lchnk )
   call outfld( 'CMFDQ  ', ptend_loc%q(1,1,1)        , pcols   , lchnk )
   call outfld( 'CMFDICE', ptend_loc%q(1,1,ixcldice) , pcols   , lchnk )
   call outfld( 'CMFDLIQ', ptend_loc%q(1,1,ixcldliq) , pcols   , lchnk )
   call outfld( 'CMFMC'  , cmfmc                     , pcols   , lchnk )
   call outfld( 'QC'     , qc2                       , pcols   , lchnk )
   call outfld( 'CMFDQR' , rprdsh                    , pcols   , lchnk )
   call outfld( 'CMFSL'  , cmfsl                     , pcols   , lchnk )
   call outfld( 'CMFLQ'  , cmflq                     , pcols   , lchnk )
   call outfld( 'DQP'    , qc2                       , pcols   , lchnk )
   call outfld( 'CBMF'   , cbmf                      , pcols   , lchnk )
   call outfld( 'CLDTOP' , cnt                       , pcols   , lchnk )
   call outfld( 'CLDBOT' , cnb                       , pcols   , lchnk )
   call outfld( 'PCLDTOP', pcnt                      , pcols   , lchnk )
   call outfld( 'PCLDBOT', pcnb                      , pcols   , lchnk )  
   call outfld( 'FREQSH' , freqsh                    , pcols   , lchnk )

   if( shallow_scheme .eq. 'UW' ) then
      call outfld( 'UWFLXPRC', flxprec                  , pcols   , lchnk )  
      call outfld( 'UWFLXSNW' , flxsnow                 , pcols   , lchnk )
   endif

   ! ---------------------------------------------------------------- !
   ! Add tendency from this process to tend from other processes here !
   ! ---------------------------------------------------------------- !

   call physics_ptend_init(ptend_all, state1%psetcols, 'convect_shallow')
   call physics_ptend_sum( ptend_loc, ptend_all, ncol )

   ! ----------------------------------------------------------------------------- !
   ! For diagnostic purpose, print out 'QT,SL,SLV,T,RH' just before cumulus scheme !
   ! ----------------------------------------------------------------------------- !

   sl_preCu(:ncol,:pver)  = state1%s(:ncol,:pver) -   latvap           * state1%q(:ncol,:pver,ixcldliq) &
                                                  - ( latvap + latice) * state1%q(:ncol,:pver,ixcldice)
   qt_preCu(:ncol,:pver)  = state1%q(:ncol,:pver,1) + state1%q(:ncol,:pver,ixcldliq) &
                                                    + state1%q(:ncol,:pver,ixcldice)
   slv_preCu(:ncol,:pver) = sl_preCu(:ncol,:pver) * ( 1._r8 + zvir * qt_preCu(:ncol,:pver) )

   t_preCu(:ncol,:)       = state1%t(:ncol,:pver)
   call qsat(state1%t(:ncol,:), state1%pmid(:ncol,:), &
        tem2(:ncol,:), ftem(:ncol,:))
   ftem_preCu(:ncol,:)    = state1%q(:ncol,:,1) / ftem(:ncol,:) * 100._r8

   call outfld( 'qt_pre_Cu      ', qt_preCu               , pcols, lchnk )
   call outfld( 'sl_pre_Cu      ', sl_preCu               , pcols, lchnk )
   call outfld( 'slv_pre_Cu     ', slv_preCu              , pcols, lchnk )
   call outfld( 'u_pre_Cu       ', state1%u               , pcols, lchnk )
   call outfld( 'v_pre_Cu       ', state1%v               , pcols, lchnk )
   call outfld( 'qv_pre_Cu      ', state1%q(:,:,1)        , pcols, lchnk )
   call outfld( 'ql_pre_Cu      ', state1%q(:,:,ixcldliq) , pcols, lchnk )
   call outfld( 'qi_pre_Cu      ', state1%q(:,:,ixcldice) , pcols, lchnk )
   call outfld( 't_pre_Cu       ', state1%t               , pcols, lchnk )
   call outfld( 'rh_pre_Cu      ', ftem_preCu             , pcols, lchnk )

   ! ----------------------------------------------- ! 
   ! Update physics state type state1 with ptend_loc ! 
   ! ----------------------------------------------- !

   call physics_update( state1, ptend_loc, ztodt )

   ! ----------------------------------------------------------------------------- !
   ! For diagnostic purpose, print out 'QT,SL,SLV,t,RH' just after cumulus scheme  !
   ! ----------------------------------------------------------------------------- !

   sl(:ncol,:pver)  = state1%s(:ncol,:pver) -   latvap           * state1%q(:ncol,:pver,ixcldliq) &
                                            - ( latvap + latice) * state1%q(:ncol,:pver,ixcldice)
   qt(:ncol,:pver)  = state1%q(:ncol,:pver,1) + state1%q(:ncol,:pver,ixcldliq) &
                                              + state1%q(:ncol,:pver,ixcldice)
   slv(:ncol,:pver) = sl(:ncol,:pver) * ( 1._r8 + zvir * qt(:ncol,:pver) )

   call qsat(state1%t(:ncol,:), state1%pmid(:ncol,:), &
        tem2(:ncol,:), ftem(:ncol,:))
   ftem(:ncol,:)    = state1%q(:ncol,:,1) / ftem(:ncol,:) * 100._r8

   call outfld( 'qt_aft_Cu      ', qt                     , pcols, lchnk )
   call outfld( 'sl_aft_Cu      ', sl                     , pcols, lchnk )
   call outfld( 'slv_aft_Cu     ', slv                    , pcols, lchnk )
   call outfld( 'u_aft_Cu       ', state1%u               , pcols, lchnk )
   call outfld( 'v_aft_Cu       ', state1%v               , pcols, lchnk )
   call outfld( 'qv_aft_Cu      ', state1%q(:,:,1)        , pcols, lchnk )
   call outfld( 'ql_aft_Cu      ', state1%q(:,:,ixcldliq) , pcols, lchnk )
   call outfld( 'qi_aft_Cu      ', state1%q(:,:,ixcldice) , pcols, lchnk )
   call outfld( 't_aft_Cu       ', state1%t               , pcols, lchnk )
   call outfld( 'rh_aft_Cu      ', ftem                   , pcols, lchnk )

   tten(:ncol,:)  = ( state1%t(:ncol,:pver) - t_preCu(:ncol,:) ) / ztodt 
   rhten(:ncol,:) = ( ftem(:ncol,:) - ftem_preCu(:ncol,:) ) / ztodt 

   call outfld( 'tten_Cu        ', tten                           , pcols, lchnk )
   call outfld( 'rhten_Cu       ', rhten                          , pcols, lchnk )


   ! ------------------------------------------------------------------------ !
   ! UW-Shallow Cumulus scheme includes                                       !
   ! evaporation physics inside in it. So when 'shallow_scheme = UW', we must !
   ! NOT perform below 'zm_conv_evap'.                                        !
   ! ------------------------------------------------------------------------ !

   if( shallow_scheme .eq. 'Hack' ) then

   ! ------------------------------------------------------------------------------- !
   ! Determine the phase of the precipitation produced and add latent heat of fusion !
   ! Evaporate some of the precip directly into the environment (Sundqvist)          !
   ! Allow this to use the updated state1 and a fresh ptend_loc type                 !
   ! Heating and specific humidity tendencies produced                               !
   ! ------------------------------------------------------------------------------- !

   ! --------------------------------- !
   ! initialize ptend for next process !
   ! --------------------------------- !

    lq(1) = .TRUE.
    lq(2:) = .FALSE.
    call physics_ptend_init(ptend_loc, state1%psetcols, 'zm_conv_evap', ls=.true., lq=lq)

    call pbuf_get_field(pbuf, sh_flxprc_idx, flxprec    )
    call pbuf_get_field(pbuf, sh_flxsnw_idx, flxsnow    )
    call pbuf_get_field(pbuf, sh_cldliq_idx, sh_cldliq  )
    call pbuf_get_field(pbuf, sh_cldice_idx, sh_cldice  )

    !! clouds have no water... :)
    sh_cldliq(:ncol,:) = 0._r8
    sh_cldice(:ncol,:) = 0._r8

    call zm_conv_evap( state1%ncol, state1%lchnk,                                    &
                       state1%t, state1%pmid, state1%pdel, state1%q(:pcols,:pver,1), &
                       ptend_loc%s, tend_s_snwprd, tend_s_snwevmlt,                  & 
                       ptend_loc%q(:pcols,:pver,1),                                  &
                       rprdsh, cld, ztodt,                                           &
                       precc, snow, ntprprd, ntsnprd , flxprec, flxsnow )

   ! ------------------------------------------ !
   ! record history variables from zm_conv_evap !
   ! ------------------------------------------ !

   evapcsh(:ncol,:pver) = ptend_loc%q(:ncol,:pver,1)

   ftem(:ncol,:pver) = ptend_loc%s(:ncol,:pver) / cpair
   call outfld( 'EVAPTCM '       , ftem                           , pcols, lchnk )
   ftem(:ncol,:pver) = tend_s_snwprd(:ncol,:pver) / cpair
   call outfld( 'FZSNTCM '       , ftem                           , pcols, lchnk )
   ftem(:ncol,:pver) = tend_s_snwevmlt(:ncol,:pver) / cpair
   call outfld( 'EVSNTCM '       , ftem                           , pcols, lchnk )
   call outfld( 'EVAPQCM '       , ptend_loc%q(1,1,1)             , pcols, lchnk )
   call outfld( 'PRECSH  '       , precc                          , pcols, lchnk )
   call outfld( 'HKFLXPRC'       , flxprec                        , pcols, lchnk )
   call outfld( 'HKFLXSNW'       , flxsnow                        , pcols, lchnk )
   call outfld( 'HKNTPRPD'       , ntprprd                        , pcols, lchnk )
   call outfld( 'HKNTSNPD'       , ntsnprd                        , pcols, lchnk )
   call outfld( 'HKEIHEAT'       , ptend_loc%s                    , pcols, lchnk )

   ! ---------------------------------------------------------------- !      
   ! Add tendency from this process to tend from other processes here !
   ! ---------------------------------------------------------------- !

   call physics_ptend_sum( ptend_loc, ptend_all, ncol )
   call physics_ptend_dealloc(ptend_loc)

   ! -------------------------------------------- !
   ! Do not perform evaporation process for UW-Cu !
   ! -------------------------------------------- !

   end if

   ! ------------------------------------------------------------- !
   ! Update name of parameterization tendencies to send to tphysbc !
   ! ------------------------------------------------------------- !

   call physics_state_dealloc(state1)

  end subroutine convect_shallow_tend

  end module convect_shallow
