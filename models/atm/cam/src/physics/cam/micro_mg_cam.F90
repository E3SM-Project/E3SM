module micro_mg_cam

!---------------------------------------------------------------------------------
!
!  CAM Interfaces for MG microphysics
!
!---------------------------------------------------------------------------------

use shr_kind_mod,   only: r8=>shr_kind_r8
use spmd_utils,     only: masterproc
use ppgrid,         only: pver, pverp, psubcols
use physconst,      only: gravit, rair, tmelt, cpair, rh2o, rhoh2o, &
     latvap, latice, mwdry
use phys_control,   only: phys_getopts


use physics_types,  only: physics_state, physics_ptend, physics_ptend_init, &
     physics_state_copy, physics_ptend_copy, &
     physics_update, physics_state_dealloc, &
     physics_ptend_sum
use physics_buffer, only: physics_buffer_desc, pbuf_add_field, dyn_time_lvls, &
     pbuf_old_tim_idx, pbuf_get_index, dtype_r8, dtype_i4, &
     pbuf_get_field, pbuf_set_field, col_type_subcol, pbuf_register_subcol
use constituents,   only: cnst_add, cnst_get_ind, &
     cnst_name, cnst_longname, sflxnam, apcnst, bpcnst, pcnst

use cldwat2m_macro, only: rhmini

use cam_history,    only: addfld, add_default, phys_decomp, outfld

use cam_logfile,    only: iulog
use abortutils,     only: endrun
use error_messages, only: handle_errmsg
use ref_pres,       only: top_lev=>trop_cloud_top_lev

use subcol_utils,   only: subcol_get_scheme

implicit none
private
save

public :: &
     micro_mg_cam_readnl,          &
     micro_mg_cam_register,        &
     micro_mg_cam_init_cnst,       &
     micro_mg_cam_implements_cnst, &
     micro_mg_cam_init,            &
     micro_mg_cam_tend

integer :: micro_mg_version     = 1      ! Version number for MG.
integer :: micro_mg_sub_version = 0      ! Second part of version number.

logical :: microp_uniform

logical, public :: do_cldliq ! Prognose cldliq flag
logical, public :: do_cldice ! Prognose cldice flag

integer, parameter :: ncnst = 4       ! Number of constituents
character(len=8), parameter :: &      ! Constituent names
     cnst_names(ncnst) = (/'CLDLIQ', 'CLDICE','NUMLIQ','NUMICE'/)

integer :: &
     ixcldliq,      &! cloud liquid amount index
     ixcldice,      &! cloud ice amount index
     ixnumliq,      &! cloud liquid number index
     ixnumice        ! cloud ice water index

! Physics buffer indices for fields registered by this module
integer :: &
     cldo_idx,           &
     qme_idx,            &
     prain_idx,          &
     nevapr_idx,         &
     wsedl_idx,          &
     rei_idx,            &
     rel_idx,            &
     dei_idx,            &
     mu_idx,             &
     lambdac_idx,        &
     iciwpst_idx,        &
     iclwpst_idx,        &
     des_idx,            &
     icswp_idx,          &
     cldfsnow_idx,       &
     rate1_cw2pr_st_idx = -1, &
     ls_flxprc_idx,      &
     ls_flxsnw_idx,      &
     relvar_idx,         &
     cmeliq_idx,         &
     accre_enhan_idx

! Fields needed as inputs to COSP
integer :: &
     ls_mrprc_idx,    ls_mrsnw_idx,    &
     ls_reffrain_idx, ls_reffsnow_idx, &
     cv_reffliq_idx,  cv_reffice_idx

! Fields needed by Park macrophysics
integer :: &
     cc_t_idx,  cc_qv_idx, &
     cc_ql_idx, cc_qi_idx, &
     cc_nl_idx, cc_ni_idx, &
     cc_qlst_idx

! Used to replace aspects of MG microphysics
! (e.g. by CARMA)
integer :: tnd_qsnow_idx, tnd_nsnow_idx, re_ice_idx

! Index fields for precipitation efficiency.
integer :: acpr_idx, acgcme_idx, acnum_idx

! Physics buffer indices for fields registered by other modules
integer :: &
     ast_idx = -1,            &
     aist_idx = -1,           &
     alst_idx = -1,           &
     cld_idx = -1,            &
     concld_idx = -1

! Pbuf fields needed for subcol_SILHS
integer :: &
   qrain_idx=-1, qsnow_idx=-1,    &
   nrain_idx=-1, nsnow_idx=-1

integer :: &
     naai_idx = -1,           &
     naai_hom_idx = -1,       &
     npccn_idx = -1,          &
     rndst_idx = -1,          &
     nacon_idx = -1,          &
     prec_str_idx = -1,       &
     snow_str_idx = -1,       &
     prec_pcw_idx = -1,       &
     snow_pcw_idx = -1,       &
     prec_sed_idx = -1,       &
     snow_sed_idx = -1


!===============================================================================
contains
!===============================================================================

subroutine micro_mg_cam_readnl(nlfile)

  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use mpishorthand

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  ! Namelist variables
  logical :: micro_mg_do_cldice   = .true. ! do_cldice = .true., MG microphysics is prognosing cldice
  logical :: micro_mg_do_cldliq   = .true. ! do_cldliq = .true., MG microphysics is prognosing cldliq

  ! Local variables
  integer :: unitn, ierr
  character(len=*), parameter :: subname = 'micro_mg_cam_readnl'

  namelist /micro_mg_nl/ micro_mg_version, micro_mg_sub_version, &
       micro_mg_do_cldice, micro_mg_do_cldliq, microp_uniform

  !-----------------------------------------------------------------------------

  if (masterproc) then
     unitn = getunit()
     open( unitn, file=trim(nlfile), status='old' )
     call find_group_name(unitn, 'micro_mg_nl', status=ierr)
     if (ierr == 0) then
        read(unitn, micro_mg_nl, iostat=ierr)
        if (ierr /= 0) then
           call endrun(subname // ':: ERROR reading namelist')
        end if
     end if
     close(unitn)
     call freeunit(unitn)

     ! set local variables
     do_cldice  = micro_mg_do_cldice
     do_cldliq  = micro_mg_do_cldliq

     ! Verify that version numbers are valid.
     select case (micro_mg_version)
     case (1)
        select case (micro_mg_sub_version)
        case(0)
           ! MG version 1.0
        case(5)
           ! MG version 1.5 - MG2 development
        case default
           call bad_version_endrun()
        end select
     case default
        call bad_version_endrun()
     end select

  end if

#ifdef SPMD
  ! Broadcast namelist variables
  call mpibcast(micro_mg_version,     1, mpiint, 0, mpicom)
  call mpibcast(micro_mg_sub_version, 1, mpiint, 0, mpicom)
  call mpibcast(do_cldice,            1, mpilog, 0, mpicom)
  call mpibcast(do_cldliq,            1, mpilog, 0, mpicom)
  call mpibcast(microp_uniform,       1, mpilog, 0, mpicom)
#endif

contains

  subroutine bad_version_endrun
    ! Endrun wrapper with a more useful error message.
    character(len=128) :: errstring
    write(errstring,*) "Invalid version number specified for MG microphysics: ", &
         micro_mg_version,".",micro_mg_sub_version
    call endrun(errstring)
  end subroutine bad_version_endrun

end subroutine micro_mg_cam_readnl

!================================================================================================

subroutine micro_mg_cam_register

  ! Register microphysics constituents and fields in the physics buffer.
  !-----------------------------------------------------------------------

  use ppgrid,          only: pcols

  logical :: prog_modal_aero
  logical :: use_subcol_microp  ! If true, then are using subcolumns in microphysics
  logical :: save_subcol_microp ! If true, then need to store sub-columnized fields in pbuf

  call phys_getopts(use_subcol_microp_out = use_subcol_microp, &
                    prog_modal_aero_out   = prog_modal_aero )

  ! Register microphysics constituents and save indices.

  call cnst_add(cnst_names(1), mwdry, cpair, 0._r8, ixcldliq, &
       longname='Grid box averaged cloud liquid amount', is_convtran1=.true.)
  call cnst_add(cnst_names(2), mwdry, cpair, 0._r8, ixcldice, &
       longname='Grid box averaged cloud ice amount', is_convtran1=.true.)
  ! The next statements should have "is_convtran1=.true.", but this would change
  ! answers.
  call cnst_add(cnst_names(3), mwdry, cpair, 0._r8, ixnumliq, &
       longname='Grid box averaged cloud liquid number', is_convtran1=.false.)
  call cnst_add(cnst_names(4), mwdry, cpair, 0._r8, ixnumice, &
       longname='Grid box averaged cloud ice number', is_convtran1=.false.)

  ! Request physics buffer space for fields that persist across timesteps.

  call pbuf_add_field('CLDO','global',dtype_r8,(/pcols,pver,dyn_time_lvls/), cldo_idx)

  ! Physics buffer variables for convective cloud properties.

  call pbuf_add_field('QME',        'physpkg',dtype_r8,(/pcols,pver/), qme_idx)
  call pbuf_add_field('PRAIN',      'physpkg',dtype_r8,(/pcols,pver/), prain_idx)
  call pbuf_add_field('NEVAPR',     'physpkg',dtype_r8,(/pcols,pver/), nevapr_idx)

  call pbuf_add_field('WSEDL',      'physpkg',dtype_r8,(/pcols,pver/), wsedl_idx)

  call pbuf_add_field('REI',        'physpkg',dtype_r8,(/pcols,pver/), rei_idx)
  call pbuf_add_field('REL',        'physpkg',dtype_r8,(/pcols,pver/), rel_idx)

  ! Mitchell ice effective diameter for radiation
  call pbuf_add_field('DEI',        'physpkg',dtype_r8,(/pcols,pver/), dei_idx)
  ! Size distribution shape parameter for radiation
  call pbuf_add_field('MU',         'physpkg',dtype_r8,(/pcols,pver/), mu_idx)
  ! Size distribution shape parameter for radiation
  call pbuf_add_field('LAMBDAC',    'physpkg',dtype_r8,(/pcols,pver/), lambdac_idx)

  ! Stratiform only in cloud ice water path for radiation
  call pbuf_add_field('ICIWPST',    'physpkg',dtype_r8,(/pcols,pver/), iciwpst_idx)
  ! Stratiform in cloud liquid water path for radiation
  call pbuf_add_field('ICLWPST',    'physpkg',dtype_r8,(/pcols,pver/), iclwpst_idx)

  ! Snow effective diameter for radiation
  call pbuf_add_field('DES',        'physpkg',dtype_r8,(/pcols,pver/), des_idx)
  ! In cloud snow water path for radiation
  call pbuf_add_field('ICSWP',      'physpkg',dtype_r8,(/pcols,pver/), icswp_idx)
  ! Cloud fraction for liquid drops + snow
  call pbuf_add_field('CLDFSNOW ',  'physpkg',dtype_r8,(/pcols,pver,dyn_time_lvls/), cldfsnow_idx)

  if (prog_modal_aero) then
     call pbuf_add_field('RATE1_CW2PR_ST','physpkg',dtype_r8,(/pcols,pver/), rate1_cw2pr_st_idx)
  endif

  call pbuf_add_field('LS_FLXPRC',  'physpkg',dtype_r8,(/pcols,pverp/), ls_flxprc_idx)
  call pbuf_add_field('LS_FLXSNW',  'physpkg',dtype_r8,(/pcols,pverp/), ls_flxsnw_idx)


  ! Fields needed as inputs to COSP
  call pbuf_add_field('LS_MRPRC',   'physpkg',dtype_r8,(/pcols,pver/), ls_mrprc_idx)
  call pbuf_add_field('LS_MRSNW',   'physpkg',dtype_r8,(/pcols,pver/), ls_mrsnw_idx)
  call pbuf_add_field('LS_REFFRAIN','physpkg',dtype_r8,(/pcols,pver/), ls_reffrain_idx)
  call pbuf_add_field('LS_REFFSNOW','physpkg',dtype_r8,(/pcols,pver/), ls_reffsnow_idx)
  call pbuf_add_field('CV_REFFLIQ', 'physpkg',dtype_r8,(/pcols,pver/), cv_reffliq_idx)
  call pbuf_add_field('CV_REFFICE', 'physpkg',dtype_r8,(/pcols,pver/), cv_reffice_idx)

  ! CC_* Fields needed by Park macrophysics
  call pbuf_add_field('CC_T',     'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_t_idx)
  call pbuf_add_field('CC_qv',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_qv_idx)
  call pbuf_add_field('CC_ql',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_ql_idx)
  call pbuf_add_field('CC_qi',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_qi_idx)
  call pbuf_add_field('CC_nl',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_nl_idx)
  call pbuf_add_field('CC_ni',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_ni_idx)
  call pbuf_add_field('CC_qlst',  'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_qlst_idx)

  ! Register subcolumn pbuf fields
  if (use_subcol_microp) then
    ! Global pbuf fields
    call pbuf_register_subcol('CLDO',        'micro_mg_cam_register', cldo_idx)

    ! CC_* Fields needed by Park macrophysics
    call pbuf_register_subcol('CC_T',        'micro_mg_cam_register', cc_t_idx)
    call pbuf_register_subcol('CC_qv',       'micro_mg_cam_register', cc_qv_idx)
    call pbuf_register_subcol('CC_ql',       'micro_mg_cam_register', cc_ql_idx)
    call pbuf_register_subcol('CC_qi',       'micro_mg_cam_register', cc_qi_idx)
    call pbuf_register_subcol('CC_nl',       'micro_mg_cam_register', cc_nl_idx)
    call pbuf_register_subcol('CC_ni',       'micro_mg_cam_register', cc_ni_idx)
    call pbuf_register_subcol('CC_qlst',     'micro_mg_cam_register', cc_qlst_idx)

    ! Physpkg pbuf fields
    ! Physics buffer variables for convective cloud properties.

    call pbuf_register_subcol('QME',         'micro_mg_cam_register', qme_idx)
    call pbuf_register_subcol('PRAIN',       'micro_mg_cam_register', prain_idx)
    call pbuf_register_subcol('NEVAPR',      'micro_mg_cam_register', nevapr_idx)

    call pbuf_register_subcol('WSEDL',       'micro_mg_cam_register', wsedl_idx)

    call pbuf_register_subcol('REI',         'micro_mg_cam_register', rei_idx)
    call pbuf_register_subcol('REL',         'micro_mg_cam_register', rel_idx)

    ! Mitchell ice effective diameter for radiation
    call pbuf_register_subcol('DEI',         'micro_mg_cam_register', dei_idx)
    ! Size distribution shape parameter for radiation
    call pbuf_register_subcol('MU',          'micro_mg_cam_register', mu_idx)
    ! Size distribution shape parameter for radiation
    call pbuf_register_subcol('LAMBDAC',     'micro_mg_cam_register', lambdac_idx)

    ! Stratiform only in cloud ice water path for radiation
    call pbuf_register_subcol('ICIWPST',     'micro_mg_cam_register', iciwpst_idx)
    ! Stratiform in cloud liquid water path for radiation
    call pbuf_register_subcol('ICLWPST',     'micro_mg_cam_register', iclwpst_idx)

    ! Snow effective diameter for radiation
    call pbuf_register_subcol('DES',         'micro_mg_cam_register', des_idx)
    ! In cloud snow water path for radiation
    call pbuf_register_subcol('ICSWP',       'micro_mg_cam_register', icswp_idx)
    ! Cloud fraction for liquid drops + snow
    call pbuf_register_subcol('CLDFSNOW ',   'micro_mg_cam_register', cldfsnow_idx)

    if (prog_modal_aero) then
      call pbuf_register_subcol('RATE1_CW2PR_ST', 'micro_mg_cam_register', rate1_cw2pr_st_idx)
    end if

    call pbuf_register_subcol('LS_FLXPRC',   'micro_mg_cam_register', ls_flxprc_idx)
    call pbuf_register_subcol('LS_FLXSNW',   'micro_mg_cam_register', ls_flxsnw_idx)

    ! Fields needed as inputs to COSP
    call pbuf_register_subcol('LS_MRPRC',    'micro_mg_cam_register', ls_mrprc_idx)
    call pbuf_register_subcol('LS_MRSNW',    'micro_mg_cam_register', ls_mrsnw_idx)
    call pbuf_register_subcol('LS_REFFRAIN', 'micro_mg_cam_register', ls_reffrain_idx)
    call pbuf_register_subcol('LS_REFFSNOW', 'micro_mg_cam_register', ls_reffsnow_idx)
    call pbuf_register_subcol('CV_REFFLIQ',  'micro_mg_cam_register', cv_reffliq_idx)
    call pbuf_register_subcol('CV_REFFICE',  'micro_mg_cam_register', cv_reffice_idx)
  end if

  ! Additional pbuf for CARMA interface
  call pbuf_add_field('TND_QSNOW',  'physpkg',dtype_r8,(/pcols,pver/), tnd_qsnow_idx)
  call pbuf_add_field('TND_NSNOW',  'physpkg',dtype_r8,(/pcols,pver/), tnd_nsnow_idx)
  call pbuf_add_field('RE_ICE',     'physpkg',dtype_r8,(/pcols,pver/), re_ice_idx)

  ! Precipitation efficiency fields across timesteps.
  call pbuf_add_field('ACPRECL',    'global',dtype_r8,(/pcols/), acpr_idx)   ! accumulated precip
  call pbuf_add_field('ACGCME',     'global',dtype_r8,(/pcols/), acgcme_idx) ! accumulated condensation
  call pbuf_add_field('ACNUM',      'global',dtype_i4,(/pcols/), acnum_idx)  ! counter for accumulated # timesteps

  ! SGS variability  -- These could be reset by CLUBB so they need to be grid only
  call pbuf_add_field('RELVAR',     'global',dtype_r8,(/pcols,pver/), relvar_idx)
  call pbuf_add_field('ACCRE_ENHAN','global',dtype_r8,(/pcols,pver/), accre_enhan_idx)

  ! Diagnostic fields needed for subcol_SILHS, need to be grid-only
  if (subcol_get_scheme() == 'SILHS') then
     call pbuf_add_field('QRAIN',   'global',dtype_r8,(/pcols,pver/), qrain_idx)
     call pbuf_add_field('QSNOW',   'global',dtype_r8,(/pcols,pver/), qsnow_idx)
     call pbuf_add_field('NRAIN',   'global',dtype_r8,(/pcols,pver/), nrain_idx)
     call pbuf_add_field('NSNOW',   'global',dtype_r8,(/pcols,pver/), nsnow_idx)
  end if
   


end subroutine micro_mg_cam_register

!===============================================================================

function micro_mg_cam_implements_cnst(name)

  ! Return true if specified constituent is implemented by the
  ! microphysics package

  character(len=*), intent(in) :: name        ! constituent name
  logical :: micro_mg_cam_implements_cnst    ! return value

  ! Local workspace
  integer :: m
  !-----------------------------------------------------------------------

  micro_mg_cam_implements_cnst = any(name == cnst_names)

end function micro_mg_cam_implements_cnst

!===============================================================================

subroutine micro_mg_cam_init_cnst(name, q, gcid)

  ! Initialize the microphysics constituents, if they are
  ! not read from the initial file.

  character(len=*), intent(in)  :: name     ! constituent name
  real(r8),         intent(out) :: q(:,:)   ! mass mixing ratio (gcol, plev)
  integer,          intent(in)  :: gcid(:)  ! global column id
  !-----------------------------------------------------------------------

  if (micro_mg_cam_implements_cnst(name)) q = 0.0_r8

end subroutine micro_mg_cam_init_cnst

!===============================================================================

subroutine micro_mg_cam_init(pbuf2d)
  use time_manager, only: is_first_step
  use micro_mg_utils, only: micro_mg_utils_init
  use micro_mg1_0, only: micro_mg_init1_0 => micro_mg_init
  use micro_mg1_5, only: micro_mg_init1_5 => micro_mg_init

  !-----------------------------------------------------------------------
  !
  ! Initialization for MG microphysics
  !
  !-----------------------------------------------------------------------

  type(physics_buffer_desc), pointer :: pbuf2d(:,:)

  integer :: m, mm
  logical :: history_amwg         ! output the variables used by the AMWG diag package
  logical :: history_budget       ! Output tendencies and state variables for CAM4
                                  ! temperature, water vapor, cloud ice and cloud
                                  ! liquid budgets.
  logical :: use_subcol_microp
  integer :: budget_histfile      ! output history file number for budget fields
  integer :: ierr

  character(128) :: errstring     ! return status (non-blank for error return)

  !-----------------------------------------------------------------------

  call phys_getopts(use_subcol_microp_out = use_subcol_microp)

  if (masterproc) then
     write(iulog,"(A,I2,A,I2)") "Initializing MG version ",micro_mg_version,".",micro_mg_sub_version
     if (.not. do_cldliq) &
          write(iulog,*) "MG prognostic cloud liquid has been turned off via namelist."
     if (.not. do_cldice) &
          write(iulog,*) "MG prognostic cloud ice has been turned off via namelist."
  end if

  select case (micro_mg_version)
  case (1)
     ! MG 1 does not initialize micro_mg_utils, so have to do it here.
     call micro_mg_utils_init(r8, rh2o, cpair, tmelt, latvap, latice, &
          errstring)
     call handle_errmsg(errstring, subname="micro_mg_utils_init")

     select case (micro_mg_sub_version)
     case (0)
        call micro_mg_init1_0( &
             r8, gravit, rair, rh2o, cpair, &
             rhoh2o, tmelt, latvap, latice, &
             rhmini, errstring)
     case (5)
        call micro_mg_init1_5( &
             r8, gravit, rair, rh2o, cpair, &
             tmelt, latvap, latice, rhmini, &
             microp_uniform, do_cldice, errstring)
     end select
  end select

  call handle_errmsg(errstring, subname="micro_mg_init")

  ! Register history variables
  do m = 1, ncnst
     call cnst_get_ind(cnst_names(m), mm)
     if ( any(mm == (/ ixcldliq, ixcldice /)) ) then
        ! mass mixing ratios
        call addfld(cnst_name(mm), 'kg/kg   ', pver, 'A', cnst_longname(mm)                   , phys_decomp)
        call addfld(sflxnam(mm),   'kg/m2/s ',    1, 'A', trim(cnst_name(mm))//' surface flux', phys_decomp)
     else if ( any(mm == (/ ixnumliq, ixnumice /)) ) then
        ! number concentrations
        call addfld(cnst_name(mm), '1/kg    ', pver, 'A', cnst_longname(mm)                   , phys_decomp)
        call addfld(sflxnam(mm),   '1/m2/s  ',    1, 'A', trim(cnst_name(mm))//' surface flux', phys_decomp)
     else
        call endrun( "micro_mg_cam_init: &
             &Could not call addfld for constituent with unknown units.")
     endif
  end do

  call addfld(apcnst(ixcldliq), 'kg/kg   ', pver, 'A', trim(cnst_name(ixcldliq))//' after physics'  , phys_decomp)
  call addfld(apcnst(ixcldice), 'kg/kg   ', pver, 'A', trim(cnst_name(ixcldice))//' after physics'  , phys_decomp)
  call addfld(bpcnst(ixcldliq), 'kg/kg   ', pver, 'A', trim(cnst_name(ixcldliq))//' before physics' , phys_decomp)
  call addfld(bpcnst(ixcldice), 'kg/kg   ', pver, 'A', trim(cnst_name(ixcldice))//' before physics' , phys_decomp)

  call addfld ('CME      ', 'kg/kg/s ', pver, 'A', 'Rate of cond-evap within the cloud'                      ,phys_decomp)
  call addfld ('PRODPREC ', 'kg/kg/s ', pver, 'A', 'Rate of conversion of condensate to precip'              ,phys_decomp)
  call addfld ('EVAPPREC ', 'kg/kg/s ', pver, 'A', 'Rate of evaporation of falling precip'                   ,phys_decomp)
  call addfld ('EVAPSNOW ', 'kg/kg/s ', pver, 'A', 'Rate of evaporation of falling snow'                     ,phys_decomp)
  call addfld ('HPROGCLD ', 'W/kg'    , pver, 'A', 'Heating from prognostic clouds'                          ,phys_decomp)
  call addfld ('FICE     ', 'fraction', pver, 'A', 'Fractional ice content within cloud'                     ,phys_decomp)
  call addfld ('ICWMRST  ', 'kg/kg   ', pver, 'A', 'Prognostic in-stratus water mixing ratio'                ,phys_decomp)
  call addfld ('ICIMRST  ', 'kg/kg   ', pver, 'A', 'Prognostic in-stratus ice mixing ratio'                  ,phys_decomp)

  ! MG microphysics diagnostics
  call addfld ('QCSEVAP  ', 'kg/kg/s ', pver, 'A', 'Rate of evaporation of falling cloud water'              ,phys_decomp)
  call addfld ('QISEVAP  ', 'kg/kg/s ', pver, 'A', 'Rate of sublimation of falling cloud ice'                ,phys_decomp)
  call addfld ('QVRES    ', 'kg/kg/s ', pver, 'A', 'Rate of residual condensation term'                      ,phys_decomp)
  call addfld ('CMEIOUT  ', 'kg/kg/s ', pver, 'A', 'Rate of deposition/sublimation of cloud ice'             ,phys_decomp)
  call addfld ('VTRMC    ', 'm/s     ', pver, 'A', 'Mass-weighted cloud water fallspeed'                     ,phys_decomp)
  call addfld ('VTRMI    ', 'm/s     ', pver, 'A', 'Mass-weighted cloud ice fallspeed'                       ,phys_decomp)
  call addfld ('QCSEDTEN ', 'kg/kg/s ', pver, 'A', 'Cloud water mixing ratio tendency from sedimentation'    ,phys_decomp)
  call addfld ('QISEDTEN ', 'kg/kg/s ', pver, 'A', 'Cloud ice mixing ratio tendency from sedimentation'      ,phys_decomp)
  call addfld ('PRAO     ', 'kg/kg/s ', pver, 'A', 'Accretion of cloud water by rain'                        ,phys_decomp)
  call addfld ('PRCO     ', 'kg/kg/s ', pver, 'A', 'Autoconversion of cloud water'                           ,phys_decomp)
  call addfld ('MNUCCCO  ', 'kg/kg/s ', pver, 'A', 'Immersion freezing of cloud water'                       ,phys_decomp)
  call addfld ('MNUCCTO  ', 'kg/kg/s ', pver, 'A', 'Contact freezing of cloud water'                         ,phys_decomp)
  call addfld ('MNUCCDO  ', 'kg/kg/s ', pver, 'A', 'Homogeneous and heterogeneous nucleation from vapor'     ,phys_decomp)
  call addfld ('MNUCCDOhet','kg/kg/s ', pver, 'A', 'Heterogeneous nucleation from vapor'                     ,phys_decomp)
  call addfld ('MSACWIO  ', 'kg/kg/s ', pver, 'A', 'Conversion of cloud water from rime-splintering'         ,phys_decomp)
  call addfld ('PSACWSO  ', 'kg/kg/s ', pver, 'A', 'Accretion of cloud water by snow'                        ,phys_decomp)
  call addfld ('BERGSO   ', 'kg/kg/s ', pver, 'A', 'Conversion of cloud water to snow from bergeron'         ,phys_decomp)
  call addfld ('BERGO    ', 'kg/kg/s ', pver, 'A', 'Conversion of cloud water to cloud ice from bergeron'    ,phys_decomp)
  call addfld ('MELTO    ', 'kg/kg/s ', pver, 'A', 'Melting of cloud ice'                                    ,phys_decomp)
  call addfld ('HOMOO    ', 'kg/kg/s ', pver, 'A', 'Homogeneous freezing of cloud water'                     ,phys_decomp)
  call addfld ('QCRESO   ', 'kg/kg/s ', pver, 'A', 'Residual condensation term for cloud water'              ,phys_decomp)
  call addfld ('PRCIO    ', 'kg/kg/s ', pver, 'A', 'Autoconversion of cloud ice'                             ,phys_decomp)
  call addfld ('PRAIO    ', 'kg/kg/s ', pver, 'A', 'Accretion of cloud ice by rain'                          ,phys_decomp)
  call addfld ('QIRESO   ', 'kg/kg/s ', pver, 'A', 'Residual deposition term for cloud ice'                  ,phys_decomp)
  call addfld ('MNUCCRO  ', 'kg/kg/s ', pver, 'A', 'Heterogeneous freezing of rain to snow'                  ,phys_decomp)
  call addfld ('PRACSO   ', 'kg/kg/s ', pver, 'A', 'Accretion of rain by snow'                               ,phys_decomp)
  call addfld ('MELTSDT  ', 'W/kg    ', pver, 'A', 'Latent heating rate due to melting of snow'              ,phys_decomp)
  call addfld ('FRZRDT   ', 'W/kg    ', pver, 'A', 'Latent heating rate due to homogeneous freezing of rain' ,phys_decomp)

  ! History variables for CAM5 microphysics
  call addfld ('MPDT     ', 'W/kg    ', pver, 'A', 'Heating tendency - Morrison microphysics'                ,phys_decomp)
  call addfld ('MPDQ     ', 'kg/kg/s ', pver, 'A', 'Q tendency - Morrison microphysics'                      ,phys_decomp)
  call addfld ('MPDLIQ   ', 'kg/kg/s ', pver, 'A', 'CLDLIQ tendency - Morrison microphysics'                 ,phys_decomp)
  call addfld ('MPDICE   ', 'kg/kg/s ', pver, 'A', 'CLDICE tendency - Morrison microphysics'                 ,phys_decomp)
  call addfld ('MPDW2V   ', 'kg/kg/s ', pver, 'A', 'Water <--> Vapor tendency - Morrison microphysics'       ,phys_decomp)
  call addfld ('MPDW2I   ', 'kg/kg/s ', pver, 'A', 'Water <--> Ice tendency - Morrison microphysics'         ,phys_decomp)
  call addfld ('MPDW2P   ', 'kg/kg/s ', pver, 'A', 'Water <--> Precip tendency - Morrison microphysics'      ,phys_decomp)
  call addfld ('MPDI2V   ', 'kg/kg/s ', pver, 'A', 'Ice <--> Vapor tendency - Morrison microphysics'         ,phys_decomp)
  call addfld ('MPDI2W   ', 'kg/kg/s ', pver, 'A', 'Ice <--> Water tendency - Morrison microphysics'         ,phys_decomp)
  call addfld ('MPDI2P   ', 'kg/kg/s ', pver, 'A', 'Ice <--> Precip tendency - Morrison microphysics'        ,phys_decomp)
  call addfld ('ICWNC    ', 'm-3     ', pver, 'A', 'Prognostic in-cloud water number conc'                   ,phys_decomp)
  call addfld ('ICINC    ', 'm-3     ', pver, 'A', 'Prognostic in-cloud ice number conc'                     ,phys_decomp)
  call addfld ('EFFLIQ_IND','Micron  ', pver, 'A', 'Prognostic droplet effective radius (indirect effect)'   ,phys_decomp)
  call addfld ('CDNUMC   ', '1/m2    ', 1,    'A', 'Vertically-integrated droplet concentration'             ,phys_decomp)
  call addfld ('MPICLWPI ', 'kg/m2   ', 1,    'A', 'Vertically-integrated &
       &in-cloud Initial Liquid WP (Before Micro)' ,phys_decomp)
  call addfld ('MPICIWPI ', 'kg/m2   ', 1,    'A', 'Vertically-integrated &
       &in-cloud Initial Ice WP (Before Micro)'    ,phys_decomp)

  ! This is provided as an example on how to write out subcolumn output
  ! NOTE -- only 'I' should be used for sub-column fields as subc-columns could shift from time-step to time-step
  if (use_subcol_microp) then
     call addfld('FICE_SCOL', 'fraction', psubcols*pver, 'I', &
                 'Sub-column fractional ice content within cloud', phys_decomp, &
                 mdimnames=(/'psubcols','lev     '/), flag_xyfill=.true., fill_value=1.e30_r8)
  end if

  ! Averaging for cloud particle number and size
  call addfld ('AWNC     ', 'm-3     ', pver, 'A', 'Average cloud water number conc'                         ,phys_decomp)
  call addfld ('AWNI     ', 'm-3     ', pver, 'A', 'Average cloud ice number conc'                           ,phys_decomp)
  call addfld ('AREL     ', 'Micron  ', pver, 'A', 'Average droplet effective radius'                        ,phys_decomp)
  call addfld ('AREI     ', 'Micron  ', pver, 'A', 'Average ice effective radius'                            ,phys_decomp)
  ! Frequency arrays for above
  call addfld ('FREQL    ', 'fraction', pver, 'A', 'Fractional occurrence of liquid'                          ,phys_decomp)
  call addfld ('FREQI    ', 'fraction', pver, 'A', 'Fractional occurrence of ice'                             ,phys_decomp)

  ! Average cloud top particle size and number (liq, ice) and frequency
  call addfld ('ACTREL   ', 'Micron  ', 1,    'A', 'Average Cloud Top droplet effective radius'              ,phys_decomp)
  call addfld ('ACTREI   ', 'Micron  ', 1,    'A', 'Average Cloud Top ice effective radius'                  ,phys_decomp)
  call addfld ('ACTNL    ', 'Micron  ', 1,    'A', 'Average Cloud Top droplet number'                        ,phys_decomp)
  call addfld ('ACTNI    ', 'Micron  ', 1,    'A', 'Average Cloud Top ice number'                            ,phys_decomp)

  call addfld ('FCTL     ', 'fraction', 1,    'A', 'Fractional occurrence of cloud top liquid'                ,phys_decomp)
  call addfld ('FCTI     ', 'fraction', 1,    'A', 'Fractional occurrence of cloud top ice'                   ,phys_decomp)

  call addfld ('LS_FLXPRC', 'kg/m2/s', pverp, 'A', 'ls stratiform gbm interface rain+snow flux', phys_decomp)
  call addfld ('LS_FLXSNW', 'kg/m2/s', pverp, 'A', 'ls stratiform gbm interface snow flux', phys_decomp)

  call addfld ('REL', 'micron', pver, 'A', 'MG REL stratiform cloud effective radius liquid', phys_decomp)
  call addfld ('REI', 'micron', pver, 'A', 'MG REI stratiform cloud effective radius ice', phys_decomp)
  call addfld ('LS_REFFRAIN', 'micron', pver, 'A', 'ls stratiform rain effective radius', phys_decomp)
  call addfld ('LS_REFFSNOW', 'micron', pver, 'A', 'ls stratiform snow effective radius', phys_decomp)
  call addfld ('CV_REFFLIQ', 'micron', pver, 'A', 'convective cloud liq effective radius', phys_decomp)
  call addfld ('CV_REFFICE', 'micron', pver, 'A', 'convective cloud ice effective radius', phys_decomp)

  ! diagnostic precip
  call addfld ('QRAIN   ','kg/kg   ',pver, 'A','Diagnostic grid-mean rain mixing ratio'         ,phys_decomp)
  call addfld ('QSNOW   ','kg/kg   ',pver, 'A','Diagnostic grid-mean snow mixing ratio'         ,phys_decomp)
  call addfld ('NRAIN   ','m-3     ',pver, 'A','Diagnostic grid-mean rain number conc'         ,phys_decomp)
  call addfld ('NSNOW   ','m-3     ',pver, 'A','Diagnostic grid-mean snow number conc'         ,phys_decomp)

  ! size of precip
  call addfld ('RERCLD   ','m      ',pver, 'A','Diagnostic effective radius of Liquid Cloud and Rain' ,phys_decomp)
  call addfld ('DSNOW   ','m       ',pver, 'A','Diagnostic grid-mean snow diameter'         ,phys_decomp)

  ! diagnostic radar reflectivity, cloud-averaged
  call addfld ('REFL  ','DBz  ',pver, 'A','94 GHz radar reflectivity'       ,phys_decomp)
  call addfld ('AREFL  ','DBz  ',pver, 'A','Average 94 GHz radar reflectivity'       ,phys_decomp)
  call addfld ('FREFL  ','fraction  ',pver, 'A','Fractional occurrence of radar reflectivity'       ,phys_decomp)

  call addfld ('CSRFL  ','DBz  ',pver, 'A','94 GHz radar reflectivity (CloudSat thresholds)'       ,phys_decomp)
  call addfld ('ACSRFL  ','DBz  ',pver, 'A','Average 94 GHz radar reflectivity (CloudSat thresholds)'       ,phys_decomp)
  call addfld ('FCSRFL  ','fraction  ',pver, 'A','Fractional occurrence of radar reflectivity (CloudSat thresholds)' &
       ,phys_decomp)

  call addfld ('AREFLZ ','mm^6/m^3 ',pver, 'A','Average 94 GHz radar reflectivity'       ,phys_decomp)

  ! Aerosol information
  call addfld ('NCAL    ','1/m3   ',pver, 'A','Number Concentation Activated for Liquid',phys_decomp)
  call addfld ('NCAI    ','1/m3   ',pver, 'A','Number Concentation Activated for Ice',phys_decomp)

  ! Average rain and snow mixing ratio (Q), number (N) and diameter (D), with frequency
  call addfld ('AQRAIN   ','kg/kg   ',pver, 'A','Average rain mixing ratio'         ,phys_decomp)
  call addfld ('AQSNOW   ','kg/kg   ',pver, 'A','Average snow mixing ratio'         ,phys_decomp)
  call addfld ('ANRAIN   ','m-3     ',pver, 'A','Average rain number conc'         ,phys_decomp)
  call addfld ('ANSNOW   ','m-3     ',pver, 'A','Average snow number conc'         ,phys_decomp)
  call addfld ('ADRAIN   ','Micron  ',pver, 'A','Average rain effective Diameter'         ,phys_decomp)
  call addfld ('ADSNOW   ','Micron  ',pver, 'A','Average snow effective Diameter'         ,phys_decomp)
  call addfld ('FREQR  ','fraction  ',pver, 'A','Fractional occurrence of rain'       ,phys_decomp)
  call addfld ('FREQS  ','fraction  ',pver, 'A','Fractional occurrence of snow'       ,phys_decomp)

  ! precipitation efficiency & other diagnostic fields
  call addfld('PE'    , '1',       1, 'A', 'Stratiform Precipitation Efficiency  (precip/cmeliq)',       phys_decomp )
  call addfld('APRL'  , 'm/s',     1, 'A', 'Average Stratiform Precip Rate over efficiency calculation', phys_decomp )
  call addfld('PEFRAC', '1',       1, 'A', 'Fraction of timesteps precip efficiency reported',           phys_decomp )
  call addfld('VPRCO' , 'kg/kg/s', 1, 'A', 'Vertical average of autoconversion rate',                         phys_decomp )
  call addfld('VPRAO' , 'kg/kg/s', 1, 'A', 'Vertical average of accretion rate',                    phys_decomp )
  call addfld('RACAU' , 'kg/kg/s', 1, 'A', 'Accretion/autoconversion ratio from vertical average',       phys_decomp )

  ! determine the add_default fields
  call phys_getopts(history_amwg_out           = history_amwg         , &
       history_budget_out         = history_budget       , &
       history_budget_histfile_num_out = budget_histfile)

  if (history_amwg) then
     call add_default ('FICE    ', 1, ' ')
     call add_default ('AQRAIN   ', 1, ' ')
     call add_default ('AQSNOW   ', 1, ' ')
     call add_default ('ANRAIN   ', 1, ' ')
     call add_default ('ANSNOW   ', 1, ' ')
     call add_default ('AREI     ', 1, ' ')
     call add_default ('AREL     ', 1, ' ')
     call add_default ('AWNC     ', 1, ' ')
     call add_default ('AWNI     ', 1, ' ')
     call add_default ('CDNUMC   ', 1, ' ')
     call add_default ('FREQR    ', 1, ' ')
     call add_default ('FREQS    ', 1, ' ')
     call add_default ('FREQL    ', 1, ' ')
     call add_default ('FREQI    ', 1, ' ')
     do m = 1, ncnst
        call cnst_get_ind(cnst_names(m), mm)
        call add_default(cnst_name(mm), 1, ' ')
        ! call add_default(sflxnam(mm),   1, ' ')
     end do
  end if

  if ( history_budget ) then
     call add_default ('EVAPSNOW ', budget_histfile, ' ')
     call add_default ('EVAPPREC ', budget_histfile, ' ')
     call add_default ('QVRES    ', budget_histfile, ' ')
     call add_default ('QISEVAP  ', budget_histfile, ' ')
     call add_default ('QCSEVAP  ', budget_histfile, ' ')
     call add_default ('QISEDTEN ', budget_histfile, ' ')
     call add_default ('QCSEDTEN ', budget_histfile, ' ')
     call add_default ('QIRESO   ', budget_histfile, ' ')
     call add_default ('QCRESO   ', budget_histfile, ' ')
     call add_default ('PSACWSO  ', budget_histfile, ' ')
     call add_default ('PRCO     ', budget_histfile, ' ')
     call add_default ('PRCIO    ', budget_histfile, ' ')
     call add_default ('PRAO     ', budget_histfile, ' ')
     call add_default ('PRAIO    ', budget_histfile, ' ')
     call add_default ('PRACSO   ', budget_histfile, ' ')
     call add_default ('MSACWIO  ', budget_histfile, ' ')
     call add_default ('MPDW2V   ', budget_histfile, ' ')
     call add_default ('MPDW2P   ', budget_histfile, ' ')
     call add_default ('MPDW2I   ', budget_histfile, ' ')
     call add_default ('MPDT     ', budget_histfile, ' ')
     call add_default ('MPDQ     ', budget_histfile, ' ')
     call add_default ('MPDLIQ   ', budget_histfile, ' ')
     call add_default ('MPDICE   ', budget_histfile, ' ')
     call add_default ('MPDI2W   ', budget_histfile, ' ')
     call add_default ('MPDI2V   ', budget_histfile, ' ')
     call add_default ('MPDI2P   ', budget_histfile, ' ')
     call add_default ('MNUCCTO  ', budget_histfile, ' ')
     call add_default ('MNUCCRO  ', budget_histfile, ' ')
     call add_default ('MNUCCCO  ', budget_histfile, ' ')
     call add_default ('MELTSDT  ', budget_histfile, ' ')
     call add_default ('MELTO    ', budget_histfile, ' ')
     call add_default ('HOMOO    ', budget_histfile, ' ')
     call add_default ('FRZRDT   ', budget_histfile, ' ')
     call add_default ('CMEIOUT  ', budget_histfile, ' ')
     call add_default ('BERGSO   ', budget_histfile, ' ')
     call add_default ('BERGO    ', budget_histfile, ' ')

     call add_default(cnst_name(ixcldliq), budget_histfile, ' ')
     call add_default(cnst_name(ixcldice), budget_histfile, ' ')
     call add_default(apcnst   (ixcldliq), budget_histfile, ' ')
     call add_default(apcnst   (ixcldice), budget_histfile, ' ')
     call add_default(bpcnst   (ixcldliq), budget_histfile, ' ')
     call add_default(bpcnst   (ixcldice), budget_histfile, ' ')

  end if

  ! physics buffer indices
  ast_idx      = pbuf_get_index('AST')
  aist_idx     = pbuf_get_index('AIST')
  alst_idx     = pbuf_get_index('ALST')
  cld_idx      = pbuf_get_index('CLD')
  concld_idx   = pbuf_get_index('CONCLD')

  naai_idx     = pbuf_get_index('NAAI')
  naai_hom_idx = pbuf_get_index('NAAI_HOM')
  npccn_idx    = pbuf_get_index('NPCCN')
  rndst_idx    = pbuf_get_index('RNDST')
  nacon_idx    = pbuf_get_index('NACON')

  prec_str_idx = pbuf_get_index('PREC_STR')
  snow_str_idx = pbuf_get_index('SNOW_STR')
  prec_sed_idx = pbuf_get_index('PREC_SED')
  snow_sed_idx = pbuf_get_index('SNOW_SED')
  prec_pcw_idx = pbuf_get_index('PREC_PCW')
  snow_pcw_idx = pbuf_get_index('SNOW_PCW')

  cmeliq_idx   = pbuf_get_index('CMELIQ')

  ! These fields may have been added, so don't abort if they have not been
  qrain_idx    = pbuf_get_index('QRAIN', ierr)
  qsnow_idx    = pbuf_get_index('QSNOW', ierr)
  nrain_idx    = pbuf_get_index('NRAIN', ierr)
  nsnow_idx    = pbuf_get_index('NSNOW', ierr)

  ! Initialize physics buffer grid fields for accumulating precip and condensation
  if (is_first_step()) then
     call pbuf_set_field(pbuf2d, cldo_idx,   0._r8)
     call pbuf_set_field(pbuf2d, cc_t_idx,   0._r8)
     call pbuf_set_field(pbuf2d, cc_qv_idx,  0._r8)
     call pbuf_set_field(pbuf2d, cc_ql_idx,  0._r8)
     call pbuf_set_field(pbuf2d, cc_qi_idx,  0._r8)
     call pbuf_set_field(pbuf2d, cc_nl_idx,  0._r8)
     call pbuf_set_field(pbuf2d, cc_ni_idx,  0._r8)
     call pbuf_set_field(pbuf2d, cc_qlst_idx,0._r8)
     call pbuf_set_field(pbuf2d, acpr_idx,   0._r8)
     call pbuf_set_field(pbuf2d, acgcme_idx, 0._r8)
     call pbuf_set_field(pbuf2d, acnum_idx,  0)
     call pbuf_set_field(pbuf2d, relvar_idx, 2._r8)
     call pbuf_set_field(pbuf2d, accre_enhan_idx, 1._r8)

     if (qrain_idx > 0)   call pbuf_set_field(pbuf2d, qrain_idx, 0._r8)
     if (qsnow_idx > 0)   call pbuf_set_field(pbuf2d, qsnow_idx, 0._r8)
     if (nrain_idx > 0)   call pbuf_set_field(pbuf2d, nrain_idx, 0._r8)
     if (nsnow_idx > 0)   call pbuf_set_field(pbuf2d, nsnow_idx, 0._r8)

     ! If sub-columns turned on, need to set the sub-column fields as well
     if (use_subcol_microp) then
        call pbuf_set_field(pbuf2d, cldo_idx,   0._r8, col_type=col_type_subcol)
        call pbuf_set_field(pbuf2d, cc_t_idx,   0._r8, col_type=col_type_subcol)
        call pbuf_set_field(pbuf2d, cc_qv_idx,  0._r8, col_type=col_type_subcol)
        call pbuf_set_field(pbuf2d, cc_ql_idx,  0._r8, col_type=col_type_subcol)
        call pbuf_set_field(pbuf2d, cc_qi_idx,  0._r8, col_type=col_type_subcol)
        call pbuf_set_field(pbuf2d, cc_nl_idx,  0._r8, col_type=col_type_subcol)
        call pbuf_set_field(pbuf2d, cc_ni_idx,  0._r8, col_type=col_type_subcol)
        call pbuf_set_field(pbuf2d, cc_qlst_idx,0._r8, col_type=col_type_subcol)
     end if

  end if

end subroutine micro_mg_cam_init


!===============================================================================

subroutine micro_mg_cam_tend(state, ptend, dtime, pbuf)

  use micro_mg_utils, only: size_dist_param_basic, size_dist_param_liq, &
       mg_liq_props, mg_ice_props, avg_diameter, rhoi, rhosn, rhow, rhows, &
       qsmall, mincld
  use micro_mg1_0, only: micro_mg_tend1_0 => micro_mg_tend
  use micro_mg1_5, only: micro_mg_tend1_5 => micro_mg_tend, &
       micro_mg_get_cols1_5 => micro_mg_get_cols

  use ppgrid,          only: pcols
  use physics_buffer,  only: pbuf_col_type_index
  use subcol,          only: subcol_field_avg

  type(physics_state),         intent(in)    :: state
  type(physics_ptend),         intent(out)   :: ptend
  real(r8),                    intent(in)    :: dtime
  type(physics_buffer_desc),   pointer       :: pbuf(:)

  ! Local variables
  logical :: microp_uniform = .false. ! True = configure microphysics for sub-columns
  ! False = use in regular mode w/o sub-columns
  integer :: lchnk, ncol, psetcols, ngrdcol

  integer :: i, k, itim_old, it

  real(r8), pointer :: naai(:,:)      ! ice nucleation number
  real(r8), pointer :: naai_hom(:,:)  ! ice nucleation number (homogeneous)
  real(r8), pointer :: npccn(:,:)     ! liquid activation number tendency
  real(r8), pointer :: rndst(:,:,:)
  real(r8), pointer :: nacon(:,:,:)

  real(r8), pointer :: prec_str(:)          ! [Total] Sfc flux of precip from stratiform [ m/s ]
  real(r8), pointer :: snow_str(:)          ! [Total] Sfc flux of snow from stratiform   [ m/s ]
  real(r8), pointer :: prec_sed(:)          ! Surface flux of total cloud water from sedimentation
  real(r8), pointer :: snow_sed(:)          ! Surface flux of cloud ice from sedimentation
  real(r8), pointer :: prec_pcw(:)          ! Sfc flux of precip from microphysics [ m/s ]
  real(r8), pointer :: snow_pcw(:)          ! Sfc flux of snow from microphysics [ m/s ]

  real(r8), pointer :: ast(:,:)          ! Relative humidity cloud fraction
  real(r8), pointer :: alst_mic(:,:)
  real(r8), pointer :: aist_mic(:,:)
  real(r8), pointer :: cldo(:,:)         ! Old cloud fraction
  real(r8), pointer :: nevapr(:,:)       ! Evaporation of total precipitation (rain + snow)
  real(r8), pointer :: relvar(:,:)       ! relative variance of cloud water
  real(r8), pointer :: accre_enhan(:,:)  ! optional accretion enhancement for experimentation
  real(r8), pointer :: prain(:,:)        ! Total precipitation (rain + snow)
  real(r8), pointer :: dei(:,:)          ! Ice effective diameter (meters) (AG: microns?)
  real(r8), pointer :: mu(:,:)           ! Size distribution shape parameter for radiation
  real(r8), pointer :: lambdac(:,:)      ! Size distribution slope parameter for radiation
  real(r8), pointer :: des(:,:)          ! Snow effective diameter (m)

  real(r8) :: rho(state%psetcols,pver)
  real(r8) :: ncic(state%psetcols,pver)
  real(r8) :: niic(state%psetcols,pver)

  real(r8) :: rate1cld(state%psetcols,pver) ! array to hold rate1ord_cw2pr_st from microphysics

  real(r8) :: tlat(state%psetcols,pver)
  real(r8) :: qvlat(state%psetcols,pver)
  real(r8) :: qcten(state%psetcols,pver)
  real(r8) :: qiten(state%psetcols,pver)
  real(r8) :: ncten(state%psetcols,pver)
  real(r8) :: niten(state%psetcols,pver)
  real(r8) :: prect(state%psetcols)
  real(r8) :: preci(state%psetcols)


  real(r8) :: evapsnow(state%psetcols,pver)                    ! Local evaporation of snow
  real(r8) :: prodsnow(state%psetcols,pver)                    ! Local production of snow
  real(r8) :: cmeice(state%psetcols,pver)                      ! Rate of cond-evap of ice within the cloud
  real(r8) :: qsout(state%psetcols,pver)                       ! Snow mixing ratio
  real(r8) :: rflx(state%psetcols,pver+1)                   ! grid-box average rain flux (kg m^-2 s^-1)
  real(r8) :: sflx(state%psetcols,pver+1)                   ! grid-box average snow flux (kg m^-2 s^-1)
  real(r8) :: qrout(state%psetcols,pver)                     ! Rain mixing ratio
  real(r8) :: reff_rain(state%psetcols,pver)                ! rain effective radius (um)
  real(r8) :: reff_snow(state%psetcols,pver)                ! snow effective radius (um)
  real(r8) :: qcsevap(state%psetcols,pver)                     ! Evaporation of falling cloud water
  real(r8) :: qisevap(state%psetcols,pver)                     ! Sublimation of falling cloud ice
  real(r8) :: qvres(state%psetcols,pver)                       ! Residual condensation term to remove excess saturation
  real(r8) :: cmeiout(state%psetcols,pver)                     ! Deposition/sublimation rate of cloud ice
  real(r8) :: vtrmc(state%psetcols,pver)                       ! Mass-weighted cloud water fallspeed
  real(r8) :: vtrmi(state%psetcols,pver)                       ! Mass-weighted cloud ice fallspeed
  real(r8) :: qcsedten(state%psetcols,pver)                    ! Cloud water mixing ratio tendency from sedimentation
  real(r8) :: qisedten(state%psetcols,pver)                    ! Cloud ice mixing ratio tendency from sedimentation
  real(r8) :: prao(state%psetcols,pver)
  real(r8) :: prco(state%psetcols,pver)
  real(r8) :: mnuccco(state%psetcols,pver)
  real(r8) :: mnuccto(state%psetcols,pver)
  real(r8) :: msacwio(state%psetcols,pver)
  real(r8) :: psacwso(state%psetcols,pver)
  real(r8) :: bergso(state%psetcols,pver)
  real(r8) :: bergo(state%psetcols,pver)
  real(r8) :: melto(state%psetcols,pver)
  real(r8) :: homoo(state%psetcols,pver)
  real(r8) :: qcreso(state%psetcols,pver)
  real(r8) :: prcio(state%psetcols,pver)
  real(r8) :: praio(state%psetcols,pver)
  real(r8) :: qireso(state%psetcols,pver)
  real(r8) :: mnuccro(state%psetcols,pver)
  real(r8) :: pracso (state%psetcols,pver)
  real(r8) :: meltsdt(state%psetcols,pver)
  real(r8) :: frzrdt (state%psetcols,pver)
  real(r8) :: mnuccdo(state%psetcols,pver)
  real(r8) :: nrout(state%psetcols,pver)
  real(r8) :: nsout(state%psetcols,pver)
  real(r8) :: refl(state%psetcols,pver)   ! analytic radar reflectivity
  real(r8) :: arefl(state%psetcols,pver)  !average reflectivity will zero points outside valid range
  real(r8) :: areflz(state%psetcols,pver) !average reflectivity in z.
  real(r8) :: frefl(state%psetcols,pver)
  real(r8) :: csrfl(state%psetcols,pver)  !cloudsat reflectivity
  real(r8) :: acsrfl(state%psetcols,pver) !cloudsat average
  real(r8) :: fcsrfl(state%psetcols,pver)
  real(r8) :: rercld(state%psetcols,pver) ! effective radius calculation for rain + cloud
  real(r8) :: ncai(state%psetcols,pver)   ! output number conc of ice nuclei available (1/m3)
  real(r8) :: ncal(state%psetcols,pver)   ! output number conc of CCN (1/m3)
  real(r8) :: qrout2(state%psetcols,pver)
  real(r8) :: qsout2(state%psetcols,pver)
  real(r8) :: nrout2(state%psetcols,pver)
  real(r8) :: nsout2(state%psetcols,pver)
  real(r8) :: drout2(state%psetcols,pver) ! mean rain particle diameter (m)
  real(r8) :: dsout2(state%psetcols,pver) ! mean snow particle diameter (m)
  real(r8) :: freqs(state%psetcols,pver)
  real(r8) :: freqr(state%psetcols,pver)
  real(r8) :: nfice(state%psetcols,pver)

  real(r8) :: mnuccdohet(state%psetcols,pver)

  ! physics buffer fields for COSP simulator
  real(r8), pointer :: mgflxprc(:,:)     ! MG grid-box mean flux_large_scale_cloud_rain+snow at interfaces (kg/m2/s)
  real(r8), pointer :: mgflxsnw(:,:)     ! MG grid-box mean flux_large_scale_cloud_snow at interfaces (kg/m2/s)
  real(r8), pointer :: mgmrprc(:,:)      ! MG grid-box mean mixingratio_large_scale_cloud_rain+snow at interfaces (kg/kg)
  real(r8), pointer :: mgmrsnw(:,:)      ! MG grid-box mean mixingratio_large_scale_cloud_snow at interfaces (kg/kg)
  real(r8), pointer :: mgreffrain_grid(:,:)   ! MG diagnostic rain effective radius (um)
  real(r8), pointer :: mgreffsnow_grid(:,:)   ! MG diagnostic snow effective radius (um)
  real(r8), pointer :: cvreffliq(:,:)    ! convective cloud liquid effective radius (um)
  real(r8), pointer :: cvreffice(:,:)    ! convective cloud ice effective radius (um)

  ! physics buffer fields used with CARMA
  real(r8), pointer, dimension(:,:) :: tnd_qsnow    ! external tendency on snow mass (kg/kg/s)
  real(r8), pointer, dimension(:,:) :: tnd_nsnow    ! external tendency on snow number(#/kg/s)
  real(r8), pointer, dimension(:,:) :: re_ice       ! ice effective radius (m)

  real(r8), pointer :: rate1ord_cw2pr_st(:,:) ! 1st order rate for direct conversion of
  ! strat. cloud water to precip (1/s)    ! rce 2010/05/01
  real(r8), pointer :: wsedl(:,:)        ! Sedimentation velocity of liquid stratus cloud droplet [ m/s ]


  real(r8), pointer :: CC_T(:,:)         ! Grid-mean microphysical tendency
  real(r8), pointer :: CC_qv(:,:)        ! Grid-mean microphysical tendency
  real(r8), pointer :: CC_ql(:,:)        ! Grid-mean microphysical tendency
  real(r8), pointer :: CC_qi(:,:)        ! Grid-mean microphysical tendency
  real(r8), pointer :: CC_nl(:,:)        ! Grid-mean microphysical tendency
  real(r8), pointer :: CC_ni(:,:)        ! Grid-mean microphysical tendency
  real(r8), pointer :: CC_qlst(:,:)      ! In-liquid stratus microphysical tendency

  real(r8), pointer :: qme(:,:)

  ! A local copy of state is used for diagnostic calculations
  type(physics_state) :: state_loc
  type(physics_ptend) :: ptend_loc

  real(r8) :: icecldf(state%psetcols,pver)                     ! Ice cloud fraction
  real(r8) :: liqcldf(state%psetcols,pver)                     ! Liquid cloud fraction (combined into cloud)

  real(r8), pointer :: rel(:,:)          ! Liquid effective drop radius (microns)
  real(r8), pointer :: rei(:,:)          ! Ice effective drop size (microns)
  real(r8) :: rel_fn(state%psetcols,pver)         ! Ice effective drop size at fixed number (indirect effect) (microns)

  ! in-cloud water quantities adjusted for convective water
  real(r8) :: allcld_ice(state%psetcols,pver)                 ! All-cloud cloud ice
  real(r8) :: allcld_liq(state%psetcols,pver)                 ! All-cloud liquid

  real(r8), pointer :: cmeliq(:,:)

  real(r8), pointer :: cld(:,:)          ! Total cloud fraction
  real(r8), pointer :: concld(:,:)       ! Convective cloud fraction
  real(r8), pointer :: iciwpst(:,:)      ! Stratiform in-cloud ice water path for radiation
  real(r8), pointer :: iclwpst(:,:)      ! Stratiform in-cloud liquid water path for radiation
  real(r8), pointer :: cldfsnow(:,:)     ! Cloud fraction for liquid+snow
  real(r8), pointer :: icswp(:,:)        ! In-cloud snow water path

  real(r8) :: icimrst(state%psetcols,pver)                     ! In stratus ice mixing ratio
  real(r8) :: icwmrst(state%psetcols,pver)                     ! In stratus water mixing ratio
  real(r8) :: icinc(state%psetcols,pver)                       ! In cloud ice number conc
  real(r8) :: icwnc(state%psetcols,pver)                       ! In cloud water number conc

  real(r8) :: iclwpi(state%psetcols)                           ! Vertically-integrated in-cloud Liquid WP before microphysics
  real(r8) :: iciwpi(state%psetcols)                           ! Vertically-integrated in-cloud Ice WP before microphysics

  ! Averaging arrays for effective radius and number....
  real(r8) :: efiout_grid(pcols,pver)
  real(r8) :: efcout_grid(pcols,pver)
  real(r8) :: ncout_grid(pcols,pver)
  real(r8) :: niout_grid(pcols,pver)
  real(r8) :: freqi_grid(pcols,pver)
  real(r8) :: freql_grid(pcols,pver)

  real(r8) :: cdnumc_grid(pcols)                           ! Vertically-integrated droplet concentration
  real(r8) :: icecldf_grid_out(pcols,pver)                 ! Ice cloud fraction
  real(r8) :: liqcldf_grid_out(pcols,pver)                 ! Liquid cloud fraction (combined into cloud)
  real(r8) :: icimrst_grid_out(pcols,pver)                 ! In stratus ice mixing ratio
  real(r8) :: icwmrst_grid_out(pcols,pver)                 ! In stratus water mixing ratio

  ! Average cloud top radius & number
  real(r8) :: ctrel_grid(pcols)
  real(r8) :: ctrei_grid(pcols)
  real(r8) :: ctnl_grid(pcols)
  real(r8) :: ctni_grid(pcols)
  real(r8) :: fcti_grid(pcols)
  real(r8) :: fctl_grid(pcols)

  real(r8) :: ftem_grid(pcols,pver)

  ! Variables for precip efficiency calculation
  real(r8) :: minlwp        ! LWP threshold

  real(r8), pointer, dimension(:) :: acprecl_grid ! accumulated precip across timesteps
  real(r8), pointer, dimension(:) :: acgcme_grid  ! accumulated condensation across timesteps
  integer,  pointer, dimension(:) :: acnum_grid   ! counter for # timesteps accumulated

  ! Variables for liquid water path and column condensation
  real(r8) :: tgliqwp_grid(pcols)   ! column liquid
  real(r8) :: tgcmeliq_grid(pcols)  ! column condensation rate (units)

  real(r8) :: pe_grid(pcols)        ! precip efficiency for output
  real(r8) :: pefrac_grid(pcols)    ! fraction of time precip efficiency is written out
  real(r8) :: tpr_grid(pcols)       ! average accumulated precipitation rate in pe calculation

  ! variables for autoconversion and accretion vertical averages
  real(r8) :: vprco_grid(pcols)     ! vertical average autoconversion
  real(r8) :: vprao_grid(pcols)     ! vertical average accretion
  real(r8) :: racau_grid(pcols)     ! ratio of vertical averages
  integer  :: cnt_grid(pcols)       ! counters
  logical  :: lq(pcnst)

  real(r8) :: qc(state%psetcols,pver)    ! cloud water mixing ratio (kg/kg)
  real(r8) :: qi(state%psetcols,pver)    ! cloud ice mixing ratio (kg/kg)
  real(r8) :: nc(state%psetcols,pver)    ! cloud water number conc (1/kg)
  real(r8) :: ni(state%psetcols,pver)    ! cloud ice number conc (1/kg)

  real(r8) :: icimrst_grid(pcols,pver) ! stratus ice mixing ratio - on grid
  real(r8) :: icwmrst_grid(pcols,pver) ! stratus water mixing ratio - on grid

  real(r8),pointer :: lambdac_grid(:,:)
  real(r8),pointer :: mu_grid(:,:)
  real(r8),pointer :: rel_grid(:,:)
  real(r8),pointer :: rei_grid(:,:)
  real(r8),pointer :: dei_grid(:,:)
  real(r8),pointer :: des_grid(:,:)
  real(r8),pointer :: iclwpst_grid(:,:)

  real(r8) :: rho_grid(pcols,pver)
  real(r8) :: liqcldf_grid(pcols,pver)
  real(r8) :: qsout_grid(pcols,pver)
  real(r8) :: ncic_grid(pcols,pver)
  real(r8) :: niic_grid(pcols,pver)
  real(r8) :: rel_fn_grid(pcols,pver)    ! Ice effective drop size at fixed number (indirect effect) (microns) - on grid
  real(r8) :: qrout_grid(pcols,pver)
  real(r8) :: drout2_grid(pcols,pver)
  real(r8) :: dsout2_grid(pcols,pver)
  real(r8) :: nsout_grid(pcols,pver)
  real(r8) :: nrout_grid(pcols,pver)
  real(r8) :: reff_rain_grid(pcols,pver)
  real(r8) :: reff_snow_grid(pcols,pver)
  real(r8) :: cld_grid(pcols,pver)
  real(r8) :: pdel_grid(pcols,pver)
  real(r8) :: prco_grid(pcols,pver)
  real(r8) :: prao_grid(pcols,pver)
  real(r8) :: q_ixnumliq_grid(pcols,pver)
  real(r8) :: icecldf_grid(pcols,pver)
  real(r8) :: icwnc_grid(pcols,pver)
  real(r8) :: icinc_grid(pcols,pver)
  real(r8) :: qcreso_grid(pcols,pver)
  real(r8) :: melto_grid(pcols,pver)
  real(r8) :: mnuccco_grid(pcols,pver)
  real(r8) :: mnuccto_grid(pcols,pver)
  real(r8) :: bergo_grid(pcols,pver)
  real(r8) :: homoo_grid(pcols,pver)
  real(r8) :: msacwio_grid(pcols,pver)
  real(r8) :: psacwso_grid(pcols,pver)
  real(r8) :: bergso_grid(pcols,pver)
  real(r8) :: cmeiout_grid(pcols,pver)
  real(r8) :: qireso_grid(pcols,pver)
  real(r8) :: prcio_grid(pcols,pver)
  real(r8) :: praio_grid(pcols,pver)

  real(r8),pointer :: cmeliq_grid(:,:)

  real(r8),pointer :: prec_str_grid(:)
  real(r8),pointer :: snow_str_grid(:)
  real(r8),pointer :: prec_pcw_grid(:)
  real(r8),pointer :: snow_pcw_grid(:)
  real(r8),pointer :: prec_sed_grid(:)
  real(r8),pointer :: snow_sed_grid(:)
  real(r8),pointer :: cldo_grid(:,:)
  real(r8),pointer :: nevapr_grid(:,:)
  real(r8),pointer :: prain_grid(:,:)
  real(r8),pointer :: mgflxprc_grid(:,:)
  real(r8),pointer :: mgflxsnw_grid(:,:)
  real(r8),pointer :: mgmrprc_grid(:,:)
  real(r8),pointer :: mgmrsnw_grid(:,:)
  real(r8),pointer :: cvreffliq_grid(:,:)
  real(r8),pointer :: cvreffice_grid(:,:)
  real(r8),pointer :: rate1ord_cw2pr_st_grid(:,:)
  real(r8),pointer :: wsedl_grid(:,:)
  real(r8),pointer :: CC_t_grid(:,:)
  real(r8),pointer :: CC_qv_grid(:,:)
  real(r8),pointer :: CC_ql_grid(:,:)
  real(r8),pointer :: CC_qi_grid(:,:)
  real(r8),pointer :: CC_nl_grid(:,:)
  real(r8),pointer :: CC_ni_grid(:,:)
  real(r8),pointer :: CC_qlst_grid(:,:)
  real(r8),pointer :: qme_grid(:,:)
  real(r8),pointer :: iciwpst_grid(:,:)
  real(r8),pointer :: icswp_grid(:,:)
  real(r8),pointer :: ast_grid(:,:)
  real(r8),pointer :: cldfsnow_grid(:,:)

  real(r8),pointer :: qrout_grid_ptr(:,:)
  real(r8),pointer :: qsout_grid_ptr(:,:)
  real(r8),pointer :: nrout_grid_ptr(:,:)
  real(r8),pointer :: nsout_grid_ptr(:,:)

  
  integer :: nlev   ! number of levels where cloud physics is done
  integer :: mgncol ! size of mgcols
  integer :: col_type               ! Flag to store whether accessing grid or sub-columns in pbuf_get_field
  integer, allocatable :: mgcols(:) ! Columns with microphysics performed

  logical :: use_subcol_microp

  character(128) :: errstring   ! return status (non-blank for error return)

  ! For rrtmg optics specified distribution.
  real(r8), parameter :: dcon   = 25.e-6_r8         ! Convective size distribution effective radius (meters)
  real(r8), parameter :: mucon  = 5.3_r8            ! Convective size distribution shape parameter
  real(r8), parameter :: deicon = 50._r8            ! Convective ice effective diameter (meters)

  !-------------------------------------------------------------------------------

  ! Find the number of levels used in the microphysics.
  nlev     = pver - top_lev + 1

  lchnk    = state%lchnk
  ncol     = state%ncol
  psetcols = state%psetcols
  ngrdcol  = state%ngrdcol

  itim_old = pbuf_old_tim_idx()

  call phys_getopts(use_subcol_microp_out = use_subcol_microp)

  ! Set the col_type flag to grid or subcolumn dependent on the value of use_subcol_microp
  call pbuf_col_type_index(use_subcol_microp, col_type=col_type)

  !-----------------------
  ! These physics buffer fields are read only and not set in this parameterization
  ! If these fields do not have subcolumn data, copy the grid to the subcolumn if subcolumns is turned on
  ! If subcolumns is not turned on, then these fields will be grid data

  call pbuf_get_field(pbuf, naai_idx,        naai,        col_type=col_type, copy_if_needed=use_subcol_microp)
  call pbuf_get_field(pbuf, naai_hom_idx,    naai_hom,    col_type=col_type, copy_if_needed=use_subcol_microp)
  call pbuf_get_field(pbuf, npccn_idx,       npccn,       col_type=col_type, copy_if_needed=use_subcol_microp)
  call pbuf_get_field(pbuf, rndst_idx,       rndst,       col_type=col_type, copy_if_needed=use_subcol_microp)
  call pbuf_get_field(pbuf, nacon_idx,       nacon,       col_type=col_type, copy_if_needed=use_subcol_microp)
  call pbuf_get_field(pbuf, tnd_qsnow_idx,   tnd_qsnow,   col_type=col_type, copy_if_needed=use_subcol_microp)
  call pbuf_get_field(pbuf, tnd_nsnow_idx,   tnd_nsnow,   col_type=col_type, copy_if_needed=use_subcol_microp)
  call pbuf_get_field(pbuf, re_ice_idx,      re_ice,      col_type=col_type, copy_if_needed=use_subcol_microp)
  call pbuf_get_field(pbuf, relvar_idx,      relvar,      col_type=col_type, copy_if_needed=use_subcol_microp)
  call pbuf_get_field(pbuf, accre_enhan_idx, accre_enhan, col_type=col_type, copy_if_needed=use_subcol_microp)
  call pbuf_get_field(pbuf, cmeliq_idx,      cmeliq,      col_type=col_type, copy_if_needed=use_subcol_microp)

  call pbuf_get_field(pbuf, cld_idx,         cld,     start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), &
                                                          col_type=col_type, copy_if_needed=use_subcol_microp)
  call pbuf_get_field(pbuf, concld_idx,      concld,  start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), &
                                                          col_type=col_type, copy_if_needed=use_subcol_microp)
  call pbuf_get_field(pbuf, ast_idx,         ast,     start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), &
                                                          col_type=col_type, copy_if_needed=use_subcol_microp)

  !-----------------------
  ! These physics buffer fields are calculated and set in this parameterization
  ! If subcolumns is turned on, then these fields will be calculated on a subcolumn grid, otherwise they will be a normal grid

  call pbuf_get_field(pbuf, prec_str_idx,    prec_str,    col_type=col_type)
  call pbuf_get_field(pbuf, snow_str_idx,    snow_str,    col_type=col_type)
  call pbuf_get_field(pbuf, prec_pcw_idx,    prec_pcw,    col_type=col_type)
  call pbuf_get_field(pbuf, snow_pcw_idx,    snow_pcw,    col_type=col_type)
  call pbuf_get_field(pbuf, prec_sed_idx,    prec_sed,    col_type=col_type)
  call pbuf_get_field(pbuf, snow_sed_idx,    snow_sed,    col_type=col_type)
  call pbuf_get_field(pbuf, nevapr_idx,      nevapr,      col_type=col_type)
  call pbuf_get_field(pbuf, prain_idx,       prain,       col_type=col_type)
  call pbuf_get_field(pbuf, dei_idx,         dei,         col_type=col_type)
  call pbuf_get_field(pbuf, mu_idx,          mu,          col_type=col_type)
  call pbuf_get_field(pbuf, lambdac_idx,     lambdac,     col_type=col_type)
  call pbuf_get_field(pbuf, des_idx,         des,         col_type=col_type)
  call pbuf_get_field(pbuf, ls_flxprc_idx,   mgflxprc,    col_type=col_type)
  call pbuf_get_field(pbuf, ls_flxsnw_idx,   mgflxsnw,    col_type=col_type)
  call pbuf_get_field(pbuf, ls_mrprc_idx,    mgmrprc,     col_type=col_type)
  call pbuf_get_field(pbuf, ls_mrsnw_idx,    mgmrsnw,     col_type=col_type)
  call pbuf_get_field(pbuf, cv_reffliq_idx,  cvreffliq,   col_type=col_type)
  call pbuf_get_field(pbuf, cv_reffice_idx,  cvreffice,   col_type=col_type)
  call pbuf_get_field(pbuf, iciwpst_idx,     iciwpst,     col_type=col_type)
  call pbuf_get_field(pbuf, iclwpst_idx,     iclwpst,     col_type=col_type)
  call pbuf_get_field(pbuf, icswp_idx,       icswp,       col_type=col_type)
  call pbuf_get_field(pbuf, rel_idx,         rel,         col_type=col_type)
  call pbuf_get_field(pbuf, rei_idx,         rei,         col_type=col_type)
  call pbuf_get_field(pbuf, wsedl_idx,       wsedl,       col_type=col_type)
  call pbuf_get_field(pbuf, qme_idx,         qme,         col_type=col_type)

  call pbuf_get_field(pbuf, cldo_idx,        cldo,     start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
  call pbuf_get_field(pbuf, cldfsnow_idx,    cldfsnow, start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
  call pbuf_get_field(pbuf, cc_t_idx,        CC_t,     start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type )
  call pbuf_get_field(pbuf, cc_qv_idx,       CC_qv,    start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type )
  call pbuf_get_field(pbuf, cc_ql_idx,       CC_ql,    start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type )
  call pbuf_get_field(pbuf, cc_qi_idx,       CC_qi,    start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type )
  call pbuf_get_field(pbuf, cc_nl_idx,       CC_nl,    start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type )
  call pbuf_get_field(pbuf, cc_ni_idx,       CC_ni,    start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type )
  call pbuf_get_field(pbuf, cc_qlst_idx,     CC_qlst,  start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type )

  if (rate1_cw2pr_st_idx > 0) then
     call pbuf_get_field(pbuf, rate1_cw2pr_st_idx, rate1ord_cw2pr_st, col_type=col_type)
  end if

  if (qrain_idx > 0) call pbuf_get_field(pbuf, qrain_idx, qrout_grid_ptr)
  if (qsnow_idx > 0) call pbuf_get_field(pbuf, qsnow_idx, qsout_grid_ptr)
  if (nrain_idx > 0) call pbuf_get_field(pbuf, nrain_idx, nrout_grid_ptr)
  if (nsnow_idx > 0) call pbuf_get_field(pbuf, nsnow_idx, nsout_grid_ptr)

  !-----------------------
  ! If subcolumns is turned on, all calculated fields which are on subcolumns 
  ! need to be retrieved on the grid as well for storing averaged values

  if (use_subcol_microp) then
     call pbuf_get_field(pbuf, prec_str_idx,    prec_str_grid)
     call pbuf_get_field(pbuf, snow_str_idx,    snow_str_grid)
     call pbuf_get_field(pbuf, prec_pcw_idx,    prec_pcw_grid)
     call pbuf_get_field(pbuf, snow_pcw_idx,    snow_pcw_grid)
     call pbuf_get_field(pbuf, prec_sed_idx,    prec_sed_grid)
     call pbuf_get_field(pbuf, snow_sed_idx,    snow_sed_grid)
     call pbuf_get_field(pbuf, nevapr_idx,      nevapr_grid)
     call pbuf_get_field(pbuf, prain_idx,       prain_grid)
     call pbuf_get_field(pbuf, dei_idx,         dei_grid)
     call pbuf_get_field(pbuf, mu_idx,          mu_grid)
     call pbuf_get_field(pbuf, lambdac_idx,     lambdac_grid)
     call pbuf_get_field(pbuf, des_idx,         des_grid)
     call pbuf_get_field(pbuf, ls_flxprc_idx,   mgflxprc_grid)
     call pbuf_get_field(pbuf, ls_flxsnw_idx,   mgflxsnw_grid)
     call pbuf_get_field(pbuf, ls_mrprc_idx,    mgmrprc_grid)
     call pbuf_get_field(pbuf, ls_mrsnw_idx,    mgmrsnw_grid)
     call pbuf_get_field(pbuf, cv_reffliq_idx,  cvreffliq_grid)
     call pbuf_get_field(pbuf, cv_reffice_idx,  cvreffice_grid)
     call pbuf_get_field(pbuf, iciwpst_idx,     iciwpst_grid)
     call pbuf_get_field(pbuf, iclwpst_idx,     iclwpst_grid)
     call pbuf_get_field(pbuf, icswp_idx,       icswp_grid)
     call pbuf_get_field(pbuf, rel_idx,         rel_grid)
     call pbuf_get_field(pbuf, rei_idx,         rei_grid)
     call pbuf_get_field(pbuf, wsedl_idx,       wsedl_grid)
     call pbuf_get_field(pbuf, qme_idx,         qme_grid)

     call pbuf_get_field(pbuf, cldo_idx,     cldo_grid,     start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
     call pbuf_get_field(pbuf, cldfsnow_idx, cldfsnow_grid, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
     call pbuf_get_field(pbuf, cc_t_idx,     CC_t_grid,     start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
     call pbuf_get_field(pbuf, cc_qv_idx,    CC_qv_grid,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
     call pbuf_get_field(pbuf, cc_ql_idx,    CC_ql_grid,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
     call pbuf_get_field(pbuf, cc_qi_idx,    CC_qi_grid,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
     call pbuf_get_field(pbuf, cc_nl_idx,    CC_nl_grid,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
     call pbuf_get_field(pbuf, cc_ni_idx,    CC_ni_grid,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
     call pbuf_get_field(pbuf, cc_qlst_idx,  CC_qlst_grid,  start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

     if (rate1_cw2pr_st_idx > 0) then
        call pbuf_get_field(pbuf, rate1_cw2pr_st_idx, rate1ord_cw2pr_st_grid)
     end if

  end if

  !-----------------------
  ! These are only on the grid regardless of whether subcolumns are turned on or not
  call pbuf_get_field(pbuf, ls_reffrain_idx, mgreffrain_grid)
  call pbuf_get_field(pbuf, ls_reffsnow_idx, mgreffsnow_grid)
  call pbuf_get_field(pbuf, acpr_idx,        acprecl_grid)
  call pbuf_get_field(pbuf, acgcme_idx,      acgcme_grid)
  call pbuf_get_field(pbuf, acnum_idx,       acnum_grid)
  call pbuf_get_field(pbuf, cmeliq_idx,      cmeliq_grid)
  call pbuf_get_field(pbuf, ast_idx,         ast_grid, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))


  !-------------------------------------------------------------------------------------
  ! Microphysics assumes 'liquid stratus frac = ice stratus frac
  !                      = max( liquid stratus frac, ice stratus frac )'.
  alst_mic => ast
  aist_mic => ast

  ! Output initial in-cloud LWP (before microphysics)

  iclwpi = 0._r8
  iciwpi = 0._r8

  do i = 1, ncol
     do k = top_lev, pver
        iclwpi(i) = iclwpi(i) + &
             min(state%q(i,k,ixcldliq) / max(mincld,ast(i,k)),0.005_r8) &
             * state%pdel(i,k) / gravit
        iciwpi(i) = iciwpi(i) + &
             min(state%q(i,k,ixcldice) / max(mincld,ast(i,k)),0.005_r8) &
             * state%pdel(i,k) / gravit
     end do
  end do

  cldo(:ncol,top_lev:pver)=ast(:ncol,top_lev:pver)

  ! Initialize local state from input.
  call physics_state_copy(state, state_loc)

  ! Initialize ptend for output.
  lq = .false.
  lq(1) = .true.
  lq(ixcldliq) = .true.
  lq(ixcldice) = .true.
  lq(ixnumliq) = .true.
  lq(ixnumice) = .true.

  ! the name 'cldwat' triggers special tests on cldliq
  ! and cldice in physics_update
  call physics_ptend_init(ptend, psetcols, "cldwat", ls=.true., lq=lq)

  select case (micro_mg_version)
  case (1)
     select case (micro_mg_sub_version)
     case (0)

        qc = state_loc%q(:,:,ixcldliq)
        qi = state_loc%q(:,:,ixcldice)
        nc = state_loc%q(:,:,ixnumliq)
        ni = state_loc%q(:,:,ixnumice)

        call micro_mg_tend1_0( &
             microp_uniform, psetcols, pver, ncol, top_lev, dtime, &
             state_loc%t, state_loc%q(:,:,1), qc, qi, nc,     &
             ni, state_loc%pmid, state_loc%pdel, ast, alst_mic,&
             relvar, accre_enhan,                             &
             aist_mic, rate1cld, naai, npccn,                 &
             rndst, nacon, tlat, qvlat, qcten,                &
             qiten, ncten, niten, rel, rel_fn,                &
             rei, prect, preci, nevapr, evapsnow,             &
             prain, prodsnow, cmeice, dei, mu,                &
             lambdac, qsout, des, rflx, sflx,                 &
             qrout, reff_rain, reff_snow, qcsevap, qisevap,   &
             qvres, cmeiout, vtrmc, vtrmi, qcsedten,          &
             qisedten, prao, prco, mnuccco, mnuccto,          &
             msacwio, psacwso, bergso, bergo, melto,          &
             homoo, qcreso, prcio, praio, qireso,             &
             mnuccro, pracso, meltsdt, frzrdt, mnuccdo,       &
             nrout, nsout, refl, arefl, areflz,               &
             frefl, csrfl, acsrfl, fcsrfl, rercld,            &
             ncai, ncal, qrout2, qsout2, nrout2,              &
             nsout2, drout2, dsout2, freqs, freqr,            &
             nfice, do_cldice, tnd_qsnow,                     &
             tnd_nsnow, re_ice, errstring)


     case (5)

        call micro_mg_get_cols1_5(ncol, nlev, top_lev, state%q(:,:,ixcldliq), &
             state%q(:,:,ixcldice), mgncol, mgcols)

        call micro_mg_tend1_5( &
             mgncol,   mgcols,   nlev,     top_lev,  dtime,              &
             state_loc%t,        state_loc%q(:,:,1),                     &
             state_loc%q(:,:,ixcldliq),    state_loc%q(:,:,ixcldice),    &
             state_loc%q(:,:,ixnumliq),    state_loc%q(:,:,ixnumice),    &
             relvar,             accre_enhan,                            &
             state_loc%pmid,     state_loc%pdel,     state_loc%pint,     &
             ast,                alst_mic,           aist_mic,           &
             rate1cld,           naai,     npccn,    rndst,    nacon,    &
             tlat,     qvlat,    qcten,    qiten,    ncten,    niten,    &
             rel,     rel_fn,  rei,               prect,    preci,    &
             nevapr,   evapsnow, prain,    prodsnow, cmeice,   dei,      &
             mu,       lambdac,  qsout,    des,      rflx,     sflx,     &
             qrout,              reff_rain,          reff_snow,          &
             qcsevap,  qisevap,  qvres,    cmeiout,  vtrmc,    vtrmi,    &
             qcsedten, qisedten, prao,     prco,     mnuccco,  mnuccto,  &
             msacwio,  psacwso,  bergso,   bergo,    melto,    homoo,    &
             qcreso,             prcio,    praio,    qireso,             &
             mnuccro,  pracso,   meltsdt,  frzrdt,   mnuccdo,            &
             nrout,    nsout,    refl,     arefl,    areflz,   frefl,    &
             csrfl,    acsrfl,   fcsrfl,             rercld,             &
             ncai,     ncal,     qrout2,   qsout2,   nrout2,   nsout2,   &
             drout2,   dsout2,   freqs,    freqr,    nfice,              &
             tnd_qsnow,          tnd_nsnow,          re_ice,             &
             errstring)

        call handle_errmsg(errstring, subname="micro_mg_tend1_5")
     end select
  end select

  call handle_errmsg(errstring, subname="micro_mg_tend")

  call physics_ptend_init(ptend_loc, psetcols, "micro_mg", ls=.true., lq=lq)

  ! Set local tendency.
  ptend_loc%s(:ncol,top_lev:pver)          =  tlat(:ncol,top_lev:pver)
  ptend_loc%q(:ncol,top_lev:pver,1)        = qvlat(:ncol,top_lev:pver)
  ptend_loc%q(:ncol,top_lev:pver,ixcldliq) = qcten(:ncol,top_lev:pver)
  ptend_loc%q(:ncol,top_lev:pver,ixcldice) = qiten(:ncol,top_lev:pver)
  ptend_loc%q(:ncol,top_lev:pver,ixnumliq) = ncten(:ncol,top_lev:pver)
  ptend_loc%q(:ncol,top_lev:pver,ixnumice) = niten(:ncol,top_lev:pver)

  ! Sum into overall ptend
  call physics_ptend_sum(ptend_loc, ptend, ncol)

  ! Update local state
  call physics_update(state_loc, ptend_loc, dtime)

  ! Check to make sure that the microphysics code is respecting the flags that control
  ! whether MG should be prognosing cloud ice and cloud liquid or not.
  if (.not. do_cldice) then
     if (any(ptend%q(:ncol,top_lev:pver,ixcldice) /= 0.0_r8)) &
          call endrun("micro_mg_cam:ERROR - MG microphysics is configured not to prognose cloud ice,"// &
          " but micro_mg_tend has ice mass tendencies.")
     if (any(ptend%q(:ncol,top_lev:pver,ixnumice) /= 0.0_r8)) &
          call endrun("micro_mg_cam:ERROR - MG microphysics is configured not to prognose cloud ice,"// &
          " but micro_mg_tend has ice number tendencies.")
  end if
  if (.not. do_cldliq) then
     if (any(ptend%q(:ncol,top_lev:pver,ixcldliq) /= 0.0_r8)) &
          call endrun("micro_mg_cam:ERROR - MG microphysics is configured not to prognose cloud liquid,"// &
          " but micro_mg_tend has liquid mass tendencies.")
     if (any(ptend%q(:ncol,top_lev:pver,ixnumliq) /= 0.0_r8)) &
          call endrun("micro_mg_cam:ERROR - MG microphysics is configured not to prognose cloud liquid,"// &
          " but micro_mg_tend has liquid number tendencies.")
  end if


  mnuccdohet = 0._r8
  do k=top_lev,pver
     do i=1,ncol
        if (naai(i,k) > 0._r8) then
           mnuccdohet(i,k) = mnuccdo(i,k) - (naai_hom(i,k)/naai(i,k))*mnuccdo(i,k)
        end if
     end do
  end do

  mgflxprc(:ncol,top_lev:pverp) = rflx(:ncol,top_lev:pverp) + sflx(:ncol,top_lev:pverp)
  mgflxsnw(:ncol,top_lev:pverp) = sflx(:ncol,top_lev:pverp)

  mgmrprc(:ncol,top_lev:pver) = qrout(:ncol,top_lev:pver) + qsout(:ncol,top_lev:pver)
  mgmrsnw(:ncol,top_lev:pver) = qsout(:ncol,top_lev:pver)

  !! calculate effective radius of convective liquid and ice using dcon and deicon (not used by code, not useful for COSP)
  !! hard-coded as average of hard-coded values used for deep/shallow convective detrainment (near line 1502/1505)
  cvreffliq(:ncol,top_lev:pver) = 9.0_r8
  cvreffice(:ncol,top_lev:pver) = 37.0_r8


  ! Reassign rate1 if modal aerosols
  if (rate1_cw2pr_st_idx > 0) then
     rate1ord_cw2pr_st(:ncol,top_lev:pver) = rate1cld(:ncol,top_lev:pver)
  end if

  ! Sedimentation velocity for liquid stratus cloud droplet
  wsedl(:ncol,top_lev:pver) = vtrmc(:ncol,top_lev:pver)

  ! Microphysical tendencies for use in the macrophysics at the next time step
  CC_T(:ncol,top_lev:pver)    =  tlat(:ncol,top_lev:pver)/cpair
  CC_qv(:ncol,top_lev:pver)   = qvlat(:ncol,top_lev:pver)
  CC_ql(:ncol,top_lev:pver)   = qcten(:ncol,top_lev:pver)
  CC_qi(:ncol,top_lev:pver)   = qiten(:ncol,top_lev:pver)
  CC_nl(:ncol,top_lev:pver)   = ncten(:ncol,top_lev:pver)
  CC_ni(:ncol,top_lev:pver)   = niten(:ncol,top_lev:pver)
  CC_qlst(:ncol,top_lev:pver) = qcten(:ncol,top_lev:pver)/max(0.01_r8,alst_mic(:ncol,top_lev:pver))

  ! Net micro_mg_cam condensation rate
  qme(:ncol,top_lev:pver) = cmeliq(:ncol,top_lev:pver) + cmeiout(:ncol,top_lev:pver)

  ! For precip, accumulate only total precip in prec_pcw and snow_pcw variables.
  ! Other precip output variables are set to 0
  prec_pcw(:ncol) = prect(:ncol)
  snow_pcw(:ncol) = preci(:ncol)
  prec_sed(:ncol) = 0._r8
  snow_sed(:ncol) = 0._r8
  prec_str(:ncol) = prec_pcw(:ncol) + prec_sed(:ncol)
  snow_str(:ncol) = snow_pcw(:ncol) + snow_sed(:ncol)

  icecldf(:ncol,top_lev:pver) = ast(:ncol,top_lev:pver)
  liqcldf(:ncol,top_lev:pver) = ast(:ncol,top_lev:pver)


  ! ------------------------------------------------------------ !
  ! Compute in cloud ice and liquid mixing ratios                !
  ! Note that 'iclwp, iciwp' are used for radiation computation. !
  ! ------------------------------------------------------------ !


  icinc = 0._r8
  icwnc = 0._r8
  iciwpst = 0._r8
  iclwpst = 0._r8
  icswp = 0._r8
  cldfsnow = 0._r8

  do k = top_lev, pver
     do i = 1, ncol
        ! Limits for in-cloud mixing ratios consistent with MG microphysics
        ! in-cloud mixing ratio maximum limit of 0.005 kg/kg
        icimrst(i,k)   = min( state_loc%q(i,k,ixcldice) / max(mincld,icecldf(i,k)),0.005_r8 )
        icwmrst(i,k)   = min( state_loc%q(i,k,ixcldliq) / max(mincld,liqcldf(i,k)),0.005_r8 )
        icinc(i,k)     = state_loc%q(i,k,ixnumice) / max(mincld,icecldf(i,k)) * &
             state_loc%pmid(i,k) / (287.15_r8*state_loc%t(i,k))
        icwnc(i,k)     = state_loc%q(i,k,ixnumliq) / max(mincld,liqcldf(i,k)) * &
             state_loc%pmid(i,k) / (287.15_r8*state_loc%t(i,k))
        ! Calculate micro_mg_cam cloud water paths in each layer
        ! Note: uses stratiform cloud fraction!
        iciwpst(i,k)   = min(state_loc%q(i,k,ixcldice)/max(mincld,ast(i,k)),0.005_r8) * state_loc%pdel(i,k) / gravit
        iclwpst(i,k)   = min(state_loc%q(i,k,ixcldliq)/max(mincld,ast(i,k)),0.005_r8) * state_loc%pdel(i,k) / gravit

        ! ------------------------------ !
        ! Adjust cloud fraction for snow !
        ! ------------------------------ !
        cldfsnow(i,k) = cld(i,k)
        ! If cloud and only ice ( no convective cloud or ice ), then set to 0.
        if( ( cldfsnow(i,k) .gt. 1.e-4_r8 ) .and. &
             ( concld(i,k)   .lt. 1.e-4_r8 ) .and. &
             ( state_loc%q(i,k,ixcldliq) .lt. 1.e-10_r8 ) ) then
           cldfsnow(i,k) = 0._r8
        end if
        ! If no cloud and snow, then set to 0.25
        if( ( cldfsnow(i,k) .lt. 1.e-4_r8 ) .and. ( qsout(i,k) .gt. 1.e-6_r8 ) ) then
           cldfsnow(i,k) = 0.25_r8
        end if
        ! Calculate in-cloud snow water path
        icswp(i,k) = qsout(i,k) / max( mincld, cldfsnow(i,k) ) * state_loc%pdel(i,k) / gravit
     end do
  end do


  ! ------------------------------------------------------ !
  ! ------------------------------------------------------ !
  ! All code  from here to the end is on grid columns only !
  ! ------------------------------------------------------ !
  ! ------------------------------------------------------ !
  
  ! Average the fields which are needed later in this paramterization to be on the grid
  if (use_subcol_microp) then
     call subcol_field_avg(lambdac,   ngrdcol, lchnk, lambdac_grid)
     call subcol_field_avg(mu,        ngrdcol, lchnk, mu_grid)
     call subcol_field_avg(rel,       ngrdcol, lchnk, rel_grid)
     call subcol_field_avg(rei,       ngrdcol, lchnk, rei_grid)
     call subcol_field_avg(dei,       ngrdcol, lchnk, dei_grid)
     call subcol_field_avg(prec_str,  ngrdcol, lchnk, prec_str_grid)
     call subcol_field_avg(iclwpst,   ngrdcol, lchnk, iclwpst_grid)
     call subcol_field_avg(cvreffliq, ngrdcol, lchnk, cvreffliq_grid)
     call subcol_field_avg(cvreffice, ngrdcol, lchnk, cvreffice_grid)
     call subcol_field_avg(mgflxprc,  ngrdcol, lchnk, mgflxprc_grid)
     call subcol_field_avg(mgflxsnw,  ngrdcol, lchnk, mgflxsnw_grid)
     call subcol_field_avg(qme,       ngrdcol, lchnk, qme_grid)
     call subcol_field_avg(nevapr,    ngrdcol, lchnk, nevapr_grid)
     call subcol_field_avg(prain,     ngrdcol, lchnk, prain_grid)

     ! Average fields which are not in pbuf
     call subcol_field_avg(qrout,     ngrdcol, lchnk, qrout_grid)
     call subcol_field_avg(qsout,     ngrdcol, lchnk, qsout_grid)
     call subcol_field_avg(nsout,     ngrdcol, lchnk, nsout_grid)
     call subcol_field_avg(nrout,     ngrdcol, lchnk, nrout_grid)
     call subcol_field_avg(cld,       ngrdcol, lchnk, cld_grid)
     call subcol_field_avg(qcreso,    ngrdcol, lchnk, qcreso_grid)
     call subcol_field_avg(melto,     ngrdcol, lchnk, melto_grid)
     call subcol_field_avg(mnuccco,   ngrdcol, lchnk, mnuccco_grid)
     call subcol_field_avg(mnuccto,   ngrdcol, lchnk, mnuccto_grid)
     call subcol_field_avg(bergo,     ngrdcol, lchnk, bergo_grid)
     call subcol_field_avg(homoo,     ngrdcol, lchnk, homoo_grid)
     call subcol_field_avg(msacwio,   ngrdcol, lchnk, msacwio_grid)
     call subcol_field_avg(psacwso,   ngrdcol, lchnk, psacwso_grid)
     call subcol_field_avg(bergso,    ngrdcol, lchnk, bergso_grid)
     call subcol_field_avg(cmeiout,   ngrdcol, lchnk, cmeiout_grid)
     call subcol_field_avg(qireso,    ngrdcol, lchnk, qireso_grid)
     call subcol_field_avg(prcio,     ngrdcol, lchnk, prcio_grid)
     call subcol_field_avg(praio,     ngrdcol, lchnk, praio_grid)
     call subcol_field_avg(icwmrst,   ngrdcol, lchnk, icwmrst_grid)
     call subcol_field_avg(icimrst,   ngrdcol, lchnk, icimrst_grid)
     call subcol_field_avg(liqcldf,   ngrdcol, lchnk, liqcldf_grid)
     call subcol_field_avg(icecldf,   ngrdcol, lchnk, icecldf_grid)
     call subcol_field_avg(icwnc,     ngrdcol, lchnk, icwnc_grid)
     call subcol_field_avg(icinc,     ngrdcol, lchnk, icinc_grid)
     call subcol_field_avg(state_loc%pdel,            ngrdcol, lchnk, pdel_grid)
     call subcol_field_avg(state_loc%q(:,:,ixnumliq), ngrdcol, lchnk, q_ixnumliq_grid)
     call subcol_field_avg(prao,      ngrdcol, lchnk, prao_grid)
     call subcol_field_avg(prco,      ngrdcol, lchnk, prco_grid)

  else  ! fields already on grids, so just assign
     lambdac_grid    => lambdac
     mu_grid         => mu
     rel_grid        => rel
     rei_grid        => rei
     dei_grid        => dei
     prec_str_grid   => prec_str
     iclwpst_grid    => iclwpst
     cvreffliq_grid  => cvreffliq
     cvreffice_grid  => cvreffice
     mgflxprc_grid   => mgflxprc
     mgflxsnw_grid   => mgflxsnw
     qme_grid        => qme
     nevapr_grid     => nevapr
     prain_grid      => prain

     ! This pbuf field needs to be assigned.  There is no corresponding subcol_field_avg
     ! as it is reset before it is used and would be a needless calculation
     des_grid        => des

     qrout_grid      = qrout
     qsout_grid      = qsout
     nsout_grid      = nsout
     nrout_grid      = nrout
     cld_grid        = cld
     qcreso_grid     = qcreso
     melto_grid      = melto
     mnuccco_grid    = mnuccco
     mnuccto_grid    = mnuccto
     bergo_grid      = bergo
     homoo_grid      = homoo
     msacwio_grid    = msacwio
     psacwso_grid    = psacwso
     bergso_grid     = bergso
     cmeiout_grid    = cmeiout
     qireso_grid     = qireso
     prcio_grid      = prcio
     praio_grid      = praio
     icwmrst_grid    = icwmrst
     icimrst_grid    = icimrst
     liqcldf_grid    = liqcldf
     icecldf_grid    = icecldf
     icwnc_grid      = icwnc
     icinc_grid      = icinc
     pdel_grid       = state_loc%pdel
     q_ixnumliq_grid = state_loc%q(:,:,ixnumliq)
     prao_grid       = prao
     prco_grid       = prco

  end if

  ! If on subcolumns, average the rest of the pbuf fields which were modified on subcolumns but are not used further in 
  ! this parameterization  (no need to assign in the non-subcolumn case -- the else step)
  if (use_subcol_microp) then
     call subcol_field_avg(snow_str,    ngrdcol, lchnk, snow_str_grid)
     call subcol_field_avg(prec_pcw,    ngrdcol, lchnk, prec_pcw_grid)
     call subcol_field_avg(snow_pcw,    ngrdcol, lchnk, snow_pcw_grid)
     call subcol_field_avg(prec_sed,    ngrdcol, lchnk, prec_sed_grid)
     call subcol_field_avg(snow_sed,    ngrdcol, lchnk, snow_sed_grid)
     call subcol_field_avg(cldo,        ngrdcol, lchnk, cldo_grid)
     call subcol_field_avg(mgmrprc,     ngrdcol, lchnk, mgmrprc_grid)
     call subcol_field_avg(mgmrsnw,     ngrdcol, lchnk, mgmrsnw_grid)
     call subcol_field_avg(wsedl,       ngrdcol, lchnk, wsedl_grid)
     call subcol_field_avg(cc_t,        ngrdcol, lchnk, cc_t_grid)
     call subcol_field_avg(cc_qv,       ngrdcol, lchnk, cc_qv_grid)
     call subcol_field_avg(cc_ql,       ngrdcol, lchnk, cc_ql_grid)
     call subcol_field_avg(cc_qi,       ngrdcol, lchnk, cc_qi_grid)
     call subcol_field_avg(cc_nl,       ngrdcol, lchnk, cc_nl_grid)
     call subcol_field_avg(cc_ni,       ngrdcol, lchnk, cc_ni_grid)
     call subcol_field_avg(cc_qlst,     ngrdcol, lchnk, cc_qlst_grid)
     call subcol_field_avg(iciwpst,     ngrdcol, lchnk, iciwpst_grid)
     call subcol_field_avg(icswp,       ngrdcol, lchnk, icswp_grid)
     call subcol_field_avg(cldfsnow,    ngrdcol, lchnk, cldfsnow_grid)

     if (rate1_cw2pr_st_idx > 0) then
        call subcol_field_avg(rate1ord_cw2pr_st,    ngrdcol, lchnk, rate1ord_cw2pr_st_grid)
     end if

  end if

  ! ------------------------------------- !
  ! Size distribution calculation         !
  ! ------------------------------------- !


  ! Calculate rho (on subcolumns if turned on) for size distribution parameter calculations and average it if needed
  rho(:ncol,top_lev:) = state%pmid(:ncol,top_lev:) / &
       (rair*state%t(:ncol,top_lev:))
  if (use_subcol_microp) then
     call subcol_field_avg(rho, ngrdcol, lchnk, rho_grid)
  else
     rho_grid = rho
  end if

  ! Effective radius for cloud liquid, fixed number.
  mu_grid      = 0._r8
  lambdac_grid = 0._r8
  rel_fn_grid  = 10._r8
  ncic_grid    = 1.e8_r8

  call size_dist_param_liq(mg_liq_props, icwmrst_grid(:ngrdcol,top_lev:), &
       ncic_grid(:ngrdcol,top_lev:), rho_grid(:ngrdcol,top_lev:), &
       mu_grid(:ngrdcol,top_lev:), lambdac_grid(:ngrdcol,top_lev:))

  where (icwmrst_grid(:ngrdcol,top_lev:) > qsmall)
     rel_fn_grid(:ngrdcol,top_lev:) = &
          (mu_grid(:ngrdcol,top_lev:) + 3._r8)/ &
          lambdac_grid(:ngrdcol,top_lev:)/2._r8 * 1.e6_r8
  end where

  ! Effective radius for cloud liquid, and size parameters mu_grid and lambdac_grid.
  mu_grid      = 0._r8
  lambdac_grid = 0._r8
  rel_grid     = 10._r8

  ! Calculate ncic (on subcolumns if turned on) and average it if needed
  ncic(:ncol,top_lev:) = state_loc%q(:ncol,top_lev:,ixnumliq) / &
       max(mincld,liqcldf(:ncol,top_lev:))
  if (use_subcol_microp) then
     call subcol_field_avg(ncic, ngrdcol, lchnk, ncic_grid)
  else
     ncic_grid=ncic
  endif

  call size_dist_param_liq(mg_liq_props, icwmrst_grid(:ngrdcol,top_lev:), &
       ncic_grid(:ngrdcol,top_lev:), rho_grid(:ngrdcol,top_lev:), &
       mu_grid(:ngrdcol,top_lev:), lambdac_grid(:ngrdcol,top_lev:))

  where (icwmrst_grid(:ngrdcol,top_lev:) >= qsmall)
     rel_grid(:ngrdcol,top_lev:) = &
          (mu_grid(:ngrdcol,top_lev:) + 3._r8) / &
          lambdac_grid(:ngrdcol,top_lev:)/2._r8 * 1.e6_r8
  elsewhere
     ! Deal with the fact that size_dist_param_liq sets mu_grid to -100 wherever
     ! there is no cloud.
     mu_grid(:ngrdcol,top_lev:) = 0._r8
  end where

  ! Rain/Snow effective diameter.
  ! Note -- These five fields are calculated in micro_mg_tend but are overwritten here
  drout2_grid = 0._r8
  reff_rain_grid = 0._r8
  des_grid = 0._r8
  dsout2_grid = 0._r8
  reff_snow_grid = 0._r8

  where (qrout_grid(:ngrdcol,top_lev:) >= 1.e-7_r8)

     drout2_grid(:ngrdcol,top_lev:) = avg_diameter(qrout_grid(:ngrdcol,top_lev:), &
          nrout_grid(:ngrdcol,top_lev:) * rho_grid(:ngrdcol,top_lev:), &
          rho_grid(:ngrdcol,top_lev:), rhow)

     reff_rain_grid(:ngrdcol,top_lev:) = drout2_grid(:ngrdcol,top_lev:) * &
          1.5_r8 * 1.e6_r8

  end where

  where (qsout_grid(:ngrdcol,top_lev:) >= 1.e-7_r8)

     dsout2_grid(:ngrdcol,top_lev:) = avg_diameter( qsout_grid(:ngrdcol,top_lev:), &
          nsout_grid(:ngrdcol,top_lev:) * rho_grid(:ngrdcol,top_lev:), &
          rho_grid(:ngrdcol,top_lev:), rhosn)

     des_grid(:ngrdcol,top_lev:) = dsout2_grid(:ngrdcol,top_lev:) * 3._r8 * rhosn/rhows

     reff_snow_grid(:ngrdcol,top_lev:) = dsout2_grid(:ngrdcol,top_lev:) * &
          1.5_r8 * 1.e6_r8

  end where

  ! Effective radius and diameter for cloud ice.
  ! These must always be on the grid
  rei_grid = 25._r8

  ! Calculate niic (on subcolumns if turned on) and average it if needed
  niic(:ncol,top_lev:) = state_loc%q(:ncol,top_lev:,ixnumice) / &
       max(mincld,icecldf(:ncol,top_lev:))
  if (use_subcol_microp) then
     call subcol_field_avg(niic,    ngrdcol, lchnk, niic_grid)
  else
     niic_grid    = niic
  end if

  call size_dist_param_basic(mg_ice_props, icimrst_grid(:ngrdcol,top_lev:), &
       niic_grid(:ngrdcol,top_lev:), rei_grid(:ngrdcol,top_lev:))

  where (icimrst_grid(:ngrdcol,top_lev:) >= qsmall)
     rei_grid(:ngrdcol,top_lev:) = 1.5_r8/rei_grid(:ngrdcol,top_lev:) &
          * 1.e6_r8
  elsewhere
     rei_grid(:ngrdcol,top_lev:) = 25._r8
  end where

  dei_grid = rei_grid * rhoi/rhows * 2._r8


  ! Limiters for low cloud fraction.
  do k = top_lev, pver
     do i = 1, ngrdcol
        ! Convert snow effective diameter to microns
        des_grid(i,k) = des_grid(i,k) * 1.e6_r8
        if ( ast_grid(i,k) < 1.e-4_r8 ) then
           mu_grid(i,k) = mucon
           lambdac_grid(i,k) = (mucon + 1._r8)/dcon
           dei_grid(i,k) = deicon
        end if
     end do
  end do

  mgreffrain_grid(:ngrdcol,top_lev:pver) = reff_rain_grid(:ngrdcol,top_lev:pver)
  mgreffsnow_grid(:ngrdcol,top_lev:pver) = reff_snow_grid(:ngrdcol,top_lev:pver)


  ! ------------------------------------- !
  ! Precipitation efficiency Calculation  !
  ! ------------------------------------- !


  !-----------------------------------------------------------------------
  ! Liquid water path

  ! Compute liquid water paths, and column condensation
  tgliqwp_grid(:ngrdcol) = 0._r8
  tgcmeliq_grid(:ngrdcol) = 0._r8
  
  do k = top_lev, pver
     do i = 1, ngrdcol
        tgliqwp_grid(i)  = tgliqwp_grid(i) + iclwpst_grid(i,k)*cld_grid(i,k)

        if (cmeliq_grid(i,k) > 1.e-12_r8) then
           !convert cmeliq to right units:  kgh2o/kgair/s  *  kgair/m2  / kgh2o/m3  = m/s
           tgcmeliq_grid(i) = tgcmeliq_grid(i) + cmeliq_grid(i,k) * (pdel_grid(i,k) / gravit) / rhoh2o
        end if
     end do
  end do

  ! note: 1e-6 kgho2/kgair/s * 1000. pa / (9.81 m/s2) / 1000 kgh2o/m3 = 1e-7 m/s
  ! this is 1ppmv of h2o in 10hpa
  ! alternatively: 0.1 mm/day * 1.e-4 m/mm * 1/86400 day/s = 1.e-9

  !-----------------------------------------------------------------------
  ! precipitation efficiency calculation  (accumulate cme and precip)

  minlwp = 0.01_r8        !minimum lwp threshold (kg/m3)

  ! zero out precip efficiency and total averaged precip
  pe_grid(:ngrdcol)     = 0._r8
  tpr_grid(:ngrdcol)    = 0._r8
  pefrac_grid(:ngrdcol) = 0._r8

  ! accumulate precip and condensation
  do i = 1, ngrdcol

     acgcme_grid(i)  = acgcme_grid(i) + tgcmeliq_grid(i)
     acprecl_grid(i) = acprecl_grid(i) + prec_str_grid(i)
     acnum_grid(i)   = acnum_grid(i) + 1

     ! if LWP is zero, then 'end of cloud': calculate precip efficiency
     if (tgliqwp_grid(i) < minlwp) then
        if (acprecl_grid(i) > 5.e-8_r8) then
           tpr_grid(i) = max(acprecl_grid(i)/acnum_grid(i), 1.e-15_r8)
           if (acgcme_grid(i) > 1.e-10_r8) then
              pe_grid(i) = min(max(acprecl_grid(i)/acgcme_grid(i), 1.e-15_r8), 1.e5_r8)
              pefrac_grid(i) = 1._r8
           end if
        end if

        ! reset counters
        !        if (pe_grid(i) /= 0._r8 .and. (pe_grid(i) < 1.e-8_r8 .or. pe_grid(i) > 1.e3_r8)) then
        !           write (iulog,*) 'PE_grid:ANOMALY  pe_grid, acprecl_grid, acgcme_grid, tpr_grid, acnum_grid ',pe_grid(i),&
        !                           acprecl_grid(i), acgcme_grid(i), tpr_grid(i), acnum_grid(i)
        !        endif

        acprecl_grid(i) = 0._r8
        acgcme_grid(i)  = 0._r8
        acnum_grid(i)   = 0
     end if               ! end LWP zero conditional

     ! if never find any rain....(after 10^3 timesteps...)
     if (acnum_grid(i) > 1000) then
        acnum_grid(i)   = 0
        acprecl_grid(i) = 0._r8
        acgcme_grid(i)  = 0._r8
     end if

  end do

  !-----------------------------------------------------------------------
  ! vertical average of non-zero accretion, autoconversion and ratio.
  ! vars: vprco_grid(i),vprao_grid(i),racau_grid(i),cnt_grid

  vprao_grid = 0._r8
  cnt_grid = 0
  do k = top_lev, pver
     vprao_grid(:ngrdcol) = vprao_grid(:ngrdcol) + prao_grid(:ngrdcol,k)
     where (prao_grid(:ngrdcol,k) /= 0._r8) cnt_grid(:ngrdcol) = cnt_grid(:ngrdcol) + 1
  end do

  where (cnt_grid > 0) vprao_grid = vprao_grid/cnt_grid

  vprco_grid = 0._r8
  cnt_grid = 0
  do k = top_lev, pver
     vprco_grid(:ngrdcol) = vprco_grid(:ngrdcol) + prco_grid(:ngrdcol,k)
     where (prco_grid(:ngrdcol,k) /= 0._r8) cnt_grid(:ngrdcol) = cnt_grid(:ngrdcol) + 1
  end do

  where (cnt_grid > 0)
     vprco_grid = vprco_grid/cnt_grid
     racau_grid = vprao_grid/vprco_grid
  elsewhere
     racau_grid = 0._r8
  end where

  racau_grid = min(racau_grid, 1.e10_r8)

  ! --------------------- !
  ! History Output Fields !
  ! --------------------- !

  ! Column droplet concentration
  cdnumc_grid(:ngrdcol) = sum(q_ixnumliq_grid(:ngrdcol,top_lev:pver) * &
       pdel_grid(:ngrdcol,top_lev:pver)/gravit, dim=2)

  ! Averaging for new output fields
  efcout_grid      = 0._r8
  efiout_grid      = 0._r8
  ncout_grid       = 0._r8
  niout_grid       = 0._r8
  freql_grid       = 0._r8
  freqi_grid       = 0._r8
  liqcldf_grid_out = 0._r8
  icecldf_grid_out = 0._r8
  icwmrst_grid_out = 0._r8
  icimrst_grid_out = 0._r8

  do k = top_lev, pver
     do i = 1, ngrdcol
        if ( liqcldf_grid(i,k) > 0.01_r8 .and. icwmrst_grid(i,k) > 5.e-5_r8 ) then
           efcout_grid(i,k) = rel_grid(i,k) * liqcldf_grid(i,k)
           ncout_grid(i,k)  = icwnc_grid(i,k) * liqcldf_grid(i,k)
           freql_grid(i,k)  = liqcldf_grid(i,k)
           liqcldf_grid_out(i,k) = liqcldf_grid(i,k)
           icwmrst_grid_out(i,k) = icwmrst_grid(i,k)
        end if
        if ( icecldf_grid(i,k) > 0.01_r8 .and. icimrst_grid(i,k) > 1.e-6_r8 ) then
           efiout_grid(i,k) = rei_grid(i,k) * icecldf_grid(i,k)
           niout_grid(i,k)  = icinc_grid(i,k) * icecldf_grid(i,k)
           freqi_grid(i,k)  = icecldf_grid(i,k)
           icecldf_grid_out(i,k) = icecldf_grid(i,k)
           icimrst_grid_out(i,k) = icimrst_grid(i,k)
        end if
     end do
  end do

  ! Cloud top effective radius and number.
  fcti_grid  = 0._r8
  fctl_grid  = 0._r8
  ctrel_grid = 0._r8
  ctrei_grid = 0._r8
  ctnl_grid  = 0._r8
  ctni_grid  = 0._r8
  do i = 1, ngrdcol
     do k = top_lev, pver
        if ( liqcldf_grid(i,k) > 0.01_r8 .and. icwmrst_grid(i,k) > 1.e-7_r8 ) then
           ctrel_grid(i) = rel_grid(i,k) * liqcldf_grid(i,k)
           ctnl_grid(i)  = icwnc_grid(i,k) * liqcldf_grid(i,k)
           fctl_grid(i)  = liqcldf_grid(i,k)
           exit
        end if
        if ( icecldf_grid(i,k) > 0.01_r8 .and. icimrst_grid(i,k) > 1.e-7_r8 ) then
           ctrei_grid(i) = rei_grid(i,k) * icecldf_grid(i,k)
           ctni_grid(i)  = icinc_grid(i,k) * icecldf_grid(i,k)
           fcti_grid(i)  = icecldf_grid(i,k)
           exit
        end if
     end do
  end do


  ! Assign the values to the pbuf pointers if they exist in pbuf
  if (qrain_idx > 0)  qrout_grid_ptr = qrout_grid
  if (qsnow_idx > 0)  qsout_grid_ptr = qsout_grid
  if (nrain_idx > 0)  nrout_grid_ptr = nrout_grid
  if (nsnow_idx > 0)  nsout_grid_ptr = nsout_grid

  ! --------------------------------------------- !
  ! General outfield calls for microphysics       !
  ! --------------------------------------------- !

  ! Output a handle of variables which are calculated on the fly
  ftem_grid = 0._r8

  ftem_grid(:ngrdcol,top_lev:pver) =  qcreso_grid(:ngrdcol,top_lev:pver)
  call outfld( 'MPDW2V', ftem_grid, pcols, lchnk)

  ftem_grid(:ngrdcol,top_lev:pver) =  melto_grid(:ngrdcol,top_lev:pver) - mnuccco_grid(:ngrdcol,top_lev:pver)&
        - mnuccto_grid(:ngrdcol,top_lev:pver) -  bergo_grid(:ngrdcol,top_lev:pver) - homoo_grid(:ngrdcol,top_lev:pver)&
        - msacwio_grid(:ngrdcol,top_lev:pver)
  call outfld( 'MPDW2I', ftem_grid, pcols, lchnk)

  ftem_grid(:ngrdcol,top_lev:pver) = -prao_grid(:ngrdcol,top_lev:pver) - prco_grid(:ngrdcol,top_lev:pver)&
        - psacwso_grid(:ngrdcol,top_lev:pver) - bergso_grid(:ngrdcol,top_lev:pver)
  call outfld( 'MPDW2P', ftem_grid, pcols, lchnk)

  ftem_grid(:ngrdcol,top_lev:pver) =  cmeiout_grid(:ngrdcol,top_lev:pver) + qireso_grid(:ngrdcol,top_lev:pver)
  call outfld( 'MPDI2V', ftem_grid, pcols, lchnk)

  ftem_grid(:ngrdcol,top_lev:pver) = -melto_grid(:ngrdcol,top_lev:pver) + mnuccco_grid(:ngrdcol,top_lev:pver) &
       + mnuccto_grid(:ngrdcol,top_lev:pver) +  bergo_grid(:ngrdcol,top_lev:pver) + homoo_grid(:ngrdcol,top_lev:pver)&
       + msacwio_grid(:ngrdcol,top_lev:pver)
  call outfld( 'MPDI2W', ftem_grid, pcols, lchnk)

  ftem_grid(:ngrdcol,top_lev:pver) = -prcio_grid(:ngrdcol,top_lev:pver) - praio_grid(:ngrdcol,top_lev:pver)
  call outfld( 'MPDI2P', ftem_grid, pcols, lchnk)

  ! Output fields which have not been averaged already, averaging if use_subcol_microp is true
  call outfld('MPICLWPI',    iclwpi,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('MPICIWPI',    iciwpi,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('REFL',        refl,        psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('AREFL',       arefl,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('AREFLZ',      areflz,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('FREFL',       frefl,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('CSRFL',       csrfl,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('ACSRFL',      acsrfl,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('FCSRFL',      fcsrfl,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('RERCLD',      rercld,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('NCAL',        ncal,        psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('NCAI',        ncai,        psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('AQRAIN',      qrout2,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('AQSNOW',      qsout2,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('ANRAIN',      nrout2,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('ANSNOW',      nsout2,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('FREQR',       freqr,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('FREQS',       freqs,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('MPDT',        tlat,        psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('MPDQ',        qvlat,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('MPDLIQ',      qcten,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('MPDICE',      qiten,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('EVAPSNOW',    evapsnow,    psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('QCSEVAP',     qcsevap,     psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('QISEVAP',     qisevap,     psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('QVRES',       qvres,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('VTRMC',       vtrmc,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('VTRMI',       vtrmi,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('QCSEDTEN',    qcsedten,    psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('QISEDTEN',    qisedten,    psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('MNUCCDO',     mnuccdo,     psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('MNUCCDOhet',  mnuccdohet,  psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('MNUCCRO',     mnuccro,     psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('PRACSO',      pracso ,     psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('MELTSDT',     meltsdt,     psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('FRZRDT',      frzrdt ,     psetcols, lchnk, avg_subcol_field=use_subcol_microp)
  call outfld('FICE',        nfice,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)

  ! Example subcolumn outfld call
  if (use_subcol_microp) then
     call outfld('FICE_SCOL',   nfice,       psubcols*pcols, lchnk)
  end if


  ! Output fields which are already on the grid
  call outfld('QRAIN',       qrout_grid,       pcols, lchnk)
  call outfld('QSNOW',       qsout_grid,       pcols, lchnk)
  call outfld('NRAIN',       nrout_grid,       pcols, lchnk)
  call outfld('NSNOW',       nsout_grid,       pcols, lchnk)
  call outfld('CV_REFFLIQ',  cvreffliq_grid,   pcols, lchnk)
  call outfld('CV_REFFICE',  cvreffice_grid,   pcols, lchnk)
  call outfld('LS_FLXPRC',   mgflxprc_grid,    pcols, lchnk)
  call outfld('LS_FLXSNW',   mgflxsnw_grid,    pcols, lchnk)
  call outfld('CME',         qme_grid,         pcols, lchnk)
  call outfld('PRODPREC',    prain_grid,       pcols, lchnk)
  call outfld('EVAPPREC',    nevapr_grid,      pcols, lchnk)
  call outfld('QCRESO',      qcreso_grid,      pcols, lchnk)
  call outfld('LS_REFFRAIN', mgreffrain_grid,  pcols, lchnk)
  call outfld('LS_REFFSNOW', mgreffsnow_grid,  pcols, lchnk)
  call outfld('DSNOW',       des_grid,         pcols, lchnk)
  call outfld('ADRAIN',      drout2_grid,      pcols, lchnk)
  call outfld('ADSNOW',      dsout2_grid,      pcols, lchnk)
  call outfld('PE',          pe_grid,          pcols, lchnk)
  call outfld('PEFRAC',      pefrac_grid,      pcols, lchnk)
  call outfld('APRL',        tpr_grid,         pcols, lchnk)
  call outfld('VPRAO',       vprao_grid,       pcols, lchnk)
  call outfld('VPRCO',       vprco_grid,       pcols, lchnk)
  call outfld('RACAU',       racau_grid,       pcols, lchnk)
  call outfld('AREL',        efcout_grid,      pcols, lchnk)
  call outfld('AREI',        efiout_grid,      pcols, lchnk)
  call outfld('AWNC' ,       ncout_grid,       pcols, lchnk)
  call outfld('AWNI' ,       niout_grid,       pcols, lchnk)
  call outfld('FREQL',       freql_grid,       pcols, lchnk)
  call outfld('FREQI',       freqi_grid,       pcols, lchnk)
  call outfld('ACTREL',      ctrel_grid,       pcols, lchnk)
  call outfld('ACTREI',      ctrei_grid,       pcols, lchnk)
  call outfld('ACTNL',       ctnl_grid,        pcols, lchnk)
  call outfld('ACTNI',       ctni_grid,        pcols, lchnk)
  call outfld('FCTL',        fctl_grid,        pcols, lchnk)
  call outfld('FCTI',        fcti_grid,        pcols, lchnk)
  call outfld('ICINC',       icinc_grid,       pcols, lchnk)
  call outfld('ICWNC',       icwnc_grid,       pcols, lchnk)
  call outfld('EFFLIQ_IND',  rel_fn_grid,      pcols, lchnk)
  call outfld('CDNUMC',      cdnumc_grid,      pcols, lchnk)
  call outfld('REL',         rel_grid,         pcols, lchnk)
  call outfld('REI',         rei_grid,         pcols, lchnk)
  call outfld('ICIMRST',     icimrst_grid_out, pcols, lchnk)
  call outfld('ICWMRST',     icwmrst_grid_out, pcols, lchnk)
  call outfld('CMEIOUT',     cmeiout_grid,     pcols, lchnk)
  call outfld('PRAO',        prao_grid,        pcols, lchnk)
  call outfld('PRCO',        prco_grid,        pcols, lchnk)
  call outfld('MNUCCCO',     mnuccco_grid,     pcols, lchnk)
  call outfld('MNUCCTO',     mnuccto_grid,     pcols, lchnk)
  call outfld('MSACWIO',     msacwio_grid,     pcols, lchnk)
  call outfld('PSACWSO',     psacwso_grid,     pcols, lchnk)
  call outfld('BERGSO',      bergso_grid,      pcols, lchnk)
  call outfld('BERGO',       bergo_grid,       pcols, lchnk)
  call outfld('MELTO',       melto_grid,       pcols, lchnk)
  call outfld('HOMOO',       homoo_grid,       pcols, lchnk)
  call outfld('PRCIO',       prcio_grid,       pcols, lchnk)
  call outfld('PRAIO',       praio_grid,       pcols, lchnk)
  call outfld('QIRESO',      qireso_grid,      pcols, lchnk)

  ! ptend_loc is deallocated in physics_update above
  call physics_state_dealloc(state_loc)

end subroutine micro_mg_cam_tend

end module micro_mg_cam
