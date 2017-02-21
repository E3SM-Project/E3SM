module micro_mg_cam

!---------------------------------------------------------------------------------
!
!  CAM Interfaces for MG microphysics
!
!---------------------------------------------------------------------------------
!
! How to add new packed MG inputs to micro_mg_cam_tend:
!
! If you have an input with first dimension [psetcols, pver], the procedure
! for adding inputs is as follows:
!
! 1) In addition to any variables you need to declare for the "unpacked"
!    (CAM format) version, you must declare an allocatable or pointer array
!    for the "packed" (MG format) version.
!
! 2) After micro_mg_get_cols is called, allocate the "packed" array with
!    size [mgncol, nlev].
!
! 3) Add a call similar to the following line (look before the
!    micro_mg_tend calls to see similar lines):
!
!      packed_array = packer%pack(original_array)
!
!    The packed array can then be passed into any of the MG schemes.
!
! This same procedure will also work for 1D arrays of size psetcols, 3-D
! arrays with psetcols and pver as the first dimensions, and for arrays of
! dimension [psetcols, pverp]. You only have to modify the allocation of
! the packed array before the "pack" call.
!
!---------------------------------------------------------------------------------
!
! How to add new packed MG outputs to micro_mg_cam_tend:
!
! 1) As with inputs, in addition to the unpacked outputs you must declare
!    an allocatable or pointer array for packed data. The unpacked and
!    packed arrays must *also* be targets or pointers (but cannot be both).
!
! 2) Again as for inputs, allocate the packed array using mgncol and nlev,
!    which are set in micro_mg_get_cols.
!
! 3) Add the field to post-processing as in the following line (again,
!    there are many examples before the micro_mg_tend calls):
!
!      call post_proc%add_field(p(final_array),p(packed_array))
!
!    This registers the field for post-MG averaging, and to scatter to the
!    final, unpacked version of the array.
!
!    By default, any columns/levels that are not operated on by MG will be
!    set to 0 on output; this value can be adjusted using the "fillvalue"
!    optional argument to post_proc%add_field.
!
!    Also by default, outputs from multiple substeps will be averaged after
!    MG's substepping is complete. Passing the optional argument
!    "accum_method=accum_null" will change this behavior so that the last
!    substep is always output.
!
! This procedure works on 1-D and 2-D outputs. Note that the final,
! unpacked arrays are not set until the call to
! "post_proc%process_and_unpack", which sets every single field that was
! added with post_proc%add_field.
!
!---------------------------------------------------------------------------------

use shr_kind_mod,   only: r8=>shr_kind_r8
use spmd_utils,     only: masterproc
use ppgrid,         only: pcols, pver, pverp, psubcols
use physconst,      only: gravit, rair, tmelt, cpair, rh2o, rhoh2o, &
                          latvap, latice, mwh2o, mwdry
use phys_control,   only: phys_getopts, use_hetfrz_classnuc


use physics_types,  only: physics_state, physics_ptend, &
                          physics_ptend_init, physics_state_copy, &
                          physics_update, physics_state_dealloc, &
                          physics_ptend_sum, physics_ptend_scale

use physics_buffer, only: physics_buffer_desc, pbuf_add_field, dyn_time_lvls, &
                          pbuf_old_tim_idx, pbuf_get_index, dtype_r8, dtype_i4, &
                          pbuf_get_field, pbuf_set_field, col_type_subcol, &
                          pbuf_register_subcol
use constituents,   only: cnst_add, cnst_get_ind, &
                          cnst_name, cnst_longname, sflxnam, apcnst, bpcnst, pcnst

use cldfrc2m,       only: rhmini=>rhmini_const

use cam_history,    only: addfld, horiz_only, add_default, outfld

use cam_logfile,    only: iulog
use cam_abortutils, only: endrun
use error_messages, only: handle_errmsg
use ref_pres,       only: top_lev=>trop_cloud_top_lev

use subcol_utils,   only: subcol_get_scheme
use perf_mod,       only: t_startf, t_stopf

implicit none
private
save

public :: &
   micro_mg_cam_readnl,          &
   micro_mg_cam_register,        &
   micro_mg_cam_init_cnst,       &
   micro_mg_cam_implements_cnst, &
   micro_mg_cam_init,            &
   micro_mg_cam_tend,            &
   micro_mg_dcs_tdep,            &
   micro_mg_version

integer :: micro_mg_version     = 1      ! Version number for MG.
integer :: micro_mg_sub_version = 0      ! Second part of version number.

real(r8) :: micro_mg_dcs = -1._r8

logical :: microp_uniform

!!== KZ_DCS
logical :: micro_mg_dcs_tdep  = .false.! if set to true, use temperature dependent dcs
!!== KZ_DCS


character(len=16) :: micro_mg_precip_frac_method = 'max_overlap' ! type of precipitation fraction method

real(r8)          :: micro_mg_berg_eff_factor    = 1.0_r8        ! berg efficiency factor

real(r8)          :: ice_sed_ai                  = 700.0_r8      ! Fall speed parameter for cloud ice

logical, public :: do_cldliq ! Prognose cldliq flag
logical, public :: do_cldice ! Prognose cldice flag
logical, public :: do_icesuper ! ice supersaturation code

integer :: num_steps ! Number of MG substeps

integer :: ncnst = 4       ! Number of constituents

character(len=8), parameter :: &      ! Constituent names
   cnst_names(8) = (/'CLDLIQ', 'CLDICE','NUMLIQ','NUMICE', &
                     'RAINQM', 'SNOWQM','NUMRAI','NUMSNO'/)

integer :: &
   ixcldliq = -1,      &! cloud liquid amount index
   ixcldice = -1,      &! cloud ice amount index
   ixnumliq = -1,      &! cloud liquid number index
   ixnumice = -1,      &! cloud ice water index
   ixrain = -1,        &! rain index
   ixsnow = -1,        &! snow index
   ixnumrain = -1,     &! rain number index
   ixnumsnow = -1       ! snow number index

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
   prer_evap_idx,            &
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

! Fields for UNICON
integer :: &
     am_evp_st_idx,      &! Evaporation area of stratiform precipitation
     evprain_st_idx,     &! Evaporation rate of stratiform rain [kg/kg/s]. >= 0.
     evpsnow_st_idx       ! Evaporation rate of stratiform snow [kg/kg/s]. >= 0.

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
integer :: &
     tnd_qsnow_idx = -1, &
     tnd_nsnow_idx = -1, &
     re_ice_idx = -1

! Index fields for precipitation efficiency.
integer :: &
     acpr_idx = -1, &
     acgcme_idx = -1, &
     acnum_idx = -1

! Physics buffer indices for fields registered by other modules
integer :: &
   ast_idx = -1,            &
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

! pbuf fields for heterogeneous freezing
integer :: &
   frzimm_idx = -1, &
   frzcnt_idx = -1, &
   frzdep_idx = -1

   logical :: allow_sed_supersat  ! allow supersaturated conditions after sedimentation loop
   real(r8) :: micro_mg_accre_enhan_fac = huge(1.0_r8) !Accretion enhancement factor from namelist
   real(r8) :: prc_coef1_in             = huge(1.0_r8)
   real(r8) :: prc_exp_in               = huge(1.0_r8)
   real(r8) :: prc_exp1_in              = huge(1.0_r8)
   real(r8) :: cld_sed_in               = huge(1.0_r8) !scale fac for cloud sedimentation velocity
   logical  :: mg_prc_coeff_fix_in      = .false. !temporary variable to maintain BFB, MUST be removed
   logical  :: rrtmg_temp_fix           = .false. !temporary variable to maintain BFB, MUST be removed

interface p
   module procedure p1
   module procedure p2
end interface p


!===============================================================================
contains
!===============================================================================

subroutine micro_mg_cam_readnl(nlfile)

  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use mpishorthand

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  ! Namelist variables
  logical :: micro_mg_do_cldice = .true. ! do_cldice = .true., MG microphysics is prognosing cldice
  logical :: micro_mg_do_cldliq = .true. ! do_cldliq = .true., MG microphysics is prognosing cldliq
  integer :: micro_mg_num_steps = 1      ! Number of substepping iterations done by MG (1.5 only for now).
  logical :: ice_supersat         = .false.! ice supersaturation code


  ! Local variables
  integer :: unitn, ierr
  character(len=*), parameter :: subname = 'micro_mg_cam_readnl'

  namelist /micro_mg_nl/ micro_mg_version, micro_mg_sub_version, &
       micro_mg_do_cldice, micro_mg_do_cldliq, micro_mg_num_steps, ice_sed_ai,&
!!== KZ_DCS
       micro_mg_dcs_tdep, & 
!!== KZ_DCS
       microp_uniform, micro_mg_dcs, micro_mg_precip_frac_method, micro_mg_berg_eff_factor, &
       ice_supersat

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
     do_cldice = micro_mg_do_cldice
     do_cldliq = micro_mg_do_cldliq
     num_steps = micro_mg_num_steps
     do_icesuper = ice_supersat

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
     case (2)
        select case (micro_mg_sub_version)
        case(0)
           ! MG version 2.0
        case default
           call bad_version_endrun()
        end select
     case default
        call bad_version_endrun()
     end select

     if (micro_mg_dcs < 0._r8) call endrun( "micro_mg_cam_readnl: &
              &micro_mg_dcs has not been set to a valid value.")
  end if

#ifdef SPMD
  ! Broadcast namelist variables
  call mpibcast(micro_mg_version,            1, mpiint, 0, mpicom)
  call mpibcast(micro_mg_sub_version,        1, mpiint, 0, mpicom)
  call mpibcast(do_cldice,                   1, mpilog, 0, mpicom)
  call mpibcast(do_cldliq,                   1, mpilog, 0, mpicom)
  call mpibcast(micro_mg_dcs_tdep,           1, mpilog, 0, mpicom)
  call mpibcast(num_steps,                   1, mpiint, 0, mpicom)
  call mpibcast(microp_uniform,              1, mpilog, 0, mpicom)
  call mpibcast(micro_mg_dcs,                1, mpir8,  0, mpicom)
  call mpibcast(micro_mg_berg_eff_factor,    1, mpir8,  0, mpicom)
  call mpibcast(ice_sed_ai,                  1, mpir8,  0, mpicom)
  call mpibcast(micro_mg_precip_frac_method, 16, mpichar,0, mpicom)
  call mpibcast(do_icesuper,                 1, mpilog, 0, mpicom)

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

  logical :: prog_modal_aero
  logical :: use_subcol_microp  ! If true, then are using subcolumns in microphysics
  logical :: save_subcol_microp ! If true, then need to store sub-columnized fields in pbuf

  call phys_getopts(use_subcol_microp_out = use_subcol_microp, &
                    prog_modal_aero_out   = prog_modal_aero, &
                    micro_mg_accre_enhan_fac_out = micro_mg_accre_enhan_fac)

  ! Register microphysics constituents and save indices.

  call cnst_add(cnst_names(1), mwdry, cpair, 0._r8, ixcldliq, &
       longname='Grid box averaged cloud liquid amount', is_convtran1=.true.)
  call cnst_add(cnst_names(2), mwdry, cpair, 0._r8, ixcldice, &
       longname='Grid box averaged cloud ice amount', is_convtran1=.true.)

   ! The next statements should have "is_convtran1=.true.", but this would change
   ! answers for MG 1.0. Thus make an exception for that version only.
   if (micro_mg_version == 1 .and. micro_mg_sub_version == 0) then
      call cnst_add(cnst_names(3), mwh2o, cpair, 0._r8, ixnumliq, &
           longname='Grid box averaged cloud liquid number', is_convtran1=.false.)
      call cnst_add(cnst_names(4), mwh2o, cpair, 0._r8, ixnumice, &
           longname='Grid box averaged cloud ice number', is_convtran1=.false.)
   else
      call cnst_add(cnst_names(3), mwh2o, cpair, 0._r8, ixnumliq, &
           longname='Grid box averaged cloud liquid number', is_convtran1=.true.)
      call cnst_add(cnst_names(4), mwh2o, cpair, 0._r8, ixnumice, &
           longname='Grid box averaged cloud ice number', is_convtran1=.true.)
   end if

   ! Note is_convtran1 is set to .true.
   if (micro_mg_version > 1) then
      call cnst_add(cnst_names(5), mwh2o, cpair, 0._r8, ixrain, &
           longname='Grid box averaged rain amount', is_convtran1=.true.)
      call cnst_add(cnst_names(6), mwh2o, cpair, 0._r8, ixsnow, &
           longname='Grid box averaged snow amount', is_convtran1=.true.)
      call cnst_add(cnst_names(7), mwh2o, cpair, 0._r8, ixnumrain, &
           longname='Grid box averaged rain number', is_convtran1=.true.)
      call cnst_add(cnst_names(8), mwh2o, cpair, 0._r8, ixnumsnow, &
           longname='Grid box averaged snow number', is_convtran1=.true.)
   end if

  ! Request physics buffer space for fields that persist across timesteps.

  call pbuf_add_field('CLDO','global',dtype_r8,(/pcols,pver,dyn_time_lvls/), cldo_idx)

  ! Physics buffer variables for convective cloud properties.

  call pbuf_add_field('QME',        'physpkg',dtype_r8,(/pcols,pver/), qme_idx)
  call pbuf_add_field('PRAIN',      'physpkg',dtype_r8,(/pcols,pver/), prain_idx)
  call pbuf_add_field('NEVAPR',     'physpkg',dtype_r8,(/pcols,pver/), nevapr_idx)
  call pbuf_add_field('PRER_EVAP',  'global', dtype_r8,(/pcols,pver/), prer_evap_idx)

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

  ! Fields for UNICON
  call pbuf_add_field('am_evp_st',  'global', dtype_r8, (/pcols,pver/), am_evp_st_idx)
  call pbuf_add_field('evprain_st', 'global', dtype_r8, (/pcols,pver/), evprain_st_idx)
  call pbuf_add_field('evpsnow_st', 'global', dtype_r8, (/pcols,pver/), evpsnow_st_idx)

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
    call pbuf_register_subcol('PRER_EVAP',   'micro_mg_cam_register', prer_evap_idx)

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
   use time_manager,   only: is_first_step
   use micro_mg_utils, only: micro_mg_utils_init
   use micro_mg1_0, only: micro_mg_init1_0 => micro_mg_init
   use micro_mg1_5, only: micro_mg_init1_5 => micro_mg_init
   use micro_mg2_0, only: micro_mg_init2_0 => micro_mg_init

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
   logical :: do_clubb_sgs
   integer :: budget_histfile      ! output history file number for budget fields
   integer :: ierr
   character(128) :: errstring     ! return status (non-blank for error return)

   !-----------------------------------------------------------------------

   call phys_getopts(use_subcol_microp_out= use_subcol_microp, &
                     do_clubb_sgs_out     = do_clubb_sgs,      &
                     prc_coef1_out        = prc_coef1_in,      &
                     prc_exp_out          = prc_exp_in,        &
                     prc_exp1_out         = prc_exp1_in,       &
                     cld_sed_out          = cld_sed_in,        &
                     mg_prc_coeff_fix_out = mg_prc_coeff_fix_in, &
                     rrtmg_temp_fix_out   = rrtmg_temp_fix       )

   if (do_clubb_sgs) then
     allow_sed_supersat = .false.
   else
     allow_sed_supersat = .true.
   endif

   if (masterproc) then
      write(iulog,"(A,I2,A,I2)") "Initializing MG version ",micro_mg_version,".",micro_mg_sub_version
      if (.not. do_cldliq) &
           write(iulog,*) "MG prognostic cloud liquid has been turned off via namelist."
      if (.not. do_cldice) &
           write(iulog,*) "MG prognostic cloud ice has been turned off via namelist."
      write(iulog,*) "Number of microphysics substeps is: ",num_steps
   end if

   select case (micro_mg_version)
   case (1)
      ! Set constituent number for later loops.
      ncnst = 4

      select case (micro_mg_sub_version)
      case (0)
         ! MG 1 does not initialize micro_mg_utils, so have to do it here.
         call micro_mg_utils_init(r8, rh2o, cpair, tmelt, latvap, latice, &
              micro_mg_dcs, ice_sed_ai, errstring)
         call handle_errmsg(errstring, subname="micro_mg_utils_init")

         call micro_mg_init1_0( &
              r8, gravit, rair, rh2o, cpair, &
              rhoh2o, tmelt, latvap, latice, &
              rhmini, micro_mg_dcs, micro_mg_dcs_tdep, use_hetfrz_classnuc, &
              micro_mg_precip_frac_method, micro_mg_berg_eff_factor, errstring)
      case (5)
         ! MG 1 does not initialize micro_mg_utils, so have to do it here.
         call micro_mg_utils_init(r8, rh2o, cpair, tmelt, latvap, latice, &
              micro_mg_dcs, ice_sed_ai, errstring)
         call handle_errmsg(errstring, subname="micro_mg_utils_init")

         call micro_mg_init1_5( &
              r8, gravit, rair, rh2o, cpair, &
              tmelt, latvap, latice, rhmini, &
              micro_mg_dcs,                  &
              micro_mg_dcs_tdep,             &
              microp_uniform, do_cldice, use_hetfrz_classnuc, &
              micro_mg_precip_frac_method, micro_mg_berg_eff_factor, errstring)
      end select
   case (2)
      ! Set constituent number for later loops.
      ncnst = 8

      select case (micro_mg_sub_version)
      case (0)
         call micro_mg_init2_0( &
              r8, gravit, rair, rh2o, cpair, &
              tmelt, latvap, latice, rhmini, &
              micro_mg_dcs,                  &
              micro_mg_dcs_tdep,             &
              microp_uniform, do_cldice, use_hetfrz_classnuc, &
              micro_mg_precip_frac_method, micro_mg_berg_eff_factor, &
              allow_sed_supersat, ice_sed_ai, prc_coef1_in,prc_exp_in, &
              prc_exp1_in, cld_sed_in, mg_prc_coeff_fix_in, errstring)
      end select
   end select

   call handle_errmsg(errstring, subname="micro_mg_init")

   ! Register history variables
   do m = 1, ncnst
      call cnst_get_ind(cnst_names(m), mm)
      if ( any(mm == (/ ixcldliq, ixcldice, ixrain, ixsnow /)) ) then
         ! mass mixing ratios
         call addfld(cnst_name(mm), (/ 'lev' /), 'A', 'kg/kg', cnst_longname(mm)                   )
         call addfld(sflxnam(mm),    horiz_only, 'A',   'kg/m2/s', trim(cnst_name(mm))//' surface flux')
      else if ( any(mm == (/ ixnumliq, ixnumice, ixnumrain, ixnumsnow /)) ) then
         ! number concentrations
         call addfld(cnst_name(mm), (/ 'lev' /), 'A', '1/kg', cnst_longname(mm)                   )
         call addfld(sflxnam(mm),    horiz_only, 'A',   '1/m2/s', trim(cnst_name(mm))//' surface flux')
      else
         call endrun( "micro_mg_cam_init: &
              &Could not call addfld for constituent with unknown units.")
      endif
   end do

   call addfld(apcnst(ixcldliq), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldliq))//' after physics'  )
   call addfld(apcnst(ixcldice), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldice))//' after physics'  )
   call addfld(bpcnst(ixcldliq), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldliq))//' before physics' )
   call addfld(bpcnst(ixcldice), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldice))//' before physics' )

   if (micro_mg_version > 1) then
      call addfld(apcnst(ixrain), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixrain))//' after physics'  )
      call addfld(apcnst(ixsnow), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixsnow))//' after physics'  )
      call addfld(bpcnst(ixrain), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixrain))//' before physics' )
      call addfld(bpcnst(ixsnow), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixsnow))//' before physics' )
   end if

   call addfld ('CME', (/ 'lev' /), 'A', 'kg/kg/s', 'Rate of cond-evap within the cloud'                      )
   call addfld ('PRODPREC', (/ 'lev' /), 'A', 'kg/kg/s', 'Rate of conversion of condensate to precip'              )
   call addfld ('EVAPPREC', (/ 'lev' /), 'A', 'kg/kg/s', 'Rate of evaporation of falling precip'                   )
   call addfld ('EVAPSNOW', (/ 'lev' /), 'A', 'kg/kg/s', 'Rate of evaporation of falling snow'                     )
   call addfld ('HPROGCLD', (/ 'lev' /), 'A', 'W/kg'    , 'Heating from prognostic clouds'                          )
   call addfld ('FICE', (/ 'lev' /), 'A', 'fraction', 'Fractional ice content within cloud'                     )
   call addfld ('ICWMRST', (/ 'lev' /), 'A', 'kg/kg', 'Prognostic in-stratus water mixing ratio'                )
   call addfld ('ICIMRST', (/ 'lev' /), 'A', 'kg/kg', 'Prognostic in-stratus ice mixing ratio'                  )

   ! MG microphysics diagnostics
   call addfld ('QCSEVAP', (/ 'lev' /), 'A', 'kg/kg/s', 'Rate of evaporation of falling cloud water'              )
   call addfld ('QISEVAP', (/ 'lev' /), 'A', 'kg/kg/s', 'Rate of sublimation of falling cloud ice'                )
   call addfld ('QVRES', (/ 'lev' /), 'A', 'kg/kg/s', 'Rate of residual condensation term'                      )
   call addfld ('CMEIOUT', (/ 'lev' /), 'A', 'kg/kg/s', 'Rate of deposition/sublimation of cloud ice'             )
   call addfld ('VTRMC', (/ 'lev' /), 'A', 'm/s', 'Mass-weighted cloud water fallspeed'                     )
   call addfld ('VTRMI', (/ 'lev' /), 'A', 'm/s', 'Mass-weighted cloud ice fallspeed'                       )
   call addfld ('QCSEDTEN', (/ 'lev' /), 'A', 'kg/kg/s', 'Cloud water mixing ratio tendency from sedimentation'    )
   call addfld ('QISEDTEN', (/ 'lev' /), 'A', 'kg/kg/s', 'Cloud ice mixing ratio tendency from sedimentation'      )
   call addfld ('PRAO', (/ 'lev' /), 'A', 'kg/kg/s', 'Accretion of cloud water by rain'                        )
   call addfld ('PRCO', (/ 'lev' /), 'A', 'kg/kg/s', 'Autoconversion of cloud water'                           )
   call addfld ('MNUCCCO', (/ 'lev' /), 'A', 'kg/kg/s', 'Immersion freezing of cloud water'                       )
   call addfld ('MNUCCTO', (/ 'lev' /), 'A', 'kg/kg/s', 'Contact freezing of cloud water'                         )
   call addfld ('MNUCCDO', (/ 'lev' /), 'A', 'kg/kg/s', 'Homogeneous and heterogeneous nucleation from vapor'     )
   call addfld ('MNUCCDOhet', (/ 'lev' /), 'A','kg/kg/s', 'Heterogeneous nucleation from vapor'                     )
   call addfld ('MSACWIO', (/ 'lev' /), 'A', 'kg/kg/s', 'Conversion of cloud water from rime-splintering'         )
   call addfld ('PSACWSO', (/ 'lev' /), 'A', 'kg/kg/s', 'Accretion of cloud water by snow'                        )
   call addfld ('BERGSO', (/ 'lev' /), 'A', 'kg/kg/s', 'Conversion of cloud water to snow from bergeron'         )
   call addfld ('BERGO', (/ 'lev' /), 'A', 'kg/kg/s', 'Conversion of cloud water to cloud ice from bergeron'    )
   call addfld ('MELTO', (/ 'lev' /), 'A', 'kg/kg/s', 'Melting of cloud ice'                                    )
   call addfld ('HOMOO', (/ 'lev' /), 'A', 'kg/kg/s', 'Homogeneous freezing of cloud water'                     )
   call addfld ('QCRESO', (/ 'lev' /), 'A', 'kg/kg/s', 'Residual condensation term for cloud water'              )
   call addfld ('PRCIO', (/ 'lev' /), 'A', 'kg/kg/s', 'Autoconversion of cloud ice'                             )
   call addfld ('PRAIO', (/ 'lev' /), 'A', 'kg/kg/s', 'Accretion of cloud ice by rain'                          )
   call addfld ('QIRESO', (/ 'lev' /), 'A', 'kg/kg/s', 'Residual deposition term for cloud ice'                  )
   call addfld ('MNUCCRO', (/ 'lev' /), 'A', 'kg/kg/s', 'Heterogeneous freezing of rain to snow'                  )
   call addfld ('PRACSO', (/ 'lev' /), 'A', 'kg/kg/s', 'Accretion of rain by snow'                               )
   call addfld ('MELTSDT', (/ 'lev' /), 'A', 'W/kg', 'Latent heating rate due to melting of snow'              )
   call addfld ('FRZRDT', (/ 'lev' /), 'A', 'W/kg', 'Latent heating rate due to homogeneous freezing of rain' )
   if (micro_mg_version > 1) then
      call addfld ('QRSEDTEN', (/ 'lev' /), 'A', 'kg/kg/s', 'Rain mixing ratio tendency from sedimentation'           )
      call addfld ('QSSEDTEN', (/ 'lev' /), 'A', 'kg/kg/s', 'Snow mixing ratio tendency from sedimentation'           )
   end if

   ! History variables for CAM5 microphysics
   call addfld ('MPDT', (/ 'lev' /), 'A', 'W/kg', 'Heating tendency - Morrison microphysics'                )
   call addfld ('MPDQ', (/ 'lev' /), 'A', 'kg/kg/s', 'Q tendency - Morrison microphysics'                      )
   call addfld ('MPDLIQ', (/ 'lev' /), 'A', 'kg/kg/s', 'CLDLIQ tendency - Morrison microphysics'                 )
   call addfld ('MPDICE', (/ 'lev' /), 'A', 'kg/kg/s', 'CLDICE tendency - Morrison microphysics'                 )
   call addfld ('MPDNLIQ',    (/ 'lev' /), 'A', '1/kg/s',   'NUMLIQ tendency - Morrison microphysics'            )
   call addfld ('MPDNICE',    (/ 'lev' /), 'A', '1/kg/s',   'NUMICE tendency - Morrison microphysics'            )
   call addfld ('MPDW2V', (/ 'lev' /), 'A', 'kg/kg/s', 'Water <--> Vapor tendency - Morrison microphysics'       )
   call addfld ('MPDW2I', (/ 'lev' /), 'A', 'kg/kg/s', 'Water <--> Ice tendency - Morrison microphysics'         )
   call addfld ('MPDW2P', (/ 'lev' /), 'A', 'kg/kg/s', 'Water <--> Precip tendency - Morrison microphysics'      )
   call addfld ('MPDI2V', (/ 'lev' /), 'A', 'kg/kg/s', 'Ice <--> Vapor tendency - Morrison microphysics'         )
   call addfld ('MPDI2W', (/ 'lev' /), 'A', 'kg/kg/s', 'Ice <--> Water tendency - Morrison microphysics'         )
   call addfld ('MPDI2P', (/ 'lev' /), 'A', 'kg/kg/s', 'Ice <--> Precip tendency - Morrison microphysics'        )
   call addfld ('ICWNC', (/ 'lev' /), 'A', 'm-3', 'Prognostic in-cloud water number conc'                   )
   call addfld ('ICINC', (/ 'lev' /), 'A', 'm-3', 'Prognostic in-cloud ice number conc'                     )
   call addfld ('EFFLIQ_IND', (/ 'lev' /), 'A','Micron', 'Prognostic droplet effective radius (indirect effect)'   )
   call addfld ('CDNUMC', horiz_only,    'A', '1/m2', 'Vertically-integrated droplet concentration'             )
   call addfld ('MPICLWPI', horiz_only,    'A', 'kg/m2', 'Vertically-integrated &
        &in-cloud Initial Liquid WP (Before Micro)' )
   call addfld ('MPICIWPI', horiz_only,    'A', 'kg/m2', 'Vertically-integrated &
        &in-cloud Initial Ice WP (Before Micro)'    )

   ! This is provided as an example on how to write out subcolumn output
   ! NOTE -- only 'I' should be used for sub-column fields as subc-columns could shift from time-step to time-step
   if (use_subcol_microp) then
      call addfld('FICE_SCOL', (/'psubcols','lev     '/), 'I', 'fraction', &
           'Sub-column fractional ice content within cloud', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('MPDICE_SCOL', (/'psubcols','lev     '/), 'I', 'kg/kg/s', &
           'Sub-column CLDICE tendency - Morrison microphysics', flag_xyfill=.true., fill_value=1.e30_r8)
      call addfld('MPDLIQ_SCOL', (/'psubcols','lev     '/), 'I', 'kg/kg/s', &
           'Sub-column CLDLIQ tendency - Morrison microphysics', flag_xyfill=.true., fill_value=1.e30_r8)
   end if

   ! Averaging for cloud particle number and size
   call addfld ('AWNC', (/ 'lev' /), 'A', 'm-3', 'Average cloud water number conc'                         )
   call addfld ('AWNI', (/ 'lev' /), 'A', 'm-3', 'Average cloud ice number conc'                           )
   call addfld ('AREL', (/ 'lev' /), 'A', 'Micron', 'Average droplet effective radius'                        )
   call addfld ('AREI', (/ 'lev' /), 'A', 'Micron', 'Average ice effective radius'                            )
   ! Frequency arrays for above
   call addfld ('FREQL', (/ 'lev' /), 'A', 'fraction', 'Fractional occurrence of liquid'                          )
   call addfld ('FREQI', (/ 'lev' /), 'A', 'fraction', 'Fractional occurrence of ice'                             )

   ! Average cloud top particle size and number (liq, ice) and frequency
   call addfld ('ACTREL', horiz_only,    'A', 'Micron', 'Average Cloud Top droplet effective radius'              )
   call addfld ('ACTREI', horiz_only,    'A', 'Micron', 'Average Cloud Top ice effective radius'                  )
   call addfld ('ACTNL', horiz_only,    'A', 'Micron', 'Average Cloud Top droplet number'                        )
   call addfld ('ACTNI', horiz_only,    'A', 'Micron', 'Average Cloud Top ice number'                            )

   call addfld ('FCTL', horiz_only,    'A', 'fraction', 'Fractional occurrence of cloud top liquid'                )
   call addfld ('FCTI', horiz_only,    'A', 'fraction', 'Fractional occurrence of cloud top ice'                   )

   call addfld ('LS_FLXPRC', (/ 'ilev' /), 'A', 'kg/m2/s', 'ls stratiform gbm interface rain+snow flux')
   call addfld ('LS_FLXSNW', (/ 'ilev' /), 'A', 'kg/m2/s', 'ls stratiform gbm interface snow flux')

   call addfld ('REL', (/ 'lev' /), 'A', 'micron', 'MG REL stratiform cloud effective radius liquid')
   call addfld ('REI', (/ 'lev' /), 'A', 'micron', 'MG REI stratiform cloud effective radius ice')
   call addfld ('LS_REFFRAIN', (/ 'lev' /), 'A', 'micron', 'ls stratiform rain effective radius')
   call addfld ('LS_REFFSNOW', (/ 'lev' /), 'A', 'micron', 'ls stratiform snow effective radius')
   call addfld ('CV_REFFLIQ', (/ 'lev' /), 'A', 'micron', 'convective cloud liq effective radius')
   call addfld ('CV_REFFICE', (/ 'lev' /), 'A', 'micron', 'convective cloud ice effective radius')

!!== KZ_DCS
   call addfld ('DCST',(/ 'lev' /), 'A','m','dcs')
!!== KZ_DCS
   ! diagnostic precip
   call addfld ('QRAIN',(/ 'lev' /), 'A','kg/kg','Diagnostic grid-mean rain mixing ratio'         )
   call addfld ('QSNOW',(/ 'lev' /), 'A','kg/kg','Diagnostic grid-mean snow mixing ratio'         )
   call addfld ('NRAIN',(/ 'lev' /), 'A','m-3','Diagnostic grid-mean rain number conc'         )
   call addfld ('NSNOW',(/ 'lev' /), 'A','m-3','Diagnostic grid-mean snow number conc'         )

   ! size of precip
   call addfld ('RERCLD',(/ 'lev' /), 'A','m','Diagnostic effective radius of Liquid Cloud and Rain' )
   call addfld ('DSNOW',(/ 'lev' /), 'A','m','Diagnostic grid-mean snow diameter'         )

   ! diagnostic radar reflectivity, cloud-averaged
   call addfld ('REFL',(/ 'lev' /), 'A','DBz','94 GHz radar reflectivity'       )
   call addfld ('AREFL',(/ 'lev' /), 'A','DBz','Average 94 GHz radar reflectivity'       )
   call addfld ('FREFL',(/ 'lev' /), 'A','fraction','Fractional occurrence of radar reflectivity'       )

   call addfld ('CSRFL',(/ 'lev' /), 'A','DBz','94 GHz radar reflectivity (CloudSat thresholds)'       )
   call addfld ('ACSRFL',(/ 'lev' /), 'A','DBz','Average 94 GHz radar reflectivity (CloudSat thresholds)'       )
   call addfld ('FCSRFL',(/ 'lev' /), 'A','fraction','Fractional occurrence of radar reflectivity (CloudSat thresholds)' &
        )

   call addfld ('AREFLZ',(/ 'lev' /), 'A','mm^6/m^3','Average 94 GHz radar reflectivity'       )

   ! Aerosol information
   call addfld ('NCAL',(/ 'lev' /), 'A','1/m3','Number Concentation Activated for Liquid')
   call addfld ('NCAI',(/ 'lev' /), 'A','1/m3','Number Concentation Activated for Ice')

   ! Average rain and snow mixing ratio (Q), number (N) and diameter (D), with frequency
   call addfld ('AQRAIN',(/ 'lev' /), 'A','kg/kg','Average rain mixing ratio'         )
   call addfld ('AQSNOW',(/ 'lev' /), 'A','kg/kg','Average snow mixing ratio'         )
   call addfld ('ANRAIN',(/ 'lev' /), 'A','m-3','Average rain number conc'         )
   call addfld ('ANSNOW',(/ 'lev' /), 'A','m-3','Average snow number conc'         )
   call addfld ('ADRAIN',(/ 'lev' /), 'A','Micron','Average rain effective Diameter'         )
   call addfld ('ADSNOW',(/ 'lev' /), 'A','Micron','Average snow effective Diameter'         )
   call addfld ('FREQR',(/ 'lev' /), 'A','fraction','Fractional occurrence of rain'       )
   call addfld ('FREQS',(/ 'lev' /), 'A','fraction','Fractional occurrence of snow'       )

   ! precipitation efficiency & other diagnostic fields
   call addfld('PE'    ,       horiz_only, 'A', '1', 'Stratiform Precipitation Efficiency  (precip/cmeliq)' )
   call addfld('APRL'  ,     horiz_only, 'A', 'm/s', 'Average Stratiform Precip Rate over efficiency calculation' )
   call addfld('PEFRAC',       horiz_only, 'A', '1', 'Fraction of timesteps precip efficiency reported' )
   call addfld('VPRCO' , horiz_only, 'A', 'kg/kg/s', 'Vertical average of autoconversion rate' )
   call addfld('VPRAO' , horiz_only, 'A', 'kg/kg/s', 'Vertical average of accretion rate' )
   call addfld('RACAU' , horiz_only, 'A', 'kg/kg/s', 'Accretion/autoconversion ratio from vertical average' )

   if (micro_mg_version > 1) then
      call addfld('UMR', (/ 'lev' /), 'A',   'm/s', 'Mass-weighted rain  fallspeed'              )
      call addfld('UMS', (/ 'lev' /), 'A',   'm/s', 'Mass-weighted snow fallspeed'               )
   end if

   ! qc limiter (only output in versions 1.5 and later)
   if (.not. (micro_mg_version == 1 .and. micro_mg_sub_version == 0)) then
      call addfld('QCRAT', (/ 'lev' /), 'A', 'fraction', 'Qc Limiter: Fraction of qc tendency applied')
   end if

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
      call add_default ('ADRAIN   ', 1, ' ')
      call add_default ('ADSNOW   ', 1, ' ')
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
      if (micro_mg_version > 1) then
         call add_default ('QRSEDTEN ', budget_histfile, ' ')
         call add_default ('QSSEDTEN ', budget_histfile, ' ')
      end if
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
      if (micro_mg_version > 1) then
         call add_default(cnst_name(ixrain), budget_histfile, ' ')
         call add_default(cnst_name(ixsnow), budget_histfile, ' ')
         call add_default(apcnst   (ixrain), budget_histfile, ' ')
         call add_default(apcnst   (ixsnow), budget_histfile, ' ')
         call add_default(bpcnst   (ixrain), budget_histfile, ' ')
         call add_default(bpcnst   (ixsnow), budget_histfile, ' ')
      end if

   end if

   ! physics buffer indices
   ast_idx      = pbuf_get_index('AST')
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

   cmeliq_idx = pbuf_get_index('CMELIQ')

   ! These fields may have been added, so don't abort if they have not been
   qrain_idx    = pbuf_get_index('QRAIN', ierr)
   qsnow_idx    = pbuf_get_index('QSNOW', ierr)
   nrain_idx    = pbuf_get_index('NRAIN', ierr)
   nsnow_idx    = pbuf_get_index('NSNOW', ierr)

  ! fields for heterogeneous freezing
  frzimm_idx = pbuf_get_index('FRZIMM', ierr)
  frzcnt_idx = pbuf_get_index('FRZCNT', ierr)
  frzdep_idx = pbuf_get_index('FRZDEP', ierr)

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
      call pbuf_set_field(pbuf2d, accre_enhan_idx, micro_mg_accre_enhan_fac)
      call pbuf_set_field(pbuf2d, am_evp_st_idx,  0._r8)
      call pbuf_set_field(pbuf2d, evprain_st_idx, 0._r8)
      call pbuf_set_field(pbuf2d, evpsnow_st_idx, 0._r8)
      call pbuf_set_field(pbuf2d, prer_evap_idx,  0._r8)

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

   use micro_mg_data, only: MGPacker, MGPostProc, accum_null, accum_mean

   use micro_mg1_0, only: micro_mg_tend1_0 => micro_mg_tend, &
        micro_mg_get_cols1_0 => micro_mg_get_cols
   use micro_mg1_5, only: micro_mg_tend1_5 => micro_mg_tend, &
        micro_mg_get_cols1_5 => micro_mg_get_cols
   use micro_mg2_0, only: micro_mg_tend2_0 => micro_mg_tend, &
        micro_mg_get_cols2_0 => micro_mg_get_cols

   use physics_buffer,  only: pbuf_col_type_index
   use subcol,          only: subcol_field_avg

   type(physics_state),         intent(in)    :: state
   type(physics_ptend),         intent(out)   :: ptend
   real(r8),                    intent(in)    :: dtime
   type(physics_buffer_desc),   pointer       :: pbuf(:)

   ! Local variables
   integer :: lchnk, ncol, psetcols, ngrdcol

   integer :: i, k, itim_old, it

   real(r8), pointer :: naai(:,:)      ! ice nucleation number
   real(r8), pointer :: naai_hom(:,:)  ! ice nucleation number (homogeneous)
   real(r8), pointer :: npccn(:,:)     ! liquid activation number tendency
   real(r8), pointer :: rndst(:,:,:)
   real(r8), pointer :: nacon(:,:,:)
   real(r8), pointer :: am_evp_st_grid(:,:)    ! Evaporation area of stratiform precipitation. 0<= am_evp_st <=1.
   real(r8), pointer :: evprain_st_grid(:,:)   ! Evaporation rate of stratiform rain [kg/kg/s]
   real(r8), pointer :: evpsnow_st_grid(:,:)   ! Evaporation rate of stratiform snow [kg/kg/s]

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
   real(r8), pointer :: prer_evap(:,:)    ! precipitation evaporation rate
   real(r8), pointer :: relvar(:,:)       ! relative variance of cloud water
   real(r8), pointer :: accre_enhan(:,:)  ! optional accretion enhancement for experimentation
   real(r8), pointer :: prain(:,:)        ! Total precipitation (rain + snow)
   real(r8), pointer :: dei(:,:)          ! Ice effective diameter (meters) (AG: microns?)
   real(r8), pointer :: mu(:,:)           ! Size distribution shape parameter for radiation
   real(r8), pointer :: lambdac(:,:)      ! Size distribution slope parameter for radiation
   real(r8), pointer :: des(:,:)          ! Snow effective diameter (m)

   real(r8) :: rho(state%psetcols,pver)
   real(r8) :: cldmax(state%psetcols,pver)

   real(r8), target :: rate1cld(state%psetcols,pver) ! array to hold rate1ord_cw2pr_st from microphysics

   real(r8), target :: tlat(state%psetcols,pver)
   real(r8), target :: qvlat(state%psetcols,pver)
   real(r8), target :: qcten(state%psetcols,pver)
   real(r8), target :: qiten(state%psetcols,pver)
   real(r8), target :: ncten(state%psetcols,pver)
   real(r8), target :: niten(state%psetcols,pver)

   real(r8), target :: qrten(state%psetcols,pver)
   real(r8), target :: qsten(state%psetcols,pver)
   real(r8), target :: nrten(state%psetcols,pver)
   real(r8), target :: nsten(state%psetcols,pver)

   real(r8), target :: prect(state%psetcols)
   real(r8), target :: preci(state%psetcols)
   real(r8), target :: am_evp_st(state%psetcols,pver)  ! Area over which precip evaporates
   real(r8), target :: evapsnow(state%psetcols,pver)   ! Local evaporation of snow
   real(r8), target :: prodsnow(state%psetcols,pver)   ! Local production of snow
   real(r8), target :: cmeice(state%psetcols,pver)     ! Rate of cond-evap of ice within the cloud
   real(r8), target :: qsout(state%psetcols,pver)      ! Snow mixing ratio
   real(r8), target :: rflx(state%psetcols,pverp)      ! grid-box average rain flux (kg m^-2 s^-1)
   real(r8), target :: sflx(state%psetcols,pverp)      ! grid-box average snow flux (kg m^-2 s^-1)
   real(r8), target :: qrout(state%psetcols,pver)      ! Rain mixing ratio
   real(r8), target :: qcsevap(state%psetcols,pver)    ! Evaporation of falling cloud water
   real(r8), target :: qisevap(state%psetcols,pver)    ! Sublimation of falling cloud ice
   real(r8), target :: qvres(state%psetcols,pver)      ! Residual condensation term to remove excess saturation
   real(r8), target :: cmeiout(state%psetcols,pver)    ! Deposition/sublimation rate of cloud ice
   real(r8), target :: vtrmc(state%psetcols,pver)      ! Mass-weighted cloud water fallspeed
   real(r8), target :: vtrmi(state%psetcols,pver)      ! Mass-weighted cloud ice fallspeed
   real(r8), target :: umr(state%psetcols,pver)        ! Mass-weighted rain fallspeed
   real(r8), target :: ums(state%psetcols,pver)        ! Mass-weighted snow fallspeed
   real(r8), target :: qcsedten(state%psetcols,pver)   ! Cloud water mixing ratio tendency from sedimentation
   real(r8), target :: qisedten(state%psetcols,pver)   ! Cloud ice mixing ratio tendency from sedimentation
   real(r8), target :: qrsedten(state%psetcols,pver)   ! Rain mixing ratio tendency from sedimentation
   real(r8), target :: qssedten(state%psetcols,pver)   ! Snow mixing ratio tendency from sedimentation

   real(r8), target :: prao(state%psetcols,pver)
   real(r8), target :: prco(state%psetcols,pver)
   real(r8), target :: mnuccco(state%psetcols,pver)
   real(r8), target :: mnuccto(state%psetcols,pver)
   real(r8), target :: msacwio(state%psetcols,pver)
   real(r8), target :: psacwso(state%psetcols,pver)
   real(r8), target :: bergso(state%psetcols,pver)
   real(r8), target :: bergo(state%psetcols,pver)
   real(r8), target :: melto(state%psetcols,pver)
   real(r8), target :: homoo(state%psetcols,pver)
   real(r8), target :: qcreso(state%psetcols,pver)
   real(r8), target :: prcio(state%psetcols,pver)
   real(r8), target :: praio(state%psetcols,pver)
   real(r8), target :: qireso(state%psetcols,pver)
   real(r8), target :: mnuccro(state%psetcols,pver)
   real(r8), target :: pracso (state%psetcols,pver)
   real(r8), target :: meltsdt(state%psetcols,pver)
   real(r8), target :: frzrdt (state%psetcols,pver)
   real(r8), target :: mnuccdo(state%psetcols,pver)
   real(r8), target :: nrout(state%psetcols,pver)
   real(r8), target :: nsout(state%psetcols,pver)
   real(r8), target :: refl(state%psetcols,pver)    ! analytic radar reflectivity
   real(r8), target :: arefl(state%psetcols,pver)   ! average reflectivity will zero points outside valid range
   real(r8), target :: areflz(state%psetcols,pver)  ! average reflectivity in z.
   real(r8), target :: frefl(state%psetcols,pver)
   real(r8), target :: csrfl(state%psetcols,pver)   ! cloudsat reflectivity
   real(r8), target :: acsrfl(state%psetcols,pver)  ! cloudsat average
   real(r8), target :: fcsrfl(state%psetcols,pver)
   real(r8), target :: rercld(state%psetcols,pver)  ! effective radius calculation for rain + cloud
   real(r8), target :: ncai(state%psetcols,pver)    ! output number conc of ice nuclei available (1/m3)
   real(r8), target :: ncal(state%psetcols,pver)    ! output number conc of CCN (1/m3)
   real(r8), target :: qrout2(state%psetcols,pver)
   real(r8), target :: qsout2(state%psetcols,pver)
   real(r8), target :: nrout2(state%psetcols,pver)
   real(r8), target :: nsout2(state%psetcols,pver)
   real(r8), target :: freqs(state%psetcols,pver)
   real(r8), target :: freqr(state%psetcols,pver)
   real(r8), target :: nfice(state%psetcols,pver)
   real(r8), target :: qcrat(state%psetcols,pver)   ! qc limiter ratio (1=no limit)

   ! Object that packs columns with clouds/precip.
   type(MGPacker) :: packer

   ! Packed versions of inputs.
   real(r8), allocatable :: packed_t(:,:)
   real(r8), allocatable :: packed_q(:,:)
   real(r8), allocatable :: packed_qc(:,:)
   real(r8), allocatable :: packed_nc(:,:)
   real(r8), allocatable :: packed_qi(:,:)
   real(r8), allocatable :: packed_ni(:,:)
   real(r8), allocatable :: packed_qr(:,:)
   real(r8), allocatable :: packed_nr(:,:)
   real(r8), allocatable :: packed_qs(:,:)
   real(r8), allocatable :: packed_ns(:,:)

   real(r8), allocatable :: packed_relvar(:,:)
   real(r8), allocatable :: packed_accre_enhan(:,:)

   real(r8), allocatable :: packed_p(:,:)
   real(r8), allocatable :: packed_pdel(:,:)

   ! This is only needed for MG1.5, and can be removed when support for
   ! that version is dropped.
   real(r8), allocatable :: packed_pint(:,:)

   real(r8), allocatable :: packed_cldn(:,:)
   real(r8), allocatable :: packed_liqcldf(:,:)
   real(r8), allocatable :: packed_icecldf(:,:)

   real(r8), allocatable :: packed_naai(:,:)
   real(r8), allocatable :: packed_npccn(:,:)

   real(r8), allocatable :: packed_rndst(:,:,:)
   real(r8), allocatable :: packed_nacon(:,:,:)

   ! Optional outputs.
   real(r8), pointer :: packed_tnd_qsnow(:,:)
   real(r8), pointer :: packed_tnd_nsnow(:,:)
   real(r8), pointer :: packed_re_ice(:,:)

   real(r8), pointer :: packed_frzimm(:,:)
   real(r8), pointer :: packed_frzcnt(:,:)
   real(r8), pointer :: packed_frzdep(:,:)

   ! Output field post-processing.
   type(MGPostProc) :: post_proc

   ! Packed versions of outputs.
   real(r8), allocatable, target :: packed_rate1ord_cw2pr_st(:,:)
   real(r8), allocatable, target :: packed_tlat(:,:)
   real(r8), allocatable, target :: packed_qvlat(:,:)
   real(r8), allocatable, target :: packed_qctend(:,:)
   real(r8), allocatable, target :: packed_qitend(:,:)
   real(r8), allocatable, target :: packed_nctend(:,:)
   real(r8), allocatable, target :: packed_nitend(:,:)

   real(r8), allocatable, target :: packed_qrtend(:,:)
   real(r8), allocatable, target :: packed_qstend(:,:)
   real(r8), allocatable, target :: packed_nrtend(:,:)
   real(r8), allocatable, target :: packed_nstend(:,:)

   real(r8), allocatable, target :: packed_prect(:)
   real(r8), allocatable, target :: packed_preci(:)
   real(r8), allocatable, target :: packed_nevapr(:,:)
   real(r8), allocatable, target :: packed_am_evp_st(:,:)
   real(r8), allocatable, target :: packed_evapsnow(:,:)
   real(r8), allocatable, target :: packed_prain(:,:)
   real(r8), allocatable, target :: packed_prodsnow(:,:)
   real(r8), allocatable, target :: packed_cmeout(:,:)
   real(r8), allocatable, target :: packed_qsout(:,:)
   real(r8), allocatable, target :: packed_rflx(:,:)
   real(r8), allocatable, target :: packed_sflx(:,:)
   real(r8), allocatable, target :: packed_qrout(:,:)
   real(r8), allocatable, target :: packed_qcsevap(:,:)
   real(r8), allocatable, target :: packed_qisevap(:,:)
   real(r8), allocatable, target :: packed_qvres(:,:)
   real(r8), allocatable, target :: packed_cmei(:,:)
   real(r8), allocatable, target :: packed_vtrmc(:,:)
   real(r8), allocatable, target :: packed_vtrmi(:,:)
   real(r8), allocatable, target :: packed_qcsedten(:,:)
   real(r8), allocatable, target :: packed_qisedten(:,:)
   real(r8), allocatable, target :: packed_qrsedten(:,:)
   real(r8), allocatable, target :: packed_qssedten(:,:)
   real(r8), allocatable, target :: packed_umr(:,:)
   real(r8), allocatable, target :: packed_ums(:,:)
   real(r8), allocatable, target :: packed_pra(:,:)
   real(r8), allocatable, target :: packed_prc(:,:)
   real(r8), allocatable, target :: packed_mnuccc(:,:)
   real(r8), allocatable, target :: packed_mnucct(:,:)
   real(r8), allocatable, target :: packed_msacwi(:,:)
   real(r8), allocatable, target :: packed_psacws(:,:)
   real(r8), allocatable, target :: packed_bergs(:,:)
   real(r8), allocatable, target :: packed_berg(:,:)
   real(r8), allocatable, target :: packed_melt(:,:)
   real(r8), allocatable, target :: packed_homo(:,:)
   real(r8), allocatable, target :: packed_qcres(:,:)
   real(r8), allocatable, target :: packed_prci(:,:)
   real(r8), allocatable, target :: packed_prai(:,:)
   real(r8), allocatable, target :: packed_qires(:,:)
   real(r8), allocatable, target :: packed_mnuccr(:,:)
   real(r8), allocatable, target :: packed_pracs(:,:)
   real(r8), allocatable, target :: packed_meltsdt(:,:)
   real(r8), allocatable, target :: packed_frzrdt(:,:)
   real(r8), allocatable, target :: packed_mnuccd(:,:)
   real(r8), allocatable, target :: packed_nrout(:,:)
   real(r8), allocatable, target :: packed_nsout(:,:)
   real(r8), allocatable, target :: packed_refl(:,:)
   real(r8), allocatable, target :: packed_arefl(:,:)
   real(r8), allocatable, target :: packed_areflz(:,:)
   real(r8), allocatable, target :: packed_frefl(:,:)
   real(r8), allocatable, target :: packed_csrfl(:,:)
   real(r8), allocatable, target :: packed_acsrfl(:,:)
   real(r8), allocatable, target :: packed_fcsrfl(:,:)
   real(r8), allocatable, target :: packed_rercld(:,:)
   real(r8), allocatable, target :: packed_ncai(:,:)
   real(r8), allocatable, target :: packed_ncal(:,:)
   real(r8), allocatable, target :: packed_qrout2(:,:)
   real(r8), allocatable, target :: packed_qsout2(:,:)
   real(r8), allocatable, target :: packed_nrout2(:,:)
   real(r8), allocatable, target :: packed_nsout2(:,:)
   real(r8), allocatable, target :: packed_freqs(:,:)
   real(r8), allocatable, target :: packed_freqr(:,:)
   real(r8), allocatable, target :: packed_nfice(:,:)
   real(r8), allocatable, target :: packed_prer_evap(:,:)
   real(r8), allocatable, target :: packed_qcrat(:,:)

   real(r8), allocatable, target :: packed_rel(:,:)
   real(r8), allocatable, target :: packed_rei(:,:)
   real(r8), allocatable, target :: packed_lambdac(:,:)
   real(r8), allocatable, target :: packed_mu(:,:)
   real(r8), allocatable, target :: packed_des(:,:)
   real(r8), allocatable, target :: packed_dei(:,:)

   ! Dummy arrays for cases where we throw away the MG version and
   ! recalculate sizes on the CAM grid to avoid time/subcolumn averaging
   ! issues.
   real(r8), allocatable :: rel_fn_dum(:,:)
   real(r8), allocatable :: dsout2_dum(:,:)
   real(r8), allocatable :: drout_dum(:,:)
   real(r8), allocatable :: reff_rain_dum(:,:)
   real(r8), allocatable :: reff_snow_dum(:,:)

   ! Heterogeneous-only version of mnuccdo.
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

  ! variables for heterogeneous freezing
  real(r8), pointer :: frzimm(:,:)
  real(r8), pointer :: frzcnt(:,:)
  real(r8), pointer :: frzdep(:,:)

   real(r8), pointer :: qme(:,:)

   ! A local copy of state is used for diagnostic calculations
   type(physics_state) :: state_loc
   type(physics_ptend) :: ptend_loc

   real(r8) :: icecldf(state%psetcols,pver) ! Ice cloud fraction
   real(r8) :: liqcldf(state%psetcols,pver) ! Liquid cloud fraction (combined into cloud)

   real(r8), pointer :: rel(:,:)          ! Liquid effective drop radius (microns)
   real(r8), pointer :: rei(:,:)          ! Ice effective drop size (microns)

   real(r8), pointer :: cmeliq(:,:)

   real(r8), pointer :: cld(:,:)          ! Total cloud fraction
   real(r8), pointer :: concld(:,:)       ! Convective cloud fraction
   real(r8), pointer :: iciwpst(:,:)      ! Stratiform in-cloud ice water path for radiation
   real(r8), pointer :: iclwpst(:,:)      ! Stratiform in-cloud liquid water path for radiation
   real(r8), pointer :: cldfsnow(:,:)     ! Cloud fraction for liquid+snow
   real(r8), pointer :: icswp(:,:)        ! In-cloud snow water path

   real(r8) :: icimrst(state%psetcols,pver) ! In stratus ice mixing ratio
   real(r8) :: icwmrst(state%psetcols,pver) ! In stratus water mixing ratio
   real(r8) :: icinc(state%psetcols,pver)   ! In cloud ice number conc
   real(r8) :: icwnc(state%psetcols,pver)   ! In cloud water number conc

   real(r8) :: iclwpi(state%psetcols)       ! Vertically-integrated in-cloud Liquid WP before microphysics
   real(r8) :: iciwpi(state%psetcols)       ! Vertically-integrated in-cloud Ice WP before microphysics

   ! Averaging arrays for effective radius and number....
   real(r8) :: efiout_grid(pcols,pver)
   real(r8) :: efcout_grid(pcols,pver)
   real(r8) :: ncout_grid(pcols,pver)
   real(r8) :: niout_grid(pcols,pver)
   real(r8) :: freqi_grid(pcols,pver)
   real(r8) :: freql_grid(pcols,pver)

   real(r8) :: cdnumc_grid(pcols)           ! Vertically-integrated droplet concentration
   real(r8) :: icimrst_grid_out(pcols,pver) ! In stratus ice mixing ratio
   real(r8) :: icwmrst_grid_out(pcols,pver) ! In stratus water mixing ratio

   ! Cloud fraction used for precipitation.
   real(r8) :: cldmax_grid(pcols,pver)

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

   real(r8) :: icimrst_grid(pcols,pver) ! stratus ice mixing ratio - on grid
   real(r8) :: icwmrst_grid(pcols,pver) ! stratus water mixing ratio - on grid

   real(r8), pointer :: lambdac_grid(:,:)
   real(r8), pointer :: mu_grid(:,:)
   real(r8), pointer :: rel_grid(:,:)
   real(r8), pointer :: rei_grid(:,:)
   real(r8), pointer :: dei_grid(:,:)
   real(r8), pointer :: des_grid(:,:)
   real(r8), pointer :: iclwpst_grid(:,:)

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

   real(r8) :: nc_grid(pcols,pver)
   real(r8) :: ni_grid(pcols,pver)
   real(r8) :: qr_grid(pcols,pver)
   real(r8) :: nr_grid(pcols,pver)
   real(r8) :: qs_grid(pcols,pver)
   real(r8) :: ns_grid(pcols,pver)

   real(r8), pointer :: cmeliq_grid(:,:)

   real(r8), pointer :: prec_str_grid(:)
   real(r8), pointer :: snow_str_grid(:)
   real(r8), pointer :: prec_pcw_grid(:)
   real(r8), pointer :: snow_pcw_grid(:)
   real(r8), pointer :: prec_sed_grid(:)
   real(r8), pointer :: snow_sed_grid(:)
   real(r8), pointer :: cldo_grid(:,:)
   real(r8), pointer :: nevapr_grid(:,:)
   real(r8), pointer :: prain_grid(:,:)
   real(r8), pointer :: mgflxprc_grid(:,:)
   real(r8), pointer :: mgflxsnw_grid(:,:)
   real(r8), pointer :: mgmrprc_grid(:,:)
   real(r8), pointer :: mgmrsnw_grid(:,:)
   real(r8), pointer :: cvreffliq_grid(:,:)
   real(r8), pointer :: cvreffice_grid(:,:)
   real(r8), pointer :: rate1ord_cw2pr_st_grid(:,:)
   real(r8), pointer :: wsedl_grid(:,:)
   real(r8), pointer :: CC_t_grid(:,:)
   real(r8), pointer :: CC_qv_grid(:,:)
   real(r8), pointer :: CC_ql_grid(:,:)
   real(r8), pointer :: CC_qi_grid(:,:)
   real(r8), pointer :: CC_nl_grid(:,:)
   real(r8), pointer :: CC_ni_grid(:,:)
   real(r8), pointer :: CC_qlst_grid(:,:)
   real(r8), pointer :: qme_grid(:,:)
   real(r8), pointer :: iciwpst_grid(:,:)
   real(r8), pointer :: icswp_grid(:,:)
   real(r8), pointer :: ast_grid(:,:)
   real(r8), pointer :: cldfsnow_grid(:,:)

   real(r8), pointer :: qrout_grid_ptr(:,:)
   real(r8), pointer :: qsout_grid_ptr(:,:)
   real(r8), pointer :: nrout_grid_ptr(:,:)
   real(r8), pointer :: nsout_grid_ptr(:,:)

   integer :: nlev   ! number of levels where cloud physics is done
   integer :: mgncol ! size of mgcols
   integer, allocatable :: mgcols(:) ! Columns with microphysics performed

   logical :: use_subcol_microp
   integer :: col_type ! Flag to store whether accessing grid or sub-columns in pbuf_get_field

   character(128) :: errstring   ! return status (non-blank for error return)

   ! For rrtmg optics. specified distribution.
   real(r8), parameter :: dcon   = 25.e-6_r8         ! Convective size distribution effective radius (meters)
   real(r8), parameter :: mucon  = 5.3_r8            ! Convective size distribution shape parameter
   real(r8), parameter :: deicon = 50._r8            ! Convective ice effective diameter (meters)

   real(r8), pointer :: pckdptr(:,:)

   !-------------------------------------------------------------------------------

   call t_startf('micro_mg_cam_tend_init')

   ! Find the number of levels used in the microphysics.
   nlev  = pver - top_lev + 1

   lchnk = state%lchnk
   ncol  = state%ncol
   psetcols = state%psetcols
   ngrdcol  = state%ngrdcol

   itim_old = pbuf_old_tim_idx()

   call phys_getopts(use_subcol_microp_out=use_subcol_microp)

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
   call pbuf_get_field(pbuf, relvar_idx,      relvar,      col_type=col_type, copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, accre_enhan_idx, accre_enhan, col_type=col_type, copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, cmeliq_idx,      cmeliq,      col_type=col_type, copy_if_needed=use_subcol_microp)

   call pbuf_get_field(pbuf, cld_idx,         cld,     start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), &
        col_type=col_type, copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, concld_idx,      concld,  start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), &
        col_type=col_type, copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, ast_idx,         ast,     start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), &
        col_type=col_type, copy_if_needed=use_subcol_microp)

   if (.not. do_cldice) then
      call pbuf_get_field(pbuf, tnd_qsnow_idx,   tnd_qsnow,   col_type=col_type, copy_if_needed=use_subcol_microp)
      call pbuf_get_field(pbuf, tnd_nsnow_idx,   tnd_nsnow,   col_type=col_type, copy_if_needed=use_subcol_microp)
      call pbuf_get_field(pbuf, re_ice_idx,      re_ice,      col_type=col_type, copy_if_needed=use_subcol_microp)
   end if

   if (use_hetfrz_classnuc) then
      call pbuf_get_field(pbuf, frzimm_idx, frzimm, col_type=col_type, copy_if_needed=use_subcol_microp)
      call pbuf_get_field(pbuf, frzcnt_idx, frzcnt, col_type=col_type, copy_if_needed=use_subcol_microp)
      call pbuf_get_field(pbuf, frzdep_idx, frzdep, col_type=col_type, copy_if_needed=use_subcol_microp)
   end if

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
   call pbuf_get_field(pbuf, prer_evap_idx,   prer_evap,   col_type=col_type)
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
   call pbuf_get_field(pbuf, cc_t_idx,        CC_t,     start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cc_qv_idx,       CC_qv,    start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cc_ql_idx,       CC_ql,    start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cc_qi_idx,       CC_qi,    start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cc_nl_idx,       CC_nl,    start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cc_ni_idx,       CC_ni,    start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)
   call pbuf_get_field(pbuf, cc_qlst_idx,     CC_qlst,  start=(/1,1,itim_old/), kount=(/psetcols,pver,1/), col_type=col_type)

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

   call pbuf_get_field(pbuf, evprain_st_idx,  evprain_st_grid)
   call pbuf_get_field(pbuf, evpsnow_st_idx,  evpsnow_st_grid)

   ! Only MG 1 defines this field so far.
   if (micro_mg_version == 1 .and. micro_mg_sub_version == 0) then
      call pbuf_get_field(pbuf, am_evp_st_idx,   am_evp_st_grid)
   end if
   
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
   if (micro_mg_version > 1) then
      lq(ixrain) = .true.
      lq(ixsnow) = .true.
      lq(ixnumrain) = .true.
      lq(ixnumsnow) = .true.
   end if

   ! the name 'cldwat' triggers special tests on cldliq
   ! and cldice in physics_update
   call physics_ptend_init(ptend, psetcols, "cldwat", ls=.true., lq=lq)

   select case (micro_mg_version)
   case (1)
      select case (micro_mg_sub_version)
      case (0)
         call micro_mg_get_cols1_0(ncol, nlev, top_lev, state%q(:,:,ixcldliq), &
              state%q(:,:,ixcldice), mgncol, mgcols)
      case (5)
         call micro_mg_get_cols1_5(ncol, nlev, top_lev, state%q(:,:,ixcldliq), &
              state%q(:,:,ixcldice), mgncol, mgcols)
      end select
   case (2)
      call micro_mg_get_cols2_0(ncol, nlev, top_lev, state%q(:,:,ixcldliq), &
           state%q(:,:,ixcldice), state%q(:,:,ixrain), state%q(:,:,ixsnow), &
           mgncol, mgcols)
   end select

   packer = MGPacker(psetcols, pver, mgcols, top_lev)
   post_proc = MGPostProc(packer)

   allocate(packed_rate1ord_cw2pr_st(mgncol,nlev))
   pckdptr => packed_rate1ord_cw2pr_st ! workaround an apparent pgi compiler bug on goldbach
   call post_proc%add_field(p(rate1cld), pckdptr)
   allocate(packed_tlat(mgncol,nlev))
   call post_proc%add_field(p(tlat), p(packed_tlat))
   allocate(packed_qvlat(mgncol,nlev))
   call post_proc%add_field(p(qvlat), p(packed_qvlat))
   allocate(packed_qctend(mgncol,nlev))
   call post_proc%add_field(p(qcten), p(packed_qctend))
   allocate(packed_qitend(mgncol,nlev))
   call post_proc%add_field(p(qiten), p(packed_qitend))
   allocate(packed_nctend(mgncol,nlev))
   call post_proc%add_field(p(ncten), p(packed_nctend))
   allocate(packed_nitend(mgncol,nlev))
   call post_proc%add_field(p(niten), p(packed_nitend))

   if (micro_mg_version > 1) then
      allocate(packed_qrtend(mgncol,nlev))
      call post_proc%add_field(p(qrten), p(packed_qrtend))
      allocate(packed_qstend(mgncol,nlev))
      call post_proc%add_field(p(qsten), p(packed_qstend))
      allocate(packed_nrtend(mgncol,nlev))
      call post_proc%add_field(p(nrten), p(packed_nrtend))
      allocate(packed_nstend(mgncol,nlev))
      call post_proc%add_field(p(nsten), p(packed_nstend))
      allocate(packed_umr(mgncol,nlev))
      call post_proc%add_field(p(umr), p(packed_umr))
      allocate(packed_ums(mgncol,nlev))
      call post_proc%add_field(p(ums), p(packed_ums))
   else if (micro_mg_sub_version == 0) then
      allocate(packed_am_evp_st(mgncol,nlev))
      call post_proc%add_field(p(am_evp_st), p(packed_am_evp_st))
   end if

   allocate(packed_prect(mgncol))
   call post_proc%add_field(p(prect), p(packed_prect))
   allocate(packed_preci(mgncol))
   call post_proc%add_field(p(preci), p(packed_preci))
   allocate(packed_nevapr(mgncol,nlev))
   call post_proc%add_field(p(nevapr), p(packed_nevapr))
   allocate(packed_evapsnow(mgncol,nlev))
   call post_proc%add_field(p(evapsnow), p(packed_evapsnow))
   allocate(packed_prain(mgncol,nlev))
   call post_proc%add_field(p(prain), p(packed_prain))
   allocate(packed_prodsnow(mgncol,nlev))
   call post_proc%add_field(p(prodsnow), p(packed_prodsnow))
   allocate(packed_cmeout(mgncol,nlev))
   call post_proc%add_field(p(cmeice), p(packed_cmeout))
   allocate(packed_qsout(mgncol,nlev))
   call post_proc%add_field(p(qsout), p(packed_qsout))
   allocate(packed_rflx(mgncol,nlev+1))
   call post_proc%add_field(p(rflx), p(packed_rflx))
   allocate(packed_sflx(mgncol,nlev+1))
   call post_proc%add_field(p(sflx), p(packed_sflx))
   allocate(packed_qrout(mgncol,nlev))
   call post_proc%add_field(p(qrout), p(packed_qrout))
   allocate(packed_qcsevap(mgncol,nlev))
   call post_proc%add_field(p(qcsevap), p(packed_qcsevap))
   allocate(packed_qisevap(mgncol,nlev))
   call post_proc%add_field(p(qisevap), p(packed_qisevap))
   allocate(packed_qvres(mgncol,nlev))
   call post_proc%add_field(p(qvres), p(packed_qvres))
   allocate(packed_cmei(mgncol,nlev))
   call post_proc%add_field(p(cmeiout), p(packed_cmei))
   allocate(packed_vtrmc(mgncol,nlev))
   call post_proc%add_field(p(vtrmc), p(packed_vtrmc))
   allocate(packed_vtrmi(mgncol,nlev))
   call post_proc%add_field(p(vtrmi), p(packed_vtrmi))
   allocate(packed_qcsedten(mgncol,nlev))
   call post_proc%add_field(p(qcsedten), p(packed_qcsedten))
   allocate(packed_qisedten(mgncol,nlev))
   call post_proc%add_field(p(qisedten), p(packed_qisedten))
   if (micro_mg_version > 1) then
      allocate(packed_qrsedten(mgncol,nlev))
      call post_proc%add_field(p(qrsedten), p(packed_qrsedten))
      allocate(packed_qssedten(mgncol,nlev))
      call post_proc%add_field(p(qssedten), p(packed_qssedten))
   end if

   allocate(packed_pra(mgncol,nlev))
   call post_proc%add_field(p(prao), p(packed_pra))
   allocate(packed_prc(mgncol,nlev))
   call post_proc%add_field(p(prco), p(packed_prc))
   allocate(packed_mnuccc(mgncol,nlev))
   call post_proc%add_field(p(mnuccco), p(packed_mnuccc))
   allocate(packed_mnucct(mgncol,nlev))
   call post_proc%add_field(p(mnuccto), p(packed_mnucct))
   allocate(packed_msacwi(mgncol,nlev))
   call post_proc%add_field(p(msacwio), p(packed_msacwi))
   allocate(packed_psacws(mgncol,nlev))
   call post_proc%add_field(p(psacwso), p(packed_psacws))
   allocate(packed_bergs(mgncol,nlev))
   call post_proc%add_field(p(bergso), p(packed_bergs))
   allocate(packed_berg(mgncol,nlev))
   call post_proc%add_field(p(bergo), p(packed_berg))
   allocate(packed_melt(mgncol,nlev))
   call post_proc%add_field(p(melto), p(packed_melt))
   allocate(packed_homo(mgncol,nlev))
   call post_proc%add_field(p(homoo), p(packed_homo))
   allocate(packed_qcres(mgncol,nlev))
   call post_proc%add_field(p(qcreso), p(packed_qcres))
   allocate(packed_prci(mgncol,nlev))
   call post_proc%add_field(p(prcio), p(packed_prci))
   allocate(packed_prai(mgncol,nlev))
   call post_proc%add_field(p(praio), p(packed_prai))
   allocate(packed_qires(mgncol,nlev))
   call post_proc%add_field(p(qireso), p(packed_qires))
   allocate(packed_mnuccr(mgncol,nlev))
   call post_proc%add_field(p(mnuccro), p(packed_mnuccr))
   allocate(packed_pracs(mgncol,nlev))
   call post_proc%add_field(p(pracso), p(packed_pracs))
   allocate(packed_meltsdt(mgncol,nlev))
   call post_proc%add_field(p(meltsdt), p(packed_meltsdt))
   allocate(packed_frzrdt(mgncol,nlev))
   call post_proc%add_field(p(frzrdt), p(packed_frzrdt))
   allocate(packed_mnuccd(mgncol,nlev))
   call post_proc%add_field(p(mnuccdo), p(packed_mnuccd))
   allocate(packed_nrout(mgncol,nlev))
   call post_proc%add_field(p(nrout), p(packed_nrout))
   allocate(packed_nsout(mgncol,nlev))
   call post_proc%add_field(p(nsout), p(packed_nsout))

   allocate(packed_refl(mgncol,nlev))
   call post_proc%add_field(p(refl), p(packed_refl), fillvalue=-9999._r8)
   allocate(packed_arefl(mgncol,nlev))
   call post_proc%add_field(p(arefl), p(packed_arefl))
   allocate(packed_areflz(mgncol,nlev))
   call post_proc%add_field(p(areflz), p(packed_areflz))
   allocate(packed_frefl(mgncol,nlev))
   call post_proc%add_field(p(frefl), p(packed_frefl))
   allocate(packed_csrfl(mgncol,nlev))
   call post_proc%add_field(p(csrfl), p(packed_csrfl), fillvalue=-9999._r8)
   allocate(packed_acsrfl(mgncol,nlev))
   call post_proc%add_field(p(acsrfl), p(packed_acsrfl))
   allocate(packed_fcsrfl(mgncol,nlev))
   call post_proc%add_field(p(fcsrfl), p(packed_fcsrfl))

   allocate(packed_rercld(mgncol,nlev))
   call post_proc%add_field(p(rercld), p(packed_rercld))
   allocate(packed_ncai(mgncol,nlev))
   call post_proc%add_field(p(ncai), p(packed_ncai))
   allocate(packed_ncal(mgncol,nlev))
   call post_proc%add_field(p(ncal), p(packed_ncal))
   allocate(packed_qrout2(mgncol,nlev))
   call post_proc%add_field(p(qrout2), p(packed_qrout2))
   allocate(packed_qsout2(mgncol,nlev))
   call post_proc%add_field(p(qsout2), p(packed_qsout2))
   allocate(packed_nrout2(mgncol,nlev))
   call post_proc%add_field(p(nrout2), p(packed_nrout2))
   allocate(packed_nsout2(mgncol,nlev))
   call post_proc%add_field(p(nsout2), p(packed_nsout2))
   allocate(packed_freqs(mgncol,nlev))
   call post_proc%add_field(p(freqs), p(packed_freqs))
   allocate(packed_freqr(mgncol,nlev))
   call post_proc%add_field(p(freqr), p(packed_freqr))
   allocate(packed_nfice(mgncol,nlev))
   call post_proc%add_field(p(nfice), p(packed_nfice))
   if (micro_mg_version /= 1 .or. micro_mg_sub_version /= 0) then
      allocate(packed_qcrat(mgncol,nlev))
      call post_proc%add_field(p(qcrat), p(packed_qcrat), fillvalue=1._r8)
   end if

   ! The following are all variables related to sizes, where it does not
   ! necessarily make sense to average over time steps. Instead, we keep
   ! the value from the last substep, which is what "accum_null" does.
   allocate(packed_rel(mgncol,nlev))
   call post_proc%add_field(p(rel), p(packed_rel), &
        fillvalue=10._r8, accum_method=accum_null)
   allocate(packed_rei(mgncol,nlev))
   call post_proc%add_field(p(rei), p(packed_rei), &
        fillvalue=25._r8, accum_method=accum_null)
   allocate(packed_lambdac(mgncol,nlev))
   call post_proc%add_field(p(lambdac), p(packed_lambdac), &
        accum_method=accum_null)
   allocate(packed_mu(mgncol,nlev))
   call post_proc%add_field(p(mu), p(packed_mu), &
        accum_method=accum_null)
   allocate(packed_des(mgncol,nlev))
   call post_proc%add_field(p(des), p(packed_des), &
        accum_method=accum_null)
   allocate(packed_dei(mgncol,nlev))
   call post_proc%add_field(p(dei), p(packed_dei), &
        accum_method=accum_null)
   allocate(packed_prer_evap(mgncol,nlev))
   call post_proc%add_field(p(prer_evap), p(packed_prer_evap), &
        accum_method=accum_null)

   ! Allocate all the dummies with MG sizes.
   allocate(rel_fn_dum(mgncol,nlev))
   allocate(dsout2_dum(mgncol,nlev))
   allocate(drout_dum(mgncol,nlev))
   allocate(reff_rain_dum(mgncol,nlev))
   allocate(reff_snow_dum(mgncol,nlev))

   ! Pack input variables that are not updated during substeps.
   allocate(packed_relvar(mgncol,nlev))
   packed_relvar = packer%pack(relvar)
   allocate(packed_accre_enhan(mgncol,nlev))
   packed_accre_enhan = packer%pack(accre_enhan)

   allocate(packed_p(mgncol,nlev))
   packed_p = packer%pack(state_loc%pmid)
   allocate(packed_pdel(mgncol,nlev))
   packed_pdel = packer%pack(state_loc%pdel)

   allocate(packed_pint(mgncol,nlev+1))
   packed_pint = packer%pack_interface(state_loc%pint)

   allocate(packed_cldn(mgncol,nlev))
   packed_cldn = packer%pack(ast)
   allocate(packed_liqcldf(mgncol,nlev))
   packed_liqcldf = packer%pack(alst_mic)
   allocate(packed_icecldf(mgncol,nlev))
   packed_icecldf = packer%pack(aist_mic)

   allocate(packed_naai(mgncol,nlev))
   packed_naai = packer%pack(naai)
   allocate(packed_npccn(mgncol,nlev))
   packed_npccn = packer%pack(npccn)

   allocate(packed_rndst(mgncol,nlev,size(rndst, 3)))
   packed_rndst = packer%pack(rndst)
   allocate(packed_nacon(mgncol,nlev,size(nacon, 3)))
   packed_nacon = packer%pack(nacon)

   if (.not. do_cldice) then
      allocate(packed_tnd_qsnow(mgncol,nlev))
      packed_tnd_qsnow = packer%pack(tnd_qsnow)
      allocate(packed_tnd_nsnow(mgncol,nlev))
      packed_tnd_nsnow = packer%pack(tnd_nsnow)
      allocate(packed_re_ice(mgncol,nlev))
      packed_re_ice = packer%pack(re_ice)
   else
      nullify(packed_tnd_qsnow)
      nullify(packed_tnd_nsnow)
      nullify(packed_re_ice)
   end if

   if (use_hetfrz_classnuc) then
      allocate(packed_frzimm(mgncol,nlev))
      packed_frzimm = packer%pack(frzimm)
      allocate(packed_frzcnt(mgncol,nlev))
      packed_frzcnt = packer%pack(frzcnt)
      allocate(packed_frzdep(mgncol,nlev))
      packed_frzdep = packer%pack(frzdep)
   else
      nullify(packed_frzimm)
      nullify(packed_frzcnt)
      nullify(packed_frzdep)
   end if

   ! Allocate input variables that are updated during substeps.
   allocate(packed_t(mgncol,nlev))
   allocate(packed_q(mgncol,nlev))
   allocate(packed_qc(mgncol,nlev))
   allocate(packed_nc(mgncol,nlev))
   allocate(packed_qi(mgncol,nlev))
   allocate(packed_ni(mgncol,nlev))
   if (micro_mg_version > 1) then
      allocate(packed_qr(mgncol,nlev))
      allocate(packed_nr(mgncol,nlev))
      allocate(packed_qs(mgncol,nlev))
      allocate(packed_ns(mgncol,nlev))
   end if
   call t_stopf('micro_mg_cam_tend_init')

   call t_startf('micro_mg_cam_tend_loop')
   do it = 1, num_steps

      ! Pack input variables that are updated during substeps.
      packed_t = packer%pack(state_loc%t)
      packed_q = packer%pack(state_loc%q(:,:,1))
      packed_qc = packer%pack(state_loc%q(:,:,ixcldliq))
      packed_nc = packer%pack(state_loc%q(:,:,ixnumliq))
      packed_qi = packer%pack(state_loc%q(:,:,ixcldice))
      packed_ni = packer%pack(state_loc%q(:,:,ixnumice))
      if (micro_mg_version > 1) then
         packed_qr = packer%pack(state_loc%q(:,:,ixrain))
         packed_nr = packer%pack(state_loc%q(:,:,ixnumrain))
         packed_qs = packer%pack(state_loc%q(:,:,ixsnow))
         packed_ns = packer%pack(state_loc%q(:,:,ixnumsnow))
      end if

      select case (micro_mg_version)
      case (1)
         select case (micro_mg_sub_version)
         case (0)

            call t_startf('micro_mg_tend1')
            call micro_mg_tend1_0( &
                 microp_uniform, mgncol, nlev, mgncol, 1, dtime/num_steps, &
                 packed_t, packed_q, packed_qc, packed_qi, packed_nc,     &
                 packed_ni, packed_p, packed_pdel, packed_cldn, packed_liqcldf,&
                 packed_relvar, packed_accre_enhan,                             &
                 packed_icecldf, packed_rate1ord_cw2pr_st, packed_naai, packed_npccn,                 &
                 packed_rndst, packed_nacon, packed_tlat, packed_qvlat, packed_qctend,                &
                 packed_qitend, packed_nctend, packed_nitend, packed_rel, rel_fn_dum,      &
                 packed_rei, packed_prect, packed_preci, packed_nevapr, packed_evapsnow, packed_am_evp_st, &
                 packed_prain, packed_prodsnow, packed_cmeout, packed_dei, packed_mu,                &
                 packed_lambdac, packed_qsout, packed_des, packed_rflx, packed_sflx,                 &
                 packed_qrout, reff_rain_dum, reff_snow_dum, packed_qcsevap, packed_qisevap,   &
                 packed_qvres, packed_cmei, packed_vtrmc, packed_vtrmi, packed_qcsedten,          &
                 packed_qisedten, packed_pra, packed_prc, packed_mnuccc, packed_mnucct,          &
                 packed_msacwi, packed_psacws, packed_bergs, packed_berg, packed_melt,          &
                 packed_homo, packed_qcres, packed_prci, packed_prai, packed_qires,             &
                 packed_mnuccr, packed_pracs, packed_meltsdt, packed_frzrdt, packed_mnuccd,       &
                 packed_nrout, packed_nsout, packed_refl, packed_arefl, packed_areflz,               &
                 packed_frefl, packed_csrfl, packed_acsrfl, packed_fcsrfl, packed_rercld,            &
                 packed_ncai, packed_ncal, packed_qrout2, packed_qsout2, packed_nrout2,              &
                 packed_nsout2, drout_dum, dsout2_dum, packed_freqs,packed_freqr,            &
                 packed_nfice, packed_prer_evap, do_cldice, errstring,                      &
		 packed_tnd_qsnow, packed_tnd_nsnow, packed_re_ice,             &
                 packed_frzimm, packed_frzcnt, packed_frzdep)
            call t_stopf('micro_mg_tend1')

         case (5)

            call t_startf('micro_mg_tend1_5')
            call micro_mg_tend1_5( &
                 mgncol,   nlev,     dtime/num_steps,    &
                 packed_t,       packed_q,                     &
                 packed_qc,      packed_qi,    &
                 packed_nc,      packed_ni,    &
                 packed_relvar,             packed_accre_enhan,                            &
                 packed_p,     packed_pdel,     packed_pint,     &
                 packed_cldn,                packed_liqcldf,           packed_icecldf,           &
                 packed_rate1ord_cw2pr_st,           packed_naai,     packed_npccn,    packed_rndst,    packed_nacon,    &
                 packed_tlat,     packed_qvlat,    packed_qctend,    packed_qitend,    packed_nctend,    packed_nitend,    &
                 packed_rel,      rel_fn_dum,   packed_rei,                packed_prect,    packed_preci,    &
                 packed_nevapr,   packed_evapsnow, packed_prain,    packed_prodsnow, packed_cmeout,   packed_dei,      &
                 packed_mu,       packed_lambdac,  packed_qsout,    packed_des,      packed_rflx,     packed_sflx,     &
                 packed_qrout,              reff_rain_dum,          reff_snow_dum,          &
                 packed_qcsevap,  packed_qisevap,  packed_qvres,    packed_cmei,  packed_vtrmc,   packed_vtrmi,    &
                 packed_qcsedten, packed_qisedten, packed_pra,     packed_prc,     packed_mnuccc,  packed_mnucct,  &
                 packed_msacwi,  packed_psacws,  packed_bergs,   packed_berg,    packed_melt,    packed_homo,    &
                 packed_qcres,             packed_prci,    packed_prai,    packed_qires,             &
                 packed_mnuccr,  packed_pracs,   packed_meltsdt,  packed_frzrdt,   packed_mnuccd,            &
                 packed_nrout,   packed_nsout,    packed_refl,     packed_arefl,    packed_areflz,   packed_frefl,    &
                 packed_csrfl,    packed_acsrfl,   packed_fcsrfl,             packed_rercld,             &
                 packed_ncai,     packed_ncal,     packed_qrout2,   packed_qsout2,   packed_nrout2,   packed_nsout2,   &
                 drout_dum,   dsout2_dum,   packed_freqs,    packed_freqr,    packed_nfice,    packed_qcrat,    &
                 errstring, &
                 packed_tnd_qsnow,          packed_tnd_nsnow,          packed_re_ice, packed_prer_evap,             &
                 packed_frzimm, packed_frzcnt, packed_frzdep)
            call t_stopf('micro_mg_tend1_5')

         end select
      case(2)
         select case (micro_mg_sub_version)
         case (0)

            call t_startf('micro_mg_tend2')
            call micro_mg_tend2_0( &
                 mgncol,         nlev,           dtime/num_steps,&
                 packed_t,               packed_q,               &
                 packed_qc,              packed_qi,              &
                 packed_nc,              packed_ni,              &
                 packed_qr,              packed_qs,              &
                 packed_nr,              packed_ns,              &
                 packed_relvar,          packed_accre_enhan,     &
                 packed_p,               packed_pdel,            &
                 packed_cldn,    packed_liqcldf, packed_icecldf, &
                 packed_rate1ord_cw2pr_st,                       &
                 packed_naai,            packed_npccn,           &
                 packed_rndst,           packed_nacon,           &
                 packed_tlat,            packed_qvlat,           &
                 packed_qctend,          packed_qitend,          &
                 packed_nctend,          packed_nitend,          &
                 packed_qrtend,          packed_qstend,          &
                 packed_nrtend,          packed_nstend,          &
                 packed_rel,     rel_fn_dum,     packed_rei,     &
                 packed_prect,           packed_preci,           &
                 packed_nevapr,          packed_evapsnow,        &
                 packed_prain,           packed_prodsnow,        &
                 packed_cmeout,          packed_dei,             &
                 packed_mu,              packed_lambdac,         &
                 packed_qsout,           packed_des,             &
                 packed_rflx,    packed_sflx,    packed_qrout,   &
                 reff_rain_dum,          reff_snow_dum,          &
                 packed_qcsevap, packed_qisevap, packed_qvres,   &
                 packed_cmei,    packed_vtrmc,   packed_vtrmi,   &
                 packed_umr,             packed_ums,             &
                 packed_qcsedten,        packed_qisedten,        &
                 packed_qrsedten,        packed_qssedten,        &
                 packed_pra,             packed_prc,             &
                 packed_mnuccc,  packed_mnucct,  packed_msacwi,  &
                 packed_psacws,  packed_bergs,   packed_berg,    &
                 packed_melt,            packed_homo,            &
                 packed_qcres,   packed_prci,    packed_prai,    &
                 packed_qires,   packed_mnuccr,  packed_pracs,   &
                 packed_meltsdt, packed_frzrdt,  packed_mnuccd,  &
                 packed_nrout,           packed_nsout,           &
                 packed_refl,    packed_arefl,   packed_areflz,  &
                 packed_frefl,   packed_csrfl,   packed_acsrfl,  &
                 packed_fcsrfl,          packed_rercld,          &
                 packed_ncai,            packed_ncal,            &
                 packed_qrout2,          packed_qsout2,          &
                 packed_nrout2,          packed_nsout2,          &
                 drout_dum,              dsout2_dum,             &
                 packed_freqs,           packed_freqr,           &
                 packed_nfice,           packed_qcrat,           &
                 errstring, &
                 packed_tnd_qsnow,packed_tnd_nsnow,packed_re_ice,&
		 packed_prer_evap,                                     &
                 packed_frzimm,  packed_frzcnt,  packed_frzdep   )
            call t_stopf('micro_mg_tend2')
         end select
      end select

      call handle_errmsg(errstring, subname="micro_mg_tend")

      call physics_ptend_init(ptend_loc, psetcols, "micro_mg", &
                              ls=.true., lq=lq)

      ! Set local tendency.
      ptend_loc%s               = packer%unpack(packed_tlat, 0._r8)
      ptend_loc%q(:,:,1)        = packer%unpack(packed_qvlat, 0._r8)
      ptend_loc%q(:,:,ixcldliq) = packer%unpack(packed_qctend, 0._r8)
      ptend_loc%q(:,:,ixcldice) = packer%unpack(packed_qitend, 0._r8)
      ptend_loc%q(:,:,ixnumliq) = packer%unpack(packed_nctend, &
           -state_loc%q(:,:,ixnumliq)/(dtime/num_steps))
      if (do_cldice) then
         ptend_loc%q(:,:,ixnumice) = packer%unpack(packed_nitend, &
              -state_loc%q(:,:,ixnumice)/(dtime/num_steps))
      else
         ! In this case, the tendency should be all 0.
         if (any(packed_nitend /= 0._r8)) &
              call endrun("micro_mg_cam:ERROR - MG microphysics is configured not to prognose cloud ice,"// &
              " but micro_mg_tend has ice number tendencies.")
         ptend_loc%q(:,:,ixnumice) = 0._r8
      end if

      if (micro_mg_version > 1) then
         ptend_loc%q(:,:,ixrain)    = packer%unpack(packed_qrtend, 0._r8)
         ptend_loc%q(:,:,ixsnow)    = packer%unpack(packed_qstend, 0._r8)
         ptend_loc%q(:,:,ixnumrain) = packer%unpack(packed_nrtend, &
              -state_loc%q(:,:,ixnumrain)/(dtime/num_steps))
         ptend_loc%q(:,:,ixnumsnow) = packer%unpack(packed_nstend, &
              -state_loc%q(:,:,ixnumsnow)/(dtime/num_steps))
      end if

      ! Sum into overall ptend
      call physics_ptend_sum(ptend_loc, ptend, ncol)

      ! Update local state
      call physics_update(state_loc, ptend_loc, dtime/num_steps)

      ! Sum all outputs for averaging.
      call post_proc%accumulate()

   end do
   call t_stopf('micro_mg_cam_tend_loop')

   call t_startf('micro_mg_cam_tend_fini')
   ! Divide ptend by substeps.
   call physics_ptend_scale(ptend, 1._r8/num_steps, ncol)

   ! Use summed outputs to produce averages
   call post_proc%process_and_unpack()

   call post_proc%finalize()

   if (associated(packed_tnd_qsnow)) deallocate(packed_tnd_qsnow)
   if (associated(packed_tnd_nsnow)) deallocate(packed_tnd_nsnow)
   if (associated(packed_re_ice)) deallocate(packed_re_ice)
   if (associated(packed_frzimm)) deallocate(packed_frzimm)
   if (associated(packed_frzcnt)) deallocate(packed_frzcnt)
   if (associated(packed_frzdep)) deallocate(packed_frzdep)

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
   ! Do not subscript by ncol here, because in physpkg we divide the whole
   ! array and need to avoid an FPE due to uninitialized data.
   prec_pcw = prect
   snow_pcw = preci
   prec_sed = 0._r8
   snow_sed = 0._r8
   prec_str = prec_pcw + prec_sed
   snow_str = snow_pcw + snow_sed

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
         !BSINGH- Following code MUST be reworked. This is TEMPORARY solution to maintain BFB
         !PMA: .lt. is replaced by .le. following part of NCAR RRTMG bug fix         
         if(rrtmg_temp_fix ) then
            if( (cldfsnow(i,k) .le. 1.e-4_r8) .and. (qsout(i,k) .gt. 1.e-6_r8) ) then
               cldfsnow(i,k) = 0.25_r8
            endif
         else
            if( (cldfsnow(i,k) .lt. 1.e-4_r8) .and. (qsout(i,k) .gt. 1.e-6_r8) ) then
               cldfsnow(i,k) = 0.25_r8
            endif
         endif
         ! Calculate in-cloud snow water path
         icswp(i,k) = qsout(i,k) / max( mincld, cldfsnow(i,k) ) * state_loc%pdel(i,k) / gravit
      end do
   end do

   ! Calculate cloud fraction for prognostic precip sizes.
   if (micro_mg_version > 1) then
      ! Cloud fraction for purposes of precipitation is maximum cloud
      ! fraction out of all the layers that the precipitation may be
      ! falling down from.
      cldmax = max(mincld, ast)
      do k = top_lev+1, pver
         where (state_loc%q(:ncol,k-1,ixrain) >= qsmall .or. &
              state_loc%q(:ncol,k-1,ixsnow) >= qsmall)
            cldmax(:ncol,k) = max(cldmax(:ncol,k-1), cldmax(:ncol,k))
         end where
      end do
   end if

   ! ------------------------------------------------------ !
   ! ------------------------------------------------------ !
   ! All code from here to the end is on grid columns only  !
   ! ------------------------------------------------------ !
   ! ------------------------------------------------------ !

   ! Average the fields which are needed later in this paramterization to be on the grid
   if (use_subcol_microp) then
      call subcol_field_avg(prec_str,  ngrdcol, lchnk, prec_str_grid)
      call subcol_field_avg(iclwpst,   ngrdcol, lchnk, iclwpst_grid)
      call subcol_field_avg(cvreffliq, ngrdcol, lchnk, cvreffliq_grid)
      call subcol_field_avg(cvreffice, ngrdcol, lchnk, cvreffice_grid)
      call subcol_field_avg(mgflxprc,  ngrdcol, lchnk, mgflxprc_grid)
      call subcol_field_avg(mgflxsnw,  ngrdcol, lchnk, mgflxsnw_grid)
      call subcol_field_avg(qme,       ngrdcol, lchnk, qme_grid)
      call subcol_field_avg(nevapr,    ngrdcol, lchnk, nevapr_grid)
      call subcol_field_avg(prain,     ngrdcol, lchnk, prain_grid)
      call subcol_field_avg(evapsnow,  ngrdcol, lchnk, evpsnow_st_grid)

      if (micro_mg_version == 1 .and. micro_mg_sub_version == 0) then
         call subcol_field_avg(am_evp_st, ngrdcol, lchnk, am_evp_st_grid)
      end if

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
      call subcol_field_avg(prao,      ngrdcol, lchnk, prao_grid)
      call subcol_field_avg(prco,      ngrdcol, lchnk, prco_grid)

      call subcol_field_avg(state_loc%q(:,:,ixnumliq), ngrdcol, lchnk, nc_grid)
      call subcol_field_avg(state_loc%q(:,:,ixnumice), ngrdcol, lchnk, ni_grid)

      if (micro_mg_version > 1) then
         call subcol_field_avg(cldmax,    ngrdcol, lchnk, cldmax_grid)

         call subcol_field_avg(state_loc%q(:,:,ixrain),    ngrdcol, lchnk, qr_grid)
         call subcol_field_avg(state_loc%q(:,:,ixnumrain), ngrdcol, lchnk, nr_grid)
         call subcol_field_avg(state_loc%q(:,:,ixsnow),    ngrdcol, lchnk, qs_grid)
         call subcol_field_avg(state_loc%q(:,:,ixnumsnow), ngrdcol, lchnk, ns_grid)
      end if

   else
      ! These pbuf fields need to be assigned.  There is no corresponding subcol_field_avg
      ! as they are reset before being used, so it would be a needless calculation
      lambdac_grid    => lambdac
      mu_grid         => mu
      rel_grid        => rel
      rei_grid        => rei
      dei_grid        => dei
      des_grid        => des

      ! fields already on grids, so just assign
      prec_str_grid   => prec_str
      iclwpst_grid    => iclwpst
      cvreffliq_grid  => cvreffliq
      cvreffice_grid  => cvreffice
      mgflxprc_grid   => mgflxprc
      mgflxsnw_grid   => mgflxsnw
      qme_grid        => qme
      nevapr_grid     => nevapr
      prain_grid      => prain

      if (micro_mg_version == 1 .and. micro_mg_sub_version == 0) then
         am_evp_st_grid  = am_evp_st
      end if

      evpsnow_st_grid = evapsnow
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
      prao_grid       = prao
      prco_grid       = prco

      nc_grid = state_loc%q(:,:,ixnumliq)
      ni_grid = state_loc%q(:,:,ixnumice)

      if (micro_mg_version > 1) then
         cldmax_grid = cldmax

         qr_grid = state_loc%q(:,:,ixrain)
         nr_grid = state_loc%q(:,:,ixnumrain)
         qs_grid = state_loc%q(:,:,ixsnow)
         ns_grid = state_loc%q(:,:,ixnumsnow)
      end if

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

   ! Calculate rho (on subcolumns if turned on) for size distribution
   ! parameter calculations and average it if needed
   !
   ! State instead of state_loc to preserve answers for MG1 (and in any
   ! case, it is unlikely to make much difference).
   rho(:ncol,top_lev:) = state%pmid(:ncol,top_lev:) / &
        (rair*state%t(:ncol,top_lev:))
   if (use_subcol_microp) then
      call subcol_field_avg(rho, ngrdcol, lchnk, rho_grid)
   else
      rho_grid = rho
   end if

   ! Effective radius for cloud liquid, fixed number.
   mu_grid = 0._r8
   lambdac_grid = 0._r8
   rel_fn_grid = 10._r8

   ncic_grid = 1.e8_r8

   call size_dist_param_liq(mg_liq_props, icwmrst_grid(:ngrdcol,top_lev:), &
        ncic_grid(:ngrdcol,top_lev:), rho_grid(:ngrdcol,top_lev:), &
        mu_grid(:ngrdcol,top_lev:), lambdac_grid(:ngrdcol,top_lev:))

   where (icwmrst_grid(:ngrdcol,top_lev:) > qsmall)
      rel_fn_grid(:ngrdcol,top_lev:) = &
           (mu_grid(:ngrdcol,top_lev:) + 3._r8)/ &
           lambdac_grid(:ngrdcol,top_lev:)/2._r8 * 1.e6_r8
   end where

   ! Effective radius for cloud liquid, and size parameters
   ! mu_grid and lambdac_grid.
   mu_grid = 0._r8
   lambdac_grid = 0._r8
   rel_grid = 10._r8

   ! Calculate ncic on the grid
   ncic_grid(:ngrdcol,top_lev:) = nc_grid(:ngrdcol,top_lev:) / &
        max(mincld,liqcldf_grid(:ngrdcol,top_lev:))

   call size_dist_param_liq(mg_liq_props, icwmrst_grid(:ngrdcol,top_lev:), &
        ncic_grid(:ngrdcol,top_lev:), rho_grid(:ngrdcol,top_lev:), &
        mu_grid(:ngrdcol,top_lev:), lambdac_grid(:ngrdcol,top_lev:))

   where (icwmrst_grid(:ngrdcol,top_lev:) >= qsmall)
      rel_grid(:ngrdcol,top_lev:) = &
           (mu_grid(:ngrdcol,top_lev:) + 3._r8) / &
           lambdac_grid(:ngrdcol,top_lev:)/2._r8 * 1.e6_r8
   elsewhere
      ! Deal with the fact that size_dist_param_liq sets mu_grid to -100
      ! wherever there is no cloud.
      mu_grid(:ngrdcol,top_lev:) = 0._r8
   end where

   ! Rain/Snow effective diameter.
   drout2_grid = 0._r8
   reff_rain_grid = 0._r8
   des_grid = 0._r8
   dsout2_grid = 0._r8
   reff_snow_grid = 0._r8

   if (micro_mg_version > 1) then
      ! Prognostic precipitation

      where (qr_grid(:ngrdcol,top_lev:) >= 1.e-7_r8)
         drout2_grid(:ngrdcol,top_lev:) = avg_diameter( &
              qr_grid(:ngrdcol,top_lev:), &
              nr_grid(:ngrdcol,top_lev:) * rho_grid(:ngrdcol,top_lev:), &
              rho_grid(:ngrdcol,top_lev:), rhow)

         reff_rain_grid(:ngrdcol,top_lev:) = drout2_grid(:ngrdcol,top_lev:) * &
              1.5_r8 * 1.e6_r8
      end where

      where (qs_grid(:ngrdcol,top_lev:) >= 1.e-7_r8)
         dsout2_grid(:ngrdcol,top_lev:) = avg_diameter( &
              qs_grid(:ngrdcol,top_lev:), &
              ns_grid(:ngrdcol,top_lev:) * rho_grid(:ngrdcol,top_lev:), &
              rho_grid(:ngrdcol,top_lev:), rhosn)

         des_grid(:ngrdcol,top_lev:) = dsout2_grid(:ngrdcol,top_lev:) *&
              3._r8 * rhosn/rhows

         reff_snow_grid(:ngrdcol,top_lev:) = dsout2_grid(:ngrdcol,top_lev:) * &
              1.5_r8 * 1.e6_r8
      end where

   else
      ! Diagnostic precipitation

      where (qrout_grid(:ngrdcol,top_lev:) >= 1.e-7_r8)
         drout2_grid(:ngrdcol,top_lev:) = avg_diameter( &
              qrout_grid(:ngrdcol,top_lev:), &
              nrout_grid(:ngrdcol,top_lev:) * rho_grid(:ngrdcol,top_lev:), &
              rho_grid(:ngrdcol,top_lev:), rhow)

         reff_rain_grid(:ngrdcol,top_lev:) = drout2_grid(:ngrdcol,top_lev:) * &
              1.5_r8 * 1.e6_r8
      end where

      where (qsout_grid(:ngrdcol,top_lev:) >= 1.e-7_r8)
         dsout2_grid(:ngrdcol,top_lev:) = avg_diameter( &
              qsout_grid(:ngrdcol,top_lev:), &
              nsout_grid(:ngrdcol,top_lev:) * rho_grid(:ngrdcol,top_lev:), &
              rho_grid(:ngrdcol,top_lev:), rhosn)

         des_grid(:ngrdcol,top_lev:) = dsout2_grid(:ngrdcol,top_lev:) &
              * 3._r8 * rhosn/rhows

         reff_snow_grid(:ngrdcol,top_lev:) = &
              dsout2_grid(:ngrdcol,top_lev:) * 1.5_r8 * 1.e6_r8
      end where

   end if

   ! Effective radius and diameter for cloud ice.
   rei_grid = 25._r8

   niic_grid(:ngrdcol,top_lev:) = ni_grid(:ngrdcol,top_lev:) / &
        max(mincld,icecldf_grid(:ngrdcol,top_lev:))

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
            tgcmeliq_grid(i) = tgcmeliq_grid(i) + cmeliq_grid(i,k) * &
                 (pdel_grid(i,k) / gravit) / rhoh2o
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
!           write (iulog,*) 'PE_grid:ANOMALY  pe_grid, acprecl_grid, acgcme_grid, tpr_grid, acnum_grid ', &
!                           pe_grid(i),acprecl_grid(i), acgcme_grid(i), tpr_grid(i), acnum_grid(i)
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
   cdnumc_grid(:ngrdcol) = sum(nc_grid(:ngrdcol,top_lev:pver) * &
        pdel_grid(:ngrdcol,top_lev:pver)/gravit, dim=2)

   ! Averaging for new output fields
   efcout_grid      = 0._r8
   efiout_grid      = 0._r8
   ncout_grid       = 0._r8
   niout_grid       = 0._r8
   freql_grid       = 0._r8
   freqi_grid       = 0._r8
   icwmrst_grid_out = 0._r8
   icimrst_grid_out = 0._r8

   do k = top_lev, pver
      do i = 1, ngrdcol
         if ( liqcldf_grid(i,k) > 0.01_r8 .and. icwmrst_grid(i,k) > 5.e-5_r8 ) then
            efcout_grid(i,k) = rel_grid(i,k) * liqcldf_grid(i,k)
            ncout_grid(i,k)  = icwnc_grid(i,k) * liqcldf_grid(i,k)
            freql_grid(i,k)  = liqcldf_grid(i,k)
            icwmrst_grid_out(i,k) = icwmrst_grid(i,k)
         end if
         if ( icecldf_grid(i,k) > 0.01_r8 .and. icimrst_grid(i,k) > 1.e-6_r8 ) then
            efiout_grid(i,k) = rei_grid(i,k) * icecldf_grid(i,k)
            niout_grid(i,k)  = icinc_grid(i,k) * icecldf_grid(i,k)
            freqi_grid(i,k)  = icecldf_grid(i,k)
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

   ! Evaporation of stratiform precipitation fields for UNICON
   evprain_st_grid(:ngrdcol,:pver) = nevapr_grid(:ngrdcol,:pver) - evpsnow_st_grid(:ngrdcol,:pver)
   do k = top_lev, pver
      do i = 1, ngrdcol
         evprain_st_grid(i,k) = max(evprain_st_grid(i,k), 0._r8)
         evpsnow_st_grid(i,k) = max(evpsnow_st_grid(i,k), 0._r8)
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
   call outfld('MPDNLIQ',     ncten,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('MPDNICE',     niten,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('EVAPSNOW',    evapsnow,    psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('QCSEVAP',     qcsevap,     psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('QISEVAP',     qisevap,     psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('QVRES',       qvres,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('VTRMC',       vtrmc,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('VTRMI',       vtrmi,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('QCSEDTEN',    qcsedten,    psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('QISEDTEN',    qisedten,    psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   if (micro_mg_version > 1) then
      call outfld('QRSEDTEN',    qrsedten,    psetcols, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('QSSEDTEN',    qssedten,    psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   end if
   call outfld('MNUCCDO',     mnuccdo,     psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('MNUCCDOhet',  mnuccdohet,  psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('MNUCCRO',     mnuccro,     psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('PRACSO',      pracso ,     psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('MELTSDT',     meltsdt,     psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('FRZRDT',      frzrdt ,     psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   call outfld('FICE',        nfice,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)

   if (micro_mg_version > 1) then
      call outfld('UMR',      umr,         psetcols, lchnk, avg_subcol_field=use_subcol_microp)
      call outfld('UMS',      ums,         psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   end if

   if (.not. (micro_mg_version == 1 .and. micro_mg_sub_version == 0)) then
      call outfld('QCRAT',    qcrat,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
   end if

   ! Example subcolumn outfld call
   if (use_subcol_microp) then
      call outfld('FICE_SCOL',   nfice,       psubcols*pcols, lchnk)
      call outfld('MPDLIQ_SCOL', qcten,       psubcols*pcols, lchnk)
      call outfld('MPDICE_SCOL', qiten,       psubcols*pcols, lchnk)
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
   call t_stopf('micro_mg_cam_tend_fini')

end subroutine micro_mg_cam_tend

function p1(tin) result(pout)
  real(r8), target, intent(in) :: tin(:)
  real(r8), pointer :: pout(:)
  pout => tin
end function p1

function p2(tin) result(pout)
  real(r8), target, intent(in) :: tin(:,:)
  real(r8), pointer :: pout(:,:)
  pout => tin
end function p2

end module micro_mg_cam
