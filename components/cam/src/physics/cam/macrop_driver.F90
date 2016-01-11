
module macrop_driver

  !-------------------------------------------------------------------------------------------------------
  ! Purpose:
  !
  ! Provides the CAM interface to the prognostic cloud macrophysics
  !
  ! Author: Andrew Gettelman, Cheryl Craig October 2010
  ! Origin: modified from stratiform.F90 elements 
  !    (Boville 2002, Coleman 2004, Park 2009, Kay 2010)
  !-------------------------------------------------------------------------------------------------------

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use spmd_utils,    only: masterproc
  use ppgrid,        only: pcols, pver, pverp
  use physconst,     only: latice, latvap
  use phys_control,  only: phys_getopts
  use constituents,  only: cnst_get_ind, pcnst
  use physics_buffer,    only: physics_buffer_desc, pbuf_set_field, pbuf_get_field, pbuf_old_tim_idx
  use time_manager,      only: is_first_step
  use cldwat2m_macro,    only: ini_macro
  use perf_mod,          only: t_startf, t_stopf
  use cam_logfile,       only: iulog
  use cam_abortutils,    only: endrun

  implicit none
  private
  save

  public :: macrop_driver_readnl
  public :: macrop_driver_register
  public :: macrop_driver_init
  public :: macrop_driver_tend
  public :: ice_macro_tend

  logical, public :: do_cldice             ! .true., park macrophysics is prognosing cldice
  logical, public :: do_cldliq             ! .true., park macrophysics is prognosing cldliq
  logical, public :: do_detrain            ! .true., park macrophysics is detraining ice into stratiform

  ! ------------------------- !
  ! Private Module Parameters !
  ! ------------------------- !

  ! 'cu_det_st' : If .true. (.false.), detrain cumulus liquid condensate into the pre-existing liquid stratus 
  !               (environment) without (with) macrophysical evaporation. If there is no pre-esisting stratus, 
  !               evaporate cumulus liquid condensate. This option only influences the treatment of cumulus
  !               liquid condensate, not cumulus ice condensate.

  logical, parameter :: cu_det_st  = .false.  

  integer,          parameter :: i_adet     =  1   ! If '0' ('1'), 'a_det' is maximally overlapped with 'ast' (simply added).
                                                   ! BE CAREFUL : The same parameter should be consistently specified within 'cldwat2m_macro.F90'.

  logical :: micro_do_icesupersat

  ! Parameters used for selecting generalized critical RH for liquid and ice stratus
  integer :: rhminl_opt = 0
  integer :: rhmini_opt = 0


  character(len=16) :: shallow_scheme
  logical           :: use_shfrc                       ! Local copy of flag from convect_shallow_use_shfrc

  integer :: &
    ixcldliq,     &! cloud liquid amount index
    ixcldice,     &! cloud ice amount index
    ixnumliq,     &! cloud liquid number index
    ixnumice,     &! cloud ice water index
    qcwat_idx,    &! qcwat index in physics buffer
    lcwat_idx,    &! lcwat index in physics buffer
    iccwat_idx,   &! iccwat index in physics buffer
    nlwat_idx,    &! nlwat index in physics buffer
    niwat_idx,    &! niwat index in physics buffer
    tcwat_idx,    &! tcwat index in physics buffer
    CC_T_idx,     &!
    CC_qv_idx,    &!
    CC_ql_idx,    &!
    CC_qi_idx,    &!
    CC_nl_idx,    &!
    CC_ni_idx,    &!
    CC_qlst_idx,  &!
    cld_idx,      &! cld index in physics buffer
    ast_idx,      &! stratiform cloud fraction index in physics buffer
    aist_idx,     &! ice stratiform cloud fraction index in physics buffer
    alst_idx,     &! liquid stratiform cloud fraction index in physics buffer
    qist_idx,     &! ice stratiform in-cloud IWC 
    qlst_idx,     &! liquid stratiform in-cloud LWC  
    concld_idx,   &! concld index in physics buffer
    fice_idx,     &  
    cmeliq_idx,   &  
    shfrc_idx,    &
    naai_idx 

  logical :: liqcf_fix

  integer :: &
    tke_idx = -1,       &! tke defined at the model interfaces
    qtl_flx_idx = -1,   &! overbar(w'qtl' where qtl = qv + ql) from the PBL scheme
    qti_flx_idx = -1,   &! overbar(w'qti' where qti = qv + qi) from the PBL scheme
    cmfr_det_idx = -1,  &! detrained convective mass flux from UNICON
    qlr_det_idx = -1,   &! detrained convective ql from UNICON  
    qir_det_idx = -1     ! detrained convective qi from UNICON  
    integer :: rqcr_l_idx = -1    ! mass tendency by detrained convective liquid condensate [ kg/kg/s ] >= 0.
    integer :: rqcr_i_idx = -1    ! mass tendency by detrained convective ice    condensate [ kg/kg/s ] >= 0.
    integer :: rncr_l_idx = -1    ! number tendency by detrained convective liquid condensate [ #/kg/s ] >= 0.
    integer :: rncr_i_idx = -1    ! number tendency by detrained convective ice    condensate [ #/kg/s ] >= 0.

  contains

  ! ===============================================================================
  subroutine macrop_driver_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Namelist variables
   logical  :: macro_park_do_cldice  = .true.   ! do_cldice = .true., park macrophysics is prognosing cldice
   logical  :: macro_park_do_cldliq  = .true.   ! do_cldliq = .true., park macrophysics is prognosing cldliq
   logical  :: macro_park_do_detrain = .true.   ! do_detrain = .true., park macrophysics is detraining ice into stratiform

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'macrop_driver_readnl'

   namelist /macro_park_nl/ macro_park_do_cldice, macro_park_do_cldliq, macro_park_do_detrain
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'macro_park_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, macro_park_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)

      ! set local variables

      do_cldice  = macro_park_do_cldice
      do_cldliq  = macro_park_do_cldliq
      do_detrain = macro_park_do_detrain

   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(do_cldice,             1, mpilog, 0, mpicom)
   call mpibcast(do_cldliq,             1, mpilog, 0, mpicom)
   call mpibcast(do_detrain,            1, mpilog, 0, mpicom)
#endif

end subroutine macrop_driver_readnl

  !================================================================================================

  subroutine macrop_driver_register

  !---------------------------------------------------------------------- !
  !                                                                       !
  ! Register the constituents (cloud liquid and cloud ice) and the fields !
  ! in the physics buffer.                                                !
  !                                                                       !
  !---------------------------------------------------------------------- !

   
   use physics_buffer, only : pbuf_add_field, dtype_r8, dyn_time_lvls

  !-----------------------------------------------------------------------

    call phys_getopts(shallow_scheme_out=shallow_scheme)

    call pbuf_add_field('AST',      'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), ast_idx)
    call pbuf_add_field('AIST',     'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), aist_idx)
    call pbuf_add_field('ALST',     'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), alst_idx)
    call pbuf_add_field('QIST',     'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), qist_idx)
    call pbuf_add_field('QLST',     'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), qlst_idx)
    call pbuf_add_field('CLD',      'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cld_idx)
    call pbuf_add_field('CONCLD',   'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), concld_idx)

    call pbuf_add_field('QCWAT',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), qcwat_idx)
    call pbuf_add_field('LCWAT',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), lcwat_idx)
    call pbuf_add_field('ICCWAT',   'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), iccwat_idx)
    call pbuf_add_field('NLWAT',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), nlwat_idx)
    call pbuf_add_field('NIWAT',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), niwat_idx)
    call pbuf_add_field('TCWAT',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), tcwat_idx)

    call pbuf_add_field('FICE',     'physpkg', dtype_r8, (/pcols,pver/), fice_idx)

    call pbuf_add_field('CMELIQ',   'physpkg', dtype_r8, (/pcols,pver/), cmeliq_idx)

! JHYoon
! These are created in unicon_cam.F90. I'm adding this to make CAM5 default run
    if( shallow_scheme .ne. 'UNICON' ) then
       call pbuf_add_field('cmfr_det',      'global', dtype_r8, (/pcols,pver/), cmfr_det_idx)
       call pbuf_add_field('qlr_det',      'global', dtype_r8, (/pcols,pver/), qlr_det_idx)
       call pbuf_add_field('qir_det',      'global', dtype_r8, (/pcols,pver/), qir_det_idx)
       call pbuf_add_field('rqcr_l',      'global', dtype_r8, (/pcols,pver/), rqcr_l_idx)
       call pbuf_add_field('rqcr_i',      'global', dtype_r8, (/pcols,pver/), rqcr_i_idx)
       call pbuf_add_field('rncr_l',      'global', dtype_r8, (/pcols,pver/), rncr_l_idx)
       call pbuf_add_field('rncr_i',      'global', dtype_r8, (/pcols,pver/), rncr_i_idx)
    end if
! JHYoon

  end subroutine macrop_driver_register

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine macrop_driver_init(pbuf2d)

  !-------------------------------------------- !
  !                                             !
  ! Initialize the cloud water parameterization !
  !                                             ! 
  !-------------------------------------------- !
    use physics_buffer, only : pbuf_get_index
    use cam_history,     only: addfld, add_default
    use convect_shallow, only: convect_shallow_use_shfrc
    
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    logical              :: history_aerosol      ! Output the MAM aerosol tendencies
    logical              :: history_budget       ! Output tendencies and state variables for CAM4
                                                 ! temperature, water vapor, cloud ice and cloud
                                                 ! liquid budgets.
    integer              :: history_budget_histfile_num ! output history file number for budget fields
    integer :: istat
    character(len=*), parameter :: subname = 'macrop_driver_init'
    !-----------------------------------------------------------------------

    ! Initialization routine for cloud macrophysics
    if (shallow_scheme .eq. 'UNICON') rhminl_opt = 1
    call ini_macro(rhminl_opt, rhmini_opt)

    call phys_getopts(history_aerosol_out              = history_aerosol      , &
                      history_budget_out               = history_budget       , &
                      history_budget_histfile_num_out  = history_budget_histfile_num, &
                      micro_do_icesupersat_out         = micro_do_icesupersat, &
                      liqcf_fix_out                    = liqcf_fix )

  ! Find out whether shfrc from convect_shallow will be used in cldfrc

    if( convect_shallow_use_shfrc() ) then
        use_shfrc = .true.
        shfrc_idx = pbuf_get_index('shfrc')
    else 
        use_shfrc = .false.
    endif

  ! Sungsu added for use in the macrophysics

    if( shallow_scheme .eq. 'UNICON' ) then
        tke_idx      = pbuf_get_index('tke')
        qtl_flx_idx  = pbuf_get_index('qtl_flx')
        qti_flx_idx  = pbuf_get_index('qti_flx')
        cmfr_det_idx = pbuf_get_index('cmfr_det')
        qlr_det_idx  = pbuf_get_index('qlr_det')
        qir_det_idx  = pbuf_get_index('qir_det')
        rqcr_l_idx   = pbuf_get_index('rqcr_l')
        rqcr_i_idx   = pbuf_get_index('rqcr_i')
        rncr_l_idx   = pbuf_get_index('rncr_l')
        rncr_i_idx   = pbuf_get_index('rncr_i')
    endif

  ! Sungsu added for use in the macrophysics

    call addfld ('DPDLFLIQ', (/ 'lev' /), 'A', 'kg/kg/s', 'Detrained liquid water from deep convection'             )
    call addfld ('DPDLFICE', (/ 'lev' /), 'A', 'kg/kg/s', 'Detrained ice from deep convection'                      )
    call addfld ('SHDLFLIQ', (/ 'lev' /), 'A', 'kg/kg/s', 'Detrained liquid water from shallow convection'          )
    call addfld ('SHDLFICE', (/ 'lev' /), 'A', 'kg/kg/s', 'Detrained ice from shallow convection'                   )
    call addfld ('DPDLFT', (/ 'lev' /), 'A', 'K/s', 'T-tendency due to deep convective detrainment'           )
    call addfld ('SHDLFT', (/ 'lev' /), 'A', 'K/s', 'T-tendency due to shallow convective detrainment'        )

    call addfld ('ZMDLF', (/ 'lev' /), 'A', 'kg/kg/s', 'Detrained liquid water from ZM convection'               )

    call addfld ('MACPDT', (/ 'lev' /), 'A', 'W/kg', 'Heating tendency - Revised  macrophysics'                )
    call addfld ('MACPDQ', (/ 'lev' /), 'A', 'kg/kg/s', 'Q tendency - Revised macrophysics'                       )
    call addfld ('MACPDLIQ', (/ 'lev' /), 'A', 'kg/kg/s', 'CLDLIQ tendency - Revised macrophysics'                  )
    call addfld ('MACPDICE', (/ 'lev' /), 'A', 'kg/kg/s', 'CLDICE tendency - Revised macrophysics'                  )

    call addfld ('CLDVAPADJ', (/ 'lev' /), 'A', 'kg/kg/s', &
         'Q tendency associated with liq/ice adjustment - Revised macrophysics' )
    call addfld ('CLDLIQADJ', (/ 'lev' /), 'A', 'kg/kg/s', 'CLDLIQ adjustment tendency - Revised macrophysics'       )
    call addfld ('CLDICEADJ', (/ 'lev' /), 'A', 'kg/kg/s', 'CLDICE adjustment tendency - Revised macrophysics'       )
    call addfld ('CLDLIQDET', (/ 'lev' /), 'A', 'kg/kg/s', &
         'Detrainment of conv cld liq into envrionment  - Revised macrophysics' )
    call addfld ('CLDICEDET', (/ 'lev' /), 'A', 'kg/kg/s', &
         'Detrainment of conv cld ice into envrionment  - Revised macrophysics' )
    call addfld ('CLDLIQLIM', (/ 'lev' /), 'A', 'kg/kg/s', 'CLDLIQ limiting tendency - Revised macrophysics'         )
    call addfld ('CLDICELIM', (/ 'lev' /), 'A', 'kg/kg/s', 'CLDICE limiting tendency - Revised macrophysics'         )

    call addfld ('AST',        (/ 'lev' /), 'A',       '1', 'Stratus cloud fraction')
    call addfld ('LIQCLDF',        (/ 'lev' /), 'A',   '1', 'Stratus Liquid cloud fraction')
    call add_default('LIQCLDF', 1, ' ')
    call addfld ('ICECLDF',        (/ 'lev' /), 'A',   '1', 'Stratus ICE cloud fraction')

    call addfld ('CLDST', (/ 'lev' /), 'A', 'fraction', 'Stratus cloud fraction'                                  )
    call addfld ('CONCLD', (/ 'lev' /), 'A', 'fraction', 'Convective cloud cover'                                  )
 

    call addfld ('CLR_LIQ', (/ 'lev' /), 'A', 'fraction', 'Clear sky fraction for liquid stratus'           )
    call addfld ('CLR_ICE', (/ 'lev' /), 'A', 'fraction', 'Clear sky fraction for ice stratus'              )

    call addfld ('CLDLIQSTR', (/ 'lev' /), 'A', 'kg/kg', 'Stratiform CLDLIQ'                                  )
    call addfld ('CLDICESTR', (/ 'lev' /), 'A', 'kg/kg', 'Stratiform CLDICE'                                  )
    call addfld ('CLDLIQCON', (/ 'lev' /), 'A', 'kg/kg', 'Convective CLDLIQ'                                  )
    call addfld ('CLDICECON', (/ 'lev' /), 'A', 'kg/kg', 'Convective CLDICE'                                  )

    call addfld ('CLDSICE', (/ 'lev' /), 'A', 'kg/kg', 'CloudSat equivalent ice mass mixing ratio'               )

    call addfld ('CMELIQ', (/ 'lev' /), 'A', 'kg/kg/s', 'Rate of cond-evap of liq within the cloud'               )

    call addfld ('TTENDICE', (/ 'lev' /), 'A',      'K/s', 'T tendency from Ice Saturation Adjustment'       )
    call addfld ('QVTENDICE', (/ 'lev' /), 'A', 'kg/kg/s', 'Q tendency from Ice Saturation Adjustment'       )
    call addfld ('QITENDICE', (/ 'lev' /), 'A', 'kg/kg/s', 'CLDICE tendency from Ice Saturation Adjustment'       )
    call addfld ('NITENDICE', (/ 'lev' /), 'A', 'kg/kg/s', 'NUMICE tendency from Ice Saturation Adjustment'       )
    if ( history_budget ) then

          call add_default ('DPDLFLIQ ', history_budget_histfile_num, ' ')
          call add_default ('DPDLFICE ', history_budget_histfile_num, ' ')
          call add_default ('SHDLFLIQ ', history_budget_histfile_num, ' ')
          call add_default ('SHDLFICE ', history_budget_histfile_num, ' ')
          call add_default ('DPDLFT   ', history_budget_histfile_num, ' ')
          call add_default ('SHDLFT   ', history_budget_histfile_num, ' ')
          call add_default ('ZMDLF    ', history_budget_histfile_num, ' ')

          call add_default ('MACPDT   ', history_budget_histfile_num, ' ')
          call add_default ('MACPDQ   ', history_budget_histfile_num, ' ')
          call add_default ('MACPDLIQ ', history_budget_histfile_num, ' ')
          call add_default ('MACPDICE ', history_budget_histfile_num, ' ')
 
          call add_default ('CLDVAPADJ', history_budget_histfile_num, ' ')
          call add_default ('CLDLIQLIM', history_budget_histfile_num, ' ')
          call add_default ('CLDLIQDET', history_budget_histfile_num, ' ')
          call add_default ('CLDLIQADJ', history_budget_histfile_num, ' ')
          call add_default ('CLDICELIM', history_budget_histfile_num, ' ')
          call add_default ('CLDICEDET', history_budget_histfile_num, ' ')
          call add_default ('CLDICEADJ', history_budget_histfile_num, ' ')

          call add_default ('CMELIQ   ', history_budget_histfile_num, ' ')

    end if

    ! Get constituent indices
    call cnst_get_ind('CLDLIQ', ixcldliq)
    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('NUMLIQ', ixnumliq)
    call cnst_get_ind('NUMICE', ixnumice)

    ! Get physics buffer indices
    CC_T_idx    = pbuf_get_index('CC_T')
    CC_qv_idx   = pbuf_get_index('CC_qv')
    CC_ql_idx   = pbuf_get_index('CC_ql')
    CC_qi_idx   = pbuf_get_index('CC_qi')
    CC_nl_idx   = pbuf_get_index('CC_nl')
    CC_ni_idx   = pbuf_get_index('CC_ni')
    CC_qlst_idx = pbuf_get_index('CC_qlst')

    if (micro_do_icesupersat) then 
       naai_idx      = pbuf_get_index('NAAI')
    endif

    if (rhminl_opt > 0 .or. rhmini_opt > 0) then
       cmfr_det_idx = pbuf_get_index('cmfr_det', istat)
       if (istat < 0) call endrun(subname//': macrop option requires cmfr_det in pbuf')
       if (rhminl_opt > 0) then
          qlr_det_idx  = pbuf_get_index('qlr_det', istat)
          if (istat < 0) call endrun(subname//': macrop option requires qlr_det in pbuf')
       end if
       if (rhmini_opt > 0) then
          qir_det_idx  = pbuf_get_index('qir_det', istat)
          if (istat < 0) call endrun(subname//': macrop option requires qir_det in pbuf')
       end if
    end if

    if (rhminl_opt == 2 .or. rhmini_opt == 2) then
       tke_idx = pbuf_get_index('tke')
       if (rhminl_opt == 2) then
          qtl_flx_idx = pbuf_get_index('qtl_flx', istat)
          if (istat < 0) call endrun(subname//': macrop option requires qtl_flx in pbuf')
       end if
       if (rhmini_opt == 2) then
          qti_flx_idx = pbuf_get_index('qti_flx', istat)
          if (istat < 0) call endrun(subname//': macrop option requires qti_flx in pbuf')
       end if
    end if

    ! Init pbuf fields.  Note that the fields CLD, CONCLD, QCWAT, LCWAT, 
    ! ICCWAT, and TCWAT are initialized in phys_inidat.
    if (is_first_step()) then
       call pbuf_set_field(pbuf2d, ast_idx,    0._r8)
       call pbuf_set_field(pbuf2d, aist_idx,   0._r8)
       call pbuf_set_field(pbuf2d, alst_idx,   0._r8)
       call pbuf_set_field(pbuf2d, qist_idx,   0._r8)
       call pbuf_set_field(pbuf2d, qlst_idx,   0._r8)
       call pbuf_set_field(pbuf2d, nlwat_idx,  0._r8)
       call pbuf_set_field(pbuf2d, niwat_idx,  0._r8)
       call pbuf_set_field(pbuf2d, fice_idx,   0._r8)
       call pbuf_set_field(pbuf2d, cmeliq_idx, 0._r8)
    end if

  end subroutine macrop_driver_init

  !============================================================================ !
  !                                                                             !
  !============================================================================ !


  subroutine macrop_driver_tend(                             &
             state, ptend, dtime, landfrac,  &
             ocnfrac,  snowh,                       &
             dlf, dlf2, cmfmc, cmfmc2, ts,          &
             sst, zdu,       &
             pbuf, &
             det_s, det_ice, alst_o)

  !-------------------------------------------------------- !  
  !                                                         ! 
  ! Purpose:                                                !
  !                                                         !
  ! Interface to detrain, cloud fraction and                !
  !     cloud macrophysics subroutines                      !
  !                                                         ! 
  ! Author: A. Gettelman, C. Craig, Oct 2010                !
  ! based on stratiform_tend by D.B. Coleman 4/2010         !
  !                                                         !
  !-------------------------------------------------------- !

  use cloud_fraction,   only: cldfrc, cldfrc_fice
  use physics_types,    only: physics_state, physics_ptend
  use physics_types,    only: physics_ptend_init, physics_update
  use physics_types,    only: physics_ptend_sum,  physics_state_copy
  use physics_types,    only: physics_state_dealloc
  use cam_history,      only: outfld
  use constituents,     only: cnst_get_ind, pcnst
  use cldwat2m_macro,   only: mmacro_pcond
  use physconst,        only: cpair, tmelt, gravit
  use time_manager,     only: get_nstep, is_first_step 

  use ref_pres,         only: top_lev => trop_cloud_top_lev

  !
  ! Input arguments
  !

  type(physics_state), intent(in)    :: state       ! State variables
  type(physics_ptend), intent(out)   :: ptend       ! macrophysics parameterization tendencies
  type(physics_buffer_desc), pointer :: pbuf(:)     ! Physics buffer

  real(r8), intent(in)  :: dtime                    ! Timestep
  real(r8), intent(in)  :: landfrac(pcols)          ! Land fraction (fraction)
  real(r8), intent(in)  :: ocnfrac (pcols)          ! Ocean fraction (fraction)
  real(r8), intent(in)  :: snowh(pcols)             ! Snow depth over land, water equivalent (m)
  real(r8), intent(in)  :: dlf(pcols,pver)          ! Detrained water from convection schemes
  real(r8), intent(in)  :: dlf2(pcols,pver)         ! Detrained water from shallow convection scheme
  real(r8), intent(in)  :: cmfmc(pcols,pverp)       ! Deep + Shallow Convective mass flux [ kg /s/m^2 ]
  real(r8), intent(in)  :: cmfmc2(pcols,pverp)      ! Shallow convective mass flux [ kg/s/m^2 ]

  real(r8), intent(in)  :: ts(pcols)                ! Surface temperature
  real(r8), intent(in)  :: sst(pcols)               ! Sea surface temperature
  real(r8), intent(in)  :: zdu(pcols,pver)          ! Detrainment rate from deep convection


  ! These two variables are needed for energy check    
  real(r8), intent(out) :: det_s(pcols)             ! Integral of detrained static energy from ice
  real(r8), intent(out) :: det_ice(pcols)           ! Integral of detrained ice for energy check

  real(r8), intent(out) :: alst_o(pcols,pver)  ! HW: for old liquid status fraction 

  !
  ! Local variables
  !

  type(physics_state) :: state_loc                  ! Local copy of the state variable
  type(physics_ptend) :: ptend_loc                  ! Local parameterization tendencies

  integer i,k
  integer :: lchnk                                  ! Chunk identifier
  integer :: ncol                                   ! Number of atmospheric columns

  ! Physics buffer fields

  integer itim_old
  real(r8), pointer, dimension(:,:) :: qcwat        ! Cloud water old q
  real(r8), pointer, dimension(:,:) :: tcwat        ! Cloud water old temperature
  real(r8), pointer, dimension(:,:) :: lcwat        ! Cloud liquid water old q
  real(r8), pointer, dimension(:,:) :: iccwat       ! Cloud ice water old q
  real(r8), pointer, dimension(:,:) :: nlwat        ! Cloud liquid droplet number condentration. old.
  real(r8), pointer, dimension(:,:) :: niwat        ! Cloud ice    droplet number condentration. old.
  real(r8), pointer, dimension(:,:) :: CC_T         ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_qv        ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_ql        ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_qi        ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_nl        ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_ni        ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_qlst      ! In-liquid stratus microphysical tendency
  real(r8), pointer, dimension(:,:) :: cld          ! Total cloud fraction
  real(r8), pointer, dimension(:,:) :: ast          ! Relative humidity cloud fraction
  real(r8), pointer, dimension(:,:) :: aist         ! Physical ice stratus fraction
  real(r8), pointer, dimension(:,:) :: alst         ! Physical liquid stratus fraction
  real(r8), pointer, dimension(:,:) :: qist         ! Physical in-cloud IWC
  real(r8), pointer, dimension(:,:) :: qlst         ! Physical in-cloud LWC
  real(r8), pointer, dimension(:,:) :: concld       ! Convective cloud fraction

  real(r8), pointer, dimension(:,:) :: shfrc        ! Cloud fraction from shallow convection scheme

  real(r8), pointer, dimension(:,:) :: cmeliq

  real(r8), pointer, dimension(:,:) :: tke
  real(r8), pointer, dimension(:,:) :: qtl_flx
  real(r8), pointer, dimension(:,:) :: qti_flx
  real(r8), pointer, dimension(:,:) :: cmfr_det
  real(r8), pointer, dimension(:,:) :: qlr_det
  real(r8), pointer, dimension(:,:) :: qir_det
  real(r8), pointer, dimension(:,:) :: rqcr_l
  real(r8), pointer, dimension(:,:) :: rqcr_i
  real(r8), pointer, dimension(:,:) :: rncr_l
  real(r8), pointer, dimension(:,:) :: rncr_i

  ! Convective cloud to the physics buffer for purposes of ql contrib. to radn.

  real(r8), pointer, dimension(:,:) :: fice_ql      ! Cloud ice/water partitioning ratio.

  real(r8), pointer, dimension(:,:) :: naai         ! Number concentration of activated ice nuclei
 
  real(r8) :: latsub

  ! tendencies for ice saturation adjustment
  real(r8)  :: stend(pcols,pver)
  real(r8)  :: qvtend(pcols,pver)
  real(r8)  :: qitend(pcols,pver)
  real(r8)  :: initend(pcols,pver)

  ! Local variables for cldfrc

  real(r8)  cldst(pcols,pver)                       ! Stratus cloud fraction
  real(r8)  rhcloud(pcols,pver)                     ! Relative humidity cloud (last timestep)
  real(r8)  clc(pcols)                              ! Column convective cloud amount
  real(r8)  rhu00(pcols,pver)                       ! RH threshold for cloud
  real(r8)  icecldf(pcols,pver)                     ! Ice cloud fraction
  real(r8)  liqcldf(pcols,pver)                     ! Liquid cloud fraction (combined into cloud)
  real(r8)  relhum(pcols,pver)                      ! RH, output to determine drh/da
  real(r8)  adet(pcols,pver)                        ! Cloud fraction generated by detrained convective condensate

  ! Local variables for macrophysics

  real(r8)  rdtime                                  ! 1./dtime
  real(r8)  qtend(pcols,pver)                       ! Moisture tendencies
  real(r8)  ttend(pcols,pver)                       ! Temperature tendencies
  real(r8)  ltend(pcols,pver)                       ! Cloud liquid water tendencies
  real(r8)  fice(pcols,pver)                        ! Fractional ice content within cloud
  real(r8)  fsnow(pcols,pver)                       ! Fractional snow production
  real(r8)  homoo(pcols,pver)  
  real(r8)  qcreso(pcols,pver)  
  real(r8)  prcio(pcols,pver)  
  real(r8)  praio(pcols,pver)  
  real(r8)  qireso(pcols,pver)
  real(r8)  ftem(pcols,pver)
  real(r8)  pracso (pcols,pver) 
  real(r8)  dpdlfliq(pcols,pver)
  real(r8)  dpdlfice(pcols,pver)
  real(r8)  shdlfliq(pcols,pver)
  real(r8)  shdlfice(pcols,pver)
  real(r8)  dpdlft  (pcols,pver)
  real(r8)  shdlft  (pcols,pver)

  real(r8)  dum1
  real(r8)  qc(pcols,pver)
  real(r8)  qi(pcols,pver)
  real(r8)  nc(pcols,pver)
  real(r8)  ni(pcols,pver)

  logical   lq(pcnst)

  ! Output from mmacro_pcond

  real(r8)  tlat(pcols,pver)
  real(r8)  qvlat(pcols,pver)
  real(r8)  qcten(pcols,pver)
  real(r8)  qiten(pcols,pver)
  real(r8)  ncten(pcols,pver)
  real(r8)  niten(pcols,pver)

  ! Output from mmacro_pcond

  real(r8)  qvadj(pcols,pver)                       ! Macro-physics adjustment tendency from "positive_moisture" call (vapor)
  real(r8)  qladj(pcols,pver)                       ! Macro-physics adjustment tendency from "positive_moisture" call (liquid)
  real(r8)  qiadj(pcols,pver)                       ! Macro-physics adjustment tendency from "positive_moisture" call (ice)
  real(r8)  qllim(pcols,pver)                       ! Macro-physics tendency from "instratus_condensate" call (liquid)
  real(r8)  qilim(pcols,pver)                       ! Macro-physics tendency from "instratus_condensate" call (ice)

  ! For revised macophysics, mmacro_pcond

  real(r8)  itend(pcols,pver)
  real(r8)  lmitend(pcols,pver)
  real(r8)  zeros(pcols,pver)
  real(r8)  t_inout(pcols,pver)
  real(r8)  qv_inout(pcols,pver)
  real(r8)  ql_inout(pcols,pver)
  real(r8)  qi_inout(pcols,pver)
  real(r8)  concld_old(pcols,pver)

  ! Note that below 'clr_old' is defined using 'alst_old' not 'ast_old' for full consistency with the 
  ! liquid condensation process which is using 'alst' not 'ast'. 
  ! For microconsistency use 'concld_old', since 'alst_old' was computed using 'concld_old'.
  ! Since convective updraft fractional area is small, it does not matter whether 'concld' or 'concld_old' is used.
  ! Note also that 'clri_old' is defined using 'ast_old' since current microphysics is operating on 'ast_old' 
  real(r8)  clrw_old(pcols,pver) ! (1 - concld_old - alst_old)
  real(r8)  clri_old(pcols,pver) ! (1 - concld_old -  ast_old)

  real(r8)  nl_inout(pcols,pver)
  real(r8)  ni_inout(pcols,pver)

  real(r8)  nltend(pcols,pver)
  real(r8)  nitend(pcols,pver)


  ! For detraining cumulus condensate into the 'stratus' without evaporation
  ! This is for use in mmacro_pcond

  real(r8)  dlf_T(pcols,pver)
  real(r8)  dlf_qv(pcols,pver)
  real(r8)  dlf_ql(pcols,pver)
  real(r8)  dlf_qi(pcols,pver)
  real(r8)  dlf_nl(pcols,pver)
  real(r8)  dlf_ni(pcols,pver)

  ! Local variables for CFMIP calculations
  real(r8) :: mr_lsliq(pcols,pver)  ! mixing_ratio_large_scale_cloud_liquid (kg/kg)
  real(r8) :: mr_lsice(pcols,pver)  ! mixing_ratio_large_scale_cloud_ice (kg/kg)
  real(r8) :: mr_ccliq(pcols,pver)  ! mixing_ratio_convective_cloud_liquid (kg/kg)
  real(r8) :: mr_ccice(pcols,pver)  ! mixing_ratio_convective_cloud_ice (kg/kg)

  ! CloudSat equivalent ice mass mixing ratio (kg/kg)
  real(r8) :: cldsice(pcols,pver)

  ! ======================================================================

  if (micro_do_icesupersat) then 
     call pbuf_get_field(pbuf, naai_idx, naai)
  endif

  lchnk = state%lchnk
  ncol  = state%ncol

  call physics_state_copy(state, state_loc)            ! Copy state to local state_loc.

  ! Associate pointers with physics buffer fields

  itim_old = pbuf_old_tim_idx()

  call pbuf_get_field(pbuf, qcwat_idx,   qcwat,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, tcwat_idx,   tcwat,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, lcwat_idx,   lcwat,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, iccwat_idx,  iccwat,  start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, nlwat_idx,   nlwat,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, niwat_idx,   niwat,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

  call pbuf_get_field(pbuf, cc_t_idx,    cc_t,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, cc_qv_idx,   cc_qv,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, cc_ql_idx,   cc_ql,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, cc_qi_idx,   cc_qi,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, cc_nl_idx,   cc_nl,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, cc_ni_idx,   cc_ni,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, cc_qlst_idx, cc_qlst, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

  call pbuf_get_field(pbuf, cld_idx,     cld,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, concld_idx,  concld, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, ast_idx,     ast,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, aist_idx,    aist,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, alst_idx,    alst,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, qist_idx,    qist,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, qlst_idx,    qlst,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

  call pbuf_get_field(pbuf, cmeliq_idx,  cmeliq)

! For purposes of convective ql.

  call pbuf_get_field(pbuf, fice_idx,     fice_ql )


  ! Initialize convective detrainment tendency

  dlf_T(:,:)  = 0._r8
  dlf_qv(:,:) = 0._r8
  dlf_ql(:,:) = 0._r8
  dlf_qi(:,:) = 0._r8
  dlf_nl(:,:) = 0._r8
  dlf_ni(:,:) = 0._r8

   ! ------------------------------------- !
   ! From here, process computation begins ! 
   ! ------------------------------------- !

   ! ----------------------------------------------------------------------------- !
   ! Detrainment of convective condensate into the environment or stratiform cloud !
   ! ----------------------------------------------------------------------------- !

   lq(:)        = .FALSE.
   lq(ixcldliq) = .TRUE.
   lq(ixcldice) = .TRUE.
   lq(ixnumliq) = .TRUE.
   lq(ixnumice) = .TRUE.
   call physics_ptend_init(ptend_loc, state%psetcols, 'pcwdetrain', ls=.true., lq=lq)   ! Initialize local physics_ptend object

     ! Procedures :
     ! (1) Partition detrained convective cloud water into liquid and ice based on T.
     !     This also involves heating.
     !     If convection scheme can handle this internally, this step is not necssary.
     ! (2) Assuming a certain effective droplet radius, computes number concentration
     !     of detrained convective cloud liquid and ice.
     ! (3) If 'cu_det_st = .true' ('false'), detrain convective cloud 'liquid' into 
     !     the pre-existing 'liquid' stratus ( mean environment ).  The former does
     !     not involve any macrophysical evaporation while the latter does. This is
     !     a kind of 'targetted' deposition. Then, force in-stratus LWC to be bounded 
     !     by qcst_min and qcst_max in mmacro_pcond.
     ! (4) In contrast to liquid, convective ice is detrained into the environment 
     !     and involved in the sublimation. Similar bounds as liquid stratus are imposed.
     ! This is the key procesure generating upper-level cirrus clouds.
     ! The unit of dlf : [ kg/kg/s ]

   det_s(:)   = 0._r8
   det_ice(:) = 0._r8

   dpdlfliq = 0._r8
   dpdlfice = 0._r8
   shdlfliq = 0._r8
   shdlfice = 0._r8
   dpdlft   = 0._r8
   shdlft   = 0._r8

   ! --------------------------------------------------------------------------------------------- !
   ! Mar.06.2015. In case of UNICON, the liquid and ice phases of detrained convective condensates !
   !              are separately treated, without converting into liquid.                          !
   !              Note that UNICON should always be used with 'cu_det_st = .true.'.                !
   ! --------------------------------------------------------------------------------------------- !

   ! ----------- !
   ! UNICON case !
   ! ----------- !

   if( shallow_scheme .eq. 'UNICON' ) then

       call pbuf_get_field(pbuf, rqcr_l_idx, rqcr_l)
       call pbuf_get_field(pbuf, rqcr_i_idx, rqcr_i)
       call pbuf_get_field(pbuf, rncr_l_idx, rncr_l)
       call pbuf_get_field(pbuf, rncr_i_idx, rncr_i)
       do k = top_lev, pver
       do i = 1, state_loc%ncol
          ptend_loc%q(i,k,ixcldliq) = rqcr_l(i,k)
          ptend_loc%q(i,k,ixcldice) = rqcr_i(i,k)
          ptend_loc%q(i,k,ixnumliq) = rncr_l(i,k)
          ptend_loc%q(i,k,ixnumice) = rncr_i(i,k)
          ptend_loc%s(i,k)          = 0._r8
          det_ice(i)  = det_ice(i) - ptend_loc%q(i,k,ixcldice)*state_loc%pdel(i,k)/gravit
          dlf_T(i,k)  = ptend_loc%s(i,k)/cpair
          dlf_qv(i,k) = 0._r8
          dlf_ql(i,k) = ptend_loc%q(i,k,ixcldliq)
          dlf_qi(i,k) = ptend_loc%q(i,k,ixcldice)
          dlf_nl(i,k) = ptend_loc%q(i,k,ixnumliq)
          dlf_ni(i,k) = ptend_loc%q(i,k,ixnumice)
          ptend_loc%q(i,k,ixcldliq) = 0._r8
          ptend_loc%q(i,k,ixcldice) = 0._r8
          ptend_loc%q(i,k,ixnumliq) = 0._r8
          ptend_loc%q(i,k,ixnumice) = 0._r8
          ptend_loc%s(i,k)          = 0._r8
          dpdlfliq(i,k)             = 0._r8
          dpdlfice(i,k)             = 0._r8
          shdlfliq(i,k)             = 0._r8
          shdlfice(i,k)             = 0._r8
          dpdlft  (i,k)             = 0._r8
          shdlft  (i,k)             = 0._r8
       end do
       end do

   else

   ! -------------------------- !
   ! THE OTHER CONVECTION CASES !
   ! -------------------------- !

   do k = top_lev, pver
   do i = 1, state_loc%ncol
      if( state_loc%t(i,k) > 268.15_r8 ) then
          dum1 = 0.0_r8
      elseif( state_loc%t(i,k) < 238.15_r8 ) then
          dum1 = 1.0_r8
      else
          dum1 = ( 268.15_r8 - state_loc%t(i,k) ) / 30._r8
      endif

     ! If detrainment was done elsewhere, still update the variables used for output
     ! assuming that the temperature split between liquid and ice is the same as assumed
     ! here.
     if (do_detrain) then
      ptend_loc%q(i,k,ixcldliq) = dlf(i,k) * ( 1._r8 - dum1 )
      ptend_loc%q(i,k,ixcldice) = dlf(i,k) * dum1
    ! dum2                      = dlf(i,k) * ( 1._r8 - dum1 )
      ptend_loc%q(i,k,ixnumliq) = 3._r8 * ( max(0._r8, ( dlf(i,k) - dlf2(i,k) )) * ( 1._r8 - dum1 ) ) / &
           (4._r8*3.14_r8* 8.e-6_r8**3*997._r8) + & ! Deep    Convection
           3._r8 * (                         dlf2(i,k)    * ( 1._r8 - dum1 ) ) / &
           (4._r8*3.14_r8*10.e-6_r8**3*997._r8)     ! Shallow Convection 
    ! dum2                      = dlf(i,k) * dum1
      ptend_loc%q(i,k,ixnumice) = 3._r8 * ( max(0._r8, ( dlf(i,k) - dlf2(i,k) )) *  dum1 ) / &
           (4._r8*3.14_r8*25.e-6_r8**3*500._r8) + & ! Deep    Convection
           3._r8 * (                         dlf2(i,k)    *  dum1 ) / &
           (4._r8*3.14_r8*50.e-6_r8**3*500._r8)     ! Shallow Convection
      ptend_loc%s(i,k)          = dlf(i,k) * dum1 * latice
     else 
        ptend_loc%q(i,k,ixcldliq) = 0._r8
        ptend_loc%q(i,k,ixcldice) = 0._r8
        ptend_loc%q(i,k,ixnumliq) = 0._r8
        ptend_loc%q(i,k,ixnumice) = 0._r8
        ptend_loc%s(i,k)          = 0._r8
     end if
    

    ! Only rliq is saved from deep convection, which is the reserved liquid.  We need to keep
    !   track of the integrals of ice and static energy that is effected from conversion to ice
    !   so that the energy checker doesn't complain.
      det_s(i)                  = det_s(i) + ptend_loc%s(i,k)*state_loc%pdel(i,k)/gravit
      det_ice(i)                = det_ice(i) - ptend_loc%q(i,k,ixcldice)*state_loc%pdel(i,k)/gravit      

    ! Targetted detrainment of convective liquid water either directly into the
    ! existing liquid stratus or into the environment. 
      if( cu_det_st ) then
          dlf_T(i,k)  = ptend_loc%s(i,k)/cpair
          dlf_qv(i,k) = 0._r8
          dlf_ql(i,k) = ptend_loc%q(i,k,ixcldliq)
          dlf_qi(i,k) = ptend_loc%q(i,k,ixcldice)
          dlf_nl(i,k) = ptend_loc%q(i,k,ixnumliq)
          dlf_ni(i,k) = ptend_loc%q(i,k,ixnumice)         
          ptend_loc%q(i,k,ixcldliq) = 0._r8
          ptend_loc%q(i,k,ixcldice) = 0._r8
          ptend_loc%q(i,k,ixnumliq) = 0._r8
          ptend_loc%q(i,k,ixnumice) = 0._r8
          ptend_loc%s(i,k)          = 0._r8
          dpdlfliq(i,k)             = 0._r8
          dpdlfice(i,k)             = 0._r8
          shdlfliq(i,k)             = 0._r8
          shdlfice(i,k)             = 0._r8
          dpdlft  (i,k)             = 0._r8
          shdlft  (i,k)             = 0._r8
       else
          dpdlfliq(i,k) = ( dlf(i,k) - dlf2(i,k) ) * ( 1._r8 - dum1 )
          dpdlfice(i,k) = ( dlf(i,k) - dlf2(i,k) ) * ( dum1 )
          shdlfliq(i,k) = dlf2(i,k) * ( 1._r8 - dum1 )
          shdlfice(i,k) = dlf2(i,k) * ( dum1 )
          dpdlft  (i,k) = ( dlf(i,k) - dlf2(i,k) ) * dum1 * latice/cpair
          shdlft  (i,k) = dlf2(i,k) * dum1 * latice/cpair
      endif
   end do
   end do

   endif

   ! ------------------------ !
   ! END OF CONVECTION CHOICE !
   ! ------------------------ !

   call outfld( 'DPDLFLIQ ', dpdlfliq, pcols, lchnk )
   call outfld( 'DPDLFICE ', dpdlfice, pcols, lchnk )
   call outfld( 'SHDLFLIQ ', shdlfliq, pcols, lchnk )
   call outfld( 'SHDLFICE ', shdlfice, pcols, lchnk )
   call outfld( 'DPDLFT   ', dpdlft  , pcols, lchnk )
   call outfld( 'SHDLFT   ', shdlft  , pcols, lchnk )

   call outfld( 'ZMDLF',     dlf     , pcols, state_loc%lchnk )

   det_ice(:ncol) = det_ice(:ncol)/1000._r8  ! divide by density of water

   ! Add the detrainment tendency to the output tendency
   call physics_ptend_init(ptend, state%psetcols, 'macrop')
   call physics_ptend_sum(ptend_loc, ptend, ncol)

   ! update local copy of state with the detrainment tendency
   ! ptend_loc is reset to zero by this call
   call physics_update(state_loc, ptend_loc, dtime)

   if (micro_do_icesupersat) then 

      ! -------------------------------------- !
      ! Ice Saturation Adjustment Computation  !
      ! -------------------------------------- !

      lq(:)        = .FALSE.

      lq(1)        = .true.
      lq(ixcldice) = .true.
      lq(ixnumice) = .true.

      latsub = latvap + latice

      call physics_ptend_init(ptend_loc, state%psetcols, 'iceadj', ls=.true., lq=lq)

      stend(:ncol,:)=0._r8
      qvtend(:ncol,:)=0._r8
      qitend(:ncol,:)=0._r8
      initend(:ncol,:)=0._r8

      call ice_macro_tend(naai(:ncol,top_lev:pver),state%t(:ncol,top_lev:pver), &
           state%pmid(:ncol,top_lev:pver),state%q(:ncol,top_lev:pver,1),state%q(:ncol,top_lev:pver,ixcldice),&
           state%q(:ncol,top_lev:pver,ixnumice),latsub,dtime,&
           stend(:ncol,top_lev:pver),qvtend(:ncol,top_lev:pver),qitend(:ncol,top_lev:pver),&
           initend(:ncol,top_lev:pver))

      ! update local copy of state with the tendencies
      ptend_loc%q(:ncol,top_lev:pver,1)=qvtend(:ncol,top_lev:pver)
      ptend_loc%q(:ncol,top_lev:pver,ixcldice)=qitend(:ncol,top_lev:pver)  
      ptend_loc%q(:ncol,top_lev:pver,ixnumice)=initend(:ncol,top_lev:pver)
      ptend_loc%s(:ncol,top_lev:pver)=stend(:ncol,top_lev:pver) 

      ! Add the ice tendency to the output tendency
      call physics_ptend_sum(ptend_loc, ptend, ncol)
 
      ! ptend_loc is reset to zero by this call
      call physics_update(state_loc, ptend_loc, dtime)

      ! Write output for tendencies:
      call outfld( 'TTENDICE',  stend/cpair, pcols, lchnk )
      call outfld( 'QVTENDICE', qvtend, pcols, lchnk )
      call outfld( 'QITENDICE', qitend, pcols, lchnk )
      call outfld( 'NITENDICE', initend, pcols, lchnk )

   endif

   ! -------------------------------------- !
   ! Computation of Various Cloud Fractions !
   ! -------------------------------------- !

   ! ----------------------------------------------------------------------------- !
   ! Treatment of cloud fraction in CAM4 and CAM5 differs                          !  
   ! (1) CAM4                                                                      !
   !     . Cumulus AMT = Deep    Cumulus AMT ( empirical fcn of mass flux ) +      !
   !                     Shallow Cumulus AMT ( empirical fcn of mass flux )        !
   !     . Stratus AMT = max( RH stratus AMT, Stability Stratus AMT )              !
   !     . Cumulus and Stratus are 'minimally' overlapped without hierarchy.       !
   !     . Cumulus LWC,IWC is assumed to be the same as Stratus LWC,IWC            !
   ! (2) CAM5                                                                      !
   !     . Cumulus AMT = Deep    Cumulus AMT ( empirical fcn of mass flux ) +      !
   !                     Shallow Cumulus AMT ( internally fcn of mass flux and w ) !
   !     . Stratus AMT = fcn of environmental-mean RH ( no Stability Stratus )     !
   !     . Cumulus and Stratus are non-overlapped with higher priority on Cumulus  !
   !     . Cumulus ( both Deep and Shallow ) has its own LWC and IWC.              !
   ! ----------------------------------------------------------------------------- ! 

   concld_old(:ncol,top_lev:pver) = concld(:ncol,top_lev:pver)

! JHYoon:  Sungsu added below for use in the cloud macrophysics
! May.06.2015. Ideally, this should be used only for UNICON block.

   nullify(tke, qtl_flx, qti_flx, cmfr_det, qlr_det, qir_det)
   if (tke_idx      > 0) call pbuf_get_field(pbuf, tke_idx, tke)
   if (qtl_flx_idx  > 0) call pbuf_get_field(pbuf, qtl_flx_idx,  qtl_flx)
   if (qti_flx_idx  > 0) call pbuf_get_field(pbuf, qti_flx_idx,  qti_flx)
   if (cmfr_det_idx > 0) call pbuf_get_field(pbuf, cmfr_det_idx, cmfr_det)
   if (qlr_det_idx  > 0) call pbuf_get_field(pbuf, qlr_det_idx,  qlr_det)
   if (qir_det_idx  > 0) call pbuf_get_field(pbuf, qir_det_idx,  qir_det)

   clrw_old(:ncol,:top_lev-1) = 0._r8
   clri_old(:ncol,:top_lev-1) = 0._r8
   do k = top_lev, pver
      do i = 1, ncol
         clrw_old(i,k) = max( 0._r8, min( 1._r8, 1._r8 - concld(i,k) - alst(i,k) ) )      
         clri_old(i,k) = max( 0._r8, min( 1._r8, 1._r8 - concld(i,k) -  ast(i,k) ) )      
      end do
   end do
! JHYoon

   if( use_shfrc ) then
       call pbuf_get_field(pbuf, shfrc_idx, shfrc )
   else 
       allocate(shfrc(pcols,pver))
       shfrc(:,:) = 0._r8
   endif

   ! CAM5 only uses 'concld' output from the below subroutine. 
   ! Stratus ('ast' = max(alst,aist)) and total cloud fraction ('cld = ast + concld')
   ! will be computed using this updated 'concld' in the stratiform macrophysics 
   ! scheme (mmacro_pcond) later below. 

   call t_startf("cldfrc")

   call cldfrc( lchnk, ncol, pbuf,                                                 &
                state_loc%pmid, state_loc%t, state_loc%q(:,:,1), state_loc%omega,  &
                state_loc%phis, shfrc, use_shfrc,                                  &
                cld, rhcloud, clc, state_loc%pdel,                                 &
                cmfmc, cmfmc2, landfrac,snowh, concld, cldst,                      &
                ts, sst, state_loc%pint(:,pverp), zdu, ocnfrac, rhu00,             &
                state_loc%q(:,:,ixcldice), icecldf, liqcldf,                       &
                relhum, 0 )

   call t_stopf("cldfrc")

   ! ---------------------------------------------- !
   ! Stratiform Cloud Macrophysics and Microphysics !
   ! ---------------------------------------------- !

   lchnk  = state_loc%lchnk
   ncol   = state_loc%ncol
   rdtime = 1._r8/dtime

 ! Define fractional amount of stratus condensate and precipitation in ice phase.
 ! This uses a ramp ( -30 ~ -10 for fice, -5 ~ 0 for fsnow ). 
 ! The ramp within convective cloud may be different

   call cldfrc_fice( ncol, state_loc%t, fice, fsnow )


   lq(:)        = .FALSE.

   lq(1)        = .true.
   lq(ixcldice) = .true.
   lq(ixcldliq) = .true.

   lq(ixnumliq) = .true.
   lq(ixnumice) = .true.

   ! Initialize local physics_ptend object again
   call physics_ptend_init(ptend_loc, state%psetcols, 'macro_park', &
        ls=.true., lq=lq )  

 ! --------------------------------- !
 ! Liquid Macrop_Driver Macrophysics !
 ! --------------------------------- !

   call t_startf('mmacro_pcond')

   zeros(:ncol,top_lev:pver)  = 0._r8
   qc(:ncol,top_lev:pver) = state_loc%q(:ncol,top_lev:pver,ixcldliq)
   qi(:ncol,top_lev:pver) = state_loc%q(:ncol,top_lev:pver,ixcldice)
   nc(:ncol,top_lev:pver) = state_loc%q(:ncol,top_lev:pver,ixnumliq)
   ni(:ncol,top_lev:pver) = state_loc%q(:ncol,top_lev:pver,ixnumice)

 ! In CAM5, 'microphysical forcing' ( CC_... ) and 'the other advective forcings' ( ttend, ... ) 
 ! are separately provided into the prognostic microp_driver macrophysics scheme. This is an
 ! attempt to resolve in-cloud and out-cloud forcings. 

   if( get_nstep() .le. 1 ) then
       tcwat(:ncol,:)   = state_loc%t(:ncol,:)
       qcwat(:ncol,:)   = state_loc%q(:ncol,:,1)
       lcwat(:ncol,:)   = qc(:ncol,:) + qi(:ncol,:)
       iccwat(:ncol,:)  = qi(:ncol,:)
       nlwat(:ncol,:)   = nc(:ncol,:)
       niwat(:ncol,:)   = ni(:ncol,:)
       ttend(:ncol,:)   = 0._r8
       qtend(:ncol,:)   = 0._r8
       ltend(:ncol,:)   = 0._r8
       itend(:ncol,:)   = 0._r8
       nltend(:ncol,:)  = 0._r8
       nitend(:ncol,:)  = 0._r8
       CC_T(:ncol,:)    = 0._r8
       CC_qv(:ncol,:)   = 0._r8
       CC_ql(:ncol,:)   = 0._r8
       CC_qi(:ncol,:)   = 0._r8
       CC_nl(:ncol,:)   = 0._r8
       CC_ni(:ncol,:)   = 0._r8
       CC_qlst(:ncol,:) = 0._r8
   else
       ttend(:ncol,top_lev:pver)   = ( state_loc%t(:ncol,top_lev:pver)   -  tcwat(:ncol,top_lev:pver)) * rdtime &
            - CC_T(:ncol,top_lev:pver) 
       qtend(:ncol,top_lev:pver)   = ( state_loc%q(:ncol,top_lev:pver,1) -  qcwat(:ncol,top_lev:pver)) * rdtime &
            - CC_qv(:ncol,top_lev:pver)
       ltend(:ncol,top_lev:pver)   = ( qc(:ncol,top_lev:pver) + qi(:ncol,top_lev:pver) - lcwat(:ncol,top_lev:pver) ) * rdtime &
            - (CC_ql(:ncol,top_lev:pver) + CC_qi(:ncol,top_lev:pver))
       itend(:ncol,top_lev:pver)   = ( qi(:ncol,top_lev:pver)         - iccwat(:ncol,top_lev:pver)) * rdtime &
            - CC_qi(:ncol,top_lev:pver)
       nltend(:ncol,top_lev:pver)  = ( nc(:ncol,top_lev:pver)         -  nlwat(:ncol,top_lev:pver)) * rdtime &
            - CC_nl(:ncol,top_lev:pver)
       nitend(:ncol,top_lev:pver)  = ( ni(:ncol,top_lev:pver)         -  niwat(:ncol,top_lev:pver)) * rdtime &
            - CC_ni(:ncol,top_lev:pver)
   endif
   lmitend(:ncol,top_lev:pver) = ltend(:ncol,top_lev:pver) - itend(:ncol,top_lev:pver)

   if(liqcf_fix) then 
      !HW: set alst to alst_o before getting updated from macro
      if(.not.is_first_step())alst_o(:ncol,:pver) = alst(:ncol,:pver)
      !HW
   endif

   t_inout(:ncol,top_lev:pver)  =  tcwat(:ncol,top_lev:pver) 
   qv_inout(:ncol,top_lev:pver) =  qcwat(:ncol,top_lev:pver)
   ql_inout(:ncol,top_lev:pver) =  lcwat(:ncol,top_lev:pver) - iccwat(:ncol,top_lev:pver)
   qi_inout(:ncol,top_lev:pver) = iccwat(:ncol,top_lev:pver)
   nl_inout(:ncol,top_lev:pver) =  nlwat(:ncol,top_lev:pver)
   ni_inout(:ncol,top_lev:pver) =  niwat(:ncol,top_lev:pver)

 ! Liquid Microp_Driver Macrophysics.
 ! The main roles of this subroutines are
 ! (1) compute net condensation rate of stratiform liquid ( cmeliq )
 ! (2) compute liquid stratus and ice stratus fractions. 
 ! Note 'ttend...' are advective tendencies except microphysical process while
 !      'CC...'    are microphysical tendencies. 

   call mmacro_pcond( lchnk, ncol, dtime, state_loc%pmid, state_loc%pdel,        &
                      t_inout, qv_inout, ql_inout, qi_inout, nl_inout, ni_inout, &                  
                      ttend, qtend, lmitend, itend, nltend, nitend,              &
                      CC_T, CC_qv, CC_ql, CC_qi, CC_nl, CC_ni, CC_qlst,          & 
                      dlf_T, dlf_qv, dlf_ql, dlf_qi, dlf_nl, dlf_ni,             &
                      concld_old, concld, clrw_old, clri_old, landfrac, snowh,   &
                      tke, qtl_flx, qti_flx, cmfr_det, qlr_det, qir_det,         &
                      tlat, qvlat, qcten, qiten, ncten, niten,                   &
                      cmeliq, qvadj, qladj, qiadj, qllim, qilim,                 &
                      cld, alst, aist, qlst, qist, adet, do_cldice ) 

   if(liqcf_fix) then 
      if(is_first_step())alst_o(:ncol,:pver) = alst(:ncol,:pver)
   endif


 ! Copy of concld/fice to put in physics buffer
 ! Below are used only for convective cloud.

   fice_ql(:ncol,:top_lev-1)     = 0._r8
   fice_ql(:ncol,top_lev:pver)   = fice(:ncol,top_lev:pver)


 ! Compute net stratus fraction using maximum over-lapping assumption
   ast(:ncol,:top_lev-1) = 0._r8
   ast(:ncol,top_lev:pver) = max( alst(:ncol,top_lev:pver), aist(:ncol,top_lev:pver) )

   if (shallow_scheme .eq. 'UNICON') then
       if( i_adet .eq. 0 ) then
         ! Below is maximum overlap.
           ast(:ncol,top_lev:pver) = max( ast(:ncol,top_lev:pver), adet(:ncol,top_lev:pver) )
       else
         ! Below is separately added.
           ast(:ncol,top_lev:pver) = ast(:ncol,top_lev:pver) + adet(:ncol,top_lev:pver)
           ast(:ncol,top_lev:pver) = max( 0._r8, min( ast(:ncol,top_lev:pver), 1._r8 - concld(:ncol,top_lev:pver) ) )
       endif
   endif

   call t_stopf('mmacro_pcond')

   do k = top_lev, pver
      do i = 1, ncol
         ptend_loc%s(i,k)          =  tlat(i,k)
         ptend_loc%q(i,k,1)        = qvlat(i,k)
         ptend_loc%q(i,k,ixcldliq) = qcten(i,k)
         ptend_loc%q(i,k,ixcldice) = qiten(i,k)
         ptend_loc%q(i,k,ixnumliq) = ncten(i,k)
         ptend_loc%q(i,k,ixnumice) = niten(i,k)

         ! Check to make sure that the macrophysics code is respecting the flags that control
         ! whether cldwat should be prognosing cloud ice and cloud liquid or not.
         if ((.not. do_cldice) .and. (qiten(i,k) /= 0.0_r8)) then 
            call endrun("macrop_driver:ERROR - "// &
                 "Cldwat is configured not to prognose cloud ice, but mmacro_pcond has ice mass tendencies.")
         end if
         if ((.not. do_cldice) .and. (niten(i,k) /= 0.0_r8)) then 
            call endrun("macrop_driver:ERROR -"// &
                 " Cldwat is configured not to prognose cloud ice, but mmacro_pcond has ice number tendencies.")
         end if

         if ((.not. do_cldliq) .and. (qcten(i,k) /= 0.0_r8)) then 
            call endrun("macrop_driver:ERROR - "// &
                 "Cldwat is configured not to prognose cloud liquid, but mmacro_pcond has liquid mass tendencies.")
         end if
         if ((.not. do_cldliq) .and. (ncten(i,k) /= 0.0_r8)) then 
            call endrun("macrop_driver:ERROR - "// &
                 "Cldwat is configured not to prognose cloud liquid, but mmacro_pcond has liquid number tendencies.")
         end if
      end do
   end do

   ! update the output tendencies with the mmacro_pcond tendencies
   call physics_ptend_sum(ptend_loc, ptend, ncol)

   ! state_loc is the equlibrium state after macrophysics
   call physics_update(state_loc, ptend_loc, dtime)

   call outfld('CLR_LIQ', clrw_old,  pcols, lchnk)
   call outfld('CLR_ICE', clri_old,  pcols, lchnk)

   call outfld( 'MACPDT   ', tlat ,  pcols, lchnk )
   call outfld( 'MACPDQ   ', qvlat,  pcols, lchnk )
   call outfld( 'MACPDLIQ ', qcten,  pcols, lchnk )
   call outfld( 'MACPDICE ', qiten,  pcols, lchnk )
   call outfld( 'CLDVAPADJ', qvadj,  pcols, lchnk )
   call outfld( 'CLDLIQADJ', qladj,  pcols, lchnk )
   call outfld( 'CLDICEADJ', qiadj,  pcols, lchnk )
   call outfld( 'CLDLIQDET', dlf_ql, pcols, lchnk )
   call outfld( 'CLDICEDET', dlf_qi, pcols, lchnk )
   call outfld( 'CLDLIQLIM', qllim,  pcols, lchnk )
   call outfld( 'CLDICELIM', qilim,  pcols, lchnk )

   call outfld( 'ICECLDF ', aist,   pcols, lchnk )
   call outfld( 'LIQCLDF ', alst,   pcols, lchnk )
   call outfld( 'AST',      ast,    pcols, lchnk )   

   call outfld( 'CONCLD  ', concld, pcols, lchnk )
   call outfld( 'CLDST   ', cldst,  pcols, lchnk )

   call outfld( 'CMELIQ'  , cmeliq, pcols, lchnk )


   ! calculations and outfld calls for CLDLIQSTR, CLDICESTR, CLDLIQCON, CLDICECON for CFMIP

   ! initialize local variables
   mr_ccliq = 0._r8   !! not seen by radiation, so setting to 0 
   mr_ccice = 0._r8   !! not seen by radiation, so setting to 0
   mr_lsliq = 0._r8
   mr_lsice = 0._r8

   do k=top_lev,pver
      do i=1,ncol
         if (cld(i,k) .gt. 0._r8) then
            mr_lsliq(i,k) = state_loc%q(i,k,ixcldliq)
            mr_lsice(i,k) = state_loc%q(i,k,ixcldice)
         else
            mr_lsliq(i,k) = 0._r8
            mr_lsice(i,k) = 0._r8
         end if
      end do
   end do

   call outfld( 'CLDLIQSTR  ', mr_lsliq,    pcols, lchnk )
   call outfld( 'CLDICESTR  ', mr_lsice,    pcols, lchnk )
   call outfld( 'CLDLIQCON  ', mr_ccliq,    pcols, lchnk )
   call outfld( 'CLDICECON  ', mr_ccice,    pcols, lchnk )

   ! ------------------------------------------------- !
   ! Save equilibrium state variables for macrophysics !        
   ! at the next time step                             !
   ! ------------------------------------------------- !
   cldsice = 0._r8
   do k = top_lev, pver
      tcwat(:ncol,k)  = state_loc%t(:ncol,k)
      qcwat(:ncol,k)  = state_loc%q(:ncol,k,1)
      lcwat(:ncol,k)  = state_loc%q(:ncol,k,ixcldliq) + state_loc%q(:ncol,k,ixcldice)
      iccwat(:ncol,k) = state_loc%q(:ncol,k,ixcldice)
      nlwat(:ncol,k)  = state_loc%q(:ncol,k,ixnumliq)
      niwat(:ncol,k)  = state_loc%q(:ncol,k,ixnumice)
      cldsice(:ncol,k) = lcwat(:ncol,k) * min(1.0_r8, max(0.0_r8, (tmelt - tcwat(:ncol,k)) / 20._r8))
   end do

   call outfld( 'CLDSICE'    , cldsice,   pcols, lchnk )

   ! ptend_loc is deallocated in physics_update above
   call physics_state_dealloc(state_loc)

end subroutine macrop_driver_tend

! Saturation adjustment for ice
! Add ice mass if supersaturated
elemental subroutine ice_macro_tend(naai,t,p,qv,qi,ni,xxls,deltat,stend,qvtend,qitend,nitend) 

  use wv_sat_methods, only: wv_sat_qsat_ice

  real(r8), intent(in)  :: naai   !Activated number of ice nuclei 
  real(r8), intent(in)  :: t      !temperature (k)
  real(r8), intent(in)  :: p      !pressure (pa0
  real(r8), intent(in)  :: qv     !water vapor mixing ratio
  real(r8), intent(in)  :: qi     !ice mixing ratio
  real(r8), intent(in)  :: ni     !ice number concentration
  real(r8), intent(in)  :: xxls   !latent heat of sublimation
  real(r8), intent(in)  :: deltat !timestep
  real(r8), intent(out) :: stend  ! 'temperature' tendency 
  real(r8), intent(out) :: qvtend !vapor tendency
  real(r8), intent(out) :: qitend !ice mass tendency
  real(r8), intent(out) :: nitend !ice number tendency  
 
  real(r8) :: ESI
  real(r8) :: QSI
  real(r8) :: tau
  logical  :: tau_constant

  tau_constant = .true.

  stend = 0._r8
  qvtend = 0._r8
  qitend = 0._r8
  nitend = 0._r8

  ! calculate qsati from t,p,q

  call wv_sat_qsat_ice(t, p, ESI, QSI)

  if (naai.gt.1.e-18_r8.and.qv.gt.QSI) then

     !optional timescale on condensation
     !tau in sections. Try 300. or tau = f(T): 300s  t> 268, 1800s for t<238
     !     
     if (.not. tau_constant) then
        if( t.gt. 268.15_r8 ) then
           tau = 300.0_r8
        elseif(t.lt.238.15_r8 ) then
           tau = 1800._r8
        else
           tau = 300._r8 + (1800._r8 - 300._r8) * ( 268.15_r8 - t ) / 30._r8
        endif
     else
         tau = 300._r8
     end if

     qitend = (qv-QSI)/deltat !* exp(-tau/deltat)
     qvtend = 0._r8 - qitend
     stend  = qitend * xxls    ! moist static energy tend...[J/kg/s] !

     ! kg(h2o)/kg(air)/s * J/kg(h2o)  = J/kg(air)/s (=W/kg)
     ! if ice exists (more than 1 L-1) and there is condensation, do not add to number (= growth), else, add 10um ice

     if (ni.lt.1.e3_r8.and.(qi+qitend*deltat).gt.1e-18_r8) then
        nitend = nitend + 3._r8 * qitend/(4._r8*3.14_r8* 10.e-6_r8**3*997._r8)
     endif

  endif

end subroutine ice_macro_tend

end module macrop_driver
